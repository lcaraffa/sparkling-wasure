/*
 * (c) Copyright 2016 Hewlett Packard Enterprise Development LP
 *
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package sparkle.graph

import org.apache.spark.{ SparkConf, SparkContext}
import org.apache.spark.graphx.{Graph, TripletFields}

/**
  * Implementation of Loopy Belief Propagation algorithm for factor graphs.
  * For details see: https://en.wikipedia.org/wiki/Belief_propagation
  * Uses libDAI factor graph file format as an input. In addition, each factor
  * must have a unique id specified by ###ID in the file.
  * For details see: https://staff.fnwi.uva.nl/j.m.mooij/libDAI/doc/fileformats.html
  *
  */
object BP {

  /**
    * Returns graph after running the BP algorithm. Values are stored in log scale.
    * @param graph initial factor graph
    * @param maxIterations maximum iterations
    * @param eps epsilon
    * @return graph
    */
  def apply(
  graph: Graph[FGVertex, Unit],
  maxIterations: Int = 50,
  eps: Double = 1e-3): Graph[FGVertex, FGEdge] = {
    // put initial messages on edges, they will be mutated every iteration
    var newGraph = graph.mapTriplets { triplet =>
      new FGEdge(triplet.srcAttr.initMessage(triplet.dstAttr.id),
        triplet.dstAttr.initMessage(triplet.srcAttr.id), false)
    }.cache()
    val numEdges = newGraph.edges.count
    var oldGraph = newGraph
    // main algorithm loop:
    var iter = 0
    var converged = false
    while (iter < maxIterations && !converged) {
      oldGraph = newGraph
      // messages to variables are merged as a product, messages to factors are merged as lists
      val newAggMessages = newGraph.aggregateMessages[List[Message]](
        triplet => {
          triplet.sendToDst(List(triplet.attr.toDst))
          triplet.sendToSrc(List(triplet.attr.toSrc))
        },
        // TODO: extract Merge into the new AggregatedMessage class and use mutable structures
        (m1, m2) => {
          if (m1(0).fromFactor && m2(0).fromFactor) {
            List(Message(m1(0).srcId, m1(0).message.compose(m2(0).message), fromFactor = true))
          } else {
            m1 ++ m2
          }
        },
        TripletFields.EdgeOnly
      )
      newGraph = newGraph.joinVertices(newAggMessages)(
        (id, attr, msg) => attr.processMessage(msg))
        .mapTriplets { triplet =>
          val toSrc = triplet.dstAttr.sendMessage(triplet.attr.toDst)
          val toDst = triplet.srcAttr.sendMessage(triplet.attr.toSrc)
          val diffSrc = toSrc.message.maxDiff(triplet.attr.toSrc.message)
          val diffDst = toDst.message.maxDiff(triplet.attr.toDst.message)
          FGEdge(toDst, toSrc, diffSrc < eps && diffDst < eps)}
        .cache()
      if (iter == 0) {
        newGraph.edges.foreachPartition(x => {})
        converged = false
      } else {
        val numConverged = newGraph.edges.aggregate(0L)((res, edge) =>
          if (edge.attr.converged) res + 1 else res, (res1, res2) =>
            res1 + res2)
        converged = numConverged == numEdges
      }
      oldGraph.unpersist(false)
      iter += 1
    }
    // TODO: return beliefs as RDD[Beliefs] that can be computed at the end as message product
    newGraph
  }

  /**
    * Main function to run BP from the command line.
    * Required arguments: [path to input] [path to output] [iterations] [epsilon]
    * Optional argument (local) runs the algorithm in the Spark local mode
    * @param args
    */
  def main(args: Array[String]): Unit = {
    if (args.length < 4) {
      throw new IllegalArgumentException("Insufficient arguments")
    }
    val conf = if (args.length == 5 && args(4) == "local") {
      new SparkConf().setAppName("Belief Propagation Application").setMaster("local")
    } else {
      new SparkConf().setAppName("Belief Propagation Application")
    }
    val sc = new SparkContext(conf)
    val inputPath = args(0)
    val outputPath = args(1)
    val numIter = args(2).toInt
    val epsilon = args(3).toDouble
    val graph = Utils.loadLibDAIToFactorGraph(sc, inputPath)
    val time = System.nanoTime()
    val beliefs = BP(graph, maxIterations = numIter, eps = epsilon)
    // TODO: output to a file instead
    val calculatedProbabilities = beliefs.vertices.flatMap { case(id, vertex) => vertex match {
      case n: NamedVariable => Seq((n.id, n.belief))
      case _ => Seq.empty[(Long, Variable)]
      }
    }
    calculatedProbabilities.saveAsTextFile(outputPath)
  }
}
