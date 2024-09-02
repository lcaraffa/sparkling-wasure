/*
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

import java.util.Calendar

import org.apache.spark.rdd.RDD
import org.apache.spark.{SparkConf, SparkContext}
import org.apache.spark.graphx.{EdgeTriplet, Graph, _}

object TwoStateBP   {

  def apply(graph: Graph[BeliefVertex, BeliefEdge],
   maxIterations: Int = 50,
   eps: Double = 1e-3) : Graph[BeliefVertex, BeliefEdge] =
  {
    val initialMessage = Array(1.0, 1.0)
    var maxDiff: Double = 1.0
    var newGraph = graph.mapVertices((vid, vdata) => vertexProgram(0, vid, vdata, initialMessage)).cache()
    var prevGraph: Graph[BeliefVertex, BeliefEdge] = null
    var graphWithNewMessages: Graph[BeliefVertex, BeliefEdge] = null
    var i = 1
    while (i <= maxIterations && maxDiff >= eps) {
      prevGraph = newGraph
      graphWithNewMessages = newGraph.mapTriplets(putMessages).cache()
      graphWithNewMessages.edges.foreachPartition( x => {})
      val newAggMessages = graphWithNewMessages.aggregateMessages[Array[Double]](
        triplet => {
          triplet.sendToDst(triplet.attr.forwardMessage)
          triplet.sendToSrc(triplet.attr.reverseMessage)
        },
        mergeMessages,
        TripletFields.EdgeOnly)
      newGraph = graphWithNewMessages.joinVertices(newAggMessages) { (id, attr, msg) =>
        vertexProgram(i, id, attr, msg) }.cache()
      // TODO: check that edges force vertices to materialize
      newGraph.edges.foreachPartition( x => {})
      prevGraph.vertices.unpersist(false)
      prevGraph.edges.unpersist(false)
      graphWithNewMessages.edges.unpersist(false)
      maxDiff = newGraph.vertices.map { case (vid, vb) => vb.diff }.max()
      i += 1
    }
    newGraph
  } // end of apply

  private def putMessages(edge: EdgeTriplet[BeliefVertex, BeliefEdge]):
  BeliefEdge = {
    // copy factor from the edge to the message
    val factor = edge.attr.factor
    val newForwardMessage = Array.ofDim[Double](2)
    val newReverseMessage = Array.ofDim[Double](2)
    val srcBelief = edge.srcAttr.belief
    val dstBelief = edge.dstAttr.belief
    val forwardMessage = Array.ofDim[Double](2)
    val reverseMessage = Array.ofDim[Double](2)
    val maxForward = math.max(edge.attr.forwardMessage(0), edge.attr.forwardMessage(1))
    val maxReverse = math.max(edge.attr.reverseMessage(0), edge.attr.reverseMessage(1))
    forwardMessage(0) = math.exp(edge.attr.forwardMessage(0) - maxForward)
    reverseMessage(0) = math.exp(edge.attr.reverseMessage(0) - maxReverse)
    forwardMessage(1) = math.exp(edge.attr.forwardMessage(1) - maxForward)
    reverseMessage(1) = math.exp(edge.attr.reverseMessage(1) - maxReverse)
    newForwardMessage(0) = math.log(
      (factor(0)(0) * srcBelief(0) / reverseMessage(0) +
        factor(1)(0) * srcBelief(1) / reverseMessage(1)
        ) / (
        factor(0)(0) * srcBelief(0) / reverseMessage(0) +
          factor(1)(0) * srcBelief(1) / reverseMessage(1) +
          factor(0)(1) * srcBelief(0) / reverseMessage(0) +
          factor(1)(1) * srcBelief(1) / reverseMessage(1)
        )
    )
    newReverseMessage(0) = math.log(
      (factor(0)(0) * dstBelief(0) / forwardMessage(0) +
        factor(0)(1) * dstBelief(1) / forwardMessage(1)
        ) / (
        factor(0)(0) * dstBelief(0) / forwardMessage(0) +
          factor(0)(1) * dstBelief(1) / forwardMessage(1) +
          factor(1)(0) * dstBelief(0) / forwardMessage(0) +
          factor(1)(1) * dstBelief(1) / forwardMessage(1)
        )
    )
    newForwardMessage(1) = math.log(
      (factor(0)(1) * srcBelief(0) / reverseMessage(0) +
        factor(1)(1) * srcBelief(1) / reverseMessage(1)
        ) / (
        factor(0)(0) * srcBelief(0) / reverseMessage(0) +
          factor(1)(0) * srcBelief(1) / reverseMessage(1) +
          factor(0)(1) * srcBelief(0) / reverseMessage(0) +
          factor(1)(1) * srcBelief(1) / reverseMessage(1)
        )
    )
    newReverseMessage(1) = math.log(
      (factor(1)(0) * dstBelief(0) / forwardMessage(0) +
        factor(1)(1) * dstBelief(1) / forwardMessage(1)
        ) / (
        factor(0)(0) * dstBelief(0) / forwardMessage(0) +
          factor(0)(1) * dstBelief(1) / forwardMessage(1) +
          factor(1)(0) * dstBelief(0) / forwardMessage(0) +
          factor(1)(1) * dstBelief(1) / forwardMessage(1)
        )
    )
    BeliefEdge(factor, newForwardMessage, newReverseMessage)
  }

  private def vertexProgram(niter: Int, id: VertexId,
                    bVertex: BeliefVertex,
                    msgProd: Array[Double]): BeliefVertex = {
    // TODO: refactor computations
    val oldBelief1 = bVertex.belief(0)
    val oldBelief2 = bVertex.belief(1)
    val tmpBelief1 = math.log(bVertex.prior(0)) + msgProd(0)
    val tmpBelief2 = math.log(bVertex.prior(1)) + msgProd(1)
    val maxBelief = math.max(tmpBelief1, tmpBelief2)
    val tmpBelief1_1 = math.exp(tmpBelief1 - maxBelief)
    val tmpBelief1_2 = math.exp(tmpBelief2 - maxBelief)
    val sum = tmpBelief1_1 + tmpBelief1_2
    val newBelief1 = tmpBelief1_1 / sum
    val newBelief2 = tmpBelief1_2 / sum
    val maxDiff = math.max(math.abs(newBelief1 - oldBelief1), math.abs(newBelief2 - oldBelief2))
    BeliefVertex(bVertex.prior, Array(newBelief1, newBelief2), maxDiff)
  }

  def mergeMessages(a: Array[Double], b: Array[Double]): Array[Double] = {
    val res = new Array[Double](2)
    res(0) = a(0) + b(0)
    res(1) = a(1) + b(1)
    res
  }


  def main(args: Array[String]): Unit = {
    if (args.length < 5) {
      throw new IllegalArgumentException("Insufficient arguments")
    }
    val conf = if (args.length == 6 && args(5) == "local") {
      new SparkConf().setAppName("Belief Propagation Application").setMaster("local")
    } else {
      new SparkConf().setAppName("Belief Propagation Application")
    }
    val sc = new SparkContext(conf)
    val vertices: RDD[String] = sc.textFile(args(0))
    val edges: RDD[String] = sc.textFile(args(1))
    val output = args(2)
    val maxIterations = args(3).toInt
    val eps = args(4).toDouble
    val vertexRDD = vertices.map { line =>
      val fields = line.split(' ')
      (fields(0).toLong,
        BeliefVertex(Array(fields(1).toDouble, fields(2).toDouble), Array(fields(1).toDouble, fields(2).toDouble), 0.0)
      )
    }
    val edgeRDD = edges.map { line =>
      val fields = line.split(' ')
      Edge(fields(0).toLong, fields(1).toLong,
        BeliefEdge(
          Array(Array(fields(2).toDouble, fields(4).toDouble), Array(fields(3).toDouble, fields(5).toDouble)),
          Array(1.0, 1.0), Array(1.0, 1.0))
      )
    }
    val graph = Graph(vertexRDD, edgeRDD)
    val time = System.nanoTime()
    val res = TwoStateBP(graph, maxIterations, eps)
    res.vertices.mapValues(x => x.belief.mkString(" ")).saveAsTextFile(output)
  }

}
/**
 *
 * @param prior prior probabilities of each state (default - 2 states)
 * @param belief current probabilities of the state (belief)
 * @param diff difference with the belief on the previous step
 */
case class BeliefVertex(prior: Array[Double], belief: Array[Double], diff: Double)

/**
 *
 * @param factor factor table: rows represent A states, columns - B states
 * @param forwardMessage message from A to B according to the factor table, rows are the states
 * @param reverseMessage message from B to A (reverse in respect of the factor table)
 */
case class BeliefEdge(factor: Array[Array[Double]], forwardMessage: Array[Double], reverseMessage: Array[Double])

