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

import org.apache.spark.rdd.RDD
import org.apache.spark.{SparkConf, SparkContext}
import org.apache.spark.graphx.{Edge, Graph, TripletFields}
import org.apache.spark.storage.StorageLevel

object PairwiseBP2   {

  /**
    * Runs Pairwise Loopy Belief Propagation
    * @param graph graph: vertexes contain prior and belief, edges contain factors
    * @param maxIterations max iterations
    * @param eps epslion
    * @return graph after running the algorithm
    */
  def apply(graph: Graph[PVertex, PEdge], maxIterations: Int = 50, eps: Double = 1e-3, slvl : StorageLevel = StorageLevel.MEMORY_AND_DISK) : Graph[PVertex, PEdge] = {
    var newGraph = graph.cache()
    var prevGraph: Graph[PVertex, PEdge] = null
    var graphWithNewMessages: Graph[PVertex, PEdge] = null
    var iter = 0
    var converged = false
    val numVertices = newGraph.vertices.count
    var cur_e = newGraph.edges
    var cur_v = newGraph.vertices
    var old_e = cur_e;
    var old_v = cur_v;

    while (iter < maxIterations && !converged) {
      val t0 = System.nanoTime()

      newGraph = Graph(cur_v,cur_e)
      newGraph.vertices.setName("NG_BP__" + iter);
      newGraph.edges.setName("NG_BP__" + iter);
      prevGraph = newGraph
      graphWithNewMessages = newGraph.mapTriplets{ edge =>
        val tmpForwardMessage = edge.srcAttr.belief.decompose(edge.attr.reverseMessage)
        val tmpReverseMessage = edge.dstAttr.belief.decompose(edge.attr.forwardMessage)
        val factor = edge.attr.factor
        val newForwardMessage = factor.compose(tmpForwardMessage, 0).marginalize(1)
        val newReverseMessage = factor.compose(tmpReverseMessage, 1).marginalize(0)
        PEdge(factor, newForwardMessage, newReverseMessage)
      }.cache()

      graphWithNewMessages.vertices.setName("NewMess_BP__" + iter);
      graphWithNewMessages.edges.setName("NewMess_BP__" + iter);
      graphWithNewMessages.edges.foreachPartition( x => {})
      val newAggMessages = graphWithNewMessages.aggregateMessages[Variable](
        triplet => {
          triplet.sendToDst(triplet.attr.forwardMessage)
          triplet.sendToSrc(triplet.attr.reverseMessage)
        },
        (v1: Variable, v2: Variable) => v1.compose(v2),
        TripletFields.EdgeOnly)
      newGraph = graphWithNewMessages.joinVertices(newAggMessages) { (id, vertex, msg) =>
        val newBelief = vertex.prior.compose(msg)
        val maxDiff = newBelief.maxDiff(vertex.belief)
        PVertex(newBelief, vertex.prior, maxDiff < eps)
      }.cache()


      newGraph.edges.foreachPartition(x => {})


      val numConverged = newGraph.vertices.aggregate(0L)((res, vertex) =>
          if (vertex._2.converged) res + 1 else res, (res1, res2) =>
        res1 + res2)

      if(iter % 20 == 0 && iter !=0 && false){
        cur_v = newGraph.vertices
        cur_e = newGraph.edges
        cur_v.checkpoint();
        cur_e.checkpoint();        
      }else{
        cur_v = newGraph.vertices.unpersist().persist(slvl).setName("RDD_BP__" + iter)
        cur_e = newGraph.edges.unpersist().persist(slvl).setName("RDD_BP_" + iter)
        cur_v.count
        cur_e.count
        if(iter > 0){
          old_v.unpersist()
          old_e.unpersist()
        }
      }



      old_v = cur_v;
      old_e = cur_e;

      converged = numConverged == numVertices
      prevGraph.vertices.unpersist(false)
      prevGraph.edges.unpersist(false)
      graphWithNewMessages.edges.unpersist(false)
      val t1 = System.nanoTime()
      println("Belief iter:" + iter + " in " + ((t1 - t0)/1000000000.0).toString + "s")
      println("nc: " + numConverged + " nv:" + numVertices)
      iter += 1
    }
    newGraph
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
    val maxIterations = args(3).toInt
    val eps = args(4).toDouble
    val sc = new SparkContext(conf)
    val graph = loadPairwiseGraph(sc, args(0), args(1))
    val output = args(2)
    val time = System.nanoTime()
    val res = PairwiseBP(graph, maxIterations, eps)
    val calculatedProbabilities = res.vertices.mapValues(x => x.belief)
    calculatedProbabilities.saveAsTextFile(output)
  }

  /**
    * Loads pairwise factor graph from vertex and edge files
    * @param sc spark context
    * @param vertexFile vertex file: each line contains variable id and prior
    * @param edgeFile edge file: each line contains edge srcId, dstId and 2d factor in column-major format
    * @return graph
    */
  def loadPairwiseGraph(sc: SparkContext, vertexFile: String, edgeFile: String): Graph[PVertex, PEdge] = {
    val vertices: RDD[String] = sc.textFile(vertexFile)
    val edges: RDD[String] = sc.textFile(edgeFile)
    val vertexRDD: RDD[(Long, PVertex)] = vertices.map { line =>
      val fields = line.split(' ')
      val id = fields(0).toLong
      val prior = Variable(fields.tail.map(x => math.log(x.toDouble)))
      val belief = Variable(prior.cloneValues)
      (id, PVertex(belief, prior))
    }
    val edgeRDD = edges.map { line =>
      val fields = line.split(' ')
      val srcId = fields(0).toLong
      val dstId = fields(1).toLong
      val factor = fields.tail.tail.map(x => math.log(x.toDouble))
      Edge(srcId, dstId, factor)
    }
    val tmpGraph = Graph(vertexRDD, edgeRDD)
    // fixing the dimensions of factor based on vertexes
    val graph = tmpGraph.mapTriplets { triplet =>
      val srcSize = triplet.srcAttr.prior.size
      val dstSize = triplet.dstAttr.prior.size
      PEdge(Factor(Array(srcSize, dstSize), triplet.attr), Variable.fill(srcSize)(0.0), Variable.fill(dstSize)(0.0))
    }
    graph
  }
}


object PairwiseBP   {

  /**
    * Runs Pairwise Loopy Belief Propagation
    * @param graph graph: vertexes contain prior and belief, edges contain factors
    * @param maxIterations max iterations
    * @param eps epslion
    * @return graph after running the algorithm
    */
  def apply(graph: Graph[PVertex, PEdge], maxIterations: Int = 50, eps: Double = 1e-3) : Graph[PVertex, PEdge] = {
    var newGraph = graph.cache()
    var prevGraph: Graph[PVertex, PEdge] = null
    var graphWithNewMessages: Graph[PVertex, PEdge] = null
    var iter = 0
    var converged = false
    val numVertices = newGraph.vertices.count
    while (iter < maxIterations && !converged) {
      prevGraph = newGraph
      graphWithNewMessages = newGraph.mapTriplets{ edge =>
        val tmpForwardMessage = edge.srcAttr.belief.decompose(edge.attr.reverseMessage)
        val tmpReverseMessage = edge.dstAttr.belief.decompose(edge.attr.forwardMessage)
        val factor = edge.attr.factor
        val newForwardMessage = factor.compose(tmpForwardMessage, 0).marginalize(1)
        val newReverseMessage = factor.compose(tmpReverseMessage, 1).marginalize(0)
        PEdge(factor, newForwardMessage, newReverseMessage)
      }.cache()
      graphWithNewMessages.edges.foreachPartition( x => {})
      val newAggMessages = graphWithNewMessages.aggregateMessages[Variable](
        triplet => {
          triplet.sendToDst(triplet.attr.forwardMessage)
          triplet.sendToSrc(triplet.attr.reverseMessage)
        },
        (v1: Variable, v2: Variable) => v1.compose(v2),
        TripletFields.EdgeOnly)
      newGraph = graphWithNewMessages.joinVertices(newAggMessages) { (id, vertex, msg) =>
        val newBelief = vertex.prior.compose(msg)
        val maxDiff = newBelief.maxDiff(vertex.belief)
        PVertex(newBelief, vertex.prior, maxDiff < eps)
      }.cache()
      newGraph.edges.foreachPartition(x => {})
      val numConverged = newGraph.vertices.aggregate(0L)((res, vertex) =>
          if (vertex._2.converged) res + 1 else res, (res1, res2) =>
          res1 + res2)
//      logInfo("%d/%d vertices converged".format(numConverged, numVertices))
      converged = numConverged == numVertices
      prevGraph.vertices.unpersist(false)
      prevGraph.edges.unpersist(false)
      graphWithNewMessages.edges.unpersist(false)
      iter += 1
    }
//    logInfo("Total %d/%d iterations completed. Inference %s with epsilon = %f".
      //format(iter, maxIterations, if (converged) "converged" else "did not converge", eps))
    newGraph
  }
}
/**
  * Edge in pairwise BP. Represents a factor
  * @param factor factor table: rows represent A states, columns - B states
  * @param forwardMessage message from A to B according to the factor table, rows are the states
  * @param reverseMessage message from B to A (reverse in respect of the factor table)
  */
case class PEdge(factor: Factor, forwardMessage: Variable, reverseMessage: Variable)

/**
  * Vertex in pairwise BP. Represents a variable
  * @param belief belief
  * @param prior prior
  * @param converged true if converged
  */
case class PVertex(belief: Variable, prior: Variable, converged: Boolean = true)
