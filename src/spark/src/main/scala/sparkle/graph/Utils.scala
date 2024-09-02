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

import java.io._

import org.apache.spark.graphx.{Edge, Graph}
import org.apache.spark.SparkContext

import scala.collection.mutable.ArrayBuffer

/**
  * Utils contains functions for factor graph loading
  */
object Utils  {

  /**
    * Load factor graph from libDAI format
    * @param sc SparkContext
    * @param path path to a file or folder with files
    * @return factor graph
    */
  def loadLibDAIToFactorGraph(sc: SparkContext, path: String): Graph[FGVertex, Unit] = {
    val partitions = sc.binaryFiles(path).count()
    val x = sc.binaryFiles(path, partitions.toInt).map { case (_, stream) =>
      val dataStream = stream.open()
      val data = loadLibDAI(dataStream)
      dataStream.close()
      data
    }
    // TODO: refactor y => y, type of edge
    val (vertices, edges) = (x.keys.flatMap(y => y).map(x => (x.id, x)),
      x.values.flatMap(y => y).map(x => Edge(x._1, x._2, ())))
    val graph = Graph(vertices, edges)
    graph
  }

  /**
    * Returns arrays of vertices and edges loaded from libDAI file
    * @param fileName graph file name
    * @return vertices and edges
    */
  def loadLibDAI(fileName: String): (Array[FGVertex], Array[(Long, Long)]) = {
    val inputStream = new FileInputStream(fileName)
    loadLibDAI(inputStream)
  }

  /**
    * Returns arrays of vertices and edges loaded from the stream of libDAI format
    * @param stream libDAI stream
    * @return vertices and edges
    */
  def loadLibDAI(stream: InputStream): (Array[FGVertex], Array[(Long, Long)]) = {
    val reader = new BufferedReader(new InputStreamReader(stream))
    // read the number of factors in the file
    var line = dropWhile(reader, x => x.length < 1 || x.startsWith("#"))
    val numFactors = line.trim.toInt
    // TODO: check that this structure is OK and size estimation is fair enough
    val factorBuffer = new ArrayBuffer[FGVertex](numFactors * 2)
    val edgeBuffer = new ArrayBuffer[(Long, Long)](numFactors * 10)
    // read factors
    var factorCounter = 0
    // factor ids are in source format as comment e.g. ###222
    while (factorCounter < numFactors) {
      // skip to the next block with factors
      line = dropWhile(reader, l => !l.startsWith("###"))
      val factorId = line.split(" ")(1).trim.toLong
      line = dropWhile(reader, l => l.startsWith("#"))
      val varNum = line.trim.toInt
      line = dropWhile(reader, l => l.startsWith("#"))
      val varIds = line.split("\\s+").map(_.toLong)
      line = dropWhile(reader, l => l.startsWith("#"))
      val varNumValues = line.split("\\s+").map(_.toInt)
      line = dropWhile(reader, l => l.startsWith("#"))
      val nonZeroNum = line.toInt
      val indexAndValues = new Array[(Int, Double)](nonZeroNum)
      var nonZeroCounter = 0
      while (nonZeroCounter < nonZeroNum) {
        line = dropWhile(reader, l => l.startsWith("#"))
        val indexAndValue = line.split("\\s+")
        // Log conversion
        var value = indexAndValue(1).toDouble
        if (value == 0 ) value = Double.MinPositiveValue
        indexAndValues(nonZeroCounter) = (indexAndValue(0).toInt, math.log(value))
        nonZeroCounter += 1
      }
      // create Factor vertex
      val namedFactor = NamedFactor(factorId, varIds, varNumValues, nonZeroNum, indexAndValues)
      // create Variable vertex if factor has only one variable and add factor there as a prior
      if (varNum == 1) {
        val namedVariable = new NamedVariable(varIds(0),
          belief = Variable(namedFactor.factor.cloneValues),
          prior = Variable(namedFactor.factor.cloneValues))
        factorBuffer += namedVariable
      } else {
        factorBuffer += namedFactor
        // create edges between Variables and Factor
        for (varId <- varIds) {
          edgeBuffer += new Tuple2(varId, factorId)
        }
      }
      // increment of factor counter
      factorCounter += 1
    }
    reader.close()
    (factorBuffer.toArray, edgeBuffer.toArray)
  }

  private def dropWhile(reader: BufferedReader, condition: String => Boolean): String  = {
    var line = reader.readLine()
    if (line == null) {
      throw new IOException("End of file reached")
    }
    while (condition(line)) {
      line = reader.readLine()
    }
    line
  }
}
