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

/**
  * Contains math for factors.
  * Math is done in the log domain to prevent (postpone) underflow in the original domain.
  * Multiplication and division are addition and subtraction in log domain.
  */
private object FactorMath {
  // precision of computations in the original domain. serves as a safe-guard for log overflow
  val precision = math.log(Double.MinPositiveValue)
  def log(x: Double): Double = math.log(x)
  def log1p(x: Double): Double = math.log1p(x)
  def exp(x: Double): Double = math.exp(x)
  def composeLog(x: Double, y: Double): Double = {
    x + y
  }
  def decomposeLog(source: Double, y: Double): Double = {
    source - y
  }
  def logNormalize(x: Array[Double]): Unit = {
    var i = 0
    var sumExp = 0.0
    val max = x.max
    while (i < x.length) {
      sumExp += exp(x(i) - max)
      i += 1
    }
    sumExp = log(sumExp)
    i = 0
    while (i < x.length) {
      x(i) = x(i) - max - sumExp
      // if (x(i) < precision) x(i) = precision
      i += 1
    }
  }
  def logSum(x: Double, y: Double): Double = {
    val (a, b) = if (x >= y) (x, y) else (y, x)
    a + log1p(exp(b - a))
  }
  def logDiff(x: Double, y: Double): Double = {
    val (a, b) = if (x >= y) (x, y) else (y, x)
    a + log1p(-exp(b - a))
  }
}

/**
 *
 * Representation of a factor
 * Example how to loop through all variables
 *  for (i <- 0 until values.length) {
 *      var product: Int = 1
 *     for (dim <- 0 until varNumValues.length - 1) {
 *     val dimValue = (i / product) % varNumValues(dim)
 *        product *= varNumValues(dim)
 *        print(dimValue)
 *      }
 *      println(i / product)
 *    }
 *
 * Accessing a value by index
 * private def value(indices: Seq[Int]): Double = {
 *   // NB: leftmost index changes the fastest
 *   // NB: Wikipedia column-major order
 *   var offset = indices.last
 *   for (i <- states.length - 1 to 1 by -1) {
 *     offset = indices(i - 1) + states(i - 1) * offset
 *   }
 *   values(offset)
 * }
 *
 * @param states number of states
 * @param values values in vector format
 */
private [graph] class Factor private (
  protected val states: Array[Int],
  protected val values: Array[Double]) extends Serializable {

  /**
   * Total length of the factor in vector representation
   */
  val length: Int = values.length

  /**
   * Returns the number of states of a variable at index
   * @param index index
   * @return number of states
   */
  def length(index: Int): Int = {
    states(index)
  }

  /**
    * Returns the composition of factor and a message.
    * @param message message
    * @param index index of the variable-receiver of a message
    * @return factor
    */
  def compose(message: Variable, index: Int): Factor = {
    require(index >= 0, "Index must be non-negative")
    require(states(index) == message.size,
      "Number of states for variable and message must be equal")
    val result = new Array[Double](length)
    val product = states.slice(0, index).product
    var i = 0
    while (i < values.length) {
      val indexInTargetState = (i / product) % states(index)
      result(i) = FactorMath.composeLog(values(i), message.state(indexInTargetState))
      i += 1
    }
    Factor(states, result)
  }

  /**
    * Returns the normalized decomposition of factor and a message
    * @param message message
    * @param index index of the variable-receiver of a message
    * @return message
    */
  def decompose(message: Variable, index: Int): Variable = {
    require(index >= 0, "Index must be non-negative")
    require(states(index) == message.size,
      "Number of states for variable and message must be equal")
    val result = Array.fill[Double](states(index))(Double.NegativeInfinity)
    val product = states.slice(0, index).product
    var i = 0
    while (i < values.length) {
      val indexInTargetState = (i / product) % states(index)
      val decomposed = FactorMath.decomposeLog(values(i), message.state(indexInTargetState))
      val logSum = FactorMath.logSum(result(indexInTargetState), decomposed)
      result(indexInTargetState) = logSum
      i += 1
    }
    FactorMath.logNormalize(result)
    Variable(result)
  }

  /**
    * Marginalize factor by a variable
    *
    * @param index index of a variable
    * @return marginal
    */
  def marginalize(index: Int): Variable = {
    require(index >= 0 && index < states.length, "Index must be non-negative && within shape")
    val result = Array.fill[Double](states(index))(Double.NegativeInfinity)
    val product = states.slice(0, index).product
    var i = 0
    while (i < values.length) {
      val indexInTargetState = (i / product) % states(index)
      result(indexInTargetState) = FactorMath.logSum(result(indexInTargetState), values(i))
      i += 1
    }
    FactorMath.logNormalize(result)
    Variable(result)
  }

  /**
   * Clone values
   * @return values
   */
  def cloneValues: Array[Double] = {
    values.clone()
  }
}

/**
 * Fabric for factor creation
 */
object Factor {

  /**
    * Returns new factor given states and values
    * @param states states
    * @param values values
    * @return factor
    */
  def apply(states: Array[Int], values: Array[Double]): Factor = {
    new Factor(states, values)
  }
}


/**
 * Represenation of an instance of a variable
 * @param values values of a variable
 */
private [graph] class Variable private (
  protected val values: Array[Double]) extends Serializable {

  /**
    * Size of the variable (or the number of states)
    */
  val size = values.length

  /**
    * Returns the state of the variable at the given index
    * @param index index
    * @return state
    */
  def state(index: Int): Double = values(index)

  /**
    * Returns normalized composition of two variable instances
    * @param other other variable
    * @return variable
    */
  def compose(other: Variable): Variable = {
    require(this.size == other.size)
    val result = new Array[Double](size)
    var i = 0
    while (i < size) {
      result(i) = FactorMath.composeLog(this.values(i), other.values(i))
      i += 1
    }
    FactorMath.logNormalize(result)
    new Variable(result)
  }

  /**
    * Returns normalized decomposition of two variable instances
    * @param other
    * @return
    */
  def decompose(other: Variable): Variable = {
    require(this.size == other.size)
    val result = new Array[Double](size)
    var i = 0
    while (i < size) {
      result(i) = FactorMath.decomposeLog(this.values(i), other.values(i))
      i += 1
    }
    FactorMath.logNormalize(result)
    new Variable(result)
  }

  /**
    * Returns text representation of a variable instance
    * @return text
    */
  override def toString(): String = {
    // TODO: return 0 instead of MinPositiveDouble
    exp().values.mkString(" ")
  }

  /**
    * Returns variable instance in the original (non-log) domain
    * @return variable instance
    */
  def exp(): Variable = {
    // TODO: find a better place for this function
    val x = values.clone()
    var i = 0
    while (i < x.length) {
      x(i) = FactorMath.exp(x(i))
      i += 1
    }
    new Variable(x)
  }

  /**
    * Returns the maximum difference between two variables in the original domain
    * @param other other variable
    * @return maximum difference
    */
  def maxDiff(other: Variable): Double = {
    require(other.size == this.size, "Variables must have the same size")
    var i = 0
    var diff = 0.0
    while (i < values.length) {
      val d = FactorMath.exp(FactorMath.logDiff(values(i), other.values(i)))
      diff = if (d > diff) d else diff
      i += 1
    }
    diff
  }

  /**
   * Clone values
   * @return values
   */
  def cloneValues: Array[Double] = {
    values.clone()
  }
}

/**
  * Fabric for variables
  */
object Variable {

  /**
    * Returns a new variable given the values
    * @param values values
    * @return variable
    */
  def apply(values: Array[Double]): Variable = {
    new Variable(values)
  }

  /**
    * Returns a new variable given the values to fill
    * @param size size
    * @param elem fill values
    * @return variable
    */
  def fill(size: Int)(elem: => Double): Variable = {
    new Variable(Array.fill[Double](size)(elem))
  }
}
