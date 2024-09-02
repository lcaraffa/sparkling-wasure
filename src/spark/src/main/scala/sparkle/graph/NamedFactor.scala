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
  * Factor Graph trait. Contains methods for initialization, sending and processing messages
  */
private [graph] trait FGVertex extends Serializable {
  /**
    * id of the vertex
    */
  val id: Long

  /**
    * Returns a new vertex that results from processing a message
    * @param aggMessage list of messages
    * @return vertex
    */
  def processMessage(aggMessage: List[Message]): FGVertex

  /**
    * Returns a message based on the incoming message by decomposition of the belief
    * @param incomingMessage incoming message
    * @return message
    */
  def sendMessage(incomingMessage: Message): Message

  /**
    * Returns the initial message
    * @param varId id of a variable-recipient of the message
    * @return message
    */
  def initMessage(varId: Long): Message
}

/**
  * Representation of a named factor.
  * Named means that it has a mapping between the variable ids and the values.
  * Named factor also contains factor belief.
  * @param id factor id
  * @param variables list of variable ids
  * @param factor factor
  * @param belief belief
  */
private [graph] class NamedFactor(
 val id: Long,
 val variables: Array[Long],
 val factor: Factor, val belief: Factor) extends FGVertex {

  /**
    * Returns variable index in the values array by its ID
    * @param varId variable ID
    * @return index
    */
  private def varIndexById(varId: Long): Int = {
    val index = variables.indexOf(varId)
    require(index >= 0, "Variable not found in factor")
    index
  }

  /**
    * Number of values for a variable
    * @param varId variable id
    * @return number of values
    */
  def length(varId: Long): Int = {
    val index = varIndexById(varId)
    factor.length(index)
  }

  override def processMessage(aggMessage: List[Message]): FGVertex = {
    assert(aggMessage.length > 1)
    // TODO: optimize (use one factor)
    var newBelief = factor
    for(message <- aggMessage) {
      val index = varIndexById(message.srcId)
      newBelief = newBelief.compose(message.message, index)
    }
    NamedFactor(id, variables, factor, newBelief)
  }

  override def sendMessage(oldMessage: Message): Message = {
    val index = varIndexById(oldMessage.srcId)
    val newMessage = belief.decompose(oldMessage.message, index)
    Message(this.id, Variable(newMessage.cloneValues), fromFactor = true)
  }

  override def initMessage(varId: Long): Message = {
    val length = this.length(varId)
    val value = 0.0
    Message(this.id, Variable.fill(length)(value), fromFactor = true)
  }
}

/**
  * Fabric for NamedFactors
  */
object NamedFactor {

  /**
    * Return new NamedFactor
    * @param id id
    * @param variables list of variables ids
    * @param factor factor
    * @param belief belief factor
    * @return named factor
    */
  def apply(id: Long, variables: Array[Long], factor: Factor, belief: Factor): NamedFactor = {
    new NamedFactor(id, variables, factor, belief)
  }

  /**
    * Create factor from the libDAI description
    * @param id unique id
    * @param variables ids of variables
    * @param states num of variables states
    * @param nonZeroNum num of nonzero values
    * @param indexAndValues index and values
    * @return factor
    */
  def apply(
             id: Long,
             variables: Array[Long],
             states: Array[Int],
             nonZeroNum: Int,
             indexAndValues: Array[(Int, Double)]): NamedFactor = {
    val values = new Array[Double](states.product)
    var i = 0
    while (i < indexAndValues.size) {
      val (index, value) = indexAndValues(i)
      values(index) = value
      i += 1
    }
    new NamedFactor(id, variables, Factor(states, values), belief = Factor(states.clone(), values.clone()))
  }
}

/**
  * Representation of a named variable instance.
  * Contains variable id, belief and prior
  * @param id id
  * @param belief belief
  * @param prior prior
  */
private [graph] case class NamedVariable (
  val id: Long,
  val belief: Variable,
  val prior: Variable) extends FGVertex {
  override def processMessage(aggMessage: List[Message]): FGVertex = {
    val newBelief = prior.compose(aggMessage(0).message)
    NamedVariable(id, newBelief, prior)
  }

  override def sendMessage(oldMessage: Message): Message = {
    val newMessage = belief.decompose(oldMessage.message)
    Message(this.id, newMessage, fromFactor = false)
  }

  override def initMessage(varId: Long): Message = {
    val length = this.belief.size
    val value = 0.0
    Message(this.id, Variable.fill(length)(value), fromFactor = false)
  }
}

// TODO: find a better name for message field
/**
  * Representation of a message
  * @param srcId source id
  * @param message message
  * @param fromFactor true if from a factor (to a variable)
  */
private [graph] case class Message(val srcId: Long, val message: Variable, val fromFactor: Boolean)

/**
  * Representation of an edge in Belief Propagation graph
  * @param toDst message from source to destination
  * @param toSrc message from destination to source
  * @param converged true if edge converged
  */
private [graph] case class FGEdge(val toDst: Message, val toSrc: Message, val converged: Boolean)
