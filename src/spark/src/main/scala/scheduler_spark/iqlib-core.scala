package spark_ddt.core

import sys.process._
import org.apache.spark._;
import org.apache.spark.graphx._;
import scala.io.Codec


import org.apache.spark.rdd.RDD;
import org.apache.spark.SparkContext._
import org.apache.spark.SparkContext
import org.apache.spark.SparkConf
import org.apache.spark.storage.StorageLevel


import org.apache.hadoop.conf.Configuration
import org.apache.hadoop.fs.FileSystem
import org.apache.hadoop.fs.Path



object IQlibCore {

  type VData = String;
  type Value = List[VData];

  type KValue = (VertexId, Value);
  type TEdge = Edge[Value];

  type TGraph = Graph[Value, Value];
  type Triplet = EdgeTriplet[Value, Value];

  implicit class StrFun(ss: VData) {
  }
  implicit class Vdat(ss: Value) {
    def value = ss;
  }

  implicit class KValuefun(tt: KValue) {
    def key = tt._1
    def value = tt._2
    def data = tt._2(0)
    def content = tt._2(0).split(" ", 5).last
  }
  

  implicit class TEdgefun(ee: TEdge) {
    def toKValue: KValue = {
      (ee.srcId, ee.attr)
    }
  }

  def get_fpath_value(vdat : VData) : Char = {
     var cc : Char = ' '
     var acc = 1;
     while (! { cc=vdat.charAt(acc); cc }.isLetter){
     acc = acc  +1
    }
   return cc
  }

  def is_data_file(vdat: VData) : Boolean = {
    return (get_fpath_value(vdat) == 'f')
  }

}


class IQlibSched(
  slvl_glob: StorageLevel = StorageLevel.NONE,
  slvl_loop: StorageLevel = StorageLevel.NONE,
  env_map : scala.collection.immutable.Map[String,String] = Map("" -> "")
) extends Serializable {

  import IQlibCore._



  def get_storage_level() : StorageLevel = {
    return slvl_glob;
  }

  def get_storage_level_loop() : StorageLevel = {
    return slvl_loop;

  }


  def time[R](block: => R): R = {
    val t0 = System.nanoTime()
    val result = block // call-by-name
    val t1 = System.nanoTime()
    println("Elapsed time: " + (t1 - t0) / 1000000000.0 + "s")
    result
  }

  // =============== Schedul stuff ========================
  def print_rdde_vdat(record: KValue, f: String => Unit) = {
    f(record._1 + " " + record._2.size + " ")
    for (e <- record._2) { f(e + " ") }
  }

  def run_pipe_fun_KValue(fun: Seq[String], rdd: RDD[KValue], txt: String = "rdd", do_dump: Boolean = false): RDD[VData] = {
    val inname = "/tmp/rdd_" + txt + "_input";
    val outname = "/tmp/rdd_" + txt + "_output";

    println("========= [Pipe function] ==========");
    val tmp = fun.map(x => print(x + " "));
    println("")

    if (do_dump) {
      println("input saved in :" + inname);
      rdd.collect();
      rdd.saveAsTextFile(inname);
    }
    val rdd2 = rdd.pipe(
      command = fun,
      env = env_map,
      printRDDElement = print_rdde_vdat);


    println("==================================");

    rdd2.setName("RDD_PIPED_" + txt)
    rdd2;
  }

  def rdd_TEdge2KVLdat(rdd: RDD[TEdge]): RDD[KValue] = {
    rdd.map(ee => (ee.srcId, ee.attr));
  }

  def run_pipe_fun_Edge(fun: Seq[VData], rdd: RDD[TEdge], txt: String = "rdd", do_dump: Boolean = false): RDD[VData] = {
    val rdd2 = rdd_TEdge2KVLdat(rdd);
    run_pipe_fun_KValue(fun, rdd2, txt, do_dump)
  }

  def run_pipe_fun_Fullriplet(fun: Seq[VData], rdd: RDD[Triplet], txt: String = "rdd", do_dump: Boolean = false): RDD[VData] = {
    val rdd2: RDD[KValue] = rdd.map(ee => (ee.srcId, ee.srcAttr ::: ee.dstAttr ::: ee.attr));
    run_pipe_fun_KValue(fun, rdd2, txt, do_dump)
  }

  def run_pipe_fun_VertsTiplet(fun: Seq[VData], rdd: RDD[Triplet], txt: String = "rdd", do_dump: Boolean = false): RDD[VData] = {
    val rdd2: RDD[KValue] = rdd.map(ee => (ee.srcId, ee.srcAttr ::: ee.dstAttr));
    run_pipe_fun_KValue(fun, rdd2, txt, do_dump)
  }

  // ============== Tiling Stuff =============

  def split_key_value2_oldgood(ss: VData): KValue = {
    val rr = (ss.split(" ", 4))
    (rr(2).toLong, List(ss))
  }

  def split_key_value2(ss: VData): KValue = {
    var rr = "";
    var pad = 0;
    if(ss.charAt(0) == ' ')
      pad = 1
    if(ss.length > 15)
      rr = (ss.charAt(4).toString + ss.charAt(5) + ss.charAt(6) + ss.charAt(7) + ss.charAt(8)
        + ss.charAt(9) + ss.charAt(10) + ss.charAt(11) + ss.charAt(12) + ss.charAt(13) + ss.charAt(14) + ss.charAt(15)).split(" ")(pad)
    else
      rr = (ss.split(" ", 4))(2+pad)
    (rr.toLong, List(ss))
  }


  def split_key_value3(ss: VData): TEdge = {
    val rr = (ss.split(" ", 5))
    new TEdge(rr(2).toLong, rr(3).toLong, List(ss))
  }

  def filter_rdd (rdd: RDD[VData], ss: String = "", txt: String = ""): RDD[VData] = {
    var res :RDD[VData] = null;
    if (!ss.isEmpty())
      res = rdd.filter(_(0) == ss(0))
    else
      res = rdd
      return res.setName("GetKVRDD_" + txt);
   }


  def get_kvrdd(rdd: RDD[VData], ss: String = "", txt: String = ""): RDD[KValue] = {
    var res :RDD[KValue] = null;
    if (!ss.isEmpty())
      res = rdd.filter(_(0) == ss(0)).map(split_key_value2(_))
    else
      res = rdd.map(split_key_value2(_));

      return res.setName("GetKVRDD_" + txt);
   }

  def get_edgrdd(rdd: RDD[VData], ss: String = ""): RDD[TEdge] = {
    if (!ss.isEmpty())
      rdd.filter(_(0) == ss(0)).map(split_key_value3(_))
    else
      rdd.map(split_key_value3(_))
  }

  // ================ Graph stuff ====================
  def aggregate_edge(graph: TGraph): VertexRDD[TEdge] = {
    
    def msgFunTvertIni(triplet: EdgeContext[Value, Value, TEdge]) {
      val l1 = Edge(triplet.srcId, triplet.dstId, triplet.attr);
      triplet.sendToDst(l1);

    }
    def reduceFunTvertIni(a: TEdge, b: TEdge): TEdge = new TEdge(b.srcId, b.dstId, b.attr ::: a.attr)
    val res = graph.aggregateMessages[TEdge](msgFunTvertIni, reduceFunTvertIni)
    res
  }

  def aggregate_edge2kvrdd(graph: TGraph): RDD[KValue] = {
    aggregate_edge(graph).mapValues(ee => ee.attr);
  }

  def aggregate_value_clique(graph: TGraph, deg: Int,rep_value : Int): RDD[KValue] = {

    def msgFunTvertIni(triplet: EdgeContext[Value, Value, Value]) {
      triplet.sendToDst(triplet.srcAttr)

    }
    def reduceFunTvertIni(a: Value, b: Value): Value = a ::: b

    val res1: VertexRDD[Value] = graph.aggregateMessages[Value](msgFunTvertIni, reduceFunTvertIni)
    deg match {
      case 0 => graph.vertices.setName("aggregated")
      case _ => graph.vertices.union(res1).reduceByKey(_ ::: _, rep_value).setName("aggregated")
    }
  }


  def aggregate_value_clique(graph: TGraph, deg: Int): RDD[KValue] = {

    def msgFunTvertIni(triplet: EdgeContext[Value, Value, Value]) {
      triplet.sendToDst(triplet.srcAttr)

    }
    def reduceFunTvertIni(a: Value, b: Value): Value = a ::: b

    val res1: VertexRDD[Value] = graph.aggregateMessages[Value](msgFunTvertIni, reduceFunTvertIni)
    deg match {
      case 0 => graph.vertices.setName("aggregated")
      case _ => graph.vertices.union(res1).reduceByKey(_ ::: _).setName("aggregated")
    }
  }




  def aggregate_edge_clique(graph: TGraph): RDD[KValue] = {

    def msgFunEdge1(triplet: EdgeContext[Value, Value, Value]) {
      val toSrc = triplet.attr;
      val toDst = triplet.attr;
      triplet.sendToDst(toDst);
      triplet.sendToSrc(toSrc);
    }

    def reduceFunEdge1(a: Value, b: Value): Value = a ::: b;
    graph.aggregateMessages[Value](msgFunEdge1, reduceFunEdge1)

    def msgFunEdge2(
      triplet: EdgeContext[Value, Value, Value]) {
      val toSrc = triplet.dstAttr;
      val toDst = triplet.srcAttr;
      triplet.sendToDst(toDst)
      triplet.sendToSrc(toSrc)
    }

    val result = graph.aggregateMessages[Value](msgFunEdge1, reduceFunEdge1);
    val defaultV = List("");
    val graph2 = Graph(result, graph.edges, defaultV);
    graph2.aggregateMessages[Value](msgFunEdge2, reduceFunEdge1);

  }

  def update_vertices(graph: TGraph, rrr: RDD[KValue]): VertexRDD[Value] = {
    graph.vertices.aggregateUsingIndex(rrr, (vd1: Value, vd2: Value) => vd2)
  }

  def get_focal_kvrdd(graph: TGraph, deg: Int): RDD[KValue] = {
    val result_tri = aggregate_value_clique(graph, deg)
    val result_edges = aggregate_edge_clique(graph).mapValues(_.distinct)
    deg match {
      case 0 => result_tri
      case _ => (result_tri union result_edges).reduceByKey(_ ::: _);
    }
  }

  def get_focal_rdd_pts(graph: TGraph, deg: Int): RDD[KValue] = {
    val result_tri = aggregate_value_clique(graph, deg);
    result_tri
  }

  def get_edge_rdd(graph: TGraph): RDD[TEdge] = {
    graph.edges
  }

  // ========== Depercated for the moment but usefull ==============
  def aggregate_edge_clique_fat(graph: TGraph): RDD[(VertexId, List[TEdge])] = {

    def exists_in_ring(rid: VertexId, vid: VertexId, one_ring: Array[(VertexId, Array[VertexId])]): Boolean = {
      one_ring.exists(x =>
        (x._1 == rid && x._2.exists(y =>
          y == vid))) || rid == vid;
    }
    def msgFunEdge1(triplet: EdgeContext[Value, Value, List[TEdge]]) {
      val l1 = List(Edge(triplet.srcId, triplet.dstId, triplet.attr));
      val l2 = List(Edge(triplet.srcId, triplet.dstId, triplet.attr));
      triplet.sendToDst(l1);
      triplet.sendToSrc(l2);
    }

    def msgFunEdge2(
      triplet: EdgeContext[List[TEdge], Value, List[TEdge]],
      one_ring: Array[(VertexId, Array[VertexId])]) {
      val res1 = triplet.srcAttr.filter(e =>
        (exists_in_ring(triplet.dstId, e.srcId, one_ring) && exists_in_ring(triplet.dstId, e.dstId, one_ring)))
      val res2 = triplet.dstAttr.filter(e =>
        (exists_in_ring(triplet.srcId, e.srcId, one_ring) && exists_in_ring(triplet.srcId, e.dstId, one_ring)))
      triplet.sendToDst(res1)
      triplet.sendToSrc(res2)
    }

    def reduceFunEdge(a: List[TEdge], b: List[TEdge]): List[TEdge] = a ::: b

    def rm_dup(ledge: List[TEdge], one_ring: Array[(VertexId, Array[VertexId])]): List[TEdge] = {
      ledge.groupBy(x => (x.srcId, x.dstId)).map(_._2.head).toList
    }

    val one_ring = graph.collectNeighborIds(EdgeDirection.Either).collect();
    val result = graph.aggregateMessages[List[TEdge]](msgFunEdge1, reduceFunEdge);
    val defaultEdge2 = List(Edge(0L, 0L, List("void")));
    val graph2 = Graph(result, graph.edges, defaultEdge2);
    val res2 = graph2.aggregateMessages[List[TEdge]](msgFunEdge2(_, one_ring), reduceFunEdge);
    res2.mapValues(rm_dup(_, one_ring))
  }

  def send_edge_attr_to_dst(graph: TGraph): RDD[KValue] = {
    def msgFunTvertIni(triplet: EdgeContext[Value, Value, Value]) {
      triplet.sendToDst(triplet.attr)
    }
    def reduceFunTvertIni(a: Value, b: Value): Value = a ::: b

    val res1: VertexRDD[Value] = graph.aggregateMessages[Value](msgFunTvertIni, reduceFunTvertIni)
    graph.vertices.union(res1).reduceByKey(_ ::: _)
  }
}

