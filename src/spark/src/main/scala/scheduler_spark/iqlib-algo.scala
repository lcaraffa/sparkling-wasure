package spark_ddt.ddt_algo

import sys.process._
import scala.io.Source
import java.io._
import scala.xml._
import java.lang.Double
import scala.concurrent._
import scala.collection.parallel._
import java.nio.file.{ Paths, Files }
import scala.collection.mutable.ListBuffer

import java.util.concurrent.Executors
import java.util.Date;
import java.text.SimpleDateFormat;
import java.util.Calendar;

import org.apache.spark._;
import org.apache.spark.graphx._;
import org.apache.spark.graphx.PartitionStrategy._
import org.apache.spark.rdd.RDD;
import org.apache.spark.sql.SaveMode
import org.apache.spark.storage.StorageLevel
import org.apache.spark.SparkConf
import org.apache.spark.HashPartitioner

import spark_ddt.core._;
import spark_ddt.util._;
import algo_stats._;
import spark_ddt.core.IQlibCore._;
import xml_parsing._;
import bash_funcs._
import strings_opt._;
import params_parser._;
import files_opt._;
import tiling._;
import geojson_export._;


object ddt_algo {

  def saveAsPly(rdd_ply : RDD[VData],output_dir : String,plot_lvl : Int = 3)
  {
    if(plot_lvl > 2)
      saveAsPly2(rdd_ply,output_dir)
    else
      saveAsPly2(rdd_ply,output_dir)
  }

  def saveAsPly1(rdd_ply : RDD[VData],output_dir : String)
  {
    rdd_ply.saveAsTextFile(output_dir)
  }

  def saveAsPly2(rdd_ply : RDD[VData],output_dir : String){
    val nbp = rdd_ply.count.toInt
    rdd_ply.zipWithIndex().map(
      x => (x._2,x._1)
    ).partitionBy(new HashPartitioner(nbp)).map(
      x => x._2.replace(';','\n')
    ).saveAsTextFile(output_dir)

  }

  def update_time(params_scala : params_map, label : String, ps : Char = '_')
  {
    val t0_init = params_scala("t0").head.toLong
    val t_curr = System.nanoTime();
    val time_string = ((t_curr - t0_init)/1000000000.0).toString
    val id = params_scala.filterKeys(x => x contains "time").size
    val ss = ("[time]"+ label).padTo(30, ps) + "[" + time_string + "]"+ps+ps;
    println(ss);
    params_scala("[time]" + "%03d".format(id) + "_" + label) = collection.mutable.Set(time_string)
  }

  def update_global_ids(kvrdd_finalized_tri : RDD[KValue],kvrdd_stats : RDD[KValue],rep_merge : Int,iq : IQlibSched , params_cpp : params_map, sc : SparkContext) : RDD[KValue] =
  {
    val stats_cum = kvrdd_simplex_id(kvrdd_stats,sc);
    val update_global_id_cmd =  set_params(params_cpp, List(("step","update_global_id"))).to_command_line
    val res_tri_gid = iq.run_pipe_fun_KValue(
      update_global_id_cmd,
      (kvrdd_finalized_tri union stats_cum).reduceByKey(_ ::: _,rep_merge), "ext_gr", do_dump = false).filter(!_.isEmpty())
    val kvrdd_gid_tri = iq.get_kvrdd(res_tri_gid)
    return kvrdd_gid_tri;
  }

  def compute_tiling(kvrdd_inputs : RDD[KValue],iq : IQlibSched,params_ddt : params_map,params_scala : params_map ) : RDD[KValue] =
  {
    val plot_lvl = params_scala("plot_lvl").head.toInt;
    val tile_cmd =  set_params(params_ddt, List(("step","tile_ply"))).to_command_line
    val dump_ply_binary_cmd =  set_params(params_ddt, List(("step","dump_ply_binary"))).to_command_line
    val tri2geojson_cmd =  set_params(params_ddt, List(("step","tri2geojson"))).to_command_line
    val max_ppt = params_scala("max_ppt").head.toInt;
    val ndtree_depth = params_scala("ndtree_depth").head.toInt;
    val output_dir = params_scala("output_dir").head;
    val dim = params_scala("dim").head.toInt
    val t0_init = params_scala("t0").head.toLong
    val nbt = new nd_tree(dim,10,iq)
    val id_pad = nbt.nb_nodes(ndtree_depth-1);

    kvrdd_inputs.count;

    val res_tiles = iq.run_pipe_fun_KValue(
      tile_cmd ++ List("--label", "pts_tile"),
      kvrdd_inputs, "tile", do_dump = false).setName("KVRDD_TILING");

    if(plot_lvl >= 3){
      if(dim == 3){
        iq.run_pipe_fun_KValue(
          dump_ply_binary_cmd ++ List("--label", "tile_pts"),
          iq.get_kvrdd(res_tiles,"z"), "tile_pts", do_dump = false).collect()
      }
    }

    res_tiles.persist(iq.get_storage_level())
    res_tiles.count
    update_time(params_scala,"tiling");
    val kvrdd_count: RDD[KValue] = iq.get_kvrdd(res_tiles, "c",txt="count").cache();
    val rddc = kvrdd_count.map(x => (x._1.toInt + id_pad,x.value(0).split(" ").last.toInt)).reduceByKey(_+_).cache()
    val max_tile = rddc.reduce((x,y) => if(x._2 > y._2) x else y);
    val tot_nbp = rddc.map(x => x._2.toLong).reduce(_ + _);
    println("total number of points : " + tot_nbp )
    params_scala("tot_nbp") =  collection.mutable.Set(tot_nbp.toString)
    if(max_tile._2 > max_ppt){
      println("======================  WARNING  =====================")
      println(" A tile contains more point than the maximum allowed  ")
      println("====> To solve it Increase the Octree depth  ")
      println("======================  WARNING  =====================")
    }


    println("Compute id map ...")
    val (res_v,lf_new) = nbt.compute_id_rdd(rddc,max_ppt);
    val rdd_idmap = res_v.reduceByKey((u,v) => u).cache()

    update_time(params_scala,"octree");
    val rdd_idmap_sorted = lf_new.map(x => (x._1,x._2._1)).sortBy(_._2,false).zipWithIndex().map( x=> (x._1._1.toLong,x._2))
    val idmap_sorted = rdd_idmap.map(x => (x._2,x._1)).cogroup(rdd_idmap_sorted).filter(
      x => ((!x._2._2.isEmpty) && (!x._2._1.isEmpty))).map(
      x => (x._1, x._2._2.head)
    ).collectAsMap


    println("===== rdd_idmap =====")
    println("Mapping keys ...")

    val nb_leaf = rdd_idmap.map(x => x._2).distinct.count
    params_scala("nb_leaf") =  collection.mutable.Set(nb_leaf.toString)
    println("Number of leafs :" + nb_leaf)
    val rep_value = nb_leaf.toInt;


    if(plot_lvl > 3){
      res_tiles.saveAsTextFile(output_dir + "/rdd_tiles")
    }

    val kvrdd_points : RDD[KValue] = iq.get_kvrdd(res_tiles,"z").map(x => ((x._1+  id_pad),x._2)).cogroup(rdd_idmap).filter(
      x => ((!x._2._2.isEmpty) && (!x._2._1.isEmpty))).map(
      x => (x._2._2.head.toLong,x._2._1.reduce(_ ::: _))).map(x => (idmap_sorted(x._1),x._2)).reduceByKey(new HashPartitioner(rep_value),(u,v) => u ::: v).setName("KVRDD_TILING");


    if(plot_lvl > 3){
      kvrdd_points.saveAsTextFile(output_dir + "/kvrdd_points")
    }

    println("Repartition done")
    kvrdd_points.persist(iq.get_storage_level()).setName("KVRDD_POINT_TILED");
    kvrdd_points.count();

    update_time(params_scala,"mapping");
    res_tiles.unpersist()

    return kvrdd_points
  }



  def algo_loop(cur_tri : RDD[KValue], cur_edges : RDD[TEdge] , iq : IQlibSched , params_cpp : params_map,params_scala : params_map) : (RDD[KValue],RDD[TEdge],ListBuffer[(Int,RDD[VData])]) = {


    val rep_loop = params_scala("rep_loop").head.toInt;
    val output_dir = params_scala("output_dir").head;
    val plot_lvl = params_scala("plot_lvl").head.toInt;
    val dim_algo = params_scala("dim").head
    val dump_mode = params_scala("dump_mode").head;

    // C++ command init
    val insert_in_triangulation_fun = set_params(params_cpp,List(("step","insert_in_triangulation"))).to_command_line
    val tri2geojson_cmd =  set_params(params_cpp, List(("step","tri2geojson"))).to_command_line

    // Loop variables
    val defaultV = (List(""));
    val kvrdd_ply_list : ListBuffer[(Int,RDD[VData])] = ListBuffer()

    var kvrdd_cur_tri = cur_tri;
    var rdd_cur_edges = cur_edges;
    var is_not_last = true
    var acc = 0
    var nb_edges = 0L;
    while (is_not_last) {
      val kvrdd_tri_with_edges = (rdd_cur_edges.map(e => (e.dstId,e.attr)) union kvrdd_cur_tri).reduceByKey(_ ::: _,rep_loop).persist(iq.get_storage_level_loop()).setName("NEW_KVRDD_WITH_EDGES_" + acc)
      kvrdd_tri_with_edges.count
      kvrdd_cur_tri.unpersist();
      rdd_cur_edges.unpersist();

      val kvrdd_insert_new_pts = kvrdd_tri_with_edges.filter(x => (x._2.size > 1 ))
      val kvrdd_insert_non_pts = kvrdd_tri_with_edges.subtractByKey(kvrdd_insert_new_pts)

      // Inset bbox points
      val res_tri_insert_nbr = iq.run_pipe_fun_KValue(
        insert_in_triangulation_fun ++ List("--extract_edg_nbrs","","--label", "tri_loop_" + acc.toString),
        kvrdd_insert_new_pts, "insert_and_neighbors_" + acc.toString, do_dump = false).persist(iq.get_storage_level_loop()).setName("NEW_KVRDD_WITH_EDGES_" + acc)
      kvrdd_cur_tri  = iq.get_kvrdd(res_tri_insert_nbr, "t",txt = "tri_loop" + acc.toString) union kvrdd_insert_non_pts;
      kvrdd_cur_tri.persist(iq.get_storage_level_loop()).setName("KVRDD_CUR_TRI_" + acc)
      kvrdd_cur_tri.count

      //Extract the edge of the merging graph
      rdd_cur_edges = iq.get_edgrdd(res_tri_insert_nbr, "e");
      rdd_cur_edges.persist(iq.get_storage_level_loop()).setName("KVRDD_CUR_EDGES_" + acc)

      nb_edges = rdd_cur_edges.count
      res_tri_insert_nbr.unpersist()

      println("NB edges  :" + nb_edges)
      kvrdd_tri_with_edges.unpersist()

      is_not_last =  nb_edges > 0
      update_time(params_scala,"loop" + acc.toString);
      params_scala("[edges]loop_" + acc.toString) = collection.mutable.Set(nb_edges.toString)
      acc = acc + 1;


      if(plot_lvl >=3 ||  (plot_lvl >=2 && !is_not_last)){
        if(dim_algo.toInt == 2){
          val rdd_json_tri = iq.run_pipe_fun_KValue(
            tri2geojson_cmd ++ List("--label", "tri_loop_" + acc.toString,"--style","tri_main.qml"),
            kvrdd_cur_tri, "tri_loop", do_dump = false)
          rdd_json_tri.collect()
        }
      }
    }

    return (kvrdd_cur_tri, rdd_cur_edges,kvrdd_ply_list)
  }



  def compute_ddt(kvrdd_points: RDD[KValue], iq : IQlibSched , params_cpp : params_map, params_scala : params_map): (TGraph,ListBuffer[RDD[KValue]],RDD[KValue])  = {

    val dim_algo = params_scala("dim").head
    val rep_loop = params_scala("rep_loop").head.toInt;
    val rep_merge = params_scala("rep_merge").head.toInt;
    val plot_lvl = params_scala("plot_lvl").head.toInt;
    val dump_mode = params_scala("dump_mode").head;
    val output_dir = params_scala("output_dir").head;

    // C++ function
    val insert_in_triangulation_fun = set_params(params_cpp,List(("step","insert_in_triangulation"))).to_command_line
    val ply2geojson_fun =  set_params(params_cpp, List(("step","ply2geojson"))).to_command_line
    val tri2geojson_cmd =  set_params(params_cpp, List(("step","tri2geojson"))).to_command_line

    // Algo variables
    val kvrdd_log_list : ListBuffer[RDD[KValue]] = ListBuffer()

    //==== Starts algo ==========
    kvrdd_points.count
    val res_tri_local = iq.run_pipe_fun_KValue(
      insert_in_triangulation_fun ++ List("--label", "tri_tile","--extract_tri_crown",""),
      kvrdd_points, "insert_in_triangulation", do_dump = false)

    res_tri_local.persist(iq.get_storage_level()).setName("KVRDD_RES_TRI_LOCAL")
    res_tri_local.count
    update_time(params_scala,"local");

    val kvrdd_pts_crown : RDD[KValue] = iq.get_kvrdd(res_tri_local, "z",txt = "tri_crown").setName("tri_crown")
    val kvrdd_pts_finalized : RDD[KValue] = iq.get_kvrdd(res_tri_local, "r",txt = "pts_finalized").setName("pts_finalized")
    kvrdd_pts_crown.coalesce(rep_loop).persist(iq.get_storage_level_loop())
    val rdd_finalized_ply_local : RDD[VData] = iq.filter_rdd(res_tri_local, "p",txt = "finalized").setName("tri_full")
    rdd_finalized_ply_local.persist(iq.get_storage_level())
    val kvrdd_bbox : RDD[KValue] = iq.get_kvrdd(res_tri_local, "q",txt = "tri_full").setName("bbox")
    val kvrdd_count : RDD[KValue] = iq.get_kvrdd(res_tri_local, "c",txt="stats");
    val ll_bbox_pts = kvrdd_bbox.map(x => x.value).reduce(_ ++ _)
    val llrdd_input_tri = kvrdd_pts_crown.map(x => (x.key, (x.value ++ ll_bbox_pts)))
    val res_tri_bbox = iq.run_pipe_fun_KValue(
      insert_in_triangulation_fun ++ List("--extract_edg_nbrs","","--label", "tri_bbox"),
      llrdd_input_tri, "insert_in_triangulation_bbox", do_dump = false).persist(iq.get_storage_level_loop())
    res_tri_bbox.count
    update_time(params_scala,"bbox");

    var kvrdd_cur_tri  : RDD[KValue] = iq.get_kvrdd(res_tri_bbox, "t",txt = "tri_bbox")
    var rdd_cur_edges = iq.get_edgrdd(res_tri_bbox, "e");
    val defaultV = (List(""));
    val (kvrdd_last_tri,kvrdd_last_edges,kvrdd_finalized_ply_list) = {
      algo_loop(kvrdd_cur_tri,rdd_cur_edges,iq, params_cpp,params_scala);
    }

    var kvrdd_finalized_tri = kvrdd_last_tri;
    var kvrdd_stats = kvrdd_last_tri;

    // Insert finalized point into global triangulation (non optimal)
    var kvrdd_b4_union =  kvrdd_last_tri;
    if(dump_mode == "NONE"){
      kvrdd_b4_union =  (kvrdd_last_tri union kvrdd_pts_finalized).reduceByKey((u,v) => u ::: v)
    }
    val res_finalized_tri = iq.run_pipe_fun_KValue(
      insert_in_triangulation_fun ++ List(
        "--finalize_tri","",
        "--send_empty_edges","",
        "--extract_edg_nbrs","",
        "--label", "ftri"),
      kvrdd_b4_union , "merge", do_dump = false);
    val rdd_finalized_ply_crown : RDD[VData] = iq.filter_rdd(res_finalized_tri, "p",txt = "finalized").setName("tri_full")
    kvrdd_finalized_tri = iq.get_kvrdd(res_finalized_tri, "t",txt = "merged_tri")
    kvrdd_stats = iq.get_kvrdd(res_finalized_tri, "s",txt="stats").cache()
    val kvrdd_finalized_edges =  iq.get_edgrdd(res_finalized_tri, "e");
    kvrdd_stats.count();
    update_time(params_scala,"ddtdone");

    // Ploting section
    if(plot_lvl >=2 && dim_algo.toInt == 2){
      val rdd_json_merged_tri = iq.run_pipe_fun_KValue(
        tri2geojson_cmd ++ List("--label", "res_merged_tri","--style","tri_main.qml"),
        kvrdd_finalized_tri, "res_merged_tri", do_dump = false)
      rdd_json_merged_tri.collect()
      update_time(params_scala,"plotingdone");
    }

    // Dumping section
    if(dump_mode != "NONE"){
      saveAsPly(rdd_finalized_ply_crown,output_dir + "/rdd_ply_finalized_last",plot_lvl)
      kvrdd_last_tri.unpersist();
      update_time(params_scala,"shareddumped");

      saveAsPly(rdd_finalized_ply_local,output_dir + "/rdd_ply_finalized_init",plot_lvl)
      update_time(params_scala,"localdumped");
      update_time(params_scala,"W:Dumpingdone");
    }

    val final_graph = Graph(kvrdd_finalized_tri, kvrdd_finalized_edges, defaultV)
    final_graph.vertices.setName("graph_final_tri");
    final_graph.edges.setName("graph_");
    return (final_graph,kvrdd_log_list,kvrdd_stats)
  }


  def compute_ddt_nograph(kvrdd_points: RDD[KValue], iq : IQlibSched , params_cpp : params_map, params_scala : params_map): (RDD[KValue],RDD[TEdge],ListBuffer[RDD[KValue]],RDD[KValue])  = {

    val dim_algo = params_scala("dim").head
    val rep_loop = params_scala("rep_loop").head.toInt;
    val rep_merge = params_scala("rep_merge").head.toInt;
    val plot_lvl = params_scala("plot_lvl").head.toInt;
    val dump_mode = params_scala("dump_mode").head;
    val output_dir = params_scala("output_dir").head;

    // C++ function
    val insert_in_triangulation_fun = set_params(params_cpp,List(("step","insert_in_triangulation"))).to_command_line
    val ply2geojson_fun =  set_params(params_cpp, List(("step","ply2geojson"))).to_command_line
    val tri2geojson_cmd =  set_params(params_cpp, List(("step","tri2geojson"))).to_command_line


    // Algo variables
    // Olg log system
    val kvrdd_log_list : ListBuffer[RDD[KValue]] = ListBuffer()


    //==== Starts algo ==========
    kvrdd_points.count
    val res_tri_local = iq.run_pipe_fun_KValue(
      insert_in_triangulation_fun ++ List("--label", "tri_tile","--extract_tri_crown",""),
      kvrdd_points, "insert_in_triangulation", do_dump = false)


    res_tri_local.persist(iq.get_storage_level()).setName("KVRDD_RES_TRI_LOCAL")
    res_tri_local.count
    update_time(params_scala,"local");

    val kvrdd_pts_crown : RDD[KValue] = iq.get_kvrdd(res_tri_local, "z",txt = "tri_crown").setName("tri_crown")
    val kvrdd_pts_finalized : RDD[KValue] = iq.get_kvrdd(res_tri_local, "r",txt = "pts_finalized").setName("pts_finalized")
    kvrdd_pts_crown.coalesce(rep_loop).persist(iq.get_storage_level_loop())
    val rdd_finalized_ply_local : RDD[VData] = iq.filter_rdd(res_tri_local, "p",txt = "finalized").setName("tri_full")
    rdd_finalized_ply_local.persist(iq.get_storage_level())
    val kvrdd_bbox : RDD[KValue] = iq.get_kvrdd(res_tri_local, "q",txt = "tri_full").setName("bbox")
    val kvrdd_count : RDD[KValue] = iq.get_kvrdd(res_tri_local, "c",txt="stats");


    val ll_bbox_pts = kvrdd_bbox.map(x => x.value).reduce(_ ++ _)
    val llrdd_input_tri = kvrdd_pts_crown.map(x => (x.key, (x.value ++ ll_bbox_pts)))
    val res_tri_bbox = iq.run_pipe_fun_KValue(
      insert_in_triangulation_fun ++ List("--extract_edg_nbrs","","--label", "tri_bbox"),
      llrdd_input_tri, "insert_in_triangulation_bbox", do_dump = false).persist(iq.get_storage_level_loop())
    res_tri_bbox.count
    update_time(params_scala,"bbox");


    // ============= DO LOOP ==============
    var kvrdd_cur_tri  : RDD[KValue] = iq.get_kvrdd(res_tri_bbox, "t",txt = "tri_bbox")
    var rdd_cur_edges = iq.get_edgrdd(res_tri_bbox, "e");
    val defaultV = (List(""));
    val (kvrdd_last_tri,kvrdd_last_edges,kvrdd_finalized_ply_list) = {
      algo_loop(kvrdd_cur_tri,rdd_cur_edges,iq, params_cpp,params_scala);
    }

    var kvrdd_finalized_tri = kvrdd_last_tri;
    var kvrdd_stats = kvrdd_last_tri;

    // Insert finalized point into global triangulation (non optimal)
    var kvrdd_b4_union =  kvrdd_last_tri;
    if(dump_mode == "NONE"){
      kvrdd_b4_union =  (kvrdd_last_tri union kvrdd_pts_finalized).reduceByKey((u,v) => u ::: v)
    }
    val res_finalized_tri = iq.run_pipe_fun_KValue(
      insert_in_triangulation_fun ++ List(
        "--finalize_tri","",
        "--send_empty_edges","",
        "--extract_edg_nbrs","",
        "--label", "ftri"),
      kvrdd_b4_union , "merge", do_dump = false);
    val rdd_finalized_ply_crown : RDD[VData] = iq.filter_rdd(res_finalized_tri, "p",txt = "finalized").setName("tri_full")
    kvrdd_finalized_tri = iq.get_kvrdd(res_finalized_tri, "t",txt = "merged_tri")
    kvrdd_stats = iq.get_kvrdd(res_finalized_tri, "s",txt="stats").cache()
    val kvrdd_finalized_edges =  iq.get_edgrdd(res_finalized_tri, "e");
    kvrdd_stats.count();

    update_time(params_scala,"ddtdone");

    // Ploting section
    if(plot_lvl >=2 && dim_algo.toInt == 2){
      val rdd_json_merged_tri = iq.run_pipe_fun_KValue(
        tri2geojson_cmd ++ List("--label", "res_merged_tri","--style","tri_main.qml"),
        kvrdd_finalized_tri, "res_merged_tri", do_dump = false)
      rdd_json_merged_tri.collect()
      update_time(params_scala,"plotingdone");
    }

    // Dumping section
    if(dump_mode != "NONE"){
      saveAsPly(rdd_finalized_ply_crown,output_dir + "/rdd_ply_finalized_last",plot_lvl)
      kvrdd_last_tri.unpersist();
      update_time(params_scala,"shareddumped");

      saveAsPly(rdd_finalized_ply_local,output_dir + "/rdd_ply_finalized_init",plot_lvl)
      update_time(params_scala,"localdumped");
      update_time(params_scala,"W:Dumpingdone");
    }

    return (kvrdd_finalized_tri, kvrdd_finalized_edges,kvrdd_log_list,kvrdd_stats)
  }


  // Ok this is not optimal but it's working
  def extract_2D_voronoi(graph_tri : TGraph,  iq : IQlibSched,params_cpp : params_map,params_scala : params_map){

    val update_global_id_cmd =  set_params(params_cpp, List(("step","update_global_id"))).to_command_line
    val extract_graph_local_cmd =  set_params(params_cpp, List(("step","extract_voronoi"),("area_processed","1"))).to_command_line
    val extract_graph_shared_cmd =  set_params(params_cpp, List(("step","extract_voronoi"),("area_processed","2"))).to_command_line
    val kvrdd_gid_tri = graph_tri.vertices //iq.get_kvrdd(res_tri_gid)

    // Build the full graph
    val graph_full = graph_tri
    val input_vertex : RDD[KValue] =  graph_full.vertices
    val input_edges : RDD[KValue] =  graph_full.convertToCanonicalEdges().triplets.map(ee => (ee.srcId,ee.srcAttr ++ ee.dstAttr))

    // Full graph local
    val full_graph_local = iq.run_pipe_fun_KValue(
      extract_graph_local_cmd,
      input_vertex, "ext_gr", do_dump = false).filter(!_.isEmpty())

    // Full graph shared
    val full_graph_shared = iq.run_pipe_fun_KValue(
      extract_graph_shared_cmd,
      input_edges, "ext_gr", do_dump = false).filter(!_.isEmpty())

    val tri_vertex = full_graph_local.filter(x => x(0) == 'v').map(
      x => x.split(" ")).map(x => x.splitAt(2)).map(cc => (cc._1(1).toLong, cc._2.map(_.toDouble)))
    val tri_simplex = full_graph_local.filter(x => x(0) == 's').filter(x => x.count(_ == 's') == 1).map(
      x => x.split(" ")).map(x => x.splitAt(2)).map(cc => (cc._1(1).toLong, cc._2.map(_.toDouble)))
    val tri_edges = (full_graph_local.filter(x => x(0) == 'e') union full_graph_shared.filter(x => x(0) == 'e')).map(
      x => x.split(" ")).map(cc => Edge(cc(1).toLong, cc(2).toLong,""))

    val g_voronoi = Graph(tri_simplex,tri_edges)

    val geojson_header= raw"""{
    "type": "FeatureCollection",
    "features": [
    """
    val geojson_footer = raw"""
     ]
     }
  """

    val geojson_fh = raw"""
     {
      "type": "Feature",
      "geometry": {
        "type": "LineString",
        "coordinates": [
  """
    val geojson_fb =   raw"""
        ]
      }
    }
  """
    val vr_str = geojson_header + g_voronoi.triplets.filter(tt => tt.dstAttr != null &&  tt.srcAttr != null).map(
      tt => geojson_fh + "[" + tt.srcAttr(0) + "," + tt.srcAttr(1) + "],[" + tt.dstAttr(0) + "," + tt.dstAttr(1) + "]" + geojson_fb ).reduce(_ + ","+_) + geojson_footer

    val bname =  params_scala("output_dir").head +"/voronoi";
    val pw1 = new PrintWriter(new File(bname + ".geojson" ))
    pw1.write(vr_str)
    pw1.close

  }
}

