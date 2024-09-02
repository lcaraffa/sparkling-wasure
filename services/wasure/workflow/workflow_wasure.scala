import sys.process._

import scala.xml._
import scala.concurrent._
import scala.collection.parallel._
import scala.collection.mutable.ListBuffer
import collection.mutable

import java.io._
import java.lang.Double
import java.util.concurrent.Executors
import java.util.Date;
import java.util.Calendar;
import java.text.SimpleDateFormat;
import java.text.DecimalFormat
import java.io.PrintWriter
import java.nio.file.{Paths, Files}
import java.nio.charset.StandardCharsets

import org.apache.hadoop.fs.FileSystem
import org.apache.hadoop.fs.Path
import org.apache.hadoop.fs.permission.FsPermission

import org.apache.spark._
import org.apache.spark.graphx._
import org.apache.spark.graphx.PartitionStrategy._
import org.apache.spark.rdd.RDD;
import org.apache.spark.sql.SaveMode
import org.apache.spark.storage.StorageLevel
import org.apache.spark.SparkConf

import spark_ddt.core._
import spark_ddt.util._
import spark_ddt.core.IQlibCore._
import spark_ddt.ddt_algo._
import spark_ddt.wasure._

import sparkle.graph._
import tiling._
import algo_stats._
import dataset_processing._
import bash_funcs._
import geojson_export._
import xml_parsing._
import strings_opt._
import params_parser._
import files_opt._


//=============================================
//==== Configuration and file sysyem init  ====
val conf = new SparkConf().setAppName("DDT")
val fs = FileSystem.get(sc.hadoopConfiguration);

// Checkpoint
val do_checkpoint = false
val checkpoint_dir_string = "/home/laurent/shared_spark/checkpoint/"
val checkpoint_dir_path = new Path(checkpoint_dir_string)
if(do_checkpoint){
  sc.setCheckpointDir(checkpoint_dir_string)
  if (fs.exists(checkpoint_dir_path))
    fs.delete(checkpoint_dir_path, true)
  fs.mkdirs(checkpoint_dir_path,new FsPermission("777"))
}

// Metadata extraction
val output_dir = get_bash_variable("OUTPUT_DATA_DIR").replaceAll("//", "/");
val input_dir = get_bash_variable("INPUT_DATA_DIR").replaceAll("//", "/");
val env_xml = get_bash_variable("PARAM_PATH");
val ddt_main_dir = get_bash_variable("DDT_MAIN_DIR");
val current_plateform = get_bash_variable("CURRENT_PLATEFORM");
val global_build_dir = get_bash_variable("GLOBAL_BUILD_DIR");
val executor_cores = get_bash_variable("MULTIVAC_EXECUTOR_CORE");
val num_executors = get_bash_variable("MULTIVAC_NUM_EXECUTORS");


// Get Params list from xml
val xml_string = sc.wholeTextFiles(env_xml).collect()(0)._2
var param_list = parse_xml_datasets_string(xml_string)
val df_par = sc.defaultParallelism;
val params_scala = param_list(0) // We only process 1 set of parameter in this workflow

// ===============================================
// ==== Scala and param initialization ===========
// Param scala is mutable, get params set the default value to the collection if it's empty
//  Se for instance the xml documentation / Algorithm params for the effect

// System params
val dim = params_scala.get_param("dim", "3").toInt
val ddt_kernel_dir = params_scala.get_param("ddt_kernel", "build-spark-Release-" + dim.toString)
val build_dir = global_build_dir + "/" + ddt_kernel_dir
val slvl_glob = StorageLevel.fromString(params_scala.get_param("StorageLevel", "DISK_ONLY"))
val slvl_loop = StorageLevel.fromString(params_scala.get_param("StorageLevelLoop", "DISK_ONLY"))

// General Algo params
val bbox = params_scala.get_param("bbox", "")
val do_profile = params_scala.get_param("do_profile", "false").toBoolean;
val do_stats = params_scala.get_param("do_stats", "false").toBoolean;
val plot_lvl = params_scala.get_param("plot_lvl", "1").toInt;
val regexp_filter = params_scala.get_param("regexp_filter", "");
val max_ppt = params_scala.get_param("max_ppt", "10000").toInt
val do_dump_debug = params_scala.get_param("dump_debug", "false").toBoolean
val ndtree_depth = params_scala.get_param("ndtree_depth", "4").toInt
val nbp =  params_scala.get_param("nbp", "10000").toInt
val datatype =  params_scala.get_param("datatype", "")
val spark_core_max = params_scala.get_param("spark_core_max", df_par.toString).toInt
val algo_seed =  params_scala.get_param("algo_seed",scala.util.Random.nextInt(100000).toString);

// Surface reconstruction prarams
val wasure_mode = params_scala.get_param("mode", "-1")
val pscale = params_scala.get_param("pscale", "0.05").toFloat
val nb_samples = params_scala.get_param("nb_samples", "3").toFloat
val rat_ray_sample = params_scala.get_param("rat_ray_sample", "1").toFloat
val min_ppt = params_scala.get_param("min_ppt", "100").toInt
val dst_scale = params_scala.get_param("dst_scale", "-1").toFloat
val lambda = params_scala.get_param("lambda", "4").toFloat
val coef_mult = params_scala.get_param("coef_mult", "1").toFloat
var max_opt_it = params_scala.get_param("max_opt_it", "30").toInt

val main_algo_opt = params_scala.get_param("algo_opt", "seg_lagrange_weight")
val stats_mod_it = params_scala.get_param("stats_mod_it", ((max_opt_it)/2).toString).toInt

params_scala("dump_mode") = collection.mutable.Set("NONE")
val fmt = new java.text.DecimalFormat("##0.##############")
val dateFormatter = new SimpleDateFormat("dd_MM_yyyy_hh_mm_ss")
val datestring = dateFormatter.format(Calendar.getInstance().getTime());
val cur_output_dir ={output_dir  + "outputs"} //sc.applicationId + "_" + datestring + "_"+ params_scala("name").head }
fs.mkdirs(new Path(cur_output_dir),new FsPermission("777"))

var env_map: scala.collection.immutable.Map[String, String] = Map("" -> "")

val cpp_exec_path = current_plateform.toLowerCase match {
  case _  => build_dir + "/bin/"
}

val python_exec_path = "python3 /nothing/to/run/in/python/"
val iq = new IQlibSched(slvl_glob,slvl_loop,env_map)


// Set the c++ command line object
val params_new = new Hash_StringSeq with mutable.MultiMap[String, String]
val params_ddt =  set_params(params_new,List(
  ("exec_path", cpp_exec_path + "ddt-stream-exe"),
  ("dim",params_scala("dim").head),
  ("bbox",params_scala("bbox").head),
  ("ech_input","1"),
  ("input_dir",input_dir),
  ("output_dir",cur_output_dir),
  ("min_ppt",params_scala("min_ppt").head),
  ("dump_mode",params_scala("dump_mode").head),
  ("seed",algo_seed)
))


val params_python =  set_params(params_new,List(
  ("exec_path", python_exec_path + "python_stream.py")
))

val params_wasure =  set_params(params_new,List(
  ("exec_path", cpp_exec_path + "wasure-stream-exe"),
  ("dim",params_scala("dim").head),
  ("bbox",params_scala("bbox").head),
  ("lambda",params_scala("lambda").head),
  ("dst_scale",params_scala("dst_scale").head),
  ("pscale",params_scala("pscale").head),
  ("rat_ray_sample",params_scala("rat_ray_sample").head),
  ("nb_samples",params_scala("nb_samples").head),
  ("mode",params_scala("mode").head),
  ("input_dir",input_dir),
  ("output_dir",cur_output_dir),
  ("seed",algo_seed)
))

if(do_dump_debug)
  params_wasure("dump_debug") = collection.mutable.Set("")

val floatFormat = new DecimalFormat("#.###")
val nbt_side = Math.min(math.pow(2,ndtree_depth),256).toInt
val tot_nbt = scala.math.pow(nbt_side,dim).toInt;
val nbp_per_tile = nbp/tot_nbt;

val rep_value = 100
var nb_leaf = tot_nbt;
params_scala("rep_value") = collection.mutable.Set(rep_value.toString)
params_ddt("output_dir") = collection.mutable.Set(cur_output_dir)
params_scala("output_dir") = collection.mutable.Set(cur_output_dir)
params_scala("ddt_main_dir") = collection.mutable.Set(ddt_main_dir)
params_scala("executor_cores") =  collection.mutable.Set(executor_cores)
params_scala("num_executors") =  collection.mutable.Set(num_executors)
params_ddt("nbt_side") =  collection.mutable.Set(nbt_side.toString)
params_wasure("nbt_side") =  collection.mutable.Set(nbt_side.toString)


println("")
println("=======================================================")
params_scala.map(x => println((x._1 + " ").padTo(15, '-') + "->  " + x._2.head))


// python
val python_cmd = params_python.to_command_line

// General c++ commands
val ply2geojson_cmd =  set_params(params_ddt, List(("step","ply2geojson"))).to_command_line
val tri2geojson_cmd =  set_params(params_ddt, List(("step","tri2geojson"))).to_command_line
val ply2dataset_cmd =  set_params(params_ddt, List(("step","ply2dataset"))).to_command_line
val extract_struct_cmd =  set_params(params_ddt, List(("step","extract_struct"))).to_command_line
val dump_ply_binary_cmd =  set_params(params_ddt, List(("step","dump_ply_binary"),("output_dir", cur_output_dir))).to_command_line
val id_cmd = List(build_dir + "/bin/identity-exe");

// Wausre surface reconstruction commands
val dim_cmd =  set_params(params_wasure, List(("step","dim"))).to_command_line
val dst_cmd =  set_params(params_wasure, List(("step","dst"))).to_command_line
val regularize_slave_focal_cmd =  set_params(params_wasure, List(("step","regularize_slave_focal"))).to_command_line
val regularize_slave_extract_cmd =  set_params(params_wasure, List(("step","regularize_slave_extract"))).to_command_line
val regularize_slave_insert_cmd =  set_params(params_wasure, List(("step","regularize_slave_insert"))).to_command_line
val extract_graph_cmd =  set_params(params_wasure, List(("step","extract_graph"))).to_command_line
val fill_graph_cmd =  set_params(params_wasure, List(("step","fill_graph"))).to_command_line
val ext_cmd =  set_params(params_wasure, List(("step","extract_surface"))).to_command_line
val tri2geojson_wasure_cmd =  set_params(params_wasure, List(("step","tri2geojson"))).to_command_line
val wasure_ply2geojson_cmd =  set_params(params_wasure, List(("step","ply2geojson"))).to_command_line

// =================================================
// ============  Parsing and init data ===========
var kvrdd_points: RDD[KValue] = sc.parallelize(List((0L,List(""))));
var kvrdd_inputs = format_data(
  params_scala,
  params_ddt,
  global_build_dir,
  ddt_main_dir,
  input_dir,
  df_par,
  sc ,
  iq
)

val t0 = System.nanoTime()
params_scala("t0") = collection.mutable.Set(t0.toString)
println("======== Tiling =============")
kvrdd_points = ddt_algo.compute_tiling(kvrdd_inputs,iq,params_ddt,params_scala);
ddt_algo.update_time(params_scala,"tiling");
nb_leaf = params_scala("nb_leaf").head.toInt;

if(nb_leaf == 1){
  max_opt_it = 2
}

var rep_merge = ((if((nb_leaf) < spark_core_max) spark_core_max else  nb_leaf));
var rep_loop = nb_leaf;
if(ndtree_depth > 8){
  rep_loop = spark_core_max*10;
  rep_merge = spark_core_max*10;
}
params_scala("rep_loop") = collection.mutable.Set(rep_loop.toString)
params_scala("rep_merge") = collection.mutable.Set(rep_merge.toString)


println("======== Dimenstionnality =============")
val res_dim = iq.run_pipe_fun_KValue(
  dim_cmd ++ List("--label", "dim"),
  kvrdd_points, "dim", do_dump = false).persist(slvl_glob)
val kvrdd_dim = iq.get_kvrdd(res_dim,"z");
val kvrdd_simp = iq.get_kvrdd(res_dim,"x").reduceByKey((u,v) => u ::: v,rep_loop);


// Test python
val res_dim_python = iq.run_pipe_fun_KValue(
  python_cmd,
  kvrdd_points, "dim", do_dump = false).persist(slvl_glob)

var input_ddt = kvrdd_points;

if(pscale > 0)
  input_ddt = kvrdd_simp;
ddt_algo.update_time(params_scala,"Dimensionality");

println("=========== Delauay triangulatin computation ==================")
val defaultV = (List(""));

val (kvrdd_tri_vertices,kvrdd_tri_edges,log_tri,stats_tri)  = ddt_algo.compute_ddt_nograph(
  kvrdd_points = input_ddt,
  iq = iq,
  params_cpp = params_ddt,
  params_scala = params_scala
);
kvrdd_tri_vertices.persist(slvl_glob);
kvrdd_tri_edges.persist(slvl_glob);
val kvrdd_tri_gid = ddt_algo.update_global_ids(kvrdd_tri_vertices,stats_tri,rep_merge,iq, params_ddt,sc)
kvrdd_tri_gid.persist(slvl_glob);
kvrdd_tri_gid.count()
kvrdd_tri_gid.unpersist()

val stats_kvrdd = kvrdd_simplex_id(stats_tri,sc)
val graph_stats = Graph(stats_kvrdd, kvrdd_tri_edges, List(""));
val input_dst = (kvrdd_tri_gid).union(kvrdd_dim).union(graph_stats.vertices).reduceByKey(_ ::: _,rep_merge).setName("input_dst");
input_dst.persist(slvl_glob);
input_dst.count()
res_dim.unpersist()

println("============= Simplex score computation ===============")
ddt_algo.update_time(params_scala,"compute_dst");
val res_dst = iq.run_pipe_fun_KValue(
  dst_cmd ++ List("--label", "dst"),
  input_dst, "dst", do_dump = false).persist(slvl_glob).setName("res_dst");
res_dst.count

input_dst.unpersist()
kvrdd_points.unpersist();
val kvrdd_dst = iq.get_kvrdd(res_dst,"t");

println("============= Regularize ===============")
ddt_algo.update_time(params_scala,"regularize");
// Message passing regularize
val res_regularize_extract = iq.run_pipe_fun_KValue(
  regularize_slave_extract_cmd ++ List("--label", "regularize_slave_extract"),
  kvrdd_dst, "regularize", do_dump = false).persist(slvl_glob).setName("res_reg");
res_regularize_extract.count()
res_dst.unpersist()
val rdd_shared_edges = iq.get_edgrdd(res_regularize_extract,"f")
val input_insert = (
  rdd_shared_edges.map(e => (e.dstId, e.attr)) union  kvrdd_dst
).reduceByKey(_ ::: _,rep_loop).persist(slvl_loop).setName("NEW_KVRDD_WITH_EDGES")
val res_regularize = iq.run_pipe_fun_KValue(
  regularize_slave_insert_cmd ++ List("--label", "regularize_slave_insert"),
  input_insert, "regularize", do_dump = false).persist(slvl_glob).setName("res_reg");

res_regularize.count()
res_regularize_extract.unpersist()
input_insert.unpersist()
val kvrdd_reg = iq.get_kvrdd(res_regularize,"t");
val kvrdd_shr = iq.get_edgrdd(res_regularize,"e")

val graph_reg = Graph(kvrdd_reg, kvrdd_tri_edges, List("")).partitionBy(EdgePartition1D,rep_merge);
graph_reg.vertices.setName("graph_reg");
graph_reg.edges.setName("graph_reg");

ddt_algo.update_time(params_scala,"optimization");
val lambda_list = params_scala("lambda").map(_.toDouble).toList.map(fmt.format(_).replace(',','.'))
var coef_mult_list = params_scala("coef_mult").map(_.toDouble).toList.sortWith(_ > _).map(fmt.format(_).replace(',','.'))
val algo_list = params_scala("algo_opt").toList


var stats_list_1 = new ListBuffer[(Int,(Float,Float))]()
var stats_list_2 = new ListBuffer[(Int,(Float,Float))]()
var stats_list_3 = new ListBuffer[(Int,(Float,Float))]()

val test_name = "coef_mult_lag"
var algo_id_acc = 0;

algo_list.foreach{ cur_algo =>
  if(cur_algo != "seg_lagrange_weight"){
    coef_mult_list = List("110000000000","110000000","110000")
  }
  lambda_list.foreach{ ll =>
    coef_mult_list.foreach{ coef_mult =>

      params_wasure("lambda") = collection.mutable.Set(ll)
      params_wasure("coef_mult") = collection.mutable.Set(coef_mult)
      val datestring = dateFormatter.format(Calendar.getInstance().getTime());
      val ext_name = test_name + "_" + cur_algo + "_" + algo_id_acc + "_ll_" + ll + "_cm_" + coef_mult + "_" + fmt.format(max_opt_it).replace(',','.') + "_" + cur_algo + "_"  + datestring;
      algo_id_acc+=1;


      val seg_cmd =  set_params(params_wasure, List(("step",cur_algo))).to_command_line;
      var input_seg2 =  (kvrdd_reg union kvrdd_shr.map(e => (e.srcId, e.attr))).reduceByKey(_ ::: _);
      input_seg2.persist(slvl_loop);
      input_seg2.count();
      var acc_loop = 0;

      val ext_cmd_vertex =  set_params(params_wasure, List(("step","extract_surface"),("area_processed","1"))).to_command_line
      val ext_cmd_edges =  set_params(params_wasure, List(("step","extract_surface"),("area_processed","2"))).to_command_line

        while (acc_loop < max_opt_it) {
          val t0_loop = System.nanoTime()
          val acc_loop_str = "%03d".format(acc_loop)
          val res_seg = iq.run_pipe_fun_KValue(
            seg_cmd ++ List("--label", "seg" + acc_loop_str) ++ (if(acc_loop == max_opt_it-1) List("--do_finalize") else List()),
            input_seg2
              ,cur_algo, do_dump = false).persist(slvl_loop)
          res_seg.count()
          val kvrdd_seg = iq.get_kvrdd(res_seg,"t");
          input_seg2.unpersist()
          
          val do_it_stats = (((acc_loop % stats_mod_it) == 0 || acc_loop == 1 || acc_loop == 3 || acc_loop == max_opt_it -1) && do_stats)

          if(acc_loop == max_opt_it -1){
            ddt_algo.update_time(params_scala,"final_before_dumping");
          }

          if(acc_loop == max_opt_it -1 || do_it_stats ){
            val rdd_local_edges = iq.get_edgrdd(res_seg,"e");
            val graph_seg = Graph(kvrdd_seg, rdd_local_edges, List("")).partitionBy(EdgePartition1D,rep_merge);

              val rdd_ply_surface = iq.run_pipe_fun_KValue(
                ext_cmd ++ List("--label","ext_seg" + ext_name + "_" + acc_loop_str),
                kvrdd_seg, "seg", do_dump = false).persist(slvl_glob)
              rdd_ply_surface.count()
              res_seg.unpersist()
              val ply_dir = cur_output_dir + "/tiles"// + "/plydist_" + ext_name + "_gc_" + algo_id_acc.toString + "_" + acc_loop_str
              ddt_algo.saveAsPly(rdd_ply_surface,ply_dir,plot_lvl);
              wasure_algo.partition2ply(cur_output_dir,algo_id_acc.toString,sc);
          }

          if(acc_loop != max_opt_it -1){
            val rdd_local_edges = iq.get_edgrdd(res_seg,"e");
            val rdd_shared_edges = iq.get_edgrdd(res_seg,"f")
            val rdd_stats  = iq.get_kvrdd(res_seg,"s")
            // Sum the stats between all the tiles
            var stats_1 = rdd_stats.map(x => x._2(0).split(" ").takeRight(2)).map(x => (x(0).toFloat,x(1).toFloat)).reduce( (x,y) => (x._1+y._1,x._2+y._2))
            var stats_2 = (0F,0F);
            var stats_3 = (0F,0F);

            input_seg2 = (
              rdd_local_edges.map(e => (e.srcId, e.attr)) union
                rdd_shared_edges.map(e => (e.dstId, e.attr)) union  kvrdd_seg
            ).reduceByKey(_ ::: _,rep_loop).persist(slvl_loop).setName("NEW_KVRDD_WITH_EDGES_" + acc_loop_str)
	    if(acc_loop % 10 ==0){
	        input_seg2.localCheckpoint()
	    }
            input_seg2.count()
            res_seg.unpersist()

            // Compute the global surface and compare
            if(do_it_stats){
              val seg_cmd_full =  set_params(params_wasure, List(("step","seg_global_extract"))).to_command_line
              val input_extract = kvrdd_seg.map(x => (0L,x._2)).reduceByKey(_ ::: _)
              val res_surface = iq.run_pipe_fun_KValue(
                seg_cmd_full ++ List("--label", "seg"),
                input_extract , "seg_lagrange", do_dump = false)
              val rdd_ply  = res_surface.filter(_(0) == 'p');
              val rdd_stats  = iq.get_kvrdd(res_surface,"s")
              val stats_collect = rdd_stats.map(x => x._2(0).split("]").tail(0).split(" ").filter(_.nonEmpty).map(_.toFloat)).collect()(0)
              stats_2 = (stats_collect(0),stats_collect(1))
              stats_3 = (stats_collect(2),stats_collect(3))
              rdd_stats.collect()
              if(acc_loop == 0){
                val ply_dir = cur_output_dir + "/tiles" // + "/plyglob_" + ext_name + "_gc_" + algo_id_acc.toString + "_" + acc_loop_str + "_global3"
                ddt_algo.saveAsPly(rdd_ply,ply_dir,plot_lvl)
                wasure_algo.partition2ply(cur_output_dir,algo_id_acc.toString,sc);
              }
            }
            val t1_loop = System.nanoTime();
            println("[it " + acc_loop_str + "]");
            ddt_algo.update_time(params_scala,"opt_loop" + algo_id_acc.toString);
            val t1 = ( acc_loop,stats_1)
            stats_list_1 += t1;
            println("   % overlap -> " + floatFormat.format(100*stats_1._1/stats_1._2.toFloat) + "% " + stats_1 );
            if(do_it_stats){
              val t2 = ( acc_loop,stats_2)
              val t3 = ( acc_loop,stats_3)
              stats_list_2 += t2;
              stats_list_3 += t3;
              println("   % GC      -> " + floatFormat.format(100*stats_2._1/stats_2._2.toFloat) + "% " + stats_2 );
              println("   % Energy  -> " + stats_3  );
            }
            println("")
          }
          acc_loop = acc_loop + 1;
        }

      if(do_stats){
        wasure_algo.dump_it_stats(cur_output_dir + "/" + ext_name + "_stats_conv.csv",stats_list_1,sc);
        wasure_algo.dump_it_stats(cur_output_dir + "/" + ext_name + "_stats_gtdiff.csv",stats_list_2,sc);
        wasure_algo.dump_it_stats(cur_output_dir + "/" + ext_name + "_stats_energy.csv",stats_list_3,sc);
      }
      stats_list_1.clear()
      stats_list_2.clear()
      stats_list_3.clear()
    }
  }
}

ddt_algo.update_time(params_scala,"final_after_dumping");
dump_json(params_scala,cur_output_dir + "/full_workflow_params.json",sc);


