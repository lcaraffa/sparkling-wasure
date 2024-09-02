package spark_ddt.util

import org.apache.spark.rdd.RDD;
import org.apache.hadoop.fs.FileSystem
import org.apache.spark.SparkContext
import scala.collection.mutable.ListBuffer
// import scala.util.parsing.json._
import scala.xml.XML
import java.nio.file.{ Paths, Files }
import java.io.File
import org.apache.hadoop.fs.Path
import java.io._
import scala.language.postfixOps

import java.io.BufferedOutputStream
import spark_ddt.core._;
import spark_ddt.core.IQlibCore._;

import java.util.Date;
import java.text.SimpleDateFormat;
import java.util.Calendar;

import org.apache.hadoop.fs.permission.FsPermission
import org.apache.spark.sql._
import org.apache.spark.sql.types._
import org.apache.spark.sql.types.DataType._
import org.apache.spark.storage.StorageLevel
import collection.mutable


object files_opt {
  def getListOfFiles(dir: String): List[File] = {
    val d = new File(dir)
    if (d.exists && d.isDirectory) {
      d.listFiles.filter(_.isFile).toList
    } else {
      List[File]()
    }
  }
}

object strings_opt {
  def get_bname(ss: String) = { ss.substring(0, ss.lastIndexOf('.')) }
  def parsable(ss: String) = { ss.replaceAll("\\s", "").replaceAll("//", "/") };
}

object bash_funcs {
  def get_bash_variable(vname: String,def_val: String = ""): String = {
    try {
      val res = sys.env(vname)
      return res
    } catch {
      case _:
          Throwable => return def_val
    }
  }


  def print_section(sec_name:String)  {
    println(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
    println("               " + sec_name + "          ")
    println(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
  }

}

object xml_parsing {

  def get_param_from_xml_old(ss: String, lab_name: String, def_val: String = ""): String = {
    if (!does_exist(ss))
      return def_val
    val xml = XML.loadFile(ss);
    val xml_params = xml \ "params";
    for (xml_param <- xml_params(0).child) {
      val label = xml_param.label;
      if (label == lab_name) {
        val value: String = xml_param.text;
        return value;
      }
    }
    return def_val;
  }

  def get_param_from_xml(ss: String, lab_name: String, def_val: String = ""): String = {
    if (!does_exist(ss))
      return def_val
    val xml = XML.loadFile(ss);

    for (xml_params <- xml(0).child) {
      for (xml_param <- xml_params(0).child) {
        val label = xml_param.label;
        if (label == lab_name) {
          val value: String = xml_param.text;
          return value;
        }
      }
    }
    return def_val;
  }


  def get_list_param_from_xml(ss: String, lab_name: String, def_val: String = ""): List[String] = {
    if (!does_exist(ss))
      return List(def_val)
    val xml = XML.loadFile(ss);
    var ll : List[String] = List()
    for (xml_params <- xml(0).child) {
      for (xml_param <- xml_params(0).child) {
        val label = xml_param.label;
        if (label == lab_name) {
          val value: String = xml_param.text;
          ll = ll :+ value;
        }
      }
    }
    return ll;
  }


  def get_list_param_from_xml_old(ss: String, lab_name: String, def_val: String = ""): List[String] = {
    if (!does_exist(ss))
      return List(def_val)
    val xml = XML.loadFile(ss);
    val xml_params = xml \ "params";
    var ll : List[String] = List()
    for (xml_param <- xml_params(0).child) {
      val label = xml_param.label;
      if (label == lab_name) {
        val value: String = xml_param.text;
        ll = ll :+ value;
      }
    }
    return ll;
  }


  def does_exist(ss: String): Boolean = {
    return Files.exists(Paths.get(ss));
  }

  def does_exist(lls: List[String]): Boolean = {
    if (lls.length > 0)
      return lls.map(does_exist(_)).reduce(_ && _);
    else false;
  }
}



object algo_stats {


  def sum_simplex_id( stats_tri_last : RDD[KValue]) : List[(Int, List[Int])] = {
    val stats_tri_last_sorted = stats_tri_last.map(x => (x.key,x.content.split(" ").map(_.toInt))).sortByKey().collect.map(x => (x._1.toInt,x._2))
    val array_id = (stats_tri_last_sorted.map(x => (x._1)) :+ -1).drop(1)
    val stats_shifted = array_id.zip(stats_tri_last_sorted).map(x => (x._1,x._2._2))
    val res =  stats_shifted.map(x => List(x)).foldLeft(List((stats_tri_last_sorted.head._1,List(0,0,0))))((res,y) => res ::: List((y(0)._1,res.last._2.zip(y.head._2).map { case (x, y) => x + y } )))
    return res
  }

  def kvrdd_simplex_id( stats_tri_last : RDD[KValue], sc : SparkContext) : RDD[KValue] = {
    val summed_id = sum_simplex_id(stats_tri_last).dropRight(1);
    return sc.parallelize(summed_id.map( x =>
      (x._1.toLong, List("s 1 " + x._1.toString +  " s " + x._2.map(_.toString).reduce((y,z) => y + " " + z)))
    ))
  }



  def compute_time(dat_stat: Array[(String, String)], label_filter : List[String]) : Float = {
    if(dat_stat.size == 0){
      return 0
    }else{
      val ll = dat_stat.filter(x => (label_filter.map( y => x._1 contains y).reduce(_ && _)))
      return ll.map( x => x._2.toFloat).reduce(_ + _)
    }
  }


  def dump_stats(kvrdd_ll : RDD[KValue],stats_filename : String, sc : SparkContext ) {
    val fs = FileSystem.get(sc.hadoopConfiguration);
    val kvrdd_time = kvrdd_ll.map(x => (0,x._2)).reduceByKey(_ ++ _).values.collect()(0).map(y => y.split(" ", 5).last.split(";").map(x => x.split(":")).map(x => (x(0),x(1)))).reduce(_ ++ _)
    val list_label = kvrdd_time.map( x => x._1.filter(!_.isDigit)).distinct;
    val output = fs.create(new Path(stats_filename));
    val os = new BufferedOutputStream(output)
    os.write((list_label.reduce(_ + "," + _) + "\n").getBytes);
    os.write(list_label.map(x => compute_time(kvrdd_time,List(x)).toString).reduce(_ ++ "," ++ _).getBytes)
    os.close()
  }


}


object dataset_processing {
  def ply2type(typename: String) : DataType = {
    typename match {
      case "uchar" => ByteType case "uint8" => ByteType
      case "char" => ByteType case "int8" => ByteType
      case "ushort" => ShortType case "uint16" => ShortType
      case "short" => ShortType case "int16" => ShortType
      case "uint" => IntegerType case "uint32" => IntegerType
      case "int" => IntegerType case "int32" => IntegerType
      case "ulong" => LongType case "uint64" => LongType
      case "long" => LongType case "int64" => LongType
      case "float" => FloatType case "float32" => FloatType
      case "double" => DoubleType case "float64" => DoubleType
      case "list" => StringType
      case other => sys.error(s"Unsupported type $other")
    }
  }


  def type2ply(datatype: DataType) : String = {
    datatype match {
      case ByteType => "uchar"
      case ShortType => "int16"
      case IntegerType => "int32"
      case LongType => "int64"
      case FloatType => "float"
      case DoubleType => "double"
      case  StringType => "list uchar int"
      case other => sys.error(s"Unsupported type $other")
    }
  }


  def get_header_schema( header_str: String,elem_pos: Int) : StructType = {
    val splited_header = header_str.split("element")(elem_pos)
    StructType(splited_header.split(";").filter( _ contains "property").map(_.split(" ")).map(elems =>
      List(StructField(elems.last, ply2type(elems(1))))).reduce(_ ::: _))
  }

}



object params_parser {

  type Hash_StringSeq = mutable.HashMap[String, mutable.Set[String]]
  type params_map = Hash_StringSeq with mutable.MultiMap[String, String]
  type Param_List = ListBuffer[params_map]

  def format_data(
    params_scala : params_map,
    params_ddt : params_map,
    global_build_dir : String,
    ddt_main_dir : String,
    input_dir : String,
    df_par : Int,
    sc : SparkContext,
    iq : IQlibSched

  ) : RDD[KValue] = {
    import files_opt._

    val fs = FileSystem.get(sc.hadoopConfiguration);
    val datatype =  params_scala.get_param("datatype", "")
    val regexp_filter = params_scala.get_param("regexp_filter", "");
    val slvl_glob = StorageLevel.fromString(params_scala.get_param("StorageLevel", "DISK_ONLY"))
    var rep_value = df_par.toLong;

    val generate_points_cmd =  set_params(params_ddt,  List(("step","generate_points_" + datatype))).to_command_line
    val ser2datastruct_cmd =  set_params(params_ddt, List(("step","serialized2datastruct"))).to_command_line
    var kvrdd_points: RDD[KValue] = sc.parallelize(List((0L,List(""))));
    var kvrdd_inputs: RDD[KValue] = sc.parallelize(List((0L,List(""))));
    var kvrdd_inputs_struct  : RDD[KValue] = sc.parallelize(List((0L,List(""))));
    println("")
    println("======== GENERATE DATA =============")
    println("datatype : " + datatype)
    datatype match {
      case "random_normal" | "random_uniform" => {
        val nbp =  params_scala.get_param("nbp", "10000").toLong
        val max_ppt = params_scala.get_param("max_ppt", "10000").toLong
        val nbp_and_depth = params_scala.get_param("nbp_and_depth", "")
        if(nbp_and_depth contains ":"){
          params_scala("nbp") =  collection.mutable.Set(nbp_and_depth.split(":")(0))
          params_scala("ndtree_depth") = collection.mutable.Set(nbp_and_depth.split(":")(1))
        }

        rep_value = nbp/max_ppt;
        var nb_kernel = ((nbp/max_ppt)*3)
        params_scala("nb_kernel") = collection.mutable.Set(nb_kernel.toString)

        println("nb kernels : " +nb_kernel)
        if(nb_kernel == 0){
          println(" ====== WARNING =====")
          println(" NB_kernel == 0, set to 1")
          nb_kernel = 1
        }
        rep_value = nbp/max_ppt
        val nbp_per_kernel = nbp/nb_kernel;
        val input_rdd_raw: RDD[KValue] = sc.parallelize(List.range(0,nb_kernel)).map(
          x => (x.toLong, List(""))).repartition(rep_value.toInt)
        val raw_inputs = iq.run_pipe_fun_KValue(
          generate_points_cmd ++ List("--label","lab2","--nbp",nbp_per_kernel.toString),
          input_rdd_raw, "generate_points", do_dump = false)
        kvrdd_inputs = iq.get_kvrdd(raw_inputs, "z",txt="pts").persist(slvl_glob).setName("KVRDD_INPUT");

      }
      case "files" => {
        println("")
        println("======== LOAD DATA  file =============")
        val ss_reg = regexp_filter.r
        val ply_input = getListOfFiles(input_dir).filter(
          x => ((x.toString endsWith ".ply") && ((ss_reg).findFirstIn(x.toString).isDefined)))
        kvrdd_inputs = iq.get_kvrdd(sc.parallelize(ply_input.map(
          fname => "p 1 " + ("([0-9]+)".r).findAllIn(fname.toString).toArray.last + " f " + fname.toString)),"p")
        kvrdd_inputs_struct = kvrdd_inputs
      }

      case "laz" => {
        println("")
        println("======== LOAD DATA  file =============")
        val ss_reg = regexp_filter.r
        val ply_input = getListOfFiles(input_dir).filter(
          x => ((x.toString endsWith ".laz") && ((ss_reg).findFirstIn(x.toString).isDefined)))
        kvrdd_inputs = iq.get_kvrdd(sc.parallelize(ply_input.map(
          fname => "p 1 " + ("([0-9]+)".r).findAllIn(fname.toString).toArray.last + " f " + fname.toString)),"p")
        kvrdd_inputs_struct = kvrdd_inputs//iq.get_kvrdd(struct_inputs)
      }

      case "plystream" => {
        println("")
        println("======== LOAD DATA filestream =============")
        val ss_reg = regexp_filter.r
        val nb_file = fs.listStatus(new Path(input_dir)).map(x => fs.listStatus(x.getPath)).reduce(_ ++ _).map(_.getPath).filter(
          x => ((x.toString endsWith ".ply")) && ((ss_reg).findFirstIn(x.toString).isDefined)
        ).size

        val ply_input = fs.listStatus(new Path(input_dir)).map(x => fs.listStatus(x.getPath)).reduce(_ ++ _).map(_.getPath).filter(
          x => ((x.toString endsWith ".ply")) && ((ss_reg).findFirstIn(x.toString).isDefined)
        )

        kvrdd_inputs_struct = sc.wholeTextFiles(input_dir + "*.ply").map(x =>
          x._2.replaceAll("[\n\r](p|f|e|c)",";$1").replace("end_header", "end_header;").replace("\n"," ")).zipWithIndex.map(
          e => (e._2.toLong,List("g 1 " + e._2.toString + " s " +  e._1.toString))
        )

      }

      case "b64stream" => {
        println("")
        println("======== LOAD DATA filestream =============")
        val ss_reg = regexp_filter.r
        val nb_file = fs.listStatus(new Path(input_dir)).map(x => fs.listStatus(x.getPath)).reduce(_ ++ _).map(_.getPath).filter(
          x => ((x.toString endsWith ".stream")) && ((ss_reg).findFirstIn(x.toString).isDefined)
        ).size

        val ply_input = fs.listStatus(new Path(input_dir)).map(x => fs.listStatus(x.getPath)).reduce(_ ++ _).map(_.getPath).filter(
          x => ((x.toString endsWith ".stream")) && ((ss_reg).findFirstIn(x.toString).isDefined)
        )

        kvrdd_inputs_struct = sc.textFile(input_dir + "*.stream").zipWithIndex.map(
          e => (e._2.toLong,List("z 1 " + e._2.toString + " z " +  e._1.toString))
        ).repartition(nb_file).setName("KVRDD_INPUT")
      }
      case _ => {
        println(" ERROR : DATATYPE " + datatype + " does not exists")
      }


    }


    return kvrdd_inputs_struct

  }


  def dump_json(pm : params_map,json_filename : String, sc : SparkContext)  = {
    println()
  }


  def get_param_map(xml_params : scala.xml.Node) : params_map = {
    val mm1 = new Hash_StringSeq with mutable.MultiMap[String, String]
    for (xml_param <- xml_params(0).child) {
      val label = xml_param.label;

      if(label != "#PCDATA"){
        val value: String = xml_param.text
        for(it <- value.split(","))
          mm1.addBinding(label,it)
      }
    }
    if(mm1.size !=0){
      mm1.addBinding("name",xml_params.label)
    }
    return mm1;
  }

  def parse_xml_datasets(ss : String) : Param_List = {
    val ll_params : Param_List  = ListBuffer()
    val xml = XML.loadFile(ss);
    val xml_datasets = xml \ "datasets";
    for (xml_params <- xml_datasets(0).child) {
      val mm1 = get_param_map(xml_params)
      if(mm1.size !=0 && mm1.get_param("do_process", "true").toBoolean ){
        val do_expand_joker = mm1.get_param("do_expand", "false").toBoolean
        if(do_expand_joker){

          val param_list_2 = new Param_List()
          val temp_map = new Hash_StringSeq with mutable.MultiMap[String, String]
          for((kk,svv) <- mm1){
            temp_map(kk) =  scala.collection.mutable.Set(svv.head)
          }
          param_list_2 += temp_map
          for((kk,svv) <- mm1){
            val param_list_3 = new Param_List()
            for(vv <- svv.tail){
              for( cur_mm1 <- param_list_2){
                val ccpm = new Hash_StringSeq with mutable.MultiMap[String, String]
                for((kk2,svv2) <- cur_mm1){
                  if(kk2 == kk)
                    ccpm(kk2) = scala.collection.mutable.Set(vv)
                  else
                    ccpm(kk2) = svv2
                }
                param_list_3 += ccpm
              }
            }
            param_list_2 ++= param_list_3
          }
          ll_params ++= param_list_2
        }else{
          ll_params += mm1
        }

      }
    }
    return ll_params;
  }




  def parse_xml_datasets_2(ss : String) : Param_List = {
    val ll_params : Param_List  = ListBuffer()
    val xml = XML.loadFile(ss);
    val xml_datasets = xml \ "datasets";
    for (xml_params <- xml_datasets(0).child) {
      val mm1 = get_param_map(xml_params)
      if(mm1.size !=0 && mm1.get_param("do_process", "true").toBoolean ){
        val do_expand_joker = mm1.get_param("do_expand", "false").toBoolean
        if(do_expand_joker){
          val sort_criterion = mm1.get_param("sort_criterion", "")
          // =============
          var param_list_2 = new Param_List()
          val temp_map = new Hash_StringSeq with mutable.MultiMap[String, String]
          for((kk,svv) <- mm1){
            temp_map(kk) =  scala.collection.mutable.Set(svv.head)
          }
          val nbp_and_depth = temp_map.get_param("nbp_and_depth", "")
          if(nbp_and_depth contains ":"){
            temp_map("nbp") =  collection.mutable.Set(nbp_and_depth.split(":")(0))
            temp_map("ndtree_depth") = collection.mutable.Set(nbp_and_depth.split(":")(1))
          }
          param_list_2 += temp_map
          for((kk,svv) <- mm1){
            val param_list_3 = new Param_List()
            for(vv <- svv.tail){
              for( cur_mm1 <- param_list_2){
                val ccpm = new Hash_StringSeq with mutable.MultiMap[String, String]
                for((kk2,svv2) <- cur_mm1){
                  if(kk2 == kk)
                    ccpm(kk2) = scala.collection.mutable.Set(vv)
                  else
                    ccpm(kk2) = svv2
                }
                val nbp_and_depth = ccpm.get_param("nbp_and_depth", "")
                if(nbp_and_depth contains ":"){
                  ccpm("nbp") =  collection.mutable.Set(nbp_and_depth.split(":")(0))
                  ccpm("ndtree_depth") = collection.mutable.Set(nbp_and_depth.split(":")(1))
                }
                param_list_3 += ccpm
              }
            }
            param_list_2 ++= param_list_3
          }

          if(!(sort_criterion isEmpty)){
            val s_cc =  sort_criterion.split(":")(0)
            val s_oo =  sort_criterion.split(":")(1)
            if(s_oo == "gt"){
              println("sort gt")
              val param_list_sorted = param_list_2.sortWith((x,y) => x(s_cc).head.toInt > y(s_cc).head.toInt)
              param_list_2 = param_list_sorted;
            }else{
              println("sort lt")
              val param_list_sorted = param_list_2.sortWith((x,y) => x(s_cc).head.toInt < y(s_cc).head.toInt)
              param_list_2 = param_list_sorted;
            }
          }

          ll_params ++= param_list_2
        }else{
          ll_params += mm1
        }

      }
    }
    return ll_params;
  }

  def parse_xml_datasets_string(ss : String) : Param_List = {
    val ll_params : Param_List  = ListBuffer()
    val xml = XML.loadString(ss);
    val xml_datasets = xml \ "datasets";
    for (xml_params <- xml_datasets(0).child) {
      val mm1 = get_param_map(xml_params)
      if(mm1.size !=0 && mm1.get_param("do_process", "true").toBoolean ){
        val do_expand_joker = mm1.get_param("do_expand", "false").toBoolean
        if(do_expand_joker){
          val sort_criterion = mm1.get_param("sort_criterion", "")
          // =============
          var param_list_2 = new Param_List()
          val temp_map = new Hash_StringSeq with mutable.MultiMap[String, String]
          for((kk,svv) <- mm1){
            temp_map(kk) =  scala.collection.mutable.Set(svv.head)
          }
          val nbp_and_depth = temp_map.get_param("nbp_and_depth", "")
          if(nbp_and_depth contains ":"){
            temp_map("nbp") =  collection.mutable.Set(nbp_and_depth.split(":")(0))
            temp_map("ndtree_depth") = collection.mutable.Set(nbp_and_depth.split(":")(1))
          }
          param_list_2 += temp_map
          for((kk,svv) <- mm1){
            val param_list_3 = new Param_List()
            for(vv <- svv.tail){
              for( cur_mm1 <- param_list_2){
                val ccpm = new Hash_StringSeq with mutable.MultiMap[String, String]
                for((kk2,svv2) <- cur_mm1){
                  if(kk2 == kk)
                    ccpm(kk2) = scala.collection.mutable.Set(vv)
                  else
                    ccpm(kk2) = svv2
                }
                val nbp_and_depth = ccpm.get_param("nbp_and_depth", "")
                if(nbp_and_depth contains ":"){
                  ccpm("nbp") =  collection.mutable.Set(nbp_and_depth.split(":")(0))
                  ccpm("ndtree_depth") = collection.mutable.Set(nbp_and_depth.split(":")(1))
                }
                param_list_3 += ccpm
              }
            }
            param_list_2 ++= param_list_3
          }

          if(!(sort_criterion isEmpty)){
            val s_cc =  sort_criterion.split(":")(0)
            val s_oo =  sort_criterion.split(":")(1)
            if(s_oo == "gt"){
              println("sort gt")
              val param_list_sorted = param_list_2.sortWith((x,y) => x(s_cc).head.toLong > y(s_cc).head.toLong)
              param_list_2 = param_list_sorted;
            }else{
              println("sort lt")
              val param_list_sorted = param_list_2.sortWith((x,y) => x(s_cc).head.toLong < y(s_cc).head.toLong)
              param_list_2 = param_list_sorted;
            }
          }

          ll_params ++= param_list_2
        }else{
          ll_params += mm1
        }

      }
    }
    return ll_params;
  }

  def parse_xml_params(ss : String) : params_map = {
    val xml = XML.loadFile(ss);
    val xml_datasets = xml \ "params";
    val mm = get_param_map(xml_datasets(0).child.head)
    return mm
  }


  def get_command_line_from_param(pmap : params_map) : List[String] = {
    pmap.map{case(k,v)=>List("--" + k, v.head)}.reduce(_ ++_)
  }


  // Function to generate a random 10-character string
  def randomString(length: Int): String = {
    val chars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
    val sb = new StringBuilder(length)
    for (_ <- 1 to length) {
      val randomNum = scala.util.Random.nextInt(chars.length)
      sb.append(chars.charAt(randomNum))
    }
    sb.toString()
  }


  def set_params( cur_map : params_map,vals : List[(String,String)]) : params_map = {
    val nmap = new Hash_StringSeq with mutable.MultiMap[String, String]
    nmap ++= cur_map
    for(cc <- vals){
      nmap.remove(cc._1)
      nmap.addBinding(cc._1,cc._2)
    }
    return nmap
  }

  implicit class cparam(pmap: params_map) {
    def to_command_line = {
      pmap("exec_path").head.split(" ").toList  ++ pmap.map{
        case(k,v)=> (if(k != "exec_path") List("--" + k, v.head) else List(""))
      }.reduce(_ ++ _)}.toList

    def get_param(ss : String,df : String) : String = {
      if(pmap.keySet.exists(_ == ss))
        return pmap(ss).head
      else {
        pmap(ss) = collection.mutable.Set(df)
        return df
      }
    }

    def exists(ll : List[String]) = { ll.map(x => pmap.keySet.exists(_ == x)).reduce(_ && _)}
  }

}

