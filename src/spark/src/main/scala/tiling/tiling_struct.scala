package tiling

import org.apache.spark.storage.StorageLevel
import scala.io.Source
import java.io._
import scala.xml._
import java.lang.Double
import sys.process._
import scala.xml.XML
import java.nio.file.{ Paths, Files }
import java.util.concurrent.Executors
import scala.concurrent._
import scala.collection.parallel._

import org.apache.spark.rdd.RDD;

import java.util.Calendar
import java.text.SimpleDateFormat;

import spark_ddt.core._;
import spark_ddt.core.IQlibCore._;
import spark_ddt.util._;
import params_parser._;



class Bbox(val D: Int) {
  var tab: Array[scala.Double] = new Array[scala.Double](D * 2);
  def this( tstring : String ) {
    this(0)
    tab = tstring.split(" +").map(_.toDouble)
  }

   def centroid() : String = {
     "["+ ((tab(0) +tab(1))/2) + "," +((tab(2) +tab(3))/2)+"]"
   }

  def toGeojson() : String =  {
    val dim = tab.size/2
    var tab_s = {
      "[["+tab(0)+","+tab(2)+"]," +
      "["+tab(0)+","+tab(3)+"]," +
      "["+tab(1)+","+tab(3)+"]," +
      "["+tab(1)+","+tab(2)+"]," +
      "["+tab(0)+","+tab(2)+"]]" 
    }

    return raw"""
    {
    "type": "Feature",
    "geometry": {
      "type": "Polygon",
      "coordinates": [ $tab_s ]
      }
    }  
     """
  }
}



class nd_tree(dim : Int,depth : Int,iq : IQlibSched) extends Serializable {
  import org.apache.spark.rdd.RDD
  import scala.collection.mutable.ListBuffer

  val d_tab = Range(0,depth).map(this.nb_nodes(_))
  def nb_nodes(dp : Int) : Int = {
    return ((math.pow(math.pow(2,dim),dp+1)-1)/(math.pow(2,dim)-1)).toInt
  }

  def get_depht(id : Int) : Int = {
    d_tab.indexWhere(_ > id)
  }

  def generate_id(dp : Int): Range = {
    return Range(nb_nodes(dp-2),nb_nodes(dp-1))
  }

  def generate_id(): Range = {
    return generate_id(depth)
  }


  def coords2id(coords : Array[Int], dp : Int) : Int =  {
    var acc=0;
    var id=0;
    val side=math.pow(2,dp).toInt;
    for( cc <- coords){
      id=id + cc*math.pow(side,acc).toInt
      acc=acc+1;
    }
    return id + nb_nodes(dp-1)
  }

  def id2coords(id : Int) : Array[Int] =  {
    val dp = this.get_depht(id);
    var n : Int = depth_id(id)
    val side=math.pow(2,dp).toInt;
    val coords = Array.range(0,dim)
    var acc=0

    while(acc<dim){
      var r=n%side;
      coords(acc)=r;
      n=(n-r)/side;
      acc=acc+1
    }
    return coords
  }

  def is_root(rid : Int, id : Int ) : Boolean = {
    var cid = id;
    while(cid != 0){
      if(cid == rid)
        return true;
      else
        cid = root_id(cid);
    }
    return false
  }

  def root_id(id : Int) : Int = {
    val dp = get_depht(id);
    val coords = id2coords(id);
    return coords2id(coords.map(_/2),dp-1)
  }

  def depth_id(id : Int) : Int = {
    val dp = get_depht(id);
    val lid = id - d_tab(dp-1)
    return  lid
  }



  def compute_id_map(l1 : RDD[(Int, Int)], lim : Int, min_nbt : Int = 1) : scala.collection.Map[Int,Int] = {
    return compute_id_map_break_lineage(l1,lim,min_nbt);

  }


  def compute_id_rdd(l1 : RDD[(Int, Int)], lim : Int, min_nbt : Int = 1) : (RDD[(Long, Long)],RDD[(Int, (Int, List[Int]))]) = {
    var li_new = l1.map(x=> (x._1,(x._2,List(x._1))));
    li_new.cache().count

    var acc = 0;
    var lf_new = li_new.filter(x => x._2._1 > lim).cache()
    lf_new.count
    println("Distributed Nd tree construction ...")
    while(li_new.count > 1){
      val root_li_new = li_new.map(x => (root_id(x._1),x._2)).reduceByKey((x,y) => (x._1+y._1,x._2 ::: y._2)).cache().setName("RDD_ROOT_LI_NEW" + acc)
      val root_li_stage = li_new.map(x => (root_id(x._1),(x._2._1,List(x._1)))).reduceByKey((x,y) => (x._1+y._1,x._2 ::: y._2)).flatMap( x=> x._2._2.map(y => (y , x._2._1)))
      val li_new_filt = li_new.cogroup(root_li_stage).filter(
        cc => ((!cc._2._2.isEmpty) && (!cc._2._1.isEmpty))).filter(
        cc => cc._2._2.head > lim && cc._2._1.head._1 <=lim).map(cc => (cc._1,cc._2._1.head));
      val unn = (lf_new union li_new_filt).cache.setName("RDD_UNN_" + acc)
      unn.count
      lf_new.unpersist()
      lf_new = unn;

      root_li_new.count
      li_new.unpersist()
      li_new = root_li_new
      acc = acc + 1;
    }
    if(lf_new.count() == 0)
      lf_new = li_new.map(x => (1,x._2)).reduceByKey((x,y) => (x._1+y._1,x._2 ::: y._2));
    val fff_map = lf_new.flatMap( x=> x._2._2.map(y => (y.toLong,x._1.toLong))).persist(iq.get_storage_level()).setName("KDTREE_MAP")
    lf_new.cache().count
    fff_map.count
    li_new.unpersist()


    return (fff_map,lf_new)
  }




  def compute_id_map_break_lineage(l1 : RDD[(Int, Int)], lim : Int, min_nbt : Int = 1) : scala.collection.Map[Int,Int] = {
    var li_new = l1.map(x=> (x._1,(x._2,List(x._1))));
    li_new.persist(iq.get_storage_level()).setName("ND_TREE_INIT")

    var lf_new = li_new.filter(x => x._2._1 > lim).flatMap( x=> x._2._2.map(y => (y,x._1))).collectAsMap()
    println("Nd tree construction ...")
    val list_map : ListBuffer[scala.collection.Map[Int,Int]] = ListBuffer()
    while(li_new.count > 1){
      val root_li_new = li_new.map(x => (root_id(x._1),x._2)).reduceByKey((x,y) => (x._1+y._1,x._2 ::: y._2),1000).persist(iq.get_storage_level()).setName("ND_TREE_LOOP")
      val root_map = root_li_new.collect().toMap
      list_map += li_new.filter(x => root_map(root_id(x._1))._1 > lim && x._2._1 <=lim).flatMap( x=> x._2._2.map(y => (y,x._1))).collectAsMap()
      li_new = root_li_new

    }
    println("Constructing the nd tree id map")
    val fff_map = list_map.reduce(_ ++ _) ++ lf_new
    return fff_map
  }
}



object geojson_export{


def export_graph(kvrdd_bbox : RDD[KValue],gg : TGraph,params : params_map,prefix : String = "")  {
  val geojson_header= raw"""{
    "type": "FeatureCollection",
    "features": [
    """
  val geojson_footer = raw"""
     ]
     }
  """

   // Bbox
  var geojson_bbox : String =geojson_header;
  val b_list=kvrdd_bbox.collect().map(x => (x._1,new Bbox(x.content.toString)));
  val b_map = b_list.toMap
  val str_bbox =   "["+gg.edges.collect().map(x => "["+ b_map(x.srcId).centroid +  "," + b_map(x.dstId).centroid +"]").reduce(_+","+_)+"]";
  geojson_bbox += b_list.map(x => x._2.toGeojson).reduce(_ + ","+_)
  geojson_bbox+=geojson_footer;

  // Graph txt
  var geojson_graph : String =geojson_header;
  geojson_graph += raw"""
    {
    "type": "Feature",
    "geometry": {
      "type": "MultiLineString",
      "coordinates":  $str_bbox 
      }
    }  
     """
  geojson_graph+=geojson_footer;


  val bname =  params("output_dir").head +"/graph" + prefix;
  val pw1 = new PrintWriter(new File(bname + "_bbox.geojson" ))
    ("cp " + params("ddt_main_dir").head +"/src/qgis/graph_polygon.qml " +  (bname + "_bbox.qml")).!
    pw1.write(geojson_bbox)
  pw1.close
  val pw2 = new PrintWriter(new File(bname + "_graph.geojson" ))
    ("cp " + params("ddt_main_dir").head +"/src/qgis/graph_linestring.qml " +  (bname + "_graph.qml")).!
    pw2.write(geojson_graph)
  pw2.close
}

}
