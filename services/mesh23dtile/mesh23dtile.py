import sys
import os
import re
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt
import random
import json
import string
from pathlib import Path
from datetime import datetime
import pymeshlab
from octree import *
import pyproj
import trimesh

num_colors = 20
cmap = plt.cm.get_cmap('tab20c', num_colors)
size_s = 10  # number of characters in the string.
target_face_num=100000

max_depth = 5
geom_error = [20,10,5,2,1,0]
transformer = pyproj.Transformer.from_crs("epsg:2154", "epsg:4979")
transformer_glob = pyproj.Transformer.from_crs("epsg:2154", "epsg:4978")

def is_inside_bbox(bbox,pts) :
    return ((pts[0] > bbox[0] and pts[0] < bbox[1]) and
            (pts[1] > bbox[2] and pts[1] < bbox[3]) and
            (pts[2] > bbox[4] and pts[2] < bbox[5]))

def get_bb_center(bb) :
    return ((bb[1]+bb[0])/2,(bb[3]+bb[2])/2,(bb[5]+bb[4])/2)

def bbox_union(bb1,bb2) :
    if len(bb1) == 0 :
        return bb2
    if len(bb2) == 0 :
        return bb1
    return [min(bb1[0],bb2[0]),max(bb1[1],bb2[1]),
            min(bb1[2],bb2[2]),max(bb1[3],bb2[3]),
            min(bb1[4],bb2[4]),max(bb1[5],bb2[5])]

def trans_bbox(bb,coords) :
    x_shift = float(coords.split("x")[0])
    y_shift = float(coords.split("x")[1])
    # Loop through all OBJ files in the input directory
    xmin_new, ymin_new, zmin_new = transformer.transform(bb[0]+x_shift,bb[2]+y_shift,bb[4])
    xmax_new, ymax_new, zmax_new = transformer.transform(bb[1]+x_shift,bb[3]+y_shift,bb[5])
    return [xmin_new, xmax_new, ymin_new,ymax_new, zmin_new, zmax_new]

def bbox2region(bb,coords) :
    x_shift = float(coords.split("x")[0])
    y_shift = float(coords.split("x")[1])
    xmin_new, ymin_new, zmin_new = transformer.transform(bb[0]+x_shift,bb[2]+y_shift,bb[4])
    xmax_new, ymax_new, zmax_new = transformer.transform(bb[1]+x_shift,bb[3]+y_shift,bb[5])
    return list([ymin_new,xmin_new,ymax_new,xmax_new,zmin_new,zmax_new])


def bbox23Dbox(bb) :
    cc = list(get_bb_center(bb))
    cc += [abs(bb[1]-bb[0])/2,0,0,
           0,abs(bb[3]-bb[2])/2,0,
           0,0,abs(bb[5]-bb[4])/2]
    return np.around(cc,5).tolist()


class tree_obj(object):
    def __init__(self, name, position,bbox):
        self.name = name
        self.position = position
        self.bbox = bbox
    def get_name(self):
         return u"{0}".format(self.name)
    def __str__(self):
        return u"{0}".format(self.name)

def print_node(depth,str_node,list_sons) :
    print("node depth " + str(depth) + " : " + str_node + " <= " + ' '.join(list_sons))    

#def depth2geomerr(depth) :
    
    
def node2dict(bbox,name,children,depth,coords) :
    leaf_dict = {}

    #import pdb; pdb.set_trace() 
    #leaf_dict["boundingVolume"] = { "region": bbox }
    if coords :
        leaf_dict["boundingVolume"] = { "region": bbox2region(bbox,coords) }
    else :
        leaf_dict["boundingVolume"] = { "box": bbox23Dbox(bbox) }

    leaf_dict["raw_bbox"] = { "box" :  bbox }        
    leaf_dict["content"] =  { "uri" :  str(name)  }
    leaf_dict["geometricError"] = str(geom_error[depth])
    leaf_dict["refine"] = "REPLACE"
    leaf_dict["children"] = children
    return leaf_dict


def merge_subtree(node_tt,depth,output_dir,coords) :
    joint_string = []
    node_name_obj = "tiles/" + str(depth) + "_" + ''.join(random.choices(string.ascii_uppercase + string.digits, k = size_s))  + '.obj'
    node_name_ply = "tiles/" + str(depth) + "_" + ''.join(random.choices(string.ascii_uppercase + string.digits, k = size_s))  + '.ply'
    file_name = output_dir +  node_name_obj
    full_bbox = []
    children_dict = []
    if node_tt.isLeafNode :
        for x in node_tt.data :
            #leaf_dict = node2dict(x.bbox,x.name,[])
            full_bbox = bbox_union(full_bbox,x.bbox)
            joint_string += [x.name]
            #children_dict += [leaf_dict]


    for bb in node_tt.branches :
        if bb is None :
            continue
        (sub_bbox,sub_string,sub_dict) = merge_subtree(bb,depth+1,output_dir,coords)
        children_dict+= [sub_dict]
        joint_string += sub_string
        full_bbox = bbox_union(full_bbox,sub_bbox)

    if inputs["meshlab_mode"] == "python" : 
        ms = pymeshlab.MeshSet()
        for ff in joint_string :
            ms.load_new_mesh(ff)
        ms.flatten_visible_layers()
        if not node_tt.isLeafNode :
            ms.simplification_quadric_edge_collapse_decimation(targetfacenum = target_face_num,preserveboundary = True)
        cc = cmap(random.randrange(num_colors))
        #ms.per_face_color_function(r=str(cc[0]*255),g=str(cc[1]*255),b=str(cc[2]*255))
        ms.save_current_mesh(file_name)
        ms.save_current_mesh(output_dir +  node_name_ply)

    node_dict = node2dict(full_bbox,node_name_obj,children_dict,depth,coords)
    print_node(depth,file_name,joint_string)            
    return  (full_bbox,[file_name],node_dict);


def extract_bbox_form_header(full_path) : 
    header_string = os.popen("head -n 5 " + full_path).read()
    return  [float(i) for i in re.split(r'\s{1,}',list(filter(lambda x: "comment bbox" in x , header_string.split("\n")))[0])[2:][:-1]]


def build_3DT(inputs) :
    full_bbox = []
    if False :
        full_bbox = [float(i) for i in inputs["bbox"].split(" ")]
    else : 
        for ff in os.listdir(inputs["input_dir"]):
            if ff.endswith(".ply"):
                full_path = os.path.join(inputs["input_dir"], ff)
                ply_bbox = extract_bbox_form_header(full_path)
                full_bbox = bbox_union(full_bbox,ply_bbox)

        bbox_len = max(full_bbox[1] - full_bbox[0],
                       full_bbox[3] - full_bbox[2],
                       full_bbox[5] - full_bbox[4])
        full_bbox[3] = full_bbox[2] + bbox_len
        full_bbox[5] = full_bbox[4] + bbox_len
    print("full bbox :" + str(full_bbox))

    origin = get_bb_center(full_bbox)
    myTree = Octree(
            full_bbox[1] - full_bbox[0],
            origin,
            max_type="depth",
            max_value=max_depth
    )

            
    for ff in os.listdir(inputs["input_dir"]):
        if ff.endswith(".ply"):
            full_path = os.path.join(inputs["input_dir"], ff)
            bb2 = extract_bbox_form_header(full_path)
            tile_center = get_bb_center(bb2)
            if (len(bb2) == 0) :
                continue;
            do_keep = True
            if inputs["mode"] == "intersect" :
                do_keep = False
            for x in range(2):
                for y in range(2):
                    for z in range(2):
                        pb = [bb2[x],bb2[2+y],bb2[4+z]]
                        is_in = is_inside_bbox(full_bbox,pb)
                        if inputs["mode"] == "intersect" : 
                            do_keep = do_keep or is_in
                        if inputs["mode"] == "strict" :
                            do_keep = do_keep and is_in
            if  do_keep :
                new_ob=tree_obj(full_path,tile_center,bb2)
                myTree.insertNode(tile_center,new_ob)
        

    tile_dic=[]
    tile_output_dir = inputs["output_dir"] + "/"
    Path(tile_output_dir + "/tiles" ).mkdir(parents=True, exist_ok=True)

    if inputs["mode_proj"] == '0' :
        coords  = inputs["coords"]
        inputs["coords"]=''
        x_shift = float(coords.split("x")[0])
        y_shift = float(coords.split("x")[1])
        xmin_new, ymin_new, zmin_new = transformer_glob.transform(x_shift,y_shift,0)
    else :
        xmin_new=0
        ymin_new=0
        zmin_new=0
        
    (sub_bbox,sub_string,sub_dict) = merge_subtree(myTree.root,0,tile_output_dir,inputs["coords"])

    final_dict = {}
    
    sub_dict["transform"]= [
        1,
        0,
        0,
        0,
        0,
        1,
        0,
        0,
        0,
        0,
        1,
        0,
      xmin_new,	
        ymin_new,
        zmin_new,
      1
    ]

    final_dict["asset"] = { "version" : "1.0" }
    final_dict["geometricError"] = "20"
    final_dict["root"] = sub_dict
    json_name = inputs["output_dir"] +  "/tileset_tmp.json"    
    with open(json_name, 'w') as fp:
        json.dump(final_dict, fp)

    
if __name__ == '__main__':
    random.seed(10)
    ###### Input Param parsing / setting =============
    parser = argparse.ArgumentParser(description='conv_ori')
    parser.add_argument('--input_dir', default='',
                        help='give the input ply dir')
    parser.add_argument('--output_dir', default='',
                        help='give the output  dir')
    parser.add_argument('--bbox', default='',
                        help='give bbox')
    parser.add_argument('--meshlab_mode',  default='python',
                        help='print command instead of merging with meshlab')
    parser.add_argument('--mode', default='intersect',
                        help='give the mode, "intersect" => just touchingt, "strict" => the bbox is included')
    parser.add_argument('--coords', default='',
                        help='translation coords')
    parser.add_argument('--mode_proj', default='0',
                        help='mode_proj')
    

    
    args=parser.parse_args()
    inputs=vars(args)


    if  args.input_dir and args.output_dir   :
        inputs["input_dir"]= os.path.expanduser(args.input_dir)
        inputs["output_dir"] = args.output_dir
    else :
        print("Error, ")
        print("--input_dir is mandatory")
        quit()

    if  args.bbox :
        inputs["bbox"]=args.bbox
        inputs["mode"]=args.mode
    if  args.bbox :
        inputs["bbox"]=args.bbox
        inputs["mode"]=args.mode
        inputs["meshlab_mode"]=args.meshlab_mode

    
    print("\n=== Params ===  \n" + "\n".join("{} ==> {}".format(k, v) for k, v in inputs.items()))
    build_3DT(inputs)


