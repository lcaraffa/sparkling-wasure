import convertcloud as cvc
import numpy as np
import argparse
import os
import io

def convert_ply(inputs) :
    list_name = []
    cur_dir = os.path.basename(os.path.normpath(inputs["input_dir"])) 
    
    for ff in os.listdir(inputs["input_dir"]):
        if ff.endswith(".npts"):
            full_path = os.path.join(inputs["input_dir"], ff)
            #bname=os.path.basename(ff)
            bname=os.path.splitext(ff)[0]
            print(bname)
            output_path_dir=inputs["output_dir"] + "/" + bname 
            output_path_file=output_path_dir + "/" + bname + ".ply" 
            if not os.path.isdir(output_path_dir) : 
                os.mkdir(output_path_dir)
            conv = cvc.Converter()
            conv.load_points(full_path)
            conv.set_ori_dist(5);
            conv.convert(output_path_file)




if __name__ == '__main__':
    ###### Input Param parsing / setting =============
    parser = argparse.ArgumentParser(description='conv_ori')
    parser.add_argument('--input_dir', default='',
                        help='give the input ply dir')
    parser.add_argument('--output_dir', default='',
                        help='give the output  dir')

    args=parser.parse_args()
    inputs=vars(args)

    
    if  args.input_dir   :
        inputs["input_dir"]= os.path.expanduser(args.input_dir)
    else :
        print("Error, ")
        print("--input_dir is mandatory")
        quit()

    if  args.output_dir :
        inputs["output_dir"] = args.output_dir
    else :
        inputs["output_dir"] = args.input_dir

    
    print("\n=== Params ===  \n" + "\n".join("{} ==> {}".format(k, v) for k, v in inputs.items()))
    convert_ply(inputs)

            




# print(str(is_inside_bbox(bb1,[-29,-29,29])))


## "<bbox>3800x4800:4500x5500:0x1000</bbox>"
# 4000 4500 4750 5250 0x1000
