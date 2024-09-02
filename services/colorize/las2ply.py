import argparse
import laspy
import os
import numpy as np
import pymeshlab
# Ouvrez le fichier LAS en lecture


def process(inputs) :

        bbox = [s.split("x") for s in inputs["bbox"].split(":")]
        inFile = laspy.read(inputs["input_file"])
        outFile = open(inputs["output_file"], "w")

        outFile.write("ply\n")
        outFile.write("format ascii 1.0\n")
        outFile.write("element vertex {}\n".format(len(inFile)))
        outFile.write("property float x\n")
        outFile.write("property float y\n")
        outFile.write("property float z\n")
        outFile.write("property uchar red\n")
        outFile.write("property uchar green\n")
        outFile.write("property uchar blue\n")
        outFile.write("end_header\n")


        vv1 = np.vstack([inFile.x,inFile.y,inFile.z]).T
        vv2 = np.vstack([(inFile.red/256).astype(int) ,(inFile.green/256).astype(int),(inFile.blue/256).astype(int)]).T
        for jj in range(0,len(vv1)) :
                pp1 = vv1[jj]
                pp2 = vv2[jj]
                #
                for dd in range(0,3) :
                        #import pdb; pdb.set_trace()
                        outFile.write(str(pp1[dd] - int(float(bbox[dd][0])))+ " ")
                outFile.write(" ".join(str(ii) for ii in pp2)+ "\n")
        outFile.close()





if __name__ == '__main__':

    ###### Input Param parsing / setting =============
    parser = argparse.ArgumentParser(description='conv_ori')
    parser.add_argument('--input_file', default='',
                        help='inpuy laz')
    parser.add_argument('--output_file', default='',
                        help='output ply')
    parser.add_argument('--bbox', default='',
                        help='give bbox')
    parser.add_argument('--meshlab_mode',  default='python',
                        help='print command instead of merging with meshlab')
    
    
    args=parser.parse_args()
    inputs=vars(args)


    if  args.input_file and args.output_file   :
        inputs["input_file"]= os.path.expanduser(args.input_file)
        inputs["output_file"] = args.output_file
    else :
        print("Error, ")
        print("--input_file is mandatory")
        quit()

    if  args.bbox :
        inputs["bbox"]=args.bbox



        
        
    print("\n=== Params ===  \n" + "\n".join("{} ==> {}".format(k, v) for k, v in inputs.items()))
    process(inputs)



