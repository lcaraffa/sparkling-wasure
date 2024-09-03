import pymeshlab
import argparse
import subprocess
import os


# Check for intersection
def do_boxes_intersect(min1, max1, min2, max2):
    return not (min1[0] > max2[0] or max1[0] < min2[0] or
                min1[1] > max2[1] or max1[1] < min2[1] or
                min1[2] > max2[2] or max1[2] < min2[2])


def bounding_boxes_intersect(bbox1, bary):
    """ Check if two bounding boxes intersect """
    for dim in range(3):  # Check each dimension (x, y, z)
        if bbox1.min()[dim] > bary[dim] or bbox1.max()[dim] < bary[dim]:
            return False
    return True

def transfer_colors(input_mesh, input_pts_dir, output_colorized_dir):
    ms = pymeshlab.MeshSet()
    ms.load_new_mesh(input_mesh)
    bb1 = ms.mesh(0).bounding_box()
    acc = 1
    glob_intersect = False
    for filename in os.listdir(input_pts_dir):
        if filename.endswith("ply"):
            name, extension = os.path.splitext(filename)
            file_path = os.path.join(input_pts_dir, filename)
            ms2 = pymeshlab.MeshSet()
            ms2.load_new_mesh(file_path)        
            bb2 = ms2.mesh(0).bounding_box()
            intersect = do_boxes_intersect(bb1.min(),bb1.max(),bb2.min(),bb2.max())
            if intersect :
                print("[intersect] " + input_mesh + " " + filename)
                ms.load_new_mesh(file_path)
                ms.vertex_attribute_transfer(sourcemesh=acc,targetmesh=0,colortransfer=True,upperbound=pymeshlab.Percentage(10))
                glob_intersect = True
                acc = acc+1
            else :
                print("[         ] " + input_mesh + " " + filename)
            
    if glob_intersect :
        filename=input_mesh.split("/")[-1]
        new_filename = f'{output_colorized_dir}/{filename}'
        ms.set_current_mesh(0)
        ms.save_current_mesh(new_filename)
        command = f"sed -n '3p' {input_mesh} | sed -i '3r /dev/stdin'  {new_filename}"
        subprocess.run(command, shell=True, check=True)

def main():
    parser = argparse.ArgumentParser(description="Transférer les couleurs d'un nuage de points vers un maillage.")
    parser.add_argument('--input_mesh', required=True, help="Chemin vers le fichier du nuage de points avec couleurs (PLY).")
    parser.add_argument('--input_pts_dir', required=True, help="Chemin vers le fichier du maillage sans couleur (PLY).")
    parser.add_argument('--output_colorized_mesh', required=True, help="Chemin vers le fichier de sortie du maillage coloré (PLY).")
    args = parser.parse_args()
    transfer_colors(args.input_mesh, args.input_pts_dir, args.output_colorized_mesh)

if __name__ == "__main__":
    main()
