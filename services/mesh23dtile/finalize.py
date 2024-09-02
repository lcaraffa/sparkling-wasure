import json
import sys
import os
import numpy as np

geom_error = [20,10,5,2,1,0]
def replace_box_values(dd,dpath,depth,list_trans,bbv):
    glob_bbv =  bbv
    list_trans_rec = list_trans.copy()

    if isinstance(dd, dict):
        if 'content' in dd :
            bname = os.path.basename(dd['content']['uri'])
            print(bname)
            fname = os.path.splitext(bname)[0]
            full_path = dpath + "/" + fname + "/tileset.json"
            fo = open(full_path, 'r')
            data_tile = json.load(fo)
            loc_bbox = dd['raw_bbox']['box']

            if depth == 1 :
                dd['transform'] = data_tile['root']['children'][0]['transform']
                dd['boundingVolume'] = data_tile['root']['children'][0]['boundingVolume']
                glob_trans = np.array(data_tile['root']['children'][0]['transform']).reshape([4,4]).T
                glob_bbv = np.array(data_tile['root']['children'][0]['boundingVolume']['box'])
                list_trans_rec.append(glob_trans)
            else :

                local_trans =  np.array(data_tile['root']['children'][0]['transform']).reshape([4,4]).T
                glob_trans_loc = local_trans
                for mat in list_trans_rec :
                    glob_trans_loc = np.dot(glob_trans_loc, np.linalg.inv(mat))
                list_trans_rec.append(glob_trans_loc)
                local_bbv =  np.array(data_tile['root']['children'][0]['boundingVolume']['box'])
                dd['boundingVolume'] = data_tile['root']['children'][0]['boundingVolume']
                dd['transform'] = glob_trans_loc.reshape([4,4]).T.reshape(-1).tolist()
                dd['geometricError'] = geom_error[int(depth/2)]
                dd["refine"] = 'REPLACE'
            dd['content'] = {'uri': fname + "/" + data_tile['root']['children'][0]['content']['uri']}

        for key, value in dd.items():
            if isinstance(value, (dict, list)):
                replace_box_values(value,dpath,depth+1,list_trans_rec,glob_bbv)
    if isinstance(dd, list):
        for item in dd:
            if isinstance(item, (dict, list)):
                replace_box_values(item,dpath,depth+1,list_trans_rec,glob_bbv)

def main(tile_path):
    try:
        with open(tile_path + '/tileset_tmp.json', 'r') as file:
            data = json.load(file)
        replace_box_values(data,tile_path,0,[],[])
        with open(tile_path + '/tileset.json', 'w') as file:
            json.dump(data, file, indent=4)
        
    except FileNotFoundError as err:
        print(err)
    except json.JSONDecodeError:
        print(f"Invalid JSON format in file: {input_file_path}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py dir")
    else:
        main(sys.argv[1])
