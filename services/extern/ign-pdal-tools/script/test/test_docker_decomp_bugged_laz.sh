rm -rf tmp/*
mkdir tmp/gpao_output_colored

# Todo: Ajouter le fichier dans GIT avec Git LFS
INPUT_DIR=/media/data/Bug_ouverture_laz/one
OUPUT_DIR=`pwd`/tmp/gpao_output_colored

docker run -e http_proxy=$http_proxy -e https_proxy=$https_proxy --rm --network host \
-v $INPUT_DIR:/input \
-v $OUPUT_DIR:/output \
lidar_hd/pdal_tools \
python -u -m pdaltools.color \
-i /input/436000_6469000.laz \
-o /output/436000_6469000.las \
-r 0.1 -t 5 --rvb --ir

