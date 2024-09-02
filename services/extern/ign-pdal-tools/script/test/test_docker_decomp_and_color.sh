rm -rf tmp/*
mkdir tmp/gpao_output_colored

# Todo: Ajouter le fichier en local
INPUT_DIR=`pwd`/../lidarExpress/data/one_micro_laz
OUPUT_DIR=`pwd`/tmp/gpao_output_colored

docker run -e http_proxy=$http_proxy -e https_proxy=$https_proxy --rm --network host \
-v $INPUT_DIR:/input \
-v $OUPUT_DIR:/output \
lidar_hd/pdal_tools \
python -m pdaltools.color \
-i /input/Semis_2021_0785_6378_LA93_IGN69_light.laz \
-o /output/Semis_2021_0785_6378_LA93_IGN69_light.las \
-r 0.1 -t 5 --rvb --ir
