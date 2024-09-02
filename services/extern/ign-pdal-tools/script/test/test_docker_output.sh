rm -rf tmp/*
mkdir tmp

INPUT_DIR=`pwd`/test/data/
OUPUT_DIR=`pwd`/tmp

docker run -e http_proxy=$http_proxy -e https_proxy=$https_proxy --rm --network host \
-v $INPUT_DIR:/input \
-v $OUPUT_DIR:/output \
lidar_hd/pdal_tools:latest \
python -u -m pdaltools.color \
-i /input/test_data_77050_627755_LA93_IGN69_ground.las \
-o /output/test_docker_capture_output.las \
-r 5 -t 10 --rvb --ir
