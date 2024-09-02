# TODO Tester avec un fichier qui existe !
docker run \
-v /var/data/store-lidarhd/developpement/lidarexpress/tests/test0/input_laz:/input \
-v /tmp:/output lidar_hd/pdal_tools \
python -u -m pdaltools.color -i /input/923000_6308000.laz -o /output/923000_6308000.las \
-r 0.1 -t 5 --rvb --ir
