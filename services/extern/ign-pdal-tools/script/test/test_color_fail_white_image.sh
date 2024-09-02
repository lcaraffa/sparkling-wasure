mkdir tmp

echo
echo First test should faild with white image error.
echo

python -m pdaltools.color \
-i ./test/data/sample_lareunion_epsg2975.laz \
-o ./tmp/output.tif \
--rvb -c

echo 
echo Second test should succeed, because we use RGB stream of 20 cm resolution.
echo

python -m pdaltools.color \
-i ./test/data/sample_lareunion_epsg2975.laz \
-o ./tmp/output.tif \
--rvb -c --stream-RGB HR.ORTHOIMAGERY.ORTHOPHOTOS
