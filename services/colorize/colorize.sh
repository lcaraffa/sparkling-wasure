# Parse command-line options

export CUR_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}")/" && pwd )"
#source ${CUR_DIR}/../../algo-env.sh
source /opt/conda/etc/profile.d/conda.sh

while [[ "$#" -gt 0 ]]; do
  case $1 in
      --input_dir) INPUT_DIR="${2%/}/"; shift ;;
      --output_dir) OUTPUT_DIR="${2%/}/"; shift ;;
    *) echo "Unknown parameter passed: $1"; usage ;;
  esac
  shift
done

mkdir -p ${OUTPUT_DIR}/{colorized_tiles,colorized_ply,colorized_laz}

# If you want to run this script in the current directory
BBOX=$(cat ${OUTPUT_DIR}/wasure_metadata_3d_gen.xml | grep "bbox_ori" | sed 's/<[^>]*>//g')
num_jobs="\j"

conda deactivate
conda activate pdaltools
for lazfile in ${INPUT_DIR}/*.laz ${INPUT_DIR}/*.las; do
    if [ -f "${lazfile}" ]; then
	while (( ${num_jobs@P} >= ${NUM_PROCESS} )); do
            wait -n
	done
	bname=$(basename ${lazfile})
	output_name="${bname%.laz}.las"
	if [ ! -f "${OUTPUT_DIR}/colorized_laz/${output_name}" ]; then
	    python ${CUR_DIR}/../extern/ign-pdal-tools/pdaltools/color.py --input ${lazfile} --output ${OUTPUT_DIR}/colorized_laz/${output_name} --rvb --ir --check-images &
	fi
    fi
done
wait

conda deactivate
conda activate mesh23Dtile
for lazfile in ${OUTPUT_DIR}/colorized_laz/*.laz ${OUTPUT_DIR}/colorized_laz/*.las; do
    if [ -f "${lazfile}" ]; then
	while (( ${num_jobs@P} >= ${NUM_PROCESS} )); do
            wait -n
	done
	bname=$(basename ${lazfile})
	output_name="${bname%.laz}.ply"
	python ${CUR_DIR}/las2ply.py --input_file ${lazfile} --output_file ${OUTPUT_DIR}/colorized_ply/${output_name} --bbox ${BBOX} &	
    fi
done
wait


for meshfile in ${OUTPUT_DIR}/outputs/tiles/*.ply ; do
    while (( ${num_jobs@P} >= ${NUM_PROCESS} )); do
        wait -n
    done
    
    python  ${CUR_DIR}/colorize.py --input_pts_dir ${OUTPUT_DIR}/colorized_ply/ --input_mesh ${meshfile} --output_colorized_mesh ${OUTPUT_DIR}/colorized_tiles &
done


wait
