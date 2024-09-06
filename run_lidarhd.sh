#!/bin/bash
source ./algo-env.sh

usage() {
  echo "Usage: $0 --input_dir <input_dir>  --output_dir <output_dir> [--params <xml_file_path>] [--colorize]"
  exit 1
}

while [[ "$#" -gt 0 ]]; do
  case $1 in
      --input_dir) INPUT_DIR="$2"; shift ;;      
      --output_dir) OUTPUT_DIR="$2"; shift ;;
      --params) PARAMS="${2%/}"; shift ;;
      --colorize) DO_COLORIZE="TRUE" ;;
    *) echo "Unknown parameter passed: $1"; usage ;;
  esac
  shift
done

[ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ] && { echo "Error." ; usage ; }


if [ -e "${PARAMS}" ]; then
    PARAMS_CMD=" --params ${PARAMS}"
fi
if [ -n "${DO_COLORIZE}" ]; then
    COLORIZE_FLAG=" --colorize"
fi

CMD="${APP_DIR}/run_workflow.sh --input_dir ${INPUT_DIR} --output_dir ${OUTPUT_DIR} ${PARAMS_CMD} ${COLORIZE_FLAG}"
run_cmd_container



 
