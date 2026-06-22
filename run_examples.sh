#!/bin/bash
source ./algo-env.sh

export COLORIZE_FLAG="--colorize"

# Example
INPUT_DIR=${PWD}/datas/maison/
PARAMS=${INPUT_DIR}/wasure_metadata.xml
OUTPUT_DIR=${PWD}/output_example/maison/
CMD="${APP_DIR}/run_workflow.sh --input_dir ${INPUT_DIR} --output_dir ${OUTPUT_DIR} ${COLORIZE_FLAG}"
run_cmd_container
