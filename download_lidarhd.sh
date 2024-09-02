#!/bin/bash
source ./algo-env.sh

usage() {
  echo "Usage: $0 --list_files <txt_file>  --project_path <project_path>"
  exit 1
}

while [[ "$#" -gt 0 ]]; do
  case $1 in
    --list_files) LIST_FILES="$2"; shift ;;
    --project_path) PROJECT_PATH="$2"; shift ;;
    *) echo "Unknown parameter passed: $1"; usage ;;
  esac
  shift
done

[ -z "$LIST_FILES" ] || [ -z "$PROJECT_PATH" ] && { echo "Error." ; usage ; }

mkdir -p ${PROJECT_PATH}/{inputs,outputs}

while IFS= read -r line; do
    FILENAME=$(basename "${line}")
    FILEPATH=${PROJECT_PATH}/inputs/${FILENAME}
    if [ ! -e "${FILEPATH}" ]; then
	wget -O ${FILEPATH}  ${line}
    fi
done < "${LIST_FILES}"




 
