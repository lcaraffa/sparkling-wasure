#!/bin/bash


COMPILE_TYPE=Release

CUR_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
MAIN_PROJECT_DIR="$(dirname ${CUR_DIR})"
echo "MAIN_PROJECT_DIR ==> ${MAIN_PROJECT_DIR}"

function build {

    while getopts "d:t:f:j:" OPTION
    do
	case $OPTION in
	    d)
		COMPILE_DIR="${OPTARG}"
		;;
	    t)
		COMPILE_TYPE="${OPTARG}"
		;;
	    f)
		COMPILE_FLAG="${OPTARG}"
		;;
	    j)
		NB_PROC_FLAG="${OPTARG}"
		;;
	esac
    done

    echo ""
    echo "------------------------------------------"
    echo "======> compiling : ${COMPILE_DIR} ..."
    echo "------------------------------------------"
    echo ""

    echo "[$0] Installing/Updating ${COMPILE_DIR} ..."
    echo "compile type = $COMPILE_TYPE"
    mkdir -p ${COMPILE_DIR}/build-$COMPILE_TYPE
    cd ${COMPILE_DIR}/build-$COMPILE_TYPE
    cmake -DMAIN_PROJECT_DIR=${MAIN_PROJECT_DIR} -DCMAKE_BUILD_TYPE=$COMPILE_TYPE ..
    make ${NB_PROC_FLAG}
    make install ${NB_PROC_FLAG}
    rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
    cd -

}

$@
