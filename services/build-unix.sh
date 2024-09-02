#!/bin/bash

libs=("./extern/graphcut" "./extern/LAStools" "./extern/QPBO" "./extern/ann_1.1.2" "./extern/double-conv" "./extern/tinyply" "./ddt"  "./wasure")
CUR_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
COMPILE_TYPE=Release
DDT_TRAITS="D2"

function build {
    while getopts "c:f:j:b:t:" OPTION
    do
	case $OPTION in
	    c)
		COMPILE_TYPE="${OPTARG}"
		;;
	    j)
		NB_PROC_FLAG="${OPTARG}"
		;;
	    b)
		GLOBAL_BUILD_DIR="${OPTARG}"
		;;
	    t)
		DDT_TRAITS="${OPTARG}"
		;;
	esac
    done

    if [[ -z ${GLOBAL_BUILD_DIR} ]];
    then
	echo "---- Err in $0 : bad args -----"
	echo "GLOBAL_BUILD_DIR not set"
	exit 1

    fi

    export EXTERN_PROJECT_SRC_DIR=${CUR_DIR}/extern/
    export CMAKE_FIND_DIR=${CUR_DIR}/cmake/
    export GLOBAL_BUILD_DIR=${GLOBAL_BUILD_DIR}
    export GLOBAL_LIBS_DIR=${GLOBAL_BUILD_DIR}/libs/
    export GLOBAL_EXE_DIR=${GLOBAL_BUILD_DIR}/bin/


    echo ""
    echo "------------------------------------------"
    echo "globs var"
    echo "EXTERN_PROJECT_SRC_DIR:$EXTERN_PROJECT_SRC_DIR"
    echo "CMAKE_FIND_DIR:$CMAKE_FIND_DIR"
    echo "GLOBAL_LIBS_DIR:$GLOBAL_LIBS_DIR"
    echo "GLOBAL_EXE_DIR:$GLOBAL_EXE_DIR"
    

    
    mkdir -p ${GLOBAL_LIBS_DIR}
    mkdir -p ${GLOBAL_EXE_DIR}
    for i in ${!libs[@]}; do
	CURRENT_PROJECT_SRC_DIR=${CUR_DIR}${libs[i]}
	CURRENT_PROJECT_BUILD_DIR=${GLOBAL_BUILD_DIR}${libs[i]}
	echo ""
	 echo "$(tput setaf 7)$(tput setaf 1) ======> compiling/Installing : ${libs[i]} ..."
	 echo "CURRENT_PROJECT_SRC_DIR:${CUR_DIR}${libs[i]}"
	 echo "CURRENT_PROJECT_BUILD_DIR:${GLOBAL_BUILD_DIR}${libs[i]}"
	 echo "------------------------------------------"


	 mkdir -p ${CURRENT_PROJECT_BUILD_DIR}
	if [[ -f ${CURRENT_PROJECT_SRC_DIR}/build-unix.sh ]];
	then
	    cd ${CURRENT_PROJECT_SRC_DIR};
	    ./build-unix.sh build -b ${CURRENT_PROJECT_BUILD_DIR}
	    cd -
	elif [[ -f ${CURRENT_PROJECT_SRC_DIR}/CMakeLists.txt ]];
	then
	    ${DDT_MAIN_DIR}/src/common/build-unix-cmake.sh build -d ${CURRENT_PROJECT_SRC_DIR} -b ${CURRENT_PROJECT_BUILD_DIR} -c ${COMPILE_TYPE} -t ${DDT_TRAITS};
	    rc=$?; if [[ $rc != 0 ]]; then return $rc; fi
	else
	    echo "error no build-unix.sh or CMakelists.txt found in ${CURRENT_PROJECT_SRC_DIR} "
	    exit 1
	fi
    done
}

$@

