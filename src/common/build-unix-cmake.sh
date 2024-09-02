#!/bin/bash


function build {

    COMPILE_TYPE=Release
    DDT_TRAITS="D2"
    while getopts "d:t:c:f:b:m:" OPTION
    do
	echo " $0 : ${OPTION} ==> ${OPTARG}"
	case $OPTION in
	    d)
		CURRENT_PROJECT_SRC_DIR="${OPTARG}"
		;;
	    c)
		COMPILE_TYPE="${OPTARG}"
		;;
	    f)
		CMAKE_FLAGS="${OPTARG}"
		;;
	    b)
		CURRENT_PROJECT_BUILD_DIR="${OPTARG}"
		;;
	    t)
		DDT_TRAITS="${OPTARG}"
		;;
	    m)
		MAKE_FLAGS="${OPTARG}"
		;;
	esac
    done

    # echo ""
    # echo "------------------------------------------"
    # echo "======> compiling/Installing : ${CURRENT_PROJECT_SRC_DIR} ..."
    # echo "------------------------------------------"
    # echo ""
    if [[ -z ${CURRENT_PROJECT_SRC_DIR} || -z ${CURRENT_PROJECT_BUILD_DIR} ]];
    then
	echo "---- Err in $0 : bad args -----"
	echo "$0 build -d project_build (where CMakelist.txt)"
	exit 1;
    fi

    
    echo "[$0] Installing/Updating ${CURRENT_PROJECT_SRC_DIR} ..."
    echo "compile type = $COMPILE_TYPE"
    mkdir -p ${CURRENT_PROJECT_BUILD_DIR}
    cd ${CURRENT_PROJECT_BUILD_DIR}
    echo "hihi $PWD ===  ${CURRENT_PROJECT_BUILD_DIR} === ${CMAKE_FLAGS} === ${CURRENT_PROJECT_SRC_DIR}"
    cmake  -DCMAKE_BUILD_TYPE=$COMPILE_TYPE ${CMAKE_FLAGS} ${CURRENT_PROJECT_SRC_DIR} -DDDT_TRAITS=${DDT_TRAITS} && make -j 2
    rc=$?;
    echo "rc ===> $rc"
    if [[ $rc != 0 ]]; then return $rc; fi
    make install
    cd -

}

$@
