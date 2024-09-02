#!/bin/bash

function run_fun () {

    ## If inside the docker
    if [ -f /.dockerenv ] ;
    then
	eval ${BASH_CMD}
    else
	if [[ "$(docker images -q $NAME_IMG 2> /dev/null)" == "" ]]; then
	    echo "image $NAME_IMG does not exists... "
	    echo "do ./docker.sh build_compile"
	    exit 1
	else
	    if [ -z "$CONTAINER_NAME" ];
	    then
		CONTAINER_NAME=$(cat /dev/urandom | tr -cd 'a-f0-9' | head -c 8)
	    fi

	    if [ ! -z "$http_proxy" ];
	    then
	    	PROXY_LINE="-e http_proxy='$http_proxy' -e http_proxy_port='$http_proxy_port'  -e http_proxy_ip='$http_proxy_ip'"
	    fi
	    CMD="docker run  -d $MOUNT_CMD -u 0 --cap-add SYS_ADMIN  --privileged --net host $PROXY_LINE -e DDT_MAIN_DIR='$DDT_MAIN_DIR_DOCKER' -e COLUMNS="`tput cols`" -e LINES="`tput lines`" --name $CONTAINER_NAME  -ti ${NAME_IMG}"
	    eval $CMD 
	    container_ip=$(docker inspect $CONTAINER_NAME | grep IPAddress | sed 's/[^0-9.]*//g')
	    echo "ip:$container_ip"
	    echo ""
	    echo "==========================================================="
	    echo "====> DOCKER EXEC :"
	    echo "$CMD"
	    if [ -z "$BASH_CMD" ] && [ -z "${DEBUG_MODE}" ];
	    then
		CMD="docker exec -i -t -u 0 $CONTAINER_NAME /bin/bash"
		echo "==> $CMD"
		eval $CMD
	    else
		case "${DEBUG_MODE}" in
		    bash | shell)
			DOCKER_EXE="docker exec -i -t -u 0 $CONTAINER_NAME /bin/bash"
			;;
		    scala)
			DOCKER_EXE="docker exec -i -t  -u 0  $CONTAINER_NAME  bash -c \"${BASH_CMD}\""    
			;;
		    *)
			DOCKER_EXE="docker exec $CONTAINER_NAME  bash -c \"${BASH_CMD}\""
		esac

		echo "====> DOCKER RUN :"
		echo "$DOCKER_EXE"
		echo ""
		echo "## =======> TYPE THE FOLLOWING BASH CMD TO START <======="
		echo "${RED}$BASH_CMD ${NC}"
		echo ""
		eval $DOCKER_EXE
		rc=$?;
		if [[ $rc != 0 ]];
		then
		    exit $rc;
		else
		    return 0;
		fi
	    fi
	    if [ -z "$DETACHED_TRUE" ];
	    then
		docker rm -f $CONTAINER_NAME
	    else
		echo  " ======>> Container $CONTAINER_NAME detached  "
		echo  " kill it with: docker rm -f $CONTAINER_NAME"
	    fi
	fi
    fi

}


while getopts "l:m:i:c:d:z" OPTION
do
    case $OPTION in
	l)
	    BASH_CMD="${OPTARG}"
	    ;;
	m)
	    MOUNT_CMD="${OPTARG}"
	    ;;
	i)
	    NAME_IMG="${OPTARG}"
	    ;;
	c)
	    CONTAINER_NAME="${OPTARG}"
	    ;;
	z)
	    DETACHED_TRUE="TRUE"
	    ;;
	d)
	    DEBUG_MODE="${OPTARG}"
	    ;;
    esac
    
done


if [[ -z ${NAME_IMG} ]] ;
then
    echo "---- Err : bad args -----"
    echo "$0 -i name_img [-l bash_cmd -m mount_cmd]"
    exit 1;
fi
run_fun

exit 0;
