if ! [ -f /.dockerenv ]; then
    echo "Error : should be run inside docker";
    exit 1
fi
cd ${DDT_MAIN_DIR}

INPUT_SCRIPT=${DDT_MAIN_DIR}/src/scala/ddt_stream.scala
#source ${DDT_MAIN_DIR}/algo-env.sh


function export_params {
    echo "export SPARK_MASTER_HOST=${1}" >> /usr/local/bin/spark-3.5.0-bin-hadoop3-scala2.13/conf/spark-env.sh
    echo "export LOCAL_DIRS=${SPARK_TMP_DIR}" >> /usr/local/bin/spark-3.5.0-bin-hadoop3-scala2.13/conf/spark-env.sh
    echo "export SPARK_LOCAL_DIRS=${SPARK_TMP_DIR}" >> /usr/local/bin/spark-3.5.0-bin-hadoop3-scala2.13/conf/spark-env.sh
    echo "export SPARK_EXECUTOR_MEMORY=${SPARK_EXECUTOR_MEMORY}" >> /usr/local/bin/spark-3.5.0-bin-hadoop3-scala2.13/conf/spark-env.sh
    echo "export SPARK_DRIVER_MEMORY=${SPARK_DRIVER_MEMORY}" >> /usr/local/bin/spark-3.5.0-bin-hadoop3-scala2.13/conf/spark-env.sh
    echo "export SPARK_WORKER_CORES=${SPARK_WORKER_CORES}" >> /usr/local/bin/spark-3.5.0-bin-hadoop3-scala2.13/conf/spark-env.sh
    echo "export SPARK_WORKER_MEMORY=${SPARK_WORKER_MEMORY}" >> /usr/local/bin/spark-3.5.0-bin-hadoop3-scala2.13/conf/spark-env.sh
    echo "export YARN_LOG_DIR=${SPARK_TMP_DIR}" >> /usr/local/bin/spark-3.5.0-bin-hadoop3-scala2.13/conf/spark-env.sh
    echo "======== PARAMS ========="
    cat /usr/local/bin/spark-3.5.0-bin-hadoop3-scala2.13/conf/spark-env.sh
    echo "======== PARAMS ========="
}

function run_local (){
    ${DDT_MAIN_DIR}/src/spark/spark.sh start_all ${MASTER_IP}
    /usr/local/bin/spark-3.5.0-bin-hadoop3-scala2.13/bin/spark-shell \
        -i ${INPUT_SCRIPT}  \
        --jars ${GLOBAL_BUILD_DIR}/spark/target/scala-2.13/iqlib-spark_2.13-1.0.jar  \
        --master spark://localhost:7077  -Dspark.executor.memory=1g -Dspark.driver.memory=1g
}


function run_master (){

    export_params ${MASTER_IP}


    ${DDT_MAIN_DIR}/src/spark/spark.sh start_all ${MASTER_IP}


    echo ""
    echo "//  ===============   INFO ======================="
    echo ""


    CMD="/usr/local/bin/spark-3.5.0-bin-hadoop3-scala2.13/bin/spark-shell  \
					 --jars ${GLOBAL_BUILD_DIR}/spark/target/scala-2.13/iqlib-spark_2.13-1.0.jar  \
					 --master spark://${MASTER_IP}:7077  \
					 --conf \"spark.local.dirs=${SPARK_TMP_DIR}\"  \
					 --conf \"spark.rdd.compress=true\"  \
					 --conf \"spark.eventLog.enabled=true\"  \
					 --conf \"spark.memory.fraction=0.2\" \
					 --conf \"spark.driver.allowMultipleContexts=true\" \
					 --conf \"spark.memory.storageFraction=0.8\" \
					 --conf \"spark.worker.cleanup.enabled=true\" \
					 --conf \"spark.worker.cleanup.interva=350\" \
					 --conf \"spark.memory.offHeap.enabled=true\" \
					 --conf \"spark.memory.offHeap.size=10g\" \
					 --conf \"spark.network.timeout=10000000\" \
					 --conf \"spark.history.fs.logDirectory=${SPARK_HISTORY_DIR}\" \
					 --conf \"spark.serializer=org.apache.spark.serializer.KryoSerializer\"  \
					 --conf \"yarn.nodemanager.log-dirs=${SPARK_HISTORY_DIR}\" \
					 -Djava.io.tmpdir=\"${SPARK_TMP_DIR}\" \
					 -Dspark.executor.memory=${SPARK_EXECUTOR_MEMORY} \
					 -Dspark.driver.memory=${SPARK_DRIVER_MEMORY} "


    if [ "$DEBUG_MODE" = true ] ; then
	echo ""
	echo "//  ===> to run the script, do:"
	echo ":load ${INPUT_SCRIPT}"
	echo ""
	echo ""
	eval ${CMD}
    else
     	echo ":load ${INPUT_SCRIPT}" | eval ${CMD}
    fi

	}

function run_slave (){
    export_params ${MASTER_IP}
    ${DDT_MAIN_DIR}/src/spark/spark.sh start_spark_slave ${MASTER_IP}
    ACC=0
    while true; do
	echo "slave runing on master : ${MASTER_IP}"
	MEMORY=$(du -h -d 1 ${SPARK_TMP_DIR} 2> /dev/null  |  tail -n 1 | awk '{print $1;}') 
	echo "loop:$ACC memory ${SPARK_TMP_DIR} :$MEMORY"
	sleep 5;
	ACC=$((ACC+1))
	sleep 2s
    done
}


function help (){
    echo "$0 -i -o"
}

while getopts "i:o:p:s:f:m:b:r:c:d" OPTION
do
    case $OPTION in
        i)
            export INPUT_DATA_DIR="${OPTARG}"
            ;;
        o)
            export OUTPUT_DATA_DIR="${OPTARG}"
            ;;
	p)
            export PARAM_PATH="${OPTARG}"
            ;;
        s)
            SPARK_CONF="${OPTARG,,}"
            ;;
	f)
            INPUT_SCRIPT="${OPTARG}"
            ;;
	m)
            MASTER_IP="${OPTARG}"
            ;;
	c)
	    SPARK_WORKER_CORES="${OPTARG}"
            ;;
	d)
	    DEBUG_MODE=true
	    ;;
	b)
	    export GLOBAL_BUILD_DIR="${OPTARG}"
	    ;;
	r)
            export ALGO_SEED="${OPTARG}"
            ;;
    esac

done

if [[ ! -d $INPUT_DATA_DIR || ! -d $OUTPUT_DATA_DIR  || -z $PARAM_PATH ]];
then
    echo "---- Err : bad args -----"
    echo "INPUT_DATA_DIR=\"${INPUT_DATA_DIR}\" or OUTPUT_DATA_DIR=\"${OUTPUT_DATA_DIR}\" or PARAM_PATH=\"${PARAM_PATH}\" does not exists" 
    exit 1;
fi



if [ -z "$SPARK_CONF" ]
then
    SPARK_CONF="local"
fi



echo "SPARK_CONF ---- $SPARK_CONF"
echo "MASTER IP ----- $MASTER_IP"



case "$SPARK_CONF" in
    "local")
	MASTER_IP="localhost"
	run_local
	;;
    "master")
	if [[ -z ${MASTER_IP} ]] ;
        then
            echo "---- Error -----"
            echo "if you are master you should give your ip \"-m\""
            exit 1;
        fi
	echo "===================="

	run_master
        ;;
    "slave")
        if [[ -z ${MASTER_IP} ]] ;
        then
            echo "---- Error -----"
            echo "if you are slave you should give the master ip with \"-m\""
            exit 1;
        fi
        run_slave
        ;;
    *)
        echo "error, spark_conf = [local,master,slave]"
        exit 1
        ;;
esac





