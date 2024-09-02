echo "start $0"
echo "----------------"
export SCRIPT_DIR="$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd  )/"
#source ${SCRIPT_DIR}/../../algo-env.sh

function start_spark_master {
    if [[ ! -z ${1} ]];
    then
        echo ""
        echo "====================== SPARK MASTER ======================="
        echo "Try to start master on : ${1}"
        echo "spark monitoring - http://${1}:8080/"
        /usr/local/bin/spark-3.5.0-bin-hadoop3-scala2.13/sbin/start-master.sh
    else
        echo "SPARK_MASTER_HOST not set "
        exit 1;
    fi
}

function start_spark_slave {
    if [ -z "$1" ]
    then
        echo "No master ip supplied : do"
        echo "$0 ip.of.the.master"
    fi
    echo "======================= SPARK SLAVE ============================"
    echo "Start slave on : ${1}"
    /usr/local/bin/spark-3.5.0-bin-hadoop3-scala2.13/sbin/start-slave.sh spark://${1}:7077 --work-dir ${TMP_DIR}
}

function start_spark_history {
    echo "======================= SPARK History ============================"
    echo "Start history server"
    /usr/local/bin/spark-3.5.0-bin-hadoop3-scala2.13/sbin/start-history-server.sh 
}

function start_all {
    start_spark_master ${1}
    start_spark_slave ${1}
    start_spark_history
}

function stop_all {
     /usr/local/bin/spark-3.5.0-bin-hadoop3-scala2.13/sbin/stop-all.sh
}

if [ $# -eq 0 ]; then
    echo "No arguments provided => start master and slave"
    start_all
    exit 0
else
    $@
fi
