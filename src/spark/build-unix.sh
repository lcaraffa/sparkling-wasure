echo "==================="
echo "run  $0 ... " 

PWD=$(pwd)
if [ ! -z "$http_proxy" ];
then
    export JAVA_OPTS="$JAVA_OPTS -Dhttp.proxySet=true -Dhttp.proxyHost=$http_proxy_ip -Dhttp.proxyPort=$http_proxy_port -Dhttps.proxySet=true -Dhttps.proxyHost=$http_proxy_ip -Dhttps.proxyPort=$http_proxy_port"
    echo "JAVA_OPTS  ==> $JAVA_OPTS"
else
    echo " ============= WARNING if you are behind a proxy =========="
    echo "http_proxy not set, please set it if you are behind a proxy"
    echo "or sbt may fail"
    echo " ============= WARNING if you are behind a proxy =========="
fi
/usr/local/bin/sbt/bin/sbt package  -ivy /tmp/ivy2 -java-home /usr/lib/jvm/java-8-openjdk-amd64/
rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
exit 0;

