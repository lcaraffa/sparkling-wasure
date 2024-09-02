parallel (

	"conda":{
	node('linux_conda') {

		stage('init') {
		gitlabCommitStatus("init") {
			checkout scm
		}
		}

		stage('mamba') {
		gitlabCommitStatus("mamba") {
			sh "make mamba-env-update"
		}
		}

		stage('test') {
		gitlabCommitStatus("test") {
			sh "mamba run -n pdaltools make testing"
		}
		}

		stage('build') {
		gitlabCommitStatus("build") {
			if (env.BRANCH_NAME == 'master') {
				sh "mamba run -n pdaltools make build"
			} else {
				echo "Nothing to do, because branch is not master"
			}
		}
		}

		stage('deploy') {
		gitlabCommitStatus("deploy") {
			if (env.BRANCH_NAME == 'master') {
				sh "mamba run -n pdaltools make check"
				withCredentials([usernamePassword(credentialsId: 'pypi', usernameVariable: 'USERNAME', passwordVariable: 'PASSWORD')]) {
					sh "mamba run -n pdaltools twine upload dist/* -u ${USERNAME} -p ${PASSWORD}"
				}
			} else {
				echo "Nothing to do, because branch is not master"
			}
	
		}
		}

	}
	},

	"docker":{	
	node('DOCKER') {

		stage('build-docker-image') {
			gitlabCommitStatus("build-docker-image") {
				if (env.BRANCH_NAME == 'master') {
					checkout scm
					sh "make docker-build"
				} else {
					echo "Nothing to do, because branch is not master"
				}
			}
		}

		stage('test-docker-images') {
			gitlabCommitStatus("test-docker-images") {
				if (env.BRANCH_NAME == 'master') {
					sh "./script/test/test_docker_output.sh"
				} else {
					echo "Nothing to do, because branch is not master"
				}
			}
		}

        stage('deploy-docker-image') {
			gitlabCommitStatus("build-docker-image") {

				if (env.BRANCH_NAME == 'master') {
					withCredentials([string(credentialsId: 'svc_lidarhd', variable: 'svc_lidarhd')]) {
						sh "docker/deploy-jenkins.sh ${svc_lidarhd}"
					}
				} else {
					echo "Nothing to do, because branch is not master"
				}

			}
		}

	}
	}
)