# Makefile to manage main tasks
# cf. https://blog.ianpreston.ca/conda/python/bash/2020/05/13/conda_envs.html#makefile

# Oneshell means I can run multiple lines in a recipe in the same shell, so I don't have to
# chain commands together with semicolon
.ONESHELL:
SHELL = /bin/bash


##############################
# Install
##############################
install:
	mamba env update -n pdaltools -f environment.yml


##############################
# Dev/Contrib tools
##############################

testing:
	python -m pytest ./test -s --log-cli-level DEBUG -m "not geopf"

testing_full:
	python -m pytest ./test -s --log-cli-level DEBUG

install-precommit:
	pre-commit install


##############################
# Build/deploy pip lib
##############################

deploy: check
	twine upload dist/*

check: dist/ign-pdal-tool*.tar.gz
	twine check dist/*

dist/ign-pdal-tool*.tar.gz:
	python -m build

build: clean
	python -m build

clean:
	rm -rf tmp
	rm -rf ign_pdal_tools.egg-info
	rm -rf dist

##############################
# Build/deploy Docker image
##############################

REGISTRY=ghcr.io
NAMESPACE=ignf
IMAGE_NAME=ign-pdal-tools
VERSION=`python -m pdaltools._version`
FULL_IMAGE_NAME=${REGISTRY}/${NAMESPACE}/${IMAGE_NAME}:${VERSION}

docker-build: clean
	docker build --no-cache -t ${IMAGE_NAME}:${VERSION} -f Dockerfile .

docker-test:
	docker run --rm -it ${IMAGE_NAME}:${VERSION} python -m pytest -s

docker-ssh:
	docker run --rm -it ${IMAGE_NAME}:${VERSION} 

docker-remove:
	docker rmi -f `docker images | grep ${IMAGE_NAME}:${VERSION} | tr -s ' ' | cut -d ' ' -f 3`

docker-deploy:
	docker tag ${IMAGE_NAME}:${VERSION} ${FULL_IMAGE_NAME}
	docker push ${FULL_IMAGE_NAME}
