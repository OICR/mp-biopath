DATA_DIR=../PathwayAnalysis
MP_BIOPATH_RELEASE="1.0.4"

.PHONY: build-env
build-env:
	docker build --no-cache -t oicr/mpbiopath-env:latest -f docker/Dockerfile.env .

.PHONY: build-app
build-app:
	docker build --no-cache -t oicr/mpbiopath:latest -f docker/Dockerfile.app .

.PHONY: run-bash
run-bash:
	docker run -v ${PWD}/${DATA_DIR}:/data -v ${PWD}:/app -w /app -it oicr/mpbiopath-env:latest /bin/bash

.PHONY: run-tests
run-tests:
	docker run -v ${PWD}:/app -w /app -it oicr/mpbiopath-env:latest julia tests/run.jl

.PHONY: build-app-release
build-app-release:
	docker build --no-cache -t oicr/mpbiopath:${MP_BIOPATH_RELEASE} -f docker/Dockerfile.app .

.PHONY: run-bash-release
run-bash-release:
	docker run -v ${PWD}/${DATA_DIR}:/data -v ${PWD}:/app -w /app -it oicr/${MP_BIOPATH_RELEASE} /bin/bash

.PHONY: run-tests-release
run-tests-release:
	docker run -v ${PWD}:/app -w /app -it oicr/${MP_BIOPATH_RELEASE} julia tests/run.jl
