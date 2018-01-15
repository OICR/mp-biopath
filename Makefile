DATA_DIR=../PathwayAnalysis

.PHONY: build-env
build-env:
	docker build -t oicr/mpbiopath-env:1.0.1-SNAPSHOT -f docker/Dockerfile.env .

.PHONY: build-app
build-app:
	docker build -t oicr/mpbiopath:1.0.1 -f docker/Dockerfile.app .

.PHONY: run-bash
run-bash:
	docker run -v ${PWD}/${DATA_DIR}:/data -v ${PWD}:/app -it oicr/mpbiopath-env:1.0.1-SNAPSHOT /bin/bash
