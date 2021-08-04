DATA_DIR=../PathwayAnalysis

.PHONY: build-env
build-env:
	docker build --no-cache -t oicr/mpbiopath-env:1.0.4 -f docker/Dockerfile.env .

.PHONY: build-app
build-app:
	docker build --no-cache -t oicr/mpbiopath:1.0.4 -f docker/Dockerfile.app .

.PHONY: run-bash
run-bash:
	docker run -v ${PWD}/${DATA_DIR}:/data -v ${PWD}:/app -w /app -it oicr/mpbiopath-env:1.0.4 /bin/bash

.PHONY: get-test-expression
get-test-expression:
	curl https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-2836/resources/ExperimentDownloadSupplier.RnaSeqBaseline/tpms.tsv -o tests/files/E-MTAB-2836.tsv

.PHONY: run-tests
run-tests:
	docker run -v ${PWD}:/app -w /app -it oicr/mpbiopath-env:1.0.4 julia tests/run.jl
