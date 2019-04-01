DATA_DIR=../PathwayAnalysis

.PHONY: build-env
build-env:
	docker build -t oicr/mpbiopath-env:1.0.3-SNAPSHOT -f docker/Dockerfile.env .

.PHONY: build-app
build-app:
	docker build -t oicr/mpbiopath:1.0.3-SNAPSHOT -f docker/Dockerfile.app .

.PHONY: run-bash
run-bash:
	docker run -v ${PWD}/${DATA_DIR}:/data -v ${PWD}:/app -it oicr/mpbiopath-env:1.0.3-SNAPSHOT /bin/bash

.PHONY: get-test-expression
get-test-expression:
	curl https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-2836/resources/ExperimentDownloadSupplier.RnaSeqBaseline/tpms.tsv -o test/files/E-MTAB-2836.tsv
	tail -n +5 test/files/E-MTAB-2836.tsv  > test/files/E-MTAB-2836-no-header.tsv

.PHONY: run-tests
run-tests: get-test-expression
	julia test/runtests.jl
