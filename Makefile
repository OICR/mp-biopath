.PHONY: build-env
build-env:
	docker build -t oicr/mpbiopath-env:0.0.4 -f docker/Dockerfile.env .
