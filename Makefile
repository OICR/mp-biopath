.PHONY: build-env
build-env:
	docker build -t oicr/mpbiopath-env:0.0.4-SNAPSHOT -f docker/Dockerfile.env .

.PHONY: build-app
build-app:
	docker build -t oicr/mpbiopath:0.0.4-SNAPSHOT -f docker/Dockerfile.app .
