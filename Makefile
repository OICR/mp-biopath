.PHONY: build-env
build-env:
	docker build -t oicr/mpbiopath-env:1.0.0-SNAPSHOT -f docker/Dockerfile.env .

.PHONY: build-app
build-app:
	docker build -t oicr/mpbiopath:1.0.0-SNAPSHOT -f docker/Dockerfile.app .
