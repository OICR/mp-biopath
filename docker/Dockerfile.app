FROM oicr/mpbiopath-env:0.0.4-SNAPSHOT

RUN mkdir /app

COPY . /app/

WORKDIR /app

CMD ["julia bin/runInference.jl --help"]
