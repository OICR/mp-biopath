FROM mpbiopath-env:0.0.1

RUN mkdir /app

COPY . /app/

WORKDIR /app

CMD ["julia bin/runInference.jl --help"]
