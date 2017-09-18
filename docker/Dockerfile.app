FROM oicr/mpbiopath-env:0.0.4-SNAPSHOT

RUN mkdir /app

COPY . /app/

WORKDIR /app

ENV JULIA_LOAD_PATH="/app/src/"

CMD chmod u+x /app/bin/mp-biopath

CMD ["./bin/mpbiopath --help"]
