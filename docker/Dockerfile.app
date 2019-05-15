FROM oicr/mpbiopath-env:latest

RUN mkdir /app

COPY . /app/

WORKDIR /app

ENV JULIA_LOAD_PATH=/app/src/:

CMD chmod u+x /app/bin/mp-biopath

CMD ["./bin/mpbiopath --help"]
