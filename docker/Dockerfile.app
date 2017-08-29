FROM oicr/mpbiopath-env:0.0.4-SNAPSHOT

RUN mkdir /app

COPY . /app/

WORKDIR /app

CMD chmod u+x /app/bin/mp-biopath

CMD ["./bin/mpbiopath --help"]
