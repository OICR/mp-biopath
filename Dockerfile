FROM julia:0.5.1

RUN mkdir /app

COPY . /app/

RUN cd /app && julia bin/INSTALL.jl;

WORKDIR /app

CMD ["julia bin/runInference.jl --help"]
