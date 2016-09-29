# MP program for performing inference on Reactome Pathways

#Installation

##Gurobi

Got to http://www.gurobi.com/ for instructions on installing Gurobi

http://www.gurobi.com/documentation/6.5/quickstart_linux.pdf

##Julia (http://julialang.org/)

on Linux type

```bash
apt-get install julia
```

##Julia dependencies 

use the Julia console to install Julia external modules

###Jump
https://jump.readthedocs.io/en/latest/index.html
```julia
Pkg.add("JuMP")
```

###Gurobi

```julia
Pkg.add("Gurobi")
```

##Running the package

```bash
julia analyzeBioNetwork.jl <pi filename>
```

