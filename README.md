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

##Essential Genes source
http://ogee.medgenius.info/browse/

##Convert gene names to HUGO

http://useast.ensembl.org/Help/Faq?id=125

#Got list of unique genes with following command

```bash
cut -f 3 essential_9606_all.txt | sort | uniq > essential_9606_all_gene_ids.txt
```
