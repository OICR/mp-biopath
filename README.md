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

###Gurobi

```julia
Pkg.add("Gurobi")
```

##Running the package

```bash
julia bin/analyzeBioNetwork.jl --help
```

##Essential Genes source
http://ogee.medgenius.info/browse/

##Convert gene names to HUGO

http://useast.ensembl.org/Help/Faq?id=125

Used biomart central to download a list of HUGO gene symbols to match the ensemble geneIDs. The exact query can be found in the script "get_gene_id_mapping_from_biomart_central.pl" that was generated by Biomart Central.
 
#Got list of unique genes with following command

```bash
cut -f 3 essential_9606_all.txt | sort | uniq > essential_9606_all_gene_ids.txt
```
##Synthetically lethal 

The two files sl_human and sdl_human were downloaded from: http://histone.sce.ntu.edu.sg/SynLethDB/downloadPage.php

And more information about the data can be found in the following articles:
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4445436/ "The authors also analyzed synthetic dosage lethality (SDL), which is the overexpression of a gene that leads to the essentiality of a second gene, making a SDL pair"
- http://nar.oxfordjournals.org/content/early/2015/10/29/nar.gkv1108.full 

##Coexpression Data

http://coxpresdb.jp/download.shtml

##Coespression DB API

http://coxpresdb.jp/help/API.shtml

##test 

```bash
julia bin/analyzeBioNetwork.jl ~/git/PGM/Pathways\ where\ negative\ OR\ is\ AND/PIP3_activates_AKT_signaling_Sept23_2016_sorted_checked_patched_1_NegativeORtoAND.txt ~/git/PGM/InputDataForTesting/Mock_Short_Input_Data_for_PIP3_Activates_AKT_Signaling_Text.txt --key-outputs
```bash
