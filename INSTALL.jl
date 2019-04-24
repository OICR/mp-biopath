import Pkg
Pkg.update()

Pkg.add("JuMP")
#Pkg.add("CoinOptServices")
##Pkg.build("CoinOptServices")
Pkg.add("ArgParse")
Pkg.add("Nullables")
#
Pkg.add("Ipopt")
##Pkg.add("Cbc")
##Pkg.add("NLopt")
#
##Pkg.add("AmplNLWriter")
#Pkg.add("Requests")
#Pkg.add("FactCheck")
Pkg.add("CSV")
Pkg.add("YAML")
Pkg.add("JSON")

#Pkg.add("CategoricalArrays")
Pkg.add("DocOpt")
#Pkg.add("DataTables")
Pkg.add("DataFrames")

##for tsne
##Pkg.add("Gadfly")
##Pkg.add("Cairo")
##Pkg.add("Fontconfig")
