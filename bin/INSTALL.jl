Pkg.add("JuMP")
Pkg.add("ArgParse")
Pkg.add("CoinOptServices")
Pkg.add("AmplNLWriter")
Pkg.add("DataFrames")
Pkg.add("DataArrays")
Pkg.add("ExcelReaders")
Pkg.add("Conda")
ENV["PYTHON"]=""
Pkg.build("PyCall")

Pkg.add("DataFrames")
