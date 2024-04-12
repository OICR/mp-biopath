using Documenter
push!(LOAD_PATH,"../src/")

makedocs(
    sitename = "mp-biopath",
    doctest = false,
    pages = [
        "Home" => "index.md",
        "Introduction" => "introduction.md",
        "Getting Started" => "getting-started.md",
        "Pathways" => "pathways.md",
        "Results" => "results.md",
        "Expression" => "expression.md"
    ]
) 
#include("/Users/tbarbazuk/mp-biopath/src/expression.jl")
