using Documenter
using HaloArrays

makedocs(;
    modules=[HaloArrays],
    sitename="HaloArrays.jl",
    authors="Davide Lasagna and contributors",
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Getting started" => "manual/getting_started.md",
            "Halo indexing" => "manual/indexing.md",
            "Halo exchange" => "manual/exchange.md",
            "Internals" => "manual/internals.md",
        ],
        "API reference" => "api.md",
    ],
    checkdocs=:none,
)

deploydocs(;
    repo="github.com/Davide-Lasagna-s-Lab/HaloArrays.jl.git",
    devbranch="master",
)
