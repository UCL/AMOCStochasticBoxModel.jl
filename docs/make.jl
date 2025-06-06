using AMOCStochasticBoxModel
using Documenter

DocMeta.setdocmeta!(AMOCStochasticBoxModel, :DocTestSetup, :(using AMOCStochasticBoxModel); recursive=true)

makedocs(;
    modules=[AMOCStochasticBoxModel],
    authors="Matt Graham <matthew.m.graham@gmail.com> and contributors",
    sitename="AMOCStochasticBoxModel.jl",
    format=Documenter.HTML(;
        canonical="https://matt-graham.github.io/AMOCStochasticBoxModel.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/matt-graham/AMOCStochasticBoxModel.jl",
    devbranch="main",
)
