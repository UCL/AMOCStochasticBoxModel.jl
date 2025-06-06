using AMOCStochasticBoxModel
using Documenter

DocMeta.setdocmeta!(AMOCStochasticBoxModel, :DocTestSetup, :(using AMOCStochasticBoxModel); recursive=true)

makedocs(;
    modules=[AMOCStochasticBoxModel],
    authors="Matt Graham and contributors",
    sitename="AMOCStochasticBoxModel.jl",
    format=Documenter.HTML(;
        canonical="https://github-pages.ucl.ac.uk/AMOCStochasticBoxModel.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/UCL/AMOCStochasticBoxModel.jl",
    devbranch="main",
)
