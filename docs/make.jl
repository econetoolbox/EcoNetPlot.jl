using Documenter
using EcologicalNetworksDynamics
using EcoNetPlot

DocMeta.setdocmeta!(
    EcoNetPostProcessing,
    :DocTestSetup,
    :(using EcoNetPlot, EcologicalNetworksDynamics);
    recursive=true,
)

makedocs(;
    authors="Ismaël Lajaaiti",
    pages=[
        "Home" => "index.md",
        "Functions" => "docstrings.md",
    ],
    sitename="EcoNetPlot.jl",
    repo="https://github.com/econetoolbox/EcoNetPlot.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    modules=[EcoNetPlot],
)

deploydocs(;
    repo="github.com/econetoolbox/EcoNetPlot.jl",
    devbranch="main",
    branch="gh-pages",
)
