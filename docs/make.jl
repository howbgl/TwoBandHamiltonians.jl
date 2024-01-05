using TwoBandHamiltonians
using Documenter
using Unitful

DocMeta.setdocmeta!(TwoBandHamiltonians, :DocTestSetup, :(using TwoBandHamiltonians); recursive=true)

makedocs(;
    modules=[TwoBandHamiltonians],
    authors="Wolfgang Hogger",
    repo="https://github.com/howbgl/TwoBandHamiltonians.jl/blob/{commit}{path}#{line}",
    sitename="TwoBandHamiltonians.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://howbgl.github.io/TwoBandHamiltonians.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Reference" => "reference.md"
    ],
)

deploydocs(;
    repo="github.com/howbgl/TwoBandHamiltonians.jl",
    devbranch="main",
)
