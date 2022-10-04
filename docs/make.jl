using SingularVectors
using Documenter

DocMeta.setdocmeta!(SingularVectors, :DocTestSetup, :(using SingularVectors); recursive=true)

makedocs(;
    modules=[SingularVectors],
    authors="rex <1073853456@qq.com> and contributors",
    repo="https://github.com/RexWzh/SingularVectors.jl/blob/{commit}{path}#{line}",
    sitename="SingularVectors.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://RexWzh.github.io/SingularVectors.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/RexWzh/SingularVectors.jl",
    devbranch="main",
)
