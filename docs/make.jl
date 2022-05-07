using SingularVectors
using Documenter

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true"
)

makedocs(;
    modules=[SingularVectors],
    sitename="SingularVectors.jl",
    format=format,
    pages=[
        "Home" => "index.md",
        "Reference" => "reference.md",
    ],
)

deploydocs(;
    repo="github.com/RexWzh/SingularVectors.jl",
    devbranch="master",
)
