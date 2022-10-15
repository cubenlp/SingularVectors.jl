using SingularVectors
using Documenter

DocMeta.setdocmeta!(SingularVectors, :DocTestSetup, :(using SingularVectors); recursive=true)

makedocs(;
    modules=[SingularVectors],
    sitename="SingularVectors.jl"
)

deploydocs(;
    repo="github.com/RexWzh/SingularVectors.jl"
)
