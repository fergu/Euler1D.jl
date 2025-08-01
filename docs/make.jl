using Documenter
using Euler1D

makedocs(;
    sitename="Euler1D",
    modules=[Euler1D],
    pages=[
        "index.md",
        "Examples.md",
        "Methodology.md",
        "FunctionReference.md"
    ],
    checkdocs=:exports,
)

deploydocs(
    repo = "github.com/fergu/Euler1D.jl.git",
)
