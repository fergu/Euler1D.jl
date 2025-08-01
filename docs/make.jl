using Documenter
using Euler1D

makedocs(;
    sitename="Euler1D",
    modules=[Euler1D],
    format=Documenter.HTML(
        prettyurls=false,
    ),
    pages=[
        "index.md",
        "Examples.md",
        "Methodology.md",
        "FunctionReference.md"
    ],
    checkdocs=:exports,
)
