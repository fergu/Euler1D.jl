using Documenter
using Euler1D

makedocs(;
    sitename="Euler1D",
    modules=[Euler1D],
    format=Documenter.HTML(
        prettyurls=false,
    ),
    checkdocs=:exports,
)
