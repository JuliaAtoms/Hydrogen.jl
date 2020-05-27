using Hydrogen
using Documenter

makedocs(;
    modules=[Hydrogen],
    authors="Stefanos Carlström <stefanos.carlstrom@gmail.com> and contributors",
    repo="https://github.com/jagot/Hydrogen.jl/blob/{commit}{path}#L{line}",
    sitename="Hydrogen.jl",
    format=Documenter.HTML(;
                           prettyurls=get(ENV, "CI", "false") == "true",
                           canonical="https://jagot.github.io/Hydrogen.jl",
                           assets=String["assets/latex.js"],
                           mathengine = Documenter.MathJax()
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jagot/Hydrogen.jl",
)
