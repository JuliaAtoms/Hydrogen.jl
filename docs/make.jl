using Hydrogen
using Documenter

DocMeta.setdocmeta!(Hydrogen, :DocTestSetup, :(using Hydrogen, AtomicLevels, HalfIntegers);
                    recursive=true)
makedocs(;
    modules=[Hydrogen],
    authors="Stefanos Carlström <stefanos.carlstrom@gmail.com> and contributors",
    repo="https://github.com/JuliaAtoms/Hydrogen.jl/blob/{commit}{path}#L{line}",
    sitename="Hydrogen.jl",
    format=Documenter.HTML(;
                           prettyurls=get(ENV, "CI", "false") == "true",
                           canonical="https://juliatatoms.github.io/Hydrogen.jl",
                           assets=String["assets/latex.js"],
                           mathengine = Documenter.MathJax()
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaAtoms/Hydrogen.jl",
)
