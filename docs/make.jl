using Pkg: Pkg
Pkg.activate(@__DIR__)
using Literate, Documenter, DocumenterCitations, TicraUtilities

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:numeric
)

literate_list = ["Contents", "index", "Tutorial", "autodocs"]
for file in literate_list
    Literate.markdown(joinpath("literate", file*".jl"), "src")
end
makedocs(;
    format = Documenter.HTML(
        assets=String["assets/citations.css"],
    ),
    pages = [
        "Contents" => "Contents.md",
        "Introduction" => "index.md",
        "Tutorial" => "Tutorial.md",
        "API Reference" => "autodocs.md",
        "References" => "references.md",
    ],
    sitename = "TicraUtilities",
    plugins = [bib,],
)