using Pkg: Pkg
Pkg.activate(@__DIR__)
using Literate, Documenter, TicraUtilities

literate_list = ["index", "Tutorial", "autodocs"]
for file in literate_list
    Literate.markdown(joinpath("literate", file*".jl"), "src")
end
makedocs(
    pages = [
        "Introduction" => "index.md",
        "Tutorial" => "Tutorial.md",
        "API Reference" => "autodocs.md",
    ],
    sitename = "TicraUtilities"
)