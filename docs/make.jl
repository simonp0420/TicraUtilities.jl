using Pkg: Pkg
Pkg.activate(@__DIR__)
using Literate, Documenter, DocumenterCitations, TicraUtilities
olddir = pwd()
cd(@__DIR__)

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:numeric
)

literate_list = ["Contents", "index", "Tutorial", "autodocs"]
for file in literate_list
    Literate.markdown(joinpath("literate", file*".jl"), "src")
end

#Literate.notebook(joinpath("literate", "Tutorial.jl"), "notebooks"; execute=false)

# Define the local path to package source
#package_root = abspath(joinpath(@__DIR__, ".."))

# Define the remote repository information using a Documenter.Remotes type (e.g., 
#remote_repo = Documenter.Remotes.GitHub("simonp0420", "TicraUtilities.jl") 


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
#    repo = remote_repo,
#    remotes = nothing # Dict(package_root => remote_repo),
)

deploydocs(
    repo = "github.com/simonp0420/TicraUtilities.jl.git",
)

cd(olddir)