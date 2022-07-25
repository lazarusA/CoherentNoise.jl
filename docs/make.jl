using CoherentNoise
using Images, FileIO, Chain
using Documenter

DocMeta.setdocmeta!(CoherentNoise, :DocTestSetup, :(using CoherentNoise); recursive=true)

makedocs(;
    modules=[CoherentNoise],
    authors="Michael Fiano <mail@mfiano.net> and contributors",
    repo="https://github.com/mfiano/CoherentNoise.jl/blob/{commit}{path}#{line}",
    sitename="CoherentNoise.jl",
    #= linkcheck=true, =#
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mfiano.github.io/CoherentNoise.jl",
        edit_link="main",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
        "Overview" => "overview.md",
        "Getting Started" => "getting_started.md",
        "Tutorial" => "tutorial.md",
        "Gallery" => "gallery.md",
        "API Reference" => "reference.md"
    ]
)

deploydocs(;
    repo="github.com/mfiano/CoherentNoise.jl",
    devbranch="main"
)
