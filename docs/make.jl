using CoherentNoise
using FileIO, Chain
using Documenter
using DocumenterVitepress

DocMeta.setdocmeta!(CoherentNoise, :DocTestSetup, :(using CoherentNoise); recursive=true)

makedocs(;
    modules=[CoherentNoise],
    authors="Michael Fiano <mail@mfiano.net> and contributors",
    sitename="CoherentNoise.jl",
    clean=true,
    doctest=true,
    linkcheck=false,
    warnonly=[:missing_docs],
    format=DocumenterVitepress.MarkdownVitepress(;
        repo="github.com/lazarusA/CoherentNoise.jl", # this must be the full URL!
        devbranch="main",
        devurl="dev",
    ),
    draft=false,
    source="src",
    build="build",
    pages=[
        "Home" => "index.md",
        "Overview" => "overview.md",
        "Algorithms" => "algorithms.md",
        "Getting Started" => "getting_started.md",
        "Tutorial" => "tutorial.md",
        "Gallery" => "gallery.md",
        "API Reference" => "reference.md"
    ]
)

deploydocs(;
    repo="github.com/lazarusA/CoherentNoise.jl",
    target = "build", # this is where Vitepress stores its output
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true
)