```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "CoherentNoise.jl"
  text: "All the noise you need, and more!"
  tagline: A comprehensive suite of coherent noise algorithms
  image:
    src: /logo.png
    alt: CoherentNoise
  actions:
    - theme: brand
      text: Getting Started
      link: /getting_started
    - theme: alt
      text: View on Github
      link: https://github.com/lazarusA/CoherentNoise.jl
    - theme: alt
      text: API
      link: /reference
features:
  - title: Coherent noise?
    details: What is it?
    link: /overview
  - title: Algorithms
    details: Explore the ones available.
    link: /algorithms
  - title: Gallery
    details: User's creations. Feel free to submit your own.
    link: /gallery
---
```

How to Install CoherentNoise.jl?

Its easy to install CoherentNoise.jl. Since is registered in the Julia General registry, you can simply run the following command in the Julia REPL:

```julia
julia> using Pkg
julia> Pkg.add("CoherentNoise")
```