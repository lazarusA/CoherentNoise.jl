```@setup tutorial
using CoherentNoise, ImageIO, FileIO, ColorSchemes, Chain
function generate(file, sampler; xbounds=(-1.0, 1.0), ybounds=(-1.0, 1.0), colorscheme=nothing)
    img = gen_image(
        sampler, w=320, h=320, xbounds=xbounds, ybounds=ybounds, colorscheme=colorscheme)
    save(joinpath(@__DIR__, "assets/tutorial/$(file).png"), img)
end
```

## [First steps](@id first_steps)

Before following along with this tutorial, it is a good idea to get familiar with how to work with
samplers. See [Getting Started](@ref installation).

### Visualizing results

CoherentNoise provides basic support for generating 2-dimensional slices of noise data, which can be
conveniently written to disk as an image file, or rendered with tools such as
[ImageView.jl](https://github.com/JuliaImages/ImageView.jl).

To generate an image that samples from a sampler for each pixel, it is as simple as the following,
with `sampler` being any instance of the supported sampler types:

```@example tutorial
sampler = opensimplex2_2d(seed=1)
gen_image(sampler)
generate("tutorial01", sampler) # hide
```

This produces an array of RGB values using
[ColorTypes.jl](https://github.com/JuliaGraphics/ColorTypes.jl). You can then write it out as an
image file with something like:

```@julia
using FileIO
img = gen_image(sampler)
save(img, "/tmp/noise.png")
```

Which, for this example, should produce something similar to the following image:

![](assets/tutorial/tutorial01.png)

### Modifiers

Modifiers are the backbone of producing interesting noise results. Using a modifier, we can alter
the input coordinates before sampling, or the output value after sampling. Modifiers are just
samplers, so they have no effect until `sample` is called on them. We can therefor compose long
pipelines of noises generators, and finally evaluate the graph, as the end result is always just an
ordinary sampler.

For example, taking our previous example from [Visualizing results](#Visualizing-results), we can
start modifying it with a simple modifier; unary subtraction (negation) of the output values.

```@example tutorial
sampler = -opensimplex2_2d(seed=1)
gen_image(sampler)
generate("tutorial02", sampler) # hide
```

![](assets/tutorial/tutorial02.png)

Notice the `-` operator in front of the sampler function call; that was the only change. This
modifier simply negates the output values.

We can get more advanced by modifying the output values *according* to the outputs of other
samplers. One such modifier is one created by the `mix()` function. This function accepts 3 samplers
as input, and linearly interpolates between the outputs of the first two samplers, using the third
sampler's output as a blending factor.

Let's mix our original sampler together with another new sampler, with interpolation controlled by
yet another new sampler. But first, let's see what each of these look like individually:

```@example tutorial
sampler1 = opensimplex2_2d(seed=1)
sampler2 = opensimplex2_2d(seed=2)
sampler3 = perlin_3d(seed=1)
gen_image(sampler1)
gen_image(sampler2)
gen_image(sampler3)
generate("tutorial03", sampler1) # hide
generate("tutorial04", sampler2) # hide
generate("tutorial05", sampler3) # hide
```

sampler1                            |  sampler2                           | sampler3
:----------------------------------:|:-----------------------------------:|:----------------------------------:
![](assets/tutorial/tutorial03.png) | ![](assets/tutorial/tutorial04.png) | ![](assets/tutorial/tutorial05.png)

Interpolating these with mix is done as follows:

```@example tutorial
sampler1 = opensimplex2_2d(seed=1)
sampler2 = opensimplex2_2d(seed=2)
sampler3 = perlin_2d(seed=1)
mixed = mix(sampler1, sampler2, sampler3)
gen_image(mixed)
generate("tutorial06", mixed) # hide
```

![](assets/tutorial/tutorial06.png)

It turns out, this doesn't look very interesting. This is mostly due to the fact that we are zoomed
in very close to the generated noise. By default, `gen_image()` samples from coordinate [-1.0, 1.0]
on both the X and Y axes. This can be tuned with the `xbounds` and `ybounds` parameters, but rather
than doing that, let's use the `scale` modifier instead.

The `scale` modifier modifies the input coordinates of a source sampler before it is sampled from.
So far we have only seen an output modifying sampler, but one can also modify input coordinates with
other modifiers. Here is how we can scale down our previous example. Throughout this tutorial, the
full source will be included with each example, to make it easier to copy and paste to your REPL to
try things out on your own.

```@example tutorial
sampler1 = opensimplex2_2d(seed=1)
sampler2 = opensimplex2_2d(seed=2)
sampler3 = perlin_2d(seed=1)
mixed = mix(sampler1, sampler2, sampler3)
scaled = scale(mixed, 0.5)
gen_image(scaled)
generate("tutorial07", scaled) # hide
```

![](assets/tutorial/tutorial07.png)

We decreased the scale by half, making the result look slightly more interesting.

There are many more modifiers. It is recommended to check out the [API Reference](@ref) for a full
list of modifiers (and other samplers). Also, be sure to check out the [Gallery](@ref
user_creations) for more examples.

### Pipelines

Once you start composing different noises together with modifiers, the code gets a little
un-manageable and harder to read. We can leverage Julia's builtin pipe operator, or more
conveniently, the [Chain.jl](https://github.com/jkrumbiegel/Chain.jl) package's `@chain` macro.

Let's rewrite our last example using `@chain`

```@example tutorial
@chain mix(opensimplex2_2d(seed=1), opensimplex2_2d(seed=2), perlin_2d(seed=1)) begin
    scale(0.5)
    gen_image
end
@chain mix(opensimplex2_2d(seed=1), opensimplex2_2d(seed=2), perlin_2d(seed=1)) begin # hide
    scale(0.5) # hide
    generate("tutorial08", _) # hide
end # hide
```

![](assets/tutorial/tutorial08.png)

Admittedly, this doesn't look that much better, but the more your noise composition grows, the more
this variable-free clear flow of control will be appealing. The remainder of this tutorial will use
@chain, but it's not a necessary part of the workflow. You can also write the pipeline code manually
without it, if this style is not for you.

### Fractals

A fractal is another type of sampler that applies multiple octaves of a given source sampler,
varying their frequency. They make the noise much less regular looking, and offers a good starting
point for transforming noises into the desired texture you are trying to simulate, should you be
using this library for procedural image generation.

CoherentNoise has a few different fractal sampler types, some of which will be demonstrated below.

We will start with the following regular noise, and apply different fractal samplers to it:

```@example tutorial
@chain opensimplex2_3d(seed=1) begin
    gen_image
end
@chain opensimplex2_3d(seed=1) begin # hide
    generate("tutorial09", _) # hide
end # hide
```

![](assets/tutorial/tutorial09.png)

#### fBm (fractional Brownian motion)

The most basic example is that with only a `source` argument supplied.

```@example tutorial
@chain opensimplex2_3d(seed=1) begin
    fbm_fractal_3d(source=_)
    gen_image
end
@chain opensimplex2_3d(seed=1) begin # hide
    fbm_fractal_3d(source=_) # hide
    generate("tutorial10", _) # hide
end # hide
```

![](assets/tutorial/tutorial10.png)

Already, this is starting to look more natural.

All fractal samplers have various knobs to tweak to control the output. Please refer to the [API
Reference](@ref) for details, but here are some more examples of fBm using our same source
sampler:

```@example tutorial
@chain opensimplex2_3d(seed=1) begin
    fbm_fractal_3d(seed=1, source=_, octaves=3)
    gen_image
end
@chain opensimplex2_3d(seed=1) begin # hide
    fbm_fractal_3d(seed=1, source=_, octaves=3) # hide
    generate("tutorial11", _) # hide
end # hide
```

![](assets/tutorial/tutorial11.png)

```@example tutorial
@chain opensimplex2_3d(seed=1) begin
    fbm_fractal_3d(seed=1, source=_, frequency=3)
    gen_image
end
@chain opensimplex2_3d(seed=1) begin # hide
    fbm_fractal_3d(seed=1, source=_, frequency=3) # hide
    generate("tutorial12", _) # hide
end # hide
```

![](assets/tutorial/tutorial12.png)

```@example tutorial
@chain opensimplex2_3d(seed=1) begin
    fbm_fractal_3d(seed=1, source=_, frequency=3, persistence=0.7)
    gen_image
end
@chain opensimplex2_3d(seed=1) begin # hide
    fbm_fractal_3d(seed=1, source=_, frequency=3, persistence=0.7) # hide
    generate("tutorial13", _) # hide
end # hide
```

![](assets/tutorial/tutorial13.png)


```@example tutorial
@chain opensimplex2_3d(seed=1) begin
    fbm_fractal_3d(seed=1, source=_, frequency=3, persistence=0.7, lacunarity=1.2)
    gen_image
end
@chain opensimplex2_3d(seed=1) begin # hide
    fbm_fractal_3d(seed=1, source=_, frequency=3, persistence=0.7, lacunarity=1.5) # hide
    generate("tutorial14", _) # hide
end # hide
```

![](assets/tutorial/tutorial14.png)

#### Billow

`Billow` is another fractal sampler type that produces a billowy appearance. It could be used as the
basis for a clouds texture, among many other uses. It supports all of the same parameters as fBm.
Here is an example:

```@example tutorial
@chain opensimplex2_3d(seed=1) begin
    billow_fractal_3d(seed=1, source=_, frequency=2.5, persistence=0.55, lacunarity=1.7)
    gen_image
end
@chain opensimplex2_3d(seed=1) begin # hide
    billow_fractal_3d(seed=1, source=_, frequency=2.5, persistence=0.55, lacunarity=1.7) # hide
    generate("tutorial15", _) # hide
end # hide
```

![](assets/tutorial/tutorial15.png)

#### Ridged

`Ridged` is a fractal (actually a [multifractal](https://en.wikipedia.org/wiki/Multifractal_system))
sampler type that produces ridges. It could be used as the basis for a mountainous terrain texture,
among many other uses. It supports all of the same parameters as all other fractals, with the
addition of an `attenuation` parameter. Here are some examples:

```@example tutorial
@chain opensimplex2_3d(seed=1) begin
    ridged_fractal_3d(seed=1, source=_, frequency=2.5, persistence=0.4, attenuation=1)
    gen_image
end
@chain opensimplex2_3d(seed=1) begin # hide
    ridged_fractal_3d(seed=1, source=_, frequency=2.5, persistence=0.4, attenuation=1) # hide
    generate("tutorial16", _) # hide
end # hide
```

![](assets/tutorial/tutorial16.png)

```@example tutorial
@chain opensimplex2s_3d(seed=1) begin
    ridged_fractal_3d(seed=1, source=_, octaves=6, frequency=1.8, persistence=0.48)
    gen_image
end
@chain opensimplex2s_3d(seed=1) begin # hide
    ridged_fractal_3d(seed=1, source=_, octaves=6, frequency=1.8, persistence=0.48) # hide
    generate("tutorial17", _) # hide
end # hide
```

![](assets/tutorial/tutorial17.png)

There are a couple other fractal types that we didn't go over, but feel free to experiment with
them. Their interface is the same as any other fractal type.

### Adding color

CoherentNoise works in combination with
[ColorSchemes.jl](https://github.com/JuliaGraphics/ColorSchemes.jl) to add color to generated
images. To do so, just pass the `colorscheme` keyword argument to `gen_image` with any color scheme
provided by ColorSchemes.jl, or [create your
own](https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#Make-your-own-colorscheme).
Adding some color to the last Ridged fractal example looks like this:

```@example tutorial
@chain opensimplex2s_3d(seed=1) begin
    ridged_fractal_3d(seed=1, source=_, octaves=6, frequency=1.2, attenuation=1.3, lacunarity=4.5, persistence=0.58)
    gen_image(colorscheme=ColorSchemes.terrain)
end
@chain opensimplex2s_3d(seed=1) begin # hide
    ridged_fractal_3d(seed=1, source=_, octaves=6, frequency=1.2, attenuation=1.3, lacunarity=4.5, persistence=0.58) # hide
    generate("tutorial18", _, colorscheme=ColorSchemes.terrain) # hide
end # hide
```

![](assets/tutorial/tutorial18.png)

Here is the same thing with the noise output negated with the unary `-` modifier:

```@example tutorial
@chain opensimplex2s_3d(seed=1) begin
    ridged_fractal_3d(seed=1, source=_, octaves=6, frequency=1.2, attenuation=1.3, lacunarity=4.5, persistence=0.58)
    -
    gen_image(colorscheme=ColorSchemes.terrain)
end
@chain opensimplex2s_3d(seed=1) begin # hide
    ridged_fractal_3d(seed=1, source=_, octaves=6, frequency=1.2, attenuation=1.3, lacunarity=4.5, persistence=0.58) # hide
    - # hide
    generate("tutorial19", _, colorscheme=ColorSchemes.terrain) # hide
end # hide
```

![](assets/tutorial/tutorial19.png)

That looks pretty interesting. With a bit more work, you could get it to look much more realistic
and terrain-like. For example, to make it more arid, you can simply add an `abs` modifier:

```@example tutorial
@chain opensimplex2s_3d(seed=1) begin
    ridged_fractal_3d(seed=1, source=_, octaves=6, frequency=1.2, attenuation=1.3, lacunarity=4.5, persistence=0.58)
    -
    abs
    gen_image(colorscheme=ColorSchemes.terrain)
end
@chain opensimplex2s_3d(seed=1) begin # hide
    ridged_fractal_3d(seed=1, source=_, octaves=6, frequency=1.2, attenuation=1.3, lacunarity=4.5, persistence=0.58) # hide
    - # hide
    abs # hide
    generate("tutorial20", _, colorscheme=ColorSchemes.terrain) # hide
end # hide
```

![](assets/tutorial/tutorial20.png)

That concludes the basic tutorial. Feel free to experiment with complex noise pipelines, and submit
your results as a pull request to be included in the [Gallery](@ref user_creations) section!
