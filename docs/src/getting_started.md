```@setup getting-started
using CoherentNoise
```
## [Installation](@id installation)

CoherentNoise.jl exists in the Julia's General registry. In Julia ≥ 1.0, you can add it using the
package manager with:

```julia-repl
pkg> add CoherentNoise
```

## Basic Usage

The following is only a brief explanation of basic usage. For more advanced usage examples, please
see the [Tutorial](@ref first_steps) section.

To get started, let's get a feel for how to create a sampler for one of the supported noise
algorithms, and sample from it. Perlin Improved noise is a well-known algorithm, so let's create a
2-dimensional Perlin Improved noise sampler and sample from it:

```@example getting-started
using CoherentNoise
sampler = perlin_2d()
sample(sampler, 120.2, 42.8)
```

To create a sampler, we simply call the function corresponding to the noise algorithm and
dimensionality. We can then sample from it using the `sample` function, passing in as arguments a
sampler, and multiple `Real` arguments corresponding to the coordinates in the sampler's noise
space. In this example, we have a 2-dimensional sampler, so we passed 2 numbers; one for the
coordinate along the X axis, and one for the coordinate along the Y axis.

!!! note
    It is strongly recommended to use floating-point numbers with a fractional component (non-zero
    value after the decimal point), as some algorithms (notably older gradient noises like
    Perlin Noise or Perlin "Improved Noise"), will return zero for integral coordinates.

All samplers have their own distinct random number generator, and it is how we can retrieve
different results with the same input coordinates. By default, a sampler's seed is set to `nothing`,
which corresponds to using your machine's hardware random number generator to create a seed. This
will result in non-deterministic results, as each instance of a sampler will produce different
results when sampled from. We can change the seed on a per-sampler basis, by passing the `seed`
keyword argument, in order to have reproducible results:

```@example getting-started
sampler = perlin_2d(seed=42)
sample(sampler, 120.2, 42.8)
```

In summary, one creates a sampler by calling a function for the desired algorithm and
dimensionality, and we can alter the the sampler's output by changing its seed.

Of particular note is that all samplers accept a `seed` keyword argument; even those that don't make
use of any random numbers in their implementation. This is required for the composition pipeline
feature described in the [Tutorial](@ref first_steps).
