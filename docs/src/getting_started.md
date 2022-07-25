## Installation

CoherentNoise.jl exists in the Julia's General registry. In Julia ≥ 1.0, you can add it using the
package manager with:

```julia-repl
pkg> add CoherentNoise
```

## Basic Usage

The following is only a brief explanation of basic usage. For more advanced usage examples, please
see the [Tutorial](tutorial.html) section.

To get started, let's get a feel for how to create a sampler and sample from it for one of the
supported noise algorithms. Perlin noise is a well-known algorithm so let's create a 2-dimensional
Perlin noise sample and sample from it:

```@example 1
using CoherentNoise
sampler = Perlin{2}()
sample(sampler, 120.2, 42.8)
```

To create a sampler we simply invoke the constructor for its type, passing in one of the supported
dimensions as a type parameter. We can then sample from it using the `sample` function, passing in
as arguments, a sampler, and multiple `Real` arguments corresponding to the coordinates in the
sampler's noise space. In this example, we have a 2-dimensional sampler, so we passed 2 numbers; one
for the coordinate along the X axis, and one for the coordinate along the Y axis.

!!! note
    It is strongly recommended to use floating-point numbers with a fractional component (non-zero
    value after the decimal point), as some algorithms (noteably Perlin noise), will return zero for
    integral coordinates.

There is one potential problem with the above example though: every time we create a sampler and
sample from it, it will produce a different output:

```@example 1
sampler = Perlin{2}()
sample(sampler, 120.2, 42.8)
```

Doesn't that conflict with the first property of coherent noise mentioned in "[What is coherent
noise?](@ref)"? No, it doesn't. By not supplying any arguments to the sampler constructor, we are
telling CoherentNoise to seed this sampler from your machine's hardware random number generator. All
samplers have their own distinct random number generator, and it is how we can retrieve different
results with the same input. We can change this behavior on a per-sampler basis, by supplying our
own seed:

```@example 2
using CoherentNoise
sampler = Perlin{2}(seed=42)
sample(sampler, 120.2, 42.8)
```

```@example 2
sampler = Perlin{2}(seed=42)
sample(sampler, 120.2, 42.8)
```

Now, when we run this multiple times, we will always get the same result, and this is guaranteed to
be reproducible across different machines.

!!! note
    Julia is free to change the implementation of its random number generator algorithms, even
    across patch versions, so this reproducibility may not hold across Julia versions. In the
    future, it might be worth considering to add support for CoherentNoise to use
    [StableRNGs.jl](https://github.com/JuliaRandom/StableRNGs.jl) or
    [RandomNumbers.jl](https://github.com/JuliaRandom/RandomNumbers.jl) for truly portable
    reproducibility, but this is not a high priority.

In summary, one creates a sampler by calling a constructor for the desired algorithm and
dimensionality, and reproducibility can be controlled by supplying a `seed` keyword argument on a
per-sampler basis.

Of particular note is that all samplers accept a `seed` keyword argument; even those that don't make
use of any random numbers in their implementation. This is required for the composition pipeline
feature described in the [Tutorial](tutorial.html).
