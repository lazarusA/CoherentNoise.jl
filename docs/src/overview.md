## What is coherent noise?

Unlike its more random cousin, coherent noise has interesting properties that give it a smooth
appearance. A coherent noise function is a function that emits a single scalar value as its output.
All coherent noise functions exhibit the following three properties:

1. The same input arguments will always result in the same output.
2. A small change in input will result in a small change in the output.
3. A large change in input will result in a random output.

The dimensionality of a coherent noise function is determined by how many input arguments it
accepts. Since inputs to a coherent noise function could be driven from any source of data, and
output values are spatially coherent in each  dimension, one could feed image coordinates into a
3-dimensional coherent noise function's first two arguments, and an incrementing time value as the
third argument, in order to smoothly animate an image effect such as clouds, fire, landscapes, or
anything else creativity dreams up. Similarly, one could render a 3-dimensional volume, and animate
it along a 4th dimension.

Coherent noise functions are especially ubiquitous in applications where memory is scarce, and in
game development where one wants to generate 2-dimensional textures or 3-dimensional objects without
storing large chunks of image or geometry data on disk.

## History

CoherentNoise started out as a [Common Lisp library](https://git.mfiano.net/mfiano/cricket)
providing a collection of noise algorithms for use in game development, but evolved over time into a
more general form we thought could be useful to others outside of that domain. Additionally, the
author has since migrated from Common Lisp to Julia, and here this port with even more functionality
was born. The implementation borrows a lot of its ideas from [C++'s
libnoise](http://libnoise.sourceforge.net/) and other noise libraries, but adds numerous
enhancements, corrections, and additional noise algorithms.

As mentioned, the main author, Michael Fiano, has recently moved to Julia after nearly 20 years of
exclusive Common Lisp development. CoherentNoise is a realization of their first Julia library, and
as such, code may not be the most idiomatic in nature. However, it was a goal to produce something
as Julian as possible, and with a focus on performance and precision, than to directly port the
original Common Lisp implementation. If something can be improved in this regard, please open an
issue or pull request -- it would be wonderful if I could stengthen my Julia knowledge by discussion
or contribution by more experienced Julia developers that find value in CoherentNoise.

## What's different from other libraries?

- Many more noise algorithms than libnoise supported, and well-tested.

- Algorithms are consistently implemented in 2, 3, and 4 dimensions, whereas libnoise only has
  3-dimensional gradient noise.

- Some terminology corrections. For example, libnoise refers to multiple octaves of gradient noise
  as "Perlin" noise. Perlin is an implementation of gradient noise, and is not defined in terms of a
  fractal; infact, any noise algorithm can be "fractalized".

- Fix quite a few fundamental mathematical errors. It was quite the task re-implementing some
  algorithms due to incorrect math, or the inclusion of magic numbers that differed among
  implementations for various reasons.

- Ensure that all samplers emit values in the [-1.0, 1.0] range; something that becomes very
  important when working with fractal-based noises, composition of complex noises, and in general is
  something just expected when working with coherent noise.

- The `PerlinImproved` noise implementation is [Ken Perlin's "Improved Noise"
  algorithm](https://cs.nyu.edu/~perlin/noise/), unlike his classic Perlin noise algorithm found in
  libnoise, which is sometimes just called gradient noise.

- Fractal samplers can be driven by any source sampler, not just gradient noise -- even other
  fractal samplers.

- Additional modifier samplers.

### How is it used?

In summary, you call a function corresponding to the coherent noise algorithm and dimensionality you
want to generate, in order to create a "sampler". A sampler can then be sampled from with the
`sample` function, by passing to it a sampler instance, along with input coordinates for each of the
sampler's dimensions.

In addition to noise algorithm samplers, there are simple pattern samplers, and higher-order fractal
and modifier samplers.

Pattern samplers create simple patterns, like checkerboards, and are useful in noise composition
pipelines.

Higher-order samplers include both fractal samplers and modifier samplers.

Fractal sampler functions accept another sampler as an input source, and produce much more
interesting results by combining the results of multiple octaves of the source.

Modifier sampler functions also accept other sampler instances as input; and modify their input
coordinates or output values in interesting ways. Modifiers are the crux of CoherentNoise's
composition pipeline support, and allow merging different samplers together, even other modifiers,
into a composition of potentially very complex and interesting results.

Finally, there are tools available to conveniently work with 2-dimensional slices of
multi-dimensional noise spaces in order to generate entire images effortlessly.

## Available samplers

The following noise algorithms and other samplers are implemented. If you would like to see support
for other algorithms, please file an issue or submit a pull request.

### Noise samplers

- Value: (1D, 2D, 3D, 4D)
- Perlin: (1D, 2D, 3D, 4D)
- Simplex: (1D, 2D, 3D, 4D)
- OpenSimplex: (2D, 3D, 4D)
- OpenSimplex2: (2D, 3D, 4D)
- OpenSimplex2S: (2D, 3D, 4D)
- Worley: (1D, 2D, 3D, 4D)

### Pattern samplers

- Constant: (1D)
- Checkered: (2D)
- Cylinders: (2D)
- Spheres: (3D)

### Fractal samplers

- fBm: (fractional Brownian motion; 1D, 2D, 3D, 4D)
- Billow: (1D, 2D, 3D, 4D)
- Multifractal: (1D, 2D, 3D, 4D)
- Hybrid: (1D, 2D, 3D, 4D)
- Ridged: (1D, 2D, 3D, 4D)

### Modifier samplers

- `+`: binary addition of two samplers
- `+`: binary addition of a sampler and a scalar
- `-`: binary subtraction of two samplers
- `-`: binary subtraction of a sampler and a scalar
- `-`: unary subtraction of a sampler
- `*`: binary multiplication of two samplers
- `*`: binary multiplication of a sampler and a scalar
- `/`: binary division of two samplers
- `/`: binary division of a sampler and a scalar
- `^`: binary exponentiation of two samplers
- `abs`: absolute value of a sampler
- `cache`: cache the last output of a sampler
- `clamp`: restrict the range of a sampler by the bounds specified by the output of two other
  samplers
- `clamp`: restrict the range of a sampler by the bounds specified by two scalars
- `copysign`: change the sign of a sampler to that of the output of another sampler
- `curve`: maps the output of a sampler onto an arbitrary user-defined curve
- `max`: take the maximum value of two samplers
- `min`: take the minimum value of two samplers
- `mix`: linearly interpolate the result of two samplers by the result of a third alpha sampler
- `mix`: linearly interpolate the result of two samplers by a scalar
- `muladd`: multiplication followed by addition of a sampler by two scalars
- `rotate`: rotate the noise space of a sampler around its origin by scalars denoting axis angles of
  rotation in radians
- `scale`: scale the input coordinates of a sampler by scalars for the specified axes
- `scale`: uniformly scale the input coordinates of a sampler along all axes
- `select`: select the result of one of two samplers decided upon by the output of a third sampler
- `terrace`: map the output of a sampler onto a terrace-forming curve
- `translate`: translate (move) the noise space of a sampler along the specified axes
- `turbulence`: randomly displace the input coordinates of a sampler before sampling from it
- `warp`: domain-warps the sampler by modifying its input coordinates to be summed with the result
  of another supplied sampler for each coordinate
