![Splash Image](docs/src/assets/splash.png)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://lazarusa.github.io/CoherentNoise.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://lazarusa.github.io/CoherentNoise.jl/dev/)
[![Build Status](https://github.com/lazarusa/CoherentNoise.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lazarusa/CoherentNoise.jl/actions/workflows/CI.yml?query=branch%3Amain)
![License](https://img.shields.io/github/license/lazarusa/CoherentNoise.jl)


## What's different from other libraries (eg. libnoise)?

- Many more noise algorithms than libnoise supported, and well-tested.

- Algorithms are consistently implemented in multiple dimensions, whereas libnoise only has
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

- The Perlin noise implementation is [Ken Perlin's "Improved Noise"
  algorithm](https://cs.nyu.edu/~perlin/noise/), unlike his classic Perlin noise algorithm found in
  libnoise, which is sometimes just called gradient noise.

- Fractal samplers can be driven by any source sampler, not just gradient noise -- even other
  fractal samplers.

- Additional modifier samplers.

- Traditional Simplex noise variants, even the original Simplex by Ken Perlin, all overshoot the
  radial extent used for the signal reconstruction kernel. This results in unwanted artifacts when
  the noise is used for certain applications. CoherentNoise adds an option to use a more correct
  kernel to remove these artifacts. This is disabled by default. See the documentation for the
  respective Simplex noise variant for more details. To the best of my knowledge, CoherentNoise is
  the only library that offers this option for all Simplex noise variants.


## Available samplers

The following noise algorithms and other samplers are implemented. If you would like to see support
for other algorithms, please file an issue or submit a pull request.

### Noise samplers

In the order of "recommended to try first":

- OpenSimplex2: (2D, 3D, 4D)
- Simplex: (1D, 2D, 3D, 4D)
- OpenSimplex2S: (2D, 3D, 4D)
- OpenSimplex: (2D, 3D, 4D)
- Perlin: (1D, 2D, 3D, 4D)
- Value cubic: (1D, 2D, 3D, 4D)
- Value: (1D, 2D, 3D, 4D)
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

### [Modifier samplers](https://mfiano.github.io/CoherentNoise.jl/stable/overview/#Modifier-samplers)

