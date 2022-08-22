module CoherentNoise

### Dependencies

# Language dependencies
using Random: shuffle

# Third-party dependencies
using CircularArrays: CircularVector
using ColorSchemes: ColorScheme
using ColorTypes: RGB
using Distributions: Uniform
using FastPow
using RandomNumbers.Xorshifts: Xoroshiro128Star

### Includes

include("docs.jl")
include("common.jl")
include("pattern.jl")
include("noise.jl")
include("fractal.jl")
include("modifier.jl")

### Exports

# Common
export AbstractSampler, sample, gen_image

# Noise
export NoiseSampler
export value_1d, value_2d, value_3d, value_4d
export perlin_1d, perlin_2d, perlin_3d, perlin_4d
export simplex_1d, simplex_2d, simplex_3d, simplex_4d
export opensimplex_2d, opensimplex_3d, opensimplex_4d
export opensimplex2_2d, opensimplex2_3d, opensimplex2_4d
export opensimplex2s_2d, opensimplex2s_3d, opensimplex2s_4d
export worley_1d, worley_2d, worley_3d, worley_4d

# Fractals
export FractalSampler
export fbm_fractal_1d, fbm_fractal_2d, fbm_fractal_3d, fbm_fractal_4d
export billow_fractal_1d, billow_fractal_2d, billow_fractal_3d, billow_fractal_4d
export multi_fractal_1d, multi_fractal_2d, multi_fractal_3d, multi_fractal_4d
export hybrid_fractal_1d, hybrid_fractal_2d, hybrid_fractal_3d, hybrid_fractal_4d
export ridged_fractal_1d, ridged_fractal_2d, ridged_fractal_3d, ridged_fractal_4d

# Patterns
export PatternSampler
export constant_1d, checkered_2d, cylinders_2d, spheres_3d

# Modifiers
export ModifierSampler
export cache, curve, mix, rotate, scale, select, terrace
export translate, turbulence, warp

end
