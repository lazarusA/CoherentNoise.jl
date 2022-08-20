module CoherentNoise

using Random: Xoshiro, RandomDevice, seed!, shuffle

using CircularArrays: CircularVector
using Distributions: Uniform
using FastPow
using ColorTypes: RGB
using ColorSchemes: ColorScheme

include("common.jl")
include("pattern_constant_1d.jl")
include("pattern_checkered_2d.jl")
include("pattern_cylinders_2d.jl")
include("pattern_spheres_3d.jl")
include("noise_value.jl")
include("noise_perlin_improved.jl")
include("noise_simplex.jl")
include("noise_opensimplex.jl")
include("noise_opensimplex2.jl")
include("noise_opensimplex2s.jl")
include("noise_worley.jl")
include("fractal.jl")
include("fractal_fbm.jl")
include("fractal_billow.jl")
include("fractal_multi.jl")
include("fractal_hybrid.jl")
include("fractal_ridged.jl")
include("mod_abs.jl")
include("mod_add.jl")
include("mod_cache.jl")
include("mod_clamp.jl")
include("mod_copysign.jl")
include("mod_curve.jl")
include("mod_div.jl")
include("mod_max.jl")
include("mod_min.jl")
include("mod_mix.jl")
include("mod_mul.jl")
include("mod_muladd.jl")
include("mod_pow.jl")
include("mod_rotate.jl")
include("mod_scale.jl")
include("mod_select.jl")
include("mod_sub.jl")
include("mod_terrace.jl")
include("mod_translate.jl")
include("mod_turbulence.jl")
include("mod_warp.jl")

# Common
export AbstractSampler, sample, gen_image

# Noise
export NoiseSampler
export value_2d, value_3d, value_4d
export perlin_improved_2d, perlin_improved_3d, perlin_improved_4d
export simplex_2d, simplex_3d, simplex_4d
export opensimplex_2d, opensimplex_3d, opensimplex_4d
export opensimplex2_2d, opensimplex2_3d, opensimplex2_4d
export opensimplex2s_2d, opensimplex2s_3d, opensimplex2s_4d
export worley_2d, worley_3d, worley_4d

# Fractals
export FractalSampler
export fbm_fractal_2d, fbm_fractal_3d, fbm_fractal_4d
export billow_fractal_2d, billow_fractal_3d, billow_fractal_4d
export multi_fractal_2d, multi_fractal_3d, multi_fractal_4d
export hybrid_fractal_2d, hybrid_fractal_3d, hybrid_fractal_4d
export ridged_fractal_2d, ridged_fractal_3d, ridged_fractal_4d

# Patterns
export PatternSampler
export constant_1d, checkered_2d, cylinders_2d, spheres_3d

# Modifiers
export ModifierSampler
export cache, curve, mix, rotate, scale, select, terrace
export translate, turbulence, warp

end
