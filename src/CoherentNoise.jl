module CoherentNoise

include("common.jl")
include("generators.jl")

using Reexport

@reexport using .Common: AbstractSampler, sample, gen_image
@reexport using .Noise: NoiseSampler
@reexport using .Noise.ValueNoise: value_2d, value_3d, value_4d
@reexport using .Noise.PerlinImprovedNoise:
    perlin_improved_2d, perlin_improved_3d, perlin_improved_4d
@reexport using .Noise.SimplexNoise: simplex_2d, simplex_3d, simplex_4d
@reexport using .Noise.OpenSimplexNoise: opensimplex_2d, opensimplex_3d, opensimplex_4d
@reexport using .Noise.OpenSimplex2Noise: opensimplex2_2d, opensimplex2_3d, opensimplex2_4d
@reexport using .Noise.OpenSimplex2SNoise: opensimplex2s_2d, opensimplex2s_3d, opensimplex2s_4d
@reexport using .Noise.WorleyNoise: worley_2d, worley_3d, worley_4d
@reexport using .Fractals: FractalSampler
@reexport using .Fractals: fbm_fractal_2d, fbm_fractal_3d, fbm_fractal_4d
@reexport using .Fractals: billow_fractal_2d, billow_fractal_3d, billow_fractal_4d
@reexport using .Fractals: multi_fractal_2d, multi_fractal_3d, multi_fractal_4d
@reexport using .Fractals: hybrid_fractal_2d, hybrid_fractal_3d, hybrid_fractal_4d
@reexport using .Fractals: ridged_fractal_2d, ridged_fractal_3d, ridged_fractal_4d
@reexport using .Patterns: PatternSampler, constant_1d, checkered_2d, cylinders_2d, spheres_3d
@reexport using .Modifiers: ModifierSampler, cache, curve, mix, rotate, scale, select, terrace
@reexport using .Modifiers: translate, turbulence, warp

end
