module CoherentNoise

include("common.jl")
include("generators.jl")

using Reexport

@reexport using .Common: AbstractSampler, sample, gen_image
@reexport using .Noise: NoiseSampler
@reexport using .Noise.ValueNoise: Value
@reexport using .Noise.PerlinNoise: Perlin
@reexport using .Noise.SimplexNoise: Simplex
@reexport using .Noise.OpenSimplexNoise: OpenSimplex
@reexport using .Noise.OpenSimplex2Noise: OpenSimplex2
@reexport using .Noise.OpenSimplex2Noise: Standard, ImproveX, ImproveXY, ImproveXZ, ImproveXYZ
@reexport using .Noise.OpenSimplex2SNoise: OpenSimplex2S
@reexport using .Noise.CellularNoise: Cellular
@reexport using .Noise.CellularNoise: Manhattan, Euclidean, Euclidean², Chebyshev, Minkowski4
@reexport using .Noise.CellularNoise: Distance1, Distance2, DistanceAdd, DistanceSub, DistanceMul
@reexport using .Noise.CellularNoise: DistanceDiv, CellValue
@reexport using .Fractals: FractalSampler, FBM, Billow, Multifractal, Hybrid, Ridged
@reexport using .Patterns: PatternSampler, Constant, Checkered, Cylinders, Spheres
@reexport using .Modifiers: ModifierSampler, cache, curve, mix, rotate, scale, select, terrace
@reexport using .Modifiers: translate, turbulence, warp

end
