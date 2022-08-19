module Patterns
include("pattern.jl")
include("pattern_constant_1d.jl")
include("pattern_checkered_2d.jl")
include("pattern_cylinders_2d.jl")
include("pattern_spheres_3d.jl")
end

module Noise
include("noise.jl")
include("noise_value.jl")
include("noise_perlin_improved.jl")
include("noise_simplex.jl")
include("noise_opensimplex.jl")
include("noise_opensimplex2.jl")
include("noise_opensimplex2s.jl")
include("noise_worley.jl")
end

module Fractals
include("fractal.jl")
include("fractal_fbm.jl")
include("fractal_billow.jl")
include("fractal_multi.jl")
include("fractal_hybrid.jl")
include("fractal_ridged.jl")
end

module Modifiers
include("modifiers.jl")
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
end
