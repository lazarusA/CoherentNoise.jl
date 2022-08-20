struct PerlinImproved{N} <: NoiseSampler{N}
    random_state::RandomState
    state::PerlinState
end

@inline function perlin_improved(dims, seed)
    rs = RandomState(seed)
    PerlinImproved{dims}(rs, PerlinState(rs))
end

HashTrait(::Type{<:PerlinImproved}) = IsPerlinHashed()

include("noise_perlin_improved_2d.jl")
include("noise_perlin_improved_3d.jl")
include("noise_perlin_improved_4d.jl")
