module PerlinNoise

using ...Common: RandomState, PerlinState, Seed, IsPerlinHashed, lerp, hash_coords, curve5
import ...Common: HashTrait, sample
using ..Noise: NoiseSampler

"""
An `N`-dimensional sampler that generates Perlin "Improved" noise.

`N` must be an integer in the range [2, 4].
"""
struct Perlin{N} <: NoiseSampler{N}
    random_state::RandomState
    state::PerlinState
end

"""
    Perlin{N}(; seed=nothing)

Construct a sampler that outputs `N`-dimensional Perlin "Improved" noise when it is sampler from.

# Arguments

  - `seed::Union{Int,Nothing}=nothing`: An integer used to seed the random number generator for this
    sampler, or `nothing`. If a seed is not supplied, one will be generated automatically which will
    negatively affect reproducibility.
"""
function Perlin{N}(; seed::Seed=nothing) where {N}
    rs = RandomState(seed)
    Perlin{N}(rs, PerlinState(rs))
end

HashTrait(::Type{<:Perlin}) = IsPerlinHashed()

include("noise_perlin_2d.jl")
include("noise_perlin_3d.jl")
include("noise_perlin_4d.jl")

end
