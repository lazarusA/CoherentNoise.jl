module OpenSimplexNoise

using Random: shuffle
using FastPow
using ...Common: RandomState, PerlinState, Seed
import ...Common: sample
using ..Noise: NoiseSampler

"""
An `N`-dimensional sampler that generates OpenSimplex (legacy) noise.

`N` must be an integer in the range [2, 4].

Note: This version of OpenSimplex is considered "legacy" by its original author. Consider using one
of the newer algorithms, [`OpenSimplex2`](@ref Main.CoherentNoise.Noise.OpenSimplex2Noise.OpenSimplex2) 
or [`OpenSimplex2S`](@ref Main.CoherentNoise.Noise.OpenSimplex2SNoise.OpenSimplex2S).
"""
struct OpenSimplex{N} <: NoiseSampler{N}
    random_state::RandomState
    state::PerlinState
    table::Vector{UInt8}
end

"""
    OpenSimplex{N}(; seed=nothing)

Construct a sampler that outputs `N`-dimensional OpenSimplex noise when it is sampler from.

# Arguments

  - `seed::Union{Int,Nothing}=nothing`: An integer used to seed the random number generator for this
    sampler, or `nothing`. If a seed is not supplied, one will be generated automatically which will
    negatively affect reproducibility.
"""
function OpenSimplex{N}(; seed::Seed=nothing) where {N}
    rs = RandomState(seed)
    table = shuffle(rs.rng, UInt8.([repeat(0:3:45, 11)..., repeat(48:3:69, 10)...]))
    OpenSimplex{N}(rs, PerlinState(rs), table)
end

include("noise_opensimplex_2d.jl")
include("noise_opensimplex_3d.jl")
include("noise_opensimplex_4d.jl")

end
