module SimplexNoise

using FastPow
using ...Common: RandomState, PerlinState, Seed, IsPerlinHashed, hash_coords
import ...Common: HashTrait, sample
using ..Noise: NoiseSampler

"""
An `N`-dimensional sampler that generates Simplex noise.

`N` must be an integer in the range [2, 4].

Note: ≥ 3 dimensions is patent-protected for certain applications. Consider using one of the
OpenSimplex implementations, instead.
"""
struct Simplex{N} <: NoiseSampler{N}
    random_state::RandomState
    state::PerlinState
end

"""
    Simplex{N}(; seed=nothing)

Construct a sampler that outputs `N`-dimensional Simplex noise when it is sampler from.

# Arguments

  - `seed::Union{Int,Nothing}=nothing`: An integer used to seed the random number generator for this
    sampler, or `nothing`. If a seed is not supplied, one will be generated automatically which will
    negatively affect reproducibility.
"""
function Simplex{N}(; seed::Seed=nothing) where {N}
    rs = RandomState(seed)
    Simplex{N}(rs, PerlinState(rs))
end

HashTrait(::Type{<:Simplex}) = IsPerlinHashed()

include("noise_simplex_2d.jl")
include("noise_simplex_3d.jl")
include("noise_simplex_4d.jl")

end
