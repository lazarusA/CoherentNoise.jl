module ValueNoise

using ...Common: RandomState, Seed, IsValueHashed, hash_coords, lerp, curve3
using ...Common: PRIME_X, PRIME_Y, PRIME_Z, PRIME_W, get_seed
import ...Common: HashTrait, sample
using ..Noise: NoiseSampler

"""
An `N`-dimensional sampler that generates value noise.

`N` must be an integer in the range [2, 4].
"""
struct Value{N} <: NoiseSampler{N}
    random_state::RandomState
end

"""
    Value{N}(; seed=nothing)

Construct a sampler that outputs `N`-dimensional value noise when it is sampler from.

# Arguments

  - `seed::Union{Int,Nothing}=nothing`: An integer used to seed the random number generator for this
    sampler, or `nothing`. If a seed is not supplied, one will be generated automatically which will
    negatively affect reproducibility.
"""
Value{N}(; seed::Seed=nothing) where {N} = Value{N}(RandomState(seed))

HashTrait(::Type{<:Value}) = IsValueHashed()

include("noise_value_2d.jl")
include("noise_value_3d.jl")
include("noise_value_4d.jl")

end
