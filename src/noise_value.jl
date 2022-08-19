module ValueNoise

using ...Common: RandomState, Seed, IsValueHashed, hash_coords, lerp, curve3
using ...Common: PRIME_X, PRIME_Y, PRIME_Z, PRIME_W, get_seed
import ...Common: HashTrait, sample
using ..Noise: NoiseSampler

struct Value{N} <: NoiseSampler{N}
    random_state::RandomState
end

@inline value(dims, seed::Seed) = Value{dims}(RandomState(seed))

HashTrait(::Type{<:Value}) = IsValueHashed()

include("noise_value_2d.jl")
include("noise_value_3d.jl")
include("noise_value_4d.jl")

end
