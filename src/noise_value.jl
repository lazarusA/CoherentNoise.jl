struct Value{N} <: NoiseSampler{N}
    random_state::RandomState
end

@inline value(dims, seed=0) = Value{dims}(RandomState(seed))

HashTrait(::Type{<:Value}) = IsValueHashed()

include("noise_value_2d.jl")
include("noise_value_3d.jl")
include("noise_value_4d.jl")
