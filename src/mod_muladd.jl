struct Muladd{N,S} <: ModifierSampler{N}
    random_state::RandomState
    source::S
    strength::Float64
    bias::Float64
end

"""
    muladd(x::AbstractSampler, strength::Real=1.0, bias::Real=0.0)

Construct a modifier sampler that performs multiplies the output of sampler `x` by the scalar
`strength`, followed by adding the scalar `bias`. sampler `y`.
"""
function Base.muladd(x::S, strength::Real=1.0, bias::Real=0.0) where {N,S<:AbstractSampler{N}}
    Muladd{N,S}(random_state(x), x, Float64(strength), Float64(bias))
end

@inline function sample(
    sampler::Muladd{N,S},
    coords::Vararg{Real,N},
) where {N,S<:AbstractSampler{N}}
    x = sample(sampler.source, coords...)
    muladd(x, sampler.strength, sampler.bias)
end
