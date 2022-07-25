struct Pow{N,S} <: ModifierSampler{N}
    random_state::RandomState
    source::S
    power::Float64
end

"""
    ^(x::AbstractSampler, y::Real)

Construct a modifier sampler that raises the output of sampler `x` to the power of the scalar `y`.
"""
Base.:^(x::S, y::Real) where {N,S<:AbstractSampler{N}} = Pow{N,S}(random_state(x), x, Float64(y))

@inline function sample(sampler::Pow{N}, coords::Vararg{Real,N}) where {N}
    x = sample(sampler.source, coords...)
    abs((x + 1) * 0.5)^sampler.power * 2 - 1
end
