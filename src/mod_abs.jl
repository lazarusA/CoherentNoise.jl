struct Abs{N,S} <: ModifierSampler{N}
    random_state::RandomState
    source::S
end

"""
    abs(x::AbstractSampler)

Construct a modifier sampler that outputs the absolute value of its source when it is sampled from.
"""
function Base.abs(x::S) where {N,S<:AbstractSampler{N}}
    Abs{N,S}(x.random_state, x)
end

@inline function sample(sampler::Abs{N}, coords::Vararg{Real,N}) where {N}
    abs(sample(sampler.source, coords...))
end
