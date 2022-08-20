mutable struct Cache{N,S} <: ModifierSampler{N}
    random_state::RandomState
    source::S
    is_cached::Bool
    coords::NTuple{N,Float64}
    value::Float64
end

"""
    cache(x::AbstractSampler)

Construct a modifier sampler that caches the set of the input coordinates and their corresponding
output value of its source sampler. If the input coordinates differs from the previously cached
output, the cache is invalidated and the new output is cached.

Caching is useful if a sampler is used as a source for multiple modifiers. Without caching, the
duplicated input sources would redundantly compute the same outputs, which would be expensive,
especially if long pipelines share a long subgraph.
"""
function cache(x::S) where {N,S<:AbstractSampler{N}}
    Cache{N,S}(x.random_state, x, false, ntuple(i -> 0.0, N), 0.0)
end

@inline function sample(sampler::Cache{N}, coords::Vararg{Real,N}) where {N}
    if !sampler.is_cached || coords != sampler.coords
        sampler.is_cached = true
        sampler.coords = coords
        sampler.value = sample(sampler.source, coords...)
    end
    sampler.value
end
