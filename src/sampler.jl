"""
Supertype for all samplers.
"""
abstract type AbstractSampler{N} end

"""
    sample(sampler::AbstractSampler, x::Real)
    sample(sampler::AbstractSampler, x::Real, y::Real)
    sample(sampler::AbstractSampler, x::Real, y::Real, z::Real)
    sample(sampler::AbstractSampler, x::Real, y::Real, z::Real, w::Real)

Sample from `sampler` with the supplied coordinates. The number of coordinates should match the
dimensionality of the sampler type.
"""
function sample end
