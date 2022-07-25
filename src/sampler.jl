"""
Supertype for all `N`-dimensional samplers.
"""
abstract type AbstractSampler{N} end

"""
    sample(sampler::AbstractSampler, x)
    sample(sampler::AbstractSampler, x, y)
    sample(sampler::AbstractSampler, x, y, z)
    sample(sampler::AbstractSampler, x, y, z, w)

Sample from `sampler` with the supplied coordinates. The number of coordinates should match the
dimensionality of the sampler type.
"""
function sample end
