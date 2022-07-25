"""
A 3-dimensional sampler that generates a pattern of concentric spheres.
"""
struct Spheres <: PatternSampler{3}
    random_state::RandomState
    frequency::Float64
end

"""
    Spheres(; seed=nothing, frequency=1.0)

Construct a sampler that outputs values that form a pattern representing concentric spheres when it
is sampled from.

# Arguments

  - `seed::Union{Int,Nothing}=nothing`: An integer used to seed the random number generator for this
    sampler, or `nothing`. If a seed is not supplied, one will be generated automatically which will
    negatively affect the reproducibility of any samplers further down the pipeline.

  - `frequency::Real=1.0`: The frequency of the signal, which controls how small or large the
    spheres are.
"""
Spheres(; seed::Seed=nothing, frequency::Real=1.0) = Spheres(RandomState(seed), frequency)

function sample(sampler::S, x::T, y::T, z::T) where {S<:Spheres,T<:Real}
    x, y, z = (x, y, z) .* sampler.frequency
    distance_center = sqrt(x^2 + y^2 + z^2)
    distance_small = distance_center - floor(Int, distance_center)
    distance_large = 1 - distance_small
    nearest = min(distance_small, distance_large)
    1 - 4nearest
end
