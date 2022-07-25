"""
A 2-dimensional sampler that generates a pattern of concentric cylinders.
"""
struct Cylinders <: PatternSampler{2}
    random_state::RandomState
    frequency::Float64
end

"""
    Cylinders(; seed=nothing, frequency=1.0)

Construct a sampler that outputs values that form a pattern representing concentric cylinders when
it is sampled from.

# Arguments

  - `seed::Union{Int,Nothing}=nothing`: An integer used to seed the random number generator for this
    sampler, or `nothing`. If a seed is not supplied, one will be generated automatically which will
    negatively affect the reproducibility of any samplers further down the pipeline.

  - `frequency::Real=1.0`: The frequency of the signal, which controls how small or large the
    cylinders are.
"""
Cylinders(; seed::Seed=nothing, frequency::Real=1.0) = Cylinders(RandomState(seed), frequency)

function sample(sampler::S, x::T, y::T) where {S<:Cylinders,T<:Real}
    x, y = (x, y) .* sampler.frequency
    distance_center = sqrt(x^2 + y^2)
    distance_small = distance_center - floor(Int, distance_center)
    distance_large = 1 - distance_small
    nearest = min(distance_small, distance_large)
    1 - 4nearest
end
