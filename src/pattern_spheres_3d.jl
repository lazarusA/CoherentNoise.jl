struct Spheres <: PatternSampler{3}
    random_state::RandomState
    frequency::Float64
end

"""
    spheres_3d(; kwargs...)

Construct a sampler that outputs values that form a pattern representing concentric spheres when it
is sampled from.

# Arguments

  - `seed=0`: An integer used to seed the random number generator for this sampler.

  - `frequency=1.0`: The frequency of the signal, which controls how small or large the spheres are.
"""
spheres_3d(; seed=0, frequency=1.0) = Spheres(RandomState(seed), frequency)

function sample(sampler::S, x::T, y::T, z::T) where {S<:Spheres,T<:Real}
    x, y, z = (x, y, z) .* sampler.frequency
    distance_center = sqrt(x^2 + y^2 + z^2)
    distance_small = distance_center - floor(Int, distance_center)
    distance_large = 1 - distance_small
    nearest = min(distance_small, distance_large)
    1 - 4nearest
end
