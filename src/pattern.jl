### Constant

struct Constant <: PatternSampler{1}
    random_state::RandomState
    value::Float64
end

"""
    constant_1d(; seed=nothing, value=0.0)

Construct a sampler that constantly outputs `value` each time it is sampled from.

This is useful for debugging and applications where you need to combine a constant value.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

  - `value`: A constant value to emit each time this sampler is sampled from.
"""
constant_1d(; seed=nothing, value=0.0) = Constant(RandomState(seed), value)

sample(sampler::S, _) where {S<:Constant} = sampler.value

### Checkered

struct Checkered <: PatternSampler{2}
    random_state::RandomState
end

"""
    checkered_2d(; seed=nothing)

Construct a sampler that outputs values in a checkerboard-like pattern when it is sampled from.

That is, output values will only ever be -1.0 or 1.0.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.
"""
checkered_2d(; seed=nothing) = Checkered(RandomState(seed))

function sample(::S, x::T, y::T) where {S<:Checkered,T<:Real}
    iszero((floor(Int, x) & 1) ⊻ (floor(Int, y) & 1)) ? 1.0 : -1.0
end

### Cylinders

struct Cylinders <: PatternSampler{2}
    random_state::RandomState
    frequency::Float64
end

"""
    cylinders_2d(; seed=nothing, frequency=1.0)

Construct a sampler that outputs values that form a pattern representing concentric cylinders when
it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

  - `frequency`: The frequency of the signal, which controls how small or large the cylinders
    are.
"""
cylinders_2d(; seed=nothing, frequency=1.0) = Cylinders(RandomState(seed), frequency)

function sample(sampler::S, x::T, y::T) where {S<:Cylinders,T<:Real}
    x, y = (x, y) .* sampler.frequency
    distance_center = sqrt(x^2 + y^2)
    distance_small = distance_center - floor(Int, distance_center)
    distance_large = 1 - distance_small
    nearest = min(distance_small, distance_large)
    1 - 4nearest
end

### Spheres

struct Spheres <: PatternSampler{3}
    random_state::RandomState
    frequency::Float64
end

"""
    spheres_3d(; seed=nothing, frequency=1.0)

Construct a sampler that outputs values that form a pattern representing concentric spheres when it
is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

  - `frequency`: The frequency of the signal, which controls how small or large the spheres are.
"""
spheres_3d(; seed=nothing, frequency=1.0) = Spheres(RandomState(seed), frequency)

function sample(sampler::S, x::T, y::T, z::T) where {S<:Spheres,T<:Real}
    x, y, z = (x, y, z) .* sampler.frequency
    distance_center = sqrt(x^2 + y^2 + z^2)
    distance_small = distance_center - floor(Int, distance_center)
    distance_large = 1 - distance_small
    nearest = min(distance_small, distance_large)
    1 - 4nearest
end
