### Constant

struct Constant <: PatternSampler{1}
    random_state::RandomState
    value::Float64
end

@doc doc_constant_1d
constant_1d(; seed=0, value=0.0) = Constant(RandomState(seed), value)

sample(sampler::S, _) where {S<:Constant} = sampler.value

### Checkered

struct Checkered <: PatternSampler{2}
    random_state::RandomState
end

@doc doc_checkered_2d
checkered_2d(; seed=0) = Checkered(RandomState(seed))

function sample(::S, x::T, y::T) where {S<:Checkered,T<:Real}
    iszero((floor(Int, x) & 1) ⊻ (floor(Int, y) & 1)) ? 1.0 : -1.0
end

### Cylinders

struct Cylinders <: PatternSampler{2}
    random_state::RandomState
    frequency::Float64
end

@doc doc_cylinders_2d
cylinders_2d(; seed=0, frequency=1.0) = Cylinders(RandomState(seed), frequency)

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

@doc doc_spheres_3d
spheres_3d(; seed=0, frequency=1.0) = Spheres(RandomState(seed), frequency)

function sample(sampler::S, x::T, y::T, z::T) where {S<:Spheres,T<:Real}
    x, y, z = (x, y, z) .* sampler.frequency
    distance_center = sqrt(x^2 + y^2 + z^2)
    distance_small = distance_center - floor(Int, distance_center)
    distance_large = 1 - distance_small
    nearest = min(distance_small, distance_large)
    1 - 4nearest
end
