struct Value{N} <: NoiseSampler{N}
    random_state::RandomState
end

@inline _value(dims, seed=0) = Value{dims}(RandomState(seed))

HashTrait(::Type{<:Value}) = IsValueHashed()

### 2D

@doc doc_value_2d
value_2d(; seed=0) = _value(2, seed)

function sample(sampler::S, x::T, y::T) where {S<:Value{2},T<:Real}
    seed = sampler.random_state.seed
    primes = (PRIME_X, PRIME_Y)
    X, Y = floor.(Int, (x, y))
    X1, Y1 = (X, Y) .* primes
    X2, Y2 = (X1, Y1) .+ primes
    x1, y1 = curve3.((x, y) .- (X, Y))
    p1 = lerp(hash_coords(S, seed, X1, Y1), hash_coords(S, seed, X2, Y1), x1)
    p2 = lerp(hash_coords(S, seed, X1, Y2), hash_coords(S, seed, X2, Y2), x1)
    lerp(p1, p2, y1) - 1
end

### 3D

@doc doc_value_3d
value_3d(; seed=0) = _value(3, seed)

function sample(sampler::S, x::T, y::T, z::T) where {S<:Value{3},T<:Real}
    seed = sampler.random_state.seed
    primes = (PRIME_X, PRIME_Y, PRIME_Z)
    X, Y, Z = floor.(Int, (x, y, z))
    X1, Y1, Z1 = (X, Y, Z) .* primes
    X2, Y2, Z2 = (X1, Y1, Z1) .+ primes
    x1, y1, z1 = curve3.((x, y, z) .- (X, Y, Z))
    p1 = lerp(hash_coords(S, seed, X1, Y1, Z1), hash_coords(S, seed, X2, Y1, Z1), x1)
    p2 = lerp(hash_coords(S, seed, X1, Y2, Z1), hash_coords(S, seed, X2, Y2, Z1), x1)
    p3 = lerp(hash_coords(S, seed, X1, Y1, Z2), hash_coords(S, seed, X2, Y1, Z2), x1)
    p4 = lerp(hash_coords(S, seed, X1, Y2, Z2), hash_coords(S, seed, X2, Y2, Z2), x1)
    lerp(lerp(p1, p2, y1), lerp(p3, p4, y1), z1) - 1
end

### 4D

@doc doc_value_4d
value_4d(; seed=0) = _value(4, seed)

function sample(sampler::S, x::T, y::T, z::T, w::T) where {S<:Value{4},T<:Real}
    seed = sampler.random_state.seed
    primes = (PRIME_X, PRIME_Y, PRIME_Z, PRIME_W)
    X, Y, Z, W = floor.(Int, (x, y, z, w))
    X1, Y1, Z1, W1 = (X, Y, Z, W) .* primes
    X2, Y2, Z2, W2 = (X1, Y1, Z1, W1) .+ primes
    x1, y1, z1, w1 = curve3.((x, y, z, w) .- (X, Y, Z, W))
    p1 = lerp(hash_coords(S, seed, X1, Y1, Z1, W1), hash_coords(S, seed, X2, Y1, Z1, W1), x1)
    p2 = lerp(hash_coords(S, seed, X1, Y2, Z1, W1), hash_coords(S, seed, X2, Y2, Z1, W1), x1)
    p3 = lerp(hash_coords(S, seed, X1, Y1, Z2, W1), hash_coords(S, seed, X2, Y1, Z2, W1), x1)
    p4 = lerp(hash_coords(S, seed, X1, Y2, Z2, W1), hash_coords(S, seed, X2, Y2, Z2, W1), x1)
    p5 = lerp(hash_coords(S, seed, X1, Y1, Z1, W2), hash_coords(S, seed, X2, Y1, Z1, W2), x1)
    p6 = lerp(hash_coords(S, seed, X1, Y2, Z1, W2), hash_coords(S, seed, X2, Y2, Z1, W2), x1)
    p7 = lerp(hash_coords(S, seed, X1, Y1, Z2, W2), hash_coords(S, seed, X2, Y1, Z2, W2), x1)
    p8 = lerp(hash_coords(S, seed, X1, Y2, Z2, W2), hash_coords(S, seed, X2, Y2, Z2, W2), x1)
    p9 = lerp(lerp(p1, p2, y1), lerp(p3, p4, y1), z1)
    p10 = lerp(lerp(p5, p6, y1), lerp(p7, p8, y1), z1)
    lerp(p9, p10, w1) - 1
end
