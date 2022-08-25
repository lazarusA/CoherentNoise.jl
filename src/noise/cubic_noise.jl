struct Cubic{N} <: NoiseSampler{N}
    random_state::RandomState
end

@inline _cubic(dims, seed=0) = Cubic{dims}(RandomState(seed))

HashTrait(::Type{<:Cubic}) = IsValueHashed()

@inline function _hash(S::Type{<:Cubic}, seed, t, c1, c2, c3, c4)
    h1 = hash_coords(S, seed, c1...)
    h2 = hash_coords(S, seed, c2...)
    h3 = hash_coords(S, seed, c3...)
    h4 = hash_coords(S, seed, c4...)
    cubic_interpolate(h1, h2, h3, h4, t)
end

### 1D

@doc doc_cubic_1d
cubic_1d(; seed=0) = _cubic(1, seed)

function sample(sampler::S, x::Real) where {S<:Cubic{1}}
    seed = sampler.random_state.seed
    X = floor.(Int, x)
    X1 = X .* PRIME_X
    X2 = X1 + PRIME_X
    (_hash(S, seed, x - X, X1 - PRIME_X, X1, X2, X2 + PRIME_X * 2) - 1) / 1.5
end

### 2D

@doc doc_cubic_2d
cubic_2d(; seed=0) = _cubic(2, seed)

function sample(sampler::S, x::T, y::T) where {S<:Cubic{2},T<:Real}
    seed = sampler.random_state.seed
    primes = (PRIME_X, PRIME_Y)
    XY = floor.(Int, (x, y))
    X2, Y2 = XY .* primes
    X1, Y1 = (X2, Y2) .- primes
    X3, Y3 = (X2, Y2) .+ primes
    X4, Y4 = (X2, Y2) .+ (primes .* 2)
    x1, y1 = (x, y) .- XY
    c1 = _hash(S, seed, x1, (X1, Y1), (X2, Y1), (X3, Y1), (X4, Y1))
    c2 = _hash(S, seed, x1, (X1, Y2), (X2, Y2), (X3, Y2), (X4, Y2))
    c3 = _hash(S, seed, x1, (X1, Y3), (X2, Y3), (X3, Y3), (X4, Y3))
    c4 = _hash(S, seed, x1, (X1, Y4), (X2, Y4), (X3, Y4), (X4, Y4))
    (cubic_interpolate(c1, c2, c3, c4, y1) - 1) / 2.25
end

### 3D

@doc doc_cubic_3d
cubic_3d(; seed=0) = _cubic(3, seed)

function sample(sampler::S, x::T, y::T, z::T) where {S<:Cubic{3},T<:Real}
    seed = sampler.random_state.seed
    primes = (PRIME_X, PRIME_Y, PRIME_Z)
    XYZ = floor.(Int, (x, y, z))
    X2, Y2, Z2 = XYZ .* primes
    X1, Y1, Z1 = (X2, Y2, Z2) .- primes
    X3, Y3, Z3 = (X2, Y2, Z2) .+ primes
    X4, Y4, Z4 = (X2, Y2, Z2) .+ primes .* 2
    x1, y1, z1 = (x, y, z) .- XYZ
    c1 = _hash(S, seed, x1, (X1, Y1, Z1), (X2, Y1, Z1), (X3, Y1, Z1), (X4, Y1, Z1))
    c2 = _hash(S, seed, x1, (X1, Y2, Z1), (X2, Y2, Z1), (X3, Y2, Z1), (X4, Y2, Z1))
    c3 = _hash(S, seed, x1, (X1, Y3, Z1), (X2, Y3, Z1), (X3, Y3, Z1), (X4, Y3, Z1))
    c4 = _hash(S, seed, x1, (X1, Y4, Z1), (X2, Y4, Z1), (X3, Y4, Z1), (X4, Y4, Z1))
    c5 = _hash(S, seed, x1, (X1, Y1, Z2), (X2, Y1, Z2), (X3, Y1, Z2), (X4, Y1, Z2))
    c6 = _hash(S, seed, x1, (X1, Y2, Z2), (X2, Y2, Z2), (X3, Y2, Z2), (X4, Y2, Z2))
    c7 = _hash(S, seed, x1, (X1, Y3, Z2), (X2, Y3, Z2), (X3, Y3, Z2), (X4, Y3, Z2))
    c8 = _hash(S, seed, x1, (X1, Y4, Z2), (X2, Y4, Z2), (X3, Y4, Z2), (X4, Y4, Z2))
    c9 = _hash(S, seed, x1, (X1, Y1, Z3), (X2, Y1, Z3), (X3, Y1, Z3), (X4, Y1, Z3))
    c10 = _hash(S, seed, x1, (X1, Y2, Z3), (X2, Y2, Z3), (X3, Y2, Z3), (X4, Y2, Z3))
    c11 = _hash(S, seed, x1, (X1, Y3, Z3), (X2, Y3, Z3), (X3, Y3, Z3), (X4, Y3, Z3))
    c12 = _hash(S, seed, x1, (X1, Y4, Z3), (X2, Y4, Z3), (X3, Y4, Z3), (X4, Y4, Z3))
    c13 = _hash(S, seed, x1, (X1, Y1, Z4), (X2, Y1, Z4), (X3, Y1, Z4), (X4, Y1, Z4))
    c14 = _hash(S, seed, x1, (X1, Y2, Z4), (X2, Y2, Z4), (X3, Y2, Z4), (X4, Y2, Z4))
    c15 = _hash(S, seed, x1, (X1, Y3, Z4), (X2, Y3, Z4), (X3, Y3, Z4), (X4, Y3, Z4))
    c16 = _hash(S, seed, x1, (X1, Y4, Z4), (X2, Y4, Z4), (X3, Y4, Z4), (X4, Y4, Z4))
    r1 = cubic_interpolate(c1, c2, c3, c4, y1)
    r2 = cubic_interpolate(c5, c6, c7, c8, y1)
    r3 = cubic_interpolate(c9, c10, c11, c12, y1)
    r4 = cubic_interpolate(c13, c14, c15, c16, y1)
    (cubic_interpolate(r1, r2, r3, r4, z1) - 1) / 2.25
end

### 4D

@doc doc_cubic_4d
cubic_4d(; seed=0) = _cubic(4, seed)

function sample(sampler::S, x::T, y::T, z::T, w::T) where {S<:Cubic{4},T<:Real}
    seed = sampler.random_state.seed
    primes = (PRIME_X, PRIME_Y, PRIME_Z, PRIME_W)
    XYZW = floor.(Int, (x, y, z, w))
    X2, Y2, Z2, W2 = XYZW .* primes
    X1, Y1, Z1, W1 = (X2, Y2, Z2, W2) .- primes
    X3, Y3, Z3, W3 = (X2, Y2, Z2, W2) .+ primes
    X4, Y4, Z4, W4 = (X2, Y2, Z2, W2) .+ primes .* 2
    x1, y1, z1, w1 = (x, y, z, w) .- XYZW
    c1 = _hash(S, seed, x1, (X1, Y1, Z1, W1), (X2, Y1, Z1, W1), (X3, Y1, Z1, W1), (X4, Y1, Z1, W1))
    c2 = _hash(S, seed, x1, (X1, Y2, Z1, W1), (X2, Y2, Z1, W1), (X3, Y2, Z1, W1), (X4, Y2, Z1, W1))
    c3 = _hash(S, seed, x1, (X1, Y3, Z1, W1), (X2, Y3, Z1, W1), (X3, Y3, Z1, W1), (X4, Y3, Z1, W1))
    c4 = _hash(S, seed, x1, (X1, Y4, Z1, W1), (X2, Y4, Z1, W1), (X3, Y4, Z1, W1), (X4, Y4, Z1, W1))
    c5 = _hash(S, seed, x1, (X1, Y1, Z2, W1), (X2, Y1, Z2, W1), (X3, Y1, Z2, W1), (X4, Y1, Z2, W1))
    c6 = _hash(S, seed, x1, (X1, Y2, Z2, W1), (X2, Y2, Z2, W1), (X3, Y2, Z2, W1), (X4, Y2, Z2, W1))
    c7 = _hash(S, seed, x1, (X1, Y3, Z2, W1), (X2, Y3, Z2, W1), (X3, Y3, Z2, W1), (X4, Y3, Z2, W1))
    c8 = _hash(S, seed, x1, (X1, Y4, Z2, W1), (X2, Y4, Z2, W1), (X3, Y4, Z2, W1), (X4, Y4, Z2, W1))
    c9 = _hash(S, seed, x1, (X1, Y1, Z3, W1), (X2, Y1, Z3, W1), (X3, Y1, Z3, W1), (X4, Y1, Z3, W1))
    c10 = _hash(S, seed, x1, (X1, Y2, Z3, W1), (X2, Y2, Z3, W1), (X3, Y2, Z3, W1), (X4, Y2, Z3, W1))
    c11 = _hash(S, seed, x1, (X1, Y3, Z3, W1), (X2, Y3, Z3, W1), (X3, Y3, Z3, W1), (X4, Y3, Z3, W1))
    c12 = _hash(S, seed, x1, (X1, Y4, Z3, W1), (X2, Y4, Z3, W1), (X3, Y4, Z3, W1), (X4, Y4, Z3, W1))
    c13 = _hash(S, seed, x1, (X1, Y1, Z4, W1), (X2, Y1, Z4, W1), (X3, Y1, Z4, W1), (X4, Y1, Z4, W1))
    c14 = _hash(S, seed, x1, (X1, Y2, Z4, W1), (X2, Y2, Z4, W1), (X3, Y2, Z4, W1), (X4, Y2, Z4, W1))
    c15 = _hash(S, seed, x1, (X1, Y3, Z4, W1), (X2, Y3, Z4, W1), (X3, Y3, Z4, W1), (X4, Y3, Z4, W1))
    c16 = _hash(S, seed, x1, (X1, Y4, Z4, W1), (X2, Y4, Z4, W1), (X3, Y4, Z4, W1), (X4, Y4, Z4, W1))
    c17 = _hash(S, seed, x1, (X1, Y1, Z1, W2), (X2, Y1, Z1, W2), (X3, Y1, Z1, W2), (X4, Y1, Z1, W2))
    c18 = _hash(S, seed, x1, (X1, Y2, Z1, W2), (X2, Y2, Z1, W2), (X3, Y2, Z1, W2), (X4, Y2, Z1, W2))
    c19 = _hash(S, seed, x1, (X1, Y3, Z1, W2), (X2, Y3, Z1, W2), (X3, Y3, Z1, W2), (X4, Y3, Z1, W2))
    c20 = _hash(S, seed, x1, (X1, Y4, Z1, W2), (X2, Y4, Z1, W2), (X3, Y4, Z1, W2), (X4, Y4, Z1, W2))
    c21 = _hash(S, seed, x1, (X1, Y1, Z2, W2), (X2, Y1, Z2, W2), (X3, Y1, Z2, W2), (X4, Y1, Z2, W2))
    c22 = _hash(S, seed, x1, (X1, Y2, Z2, W2), (X2, Y2, Z2, W2), (X3, Y2, Z2, W2), (X4, Y2, Z2, W2))
    c23 = _hash(S, seed, x1, (X1, Y3, Z2, W2), (X2, Y3, Z2, W2), (X3, Y3, Z2, W2), (X4, Y3, Z2, W2))
    c24 = _hash(S, seed, x1, (X1, Y4, Z2, W2), (X2, Y4, Z2, W2), (X3, Y4, Z2, W2), (X4, Y4, Z2, W2))
    c25 = _hash(S, seed, x1, (X1, Y1, Z3, W2), (X2, Y1, Z3, W2), (X3, Y1, Z3, W2), (X4, Y1, Z3, W2))
    c26 = _hash(S, seed, x1, (X1, Y2, Z3, W2), (X2, Y2, Z3, W2), (X3, Y2, Z3, W2), (X4, Y2, Z3, W2))
    c27 = _hash(S, seed, x1, (X1, Y3, Z3, W2), (X2, Y3, Z3, W2), (X3, Y3, Z3, W2), (X4, Y3, Z3, W2))
    c28 = _hash(S, seed, x1, (X1, Y4, Z3, W2), (X2, Y4, Z3, W2), (X3, Y4, Z3, W2), (X4, Y4, Z3, W2))
    c29 = _hash(S, seed, x1, (X1, Y1, Z4, W2), (X2, Y1, Z4, W2), (X3, Y1, Z4, W2), (X4, Y1, Z4, W2))
    c30 = _hash(S, seed, x1, (X1, Y2, Z4, W2), (X2, Y2, Z4, W2), (X3, Y2, Z4, W2), (X4, Y2, Z4, W2))
    c31 = _hash(S, seed, x1, (X1, Y3, Z4, W2), (X2, Y3, Z4, W2), (X3, Y3, Z4, W2), (X4, Y3, Z4, W2))
    c32 = _hash(S, seed, x1, (X1, Y4, Z4, W2), (X2, Y4, Z4, W2), (X3, Y4, Z4, W2), (X4, Y4, Z4, W2))
    c33 = _hash(S, seed, x1, (X1, Y1, Z1, W3), (X2, Y1, Z1, W3), (X3, Y1, Z1, W3), (X4, Y1, Z1, W3))
    c34 = _hash(S, seed, x1, (X1, Y2, Z1, W3), (X2, Y2, Z1, W3), (X3, Y2, Z1, W3), (X4, Y2, Z1, W3))
    c35 = _hash(S, seed, x1, (X1, Y3, Z1, W3), (X2, Y3, Z1, W3), (X3, Y3, Z1, W3), (X4, Y3, Z1, W3))
    c36 = _hash(S, seed, x1, (X1, Y4, Z1, W3), (X2, Y4, Z1, W3), (X3, Y4, Z1, W3), (X4, Y4, Z1, W3))
    c37 = _hash(S, seed, x1, (X1, Y1, Z2, W3), (X2, Y1, Z2, W3), (X3, Y1, Z2, W3), (X4, Y1, Z2, W3))
    c38 = _hash(S, seed, x1, (X1, Y2, Z2, W3), (X2, Y2, Z2, W3), (X3, Y2, Z2, W3), (X4, Y2, Z2, W3))
    c39 = _hash(S, seed, x1, (X1, Y3, Z2, W3), (X2, Y3, Z2, W3), (X3, Y3, Z2, W3), (X4, Y3, Z2, W3))
    c40 = _hash(S, seed, x1, (X1, Y4, Z2, W3), (X2, Y4, Z2, W3), (X3, Y4, Z2, W3), (X4, Y4, Z2, W3))
    c41 = _hash(S, seed, x1, (X1, Y1, Z3, W3), (X2, Y1, Z3, W3), (X3, Y1, Z3, W3), (X4, Y1, Z3, W3))
    c42 = _hash(S, seed, x1, (X1, Y2, Z3, W3), (X2, Y2, Z3, W3), (X3, Y2, Z3, W3), (X4, Y2, Z3, W3))
    c43 = _hash(S, seed, x1, (X1, Y3, Z3, W3), (X2, Y3, Z3, W3), (X3, Y3, Z3, W3), (X4, Y3, Z3, W3))
    c44 = _hash(S, seed, x1, (X1, Y4, Z3, W3), (X2, Y4, Z3, W3), (X3, Y4, Z3, W3), (X4, Y4, Z3, W3))
    c45 = _hash(S, seed, x1, (X1, Y1, Z4, W3), (X2, Y1, Z4, W3), (X3, Y1, Z4, W3), (X4, Y1, Z4, W3))
    c46 = _hash(S, seed, x1, (X1, Y2, Z4, W3), (X2, Y2, Z4, W3), (X3, Y2, Z4, W3), (X4, Y2, Z4, W3))
    c47 = _hash(S, seed, x1, (X1, Y3, Z4, W3), (X2, Y3, Z4, W3), (X3, Y3, Z4, W3), (X4, Y3, Z4, W3))
    c48 = _hash(S, seed, x1, (X1, Y4, Z4, W3), (X2, Y4, Z4, W3), (X3, Y4, Z4, W3), (X4, Y4, Z4, W3))
    c49 = _hash(S, seed, x1, (X1, Y1, Z1, W4), (X2, Y1, Z1, W4), (X3, Y1, Z1, W4), (X4, Y1, Z1, W4))
    c50 = _hash(S, seed, x1, (X1, Y2, Z1, W4), (X2, Y2, Z1, W4), (X3, Y2, Z1, W4), (X4, Y2, Z1, W4))
    c51 = _hash(S, seed, x1, (X1, Y3, Z1, W4), (X2, Y3, Z1, W4), (X3, Y3, Z1, W4), (X4, Y3, Z1, W4))
    c52 = _hash(S, seed, x1, (X1, Y4, Z1, W4), (X2, Y4, Z1, W4), (X3, Y4, Z1, W4), (X4, Y4, Z1, W4))
    c53 = _hash(S, seed, x1, (X1, Y1, Z2, W4), (X2, Y1, Z2, W4), (X3, Y1, Z2, W4), (X4, Y1, Z2, W4))
    c54 = _hash(S, seed, x1, (X1, Y2, Z2, W4), (X2, Y2, Z2, W4), (X3, Y2, Z2, W4), (X4, Y2, Z2, W4))
    c55 = _hash(S, seed, x1, (X1, Y3, Z2, W4), (X2, Y3, Z2, W4), (X3, Y3, Z2, W4), (X4, Y3, Z2, W4))
    c56 = _hash(S, seed, x1, (X1, Y4, Z2, W4), (X2, Y4, Z2, W4), (X3, Y4, Z2, W4), (X4, Y4, Z2, W4))
    c57 = _hash(S, seed, x1, (X1, Y1, Z3, W4), (X2, Y1, Z3, W4), (X3, Y1, Z3, W4), (X4, Y1, Z3, W4))
    c58 = _hash(S, seed, x1, (X1, Y2, Z3, W4), (X2, Y2, Z3, W4), (X3, Y2, Z3, W4), (X4, Y2, Z3, W4))
    c59 = _hash(S, seed, x1, (X1, Y3, Z3, W4), (X2, Y3, Z3, W4), (X3, Y3, Z3, W4), (X4, Y3, Z3, W4))
    c60 = _hash(S, seed, x1, (X1, Y4, Z3, W4), (X2, Y4, Z3, W4), (X3, Y4, Z3, W4), (X4, Y4, Z3, W4))
    c61 = _hash(S, seed, x1, (X1, Y1, Z4, W4), (X2, Y1, Z4, W4), (X3, Y1, Z4, W4), (X4, Y1, Z4, W4))
    c62 = _hash(S, seed, x1, (X1, Y2, Z4, W4), (X2, Y2, Z4, W4), (X3, Y2, Z4, W4), (X4, Y2, Z4, W4))
    c63 = _hash(S, seed, x1, (X1, Y3, Z4, W4), (X2, Y3, Z4, W4), (X3, Y3, Z4, W4), (X4, Y3, Z4, W4))
    c64 = _hash(S, seed, x1, (X1, Y4, Z4, W4), (X2, Y4, Z4, W4), (X3, Y4, Z4, W4), (X4, Y4, Z4, W4))
    r1a = cubic_interpolate(c1, c2, c3, c4, y1)
    r1b = cubic_interpolate(c5, c6, c7, c8, y1)
    r1c = cubic_interpolate(c9, c10, c11, c12, y1)
    r1d = cubic_interpolate(c13, c14, c15, c16, y1)
    r2a = cubic_interpolate(c17, c18, c19, c20, y1)
    r2b = cubic_interpolate(c21, c22, c23, c24, y1)
    r2c = cubic_interpolate(c25, c26, c27, c28, y1)
    r2d = cubic_interpolate(c29, c30, c31, c32, y1)
    r3a = cubic_interpolate(c33, c34, c35, c36, y1)
    r3b = cubic_interpolate(c37, c38, c39, c40, y1)
    r3c = cubic_interpolate(c41, c42, c43, c44, y1)
    r3d = cubic_interpolate(c45, c46, c47, c48, y1)
    r4a = cubic_interpolate(c49, c50, c51, c52, y1)
    r4b = cubic_interpolate(c53, c54, c55, c56, y1)
    r4c = cubic_interpolate(c57, c58, c59, c60, y1)
    r4d = cubic_interpolate(c61, c62, c63, c64, y1)
    r1 = cubic_interpolate(r1a, r1b, r1c, r1d, z1)
    r2 = cubic_interpolate(r2a, r2b, r2c, r2d, z1)
    r3 = cubic_interpolate(r3a, r3b, r3c, r3d, z1)
    r4 = cubic_interpolate(r4a, r4b, r4c, r4d, z1)
    (cubic_interpolate(r1, r2, r3, r4, w1) - 1) / 2.25
end
