struct Simplex{N} <: NoiseSampler{N}
    random_state::RandomState
    perlin_state::PerlinState
    simplex_state::SimplexState
end

@inline function _simplex(dims, seed, smooth)
    T = Simplex{dims}
    rs = RandomState(seed)
    T(rs, PerlinState(rs), SimplexState(T, Val(smooth)))
end

HashTrait(::Type{<:Simplex}) = IsPerlinHashed()

SimplexState(::Type{<:Simplex{1}}, ::Val) = SimplexState(1.0, 0.395)
SimplexState(::Type{<:Simplex{2}}, ::Val) = SimplexState(0.5, 45.23065)
SimplexState(::Type{<:Simplex{3}}, ::Val{true}) = SimplexState(0.5, 76.88075002223152)
SimplexState(::Type{<:Simplex{3}}, ::Val{false}) = SimplexState(0.6, 32.69428328944204)
SimplexState(::Type{<:Simplex{4}}, ::Val{true}) = SimplexState(0.5, 62.0)
SimplexState(::Type{<:Simplex{4}}, ::Val{false}) = SimplexState(0.6, 27.0)

# 1D

"""
    simplex_1d(; seed=nothing)

Construct a sampler that outputs 1-dimensional Perlin Simplex noise when it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.
"""
simplex_1d(; seed=nothing) = _simplex(1, seed, false)

@inline function grad(::Type{Simplex{1}}, falloff, hash, x)
    s = falloff - x^2
    h = hash & 15
    u = (h & 7) + 1
    g = iszero(h & 8) ? -u * x : u * x
    pow4(s) * g
end

function sample(sampler::S, x::Real) where {S<:Simplex{1}}
    t = sampler.perlin_state.table
    state = sampler.simplex_state
    falloff = state.falloff
    X = floor(Int, x)
    X1 = X & 255 + 1
    x1 = x - X
    @inbounds begin
        p1 = grad(S, falloff, t[X1], x1)
        p2 = grad(S, falloff, t[X1+1], x1 - 1)
    end
    (p1 + p2) * state.scale_factor
end

# 2D

const SIMPLEX_SKEW_2D = (sqrt(3) - 1) / 2
const SIMPLEX_UNSKEW_2D = (3 - sqrt(3)) / 6

"""
    simplex_2d(; seed=nothing)

Construct a sampler that outputs 2-dimensional Perlin Simplex noise when it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.
"""
simplex_2d(; seed=nothing) = _simplex(2, seed, false)

@inline function grad(::Type{Simplex{2}}, falloff, hash, x, y)
    s = falloff - x^2 - y^2
    h = hash & 7
    u, v = h < 4 ? (x, y) : (y, x)
    g = (iszero(h & 1) ? -u : u) + (iszero(h & 2) ? -2v : 2v)
    s > 0 ? pow4(s) * g : 0.0
end

@inline get_simplex(::Type{Simplex{2}}, x, y) = x > y ? (1, 0) : (0, 1)

function sample(sampler::S, x::T, y::T) where {S<:Simplex{2},T<:Real}
    t = sampler.perlin_state.table
    state = sampler.simplex_state
    falloff = state.falloff
    s = (x + y) * SIMPLEX_SKEW_2D
    X, Y = floor.(Int, (x, y) .+ s)
    X1, Y1 = (X, Y) .& 255 .+ 1
    tx = (X + Y) * SIMPLEX_UNSKEW_2D
    xy = (x, y) .- (X, Y) .+ tx
    X2, Y2 = get_simplex(S, xy...)
    @inbounds begin
        p1 = grad(S, falloff, t[t[X1]+Y1], xy...)
        p2 = grad(S, falloff, t[t[X1+X2]+Y1+Y2], xy .- (X2, Y2) .+ SIMPLEX_UNSKEW_2D...)
        p3 = grad(S, falloff, t[t[X1+1]+Y1+1], xy .- 1 .+ 2SIMPLEX_UNSKEW_2D...)
    end
    (p1 + p2 + p3) * state.scale_factor
end

# 3D

const SIMPLEX_SKEW_3D = 1 / 3
const SIMPLEX_UNSKEW_3D = 1 / 6

"""
    simplex_3d(; seed=nothing, smooth=false)

Construct a sampler that outputs 3-dimensional Perlin Simplex noise when it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

  - `smooth`: Specify whether to have continuous gradients.
    Simplex variants, even the original Simplex noise by Ken Perlin, overshoot the radial extent for
    the signal reconstruction kernel in order to improve the visual of the noise. Normally this is
    okay, especially if layering multiple octaves of the noise. However, in some applications, such
    as creating height or bump maps, this will produce discontinuities visually identified by
    jarring creases in the generated noise.

    This option changes the falloff in order to produce smooth continuous noise, however, the
    resulting noise may look quite different than the non-smooth option, depending on the Simplex
    variant.

    The default value is `false`, in order to be true to the original implementation.
"""
simplex_3d(; seed=nothing, smooth=false) = _simplex(3, seed, smooth)

@inline function grad(S::Type{Simplex{3}}, falloff, hash, x, y, z)
    s = falloff - x^2 - y^2 - z^2
    h = hash & 15
    u = h < 8 ? x : y
    v = h < 4 ? y : h == 12 || h == 14 ? x : z
    g = hash_coords(S, h, u, v)
    s > 0 ? pow4(s) * g : 0.0
end

@inline function get_simplex(::Type{Simplex{3}}, x, y, z)
    if x ≥ y
        y ≥ z ? (1, 0, 0, 1, 1, 0) : x ≥ z ? (1, 0, 0, 1, 0, 1) : (0, 0, 1, 1, 0, 1)
    else
        y < z ? (0, 0, 1, 0, 1, 1) : x < z ? (0, 1, 0, 0, 1, 1) : (0, 1, 0, 1, 1, 0)
    end
end

function sample(sampler::S, x::T, y::T, z::T) where {S<:Simplex{3},T<:Real}
    t = sampler.perlin_state.table
    state = sampler.simplex_state
    falloff = state.falloff
    s = (x + y + z) * SIMPLEX_SKEW_3D
    X, Y, Z = floor.(Int, (x, y, z) .+ s)
    X1, Y1, Z1 = (X, Y, Z) .& 255 .+ 1
    tx = (X + Y + Z) * SIMPLEX_UNSKEW_3D
    xyz = (x, y, z) .- (X, Y, Z) .+ tx
    X2, Y2, Z2, X3, Y3, Z3 = get_simplex(S, xyz...)
    @inbounds begin
        hash1 = t[t[t[Z1]+Y1]+X1]
        hash2 = t[t[t[Z1+Z2]+Y1+Y2]+X1+X2]
        hash3 = t[t[t[Z1+Z3]+Y1+Y3]+X1+X3]
        hash4 = t[t[t[Z1+1]+Y1+1]+X1+1]
    end
    p1 = grad(S, falloff, hash1, xyz...)
    p2 = grad(S, falloff, hash2, xyz .- (X2, Y2, Z2) .+ SIMPLEX_UNSKEW_3D...)
    p3 = grad(S, falloff, hash3, xyz .- (X3, Y3, Z3) .+ 2SIMPLEX_UNSKEW_3D...)
    p4 = grad(S, falloff, hash4, xyz .- 1 .+ 3SIMPLEX_UNSKEW_3D...)
    (p1 + p2 + p3 + p4) * state.scale_factor
end

# 4D

const SIMPLEX_SKEW_4D = (sqrt(5) - 1) / 4
const SIMPLEX_UNSKEW_4D = (5 - sqrt(5)) / 20
const SIMPLEX_GRADIENTS_4D = [
    0x0, 0x1, 0x2, 0x3, 0x0, 0x1, 0x3, 0x2, 0x0, 0x0, 0x0, 0x0, 0x0, 0x2, 0x3, 0x1, 0x0, 0x0,
    0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x1, 0x2, 0x3, 0x0, 0x0, 0x2, 0x1, 0x3,
    0x0, 0x0, 0x0, 0x0, 0x0, 0x3, 0x1, 0x2, 0x0, 0x3, 0x2, 0x1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
    0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x1, 0x3, 0x2, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
    0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
    0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x1, 0x2, 0x0, 0x3, 0x0, 0x0, 0x0, 0x0, 0x1, 0x3, 0x0, 0x2,
    0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x2, 0x3, 0x0, 0x1, 0x2, 0x3,
    0x1, 0x0, 0x1, 0x0, 0x2, 0x3, 0x1, 0x0, 0x3, 0x2, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
    0x0, 0x0, 0x0, 0x0, 0x2, 0x0, 0x3, 0x1, 0x0, 0x0, 0x0, 0x0, 0x2, 0x1, 0x3, 0x0, 0x0, 0x0,
    0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
    0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x2, 0x0, 0x1, 0x3, 0x0, 0x0,
    0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x3, 0x0, 0x1, 0x2, 0x3, 0x0, 0x2, 0x1,
    0x0, 0x0, 0x0, 0x0, 0x3, 0x1, 0x2, 0x0, 0x2, 0x1, 0x0, 0x3, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
    0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x3, 0x1, 0x0, 0x2, 0x0, 0x0, 0x0, 0x0, 0x3, 0x2, 0x0, 0x1,
    0x3, 0x2, 0x1, 0x0]

"""
    simplex_4d(; seed=nothing, smooth=false)

Construct a sampler that outputs 4-dimensional Perlin Simplex noise when it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

  - `smooth`: Specify whether to have continuous gradients.
    Simplex variants, even the original Simplex noise by Ken Perlin, overshoot the radial extent for
    the signal reconstruction kernel in order to improve the visual of the noise. Normally this is
    okay, especially if layering multiple octaves of the noise. However, in some applications, such
    as creating height or bump maps, this will produce discontinuities visually identified by
    jarring creases in the generated noise.

    This option changes the falloff in order to produce smooth continuous noise, however, the
    resulting noise may look quite different than the non-smooth option, depending on the Simplex
    variant.

    The default value is `false`, in order to be true to the original implementation.
"""
simplex_4d(; seed=nothing, smooth=false) = _simplex(4, seed, smooth)

@inline function grad(S::Type{Simplex{4}}, falloff, hash, x, y, z, w)
    s = falloff - x^2 - y^2 - z^2 - w^2
    h = hash & 31
    u = h < 24 ? x : y
    v = h < 16 ? y : z
    w = h < 8 ? z : w
    g = hash_coords(S, h, u, v, w)
    s > 0 ? pow4(s) * g : 0.0
end

@inline function get_simplex(::Type{Simplex{4}}, x, y, z, w)
    t = SIMPLEX_GRADIENTS_4D
    c1 = x > y ? 32 : 0
    c2 = x > z ? 16 : 0
    c3 = y > z ? 8 : 0
    c4 = x > w ? 4 : 0
    c5 = y > w ? 2 : 0
    c6 = z > w ? 1 : 0
    c = (c1 + c2 + c3 + c4 + c5 + c6) * 4
    a1, a2, a3, a4 = @inbounds (t[c+i] for i in 1:4)
    (
        a1 ≥ 3 ? 1 : 0, a2 ≥ 3 ? 1 : 0, a3 ≥ 3 ? 1 : 0, a4 ≥ 3 ? 1 : 0,
        a1 ≥ 2 ? 1 : 0, a2 ≥ 2 ? 1 : 0, a3 ≥ 2 ? 1 : 0, a4 ≥ 2 ? 1 : 0,
        a1 ≥ 1 ? 1 : 0, a2 ≥ 1 ? 1 : 0, a3 ≥ 1 ? 1 : 0, a4 ≥ 1 ? 1 : 0,
    )
end

function sample(sampler::S, x::T, y::T, z::T, w::T) where {S<:Simplex{4},T<:Real}
    t = sampler.perlin_state.table
    state = sampler.simplex_state
    falloff = state.falloff
    s = (x + y + z + w) * SIMPLEX_SKEW_4D
    X, Y, Z, W = floor.(Int, (x, y, z, w) .+ s)
    X1, Y1, Z1, W1 = (X, Y, Z, W) .& 255 .+ 1
    tx = (X + Y + Z + W) * SIMPLEX_UNSKEW_4D
    v1 = (x, y, z, w) .- (X, Y, Z, W) .+ tx
    X2, Y2, Z2, W2, X3, Y3, Z3, W3, X4, Y4, Z4, W4 = get_simplex(S, v1...)
    v2 = v1 .- (X2, Y2, Z2, W2) .+ SIMPLEX_UNSKEW_4D
    v3 = v1 .- (X3, Y3, Z3, W3) .+ 2SIMPLEX_UNSKEW_4D
    v4 = v1 .- (X4, Y4, Z4, W4) .+ 3SIMPLEX_UNSKEW_4D
    v5 = v1 .- 1 .+ 4SIMPLEX_UNSKEW_4D
    p1 = grad(S, falloff, t[t[t[t[W1]+Z1]+Y1]+X1], v1...)
    p2 = grad(S, falloff, t[t[t[t[W1+W2]+Z1+Z2]+Y1+Y2]+X1+X2], v2...)
    p3 = grad(S, falloff, t[t[t[t[W1+W3]+Z1+Z3]+Y1+Y3]+X1+X3], v3...)
    p4 = grad(S, falloff, t[t[t[t[W1+W4]+Z1+Z4]+Y1+Y4]+X1+X4], v4...)
    p5 = grad(S, falloff, t[t[t[t[W1+1]+Z1+1]+Y1+1]+X1+1], v5...)
    (p1 + p2 + p3 + p4 + p5) * state.scale_factor
end
