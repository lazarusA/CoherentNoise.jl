struct Simplex{N} <: NoiseSampler{N}
    random_state::RandomState
    state::PerlinState
end

@inline function _simplex(dims, seed)
    rs = RandomState(seed)
    Simplex{dims}(rs, PerlinState(rs))
end

HashTrait(::Type{<:Simplex}) = IsPerlinHashed()

# 1D

@doc doc_simplex_1d
simplex_1d(; seed=0) = _simplex(1, seed)

@inline function grad(::Type{Simplex{1}}, hash, x)
    s = 1 - x^2
    h = hash & 15
    u = (h & 7) + 1
    g = iszero(h & 8) ? u * x : -u * x
    s > 0 ? pow4(s) * g : 0.0
end

function sample(sampler::S, x::Real) where {S<:Simplex{1}}
    t = sampler.state.table
    X = floor(Int, x)
    x1 = x - X
    p1 = grad(S, t[X], x1)
    p2 = grad(S, t[X+1], x1 - 1)
    (p1 + p2) * 0.395
end

# 2D

const SIMPLEX_SKEW_2D = (sqrt(3) - 1) / 2
const SIMPLEX_UNSKEW_2D = (3 - sqrt(3)) / 6

@doc doc_simplex_2d
simplex_2d(; seed=0) = _simplex(2, seed)

@inline function grad(::Type{Simplex{2}}, hash, x, y)
    s = 0.5 - x^2 - y^2
    h = hash & 7
    u, v = h < 4 ? (x, y) : (y, x)
    g = (iszero(h & 1) ? u : -u) + (iszero(h & 2) ? 2v : -2v)
    s > 0 ? pow4(s) * g : 0.0
end

@inline get_simplex(::Type{Simplex{2}}, x, y) = x > y ? (1, 0) : (0, 1)

function sample(sampler::S, x::T, y::T) where {S<:Simplex{2},T<:Real}
    t = sampler.state.table
    s = (x + y) * SIMPLEX_SKEW_2D
    X, Y = floor.(Int, (x, y) .+ s)
    tx = (X + Y) * SIMPLEX_UNSKEW_2D
    xy = (x, y) .- (X, Y) .+ tx
    X1, Y1 = get_simplex(S, xy...)
    p1 = grad(S, t[t[X]+Y], xy...)
    p2 = grad(S, t[t[X+X1]+Y+Y1], xy .- (X1, Y1) .+ SIMPLEX_UNSKEW_2D...)
    p3 = grad(S, t[t[X+1]+Y+1], xy .- 1 .+ 2SIMPLEX_UNSKEW_2D...)
    (p1 + p2 + p3) * 45.23065
end

# 3D

const SIMPLEX_SKEW_3D = 1 / 3
const SIMPLEX_UNSKEW_3D = 1 / 6

@doc doc_simplex_3d
simplex_3d(; seed=0) = _simplex(3, seed)

@inline function grad(S::Type{Simplex{3}}, hash, x, y, z)
    s = 0.6 - x^2 - y^2 - z^2
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
    t = sampler.state.table
    s = (x + y + z) * SIMPLEX_SKEW_3D
    X, Y, Z = floor.(Int, (x, y, z) .+ s)
    tx = (X + Y + Z) * SIMPLEX_UNSKEW_3D
    xyz = (x, y, z) .- (X, Y, Z) .+ tx
    X1, Y1, Z1, X2, Y2, Z2 = get_simplex(S, xyz...)
    p1 = grad(S, t[t[t[Z]+Y]+X], xyz...)
    p2 = grad(S, t[t[t[Z+Z1]+Y+Y1]+X+X1], xyz .- (X1, Y1, Z1) .+ SIMPLEX_UNSKEW_3D...)
    p3 = grad(S, t[t[t[Z+Z2]+Y+Y2]+X+X2], xyz .- (X2, Y2, Z2) .+ 2SIMPLEX_UNSKEW_3D...)
    p4 = grad(S, t[t[t[Z+1]+Y+1]+X+1], xyz .- 1 .+ 3SIMPLEX_UNSKEW_3D...)
    (p1 + p2 + p3 + p4) * 32
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

@doc doc_simplex_4d
simplex_4d(; seed=0) = _simplex(4, seed)

@inline function grad(S::Type{Simplex{4}}, hash, x, y, z, w)
    s = 0.6 - x^2 - y^2 - z^2 - w^2
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
    t = sampler.state.table
    s = (x + y + z + w) * SIMPLEX_SKEW_4D
    X, Y, Z, W = floor.(Int, (x, y, z, w) .+ s)
    tx = (X + Y + Z + W) * SIMPLEX_UNSKEW_4D
    v1 = (x, y, z, w) .- (X, Y, Z, W) .+ tx
    X1, Y1, Z1, W1, X2, Y2, Z2, W2, X3, Y3, Z3, W3 = get_simplex(S, v1...)
    v2 = v1 .- (X1, Y1, Z1, W1) .+ SIMPLEX_UNSKEW_4D
    v3 = v1 .- (X2, Y2, Z2, W2) .+ 2SIMPLEX_UNSKEW_4D
    v4 = v1 .- (X3, Y3, Z3, W3) .+ 3SIMPLEX_UNSKEW_4D
    v5 = v1 .- 1 .+ 4SIMPLEX_UNSKEW_4D
    p1 = grad(S, t[t[t[t[W]+Z]+Y]+X], v1...)
    p2 = grad(S, t[t[t[t[W+W1]+Z+Z1]+Y+Y1]+X+X1], v2...)
    p3 = grad(S, t[t[t[t[W+W2]+Z+Z2]+Y+Y2]+X+X2], v3...)
    p4 = grad(S, t[t[t[t[W+W3]+Z+Z3]+Y+Y3]+X+X3], v4...)
    p5 = grad(S, t[t[t[t[W+1]+Z+1]+Y+1]+X+1], v5...)
    (p1 + p2 + p3 + p4 + p5) * 27
end
