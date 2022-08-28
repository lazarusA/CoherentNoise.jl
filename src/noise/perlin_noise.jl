struct Perlin{N} <: NoiseSampler{N}
    random_state::RandomState
    state::PerlinState
end

@inline function _perlin(dims, seed)
    rs = RandomState(seed)
    Perlin{dims}(rs, PerlinState(rs))
end

HashTrait(::Type{<:Perlin}) = IsPerlinHashed()

### 1D

@doc doc_perlin_1d
perlin_1d(; seed=nothing) = _perlin(1, seed)
@deprecate perlin_improved_1d(; kwargs...) perlin_1d(; kwargs...)

@inline function grad(S::Type{Perlin{1}}, hash, x)
    hash_coords(S, hash, x)
end

function sample(sampler::S, x::T) where {S<:Perlin{1},T<:Real}
    t = sampler.state.table
    X = floor(Int, x)
    x1 = x - X
    x2 = x1 - 1
    fx = curve5(x1)
    lerp(grad(S, t[t[X]], x1), grad(S, t[t[X+1]], x2), fx)
end

### 2D

@doc doc_perlin_2d
perlin_2d(; seed=nothing) = _perlin(2, seed)
@deprecate perlin_improved_2d(; kwargs...) perlin_2d(; kwargs...)

@inline function grad(S::Type{Perlin{2}}, hash, x, y)
    h = hash & 7
    u, v = h < 4 ? (x, y) : (y, x)
    hash_coords(S, h, u, v)
end

function sample(sampler::S, x::T, y::T) where {S<:Perlin{2},T<:Real}
    t = sampler.state.table
    X, Y = floor.(Int, (x, y))
    x1, y1 = (x, y) .- (X, Y)
    x2, y2 = (x1, y1) .- 1
    fx, fy = curve5.((x1, y1))
    a, b = (t[X] + Y, t[X+1] + Y)
    c1, c2, c3, c4 = (t[a], t[a+1], t[b], t[b+1])
    p1 = lerp(grad(S, t[c1], x1, y1), grad(S, t[c3], x2, y1), fx)
    p2 = lerp(grad(S, t[c2], x1, y2), grad(S, t[c4], x2, y2), fx)
    lerp(p1, p2, fy)
end

### 3D

@doc doc_perlin_3d
perlin_3d(; seed=nothing) = _perlin(3, seed)
@deprecate perlin_improved_3d(; kwargs...) perlin_3d(; kwargs...)

@inline function grad(S::Type{Perlin{3}}, hash, x, y, z)
    h = hash & 15
    u = h < 8 ? x : y
    v = h < 4 ? y : h == 12 || h == 14 ? x : z
    hash_coords(S, h, u, v)
end

function sample(sampler::S, x::T, y::T, z::T) where {S<:Perlin{3},T<:Real}
    t = sampler.state.table
    X, Y, Z = floor.(Int, (x, y, z))
    x1, y1, z1 = (x, y, z) .- (X, Y, Z)
    x2, y2, z2 = (x1, y1, z1) .- 1
    fx, fy, fz = curve5.((x1, y1, z1))
    a, b = (t[X] + Y, t[X+1] + Y)
    c1, c2, c3, c4 = (t[a] + Z, t[a+1] + Z, t[b] + Z, t[b+1] + Z)
    p1 = lerp(grad(S, t[c1], x1, y1, z1), grad(S, t[c3], x2, y1, z1), fx)
    p2 = lerp(grad(S, t[c2], x1, y2, z1), grad(S, t[c4], x2, y2, z1), fx)
    p3 = lerp(grad(S, t[c1+1], x1, y1, z2), grad(S, t[c3+1], x2, y1, z2), fx)
    p4 = lerp(grad(S, t[c2+1], x1, y2, z2), grad(S, t[c4+1], x2, y2, z2), fx)
    lerp(lerp(p1, p2, fy), lerp(p3, p4, fy), fz) * 0.9649214285521897
end

### 4D
#
@doc doc_perlin_4d
perlin_4d(; seed=nothing) = _perlin(4, seed)
@deprecate perlin_improved_4d(; kwargs...) perlin_4d(; kwargs...)

@inline function grad(S::Type{Perlin{4}}, hash, x, y, z, w)
    h1 = hash & 31
    h2 = h1 >> 3
    if h2 == 1
        hash_coords(S, h1, w, x, y)
    elseif h2 == 2
        hash_coords(S, h1, z, w, x)
    else
        hash_coords(S, h1, y, z, w)
    end
end

function sample(sampler::S, x::T, y::T, z::T, w::T) where {S<:Perlin{4},T<:Real}
    t = sampler.state.table
    X, Y, Z, W = floor.(Int, (x, y, z, w))
    x1, y1, z1, w1 = (x, y, z, w) .- (X, Y, Z, W)
    x2, y2, z2, w2 = (x1, y1, z1, w1) .- 1
    fx, fy, fz, fw = curve5.((x1, y1, z1, w1))
    a, b = (t[X] + Y, t[X+1] + Y)
    aa, ab, ba, bb = (t[a] + Z, t[a+1] + Z, t[b] + Z, t[b+1] + Z)
    c1, c2, c3, c4 = (t[aa] + W, t[ab] + W, t[aa+1] + W, t[ab+1] + W)
    c5, c6, c7, c8 = (t[ba] + W, t[bb] + W, t[ba+1] + W, t[bb+1] + W)
    p1 = lerp(grad(S, t[c1], x1, y1, z1, w1), grad(S, t[c5], x2, y1, z1, w1), fx)
    p2 = lerp(grad(S, t[c2], x1, y2, z1, w1), grad(S, t[c6], x2, y2, z1, w1), fx)
    p3 = lerp(grad(S, t[c3], x1, y1, z2, w1), grad(S, t[c7], x2, y1, z2, w1), fx)
    p4 = lerp(grad(S, t[c4], x1, y2, z2, w1), grad(S, t[c8], x2, y2, z2, w1), fx)
    p5 = lerp(grad(S, t[c1+1], x1, y1, z1, w2), grad(S, t[c5+1], x2, y1, z1, w2), fx)
    p6 = lerp(grad(S, t[c2+1], x1, y2, z1, w2), grad(S, t[c6+1], x2, y2, z1, w2), fx)
    p7 = lerp(grad(S, t[c3+1], x1, y1, z2, w2), grad(S, t[c7+1], x2, y1, z2, w2), fx)
    p8 = lerp(grad(S, t[c4+1], x1, y2, z2, w2), grad(S, t[c8+1], x2, y2, z2, w2), fx)
    p9 = lerp(lerp(p1, p2, fy), lerp(p3, p4, fy), fz)
    p10 = lerp(lerp(p5, p6, fy), lerp(p7, p8, fy), fz)
    lerp(p9, p10, fw) * 0.84
end
