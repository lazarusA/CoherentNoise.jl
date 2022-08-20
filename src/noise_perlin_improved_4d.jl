"""
    perlin_improved_4d(; kwargs...)

Construct a sampler that outputs 4-dimensional Perlin "Improved" noise when it is sampled from.

# Arguments

  - `seed=0`: An integer used to seed the random number generator for this sampler.
"""
perlin_improved_4d(; seed=0) = perlin_improved(4, seed)

@inline function grad(S::Type{PerlinImproved{4}}, hash, x, y, z, w)
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

function sample(sampler::S, x::T, y::T, z::T, w::T) where {S<:PerlinImproved{4},T<:Real}
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
