const SCALE_3D = 0.9649214285521897

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
    lerp(lerp(p1, p2, fy), lerp(p3, p4, fy), fz) * SCALE_3D
end
