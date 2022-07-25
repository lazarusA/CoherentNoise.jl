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
