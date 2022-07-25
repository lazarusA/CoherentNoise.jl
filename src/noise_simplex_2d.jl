const SKEW_2D = (sqrt(3) - 1) / 2

const UNSKEW_2D = (3 - sqrt(3)) / 6

const SCALE_2D = 45.23065

@inline function grad(::Type{Simplex{2}}, hash, x, y)
    s = 0.5 - x^2 - y^2
    h = hash & 7
    u, v = h < 4 ? (x, y) : (y, x)
    g = (iszero(h & 1) ? u : -u) + (iszero(h & 2) ? 2v : -2v)
    @fastpow s > 0 ? s^4 * g : 0.0
end

@inline get_simplex(::Type{Simplex{2}}, x, y) = x > y ? (1, 0) : (0, 1)

function sample(sampler::S, x::T, y::T) where {S<:Simplex{2},T<:Real}
    t = sampler.state.table
    s = (x + y) * SKEW_2D
    X, Y = floor.(Int, (x, y) .+ s)
    tx = (X + Y) * UNSKEW_2D
    xy = (x, y) .- (X, Y) .+ tx
    X1, Y1 = get_simplex(S, xy...)
    p1 = grad(S, t[t[X]+Y], xy...)
    p2 = grad(S, t[t[X+X1]+Y+Y1], xy .- (X1, Y1) .+ UNSKEW_2D...)
    p3 = grad(S, t[t[X+1]+Y+1], xy .- 1 .+ 2UNSKEW_2D...)
    (p1 + p2 + p3) * SCALE_2D
end
