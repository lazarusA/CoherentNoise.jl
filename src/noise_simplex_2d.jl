const SIMPLEX_SKEW_2D = (sqrt(3) - 1) / 2
const SIMPLEX_UNSKEW_2D = (3 - sqrt(3)) / 6

"""
    simplex_2d(; kwargs...)

Construct a sampler that outputs 2-dimensional Perlin Simplex noise when it is sampled from.

# Arguments

  - `seed=0`: An integer used to seed the random number generator for this sampler.
"""
simplex_2d(; seed=0) = simplex(2, seed)

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
