"""
    value_3d(; kwargs...)

Construct a sampler that outputs 3-dimensonal value noise when it is sampled from.

# Arguments

  - `seed=nothing`: An integer used to seed the random number generator for this sampler, or
    `nothing`. If a seed is not supplied, one will be generated automatically which will negatively
    affect reproducibility.
"""
value_3d(; seed=nothing) = value(3, seed)

function sample(sampler::S, x::T, y::T, z::T) where {S<:Value{3},T<:Real}
    seed = get_seed(sampler)
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
