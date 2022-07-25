function sample(sampler::S, x::T, y::T, z::T, w::T) where {S<:Value{4},T<:Real}
    seed = get_seed(sampler)
    primes = (PRIME_X, PRIME_Y, PRIME_Z, PRIME_W)
    X, Y, Z, W = floor.(Int, (x, y, z, w))
    X1, Y1, Z1, W1 = (X, Y, Z, W) .* primes
    X2, Y2, Z2, W2 = (X1, Y1, Z1, W1) .+ primes
    x1, y1, z1, w1 = curve3.((x, y, z, w) .- (X, Y, Z, W))
    p1 = lerp(hash_coords(S, seed, X1, Y1, Z1, W1), hash_coords(S, seed, X2, Y1, Z1, W1), x1)
    p2 = lerp(hash_coords(S, seed, X1, Y2, Z1, W1), hash_coords(S, seed, X2, Y2, Z1, W1), x1)
    p3 = lerp(hash_coords(S, seed, X1, Y1, Z2, W1), hash_coords(S, seed, X2, Y1, Z2, W1), x1)
    p4 = lerp(hash_coords(S, seed, X1, Y2, Z2, W1), hash_coords(S, seed, X2, Y2, Z2, W1), x1)
    p5 = lerp(hash_coords(S, seed, X1, Y1, Z1, W2), hash_coords(S, seed, X2, Y1, Z1, W2), x1)
    p6 = lerp(hash_coords(S, seed, X1, Y2, Z1, W2), hash_coords(S, seed, X2, Y2, Z1, W2), x1)
    p7 = lerp(hash_coords(S, seed, X1, Y1, Z2, W2), hash_coords(S, seed, X2, Y1, Z2, W2), x1)
    p8 = lerp(hash_coords(S, seed, X1, Y2, Z2, W2), hash_coords(S, seed, X2, Y2, Z2, W2), x1)
    p9 = lerp(lerp(p1, p2, y1), lerp(p3, p4, y1), z1)
    p10 = lerp(lerp(p5, p6, y1), lerp(p7, p8, y1), z1)
    lerp(p9, p10, w1) - 1
end
