function sample(sampler::S, x::T, y::T) where {S<:Value{2},T<:Real}
    seed = get_seed(sampler)
    primes = (PRIME_X, PRIME_Y)
    X, Y = floor.(Int, (x, y))
    X1, Y1 = (X, Y) .* primes
    X2, Y2 = (X1, Y1) .+ primes
    x1, y1 = curve3.((x, y) .- (X, Y))
    p1 = lerp(hash_coords(S, seed, X1, Y1), hash_coords(S, seed, X2, Y1), x1)
    p2 = lerp(hash_coords(S, seed, X1, Y2), hash_coords(S, seed, X2, Y2), x1)
    lerp(p1, p2, y1) - 1
end
