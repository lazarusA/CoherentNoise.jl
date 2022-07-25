abstract type HashTrait end
struct IsPerlinHashed <: HashTrait end
struct IsValueHashed <: HashTrait end

@inline hash_coords(sampler::S, args...) where {S} = hash_coords(HashTrait(sampler), args...)

@inline function hash_coords(::IsPerlinHashed, hash, u, v)
    (iszero(hash & 1) ? u : -u) + (iszero(hash & 2) ? v : -v)
end

@inline function hash_coords(::IsPerlinHashed, hash, u, v, w)
    (iszero(hash & 1) ? u : -u) + (iszero(hash & 2) ? v : -v) + (iszero(hash & 4) ? w : -w)
end

@inline function hash_coords(::IsValueHashed, seed, coords...)
    hash = HASH1 * ⊻(seed, coords...)^2
    (hash ⊻ hash >> 19) % UInt32 / HASH2
end
