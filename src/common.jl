### Abstract sampler types

"Supertype for all samplers."
abstract type AbstractSampler{N} end

"Supertype for all pattern samplers."
abstract type PatternSampler{N} <: AbstractSampler{N} end

"Supertype for all noise algorithm samplers."
abstract type NoiseSampler{N} <: AbstractSampler{N} end

"Supertype for all fractal samplers."
abstract type FractalSampler{N} <: AbstractSampler{N} end

"Supertype for all modifier samplers."
abstract type ModifierSampler{N} <: AbstractSampler{N} end

### Constants shared by multiple samplers

const HASH1 = 668_265_261
const HASH2 = 2_147_483_648
const HASH_MULTIPLIER = 0x53a3f72deec546f5
const PRIME_X = 0x5205402b9270c86f
const PRIME_Y = 0x598cd327003817b5
const PRIME_Z = 0x5bcc226e9fa0bacb
const PRIME_W = 0x56cc5227e58f554b
const ROOT_2_OVER_2 = 0.7071067811865476
const ROOT_3_OVER_3 = 0.577350269189626

### Random number generation

struct RandomState
    seed::UInt64
    rng::Xoroshiro128Star
end

@inline RandomState(seed) = RandomState(seed, Xoroshiro128Star(seed))

### State for Perlin-based samplers

struct PerlinState
    table::CircularVector{UInt8,Vector{UInt8}}
    PerlinState(rs) = new(shuffle(rs.rng, 0x00:0xff) |> CircularVector)
end

### State for Simplex-based samplers

struct SimplexState
    falloff::Float64
    scale_factor::Float64
end

### Trait to decide how some samplers should generate hash values

abstract type HashTrait end
struct IsPerlinHashed <: HashTrait end
struct IsValueHashed <: HashTrait end

@inline hash_coords(sampler::S, args...) where {S} = hash_coords(HashTrait(sampler), args...)

@inline function hash_coords(::IsPerlinHashed, hash, u)
    iszero(hash & 1) ? -u : u
end

@inline function hash_coords(::IsPerlinHashed, hash, u, v)
    (iszero(hash & 1) ? -u : u) + (iszero(hash & 2) ? -v : v)
end

@inline function hash_coords(::IsPerlinHashed, hash, u, v, w)
    (iszero(hash & 1) ? -u : u) + (iszero(hash & 2) ? -v : v) + (iszero(hash & 4) ? -w : w)
end

@inline function hash_coords(::IsValueHashed, seed, coords...)
    hash = HASH1 * ⊻(seed, coords...)^2
    (hash ⊻ hash >> 19) % UInt32 / HASH2
end

### Utility functions

@inline pow3(x) = x * x * x

@inline pow4(x) = x * x * x * x

@inline curve3(t) = t^2 * (3 - 2t)

@inline curve5(t) = pow3(t) * (t * (6t - 15) + 10)

@inline lerp(a, b, t) = a + t * (b - a)

@inline function cubic_interpolate(a, b, c, d, t)
    x = (d - c) - (a - b)
    y = (a - b) - x
    z = c - a
    x * pow3(t) + y * t^2 + z * t + t
end

### Public interface functions

@doc doc_sampler
function sample end

@doc doc_gen_image
function gen_image(
    sampler::S;
    w::Integer=1024,
    h::Integer=1024,
    xbounds::NTuple{2,Float64}=(-1.0, 1.0),
    ybounds::NTuple{2,Float64}=(-1.0, 1.0),
    colorscheme::Union{ColorScheme,Nothing}=nothing,
) where {N,S<:AbstractSampler{N}}
    x1, x2 = xbounds
    y1, y2 = ybounds
    xd = (x2 - x1) / w
    yd = (y2 - y1) / h
    img = Array{RGB{Float64}}(undef, h, w)
    zw = rand(sampler.random_state.rng, Float64, N - 2) * 1000
    Threads.@threads for x in 1:h
        cx = x * xd + x1
        for y in 1:w
            cy = y * yd + y1
            value = clamp(sample(sampler, cx, cy, zw...) * 0.5 + 0.5, 0, 1)
            img[x, y] = isnothing(colorscheme) ? value : colorscheme[value]
        end
    end
    img
end
