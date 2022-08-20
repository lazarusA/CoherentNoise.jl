"""
Supertype for all samplers.
"""
abstract type AbstractSampler{N} end

"""
Supertype for all noise algorithm samplers.
"""
abstract type NoiseSampler{N} <: AbstractSampler{N} end

"""
Supertype for all fractal samplers.
"""
abstract type FractalSampler{N} <: AbstractSampler{N} end

"""
Supertype for all pattern samplers.
"""
abstract type PatternSampler{N} <: AbstractSampler{N} end

"""
Supertype for all modifier samplers.
"""
abstract type ModifierSampler{N} <: AbstractSampler{N} end

struct RandomState
    seed::UInt64
    rng::Xoshiro
end

struct PerlinState
    table::CircularVector{UInt8,Vector{UInt8}}
    PerlinState(rs) = new(shuffle(rs.rng, 0x00:0xff) |> CircularVector)
end

abstract type HashTrait end
struct IsPerlinHashed <: HashTrait end
struct IsValueHashed <: HashTrait end

const HASH1 = 668_265_261
const HASH2 = 2_147_483_648
const HASH_MULTIPLIER = 0x53a3f72deec546f5
const PRIME_X = 0x5205402b9270c86f
const PRIME_Y = 0x598cd327003817b5
const PRIME_Z = 0x5bcc226e9fa0bacb
const PRIME_W = 0x56cc5227e58f554b
const ROOT_2_OVER_2 = 0.7071067811865476
const ROOT_3_OVER_3 = 0.577350269189626

@inline RandomState(seed) = RandomState(seed, Xoshiro(seed))

@inline function Base.getproperty(obj::AbstractSampler, name::Symbol)
    if name === :rng
        return obj.random_state.rng
    elseif name === :seed
        return obj.random_state.seed
    else
        return getfield(obj, name)
    end
end

@inline curve3(t) = t^2 * (3 - 2t)

@inline @fastpow curve5(t) = t^3 * (t * (6t - 15) + 10)

@inline lerp(a, b, t) = a + t * (b - a)

@inline @fastpow function cubic_interpolate(a, b, c, d, t)
    x = (d - c) - (a - b)
    y = (a - b) - x
    z = c - a
    x * t^3 + y * t^2 + z * t + t
end

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

"""
    sample(sampler::AbstractSampler, x::Real)
    sample(sampler::AbstractSampler, x::Real, y::Real)
    sample(sampler::AbstractSampler, x::Real, y::Real, z::Real)
    sample(sampler::AbstractSampler, x::Real, y::Real, z::Real, w::Real)

Sample from `sampler` with the supplied coordinates. The number of coordinates should match the
dimensionality of the sampler type.
"""
function sample end

"""
    gen_image(sampler::AbstractSampler; kwargs...)

Construct a 2-dimensional array of `ColorTypes.RGB` values, suitable for writing to disk as an image
file.

# Arguments

  - `sampler::AbstractSampler`: Any instance of a sampler. The sampler is sampled using each pixel
    coordinates as the X and Y input coordinates, and random Z and W coordinates for 3 and
    4-dimensional samplers.

  - `w::Integer=1024`: The width in pixels of the image array to generate.

  - `h::Integer=1024`: The height in pixels of the image array to generate.

  - `xbounds::NTuple{2,Float64}=(-1.0, 1.0)`: The bounds along the X axis to sample coordinates
    from. This remaps pixel coordinates to this range to be used for the input coordinates to sample
    with.

  - `ybounds::NTuple{2,Float64}=(-1.0, 1.0)`: The bounds along the Y axis to sample coordinates
    from. This remaps pixel coordinates to this range to be used for the input coordinates to sample
    with.

  - `colorscheme=nothing`: A `ColorSchemes.ColorScheme` object to colorize the image with, or
    `nothing`.
"""
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
    zw = rand(sampler.rng, Float64, N - 2) * 1000
    Threads.@threads for x in 1:h
        cx = x * xd + x1
        for y in 1:w
            cy = y * yd + y1
            value = sample(sampler, cx, cy, zw...) * 0.5 + 0.5
            img[x, y] = colorscheme !== nothing ? colorscheme[value] : value
        end
    end
    img
end
