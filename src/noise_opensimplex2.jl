module OpenSimplex2Noise

using CircularArrays: CircularVector
using FastPow
using ...Common: RandomState, Seed, get_seed
import ...Common: sample
using ..Noise: NoiseSampler

"""
Supertype of all orientations that can be applied to the [`OpenSimplex2`](@ref) or
[`OpenSimplex2S`](@ref Main.CoherentNoise.Noise.OpenSimplex2SNoise.OpenSimplex2S) sampler types to re-orient the noise space.
"""
abstract type Orientation end

"""
The standard orientation for [`OpenSimplex2`](@ref) and [`OpenSimplex2S`](@ref Main.CoherentNoise.Noise.OpenSimplex2SNoise.OpenSimplex2S).
"""
struct Standard <: Orientation end
"""
Re-orient the noise space with the Y axis pointing down the main diagonal.
"""
struct ImproveX <: Orientation end
"""
Re-orient the noise space to have better visual isotropy in the XY plane.
"""
struct ImproveXY <: Orientation end
"""
Re-orient the noise space to have better visual isotropy in the XZ plane.
"""
struct ImproveXZ <: Orientation end
"""
Re-orient a 4-dimensional noise space to be better suited for time-varied animations, where the W
axis is time.
"""
struct ImproveXYZ <: Orientation end

"""
An `N`-dimensional sampler that generates OpenSimplex2 noise, with an orientation of `O`.

OpenSimplex2 is not as smooth as [`OpenSimplex2S`](@ref Main.CoherentNoise.Noise.OpenSimplex2SNoise.OpenSimplex2S), although more performant.

`N` must be an integer in the range [2, 4].

`O` must be an [`Orientation`](@ref) type.
"""
struct OpenSimplex2{N,O<:Orientation} <: NoiseSampler{N}
    random_state::RandomState
end

"""
    OpenSimplex2{N,O}(; seed=nothing)

Construct a sampler that outputs `N`-dimensional OpenSimplex2 noise when it is sampler from.

The noise space is re-oriented according by the supplied [`Orientation`](@ref) type `O`.

# Arguments

  - `seed::Union{Int,Nothing}=nothing`: An integer used to seed the random number generator for this
    sampler, or `nothing`. If a seed is not supplied, one will be generated automatically which will
    negatively affect reproducibility.
"""
function OpenSimplex2{N,O}(; seed::Seed=nothing) where {N,O}
    rs = RandomState(seed)
    OpenSimplex2{N,O}(rs)
end

"""
    OpenSimplex2{N}(; seed=nothing)

Construct a sampler that outputs `N`-dimensional OpenSimplex2 noise when it is sampler from.

The noise space will be re-oriented using [`Standard`](@ref) orientation.

# Arguments

  - `seed::Union{Int,Nothing}=nothing`: An integer used to seed the random number generator for this
    sampler, or `nothing`. If a seed is not supplied, one will be generated automatically which will
    negatively affect reproducibility.
"""
OpenSimplex2{N}(; seed::Seed=nothing) where {N} = OpenSimplex2{N,Standard}(seed=seed)

const PRIME_X = 0x5205402b9270c86f

const PRIME_Y = 0x598cd327003817b5

const PRIME_Z = 0x5bcc226e9fa0bacb

const PRIME_W = 0x56cc5227e58f554b

const HASH_MULTIPLIER = 0x53a3f72deec546f5

const ROOT_2_OVER_2 = 0.7071067811865476

const ROOT_3_OVER_3 = 0.577350269189626

include("noise_opensimplex2_2d.jl")
include("noise_opensimplex2_3d.jl")
include("noise_opensimplex2_4d.jl")

end
