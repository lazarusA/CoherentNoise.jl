module OpenSimplex2SNoise

using CircularArrays: CircularVector
using FastPow
using ...Common: RandomState, Seed, get_seed
import ...Common: sample
using ..Noise: NoiseSampler
using ..OpenSimplex2Noise: PRIME_X, PRIME_Y, PRIME_Z, PRIME_W, ROOT_2_OVER_2, ROOT_3_OVER_3, grad
using ..OpenSimplex2Noise: Orientation, Standard, ImproveX, ImproveXY, ImproveXZ, ImproveXYZ

"""
An `N`-dimensional sampler that generates OpenSimplex2S noise, with an orientation of `O`.

OpenSimplex2S is smoother than [`OpenSimplex2`](@ref Main.CoherentNoise.Noise.OpenSimplex2Noise.OpenSimplex2), at the expense of being less performant.

`N` must be an integer in the range [2, 4].

`O` must be an [`Orientation`](@ref) type.
"""
struct OpenSimplex2S{N,O<:Orientation} <: NoiseSampler{N}
    random_state::RandomState
end

"""
    OpenSimplex2S{N,O}(; seed=nothing)

Construct a sampler that outputs `N`-dimensional OpenSimplex2S noise when it is sampler from.

The noise space is re-oriented according by the supplied [`Orientation`](@ref) type `O`.

# Arguments

  - `seed::Union{Int,Nothing}=nothing`: An integer used to seed the random number generator for this
    sampler, or `nothing`. If a seed is not supplied, one will be generated automatically which will
    negatively affect reproducibility.
"""
function OpenSimplex2S{N,O}(; seed::Seed=nothing) where {N,O}
    rs = RandomState(seed)
    OpenSimplex2S{N,O}(rs)
end

"""
    OpenSimplex2S{N}(; seed=nothing)

Construct a sampler that outputs `N`-dimensional OpenSimplex2S noise when it is sampler from.

The noise space will be re-oriented using [`Standard`](@ref) orientation.

# Arguments

  - `seed::Union{Int,Nothing}=nothing`: An integer used to seed the random number generator for this
    sampler, or `nothing`. If a seed is not supplied, one will be generated automatically which will
    negatively affect reproducibility.
"""
OpenSimplex2S{N}(; seed::Seed=nothing) where {N} = OpenSimplex2S{N,Standard}(seed=seed)

include("noise_opensimplex2s_2d.jl")
include("noise_opensimplex2s_3d.jl")
include("noise_opensimplex2s_4d.jl")

end
