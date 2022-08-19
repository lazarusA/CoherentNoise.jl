module OpenSimplex2Noise

using CircularArrays: CircularVector
using FastPow
using ...Common: RandomState, Seed, get_seed
import ...Common: sample
using ..Noise: NoiseSampler

abstract type Orientation end
struct Standard <: Orientation end
struct ImproveX <: Orientation end
struct ImproveXY <: Orientation end
struct ImproveXZ <: Orientation end
struct ImproveXYZ <: Orientation end

struct OpenSimplex2{N,O<:Orientation} <: NoiseSampler{N}
    random_state::RandomState
end

const PRIME_X = 0x5205402b9270c86f
const PRIME_Y = 0x598cd327003817b5
const PRIME_Z = 0x5bcc226e9fa0bacb
const PRIME_W = 0x56cc5227e58f554b
const HASH_MULTIPLIER = 0x53a3f72deec546f5
const ROOT_2_OVER_2 = 0.7071067811865476
const ROOT_3_OVER_3 = 0.577350269189626

@inline function opensimplex2(dims, seed, orientation)
    rs = RandomState(seed)
    orientation = orientation_type(Val(orientation))
    OpenSimplex2{dims,orientation}(rs)
end

@inline orientation_type(::Val{nothing}) = Standard
@inline orientation_type(::Val{:x}) = ImproveX
@inline orientation_type(::Val{:xy}) = ImproveXY
@inline orientation_type(::Val{:xz}) = ImproveXZ
@inline orientation_type(::Val{:xyz}) = ImproveXYZ

include("noise_opensimplex2_2d.jl")
include("noise_opensimplex2_3d.jl")
include("noise_opensimplex2_4d.jl")

end
