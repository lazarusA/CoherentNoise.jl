module OpenSimplex2SNoise

using CircularArrays: CircularVector
using FastPow
using ...Common: RandomState, Seed, get_seed
import ...Common: sample
using ..Noise: NoiseSampler
using ..OpenSimplex2Noise: PRIME_X, PRIME_Y, PRIME_Z, PRIME_W, ROOT_2_OVER_2, ROOT_3_OVER_3, grad
using ..OpenSimplex2Noise: Orientation, Standard, ImproveX, ImproveXY, ImproveXZ, ImproveXYZ
using ..OpenSimplex2Noise: orientation_type

struct OpenSimplex2S{N,O<:Orientation} <: NoiseSampler{N}
    random_state::RandomState
end

@inline function opensimplex2s(dims, seed, orientation)
    rs = RandomState(seed)
    orientation = orientation_type(Val(orientation))
    OpenSimplex2S{dims,orientation}(rs)
end

include("noise_opensimplex2s_2d.jl")
include("noise_opensimplex2s_3d.jl")
include("noise_opensimplex2s_4d.jl")

end
