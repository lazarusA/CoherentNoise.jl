struct OpenSimplex2S{N,O<:Orientation} <: NoiseSampler{N}
    random_state::RandomState
end

@inline function opensimplex2s(dims, seed, orientation)
    rs = RandomState(seed)
    orientation = os2_orientation_type(Val(orientation))
    OpenSimplex2S{dims,orientation}(rs)
end

include("noise_opensimplex2s_2d.jl")
include("noise_opensimplex2s_3d.jl")
include("noise_opensimplex2s_4d.jl")
