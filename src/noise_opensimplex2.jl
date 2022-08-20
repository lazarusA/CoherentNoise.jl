abstract type Orientation end
struct OrientStandard <: Orientation end
struct OrientX <: Orientation end
struct OrientXY <: Orientation end
struct OrientXZ <: Orientation end
struct OrientXYZ <: Orientation end

struct OpenSimplex2{N,O<:Orientation} <: NoiseSampler{N}
    random_state::RandomState
end

@inline function opensimplex2(dims, seed, orientation)
    rs = RandomState(seed)
    orientation = os2_orientation_type(Val(orientation))
    OpenSimplex2{dims,orientation}(rs)
end

@inline os2_orientation_type(::Val{nothing}) = OrientStandard
@inline os2_orientation_type(::Val{:x}) = OrientX
@inline os2_orientation_type(::Val{:xy}) = OrientXY
@inline os2_orientation_type(::Val{:xz}) = OrientXZ
@inline os2_orientation_type(::Val{:xyz}) = OrientXYZ

include("noise_opensimplex2_2d.jl")
include("noise_opensimplex2_3d.jl")
include("noise_opensimplex2_4d.jl")
