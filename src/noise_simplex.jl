struct Simplex{N} <: NoiseSampler{N}
    random_state::RandomState
    state::PerlinState
end

@inline function simplex(dims, seed)
    rs = RandomState(seed)
    Simplex{dims}(rs, PerlinState(rs))
end

HashTrait(::Type{<:Simplex}) = IsPerlinHashed()

include("noise_simplex_2d.jl")
include("noise_simplex_3d.jl")
include("noise_simplex_4d.jl")
