module OpenSimplexNoise

using Random: shuffle
using FastPow
using ...Common: RandomState, PerlinState, Seed
import ...Common: sample
using ..Noise: NoiseSampler

struct OpenSimplex{N} <: NoiseSampler{N}
    random_state::RandomState
    state::PerlinState
    table::Vector{UInt8}
end

@inline function opensimplex(dims, seed::Seed)
    rs = RandomState(seed)
    table = shuffle(rs.rng, UInt8.([repeat(0:3:45, 11)..., repeat(48:3:69, 10)...]))
    OpenSimplex{dims}(rs, PerlinState(rs), table)
end

include("noise_opensimplex_2d.jl")
include("noise_opensimplex_3d.jl")
include("noise_opensimplex_4d.jl")

end
