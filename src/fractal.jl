using Random: seed!
using ..Common: AbstractSampler, RandomState, Seed
import ..Common: sample, random_state, rng, get_seed
using ..Noise: NoiseSampler
using ..Noise.OpenSimplex2Noise: opensimplex2_2d, opensimplex2_3d, opensimplex2_4d
using ..Noise.OpenSimplex2SNoise: opensimplex2s_2d, opensimplex2s_3d, opensimplex2s_4d

"""
Supertype for all fractal samplers.
"""
abstract type FractalSampler{N} <: AbstractSampler{N} end

struct State{N,S,O}
    random_state::RandomState
    sources::NTuple{O,S}
    scale::Float64
    frequency::Float64
    lacunarity::Float64
    persistence::Float64
    attenuation::Float64
end

@inline function State{N,F,O}(
    seed::Seed,
    source::S,
    frequency::Real,
    lacunarity::Real,
    persistence::Real,
    attenuation::Real,
) where {N,F<:FractalSampler,O,S<:AbstractSampler{N}}
    rs = RandomState(seed)
    sources = ntuple(_ -> deepcopy(source), O)
    scale = scale_factor(F, O, persistence, attenuation)
    State{N,S,O}(rs, sources, scale, frequency, lacunarity, persistence, attenuation)
end

@inline random_state(sampler::FractalSampler) = sampler.state.random_state
@inline rng(sampler::FractalSampler) = sampler.state.random_state.rng
@inline get_seed(sampler::FractalSampler) = sampler.state.random_state.seed
