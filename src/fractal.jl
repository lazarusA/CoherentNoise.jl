using ..Common: AbstractSampler, RandomState, Seed
import ..Common: sample, random_state, rng, get_seed
using ..Noise: NoiseSampler
using ..Noise.OpenSimplex2Noise: OpenSimplex2
using ..Noise.OpenSimplex2SNoise: OpenSimplex2S

"""
Supertype for all `N`-dimensional fractal samplers.
"""
abstract type FractalSampler{N} <: AbstractSampler{N} end

const Source = Union{Type{<:AbstractSampler},AbstractSampler}

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
    S::Source,
    frequency::Real,
    lacunarity::Real,
    persistence::Real,
    attenuation::Real;
    args...,
) where {N,F<:FractalSampler,O}
    if !(O in 1:32)
        throw(ArgumentError("fractal octaves must be in the range [1,32]"))
    end
    rs = RandomState(seed)
    sources, ST = build_sources(rs, S, N, O; args...)
    scale = scale_factor(F, O, persistence, attenuation)
    (State{N,ST,O}(rs, sources, scale, frequency, lacunarity, persistence, attenuation), ST)
end

function check_source_dims(fdims, sdims)
    if fdims !== sdims
        throw(ArgumentError("source dimensions do not match fractal"))
    end
end

function check_source_type(source_type)
    if !hasmethod(source_type, Tuple{}, (:seed,))
        throw(ArgumentError("invalid fractal source type: $(source_type)"))
    end
end

function build_sources(
    rs,
    source_type::Type{<:AbstractSampler{N}},
    dims,
    octaves;
    args...,
) where {N}
    check_source_dims(N, dims)
    check_source_type(source_type)
    T = typeof(source_type())
    (ntuple(i -> T(seed=rand(rs.rng, UInt64); args...), octaves), T)
end

function build_sources(rs, source_type::Type{<:AbstractSampler}, dims, octaves; args...)
    build_sources(rs, source_type{dims}, dims, octaves; args...)
end

function build_sources(rs, source::AbstractSampler{N}, dims, octaves; args...) where {N}
    check_source_dims(N, dims)
    (ntuple(i -> source, octaves), typeof(source))
end

@inline random_state(sampler::FractalSampler) = sampler.state.random_state
@inline rng(sampler::FractalSampler) = sampler.state.random_state.rng
@inline get_seed(sampler::FractalSampler) = sampler.state.random_state.seed
