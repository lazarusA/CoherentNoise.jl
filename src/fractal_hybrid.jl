"""
An `N`-dimensional sampler that generates a hybrid multifractal.

`N` must be an integer in the range [2, 4].

`S` may be an instance of any [`AbstractSampler`](@ref), or an un-parameterized concrete type of
one.

`O` denotes the number of octaves to apply for this fractal.
"""
struct Hybrid{N,S,O} <: FractalSampler{N}
    state::State{N,S,O}
end

"""
    Hybrid{N}(; <kwargs>)

Construct a sampler that outputs an `N`-dimensional hybrid multifractal noise when it is sampler
from.

# Arguments

  - `seed::Union{Int,Nothing}=nothing`: An integer used to seed the random number generator for this
    sampler, or `nothing`. If a seed is not supplied, one will be generated automatically which will
    negatively affect reproducibility.
  - `source::Union{Type{<:AbstractSampler},AbstractSampler}=OpenSimplex2S`: A sampler type or
    sampler instance to use as the source of the fractal. If this is an instance, its dimensionality
    must match `N`.
  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.
  - `frequency=1.0`: The frequency of the first octave's signal.
  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.
  - `persistence=0.25`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
  - `args...`: Additional arguments to pass to the `source` if it is a type instead of an instance.
"""
function Hybrid{N}(;
    seed::Seed=nothing,
    source::Source=OpenSimplex2S,
    octaves::Int=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.25,
    args...,
) where {N}
    O = octaves
    fs, S = State{N,Hybrid,O}(seed, source, frequency, lacunarity, persistence, 1.0; args...)
    Hybrid{N,S,O}(fs)
end

function scale_factor(::Type{Hybrid}, octaves, persistence, _)
    amplitude = persistence
    weight = persistence^2
    result = persistence + weight
    for _ in 1:octaves-2
        amplitude *= persistence
        weight = max(weight, 1)
        weight *= amplitude
        result += weight
    end
    result
end

function sample(sampler::Hybrid{N}, coords::Vararg{Real,N}) where {N}
    state = sampler.state
    sources = state.sources
    persistence = state.persistence
    amplitude = persistence
    coords = coords .* state.frequency
    result = sample(sources[1], coords...) * amplitude
    coords = coords .* state.lacunarity
    weight = sample(sources[2], coords...) * amplitude * result
    result += weight
    for source in sources[3:end]
        amplitude *= persistence
        weight = max(weight, 1)
        weight *= sample(source, coords...) * amplitude
        result += weight
    end
    result / state.scale
end
