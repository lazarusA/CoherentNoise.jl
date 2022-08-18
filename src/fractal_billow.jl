"""
An `N`-dimensional sampler that generates a billowy fractal.

`N` must be an integer in the range [2, 4].

`S` may be an instance of any [`AbstractSampler`](@ref), or an un-parameterized concrete type of
one.

`O` denotes the number of octaves to apply for this fractal.
"""
struct Billow{N,S,O} <: FractalSampler{N}
    state::State{N,S,O}
end

"""
    Billow{N}(; kwargs...)

Construct a sampler that outputs an `N`-dimensional billowy fractal noise when it is sampler from.

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
  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
  - `args...`: Additional arguments to pass to the `source` if it is a type instead of an instance.
"""
function Billow{N}(;
    seed::Seed=nothing,
    source::Source=OpenSimplex2S,
    octaves::Int=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
    args...,
) where {N}
    O = octaves
    fs, S = State{N,Billow,O}(seed, source, frequency, lacunarity, persistence, 1.0; args...)
    Billow{N,S,O}(fs)
end

function scale_factor(::Type{Billow}, octaves, persistence, _)
    sum(persistence .^ (0:octaves-1))
end

function sample(sampler::Billow{N}, coords::Vararg{Real,N}) where {N}
    state = sampler.state
    persistence = state.persistence
    lacunarity = state.lacunarity
    coords = coords .* state.frequency
    amplitude = 1.0
    result = 0.0
    for source in state.sources
        result += (abs(sample(source, coords...)) * 2 - 1) * amplitude
        amplitude *= persistence
        coords = coords .* lacunarity
    end
    result / state.scale
end
