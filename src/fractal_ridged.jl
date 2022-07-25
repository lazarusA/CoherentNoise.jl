"""
An `N`-dimensional sampler that generates a ridged multifractal.

`N` must be an integer in the range [2, 4].

`S` may be an instance of any [`AbstractSampler`](@ref), or an un-parameterized concrete type of
one.

`O` denotes the number of octaves to apply for this fractal.
"""
struct Ridged{N,S,O} <: FractalSampler{N}
    state::State{N,S,O}
end

"""
    Ridged{N}(; <kwargs>)

Construct a sampler that outputs an `N`-dimensional ridged multifractal noise when it is sampler
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
  - `persistence=1.0`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
  - `attenuation=2.0`: The attenuation to apply to the weight of each octave.
  - `args...`: Additional arguments to pass to the `source` if it is a type instead of an instance.
"""
function Ridged{N}(;
    seed::Seed=nothing,
    source::Source=OpenSimplex2S,
    octaves::Int=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=1.0,
    attenuation=2.0,
    args...,
) where {N}
    O = octaves
    fs, S = State{N,Ridged,O}(
        seed,
        source,
        frequency,
        lacunarity,
        persistence,
        attenuation;
        args...)
    Ridged{N,S,O}(fs)
end

function scale_factor(::Type{Ridged}, octaves, persistence, attenuation)
    amplitude = 1.0
    weight = 1.0
    result = 0.0
    for _ in 1:octaves
        sample = weight * amplitude
        amplitude *= persistence
        weight = clamp(sample / attenuation, 0, 1)
        result += sample
    end
    2 / result
end

function sample(sampler::Ridged{N}, coords::Vararg{Real,N}) where {N}
    state = sampler.state
    persistence = state.persistence
    lacunarity = state.lacunarity
    attenuation = state.attenuation
    coords = coords .* state.frequency
    amplitude = 1.0
    weight = 1.0
    result = 0.0
    for source in state.sources
        temp = (1 - abs(sample(source, coords...)))^2 * weight * amplitude
        amplitude *= persistence
        coords = coords .* lacunarity
        weight = clamp(temp / attenuation, 0, 1)
        result += temp
    end
    result * state.scale - 1
end
