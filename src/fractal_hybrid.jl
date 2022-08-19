struct HybridFractal{N,S,O} <: FractalSampler{N}
    state::State{N,S,O}
end

@inline function HybridFractal{N}(
    seed::Seed,
    source::S,
    octaves::Int,
    frequency,
    lacunarity,
    persistence,
) where {N,S<:AbstractSampler{N}}
    O = octaves
    fs = State{N,HybridFractal,O}(seed, source, frequency, lacunarity, persistence, 1.0)
    HybridFractal{N,S,O}(fs)
end

"""
    hybrid_fractal_2d(; kwargs...)

Construct a sampler that outputs a 2-dimensional hybrid multifractal noise when it is sampled from.

# Arguments

  - `seed=nothing`: An integer used to seed the random number generator for this sampler, or
    `nothing`. If a seed is not supplied, one will be generated automatically which will negatively
    affect reproducibility.

  - `source::AbstractSampler=opensimplex2s_2d()`: A 2-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.25`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""
function hybrid_fractal_2d(;
    seed=nothing,
    source=opensimplex2s_2d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.25,
)
    HybridFractal{2}(seed, source, octaves, frequency, lacunarity, persistence)
end

"""
    hybrid_fractal_3d(; kwargs...)

Construct a sampler that outputs a 3-dimensional hybrid multifractal noise when it is sampled from.

# Arguments

  - `seed=nothing`: An integer used to seed the random number generator for this sampler, or
    `nothing`. If a seed is not supplied, one will be generated automatically which will negatively
    affect reproducibility.

  - `source::AbstractSampler=opensimplex2s_3d()`: A 3-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.25`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""
function hybrid_fractal_3d(;
    seed=nothing,
    source=opensimplex2s_3d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.25,
)
    HybridFractal{3}(seed, source, octaves, frequency, lacunarity, persistence)
end

"""
    hybrid_fractal_4d(; kwargs...)

Construct a sampler that outputs a 4-dimensional hybrid multifractal noise when it is sampled from.

# Arguments

  - `seed=nothing`: An integer used to seed the random number generator for this sampler, or
    `nothing`. If a seed is not supplied, one will be generated automatically which will negatively
    affect reproducibility.

  - `source::AbstractSampler=opensimplex2s_4d()`: A 4-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.25`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""
function hybrid_fractal_4d(;
    seed=nothing,
    source=opensimplex2s_4d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.25,
)
    HybridFractal{4}(seed, source, octaves, frequency, lacunarity, persistence)
end

function scale_factor(::Type{HybridFractal}, octaves, persistence, _)
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

function sample(sampler::HybridFractal{N}, coords::Vararg{Real,N}) where {N}
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
