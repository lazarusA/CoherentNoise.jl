struct FBMFractal{N,S,O} <: FractalSampler{N}
    state::State{N,S,O}
end

@inline function FBMFractal{N}(
    seed::Seed,
    source::S,
    octaves::Int,
    frequency,
    lacunarity,
    persistence,
) where {N,S<:AbstractSampler{N}}
    O = octaves
    fs = State{N,FBMFractal,O}(seed, source, frequency, lacunarity, persistence, 1.0)
    FBMFractal{N,S,O}(fs)
end

"""
    fbm_fractal_2d(; kwargs...)

Construct a sampler that outputs a 2-dimensional fractional Brownian motion fractal noise when it
is sampled from.

# Arguments

  - `seed=nothing`: An integer used to seed the random number generator for this sampler, or
    `nothing`. If a seed is not supplied, one will be generated automatically which will negatively
    affect reproducibility.

  - `source::AbstractSampler=opensimplex2_2d()`: A 2-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""
function fbm_fractal_2d(;
    seed=nothing,
    source=opensimplex2_2d(seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    FBMFractal{2}(seed, source, octaves, frequency, lacunarity, persistence)
end

"""
    fbm_fractal_3d(; kwargs...)

Construct a sampler that outputs a 3-dimensional fractional Brownian motion fractal noise when it
is sampled from.

# Arguments

  - `seed=nothing`: An integer used to seed the random number generator for this sampler, or
    `nothing`. If a seed is not supplied, one will be generated automatically which will negatively
    affect reproducibility.

  - `source::AbstractSampler=opensimplex2_3d()`: A 3-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""
function fbm_fractal_3d(;
    seed=nothing,
    source=opensimplex2_3d(seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    FBMFractal{3}(seed, source, octaves, frequency, lacunarity, persistence)
end

"""
    fbm_fractal_4d(; kwargs...)

Construct a sampler that outputs a 4-dimensional fractional Brownian motion fractal noise when it
is sampled from.

# Arguments

  - `seed=nothing`: An integer used to seed the random number generator for this sampler, or
    `nothing`. If a seed is not supplied, one will be generated automatically which will negatively
    affect reproducibility.

  - `source::AbstractSampler=opensimplex2_4d()`: A 4-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""
function fbm_fractal_4d(;
    seed=nothing,
    source=opensimplex2_4d(seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    FBMFractal{4}(seed, source, octaves, frequency, lacunarity, persistence)
end

@inline function scale_factor(::Type{FBMFractal}, octaves, persistence, _)
    sum(persistence .^ (0:octaves-1))
end

function sample(sampler::FBMFractal{N}, coords::Vararg{Real,N}) where {N}
    state = sampler.state
    persistence = state.persistence
    lacunarity = state.lacunarity
    coords = coords .* state.frequency
    amplitude = 1.0
    result = 0.0
    for source in state.sources
        result += sample(source, coords...) * amplitude
        amplitude *= persistence
        coords = coords .* lacunarity
    end
    result / state.scale
end
