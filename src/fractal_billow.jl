struct BillowFractal{N,S,O} <: FractalSampler{N}
    state::FractalState{N,S,O}
end

@inline function BillowFractal{N}(
    seed,
    source::S,
    octaves,
    frequency,
    lacunarity,
    persistence,
) where {N,S<:AbstractSampler{N}}
    O = octaves
    fs = FractalState{N,BillowFractal,O}(seed, source, frequency, lacunarity, persistence, 1.0)
    BillowFractal{N,S,O}(fs)
end

"""
    billow_fractal_2d(; kwargs...)

Construct a sampler that outputs a 2-dimensional billow fractal noise when it is sampled from.

# Arguments

  - `seed=0`: An integer used to seed the random number generator for this sampler.

  - `source::AbstractSampler=opensimplex2s_2d()`: A 2-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""
function billow_fractal_2d(;
    seed=0,
    source=opensimplex2s_2d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    BillowFractal{2}(seed, source, octaves, frequency, lacunarity, persistence)
end

"""
    billow_fractal_3d(; kwargs...)

Construct a sampler that outputs a 3-dimensional billow fractal noise when it is sampled from.

# Arguments

  - `seed=0`: An integer used to seed the random number generator for this sampler.

  - `source::AbstractSampler=opensimplex2s_3d()`: A 3-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""
function billow_fractal_3d(;
    seed=0,
    source=opensimplex2s_3d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    BillowFractal{3}(seed, source, octaves, frequency, lacunarity, persistence)
end

"""
    billow_fractal_4d(; kwargs...)

Construct a sampler that outputs a 4-dimensional billow fractal noise when it is sampled from.

# Arguments

  - `seed=0`: An integer used to seed the random number generator for this sampler.

  - `source::AbstractSampler=opensimplex2s_4d()`: A 4-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""
function billow_fractal_4d(;
    seed=0,
    source=opensimplex2s_4d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    BillowFractal{4}(seed, source, octaves, frequency, lacunarity, persistence)
end

@inline function scale_factor(::Type{BillowFractal}, octaves, persistence, _)
    sum(persistence .^ (0:octaves-1))
end

function sample(sampler::BillowFractal{N}, coords::Vararg{Real,N}) where {N}
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
