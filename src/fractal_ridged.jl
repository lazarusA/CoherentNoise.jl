struct RidgedFractal{N,S,O} <: FractalSampler{N}
    state::State{N,S,O}
end

@inline function RidgedFractal{N}(
    seed::Seed,
    source::S,
    octaves::Int,
    frequency,
    lacunarity,
    persistence,
    attenuation,
) where {N,S<:AbstractSampler{N}}
    O = octaves
    fs = State{N,RidgedFractal,O}(seed, source, frequency, lacunarity, persistence, attenuation)
    RidgedFractal{N,S,O}(fs)
end

"""
    ridged_fractal_2d(; kwargs...)

Construct a sampler that outputs a 2-dimensional ridged multifractal noise when it is sampled from.

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

  - `persistence=1.0`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.

  - `attenuation=2.0`: The attenuation to apply to the weight of each octave.
"""
function ridged_fractal_2d(;
    seed=nothing,
    source=opensimplex2s_2d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=1.0,
    attenuation=2.0,
)
    RidgedFractal{2}(seed, source, octaves, frequency, lacunarity, persistence, attenuation)
end

"""
    ridged_fractal_3d(; kwargs...)

Construct a sampler that outputs a 3-dimensional ridged multifractal noise when it is sampled from.

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

  - `persistence=1.0`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.

  - `attenuation=2.0`: The attenuation to apply to the weight of each octave.
"""
function ridged_fractal_3d(;
    seed=nothing,
    source=opensimplex2s_3d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=1.0,
    attenuation=2.0,
)
    RidgedFractal{3}(seed, source, octaves, frequency, lacunarity, persistence, attenuation)
end

"""
    ridged_fractal_4d(; kwargs...)

Construct a sampler that outputs a 4-dimensional ridged multifractal noise when it is sampled from.

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

  - `persistence=1.0`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.

  - `attenuation=2.0`: The attenuation to apply to the weight of each octave.
"""
function ridged_fractal_4d(;
    seed=nothing,
    source=opensimplex2s_4d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=1.0,
    attenuation=2.0,
)
    RidgedFractal{4}(seed, source, octaves, frequency, lacunarity, persistence, attenuation)
end

@inline function scale_factor(::Type{RidgedFractal}, octaves, persistence, attenuation)
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

function sample(sampler::RidgedFractal{N}, coords::Vararg{Real,N}) where {N}
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
