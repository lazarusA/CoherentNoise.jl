struct MultiFractal{N,S,O} <: FractalSampler{N}
    state::State{N,S,O}
end

@inline function MultiFractal{N}(
    seed::Seed,
    source::S,
    octaves::Int,
    frequency,
    lacunarity,
    persistence,
) where {N,S<:AbstractSampler{N}}
    O = octaves
    fs = State{N,MultiFractal,O}(seed, source, frequency, lacunarity, persistence, 1.0)
    MultiFractal{N,S,O}(fs)
end

"""
    multi_fractal_2d(; kwargs...)

Construct a sampler that outputs a 2-dimensional multifractal noise when it is sampled from.

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

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""
function multi_fractal_2d(;
    seed=nothing,
    source=opensimplex2_2d(seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    MultiFractal{2}(seed, source, octaves, frequency, lacunarity, persistence)
end

"""
    multi_fractal_3d(; kwargs...)

Construct a sampler that outputs a 3-dimensional multifractal noise when it is sampled from.

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

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""
function multi_fractal_3d(;
    seed=nothing,
    source=opensimplex2_3d(seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    MultiFractal{3}(seed, source, octaves, frequency, lacunarity, persistence)
end

"""
    multi_fractal_4d(; kwargs...)

Construct a sampler that outputs a 4-dimensional multifractal noise when it is sampled from.

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

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""
function multi_fractal_4d(;
    seed=nothing,
    source=opensimplex2_4d(seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    MultiFractal{4}(seed, source, octaves, frequency, lacunarity, persistence)
end

@inline function scale_factor(::Type{MultiFractal}, octaves, persistence, _)
    reduce(1:octaves-1, init=1) do result, i
        result += result * persistence^i
    end
end

function sample(sampler::MultiFractal{N}, coords::Vararg{Real,N}) where {N}
    state = sampler.state
    sources = state.sources
    persistence = state.persistence
    lacunarity = state.lacunarity
    coords = coords .* state.frequency
    amplitude = 1.0
    result = sample(sources[1], coords...)
    for source in sources[2:end]
        amplitude *= persistence
        coords = coords .* lacunarity
        result += sample(source, coords...) * result * amplitude
    end
    result / state.scale
end
