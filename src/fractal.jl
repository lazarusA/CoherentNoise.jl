### State composed into all fractal sampler types

struct FractalState{N,S,O}
    random_state::RandomState
    offset::NTuple{N,Float64}
    sources::NTuple{O,S}
    scale::Float64
    frequency::Float64
    lacunarity::Float64
    persistence::Float64
    attenuation::Float64
end

@inline function FractalState{N,F,O}(
    seed,
    source::S,
    frequency,
    lacunarity,
    persistence,
    attenuation,
) where {N,F<:FractalSampler,O,S<:AbstractSampler{N}}
    rs = RandomState(seed)
    offset = ntuple(i -> rand(rs.rng), N)
    sources = ntuple(_ -> deepcopy(source), O)
    scale = scale_factor(F, O, persistence, attenuation)
    FractalState{N,S,O}(rs, offset, sources, scale, frequency, lacunarity, persistence, attenuation)
end

@inline function Base.getproperty(obj::FractalSampler, name::Symbol)
    name === :random_state ? getfield(obj, :state).random_state : getfield(obj, name)
end

### fBm

struct FBM{N,S,O} <: FractalSampler{N}
    state::FractalState{N,S,O}
end

@inline function FBM{N}(
    seed,
    source::S,
    octaves,
    frequency,
    lacunarity,
    persistence,
) where {N,S<:AbstractSampler{N}}
    O = octaves
    fs = FractalState{N,FBM,O}(seed, source, frequency, lacunarity, persistence, 1.0)
    FBM{N,S,O}(fs)
end

@inline function scale_factor(::Type{FBM}, octaves, persistence, _)
    sum(persistence .^ (0:octaves-1))
end

"""
    fbm_fractal_1d(; kwargs...)

Construct a sampler that outputs a 1-dimensional fractional Brownian motion fractal noise when it
is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

  - `source::AbstractSampler=simplex_1d()`: A 1-dimensional sampler instance to use as the source of
    the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""
function fbm_fractal_1d(;
    seed=nothing,
    source=simplex_1d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    FBM{1}(seed, source, octaves, frequency, lacunarity, persistence)
end

"""
    fbm_fractal_2d(; kwargs...)

Construct a sampler that outputs a 2-dimensional fractional Brownian motion fractal noise when it
is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

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
    source=opensimplex2_2d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    FBM{2}(seed, source, octaves, frequency, lacunarity, persistence)
end

"""
    fbm_fractal_3d(; kwargs...)

Construct a sampler that outputs a 3-dimensional fractional Brownian motion fractal noise when it
is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

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
    source=opensimplex2_3d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    FBM{3}(seed, source, octaves, frequency, lacunarity, persistence)
end

"""
    fbm_fractal_4d(; kwargs...)

Construct a sampler that outputs a 4-dimensional fractional Brownian motion fractal noise when it
is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

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
    source=opensimplex2_4d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    FBM{4}(seed, source, octaves, frequency, lacunarity, persistence)
end

function sample(sampler::FBM{N}, coords::Vararg{Real,N}) where {N}
    state = sampler.state
    offset = state.offset
    persistence = state.persistence
    lacunarity = state.lacunarity
    coords = coords .* state.frequency
    amplitude = 1.0
    result = 0.0
    for source in state.sources
        result += sample(source, coords...) * amplitude
        amplitude *= persistence
        coords = coords .* lacunarity .+ offset
    end
    result / state.scale
end

### Billow

struct Billow{N,S,O} <: FractalSampler{N}
    state::FractalState{N,S,O}
end

@inline function Billow{N}(
    seed,
    source::S,
    octaves,
    frequency,
    lacunarity,
    persistence,
) where {N,S<:AbstractSampler{N}}
    O = octaves
    fs = FractalState{N,Billow,O}(seed, source, frequency, lacunarity, persistence, 1.0)
    Billow{N,S,O}(fs)
end

@inline function scale_factor(::Type{Billow}, octaves, persistence, _)
    sum(persistence .^ (0:octaves-1))
end

"""
    billow_fractal_1d(; kwargs...)

Construct a sampler that outputs a 1-dimensional billow fractal noise when it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

  - `source::AbstractSampler=simplex_1d()`: A 1-dimensional sampler instance to use as the source of
    the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""
function billow_fractal_1d(;
    seed=nothing,
    source=simplex_1d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    Billow{1}(seed, source, octaves, frequency, lacunarity, persistence)
end

"""
    billow_fractal_2d(; kwargs...)

Construct a sampler that outputs a 2-dimensional billow fractal noise when it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

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
    seed=nothing,
    source=opensimplex2s_2d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    Billow{2}(seed, source, octaves, frequency, lacunarity, persistence)
end

"""
    billow_fractal_3d(; kwargs...)

Construct a sampler that outputs a 3-dimensional billow fractal noise when it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

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
    seed=nothing,
    source=opensimplex2s_3d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    Billow{3}(seed, source, octaves, frequency, lacunarity, persistence)
end

"""
    billow_fractal_4d(; kwargs...)

Construct a sampler that outputs a 4-dimensional billow fractal noise when it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

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
    seed=nothing,
    source=opensimplex2s_4d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    Billow{4}(seed, source, octaves, frequency, lacunarity, persistence)
end

function sample(sampler::Billow{N}, coords::Vararg{Real,N}) where {N}
    state = sampler.state
    offset = state.offset
    persistence = state.persistence
    lacunarity = state.lacunarity
    coords = coords .* state.frequency
    amplitude = 1.0
    result = 0.0
    for source in state.sources
        result += (abs(sample(source, coords...)) * 2 - 1) * amplitude
        amplitude *= persistence
        coords = coords .* lacunarity .+ offset
    end
    result / state.scale
end

### Multifractal

struct Multi{N,S,O} <: FractalSampler{N}
    state::FractalState{N,S,O}
end

@inline function Multi{N}(
    seed,
    source::S,
    octaves::Int,
    frequency,
    lacunarity,
    persistence,
) where {N,S<:AbstractSampler{N}}
    O = octaves
    fs = FractalState{N,Multi,O}(seed, source, frequency, lacunarity, persistence, 1.0)
    Multi{N,S,O}(fs)
end

@inline function scale_factor(::Type{Multi}, octaves, persistence, _)
    reduce(1:octaves-1, init=1) do result, i
        result += result * persistence^i
    end
end

"""
    multi_fractal_1d(; kwargs...)

Construct a sampler that outputs a 1-dimensional multifractal noise when it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

  - `source::AbstractSampler=simplex_1d()`: A 1-dimensional sampler instance to use as the source of
    the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""
function multi_fractal_1d(;
    seed=nothing,
    source=simplex_1d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    Multi{1}(seed, source, octaves, frequency, lacunarity, persistence)
end

"""
    multi_fractal_2d(; kwargs...)

Construct a sampler that outputs a 2-dimensional multifractal noise when it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

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
    source=opensimplex2s_2d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    Multi{2}(seed, source, octaves, frequency, lacunarity, persistence)
end

"""
    multi_fractal_3d(; kwargs...)

Construct a sampler that outputs a 3-dimensional multifractal noise when it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

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
    source=opensimplex2s_3d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    Multi{3}(seed, source, octaves, frequency, lacunarity, persistence)
end

"""
    multi_fractal_4d(; kwargs...)

Construct a sampler that outputs a 4-dimensional multifractal noise when it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

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
    source=opensimplex2s_4d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    Multi{4}(seed, source, octaves, frequency, lacunarity, persistence)
end

function sample(sampler::Multi{N}, coords::Vararg{Real,N}) where {N}
    state = sampler.state
    offset = state.offset
    sources = state.sources
    persistence = state.persistence
    lacunarity = state.lacunarity
    coords = coords .* state.frequency
    amplitude = 1.0
    result = sample(sources[1], coords...)
    for source in sources[2:end]
        amplitude *= persistence
        coords = coords .* lacunarity .+ offset
        result += sample(source, coords...) * result * amplitude
    end
    result / state.scale
end

### Hybrid

struct Hybrid{N,S,O} <: FractalSampler{N}
    state::FractalState{N,S,O}
end

@inline function Hybrid{N}(
    seed,
    source::S,
    octaves::Int,
    frequency,
    lacunarity,
    persistence,
) where {N,S<:AbstractSampler{N}}
    O = octaves
    fs = FractalState{N,Hybrid,O}(seed, source, frequency, lacunarity, persistence, 1.0)
    Hybrid{N,S,O}(fs)
end

@inline function scale_factor(::Type{Hybrid}, octaves, persistence, _)
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

"""
    hybrid_fractal_1d(; kwargs...)

Construct a sampler that outputs a 1-dimensional hybrid multifractal noise when it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

  - `source::AbstractSampler=simplex_1d()`: A 1-dimensional sampler instance to use as the source of
    the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.25`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""
function hybrid_fractal_1d(;
    seed=nothing,
    source=simplex_1d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.25,
)
    Hybrid{1}(seed, source, octaves, frequency, lacunarity, persistence)
end

"""
    hybrid_fractal_2d(; kwargs...)

Construct a sampler that outputs a 2-dimensional hybrid multifractal noise when it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

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
    Hybrid{2}(seed, source, octaves, frequency, lacunarity, persistence)
end

"""
    hybrid_fractal_3d(; kwargs...)

Construct a sampler that outputs a 3-dimensional hybrid multifractal noise when it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

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
    Hybrid{3}(seed, source, octaves, frequency, lacunarity, persistence)
end

"""
    hybrid_fractal_4d(; kwargs...)

Construct a sampler that outputs a 4-dimensional hybrid multifractal noise when it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

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
    Hybrid{4}(seed, source, octaves, frequency, lacunarity, persistence)
end

function sample(sampler::Hybrid{N}, coords::Vararg{Real,N}) where {N}
    state = sampler.state
    offset = state.offset
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
        weight *= sample(source, (coords .+ offset)...) * amplitude
        result += weight
    end
    result / state.scale
end

### Ridged

struct Ridged{N,S,O} <: FractalSampler{N}
    state::FractalState{N,S,O}
end

@inline function Ridged{N}(
    seed,
    source::S,
    octaves::Int,
    frequency,
    lacunarity,
    persistence,
    attenuation,
) where {N,S<:AbstractSampler{N}}
    O = octaves
    fs = FractalState{N,Ridged,O}(seed, source, frequency, lacunarity, persistence, attenuation)
    Ridged{N,S,O}(fs)
end

@inline function scale_factor(::Type{Ridged}, octaves, persistence, attenuation)
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

"""
    ridged_fractal_1d(; kwargs...)

Construct a sampler that outputs a 1-dimensional ridged multifractal noise when it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

  - `source::AbstractSampler=simplex_1d()`: A 1-dimensional sampler instance to use as the source of
    the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=1.0`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.

   - `attenuation=2.0`: The attenuation to apply to the weight of each octave.
"""
function ridged_fractal_1d(;
    seed=nothing,
    source=simplex_1d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=1.0,
    attenuation=2.0,
)
    Ridged{1}(seed, source, octaves, frequency, lacunarity, persistence, attenuation)
end

"""
    ridged_fractal_2d(; kwargs...)

Construct a sampler that outputs a 2-dimensional ridged multifractal noise when it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

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
    Ridged{2}(seed, source, octaves, frequency, lacunarity, persistence, attenuation)
end

"""
    ridged_fractal_3d(; kwargs...)

Construct a sampler that outputs a 3-dimensional ridged multifractal noise when it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

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
    Ridged{3}(seed, source, octaves, frequency, lacunarity, persistence, attenuation)
end

"""
    ridged_fractal_4d(; kwargs...)

Construct a sampler that outputs a 4-dimensional ridged multifractal noise when it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

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
    Ridged{4}(seed, source, octaves, frequency, lacunarity, persistence, attenuation)
end

function sample(sampler::Ridged{N}, coords::Vararg{Real,N}) where {N}
    state = sampler.state
    offset = state.offset
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
        coords = coords .* lacunarity .+ offset
        weight = clamp(temp / attenuation, 0, 1)
        result += temp
    end
    result * state.scale - 1
end
