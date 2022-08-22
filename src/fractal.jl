### State composed into all fractal sampler types

struct FractalState{N,S,O}
    sources::NTuple{O,S}
    scale::Float64
    frequency::Float64
    lacunarity::Float64
    persistence::Float64
    attenuation::Float64
end

@inline function FractalState{N,F,O}(
    source::S,
    frequency,
    lacunarity,
    persistence,
    attenuation,
) where {N,F<:FractalSampler,O,S<:AbstractSampler{N}}
    sources = ntuple(_ -> deepcopy(source), O)
    scale = scale_factor(F, O, persistence, attenuation)
    FractalState{N,S,O}(sources, scale, frequency, lacunarity, persistence, attenuation)
end

### fBm

struct FBMFractal{N,S,O} <: FractalSampler{N}
    random_state::RandomState
    state::FractalState{N,S,O}
end

@inline function FBMFractal{N}(
    seed,
    source::S,
    octaves,
    frequency,
    lacunarity,
    persistence,
) where {N,S<:AbstractSampler{N}}
    O = octaves
    rs = RandomState(seed)
    fs = FractalState{N,FBMFractal,O}(source, frequency, lacunarity, persistence, 1.0)
    FBMFractal{N,S,O}(rs, fs)
end

@inline function scale_factor(::Type{FBMFractal}, octaves, persistence, _)
    sum(persistence .^ (0:octaves-1))
end

@doc doc_fbm_fractal_2d
function fbm_fractal_2d(;
    seed=0,
    source=opensimplex2s_2d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    FBMFractal{2}(seed, source, octaves, frequency, lacunarity, persistence)
end

@doc doc_fbm_fractal_3d
function fbm_fractal_3d(;
    seed=0,
    source=opensimplex2s_3d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    FBMFractal{3}(seed, source, octaves, frequency, lacunarity, persistence)
end

@doc doc_fbm_fractal_4d
function fbm_fractal_4d(;
    seed=0,
    source=opensimplex2s_4d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    FBMFractal{4}(seed, source, octaves, frequency, lacunarity, persistence)
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

### Billow

struct BillowFractal{N,S,O} <: FractalSampler{N}
    random_state::RandomState
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
    rs = RandomState(seed)
    fs = FractalState{N,BillowFractal,O}(source, frequency, lacunarity, persistence, 1.0)
    BillowFractal{N,S,O}(rs, fs)
end

@inline function scale_factor(::Type{BillowFractal}, octaves, persistence, _)
    sum(persistence .^ (0:octaves-1))
end

@doc doc_billow_fractal_2d
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

@doc doc_billow_fractal_3d
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

@doc doc_billow_fractal_4d
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

### Multifractal

struct MultiFractal{N,S,O} <: FractalSampler{N}
    random_state::RandomState
    state::FractalState{N,S,O}
end

@inline function MultiFractal{N}(
    seed,
    source::S,
    octaves::Int,
    frequency,
    lacunarity,
    persistence,
) where {N,S<:AbstractSampler{N}}
    O = octaves
    rs = RandomState(seed)
    fs = FractalState{N,MultiFractal,O}(source, frequency, lacunarity, persistence, 1.0)
    MultiFractal{N,S,O}(rs, fs)
end

@inline function scale_factor(::Type{MultiFractal}, octaves, persistence, _)
    reduce(1:octaves-1, init=1) do result, i
        result += result * persistence^i
    end
end

@doc doc_multi_fractal_2d
function multi_fractal_2d(;
    seed=0,
    source=opensimplex2_2d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    MultiFractal{2}(seed, source, octaves, frequency, lacunarity, persistence)
end

@doc doc_multi_fractal_3d
function multi_fractal_3d(;
    seed=0,
    source=opensimplex2_3d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    MultiFractal{3}(seed, source, octaves, frequency, lacunarity, persistence)
end

@doc doc_multi_fractal_4d
function multi_fractal_4d(;
    seed=0,
    source=opensimplex2_4d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.5,
)
    MultiFractal{4}(seed, source, octaves, frequency, lacunarity, persistence)
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

### Hybrid

struct HybridFractal{N,S,O} <: FractalSampler{N}
    random_state::RandomState
    state::FractalState{N,S,O}
end

@inline function HybridFractal{N}(
    seed,
    source::S,
    octaves::Int,
    frequency,
    lacunarity,
    persistence,
) where {N,S<:AbstractSampler{N}}
    O = octaves
    rs = RandomState(seed)
    fs = FractalState{N,HybridFractal,O}(source, frequency, lacunarity, persistence, 1.0)
    HybridFractal{N,S,O}(rs, fs)
end

@inline function scale_factor(::Type{HybridFractal}, octaves, persistence, _)
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

@doc doc_hybrid_fractal_2d
function hybrid_fractal_2d(;
    seed=0,
    source=opensimplex2s_2d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.25,
)
    HybridFractal{2}(seed, source, octaves, frequency, lacunarity, persistence)
end

@doc doc_hybrid_fractal_3d
function hybrid_fractal_3d(;
    seed=0,
    source=opensimplex2s_3d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.25,
)
    HybridFractal{3}(seed, source, octaves, frequency, lacunarity, persistence)
end

@doc doc_hybrid_fractal_4d
function hybrid_fractal_4d(;
    seed=0,
    source=opensimplex2s_4d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=0.25,
)
    HybridFractal{4}(seed, source, octaves, frequency, lacunarity, persistence)
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

### Ridged

struct RidgedFractal{N,S,O} <: FractalSampler{N}
    random_state::RandomState
    state::FractalState{N,S,O}
end

@inline function RidgedFractal{N}(
    seed,
    source::S,
    octaves::Int,
    frequency,
    lacunarity,
    persistence,
    attenuation,
) where {N,S<:AbstractSampler{N}}
    O = octaves
    rs = RandomState(seed)
    fs = FractalState{N,RidgedFractal,O}(source, frequency, lacunarity, persistence, attenuation)
    RidgedFractal{N,S,O}(rs, fs)
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

@doc doc_ridged_fractal_2d
#  - `attenuation=2.0`: The attenuation to apply to the weight of each octave.
function ridged_fractal_2d(;
    seed=0,
    source=opensimplex2s_2d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=1.0,
    attenuation=2.0,
)
    RidgedFractal{2}(seed, source, octaves, frequency, lacunarity, persistence, attenuation)
end

@doc doc_ridged_fractal_3d
function ridged_fractal_3d(;
    seed=0,
    source=opensimplex2s_3d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=1.0,
    attenuation=2.0,
)
    RidgedFractal{3}(seed, source, octaves, frequency, lacunarity, persistence, attenuation)
end

@doc doc_ridged_fractal_4d
function ridged_fractal_4d(;
    seed=0,
    source=opensimplex2s_4d(seed=seed),
    octaves=4,
    frequency=1.0,
    lacunarity=2.0,
    persistence=1.0,
    attenuation=2.0,
)
    RidgedFractal{4}(seed, source, octaves, frequency, lacunarity, persistence, attenuation)
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
