struct FractalState{N,S,O}
    random_state::RandomState
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
    sources = ntuple(_ -> deepcopy(source), O)
    scale = scale_factor(F, O, persistence, attenuation)
    FractalState{N,S,O}(rs, sources, scale, frequency, lacunarity, persistence, attenuation)
end

@inline function Base.getproperty(obj::FractalSampler, name::Symbol)
    if name === :random_state
        obj.state.random_state
    elseif name === :rng
        obj.state.random_state.rng
    elseif name === :seed
        obj.state.random_state.seed
    else
        getfield(obj, name)
    end
end
