struct Clamp{N,S,Lo,Hi} <: ModifierSampler{N}
    random_state::RandomState
    source::S
    lo::Lo
    hi::Hi
end

"""
    clamp(x::AbstractSampler, lo::AbstractSampler, hi::AbstractSampler)

Construct a modifier sampler that clamps the output of sampler `x` to be within the range of of
output values from samplers `lo` and `hi`.
"""
function Base.clamp(
    x::S,
    lo::L,
    hi::H,
) where {N,N1,N2,S<:AbstractSampler{N},L<:AbstractSampler{N1},H<:AbstractSampler{N2}}
    Clamp{N,S,L,H}(random_state(x), x, lo, hi)
end

"""
    clamp(x::AbstractSampler, lo::Real=-1.0, hi::Real=1.0)

Construct a modifier sampler that clamps the output of sampler `x` to be within the range of of the
scalars `lo` and `hi`.
"""
function Base.clamp(x::S, lo::Real=-1.0, hi::Real=1.0) where {N,S<:AbstractSampler{N}}
    Clamp{N,S,Float64,Float64}(random_state(x), x, Float64(lo), Float64(hi))
end

@inline function sample(
    sampler::Clamp{N,S,L,H},
    coords::Vararg{Real,N},
) where {N,N1,N2,S,L<:AbstractSampler{N1},H<:AbstractSampler{N2}}
    lo = sample(sampler.lo, coords[1:min(N, N1)]..., ntuple(i -> 0.0, max(0, N1 - N))...)
    hi = sample(sampler.hi, coords[1:min(N, N2)]..., ntuple(i -> 0.0, max(0, N2 - N))...)
    clamp(sample(sampler.source, coords...), lo, hi)
end

@inline function sample(sampler::Clamp{N,S,Real,Real}, coords::Vararg{Real,N}) where {N,S}
    clamp(sample(sampler.source, coords...), sampler.lo, sampler.hi)
end
