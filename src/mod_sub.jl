struct Sub{N,S1,S2} <: ModifierSampler{N}
    random_state::RandomState
    source1::S1
    source2::S2
end

"""
    -(x::AbstractSampler, y::AbstractSampler)

Construct a modifier sampler that outputs the difference of the outputs of samplers `x` and `y`.
"""
function Base.:-(x::S1, y::S2) where {N1,N2,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2}}
    N = max(N1, N2)
    Sub{N,S1,S2}(random_state(x), x, y)
end

"""
    -(x::AbstractSampler, y::Real)

Construct a modifier sampler that outputs the difference of the output of sampler `x` and the scalar
`y`.
"""
function Base.:-(x::S, y::Real) where {N,S<:AbstractSampler{N}}
    Sub{N,S,Float64}(random_state(x), x, Float64(y))
end

"""
    -(x::AbstractSampler)

Construct a modifier sampler that outputs the negated output of sampler `x`.
"""
Base.:-(x::S) where {N,S<:AbstractSampler{N}} = Sub{N,S,Nothing}(random_state(x), x, nothing)

@inline function sample(
    sampler::Sub{N,S1,S2},
    coords::Vararg{Real,N},
) where {N,N1,N2,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2}}
    x = sample(sampler.source1, coords[1:N1]...)
    y = sample(sampler.source2, coords[1:N2]...)
    x - y
end

@inline function sample(
    sampler::Sub{N,S,<:Real},
    coords::Vararg{Real,N},
) where {N,S<:AbstractSampler{N}}
    sample(sampler.source1, coords...) - sampler.source2
end

@inline function sample(
    sampler::Sub{N,S,Nothing},
    coords::Vararg{Real,N},
) where {N,S<:AbstractSampler{N}}
    -sample(sampler.source1, coords...)
end
