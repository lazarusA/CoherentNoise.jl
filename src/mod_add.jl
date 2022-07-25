struct Add{N,S1,S2} <: ModifierSampler{N}
    random_state::RandomState
    source1::S1
    source2::S2
end

"""
    +(x::AbstractSampler, y::AbstractSampler)

Construct a modifier sampler that outputs the sum of the outputs of samplers `x` and `y`.
"""
function Base.:+(x::S1, y::S2) where {N1,N2,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2}}
    N = max(N1, N2)
    Add{N,S1,S2}(random_state(x), x, y)
end

"""
    +(x::AbstractSampler, y::Real)

Construct a modifier sampler that outputs the sum of the output of sampler `x` and the scalar `y`.
"""
function Base.:+(x::S, y::Real) where {N,S<:AbstractSampler{N}}
    Add{N,S,Float64}(random_state(x), x, Float64(y))
end

@inline function sample(
    sampler::Add{N,S1,S2},
    coords::Vararg{Real,N},
) where {N,N1,N2,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2}}
    x = sample(sampler.source1, coords[1:N1]...)
    y = sample(sampler.source2, coords[1:N2]...)
    x + y
end

@inline function sample(
    sampler::Add{N,S,<:Real},
    coords::Vararg{Real,N},
) where {N,S<:AbstractSampler{N}}
    sample(sampler.source1, coords...) + sampler.source2
end
