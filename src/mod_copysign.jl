struct CopySign{N,S1,S2} <: ModifierSampler{N}
    random_state::RandomState
    source1::S1
    source2::S2
end

"""
    copysign(x::AbstractSampler, y::AbstractSampler)

Construct a modifier sampler that outputs the value of sampling from `x` with the sign copied from
the value of sampling from `y`.
"""
function Base.copysign(x::S1, y::S2) where {N1,N2,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2}}
    N = max(N1, N2)
    CopySign{N,S1,S2}(x.random_state, x, y)
end

@inline function sample(
    sampler::CopySign{N,S1,S2},
    coords::Vararg{Real,N},
) where {N,N1,N2,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2}}
    x = sample(sampler.source1, coords[1:N1]...)
    y = sample(sampler.source2, coords[1:N2]...)
    copysign(x, y)
end
