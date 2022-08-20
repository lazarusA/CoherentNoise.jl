struct Mix{N,S1,S2,C} <: ModifierSampler{N}
    random_state::RandomState
    source1::S1
    source2::S2
    control::C
end

"""
    mix(x::AbstractSampler, y::AbstractSampler, t::AbstractSampler)

Construct a modifier sampler that outputs the result of linearly interpolating the output of
samplers `x` and `y` by the output of sampler `t`.
"""
function mix(
    a::S1,
    b::S2,
    t::C,
) where {N1,N2,N3,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2},C<:AbstractSampler{N3}}
    N = max(N1, N2, N3)
    Mix{N,S1,S2,C}(a.random_state, a, b, t)
end

"""
    mix(x::AbstractSampler, y::AbstractSampler, t::Scalar)

Construct a modifier sampler that outputs the result of linearly interpolating the output of
samplers `x` and `y` by the scalar `t`.
"""
function mix(a::S1, b::S2, t::Real) where {N1,N2,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2}}
    N = max(N1, N2)
    Mix{N,S1,S2,Float64}(a.random_state, a, b, clamp(t, 0.0, 1.0))
end

@inline function sample(
    sampler::Mix{N,S1,S2,C},
    coords::Vararg{Real,N},
) where {N,N1,N2,N3,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2},C<:AbstractSampler{N3}}
    a = sample(sampler.source1, coords[1:N1]...)
    b = sample(sampler.source2, coords[1:N2]...)
    t = (sample(sampler.control, coords[1:N3]...) + 1) * 0.5
    lerp(a, b, t)
end

@inline function sample(
    sampler::Mix{N,S1,S2,C},
    coords::Vararg{Real,N},
) where {N,N1,N2,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2},C<:Real}
    a = sample(sampler.source1, coords[1:N1]...)
    b = sample(sampler.source2, coords[1:N2]...)
    lerp(a, b, sampler.control)
end
