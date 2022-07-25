struct Scale{N,S} <: ModifierSampler{N}
    random_state::RandomState
    source::S
    axes::NTuple{N,Float64}
end

"""
    scale(source::AbstractSampler; x::Real=1.0, y::Real=1.0, z::Real=1.0, w::Real=1.0)

Construct a modifier sampler that scales the input coordinates of sampler `source` before sampling
from it.

Each axis can be scaled independently with `x`, `y`, `z`, or `w`.
"""
function scale(
    source::S;
    x::Real=1.0,
    y::Real=1.0,
    z::Real=1.0,
    w::Real=1.0,
) where {N,S<:AbstractSampler{N}}
    Scale{N,S}(random_state(source), source, Float64.((x, y, z, w))[1:N])
end

"""
    scale(source::AbstractSampler, scale::Real=1.0)

Construct a modifier sampler that uniformly scales the input coordinates of sampler `source` by the
scalar `scale` before sampling from it.
"""
function scale(source::S, scale::Real) where {N,S<:AbstractSampler{N}}
    Scale{N,S}(random_state(source), source, ntuple(i -> Float64(scale), N))
end

@inline function sample(sampler::Scale{N,S}, coords::Vararg{Real,N}) where {N,S<:AbstractSampler{N}}
    sample(sampler.source, (coords ./ sampler.axes)...)
end
