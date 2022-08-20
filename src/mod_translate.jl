struct Translate{N,S} <: ModifierSampler{N}
    random_state::RandomState
    source::S
    axes::NTuple{N,Float64}
end

"""
    translate(source::AbstractSampler; x::Real=0.0, y::Real=0.0, z::Real=0.0, w::Real=0.0)

Construct a modifier sampler that translates the input coordinates of sampler `source` along a
specified vector, with a direction and magnitude given by the coordinates `x`, `y`, `z`, and `w`.
"""
function translate(source::S; x=0.0, y=0.0, z=0.0, w=0.0) where {N,S<:AbstractSampler{N}}
    Translate{N,S}(source.random_state, source, Float64.((x, y, z, w))[1:N])
end

@inline function sample(
    sampler::Translate{N,S},
    coords::Vararg{Real,N},
) where {N,S<:AbstractSampler{N}}
    sample(sampler.source, (coords .+ sampler.axes)...)
end
