struct Rotate{N,S} <: ModifierSampler{N}
    random_state::RandomState
    source::S
    x::NTuple{3,Float64}
    y::NTuple{3,Float64}
    z::NTuple{3,Float64}
end

"""
    rotate(source::AbstractSampler; x::Real=0.0, y::Real=0.0, z::Real=0.0)

Construct a modifier sampler that rotates the input coordinates of sampler `source` around the
origin before sampling from it.

The coordinate system is assumed to be left-handed.

The angle of rotation is specified in radians for the corresponding axis given by `x`, `y`, and `z`.
"""
function rotate(source::S; x=0.0, y=0.0, z=0.0) where {N,S<:AbstractSampler{N}}
    cx, cy, cz = cos.((x, y, z))
    sx, sy, sz = sin.((x, y, z))
    rx = (sx * sy * sz + cy * cz, cx * sz, sy * cz - cy * sx * sz)
    ry = (sx * sy * cz - cy * sz, cx * cz, -cy * sx * cz - sy * sz)
    rz = (-sy * cx, sx, cy * cx)
    Rotate{N,S}(source.random_state, source, rx, ry, rz)
end

@inline function sample(sampler::Rotate{2,S}, coords::Vararg{Real,2}) where {S<:AbstractSampler{2}}
    x = sum((coords .* sampler.x[1:2]))
    y = sum((coords .* sampler.y[1:2]))
    sample(sampler.source, x, y)
end

@inline function sample(sampler::Rotate{3,S}, coords::Vararg{Real,3}) where {S<:AbstractSampler{3}}
    x = sum((coords .* sampler.x))
    y = sum((coords .* sampler.y))
    z = sum((coords .* sampler.z))
    sample(sampler.source, x, y, z)
end

@inline function sample(sampler::Rotate{4,S}, coords::Vararg{Real,4}) where {S<:AbstractSampler{4}}
    x = sum((coords[1:3] .* sampler.x))
    y = sum((coords[1:3] .* sampler.y))
    z = sum((coords[1:3] .* sampler.z))
    sample(sampler.source, x, y, z, coords[4])
end
