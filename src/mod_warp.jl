struct Warp{N,N1,N2,N3,N4,S,SX,SY,SZ,SW} <: ModifierSampler{N}
    random_state::RandomState
    source::S
    displacement::Tuple{SX,SY,SZ,SW}
end

"""
    warp(source::AbstractSampler; kwargs...)

Construct a modifier sampler that performs domain warping of the sampler `source` before sampling
from it.

Domain warping feeds the output of other samplers to the input of a sampler. For this modifier, each
input coordinate can specify a different sampler to warp with.

If a sampler is not supplied for `x`, `y`, `z`, or `w`, a sampler that outputs a constant zero will
be used instead.

# Arguments

  - `x::AbstractSampler=Constant()`: A sampler to warp the X axis by.
  - `y::AbstractSampler=Constant()`: A sampler to warp the Y axis by.
  - `z::AbstractSampler=Constant()`: A sampler to warp the Z axis by.
  - `w::AbstractSampler=Constant()`: A sampler to warp the W axis by.
"""
function warp(
    source::S;
    x::SX=Constant(),
    y::SY=Constant(),
    z::SZ=Constant(),
    w::SW=Constant(),
) where {
    N,N1,N2,N3,N4,
    S<:AbstractSampler{N},
    SX<:AbstractSampler{N1},
    SY<:AbstractSampler{N2},
    SZ<:AbstractSampler{N3},
    SW<:AbstractSampler{N4},
}
    Warp{N,N1,N2,N3,N4,S,SX,SY,SZ,SW}(random_state(source), source, (x, y, z, w))
end

function sample(sampler::Warp{N,N1,N2,N3,N4}, coords::Vararg{Real,N}) where {N,N1,N2,N3,N4}
    x, y, z, w = sampler.displacement
    dx = sample(x, coords[1:min(N, N1)]..., ntuple(i -> 0.0, max(0, N1 - N))...)
    dy = sample(y, coords[1:min(N, N2)]..., ntuple(i -> 0.0, max(0, N2 - N))...)
    dz = sample(z, coords[1:min(N, N3)]..., ntuple(i -> 0.0, max(0, N3 - N))...)
    dw = sample(w, coords[1:min(N, N4)]..., ntuple(i -> 0.0, max(0, N4 - N))...)
    sample(sampler.source, (coords .+ (dx, dy, dz, dw)[1:N])...)
end
