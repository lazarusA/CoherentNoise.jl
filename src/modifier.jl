# abs

struct Abs{N,S} <: ModifierSampler{N}
    random_state::RandomState
    source::S
end

"""
    abs(x::AbstractSampler)

Construct a modifier sampler that outputs the absolute value of its source when it is sampled from.
"""
function Base.abs(x::S) where {N,S<:AbstractSampler{N}}
    Abs{N,S}(x.random_state, x)
end

@inline function sample(sampler::Abs{N}, coords::Vararg{Real,N}) where {N}
    abs(sample(sampler.source, coords...))
end

# add

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
    Add{N,S1,S2}(x.random_state, x, y)
end

@inline function sample(
    sampler::Add{N,S1,S2},
    coords::Vararg{Real,N},
) where {N,N1,N2,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2}}
    x = sample(sampler.source1, coords[1:N1]...)
    y = sample(sampler.source2, coords[1:N2]...)
    x + y
end

"""
    +(x::AbstractSampler, y::Real)

Construct a modifier sampler that outputs the sum of the output of sampler `x` and the scalar `y`.
"""
function Base.:+(x::S, y::Real) where {N,S<:AbstractSampler{N}}
    Add{N,S,Float64}(x.random_state, x, Float64(y))
end

@inline function sample(
    sampler::Add{N,S,<:Real},
    coords::Vararg{Real,N},
) where {N,S<:AbstractSampler{N}}
    sample(sampler.source1, coords...) + sampler.source2
end

# cache

mutable struct Cache{N,S} <: ModifierSampler{N}
    random_state::RandomState
    source::S
    is_cached::Bool
    coords::NTuple{N,Float64}
    value::Float64
end

"""
    cache(x::AbstractSampler)

Construct a modifier sampler that caches the set of the input coordinates and their corresponding
output value of its source sampler. If the input coordinates differs from the previously cached
output, the cache is invalidated and the new output is cached.

Caching is useful if a sampler is used as a source for multiple modifiers. Without caching, the
duplicated input sources would redundantly compute the same outputs, which would be expensive,
especially if long pipelines share a long subgraph.
"""
function cache(x::S) where {N,S<:AbstractSampler{N}}
    Cache{N,S}(x.random_state, x, false, ntuple(i -> 0.0, N), 0.0)
end

@inline function sample(sampler::Cache{N}, coords::Vararg{Real,N}) where {N}
    if !sampler.is_cached || coords != sampler.coords
        sampler.is_cached = true
        sampler.coords = coords
        sampler.value = sample(sampler.source, coords...)
    end
    sampler.value
end

# clamp

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
    Clamp{N,S,L,H}(x.random_state, x, lo, hi)
end

@inline function sample(
    sampler::Clamp{N,S,L,H},
    coords::Vararg{Real,N},
) where {N,N1,N2,S,L<:AbstractSampler{N1},H<:AbstractSampler{N2}}
    lo = sample(sampler.lo, coords[1:min(N, N1)]..., ntuple(i -> 0.0, max(0, N1 - N))...)
    hi = sample(sampler.hi, coords[1:min(N, N2)]..., ntuple(i -> 0.0, max(0, N2 - N))...)
    clamp(sample(sampler.source, coords...), lo, hi)
end

"""
    clamp(x::AbstractSampler, lo=-1.0, hi=1.0)

Construct a modifier sampler that clamps the output of sampler `x` to be within the range of of the
scalars `lo` and `hi`.
"""
function Base.clamp(x::S, lo=-1.0, hi=1.0) where {N,S<:AbstractSampler{N}}
    Clamp{N,S,Float64,Float64}(x.random_state, x, Float64(lo), Float64(hi))
end

@inline function sample(sampler::Clamp{N,S,Real,Real}, coords::Vararg{Real,N}) where {N,S}
    clamp(sample(sampler.source, coords...), sampler.lo, sampler.hi)
end

# copysign

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

# curve

struct Curve{N,S} <: ModifierSampler{N}
    random_state::RandomState
    source::S
    points::Vector{Pair{Float64,Float64}}
end

"""
    curve(x::AbstractSampler, points::Vector{Pair{Float64,Float64}})

Construct a modifier sampler that outputs the result of sampling from `x` after remapping its output
to an arbitrary curve.

The curve is defined by a `Vector` of `Pair`s given by `points`. Each pair of points represents an
input and output number. The curve is a cubic spline, and so `points` must contain a list of four
point pairs at a minimum. Additionally, no two point pairs can contain the same input point value.

When sampling from sampler `x`, the output is evaluated using the curve data, and maps it to a new
output value.
"""
function curve(x::S, points::Vector{Pair{Float64,Float64}}) where {N,S<:AbstractSampler{N}}
    Curve{N,S}(x.random_state, x, sort(points, by=(p) -> p.first))
end

function sample(sampler::Curve{N}, coords::Vararg{Real,N}) where {N}
    points = sampler.points
    len = length(points)
    x = sample(sampler.source, coords...)
    p = findfirst(>(x) ∘ first, points)
    i = isnothing(p) ? len : clamp(p, 3:len)
    i1 = clamp(i - 2, 1:len)
    i2 = clamp(i - 1, 1:len)
    i3 = clamp(i, 1:len)
    i4 = clamp(i + 1, 1:len)
    if i2 == i3
        points[i2].second
    else
        in1 = points[i2].first
        in2 = points[i3].first
        a, b, c, d = (x -> x.second).((points[i1], points[i2], points[i3], points[i4]))
        cubic_interpolate(a, b, c, d, x - in1 / in2 - in1)
    end
end

# div

struct Div{N,S1,S2} <: ModifierSampler{N}
    random_state::RandomState
    source1::S1
    source2::S2
end

"""
    /(x::AbstractSampler, y::AbstractSampler)

Construct a modifier sampler that performs division of the output of sampler `x` by the output of
sampler `y`.
"""
function Base.:/(x::S1, y::S2) where {N1,N2,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2}}
    N = max(N1, N2)
    Div{N,S1,S2}(x.random_state, x, y)
end

@inline function sample(
    sampler::Div{N,S1,S2},
    coords::Vararg{Real,N},
) where {N,N1,N2,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2}}
    x = sample(sampler.source1, coords[1:N1]...)
    y = sample(sampler.source2, coords[1:N2]...)
    iszero(y) ? 0.0 : x / y
end

"""
    /(x::AbstractSampler, y::Real)

Construct a modifier sampler that performs division of the output of sampler `x` by the scalar `y`.
"""
function Base.:/(x::S, y::Real) where {N,S<:AbstractSampler{N}}
    Div{N,S,Float64}(x.random_state, x, Float64(y))
end

@inline function sample(
    sampler::Div{N,S,<:Real},
    coords::Vararg{Real,N},
) where {N,S<:AbstractSampler{N}}
    sample(sampler.source1, coords...) / sampler.source2
end

# max

struct Max{N,S1,S2} <: ModifierSampler{N}
    random_state::RandomState
    source1::S1
    source2::S2
end

"""
    max(x::AbstractSampler, y::AbstractSampler)

Construct a modifier sampler that outputs the maximum value of the outputs of sampling from samplers
`x` and `y`.
"""
function Base.max(x::S1, y::S2) where {N1,N2,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2}}
    N = max(N1, N2)
    Max{N,S1,S2}(x.random_state, x, y)
end

@inline function sample(
    sampler::Max{N,S1,S2},
    coords::Vararg{Real,N},
) where {N,N1,N2,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2}}
    x = sample(sampler.source1, coords[1:N1]...)
    y = sample(sampler.source2, coords[1:N2]...)
    max(x, y)
end

# min

struct Min{N,S1,S2} <: ModifierSampler{N}
    random_state::RandomState
    source1::S1
    source2::S2
end

"""
    min(x::AbstractSampler, y::AbstractSampler)

Construct a modifier sampler that outputs the minimum value of the outputs of sampling from samplers
`x` and `y`.
"""
function Base.min(x::S1, y::S2) where {N1,N2,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2}}
    N = max(N1, N2)
    Min{N,S1,S2}(x.random_state, x, y)
end

@inline function sample(
    sampler::Min{N,S1,S2},
    coords::Vararg{Real,N},
) where {N,N1,N2,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2}}
    x = sample(sampler.source1, coords[1:N1]...)
    y = sample(sampler.source2, coords[1:N2]...)
    min(x, y)
end

# mix

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
    mix(x::AbstractSampler, y::AbstractSampler, t::Real)

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

# mul

struct Mul{N,S1,S2} <: ModifierSampler{N}
    random_state::RandomState
    source1::S1
    source2::S2
end

"""
    *(x::AbstractSampler, y::AbstractSampler)

Construct a modifier sampler that outputs the product of the outputs samplers `x` and `y`.
"""
function Base.:*(x::S1, y::S2) where {N1,N2,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2}}
    N = max(N1, N2)
    Mul{N,S1,S2}(x.random_state, x, y)
end

@inline function sample(
    sampler::Mul{N,S1,S2},
    coords::Vararg{Real,N},
) where {N,N1,N2,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2}}
    x = sample(sampler.source1, coords[1:N1]...)
    y = sample(sampler.source2, coords[1:N2]...)
    x * y
end

"""
    *(x::AbstractSampler, y::Real)

Construct a modifier sampler that outputs the product of the output of sampler `x` by the scalar
`y`.
"""
function Base.:*(x::S, y::Real) where {N,S<:AbstractSampler{N}}
    Mul{N,S,Float64}(x.random_state, x, Float64(y))
end

@inline function sample(
    sampler::Mul{N,S,<:Real},
    coords::Vararg{Real,N},
) where {N,S<:AbstractSampler{N}}
    sample(sampler.source1, coords...) * sampler.source2
end

# muladd

struct Muladd{N,S} <: ModifierSampler{N}
    random_state::RandomState
    source::S
    strength::Float64
    bias::Float64
end

"""
    muladd(x::AbstractSampler, strength=1.0, bias=0.0)

Construct a modifier sampler that performs multiplies the output of sampler `x` by the scalar
`strength`, followed by adding the scalar `bias`. sampler `y`.
"""
function Base.muladd(x::S, strength=1.0, bias=0.0) where {N,S<:AbstractSampler{N}}
    Muladd{N,S}(x.random_state, x, Float64(strength), Float64(bias))
end

@inline function sample(
    sampler::Muladd{N,S},
    coords::Vararg{Real,N},
) where {N,S<:AbstractSampler{N}}
    x = sample(sampler.source, coords...)
    muladd(x, sampler.strength, sampler.bias)
end

# pow

struct Pow{N,S} <: ModifierSampler{N}
    random_state::RandomState
    source::S
    power::Float64
end

"""
    ^(x::AbstractSampler, y::Real)

Construct a modifier sampler that raises the output of sampler `x` to the power of the scalar `y`.
"""
Base.:^(x::S, y::Real) where {N,S<:AbstractSampler{N}} = Pow{N,S}(x.random_state, x, Float64(y))

@inline function sample(sampler::Pow{N}, coords::Vararg{Real,N}) where {N}
    x = sample(sampler.source, coords...)
    abs((x + 1) * 0.5)^sampler.power * 2 - 1
end

# rotate

struct Rotate{N,S} <: ModifierSampler{N}
    random_state::RandomState
    source::S
    x::NTuple{3,Float64}
    y::NTuple{3,Float64}
    z::NTuple{3,Float64}
end

"""
    rotate(source::AbstractSampler; x=0.0, y=0.0, z=0.0)

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

# scale

struct Scale{N,S} <: ModifierSampler{N}
    random_state::RandomState
    source::S
    axes::NTuple{N,Float64}
end

"""
    scale(source::AbstractSampler; x=1.0, y=1.0, z=1.0, w=1.0)

Construct a modifier sampler that scales the input coordinates of sampler `source` before sampling
from it.

Each axis can be scaled independently with `x`, `y`, `z`, or `w`.
"""
function scale(source::S; x=1.0, y=1.0, z=1.0, w=1.0) where {N,S<:AbstractSampler{N}}
    Scale{N,S}(source.random_state, source, Float64.((x, y, z, w))[1:N])
end

"""
    scale(source::AbstractSampler, scale::Real)

Construct a modifier sampler that uniformly scales the input coordinates of sampler `source` by the
scalar `scale` before sampling from it.
"""
function scale(source::S, scale::Real) where {N,S<:AbstractSampler{N}}
    Scale{N,S}(source.random_state, source, ntuple(i -> Float64(scale), N))
end

@inline function sample(sampler::Scale{N,S}, coords::Vararg{Real,N}) where {N,S<:AbstractSampler{N}}
    sample(sampler.source, (coords ./ sampler.axes)...)
end

# select

struct Select{N,S1,S2,C} <: ModifierSampler{N}
    random_state::RandomState
    source1::S1
    source2::S2
    control::C
    min::Float64
    max::Float64
    falloff::Float64
end

"""
    select(x::AbstractSampler, y::AbstractSampler; kwargs...)

Construct a modifier sampler that outputs either the out of sampler `x` or `y`, depending on the
output of sampler `z`.

If the output of sampler `z` is within the range denoted by `min` and `max`, the output of sampler
`y` is chosen. If the output of sampler `z` is outside of this range, the output of sampler `x` is
chosen.

# Arguments

  - `min`: A real number between -1.0 and 1.0 defining the lower bound of the selection range.

  - `max`: A real number between -1.0 and 1.0 defining the upper bound of the selection range.

  - `falloff`: A real number between 0.0 and 1.0 specifying the smoothness of the transition.
"""
function select(
    x::S1,
    y::S2,
    z::C;
    min=-1.0,
    max=1.0,
    falloff=0.0,
) where {N1,N2,N3,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2},C<:AbstractSampler{N3}}
    N = Base.max(N1, N2, N3)
    Select{N,S1,S2,C}(x.random_state, x, y, z, Float64(min), Float64(max), Float64(falloff))
end

@inline function sample(
    sampler::Select{N,S1,S2,C},
    coords::Vararg{Real,N},
) where {N,N1,N2,N3,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2},C<:AbstractSampler{N3}}
    min = sampler.min
    max = sampler.max
    falloff = Base.min(sampler.falloff, (max - min) * 0.5)
    s1 = sample(sampler.source1, coords[1:N1]...)
    s2 = sample(sampler.source2, coords[1:N2]...)
    c = sample(sampler.control, coords[1:N3]...)
    if falloff > 0
        minmf, maxmf = min - falloff, max - falloff
        minpf, maxpf = min + falloff, max + falloff
        if c < minmf
            s1
        elseif c < minpf
            lerp(s1, s2, curve3((c - minmf) / (minpf - minmf)))
        elseif c < maxpf
            lerp(s1, s2, curve3((c - maxmf) / (maxpf - maxmf)))
        else
            s1
        end
    else
        c < min || c > max ? s1 : s2
    end
end

# sub

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
    Sub{N,S1,S2}(x.random_state, x, y)
end

@inline function sample(
    sampler::Sub{N,S1,S2},
    coords::Vararg{Real,N},
) where {N,N1,N2,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2}}
    x = sample(sampler.source1, coords[1:N1]...)
    y = sample(sampler.source2, coords[1:N2]...)
    x - y
end

"""
    -(x::AbstractSampler, y::Real)

Construct a modifier sampler that outputs the difference of the output of sampler `x` and the scalar
`y`.
"""
function Base.:-(x::S, y::Real) where {N,S<:AbstractSampler{N}}
    Sub{N,S,Float64}(x.random_state, x, Float64(y))
end

@inline function sample(
    sampler::Sub{N,S,<:Real},
    coords::Vararg{Real,N},
) where {N,S<:AbstractSampler{N}}
    sample(sampler.source1, coords...) - sampler.source2
end

"""
    -(x::AbstractSampler)

Construct a modifier sampler that outputs the negated output of sampler `x`.
"""
Base.:-(x::S) where {N,S<:AbstractSampler{N}} = Sub{N,S,Nothing}(x.random_state, x, nothing)

@inline function sample(
    sampler::Sub{N,S,Nothing},
    coords::Vararg{Real,N},
) where {N,S<:AbstractSampler{N}}
    -sample(sampler.source1, coords...)
end

# terrace

struct Terrace{N,S} <: ModifierSampler{N}
    random_state::RandomState
    source::S
    points::Vector{Float64}
    invert::Bool
end

"""
    terrace(x::AbstractSampler, points::Vector{Pair{Float64,Float64}}; invert=false)

Construct a modifier sampler that outputs the result of sampling from `x` after remapping its output
to a terrace-forming curve.

The curve is defined by a `Vector` of `Float64`s given by `points`. Each point represents an input
and output number.

When sampling from sampler `x`, the output is evaluated using the curve data, and maps it to a new
output value.

# Arguments

  - `invert`: Specify whether the curve is inverted between control points.
"""
function terrace(
    x::S,
    points::Vector{Float64};
    invert=false,
) where {N,S<:AbstractSampler{N}}
    Terrace{N,S}(x.random_state, x, sort(points), invert)
end

function sample(sampler::Terrace{N}, coords::Vararg{Real,N}) where {N}
    points = sampler.points
    len = length(points)
    x = sample(sampler.source, coords...)
    p = findfirst(≥(x), points)
    i = isnothing(p) ? len : p
    i1 = clamp(i - 1, 1:len)
    i2 = clamp(i, 1:len)
    if i1 == i2
        points[i2]
    else
        in1 = points[i1]
        in2 = points[i2]
        t = (x - in1) / (in2 - in1)
        if sampler.invert
            t = 1 - t
            in1, in2 = in2, in1
        end
        lerp(in1, in2, t^2)
    end
end

# translate

struct Translate{N,S} <: ModifierSampler{N}
    random_state::RandomState
    source::S
    axes::NTuple{N,Float64}
end

"""
    translate(source::AbstractSampler; x=0.0, y=0.0, z=0.0, w=0.0)

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

# turbulence

struct Turbulence{N,N1,N2,S1,S2} <: ModifierSampler{N}
    random_state::RandomState
    source1::S1
    source2::S2
    power::Float64
    x::NTuple{N2,Float64}
    y::NTuple{N2,Float64}
    z::NTuple{N2,Float64}
    w::NTuple{N2,Float64}
end

"""
    turbulence(s1::AbstractSampler; s2::AbstractSampler; kwargs...)

Construct a modifier sampler that displaces the input coordinates of sampler `s1` by the output of
sampler `s2` with a fractional Brownian motion fractal applied to it.

Sampler `s2`'s input coordinates are randomly generated using the seed of sampler `s1`.

# Arguments

  - `frequency=1.0`: The frequency of the fractal signal to apply to sampler `s2`.

  - `roughness=3`: The number of octaves of the fractal to apply to sampler `s2`.

  - `power=1.0`: A scaling factor that is applied to the displaced result before sampling from
    sampler `s1`.
"""
function turbulence(
    s1::S1,
    s2::S2;
    frequency=1.0,
    roughness=3,
    power=1.0,
) where {N,N1,S1<:AbstractSampler{N},S2<:AbstractSampler{N1}}
    rng = s1.random_state.rng
    seed = rand(rng, UInt64)
    N2 = min(N, N1)
    x = ntuple(i -> rand(rng, Float64), N2)
    y = ntuple(i -> rand(rng, Float64), N2)
    z = ntuple(i -> rand(rng, Float64), N2)
    w = ntuple(i -> rand(rng, Float64), N2)
    s3 = FBM{N1}(seed, s2, roughness, frequency, 2.0, 0.5)
    S3 = typeof(s3)
    Turbulence{N,N1,N2,S1,S3}(s1.random_state, s1, s3, Float64(power), x, y, z, w)
end

function sample(sampler::Turbulence{N,N1,N2}, coords::Vararg{Real,N}) where {N,N1,N2}
    s2 = sampler.source2
    power = sampler.power
    zeros = ntuple(i -> 0.0, max(0, N1 - N))
    tx = sample(s2, (sampler.x .+ coords[1:N2])..., zeros...) * power
    ty = sample(s2, (sampler.y .+ coords[1:N2])..., zeros...) * power
    tz = sample(s2, (sampler.z .+ coords[1:N2])..., zeros...) * power
    tw = sample(s2, (sampler.w .+ coords[1:N2])..., zeros...) * power
    sample(sampler.source1, (coords .+ (tx, ty, tz, tw)[1:N])...)
end

# warp

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

  - `x::AbstractSampler=constant_1d()`: A sampler to warp the X axis by.

  - `y::AbstractSampler=constant_1d()`: A sampler to warp the Y axis by.

  - `z::AbstractSampler=constant_1d()`: A sampler to warp the Z axis by.

  - `w::AbstractSampler=constant_1d()`: A sampler to warp the W axis by.
"""
function warp(
    source::S;
    x::SX=constant_1d(),
    y::SY=constant_1d(),
    z::SZ=constant_1d(),
    w::SW=constant_1d(),
) where {
    N,N1,N2,N3,N4,
    S<:AbstractSampler{N},
    SX<:AbstractSampler{N1},
    SY<:AbstractSampler{N2},
    SZ<:AbstractSampler{N3},
    SW<:AbstractSampler{N4},
}
    Warp{N,N1,N2,N3,N4,S,SX,SY,SZ,SW}(source.random_state, source, (x, y, z, w))
end

function sample(sampler::Warp{N,N1,N2,N3,N4}, coords::Vararg{Real,N}) where {N,N1,N2,N3,N4}
    x, y, z, w = sampler.displacement
    dx = sample(x, coords[1:min(N, N1)]..., ntuple(i -> 0.0, max(0, N1 - N))...)
    dy = sample(y, coords[1:min(N, N2)]..., ntuple(i -> 0.0, max(0, N2 - N))...)
    dz = sample(z, coords[1:min(N, N3)]..., ntuple(i -> 0.0, max(0, N3 - N))...)
    dw = sample(w, coords[1:min(N, N4)]..., ntuple(i -> 0.0, max(0, N4 - N))...)
    sample(sampler.source, (coords .+ (dx, dy, dz, dw)[1:N])...)
end
