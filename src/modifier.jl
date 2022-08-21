# abs

struct Abs{N,S} <: ModifierSampler{N}
    random_state::RandomState
    source::S
end

@doc doc_mod_abs
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

@doc doc_mod_add
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

@doc doc_mod_add_scalar
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

@doc doc_mod_cache
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

@doc doc_mod_clamp
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

@doc doc_mod_clamp_scalar
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

@doc doc_mod_copysign
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

@doc doc_mod_curve
function curve(x::S, points::Vector{Pair{Float64,Float64}}) where {N,S<:AbstractSampler{N}}
    Curve{N,S}(x.random_state, x, sort(points, by=(p) -> p.first))
end

function sample(sampler::Curve{N}, coords::Vararg{Real,N}) where {N}
    points = sampler.points
    len = length(points)
    x = sample(sampler.source, coords...)
    p = findfirst(>(x) ∘ first, points)
    i = p !== nothing ? clamp(p, 3:len) : len
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

@doc doc_mod_div
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

@doc doc_mod_div_scalar
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

@doc doc_mod_max
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

@doc doc_mod_min
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

@doc doc_mod_mix
function mix(
    a::S1,
    b::S2,
    t::C,
) where {N1,N2,N3,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2},C<:AbstractSampler{N3}}
    N = max(N1, N2, N3)
    Mix{N,S1,S2,C}(a.random_state, a, b, t)
end

@doc doc_mod_mix_scalar
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

@doc doc_mod_mul
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

@doc doc_mod_mul_scalar
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

@doc doc_mod_muladd
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

@doc doc_mod_pow
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

@doc doc_mod_rotate
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

@doc doc_mod_scale
function scale(source::S; x=1.0, y=1.0, z=1.0, w=1.0) where {N,S<:AbstractSampler{N}}
    Scale{N,S}(source.random_state, source, Float64.((x, y, z, w))[1:N])
end

@doc doc_mod_scale_scalar
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

@doc doc_mod_select
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

@doc doc_mod_sub
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

@doc doc_mod_sub_scalar
function Base.:-(x::S, y::Real) where {N,S<:AbstractSampler{N}}
    Sub{N,S,Float64}(x.random_state, x, Float64(y))
end

@inline function sample(
    sampler::Sub{N,S,<:Real},
    coords::Vararg{Real,N},
) where {N,S<:AbstractSampler{N}}
    sample(sampler.source1, coords...) - sampler.source2
end

@doc doc_mod_sub_unary
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

@doc doc_mod_terrace
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
    i = p !== nothing ? p : len
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

@doc doc_mod_translate
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

@doc doc_mod_turbulence
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
    s3 = FBMFractal{N1}(seed=seed, source=s2, octaves=roughness, frequency=frequency)
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

@doc doc_mod_warp
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
