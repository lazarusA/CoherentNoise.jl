struct OpenSimplex2S{N,O<:Orientation} <: NoiseSampler{N}
    random_state::RandomState
    simplex_state::SimplexState
    table::Vector{Float64}
end

@inline function _opensimplex2s(dims, seed, table_size, gradients, orientation, smooth)
    rs = RandomState(seed)
    orientation = os2_orientation_type(Val(orientation))
    table = Iterators.take(Iterators.cycle(gradients), table_size) |> collect
    T = OpenSimplex2S{dims,orientation}
    T(rs, SimplexState(T, Val(smooth)), table)
end

SimplexState(::Type{<:OpenSimplex2S{2}}, ::Val) = SimplexState(2 / 3, 1.0)
SimplexState(::Type{<:OpenSimplex2S{3}}, ::Val{true}) = SimplexState(2 / 3, 1.9139664641132217)
SimplexState(::Type{<:OpenSimplex2S{3}}, ::Val{false}) = SimplexState(3 / 4, 1.0)
SimplexState(::Type{<:OpenSimplex2S{4}}, ::Val{true}) = SimplexState(2 / 3, 2.7405223931658242)
SimplexState(::Type{<:OpenSimplex2S{4}}, ::Val{false}) = SimplexState(4 / 5, 1.0)

# 2D

"""
    opensimplex2s_2d(; seed=nothing, orient=nothing)

Construct a sampler that outputs 2-dimensional OpenSimplex2S noise when it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

  - `orient`: Either the symbol `:x` or the value `nothing`:

      + `nothing`: Use the standard orientation.
      + `:x`: The noise space will be re-oriented with the Y axis pointing down the main diagonal to
        improve visual isotropy.
"""
function opensimplex2s_2d(; seed=nothing, orient=nothing)
    size = OS2_NUM_GRADIENTS_2D * 2
    gradients = OS2_GRADIENTS_NORMALIZED_2D ./ 0.05481866495625118
    _opensimplex2s(2, seed, size, gradients, orient, false)
end

@inline orient(::Type{OpenSimplex2S{2,OrientStandard}}, x, y) = (x, y) .+ OS2_SKEW_2D .* (x + y)

@inline function orient(::Type{OpenSimplex2S{2,OrientX}}, x, y)
    xx = x * ROOT_2_OVER_2
    yy = y * ROOT_2_OVER_2 * (2OS2_SKEW_2D + 1)
    (yy + xx, yy - xx)
end

@inline function contribute(seed, table, falloff, X, Y, x, y)
    a = falloff - x^2 - y^2
    a > 0 ? pow4(a) * grad(table, seed, X, Y, x, y) : 0.0
end

function sample(sampler::S, x::T, y::T) where {O,S<:OpenSimplex2S{2,O},T<:Real}
    seed = sampler.random_state.seed
    table = sampler.table
    state = sampler.simplex_state
    falloff = state.falloff
    primes = (PRIME_X, PRIME_Y)
    tr = orient(S, x, y)
    XY = floor.(Int, tr)
    xs, ys = tr .- XY
    d = xs - ys
    t = (xs + ys) * OS2_UNSKEW_2D
    us = OS2_UNSKEW_2D
    us2p1 = 2us + 1
    X1, Y1 = XY .* primes
    X2, Y2 = (X1, Y1) .+ primes
    x1, y1 = (xs, ys) .+ t
    x2, y2 = (x1, y1) .- us2p1
    a1 = falloff - x1^2 - y1^2
    c1 = pow4(a1) * grad(table, seed, X1, Y1, x1, y1)
    a2 = 2us2p1 * (1 / us + 2) * t + -2us2p1^2 + a1
    c2 = pow4(a2) * grad(table, seed, X2, Y2, x2, y2)
    result = c1 + c2
    if t < OS2_UNSKEW_2D
        X3, Y3 = (X1, Y1) .+ (primes .<< 1)
        if xs + d > 1
            result += contribute(seed, table, falloff, X3, Y2, x1 - 3us - 2, y1 - 3us - 1)
        else
            result += contribute(seed, table, falloff, X1, Y2, x1 - us, y1 - us - 1)
        end
        if ys - d > 1
            result += contribute(seed, table, falloff, X2, Y3, x1 - 3us - 1, y1 - 3us - 2)
        else
            result += contribute(seed, table, falloff, X2, Y1, x1 - us - 1, y1 - us)
        end
    else
        X4, Y4 = (X1, Y1) .- primes
        if xs + d < 0
            result += contribute(seed, table, falloff, X4, Y1, x1 + us + 1, y1 + us)
        else
            result += contribute(seed, table, falloff, X2, Y1, x1 - us - 1, y1 - us)
        end
        if ys < d
            result += contribute(seed, table, falloff, X1, Y4, x1 + us, y1 + us + 1)
        else
            result += contribute(seed, table, falloff, X1, Y2, x1 - us, y1 - us - 1)
        end
    end
    result * state.scale_factor
end

# 3D

"""
    opensimplex2s_3d(; seed=nothing, smooth=false, orient=nothing)

Construct a sampler that outputs 3-dimensional OpenSimplex2S noise when it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
    `nothing` for non-deterministic results.

  - `smooth`: Specify whether to have continuous gradients.
    Simplex variants, even the original Simplex noise by Ken Perlin, overshoot the radial extent for
    the signal reconstruction kernel in order to improve the visual of the noise. Normally this is
    okay, especially if layering multiple octaves of the noise. However, in some applications, such
    as creating height or bump maps, this will produce discontinuities visually identified by
    jarring creases in the generated noise.

    This option changes the falloff in order to produce smooth continuous noise, however, the
    resulting noise may look quite different than the non-smooth option, depending on the Simplex
    variant.

    The default value is `false`, in order to be true to the original implementation.

  - `orient`: Either the symbol `:x` or the value `nothing`:

      + `nothing`: Use the standard orientation.
      + `:x`: The noise space will be re-oriented with the Y axis pointing down the main diagonal to
        improve visual isotropy.
      + `:xy`: Re-orient the noise space to have better visual isotropy in the XY plane.
      + `:xz`: Re-orient the noise space to have better visual isotropy in the XZ plane.
"""
function opensimplex2s_3d(; seed=nothing, orient=nothing, smooth=false)
    size = OS2_NUM_GRADIENTS_3D * 4
    gradients = OS2_GRADIENTS_NORMALIZED_3D ./ 0.2781926117527186
    _opensimplex2s(3, seed, size, gradients, orient, smooth)
end

@inline function orient(::Type{OpenSimplex2S{3,OrientStandard}}, x, y, z)
    OS2_FALLBACK_ROTATE_3D * (x + y + z) .- (x, y, z)
end

@inline function orient(::Type{OpenSimplex2S{3,OrientXY}}, x, y, z)
    xy = x + y
    zz = z * ROOT_3_OVER_3
    xr, yr = (x, y) .+ xy .* OS2_ROTATE_3D_ORTHONORMALIZER .+ zz
    zr = xy * -ROOT_3_OVER_3 + zz
    (xr, yr, zr)
end

@inline function orient(::Type{OpenSimplex2S{3,OrientXZ}}, x, y, z)
    xz = x + z
    yy = y * ROOT_3_OVER_3
    xr, zr = (x, z) .+ xz .* -0.211324865405187 .+ yy
    yr = xz * -ROOT_3_OVER_3 + yy
    (xr, yr, zr)
end

@inline os2s_contribute1(seed, table, a, args...) = pow4(a) * grad(table, seed, args...)

@inline function os2s_contribute2(seed, table, a, args...)
    a > 0 ? os2s_contribute1(seed, table, a, args...) : 0.0
end

function sample(sampler::S, x::T, y::T, z::T) where {O,S<:OpenSimplex2S{3,O},T<:Real}
    seed = sampler.random_state.seed
    table = sampler.table
    seed2 = seed ⊻ -OS2_SEED_FLIP_3D
    state = sampler.simplex_state
    falloff = state.falloff
    primes = (PRIME_X, PRIME_Y, PRIME_Z)
    tr = orient(S, x, y, z)
    V = floor.(Int, tr)
    Vr = tr .- V
    Vp = V .* primes
    mask1 = trunc.(Int, -0.5 .- Vr)
    mask2 = mask1 .| 1
    X1, Y1, Z1 = Vp .+ mask1 .& primes
    X2, Y2, Z2 = Vp .+ primes
    X3, Y3, Z3 = Vp .+ .~mask1 .& primes
    X4, Y4, Z4 = Vp .+ mask1 .& primes .<< 1
    @inbounds X5 = Vp[1] + mask1[1] & PRIME_X * 2
    x1, y1, z1 = Vr .+ mask1
    x2, y2, z2 = Vr .- 0.5
    x3, y3, z3 = (x1, y1, z1) .- mask2
    x4, y4, z4 = (x2, y2, z2) .+ mask2
    a0 = falloff - x1^2 - y1^2 - z1^2
    a1 = falloff - x2^2 - y2^2 - z2^2
    xf1, yf1, zf1 = (x2, y2, z2) .* mask2 .<< 1
    xf2, yf2, zf2 = (x2, y2, z2) .* (-2 .- mask1 .<< 2) .- 1
    xf3, yf3, zf3 = (xf1, yf1, zf1) .+ a0
    xf4, yf4, zf4 = (xf2, yf2, zf2) .+ a1
    skip1, skip2, skip3 = false, false, false
    result = os2s_contribute1(seed, table, a0, X1, Y1, Z1, x1, y1, z1)
    result += os2s_contribute1(seed2, table, a1, X2, Y2, Z2, x2, y2, z2)
    if xf3 > 0
        result += os2s_contribute1(seed, table, xf3, X3, Y1, Z1, x3, y1, z1)
    else
        result += os2s_contribute2(seed, table, yf1 + zf1 + a0, X1, Y3, Z3, x1, y3, z3)
        if xf4 > 0
            result += os2s_contribute1(seed2, table, xf4, X5, Y2, Z2, x4, y2, z2)
            skip1 = true
        end
    end
    if yf3 > 0
        result += os2s_contribute1(seed, table, yf3, X1, Y3, Z1, x1, y3, z1)
    else
        result += os2s_contribute2(seed, table, xf1 + zf1 + a0, X3, Y1, Z3, x3, y1, z3)
        if yf4 > 0
            result += os2s_contribute1(seed2, table, yf4, X2, Y4, Z2, x2, y4, z2)
            skip2 = true
        end
    end
    if zf3 > 0
        result += os2s_contribute1(seed, table, zf3, X1, Y1, Z3, x1, y1, z3)
    else
        result += os2s_contribute2(seed, table, xf1 + yf1 + a0, X3, Y3, Z1, x3, y3, z1)
        if zf4 > 0
            result += os2s_contribute1(seed2, table, zf4, X2, Y2, Z4, x2, y2, z4)
            skip3 = true
        end
    end
    if !skip1
        result += os2s_contribute2(seed2, table, yf2 + zf2 + a1, X2, Y4, Z4, x2, y4, z4)
    end
    if !skip2
        result += os2s_contribute2(seed2, table, xf2 + zf2 + a1, X5, Y2, Z4, x4, y2, z4)
    end
    if !skip3
        result += os2s_contribute2(seed2, table, xf2 + yf2 + a1, X4, Y4, Z2, x4, y4, z2)
    end
    result * state.scale_factor
end

# 4D

struct OS2S_Vertex4D
    XYZW::NTuple{4,UInt64}
    xyzw::NTuple{4,Float64}
    function OS2S_Vertex4D(x::T, y::T, z::T, w::T) where {T<:Int}
        s = (x + y + z + w) * OS2S_UNSKEW_4D
        XYZW = (x, y, z, w) .* (PRIME_X, PRIME_Y, PRIME_Z, PRIME_W)
        xyzw = .-((x, y, z, w)) .- s
        new(XYZW, xyzw)
    end
end

const OS2S_SKEW_4D = 0.309016994374947
const OS2S_UNSKEW_4D = -0.138196601125011
const OS2S_VERTICES_4D = [
    [0x15 0x45 0x51 0x54 0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa],
    [0x15 0x45 0x51 0x55 0x56 0x59 0x5a 0x65 0x66 0x6a 0x95 0x96 0x9a 0xa6 0xaa],
    [0x01 0x05 0x11 0x15 0x41 0x45 0x51 0x55 0x56 0x5a 0x66 0x6a 0x96 0x9a 0xa6 0xaa],
    [0x01 0x15 0x16 0x45 0x46 0x51 0x52 0x55 0x56 0x5a 0x66 0x6a 0x96 0x9a 0xa6 0xaa 0xab],
    [0x15 0x45 0x54 0x55 0x56 0x59 0x5a 0x65 0x69 0x6a 0x95 0x99 0x9a 0xa9 0xaa],
    [0x05 0x15 0x45 0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0x95 0x96 0x99 0x9a 0xaa],
    [0x05 0x15 0x45 0x55 0x56 0x59 0x5a 0x66 0x6a 0x96 0x9a 0xaa],
    [0x05 0x15 0x16 0x45 0x46 0x55 0x56 0x59 0x5a 0x66 0x6a 0x96 0x9a 0xaa 0xab],
    [0x04 0x05 0x14 0x15 0x44 0x45 0x54 0x55 0x59 0x5a 0x69 0x6a 0x99 0x9a 0xa9 0xaa],
    [0x05 0x15 0x45 0x55 0x56 0x59 0x5a 0x69 0x6a 0x99 0x9a 0xaa],
    [0x05 0x15 0x45 0x55 0x56 0x59 0x5a 0x6a 0x9a 0xaa],
    [0x05 0x15 0x16 0x45 0x46 0x55 0x56 0x59 0x5a 0x5b 0x6a 0x9a 0xaa 0xab],
    [0x04 0x15 0x19 0x45 0x49 0x54 0x55 0x58 0x59 0x5a 0x69 0x6a 0x99 0x9a 0xa9 0xaa 0xae],
    [0x05 0x15 0x19 0x45 0x49 0x55 0x56 0x59 0x5a 0x69 0x6a 0x99 0x9a 0xaa 0xae],
    [0x05 0x15 0x19 0x45 0x49 0x55 0x56 0x59 0x5a 0x5e 0x6a 0x9a 0xaa 0xae],
    [0x05 0x15 0x1a 0x45 0x4a 0x55 0x56 0x59 0x5a 0x5b 0x5e 0x6a 0x9a 0xaa 0xab 0xae 0xaf],
    [0x15 0x51 0x54 0x55 0x56 0x59 0x65 0x66 0x69 0x6a 0x95 0xa5 0xa6 0xa9 0xaa],
    [0x11 0x15 0x51 0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0x95 0x96 0xa5 0xa6 0xaa],
    [0x11 0x15 0x51 0x55 0x56 0x5a 0x65 0x66 0x6a 0x96 0xa6 0xaa],
    [0x11 0x15 0x16 0x51 0x52 0x55 0x56 0x5a 0x65 0x66 0x6a 0x96 0xa6 0xaa 0xab],
    [0x14 0x15 0x54 0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0x95 0x99 0xa5 0xa9 0xaa],
    [0x15 0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0x95 0x9a 0xa6 0xa9 0xaa],
    [0x15 0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0x96 0x9a 0xa6 0xaa 0xab],
    [0x15 0x16 0x55 0x56 0x5a 0x66 0x6a 0x6b 0x96 0x9a 0xa6 0xaa 0xab],
    [0x14 0x15 0x54 0x55 0x59 0x5a 0x65 0x69 0x6a 0x99 0xa9 0xaa],
    [0x15 0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0x99 0x9a 0xa9 0xaa 0xae],
    [0x15 0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0x9a 0xaa],
    [0x15 0x16 0x55 0x56 0x59 0x5a 0x66 0x6a 0x6b 0x9a 0xaa 0xab],
    [0x14 0x15 0x19 0x54 0x55 0x58 0x59 0x5a 0x65 0x69 0x6a 0x99 0xa9 0xaa 0xae],
    [0x15 0x19 0x55 0x59 0x5a 0x69 0x6a 0x6e 0x99 0x9a 0xa9 0xaa 0xae],
    [0x15 0x19 0x55 0x56 0x59 0x5a 0x69 0x6a 0x6e 0x9a 0xaa 0xae],
    [0x15 0x1a 0x55 0x56 0x59 0x5a 0x6a 0x6b 0x6e 0x9a 0xaa 0xab 0xae 0xaf],
    [0x10 0x11 0x14 0x15 0x50 0x51 0x54 0x55 0x65 0x66 0x69 0x6a 0xa5 0xa6 0xa9 0xaa],
    [0x11 0x15 0x51 0x55 0x56 0x65 0x66 0x69 0x6a 0xa5 0xa6 0xaa],
    [0x11 0x15 0x51 0x55 0x56 0x65 0x66 0x6a 0xa6 0xaa],
    [0x11 0x15 0x16 0x51 0x52 0x55 0x56 0x65 0x66 0x67 0x6a 0xa6 0xaa 0xab],
    [0x14 0x15 0x54 0x55 0x59 0x65 0x66 0x69 0x6a 0xa5 0xa9 0xaa],
    [0x15 0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0xa5 0xa6 0xa9 0xaa 0xba],
    [0x15 0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0xa6 0xaa],
    [0x15 0x16 0x55 0x56 0x5a 0x65 0x66 0x6a 0x6b 0xa6 0xaa 0xab],
    [0x14 0x15 0x54 0x55 0x59 0x65 0x69 0x6a 0xa9 0xaa],
    [0x15 0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0xa9 0xaa],
    [0x15 0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0xaa],
    [0x15 0x16 0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0x6b 0xaa 0xab],
    [0x14 0x15 0x19 0x54 0x55 0x58 0x59 0x65 0x69 0x6a 0x6d 0xa9 0xaa 0xae],
    [0x15 0x19 0x55 0x59 0x5a 0x65 0x69 0x6a 0x6e 0xa9 0xaa 0xae],
    [0x15 0x19 0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0x6e 0xaa 0xae],
    [0x15 0x55 0x56 0x59 0x5a 0x66 0x69 0x6a 0x6b 0x6e 0x9a 0xaa 0xab 0xae 0xaf],
    [0x10 0x15 0x25 0x51 0x54 0x55 0x61 0x64 0x65 0x66 0x69 0x6a 0xa5 0xa6 0xa9 0xaa 0xba],
    [0x11 0x15 0x25 0x51 0x55 0x56 0x61 0x65 0x66 0x69 0x6a 0xa5 0xa6 0xaa 0xba],
    [0x11 0x15 0x25 0x51 0x55 0x56 0x61 0x65 0x66 0x6a 0x76 0xa6 0xaa 0xba],
    [0x11 0x15 0x26 0x51 0x55 0x56 0x62 0x65 0x66 0x67 0x6a 0x76 0xa6 0xaa 0xab 0xba 0xbb],
    [0x14 0x15 0x25 0x54 0x55 0x59 0x64 0x65 0x66 0x69 0x6a 0xa5 0xa9 0xaa 0xba],
    [0x15 0x25 0x55 0x65 0x66 0x69 0x6a 0x7a 0xa5 0xa6 0xa9 0xaa 0xba],
    [0x15 0x25 0x55 0x56 0x65 0x66 0x69 0x6a 0x7a 0xa6 0xaa 0xba],
    [0x15 0x26 0x55 0x56 0x65 0x66 0x6a 0x6b 0x7a 0xa6 0xaa 0xab 0xba 0xbb],
    [0x14 0x15 0x25 0x54 0x55 0x59 0x64 0x65 0x69 0x6a 0x79 0xa9 0xaa 0xba],
    [0x15 0x25 0x55 0x59 0x65 0x66 0x69 0x6a 0x7a 0xa9 0xaa 0xba],
    [0x15 0x25 0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0x7a 0xaa 0xba],
    [0x15 0x55 0x56 0x5a 0x65 0x66 0x69 0x6a 0x6b 0x7a 0xa6 0xaa 0xab 0xba 0xbb],
    [0x14 0x15 0x29 0x54 0x55 0x59 0x65 0x68 0x69 0x6a 0x6d 0x79 0xa9 0xaa 0xae 0xba 0xbe],
    [0x15 0x29 0x55 0x59 0x65 0x69 0x6a 0x6e 0x7a 0xa9 0xaa 0xae 0xba 0xbe],
    [0x15 0x55 0x59 0x5a 0x65 0x66 0x69 0x6a 0x6e 0x7a 0xa9 0xaa 0xae 0xba 0xbe],
    [0x15 0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0x6b 0x6e 0x7a 0xaa 0xab 0xae 0xba 0xbf],
    [0x45 0x51 0x54 0x55 0x56 0x59 0x65 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa],
    [0x41 0x45 0x51 0x55 0x56 0x59 0x5a 0x65 0x66 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xaa],
    [0x41 0x45 0x51 0x55 0x56 0x5a 0x66 0x95 0x96 0x9a 0xa6 0xaa],
    [0x41 0x45 0x46 0x51 0x52 0x55 0x56 0x5a 0x66 0x95 0x96 0x9a 0xa6 0xaa 0xab],
    [0x44 0x45 0x54 0x55 0x56 0x59 0x5a 0x65 0x69 0x95 0x96 0x99 0x9a 0xa5 0xa9 0xaa],
    [0x45 0x55 0x56 0x59 0x5a 0x65 0x6a 0x95 0x96 0x99 0x9a 0xa6 0xa9 0xaa],
    [0x45 0x55 0x56 0x59 0x5a 0x66 0x6a 0x95 0x96 0x99 0x9a 0xa6 0xaa 0xab],
    [0x45 0x46 0x55 0x56 0x5a 0x66 0x6a 0x96 0x9a 0x9b 0xa6 0xaa 0xab],
    [0x44 0x45 0x54 0x55 0x59 0x5a 0x69 0x95 0x99 0x9a 0xa9 0xaa],
    [0x45 0x55 0x56 0x59 0x5a 0x69 0x6a 0x95 0x96 0x99 0x9a 0xa9 0xaa 0xae],
    [0x45 0x55 0x56 0x59 0x5a 0x6a 0x95 0x96 0x99 0x9a 0xaa],
    [0x45 0x46 0x55 0x56 0x59 0x5a 0x6a 0x96 0x9a 0x9b 0xaa 0xab],
    [0x44 0x45 0x49 0x54 0x55 0x58 0x59 0x5a 0x69 0x95 0x99 0x9a 0xa9 0xaa 0xae],
    [0x45 0x49 0x55 0x59 0x5a 0x69 0x6a 0x99 0x9a 0x9e 0xa9 0xaa 0xae],
    [0x45 0x49 0x55 0x56 0x59 0x5a 0x6a 0x99 0x9a 0x9e 0xaa 0xae],
    [0x45 0x4a 0x55 0x56 0x59 0x5a 0x6a 0x9a 0x9b 0x9e 0xaa 0xab 0xae 0xaf],
    [0x50 0x51 0x54 0x55 0x56 0x59 0x65 0x66 0x69 0x95 0x96 0x99 0xa5 0xa6 0xa9 0xaa],
    [0x51 0x55 0x56 0x59 0x65 0x66 0x6a 0x95 0x96 0x9a 0xa5 0xa6 0xa9 0xaa],
    [0x51 0x55 0x56 0x5a 0x65 0x66 0x6a 0x95 0x96 0x9a 0xa5 0xa6 0xaa 0xab],
    [0x51 0x52 0x55 0x56 0x5a 0x66 0x6a 0x96 0x9a 0xa6 0xa7 0xaa 0xab],
    [0x54 0x55 0x56 0x59 0x65 0x69 0x6a 0x95 0x99 0x9a 0xa5 0xa6 0xa9 0xaa],
    [0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa],
    [0x15 0x45 0x51 0x55 0x56 0x59 0x5a 0x65 0x66 0x6a 0x95 0x96 0x9a 0xa6 0xaa 0xab],
    [0x55 0x56 0x5a 0x66 0x6a 0x96 0x9a 0xa6 0xaa 0xab],
    [0x54 0x55 0x59 0x5a 0x65 0x69 0x6a 0x95 0x99 0x9a 0xa5 0xa9 0xaa 0xae],
    [0x15 0x45 0x54 0x55 0x56 0x59 0x5a 0x65 0x69 0x6a 0x95 0x99 0x9a 0xa9 0xaa 0xae],
    [0x15 0x45 0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0x95 0x96 0x99 0x9a 0xa6 0xa9 0xaa 0xab 0xae],
    [0x55 0x56 0x59 0x5a 0x66 0x6a 0x96 0x9a 0xa6 0xaa 0xab],
    [0x54 0x55 0x58 0x59 0x5a 0x69 0x6a 0x99 0x9a 0xa9 0xaa 0xad 0xae],
    [0x55 0x59 0x5a 0x69 0x6a 0x99 0x9a 0xa9 0xaa 0xae],
    [0x55 0x56 0x59 0x5a 0x69 0x6a 0x99 0x9a 0xa9 0xaa 0xae],
    [0x55 0x56 0x59 0x5a 0x6a 0x9a 0xaa 0xab 0xae 0xaf],
    [0x50 0x51 0x54 0x55 0x65 0x66 0x69 0x95 0xa5 0xa6 0xa9 0xaa],
    [0x51 0x55 0x56 0x65 0x66 0x69 0x6a 0x95 0x96 0xa5 0xa6 0xa9 0xaa 0xba],
    [0x51 0x55 0x56 0x65 0x66 0x6a 0x95 0x96 0xa5 0xa6 0xaa],
    [0x51 0x52 0x55 0x56 0x65 0x66 0x6a 0x96 0xa6 0xa7 0xaa 0xab],
    [0x54 0x55 0x59 0x65 0x66 0x69 0x6a 0x95 0x99 0xa5 0xa6 0xa9 0xaa 0xba],
    [0x15 0x51 0x54 0x55 0x56 0x59 0x65 0x66 0x69 0x6a 0x95 0xa5 0xa6 0xa9 0xaa 0xba],
    [0x15 0x51 0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0x95 0x96 0x9a 0xa5 0xa6 0xa9 0xaa 0xab 0xba],
    [0x55 0x56 0x5a 0x65 0x66 0x6a 0x96 0x9a 0xa6 0xaa 0xab],
    [0x54 0x55 0x59 0x65 0x69 0x6a 0x95 0x99 0xa5 0xa9 0xaa],
    [0x15 0x54 0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0x95 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xae 0xba],
    [0x15 0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0x9a 0xa6 0xa9 0xaa],
    [0x15 0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0x96 0x9a 0xa6 0xaa 0xab],
    [0x54 0x55 0x58 0x59 0x65 0x69 0x6a 0x99 0xa9 0xaa 0xad 0xae],
    [0x55 0x59 0x5a 0x65 0x69 0x6a 0x99 0x9a 0xa9 0xaa 0xae],
    [0x15 0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0x99 0x9a 0xa9 0xaa 0xae],
    [0x15 0x55 0x56 0x59 0x5a 0x66 0x69 0x6a 0x9a 0xaa 0xab 0xae 0xaf],
    [0x50 0x51 0x54 0x55 0x61 0x64 0x65 0x66 0x69 0x95 0xa5 0xa6 0xa9 0xaa 0xba],
    [0x51 0x55 0x61 0x65 0x66 0x69 0x6a 0xa5 0xa6 0xa9 0xaa 0xb6 0xba],
    [0x51 0x55 0x56 0x61 0x65 0x66 0x6a 0xa5 0xa6 0xaa 0xb6 0xba],
    [0x51 0x55 0x56 0x62 0x65 0x66 0x6a 0xa6 0xa7 0xaa 0xab 0xb6 0xba 0xbb],
    [0x54 0x55 0x64 0x65 0x66 0x69 0x6a 0xa5 0xa6 0xa9 0xaa 0xb9 0xba],
    [0x55 0x65 0x66 0x69 0x6a 0xa5 0xa6 0xa9 0xaa 0xba],
    [0x55 0x56 0x65 0x66 0x69 0x6a 0xa5 0xa6 0xa9 0xaa 0xba],
    [0x55 0x56 0x65 0x66 0x6a 0xa6 0xaa 0xab 0xba 0xbb],
    [0x54 0x55 0x59 0x64 0x65 0x69 0x6a 0xa5 0xa9 0xaa 0xb9 0xba],
    [0x55 0x59 0x65 0x66 0x69 0x6a 0xa5 0xa6 0xa9 0xaa 0xba],
    [0x15 0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0xa5 0xa6 0xa9 0xaa 0xba],
    [0x15 0x55 0x56 0x5a 0x65 0x66 0x69 0x6a 0xa6 0xaa 0xab 0xba 0xbb],
    [0x54 0x55 0x59 0x65 0x68 0x69 0x6a 0xa9 0xaa 0xad 0xae 0xb9 0xba 0xbe],
    [0x55 0x59 0x65 0x69 0x6a 0xa9 0xaa 0xae 0xba 0xbe],
    [0x15 0x55 0x59 0x5a 0x65 0x66 0x69 0x6a 0xa9 0xaa 0xae 0xba 0xbe],
    [0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0xaa 0xab 0xae 0xba 0xbf],
    [0x40 0x41 0x44 0x45 0x50 0x51 0x54 0x55 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa],
    [0x41 0x45 0x51 0x55 0x56 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xaa],
    [0x41 0x45 0x51 0x55 0x56 0x95 0x96 0x9a 0xa6 0xaa],
    [0x41 0x45 0x46 0x51 0x52 0x55 0x56 0x95 0x96 0x97 0x9a 0xa6 0xaa 0xab],
    [0x44 0x45 0x54 0x55 0x59 0x95 0x96 0x99 0x9a 0xa5 0xa9 0xaa],
    [0x45 0x55 0x56 0x59 0x5a 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xea],
    [0x45 0x55 0x56 0x59 0x5a 0x95 0x96 0x99 0x9a 0xa6 0xaa],
    [0x45 0x46 0x55 0x56 0x5a 0x95 0x96 0x9a 0x9b 0xa6 0xaa 0xab],
    [0x44 0x45 0x54 0x55 0x59 0x95 0x99 0x9a 0xa9 0xaa],
    [0x45 0x55 0x56 0x59 0x5a 0x95 0x96 0x99 0x9a 0xa9 0xaa],
    [0x45 0x55 0x56 0x59 0x5a 0x95 0x96 0x99 0x9a 0xaa],
    [0x45 0x46 0x55 0x56 0x59 0x5a 0x95 0x96 0x99 0x9a 0x9b 0xaa 0xab],
    [0x44 0x45 0x49 0x54 0x55 0x58 0x59 0x95 0x99 0x9a 0x9d 0xa9 0xaa 0xae],
    [0x45 0x49 0x55 0x59 0x5a 0x95 0x99 0x9a 0x9e 0xa9 0xaa 0xae],
    [0x45 0x49 0x55 0x56 0x59 0x5a 0x95 0x96 0x99 0x9a 0x9e 0xaa 0xae],
    [0x45 0x55 0x56 0x59 0x5a 0x6a 0x96 0x99 0x9a 0x9b 0x9e 0xaa 0xab 0xae 0xaf],
    [0x50 0x51 0x54 0x55 0x65 0x95 0x96 0x99 0xa5 0xa6 0xa9 0xaa],
    [0x51 0x55 0x56 0x65 0x66 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xea],
    [0x51 0x55 0x56 0x65 0x66 0x95 0x96 0x9a 0xa5 0xa6 0xaa],
    [0x51 0x52 0x55 0x56 0x66 0x95 0x96 0x9a 0xa6 0xa7 0xaa 0xab],
    [0x54 0x55 0x59 0x65 0x69 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xea],
    [0x45 0x51 0x54 0x55 0x56 0x59 0x65 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xea],
    [0x45 0x51 0x55 0x56 0x59 0x5a 0x65 0x66 0x6a 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xab 0xea],
    [0x55 0x56 0x5a 0x66 0x6a 0x95 0x96 0x9a 0xa6 0xaa 0xab],
    [0x54 0x55 0x59 0x65 0x69 0x95 0x99 0x9a 0xa5 0xa9 0xaa],
    [0x45 0x54 0x55 0x56 0x59 0x5a 0x65 0x69 0x6a 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xae 0xea],
    [0x45 0x55 0x56 0x59 0x5a 0x6a 0x95 0x96 0x99 0x9a 0xa6 0xa9 0xaa],
    [0x45 0x55 0x56 0x59 0x5a 0x66 0x6a 0x95 0x96 0x99 0x9a 0xa6 0xaa 0xab],
    [0x54 0x55 0x58 0x59 0x69 0x95 0x99 0x9a 0xa9 0xaa 0xad 0xae],
    [0x55 0x59 0x5a 0x69 0x6a 0x95 0x99 0x9a 0xa9 0xaa 0xae],
    [0x45 0x55 0x56 0x59 0x5a 0x69 0x6a 0x95 0x96 0x99 0x9a 0xa9 0xaa 0xae],
    [0x45 0x55 0x56 0x59 0x5a 0x6a 0x96 0x99 0x9a 0xaa 0xab 0xae 0xaf],
    [0x50 0x51 0x54 0x55 0x65 0x95 0xa5 0xa6 0xa9 0xaa],
    [0x51 0x55 0x56 0x65 0x66 0x95 0x96 0xa5 0xa6 0xa9 0xaa],
    [0x51 0x55 0x56 0x65 0x66 0x95 0x96 0xa5 0xa6 0xaa],
    [0x51 0x52 0x55 0x56 0x65 0x66 0x95 0x96 0xa5 0xa6 0xa7 0xaa 0xab],
    [0x54 0x55 0x59 0x65 0x69 0x95 0x99 0xa5 0xa6 0xa9 0xaa],
    [0x51 0x54 0x55 0x56 0x59 0x65 0x66 0x69 0x6a 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xba 0xea],
    [0x51 0x55 0x56 0x65 0x66 0x6a 0x95 0x96 0x9a 0xa5 0xa6 0xa9 0xaa],
    [0x51 0x55 0x56 0x5a 0x65 0x66 0x6a 0x95 0x96 0x9a 0xa5 0xa6 0xaa 0xab],
    [0x54 0x55 0x59 0x65 0x69 0x95 0x99 0xa5 0xa9 0xaa],
    [0x54 0x55 0x59 0x65 0x69 0x6a 0x95 0x99 0x9a 0xa5 0xa6 0xa9 0xaa],
    [0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa],
    [0x55 0x56 0x59 0x5a 0x65 0x66 0x6a 0x95 0x96 0x9a 0xa6 0xa9 0xaa 0xab],
    [0x54 0x55 0x58 0x59 0x65 0x69 0x95 0x99 0xa5 0xa9 0xaa 0xad 0xae],
    [0x54 0x55 0x59 0x5a 0x65 0x69 0x6a 0x95 0x99 0x9a 0xa5 0xa9 0xaa 0xae],
    [0x55 0x56 0x59 0x5a 0x65 0x69 0x6a 0x95 0x99 0x9a 0xa6 0xa9 0xaa 0xae],
    [0x55 0x56 0x59 0x5a 0x66 0x69 0x6a 0x96 0x99 0x9a 0xa6 0xa9 0xaa 0xab 0xae 0xaf],
    [0x50 0x51 0x54 0x55 0x61 0x64 0x65 0x95 0xa5 0xa6 0xa9 0xaa 0xb5 0xba],
    [0x51 0x55 0x61 0x65 0x66 0x95 0xa5 0xa6 0xa9 0xaa 0xb6 0xba],
    [0x51 0x55 0x56 0x61 0x65 0x66 0x95 0x96 0xa5 0xa6 0xaa 0xb6 0xba],
    [0x51 0x55 0x56 0x65 0x66 0x6a 0x96 0xa5 0xa6 0xa7 0xaa 0xab 0xb6 0xba 0xbb],
    [0x54 0x55 0x64 0x65 0x69 0x95 0xa5 0xa6 0xa9 0xaa 0xb9 0xba],
    [0x55 0x65 0x66 0x69 0x6a 0x95 0xa5 0xa6 0xa9 0xaa 0xba],
    [0x51 0x55 0x56 0x65 0x66 0x69 0x6a 0x95 0x96 0xa5 0xa6 0xa9 0xaa 0xba],
    [0x51 0x55 0x56 0x65 0x66 0x6a 0x96 0xa5 0xa6 0xaa 0xab 0xba 0xbb],
    [0x54 0x55 0x59 0x64 0x65 0x69 0x95 0x99 0xa5 0xa9 0xaa 0xb9 0xba],
    [0x54 0x55 0x59 0x65 0x66 0x69 0x6a 0x95 0x99 0xa5 0xa6 0xa9 0xaa 0xba],
    [0x55 0x56 0x59 0x65 0x66 0x69 0x6a 0x95 0x9a 0xa5 0xa6 0xa9 0xaa 0xba],
    [0x55 0x56 0x5a 0x65 0x66 0x69 0x6a 0x96 0x9a 0xa5 0xa6 0xa9 0xaa 0xab 0xba 0xbb],
    [0x54 0x55 0x59 0x65 0x69 0x6a 0x99 0xa5 0xa9 0xaa 0xad 0xae 0xb9 0xba 0xbe],
    [0x54 0x55 0x59 0x65 0x69 0x6a 0x99 0xa5 0xa9 0xaa 0xae 0xba 0xbe],
    [0x55 0x59 0x5a 0x65 0x66 0x69 0x6a 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xae 0xba 0xbe],
    [0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0x9a 0xa6 0xa9 0xaa 0xab 0xae 0xba],
    [0x40 0x45 0x51 0x54 0x55 0x85 0x91 0x94 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xea],
    [0x41 0x45 0x51 0x55 0x56 0x85 0x91 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xaa 0xea],
    [0x41 0x45 0x51 0x55 0x56 0x85 0x91 0x95 0x96 0x9a 0xa6 0xaa 0xd6 0xea],
    [0x41 0x45 0x51 0x55 0x56 0x86 0x92 0x95 0x96 0x97 0x9a 0xa6 0xaa 0xab 0xd6 0xea 0xeb],
    [0x44 0x45 0x54 0x55 0x59 0x85 0x94 0x95 0x96 0x99 0x9a 0xa5 0xa9 0xaa 0xea],
    [0x45 0x55 0x85 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xda 0xea],
    [0x45 0x55 0x56 0x85 0x95 0x96 0x99 0x9a 0xa6 0xaa 0xda 0xea],
    [0x45 0x55 0x56 0x86 0x95 0x96 0x9a 0x9b 0xa6 0xaa 0xab 0xda 0xea 0xeb],
    [0x44 0x45 0x54 0x55 0x59 0x85 0x94 0x95 0x99 0x9a 0xa9 0xaa 0xd9 0xea],
    [0x45 0x55 0x59 0x85 0x95 0x96 0x99 0x9a 0xa9 0xaa 0xda 0xea],
    [0x45 0x55 0x56 0x59 0x5a 0x85 0x95 0x96 0x99 0x9a 0xaa 0xda 0xea],
    [0x45 0x55 0x56 0x5a 0x95 0x96 0x99 0x9a 0x9b 0xa6 0xaa 0xab 0xda 0xea 0xeb],
    [0x44 0x45 0x54 0x55 0x59 0x89 0x95 0x98 0x99 0x9a 0x9d 0xa9 0xaa 0xae 0xd9 0xea 0xee],
    [0x45 0x55 0x59 0x89 0x95 0x99 0x9a 0x9e 0xa9 0xaa 0xae 0xda 0xea 0xee],
    [0x45 0x55 0x59 0x5a 0x95 0x96 0x99 0x9a 0x9e 0xa9 0xaa 0xae 0xda 0xea 0xee],
    [0x45 0x55 0x56 0x59 0x5a 0x95 0x96 0x99 0x9a 0x9b 0x9e 0xaa 0xab 0xae 0xda 0xea 0xef],
    [0x50 0x51 0x54 0x55 0x65 0x91 0x94 0x95 0x96 0x99 0xa5 0xa6 0xa9 0xaa 0xea],
    [0x51 0x55 0x91 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xe6 0xea],
    [0x51 0x55 0x56 0x91 0x95 0x96 0x9a 0xa5 0xa6 0xaa 0xe6 0xea],
    [0x51 0x55 0x56 0x92 0x95 0x96 0x9a 0xa6 0xa7 0xaa 0xab 0xe6 0xea 0xeb],
    [0x54 0x55 0x94 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xe9 0xea],
    [0x55 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xea],
    [0x55 0x56 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xea],
    [0x55 0x56 0x95 0x96 0x9a 0xa6 0xaa 0xab 0xea 0xeb],
    [0x54 0x55 0x59 0x94 0x95 0x99 0x9a 0xa5 0xa9 0xaa 0xe9 0xea],
    [0x55 0x59 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xea],
    [0x45 0x55 0x56 0x59 0x5a 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xea],
    [0x45 0x55 0x56 0x5a 0x95 0x96 0x99 0x9a 0xa6 0xaa 0xab 0xea 0xeb],
    [0x54 0x55 0x59 0x95 0x98 0x99 0x9a 0xa9 0xaa 0xad 0xae 0xe9 0xea 0xee],
    [0x55 0x59 0x95 0x99 0x9a 0xa9 0xaa 0xae 0xea 0xee],
    [0x45 0x55 0x59 0x5a 0x95 0x96 0x99 0x9a 0xa9 0xaa 0xae 0xea 0xee],
    [0x55 0x56 0x59 0x5a 0x95 0x96 0x99 0x9a 0xaa 0xab 0xae 0xea 0xef],
    [0x50 0x51 0x54 0x55 0x65 0x91 0x94 0x95 0xa5 0xa6 0xa9 0xaa 0xe5 0xea],
    [0x51 0x55 0x65 0x91 0x95 0x96 0xa5 0xa6 0xa9 0xaa 0xe6 0xea],
    [0x51 0x55 0x56 0x65 0x66 0x91 0x95 0x96 0xa5 0xa6 0xaa 0xe6 0xea],
    [0x51 0x55 0x56 0x66 0x95 0x96 0x9a 0xa5 0xa6 0xa7 0xaa 0xab 0xe6 0xea 0xeb],
    [0x54 0x55 0x65 0x94 0x95 0x99 0xa5 0xa6 0xa9 0xaa 0xe9 0xea],
    [0x55 0x65 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xea],
    [0x51 0x55 0x56 0x65 0x66 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xea],
    [0x51 0x55 0x56 0x66 0x95 0x96 0x9a 0xa5 0xa6 0xaa 0xab 0xea 0xeb],
    [0x54 0x55 0x59 0x65 0x69 0x94 0x95 0x99 0xa5 0xa9 0xaa 0xe9 0xea],
    [0x54 0x55 0x59 0x65 0x69 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xea],
    [0x55 0x56 0x59 0x65 0x6a 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xea],
    [0x55 0x56 0x5a 0x66 0x6a 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xab 0xea 0xeb],
    [0x54 0x55 0x59 0x69 0x95 0x99 0x9a 0xa5 0xa9 0xaa 0xad 0xae 0xe9 0xea 0xee],
    [0x54 0x55 0x59 0x69 0x95 0x99 0x9a 0xa5 0xa9 0xaa 0xae 0xea 0xee],
    [0x55 0x59 0x5a 0x69 0x6a 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xae 0xea 0xee],
    [0x55 0x56 0x59 0x5a 0x6a 0x95 0x96 0x99 0x9a 0xa6 0xa9 0xaa 0xab 0xae 0xea],
    [0x50 0x51 0x54 0x55 0x65 0x95 0xa1 0xa4 0xa5 0xa6 0xa9 0xaa 0xb5 0xba 0xe5 0xea 0xfa],
    [0x51 0x55 0x65 0x95 0xa1 0xa5 0xa6 0xa9 0xaa 0xb6 0xba 0xe6 0xea 0xfa],
    [0x51 0x55 0x65 0x66 0x95 0x96 0xa5 0xa6 0xa9 0xaa 0xb6 0xba 0xe6 0xea 0xfa],
    [0x51 0x55 0x56 0x65 0x66 0x95 0x96 0xa5 0xa6 0xa7 0xaa 0xab 0xb6 0xba 0xe6 0xea 0xfb],
    [0x54 0x55 0x65 0x95 0xa4 0xa5 0xa6 0xa9 0xaa 0xb9 0xba 0xe9 0xea 0xfa],
    [0x55 0x65 0x95 0xa5 0xa6 0xa9 0xaa 0xba 0xea 0xfa],
    [0x51 0x55 0x65 0x66 0x95 0x96 0xa5 0xa6 0xa9 0xaa 0xba 0xea 0xfa],
    [0x55 0x56 0x65 0x66 0x95 0x96 0xa5 0xa6 0xaa 0xab 0xba 0xea 0xfb],
    [0x54 0x55 0x65 0x69 0x95 0x99 0xa5 0xa6 0xa9 0xaa 0xb9 0xba 0xe9 0xea 0xfa],
    [0x54 0x55 0x65 0x69 0x95 0x99 0xa5 0xa6 0xa9 0xaa 0xba 0xea 0xfa],
    [0x55 0x65 0x66 0x69 0x6a 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xba 0xea 0xfa],
    [0x55 0x56 0x65 0x66 0x6a 0x95 0x96 0x9a 0xa5 0xa6 0xa9 0xaa 0xab 0xba 0xea],
    [0x54 0x55 0x59 0x65 0x69 0x95 0x99 0xa5 0xa9 0xaa 0xad 0xae 0xb9 0xba 0xe9 0xea 0xfe],
    [0x55 0x59 0x65 0x69 0x95 0x99 0xa5 0xa9 0xaa 0xae 0xba 0xea 0xfe],
    [0x55 0x59 0x65 0x69 0x6a 0x95 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xae 0xba 0xea],
    [0x55 0x56 0x59 0x5a 0x65 0x66 0x69 0x6a 0x95 0x96 0x99 0x9a 0xa5 0xa6 0xa9 0xaa 0xab 0xae 0xba 0xea],
]

const OS2S_LATTICE_VERTICES_4D = let
    vertices = Vector{OS2S_Vertex4D}(undef, 256)
    for i in 0:255
        cx = ((i >> 0) & 3) - 1
        cy = ((i >> 2) & 3) - 1
        cz = ((i >> 4) & 3) - 1
        cw = ((i >> 6) & 3) - 1
        vertices[i+1] = OS2S_Vertex4D(cx, cy, cz, cw)
    end
    vertices
end

const OS2S_LOOKUP_4D_A, OS2S_LOOKUP_4D_B = let
    a = Vector{NTuple{2,Int16}}(undef, 256)
    b = Vector{OS2S_Vertex4D}(undef, 3476)
    j = 1
    for i in 1:256
        len = length(OS2S_VERTICES_4D[i])
        a[i] = (j, j + len - 1)
        for k in 1:len
            b[j] = OS2S_LATTICE_VERTICES_4D[OS2S_VERTICES_4D[i][k]+1]
            j += 1
        end
    end
    (a, b)
end

"""
    opensimplex2s_4d(; seed=nothing, smooth=false, orient=nothing)

Construct a sampler that outputs 4-dimensional OpenSimplex2S noise when it is sampled from.

# Arguments

  - `seed`: An unsigned integer used to seed the random number generator for this sampler, or
  `nothing` for non-deterministic results.

  - `smooth`: Specify whether to have continuous gradients.
    Simplex variants, even the original Simplex noise by Ken Perlin, overshoot the radial extent for
    the signal reconstruction kernel in order to improve the visual of the noise. Normally this is
    okay, especially if layering multiple octaves of the noise. However, in some applications, such
    as creating height or bump maps, this will produce discontinuities visually identified by
    jarring creases in the generated noise.

    This option changes the falloff in order to produce smooth continuous noise, however, the
    resulting noise may look quite different than the non-smooth option, depending on the Simplex
    variant.

    The default value is `false`, in order to be true to the original implementation.

  - `orient`: Either the symbol `:x` or the value `nothing`:

      + `nothing`: Use the standard orientation.
      + `:x`: The noise space will be re-oriented with the Y axis pointing down the main diagonal to
        improve visual isotropy.
      + `:xy`: Re-orient the noise space to have better visual isotropy in the XY plane.
      + `:xz`: Re-orient the noise space to have better visual isotropy in the XZ plane.
      + `:xyz`: Re-orient the noise space to be better suited for time-varied animations, where
        the W axis is time.
"""
function opensimplex2s_4d(; seed=nothing, orient=nothing, smooth=false)
    size = OS2_NUM_GRADIENTS_4D * 4
    gradients = OS2_GRADIENTS_NORMALIZED_4D ./ 0.11127401889945551
    _opensimplex2s(4, seed, size, gradients, orient, smooth)
end

@inline function orient(::Type{OpenSimplex2S{4,OrientStandard}}, x, y, z, w)
    (x, y, z, w) .+ OS2S_SKEW_4D .* (x + y + z + w)
end

@inline function orient(::Type{OpenSimplex2S{4,OrientXY}}, x, y, z, w)
    xy = x + y
    ww = w * 1.118033988749894
    zw = z * 0.28867513459481294226 + ww
    xr, yr = (x, y) .+ zw .+ xy .* -0.21132486540518699998
    zr = xy * -0.57735026918962599998 + zw
    wr = z * -0.866025403784439 + ww
    (xr, yr, zr, wr)
end

@inline function orient(::Type{OpenSimplex2S{4,OrientXZ}}, x, y, z, w)
    orient(OpenSimplex2S{4,OrientXY}, x, z, y, w)
end

@inline function orient(::Type{OpenSimplex2S{4,OrientXYZ}}, x, y, z, w)
    xyz = -(x + y + z)
    ww = w * 1.118033988749894
    s = xyz / 6 + ww
    xs, ys, zs = (x, y, z) .+ s
    ws = xyz * 0.5 + ww
    (xs, ys, zs, ws)
end

function sample(sampler::S, x::T, y::T, z::T, w::T) where {O,S<:OpenSimplex2S{4,O},T<:Real}
    seed = sampler.random_state.seed
    table = sampler.table
    state = sampler.simplex_state
    falloff = state.falloff
    primes = (PRIME_X, PRIME_Y, PRIME_Z, PRIME_W)
    tr = orient(S, x, y, z, w)
    XYZW = floor.(Int, tr)
    XYZWp = XYZW .* primes
    x1, y1, z1, w1 = tr .- XYZW
    vs = (x1, y1, z1, w1) .+ (x1 + y1 + z1 + w1) * OS2S_UNSKEW_4D
    ix, iy, iz, iw = floor.(Int, (tr .* 4)) .& 3
    index = ix << 0 | iy << 2 | iz << 4 | iw << 6
    result = 0.0
    @inbounds start, stop = OS2S_LOOKUP_4D_A[index+1]
    @inbounds for i in start:stop
        c = OS2S_LOOKUP_4D_B[i]
        V = XYZWp .+ c.XYZW
        x, y, z, w = vs .+ c.xyzw
        a = (x^2 + y^2) + (z^2 + w^2)
        if a < falloff
            result += pow4(a - falloff) * grad(table, seed, V..., x, y, z, w)
        end
    end
    result * state.scale_factor
end
