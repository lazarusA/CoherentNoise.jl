const OS2S_R²3D = 3 / 4
const OS2S_GRADIENTS_3D = OS2_GRADIENTS_NORMALIZED_3D ./ 0.2781926117527186 |> CircularVector

"""
    opensimplex2s_3d(; kwargs...)

Construct a sampler that outputs 3-dimensional OpenSimplex2S noise when it is sampled from.

# Arguments

  - `seed=0`: An integer used to seed the random number generator for this sampler.

  - `orient=nothing`: One of the following symbols or the value `nothing`:

      + `:x`: The noise space will be re-oriented with the Y axis pointing down the main diagonal to
        improve visual isotropy.

      + `:xy`: Re-orient the noise space to have better visual isotropy in the XY plane.

      + `:xz`: Re-orient the noise space to have better visual isotropy in the XZ plane.

      + `nothing`: Use the standard orientation.
"""
opensimplex2s_3d(; seed=0, orient=nothing) = opensimplex2s(3, seed, orient)

@inline @fastpow os2s_contribute1(seed, a, args...) = a^4 * grad(OS2S_GRADIENTS_3D, seed, args...)

@inline os2s_contribute2(seed, a, args...) = a > 0 ? os2s_contribute1(seed, a, args...) : 0.0

@inline function transform(::OpenSimplex2S{3,OrientStandard}, x, y, z)
    OS2_FALLBACK_ROTATE_3D * (x + y + z) .- (x, y, z)
end

@inline function transform(::OpenSimplex2S{3,OrientXY}, x, y, z)
    xy = x + y
    zz = z * ROOT_3_OVER_3
    xr, yr = (x, y) .+ xy .* OS2_ROTATE_3D_ORTHONORMALIZER .+ zz
    zr = xy * -ROOT_3_OVER_3 + zz
    (xr, yr, zr)
end

@inline function transform(::OpenSimplex2S{3,OrientXZ}, x, y, z)
    xz = x + z
    yy = y * ROOT_3_OVER_3
    xr, zr = (x, z) .+ xz .* -0.211324865405187 .+ yy
    yr = xz * -ROOT_3_OVER_3 + yy
    (xr, yr, zr)
end

@fastpow function sample(sampler::OpenSimplex2S{3}, x::T, y::T, z::T) where {T<:Real}
    seed = sampler.seed
    seed2 = seed ⊻ -OS2_SEED_FLIP_3D
    primes = (PRIME_X, PRIME_Y, PRIME_Z)
    tr = transform(sampler, x, y, z)
    V = floor.(Int, tr)
    Vr = tr .- V
    Vp = V .* primes
    mask1 = trunc.(Int, -0.5 .- Vr)
    mask2 = mask1 .| 1
    X1, Y1, Z1 = Vp .+ mask1 .& primes
    X2, Y2, Z2 = Vp .+ primes
    X3, Y3, Z3 = Vp .+ .~mask1 .& primes
    X4, Y4, Z4 = Vp .+ mask1 .& primes .<< 1
    X5 = Vp[1] + mask1[1] & PRIME_X * 2
    x1, y1, z1 = Vr .+ mask1
    x2, y2, z2 = Vr .- 0.5
    x3, y3, z3 = (x1, y1, z1) .- mask2
    x4, y4, z4 = (x2, y2, z2) .+ mask2
    a0 = OS2S_R²3D - x1^2 - y1^2 - z1^2
    a1 = OS2S_R²3D - x2^2 - y2^2 - z2^2
    xf1, yf1, zf1 = (x2, y2, z2) .* mask2 .<< 1
    xf2, yf2, zf2 = (x2, y2, z2) .* (-2 .- mask1 .<< 2) .- 1
    xf3, yf3, zf3 = (xf1, yf1, zf1) .+ a0
    xf4, yf4, zf4 = (xf2, yf2, zf2) .+ a1
    skip1, skip2, skip3 = false, false, false
    result = os2s_contribute1(seed, a0, X1, Y1, Z1, x1, y1, z1)
    result += os2s_contribute1(seed2, a1, X2, Y2, Z2, x2, y2, z2)
    if xf3 > 0
        result += os2s_contribute1(seed, xf3, X3, Y1, Z1, x3, y1, z1)
    else
        os2s_contribute2(seed, yf1 + zf1 + a0, X1, Y3, Z3, x1, y3, z3)
        if xf4 > 0
            result += os2s_contribute1(seed2, xf4, X5, Y2, Z2, x4, y2, z2)
            skip1 = true
        end
    end
    if yf3 > 0
        result += os2s_contribute1(seed, yf3, X1, Y3, Z1, x1, y3, z1)
    else
        result += os2s_contribute2(seed, xf1 + zf1 + a0, X3, Y1, Z3, x3, y1, z3)
        if yf4 > 0
            result += os2s_contribute1(seed2, yf4, X2, Y4, Z2, x2, y4, z2)
            skip2 = true
        end
    end
    if zf3 > 0
        result += os2s_contribute1(seed, zf3, X1, Y1, Z3, x1, y1, z3)
    else
        result += os2s_contribute2(seed, xf1 + yf1 + a0, X3, Y3, Z1, x3, y3, z1)
        if zf4 > 0
            result += os2s_contribute1(seed2, zf4, X2, Y2, Z4, x2, y2, z4)
            skip3 = true
        end
    end
    if !skip1
        result += os2s_contribute2(seed2, yf2 + zf2 + a1, X2, Y4, Z4, x2, y4, z4)
    end
    if !skip2
        result += os2s_contribute2(seed2, xf2 + zf2 + a1, X5, Y2, Z4, x4, y2, z4)
    end
    if !skip3
        result += os2s_contribute2(seed2, xf2 + yf2 + a1, X4, Y4, Z2, x4, y4, z2)
    end
    result
end
