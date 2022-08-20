const OS2_SKEW_2D = 0.366025403784439
const OS2_UNSKEW_2D = -0.21132486540518713
const OS2_R²2D = 0.5
const OS2_NUM_GRADIENTS_EXP_2D = 7
const OS2_NUM_GRADIENTS_2D = 1 << OS2_NUM_GRADIENTS_EXP_2D
const OS2_GRADIENTS_NORMALIZED_2D = [
    0.38268343236509, 0.923879532511287,
    0.923879532511287, 0.38268343236509,
    0.923879532511287, -0.38268343236509,
    0.38268343236509, -0.923879532511287,
    -0.38268343236509, -0.923879532511287,
    -0.923879532511287, -0.38268343236509,
    -0.923879532511287, 0.38268343236509,
    -0.38268343236509, 0.923879532511287,
    0.130526192220052, 0.99144486137381,
    0.608761429008721, 0.793353340291235,
    0.793353340291235, 0.608761429008721,
    0.99144486137381, 0.130526192220051,
    0.99144486137381, -0.130526192220051,
    0.793353340291235, -0.60876142900872,
    0.608761429008721, -0.793353340291235,
    0.130526192220052, -0.99144486137381,
    -0.130526192220052, -0.99144486137381,
    -0.608761429008721, -0.793353340291235,
    -0.793353340291235, -0.608761429008721,
    -0.99144486137381, -0.130526192220052,
    -0.99144486137381, 0.130526192220051,
    -0.793353340291235, 0.608761429008721,
    -0.608761429008721, 0.793353340291235,
    -0.130526192220052, 0.99144486137381]
const OS2_GRADIENTS_2D = OS2_GRADIENTS_NORMALIZED_2D ./ 0.01001634121365712 |> CircularVector

"""
    opensimplex2_2d(; kwargs...)

Construct a sampler that outputs 2-dimensional OpenSimplex2 noise when it is sampled from.

# Arguments

  - `seed=0`: An integer used to seed the random number generator for this sampler.

  - `orient=nothing`: Either the symbol `:x` or the value `nothing`:

      + `:x`: The noise space will be re-oriented with the Y axis pointing down the main diagonal to
        improve visual isotropy.

      + `nothing`: Use the standard orientation.
"""
opensimplex2_2d(; seed=0, orient=nothing) = opensimplex2(2, seed, orient)

@inline function grad(table, seed, X, Y, x, y)
    hash = (seed ⊻ X ⊻ Y) * HASH_MULTIPLIER
    hash ⊻= hash >> (64 - OS2_NUM_GRADIENTS_EXP_2D + 1)
    i = trunc(hash) & ((OS2_NUM_GRADIENTS_2D - 1) << 1)
    t = (table[i+1], table[(i|1)+1])
    sum((t .* (x, y)))
end

@inline transform(::OpenSimplex2{2,OrientStandard}, x, y) = (x, y) .+ OS2_SKEW_2D .* (x + y)

@inline function transform(::OpenSimplex2{2,OrientX}, x, y)
    xx = x * ROOT_2_OVER_2
    yy = y * ROOT_2_OVER_2 * (2OS2_SKEW_2D + 1)
    (yy + xx, yy - xx)
end

@fastpow function sample(sampler::OpenSimplex2{2}, x::T, y::T) where {T<:Real}
    seed = sampler.seed
    primes = (PRIME_X, PRIME_Y)
    tr = transform(sampler, x, y)
    XY = floor.(Int, tr)
    vtr = tr .- XY
    t = sum(vtr) * OS2_UNSKEW_2D
    X1, Y1 = XY .* primes
    X2, Y2 = (X1, Y1) .+ primes
    x1, y1 = vtr .+ t
    us1 = 2OS2_UNSKEW_2D + 1
    result = 0.0
    a1 = OS2_R²2D - x1^2 - y1^2
    if a1 > 0
        result += a1^4 * grad(OS2_GRADIENTS_2D, seed, X1, Y1, x1, y1)
    end
    a2 = 2us1 * (1 / OS2_UNSKEW_2D + 2) * t + -2us1^2 + a1
    if a2 > 0
        x, y = (x1, y1) .- 2OS2_UNSKEW_2D .- 1
        result += a2^4 * grad(OS2_GRADIENTS_2D, seed, X2, Y2, x, y)
    end
    if y1 > x1
        x = x1 - OS2_UNSKEW_2D
        y = y1 - OS2_UNSKEW_2D - 1
        a3 = OS2_R²2D - x^2 - y^2
        if a3 > 0
            result += a3^4 * grad(OS2_GRADIENTS_2D, seed, X1, Y2, x, y)
        end
    else
        x = x1 - OS2_UNSKEW_2D - 1
        y = y1 - OS2_UNSKEW_2D
        a4 = OS2_R²2D - x^2 - y^2
        if a4 > 0
            result += a4^4 * grad(OS2_GRADIENTS_2D, seed, X2, Y1, x, y)
        end
    end
    result
end
