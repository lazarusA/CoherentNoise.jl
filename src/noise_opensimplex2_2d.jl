const SKEW_2D = 0.366025403784439

const UNSKEW_2D = -0.21132486540518713

const R²2D = 0.5

const NUM_GRADIENTS_EXP_2D = 7

const NUM_GRADIENTS_2D = 1 << NUM_GRADIENTS_EXP_2D

const GRADIENTS_NORMALIZED_2D = [
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

const GRADIENTS_2D = GRADIENTS_NORMALIZED_2D ./ 0.01001634121365712 |> CircularVector

@inline function grad(table, seed, X, Y, x, y)
    hash = (seed ⊻ X ⊻ Y) * HASH_MULTIPLIER
    hash ⊻= hash >> (64 - NUM_GRADIENTS_EXP_2D + 1)
    i = trunc(hash) & ((NUM_GRADIENTS_2D - 1) << 1)
    t = (table[i+1], table[(i|1)+1])
    sum((t .* (x, y)))
end

@inline transform(::Type{Standard}, x, y) = (x, y) .+ SKEW_2D .* (x + y)

@inline function transform(::Type{ImproveX}, x, y)
    xx = x * ROOT_2_OVER_2
    yy = y * ROOT_2_OVER_2 * (2SKEW_2D + 1)
    (yy + xx, yy - xx)
end

@fastpow function sample(sampler::OpenSimplex2{2,O}, x::T, y::T) where {O,T<:Real}
    seed = get_seed(sampler)
    primes = (PRIME_X, PRIME_Y)
    tr = transform(O, x, y)
    XY = floor.(Int, tr)
    vtr = tr .- XY
    t = sum(vtr) * UNSKEW_2D
    X1, Y1 = XY .* primes
    X2, Y2 = (X1, Y1) .+ primes
    x1, y1 = vtr .+ t
    us1 = 2UNSKEW_2D + 1
    result = 0.0
    a1 = R²2D - x1^2 - y1^2
    if a1 > 0
        result += a1^4 * grad(GRADIENTS_2D, seed, X1, Y1, x1, y1)
    end
    a2 = 2us1 * (1 / UNSKEW_2D + 2) * t + -2us1^2 + a1
    if a2 > 0
        x, y = (x1, y1) .- 2UNSKEW_2D .- 1
        result += a2^4 * grad(GRADIENTS_2D, seed, X2, Y2, x, y)
    end
    if y1 > x1
        x = x1 - UNSKEW_2D
        y = y1 - UNSKEW_2D - 1
        a3 = R²2D - x^2 - y^2
        if a3 > 0
            result += a3^4 * grad(GRADIENTS_2D, seed, X1, Y2, x, y)
        end
    else
        x = x1 - UNSKEW_2D - 1
        y = y1 - UNSKEW_2D
        a4 = R²2D - x^2 - y^2
        if a4 > 0
            result += a4^4 * grad(GRADIENTS_2D, seed, X2, Y1, x, y)
        end
    end
    result
end
