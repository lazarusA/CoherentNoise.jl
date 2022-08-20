const OS2S_R²2D = 2 / 3
const OS2S_GRADIENTS_2D = OS2_GRADIENTS_NORMALIZED_2D ./ 0.05481866495625118 |> CircularVector

"""
    opensimplex2s_2d(; kwargs...)

Construct a sampler that outputs 2-dimensional OpenSimplex2S noise when it is sampled from.

# Arguments

  - `seed=0`: An integer used to seed the random number generator for this sampler.

  - `orient=nothing`: Either the symbol `:x` or the value `nothing`:

      + `:x`: The noise space will be re-oriented with the Y axis pointing down the main diagonal
        to improve visual isotropy.

      + `nothing`: Use the standard orientation.
"""
opensimplex2s_2d(; seed=0, orient=nothing) = opensimplex2s(2, seed, orient)

@inline @fastpow function contribute(seed, X, Y, x, y)
    a = OS2S_R²2D - x^2 - y^2
    a > 0 ? a^4 * grad(OS2S_GRADIENTS_2D, seed, X, Y, x, y) : 0.0
end

@inline transform(::OpenSimplex2S{2,OrientStandard}, x, y) = (x, y) .+ OS2_SKEW_2D .* (x + y)

@inline function transform(::OpenSimplex2S{2,OrientX}, x, y)
    xx = x * ROOT_2_OVER_2
    yy = y * ROOT_2_OVER_2 * (2OS2_SKEW_2D + 1)
    (yy + xx, yy - xx)
end

@fastpow function sample(sampler::OpenSimplex2S{2}, x::T, y::T) where {T<:Real}
    seed = sampler.seed
    primes = (PRIME_X, PRIME_Y)
    tr = transform(sampler, x, y)
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
    a1 = OS2S_R²2D - x1^2 - y1^2
    c1 = a1^4 * grad(OS2S_GRADIENTS_2D, seed, X1, Y1, x1, y1)
    a2 = 2us2p1 * (1 / us + 2) * t + -2us2p1^2 + a1
    c2 = a2^4 * grad(OS2S_GRADIENTS_2D, seed, X2, Y2, x2, y2)
    result = c1 + c2
    if t < OS2_UNSKEW_2D
        X3, Y3 = (X1, Y1) .+ (primes .<< 1)
        if xs + d > 1
            result += contribute(seed, X3, Y2, x1 - 3us - 2, y1 - 3us - 1)
        else
            result += contribute(seed, X1, Y2, x1 - us, y1 - us - 1)
        end
        if ys - d > 1
            result += contribute(seed, X2, Y3, x1 - 3us - 1, y1 - 3us - 2)
        else
            result += contribute(seed, X2, Y1, x1 - us - 1, y1 - us)
        end
    else
        X4, Y4 = (X1, Y1) .- primes
        if xs + d < 0
            result += contribute(seed, X4, Y1, x1 + us + 1, y1 + us)
        else
            result += contribute(seed, X2, Y1, x1 - us - 1, y1 - us)
        end
        if ys < d
            result += contribute(seed, X1, Y4, x1 + us, y1 + us + 1)
        else
            result += contribute(seed, X1, Y2, x1 - us, y1 - us - 1)
        end
    end
    result
end
