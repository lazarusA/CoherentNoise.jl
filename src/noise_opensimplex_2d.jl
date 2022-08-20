const OS_GRADIENTS_2D = Int8.((5, 2, 2, 5, -5, 2, -2, 5, 5, -2, 2, -5, -5, -2, -2, -5))
const OS_STRETCH_2D = (1 / sqrt(3) - 1) / 2
const OS_SQUISH_2D = (sqrt(3) - 1) / 2
const OS_SCALE_2D = 1 / 40.7

"""
    opensimplex_2d(; kwargs...)

Construct a sampler that outputs 2-dimensional legacy OpenSimplex noise when it is sampled from.

# Arguments

  - `seed=0`: An integer used to seed the random number generator for this sampler.

# See also:

  - [`OpenSimplex2`](@ref opensimplex2_2d)
  - [`OpenSimplex2S`](@ref opensimplex2s_2d)
"""
opensimplex_2d(; seed=0) = opensimplex(2, seed)

@inline function contribute(sampler::OpenSimplex{2}, X, Y, x, y)
    t = sampler.state.table
    g = OS_GRADIENTS_2D
    @inbounds i = (t[t[X]+Y] & 14) + 1
    a = 2 - x^2 - y^2
    @inbounds @fastpow a > 0 ? a^4 * (g[i] * x + g[i+1] * y) : 0.0
end

function sample(sampler::OpenSimplex{2}, x::T, y::T) where {T<:Real}
    s = (x, y) .+ ((x + y) * OS_STRETCH_2D)
    X1, Y1 = floor.(Int, s)
    x1, y1 = (x, y) .- ((X1 + Y1) * OS_SQUISH_2D) .- (X1, Y1)
    y2, x3 = (y1, x1) .- OS_SQUISH_2D
    x2, y3 = (x1, y1) .- 1 .- OS_SQUISH_2D
    Xs, Ys = s .- (X1, Y1)
    XYs = Xs + Ys
    c1 = contribute(sampler, X1 + 1, Y1, x2, y2)
    c2 = contribute(sampler, X1, Y1 + 1, x3, y3)
    result = c1 + c2
    if XYs ≤ 1
        Zs = 1 - XYs
        if Zs > Xs || Zs > Ys
            if Xs > Ys
                X2, y4 = (X1, y1) .+ 1
                Y2, x4 = (Y1, x1) .- 1
            else
                X2, y4 = (X1, y1) .- 1
                Y2, x4 = (Y1, x1) .+ 1
            end
        else
            X2, Y2 = X1 + 1, Y1 + 1
            x4, y4 = (x1, y1) .- 1 .- 2OS_SQUISH_2D
        end
    else
        Zs = 2 - XYs
        if Zs < Xs || Zs < Ys
            if Xs > Ys
                X2, Y2 = X1 + 2, Y1
                x4, y4 = (x1, y1) .- 2 .- 2OS_SQUISH_2D
            else
                X2, Y2 = X1, Y1 + 1
                x4 = x1 - 2OS_SQUISH_2D
                y4 = y1 - 2 - 2OS_SQUISH_2D
            end
        else
            X2, Y2 = X1, Y1
            x4, y4 = x1, y1
        end
        X1, Y1 = (X1, Y1) .+ 1
        x1, y1 = (x1, y1) .- 2OS_SQUISH_2D .- 1
    end
    c1 = contribute(sampler, X1, Y1, x1, y1)
    c2 = contribute(sampler, X2, Y2, x4, y4)
    (result + c1 + c2) * SCALE_2D
end
