const GRADIENTS_2D = Int8.((5, 2, 2, 5, -5, 2, -2, 5, 5, -2, 2, -5, -5, -2, -2, -5))

const STRETCH_2D = (1 / sqrt(3) - 1) / 2

const SQUISH_2D = (sqrt(3) - 1) / 2

const SCALE_2D = 1 / 40.7

@inline function contribute(sampler::OpenSimplex{2}, X, Y, x, y)
    t = sampler.state.table
    g = GRADIENTS_2D
    @inbounds i = (t[t[X]+Y] & 14) + 1
    a = 2 - x^2 - y^2
    @inbounds @fastpow a > 0 ? a^4 * (g[i] * x + g[i+1] * y) : 0.0
end

function sample(sampler::OpenSimplex{2}, x::T, y::T) where {T<:Real}
    s = (x, y) .+ ((x + y) * STRETCH_2D)
    X1, Y1 = floor.(Int, s)
    x1, y1 = (x, y) .- ((X1 + Y1) * SQUISH_2D) .- (X1, Y1)
    y2, x3 = (y1, x1) .- SQUISH_2D
    x2, y3 = (x1, y1) .- 1 .- SQUISH_2D
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
            x4, y4 = (x1, y1) .- 1 .- 2SQUISH_2D
        end
    else
        Zs = 2 - XYs
        if Zs < Xs || Zs < Ys
            if Xs > Ys
                X2, Y2 = X1 + 2, Y1
                x4, y4 = (x1, y1) .- 2 .- 2SQUISH_2D
            else
                X2, Y2 = X1, Y1 + 1
                x4 = x1 - 2SQUISH_2D
                y4 = y1 - 2 - 2SQUISH_2D
            end
        else
            X2, Y2 = X1, Y1
            x4, y4 = x1, y1
        end
        X1, Y1 = (X1, Y1) .+ 1
        x1, y1 = (x1, y1) .- 2SQUISH_2D .- 1
    end
    c1 = contribute(sampler, X1, Y1, x1, y1)
    c2 = contribute(sampler, X2, Y2, x4, y4)
    (result + c1 + c2) * SCALE_2D
end
