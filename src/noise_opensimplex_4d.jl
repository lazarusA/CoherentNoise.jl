const STRETCH_4D = (1 / sqrt(5) - 1) / 4

const SQUISH_4D = (sqrt(5) - 1) / 4

const SCALE_4D = 1 / 30

const GRADIENTS_4D = [
    3, 1, 1, 1, 1, 3, 1, 1, 1, 1, 3, 1, 1, 1, 1, 3, -3, 1, 1, 1, -1, 3, 1, 1, -1, 1, 3, 1, -1,
    1, 1, 3, 3, -1, 1, 1, 1, -3, 1, 1, 1, -1, 3, 1, 1, -1, 1, 3, -3, -1, 1, 1, -1, -3, 1, 1, -1,
    -1, 3, 1, -1, -1, 1, 3, 3, 1, -1, 1, 1, 3, -1, 1, 1, 1, -3, 1, 1, 1, -1, 3, -3, 1, 1, 1, -1,
    3, -1, 1, -1, 1, -3, 1, -1, 1, -1, 3, 3, -1, -1, 1, 1, -3, -1, 1, 1, -1, -3, 1, 1, -1, -1,
    3, -3, -1, -1, 1, -1, -3, -1, 1, -1, -1, -3, 1, -1, -1, -1, 3, 3, 1, 1, -1, 1, 3, 1, -1, 1,
    1, 3, -1, 1, 1, 1, -3, -3, 1, 1, -1, -1, 3, 1, -1, -1, 1, 3, -1, -1, 1, 1, -3, 3, -1, 1, -1,
    1, -3, 1, -1, 1, -1, 3, -1, 1, -1, 1, -3, -3, -1, 1, -1, -1, 3, 1, -1, -1, -1, 3, -1, -1,
    -1, 1, -3, 3, 1, -1, -1, 1, 3, -1, -1, 1, 1, -3, -1, 1, 1, -1, -3, -3, 1, -1, -1, -1, 3, -1,
    -1, -1, 1, -3, -1, -1, 1, -1, -3, 3, -1, -1, -1, 1, -3, -1, -1, 1, -1, -3, -1, 1, -1, -1,
    -3, -3, -1, -1, -1, -1, -3, -1, -1, -1, -1, 3, -1, -1, -1, -1, -3]

@inline function contribute(sampler::OpenSimplex{4}, X, Y, Z, W, x, y, z, w)
    t = sampler.state.table
    g = GRADIENTS_4D
    @inbounds i = t[(W+t[(Z+t[(Y+t[(X&0xff+1)])&0xff+1])&0xff+1]&0xff+1)] & 252 + 1
    a = 2 - x^2 - y^2 - z^2 - w^2
    @inbounds @fastpow a > 0 ? a^4 * (g[i] * x + g[i+1] * y + g[i+2] * z + g[i+3] * w) : 0.0
end

function sample(sampler::OpenSimplex{4}, x::T, y::T, z::T, w::T) where {T<:Real}
    s = (x, y, z, w) .+ ((x + y + z + w) * STRETCH_4D)
    X1, Y1, Z1, W1 = floor.(Int, s)
    x1, y1, z1, w1 = (x, y, z, w) .- ((X1 + Y1 + Z1 + W1) * SQUISH_4D) .- (X1, Y1, Z1, W1)
    Xs, Ys, Zs, Ws = s .- (X1, Y1, Z1, W1)
    XYZWs = Xs + Ys + Zs + Ws
    result = 0.0
    if XYZWs ≤ 1
        p1, p2, s1, s2, Us = 1, 2, Xs, Ys, 1 - XYZWs
        if s1 ≥ s2 && Zs > s2
            s2, p2 = Zs, 4
        elseif s1 < s2 && Zs > s1
            s1, p1 = Zs, 4
        end
        if s1 ≥ s2 && Ws > s2
            s2, p2 = Ws, 8
        elseif s1 < s2 && Ws > s1
            s1, p1 = Ws, 8
        end
        if Us > s1 || Us > s2
            c = s2 > s1 ? p2 : p1
            if iszero(c & 1)
                X2, X3, X4, x12, x13, x14 = X1 - 1, X1, X1, x1 + 1, x1, x1
            else
                X, = X1 + 1
                x = x1 - 1
                X2, X3, X4, x12, x13, x14 = X, X, X, x, x, x
            end
            if iszero(c & 2)
                Y2, Y3, Y4, y12, y13, y14 = Y1, Y1, Y1, y1, y1, y1
                if (c & 1) == 1
                    Y2 -= 1
                    y12 += 1
                else
                    Y3 -= 1
                    y13 += 1
                end
            else
                _Y1, _y = Y1 + 1, y1 - 1
                Y2, Y3, Y4, y12, y13, y14 = _Y1, _Y1, _Y1, _y, _y, _y
            end
            if iszero(c & 4)
                Z12, Z13, Z14, z12, z13, z14 = Z1, Z1, Z1, z1, z1, z1
                if !iszero(c & 3)
                    if (c & 3) !== 3
                        Z13 -= 1
                        z13 += 1
                    end
                else
                    Z14 -= 1
                    z14 += 1
                end
            else
                _Z1, _z = Z1 + 1, z1 - 1
                Z12, Z13, Z14, z12, z13, z14 = _Z1, _Z1, _Z1, _z, _z, _z
            end
            if iszero(c & 8)
                W12, W13, W14, w12, w13, w14 = W1, W1, W1 - 1, w1, w1, w1 + 1
            else
                _W1, _w = W1 + 1, w1 - 1
                W12, W13, W14, w12, w13, w14 = _W1, _W1, _W1, _w, _w, _w
            end
        else
            c = p1 | p2
            if iszero(c & 1)
                X2, X3, X4 = X1, X1 - 1, X1
                x12 = x1 - 2SQUISH_4D
                x13 = x1 + 1 - SQUISH_4D
                x14 = x1 - SQUISH_4D
            else
                X, x = X1 + 1, x1 - 1 - SQUISH_4D
                X2, X3, X4, x13, x14 = X, X, X, x, x
                x12 = x1 - 1 - 2SQUISH_4D
            end
            if iszero(c & 2)
                _y = y1 - SQUISH_4D
                Y2, Y3, Y4, y13, y14 = Y1, Y1, Y1, _y, _y
                y12 = y1 - 2SQUISH_4D
                if (c & 1) == 1
                    Y3 -= 1
                    y13 += 1
                else
                    Y4 -= 1
                    y14 += 1
                end
            else
                _Y1 = Y1 + 1
                _y = y1 - 1 - SQUISH_4D
                Y2, Y3, Y4, y13, y14 = _Y1, _Y1, _Y1, _y, _y
                y12 = y1 - 1 - 2SQUISH_4D
            end
            if iszero(c & 4)
                _z = z1 - SQUISH_4D
                Z12, Z13, Z14, z13, z14 = Z1, Z1, Z1, _z, _z
                z12 = z1 - 2SQUISH_4D
                if (c & 3) == 3
                    Z13 -= 1
                    z13 += 1
                else
                    Z14 -= 1
                    z14 += 1
                end
            else
                _Z1 = Z1 + 1
                _z = z1 - 1 - SQUISH_4D
                Z12, Z13, Z14, z13, z14 = _Z1, _Z1, _Z1, _z, _z
                z12 = z1 - 1 - 2SQUISH_4D
            end
            if iszero(c & 8)
                W12, W13, W14 = W1, W1, W1 - 1
                w12 = w1 - 2SQUISH_4D
                w13 = w1 - SQUISH_4D
                w14 = w1 + 1 - SQUISH_4D
            else
                _W1 = W1 + 1
                _w = w1 - 1 - SQUISH_4D
                W12, W13, W14, w13, w14 = _W1, _W1, _W1, _w, _w
                w12 = w1 - 1 - 2SQUISH_4D
            end
        end
        y2, z2, w2, x3 = (y1, z1, w1, x1) .- SQUISH_4D
        x2, y3, z4, w5 = (x1, y1, z1, w1) .- 1 .- SQUISH_4D
        z3, w3, x4, y4, w4, x5, y5, z5 = (z2, w2, x3, y2, w2, x3, y2, z2)
        c1 = contribute(sampler, X1, Y1, Z1, W1, x1, y1, z1, w1)
        c2 = contribute(sampler, X1 + 1, Y1, Z1, W1, x2, y2, z2, w2)
        c3 = contribute(sampler, X1, Y1 + 1, Z1, W1, x3, y3, z3, w3)
        c4 = contribute(sampler, X1, Y1, Z1 + 1, W1, x4, y4, z4, w4)
        c5 = contribute(sampler, X1, Y1, Z1, W1 + 1, x5, y5, z5, w5)
        result += c1 + c2 + c3 + c4 + c5
    elseif XYZWs ≥ 3
        p1, p2, s1, s2, Us = 14, 13, Xs, Ys, 4 - XYZWs
        if s1 ≤ s2 && Zs < s2
            s2, p2 = Zs, 11
        elseif s1 > s2 && Zs < s1
            s1, p1 = Zs, 11
        end
        if s1 ≤ s2 && Ws < s2
            s2, p2 = Ws, 7
        elseif s1 > s2 && Ws < s1
            s1, p1 = Ws, 7
        end
        if Us < s1 || Us < s2
            c = s2 < s1 ? p2 : p1
            if !iszero(c & 1)
                X = X1 + 1
                x = x1 - 1 - 4SQUISH_4D
                X2 = X1 + 2
                x12 = x1 - 2 - 4SQUISH_4D
                X3, X4, x13, x14 = X, X, x, x
            else
                x = x1 - 4SQUISH_4D
                X2, X3, X4, x12, x13, x14 = X1, X1, X1, x, x, x
            end
            if !iszero(c & 2)
                _Y1 = Y1 + 1
                _y = y1 - 1 - 4SQUISH_4D
                Y2, Y3, Y4, y12, y13, y14 = _Y1, _Y1, _Y1, _y, _y, _y
                if !iszero(c & 1)
                    Y3 += 1
                    y13 -= 1
                else
                    Y2 += 1
                    y12 -= 1
                end
            else
                _y = y1 - 4SQUISH_4D
                Y2, Y3, Y4, y12, y13, y14 = Y1, Y1, Y1, _y, _y, _y
            end
            if !iszero(c & 4)
                _Z1 = Z1 + 1
                _z = z1 - 1 - 4SQUISH_4D
                Z12, Z13, Z14, z12, z13, z14 = _Z1, _Z1, _Z1, _z, _z, _z
                if (c & 3) !== 3
                    if !iszero(c & 3)
                        Z13 += 1
                        z13 -= 1
                    end
                else
                    Z14 += 1
                    z14 -= 1
                end
            else
                _z = z1 - 4SQUISH_4D
                Z12, Z13, Z14, z12, z13, z14 = Z1, Z1, Z1, _z, _z, _z
            end
            if !iszero(c & 8)
                _W1 = W1 + 1
                _w = w1 - 1 - 4SQUISH_4D
                W12, W13, w12, w13 = _W1, _W1, _w, _w
                W14 = W1 + 2
                w14 = w1 - 2 - 4SQUISH_4D
            else
                _w = w1 - 4SQUISH_4D
                W12, W13, W14, w12, w13, w14 = W1, W1, W1, _w, _w, _w
            end
        else
            c = p1 & p2
            if !iszero(c & 1)
                X = X1 + 1
                X2, X3, X4 = X, X1 + 2, X
                x12 = x1 - 1 - 2SQUISH_4D
                x13 = x1 - 2 - 3SQUISH_4D
                x14 = x1 - 1 - 3SQUISH_4D
            else
                x = x1 - 3SQUISH_4D
                X2, X3, X4, x13, x14 = X1, X1, X1, x, x
                x12 = x1 - 2SQUISH_4D
            end
            if !iszero(c & 2)
                _Y1 = Y1 + 1
                _y = y1 - 1 - 3SQUISH_4D
                Y2, Y3, Y4, y13, y14 = _Y1, _Y1, _Y1, _y, _y
                y12 = y1 - 1 - 2SQUISH_4D
                if !iszero(c & 1)
                    Y4 += 1
                    y14 -= 1
                else
                    Y3 += 1
                    y13 -= 1
                end
            else
                _y = y1 - 3SQUISH_4D
                Y2, Y3, Y4, y13, y14 = Y1, Y1, Y1, _y, _y
                y12 = y1 - 2SQUISH_4D
            end
            if !iszero(c & 4)
                _Z1 = Z1 + 1
                _z = z1 - 1 - 3SQUISH_4D
                Z12, Z13, Z14, z13, z14 = _Z1, _Z1, _Z1, _z, _z
                z12 = z1 - 1 - 2SQUISH_4D
                if !iszero(c & 3)
                    Z14 += 1
                    z14 -= 1
                else
                    Z13 += 1
                    z13 -= 1
                end
            else
                _z = z1 - 3SQUISH_4D
                Z12, Z13, Z14, z13, z14 = Z1, Z1, Z1, _z, _z
                z12 = z1 - 2SQUISH_4D
            end
            if !iszero(c & 8)
                _W1 = W1 + 1
                W12, W13, W14 = _W1, _W1, W1 + 2
                w12 = w1 - 1 - 2SQUISH_4D
                w13 = w1 - 1 - 3SQUISH_4D
                w14 = w1 - 2 - 3SQUISH_4D
            else
                _w = w1 - 3SQUISH_4D
                W12, W13, W14, w13, w14 = W1, W1, W1, _w, _w
                w12 = w1 - 2SQUISH_4D
            end
        end
        w4, x5, y5, z5 = (w1, x1, y1, z1) .- 1 .- 3SQUISH_4D
        w5, z4, y3, x2 = (w1, z1, y1, x1) .- 3SQUISH_4D
        x4, y4, x3, z3, w3, y2, z2, w2 = x5, y5, x5, z5, w4, y5, z5, w4
        x1, y1, z1, w1 = (x1, y1, z1, w1) .- 1 .- 4SQUISH_4D
        c1 = contribute(sampler, X1 + 1, Y1 + 1, Z1 + 1, W1, x5, y5, z5, w5)
        c2 = contribute(sampler, X1 + 1, Y1 + 1, Z1, W1 + 1, x4, y4, z4, w4)
        c3 = contribute(sampler, X1 + 1, Y1, Z1 + 1, W1 + 1, x3, y3, z3, w3)
        c4 = contribute(sampler, X1, Y1 + 1, Z1 + 1, W1 + 1, x2, y2, z2, w2)
        c5 = contribute(sampler, X1 + 1, Y1 + 1, Z1 + 1, W1 + 1, x1, y1, z1, w1)
        result += c1 + c2 + c3 + c4 + c5
    elseif XYZWs ≤ 2
        p1, p2, s1, s2, p1_bigger, p2_bigger = 0, 0, 0.0, 0.0, true, true
        p3 = 2 - XYZWs + Xs
        p4 = 2 - XYZWs + Ys
        p5 = 2 - XYZWs + Zs
        p6 = 2 - XYZWs + Ws
        XYs, ZWs, XZs, YWs, XWs, YZs = Xs + Ys, Zs + Ws, Xs + Zs, Ys + Ws, Xs + Ws, Ys + Zs
        if XYs > ZWs
            s1, p1 = XYs, 3
        else
            s1, p1 = ZWs, 12
        end
        if XZs > YWs
            s2, p2 = XZs, 5
        else
            s2, p2 = YWs, 10
        end
        if XWs > YZs
            s = XWs
            if s1 ≥ s2 && s > s2
                s2, p2 = s, 9
            elseif s1 < s2 && s > s1
                s1, p1 = s, 9
            end
        else
            s = YZs
            if s1 ≥ s2 && s > s2
                s2, p2 = s, 6
            elseif s1 < s2 && s > s1
                s1, p1 = s, 6
            end
        end
        if s1 ≥ s2 && p3 > s2
            s2, p2, p2_bigger = p3, 1, false
        elseif s1 < s2 && p3 > s1
            s1, p1, p1_bigger = p3, 1, false
        end
        if s1 ≥ s2 && p4 > s2
            s2, p2, p2_bigger = p4, 2, false
        elseif s1 < s2 && p4 > s1
            s1, p1, p1_bigger = p4, 2, false
        end
        if s1 ≥ s2 && p5 > s2
            s2, p2, p2_bigger = p5, 4, false
        elseif s1 < s2 && p5 > s1
            s1, p1, p1_bigger = p5, 4, false
        end
        if s1 ≥ s2 && p6 > s2
            s2, p2, p2_bigger = p6, 8, false
        elseif s1 < s2 && p6 > s1
            s1, p1, p1_bigger = p6, 8, false
        end
        if p1_bigger === p2_bigger
            if p1_bigger === true
                c1 = p1 | p2
                c2 = p1 & p2
                if iszero(c1 & 1)
                    X2, X3 = X1, X1 - 1
                    x12 = x1 - 3SQUISH_4D
                    x13 = x1 + 1 - 2SQUISH_4D
                else
                    X = X1 + 1
                    X2, X3 = X, X
                    x12 = x1 - 1 - 3SQUISH_4D
                    x13 = x1 - 1 - 2SQUISH_4D
                end
                if iszero(c1 & 2)
                    Y2, Y3 = Y1, Y1 - 1
                    y12 = y1 - 3SQUISH_4D
                    y13 = y1 + 1 - 2SQUISH_4D
                else
                    _Y1 = Y1 + 1
                    Y2, Y3 = _Y1, _Y1
                    y12 = y1 - 1 - 3SQUISH_4D
                    y13 = y1 - 1 - 2SQUISH_4D
                end
                if iszero(c1 & 4)
                    Z12, Z13 = Z1, Z1 - 1
                    z12 = z1 - 3SQUISH_4D
                    z13 = z1 + 1 - 2SQUISH_4D
                else
                    _Z1 = Z1 + 1
                    Z12, Z13 = _Z1, _Z1
                    z12 = z1 - 1 - 3SQUISH_4D
                    z13 = z1 - 1 - 2SQUISH_4D
                end
                if iszero(c1 & 8)
                    W12, W13 = W1, W1 - 1
                    w12 = w1 - 3SQUISH_4D
                    w13 = w1 + 1 - 2SQUISH_4D
                else
                    _W1 = W1 + 1
                    W12, W13 = _W1, _W1
                    w12 = w1 - 1 - 3SQUISH_4D
                    w13 = w1 - 1 - 2SQUISH_4D
                end
                X4, Y4, Z14, W14 = X1, Y1, Z1, W1
                x14, y14, z14, w14 = (x1, y1, z1, w1) .- 2SQUISH_4D
                if !iszero(c2 & 1)
                    X4 += 2
                    x14 -= 2
                elseif !iszero(c2 & 2)
                    Y4 += 2
                    y14 -= 2
                elseif !iszero(c2 & 4)
                    Z14 += 2
                    z14 -= 2
                else
                    W14 += 2
                    w14 -= 2
                end
            else
                c = p1 | p2
                X4, Y4, Z14, W14 = X1, Y1, Z1, W1
                x14, y14, z14, w14 = x1, y1, z1, w1
                if iszero(c & 1)
                    X2, X3 = X1 - 1, X1
                    x12 = x1 + 1 - SQUISH_4D
                    x13 = x1 - SQUISH_4D
                else
                    X = X1 + 1
                    x = x1 - 1 - SQUISH_4D
                    X2, X3, x12, x13 = X, X, x, x
                end
                if iszero(c & 2)
                    _y = y1 - SQUISH_4D
                    Y2, Y3, y12, y13 = Y1, Y1, _y, _y
                    if (c & 1) == 1
                        Y2 -= 1
                        y12 += 1
                    else
                        Y3 -= 1
                        y13 += 1
                    end
                else
                    _Y1 = Y1 + 1
                    _y = y1 - 1 - SQUISH_4D
                    Y2, Y3, y12, y13 = _Y1, _Y1, _y, _y
                end
                if iszero(c & 4)
                    _z = z1 - SQUISH_4D
                    Z12, Z13, z12, z13 = Z1, Z1, _z, _z
                    if (c & 3) == 3
                        Z12 -= 1
                        z12 += 1
                    else
                        Z13 -= 1
                        z13 += 1
                    end
                else
                    _Z1 = Z1 + 1
                    _z = z1 - 1 - SQUISH_4D
                    Z12, Z13, z12, z13 = _Z1, _Z1, _z, _z
                end
                if iszero(c & 8)
                    W12, W13 = W1, W1 - 1
                    w12 = w1 - SQUISH_4D
                    w13 = w1 + 1 - SQUISH_4D
                else
                    _W1 = W1 + 1
                    _w = w1 - 1 - SQUISH_4D
                    W12, W13, w12, w13 = _W1, _W1, _w, _w
                end
            end
        else
            c1 = p1_bigger === true ? p1 : p2
            c2 = p1_bigger === true ? p2 : p1
            if iszero(c1 & 1)
                X2, X3 = X1 - 1, X1
                x12 = x1 + 1 - SQUISH_4D
                x13 = x1 - SQUISH_4D
            else
                X = X1 + 1
                x = x1 - 1 - SQUISH_4D
                X2, X3, x12, x13 = X, X, x, x
            end
            if iszero(c1 & 2)
                _y = y1 - SQUISH_4D
                Y2, Y3, y12, y13 = Y1, Y1, _y, _y
                if (c1 & 1) == 1
                    Y2 -= 1
                    y12 += 1
                else
                    Y3 -= 1
                    y13 += 1
                end
            else
                _Y1 = Y1 + 1
                _y = y1 - 1 - SQUISH_4D
                Y2, Y3, y12, y13 = _Y1, _Y1, _y, _y
            end
            if iszero(c1 & 4)
                _z = z1 - SQUISH_4D
                Z12, Z13, z12, z13 = Z1, Z1, _z, _z
                if (c1 & 3) == 3
                    Z12 -= 1
                    z12 += 1
                else
                    Z13 -= 1
                    z13 += 1
                end
            else
                _Z1 = Z1 + 1
                _z = z1 - 1 - SQUISH_4D
                Z12, Z13, z12, z13 = _Z1, _Z1, _z, _z
            end
            if iszero(c1 & 8)
                W12, W13 = W1, W1 - 1
                w12 = w1 - SQUISH_4D
                w13 = w1 + 1 - SQUISH_4D
            else
                _W1 = W1 + 1
                _w = w1 - 1 - SQUISH_4D
                W12, W13, w12, w13 = _W1, _W1, _w, _w
            end
            X4, Y4, Z14, W14 = X1, Y1, Z1, W1
            x14, y14, z14, w14 = (x1, y1, z1, w1) .- 2SQUISH_4D
            if !iszero(c2 & 1)
                X4 += 2
                x14 -= 2
            elseif !iszero(c2 & 2)
                Y4 += 2
                y14 -= 2
            elseif !iszero(c2 & 4)
                Z14 += 2
                z14 -= 2
            else
                W14 += 2
                w14 -= 2
            end
        end
        y2, z2, w2, x3 = (y1, z1, w1, x1) .- SQUISH_4D
        z6, w6, y7, x9 = (z1, w1, y1, x1) .- 2SQUISH_4D
        x6, y6, z7, w8 = (x1, y1, z1, w1) .- 1 .- 2SQUISH_4D
        x2, y3, z4, w5 = (x1, y1, z1, w1) .- 1 .- SQUISH_4D
        z3, w3, x4, y4, w4, x5, y5, z5 = z2, w2, x3, y2, w2, x3, y2, z2
        x7, w7, x8, y8, z8, y9, z9, w9 = x6, w6, x6, y7, z6, y6, z7, w6
        x10, y10, z10, w10, x11, y11, z11, w11 = x9, y6, z6, w8, x9, y7, z7, w8
        c1 = contribute(sampler, X1 + 1, Y1, Z1, W1, x2, y2, z2, w2)
        c2 = contribute(sampler, X1, Y1 + 1, Z1, W1, x3, y3, z3, w3)
        c3 = contribute(sampler, X1, Y1, Z1 + 1, W1, x4, y4, z4, w4)
        c4 = contribute(sampler, X1, Y1, Z1, W1 + 1, x5, y5, z5, w5)
        c5 = contribute(sampler, X1 + 1, Y1 + 1, Z1, W1, x6, y6, z6, w6)
        c6 = contribute(sampler, X1 + 1, Y1, Z1 + 1, W1, x7, y7, z7, w7)
        c7 = contribute(sampler, X1 + 1, Y1, Z1, W1 + 1, x8, y8, z8, w8)
        c8 = contribute(sampler, X1, Y1 + 1, Z1 + 1, W1, x9, y9, z9, w9)
        c9 = contribute(sampler, X1, Y1 + 1, Z1, W1 + 1, x10, y10, z10, w10)
        c10 = contribute(sampler, X1, Y1, Z1 + 1, W1 + 1, x11, y11, z11, w11)
        result += c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10
    else
        p1, p2, s1, s2, p1_bigger, p2_bigger = 0, 0, 0.0, 0.0, true, true
        p3 = 3 - XYZWs + Xs
        p4 = 3 - XYZWs + Ys
        p5 = 3 - XYZWs + Zs
        p6 = 3 - XYZWs + Ws
        XYs, ZWs, XZs, YWs, XWs, YZs = Xs + Ys, Zs + Ws, Xs + Zs, Ys + Ws, Xs + Ws, Ys + Zs
        if XYs < ZWs
            s1, p1 = XYs, 12
        else
            s1, p1 = ZWs, 3
        end
        if XZs < YWs
            s2, p2 = XZs, 10
        else
            s2, p2 = YWs, 5
        end
        if XWs < YZs
            s = XWs
            if s1 ≤ s2 && s < s2
                s2, p2 = s, 6
            elseif s1 > s2 && s < s1
                s1, p1 = s, 6
            end
        else
            s = YZs
            if s1 ≤ s2 && s < s2
                s2, p2 = s, 9
            elseif s1 > s2 && s < s1
                s1, p1 = s, 9
            end
        end
        if s1 ≤ s2 && p3 < s2
            s2, p2, p2_bigger = p3, 14, false
        elseif s1 > s2 && p3 < s1
            s1, p1, p1_bigger = p3, 14, false
        end
        if s1 ≤ s2 && p4 < s2
            s2, p2, p2_bigger = p4, 13, false
        elseif s1 > s2 && p4 < s1
            s1, p1, p1_bigger = p4, 13, false
        end
        if s1 ≤ s2 && p5 < s2
            s2, p2, p2_bigger = p5, 11, false
        elseif s1 > s2 && p5 < s1
            s1, p1, p1_bigger = p5, 11, false
        end
        if s1 ≤ s2 && p6 < s2
            s2, p2, p2_bigger = p6, 7, false
        elseif s1 > s2 && p6 < s1
            s1, p1, p1_bigger = p6, 7, false
        end
        if p1_bigger === p2_bigger
            if p1_bigger === true
                c1 = p1 & p2
                c2 = p1 | p2
                X2, X3, Y2, Y3, Z12, Z13, W12, W13 = X1, X1, Y1, Y1, Z1, Z1, W1, W1
                x12, y12, z12, w12 = (x1, y1, z1, w1) .- SQUISH_4D
                x13, y13, z13, w13 = (x1, y1, z1, w1) .- 2SQUISH_4D
                if !iszero(c1 & 1)
                    X2 += 1
                    x12 -= 1
                    X3 += 2
                    x13 -= 2
                elseif !iszero(c1 & 2)
                    Y2 += 1
                    y12 -= 1
                    Y3 += 2
                    y13 -= 2
                elseif !iszero(c1 & 4)
                    Z12 += 1
                    z12 -= 1
                    Z13 += 2
                    z13 -= 2
                else
                    W12 += 1
                    w12 -= 1
                    W13 += 2
                    w13 -= 2
                end
                X4, Y4, Z14, W14 = X1 + 1, Y1 + 1, Z1 + 1, W1 + 1
                x14, y14, z14, w14 = (x1, y1, z1, w1) .- 1 .- 2SQUISH_4D
                if iszero(c2 & 1)
                    X4 -= 2
                    x14 += 2
                elseif iszero(c2 & 2)
                    Y4 -= 2
                    y14 += 2
                elseif iszero(c2 & 4)
                    Z14 -= 2
                    z14 += 2
                else
                    W14 -= 2
                    w14 += 2
                end
            else
                c = p1 & p2
                X4, Y4, Z14, W14 = X1 + 1, Y1 + 1, Z1 + 1, W1 + 1
                x14, y14, z14, w14 = (x1, y1, z1, w1) .- 1 .- 4SQUISH_4D
                if !iszero(c & 1)
                    X2, X3 = X1 + 2, X1 + 1
                    x12 = x1 - 2 - 3SQUISH_4D
                    x13 = x1 - 1 - 3SQUISH_4D
                else
                    x = x1 - 3SQUISH_4D
                    X2, X3, x12, x13 = X1, X1, x, x
                end
                if !iszero(c & 2)
                    _Y1 = Y1 + 1
                    _y = y1 - 1 - 3SQUISH_4D
                    Y2, Y3, y12, y13 = _Y1, _Y1, _y, _y
                    if iszero(c & 1)
                        Y2 += 1
                        y12 -= 1
                    else
                        Y3 += 1
                        y13 -= 1
                    end
                else
                    _y = y1 - 3SQUISH_4D
                    Y2, Y3, y12, y13 = Y1, Y1, _y, _y
                end
                if !iszero(c & 4)
                    _Z1 = Z1 + 1
                    _z = z1 - 1 - 3SQUISH_4D
                    Z12, Z13, z12, z13 = _Z1, _Z1, _z, _z
                    if iszero(c & 3)
                        Z12 += 1
                        z12 -= 1
                    else
                        Z13 += 1
                        z13 -= 1
                    end
                else
                    _z = z1 - 3SQUISH_4D
                    Z12, Z13, z12, z13 = Z1, Z1, _z, _z
                end
                if !iszero(c & 8)
                    W12, W13 = W1 + 1, W1 + 2
                    w12 = w1 - 1 - 3SQUISH_4D
                    w13 = w1 - 2 - 3SQUISH_4D
                else
                    _w = w1 - 3SQUISH_4D
                    W12, W13, w12, w13 = W1, W1, _w, _w
                end
            end
        else
            c1 = p1_bigger === true ? p1 : p2
            c2 = p1_bigger === true ? p2 : p1
            if !iszero(c1 & 1)
                X2, X3 = X1 + 2, X1 + 1
                x12 = x1 - 2 - 3SQUISH_4D
                x13 = x1 - 1 - 3SQUISH_4D
            else
                x = x1 - 3SQUISH_4D
                X2, X3, x12, x13 = X1, X1, x, x
            end
            if !iszero(c1 & 2)
                _Y1 = Y1 + 1
                _y = y1 - 1 - 3SQUISH_4D
                Y2, Y3, y12, y13 = _Y1, _Y1, _y, _y
                if iszero(c1 & 1)
                    Y2 += 1
                    y12 -= 1
                else
                    Y3 += 1
                    y13 -= 1
                end
            else
                _y = y1 - 3SQUISH_4D
                Y2, Y3, y12, y13 = Y1, Y1, _y, _y
            end
            if !iszero(c1 & 4)
                _Z1 = Z1 + 1
                _z = z1 - 1 - 3SQUISH_4D
                Z12, Z13, z12, z13 = _Z1, _Z1, _z, _z
                if iszero(c1 & 3)
                    Z12 += 1
                    z12 -= 1
                else
                    Z13 += 1
                    z13 -= 1
                end
            else
                _z = z1 - 3SQUISH_4D
                Z12, Z13, z12, z13 = Z1, Z1, _z, _z
            end
            if !iszero(c1 & 8)
                W12, W13 = W1 + 1, W1 + 2
                w12 = w1 - 1 - 3SQUISH_4D
                w13 = w1 - 2 - 3SQUISH_4D
            else
                _w = w1 - 3SQUISH_4D
                W12, W13, w12, w13 = W1, W1, _w, _w
            end
            X4, Y4, Z14, W14 = X1 + 1, Y1 + 1, Z1 + 1, W1 + 1
            x14, y14, z14, w14 = (x1, y1, z1, w1) .- 1 .- 2SQUISH_4D
            if iszero(c2 & 1)
                X4 -= 2
                x14 += 2
            elseif iszero(c2 & 2)
                Y4 -= 2
                y14 += 2
            elseif iszero(c2 & 4)
                Z14 -= 2
                z14 += 2
            else
                W14 -= 2
                w14 += 2
            end
        end
        w4, x5, y5, z5 = (w1, x1, y1, z1) .- 1 .- 3SQUISH_4D
        x6, y6, z7, w8 = (x1, y1, z1, w1) .- 1 .- 2SQUISH_4D
        z6, w6, y7, x9 = (z1, w1, y1, x1) .- 2SQUISH_4D
        w5, z4, y3, x2 = (w1, z1, y1, x1) .- 3SQUISH_4D
        x4, y4, x3, z3, w3, y2, z2, w2 = x5, y5, x5, z5, w4, y5, z5, w4
        x7, w7, x8, y8, z8, y9, z9, w9 = x6, w6, x6, y7, z6, y6, z7, w6
        x10, y10, z10, w10, x11, y11, z11, w11 = x9, y6, z6, w8, x9, y7, z7, w8
        c1 = contribute(sampler, X1 + 1, Y1 + 1, Z1 + 1, W1, x5, y5, z5, w5)
        c2 = contribute(sampler, X1 + 1, Y1 + 1, Z1, W1 + 1, x4, y4, z4, w4)
        c3 = contribute(sampler, X1 + 1, Y1, Z1 + 1, W1 + 1, x3, y3, z3, w3)
        c4 = contribute(sampler, X1, Y1 + 1, Z1 + 1, W1 + 1, x2, y2, z2, w2)
        c5 = contribute(sampler, X1 + 1, Y1 + 1, Z1, W1, x6, y6, z6, w6)
        c6 = contribute(sampler, X1 + 1, Y1, Z1 + 1, W1, x7, y7, z7, w7)
        c7 = contribute(sampler, X1 + 1, Y1, Z1, W1 + 1, x8, y8, z8, w8)
        c8 = contribute(sampler, X1, Y1 + 1, Z1 + 1, W1, x9, y9, z9, w9)
        c9 = contribute(sampler, X1, Y1 + 1, Z1, W1 + 1, x10, y10, z10, w10)
        c10 = contribute(sampler, X1, Y1, Z1 + 1, W1 + 1, x11, y11, z11, w11)
        result += c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10
    end
    c1 = contribute(sampler, X2, Y2, Z12, W12, x12, y12, z12, w12)
    c2 = contribute(sampler, X3, Y3, Z13, W13, x13, y13, z13, w13)
    c3 = contribute(sampler, X4, Y4, Z14, W14, x14, y14, z14, w14)
    (result + c1 + c2 + c3) * SCALE_4D
end
