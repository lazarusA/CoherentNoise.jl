const STRETCH_3D = 1 / -6

const SQUISH_3D = 1 / 3

const SCALE_3D = 1 / 103

const GRADIENTS_3D =
    Int8.((
        -11, 4, 4, -4, 11, 4, -4, 4, 11, 11, 4, 4, 4, 11, 4, 4, 4, 11, -11, -4, 4, -4, -11, 4, -4,
        -4, 11, 11, -4, 4, 4, -11, 4, 4, -4, 11, -11, 4, -4, -4, 11, -4, -4, 4, -11, 11, 4, -4, 4,
        11, -4, 4, 4, -11, -11, -4, -4, -4, -11, -4, -4, -4, -11, 11, -4, -4, 4, 11, -4, 4, -4, -11,
    ))

@inline function contribute(sampler::OpenSimplex{3}, X, Y, Z, x, y, z)
    t = sampler.table
    g = GRADIENTS_3D
    @inbounds i = t[(Z+t[(Y+t[(X&0xff+1)])&0xff+1])&0xff+1] + 1
    a = 2 - x^2 - y^2 - z^2
    @inbounds @fastpow a > 0 ? a^4 * (g[i] * x + g[i+1] * y + g[i+2] * z) : 0.0
end

function sample(sampler::OpenSimplex{3}, x::T, y::T, z::T) where {T<:Real}
    s = (x, y, z) .+ ((x + y + z) * STRETCH_3D)
    X1, Y1, Z1 = floor.(Int, s)
    x1, y1, z1 = (x, y, z) .- ((X1 + Y1 + Z1) * SQUISH_3D) .- (X1, Y1, Z1)
    Xs, Ys, Zs = s .- (X1, Y1, Z1)
    XYZs = Xs + Ys + Zs
    result = 0.0
    if XYZs ≤ 1
        p1, p2, s1, s2, Ws = 1, 2, Xs, Ys, 1 - XYZs
        if s1 ≥ s2 && Zs > s2
            s2, p2 = Zs, 4
        elseif s1 < s2 && Zs > s1
            s1, p1 = Zs, 4
        end
        if Ws > s1 || Ws > s2
            c = s2 > s1 ? p2 : p1
            if iszero(c & 1)
                X2, X3, x8, x9 = X1 - 1, X1, x1 + 1, x1
            else
                X, x = X1 + 1, x1 - 1
                X2, X3, x8, x9 = X, X, x, x
            end
            if iszero(c & 2)
                Y2, Y3, y8, y9 = Y1, Y1, y1, y1
                if iszero(c & 1)
                    Y3 -= 1
                    y9 += 1
                else
                    Y2 -= 1
                    y8 += 1
                end
            else
                Y, y = Y1 + 1, y1 - 1
                Y2, Y3, y8, y9 = Y, Y, y, y
            end
            if iszero(c & 4)
                Z2, Z3, z8, z9 = Z1, Z1 - 1, z1, z1 + 1
            else
                Z, z = Z1 + 1, z1 - 1
                Z2, Z3, z8, z9 = Z, Z, z, z
            end
        else
            c = (p1 | p2) & 0xff
            if iszero(c & 1)
                X2, X3 = X1, X1 - 1
                x8 = x1 - 2SQUISH_3D
                x9 = x1 + 1 - SQUISH_3D
            else
                X = X1 + 1
                X2, X3 = X, X
                x8 = x1 - 1 - 2SQUISH_3D
                x9 = x1 - 1 - SQUISH_3D
            end
            if iszero(c & 2)
                Y2, Y3 = Y1, Y1 - 1
                y8 = y1 - 2SQUISH_3D
                y9 = y1 + 1 - SQUISH_3D
            else
                Y = Y1 + 1
                Y2, Y3 = Y, Y
                y8 = y1 - 1 - 2SQUISH_3D
                y9 = y1 - 1 - SQUISH_3D
            end
            if iszero(c & 4)
                Z2, Z3 = Z1, Z1 - 1
                z8 = z1 - 2SQUISH_3D
                z9 = z1 + 1 - SQUISH_3D
            else
                Z = Z1 + 1
                Z2, Z3 = Z, Z
                z8 = z1 - 1 - 2SQUISH_3D
                z9 = z1 - 1 - SQUISH_3D
            end
        end
        y2, z2, x3 = (y1, z1, x1) .- SQUISH_3D
        x2, y3, z4 = (x1, y1, z1) .- 1 .- SQUISH_3D
        z3, x4, y4 = z2, x3, y2
        c1 = contribute(sampler, X1, Y1, Z1, x1, y1, z1)
        c2 = contribute(sampler, X1 + 1, Y1, Z1, x2, y2, z2)
        c3 = contribute(sampler, X1, Y1 + 1, Z1, x3, y3, z3)
        c4 = contribute(sampler, X1, Y1, Z1 + 1, x4, y4, z4)
        result += c1 + c2 + c3 + c4
    elseif XYZs ≥ 2
        p1, p2, s1, s2, Ws = 6, 5, Xs, Ys, 3 - XYZs
        if s1 ≤ s2 && Zs < s2
            s2, p2 = Zs, 3
        elseif s1 > s2 && Zs < s1
            s1, p1 = Zs, 3
        end
        if Ws < s1 || Ws < s2
            c = s2 < s1 ? p2 : p1
            if !iszero(c & 1)
                X2, X3 = X1 + 2, X1 + 1
                x8 = x1 - 2 - 3SQUISH_3D
                x9 = x1 - 1 - 3SQUISH_3D
            else
                x = x1 - 3SQUISH_3D
                X2, X3, x8, x9 = X1, X1, x, x
            end
            if !iszero(c & 2)
                Y = Y1 + 1
                y = y1 - 1 - 3SQUISH_3D
                Y2, Y3, y8, y9 = Y, Y, y, y
                if !iszero(c & 1)
                    Y3 += 1
                    y9 -= 1
                else
                    Y2 += 1
                    y8 -= 1
                end
            else
                y = y1 - 3SQUISH_3D
                Y2, Y3, y8, y9 = Y1, Y1, y, y
            end
            if !iszero(c & 4)
                Z2, Z3 = Z1 + 1, Z1 + 2
                z8 = z1 - 1 - 3SQUISH_3D
                z9 = z1 - 2 - 3SQUISH_3D
            else
                z = z1 - 3SQUISH_3D
                Z2, Z3, z8, z9 = Z1, Z1, z, z
            end
        else
            c = (p1 & p2) & 0xff
            if !iszero(c & 1)
                X2, X3 = X1 + 1, X1 + 2
                x8 = x1 - 1 - SQUISH_3D
                x9 = x1 - 2 - 2SQUISH_3D
            else
                X2, X3 = X1, X1
                x8 = x1 - SQUISH_3D
                x9 = x1 - 2SQUISH_3D
            end
            if !iszero(c & 2)
                Y2, Y3 = Y1 + 1, Y1 + 2
                y8 = y1 - 1 - SQUISH_3D
                y9 = y1 - 2 - 2SQUISH_3D
            else
                Y2, Y3 = Y1, Y1
                y8 = y1 - SQUISH_3D
                y9 = y1 - 2SQUISH_3D
            end
            if !iszero(c & 4)
                Z2, Z3 = Z1 + 1, Z1 + 2
                z8 = z1 - 1 - SQUISH_3D
                z9 = z1 - 2 - 2SQUISH_3D
            else
                Z2, Z3 = Z1, Z1
                z8 = z1 - SQUISH_3D
                z9 = z1 - 2SQUISH_3D
            end
        end
        z3, x4, y4 = (z1, x1, y1) .- 1 .- 2SQUISH_3D
        x3, y2, z2 = x4, y4, z3
        z4, y3, x2 = (z1, y1, x1) .- 2SQUISH_3D
        x1, y1, z1 = (x1, y1, z1) .- 1 .- 3SQUISH_3D
        c1 = contribute(sampler, X1 + 1, Y1 + 1, Z1, x4, y4, z4)
        c2 = contribute(sampler, X1 + 1, Y1, Z1 + 1, x3, y3, z3)
        c3 = contribute(sampler, X1, Y1 + 1, Z1 + 1, x2, y2, z2)
        c4 = contribute(sampler, X1 + 1, Y1 + 1, Z1 + 1, x1, y1, z1)
        result += c1 + c2 + c3 + c4
    else
        p3, p4, p5 = Xs + Ys, Xs + Zs, Ys + Zs
        if p3 > 1
            s1, p1, p1_farthest = p3 - 1, 3, true
        else
            s1, p1, p1_farthest = 1 - p3, 4, false
        end
        if p4 > 1
            s2, p2, p2_farthest = p4 - 1, 5, true
        else
            s2, p2, p2_farthest = 1 - p4, 2, false
        end
        if p5 > 1
            s = p5 - 1
            if s1 ≤ s2 && s1 < s
                s1, p1, p1_farthest = s, 6, true
            elseif s1 > s2 && s2 < s
                s2, p2, p2_farthest = s, 6, true
            end
        else
            s = 1 - p5
            if s1 ≤ s2 && s1 < s
                s1, p1, p1_farthest = s, 1, false
            elseif s1 > s2 && s2 < s
                s2, p2, p2_farthest = s, 1, false
            end
        end
        if p1_farthest === p2_farthest
            if p1_farthest === true
                c = p1 & p2
                x8, y8, z8 = (x1, y1, z1) .- 1 .- 3SQUISH_3D
                X2, Y2, Z2 = Z1 + 1, Y1 + 1, Z1 + 1
                if !iszero(c & 1)
                    x9 = x1 - 2 - 2SQUISH_3D
                    y9, z9 = (y1, z1) .- 2SQUISH_3D
                    X3, Y3, Z3 = X1 + 2, Y1, Z1
                elseif !iszero(c & 2)
                    x9, z9 = (x1, z1) .- 2SQUISH_3D
                    y9 = y1 - 2 - 2SQUISH_3D
                    X3, Y3, Z3 = X1, Y1 + 2, Z1
                else
                    x9, y9 = (x1, y1) .- 2SQUISH_3D
                    z9 = z1 - 2 - 2SQUISH_3D
                    X3, Y3, Z3 = X1, Y1, Z1 + 2
                end
            else
                c = p1 | p2
                x8, y8, z8, X2, Y2, Z2 = x1, y1, z1, X1, Y1, Z1
                if iszero(c & 1)
                    x9 = x1 + 1 - SQUISH_3D
                    y9, z9 = (y1, z1) .- 1 .- SQUISH_3D
                    X3, Y3, Z3 = X1 - 1, Y1 + 1, Z1 + 1
                elseif iszero(c & 2)
                    x9, z9 = (x1, z1) .- 1 .- SQUISH_3D
                    y9 = y1 + 1 - SQUISH_3D
                    X3, Y3, Z3 = X1 + 1, Y1 - 1, Z1 + 1
                else
                    x9, y9 = (x1, y1) .- 1 .- SQUISH_3D
                    z9 = z1 + 1 - SQUISH_3D
                    X3, Y3, Z3 = X1 + 1, Y1 + 1, Z1 - 1
                end
            end
        else
            c1 = p1_farthest ? p1 : p2
            c2 = p1_farthest ? p2 : p1
            if iszero(c1 & 1)
                x8 = x1 + 1 - SQUISH_3D
                y8, z8 = (y1, z1) .- 1 .- SQUISH_3D
                X2, Y2, Z2 = X1 - 1, Y1 + 1, Z1 + 1
            elseif iszero(c1 & 2)
                x8, z8 = (x1, z1) .- 1 .- SQUISH_3D
                y8 = y1 + 1 - SQUISH_3D
                X2, Y2, Z2 = X1 + 1, Y1 - 1, Z1 + 1
            else
                x8, y8 = (x1, y1) .- 1 .- SQUISH_3D
                z8 = z1 + 1 - SQUISH_3D
                X2, Y2, Z2 = X1 + 1, Y1 + 1, Z1 - 1
            end
            x9, y9, z9 = (x1, y1, z1) .- 2SQUISH_3D
            X3, Y3, Z3 = X1, Y1, Z1
            if !iszero(c2 & 1)
                x9 -= 2
                X3 += 2
            elseif !iszero(c2 & 2)
                y9 -= 2
                Y3 += 2
            else
                z9 -= 2
                Z3 += 2
            end
        end
        y2, z2, x3 = (y1, z1, x1) .- SQUISH_3D
        x5, y5, z6 = (x1, y1, z1) .- 1 .- 2SQUISH_3D
        x2, y3, z4 = (x1, y1, z1) .- 1 .- SQUISH_3D
        z5, y6, x7 = (z1, y1, x1) .- 2SQUISH_3D
        z3, x4, y4, x6, y7, z7 = z2, x3, y2, x5, y5, z6
        c1 = contribute(sampler, X1 + 1, Y1, Z1, x2, y2, z2)
        c2 = contribute(sampler, X1, Y1 + 1, Z1, x3, y3, z3)
        c3 = contribute(sampler, X1, Y1, Z1 + 1, x4, y4, z4)
        c4 = contribute(sampler, X1 + 1, Y1 + 1, Z1, x5, y5, z5)
        c5 = contribute(sampler, X1 + 1, Y1, Z1 + 1, x6, y6, z6)
        c6 = contribute(sampler, X1, Y1 + 1, Z1 + 1, x7, y7, z7)
        result += c1 + c2 + c3 + c4 + c5 + c6
    end
    c1 = contribute(sampler, X2, Y2, Z2, x8, y8, z8)
    c2 = contribute(sampler, X3, Y3, Z3, x9, y9, z9)
    (result + c1 + c2) * SCALE_3D
end
