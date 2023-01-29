struct OpenSimplex{N} <: NoiseSampler{N}
    random_state::RandomState
    perlin_state::PerlinState
    simplex_state::SimplexState
    table::Vector{UInt8}
end

@inline function _opensimplex(dims, seed, smooth)
    T = OpenSimplex{dims}
    rs = RandomState(seed)
    # 0:3:45 and 48:3:69
    table = shuffle(rs.rng, [repeat(0x00:0x03:0x2d, 11); repeat(0x30:0x03:0x45, 10)])
    T(rs, PerlinState(rs), SimplexState(T, Val(smooth)), table)
end

SimplexState(::Type{<:OpenSimplex{2}}, ::Val{true}) = SimplexState(1, 0.8923035729617332)
SimplexState(::Type{<:OpenSimplex{2}}, ::Val{false}) = SimplexState(2, 1 / 40.7)
SimplexState(::Type{<:OpenSimplex{3}}, ::Val{true}) = SimplexState(1, 0.3777736183312849)
SimplexState(::Type{<:OpenSimplex{3}}, ::Val{false}) = SimplexState(2, 1 / 103)
SimplexState(::Type{<:OpenSimplex{4}}, ::Val{true}) = SimplexState(1, 1.4485981091431033)
SimplexState(::Type{<:OpenSimplex{4}}, ::Val{false}) = SimplexState(2, 1 / 30)

# 2D

const OS_GRADIENTS_2D = Int8.((5, 2, 2, 5, -5, 2, -2, 5, 5, -2, 2, -5, -5, -2, -2, -5))
const OS_STRETCH_2D = (1 / sqrt(3) - 1) / 2
const OS_SQUISH_2D = (sqrt(3) - 1) / 2

"""
    opensimplex_2d(; seed=nothing, smooth=false)

Construct a sampler that outputs 2-dimensional legacy OpenSimplex noise when it is sampled from.

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

# See also:

  - [`opensimplex2_2d`](@ref opensimplex2_2d)
  - [`openSimplex2s_2d`](@ref opensimplex2s_2d)
"""
opensimplex_2d(; seed=nothing, smooth=false) = _opensimplex(2, seed, smooth)

@inline function contribute(sampler::OpenSimplex{2}, falloff, X, Y, x, y)
    t = sampler.perlin_state.table
    g = OS_GRADIENTS_2D
    @inbounds i = (t[t[X+1]+Y+1] & 14) + 1
    a = falloff - x^2 - y^2
    @inbounds a > 0 ? pow4(a) * (g[i] * x + g[i+1] * y) : 0.0
end

function sample(sampler::OpenSimplex{2}, x::T, y::T) where {T<:Real}
    s = (x, y) .+ ((x + y) * OS_STRETCH_2D)
    state = sampler.simplex_state
    falloff = state.falloff
    X0, Y0 = floor.(Int, s)
    X1, Y1 = (X0, Y0) .& 255 .+ 1
    x1, y1 = (x, y) .- ((X0 + Y0) * OS_SQUISH_2D) .- (X0, Y0)
    y2, x3 = (y1, x1) .- OS_SQUISH_2D
    x2, y3 = (x1, y1) .- 1 .- OS_SQUISH_2D
    Xs, Ys = s .- (X0, Y0)
    XYs = Xs + Ys
    c1 = contribute(sampler, falloff, X1 + 1, Y1, x2, y2)
    c2 = contribute(sampler, falloff, X1, Y1 + 1, x3, y3)
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
    c1 = contribute(sampler, falloff, X1, Y1, x1, y1)
    c2 = contribute(sampler, falloff, X2, Y2, x4, y4)
    (result + c1 + c2) * state.scale_factor
end

# 3D

const OS_STRETCH_3D = 1 / -6
const OS_SQUISH_3D = 1 / 3
const OS_GRADIENTS_3D =
    Int8.((
        -11, 4, 4, -4, 11, 4, -4, 4, 11, 11, 4, 4, 4, 11, 4, 4, 4, 11, -11, -4, 4, -4, -11, 4, -4,
        -4, 11, 11, -4, 4, 4, -11, 4, 4, -4, 11, -11, 4, -4, -4, 11, -4, -4, 4, -11, 11, 4, -4, 4,
        11, -4, 4, 4, -11, -11, -4, -4, -4, -11, -4, -4, -4, -11, 11, -4, -4, 4, 11, -4, 4, -4, -11,
    ))

"""
    opensimplex_3d(; seed=nothing, smooth=false)

Construct a sampler that outputs 3-dimensional legacy OpenSimplex noise when it is sampled from.

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

# See also:

  - [`opensimplex2_3d`](@ref opensimplex2_3d)
  - [`opensimplex2s_3d`](@ref opensimplex2s_3d)
"""
opensimplex_3d(; seed=nothing, smooth=false) = _opensimplex(3, seed, smooth)

@inline function contribute(sampler::OpenSimplex{3}, falloff, X, Y, Z, x, y, z)
    t = sampler.table
    g = OS_GRADIENTS_3D
    @inbounds i = t[(Z+t[(Y+t[(X&0xff+1)])&0xff+1])&0xff+1] + 1
    a = falloff - x^2 - y^2 - z^2
    @inbounds a > 0 ? pow4(a) * (g[i] * x + g[i+1] * y + g[i+2] * z) : 0.0
end

function sample(sampler::OpenSimplex{3}, x::T, y::T, z::T) where {T<:Real}
    s = (x, y, z) .+ ((x + y + z) * OS_STRETCH_3D)
    state = sampler.simplex_state
    falloff = state.falloff
    X0, Y0, Z0 = floor.(Int, s)
    X1, Y1, Z1 = (X0, Y0, Z0) .& 255 .+ 1
    x1, y1, z1 = (x, y, z) .- ((X0 + Y0 + Z0) * OS_SQUISH_3D) .- (X0, Y0, Z0)
    Xs, Ys, Zs = s .- (X0, Y0, Z0)
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
                x8 = x1 - 2OS_SQUISH_3D
                x9 = x1 + 1 - OS_SQUISH_3D
            else
                X = X1 + 1
                X2, X3 = X, X
                x8 = x1 - 1 - 2OS_SQUISH_3D
                x9 = x1 - 1 - OS_SQUISH_3D
            end
            if iszero(c & 2)
                Y2, Y3 = Y1, Y1 - 1
                y8 = y1 - 2OS_SQUISH_3D
                y9 = y1 + 1 - OS_SQUISH_3D
            else
                Y = Y1 + 1
                Y2, Y3 = Y, Y
                y8 = y1 - 1 - 2OS_SQUISH_3D
                y9 = y1 - 1 - OS_SQUISH_3D
            end
            if iszero(c & 4)
                Z2, Z3 = Z1, Z1 - 1
                z8 = z1 - 2OS_SQUISH_3D
                z9 = z1 + 1 - OS_SQUISH_3D
            else
                Z = Z1 + 1
                Z2, Z3 = Z, Z
                z8 = z1 - 1 - 2OS_SQUISH_3D
                z9 = z1 - 1 - OS_SQUISH_3D
            end
        end
        y2, z2, x3 = (y1, z1, x1) .- OS_SQUISH_3D
        x2, y3, z4 = (x1, y1, z1) .- 1 .- OS_SQUISH_3D
        z3, x4, y4 = z2, x3, y2
        c1 = contribute(sampler, falloff, X1, Y1, Z1, x1, y1, z1)
        c2 = contribute(sampler, falloff, X1 + 1, Y1, Z1, x2, y2, z2)
        c3 = contribute(sampler, falloff, X1, Y1 + 1, Z1, x3, y3, z3)
        c4 = contribute(sampler, falloff, X1, Y1, Z1 + 1, x4, y4, z4)
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
                x8 = x1 - 2 - 3OS_SQUISH_3D
                x9 = x1 - 1 - 3OS_SQUISH_3D
            else
                x = x1 - 3OS_SQUISH_3D
                X2, X3, x8, x9 = X1, X1, x, x
            end
            if !iszero(c & 2)
                Y = Y1 + 1
                y = y1 - 1 - 3OS_SQUISH_3D
                Y2, Y3, y8, y9 = Y, Y, y, y
                if !iszero(c & 1)
                    Y3 += 1
                    y9 -= 1
                else
                    Y2 += 1
                    y8 -= 1
                end
            else
                y = y1 - 3OS_SQUISH_3D
                Y2, Y3, y8, y9 = Y1, Y1, y, y
            end
            if !iszero(c & 4)
                Z2, Z3 = Z1 + 1, Z1 + 2
                z8 = z1 - 1 - 3OS_SQUISH_3D
                z9 = z1 - 2 - 3OS_SQUISH_3D
            else
                z = z1 - 3OS_SQUISH_3D
                Z2, Z3, z8, z9 = Z1, Z1, z, z
            end
        else
            c = (p1 & p2) & 0xff
            if !iszero(c & 1)
                X2, X3 = X1 + 1, X1 + 2
                x8 = x1 - 1 - OS_SQUISH_3D
                x9 = x1 - 2 - 2OS_SQUISH_3D
            else
                X2, X3 = X1, X1
                x8 = x1 - OS_SQUISH_3D
                x9 = x1 - 2OS_SQUISH_3D
            end
            if !iszero(c & 2)
                Y2, Y3 = Y1 + 1, Y1 + 2
                y8 = y1 - 1 - OS_SQUISH_3D
                y9 = y1 - 2 - 2OS_SQUISH_3D
            else
                Y2, Y3 = Y1, Y1
                y8 = y1 - OS_SQUISH_3D
                y9 = y1 - 2OS_SQUISH_3D
            end
            if !iszero(c & 4)
                Z2, Z3 = Z1 + 1, Z1 + 2
                z8 = z1 - 1 - OS_SQUISH_3D
                z9 = z1 - 2 - 2OS_SQUISH_3D
            else
                Z2, Z3 = Z1, Z1
                z8 = z1 - OS_SQUISH_3D
                z9 = z1 - 2OS_SQUISH_3D
            end
        end
        z3, x4, y4 = (z1, x1, y1) .- 1 .- 2OS_SQUISH_3D
        x3, y2, z2 = x4, y4, z3
        z4, y3, x2 = (z1, y1, x1) .- 2OS_SQUISH_3D
        x1, y1, z1 = (x1, y1, z1) .- 1 .- 3OS_SQUISH_3D
        c1 = contribute(sampler, falloff, X1 + 1, Y1 + 1, Z1, x4, y4, z4)
        c2 = contribute(sampler, falloff, X1 + 1, Y1, Z1 + 1, x3, y3, z3)
        c3 = contribute(sampler, falloff, X1, Y1 + 1, Z1 + 1, x2, y2, z2)
        c4 = contribute(sampler, falloff, X1 + 1, Y1 + 1, Z1 + 1, x1, y1, z1)
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
                x8, y8, z8 = (x1, y1, z1) .- 1 .- 3OS_SQUISH_3D
                X2, Y2, Z2 = Z1 + 1, Y1 + 1, Z1 + 1
                if !iszero(c & 1)
                    x9 = x1 - 2 - 2OS_SQUISH_3D
                    y9, z9 = (y1, z1) .- 2OS_SQUISH_3D
                    X3, Y3, Z3 = X1 + 2, Y1, Z1
                elseif !iszero(c & 2)
                    x9, z9 = (x1, z1) .- 2OS_SQUISH_3D
                    y9 = y1 - 2 - 2OS_SQUISH_3D
                    X3, Y3, Z3 = X1, Y1 + 2, Z1
                else
                    x9, y9 = (x1, y1) .- 2OS_SQUISH_3D
                    z9 = z1 - 2 - 2OS_SQUISH_3D
                    X3, Y3, Z3 = X1, Y1, Z1 + 2
                end
            else
                c = p1 | p2
                x8, y8, z8, X2, Y2, Z2 = x1, y1, z1, X1, Y1, Z1
                if iszero(c & 1)
                    x9 = x1 + 1 - OS_SQUISH_3D
                    y9, z9 = (y1, z1) .- 1 .- OS_SQUISH_3D
                    X3, Y3, Z3 = X1 - 1, Y1 + 1, Z1 + 1
                elseif iszero(c & 2)
                    x9, z9 = (x1, z1) .- 1 .- OS_SQUISH_3D
                    y9 = y1 + 1 - OS_SQUISH_3D
                    X3, Y3, Z3 = X1 + 1, Y1 - 1, Z1 + 1
                else
                    x9, y9 = (x1, y1) .- 1 .- OS_SQUISH_3D
                    z9 = z1 + 1 - OS_SQUISH_3D
                    X3, Y3, Z3 = X1 + 1, Y1 + 1, Z1 - 1
                end
            end
        else
            c1 = p1_farthest ? p1 : p2
            c2 = p1_farthest ? p2 : p1
            if iszero(c1 & 1)
                x8 = x1 + 1 - OS_SQUISH_3D
                y8, z8 = (y1, z1) .- 1 .- OS_SQUISH_3D
                X2, Y2, Z2 = X1 - 1, Y1 + 1, Z1 + 1
            elseif iszero(c1 & 2)
                x8, z8 = (x1, z1) .- 1 .- OS_SQUISH_3D
                y8 = y1 + 1 - OS_SQUISH_3D
                X2, Y2, Z2 = X1 + 1, Y1 - 1, Z1 + 1
            else
                x8, y8 = (x1, y1) .- 1 .- OS_SQUISH_3D
                z8 = z1 + 1 - OS_SQUISH_3D
                X2, Y2, Z2 = X1 + 1, Y1 + 1, Z1 - 1
            end
            x9, y9, z9 = (x1, y1, z1) .- 2OS_SQUISH_3D
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
        y2, z2, x3 = (y1, z1, x1) .- OS_SQUISH_3D
        x5, y5, z6 = (x1, y1, z1) .- 1 .- 2OS_SQUISH_3D
        x2, y3, z4 = (x1, y1, z1) .- 1 .- OS_SQUISH_3D
        z5, y6, x7 = (z1, y1, x1) .- 2OS_SQUISH_3D
        z3, x4, y4, x6, y7, z7 = z2, x3, y2, x5, y5, z6
        c1 = contribute(sampler, falloff, X1 + 1, Y1, Z1, x2, y2, z2)
        c2 = contribute(sampler, falloff, X1, Y1 + 1, Z1, x3, y3, z3)
        c3 = contribute(sampler, falloff, X1, Y1, Z1 + 1, x4, y4, z4)
        c4 = contribute(sampler, falloff, X1 + 1, Y1 + 1, Z1, x5, y5, z5)
        c5 = contribute(sampler, falloff, X1 + 1, Y1, Z1 + 1, x6, y6, z6)
        c6 = contribute(sampler, falloff, X1, Y1 + 1, Z1 + 1, x7, y7, z7)
        result += c1 + c2 + c3 + c4 + c5 + c6
    end
    c1 = contribute(sampler, falloff, X2, Y2, Z2, x8, y8, z8)
    c2 = contribute(sampler, falloff, X3, Y3, Z3, x9, y9, z9)
    (result + c1 + c2) * state.scale_factor
end

# 4D

const OS_STRETCH_4D = (1 / sqrt(5) - 1) / 4
const OS_SQUISH_4D = (sqrt(5) - 1) / 4
const OS_GRADIENTS_4D = [
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

"""
    opensimplex_4d(; seed=nothing, smooth=false)

Construct a sampler that outputs 4-dimensional legacy OpenSimplex noise when it is sampled from.

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

# See also:

  - [`opensimplex2_4d`](@ref opensimplex2_4d)
  - [`opensimplex2s_4d`](@ref opensimplex2s_4d)
"""
opensimplex_4d(; seed=nothing, smooth=false) = _opensimplex(4, seed, smooth)

@inline function contribute(sampler::OpenSimplex{4}, falloff, X, Y, Z, W, x, y, z, w)
    t = sampler.perlin_state.table
    g = OS_GRADIENTS_4D
    @inbounds i = t[(W+t[(Z+t[(Y+t[(X&0xff+1)])&0xff+1])&0xff+1]&0xff+1)] & 252 + 1
    a = falloff - x^2 - y^2 - z^2 - w^2
    @inbounds a > 0 ? pow4(a) * (g[i] * x + g[i+1] * y + g[i+2] * z + g[i+3] * w) : 0.0
end

function sample(sampler::OpenSimplex{4}, x::T, y::T, z::T, w::T) where {T<:Real}
    s = (x, y, z, w) .+ ((x + y + z + w) * OS_STRETCH_4D)
    state = sampler.simplex_state
    falloff = state.falloff
    X0, Y0, Z0, W0 = floor.(Int, s)
    X1, Y1, Z1, W1 = (X0, Y0, Z0, W0) .& 255 .+ 1
    x1, y1, z1, w1 = (x, y, z, w) .- ((X0 + Y0 + Z0 + W0) * OS_SQUISH_4D) .- (X0, Y0, Z0, W0)
    Xs, Ys, Zs, Ws = s .- (X0, Y0, Z0, W0)
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
                x12 = x1 - 2OS_SQUISH_4D
                x13 = x1 + 1 - OS_SQUISH_4D
                x14 = x1 - OS_SQUISH_4D
            else
                X, x = X1 + 1, x1 - 1 - OS_SQUISH_4D
                X2, X3, X4, x13, x14 = X, X, X, x, x
                x12 = x1 - 1 - 2OS_SQUISH_4D
            end
            if iszero(c & 2)
                _y = y1 - OS_SQUISH_4D
                Y2, Y3, Y4, y13, y14 = Y1, Y1, Y1, _y, _y
                y12 = y1 - 2OS_SQUISH_4D
                if (c & 1) == 1
                    Y3 -= 1
                    y13 += 1
                else
                    Y4 -= 1
                    y14 += 1
                end
            else
                _Y1 = Y1 + 1
                _y = y1 - 1 - OS_SQUISH_4D
                Y2, Y3, Y4, y13, y14 = _Y1, _Y1, _Y1, _y, _y
                y12 = y1 - 1 - 2OS_SQUISH_4D
            end
            if iszero(c & 4)
                _z = z1 - OS_SQUISH_4D
                Z12, Z13, Z14, z13, z14 = Z1, Z1, Z1, _z, _z
                z12 = z1 - 2OS_SQUISH_4D
                if (c & 3) == 3
                    Z13 -= 1
                    z13 += 1
                else
                    Z14 -= 1
                    z14 += 1
                end
            else
                _Z1 = Z1 + 1
                _z = z1 - 1 - OS_SQUISH_4D
                Z12, Z13, Z14, z13, z14 = _Z1, _Z1, _Z1, _z, _z
                z12 = z1 - 1 - 2OS_SQUISH_4D
            end
            if iszero(c & 8)
                W12, W13, W14 = W1, W1, W1 - 1
                w12 = w1 - 2OS_SQUISH_4D
                w13 = w1 - OS_SQUISH_4D
                w14 = w1 + 1 - OS_SQUISH_4D
            else
                _W1 = W1 + 1
                _w = w1 - 1 - OS_SQUISH_4D
                W12, W13, W14, w13, w14 = _W1, _W1, _W1, _w, _w
                w12 = w1 - 1 - 2OS_SQUISH_4D
            end
        end
        y2, z2, w2, x3 = (y1, z1, w1, x1) .- OS_SQUISH_4D
        x2, y3, z4, w5 = (x1, y1, z1, w1) .- 1 .- OS_SQUISH_4D
        z3, w3, x4, y4, w4, x5, y5, z5 = (z2, w2, x3, y2, w2, x3, y2, z2)
        c1 = contribute(sampler, falloff, X1, Y1, Z1, W1, x1, y1, z1, w1)
        c2 = contribute(sampler, falloff, X1 + 1, Y1, Z1, W1, x2, y2, z2, w2)
        c3 = contribute(sampler, falloff, X1, Y1 + 1, Z1, W1, x3, y3, z3, w3)
        c4 = contribute(sampler, falloff, X1, Y1, Z1 + 1, W1, x4, y4, z4, w4)
        c5 = contribute(sampler, falloff, X1, Y1, Z1, W1 + 1, x5, y5, z5, w5)
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
                x = x1 - 1 - 4OS_SQUISH_4D
                X2 = X1 + 2
                x12 = x1 - 2 - 4OS_SQUISH_4D
                X3, X4, x13, x14 = X, X, x, x
            else
                x = x1 - 4OS_SQUISH_4D
                X2, X3, X4, x12, x13, x14 = X1, X1, X1, x, x, x
            end
            if !iszero(c & 2)
                _Y1 = Y1 + 1
                _y = y1 - 1 - 4OS_SQUISH_4D
                Y2, Y3, Y4, y12, y13, y14 = _Y1, _Y1, _Y1, _y, _y, _y
                if !iszero(c & 1)
                    Y3 += 1
                    y13 -= 1
                else
                    Y2 += 1
                    y12 -= 1
                end
            else
                _y = y1 - 4OS_SQUISH_4D
                Y2, Y3, Y4, y12, y13, y14 = Y1, Y1, Y1, _y, _y, _y
            end
            if !iszero(c & 4)
                _Z1 = Z1 + 1
                _z = z1 - 1 - 4OS_SQUISH_4D
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
                _z = z1 - 4OS_SQUISH_4D
                Z12, Z13, Z14, z12, z13, z14 = Z1, Z1, Z1, _z, _z, _z
            end
            if !iszero(c & 8)
                _W1 = W1 + 1
                _w = w1 - 1 - 4OS_SQUISH_4D
                W12, W13, w12, w13 = _W1, _W1, _w, _w
                W14 = W1 + 2
                w14 = w1 - 2 - 4OS_SQUISH_4D
            else
                _w = w1 - 4OS_SQUISH_4D
                W12, W13, W14, w12, w13, w14 = W1, W1, W1, _w, _w, _w
            end
        else
            c = p1 & p2
            if !iszero(c & 1)
                X = X1 + 1
                X2, X3, X4 = X, X1 + 2, X
                x12 = x1 - 1 - 2OS_SQUISH_4D
                x13 = x1 - 2 - 3OS_SQUISH_4D
                x14 = x1 - 1 - 3OS_SQUISH_4D
            else
                x = x1 - 3OS_SQUISH_4D
                X2, X3, X4, x13, x14 = X1, X1, X1, x, x
                x12 = x1 - 2OS_SQUISH_4D
            end
            if !iszero(c & 2)
                _Y1 = Y1 + 1
                _y = y1 - 1 - 3OS_SQUISH_4D
                Y2, Y3, Y4, y13, y14 = _Y1, _Y1, _Y1, _y, _y
                y12 = y1 - 1 - 2OS_SQUISH_4D
                if !iszero(c & 1)
                    Y4 += 1
                    y14 -= 1
                else
                    Y3 += 1
                    y13 -= 1
                end
            else
                _y = y1 - 3OS_SQUISH_4D
                Y2, Y3, Y4, y13, y14 = Y1, Y1, Y1, _y, _y
                y12 = y1 - 2OS_SQUISH_4D
            end
            if !iszero(c & 4)
                _Z1 = Z1 + 1
                _z = z1 - 1 - 3OS_SQUISH_4D
                Z12, Z13, Z14, z13, z14 = _Z1, _Z1, _Z1, _z, _z
                z12 = z1 - 1 - 2OS_SQUISH_4D
                if !iszero(c & 3)
                    Z14 += 1
                    z14 -= 1
                else
                    Z13 += 1
                    z13 -= 1
                end
            else
                _z = z1 - 3OS_SQUISH_4D
                Z12, Z13, Z14, z13, z14 = Z1, Z1, Z1, _z, _z
                z12 = z1 - 2OS_SQUISH_4D
            end
            if !iszero(c & 8)
                _W1 = W1 + 1
                W12, W13, W14 = _W1, _W1, W1 + 2
                w12 = w1 - 1 - 2OS_SQUISH_4D
                w13 = w1 - 1 - 3OS_SQUISH_4D
                w14 = w1 - 2 - 3OS_SQUISH_4D
            else
                _w = w1 - 3OS_SQUISH_4D
                W12, W13, W14, w13, w14 = W1, W1, W1, _w, _w
                w12 = w1 - 2OS_SQUISH_4D
            end
        end
        w4, x5, y5, z5 = (w1, x1, y1, z1) .- 1 .- 3OS_SQUISH_4D
        w5, z4, y3, x2 = (w1, z1, y1, x1) .- 3OS_SQUISH_4D
        x4, y4, x3, z3, w3, y2, z2, w2 = x5, y5, x5, z5, w4, y5, z5, w4
        x1, y1, z1, w1 = (x1, y1, z1, w1) .- 1 .- 4OS_SQUISH_4D
        c1 = contribute(sampler, falloff, X1 + 1, Y1 + 1, Z1 + 1, W1, x5, y5, z5, w5)
        c2 = contribute(sampler, falloff, X1 + 1, Y1 + 1, Z1, W1 + 1, x4, y4, z4, w4)
        c3 = contribute(sampler, falloff, X1 + 1, Y1, Z1 + 1, W1 + 1, x3, y3, z3, w3)
        c4 = contribute(sampler, falloff, X1, Y1 + 1, Z1 + 1, W1 + 1, x2, y2, z2, w2)
        c5 = contribute(sampler, falloff, X1 + 1, Y1 + 1, Z1 + 1, W1 + 1, x1, y1, z1, w1)
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
                    x12 = x1 - 3OS_SQUISH_4D
                    x13 = x1 + 1 - 2OS_SQUISH_4D
                else
                    X = X1 + 1
                    X2, X3 = X, X
                    x12 = x1 - 1 - 3OS_SQUISH_4D
                    x13 = x1 - 1 - 2OS_SQUISH_4D
                end
                if iszero(c1 & 2)
                    Y2, Y3 = Y1, Y1 - 1
                    y12 = y1 - 3OS_SQUISH_4D
                    y13 = y1 + 1 - 2OS_SQUISH_4D
                else
                    _Y1 = Y1 + 1
                    Y2, Y3 = _Y1, _Y1
                    y12 = y1 - 1 - 3OS_SQUISH_4D
                    y13 = y1 - 1 - 2OS_SQUISH_4D
                end
                if iszero(c1 & 4)
                    Z12, Z13 = Z1, Z1 - 1
                    z12 = z1 - 3OS_SQUISH_4D
                    z13 = z1 + 1 - 2OS_SQUISH_4D
                else
                    _Z1 = Z1 + 1
                    Z12, Z13 = _Z1, _Z1
                    z12 = z1 - 1 - 3OS_SQUISH_4D
                    z13 = z1 - 1 - 2OS_SQUISH_4D
                end
                if iszero(c1 & 8)
                    W12, W13 = W1, W1 - 1
                    w12 = w1 - 3OS_SQUISH_4D
                    w13 = w1 + 1 - 2OS_SQUISH_4D
                else
                    _W1 = W1 + 1
                    W12, W13 = _W1, _W1
                    w12 = w1 - 1 - 3OS_SQUISH_4D
                    w13 = w1 - 1 - 2OS_SQUISH_4D
                end
                X4, Y4, Z14, W14 = X1, Y1, Z1, W1
                x14, y14, z14, w14 = (x1, y1, z1, w1) .- 2OS_SQUISH_4D
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
                    x12 = x1 + 1 - OS_SQUISH_4D
                    x13 = x1 - OS_SQUISH_4D
                else
                    X = X1 + 1
                    x = x1 - 1 - OS_SQUISH_4D
                    X2, X3, x12, x13 = X, X, x, x
                end
                if iszero(c & 2)
                    _y = y1 - OS_SQUISH_4D
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
                    _y = y1 - 1 - OS_SQUISH_4D
                    Y2, Y3, y12, y13 = _Y1, _Y1, _y, _y
                end
                if iszero(c & 4)
                    _z = z1 - OS_SQUISH_4D
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
                    _z = z1 - 1 - OS_SQUISH_4D
                    Z12, Z13, z12, z13 = _Z1, _Z1, _z, _z
                end
                if iszero(c & 8)
                    W12, W13 = W1, W1 - 1
                    w12 = w1 - OS_SQUISH_4D
                    w13 = w1 + 1 - OS_SQUISH_4D
                else
                    _W1 = W1 + 1
                    _w = w1 - 1 - OS_SQUISH_4D
                    W12, W13, w12, w13 = _W1, _W1, _w, _w
                end
            end
        else
            c1 = p1_bigger === true ? p1 : p2
            c2 = p1_bigger === true ? p2 : p1
            if iszero(c1 & 1)
                X2, X3 = X1 - 1, X1
                x12 = x1 + 1 - OS_SQUISH_4D
                x13 = x1 - OS_SQUISH_4D
            else
                X = X1 + 1
                x = x1 - 1 - OS_SQUISH_4D
                X2, X3, x12, x13 = X, X, x, x
            end
            if iszero(c1 & 2)
                _y = y1 - OS_SQUISH_4D
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
                _y = y1 - 1 - OS_SQUISH_4D
                Y2, Y3, y12, y13 = _Y1, _Y1, _y, _y
            end
            if iszero(c1 & 4)
                _z = z1 - OS_SQUISH_4D
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
                _z = z1 - 1 - OS_SQUISH_4D
                Z12, Z13, z12, z13 = _Z1, _Z1, _z, _z
            end
            if iszero(c1 & 8)
                W12, W13 = W1, W1 - 1
                w12 = w1 - OS_SQUISH_4D
                w13 = w1 + 1 - OS_SQUISH_4D
            else
                _W1 = W1 + 1
                _w = w1 - 1 - OS_SQUISH_4D
                W12, W13, w12, w13 = _W1, _W1, _w, _w
            end
            X4, Y4, Z14, W14 = X1, Y1, Z1, W1
            x14, y14, z14, w14 = (x1, y1, z1, w1) .- 2OS_SQUISH_4D
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
        y2, z2, w2, x3 = (y1, z1, w1, x1) .- OS_SQUISH_4D
        z6, w6, y7, x9 = (z1, w1, y1, x1) .- 2OS_SQUISH_4D
        x6, y6, z7, w8 = (x1, y1, z1, w1) .- 1 .- 2OS_SQUISH_4D
        x2, y3, z4, w5 = (x1, y1, z1, w1) .- 1 .- OS_SQUISH_4D
        z3, w3, x4, y4, w4, x5, y5, z5 = z2, w2, x3, y2, w2, x3, y2, z2
        x7, w7, x8, y8, z8, y9, z9, w9 = x6, w6, x6, y7, z6, y6, z7, w6
        x10, y10, z10, w10, x11, y11, z11, w11 = x9, y6, z6, w8, x9, y7, z7, w8
        c1 = contribute(sampler, falloff, X1 + 1, Y1, Z1, W1, x2, y2, z2, w2)
        c2 = contribute(sampler, falloff, X1, Y1 + 1, Z1, W1, x3, y3, z3, w3)
        c3 = contribute(sampler, falloff, X1, Y1, Z1 + 1, W1, x4, y4, z4, w4)
        c4 = contribute(sampler, falloff, X1, Y1, Z1, W1 + 1, x5, y5, z5, w5)
        c5 = contribute(sampler, falloff, X1 + 1, Y1 + 1, Z1, W1, x6, y6, z6, w6)
        c6 = contribute(sampler, falloff, X1 + 1, Y1, Z1 + 1, W1, x7, y7, z7, w7)
        c7 = contribute(sampler, falloff, X1 + 1, Y1, Z1, W1 + 1, x8, y8, z8, w8)
        c8 = contribute(sampler, falloff, X1, Y1 + 1, Z1 + 1, W1, x9, y9, z9, w9)
        c9 = contribute(sampler, falloff, X1, Y1 + 1, Z1, W1 + 1, x10, y10, z10, w10)
        c10 = contribute(sampler, falloff, X1, Y1, Z1 + 1, W1 + 1, x11, y11, z11, w11)
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
                x12, y12, z12, w12 = (x1, y1, z1, w1) .- OS_SQUISH_4D
                x13, y13, z13, w13 = (x1, y1, z1, w1) .- 2OS_SQUISH_4D
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
                x14, y14, z14, w14 = (x1, y1, z1, w1) .- 1 .- 2OS_SQUISH_4D
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
                x14, y14, z14, w14 = (x1, y1, z1, w1) .- 1 .- 4OS_SQUISH_4D
                if !iszero(c & 1)
                    X2, X3 = X1 + 2, X1 + 1
                    x12 = x1 - 2 - 3OS_SQUISH_4D
                    x13 = x1 - 1 - 3OS_SQUISH_4D
                else
                    x = x1 - 3OS_SQUISH_4D
                    X2, X3, x12, x13 = X1, X1, x, x
                end
                if !iszero(c & 2)
                    _Y1 = Y1 + 1
                    _y = y1 - 1 - 3OS_SQUISH_4D
                    Y2, Y3, y12, y13 = _Y1, _Y1, _y, _y
                    if iszero(c & 1)
                        Y2 += 1
                        y12 -= 1
                    else
                        Y3 += 1
                        y13 -= 1
                    end
                else
                    _y = y1 - 3OS_SQUISH_4D
                    Y2, Y3, y12, y13 = Y1, Y1, _y, _y
                end
                if !iszero(c & 4)
                    _Z1 = Z1 + 1
                    _z = z1 - 1 - 3OS_SQUISH_4D
                    Z12, Z13, z12, z13 = _Z1, _Z1, _z, _z
                    if iszero(c & 3)
                        Z12 += 1
                        z12 -= 1
                    else
                        Z13 += 1
                        z13 -= 1
                    end
                else
                    _z = z1 - 3OS_SQUISH_4D
                    Z12, Z13, z12, z13 = Z1, Z1, _z, _z
                end
                if !iszero(c & 8)
                    W12, W13 = W1 + 1, W1 + 2
                    w12 = w1 - 1 - 3OS_SQUISH_4D
                    w13 = w1 - 2 - 3OS_SQUISH_4D
                else
                    _w = w1 - 3OS_SQUISH_4D
                    W12, W13, w12, w13 = W1, W1, _w, _w
                end
            end
        else
            c1 = p1_bigger === true ? p1 : p2
            c2 = p1_bigger === true ? p2 : p1
            if !iszero(c1 & 1)
                X2, X3 = X1 + 2, X1 + 1
                x12 = x1 - 2 - 3OS_SQUISH_4D
                x13 = x1 - 1 - 3OS_SQUISH_4D
            else
                x = x1 - 3OS_SQUISH_4D
                X2, X3, x12, x13 = X1, X1, x, x
            end
            if !iszero(c1 & 2)
                _Y1 = Y1 + 1
                _y = y1 - 1 - 3OS_SQUISH_4D
                Y2, Y3, y12, y13 = _Y1, _Y1, _y, _y
                if iszero(c1 & 1)
                    Y2 += 1
                    y12 -= 1
                else
                    Y3 += 1
                    y13 -= 1
                end
            else
                _y = y1 - 3OS_SQUISH_4D
                Y2, Y3, y12, y13 = Y1, Y1, _y, _y
            end
            if !iszero(c1 & 4)
                _Z1 = Z1 + 1
                _z = z1 - 1 - 3OS_SQUISH_4D
                Z12, Z13, z12, z13 = _Z1, _Z1, _z, _z
                if iszero(c1 & 3)
                    Z12 += 1
                    z12 -= 1
                else
                    Z13 += 1
                    z13 -= 1
                end
            else
                _z = z1 - 3OS_SQUISH_4D
                Z12, Z13, z12, z13 = Z1, Z1, _z, _z
            end
            if !iszero(c1 & 8)
                W12, W13 = W1 + 1, W1 + 2
                w12 = w1 - 1 - 3OS_SQUISH_4D
                w13 = w1 - 2 - 3OS_SQUISH_4D
            else
                _w = w1 - 3OS_SQUISH_4D
                W12, W13, w12, w13 = W1, W1, _w, _w
            end
            X4, Y4, Z14, W14 = X1 + 1, Y1 + 1, Z1 + 1, W1 + 1
            x14, y14, z14, w14 = (x1, y1, z1, w1) .- 1 .- 2OS_SQUISH_4D
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
        w4, x5, y5, z5 = (w1, x1, y1, z1) .- 1 .- 3OS_SQUISH_4D
        x6, y6, z7, w8 = (x1, y1, z1, w1) .- 1 .- 2OS_SQUISH_4D
        z6, w6, y7, x9 = (z1, w1, y1, x1) .- 2OS_SQUISH_4D
        w5, z4, y3, x2 = (w1, z1, y1, x1) .- 3OS_SQUISH_4D
        x4, y4, x3, z3, w3, y2, z2, w2 = x5, y5, x5, z5, w4, y5, z5, w4
        x7, w7, x8, y8, z8, y9, z9, w9 = x6, w6, x6, y7, z6, y6, z7, w6
        x10, y10, z10, w10, x11, y11, z11, w11 = x9, y6, z6, w8, x9, y7, z7, w8
        c1 = contribute(sampler, falloff, X1 + 1, Y1 + 1, Z1 + 1, W1, x5, y5, z5, w5)
        c2 = contribute(sampler, falloff, X1 + 1, Y1 + 1, Z1, W1 + 1, x4, y4, z4, w4)
        c3 = contribute(sampler, falloff, X1 + 1, Y1, Z1 + 1, W1 + 1, x3, y3, z3, w3)
        c4 = contribute(sampler, falloff, X1, Y1 + 1, Z1 + 1, W1 + 1, x2, y2, z2, w2)
        c5 = contribute(sampler, falloff, X1 + 1, Y1 + 1, Z1, W1, x6, y6, z6, w6)
        c6 = contribute(sampler, falloff, X1 + 1, Y1, Z1 + 1, W1, x7, y7, z7, w7)
        c7 = contribute(sampler, falloff, X1 + 1, Y1, Z1, W1 + 1, x8, y8, z8, w8)
        c8 = contribute(sampler, falloff, X1, Y1 + 1, Z1 + 1, W1, x9, y9, z9, w9)
        c9 = contribute(sampler, falloff, X1, Y1 + 1, Z1, W1 + 1, x10, y10, z10, w10)
        c10 = contribute(sampler, falloff, X1, Y1, Z1 + 1, W1 + 1, x11, y11, z11, w11)
        result += c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10
    end
    c1 = contribute(sampler, falloff, X2, Y2, Z12, W12, x12, y12, z12, w12)
    c2 = contribute(sampler, falloff, X3, Y3, Z13, W13, x13, y13, z13, w13)
    c3 = contribute(sampler, falloff, X4, Y4, Z14, W14, x14, y14, z14, w14)
    (result + c1 + c2 + c3) * state.scale_factor
end
