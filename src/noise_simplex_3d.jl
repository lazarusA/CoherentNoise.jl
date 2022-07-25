const SKEW_3D = 1 / 3

const UNSKEW_3D = 1 / 6

const SCALE_3D = 32

@inline function grad(S::Type{Simplex{3}}, hash, x, y, z)
    s = 0.6 - x^2 - y^2 - z^2
    h = hash & 15
    u = h < 8 ? x : y
    v = h < 4 ? y : h == 12 || h == 14 ? x : z
    g = hash_coords(S, h, u, v)
    @fastpow s > 0 ? s^4 * g : 0.0
end

@inline function get_simplex(::Type{Simplex{3}}, x, y, z)
    if x ≥ y
        y ≥ z ? (1, 0, 0, 1, 1, 0) : x ≥ z ? (1, 0, 0, 1, 0, 1) : (0, 0, 1, 1, 0, 1)
    else
        y < z ? (0, 0, 1, 0, 1, 1) : x < z ? (0, 1, 0, 0, 1, 1) : (0, 1, 0, 1, 1, 0)
    end
end

function sample(sampler::S, x::T, y::T, z::T) where {S<:Simplex{3},T<:Real}
    t = sampler.state.table
    s = (x + y + z) * SKEW_3D
    X, Y, Z = floor.(Int, (x, y, z) .+ s)
    tx = (X + Y + Z) * UNSKEW_3D
    xyz = (x, y, z) .- (X, Y, Z) .+ tx
    X1, Y1, Z1, X2, Y2, Z2 = get_simplex(S, xyz...)
    p1 = grad(S, t[t[t[Z]+Y]+X], xyz...)
    p2 = grad(S, t[t[t[Z+Z1]+Y+Y1]+X+X1], xyz .- (X1, Y1, Z1) .+ UNSKEW_3D...)
    p3 = grad(S, t[t[t[Z+Z2]+Y+Y2]+X+X2], xyz .- (X2, Y2, Z2) .+ 2UNSKEW_3D...)
    p4 = grad(S, t[t[t[Z+1]+Y+1]+X+1], xyz .- 1 .+ 3UNSKEW_3D...)
    (p1 + p2 + p3 + p4) * SCALE_3D
end
