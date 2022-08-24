abstract type Orientation end
struct OrientStandard <: Orientation end
struct OrientX <: Orientation end
struct OrientXY <: Orientation end
struct OrientXZ <: Orientation end
struct OrientXYZ <: Orientation end

struct OpenSimplex2{N,O<:Orientation} <: NoiseSampler{N}
    random_state::RandomState
end

@inline function _opensimplex2(dims, seed, orientation)
    rs = RandomState(seed)
    orientation = os2_orientation_type(Val(orientation))
    OpenSimplex2{dims,orientation}(rs)
end

@inline os2_orientation_type(::Val{nothing}) = OrientStandard
@inline os2_orientation_type(::Val{:x}) = OrientX
@inline os2_orientation_type(::Val{:xy}) = OrientXY
@inline os2_orientation_type(::Val{:xz}) = OrientXZ
@inline os2_orientation_type(::Val{:xyz}) = OrientXYZ

# 2D

const OS2_SKEW_2D = 0.366025403784439
const OS2_UNSKEW_2D = -0.21132486540518713
const OS2_R²2D = 0.5
const OS2_NUM_GRADIENTS_EXP_2D = 7
const OS2_NUM_GRADIENTS_2D = 1 << OS2_NUM_GRADIENTS_EXP_2D
const OS2_GRADIENTS_NORMALIZED_2D = [
    0.38268343236509, 0.923879532511287, 0.923879532511287, 0.38268343236509,
    0.923879532511287, -0.38268343236509, 0.38268343236509, -0.923879532511287,
    -0.38268343236509, -0.923879532511287, -0.923879532511287, -0.38268343236509,
    -0.923879532511287, 0.38268343236509, -0.38268343236509, 0.923879532511287,
    0.130526192220052, 0.99144486137381, 0.608761429008721, 0.793353340291235,
    0.793353340291235, 0.608761429008721, 0.99144486137381, 0.130526192220051,
    0.99144486137381, -0.130526192220051, 0.793353340291235, -0.60876142900872,
    0.608761429008721, -0.793353340291235, 0.130526192220052, -0.99144486137381,
    -0.130526192220052, -0.99144486137381, -0.608761429008721, -0.793353340291235,
    -0.793353340291235, -0.608761429008721, -0.99144486137381, -0.130526192220052,
    -0.99144486137381, 0.130526192220051, -0.793353340291235, 0.608761429008721,
    -0.608761429008721, 0.793353340291235, -0.130526192220052, 0.99144486137381]
const OS2_GRADIENTS_2D = OS2_GRADIENTS_NORMALIZED_2D ./ 0.01001634121365712 |> CircularVector

@doc doc_opensimplex2_2d
opensimplex2_2d(; seed=0, orient=nothing) = _opensimplex2(2, seed, orient)

@inline function grad(table, seed, X, Y, x, y)
    hash = (seed ⊻ X ⊻ Y) * HASH_MULTIPLIER
    hash ⊻= hash >> (64 - OS2_NUM_GRADIENTS_EXP_2D + 1)
    i = trunc(hash) & ((OS2_NUM_GRADIENTS_2D - 1) << 1)
    t = (table[i+1], table[(i|1)+1])
    sum((t .* (x, y)))
end

@inline orient(::OpenSimplex2{2,OrientStandard}, x, y) = (x, y) .+ OS2_SKEW_2D .* (x + y)

@inline function orient(::OpenSimplex2{2,OrientX}, x, y)
    xx = x * ROOT_2_OVER_2
    yy = y * ROOT_2_OVER_2 * (2OS2_SKEW_2D + 1)
    (yy + xx, yy - xx)
end

function sample(sampler::OpenSimplex2{2,O}, x::T, y::T) where {O,T<:Real}
    seed = sampler.random_state.seed
    primes = (PRIME_X, PRIME_Y)
    tr = orient(sampler, x, y)
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
        result += pow4(a1) * grad(OS2_GRADIENTS_2D, seed, X1, Y1, x1, y1)
    end
    a2 = 2us1 * (1 / OS2_UNSKEW_2D + 2) * t + -2us1^2 + a1
    if a2 > 0
        x, y = (x1, y1) .- 2OS2_UNSKEW_2D .- 1
        result += pow4(a2) * grad(OS2_GRADIENTS_2D, seed, X2, Y2, x, y)
    end
    if y1 > x1
        x = x1 - OS2_UNSKEW_2D
        y = y1 - OS2_UNSKEW_2D - 1
        a3 = OS2_R²2D - x^2 - y^2
        if a3 > 0
            result += pow4(a3) * grad(OS2_GRADIENTS_2D, seed, X1, Y2, x, y)
        end
    else
        x = x1 - OS2_UNSKEW_2D - 1
        y = y1 - OS2_UNSKEW_2D
        a4 = OS2_R²2D - x^2 - y^2
        if a4 > 0
            result += pow4(a4) * grad(OS2_GRADIENTS_2D, seed, X2, Y1, x, y)
        end
    end
    result
end

# 3D

const OS2_SEED_FLIP_3D = -0x52d547b2e96ed629
const OS2_FALLBACK_ROTATE_3D = 2 / 3
const OS2_ROTATE_3D_ORTHONORMALIZER = OS2_UNSKEW_2D
const OS2_R²3D = 0.6
const OS2_NUM_GRADIENTS_EXP_3D = 8
const OS2_NUM_GRADIENTS_3D = 1 << OS2_NUM_GRADIENTS_EXP_3D
const OS2_GRADIENTS_NORMALIZED_3D = [
    2.22474487139, 2.22474487139, -1.0, 0.0,
    2.22474487139, 2.22474487139, 1.0, 0.0,
    3.0862664687972017, 1.1721513422464978, 0.0, 0.0,
    1.1721513422464978, 3.0862664687972017, 0.0, 0.0,
    -2.22474487139, 2.22474487139, -1.0, 0.0,
    -2.22474487139, 2.22474487139, 1.0, 0.0,
    -1.1721513422464978, 3.0862664687972017, 0.0, 0.0,
    -3.0862664687972017, 1.1721513422464978, 0.0, 0.0,
    -1.0, -2.22474487139, -2.22474487139, 0.0,
    1.0, -2.22474487139, -2.22474487139, 0.0,
    0.0, -3.0862664687972017, -1.1721513422464978, 0.0,
    0.0, -1.1721513422464978, -3.0862664687972017, 0.0,
    -1.0, -2.22474487139, 2.22474487139, 0.0,
    1.0, -2.22474487139, 2.22474487139, 0.0,
    0.0, -1.1721513422464978, 3.0862664687972017, 0.0,
    0.0, -3.0862664687972017, 1.1721513422464978, 0.0,
    -2.22474487139, -2.22474487139, -1.0, 0.0,
    -2.22474487139, -2.22474487139, 1.0, 0.0,
    -3.0862664687972017, -1.1721513422464978, 0.0, 0.0,
    -1.1721513422464978, -3.0862664687972017, 0.0, 0.0,
    -2.22474487139, -1.0, -2.22474487139, 0.0,
    -2.22474487139, 1.0, -2.22474487139, 0.0,
    -1.1721513422464978, 0.0, -3.0862664687972017, 0.0,
    -3.0862664687972017, 0.0, -1.1721513422464978, 0.0,
    -2.22474487139, -1.0, 2.22474487139, 0.0,
    -2.22474487139, 1.0, 2.22474487139, 0.0,
    -3.0862664687972017, 0.0, 1.1721513422464978, 0.0,
    -1.1721513422464978, 0.0, 3.0862664687972017, 0.0,
    -1.0, 2.22474487139, -2.22474487139, 0.0,
    1.0, 2.22474487139, -2.22474487139, 0.0,
    0.0, 1.1721513422464978, -3.0862664687972017, 0.0,
    0.0, 3.0862664687972017, -1.1721513422464978, 0.0,
    -1.0, 2.22474487139, 2.22474487139, 0.0,
    1.0, 2.22474487139, 2.22474487139, 0.0,
    0.0, 3.0862664687972017, 1.1721513422464978, 0.0,
    0.0, 1.1721513422464978, 3.0862664687972017, 0.0,
    2.22474487139, -2.22474487139, -1.0, 0.0,
    2.22474487139, -2.22474487139, 1.0, 0.0,
    1.1721513422464978, -3.0862664687972017, 0.0, 0.0,
    3.0862664687972017, -1.1721513422464978, 0.0, 0.0,
    2.22474487139, -1.0, -2.22474487139, 0.0,
    2.22474487139, 1.0, -2.22474487139, 0.0,
    3.0862664687972017, 0.0, -1.1721513422464978, 0.0,
    1.1721513422464978, 0.0, -3.0862664687972017, 0.0,
    2.22474487139, -1.0, 2.22474487139, 0.0,
    2.22474487139, 1.0, 2.22474487139, 0.0,
    1.1721513422464978, 0.0, 3.0862664687972017, 0.0,
    3.0862664687972017, 0.0, 1.1721513422464978, 0.0]
const OS2_GRADIENTS_3D = OS2_GRADIENTS_NORMALIZED_3D ./ 0.07969837668935331 |> CircularVector

@doc doc_opensimplex2_3d
opensimplex2_3d(; seed=0, orient=nothing) = _opensimplex2(3, seed, orient)

@inline function grad(table, seed, X, Y, Z, x, y, z)
    hash = ((seed ⊻ X) ⊻ (Y ⊻ Z)) * HASH_MULTIPLIER
    hash ⊻= hash >> (64 - OS2_NUM_GRADIENTS_EXP_3D + 2)
    i = trunc(hash) & ((OS2_NUM_GRADIENTS_3D - 1) << 2)
    t = (table[i+1], table[(i|1)+1], table[(i|2)+1])
    sum((t .* (x, y, z)))
end

@inline function os2_contribute1(seed, a, X, Y, Z, x1, y1, z1, x2, y2, z2, xs, ys, zs)
    result = 0.0
    if a > 0
        result += pow4(a) * grad(OS2_GRADIENTS_3D, seed, X, Y, Z, x1, y1, z1)
    end
    if x2 ≥ y2 && x2 ≥ z2
        result += os2_contribute2(seed, a + 2x2, X - xs * PRIME_X, Y, Z, x1 + xs, y1, z1)
    elseif y2 ≥ x2 && y2 ≥ z2
        result += os2_contribute2(seed, a + 2y2, X, Y - ys * PRIME_Y, Z, x1, y1 + ys, z1)
    else
        result += os2_contribute2(seed, a + 2z2, X, Y, Z - zs * PRIME_Z, x1, y1, z1 + zs)
    end
    result
end

@inline function os2_contribute2(seed, a, args...)
    a > 1 ? pow4(a - 1) * grad(OS2_GRADIENTS_3D, seed, args...) : 0.0
end

@inline function orient(::OpenSimplex2{3,OrientStandard}, x, y, z)
    OS2_FALLBACK_ROTATE_3D * (x + y + z) .- (x, y, z)
end

@inline function orient(::OpenSimplex2{3,OrientXY}, x, y, z)
    xy = x + y
    zz = z * ROOT_3_OVER_3
    xr, yr = (x, y) .+ xy .* OS2_ROTATE_3D_ORTHONORMALIZER .+ zz
    zr = xy * -ROOT_3_OVER_3 + zz
    (xr, yr, zr)
end

@inline orient(::OpenSimplex2{3,OrientXZ}, x, y, z) = orient(OrientXY, x, z, y)

function sample(sampler::OpenSimplex2{3,O}, x::T, y::T, z::T) where {O,T<:Real}
    seed = sampler.random_state.seed
    primes = (PRIME_X, PRIME_Y, PRIME_Z)
    tr = orient(sampler, x, y, z)
    V = round.(Int, tr)
    XYZ = V .* primes
    x1, y1, z1 = tr .- V
    s = trunc.(Int, -1 .- (x1, y1, z1)) .| 1
    XYZ2 = XYZ .+ s .>> 1 .& primes
    xyz2 = s .* .-((x1, y1, z1))
    x4, y4, z4 = 0.5 .- xyz2
    xyz3 = s .* (x4, y4, z4)
    a1 = OS2_R²3D - x1^2 - (y1^2 + z1^2)
    c1 = os2_contribute1(seed, a1, XYZ..., x1, y1, z1, xyz2..., s...)
    a2 = a1 + 0.75 - x4 - (y4 + z4)
    c2 = os2_contribute1(seed ⊻ OS2_SEED_FLIP_3D, a2, XYZ2..., xyz3..., x4, y4, z4, .-s...)
    c1 + c2
end

# 4D

const OS2_SEED_OFFSET_4D = 0xe83dc3e0da7164d
const OS2_SKEW_4D = -0.138196601125011
const OS2_UNSKEW_4D = 0.309016994374947
const OS2_LATTICE_STEP_4D = 0.2
const OS2_R²4D = 0.6
const OS2_NUM_GRADIENTS_EXP_4D = 9
const OS2_NUM_GRADIENTS_4D = 1 << OS2_NUM_GRADIENTS_EXP_4D
const OS2_GRADIENTS_NORMALIZED_4D = [
    -0.6740059517812944, -0.3239847771997537, -0.3239847771997537, 0.5794684678643381,
    -0.7504883828755602, -0.4004672082940195, 0.15296486218853164, 0.5029860367700724,
    -0.7504883828755602, 0.15296486218853164, -0.4004672082940195, 0.5029860367700724,
    -0.8828161875373585, 0.08164729285680945, 0.08164729285680945, 0.4553054119602712,
    -0.4553054119602712, -0.08164729285680945, -0.08164729285680945, 0.8828161875373585,
    -0.5029860367700724, -0.15296486218853164, 0.4004672082940195, 0.7504883828755602,
    -0.5029860367700724, 0.4004672082940195, -0.15296486218853164, 0.7504883828755602,
    -0.5794684678643381, 0.3239847771997537, 0.3239847771997537, 0.6740059517812944,
    -0.6740059517812944, -0.3239847771997537, 0.5794684678643381, -0.3239847771997537,
    -0.7504883828755602, -0.4004672082940195, 0.5029860367700724, 0.15296486218853164,
    -0.7504883828755602, 0.15296486218853164, 0.5029860367700724, -0.4004672082940195,
    -0.8828161875373585, 0.08164729285680945, 0.4553054119602712, 0.08164729285680945,
    -0.4553054119602712, -0.08164729285680945, 0.8828161875373585, -0.08164729285680945,
    -0.5029860367700724, -0.15296486218853164, 0.7504883828755602, 0.4004672082940195,
    -0.5029860367700724, 0.4004672082940195, 0.7504883828755602, -0.15296486218853164,
    -0.5794684678643381, 0.3239847771997537, 0.6740059517812944, 0.3239847771997537,
    -0.6740059517812944, 0.5794684678643381, -0.3239847771997537, -0.3239847771997537,
    -0.7504883828755602, 0.5029860367700724, -0.4004672082940195, 0.15296486218853164,
    -0.7504883828755602, 0.5029860367700724, 0.15296486218853164, -0.4004672082940195,
    -0.8828161875373585, 0.4553054119602712, 0.08164729285680945, 0.08164729285680945,
    -0.4553054119602712, 0.8828161875373585, -0.08164729285680945, -0.08164729285680945,
    -0.5029860367700724, 0.7504883828755602, -0.15296486218853164, 0.4004672082940195,
    -0.5029860367700724, 0.7504883828755602, 0.4004672082940195, -0.15296486218853164,
    -0.5794684678643381, 0.6740059517812944, 0.3239847771997537, 0.3239847771997537,
    0.5794684678643381, -0.6740059517812944, -0.3239847771997537, -0.3239847771997537,
    0.5029860367700724, -0.7504883828755602, -0.4004672082940195, 0.15296486218853164,
    0.5029860367700724, -0.7504883828755602, 0.15296486218853164, -0.4004672082940195,
    0.4553054119602712, -0.8828161875373585, 0.08164729285680945, 0.08164729285680945,
    0.8828161875373585, -0.4553054119602712, -0.08164729285680945, -0.08164729285680945,
    0.7504883828755602, -0.5029860367700724, -0.15296486218853164, 0.4004672082940195,
    0.7504883828755602, -0.5029860367700724, 0.4004672082940195, -0.15296486218853164,
    0.6740059517812944, -0.5794684678643381, 0.3239847771997537, 0.3239847771997537,
    -0.753341017856078, -0.37968289875261624, -0.37968289875261624, -0.37968289875261624,
    -0.7821684431180708, -0.4321472685365301, -0.4321472685365301, 0.12128480194602098,
    -0.7821684431180708, -0.4321472685365301, 0.12128480194602098, -0.4321472685365301,
    -0.7821684431180708, 0.12128480194602098, -0.4321472685365301, -0.4321472685365301,
    -0.8586508742123365, -0.508629699630796, 0.044802370851755174, 0.044802370851755174,
    -0.8586508742123365, 0.044802370851755174, -0.508629699630796, 0.044802370851755174,
    -0.8586508742123365, 0.044802370851755174, 0.044802370851755174, -0.508629699630796,
    -0.9982828964265062, -0.03381941603233842, -0.03381941603233842, -0.03381941603233842,
    -0.37968289875261624, -0.753341017856078, -0.37968289875261624, -0.37968289875261624,
    -0.4321472685365301, -0.7821684431180708, -0.4321472685365301, 0.12128480194602098,
    -0.4321472685365301, -0.7821684431180708, 0.12128480194602098, -0.4321472685365301,
    0.12128480194602098, -0.7821684431180708, -0.4321472685365301, -0.4321472685365301,
    -0.508629699630796, -0.8586508742123365, 0.044802370851755174, 0.044802370851755174,
    0.044802370851755174, -0.8586508742123365, -0.508629699630796, 0.044802370851755174,
    0.044802370851755174, -0.8586508742123365, 0.044802370851755174, -0.508629699630796,
    -0.03381941603233842, -0.9982828964265062, -0.03381941603233842, -0.03381941603233842,
    -0.37968289875261624, -0.37968289875261624, -0.753341017856078, -0.37968289875261624,
    -0.4321472685365301, -0.4321472685365301, -0.7821684431180708, 0.12128480194602098,
    -0.4321472685365301, 0.12128480194602098, -0.7821684431180708, -0.4321472685365301,
    0.12128480194602098, -0.4321472685365301, -0.7821684431180708, -0.4321472685365301,
    -0.508629699630796, 0.044802370851755174, -0.8586508742123365, 0.044802370851755174,
    0.044802370851755174, -0.508629699630796, -0.8586508742123365, 0.044802370851755174,
    0.044802370851755174, 0.044802370851755174, -0.8586508742123365, -0.508629699630796,
    -0.03381941603233842, -0.03381941603233842, -0.9982828964265062, -0.03381941603233842,
    -0.37968289875261624, -0.37968289875261624, -0.37968289875261624, -0.753341017856078,
    -0.4321472685365301, -0.4321472685365301, 0.12128480194602098, -0.7821684431180708,
    -0.4321472685365301, 0.12128480194602098, -0.4321472685365301, -0.7821684431180708,
    0.12128480194602098, -0.4321472685365301, -0.4321472685365301, -0.7821684431180708,
    -0.508629699630796, 0.044802370851755174, 0.044802370851755174, -0.8586508742123365,
    0.044802370851755174, -0.508629699630796, 0.044802370851755174, -0.8586508742123365,
    0.044802370851755174, 0.044802370851755174, -0.508629699630796, -0.8586508742123365,
    -0.03381941603233842, -0.03381941603233842, -0.03381941603233842, -0.9982828964265062,
    -0.3239847771997537, -0.6740059517812944, -0.3239847771997537, 0.5794684678643381,
    -0.4004672082940195, -0.7504883828755602, 0.15296486218853164, 0.5029860367700724,
    0.15296486218853164, -0.7504883828755602, -0.4004672082940195, 0.5029860367700724,
    0.08164729285680945, -0.8828161875373585, 0.08164729285680945, 0.4553054119602712,
    -0.08164729285680945, -0.4553054119602712, -0.08164729285680945, 0.8828161875373585,
    -0.15296486218853164, -0.5029860367700724, 0.4004672082940195, 0.7504883828755602,
    0.4004672082940195, -0.5029860367700724, -0.15296486218853164, 0.7504883828755602,
    0.3239847771997537, -0.5794684678643381, 0.3239847771997537, 0.6740059517812944,
    -0.3239847771997537, -0.3239847771997537, -0.6740059517812944, 0.5794684678643381,
    -0.4004672082940195, 0.15296486218853164, -0.7504883828755602, 0.5029860367700724,
    0.15296486218853164, -0.4004672082940195, -0.7504883828755602, 0.5029860367700724,
    0.08164729285680945, 0.08164729285680945, -0.8828161875373585, 0.4553054119602712,
    -0.08164729285680945, -0.08164729285680945, -0.4553054119602712, 0.8828161875373585,
    -0.15296486218853164, 0.4004672082940195, -0.5029860367700724, 0.7504883828755602,
    0.4004672082940195, -0.15296486218853164, -0.5029860367700724, 0.7504883828755602,
    0.3239847771997537, 0.3239847771997537, -0.5794684678643381, 0.6740059517812944,
    -0.3239847771997537, -0.6740059517812944, 0.5794684678643381, -0.3239847771997537,
    -0.4004672082940195, -0.7504883828755602, 0.5029860367700724, 0.15296486218853164,
    0.15296486218853164, -0.7504883828755602, 0.5029860367700724, -0.4004672082940195,
    0.08164729285680945, -0.8828161875373585, 0.4553054119602712, 0.08164729285680945,
    -0.08164729285680945, -0.4553054119602712, 0.8828161875373585, -0.08164729285680945,
    -0.15296486218853164, -0.5029860367700724, 0.7504883828755602, 0.4004672082940195,
    0.4004672082940195, -0.5029860367700724, 0.7504883828755602, -0.15296486218853164,
    0.3239847771997537, -0.5794684678643381, 0.6740059517812944, 0.3239847771997537,
    -0.3239847771997537, -0.3239847771997537, 0.5794684678643381, -0.6740059517812944,
    -0.4004672082940195, 0.15296486218853164, 0.5029860367700724, -0.7504883828755602,
    0.15296486218853164, -0.4004672082940195, 0.5029860367700724, -0.7504883828755602,
    0.08164729285680945, 0.08164729285680945, 0.4553054119602712, -0.8828161875373585,
    -0.08164729285680945, -0.08164729285680945, 0.8828161875373585, -0.4553054119602712,
    -0.15296486218853164, 0.4004672082940195, 0.7504883828755602, -0.5029860367700724,
    0.4004672082940195, -0.15296486218853164, 0.7504883828755602, -0.5029860367700724,
    0.3239847771997537, 0.3239847771997537, 0.6740059517812944, -0.5794684678643381,
    -0.3239847771997537, 0.5794684678643381, -0.6740059517812944, -0.3239847771997537,
    -0.4004672082940195, 0.5029860367700724, -0.7504883828755602, 0.15296486218853164,
    0.15296486218853164, 0.5029860367700724, -0.7504883828755602, -0.4004672082940195,
    0.08164729285680945, 0.4553054119602712, -0.8828161875373585, 0.08164729285680945,
    -0.08164729285680945, 0.8828161875373585, -0.4553054119602712, -0.08164729285680945,
    -0.15296486218853164, 0.7504883828755602, -0.5029860367700724, 0.4004672082940195,
    0.4004672082940195, 0.7504883828755602, -0.5029860367700724, -0.15296486218853164,
    0.3239847771997537, 0.6740059517812944, -0.5794684678643381, 0.3239847771997537,
    -0.3239847771997537, 0.5794684678643381, -0.3239847771997537, -0.6740059517812944,
    -0.4004672082940195, 0.5029860367700724, 0.15296486218853164, -0.7504883828755602,
    0.15296486218853164, 0.5029860367700724, -0.4004672082940195, -0.7504883828755602,
    0.08164729285680945, 0.4553054119602712, 0.08164729285680945, -0.8828161875373585,
    -0.08164729285680945, 0.8828161875373585, -0.08164729285680945, -0.4553054119602712,
    -0.15296486218853164, 0.7504883828755602, 0.4004672082940195, -0.5029860367700724,
    0.4004672082940195, 0.7504883828755602, -0.15296486218853164, -0.5029860367700724,
    0.3239847771997537, 0.6740059517812944, 0.3239847771997537, -0.5794684678643381,
    0.5794684678643381, -0.3239847771997537, -0.6740059517812944, -0.3239847771997537,
    0.5029860367700724, -0.4004672082940195, -0.7504883828755602, 0.15296486218853164,
    0.5029860367700724, 0.15296486218853164, -0.7504883828755602, -0.4004672082940195,
    0.4553054119602712, 0.08164729285680945, -0.8828161875373585, 0.08164729285680945,
    0.8828161875373585, -0.08164729285680945, -0.4553054119602712, -0.08164729285680945,
    0.7504883828755602, -0.15296486218853164, -0.5029860367700724, 0.4004672082940195,
    0.7504883828755602, 0.4004672082940195, -0.5029860367700724, -0.15296486218853164,
    0.6740059517812944, 0.3239847771997537, -0.5794684678643381, 0.3239847771997537,
    0.5794684678643381, -0.3239847771997537, -0.3239847771997537, -0.6740059517812944,
    0.5029860367700724, -0.4004672082940195, 0.15296486218853164, -0.7504883828755602,
    0.5029860367700724, 0.15296486218853164, -0.4004672082940195, -0.7504883828755602,
    0.4553054119602712, 0.08164729285680945, 0.08164729285680945, -0.8828161875373585,
    0.8828161875373585, -0.08164729285680945, -0.08164729285680945, -0.4553054119602712,
    0.7504883828755602, -0.15296486218853164, 0.4004672082940195, -0.5029860367700724,
    0.7504883828755602, 0.4004672082940195, -0.15296486218853164, -0.5029860367700724,
    0.6740059517812944, 0.3239847771997537, 0.3239847771997537, -0.5794684678643381,
    0.03381941603233842, 0.03381941603233842, 0.03381941603233842, 0.9982828964265062,
    -0.044802370851755174, -0.044802370851755174, 0.508629699630796, 0.8586508742123365,
    -0.044802370851755174, 0.508629699630796, -0.044802370851755174, 0.8586508742123365,
    -0.12128480194602098, 0.4321472685365301, 0.4321472685365301, 0.7821684431180708,
    0.508629699630796, -0.044802370851755174, -0.044802370851755174, 0.8586508742123365,
    0.4321472685365301, -0.12128480194602098, 0.4321472685365301, 0.7821684431180708,
    0.4321472685365301, 0.4321472685365301, -0.12128480194602098, 0.7821684431180708,
    0.37968289875261624, 0.37968289875261624, 0.37968289875261624, 0.753341017856078,
    0.03381941603233842, 0.03381941603233842, 0.9982828964265062, 0.03381941603233842,
    -0.044802370851755174, 0.044802370851755174, 0.8586508742123365, 0.508629699630796,
    -0.044802370851755174, 0.508629699630796, 0.8586508742123365, -0.044802370851755174,
    -0.12128480194602098, 0.4321472685365301, 0.7821684431180708, 0.4321472685365301,
    0.508629699630796, -0.044802370851755174, 0.8586508742123365, -0.044802370851755174,
    0.4321472685365301, -0.12128480194602098, 0.7821684431180708, 0.4321472685365301,
    0.4321472685365301, 0.4321472685365301, 0.7821684431180708, -0.12128480194602098,
    0.37968289875261624, 0.37968289875261624, 0.753341017856078, 0.37968289875261624,
    0.03381941603233842, 0.9982828964265062, 0.03381941603233842, 0.03381941603233842,
    -0.044802370851755174, 0.8586508742123365, -0.044802370851755174, 0.508629699630796,
    -0.044802370851755174, 0.8586508742123365, 0.508629699630796, -0.044802370851755174,
    -0.12128480194602098, 0.7821684431180708, 0.4321472685365301, 0.4321472685365301,
    0.508629699630796, 0.8586508742123365, -0.044802370851755174, -0.044802370851755174,
    0.4321472685365301, 0.7821684431180708, -0.12128480194602098, 0.4321472685365301,
    0.4321472685365301, 0.7821684431180708, 0.4321472685365301, -0.12128480194602098,
    0.37968289875261624, 0.753341017856078, 0.37968289875261624, 0.37968289875261624,
    0.9982828964265062, 0.03381941603233842, 0.03381941603233842, 0.03381941603233842,
    0.8586508742123365, -0.044802370851755174, -0.044802370851755174, 0.508629699630796,
    0.8586508742123365, -0.044802370851755174, 0.508629699630796, -0.044802370851755174,
    0.7821684431180708, -0.12128480194602098, 0.4321472685365301, 0.4321472685365301,
    0.8586508742123365, 0.508629699630796, -0.044802370851755174, -0.044802370851755174,
    0.7821684431180708, 0.4321472685365301, -0.12128480194602098, 0.4321472685365301,
    0.7821684431180708, 0.4321472685365301, 0.4321472685365301, -0.12128480194602098,
    0.753341017856078, 0.37968289875261624, 0.37968289875261624, 0.37968289875261624]
const OS2_GRADIENTS_4D = OS2_GRADIENTS_NORMALIZED_4D ./ 0.0220065933241897 |> CircularVector

@doc doc_opensimplex2_4d
opensimplex2_4d(; seed=0, orient=nothing) = _opensimplex2(4, seed, orient)

@inline function grad(table, seed, X, Y, Z, W, x, y, z, w)
    hash = seed ⊻ (X ⊻ Y) ⊻ (Z ⊻ W) * HASH_MULTIPLIER
    hash ⊻= hash >> (64 - OS2_NUM_GRADIENTS_EXP_4D + 2)
    i = trunc(hash) & ((OS2_NUM_GRADIENTS_4D - 1) << 2)
    t = (table[i+1], table[(i|1)+1], table[(i|2)+1], table[(i|3)+1])
    sum((t .* (x, y, z, w)))
end

@inline function orient(::OpenSimplex2{4,OrientStandard}, x, y, z, w)
    (x, y, z, w) .+ OS2_SKEW_4D .* (x + y + z + w)
end

@inline function orient(::OpenSimplex2{4,OrientXY}, x, y, z, w)
    xy = x + y
    ww = w * 0.2236067977499788
    zw = z * 0.28867513459481294226 + ww
    xr, yr = (x, y) .+ zw .+ xy .* -0.21132486540518699998
    zr = xy * -0.57735026918962599998 + zw
    wr = z * -0.866025403784439 + ww
    (xr, yr, zr, wr)
end

@inline orient(::OpenSimplex2{4,OrientXZ}, x, y, z, w) = orient(OrientXY, x, z, y, w)

@inline function orient(::OpenSimplex2{4,OrientXYZ}, x, y, z, w)
    xyz = -(x + y + z)
    ww = w * 0.2236067977499788
    s = xyz / 6 + ww
    xs, ys, zs = (x, y, z) .+ s
    ws = xyz * 0.5 + ww
    (xs, ys, zs, ws)
end

function sample(sampler::OpenSimplex2{4,O}, x::T, y::T, z::T, w::T) where {O,T<:Real}
    seed = sampler.random_state.seed
    primes = (PRIME_X, PRIME_Y, PRIME_Z, PRIME_W)
    tr = orient(sampler, x, y, z, w)
    X1, Y1, Z1, W1 = floor.(Int, tr)
    X2, Y2, Z2, W2 = (X1, Y1, Z1, W1) .* primes
    v = tr .- (X1, Y1, Z1, W1)
    sv = sum(v)
    lattice = trunc(Int, sv * 1.25)
    lattice_offset = lattice * -OS2_LATTICE_STEP_4D
    x1, y1, z1, w1 = v .+ lattice_offset
    ssi = (sv + 4lattice_offset) * OS2_UNSKEW_4D
    seed += lattice * OS2_SEED_OFFSET_4D
    result = 0.0
    for i in 0:4
        score = 1 + ssi * (-1 / OS2_UNSKEW_4D)
        if x1 ≥ y1 && x1 ≥ z1 && x1 ≥ w1 && x1 ≥ score
            X2 += PRIME_X
            x1 -= 1
            ssi -= OS2_UNSKEW_4D
        elseif y1 > x1 && y1 ≥ z1 && y1 ≥ w1 && y1 ≥ score
            Y2 += PRIME_Y
            y1 -= 1
            ssi -= OS2_UNSKEW_4D
        elseif z1 > x1 && z1 > y1 && z1 ≥ w1 && z1 ≥ score
            Z2 += PRIME_Z
            z1 -= 1
            ssi -= OS2_UNSKEW_4D
        elseif w1 > x1 && w1 > y1 && w1 > z1 && w1 ≥ score
            W2 += PRIME_W
            w1 -= 1
            ssi -= OS2_UNSKEW_4D
        end
        xyzw = (x1, y1, z1, w1) .+ ssi
        a = sum(xyzw .^ 2)
        if a < OS2_R²4D
            result += pow4(a - OS2_R²4D) * grad(OS2_GRADIENTS_4D, seed, X2, Y2, Z2, W2, xyzw...)
        end
        if i !== 4
            x1, y1, z1, w1 = (x1, y1, z1, w1) .+ OS2_LATTICE_STEP_4D
            ssi += OS2_LATTICE_STEP_4D * 4OS2_UNSKEW_4D
            seed -= OS2_SEED_OFFSET_4D
            if i == lattice
                X2, Y2, Z2, W2 = (X2, Y2, Z2, W2) .- primes
                seed += 5OS2_SEED_OFFSET_4D
            end
        end
    end
    result
end
