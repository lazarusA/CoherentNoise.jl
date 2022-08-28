abstract type DistanceMetric end
struct Manhattan <: DistanceMetric end
struct Euclidean <: DistanceMetric end
struct Euclidean² <: DistanceMetric end
struct Chebyshev <: DistanceMetric end
struct Minkowski4 <: DistanceMetric end

abstract type WorleyOutput end
struct CellF1 <: WorleyOutput end
struct CellF2 <: WorleyOutput end
struct CellAdd <: WorleyOutput end
struct CellSub <: WorleyOutput end
struct CellMul <: WorleyOutput end
struct CellDiv <: WorleyOutput end
struct CellValue <: WorleyOutput end

struct Worley{N,M<:DistanceMetric,F<:WorleyOutput} <: NoiseSampler{N}
    random_state::RandomState
    table::Vector{Float64}
    jitter::Float64
end

const WORLEY_JITTER1 = 0.43701595
const WORLEY_JITTER2 = 0.39614353

@inline function Worley{N,M,F}(seed, jitter) where {N,M,F}
    rs = RandomState(seed)
    table = rand(rs.rng, Uniform(-1.0, 1.0), 2^(N - 1) * 256)
    Worley{N,M,F}(rs, table, jitter)
end

@inline function _worley(dims, seed, metric, output, jitter)
    metric = worley_metric_type(Val(metric))
    output = worley_output_type(Val(output))
    Worley{dims,metric,output}(seed, jitter)
end

@inline worley_metric_type(::Val{:manhattan}) = Manhattan
@inline worley_metric_type(::Val{:euclidean}) = Euclidean
@inline worley_metric_type(::Val{:euclidean²}) = Euclidean²
@inline worley_metric_type(::Val{:chebyshev}) = Chebyshev
@inline worley_metric_type(::Val{:minkowski4}) = Minkowski4
@inline worley_metric_type(x::Any) = @error "invalid Worley noise metric type: " x

@inline worley_output_type(::Val{:f1}) = CellF1
@inline worley_output_type(::Val{:f2}) = CellF2
@inline worley_output_type(::Val{:+}) = CellAdd
@inline worley_output_type(::Val{:-}) = CellSub
@inline worley_output_type(::Val{:*}) = CellMul
@inline worley_output_type(::Val{:/}) = CellDiv
@inline worley_output_type(::Val{:value}) = CellValue
@inline worley_output_type(x::Any) = @error "invalid Worley noise output type: " x

@inline cell_distance(::Type{Manhattan}, v...) = sum(abs, v)
@inline cell_distance(::Type{Euclidean}, v...) = cell_distance(Euclidean², v...)
@inline cell_distance(::Type{Euclidean²}, v...) = sum(v .^ 2)
@inline cell_distance(::Type{Chebyshev}, v...) = mapreduce(abs, max, v)
@inline cell_distance(::Type{Minkowski4}, v...) = sum(v .^ 2 .* v .^ 2)^0.25

@inline cell_value(::Type{CellF1}, _, min, _) = min
@inline cell_value(::Type{CellF2}, _, _, max) = max
@inline cell_value(::Type{CellAdd}, _, min, max) = (min + max) * 0.5
@inline cell_value(::Type{CellSub}, _, min, max) = max - min
@inline cell_value(::Type{CellMul}, _, min, max) = min * max
@inline cell_value(::Type{CellDiv}, _, min, max) = min / max
@inline cell_value(::Type{CellValue}, hash, _, _) = hash % UInt32 / HASH2
@inline cell_value(F, ::Type{Euclidean}, hash, min, max) = cell_value(F, hash, sqrt(min), sqrt(max))
@inline cell_value(F, ::Type{<:DistanceMetric}, args...) = cell_value(F, args...)

# 1D

@doc doc_worley_1d
function worley_1d(; seed=0, metric=:euclidean, output=:f1, jitter=1.0)
    _worley(1, seed, metric, output, jitter)
end

function sample(sampler::Worley{1,M,F}, x::Real) where {M,F}
    seed = sampler.random_state.seed
    table = sampler.table
    jitter = sampler.jitter * WORLEY_JITTER1
    r = round(Int, x) .- 1
    xr = r .- x
    xp = r * PRIME_X
    minf = floatmax(Float64)
    maxf = minf
    closest_hash::UInt32 = 0
    @inbounds for xi in 0:2
        hash = ⊻(seed, xp) * HASH1 % UInt32
        vx = table[(hash+1)&255] * jitter + xr + xi
        d = cell_distance(M, vx)
        maxf = clamp(d, minf, maxf)
        if d < minf
            minf = d
            closest_hash = hash
        end
        xp += PRIME_X
    end
    cell_value(F, M, closest_hash, minf, maxf) - 1
end

# 2D

@doc doc_worley_2d
function worley_2d(; seed=0, metric=:euclidean, output=:f1, jitter=1.0)
    _worley(2, seed, metric, output, jitter)
end

function sample(sampler::Worley{2,M,F}, x::T, y::T) where {M,F,T<:Real}
    seed = sampler.random_state.seed
    table = sampler.table
    jitter = sampler.jitter * WORLEY_JITTER1
    r = round.(Int, (x, y)) .- 1
    xr, yr = r .- (x, y)
    xp, yp_base = r .* (PRIME_X, PRIME_Y)
    minf = floatmax(Float64)
    maxf = minf
    closest_hash::UInt32 = 0
    @inbounds for xi in 0:2
        xri = xr + xi
        yp = yp_base
        for yi in 0:2
            hash = ⊻(seed, xp, yp) * HASH1 % UInt32
            vx = table[(hash+1)&511] * jitter + xri
            vy = table[((hash|1)+1)&511] * jitter + yr + yi
            d = cell_distance(M, vx, vy)
            maxf = clamp(d, minf, maxf)
            if d < minf
                minf = d
                closest_hash = hash
            end
            yp += PRIME_Y
        end
        xp += PRIME_X
    end
    cell_value(F, M, closest_hash, minf, maxf) - 1
end

# 3D

@doc doc_worley_3d
function worley_3d(; seed=0, metric=:euclidean, output=:f1, jitter=1.0)
    _worley(3, seed, metric, output, jitter)
end

function sample(sampler::Worley{3,M,F}, x::T, y::T, z::T) where {M,F,T<:Real}
    seed = sampler.random_state.seed
    table = sampler.table
    jitter = sampler.jitter * WORLEY_JITTER2
    r = round.(Int, (x, y, z)) .- 1
    xr, yr, zr = r .- (x, y, z)
    xp, yp_base, zp_base = r .* (PRIME_X, PRIME_Y, PRIME_Z)
    minf = floatmax(Float64)
    maxf = minf
    closest_hash::UInt32 = 0
    @inbounds for xi in 0:2
        xri = xr + xi
        yp = yp_base
        for yi in 0:2
            yri = yr + yi
            zp = zp_base
            for zi in 0:2
                hash = ⊻(seed, xp, yp, zp) * HASH1 % UInt32
                vx = table[(hash+1)&1023] * jitter + xri
                vy = table[((hash|1)+1)&1023] * jitter + yri
                vz = table[((hash|2)+1)&1023] * jitter + zr + zi
                d = cell_distance(M, vx, vy, vz)
                maxf = clamp(d, minf, maxf)
                if d < minf
                    minf = d
                    closest_hash = hash
                end
                zp += PRIME_Z
            end
            yp += PRIME_Y
        end
        xp += PRIME_X
    end
    cell_value(F, M, closest_hash, minf, maxf) - 1
end

# 4D

@doc doc_worley_4d
function worley_4d(; seed=0, metric=:euclidean, output=:f1, jitter=1.0)
    _worley(4, seed, metric, output, jitter)
end

function sample(sampler::Worley{4,M,F}, x::T, y::T, z::T, w::T) where {M,F,T<:Real}
    seed = sampler.random_state.seed
    table = sampler.table
    jitter = sampler.jitter * WORLEY_JITTER2
    r = round.(Int, (x, y, z, w)) .- 1
    xr, yr, zr, wr = r .- (x, y, z, w)
    xp, yp_base, zp_base, wp_base = r .* (PRIME_X, PRIME_Y, PRIME_Z, PRIME_W)
    minf = floatmax(Float64)
    maxf = minf
    closest_hash::UInt32 = 0
    @inbounds for xi in 0:2
        xri = xr + xi
        yp = yp_base
        for yi in 0:2
            yri = yr + yi
            zp = zp_base
            for zi in 0:2
                zri = zr + zi
                wp = wp_base
                for wi in 0:2
                    hash = ⊻(seed, xp, yp, zp, wp) * HASH1 % UInt32
                    vx = table[(hash+1)&2047] * jitter + xri
                    vy = table[((hash|1)+1)&2047] * jitter + yri
                    vz = table[((hash|2)+1)&2047] * jitter + zri
                    vw = table[((hash|3)+1)&2047] * jitter + wr + wi
                    d = cell_distance(M, vx, vy, vz, vw)
                    maxf = clamp(d, minf, maxf)
                    if d < minf
                        minf = d
                        closest_hash = hash
                    end
                    wp += PRIME_W
                end
                zp += PRIME_Z
            end
            yp += PRIME_Y
        end
        xp += PRIME_X
    end
    cell_value(F, M, closest_hash, minf, maxf) - 1
end
