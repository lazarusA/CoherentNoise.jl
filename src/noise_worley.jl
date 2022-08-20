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
    table = rand(rs.rng, Uniform(-1.0, 1.0), worley_table_size(Val(N)))
    Worley{N,M,F}(rs, table, jitter)
end

@inline function worley(dims, seed, metric, output, jitter)
    metric = worley_metric_type(Val(metric))
    output = worley_output_type(Val(output))
    Worley{dims,metric,output}(seed, jitter)
end

@inline worley_metric_type(::Val{:manhattan}) = Manhattan
@inline worley_metric_type(::Val{:euclidean}) = Euclidean
@inline worley_metric_type(::Val{:euclidean²}) = Euclidean²
@inline worley_metric_type(::Val{:chebyshev}) = Chebyshev
@inline worley_metric_type(::Val{:minkowski4}) = Minkowski4
@inline worley_metric_type(x::Any) = @error "invalid Worley noise metric: " x

@inline worley_output_type(::Val{:f1}) = CellF1
@inline worley_output_type(::Val{:f2}) = CellF2
@inline worley_output_type(::Val{:+}) = CellAdd
@inline worley_output_type(::Val{:-}) = CellSub
@inline worley_output_type(::Val{:*}) = CellMul
@inline worley_output_type(::Val{:/}) = CellDiv
@inline worley_output_type(::Val{:value}) = CellValue
@inline worley_output_type(x::Any) = @error "invalid Worley noise output: " x

@inline worley_table_size(::Val{2}) = 512
@inline worley_table_size(::Val{3}) = 1024
@inline worley_table_size(::Val{4}) = 2048

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

include("noise_worley_2d.jl")
include("noise_worley_3d.jl")
include("noise_worley_4d.jl")
