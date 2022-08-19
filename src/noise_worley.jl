module WorleyNoise

using Distributions: Uniform
using ...Common: RandomState, Seed, get_seed
using ...Common: HASH1, HASH2, PRIME_X, PRIME_Y, PRIME_Z, PRIME_W
import ...Common: sample
using ..Noise: NoiseSampler

abstract type Metric end
struct Manhattan <: Metric end
struct Euclidean <: Metric end
struct Euclidean² <: Metric end
struct Chebyshev <: Metric end
struct Minkowski4 <: Metric end

abstract type Output end
struct F1 <: Output end
struct F2 <: Output end
struct Add <: Output end
struct Sub <: Output end
struct Mul <: Output end
struct Div <: Output end
struct Value <: Output end

struct Worley{N,M<:Metric,F<:Output} <: NoiseSampler{N}
    random_state::RandomState
    table::Vector{Float64}
    jitter::Float64
end

const JITTER1 = 0.43701595
const JITTER2 = 0.39614353

@inline function Worley{N,M,F}(seed::Seed, jitter::Real) where {N,M,F}
    rs = RandomState(seed)
    table = rand(rs.rng, Uniform(-1.0, 1.0), table_size(Val(N)))
    Worley{N,M,F}(rs, table, jitter)
end

@inline function worley(dims, seed, metric, output, jitter)
    metric = metric_type(Val(metric))
    output = output_type(Val(output))
    Worley{dims,metric,output}(seed, jitter)
end

@inline metric_type(::Val{:manhattan}) = Manhattan
@inline metric_type(::Val{:euclidean}) = Euclidean
@inline metric_type(::Val{:euclidean²}) = Euclidean²
@inline metric_type(::Val{:chebyshev}) = Chebyshev
@inline metric_type(::Val{:minkowski4}) = Minkowski4
@inline metric_type(x::Any) = @error "invalid Worley noise metric: " x

@inline output_type(::Val{:f1}) = F1
@inline output_type(::Val{:f2}) = F2
@inline output_type(::Val{:+}) = Add
@inline output_type(::Val{:-}) = Sub
@inline output_type(::Val{:*}) = Mul
@inline output_type(::Val{:/}) = Div
@inline output_type(::Val{:value}) = Value
@inline output_type(x::Any) = @error "invalid Worley noise output: " x

@inline cell_distance(::Type{Manhattan}, v...) = sum(abs, v)
@inline cell_distance(::Type{Euclidean}, v...) = cell_distance(Euclidean², v...)
@inline cell_distance(::Type{Euclidean²}, v...) = sum(v .^ 2)
@inline cell_distance(::Type{Chebyshev}, v...) = mapreduce(abs, max, v)
@inline cell_distance(::Type{Minkowski4}, v...) = sum(v .^ 2 .* v .^ 2)^0.25

@inline cell_value(::Type{F1}, _, min, _) = min
@inline cell_value(::Type{F2}, _, _, max) = max
@inline cell_value(::Type{Add}, _, min, max) = (min + max) * 0.5
@inline cell_value(::Type{Sub}, _, min, max) = max - min
@inline cell_value(::Type{Mul}, _, min, max) = min * max
@inline cell_value(::Type{Div}, _, min, max) = min / max
@inline cell_value(::Type{Value}, hash, _, _) = hash % UInt32 / HASH2
@inline cell_value(F, ::Type{Euclidean}, hash, min, max) = cell_value(F, hash, sqrt(min), sqrt(max))
@inline cell_value(F, ::Type{<:Metric}, args...) = cell_value(F, args...)

@inline table_size(::Val{2}) = 512
@inline table_size(::Val{3}) = 1024
@inline table_size(::Val{4}) = 2048

include("noise_worley_2d.jl")
include("noise_worley_3d.jl")
include("noise_worley_4d.jl")

end
