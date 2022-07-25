module CellularNoise

using Distributions: Uniform
using ...Common: RandomState, Seed, get_seed
using ...Common: HASH1, HASH2, PRIME_X, PRIME_Y, PRIME_Z, PRIME_W
import ...Common: sample
using ..Noise: NoiseSampler

"""
Supertype of all cellular noise distance metrics.
"""
abstract type DistanceMetric end

"""
Use the Manhattan distance metric for determining the distance between two cells.
"""
struct Manhattan <: DistanceMetric end
"""
Use the Euclidean distance metric for determining the distance between two cells.
"""
struct Euclidean <: DistanceMetric end
"""
Use the squared Euclidean distance metric for determining the distance between two cells.
"""
struct Euclidean² <: DistanceMetric end
"""
Use the Chebyshev distance metric for determining the distance between two cells.
"""
struct Chebyshev <: DistanceMetric end
"""
Use the Minkowski4 (p=4) distance metric for determining the distance between two cells.
"""
struct Minkowski4 <: DistanceMetric end

"""
Supertype of all cellular noise cell functions.
"""
abstract type CellFunction end

struct Distance1 <: CellFunction end
struct Distance2 <: CellFunction end
struct DistanceAdd <: CellFunction end
struct DistanceSub <: CellFunction end
struct DistanceMul <: CellFunction end
struct DistanceDiv <: CellFunction end
struct CellValue <: CellFunction end

"""
An `N`-dimensional sampler that generates cellular (Voronoi/Worley) noise.

`N` must be an integer in the range [2, 4].
`M` must be a [`DistanceMetric`](@ref) type.
`F` must be a [`CellFunction`](@ref) type.
"""
struct Cellular{N,M<:DistanceMetric,F<:CellFunction} <: NoiseSampler{N}
    random_state::RandomState
    table::Vector{Float64}
    jitter::Float64
end

"""
    Cellular{N}(; seed=nothing)
    Cellular{N,M}(; seed=nothing)
    Cellular{N,M,F}(; seed=nothing)

Construct a sampler that outputs `N`-dimensional cellular noise when it is sampler from.

`M` if supplied, must be a [`DistanceMetric`](@ref) type, which denotes the method to use for
determining adjacent cells.

`F` if supplied, must be a [`CellFunction`](@ref) type, which denotes the "cell function" used in
the cell distance calculation.

# Arguments

  - `seed::Union{Int,Nothing}=nothing`: An integer used to seed the random number generator for this
    sampler, or `nothing`. If a seed is not supplied, one will be generated automatically which will
    negatively affect reproducibility.
  - `jitter::Real=1.0`: A `Real` number between 0.0 and 1.0, with values closer to one randomly
    distributing cells away from their grid alignment.
"""
function Cellular{N}(; seed::Seed=nothing, jitter::Real=1.0) where {N}
    Cellular{N,Euclidean,Distance1}(seed=seed, jitter=jitter)
end

function Cellular{N,M}(; seed::Seed=nothing, jitter::Real=1.0) where {N,M}
    Cellular{N,M,Distance1}(seed=seed, jitter=jitter)
end

function Cellular{N,M,F}(; seed::Seed=nothing, jitter::Real=1.0) where {N,M,F}
    rs = RandomState(seed)
    table = rand(rs.rng, Uniform(-1.0, 1.0), table_size(Val(N)))
    Cellular{N,M,F}(rs, table, jitter)
end

const JITTER1 = 0.43701595

const JITTER2 = 0.39614353

@inline cell_distance(::Type{Manhattan}, v...) = sum(abs, v)
@inline cell_distance(::Type{Euclidean}, v...) = cell_distance(Euclidean², v...)
@inline cell_distance(::Type{Euclidean²}, v...) = sum(v .^ 2)
@inline cell_distance(::Type{Chebyshev}, v...) = mapreduce(abs, max, v)
@inline cell_distance(::Type{Minkowski4}, v...) = sum(v .^ 2 .* v .^ 2)^0.25

@inline cell_value(::Type{Distance1}, _, min, _) = min
@inline cell_value(::Type{Distance2}, _, _, max) = max
@inline cell_value(::Type{DistanceAdd}, _, min, max) = (min + max) * 0.5
@inline cell_value(::Type{DistanceSub}, _, min, max) = max - min
@inline cell_value(::Type{DistanceMul}, _, min, max) = min * max
@inline cell_value(::Type{DistanceDiv}, _, min, max) = min / max
@inline cell_value(::Type{CellValue}, hash, _, _) = hash % UInt32 / HASH2
@inline cell_value(F, ::Type{Euclidean}, hash, min, max) = cell_value(F, hash, sqrt(min), sqrt(max))
@inline cell_value(F, ::Type{<:DistanceMetric}, args...) = cell_value(F, args...)

@inline table_size(::Val{2}) = 512
@inline table_size(::Val{3}) = 1024
@inline table_size(::Val{4}) = 2048

include("noise_cellular_2d.jl")
include("noise_cellular_3d.jl")
include("noise_cellular_4d.jl")

end
