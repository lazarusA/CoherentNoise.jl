struct Checkered <: PatternSampler{2}
    random_state::RandomState
end

"""
    checkered_2d(; seed=nothing)

Construct a sampler that outputs values in a checkerboard-like pattern when it is sampled from.

That is, output values will only ever be -1.0 or 1.0.

# Arguments

  - `seed=nothing`: An integer used to seed the random number generator for this sampler, or
    `nothing`. If a seed is not supplied, one will be generated automatically which will negatively
    affect the reproducibility of any samplers further down the pipeline.
"""
checkered_2d(; seed=nothing) = Checkered(RandomState(seed))

function sample(::S, x::T, y::T) where {S<:Checkered,T<:Real}
    iszero((floor(Int, x) & 1) ⊻ (floor(Int, y) & 1)) ? 1.0 : -1.0
end
