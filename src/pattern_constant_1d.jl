struct Constant <: PatternSampler{1}
    random_state::RandomState
    value::Float64
end

"""
    constant_1d(; seed=nothing, value=0.0)

Construct a sampler that constantly outputs `value` each time it is sampled from.

This is useful for debugging and applications where you need to combine a constant value.

# Arguments

  - `seed=nothing`: An integer used to seed the random number generator for this sampler, or
    `nothing`. If a seed is not supplied, one will be generated automatically which will negatively
    affect the reproducibility of any samplers further down the pipeline.
  - `value=0.0`: A constant value to emit each time this sampler is sampled from.
"""
constant_1d(; seed=nothing, value=0.0) = Constant(RandomState(seed), value)

sample(sampler::S, ::Real) where {S<:Constant} = sampler.value
