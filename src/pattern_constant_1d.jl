struct Constant <: PatternSampler{1}
    random_state::RandomState
    value::Float64
end

"""
    constant_1d(; kwargs...)

Construct a sampler that constantly outputs `value` each time it is sampled from.

This is useful for debugging and applications where you need to combine a constant value.

# Arguments

  - `seed=0`: An integer used to seed the random number generator for this sampler.

  - `value=0.0`: A constant value to emit each time this sampler is sampled from.
"""
constant_1d(; seed=0, value=0.0) = Constant(RandomState(seed), value)

sample(sampler::S, _) where {S<:Constant} = sampler.value
