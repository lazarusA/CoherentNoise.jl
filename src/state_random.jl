using Random: Xoshiro, RandomDevice, shuffle

struct RandomState
    seed::UInt64
    rng::Xoshiro
end

const Seed = Union{Integer,Nothing}

RandomState(seed::Integer) = RandomState(seed, Xoshiro(seed))
RandomState(::Nothing) = RandomState(rand(RandomDevice(), UInt64))

@inline random_state(sampler::AbstractSampler) = sampler.random_state
@inline rng(sampler::AbstractSampler) = sampler.random_state.rng
@inline get_seed(sampler::AbstractSampler) = sampler.random_state.seed
