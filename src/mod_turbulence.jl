struct Turbulence{N,N1,N2,S1,S2} <: ModifierSampler{N}
    random_state::RandomState
    source1::S1
    source2::S2
    power::Float64
    x::NTuple{N2,Float64}
    y::NTuple{N2,Float64}
    z::NTuple{N2,Float64}
    w::NTuple{N2,Float64}
end

"""
    turbulence(s1::AbstractSampler; s2::AbstractSampler; kwargs...)

Construct a modifier sampler that displaces the input coordinates of sampler `s1` by the output of
sampler `s2` with a fractional Brownian motion fractal applied to it.

Sampler `s2`'s input coordinates are randomly generated using the seed of sampler `s1`.

# Arguments

  - `frequency::Real=1.0`: The frequency of the fractal signal to apply to sampler `s2`.
  - `roughness::Real=3`: The number of octaves of the fractal to apply to sampler `s2`.
  - `power::Real=1.0`: A scaling factor that is applied to the displaced result before sampling from
    sampler `s1`.
"""
function turbulence(
    s1::S1,
    s2::S2;
    frequency::Real=1.0, roughness::Real=3, power::Real=1.0,
) where {
    N,N1,
    S1<:AbstractSampler{N},
    S2<:AbstractSampler{N1},
}
    rs = random_state(s1)
    N2 = min(N, N1)
    x = ntuple(i -> rand(rs.rng, Float64), N2)
    y = ntuple(i -> rand(rs.rng, Float64), N2)
    z = ntuple(i -> rand(rs.rng, Float64), N2)
    w = ntuple(i -> rand(rs.rng, Float64), N2)
    s3 = FBM{N1}(seed=rand(rs.rng, UInt64), source=s2, octaves=roughness, frequency=frequency)
    S3 = typeof(s3)
    Turbulence{N,N1,N2,S1,S3}(rs, s1, s3, Float64(power), x, y, z, w)
end

function sample(sampler::Turbulence{N,N1,N2}, coords::Vararg{Real,N}) where {N,N1,N2}
    s2 = sampler.source2
    power = sampler.power
    zeros = ntuple(i -> 0.0, max(0, N1 - N))
    tx = sample(s2, (sampler.x .+ coords[1:N2])..., zeros...) * power
    ty = sample(s2, (sampler.y .+ coords[1:N2])..., zeros...) * power
    tz = sample(s2, (sampler.z .+ coords[1:N2])..., zeros...) * power
    tw = sample(s2, (sampler.w .+ coords[1:N2])..., zeros...) * power
    sample(sampler.source1, (coords .+ (tx, ty, tz, tw)[1:N])...)
end
