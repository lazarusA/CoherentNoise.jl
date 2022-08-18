struct Select{N,S1,S2,C} <: ModifierSampler{N}
    random_state::RandomState
    source1::S1
    source2::S2
    control::C
    min::Float64
    max::Float64
    falloff::Float64
end

"""
    select(x::AbstractSampler, y::AbstractSampler; kwargs...)

Construct a modifier sampler that outputs either the out of sampler `x` or `y`, depending on the
output of sampler `z`.

If the output of sampler `z` is within the range denoted by `min` and `max`, the output of sampler
`y` is chosen. If the output of sampler `z` is outside of this range, the output of sampler `x` is
chosen.

# Arguments

  - `min::Real=-1.0`: A real number between -1.0 and 1.0 defining the lower bound of the selection
    range.
  - `max::Real=1.0`: A real number between -1.0 and 1.0 defining the upper bound of the selection
    range.
  - `falloff::Real=1.0`: A real number between 0.0 and 1.0 specifying the smoothness of the
    transition.
"""
function select(
    x::S1,
    y::S2,
    z::C;
    min::Real=-1.0,
    max::Real=1.0,
    falloff::Real=0.0,
) where {N1,N2,N3,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2},C<:AbstractSampler{N3}}
    N = Base.max(N1, N2, N3)
    Select{N,S1,S2,C}(random_state(x), x, y, z, Float64(min), Float64(max), Float64(falloff))
end

@inline function sample(
    sampler::Select{N,S1,S2,C},
    coords::Vararg{Real,N},
) where {N,N1,N2,N3,S1<:AbstractSampler{N1},S2<:AbstractSampler{N2},C<:AbstractSampler{N3}}
    min = sampler.min
    max = sampler.max
    falloff = Base.min(sampler.falloff, (max - min) * 0.5)
    s1 = sample(sampler.source1, coords[1:N1]...)
    s2 = sample(sampler.source2, coords[1:N2]...)
    c = sample(sampler.control, coords[1:N3]...)
    if falloff > 0
        minmf, maxmf = min - falloff, max - falloff
        minpf, maxpf = min + falloff, max + falloff
        if c < minmf
            s1
        elseif c < minpf
            lerp(s1, s2, curve3((c - minmf) / (minpf - minmf)))
        elseif c < maxpf
            lerp(s1, s2, curve3((c - maxmf) / (maxpf - maxmf)))
        else
            s1
        end
    else
        c < min || c > max ? s1 : s2
    end
end
