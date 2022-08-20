struct Terrace{N,S} <: ModifierSampler{N}
    random_state::RandomState
    source::S
    points::Vector{Float64}
    invert::Bool
end

"""
    terrace(x::AbstractSampler, points::Vector{Pair{Float64,Float64}}; invert::Bool=false)

Construct a modifier sampler that outputs the result of sampling from `x` after remapping its output
to a terrace-forming curve.

The curve is defined by a `Vector` of `Float64`s given by `points`. Each point represents an input
and output number.

When sampling from sampler `x`, the output is evaluated using the curve data, and maps it to a new
output value.

# Arguments

  - `invert::Bool=false`: Specify whether the curve is inverted between control points.
"""
function terrace(
    x::S,
    points::Vector{Float64};
    invert::Bool=false,
) where {N,S<:AbstractSampler{N}}
    Terrace{N,S}(x.random_state, x, sort(points), invert)
end

function sample(sampler::Terrace{N}, coords::Vararg{Real,N}) where {N}
    points = sampler.points
    len = length(points)
    x = sample(sampler.source, coords...)
    p = findfirst(≥(x), points)
    i = p !== nothing ? p : len
    i1 = clamp(i - 1, 1:len)
    i2 = clamp(i, 1:len)
    if i1 == i2
        points[i2]
    else
        in1 = points[i1]
        in2 = points[i2]
        t = (x - in1) / (in2 - in1)
        if sampler.invert
            t = 1 - t
            in1, in2 = in2, in1
        end
        lerp(in1, in2, t^2)
    end
end
