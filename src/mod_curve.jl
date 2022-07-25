struct Curve{N,S} <: ModifierSampler{N}
    random_state::RandomState
    source::S
    points::Vector{Pair{Float64,Float64}}
end

"""
    curve(x::AbstractSampler, points::Vector{Pair{Float64,Float64}})

Construct a modifier sampler that outputs the result of sampling from `x` after remapping its output
to an arbitrary curve.

The curve is defined by a `Vector` of `Pair`s given by `points`. Each pair of points represents an
input and output number. The curve is a cubic spline, and so `points` must contain a list of four
point pairs at a minimum. Additionally, no two point pairs can contain the same input point value.

When sampling from sampler `x`, the output is evaluated using the curve data, and maps it to a new
output value.
"""
function curve(x::S, points::Vector{Pair{Float64,Float64}}) where {N,S<:AbstractSampler{N}}
    Curve{N,S}(random_state(x), x, sort(points, by=(p) -> p.first))
end

function sample(sampler::Curve{N}, coords::Vararg{Real,N}) where {N}
    points = sampler.points
    len = length(points)
    x = sample(sampler.source, coords...)
    p = findfirst(>(x) ∘ first, points)
    i = p !== nothing ? clamp(p, 3:len) : len
    i1 = clamp(i - 2, 1:len)
    i2 = clamp(i - 1, 1:len)
    i3 = clamp(i, 1:len)
    i4 = clamp(i + 1, 1:len)
    if i2 == i3
        points[i2].second
    else
        in1 = points[i2].first
        in2 = points[i3].first
        a, b, c, d = (x -> x.second).((points[i1], points[i2], points[i3], points[i4]))
        cubic_interpolate(a, b, c, d, x - in1 / in2 - in1)
    end
end
