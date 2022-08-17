"""
    gen_image(sampler::AbstractSampler; <kwargs>)

Construct a 2-dimensional array of `Images.RGB` values, suitable for writing to disk as an image
file.

# Arguments

  - `sampler::AbstractSampler`: Any instance of a sampler. The sampler is sampled using each pixel
    coordinates as the X and Y input coordinates, and random Z and W coordinates for 3 and
    4-dimensional samplers.
  - `w::Integer=1024`: The width in pixels of the image array to generate.
  - `h::Integer=1024`: The height in pixels of the image array to generate.
  - `xbounds::NTuple{2,Float64}=(-1.0, 1.0)`: The bounds along the X axis to sample coordinates
    from. This remaps pixel coordinates to this range to be used for the input coordinates to sample
    with.
  - `ybounds::NTuple{2,Float64}=(-1.0, 1.0)`: The bounds along the Y axis to sample coordinates
    from. This remaps pixel coordinates to this range to be used for the input coordinates to sample
    with.
  - `colorscheme=nothing`: A `ColorSchemes.ColorScheme` object to colorize the image with, or
    `nothing`.
"""
function gen_image(
    sampler::S;
    w::Integer=1024,
    h::Integer=1024,
    xbounds::NTuple{2,Float64}=(-1.0, 1.0),
    ybounds::NTuple{2,Float64}=(-1.0, 1.0),
    colorscheme::Union{ColorScheme,Nothing}=nothing,
) where {N,S<:AbstractSampler{N}}
    x1, x2 = xbounds
    y1, y2 = ybounds
    xd = (x2 - x1) / w
    yd = (y2 - y1) / h
    img = Array{RGB}(undef, h, w)
    zw = rand(rng(sampler), Float64, N - 2) * 1000
    Threads.@threads for x in 1:w
        cx = x * xd + x1
        for y in 1:h
            cy = y * yd + y1
            value = sample(sampler, cx, cy, zw...) * 0.5 + 0.5
            img[y, x] = colorscheme !== nothing ? colorscheme[value] : value
        end
    end
    img
end
