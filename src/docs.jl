const _doc_seed = "An integer used to seed the random number generator for this sampler."

const doc_sampler = """
    sample(sampler::AbstractSampler, x::Real)
    sample(sampler::AbstractSampler, x::Real, y::Real)
    sample(sampler::AbstractSampler, x::Real, y::Real, z::Real)
    sample(sampler::AbstractSampler, x::Real, y::Real, z::Real, w::Real)

Sample from `sampler` with the supplied coordinates. The number of coordinates should match the
dimensionality of the sampler type.
"""

const doc_gen_image = """
    gen_image(sampler::AbstractSampler; kwargs...)

Construct a 2-dimensional array of `ColorTypes.RGB` values, suitable for writing to disk as an image
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

const doc_constant_1d = """
    constant_1d(; seed=0, value=0.0)

Construct a sampler that constantly outputs `value` each time it is sampled from.

This is useful for debugging and applications where you need to combine a constant value.

# Arguments

  - `seed`: $(_doc_seed)

  - `value`: A constant value to emit each time this sampler is sampled from.
"""

const doc_checkered_2d = """
    checkered_2d(; seed=0)

Construct a sampler that outputs values in a checkerboard-like pattern when it is sampled from.

That is, output values will only ever be -1.0 or 1.0.

# Arguments

  - `seed`: $(_doc_seed)
"""

const doc_cylinders_2d = """
    cylinders_2d(; seed=0, frequency=1.0)

Construct a sampler that outputs values that form a pattern representing concentric cylinders when
it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `frequency`: The frequency of the signal, which controls how small or large the cylinders
    are.
"""

const doc_spheres_3d = """
    spheres_3d(; seed=0, frequency=1.0)

Construct a sampler that outputs values that form a pattern representing concentric spheres when it
is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `frequency`: The frequency of the signal, which controls how small or large the spheres are.
"""

const doc_value_1d = """
    value_1d(; seed=0)

Construct a sampler that outputs 1-dimensonal value noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)
"""

const doc_value_2d = """
    value_2d(; seed=0)

Construct a sampler that outputs 2-dimensonal value noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)
"""

const doc_value_3d = """
    value_3d(; seed=0)

Construct a sampler that outputs 3-dimensonal value noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)
"""

const doc_value_4d = """
    value_4d(; seed=0)

Construct a sampler that outputs 4-dimensonal value noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)
"""

const doc_cubic_1d = """
    cubic_1d(; seed=0)

Construct a sampler that outputs 1-dimensonal cubic noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)
"""

const doc_cubic_2d = """
    cubic_2d(; seed=0)

Construct a sampler that outputs 2-dimensonal cubic noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)
"""

const doc_cubic_3d = """
    cubic_3d(; seed=0)

Construct a sampler that outputs 3-dimensonal cubic noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)
"""

const doc_cubic_4d = """
    cubic_4d(; seed=0)

Construct a sampler that outputs 4-dimensonal cubic noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)
"""

const doc_perlin_1d = """
    perlin_1d(; seed=0)
    perlin_improved_1d(; seed=0)

Construct a sampler that outputs 1-dimensional Perlin "Improved" noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

# Notes:

  - `perlin_improved_1d` is deprecated and will be removed in favor of `perlin_1d` in a future
    major version release.
"""

const doc_perlin_2d = """
    perlin_2d(; seed=0)
    perlin_improved_2d(; seed=0)

Construct a sampler that outputs 2-dimensional Perlin "Improved" noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

# Notes:

  - `perlin_improved_2d` is deprecated and will be removed in favor of `perlin_2d` in a future
    major version release.
"""

const doc_perlin_3d = """
    perlin_3d(; seed=0)
    perlin_improved_3d(; seed=0)

Construct a sampler that outputs 3-dimensional Perlin "Improved" noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

# Notes:

  - `perlin_improved_3d` is deprecated and will be removed in favor of `perlin_3d` in a future
    major version release.
"""

const doc_perlin_4d = """
    perlin_4d(; seed=0)
    perlin_improved_4d(; seed=0)

Construct a sampler that outputs 4-dimensional Perlin "Improved" noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

# Notes:

  - `perlin_improved_4d` is deprecated and will be removed in favor of `perlin_4d` in a future
    major version release.
"""

const doc_simplex_1d = """
    simplex_1d(; seed=0)

Construct a sampler that outputs 1-dimensional Perlin Simplex noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)
"""

const doc_simplex_2d = """
    simplex_2d(; seed=0)

Construct a sampler that outputs 2-dimensional Perlin Simplex noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)
"""

const doc_simplex_3d = """
    simplex_3d(; seed=0, smooth=false)

Construct a sampler that outputs 3-dimensional Perlin Simplex noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `smooth`: Specify whether to have continuous gradients.
    Simplex variants, even the original Simplex noise by Ken Perlin, overshoot the radial extent for
    the signal reconstruction kernel in order to improve the visual of the noise. Normally this is
    okay, especially if layering multiple octaves of the noise. However, in some applications, such
    as creating height or bump maps, this will produce discontinuities visually identified by
    jarring creases in the generated noise.

    This option changes the falloff in order to produce smooth continuous noise, however, the
    resulting noise may look quite different than the non-smooth option, depending on the Simplex
    variant.

    The default value is `false`, in order to be true to the original implementation.
"""

const doc_simplex_4d = """
    simplex_4d(; seed=0, smooth=false)

Construct a sampler that outputs 4-dimensional Perlin Simplex noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `smooth`: Specify whether to have continuous gradients.
    Simplex variants, even the original Simplex noise by Ken Perlin, overshoot the radial extent for
    the signal reconstruction kernel in order to improve the visual of the noise. Normally this is
    okay, especially if layering multiple octaves of the noise. However, in some applications, such
    as creating height or bump maps, this will produce discontinuities visually identified by
    jarring creases in the generated noise.

    This option changes the falloff in order to produce smooth continuous noise, however, the
    resulting noise may look quite different than the non-smooth option, depending on the Simplex
    variant.

    The default value is `false`, in order to be true to the original implementation.
"""

const doc_opensimplex_2d = """
    opensimplex_2d(; seed=0, smooth=false)

Construct a sampler that outputs 2-dimensional legacy OpenSimplex noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `smooth`: Specify whether to have continuous gradients.
    Simplex variants, even the original Simplex noise by Ken Perlin, overshoot the radial extent for
    the signal reconstruction kernel in order to improve the visual of the noise. Normally this is
    okay, especially if layering multiple octaves of the noise. However, in some applications, such
    as creating height or bump maps, this will produce discontinuities visually identified by
    jarring creases in the generated noise.

    This option changes the falloff in order to produce smooth continuous noise, however, the
    resulting noise may look quite different than the non-smooth option, depending on the Simplex
    variant.

    The default value is `false`, in order to be true to the original implementation.

# See also:

  - [`opensimplex2_2d`](@ref opensimplex2_2d)
  - [`openSimplex2s_2d`](@ref opensimplex2s_2d)
"""

const doc_opensimplex_3d = """
    opensimplex_3d(; seed=0, smooth=false)

Construct a sampler that outputs 3-dimensional legacy OpenSimplex noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `smooth`: Specify whether to have continuous gradients.
    Simplex variants, even the original Simplex noise by Ken Perlin, overshoot the radial extent for
    the signal reconstruction kernel in order to improve the visual of the noise. Normally this is
    okay, especially if layering multiple octaves of the noise. However, in some applications, such
    as creating height or bump maps, this will produce discontinuities visually identified by
    jarring creases in the generated noise.

    This option changes the falloff in order to produce smooth continuous noise, however, the
    resulting noise may look quite different than the non-smooth option, depending on the Simplex
    variant.

    The default value is `false`, in order to be true to the original implementation.

# See also:

  - [`opensimplex2_3d`](@ref opensimplex2_3d)
  - [`opensimplex2s_3d`](@ref opensimplex2s_3d)
"""

const doc_opensimplex_4d = """
    opensimplex_4d(; seed=0, smooth=false)

Construct a sampler that outputs 4-dimensional legacy OpenSimplex noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `smooth`: Specify whether to have continuous gradients.
    Simplex variants, even the original Simplex noise by Ken Perlin, overshoot the radial extent for
    the signal reconstruction kernel in order to improve the visual of the noise. Normally this is
    okay, especially if layering multiple octaves of the noise. However, in some applications, such
    as creating height or bump maps, this will produce discontinuities visually identified by
    jarring creases in the generated noise.

    This option changes the falloff in order to produce smooth continuous noise, however, the
    resulting noise may look quite different than the non-smooth option, depending on the Simplex
    variant.

    The default value is `false`, in order to be true to the original implementation.

# See also:

  - [`opensimplex2_4d`](@ref opensimplex2_4d)
  - [`opensimplex2s_4d`](@ref opensimplex2s_4d)
"""

const doc_opensimplex2_2d = """
    opensimplex2_2d(; seed=0, orient=nothing)

Construct a sampler that outputs 2-dimensional OpenSimplex2 noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `orient`: Either the symbol `:x` or the value `nothing`:

      + `nothing`: Use the standard orientation.
      + `:x`: The noise space will be re-oriented with the Y axis pointing down the main diagonal to
        improve visual isotropy.
"""

const doc_opensimplex2_3d = """
    opensimplex2_3d(; seed=0, smooth=false, orient=nothing)

Construct a sampler that outputs 3-dimensional OpenSimplex2 noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `smooth`: Specify whether to have continuous gradients.
    Simplex variants, even the original Simplex noise by Ken Perlin, overshoot the radial extent for
    the signal reconstruction kernel in order to improve the visual of the noise. Normally this is
    okay, especially if layering multiple octaves of the noise. However, in some applications, such
    as creating height or bump maps, this will produce discontinuities visually identified by
    jarring creases in the generated noise.

    This option changes the falloff in order to produce smooth continuous noise, however, the
    resulting noise may look quite different than the non-smooth option, depending on the Simplex
    variant.

    The default value is `false`, in order to be true to the original implementation.

  - `orient`: Either the symbol `:x` or the value `nothing`:

      + `nothing`: Use the standard orientation.
      + `:x`: The noise space will be re-oriented with the Y axis pointing down the main diagonal to
        improve visual isotropy.
      + `:xy`: Re-orient the noise space to have better visual isotropy in the XY plane.
      + `:xz`: Re-orient the noise space to have better visual isotropy in the XZ plane.
"""

const doc_opensimplex2_4d = """
    opensimplex2_4d(; seed=0, smooth=false, orient=nothing)

Construct a sampler that outputs 4-dimensional OpenSimplex2 noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `smooth`: Specify whether to have continuous gradients.
    Simplex variants, even the original Simplex noise by Ken Perlin, overshoot the radial extent for
    the signal reconstruction kernel in order to improve the visual of the noise. Normally this is
    okay, especially if layering multiple octaves of the noise. However, in some applications, such
    as creating height or bump maps, this will produce discontinuities visually identified by
    jarring creases in the generated noise.

    This option changes the falloff in order to produce smooth continuous noise, however, the
    resulting noise may look quite different than the non-smooth option, depending on the Simplex
    variant.

    The default value is `false`, in order to be true to the original implementation.

  - `orient`: Either the symbol `:x` or the value `nothing`:

      + `nothing`: Use the standard orientation.
      + `:x`: The noise space will be re-oriented with the Y axis pointing down the main diagonal to
        improve visual isotropy.
      + `:xy`: Re-orient the noise space to have better visual isotropy in the XY plane.
      + `:xz`: Re-orient the noise space to have better visual isotropy in the XZ plane.
      + `:xyz`: Re-orient the noise space to be better suited for time-varied animations, where
        the W axis is time.
"""

const doc_opensimplex2s_2d = """
    opensimplex2s_2d(; seed=0, orient=nothing)

Construct a sampler that outputs 2-dimensional OpenSimplex2S noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `orient`: Either the symbol `:x` or the value `nothing`:

      + `nothing`: Use the standard orientation.
      + `:x`: The noise space will be re-oriented with the Y axis pointing down the main diagonal to
        improve visual isotropy.
"""

const doc_opensimplex2s_3d = """
    opensimplex2s_3d(; seed=0, smooth=false, orient=nothing)

Construct a sampler that outputs 3-dimensional OpenSimplex2S noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `smooth`: Specify whether to have continuous gradients.
    Simplex variants, even the original Simplex noise by Ken Perlin, overshoot the radial extent for
    the signal reconstruction kernel in order to improve the visual of the noise. Normally this is
    okay, especially if layering multiple octaves of the noise. However, in some applications, such
    as creating height or bump maps, this will produce discontinuities visually identified by
    jarring creases in the generated noise.

    This option changes the falloff in order to produce smooth continuous noise, however, the
    resulting noise may look quite different than the non-smooth option, depending on the Simplex
    variant.

    The default value is `false`, in order to be true to the original implementation.

  - `orient`: Either the symbol `:x` or the value `nothing`:

      + `nothing`: Use the standard orientation.
      + `:x`: The noise space will be re-oriented with the Y axis pointing down the main diagonal to
        improve visual isotropy.
      + `:xy`: Re-orient the noise space to have better visual isotropy in the XY plane.
      + `:xz`: Re-orient the noise space to have better visual isotropy in the XZ plane.
"""

const doc_opensimplex2s_4d = """
    opensimplex2s_4d(; seed=0, smooth=false, orient=nothing)

Construct a sampler that outputs 4-dimensional OpenSimplex2S noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `smooth`: Specify whether to have continuous gradients.
    Simplex variants, even the original Simplex noise by Ken Perlin, overshoot the radial extent for
    the signal reconstruction kernel in order to improve the visual of the noise. Normally this is
    okay, especially if layering multiple octaves of the noise. However, in some applications, such
    as creating height or bump maps, this will produce discontinuities visually identified by
    jarring creases in the generated noise.

    This option changes the falloff in order to produce smooth continuous noise, however, the
    resulting noise may look quite different than the non-smooth option, depending on the Simplex
    variant.

    The default value is `false`, in order to be true to the original implementation.

  - `orient`: Either the symbol `:x` or the value `nothing`:

      + `nothing`: Use the standard orientation.
      + `:x`: The noise space will be re-oriented with the Y axis pointing down the main diagonal to
        improve visual isotropy.
      + `:xy`: Re-orient the noise space to have better visual isotropy in the XY plane.
      + `:xz`: Re-orient the noise space to have better visual isotropy in the XZ plane.
      + `:xyz`: Re-orient the noise space to be better suited for time-varied animations, where
        the W axis is time.
"""

const doc_worley_1d = """
    worley_1d(; seed=0, jitter=1.0, output=:f1)

Construct a sampler that outputs 1-dimensional Worley noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `jitter`: A `Real` number between 0.0 and 1.0, with values closer to one randomly distributing
    cells away from their grid alignment.

  - `output`: One of the following symbols:
      + `:f1`: Calculate the distance to the nearest cell as the output.
      + `:f2`: Calculate the distance to the second-nearest cell as the output.
      + `:+`: Calculate `:f1` + `:f2` as the output.
      + `:-`: Calculate `:f2` - `:f1` as the output.
      + `:*`: Calculate `:f1` * `:f2` as the output.
      + `:/`: Calculate `:f1` / `:f2` as the output.
      + `:value`: Use the cell's hash value as the output.
"""

const doc_worley_2d = """
    worley_2d(; seed=0, jitter=1.0, output=:f1, metric=:euclidean)

Construct a sampler that outputs 2-dimensional Worley noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `jitter`: A `Real` number between 0.0 and 1.0, with values closer to one randomly distributing
    cells away from their grid alignment.

  - `output`: One of the following symbols:
      + `:f1`: Calculate the distance to the nearest cell as the output.
      + `:f2`: Calculate the distance to the second-nearest cell as the output.
      + `:+`: Calculate `:f1` + `:f2` as the output.
      + `:-`: Calculate `:f2` - `:f1` as the output.
      + `:*`: Calculate `:f1` * `:f2` as the output.
      + `:/`: Calculate `:f1` / `:f2` as the output.
      + `:value`: Use the cell's hash value as the output.

  - `metric`: One of the following symbols:
      + `:manhattan`: Use the Manhattan distance to the next cell (Minkowski metric ``p = 2^0``).
      + `:euclidean`: Use the Euclidean distance to the next cell (Minkowski metric ``p = 2^1``).
      + `:euclidean²`: Same as `:euclidean` but slighter faster due to no ``\\sqrt{}``.
      + `:minkowski4`: Use Minkowski metric with ``p = 2^4`` for the distance to the next cell.
      + `:chebyshev`: Use the Chebyshev distance to the next cell (Minkowski metric ``p =
        2^\\infty``).
"""

const doc_worley_3d = """
    worley_3d(; seed=0, jitter=1.0, output=:f1, metric=:euclidean)

Construct a sampler that outputs 3-dimensional Worley noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `jitter`: A `Real` number between 0.0 and 1.0, with values closer to one randomly distributing
    cells away from their grid alignment.

  - `output`: One of the following symbols:
      + `:f1`: Calculate the distance to the nearest cell as the output.
      + `:f2`: Calculate the distance to the second-nearest cell as the output.
      + `:+`: Calculate `:f1` + `:f2` as the output.
      + `:-`: Calculate `:f2` - `:f1` as the output.
      + `:*`: Calculate `:f1` * `:f2` as the output.
      + `:/`: Calculate `:f1` / `:f2` as the output.
      + `:value`: Use the cell's hash value as the output.

  - `metric`: One of the following symbols:
      + `:manhattan`: Use the Manhattan distance to the next cell (Minkowski metric ``p = 2^0``).
      + `:euclidean`: Use the Euclidean distance to the next cell (Minkowski metric ``p = 2^1``).
      + `:euclidean²`: Same as `:euclidean` but slighter faster due to no ``\\sqrt{}``.
      + `:minkowski4`: Use Minkowski metric with ``p = 2^4`` for the distance to the next cell.
      + `:chebyshev`: Use the Chebyshev distance to the next cell (Minkowski metric ``p =
        2^\\infty``).
"""

const doc_worley_4d = """
    worley_4d(; seed=0, jitter=1.0, output=:f1, metric=:euclidean)

Construct a sampler that outputs 4-dimensional Worley noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `jitter`: A `Real` number between 0.0 and 1.0, with values closer to one randomly distributing
    cells away from their grid alignment.

  - `output`: One of the following symbols:
      + `:f1`: Calculate the distance to the nearest cell as the output.
      + `:f2`: Calculate the distance to the second-nearest cell as the output.
      + `:+`: Calculate `:f1` + `:f2` as the output.
      + `:-`: Calculate `:f2` - `:f1` as the output.
      + `:*`: Calculate `:f1` * `:f2` as the output.
      + `:/`: Calculate `:f1` / `:f2` as the output.
      + `:value`: Use the cell's hash value as the output.

  - `metric`: One of the following symbols:
      + `:manhattan`: Use the Manhattan distance to the next cell (Minkowski metric ``p = 2^0``).
      + `:euclidean`: Use the Euclidean distance to the next cell (Minkowski metric ``p = 2^1``).
      + `:euclidean²`: Same as `:euclidean` but slighter faster due to no ``\\sqrt{}``.
      + `:minkowski4`: Use Minkowski metric with ``p = 2^4`` for the distance to the next cell.
      + `:chebyshev`: Use the Chebyshev distance to the next cell (Minkowski metric ``p =
        2^\\infty``).
"""

const doc_fbm_fractal_1d = """
    fbm_fractal_1d(; kwargs...)

Construct a sampler that outputs a 1-dimensional fractional Brownian motion fractal noise when it
is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `source::AbstractSampler=simplex_1d()`: A 1-dimensional sampler instance to use as the source of
    the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""

const doc_fbm_fractal_2d = """
    fbm_fractal_2d(; kwargs...)

Construct a sampler that outputs a 2-dimensional fractional Brownian motion fractal noise when it
is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `source::AbstractSampler=opensimplex2_2d()`: A 2-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""

const doc_fbm_fractal_3d = """
    fbm_fractal_3d(; kwargs...)

Construct a sampler that outputs a 3-dimensional fractional Brownian motion fractal noise when it
is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `source::AbstractSampler=opensimplex2_3d()`: A 3-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""

const doc_fbm_fractal_4d = """
    fbm_fractal_4d(; kwargs...)

Construct a sampler that outputs a 4-dimensional fractional Brownian motion fractal noise when it
is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `source::AbstractSampler=opensimplex2_4d()`: A 4-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""

const doc_billow_fractal_1d = """
    billow_fractal_1d(; kwargs...)

Construct a sampler that outputs a 1-dimensional billow fractal noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `source::AbstractSampler=simplex_1d()`: A 1-dimensional sampler instance to use as the source of
    the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""

const doc_billow_fractal_2d = """
    billow_fractal_2d(; kwargs...)

Construct a sampler that outputs a 2-dimensional billow fractal noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `source::AbstractSampler=opensimplex2s_2d()`: A 2-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""

const doc_billow_fractal_3d = """
    billow_fractal_3d(; kwargs...)

Construct a sampler that outputs a 3-dimensional billow fractal noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `source::AbstractSampler=opensimplex2s_3d()`: A 3-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""

const doc_billow_fractal_4d = """
    billow_fractal_4d(; kwargs...)

Construct a sampler that outputs a 4-dimensional billow fractal noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `source::AbstractSampler=opensimplex2s_4d()`: A 4-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""

const doc_multi_fractal_1d = """
    multi_fractal_1d(; kwargs...)

Construct a sampler that outputs a 1-dimensional multifractal noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `source::AbstractSampler=simplex_1d()`: A 1-dimensional sampler instance to use as the source of
    the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""

const doc_multi_fractal_2d = """
    multi_fractal_2d(; kwargs...)

Construct a sampler that outputs a 2-dimensional multifractal noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `source::AbstractSampler=opensimplex2s_2d()`: A 2-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""

const doc_multi_fractal_3d = """
    multi_fractal_3d(; kwargs...)

Construct a sampler that outputs a 3-dimensional multifractal noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `source::AbstractSampler=opensimplex2s_3d()`: A 3-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""

const doc_multi_fractal_4d = """
    multi_fractal_4d(; kwargs...)

Construct a sampler that outputs a 4-dimensional multifractal noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `source::AbstractSampler=opensimplex2s_4d()`: A 4-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.5`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""

const doc_hybrid_fractal_1d = """
    hybrid_fractal_1d(; kwargs...)

Construct a sampler that outputs a 1-dimensional hybrid multifractal noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `source::AbstractSampler=simplex_1d()`: A 1-dimensional sampler instance to use as the source of
    the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.25`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""

const doc_hybrid_fractal_2d = """
    hybrid_fractal_2d(; kwargs...)

Construct a sampler that outputs a 2-dimensional hybrid multifractal noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `source::AbstractSampler=opensimplex2s_2d()`: A 2-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.25`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""

const doc_hybrid_fractal_3d = """
    hybrid_fractal_3d(; kwargs...)

Construct a sampler that outputs a 3-dimensional hybrid multifractal noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `source::AbstractSampler=opensimplex2s_3d()`: A 3-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.25`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""

const doc_hybrid_fractal_4d = """
    hybrid_fractal_4d(; kwargs...)

Construct a sampler that outputs a 4-dimensional hybrid multifractal noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `source::AbstractSampler=opensimplex2s_4d()`: A 4-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=0.25`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.
"""

const doc_ridged_fractal_1d = """
    ridged_fractal_1d(; kwargs...)

Construct a sampler that outputs a 1-dimensional ridged multifractal noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `source::AbstractSampler=simplex_1d()`: A 1-dimensional sampler instance to use as the source of
    the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=1.0`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.

   - `attenuation=2.0`: The attenuation to apply to the weight of each octave.
"""

const doc_ridged_fractal_2d = """
    ridged_fractal_2d(; kwargs...)

Construct a sampler that outputs a 2-dimensional ridged multifractal noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `source::AbstractSampler=opensimplex2s_2d()`: A 2-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=1.0`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.

   - `attenuation=2.0`: The attenuation to apply to the weight of each octave.
"""

const doc_ridged_fractal_3d = """
    ridged_fractal_3d(; kwargs...)

Construct a sampler that outputs a 3-dimensional ridged multifractal noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `source::AbstractSampler=opensimplex2s_3d()`: A 3-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=1.0`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.

   - `attenuation=2.0`: The attenuation to apply to the weight of each octave.
"""

const doc_ridged_fractal_4d = """
    ridged_fractal_4d(; kwargs...)

Construct a sampler that outputs a 4-dimensional ridged multifractal noise when it is sampled from.

# Arguments

  - `seed`: $(_doc_seed)

  - `source::AbstractSampler=opensimplex2s_4d()`: A 4-dimensional sampler instance to use as the
    source of the fractal.

  - `octaves=4`: An integer between 1 and 32, denoting the number of octaves to apply.

  - `frequency=1.0`: The frequency of the first octave's signal.

  - `lacunarity=2.0`: A multiplier that determines how quickly the frequency increases for
    successive octaves.

  - `persistence=1.0`: A multiplier that determines how quickly the amplitude diminishes for
    successive octaves.

   - `attenuation=2.0`: The attenuation to apply to the weight of each octave.
"""

const doc_mod_abs = """
    abs(x::AbstractSampler)

Construct a modifier sampler that outputs the absolute value of its source when it is sampled from.
"""

const doc_mod_add = """
    +(x::AbstractSampler, y::AbstractSampler)

Construct a modifier sampler that outputs the sum of the outputs of samplers `x` and `y`.
"""

const doc_mod_add_scalar = """
    +(x::AbstractSampler, y::Real)

Construct a modifier sampler that outputs the sum of the output of sampler `x` and the scalar `y`.
"""

const doc_mod_cache = """
    cache(x::AbstractSampler)

Construct a modifier sampler that caches the set of the input coordinates and their corresponding
output value of its source sampler. If the input coordinates differs from the previously cached
output, the cache is invalidated and the new output is cached.

Caching is useful if a sampler is used as a source for multiple modifiers. Without caching, the
duplicated input sources would redundantly compute the same outputs, which would be expensive,
especially if long pipelines share a long subgraph.
"""

const doc_mod_clamp = """
    clamp(x::AbstractSampler, lo::AbstractSampler, hi::AbstractSampler)

Construct a modifier sampler that clamps the output of sampler `x` to be within the range of of
output values from samplers `lo` and `hi`.
"""

const doc_mod_clamp_scalar = """
    clamp(x::AbstractSampler, lo=-1.0, hi=1.0)

Construct a modifier sampler that clamps the output of sampler `x` to be within the range of of the
scalars `lo` and `hi`.
"""

const doc_mod_copysign = """
    copysign(x::AbstractSampler, y::AbstractSampler)

Construct a modifier sampler that outputs the value of sampling from `x` with the sign copied from
the value of sampling from `y`.
"""

const doc_mod_curve = """
    curve(x::AbstractSampler, points::Vector{Pair{Float64,Float64}})

Construct a modifier sampler that outputs the result of sampling from `x` after remapping its output
to an arbitrary curve.

The curve is defined by a `Vector` of `Pair`s given by `points`. Each pair of points represents an
input and output number. The curve is a cubic spline, and so `points` must contain a list of four
point pairs at a minimum. Additionally, no two point pairs can contain the same input point value.

When sampling from sampler `x`, the output is evaluated using the curve data, and maps it to a new
output value.
"""

const doc_mod_div = """
    /(x::AbstractSampler, y::AbstractSampler)

Construct a modifier sampler that performs division of the output of sampler `x` by the output of
sampler `y`.
"""

const doc_mod_div_scalar = """
    /(x::AbstractSampler, y::Real)

Construct a modifier sampler that performs division of the output of sampler `x` by the scalar `y`.
"""

const doc_mod_max = """
    max(x::AbstractSampler, y::AbstractSampler)

Construct a modifier sampler that outputs the maximum value of the outputs of sampling from samplers
`x` and `y`.
"""

const doc_mod_min = """
    min(x::AbstractSampler, y::AbstractSampler)

Construct a modifier sampler that outputs the minimum value of the outputs of sampling from samplers
`x` and `y`.
"""

const doc_mod_mix = """
    mix(x::AbstractSampler, y::AbstractSampler, t::AbstractSampler)

Construct a modifier sampler that outputs the result of linearly interpolating the output of
samplers `x` and `y` by the output of sampler `t`.
"""

const doc_mod_mix_scalar = """
    mix(x::AbstractSampler, y::AbstractSampler, t::Real)

Construct a modifier sampler that outputs the result of linearly interpolating the output of
samplers `x` and `y` by the scalar `t`.
"""

const doc_mod_mul = """
    *(x::AbstractSampler, y::AbstractSampler)

Construct a modifier sampler that outputs the product of the outputs samplers `x` and `y`.
"""

const doc_mod_mul_scalar = """
    *(x::AbstractSampler, y::Real)

Construct a modifier sampler that outputs the product of the output of sampler `x` by the scalar
`y`.
"""

const doc_mod_muladd = """
    muladd(x::AbstractSampler, strength=1.0, bias=0.0)

Construct a modifier sampler that performs multiplies the output of sampler `x` by the scalar
`strength`, followed by adding the scalar `bias`. sampler `y`.
"""

const doc_mod_pow = """
    ^(x::AbstractSampler, y::Real)

Construct a modifier sampler that raises the output of sampler `x` to the power of the scalar `y`.
"""

const doc_mod_rotate = """
    rotate(source::AbstractSampler; x=0.0, y=0.0, z=0.0)

Construct a modifier sampler that rotates the input coordinates of sampler `source` around the
origin before sampling from it.

The coordinate system is assumed to be left-handed.

The angle of rotation is specified in radians for the corresponding axis given by `x`, `y`, and `z`.
"""

const doc_mod_scale = """
    scale(source::AbstractSampler; x=1.0, y=1.0, z=1.0, w=1.0)

Construct a modifier sampler that scales the input coordinates of sampler `source` before sampling
from it.

Each axis can be scaled independently with `x`, `y`, `z`, or `w`.
"""

const doc_mod_scale_scalar = """
    scale(source::AbstractSampler, scale::Real)

Construct a modifier sampler that uniformly scales the input coordinates of sampler `source` by the
scalar `scale` before sampling from it.
"""

const doc_mod_select = """
    select(x::AbstractSampler, y::AbstractSampler; kwargs...)

Construct a modifier sampler that outputs either the out of sampler `x` or `y`, depending on the
output of sampler `z`.

If the output of sampler `z` is within the range denoted by `min` and `max`, the output of sampler
`y` is chosen. If the output of sampler `z` is outside of this range, the output of sampler `x` is
chosen.

# Arguments

  - `min`: A real number between -1.0 and 1.0 defining the lower bound of the selection range.

  - `max`: A real number between -1.0 and 1.0 defining the upper bound of the selection range.

  - `falloff`: A real number between 0.0 and 1.0 specifying the smoothness of the transition.
"""

const doc_mod_sub = """
    -(x::AbstractSampler, y::AbstractSampler)

Construct a modifier sampler that outputs the difference of the outputs of samplers `x` and `y`.
"""

const doc_mod_sub_scalar = """
    -(x::AbstractSampler, y::Real)

Construct a modifier sampler that outputs the difference of the output of sampler `x` and the scalar
`y`.
"""

const doc_mod_sub_unary = """
    -(x::AbstractSampler)

Construct a modifier sampler that outputs the negated output of sampler `x`.
"""

const doc_mod_terrace = """
    terrace(x::AbstractSampler, points::Vector{Pair{Float64,Float64}}; invert=false)

Construct a modifier sampler that outputs the result of sampling from `x` after remapping its output
to a terrace-forming curve.

The curve is defined by a `Vector` of `Float64`s given by `points`. Each point represents an input
and output number.

When sampling from sampler `x`, the output is evaluated using the curve data, and maps it to a new
output value.

# Arguments

  - `invert`: Specify whether the curve is inverted between control points.
"""

const doc_mod_translate = """
    translate(source::AbstractSampler; x=0.0, y=0.0, z=0.0, w=0.0)

Construct a modifier sampler that translates the input coordinates of sampler `source` along a
specified vector, with a direction and magnitude given by the coordinates `x`, `y`, `z`, and `w`.
"""

const doc_mod_turbulence = """
    turbulence(s1::AbstractSampler; s2::AbstractSampler; kwargs...)

Construct a modifier sampler that displaces the input coordinates of sampler `s1` by the output of
sampler `s2` with a fractional Brownian motion fractal applied to it.

Sampler `s2`'s input coordinates are randomly generated using the seed of sampler `s1`.

# Arguments

  - `frequency=1.0`: The frequency of the fractal signal to apply to sampler `s2`.

  - `roughness=3`: The number of octaves of the fractal to apply to sampler `s2`.

  - `power=1.0`: A scaling factor that is applied to the displaced result before sampling from
    sampler `s1`.
"""

const doc_mod_warp = """
    warp(source::AbstractSampler; kwargs...)

Construct a modifier sampler that performs domain warping of the sampler `source` before sampling
from it.

Domain warping feeds the output of other samplers to the input of a sampler. For this modifier, each
input coordinate can specify a different sampler to warp with.

If a sampler is not supplied for `x`, `y`, `z`, or `w`, a sampler that outputs a constant zero will
be used instead.

# Arguments

  - `x::AbstractSampler=constant_1d()`: A sampler to warp the X axis by.

  - `y::AbstractSampler=constant_1d()`: A sampler to warp the Y axis by.

  - `z::AbstractSampler=constant_1d()`: A sampler to warp the Z axis by.

  - `w::AbstractSampler=constant_1d()`: A sampler to warp the W axis by.
"""
