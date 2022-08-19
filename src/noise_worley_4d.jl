"""
    worley_4d(; kwargs...)

Construct a sampler that outputs 4-dimensional Worley noise when it is sampled from.

# Arguments

  - `seed=nothing`: An integer used to seed the random number generator for this sampler, or
    `nothing`. If a seed is not supplied, one will be generated automatically which will negatively
    affect reproducibility.

  - `jitter=1.0`: A `Real` number between 0.0 and 1.0, with values closer to one randomly
    distributing cells away from their grid alignment.

  - `metric=:euclidean`: One of the following symbols:

      + `:manhattan`: Use the Manhattan distance to the next cell (Minkowski metric ``p = 2^0``).

      + `:euclidean`: Use the Euclidean distance to the next cell (Minkowski metric ``p = 2^1``).

      + `:euclidean²`: Same as `:euclidean` but slighter faster due to no ``\\sqrt{}``.

      + `:minkowski4`: Use Minkowski metric with ``p = 2^4`` for the distance to the next cell.

      + `:chebyshev`: Use the Chebyshev distance to the next cell (Minkowski metric ``p =
        2^\\infty``).

  - `output=:f1`: One of the following symbols:

      + `:f1`: Calculate the distance to the nearest cell as the output.

      + `:f2`: Calculate the distance to the second-nearest cell as the output.

      + `:+`: Calculate `:f1` + `:f2` as the output.

      + `:-`: Calculate `:f2` - `:f1` as the output.

      + `:*`: Calculate `:f1` * `:f2` as the output.

      + `:/`: Calculate `:f1` / `:f2` as the output.

      + `:value`: Use the cell's hash value as the output.
"""
function worley_4d(; seed=nothing, metric=:euclidean, output=:f1, jitter=1.0)
    worley(4, seed, metric, output, jitter)
end

function sample(sampler::Worley{4,M,F}, x::T, y::T, z::T, w::T) where {M,F,T<:Real}
    seed = get_seed(sampler)
    table = sampler.table
    jitter = sampler.jitter * JITTER2
    r = round.(Int, (x, y, z, w)) .- 1
    xr, yr, zr, wr = r .- (x, y, z, w)
    xp, yp_base, zp_base, wp_base = r .* (PRIME_X, PRIME_Y, PRIME_Z, PRIME_W)
    minf = floatmax(Float64)
    maxf = minf
    closest_hash::UInt32 = 0
    @inbounds for xi in 0:2
        xri = xr + xi
        yp = yp_base
        sxp = seed ⊻ xp * HASH1
        for yi in 0:2
            zp = zp_base
            syp = sxp ⊻ yp
            szp = syp ⊻ zp
            wp = wp_base
            hash = (szp ⊻ wp) % UInt32
            vx = table[(hash+1)&2047] * jitter + xri
            vy = table[((hash|1)+1)&2047] * jitter + yr + yi
            for zi in 0:2
                vz = table[((hash|2)+1)&2047] * jitter + zr + zi
                for wi in 0:2
                    vw = table[((hash|3)+1)&2047] * jitter + wr + wi
                    d = cell_distance(M, vx, vy, vz, vw)
                    maxf = clamp(d, minf, maxf)
                    if d < minf
                        minf = d
                        closest_hash = hash
                    end
                    wp += PRIME_W
                end
                zp += PRIME_Z
            end
            yp += PRIME_Y
        end
        xp += PRIME_X
    end
    cell_value(F, M, closest_hash, minf, maxf) - 1
end
