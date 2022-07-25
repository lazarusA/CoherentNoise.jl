function sample(sampler::Cellular{2,M,F}, x::T, y::T) where {M,F,T<:Real}
    seed = get_seed(sampler)
    table = sampler.table
    jitter = sampler.jitter * JITTER1
    r = round.(Int, (x, y)) .- 1
    xr, yr = r .- (x, y)
    xp, yp_base = r .* (PRIME_X, PRIME_Y)
    minf = floatmax(Float64)
    maxf = minf
    closest_hash::UInt32 = 0
    @inbounds for xi in 0:2
        xri = xr + xi
        yp = yp_base
        sxp = seed ⊻ xp * HASH1
        for yi in 0:2
            hash = (sxp ⊻ yp) % UInt32
            vx = table[(hash+1)&511] * jitter + xri
            vy = table[((hash|1)+1)&511] * jitter + yr + yi
            d = cell_distance(M, vx, vy)
            maxf = clamp(d, minf, maxf)
            if d < minf
                minf = d
                closest_hash = hash
            end
            yp += PRIME_Y
        end
        xp += PRIME_X
    end
    cell_value(F, M, closest_hash, minf, maxf) - 1
end
