function sample(sampler::Cellular{3,M,F}, x::T, y::T, z::T) where {M,F,T<:Real}
    seed = get_seed(sampler)
    table = sampler.table
    jitter = sampler.jitter * JITTER2
    r = round.(Int, (x, y, z)) .- 1
    xr, yr, zr = r .- (x, y, z)
    xp, yp_base, zp_base = r .* (PRIME_X, PRIME_Y, PRIME_Z)
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
            hash = (syp ⊻ zp) % UInt32
            vx = table[(hash+1)&1023] * jitter + xri
            vy = table[((hash|1)+1)&1023] * jitter + yr + yi
            for zi in 0:2
                vz = table[((hash|2)+1)&1023] * jitter + zr + zi
                d = cell_distance(M, vx, vy, vz)
                maxf = clamp(d, minf, maxf)
                if d < minf
                    minf = d
                    closest_hash = hash
                end
                zp += PRIME_Z
            end
            yp += PRIME_Y
        end
        xp += PRIME_X
    end
    cell_value(F, M, closest_hash, minf, maxf) - 1
end
