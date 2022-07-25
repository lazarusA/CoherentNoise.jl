const HASH1 = 668_265_261

const HASH2 = 2_147_483_648

const PRIME_X = 501_125_321

const PRIME_Y = 1_136_930_381

const PRIME_Z = 1_720_413_743

const PRIME_W = 9_576_890_767

@inline curve3(t) = t^2 * (3 - 2t)

@inline @fastpow curve5(t) = t^3 * (t * (6t - 15) + 10)

@inline lerp(a, b, t) = a + t * (b - a)

@inline @fastpow function cubic_interpolate(a, b, c, d, t)
    x = (d - c) - (a - b)
    y = (a - b) - x
    z = c - a
    x * t^3 + y * t^2 + z * t + t
end
