module Common

using FastPow
using Images: RGB
using ColorSchemes: grays, amp

include("sampler.jl")
include("shared.jl")
include("state_random.jl")
include("state_perlin.jl")
include("trait_hash.jl")
include("image.jl")

end
