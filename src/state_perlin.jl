using CircularArrays: CircularVector

struct PerlinState
    table::CircularVector{UInt8,Vector{UInt8}}
    PerlinState(rs) = new(shuffle(rs.rng, 0x00:0xff) |> CircularVector)
end
