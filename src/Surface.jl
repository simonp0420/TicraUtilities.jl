struct Surface
    text::String
    x::AbstractRange
    y::AbstractRange
    z::Matrix{Float64}
end

Surface() = Surface("", 0:0, 0:0, zeros(0, 0))

import Base.show
function show(io::IO, t::Surface)
    println(io, "Surface")
    println(io, "  text: ", t.text)
    println(io, "     x: ", t.x)
    println(io, "     y: ", t.y)
    println(io, "     z: ", summary(t.z))
    return nothing
end

get_x(c::Surface) = c.x
get_y(c::Surface) = c.y
get_z(c::Surface) = c.z
get_text(c::Surface) = c.text

function read_surface(fname::AbstractString)
    # Read a Ticra .sfc file.
    sfc = open(fname, "r") do fid
        text = strip(readline(fid))

        # Obtain the grid limits:
        (xs, ys, xe, ye) = parse.(Float64, split(readline(fid), t -> isspace(t) || t == ','; keepempty=false, limit=4))
        (nx, ny) = parse.(Int, split(readline(fid), t -> isspace(t) || t == ','; keepempty=false, limit=2))
        x = range(xs, xe, nx)
        y = range(ys, ye, ny)
        z = reshape(parse.(Float64, split(read(fid, String), limit=nx * ny)), ny, nx) |> permutedims # size(sfc.z) = (Nx,Ny)
        sfc = Surface(text, x, y, z)
        return sfc
    end
    return sfc
end

import Base.+
function +(s1::Surface, s2::Surface)
    (s1.x ≈ s2.x && s1.y ≈ s2.y) || error("Surfaces not defined on the same points")
    z = s1.z + s2.z
    text = string(s1.text, " + ", s2.text)
    return Surface(text, s1.x, s1.y, z)
end

import Base.-
function -(s1::Surface, s2::Surface)
    (s1.x ≈ s2.x && s1.y ≈ s2.y) || error("Surfaces not defined on the same points")
    z = s1.z - s2.z
    text = string(s1.text, " - ", s2.text)
    return Surface(text, s1.x, s1.y, z)
end

function write_surface(fname::AbstractString, sfc::Surface)
    # Write a Ticra .sfc file.
    open(fname, "w") do fid
        println(fid, sfc.text)
        @printf(fid, " %17.10e %17.10e %17.10e %17.10e\n", sfc.x[1], sfc.y[1], sfc.x[end], sfc.y[end])
        @printf(fid, "%d %d\n", length(sfc.x), length(sfc.y))
        for (i, z) in enumerate(permutedims(sfc.z))
            @printf(fid, " %17.10e", z)
            iszero(i % 4) && println(fid, "")
        end
    end
    return nothing
end

