"""
    Tfile

Ticra Field-on-stations ("T-file") data.
"""
struct Tfile
    title::String  # title string from first line of file
    nel::Int          # Number of elements in the array (or number of beams)
    nsta::Vector{Int}  # Number of stations in each partition
    npart::Int          # Number of partitions
    p1::Vector{Matrix{ComplexF64}}  # 1st pol component, element size = (nsta[kp], Nel)
    p2::Vector{Matrix{ComplexF64}}  # 2cd pol component, element size = (nsta[kp], Nel)
end

import Base.show
function show(io::IO, t::Tfile)
    println(io, "Tfile")
    println(io, "  title: ", t.title)
    println(io, "  nel\t$(t.nel)")
    for kp in 1:t.npart
        println(io, "  Partition $(lpad(kp,3)): $(t.nsta[kp]) stations")
    end
    return nothing
end


"""
    read_tfile(filename)

Return a `Tfile` object containing data read from file `filename`.
"""
function read_tfile(filename)
    open(filename, "r") do fid
        title = readline(fid)
        tline = split(readline(fid))
        (nel, npart) = (parse(Int, t) for t in tline)
        nsta = zeros(Int, npart)
        p1 = Matrix{ComplexF64}[]
        p2 = Matrix{ComplexF64}[]
        for kp = 1:npart
            nsta[kp] = parse(Int, readline(fid))
            p1kp = Matrix{ComplexF64}(undef, nsta[kp], nel)
            p2kp = Matrix{ComplexF64}(undef, nsta[kp], nel)
            for ibeam = 1:nel
                tline = split(readline(fid))    # Line containing beam ID number
                idbeam = parse(Int, tline[1])
                idbeam == ibeam || error("beam number mismatch")
                for ks = 1:nsta[kp]
                    tline = split(readline(fid))
                    (x1, x2, x3, x4) = (parse(Float64, t) for t in tline)
                    p1kp[ks, ibeam] = complex(x1, x2)
                    p2kp[ks, ibeam] = complex(x3, x4)
                end
            end
            push!(p1, p1kp)
            push!(p2, p2kp)
        end
        return Tfile(title, nel, nsta, npart, p1, p2)
    end
end

