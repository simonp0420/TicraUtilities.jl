using Printf: @printf


"""
    TicraStation 

Contains the data for a single partition of a Ticra station file.
"""
struct TicraStation
    "`npoint`: The number of stations in the partition."
    npoint::Int 

    """
    `u`: A vector of length `npoint` containing the u values (unitless) 
    of each of the optimization stations in the partition.
    """
    u::Vector{Float64}

    """
    `v`: A vector of length `npoint` containing the v values (unitless)
    of each of the optimization stations in the partition.
    """
    v::Vector{Float64}

    """
    `goal`: A vector of length `npoint` containing the optimization goals
    in dB of each of the optimization stations in the partition.
    """
    goal::Vector{Float64}

    """
    `weight`: A vector of length `npoint` containing the optimization weights
    of each of the optimization stations in the partition.
    """
    weight::Vector{Float64}

    """
    `ipol`: A vector of length `npoint` containing the polarization codes 
    for each of the stations in the partition. The codes have the following meanings:
    * 1:  Linear Ludwig 3 x ("copol") component.
    * 2:  Linear Ludwig 3 y ("cxpol") component.
    * 3:  RHCP component.
    * 4:  LHCP component.
    * 5:  Major axis of polarization ellipse.
    * 6:  Minor axis of polarization ellipse.
    * 7:  The ratio co/cx in dB.
    * 8:  The ratio cx/co in dB.
    * 9:  The ratio RHCP/LHCP in dB.
    * 10:  The ratio LHCP/RHCP in dB.
    * 11:  The ratio major/minor in dB.
    * 12:  The ratio minor/major in dB.
    """
    ipol::Vector{Int}

    """
    `rot`: A vector of length `npoint` containing the rotation in degrees
    of linear polarization reference for each of the optimization stations.
    """
    rot::Vector{Float64}

    """
    `att`: A vector of length `npoint` containing the attenuation in dB ≥ 0 
    relative to the sub-satellite point.
    """
    att::Vector{Float64}

    """
    `ID`: A vector of length `npoint` containing the ID string
     for each of the optimization stations in the partition. May be empty.
    """
    ID::Vector{String}

    idx_copol_eoc::Vector{Bool}  # Used by SHSPO only
    idx_sld::Vector{Bool}  # Used by SHSPO only
    idx_XPD_eoc::Vector{Bool}  # Used by SHSPO only
    area_weight::Vector{Float64} # Used by SHSPO only
end

TicraStation() = TicraStation(0,Float64[],Float64[],Float64[],Float64[],Int[],Float64[],Float64[],
    String[],Bool[],Bool[],Bool[],Float64[])

# Define a constructor that doesn't require the `SHSPO`-specific fields:
TicraStation(n::Int, u::AbstractVector{Float64}, v::AbstractVector{Float64}, g::AbstractVector{Float64},
    w::Vector{Float64}, ipol::AbstractVector{Int}, rot::AbstractVector{Float64},
    att::AbstractVector{Float64}, ID::AbstractVector{String}) = 
    TicraStation(n, u, v, g, w, ipol, rot, att, ID, Bool[],Bool[],Bool[],Float64[])


"""
    write_station(stationfile::AbstractString, stdat::TicraStation)
    write_station(stationfile::AbstractString, stdat::AbstractVector{TicraStation})
    
Write a Ticra POS4-compatible optimization station file.  Here, when `stdat` is a vector,
its elements are the partitions in the station file. 
"""
function write_station(stationfile::AbstractString, stadat::AbstractVector{TicraStation})
    open(stationfile, "w") do fid
        for (kpart, d) in enumerate(stadat)
            @printf(fid, "%d\n", d.npoint)
            for k = 1:d.npoint
                @printf(fid, "%9.5f%9.5f%18.9e%18.9e%3d%7.2f%7.2f %s\n",
                        d.u[k], d.v[k], d.goal[k], d.weight[k], d.ipol[k],
                        d.rot[k], d.att[k], d.ID[k])
            end
        end
        @printf(fid, "    0  GRASP8/POS4 station file terminator\n")
    end
    return
end

write_station(stafile::AbstractString, stadat::TicraStation) = write_station(stafile, TicraStation[stadat])


"""
    read_station(stationfile) -> Vector{TicraStation}
    
Read the contents of a Ticra optimization station file. The returned value
a vector of length `npart`, where `npart` is the number of partitions in the file.
"""
function read_station(stafile)
    stadat = open(stafile, "r") do fid
        kpart = 0;  # Partition counter
        stadat = TicraStation[]
        
        # Read number of points in next partition.
        tline = split(readline(fid))
        isempty(tline) && error("Unable to read number of points from first line!")
        npoint = parse(Int, first(tline))
        
        while npoint > 0
            kpart += 1  # Bump number of partitions
            # Preallocate a structure:
            d = TicraStation(npoint,
                zeros(npoint),
                zeros(npoint),
                zeros(npoint),
                zeros(npoint),
                zeros(Int,npoint),
                zeros(npoint),
                zeros(npoint),
                fill("", npoint))
            # Scan through the lines of the file:
            for k = 1:npoint   # Read each line of the partition
                tline = split(readline(fid))
                isempty(tline) && error("End of file encountered!")
                (d.u[k], d.v[k], d.goal[k], d.weight[k]) = (parse(Float64, tline[i]) for i in 1:4)
                d.ipol[k] = parse(Int, tline[5])
                (d.rot[k], d.att[k]) = (parse(Float64, tline[i]) for i in 6:7)
                length(tline) > 7 && (d.ID[k] = tline[8])
            end
            push!(stadat,d)
            
            # Read number of points in next partition.
            tline = split(readline(fid))
            if isempty(tline)
                npoint = 0
            else
                npoint = parse(Int, first(tline))
            end
        end
        return stadat
    end  # do
    return stadat
end
