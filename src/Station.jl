using Printf: @printf


"""
    Station 

A struct that contains data for a single partition of a Ticra-compatible optimization station file.


### Fields
* `npoint::Int`: The number of stations in the partition.
* `u::Vector{Float64}`: A vector of length `npoint` containing the ``u`` values (unitless direction cosine along ``x``) 
  of each of the optimization stations in the partition.
* `v::Vector{Float64}`: A vector of length `npoint` containing the ``v`` values (unitless direction cosine along ``y``) 
  of each of the optimization stations in the partition.
* `goal::Vector{Float64}`: A vector of length `npoint` containing the goal values in dB 
  for each of the optimization stations in the partition.
* `weight::Vector{Float64}`: A vector of length `npoint` containing the optimization weights
  for each of the optimization stations in the partition.
* `ipol::Vector{Int}`: A vector of length `npoint` containing the polarization codes 
  for each of the stations in the partition. The codes have the following meanings:
    - 1:  Linear Ludwig 3 x ("copol") component.
    - 2:  Linear Ludwig 3 y ("cxpol") component.
    - 3:  RHCP component.
    - 4:  LHCP component.
    - 5:  Major axis of polarization ellipse.
    - 6:  Minor axis of polarization ellipse.
    - 7:  The ratio co/cx in dB.
    - 8:  The ratio cx/co in dB.
    - 9:  The ratio RHCP/LHCP in dB.
    - 10:  The ratio LHCP/RHCP in dB.
    - 11:  The ratio major/minor in dB.
    - 12:  The ratio minor/major in dB.
    - 13:  The total power 10*log₁₀(‖E⃗‖²) in dB.
* `rot::Vector{Float64}`: A vector of length `npoint` containing the rotation in degrees
  of linear polarization reference for each of the optimization stations.
* `att::Vector{Float64}`: A vector of length `npoint` containing the attenuation in dB ≥ 0 
  relative to the sub-satellite point.
* `ID::Vector{String}`: A vector of length `npoint` containing the ID string
  for each of the optimization stations in the partition. May be empty.
"""
struct Station
    npoint::Int
    u::Vector{Float64}
    v::Vector{Float64}
    goal::Vector{Float64}
    weight::Vector{Float64}
    ipol::Vector{Int}
    rot::Vector{Float64}
    att::Vector{Float64}
    ID::Vector{String}
end

Station() = Station(0, Float64[], Float64[], Float64[], Float64[], Int[], Float64[], Float64[], String[])

"""
    get_npoint(s::Station)

Return `npoint`, the number of stations.
"""
get_npoint(s::Station) = s.npoint

"""
    get_u(s::Station)

Return `u`, the vector of unitless station direction cosines along ``x``.
"""
get_u(s::Station) = s.u

"""
    get_v(s::Station)

Return `v`, the vector of unitless station direction cosines along ``y``.
"""
get_v(s::Station) = s.v

"""
    get_goal(s::Station)

Return `goal`, the vector of optimization goals in dB.
"""
get_goal(s::Station) = s.goal

"""
    get_weight(s::Station)

Return `weight`, the vector of station optimization weights.
"""
get_weight(s::Station) = s.weight

"""
    get_ipol(s::Station)

Return `ipol`, the vector of station polarization codes.
"""
get_ipol(s::Station) = s.ipol

"""
    get_rot(s::Station)

Return `rot`, the vector of station polarization rotation angles in degrees.
"""
get_rot(s::Station) = s.rot

"""
    get_att(s::Station)

Return `att`, the vector of attenuation values wrt nadir in dB.
"""
get_att(s::Station) = s.att

"""
    get_id(s::Station)

Return `id`, the vector of station ID strings.
"""
get_id(s::Station) = s.ID


function Base.show(io::IO, sta::Station)
    print(io, "Station partition with ", length(sta.u), " stations")
end



"""
    write_stationfile(stationfile::AbstractString, stdat::Station)
    write_stationfile(stationfile::AbstractString, stdat::AbstractVector{Station})
    
Write a Ticra POS4-compatible optimization station file.  Here, when `stdat` is a vector,
its elements are the partitions in the station file. 
"""
function write_stationfile(stationfile::AbstractString, stadat::AbstractVector{Station})
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

write_stationfile(stafile::AbstractString, stadat::Station) = write_stationfile(stafile, Station[stadat])


"""
    read_stationfile(stationfile) -> Vector{Station}
    
Read the contents of a Ticra optimization station file. The returned value is
a vector of length `npart`, where `npart` is the number of partitions in the file.
"""
function read_stationfile(stafile)
    stadat = open(stafile, "r") do fid
        kpart = 0  # Partition counter
        stadat = Station[]

        # Read number of points in next partition.
        tline = split(readline(fid))
        isempty(tline) && error("Unable to read number of points from first line!")
        npoint = parse(Int, first(tline))

        while npoint > 0
            kpart += 1  # Bump number of partitions
            # Preallocate a structure:
            d = Station(npoint,
                zeros(npoint),
                zeros(npoint),
                zeros(npoint),
                zeros(npoint),
                zeros(Int, npoint),
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
            push!(stadat, d)

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
