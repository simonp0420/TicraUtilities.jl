using StaticArrays
using Dierckx: Spline1D
using DSP: unwrap

abstract type Geom end

function _check_legal_length_units(unit::AbstractString)
    u = lowercase(unit)
    @assert u == "mm" || u == "cm" || u == "m" || u == "km" || u == "in"  "Illegal length unit"
end

abstract type Rim <: Geom end

struct EllipticalRim <: Rim
    centre::      SVector{2, Float64}
    half_axis::   SVector{2, Float64}
    rotation::    Float64
    length_unit:: String

    function EllipticalRim(c,h,r,l)
        _check_legal_length_units(l)
        new(c,h,r,l)
    end
end
EllipticalRim(c,h,r) = EllipticalRim(c,h,r,"in")
EllipticalRim(c,h) = EllipticalRim(c,h,0.0,"in")

function EllipticalRim(obj::TorObj)
    @assert obj.objtype == "elliptical_rim" "Incorrect TOR object type: $(obj.objtype)"

    i = findall(==("half_axis"), obj.propname)
    @assert length(i) > 0 "elliptical_rim missing half-axis property"
    @assert length(i) == 1 "elliptical_rim contains more than one half-axis property"
    (xha, yha, xhaunit, yhaunit) = _parse_tor_xy_struct(obj.propval[i[1]])
    
    i = findall(==("centre"), obj.propname)
    if isempty(i)
        xc = yc = 0.0
        xcunit = ycunit = xhaunit
    else
        (xc,yc,xcunit,ycunit) = _parse_tor_xy_struct(obj.propval[i[1]])
        @assert xcunit == ycunit == xhaunit == yhaunit "Inconsistent units"
    end
    
    i = find(==("rotation"), obj.propname)
    if isempty(i)
        rotation = 0.0
    else
        rotation = parse(Float64, obj.propval[i[1]])
    end
    return EllipticalRim([xc,yc], [xha, yha], rotation, xcunit)
end


function Base.show(io::IO, r::EllipticalRim)
    println(io, "EllipticalRim")
    println(io, "  centre: \t$(r.centre)")
    println(io, "  half_axis: \t$(r.half_axis)")
    println(io, "  rotation: \t$(r.rotation)")
    println(io, "  length_unit: \t$(r.length_unit)")
    return nothing
end




struct SuperEllipticalRim <: Rim
    centre::      SVector{2, Float64}
    half_axis::   SVector{2, Float64}
    exponent::    Float64
    rotation::    Float64
    length_unit:: String

    function SuperEllipticalRim(c,h,e,r,l)
        @assert e > 2  "exponent must be greater than 2"
        _check_legal_length_units(l)
        new(c,h,e,r,l)
    end
end
SuperEllipticalRim(c,h,e,r) = SuperEllipticalRim(c,h,e,r,"in")
SuperEllipticalRim(c,h,e) = SuperEllipticalRim(c,h,e,0.0,"in")

function SuperEllipticalRim(obj::TorObj)
    @assert obj.objtype == "superelliptical_rim" "Incorrect TOR object type: $(obj.objtype)"

    i = findall(==("exponent"), obj.propname)
    @assert length(i) == 1 "superelliptical_rim requires exactly one exponent property"
    exponent = parse(Float64, obj.propval[i[1]])

    i = findall(==("half_axis"), obj.propname)
    @assert length(i) == 1 "superelliptical_rim requires exactly one half_axis property"
    (xha,yha,xhaunit,yhaunit) = parse_tor_xy_struct(obj.propval[i[1]])

    i = findall(==("centre"), obj.propname)
    if isempty(i)
        xc = yc = 0.0
        xcunit = ycunit = xhaunit
    else
        (xc,yc,xcunit,ycunit) = parse_tor_xy_struct(obj.propval[i[1]])
    end
    
    i = find(==("rotation"), obj.propname)
    if isempty(i)
        rotation = 0.0
    else
        rotation = parse(Float64, obj.propval[i[1]])
    end

    @assert xcunit == ycunit == xhaunit == yhaunit "Inconsistent units"
    return   SuperEllipticalRim([xc,yc], [xha,yha], exponent, rotation, xcunit)
end

import Base.show
function show(io::IO, r::SuperEllipticalRim)
    println(io, "SuperEllipticalRim")
    println(io, "  centre: \t$(r.centre)")
    println(io, "  half_axis: \t$(r.half_axis)")
    println(io, "  exponent: \t$(r.exponent)")
    println(io, "  rotation: \t$(r.rotation)")
    println(io, "  length_unit: \t$(r.length_unit)")
    return nothing
end


struct TabulatedRimXYold <: Rim
    x::                 Vector{Float64}
    y::                 Vector{Float64}
    centre::            SVector{2, Float64}
    length_unit::       String
    interpobject::      Spline1D

    function TabulatedRimXYold(x::Vector{Float64},y::Vector{Float64},centre,unit,rotation=0.0)
        _check_legal_length_units(unit)
        @assert length(x) == length(y)  "x and y have different lengths"
        # Eliminate duplicate points:
        xyp = unique([(x - centre[1]) (y - centre[2])], dims=1)
        xp = xyp[:,1]
        yp = xyp[:,2]
        # Compute rho and phi:
        rho_in = sqrt.(xp.^2 + yp.^2)
        phi_in = atan2.(yp, xp) .+ rotation
        ind = sortperm(phi_in)  # Sort into order of increasing phi
        rho_in .= rho_in[ind]
        phi_in .= phi_in[ind]
        phi_diff = diff(unwrap(phi_in))
        @assert all(>(0), phi_diff) "Duplicate phi values encountered in rim definition"
        # Ensure that supplied grid points cover [-pi,pi]:
        orig_start = 1
        if phi_in[1] > -π
            pushfirst!(phi_in, phi_in[end] - 2π)
            pushfirst!(rho_in, rho_in[end])
            orig_start = 2
        end
        if phi_in[end] < π
            push!(phi_in, phi_in[orig_start] + 2π)
            push!(rho_in, rho_in[orig_start])
        end
        rho_versus_phi_interp = Spline1D(phi_in, rho_in)
        new(x,y,centre,unit,rho_versus_phi_interp)
    end
end

TabulatedRimXYold(x,y,centre) = TabulatedRimXYold(x,y,centre,"in")

function TabulatedRimXYold(obj::TorObj)
    @assert (obj.objtype == "tabulated_rim_xy" || 
             obj.objtype == "tabulated_rim")  "Incorrect TOR object type: $(obj.objtype)"

    i = findall(==("file_form"), obj.propname)
    #isempty(i) && error("Only file_form = 'old_grasp' is allowed")
    if isempty(i)
        file_form = "old_grasp"
    else
        file_form = obj.propval[i[1]]
    end
    file_form == "old_grasp" || error("Only file_form = 'old_grasp' is allowed")
         
    # Disallowed properties:
    for prop in ["file_xy_number", "file_xy_values", 
                 "file_corner_points_id"]
        i = findall(==(prop), obj.propname)
        @assert isempty(i) "Unable to process property $prop"
    end

    i = findall(==("file_name"), obj.propname)
    file_name = obj.propval[i[1]]

    unit = "m"
    i = findall(==("unit"), obj.propname)
    isempty(i) || (unit = obj.propval[i[1]])
    _check_legal_length_units(unit)
    
    npoints = 0
    i = findall(p -> p == "number_of_points" || p == "n_points", obj.propname)
    isempty(i) || (npoints = parse(Int, obj.propval[i[1]]))

    factor = 1.0
    i = find(p -> p == "scaling_factor" || p == "factor", obj.propname)
    isempty(i) || (factor = parse(Float64, obj.propval[i[1]]))
    @assert abs(factor - 1.0) < 1e-6 "Unable to process non-unity factor"

    centre = [Inf,Inf] # Inf means to determine center automatically
    i = find(p -> p == "centre" || p == "polar_origin", obj.propname)
    if !isempty(i)
        (x, y, xunit, yunit) = _parse_tor_xy_struct(obj.propval[i[1]])
        if xunit == "automatic"
            x = Inf # Special value means to automatically determine the centroid
            xunit = yunit
        else 
            @assert xunit == yunit  "Inconsistent units"
        end
        centre = [x,y]
    end

    translation = [0.0, 0.0] 
    i = find(==("translation"), obj.propname)
    if !isempty(i)
        (x,y,xunit,yunit) = _parse_tor_xy_struct(obj.propval[i[1]])
        @assert xunit == yunit  "Inconsistent units"
        translation = [x,y]    
        @assert all(translation .== 0.0) "Nonzero translation not implemented"
    end

    rotation = 0.0
    i = findall(==("rotation"), obj.propname)
    isempty(i) || (rotation = parse(Float64, obj.propval[i[1]]))
    
    return read_tabulatedrimxyold(file_name, unit, npoints, centre, rotation)
end


Base.show(r::TabulatedRimXYold) = Base.show(STDOUT, r::TabulatedRimXYold)
function Base.show(io::IO, r::TabulatedRimXYold)
    println(io, "TabulatedRimXYold")
    println(io, "  x: \t$(summary(r.x))")
    println(io, "  y: \t$(summary(r.y))")
    println(io, "  centre: \t$(r.centre)")
    println(io, "  length_unit: \t$(r.length_unit)")
    return nothing
end


rimexponent(r::EllipticalRim) = 2.0
rimexponent(r::SuperEllipticalRim) = r.exponent

function insiderim(x::T, y::T, rim::R) where {T<:Real, R<:Union{EllipticalRim,SuperEllipticalRim}}
    (sphi, cphi) = sincosd(rim.rotation)
    (a, b) = rim.half_axis
    (x0, y0) = rim.centre
    exponent = rimexponent(rim)
    xp = (x - x0) * cphi + (y - y0) * sphi 
    yp = (yi - y0) * cphi - (xi - x0) * sphi
    inside = abs(xp/a)^exponent + abs(yp/b)^exponent ≤ 1.0
    return inside
end

function insiderim(x::T, y::T, rim::TabulatedRimXYold) where {T <: Real}
    (x0, y0) = rim.centre
    xp = x - x0
    yp = y - y0
    rhosq = xp * xp + yp * yp
    phi = atan2(yp,xp)
    rhoi = rim.interpobject[phi]
    inside = rhosq ≤ rhoi*rhoi
    return inside
end

"""
    _read_n_numbers(fid::IO, n::Integer, datatype::DataType<:Real, delim_chars)

Read `np` `Real`s from `fid`, where the number of numbers per line is not known.
Return a vector of type `datatype`.  
"""
function _read_n_numbers(io::IO, n::Integer, ::Type{T}, delim_chars) where {T <: Real}
    numbers = type[]
    while length(numbers) < n
        for word in split(readline(io), delim_chars; keepempty=false)
            push!(numbers, parse(T, word))
            length(numbers) == n && break
        end
    end
    return numbers
end

"""
read_tabulatedrimxyold(fname::AbstractString, unit::AbstractString, npoints::Int, 
                                centre::SVector{2,Float64}, rotation::Float64)

Read the data from a `grasp_old`-formatted Ticra rim file and return a `TabulatedRimXYold`
variable.
"""
function read_tabulatedrimxyold(fname::AbstractString, unit::AbstractString, npoints::Int, 
                                centre::SVector{2,Float64}, rotation::Float64)
    delim_chars = ('\t',' ',',','\n','\r','\f','\v')
    open(fname,"r") do fid
        line = readline(fid)  # Comment line
        line = strip(readline(fid))
        linesplit = split(line)
        np = parse(Int, linesplit[1])
        npoints > 0 && (np = min(np, npoints))
        # kspace == 0 -> equal spacing in phi, == 1 -> unequal spacing in phi
        kspace = length(linesplit) > 1 ? parse(Int,linesplit[2]) : 1
        # kxyrp == 1 -> points are cartesian coords, == 2 -> points are cylindrical coords
        kxyrp = length(linesplit) > 2 ?  parse(Int, linesplit[3]) : 1
        @assert np > 2  "np must be greater than 2"
        @assert 0 ≤ kspace ≤ 1 "kspace must be 0 or 1"
        @assert 1 ≤ kxyrp ≤ 2 "kxyrph must be 1 or 2"
        @assert (kspace,kxyrp) != (0,1) "(kspace,kxyrph) == (0,1) not allowed"
        if kspace == 0  
            # Points are equispaced in ϕ:
            rho = _read_n_numbers(fid, np, Float64, delim_chars)
            phi = range(0.0, 2pi * (np-1)/np, length=np)
            x = [ρ * cos(ϕ) for (ρ,ϕ) in zip(rho, phi)]
            y = [ρ * sin(ϕ) for (ρ,ϕ) in zip(rho, phi)]
        else
            if kxyrp == 1  # cartesian coordinates
                xy = reshape(_read_n_numbers(fid, 2*np, Float64, delim_chars), 2, np)
                x = @view xy[1,:]
                y = @view xy[2,:]
            else
                # cylindrical coords
                phidegrho = reshape(_read_n_numbers(fid, 2*np, Float64, delim_chars), 2, np)
                x = zeros(np)
                y = zeros(np)
                for (i, (ϕdeg, ρ)) in enumerate(eachcol(phidegrho))
                    x[i] = ρ * sind(ϕdeg)
                    y[i] = ρ * cosd(ϕdeg)
                end
            end
        end
        xcen = mean(x)
        ycen = mean(y)
        centre = [xcen, ycen]
        return TabulatedRimXYold(x,y,centre,unit,rotation)
    end
end


abstract type CoordinateSystem <: Geom end

struct CoorSys <: CoordinateSystem
    name::         String
    length_unit::  String
    origin::       SVector{3, Float64}
    x_axis::       SVector{3, Float64}
    y_axis::       SVector{3, Float64}
    base::Union{CoorSys, Nothing}
end

CoorSys(n,l,o,x,y) = CoorSys(n,l,o,x,y,nothing)

