# Read and write TEP files

using StaticArrays: SMatrix
using Unitful: Unitful

abstract type TEP end

"""
    TEPscatter <: TEP

Struct containing the data for a single frequency from a TICRA-compatible 
Tabulated Electrical Properties (TEP) file for a scattering surface (i.e., the old-style,
original definition of TEP file in use since GRASP8.)

## Fields
* `title::String`: Contains the title string.
* `theta`: An `AbstractRange` containing the θ values in degrees.
* `phi`: An `AbstractRange` containing the ϕ values in degrees.
* `sff::Array{ComplexF64, 4}`: Contains reflection coefficients for front surface incidence. 
* `sfr::Array{ComplexF64, 4}`: Contains transmission coefficients for rear surface incidence. 
* `srf::Array{ComplexF64, 4}`: Contains transmission coefficients for front surface incidence. 
* `srr::Array{ComplexF64, 4}`: Contains reflection coefficients for rear surface incidence. 

For `s` representing any of the fields `sff`, `sfr`, `srf`, and `srr`, 
`size(s) = (2, 2, length(theta), length(phi))`, and the 2×2 matrix `s[:,:,i,j]` 
is arranged in the order `[sθθ sϕθ; sθϕ sϕϕ]`.

## See Also
* `TEPperiodic`
"""
struct TEPscatter{R<:AbstractRange} <: TEP
    title::String
    theta::R
    phi::R
    sff::Array{ComplexF64, 4}
    sfr::Array{ComplexF64, 4}
    srf::Array{ComplexF64, 4}
    srr::Array{ComplexF64, 4}
end
function TEPscatter(title, theta::AbstractRange, phi::AbstractRange, sff, sfr, srf, srr)
    t, p = promote(theta, phi)
    return TEPscatter(title, t, p, sff, sfr, srf, srr)
end

TEPscatter(; title, theta, phi, sff, sfr, srf, srr) = TEPscatter(title, theta, phi, sff, sfr, srf, srr)

function Base.show(io::IO, mime::MIME"text/plain", t::TEPscatter)
    println(io, "TEPscatter")
    println(io, "  title\t$(t.title)")
    println(io, "  theta\t$(t.theta)")
    println(io, "  phi  \t$(t.phi)")
    println(io, "  sff \t$(summary(t.sff))")
    println(io, "  sfr \t$(summary(t.sfr))")
    println(io, "  srf \t$(summary(t.srf))")
    println(io, "  srr \t$(summary(t.srr))")
    return nothing
end

function Base.show(io::IO, t::TEPscatter)
    print(io, "TEPscatter with theta=$(t.theta), phi=$(t.phi)")
    return nothing
end


"""
    TEPperiodic <: TEP

Struct containing the data from a TICRA-compatible 
Tabulated Electrical Properties (TEP) file for a periodic unit cell.  Note that 
TEP files containing geometrical parameter sweeps are not yet supported.


## Fields
* `version::String`: Must be `"1.0.0"`.
* `puc_object_name::String`: Object name of the periodic unit cell.
* `puc_class_name::String`: Class name of the periodic unit cell.
* `ns`: The number of geometrical parameters which are varied.
* `theta`: An `AbstractRange` containing the θ values in units specified by field `atunit`.
* `phi`: An `AbstractRange` containing the ϕ values in units specified by field `atunit`..
* `aunit::String`: The angular units read from the TEP file.  Typically equal to `"[deg]"`.
* `freqs`: A vector of frequencies, each element of which is a `Unitful` quantity. The elements 
  may not all share the same frequency units.
* `sff::Array{ComplexF64, 5}`: Contains reflection coefficients for front surface incidence. 
* `sfr::Array{ComplexF64, 5}`: Contains transmission coefficients for rear surface incidence. 
* `srf::Array{ComplexF64, 5}`: Contains transmission coefficients for front surface incidence. 
* `srr::Array{ComplexF64, 5}`: Contains reflection coefficients for rear surface incidence. 

For `s` representing any of the fields `sff`, `sfr`, `srf`, and `srr`, 
`size(s) = (2, 2, length(theta), length(phi), length(freqs))`, and the 2×2 matrix `s[:,:,i,j,k]` 
is arranged in the order `[sTETE sTMTE; sTETM sTMTM]`.

## See Also
* `TEPscatter`
"""
struct TEPperiodic{R<:AbstractRange} <: TEP
    title::Vector{String}
    theta::R
    phi::R
    freqs::Unitful.Quantity{Float64, Unitful.𝐓^-1}
    sff::Array{ComplexF64, 5}
    sfr::Array{ComplexF64, 5}
    srf::Array{ComplexF64, 5}
    srr::Array{ComplexF64, 5}
end
function TEPperiodic(title, theta::AbstractRange, phi::AbstractRange, freqs, sff, sfr, srf, srr)
    t, p = promote(theta, phi)
    return TEPperiodic(title, t, p, freqs, sff, sfr, srf, srr)
end

TEPperiodic(; title, theta, phi, freqs, sff, sfr, srf, srr) = TEPperiodic(title, theta, phi, freqs, sff, sfr, srf, srr)

"""
    get_title(tep::TEP)

Return the title string for the `TEP` object.
"""
get_title(tep::TEP) = tep.title


"""
    get_phi(tep::TEP)

Return the phi values (degrees) stored in the `TEP` object. The returned value will be 
some type of [`AbstractRange`](https://docs.julialang.org/en/v1/base/math/#Base.range) object.
"""
get_phi(tep::TEP) = tep.phi

"""
    get_theta(tep::TEP)

Return the theta values (degrees) stored in the `TEP` object. The returned value will be 
some type of [`AbstractRange`](https://docs.julialang.org/en/v1/base/math/#Base.range) object.
"""
get_theta(tep::TEP) = tep.theta

"""
    get_sff(tep::TEP)

Return the reflection array for front incidence. The array has size `(2, 2, length(theta), length(phi)`
for a `TEPscatter` object and `(2, 2, length(theta), length(phi), length(freqs))` for a `TEPperiodic` object.
"""
get_sff(tep::TEP) = tep.sff

"""
    get_sfr(tep::TEP)

Return the transmission array for rear incidence. The array has size `(2, 2, length(theta), length(phi)`
for a `TEPscatter` object and `(2, 2, length(theta), length(phi), length(freqs))` for a `TEPperiodic` object.
"""
get_sfr(tep::TEP) = tep.sfr

"""
    get_srf(tep::TEP)

Return the transmission array for front incidence. The array has size `(2, 2, length(theta), length(phi)`
for a `TEPscatter` object and `(2, 2, length(theta), length(phi), length(freqs))` for a `TEPperiodic` object.
"""
get_srf(tep::TEP) = tep.srf

"""
    get_srr(tep::TEP)

Return the reflection array for rear incidence. The array has size `(2, 2, length(theta), length(phi)`
for a `TEPscatter` object and `(2, 2, length(theta), length(phi), length(freqs))` for a `TEPperiodic` object.
"""
get_srr(tep::TEP) = tep.srr


"""
    get_freqs(tep::TEPperiodic)

Return the vector of frequencies.  The elements of the vector will be `Unitful` quantities, each 
of which may have different frequency units.
"""
get_freqs(tep::TEPperiodic) = tep.freqs


"""
    read_tepfile(filename::AbstractString)

Read a TICRA-compatible "Tabulated Electrical Properties" (TEP) file.  The file
may be in either the original format (scattering surface) or the newer format for
periodic unit cells created by the QUPES program.  The return value depends on the 
type of TEP file encountered:
* Old-style scattering surface: For a TEP file containing data for a single frequency,
  the return value will be a scalar object of type `TEPscatter`.  If the file contains
  data for multiple frequencies, the return value will be an object of type `Vector{TEPscatter}`,
  with one element for each frequency.
* New-style periodic unit cell: The return value will be an object of type `TEPperiodic`.  Note that
  the file may not contain swept geometrical parameters.
"""
function read_tepfile(filename::AbstractString)
    line = readline(filename) # Obtain title line
    if line == "TICRA-EL_PROP-V1.0"
        return _read_tepfile_scatter(filename)
    elseif line == "[TITLE] ELECTRICAL PROPERTIES OF PERIODIC UNIT CELL"
        return _read_tepfile_periodic(filename)
    else
        error("Incorrect title line for TEP file: \"$line\"")
    end
end


function _read_tepfile_scatter(filename::AbstractString)
    teps = open(filename, "r") do fid
        filev = readline(fid)
        filev == "TICRA-EL_PROP-V1.0" || error("Incorrect title line: \"$filev\"")
        firstfreq = true
        while !eof(fid)
            title = readline(fid)

            strs = split(readline(fid), x -> isspace(x) || x == ','; keepempty=false)
            nth = parse(Int, strs[1])
            nth ≥ 4 || error("nth must be ≥ 4. Value is $nth")
            mp = parse(Int, strs[2])
            mp ≥ 1 || error("mp must exceed 0.  Value is $mp")
            phi = range(0, (mp-1)*360/mp, length = mp)
            thetamax = parse(Float64, strs[3])
            thetamax > 0 || error("thetamax must exceed 0.  Value is $thetamax")
            theta = range(0, thetamax, length = nth)
            sff, sfr, srf, srr = (zeros(ComplexF64, (2, 2, nth, mp)) for _ in 1:4)
            for j in 1:mp, i in 1:nth
                sff[:,:,i,j] .= SMatrix{2, 2, ComplexF64}(_read2lines(fid))
                srf[:,:,i,j] .= SMatrix{2, 2, ComplexF64}(_read2lines(fid))
                srr[:,:,i,j] .= SMatrix{2, 2, ComplexF64}(_read2lines(fid))
                sfr[:,:,i,j] .= SMatrix{2, 2, ComplexF64}(_read2lines(fid))
            end
            tep = TEPscatter(; title, theta, phi, sff, sfr, srf, srr)
            if firstfreq 
                teps = [tep]
                firstfreq = false
            else
                push!(teps, tep)
            end
        end # while
        return teps
    end # do closure
    if length(teps) > 1
        return teps
    else
        return only(teps)
    end
end # function

function _read2lines(fid)
    sθθ, sθϕ = _read1line(fid)
    sϕθ, sϕϕ = _read1line(fid)
    return (sθθ, sθϕ, sϕθ, sϕϕ)
end

function _read1line(fid)
    strs = split(readline(fid), x -> isspace(x) || x == ','; keepempty=false)
    length(strs) ≥ 4 || error("Not enough entries on a line of the file")
    c1 = complex(parse(Float64, strs[1]), parse(Float64, strs[2]))
    c2 = complex(parse(Float64, strs[3]), parse(Float64, strs[4]))
    return (c1, c2)
end

