using Dates: Dates

"""
    Ffd
 
Contains data for an HFSS "Far Field Incident Wave Source" object, typically obtained
from an ASCII file with extension "ffd". Note that a single `Ffd` instance 
contains the field data for a single frequency.

### Fields

* `frequency::Float64`: The frequency of the cut in Hz.  Zero if undefined.
* `theta::AbstractRange`: The theta values (in degrees) stored in the cut.
* `phi::AbstractRange`: The phi values (in degrees) stored in the cut.
* `evec`: Matrix of complex field vectors (E_θ, E_ϕ). Each column corresponds to a
  particular phi value. 
"""
@kwdef mutable struct Ffd{Tt<:AbstractRange, Tp<:AbstractRange}
    frequency::Float64 = 0.0
    theta::Tt
    phi::Tp
    evec::Matrix{SVector{2,ComplexF64}}
end

function Base.show(io::IO, mime::MIME"text/plain", t::Ffd)
    println(io, "Ffd")
    println(io, "  frequency  $(t.frequency) [Hz]")
    println(io, "  theta      $(t.theta) [°]")
    println(io, "  phi        $(t.phi) [°]")
    evec_summary = replace(summary(t.evec), "StaticArraysCore." => "")
    println(io, "  evec       " * evec_summary * " [V]")
    return nothing
end

function Base.show(io::IO, t::Ffd)
    print(io, "Ffd with frequency=$(t.frequency), theta=$(t.theta), phi=$(t.phi)")
    return nothing
end

"""
    get_theta(ffd::Ffd)

Return the theta values (in degrees) stored in the `FFd` object.  The returned value will be 
some type of [`AbstractRange`](https://docs.julialang.org/en/v1/base/math/#Base.range) object.
"""
get_theta(ffd::Ffd) = ffd.theta

"""
    get_phi(ffd::Ffd)

Return the phi values (degrees) stored in the `Ffd` object. The returned value will be 
some type of [`AbstractRange`](https://docs.julialang.org/en/v1/base/math/#Base.range) object.
"""
get_phi(ffd::Ffd) = ffd.phi

"""
    get_evec(ffd::Ffd)

Return the ntheta × nphi matrix of complex field vectors stored in the `FFD` object. Each element
of the returned matrix will be a 2-vector.
"""
get_evec(ffd::Ffd) = ffd.evec


"""
    split_ffd_line(line::AbstractString)::Vector{String}
Split an input line from a "*.ffd" file into tokens.
"""
function split_ffd_line(line::AbstractString)
    split(line, x -> isspace(x) || x == ',' || x == ';'; keepempty = false)
end

"""
    read_ffdfile(fname) -> Vector{Ffd}  

Read data from an HFSS-compatible ffd (Far Field Incident Wave Source) file.  

If `fname` contains data for only a single frequency, then return an `Ffd` struct.  If
the file is frequency-independent, then the `frequency` field of the object will be set
to 0.0.

If `fname` contains multiple frequencies, then return a vector of `Ffd` structs,
one for each frequency in the file.
"""
function read_ffdfile(fname::AbstractString)
    ffds = open(fname, "r") do fid
        ffds = Ffd[]
        lsplit = split_ffd_line(readline(fid))[1:3]
        tstart, tstop, nt = (parse(Float64, ls) for ls in lsplit)
        thetas = range(tstart, tstop, round(Int, nt))
        lsplit = split_ffd_line(readline(fid))[1:3]
        pstart, pstop, np = (parse(Float64, ls) for ls in lsplit)
        phis = range(pstart, pstop, round(Int, np))

        # 3rd line determines whether frequency-dependent or independent
        mark(fid) # So we can "unread" the line later
        lsplit = split_ffd_line(readline(fid))
        if lowercase(lsplit[1]) == "frequencies"
            fdependent = true
            nfreqs = parse(Int, lsplit[2])
        else
            fdependent = false
            nfreqs = 1
            reset(fid) # "Unread" the line
        end

        for nf in 1:nfreqs
            if fdependent
                lsplit = split_ffd_line(readline(fid))
                lowercase(lsplit[1]) == "frequency" || error("bad frequency line")
                frequency = parse(Float64, lsplit[2])
            else
                frequency = 0.0
            end

            evec = [SVector(0.0im, 0.0im) for theta in thetas, phi in phis]
            for nt in 1:length(thetas), np in 1:length(phis)
                lsplit = split_ffd_line(readline(fid))[1:4]
                etr, eti, epr, epi = (parse(Float64, es) for es in lsplit)
                evec[nt, np] = SVector(complex(etr, eti), complex(epr, epi))
            end

            ffd = Ffd(frequency=frequency, theta=thetas, phi=phis, evec=evec)
            push!(ffds, ffd)
        end
        return ffds
    end # do closure
    return isone(length(ffds)) ? only(ffds) : ffds
end # function

"""
    write_ffdfile(fname::AbstractString, ffd::Ffd")

    write_ffdfile(fname::AbstractString, ffds::AbstractVector{Ffd})

Write `Ffd` data to an HFSS-compatible "Far Field Incident Wave Source" file.
"""
write_ffdfile(fname::AbstractString, ffd::Ffd) = write_ffdfile(fname, [ffd])

function write_ffdfile(fname::AbstractString, ffds::AbstractVector{<:Ffd})
    open(fname, "w") do fid
        ts = ffds[1].theta
        println(fid, first(ts), " ", last(ts), " ", length(ts))
        ps = ffds[1].phi
        println(fid, first(ps), " ", last(ps), " ", length(ps))
        println(fid, "Frequencies ", length(ffds))
        foreach(ffd -> write_ffd_1freq(fid, ffd), ffds)
    end
    return
end

function write_ffd_1freq(fid, ffd)
    !iszero(ffd.frequency) && println(fid, "Frequency ", ffd.frequency)
    for t in 1:length(ffd.theta), p in 1:length(ffd.phi)
        e = ffd.evec[t, p]
        r1, i1 = reim(e[1])
        r2, i2 = reim(e[2])
        @printf(fid, "%21.15e %21.15e %21.15e %21.15e\n", r1, i1, r2, i2)
    end
    return
end

"""
    ffd2cut(ffd::Ffd; kwargs...) -> cut::Cut

    ffd2cut(ffds::AbstractVector{Ffd}; kwargs...) -> cuts::Vector{Cut}

    ffd2cut(ffdfile::AbstractString; kwargs...) -> cut::Cut

    ffd2cut(ffdfile::AbstractString, cutfile::AbstractString; kwargs...) -> cut::Cut

Convert an HFSS-compatible `Ffd` object to a Ticra-compatible `Cut` object.

The first positional input argument can be either a string containing the name 
of an HFSS-compatible, ASCII text ffd file, or the returned value of type `Ffd` 
(or a vector of `Ffd` objects) that results from reading such a file with `read_ffdfile`.  

The second positional argument, if present, is the name of a Ticra-compatible cut file to 
which the generated `Cut` object(s) should be written.

## Keyword Arguments
- `phi0::Bool = true`: If true, reorder the phi cuts so that they start at phi = 0. In
  this case, also remove any duplicate phi cuts.
"""
function ffd2cut(ffdfile::AbstractString; kwargs...)
    ffd = read_ffdfile(ffdfile)
    cut = ffd2cut(ffd; kwargs...)
    return cut
end

function ffd2cut(ffdfile::AbstractString, cutfile::AbstractString; kwargs...)
    ffdfile == cutfile && throw(ArgumentError("ffdfile cannot be same as cutfile"))
    cut = ffd2cut(ffdfile; kwargs...)
    write_cutfile(cutfile, cut, "Cut created by ffd2cut")
    return cut
end

ffd2cut(ffds::AbstractVector{<:Ffd}; kwargs...) = [ffd2cut(ffd; kwargs...) for ffd in ffds]

function ffd2cut(ffds::AbstractVector{<:Ffd}, cutfile::AbstractString; kwargs...)
    cut = ffd2cut(ffds; kwargs...)
    write_cutfile(cutfile, cut)
    return cut
end

function ffd2cut(ffd::Ffd; phi0 = true)
    θs_ffd = get_theta(ffd)
    ϕs_ffd = get_phi(ffd)
    evec_ffd = get_evec(ffd)
    if phi0 
        dupphi = isequal(mod(ϕs_ffd[begin], 360), mod(ϕs_ffd[end], 360))
        i0 = findfirst(iszero, ϕs_ffd)
        isnothing(i0) && error("ffd does not contain ϕ = 0")
        if isone(i0)
            if dupphi
                evec_cut = evec_ffd[:, 1:end-1]
                ϕs_cut = range(0.0, ϕs_ffd[end-1], length = length(ϕs_ffd) - 1)
            else
                evec_cut = copy(evec_ffd)
                ϕs_cut = copy(ϕs_ffd)
            end
        else
            kstart = dupphi ? 2 : 1
            inds = append!(collect(i0:length(ϕs_ffd)), kstart:i0-1)
            evec_cut = evec_ffd[:, inds]
            ϕs_cut = range(0, 360 + ϕs_ffd[i0-1], length = length(inds))
        end
    else
        evec_cut = copy(evec_ffd)
        ϕs_cut = copy(ϕs_ffd)
    end
    θs_cut = copy(θs_ffd)

    ncomp = 2 # number of polarization components
    icomp = 1 # Eθ/Eϕ decomposition
    icut = 1 # Standard polar cuts
    date, tim = split(string(Dates.now()), 'T')
    f = ffd.frequency
    if iszero(f)
        text = [string("ϕ = ", ϕ, " cut from ffd2cut on ", date, " at ", tim) for ϕ in ϕs_cut]
    else
        text = [string("f = ", f, ", ϕ = ", ϕ, " cut from ffd2cut on ", date, " at ", tim) for ϕ in ϕs_cut]
    end
    cut = Cut(; ncomp, icut, icomp, text, evec = evec_cut, theta = θs_cut, phi = ϕs_cut)
    return cut
    
end


