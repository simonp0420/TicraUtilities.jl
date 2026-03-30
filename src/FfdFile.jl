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

Base.length(::Ffd) = 1

function Base.show(io::IO, mime::MIME"text/plain", t::Ffd)
    println(io, "Ffd")
    println(io, "  frequency  $(t.frequency)")
    println(io, "  theta      $(t.theta)")
    println(io, "  phi        $(t.phi)")
    evec_summary = replace(summary(t.evec), "StaticArraysCore." => "")
    println(io, "  evec       ", evec_summary)
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
    get_frequency(ffd::Ffd)

Return the frequency (in Hz) stored in the `FFd` object.  A zero value indicates a 
frequency-independent `Ffd`.
"""
get_frequency(ffd::Ffd) = ffd.frequency

"""
    isapprox(f1::Ffd, f2::Ffd; kwargs...) -> tf::Bool

Compare two `Ffd` objects for approximate equality.

Compares most fields for perfect equality except `evec`.
The `evec` fields (`Matrix` types) are compared for approximate equality using `isapprox`,
to which any provided keyword arguments are passed.
"""
function Base.isapprox(f1::Ffd, f2::Ffd; kwargs...)
    f1.frequency == f2.frequency || return false
    f1.theta == f2.theta || return false
    f1.phi == f2.phi || return false
    return isapprox(f1.evec, f2.evec; kwargs...)
end

"""
    _split_ffd_line(line::AbstractString)::Vector{String}
Split an input line from a "*.ffd" file into tokens.
"""
function _split_ffd_line(line::AbstractString)
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
        lsplit = _split_ffd_line(readline(fid))[1:3]
        tstart, tstop, nt = (parse(Float64, ls) for ls in lsplit)
        thetas = range(tstart, tstop, round(Int, nt))
        lsplit = _split_ffd_line(readline(fid))[1:3]
        pstart, pstop, np = (parse(Float64, ls) for ls in lsplit)
        phis = range(pstart, pstop, round(Int, np))

        # 3rd line determines whether frequency-dependent or independent
        mark(fid) # So we can "unread" the line later
        lsplit = _split_ffd_line(readline(fid))
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
                lsplit = _split_ffd_line(readline(fid))
                lowercase(lsplit[1]) == "frequency" || error("bad frequency line")
                frequency = parse(Float64, lsplit[2])
            else
                frequency = 0.0
            end

            evec = [SVector(0.0im, 0.0im) for theta in thetas, phi in phis]
            for nt in 1:length(thetas), np in 1:length(phis)
                lsplit = _split_ffd_line(readline(fid))[1:4]
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

        # HFSS .ffd files are frequency-independent when there is a single frequency
        # and that frequency is zero, in which case the "Frequencies" line is omitted.
        freq_dependent = length(ffds) > 1 || !iszero(ffds[1].frequency)
        if freq_dependent
            println(fid, "Frequencies ", length(ffds))
        end

        foreach(ffd -> _write_ffd_1freq(fid, ffd), ffds)
    end
    return
end

function _write_ffd_1freq(fid, ffd)
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

"""
    ffd2sph(ffd::Ffd; keywords...) -> sph::SPHQPartition

    ffd2sph(ffds::AbstractVector{Ffd}; keywords...) -> sphs::Vector{SPHQPartition}

    ffd2sph(ffdfile::AbstractString; kwargs...) -> s::SPHQPartition

    ffd2sph(ffdfile::AbstractString, sphfile::AbstractString; kwargs...) -> s::SPHQPartition

Convert a `Ffd` object to a `SPHQPartition` using recursive FFT/IFFT methods from
the Hansen 1988 book "Spherical Near-Field Antenna Measurements.

The first positional input argument can be either a string containing the name 
of an HFSS-compatible, ASCII text .ffd file, or the returned value of type `Ffd` 
(or a vector of `Ffd` objects) that results from reading such a file with `read_ffdfile`.  
The output of this function can be passed to `write_sphfile` to create a Ticra-compatible 
file of Q-type spherical wave coefficients.  This is done automatically if the second positional
argument is provided as shown above.

If the input data extend in θ only to θ₀ < 180°, then it will be assumed that
the fields are identically zero for θ₀ < θ ≤ 180°.

## Keyword Arguments (and their default values)
* `mmax=1000`: An upper limit for the `m` (azimuthal) mode index to be included.
  The actual limit will be set to `min(Nϕ÷2, mmax)` for odd `Nϕ`, and `min(Nϕ÷2-1, mmax)`
  for even `Nϕ`, where `Nϕ` is the number of ϕ = constant polar cuts in the cut object.
* `nmax=1000`: An upper limit for the `n` (polar) mode index to be included.
  The actual limit will be the lesser of `nmax` and `Nθ-1` where `Nθ` is the number of 
  θ values included in each ϕ = constant polar cut.
* `pwrtol=1e-6`: The power tolerance.  Spherical modes are included until the excluded
  modes' power is less than `pwrtol` times the total modal power.  A zero or negative value
  precludes removal of any modes. A small, nonzero value here prevents the inclusion of many
  insignificant modes when the far-field sphere is over-sampled.
* `style = :hfss`: This argument (a `Symbol`) is only relevant for the 2-argument method of `ffd2sph` (i.e. 
  when a spherical wave file is to be written.) In that case, the default value of `:hfss` ensures
  that an HFSS-compatible spherical wave file containing frequency information is created. When
  `style = :ticra` (the only other legal value), then a standard, Ticra-compatible, Q-type, spherical
   wave expansion file will be created, one that does not contain any frequency information.
"""
function ffd2sph(ffd::Ffd; pwrtol = 1e-6, kwargs...)
    cut = ffd2cut(ffd)
    sph = cut2sph(cut; pwrtol, frequency = ffd.frequency, kwargs...)
    return sph
end

ffd2sph(ffds::AbstractVector{Ffd}; keywords...) = [ffd2sph(ffd) for ffd in ffds]

ffd2sph(ffdfile::AbstractString; kwargs...) = ffd2sph(read_ffdfile(ffdfile))

function ffd2sph(ffdfile::AbstractString, sphfile::AbstractString; style::Symbol = :ticra, kwargs...)
    ffdfile == sphfile && error("ffdfile and sphfile must be distinct")
    sph = ffd2sph(read_ffdfile(ffdfile); kwargs...)
    write_sphfile(sphfile, sph; style)
    return sph
end

"""
    cut2ffd(cut::Cut; frequency = 0.0) -> ffd::Ffd

    cut2ffd(cuts::AbstractVector{Cut}; frequency::AbstractVector{<:Real}) -> ffds::Vector{Ffd}

    cut2ffd(cutfile::AbstractString; frequency) -> ffd::Ffd

    cut2ffd(cutfile::AbstractString, ffdfile::AbstractString; frequency) -> ffd::Ffd

Convert a Ticra-compatible `Cut` object to an HFSS-compatible `Ffd` object.

The first positional input argument can be either a string containing the name 
of a Ticra-compatible, spherical polar cut file, or the returned value of type `Cut` 
(or a vector of `Cut` objects) that results from reading such a file with `read_cutfile`.  

The second positional argument, if present, is the name of an HFSS-compatible ffd file to 
which the generated `Ffd` object(s) should be written.

Note: The `Cut` object (or file) need not be in Eθ/Eϕ format (`icomp = 1`).  If needed,
polarization basis conversion to Eθ/Eϕ will be performed automatically.

## Keyword Arguments
- `frequency`: When converting a single `Cut` object, the default value is 0.0, implying that
  a frequency-independent `Ffd` object is desired. If a positive value is provided (in Hz), a
  frequency-dependent `Ffd` is created containing that single frequency.
  When the positional input argument consists of a vector of multiple `Cut` objects (or is 
  a file name containing multiple `Cut`s), then `frequency` must be an abstract vector 
  containing the same number of positive frequencies (in Hz), as needed to create the multiple-frequency
  output `Ffd` vector (or file). If the input file contains only a single frequency, then `frequency` 
  can be a scalar, or can be totally omitted, in which case a frequency-independent `Ffd` will be created.
"""
function cut2ffd end

function cut2ffd(c::Cut; frequency::Float64 = 0.0)
    frequency ≥ 0 || throw(ArgumentError("frequency must be ≥ 0"))
    get_ncomp(c) == 2 || error("Cut must contain only 2 polarization components")
    get_icut(c) == 1 || error("Cut.icut ≠ 1.  Only constant ϕ (polar) cuts are allowed.")
    c = convert_cut(c, 1) # Copy cut and ensure Eθ/Eϕ decomposition
    theta = get_theta(c)
    phi = get_phi(c)
    evec = get_evec(c)
    ffd = Ffd(; frequency, theta, phi, evec)
    return ffd
end

function cut2ffd(cuts::AbstractVector{<:Cut}; frequency::AbstractVector{<:Real})
    nf, nc = length.((frequency, cuts))
    nf == nc || error("# frequencies $nf does not equal # cuts $nc")
    all(>(0), frequency) || throw(ArgumentError("frequency must be all positive numbers"))
    ffds = [cut2ffd(cut; frequency=f) for (cut, f) in zip(cuts, frequency)]
    return ffds
end

function cut2ffd(cutfile::AbstractString; frequency = 0.0)
    cuts = read_cutfile(cutfile)
    ffds = cut2ffd(cuts; frequency)
    return ffds
end

function cut2ffd(cutfile::AbstractString, ffdfile::AbstractString; frequency = 0.0)
    ffds = cut2ffd(cutfile; frequency)
    write_ffdfile(ffdfile, ffds)
    return ffds
end


"""
    sph2ffd(sphfile::AbstractString; kwargs...) -> ffd::Ffd

    sph2ffd(sph:SPHQPartition; kwargs...) -> ffd::Ffd

    sph2ffd(sphs:Vector{SPHQPartition}; kwargs...) -> ffds::Vector{Ffd}

    sph2ffd(sphfile::AbstractString, ffdfile::AbstractString; kwargs...) -> ffd::Ffd

Convert a set of Q-type spherical wave modal coefficients to far-field electric field 
values, returned as a `Ffd` object. 

The first positional input argument can be either a string containing the name 
of a Ticra-compatible (*.sph) or HFSS-compatible (*.swef) Q-type spherical wave 
file, or or the returned value from reading such a file with `read_sphfile`.
The output of this function can be passed to `write_ffdfile` to create an HFSS-compatible 
far-field data file (.ffd file).  This is done automatically if the second positional
argument is provided as shown above.


## Keyword Arguments (all are optional, except possibly `frequency`):
- `theta`: An abstract range denoting the desired polar angles (colattitudes) 
  in degrees at which the field should be evaluated. If an empty range is provided (the default), then
  the values will be determined automatically by examining the input spherical mode content.
- `phi`: An abstract range denoting the desired azimuthal angles in degrees at 
  which the field should be evaluated.  If an empty range is provided (the default), then
  the values will be determined automatically by examining the input spherical mode content.
- `frequency`: The frequency in Hz.  If zero (the default value), then the frequency stored in the
  output will be determined from the positional input argument. However, if the `frequency` 
  argument is positive, then this value will be used in preference to the value determined from
  examining the positional argument.  For the case when `sphs` is a vector of `SPHQPartition` objects
  (or a file containing multiple partitions), then the `frequency` argument should be set to a vector
  of the same length.  Note that a nondefault value is required when the input consists of multiple 
  partitions where each partition has zero in its `frequency` field.  This case arises when the input
  partitions arise from reading Ticra-compatible spherical wave files.

## Usage Example
    ffd = sph2ffd("testfile.sph"; phi=0:5:355, theta=0:1:180)
"""
function sph2ffd(sphfile::AbstractString; kwargs...)
    sph = read_sphfile(sphfile)
    ffd = sph2ffd(sph; kwargs...)
    return ffd
end

function sph2ffd(sphfile::AbstractString, ffdfile::AbstractString; kwargs...)
    sph = read_sphfile(sphfile)
    ffd = sph2ffd(sph; kwargs...)
    write_ffdfile(ffdfile, ffd)
    return ffd
end

function sph2ffd(sphs::AbstractVector{SPHQPartition};
        frequency::AbstractVector{<:Real} = zeros(length(sphs)),
        kwargs...)

    length(sphs) == length(frequency) || throw(ArgumentError("sphs and frequency have different lengths"))
    if length(sphs) > 1
        for (sph, f) in zip(sphs, frequency)
            if iszero(sph.frequency) && iszero(f)
                error("Zero frequency not allowed in multifrequency vector{Ffd}")
            end
        end
    end
    ffds = [sph2ffd(sph; frequency = f, kwargs...) for (sph, f) in zip(sphs, frequency)]
    return ffds
end

function sph2ffd(sph::SPHQPartition;
    theta::AbstractRange=0.0:-1.0:1.0,
    phi::AbstractRange=0.0:-1.0:1.0,
    frequency = 0.0)

    ipol = 1 # Eθ/Eϕ decomposition
    cut = sph2cut(sph; theta, phi, ipol)
    if iszero(frequency)
        frequency = get_frequency(sph)
    end
    ffd = cut2ffd(cut; frequency)
    return ffd
end
