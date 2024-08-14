using Printf: @printf
using DSP: unwrap
using Dierckx: Spline1D, integrate
using DataInterpolations: CubicSpline
using QuadGK: quadgk
using StaticArrays: @SVector, @SMatrix, SVector
using LinearAlgebra: ⋅, norm


"""
`Cut` Holds data for a Ticra "Tabulated Pattern Data" object.
Note that a single `Cut` instance contains all of the cuts for a single frequency.

### Fields

* `ncomp::Int`: Number of polarization components (2 or 3)
* `icut::Int`: 1 for standard constant ϕ polar cuts, 2 for conical, constant θ cuts1
* `icomp::Int`: Polarization control parameter. 1 for Eθ and Eφ, 2 for RHCP and LHCP, 3 for Co and Cx (Ludwig 3).
* `text::Vector{String}`: Identification text for each constant angle cut.
* `theta::T<:AbstractRange`: The theta values (in degrees) stored in the cut.
* `phi::T<:AbstractRange`: The phi values (in degrees) stored in the cut.
* `evec`: Matrix of complex field vectors for the two or three polarization components.
"""
@kwdef mutable struct Cut{T <: AbstractRange, N}
    ncomp::Int
    icut::Int
    icomp::Int
    text::Vector{String}
    theta::T
    phi::T
    evec::Matrix{SVector{N,ComplexF64}}
end


import Base.show
function show(io::IO, t::Cut)
    println(io, "Cut")
    println(io, "  ncomp\t$(t.ncomp)")
    println(io, "  icut \t$(t.icut)")
    println(io, "  icomp\t$(t.icomp)")
    print(io, "  text")
    if isempty(t.text)
        println(io, "")
    elseif length(t.text) < 6
        println(io,"")
        for line in t.text
            println(io, "      \t" * line)
        end
    else
        println(io, ": ", summary(t.text))
    end
    if !isempty(t.theta)
        println(io, "  theta\t$(t.theta)")
    else
        println(io, "  theta\tEmpty range")
    end
    if !isempty(t.phi)
        println(io, "  phi  \t$(t.phi)")
    else
        println(io, "  phi  \tEmpty range")
    end
    println(io, "  evec \t$(summary(t.evec))")
    return nothing
end

"""
    get_theta(c::Cut)

Return the theta values (in degrees) stored in the cut
"""
get_theta(c::Cut) = c.theta

"""
    get_phi(c::Cut)

Return the phi values (degrees) stored in the cut.
"""
get_phi(c::Cut) = c.phi

"""
    get_evec(c::Cut)

Return the ntheta × nphi matrix of complex field vectors stored in the cut.
"""
get_evec(c::Cut) = c.evec

"""
    get_ncomp(c::Cut)

Return ncomp, the number of polarization components stored in the cut.
"""
get_ncomp(c::Cut) = c.ncomp

"""
    get_icut(c::Cut)

Return icut, the control parameter of the cut. 1 for a polar cut, 2 for a conical cut.
"""
get_icut(c::Cut) = c.icut

"""
    get_icomp(c::Cut)

Return icomp, polarization parameter. 1 for Eθ and Eφ, 2 for RHCP and LHCP, 3 for Co and Cx (Ludwig 3).
"""
get_icomp(c::Cut) = c.icomp


"""
    get_text(c::Cut)

Return a vector of strings containing the cut identification text.
"""
get_text(c::Cut) = c.text

"""
    amplitude_db(c::Cut, ipol::Int)
    amplitude_db(c::Cut, polstr::String = "copol")

Return a matrix of amplitudes in dB for some choice of polarization in the cut.
Legal values for `ipol` are 1, 2, or 3.  Legal values for `polstr` are
"copol" (the default) and "crosspol".
"""
function amplitude_db end

amplitude_db(c::Cut, ipol::Integer) = 10 * log10.(abs2.(getindex.(get_evec(c), ipol)))

function amplitude_db(c::Cut, polstr::String = "copol")
    polstr = lowercase(strip(polstr))
    minflag = maxflag = false
    if isequal(polstr, "copol")
        maxflag = true
    elseif isequal(polstr, "crosspol")
        minflag = true
    else
        throw(ArgumentError("Unknown string \"$polstr\""))
    end
    if maxflag || minflag
        p1 = abs2.(getindex.(get_evec(c), 1))
        p1maxsq = maximum(p1)
        p2 = abs2.(getindex.(get_evec(c), 2))
        p2maxsq = maximum(p2)
        if (p1maxsq > p2maxsq && maxflag) || (p1maxsq < p2maxsq && minflag)
            return 10 * log10.(p1)
        else
            return 10 * log10.(p2)
        end
    end
end


"""
    phase_deg(c::Cut, ipol::Int)
    phase_deg(c::Cut, polstr::String = "copol")

Return a matrix of phases in degrees for some choice of polarization in the cut.
Legal values for `ipol` are 1, 2, or 3.  Legal values for `polstr` are
"copol" (the default) and "crosspol".
"""
function phase_deg end

phase_deg(c::Cut, ipol::Int) = rad2deg.(angle.(getindex.(get_evec(c), ipol)))

function phase_deg(c::Cut, polstr::String = "copol")
    polstr = lowercase(polstr)
    minflag = maxflag = false
    if isequal(polstr, "copol")
        maxflag = true
    elseif isequal(polstr, "crosspol")
        minflag = true
    else
        throw(ArgumentError("Unknown string \"$polstr\""))
    end
    if maxflag || minflag
        p1maxsq = maximum(x -> abs2(x[1]), get_evec(c))
        p2maxsq = maximum(x -> abs2(x[2]), get_evec(c))
        if (p1maxsq > p2maxsq && maxflag) || (p1maxsq < p2maxsq && minflag)
            ipol = 1
        else
            ipol = 2
        end
    end
    return phase_deg(c, ipol)
end


"""
    power(cut::Cut)::Float64
    power(cut::Cut, θmax=180)
    
Compute the total radiated power in a Cut object. 
If only a single phi value is included in the cut, then assume no azimuthal variation.
The integration in the θ direction will be computed over the limits from 0 to `min(θmax, last(get_theta(cut)))`.
"""
function power(cut::Cut, θmax=180)::Float64
    # This version uses the trapezoidal rule in phi and integration of a
    # cubic spline interpolant in the theta direction.
    get_ncomp(cut) == 3 && error("Cut has 3 field components.  Only 2 allowed.")
    sym = get_theta(cut)[begin] < 0  # Symmetrical cut
    phifullmax = sym ? 180.0 : 360.0
    phi = get_phi(cut)
    nphi = length(phi)
    p = vec(sum(x -> real(x ⋅ x), get_evec(cut), dims = 2))
    if nphi > 1
        dphi = abs(phi[begin+1] - cut.phi[begin])
        # Check that the full range of phi is covered in the cut object:
        fullrange = abs(abs(phi[end] - phi[begin]) - (phifullmax - dphi)) ≤ 1e-2
        fullrange || error("Full range of phi values not provided in cut")
    else
        dphi = 360.0
    end
    phimult = deg2rad(dphi)  # Factor for trap. int. over phi
    theta = deg2rad.(get_theta(cut))
    p = p .* abs.(sin.(theta)) # sintheta weighting
    spl = Spline1D(theta, p)
    pwr = phimult * integrate(spl, theta[1], min(theta[end], deg2rad(θmax)))
    return pwr
end


"""
    read_cutfile(fname, warnflag=true)::Cut

Read the first frequency's data from a Ticra-compatible cut file.  
`warnflag`, if `true`, causes this function to issue a warning
if more than one frequency is detected in the file.
"""
function read_cutfile(fname::AbstractString, warnflag::Bool = true)
    # Read a single frequency from a (possibly multi-frequency) cut file.
    # warnflag, if true, instructs this routine to issue a warning
    # if more than one frequency is present in the file.
    t = read_cutfiles(fname)
    n = length(t)
    n > 1 && warnflag && @warn "$n frequencies found in $fname.  Returning only 1st..."
    return t[1]
end

"""
    read_cutfiles(fname)::Vector{Cut}  

Read data from a possibly multi-frequency Ticra-compatible cut file.  
Return a vector of one or more `Cut` structs.
"""
function read_cutfiles(fname::AbstractString)
    cuts = open(fname, "r") do fid
        cutphi = Float64[]
        header = String[]
        kf = 0  # Initialize frequency counter
        local dth, header, icomp, icut, kf, ncomp, nphi, nth, evecallphi, evecnext, ths, cuts
        firstcut = true
        while !eof(fid)
            textline = rstrip(readline(fid))
            str = readline(fid)
            ts = split(str)
            length(ts) == 6 && push!(ts, "2")
            (ths, dth, phi) = (parse(Float64, ts[i]) for i in (1, 2, 4))
            (nth, icomp, icut, ncomp) = (parse(Int, ts[i]) for i in (3, 5, 6, 7))
            if kf == 0 || phi ∈ cuts[end].phi # Begin a new frequency
                if kf ≠ 0 # Finish off old frequency
                    cut = cuts[end]
                    evec = reshape(evecallphi, length(cut.theta), length(cut.phi))
                    cuts[end] = Cut(cut.ncomp, cut.icut, cut.icomp, cut.text, cut.theta, cut.phi, evec)
                else
                    if ncomp == 2
                        evecnext = [@SVector[0.0+0.0im, 0.0+0.0im] for _ in 1:nth]
                    else
                        evecnext = [@SVector[0.0+0.0im, 0.0+0.0im, 0.0+0.0im] for _ in 1:nth]# Storage increment for reading field values
                    end
                end
                kf += 1
                header = String[textline]
                cutphi = phi:phi
                theta = range(start = ths, step = dth, length = nth)
                evecallphi = ncomp == 2 ? Array{SVector{2,ComplexF64},1}(undef,0) : Array{SVector{3,ComplexF64},1}(undef,0) 
                evec = ncomp == 2 ? Array{SVector{2,ComplexF64},2}(undef,0,0) : Array{SVector{3,ComplexF64},2}(undef,0,0) 
                cut = Cut(ncomp, icut, icomp, header, theta, cutphi, evec)
                if firstcut
                    cuts = [cut]
                    firstcut = false
                else
                    push!(cuts, cut)
                end
            else
                # Begin an additional phi cut at current frequency:
                cut = cuts[end]
                cutphi = first(cut.phi):(phi-last(cut.phi)):phi
                push!(cut.text, textline)
                cuts[end] = Cut(cut.ncomp, cut.icut, cut.icomp, cut.text, cut.theta, cutphi, cut.evec)
            end
            # Check consistency
            cut = cuts[end]
            @assert (icomp, icut, ncomp) == (cut.icomp, cut.icut, cut.ncomp)
            @assert (ths, dth, nth) == (first(cut.theta), cut.theta[2]-cut.theta[1], length(cut.theta)) 
            append!(evecallphi, evecnext)
            evecthisphi = @view evecallphi[end-nth+1:end]

            # Read in the field data for this phi cut
            for i = 1:nth
                str = readline(fid)
                if ncomp == 2
                    (t1, t2, t3, t4) = (parse(Float64, s) for s in split(str))
                    evecthisphi[i] = @SVector[complex(t1,t2), complex(t3,t4)]
                else
                    (t1, t2, t3, t4, t5, t6) = (parse(Float64, s) for s in split(str))
                    evecthisphi[i] = @SVector[complex(t1,t2), complex(t3,t4), complex(t5,t6)]
                end
            end
        end # while
        # Finish final cut
        cut = cuts[end]
        nphi = length(cut.phi)
        nth = length(cut.theta)
        evec = reshape(evecallphi, nth, nphi)
        cut = cuts[end]
        cuts[end] = Cut(cut.ncomp, cut.icut, cut.icomp, cut.text, cut.theta, cut.phi, evec)

        return cuts
    end # function
    return cuts
end

"""
    write_ticra_cut(fname::AbstractString, cut::Cut, title::AbstractString="Cut file created by write_ticra_cut")

    write_ticra_cut(fname::AbstractString, cuts::AbstractVector{Cut}, title::AbstractString="Cut file created by write_ticra_cut")

Write `Cut` cut data to a Ticra-compatible cut file.
"""
function write_ticra_cut(
    fname::AbstractString,
    cut::Cut{T,N},
    title::String = "Cut file created by write_ticra_cut"
) where {T,N}
    open(fname, "w") do fid
        for (n, phi) in enumerate(cut.phi)
            @printf(fid, "%s, phi = %8.3f\n", title, phi)
            @printf(fid, "%8.3f %8.3f %d %8.3f %d %d %d\n",
                cut.theta[1], cut.theta[2] - cut.theta[1],
                length(cut.theta), phi, cut.icomp, cut.icut,
                cut.ncomp)
            for m in 1:length(cut.theta)
                e = cut.evec[m, n]
                e1 = e[1]
                e2 = e[2]
                if N == 2
                    @printf(fid, " %18.10E %18.10E %18.10E %18.10E\n",
                        real(e1), imag(e1), real(e2), imag(e2))
                elseif N == 3
                    e3 = e[3]
                    @printf(fid, " %18.10E %18.10E %18.10E %18.10E %18.10E %18.10E\n",
                        real(e1), imag(e1), real(e2), imag(e2), real(e3), imag(e3))
                else
                    error("Unknown type parameter $N")
                end
            end
        end
    end
end

function write_ticra_cut(
    fname::AbstractString,
    cuts::AbstractVector{Cut{T,N}},
    title::String = "Cut file created by write_ticra_cut",
    ) where {T,N}
    open(fname, "w") do fid
        for cut in cuts
            for (n, phi) in enumerate(cut.phi)
                @printf(fid, "%s, phi = %8.3f\n", title, phi)
                @printf(fid, "%8.3f %8.3f %d %8.3f %d %d %d\n",
                    cut.theta[1], cut.theta[2] - cut.theta[1], length(cut.theta),
                    phi, cut.icomp, cut.icut, cut.ncomp)
                for m in 1:length(cut.theta)
                    e = cut.evec[m, n]
                    e1 = e[1]
                    e2 = e[2]
                    if N == 2
                        @printf(fid, " %18.10E %18.10E %18.10E %18.10E\n",
                            real(e1), imag(e1), real(e2), imag(e2))
                    elseif N == 3
                        e3 = e[3]
                        @printf(fid, " %18.10E %18.10E %18.10E %18.10E %18.10E %18.10E\n",
                            real(e1), imag(e1), real(e2), imag(e2), real(e3), imag(e3))
                    else
                            error("Unknown type parameter $N")
                    end
                end
            end
        end
    end
end


"""
    (x,y,z0,z90) = phscen(cutfile::AbstractString, fghz=11.802852677165355; min_dropoff=-10)

    (x,y,z0,z90) = phscen(cut::Cut, fghz=11.802852677165355; min_dropoff=-10)
    
Estimate phase center for a cut file or `Cut` object using NSI least-squares 
algorithm.

The three outputs are estimates of the location of the phase center relative to
the origin used in recording the data in the cut.  If `fghz` is passed in as 
an argument, the values will be expressed in units of inches.  Otherwise the 
lengths will be normalized to wavelengths. Note that `z0` and `z90` are the 
phi = 0ᵒ and phi = 90ᵒ plane estimates of the phase center location.  
In determining the phase center locations, only field values with magnitudes
in dB greater than `min_dropoff` relative to the peak field are considered.
"""
function phscen end

function phscen(cut::Cut, fghz = 11.802852677165355; min_dropoff = -10.0)
    # Determine which pol slot has main polarization:
    p1max, p2max = (maximum(x -> abs2(x[i]), get_evec(cut)) for i in 1:2)
    p = p1max > p2max ? 1 : 2
    E = getindex.(get_evec(cut), p)

    thetamin = minimum(get_theta(cut))
    x = deg2rad.(cut.theta)
    if thetamin ≥ 0.0
        x = vcat(-reverse(x), x)  # Also cover negative theta angles
    end

    # Boresight electric field:
    thetamin, index = findmin(abs, get_theta(cut))
    Ebore = E[index, begin]
    # Normalize to boresight (including phase):
    E = E ./ Ebore
    Edb = 10 * log10.(abs2.(E))
    phase = angle.(E)

    if length(cut.phi) == 1
        phi = collect(cut.phi)
    else
        phi = [0.0, 90.0]
    end
    krho = zeros(length(phi))
    kz = zeros(length(phi))
    for (k, ph) in enumerate(phi)
        index = findall(==(ph), get_phi(cut))
        isempty(index) && error("Unable to find cut phi = $(ph)!")
        length(index) > 1 && error("Too many cuts with phi = $(ph)!")
        y = phase[:, index[1]]
        cutdb = Edb[:, index]
        if thetamin ≥ 0
            index = findall(==(180 + ph), get_phi(cut))
            if isempty(index)
                index = findall(==(-180 + ph), get_phi(cut))
            end
            if isempty(index)
                if length(cut.phi) == 1
                    index = [1]
                else
                    error("Unable to find cut opposite phi = $(ph)!")
                end
            end
            y = vcat(reverse(phase[:, index[1]]), y)  # Add phases for (theta, phi+180),
            # which are the same as those for (-theta, phi).
            cutdb = vcat(reverse(Edb[:, index[1]]), cutdb)
        end

        y = unwrap(y)
        (_, kk) = findmin(abs, x)
        y .-= y[kk] # Normalize to theta = 0 value
        big_enough = findall(≥(min_dropoff), cutdb)
        (krho[k], kz[k]) = _find_phase_center(x[big_enough], y[big_enough])

    end

    wl = 11.802852677165355 / fghz  # wavelength in inches
    wn = 2π / wl

    x0 = krho[1] / wn
    if length(cut.phi) > 1
        y0 = krho[2] / wn
    else
        y0 = x0
    end
    z0 = kz ./ wn
    return (x0, y0, z0[1], z0[2])
end

function phscen(cutfile::AbstractString, fghz = 11.802852677165355; min_dropoff = -10)
    cut = read_cutfile(cutfile)
    phscen(cut, fghz; min_dropoff)
end


function _find_phase_center(theta::Vector{Float64}, phase::Vector{Float64})
    # theta and phase in radians. phase should be unwrapped.
    # Set up least squares solution:
    soln = hcat(ones(length(theta)), sin.(theta), cos.(theta)) \ phase
    phi0 = soln[1]
    krho = soln[2]
    kz = soln[3]
    return (krho, kz)
end

eval_cut(cutfile::AbstractString, fghz::Real, thetamax::Real) = eval_cut(read_cutfile(cutfile), fghz, thetamax)

"""
    (c, xn, sp, slh, et, pc, xpd) = eval_cut(cutfile, fghz, thetamax)

Determine primary pattern performance metrics for a `Cut` object or 
Ticra cut file.

## Arguments:

- `cutfile`:   A string containing the name of the cut file to evaluate,
  or a `Cut` structure as returned by `read_cutfile`.
- `fghz`: The frequency in GHz.
- `thetamax`:  The maximum theta angle in degrees for which the pattern is
  to be evaluated.

# Returned Values:

- `c`: A scalar containing copol directivity in dbi evaluated at boresight.
- `xn`: An array containing normalized crosspol depths in dB > 0 evaluated at points 
  0 ≤ theta ≤ `thetamax`.  `xn` is defined as peak copol gain in dbi - crosspol gain in dBi.
- `sp`: The spillover loss in db > 0, representing power not captured by the reflector
  because it is outside the cone defined by half-angle `thetamax`.
- `slh`: Sidelobe height.  The sum of the heights of the first few sidelobes falling within 
  the cone of interest, measured in db above their immediately preceding relative minimum 
  (in copol amplitude) in the phi = 0 cut. The value returned is 0 if no sidelobe is detected. 
- `et`:  Edge Taper.  An array containing G(0,0) - G(`thetamax`,phi) for each point in the region
  |theta| = `thetamax`. G(theta,phi) is the copol directivity in dbi evaluated at observation 
  angle (theta,phi). 
- `pc`:  The z-coordinate of the phase center location in inches wrt the
  horn aperture plane.  It is calculated as the mean of the E-plane and
  H-plane phase centers evaluated for theta values corresponding to
  field amplitudes within 10 db of the beam peak.
- `xpd`: An array containing xpd (crosspol discrimination) values in db evaluated at points
  0 <= |theta| ≤ `thetamax`.  xpd = copol gain in dbi - crosspol gain in dBi
"""
function eval_cut(cut::Cut, fghz::Real, thetamax::Real)
    pwr_tot = power(cut)  # Total power in the cut

    pol1db = 20 * log10.(norm.(getindex.(get_evec(cut), 1)))
    pol2db = 20 * log10.(norm.(getindex.(get_evec(cut), 2)))

    if maximum(pol1db) > maximum(pol2db)
        copol = pol1db
        cxpol = pol2db
    else
        copol = pol2db
        cxpol = pol1db
    end

    clamp!(cxpol, -500.0, Inf) # Eliminate negative infinities

    # Check that cut is asymmetric:
    0.0 ∈ cut.theta || error("This cut does not contain theta = 0!")
    first(cut.theta) == 0.0 || error("eval_cut requires an asymmetric cut")
    thetamax > maximum(cut.theta) && error("thetamax is too large")

    c = copol[begin, begin] # Boresight copol

    ind = filter(i -> cut.theta[i] ≤ thetamax, 1:length(cut.theta))
    indm = maximum(ind)

    theta = collect(cut.theta[ind])  # Truncate to the thetamax cone

    # Interpolate to add another point at thetamax into the retained cut data
    h1 = thetamax - cut.theta[indm]
    if h1 > 0.001 && length(cut.theta) > length(ind)
        h = h1 / (cut.theta[indm+1] - cut.theta[indm])
        push!(theta, thetamax)
        conew = zeros(1, size(cut.evec, 2))
        cxnew = zeros(1, size(cut.evec, 2))
        for k = 1:size(cut.p1, 2)
            conew[1, k] = (1 - h) * copol[indm, k] + h * copol[indm+1, k]
            cxnew[1, k] = (1 - h) * cxpol[indm, k] + h * cxpol[indm+1, k]
        end
        copol = vcat(copol[ind, :], conew, copol[indm+1:end, :])
        cxpol = vcat(cxpol[ind, :], cxnew, copol[indm+1:end, :])
        push!(ind, indm + 1)
    end

    # Edge taper:
    et = c .- vec(copol[ind[end], :])
    amp = 10 * log10.(abs2.(copol[ind, 1]))  # Phi = 0 cut
    indc = 2:length(amp)-1  # indices of central points of phi = 0 cut.

    # Locate indices of relative maxima:
    ind_lmax = 1 .+ filter(i -> amp[i] > amp[i-1] && amp[i] > amp[i+1], indc)
    # Locate indices of relative minima:
    ind_lmin = 1 .+ filter(i -> amp[i] < amp[i-1] && amp[i] < amp[i+1], indc)
    isempty(ind_lmin) && !isempty(ind_lmax) && (ind_lmin = ind_lmax)
    length(ind_lmin) < length(ind_lmax) && amp[1] < amp[ind_lmax[1]] && pushfirst!(ind_lmin, 1)

    slh = 0.0  # No sidelobe inside cone of interest.
    for k = 1:length(ind_lmax)
        slh = slh + amp[ind_lmax[k]] - amp[ind_lmin[k]]
    end

    xn = c .- vec(cxpol[ind, :])

    xpd = vec(copol[ind, :] - cxpol[ind, :])

    # Zero out the power in the cone theta <= thetamax:
    cut2theta = range(cut.theta[1], cut.theta[indm], length=length(ind))  # Truncate to the thetamax cone
    evec = get_evec(cut)[ind, :]
    cut2 = Cut(cut.ncomp, cut.icut, cut.icomp, cut.text, cut2theta, cut.phi, evec)
    pwr_beam = power(cut2)
    sp = 10 * log10(pwr_tot / pwr_beam)

    (x, y, z0, z90) = phscen(cut2, fghz)
    pc = 0.5 * (z0 + z90)
    return (c, xn, sp, slh, et, pc, xpd)
end

"""
    normalize2dir!(cut::Cut)

Normalize a Cut object so it's total power is 4π.  This results
in field magnitude squared being directivity.
"""
function normalize2dir!(cut::Cut)
    pwr = power(cut)
    c = sqrt(4π/pwr)
    @inbounds for i in eachindex(cut.evec)
        cut.evec[i] *= c
    end
    return cut
end

"""
    convert_cut!(cut::Cut, icomp::Integer)

Convert `cut` to a new polarization basis determined by `icomp`.  Legal values 
for `icomp` and their meanings:
*  1 => Eθ and Eϕ
*  2 => ERHCP and ELHCP
*  3 => Eh and Ev (Ludwig 3 co and cx)
"""
function convert_cut!(cut::Cut{Tc,N}, icomp::Integer) where {Tc,N}
    (icomp < 1 || icomp > 3) && throw(ArgumentError("icomp is not 1, 2, or 3"))
    outpol = icomp
    get_ncomp(cut) == 2 || error("Only ncomp == 2 allowed")
    (inpol = get_icomp(cut)) == outpol && return
    evec = get_evec(cut)
    @inbounds for (col, phi) in enumerate(get_phi(cut))        
        p̂s = _pol_basis_vectors(phi)
        p̂1in, p̂2in = p̂s[inpol]
        p̂1out, p̂2out = p̂s[outpol]
        if N == 2
            mat = @SMatrix [(p̂1out ⋅ p̂1in)  (p̂1out ⋅ p̂2in)
                            (p̂2out ⋅ p̂1in)  (p̂2out ⋅ p̂2in)]
        elseif N == 3
            mat = @SMatrix [(p̂1out ⋅ p̂1in)  (p̂1out ⋅ p̂2in) 0 
                            (p̂2out ⋅ p̂1in)  (p̂2out ⋅ p̂2in) 0
                                 0              0          1]
        else
            error("Unknown type parameter $N")
        end
    for row in 1:length(get_theta(cut))
            evec[row,col] = mat * evec[row,col]
        end
    end
    cut.icomp = outpol
    return
end

const iroot2 = inv(sqrt(2))

"""
    _pol_basis_vectors(ϕ::Real)
Return a tuple containing the three sets of unitary polarization basis vectors:
`((θ̂(ϕ), ϕ̂(ϕ)), (R̂(ϕ), L̂(ϕ)), (ĥ(ϕ), v̂(ϕ)))`.

The input argument `ϕ` is the azimuthal angle in *degrees*. The returned basis vectors
are of type `SVector{2, ComplexF64}` and are represented as (θ̂, ϕ̂) components.  E.g.,
in this representation `θ̂ == [1,0]` and `ϕ̂ == [0,1]`.
"""
function _pol_basis_vectors(ϕ::Real)
    sinϕ, cosϕ = sincosd(ϕ)
    θ̂ = @SVector [1.0+0im, 0.0+0.0im]
    ϕ̂ = @SVector [0.0+0.0im, 1.0+0im]
    ĥ = θ̂*cosϕ - ϕ̂*sinϕ
    v̂ = θ̂*sinϕ + ϕ̂*cosϕ
    R̂ = iroot2 * (ĥ - im*v̂)
    L̂ = iroot2 * (ĥ + im*v̂)
    return ((θ̂, ϕ̂), (R̂, L̂), (ĥ, v̂))
end

function _add_3rd_component(cut::Cut)
    ncomp_old = get_ncomp(cut)
    ncomp_old == 1 && error("Can't handle one-component cut")
    ncomp_old == 3 && return deepcopy(cut)
    evec_old = get_evec(cut)
    evec = [@SVector[t[1], t[2], rand(ComplexF64)] for t in evec_old]
    cut_new = Cut(3, get_icut(cut), get_icomp(cut), get_text(cut), 
                       get_theta(cut), get_phi(cut), evec)
    return cut_new
end

_norm²(E::SVector{N,T}) where {N,T} = sum(abs2, E[n] for n in 1:N)

"""
    sor_efficiency(cut; F, D, Oc, pol=:matched, dz=0.0) -> (ηₛ, ηᵢₗ, ηₚ, ηₓ)

Compute boresight directivity efficiency of a parabolic single-offset reflector using the feed pattern
specified by `cut`.

## Positional input arguments:
* `cut`: Either a string containing the name of a Ticra-compatible, spherical,
  polar, asymmetric cut file, or a `Cut` object as returned by the `read_cutfile`
  function.

## Keyword input arguments:
* `F`, `D`, and `Oc`:  The single-offset reflector focal length, aperture diameter, and
  center offset, respectively.  These may be expressed in any convenient length units, so 
  long as they are consistent. 
* `pol`: A `Symbol` having one of the values `:L3h`, `:L3v`, `:RHCP`, `:LHCP`, or :matched
  (any variations in terms of capitalization are acceptable).  The first two denote Ludwig 3 
  horizontal (x) and vertical (y) polarizations, the second two denote the two senses of 
  circular polarization, and `:max` (the default) uses the polarization among the 4 previously
  listed that that has the maximum field amplitude.  Polarization efficiency of the boresight 
  secondary pattern will be computed relative to the polarization specified by this argument.
* `dz`: This is a signed distance along the zfeed direction, measured in wavelengths. It allows for
  repositioning the feed in an attempt to locate the feed phase center at the reflector focal
  point. Suppose that the feed's phase center is located 0.42 wavelengths inside the horn aperture. The
  horn origin should ideally be positioned closer to the reflector, so a positive value `dz = 0.42` would
  be specified to indicate that the horn has been repositioned in this manner.  

## Return values:
* `ηₛ`: Spillover efficiency (a real number between 0 and 1).
* `ηᵢₗ`: Illumination (amplitude taper) efficiency (a real number between 0 and 1).
* `ηₚ`:  Phase error efficiency (a real number between 0 and 1).
* `ηₓ`: Polarization efficiency (a real number between 0 and 1).
"""
sor_efficiency(cut::String; kwargs...) = sor_efficiency(read_cutfile(cut); kwargs...)

function sor_efficiency(cut::Cut; pol::Symbol=:max, F::Number, D::Number, Oc::Number, dz::Real=0.0)
    iszero(first(get_theta(cut))) || error("Only asymmetric cut files are supported")
    pol = Symbol(lowercase(string(pol)))
    pol in (:lhcp, :rhcp, :l3h, :l3v, :max) || error("Illegal pol value: $pol")
    cut = deepcopy(cut) # Avoid modifying input object
    # Compute β (feed tilt angle) and θe (edge ray angle):
    xm = Oc + D/2  # Upper edge
    zm = xm^2 / 4F
    θU = atan(xm / (F - zm))  # Upper edge angle
    xm = Oc - D/2  # Lower edge
    zm = xm^2 / 4F
    θL = atan(xm / (F-zm))  # Lower edge angle
    β = 0.5 * (θU + θL)  # Bisector angle (radians)
    sβ, cβ = sincos(β)
    θe = 0.5 * (θU - θL)  # Edge angle (radians)
    
    pwr = power(cut) # Total radiated power
    θmax = θe
    θmaxdeg = rad2deg(θmax)
    θmaxdeg > last(get_theta(cut)) && error("Cut does not extend to edges of reflector!")
    pwr_cone = power(cut, θmaxdeg)


    ηₛ = pwr_cone / pwr # Spillover efficiency

    # Normalize the fields to 4*pi:
    factor = sqrt(4π / pwr)
    cut.evec .*= factor

    # Get needed trig functions
    scθ = sincosd.(get_theta(cut))
    scϕ = sincosd.(get_phi(cut))

    # Adjust phase to account for translation of dz along z-feed axis:
    if !iszero(dz)
        cfactor = @. cis(2π * dz * last(scθ))
        cut.evec .*= cfactor
    end


    # Convert fields to θ/ϕ components:
    convert_cut!(cut, 1)
    evec = get_evec(cut) # Save actual fields

    if pol == :max
        # Find the maximum norm E-field
        _, imaxnorm = findmax(_norm², evec)
        Eθϕ_maxnorm = evec[imaxnorm]
        ϕ_maxnorm = get_phi(cut)[last(Tuple(imaxnorm))]
        # Check which polarization decomposition produces largest copol:
        Eθϕ = Eθϕ_maxnorm
        (θ̂,ϕ̂), (R̂,L̂), (ĥ,v̂) = _pol_basis_vectors(ϕ_maxnorm)
        ERL = @SMatrix[R̂ ⋅ θ̂   R̂ ⋅ ϕ̂; L̂ ⋅ θ̂   L̂ ⋅ ϕ̂] * Eθϕ
        Ehv = @SMatrix[ĥ ⋅ θ̂   ĥ ⋅ ϕ̂; v̂ ⋅ θ̂   v̂ ⋅ ϕ̂] * Eθϕ
        icompm1 = argmax(maximum(abs2.(x)) for x in (ERL, Ehv))
        Emaxsq12 = abs2.((ERL, Ehv)[icompm1])
        if icompm1 == 1
            # CP components: Sense will reverse for secondary pattern!
            if Emaxsq12[1] > Emaxsq12[2] 
                pol = :lhcp     
            else
                pol = :rhcp
            end
        else 
            # Ludwig 3 components
            if Emaxsq12[1] > Emaxsq12[2] 
                pol = :l3h
            else
                pol = :l3v
            end
        end
    end

    # Compute theta and phi components due to zero phase error copol and crosspol:
    evec0 = copy(evec) # Copy far-field vectors
    cut.evec = evec0
    if pol in (:lhcp, :rhcp)
        convert_cut!(cut, 2) # CP components
    elseif pol in (:l3h, :l3v)
        convert_cut!(cut, 3) # Ludwig 3 components
    end
    for i in eachindex(evec0)
        evec0[i] = abs.(evec0[i])
    end
    convert_cut!(cut, 1) # Back to theta, phi components
    cut.evec = evec # Restore actual fields

    # Storage for ϕ integrals at each θ location
    integrand, integrand0 = (zeros(SVector{2,ComplexF64}, length(get_theta(cut))) for _ in 1:2)
    
    # Compute ϕ integrals:
    for (iθ, θ) in pairs(get_theta(cut))
        sθ, cθ = scθ[iθ]
        θ == 180 && ((sθ, cθ)  = (0.0001, -0.9999)) # Avoid irrelevant singularity at 180°
        for (iϕ, ϕ) in pairs(get_phi(cut))
            eθ, eϕ = evec[iθ, iϕ]
            e0θ, e0ϕ = evec0[iθ, iϕ]
            sϕ, cϕ = scϕ[iϕ]
            common = sθ / (1 + cθ*cβ - sθ*cϕ*sβ)^2
            f1 = cϕ * (1 + cβ*cθ) - sβ*sθ
            f2 = -sϕ * (cβ + cθ)
            integrand[iθ] += common * @SVector [f1 * eθ  + f2 * eϕ, f2 * eθ  -  f1 * eϕ]
            integrand0[iθ] +=  common * @SVector [f1 * e0θ + f2 * e0ϕ, f2 * e0θ - f1 * e0ϕ]
        end
    end
    phi = get_phi(cut)
    dphi = deg2rad(phi[2] - phi[1])
    integrand .*= dphi
    integrand0 .*= dphi

    # Perform integration over theta:
    thetarad = deg2rad.(get_theta(cut))
    spl = CubicSpline(integrand, thetarad; assume_linear_t=true)
    Ivec, errest1 = quadgk(spl, 0, θmax; rtol=1e-10)
    spl0 = CubicSpline(integrand0, thetarad; assume_linear_t=true)
    I0vec, errest2 = quadgk(spl0, 0, θmax; rtol=1e-10)

    # Compute efficiency vectors:
    factor = 2/π * F/D 
    Ivec *= factor
    I0vec *= factor # Spillover-Illum only efficiency vector
    ηₛηᵢₗ = _norm²(I0vec)  # Product of spillover and illumination effics
    ηᵢₗ = ηₛηᵢₗ / ηₛ

    # Compute polarization unitary vector:
    if pol == :l3h
        uhat = @SVector [one(ComplexF64), zero(ComplexF64)]
    elseif pol == :l3v
        uhat = @SVector [zero(ComplexF64), one(ComplexF64)]
    elseif pol == :lhcp
        uhat = @SVector [iroot2, im*iroot2]
    elseif pol == :rhcp
        uhat = @SVector [iroot2, -im*iroot2]
    else
        error("Illegal pol value: $pol")
    end 

    # Polarization efficiency:
    ηₓ = abs2(uhat ⋅  Ivec) / _norm²(Ivec)

    # Phase error efficiency:
    ηₚ = _norm²(Ivec) / ηₛηᵢₗ

    return (;ηₛ, ηᵢₗ, ηₚ, ηₓ)
end

