using Printf: @printf
using DSP: unwrap
using Dierckx: Spline1D, integrate
using StaticArrays: @SVector, @SMatrix
using LinearAlgebra: ⋅


"""
`TicraCut` Holds data for a Ticra "Tabulated Pattern Data" object.
Note that a single `TicraCut` instance contains all of the cuts for a single frequency.

### Fields

* `ncomp::Int`: Number of polarization components (2 or 3)
* `icut::Int`: 1 for standard constant ϕ polar cuts, 2 for conical, constant θ cuts1
* `icomp::Int`: Polarization control parameter. 1 for Eθ and Eφ, 2 for RHCP and LHCP, 3 for Co and Cx (Ludwig 3).
* `text::Vector{String}`: Identification text for each constant angle cut.
* `theta::T<:AbstractRange`: The theta values (in degrees) stored in the cut.
* `phi::T<:AbstractRange`: The phi values (in degrees) stored in the cut.
* `p1`, `p2`, `p3`: Matrices of complex field values for the three possible polarization components.
"""
mutable struct TicraCut{T <: AbstractRange}
    ncomp::Int
    icut::Int
    icomp::Int
    text::Vector{String}
    theta::T
    phi::T
    p1::Matrix{ComplexF64}
    p2::Matrix{ComplexF64}
    p3::Matrix{ComplexF64}
end

TicraCut() = TicraCut(
    0,
    1,
    0,
    String[],
    0.0:-1.0:10.0,
    0.0:-1.0:10.0,
    zeros(ComplexF64, 0, 0),
    zeros(ComplexF64, 0, 0),
    zeros(ComplexF64, 0, 0),
)

import Base.show
function show(io::IO, t::TicraCut)
    println(io, "TicraCut")
    println(io, "  ncomp\t$(t.ncomp)")
    println(io, "  icut \t$(t.icut)")
    println(io, "  icomp\t$(t.icomp)")
    print(io, "  text")
    if isempty(t.text)
        println(io, "")
    elseif length(t.text) < 6
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
    println(io, "  p1   \t$(summary(t.p1))")
    println(io, "  p2   \t$(summary(t.p2))")
    t.ncomp > 2 && println(io, "  p3   \t$(summary(t.p3))")
    return nothing
end

"""
    get_theta(c::TicraCut)

Return the theta values (in degrees) stored in the cut
"""
get_theta(c::TicraCut) = c.theta

"""
    get_phi(c::TicraCut)

Return the phi values (degrees) stored in the cut.
"""
get_phi(c::TicraCut) = c.phi

"""
    get_p1(c::TicraCut)

Return the ntheta × nphi matrix of complex field values stored in polarization slot 1 of the cut.
"""
get_p1(c::TicraCut) = c.p1

"""
    get_p2(c::TicraCut)

Return the ntheta × nphi matrix of complex field values stored in polarization slot 2 of the cut.
"""
get_p2(c::TicraCut) = c.p2

"""
    get_p3(c::TicraCut)

Return the ntheta × nphi matrix of complex field values stored in polarization slot 3 of the cut.
"""
get_p3(c::TicraCut) = c.p3

"""
    get_ncomp(c::TicraCut)

Return ncomp, the number of polarization components stored in the cut.
"""
get_ncomp(c::TicraCut) = c.ncomp

"""
    get_icut(c::TicraCut)

Return icut, the control parameter of the cut. 1 for a polar cut, 2 for a conical cut.
"""
get_icut(c::TicraCut) = c.icut

"""
    get_icomp(c::TicraCut)

Return icomp, polarization parameter. 1 for Eθ and Eφ, 2 for RHCP and LHCP, 3 for Co and Cx (Ludwig 3).
"""
get_icomp(c::TicraCut) = c.icomp


"""
    get_text(c::TicraCut)

Return a vector of strings containing the cut identification text.
"""
get_text(c::TicraCut) = c.text

"""
    amplitude_db(c::TicraCut, ipol::Int)
    amplitude_db(c::TicraCut, polstr::String = "copol")

Return a matrix of amplitudes in dB for some choice of polarization in the cut.
Legal values for `ipol` are 1, 2, or 3.  Legal values for `polstr` are
"copol" (the default) and "crosspol".
"""
function amplitude_db end

function amplitude_db(c::TicraCut, ipol::Int)
    if ipol == 1
        return 10 * log10.(abs2.(c.p1))
    elseif ipol == 2
        return 10 * log10.(abs2.(c.p2))
    elseif ipol == 3
        return 10 * log10.(abs2.(c.p3))
    else
        error("Illegal value: $ipol for ipol")
    end
end

function amplitude_db(c::TicraCut, polstr::String = "copol")
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
        p1 = abs2.(c.p1)
        p1maxsq = maximum(p1)
        p2 = abs2.(c.p2)
        p2maxsq = maximum(p2)
        if (p1maxsq > p2maxsq && maxflag) || (p1maxsq < p2maxsq && minflag)
            return 10 * log10.(p1)
        else
            return 10 * log10.(p2)
        end
    end
end


"""
    phase_deg(c::TicraCut, ipol::Int)
    phase_deg(c::TicraCut, polstr::String = "copol")

Return a matrix of phases in degrees for some choice of polarization in the cut.
Legal values for `ipol` are 1, 2, or 3.  Legal values for `polstr` are
"copol" (the default) and "crosspol".
"""
function phase_deg end

function phase_deg(c::TicraCut, ipol::Int)
    if ipol == 1
        return rad2deg.(angle.(c.p1))
    elseif ipol == 2
        return rad2deg.(angle.(c.p2))
    elseif ipol == 3
        return rad2deg.(angle.(c.p3))
    else
        error("Illegal value: $ipol for ipol")
    end
end

function phase_deg(c::TicraCut, polstr::String = "copol")
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
        p1maxsq = maximum(abs2, c.p1)
        p2maxsq = maximum(abs2, c.p2)
        if (p1maxsq > p2maxsq && maxflag) || (p1maxsq < p2maxsq && minflag)
            return rad2deg.(angle.(c.p1))
        else
            return rad2deg.(angle.(c.p2))
        end
    end
end


"""
    power(cut::Ticracut)::Float64
    
Compute the total radiated power in a TicraCut object. 
If only a single phi value is included in the cut, then assume no azimuthal variation.
"""
function power(cut::TicraCut)::Float64
    # This version uses the trapezoidal rule in phi and integration of a
    # cubic spline interpolant in the theta direction.
    get_ncomp(cut) == 3 && error("Cut has 3 field components.  Only 2 allowed.")
    sym = cut.theta[begin] < 0  # Symmetrical cut
    phifullmax = sym ? 180.0 : 360.0
    ntheta = length(get_theta(cut))
    phi = get_phi(cut)
    nphi = length(phi)
    p = vec(sum(abs2, get_p1(cut), dims = 2) + sum(abs2, get_p2(cut), dims = 2))
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
    pwr = phimult * integrate(spl, theta[1], theta[end])
    return pwr
end


"""
    read_ticra_cut(fname, warnflag=true)::TicraCut

Read the first frequency's data from a Ticra-compatible cut file.  
`warnflag`, if `true`, causes this function to issue a warning
if more than one frequency is detected in the file.
"""
function read_ticra_cut(fname::AbstractString, warnflag::Bool = true)
    # Read a single frequency from a (possibly multi-frequency) cut file.
    # warnflag, if true, instructs this routine to issue a warning
    # if more than one frequency is present in the file.
    t = read_ticra_cuts(fname)
    n = length(t)
    n > 1 && warnflag && @warn "$n frequencies found in $fname.  Returning only 1st..."
    return t[1]
end

"""
    read_ticra_cuts(fname)::Vector{TicraCut}  

Read data from a possibly multi-frequency Ticra-compatible cut file.  
Return a vector of one or more `TicraCut` structs.
"""
function read_ticra_cuts(fname::AbstractString)
    cuts = open(fname, "r") do fid
        cutphi = Float64[]
        header = String[]
        kf = 0  # Initialize frequency counter
        local dth, header, icomp, icut, kf, ncomp, nphi, nth, p1allphi, p2allphi, p3allphi, pnext, ths
        cuts = TicraCut[]
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
                    p1 = reshape(p1allphi, length(cut.theta), length(cut.phi))
                    p2 = reshape(p2allphi, length(cut.theta), length(cut.phi))
                    if ncomp == 3
                        p3 = reshape(p3allphi, length(cut.theta), length(cut.phi))
                    else
                        p3 = Matrix{ComplexF64}(undef, 0, 0)
                    end
                    cuts[end] = TicraCut(cut.ncomp, cut.icut, cut.icomp, cut.text, cut.theta, cut.phi, p1, p2, p3)
                else
                    pnext = zeros(ComplexF64, nth) # Storage increment for reading field values
                end
                kf += 1
                header = String[textline]
                cutphi = phi:phi
                theta = range(start = ths, step = dth, length = nth)
                p1allphi = zeros(ComplexF64, 0)
                p2allphi = zeros(ComplexF64, 0)
                ncomp == 3 && (p3allphi = zeros(ComplexF64, 0))
                p1 = p2 = p3 = zeros(ComplexF64,0,0)
                cut = TicraCut(ncomp, icut, icomp, header, theta, cutphi, p1, p2, p3)
                push!(cuts, cut)
            else
                # Begin an additional phi cut at current frequency:
                cut = cuts[end]
                cutphi = first(cut.phi):(phi-last(cut.phi)):phi
                push!(cut.text, textline)
                cuts[end] = TicraCut(cut.ncomp, cut.icut, cut.icomp, cut.text, cut.theta, cutphi, cut.p1, cut.p2, cut.p3)
            end
            # Check consistency
            cut = cuts[end]
            @assert (icomp, icut, ncomp) == (cut.icomp, cut.icut, cut.ncomp)
            @assert (ths, dth, nth) == (first(cut.theta), cut.theta[2]-cut.theta[1], length(cut.theta)) 
            append!(p1allphi, pnext)
            p1thisphi = @view p1allphi[end-nth+1:end]
            append!(p2allphi, pnext)
            p2thisphi = @view p2allphi[end-nth+1:end]
            if ncomp == 3
                append!(p3allphi, pnext)
                p3thisphi = @view p3allphi[end-nth+1:end]
            end

            # Read in the field data for this phi cut
            for i = 1:nth
                str = readline(fid)
                if ncomp == 2
                    (t1, t2, t3, t4) = (parse(Float64, s) for s in split(str))
                else
                    (t1, t2, t3, t4, t5, t6) = (parse(Float64, s) for s in split(str))
                    p3thisphi[i] = complex(t5, t6)
                end
                p1thisphi[i] = complex(t1, t2)
                p2thisphi[i] = complex(t3, t4)
            end
        end # while
        # Finish final cut
        cut = cuts[end]
        nphi = length(cut.phi)
        nth = length(cut.theta)
        p1 = reshape(p1allphi, nth, nphi)
        p2 = reshape(p2allphi, nth, nphi)
        if ncomp == 3
            p3 = reshape(p3allphi, nth, nphi)
        else
            p3 = Matrix{ComplexF64}(undef, 0, 0)
        end
        cut = cuts[end]
        cuts[end] = TicraCut(cut.ncomp, cut.icut, cut.icomp, cut.text, cut.theta, cut.phi, p1, p2, p3)

        return cuts
    end # function
    return cuts
end

"""
    write_ticra_cut(fname::AbstractString, cut::TicraCut, title::AbstractString="Cut file created by write_ticra_cut")

    write_ticra_cut(fname::AbstractString, cuts::AbstractVector{TicraCut}, title::AbstractString="Cut file created by write_ticra_cut")

Write `TicraCut` cut data to a Ticra-compatible cut file.
"""
function write_ticra_cut(
    fname::AbstractString,
    cut::TicraCut,
    title::String = "Cut file created by write_ticra_cut",
)
    open(fname, "w") do fid
        for (n, phi) in enumerate(cut.phi)
            @printf(fid, "%s, phi = %8.3f\n", title, phi)
            @printf(
                fid,
                "%8.3f %8.3f %d %8.3f %d %d %d\n",
                cut.theta[1],
                cut.theta[2] - cut.theta[1],
                length(cut.theta),
                phi,
                cut.icomp,
                cut.icut,
                cut.ncomp
            )
            for (m, theta) in enumerate(cut.theta)
                e1 = cut.p1[m, n]
                e2 = cut.p2[m, n]
                @printf(
                    fid,
                    " %18.10E %18.10E %18.10E %18.10E\n",
                    real(e1),
                    imag(e1),
                    real(e2),
                    imag(e2)
                )
            end
        end
    end
end

function write_ticra_cut(
    fname::AbstractString,
    cuts::AbstractVector{TicraCut},
    title::String = "Cut file created by write_ticra_cut",
    )
    open(fname, "w") do fid
        for cut in cuts
            for (n, phi) in enumerate(cut.phi)
                @printf(fid, "%s, phi = %8.3f\n", title, phi)
                @printf(
                    fid,
                    "%8.3f %8.3f %d %8.3f %d %d %d\n",
                    cut.theta[1],
                    cut.theta[2] - cut.theta[1],
                    length(cut.theta),
                    phi,
                    cut.icomp,
                    cut.icut,
                    cut.ncomp
                )
                for (m, theta) in enumerate(cut.theta)
                    e1 = cut.p1[m, n]
                    e2 = cut.p2[m, n]
                    @printf(
                        fid,
                        " %18.10E %18.10E %18.10E %18.10E\n",
                        real(e1),
                        imag(e1),
                        real(e2),
                        imag(e2)
                    )
                end
            end
        end
    end
end


function phscen(cutfile::AbstractString, fghz = 11.802852677165355; min_dropoff = -10)
    cut = read_ticra_cut(cutfile)
    phscen(cut, fghz; min_dropoff)
end


"""
    (x,y,z0,z90) = phscen(cutfile::AbstractString, fghz=11.802852677165355; min_dropoff=-10)

    (x,y,z0,z90) = phscen(cut::TicraCut, fghz=11.802852677165355; min_dropoff=-10)
    
Estimate phase center for a cut file or `TicraCut` object using NSI least-squares 
algorithm.

The three outputs are estimates of the location of the phase center relative to
the origin used in recording the data in the cut.  If `fghz` is passed in as 
an argument, the values will be expressed in units of inches.  Otherwise the 
lengths will be normalized to wavelengths. Note that `z0` and `z90` are the 
phi = 0ᵒ and phi = 90ᵒ plane estimates of the phase center location.  
In determining the phase center locations, only field values with magnitudes
in dB greater than `min_dropoff` relative to the peak field are considered.
"""
function phscen(cut::TicraCut, fghz = 11.802852677165355; min_dropoff = -10.0)
    # Determine which pol slot has main polarization:
    p1max = maximum(abs2, get_p1(cut)) |> sqrt
    p2max = maximum(abs2, get_p2(cut)) |> sqrt
    if p1max > p2max
        pmax = p1max
        p = 1
        E = get_p1(cut)
    else
        pmax = p2max
        p = 2
        E = get_p2(cut)
    end

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
        (krho[k], kz[k]) = find_phase_center(x[big_enough], y[big_enough])

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


function find_phase_center(theta::Vector{Float64}, phase::Vector{Float64})
    # theta and phase in radians. phase should be unwrapped.
    # Set up least squares solution:
    soln = hcat(ones(length(theta)), sin.(theta), cos.(theta)) \ phase
    phi0 = soln[1]
    krho = soln[2]
    kz = soln[3]
    return (krho, kz)
end

eval_cut(cutfile::AbstractString, fghz::Real, thetamax::Real) = eval_cut(read_ticra_cut(cutfile), fghz, thetamax)

"""
    (c, xn, sp, slh, et, pc, xpd) = eval_cut(cutfile, fghz, thetamax)

Determine primary pattern performance metrics for a `TicraCut` object or 
Ticra cut file.

## Arguments:

- `cutfile`:   A string containing the name of the cut file to evaluate,
  or a `TicraCut` structure as returned by `read_ticra_cut`.
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
function eval_cut(cut::TicraCut, fghz::Real, thetamax::Real)
    pwr_tot = power(cut)  # Total power in the cut

    pol1db = 10 * log10.(abs2.(cut.p1))
    pol2db = 10 * log10.(abs2.(cut.p2))

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
        conew = zeros(1, size(cut.p1, 2))
        cxnew = zeros(1, size(cut.p1, 2))
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
    p1 = cut.p1[ind, :]
    p2 = cut.p2[ind, :]
    cut2 = TicraCut(cut.ncomp, cut.icut, cut.icomp, cut.text, cut2theta, cut.phi, p1, p2, cut.p3)
    pwr_beam = power(cut2)
    sp = 10 * log10(pwr_tot / pwr_beam)

    (x, y, z0, z90) = phscen(cut2, fghz)
    pc = 0.5 * (z0 + z90)
    return (c, xn, sp, slh, et, pc, xpd)
end


"""
    convert_cut!(cut::TicraCut, icomp::Integer)

Convert `cut` to a new polarization basis determined by `icomp`.  Legal values 
for `icomp` and their meanings:
*  1 => Eθ and Eϕ
*  2 => ERHCP and ELHCP
*  3 => Eh and Ev (Ludwig 3 co and cx)
"""
function convert_cut!(cut, icomp)
    (icomp < 1 || icomp > 3) && throw(ArgumentError("icomp is not 1, 2, or 3"))
    outpol = icomp
    get_ncomp(cut) == 2 || error("Only ncomp == 2 allowed")
    (inpol = get_icomp(cut)) == outpol && return
    for (col, phi) in enumerate(get_phi(cut))
        sp, cp = sincosd(phi)
        p̂s = _pol_basis_vectors(sp, cp)
        for row in 1:length(get_theta(cut))
            Ein = @SVector [cut.p1[row,col], cut.p2[row,col]]
            p̂1in, p̂2in = p̂s[inpol]
            p̂1out, p̂2out = p̂s[outpol]
            mat = @SMatrix [(p̂1out ⋅ p̂1in)  (p̂1out ⋅ p̂2in)
                            (p̂2out ⋅ p̂1in)  (p̂2out ⋅ p̂2in)]
            Eout = mat * Ein
            cut.p1[row,col] = Eout[1]
            cut.p2[row,col] = Eout[2]
        end
    end
    cut.icomp = outpol
    return
end

const iroot2 = inv(sqrt(2))
function _pol_basis_vectors(sinϕ, cosϕ)
    θ̂ = @SVector [1.0+0im, 0.0+0.0im]
    ϕ̂ = @SVector [0.0+0.0im, 1.0+0im]
    ĥ = θ̂*cosϕ - ϕ̂*sinϕ
    v̂ = θ̂*sinϕ + ϕ̂*cosϕ
    R̂ = iroot2 * (ĥ - im*v̂)
    L̂ = iroot2 * (ĥ + im*v̂)
    return ((θ̂, ϕ̂), (R̂, L̂), (ĥ, v̂))
end
