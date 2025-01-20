using Printf: @printf
using DSP: unwrap
using DataInterpolations: CubicSpline
using QuadGK: quadgk
using StaticArrays: @SVector, @SMatrix, SVector
using LinearAlgebra: ⋅, norm
using RecipesBase: @recipe, @series

"""
    Cut
 
Contains data for a Ticra "Tabulated Pattern Data" object. Note that a single `Cut` instance 
contains all of the cuts for a single frequency.

### Fields

* `ncomp::Int`: Number of polarization components (2 or 3)
* `icut::Int`: 1 for standard constant ϕ polar cuts, 2 for conical, constant θ cuts.
* `icomp::Int`: Polarization control parameter. 1 for Eθ and Eφ, 2 for RHCP and LHCP, 3 for Co and Cx (Ludwig 3).
* `text::Vector{String}`: Identification text for each constant angle cut.
* `theta::Tt<:AbstractRange`: The theta values (in degrees) stored in the cut.
* `phi::Tp<:AbstractRange`: The phi values (in degrees) stored in the cut.
* `evec`: Matrix of complex field vectors for the two or three polarization components.
"""
@kwdef mutable struct Cut{Tt<:AbstractRange, Tp<:AbstractRange, N}
    ncomp::Int
    icut::Int
    icomp::Int
    text::Vector{String}
    theta::Tt
    phi::Tp
    evec::Matrix{SVector{N,ComplexF64}}
end


function Base.show(io::IO, mime::MIME"text/plain", t::Cut)
    println(io, "Cut")
    println(io, "  ncomp\t$(t.ncomp)")
    println(io, "  icut \t$(t.icut)")
    println(io, "  icomp\t$(t.icomp)")
    print(io, "  text")
    if isempty(t.text)
        println(io, "")
    elseif length(t.text) < 6
        println(io, "")
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

function Base.show(io::IO, t::Cut)
    print(io, "Cut with ncomp=$(t.ncomp), icomp=$(t.icomp), phi=$(t.phi), theta=$(t.theta)")
    return nothing
end

"""
    isapprox(c1::Cut, c2::Cut; kwargs...) -> tf::Bool

Compare two `Cut` objects for approximate equality.

Compares most fields for perfect equality except `text` and `evec`.
The `text` fields are not compared at all, and the `evec` fields 
(`Matrix` types) are compared for approximate equality using `isapprox`.
"""
function Base.isapprox(c1::Cut, c2::Cut; kwargs...)
    c1.ncomp == c2.ncomp || return false
    c1.icut == c2.icut || return false
    c1.icomp == c2.icomp || return false
    c1.theta == c2.theta || return false
    c1.phi == c2.phi || return false
    return isapprox(c1.evec, c2.evec; kwargs...)
end

"""
    maximum(cut::Cut) -> maxE

Return the maximum amplitude for any polarization component stored in the `Cut` object.
"""
Base.maximum(cut::Cut) = maximum(x -> maximum(abs, x), cut.evec)


"""
    maximum_db(cut::Cut)

Return the maximum amplitude in dB for any polarization component stored in the `Cut` object.
"""
maximum_db(cut::Cut) = 20 * log10(maximum(cut))


"""
    get_theta(c::Cut)

Return the theta values (in degrees) stored in the cut.  The returned value will be 
some type of [`AbstractRange`](https://docs.julialang.org/en/v1/base/math/#Base.range) object.
"""
get_theta(c::Cut) = c.theta

"""
    get_phi(c::Cut)

Return the phi values (degrees) stored in the cut. The returned value will be 
some type of [`AbstractRange`](https://docs.julialang.org/en/v1/base/math/#Base.range) object.
"""
get_phi(c::Cut) = c.phi

"""
    get_evec(c::Cut)

Return the ntheta × nphi matrix of complex field vectors stored in the cut. Each element
of the returned matrix will be either a 2-vector or 3-vector, depending on the number
of field components stored in the cut.
"""
get_evec(c::Cut) = c.evec

"""
    get_evec(c::Cut, ipol::Integer)

Return the ntheta × nphi matrix of complex numbers stored in polarization slot `ipol` of 
the cut. `ipol` must be positive and less than or equal to `get_ncomp(cut)`.
"""
function get_evec(c::Cut, ipol::Integer) 
    ncomp = get_ncomp(c)
    1 ≤ ipol ≤ ncomp || throw(ArgumentError("ipol is not between 1 and $(ncomp)"))
    return getindex.(get_evec(c), ipol)
end

"""
    get_ncomp(c::Cut)

Return ncomp, the number of polarization components stored in the cut.
"""
get_ncomp(c::Cut) = c.ncomp

"""
    get_icut(c::Cut)

Return icut, the control parameter of the cut. 1 for a polar cut, 2 for a conical cut.
`TicraUtilities` currently accommodates only `icut == 1`, wherein each cut is for a constant ϕ value.
"""
get_icut(c::Cut) = c.icut

"""
    get_icomp(c::Cut)

Return icomp, polarization parameter. 1 for Eθ and Eφ, 2 for RHCP and LHCP, 3 for h and v (Ludwig 3).
"""
get_icomp(c::Cut) = c.icomp


"""
    get_text(c::Cut)

Return a vector of strings containing the cut identification text.
"""
get_text(c::Cut) = c.text

"""
    amplitude_db(c::Cut, ipol::Int)
    amplitude_db(c::Cut, polsymb::Symbol = :copol)
    amplitude_db(c::Cut, polstr::String = "copol")

Return a matrix of amplitudes in dB for some choice of polarization in the cut.
Legal values for `ipol` are 1, 2, or 3, the latter only being legal if there are
three polarization components present in the cut.  Legal values for `polstr` are
"copol" (the default) and "xpol".  Capitalization is not significant. Legal
values for `polsymb` are `:copol` and `:xpol`.  Again, capitalization is not
significant.  Copol is defined as the polarization with maximum amplitude
at θ = ϕ = 0.
"""
function amplitude_db end

amplitude_db(c::Cut, ipol::Integer) = 10 * log10.(abs2.(getindex.(get_evec(c), ipol)))

amplitude_db(c::Cut, polsymb::Symbol) = amplitude_db(c, lowercase(string(polsymb)))

function amplitude_db(c::Cut, polstr::String="copol")
    polstr = lowercase(strip(polstr))
    minflag = maxflag = false
    if isequal(polstr, "copol")
        maxflag = true
    elseif isequal(polstr, "xpol")
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
    phase_deg(c::Cut, polsymb::Symbol)
    phase_deg(c::Cut, polstr::String = "copol")

Return a matrix of phases in degrees for some choice of polarization in the cut.
Legal values for `ipol` are 1, 2, or 3.  Legal values for `polstr` are
"copol" (the default) and "xpol" (capitalization is not significant).  Legal
values of `polsymb` are `:copol` and `:xpol`.  Again, capitalization is not significant.
The returned value is a `Matrix{Float64}` with `length(get_theta(cut))` rows and
`length(get_phi(cut))` columns.
"""
function phase_deg end

phase_deg(c::Cut, ipol::Int) = rad2deg.(angle.(getindex.(get_evec(c), ipol)))

phase_deg(c::Cut, polsymb::Symbol) = phase_deg(c, lowercase(string(polsymb)))


function phase_deg(c::Cut, polstr::String="copol")
    polstr = lowercase(polstr)
    minflag = maxflag = false
    if isequal(polstr, "copol")
        maxflag = true
    elseif isequal(polstr, "xpol")
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
    power(cut::Cut, θmax=180)
    
Compute the total radiated power in a Cut object. 
If only a single phi value is included in the cut, then assume no azimuthal variation.
The integration in the θ direction will be computed over the limits from 0 to `min(θmax, last(get_theta(cut)))`.
`θmax` should be specified in degrees.
"""
function power(cut::Cut, θmax=180)::Float64
    # This version uses the trapezoidal rule in phi and integration of a
    # cubic spline interpolant in the theta direction.
    get_icut(cut) == 1 || error("Not a standard polar spherical cut")
    get_ncomp(cut) == 3 && error("Cut has 3 field components.  Only 2 allowed.")
    sym = get_theta(cut)[begin] < 0  # Symmetrical cut
    phifullmax = sym ? 180.0 : 360.0
    phi = get_phi(cut)
    nphi = length(phi)
    p = vec(sum(x -> real(x ⋅ x), get_evec(cut), dims=2))
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
    p .*= abs.(sin.(theta)) # sintheta weighting
    spl = CubicSpline(p, theta; assume_linear_t=true)
    pwr = phimult * first(quadgk(spl, first(theta), min(last(theta), deg2rad(θmax)); atol=1e-12))
    return pwr
end


"""
    read_cutfile(fname) -> Vector{Cut}  

Read data from a possibly multi-frequency Ticra-compatible cut file.  

Return a single `Cut` struct or a vector `Cut` structs.  Each element of the returned vector 
corresponds to a particular operating frequency partition in the file.  If there 
is only a single partition in the file, then instead of returning a 1-element
vector, the single object of type `Cut` is returned as a scalar.
"""
function read_cutfile(fname::AbstractString)
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
                        evecnext = [@SVector[0.0 + 0.0im, 0.0 + 0.0im] for _ in 1:nth]
                    else
                        evecnext = [@SVector[0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im] for _ in 1:nth]# Storage increment for reading field values
                    end
                end
                kf += 1
                header = String[textline]
                cutphi = phi:phi
                theta = range(start=ths, step=dth, length=nth)
                evecallphi = ncomp == 2 ? Array{SVector{2,ComplexF64},1}(undef, 0) : Array{SVector{3,ComplexF64},1}(undef, 0)
                evec = ncomp == 2 ? Array{SVector{2,ComplexF64},2}(undef, 0, 0) : Array{SVector{3,ComplexF64},2}(undef, 0, 0)
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
            (icomp, icut, ncomp) == (cut.icomp, cut.icut, cut.ncomp) || error("Inconsistent cuts")
            ths == first(cut.theta) || error("Inconsistent ths")
            nth == length(cut.theta) || error("Inconsistent nth")
            isapprox(dth, cut.theta[2] - cut.theta[1], atol=1e-7) || error("Inconsistent dth")

            append!(evecallphi, evecnext)
            evecthisphi = @view evecallphi[end-nth+1:end]

            # Read in the field data for this phi cut
            for i = 1:nth
                str = readline(fid)
                if ncomp == 2
                    (t1, t2, t3, t4) = (parse(Float64, s) for s in split(str))
                    evecthisphi[i] = @SVector[complex(t1, t2), complex(t3, t4)]
                else
                    (t1, t2, t3, t4, t5, t6) = (parse(Float64, s) for s in split(str))
                    evecthisphi[i] = @SVector[complex(t1, t2), complex(t3, t4), complex(t5, t6)]
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
    end # do closure
    if length(cuts) > 1
        return cuts
    else
        return only(cuts)
    end
end # function

"""
    write_cutfile(fname::AbstractString, cut::Cut, title::AbstractString="Cut file created by write_cutfile")

    write_cutfile(fname::AbstractString, cuts::AbstractVector{Cut}, title::AbstractString="Cut file created by write_cutfile")

Write `Cut` cut data to a Ticra-compatible cut file.
"""
function write_cutfile(
    fname::AbstractString,
    cut::Cut{T1, T2, N},
    title::String="Cut file created by write_cutfile"
) where {T1<:AbstractRange, T2<:AbstractRange, N}
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

function write_cutfile(
    fname::AbstractString,
    cuts::AbstractVector{Cut{T1, T2, N}},
    title::String="Cut file created by write_cutfile",
) where {T1, T2, N}
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
    (x,y,z0,z90) = phscen(cut, fghz=11.80285268; min_dropoff=-10)
Estimate phase center for a cut file or `Cut` object using NSI least-squares algorithm.

The four outputs are estimates of the location of the phase center relative to
the origin used in recording the data in the cut.  If `fghz` is passed in as 
an argument, the values will be expressed in units of inches.  Otherwise the 
lengths will be normalized to wavelengths. Note that `z0` and `z90` are the 
phi = 0ᵒ and phi = 90ᵒ plane estimates of the phase center location.  
In determining the phase center locations, only field values with magnitudes
in dB greater than `min_dropoff` relative to the peak field are considered.
"""
function phscen end

function phscen(cut::Cut, fghz=11.802852677165355; min_dropoff=-10.0)
    cut = asym2sym(cut)  # Ensure cut is symmetric
    # Determine which pol slot has main polarization:
    p1max, p2max = (maximum(x -> abs2(x[i]), get_evec(cut)) for i in 1:2)
    p = p1max > p2max ? 1 : 2
    E = getindex.(get_evec(cut), p)

    x = deg2rad.(cut.theta)

    # Boresight electric field:
    thetamin = 0.0
    index = findfirst(iszero, get_theta(cut))
    kk = index
    Ebore = E[index, begin]
    E .*= inv(Ebore) # Normalize to boresight (including phase)
    Edb = 10 .* log10.(abs2.(E))
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
        y = unwrap(y)
        y .-= y[kk] # Normalize to theta = 0 value
        cutdb = Edb[:, index]
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

function phscen(cutfile::AbstractString, fghz=11.802852677165355; min_dropoff=-10)
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


"""
    (c, xn, sp, slh, et, pc, xpd) = eval_cut(cut, fghz, thetamax)

Determine primary pattern performance metrics for a `Cut` object or 
Ticra cut file.

## Arguments:

- `cut`:   A string containing the name of the cut file to evaluate,
  or a `Cut` object as returned by `read_cutfile`.
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
function eval_cut end

eval_cut(cutfile::AbstractString, fghz::Real, thetamax::Real) = eval_cut(read_cutfile(cutfile), fghz, thetamax)

function eval_cut(cut::Cut, fghz::Real, thetamax::Real)
    cut = sym2asym(cut) # Ensure that cut is asymmetric
    pwr_tot = power(cut)  # Total power in the cut

    pol1db = 10 .* log10.(abs2.(getindex.(get_evec(cut), 1)))
    pol2db = 10 .* log10.(abs2.(getindex.(get_evec(cut), 2)))

    if maximum(pol1db) > maximum(pol2db)
        copol = pol1db
        cxpol = pol2db
    else
        copol = pol2db
        cxpol = pol1db
    end

    clamp!(cxpol, -500.0, Inf) # Eliminate negative infinities

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

    # Zero out the power in the cone theta >= thetamax:
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
    normalize!(cut::Cut, totpower=4π)

Normalize a Cut object so it's total power is `totpower` (which defaults to `4π`).
The default value results in field magnitude squared being numerically equal to directivity.
"""
function normalize!(cut::Cut, totpower=4π)
    pwr = power(cut)
    c = sqrt(totpower / pwr)
    @inbounds for i in eachindex(cut.evec)
        cut.evec[i] *= c
    end
    return cut
end

"""
    cut2 = convert_cut(cut::Cut, icomp::Integer)

Make a copy of `cut` and convert it to a new polarization basis determined by `icomp`.  Legal values 
for `icomp` and their meanings:
*  1 => Eθ and Eϕ
*  2 => ERHCP and ELHCP
*  3 => Eh and Ev (Ludwig 3 co and cx)
"""
convert_cut(cut::Cut{Tc,N}, icomp) where {Tc,N} = convert_cut!(deepcopy(cut), icomp)


"""
    convert_cut!(cut::Cut, icomp::Integer)

Convert `cut` to a new polarization basis determined by `icomp`.  Legal values 
for `icomp` and their meanings:
*  1 => Eθ and Eϕ
*  2 => ERHCP and ELHCP
*  3 => Eh and Ev (Ludwig 3 co and cx)
"""
function convert_cut!(cut::Cut{Tct, Tcp, N}, icomp::Integer) where {Tct, Tcp, N}
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
            mat = @SMatrix [(p̂1out⋅p̂1in) (p̂1out⋅p̂2in)
                (p̂2out⋅p̂1in) (p̂2out⋅p̂2in)]
        elseif N == 3
            mat = @SMatrix [(p̂1out⋅p̂1in) (p̂1out⋅p̂2in) 0
                (p̂2out⋅p̂1in) (p̂2out⋅p̂2in) 0
                0 0 1]
        else
            error("Unknown type parameter $N")
        end
        for row in 1:length(get_theta(cut))
            evec[row, col] = mat * evec[row, col]
        end
    end
    cut.icomp = outpol
    return cut
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
    θ̂ = @SVector [1.0 + 0im, 0.0 + 0.0im]
    ϕ̂ = @SVector [0.0 + 0.0im, 1.0 + 0im]
    ĥ = θ̂ * cosϕ - ϕ̂ * sinϕ
    v̂ = θ̂ * sinϕ + ϕ̂ * cosϕ
    R̂ = iroot2 * (ĥ - im * v̂)
    L̂ = iroot2 * (ĥ + im * v̂)
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
    _find_max_pol(cut::Cut) -> (pol::Symbol, outcut::Cut)

Find the secondary polarization (either :lhcp, :rhcp, :l3h, or :l3v) that
corresponds to the boresight primary polarization of greatest magnitude.  Return
the polarization symbol, and a copy of the input `cut` that has been converted
to a polarization basis consistent with `pol` value.
"""
function _find_max_pol(cut::Cut)
    cut = deepcopy(cut) # Don't mutate original cut
    convert_cut!(cut, 1) # Change to θ,ϕ components
    evec = get_evec(cut)
    # Find the maximum norm E-field
    _, imaxnorm = findmax(_norm², evec)
    Eθϕ_maxnorm = evec[imaxnorm]
    ϕ_maxnorm = get_phi(cut)[last(Tuple(imaxnorm))]
    # Check which polarization decomposition produces largest copol:
    Eθϕ = Eθϕ_maxnorm
    (θ̂, ϕ̂), (R̂, L̂), (ĥ, v̂) = _pol_basis_vectors(ϕ_maxnorm)
    ERL = @SMatrix[R̂⋅θ̂ R̂⋅ϕ̂; L̂⋅θ̂ L̂⋅ϕ̂] * Eθϕ
    Ehv = @SMatrix[ĥ⋅θ̂ ĥ⋅ϕ̂; v̂⋅θ̂ v̂⋅ϕ̂] * Eθϕ
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
    convert_cut!(cut, icompm1 + 1)
    return (pol, cut)
end



"""
    sor_efficiency(cut; F, D, Oc, pol=:matched, dz=0.0) -> (;ηₗₒₛₛ, ηₛₚ, ηᵢₗ, ηₚₕ, ηₚₒₗ, ηₓ)

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
* `pol`: A `Symbol` having one of the values `:L3h`, `:L3v`, `:RHCP`, `:LHCP`, or :max
  (any variations in terms of capitalization are acceptable).  The first two denote Ludwig 3 
  horizontal (x) and vertical (y) polarizations, the second two denote the two senses of 
  circular polarization, and `:max` (the default) uses the polarization among the 4 previously
  listed that that has the maximum field amplitude.  Polarization efficiency and boresight polarization
  mismatch of the **secondary** pattern will be computed relative to the polarization specified by 
  this argument.  Note that for CP, this will be orthogonal to the *primary* copol.
* `dz`: This is a signed distance along the zfeed direction, measured in wavelengths. It allows for
  repositioning the feed in an attempt to locate the feed phase center at the reflector focal
  point. Suppose that the feed's phase center is located 0.42 wavelengths inside the horn aperture. The
  horn origin should ideally be positioned closer to the reflector, so a positive value `dz = 0.42` would
  be specified to indicate that the horn has been repositioned in this manner.  

## Return value:  A `NamedTuple` with the following fields:
;ηₗₒₛₛ, ηₛₚ, ηᵢₗ, ηₚₕ, ηₚₒₗ, ηₓ
* `ηₗₒₛₛ`: Feed loss efficiency (ratio of radiated power to 4π).
* `ηₛₚ`: Spillover efficiency (a real number between 0 and 1).
* `ηᵢₗ`: Illumination (amplitude taper) efficiency (a real number between 0 and 1).
* `ηₚₕ`:  Phase error efficiency (a real number between 0 and 1).
* `ηₚₒₗ`: Polarization efficiency (a real number between 0 and 1).
* `ηₓ`: Boresight polarization mismatch loss efficiency (a real number between 0 and 1)
"""
sor_efficiency(cutfile::String; kwargs...) = sor_efficiency(read_cutfile(cutfile); kwargs...)

function sor_efficiency(cut::Cut; pol::Symbol=:max, F::Number, D::Number, Oc::Number, dz::Real=0.0)
    cut = sym2asym(cut) # Ensure asymmetrical, and make copy to avoid accidental mutation
    pol = Symbol(lowercase(string(pol)))
    pol in (:lhcp, :rhcp, :l3h, :l3v, :max) || error("Illegal pol value: $pol")
    # Compute β (feed tilt angle) and θe (edge ray angle):
    xm = Oc + D / 2  # Upper edge
    zm = xm^2 / 4F
    θU = atan(xm / (F - zm))  # Upper edge angle
    xm = Oc - D / 2  # Lower edge
    zm = xm^2 / 4F
    θL = atan(xm / (F - zm))  # Lower edge angle
    β = 0.5 * (θU + θL)  # Bisector angle (radians)
    sβ, cβ = sincos(β)
    θe = 0.5 * (θU - θL)  # Edge angle (radians)

    pwr = power(cut) # Total radiated power
    ηₗₒₛₛ = clamp(pwr / 4π, 0, 1)

    θmax = θe
    θmaxdeg = rad2deg(θmax)
    θmaxdeg > last(get_theta(cut)) && error("Cut does not extend to edges of reflector!")
    pwr_cone = power(cut, θmaxdeg)
    ηₛₚ = clamp(pwr_cone / pwr, 0, 1) # Spillover efficiency

    # Normalize cut to directivity:
    cut.evec .*= inv(sqrt(ηₗₒₛₛ))

    # Get needed trig functions
    scθ = sincosd.(get_theta(cut))
    scϕ = sincosd.(get_phi(cut))

    # Adjust phase to account for translation of dz along z-feed axis:
    if !iszero(dz)
        cfactor = @. cis(2π * dz * last(scθ))
        cut.evec .*= cfactor
    end

    # Obtain new cut decomposed in components defined by pol:
    if pol == :max
        (pol, cut_copol) = _find_max_pol(cut)
    else
        cut_copol = deepcopy(cut)
        pol in (:rhcp, :lhcp) && convert_cut!(cut_copol, 2)
        pol in (:l3h, :l3v) && convert_cut!(cut_copol, 3)
    end

    # Zero out the crosspol (remember that pol refers to secondary pattern):
    if pol in (:lhcp, :l3h)
        for i in eachindex(cut_copol.evec)
            e1e2 = cut_copol.evec[i]
            cut_copol.evec[i] = SVector(first(e1e2), zero(ComplexF64))
        end
    else
        for i in eachindex(cut_copol.evec)
            e1e2 = cut_copol.evec[i]
            cut_copol.evec[i] = SVector(zero(ComplexF64), last(e1e2))
        end
    end

    # Zero out the phase variations
    cut_copol_0phase = deepcopy(cut_copol)
    for i in eachindex(cut_copol_0phase.evec)
        cut_copol_0phase.evec[i] = abs.(cut_copol_0phase.evec[i])
    end

    # Convert fields to θ/ϕ components in preparation for integration:
    foreach(x -> convert_cut!(x, 1), (cut, cut_copol, cut_copol_0phase))

    evec = get_evec(cut) # actual fields
    evec_copol = get_evec(cut_copol)
    evec_copol_0phase = get_evec(cut_copol_0phase)

    # Storage for ϕ integrals at each θ location
    itgr, itgr_copol, itgr_copol0 = (zeros(SVector{2,ComplexF64}, length(get_theta(cut))) for _ in 1:3)

    # Compute ϕ integrals:
    for (iθ, θ) in pairs(get_theta(cut))
        sθ, cθ = scθ[iθ]
        θ == 180 && ((sθ, cθ) = (0.0001, -0.9999)) # Avoid irrelevant singularity at 180°
        for (iϕ, ϕ) in pairs(get_phi(cut))
            eθ, eϕ = evec[iθ, iϕ]
            ecθ, ecϕ = evec_copol[iθ, iϕ]
            ec0θ, ec0ϕ = evec_copol_0phase[iθ, iϕ]
            sϕ, cϕ = scϕ[iϕ]
            common = sθ / (1 + cθ * cβ - sθ * cϕ * sβ)^2
            f1 = cϕ * (1 + cβ * cθ) - sβ * sθ
            f2 = -sϕ * (cβ + cθ)
            itgr[iθ] += common * @SVector [f1 * eθ + f2 * eϕ, f2 * eθ - f1 * eϕ]
            itgr_copol[iθ] += common * @SVector [f1 * ecθ + f2 * ecϕ, f2 * ecθ - f1 * ecϕ]
            itgr_copol0[iθ] += common * @SVector [f1 * ec0θ + f2 * ec0ϕ, f2 * ec0θ - f1 * ec0ϕ]
        end
    end
    phi = get_phi(cut)
    dphi = deg2rad(phi[2] - phi[1])
    itgr .*= dphi
    itgr_copol .*= dphi
    itgr_copol0 .*= dphi

    # Perform integration over theta:
    thetarad = deg2rad.(get_theta(cut))
    spl = CubicSpline(itgr, thetarad; assume_linear_t=true)
    Ivec, errest1 = quadgk(spl, 0, θmax; rtol=1e-10)
    splc = CubicSpline(itgr_copol, thetarad; assume_linear_t=true)
    Icvec, errest2 = quadgk(splc, 0, θmax; rtol=1e-10)
    splc0 = CubicSpline(itgr_copol0, thetarad; assume_linear_t=true)
    Ic0vec, errest3 = quadgk(splc0, 0, θmax; rtol=1e-10)

    # Compute efficiency vectors:
    factor = 2 / π * F / D
    Ivec *= factor
    Icvec *= factor
    Ic0vec *= factor

    ηI = clamp(_norm²(Ivec), 0, 1)
    ηₚₒₗ = clamp(ηI / _norm²(Icvec), 0, 1)  # Polarization efficiency
    ηₚₕ = clamp(_norm²(Icvec) / _norm²(Ic0vec), 0, 1)  # Phase error efficiency
    ηᵢₗ = clamp(ηI / (ηₗₒₛₛ * ηₛₚ * ηₚₒₗ * ηₚₕ), 0, 1)

    # Compute polarization unitary vector:
    if pol == :l3h
        uhat = @SVector [one(ComplexF64), zero(ComplexF64)]
    elseif pol == :l3v
        uhat = @SVector [zero(ComplexF64), one(ComplexF64)]
    elseif pol == :lhcp
        uhat = @SVector [iroot2, im * iroot2]
    elseif pol == :rhcp
        uhat = @SVector [iroot2, -im * iroot2]
    else
        error("Illegal pol value: $pol")
    end

    # Boresight polarization mismatch loss:
    ηₓ = clamp(abs2(uhat ⋅ Ivec) / _norm²(Ivec), 0, 1)

    return (; ηₗₒₛₛ, ηₛₚ, ηᵢₗ, ηₚₕ, ηₚₒₗ, ηₓ)
end

"""
    sym2asym(cut::Cut) -> cut2

Convert a symmetrical (in θ) `Cut` to a new asymmetrical `Cut`.  An asymmetrical cut
begins at θ = 0, while a symmetrical cut covers equal extents of negative and 
positive angles.  If the input `cut` is already asymmetrical, then return a copy of 
this input as the output.
"""
function sym2asym(cut::Cut)
    # Check validity of input cut
    get_icut(cut) == 1 || error("Not a standard polar spherical cut")
    phi = get_phi(cut)
    theta = get_theta(cut)
    iszero(first(theta)) && return deepcopy(cut)
    first(theta) == -last(theta) || error("Not a symmetrical cut")
    0 in theta || error("cut does not include θ=0")
    dphi = phi[2] - phi[1]
    last(phi) + dphi ≈ first(phi) + 180 || error("cut phi values not properly distributed")
    np1 = length(phi)
    nt1 = length(theta)

    cut = deepcopy(cut) # Avoid mutating the input argument
    dtheta = theta[2] - theta[1]
    icomp = get_icomp(cut)
    icut = get_icut(cut)
    ncomp = get_ncomp(cut)
    convert_cut!(cut, 1) # Convert to θ/ϕ components
    it0 = findfirst(iszero, theta)

    # Set up the new cut2, also in θ/ϕ components
    phi2 = first(phi):dphi:(360-dphi)
    np2 = length(phi2)
    theta2 = 0:dtheta:last(theta)
    nt2 = length(theta2)
    evec2 = Array{SVector{2,ComplexF64}}(undef, nt2, np2)
    text2 = ["phi = $p" for p in phi2]
    for ip1 in 1:np1
        ip2 = ip1 + np1
        for (it2, it) in enumerate(it0:nt1)
            evec2[it2, ip1] = cut.evec[it, ip1]
            evec2[it2, ip2] = -cut.evec[it0-it2+1, ip1]
        end
    end
    cut2 = Cut(ncomp, icut, 1, text2, theta2, phi2, evec2)
    convert_cut!(cut2, icomp) # Restore polarization basis
    return cut2
end

"""
    asym2sym(cut::Cut) -> cut2

Convert an asymmetrical (in θ) `Cut` to a new symmetrical `Cut`.  An asymmetrical cut
begins at θ = 0, while a symmetrical cut covers equal extents of negative and 
positive angles.  If the input `cut` is already symmetrical, then return a copy of 
this input as the output.
"""
function asym2sym(cut::Cut)
    # Check validity of input cut
    get_icut(cut) == 1 || error("Not a standard polar spherical cut")
    phi = get_phi(cut)
    theta = get_theta(cut)
    (first(theta) == -last(theta)) && return deepcopy(cut)
    iszero(first(theta)) || error("cut does begin at θ=0")
    iseven(length(phi)) || error("Number of ϕ cuts is not even")
    for ϕ in phi
        found = false
        for ϕ′ in ϕ
            if iszero(mod(ϕ - ϕ′, 360))
                found = true
                break
            end
        end
        !found && error("cut contains $ϕ but not $(ϕ±180)")
    end

    cut = deepcopy(cut) # Avoid mutating the input argument
    dtheta = theta[2] - theta[1]
    icomp = get_icomp(cut)
    icut = get_icut(cut)
    ncomp = get_ncomp(cut)
    convert_cut!(cut, 1) # Convert to θ/ϕ components
    it0 = findfirst(iszero, theta)
    nt1 = length(theta)
    np1 = length(phi)
    np1o2 = length(phi) ÷ 2

    # Set up the new cut2, also in θ/ϕ components
    phi2 = phi[1:np1÷2]
    np2 = length(phi2)
    theta2 = -last(theta):dtheta:last(theta)
    nt2 = length(theta2)
    evec2 = Array{SVector{2,ComplexF64}}(undef, nt2, np2)
    text2 = ["phi = $p" for p in phi2]
    for ip2 in 1:np2
        for i in 1:nt1
            evec2[end-i+1, ip2] = cut.evec[end-i+1, ip2]
            evec2[i, ip2] = -cut.evec[end-i+1, ip2+np2]
        end
    end
    cut2 = Cut(ncomp, icut, 1, text2, theta2, phi2, evec2)
    convert_cut!(cut2, icomp) # Restore polarization basis
    return cut2
end

"""
    eh2bor1cut(theta, fe, fh; kwargs...) -> cut::Cut

Create a [`Cut`](@ref) object for a "BOR1" horn from its E-plane and H-plane patterns.

A "BOR₁" horn is circularly symmetric and contains only TE₁ₙ and TM₁ₙ waveguide modes in its 
radiating aperture.  It's radiated far field can therefore be expressed in terms
of the E-plane and H-plane patterns it radiates when excited for linear polarization.

## Positional Input Arguments
- `theta`: A vector or range (an `AbstractVector`) of θ values (in degrees) at which the cut 
  pattern should be evaluated. The first element of `theta` must be 0, and the entries must 
  be equally spaced, as in a `range` object.
- `fe`, `fh`: The E-plane and H-plane patterns, resp.  These are either both `AbstractVector`s of the 
  same length as `theta`, or both functions which take a single input θ (in degrees) and return the 
  respective patterns evaluated at that angle.

## Keyword Arguments
- `pol`: Defines the manner in which the horn is assumed to be excited, and the polarization basis 
  selected for use in the output `Cut`.  `pol` is a `String` or `Symbol` taking one of the values
  (capitalization is not significant):
  * "l3v" or `:l3v`: (the default value) The horn is excited for "vertical" (y-directed) linear polarization 
    and the far field is expressed as Ludwig-3 components.  
  * "l3h" or `:l3h`: The horn is excited for "horizontal" (x-directed) linear polarization and the far field is
    expressed as Ludwig-3 components.
  * "rhcp" or `:rhcp`: The horn is excited for RHCP (right-hand circular polarization) and the far field is
    expressed as RHCP and LHCP components.
  * "lhcp" or `:lhcp`: The horn is excited for LHCP (left-hand circular polarization) and the far field is
    expressed as RHCP and LHCP components.
  If linear (circular) polarization is requested, then the output `Cut` object will contain eight (four)
  cuts, spaced every 45° (90°). 
- `xpd`: The crosspol level in dB < 0. Defaults to `-Inf` (negative infinity).  If finite, then
  in addition to the specified polarization, a crosspolarized contribution will be added to the cut,
  as if the horn is fed by an imperfect feed network with the specified crosspol level.
- `xpphase`: The phase (in degrees) of the crosspol contribution whose amplitude is specified by `xpd`.
"""
function eh2bor1cut(
    theta::AbstractVector,
    fe::AbstractVector, 
    fh::AbstractVector;
    pol::Union{Symbol, String} = :l3v,
    xpd::Real = -Inf,
    xpphase::Real = 0.0)

    polsymb = string(pol) |> lowercase |> Symbol
    if polsymb ∈ (:l3h, :l3v) 
        icomp = 3
        xpolsym = setdiff((:l3h, :l3v), (polsymb,)) |> only
        phi = 0:45:(360-45)
    elseif polsymb ∈ (:rhcp, :lhcp)
        icomp = 2
        xpolsym = setdiff((:rhcp, :lhcp), (polsymb,)) |> only
        phi = 0:90:(360-90)
    else 
        throw(ArgumentError("Illegal value for pol"))
    end
    xpd ≥ 0 && throw(ArgumentError("xpd must be < 0"))

    ncomp = 2 # Number of polarization components
    icut = 1 # standard polar cut
    nphi = length(phi)
    ntheta = length(theta)
    thetarange = range(first(theta), last(theta), ntheta)
    first(theta) |> iszero || throw(ArgumentError("theta does not begin at 0"))
    thetarange ≈ theta || throw(ArgumentError("theta not equally spaced"))
    evec = Matrix{SVector{2,ComplexF64}}(undef, ntheta, nphi)
    text = ["phi = " * string(p) for p in phi]
    cut = Cut(; ncomp, icut, icomp, text, theta=thetarange, phi, evec)

    for (kp, ph) in pairs(phi)
        sp, cp = sincosd(ph)
        cpsp = cp * sp
        sp², cp² = sp * sp, cp * cp
        cis2phi = cis(2 * deg2rad(ph))
        for (kt, th) in pairs(thetarange)
            if polsymb == :l3h
                cut.evec[kt, kp] = SVector(fe[kt] * cp² + fh[kt] * sp², (fe[kt] - fh[kt]) * cpsp)
            elseif polsymb == :l3v
                cut.evec[kt, kp] = SVector((fe[kt] - fh[kt]) * cpsp, fe[kt] * sp² + fh[kt] * cp²)
            elseif polsymb == :rhcp
                cut.evec[kt, kp] = SVector(0.5 * (fe[kt] + fh[kt]), 0.5 * conj(cis2phi) * (fe[kt] - fh[kt]))
            else # :lhcp
                cut.evec[kt, kp] = SVector(0.5 * cis2phi * (fe[kt] - fh[kt]), 0.5 * (fe[kt] + fh[kt]))
            end
        end
    end

    isinf(xpd) && return cut

    # Add in crosspol
    xpdpwr = 10 ^ (xpd / 10)
    a = sqrt(inv(1 + xpdpwr)) # Voltage weight for original cut
    b = a * sqrt(xpdpwr) * cis(deg2rad(xpphase)) # Voltage weight for xpol cut
    xcut = eh2bor1cut(theta, fe, fh; pol = xpolsym)
    cut.evec .= a .* cut.evec + b .* xcut.evec
    return cut
end

function eh2bor1cut(
    theta::AbstractVector,
    fe::F1, 
    fh::F2;
    kwargs...) where {F1<:Function, F2<:Function}
    fevec = [fe(t) for t in theta]
    fhvec = [fh(t) for t in theta]
    return eh2bor1cut(theta, fevec, fhvec; kwargs...)
end


"Plot recipe for Cut"
@recipe function f(cut::Cut;
    phi=get_phi(cut),
    theta=get_theta(cut),
    pol=:both, # or :copol or :xpol or 1 or 2, or "both" or "copol" or "xpol"
    quantity=:db, # or :power or :linear or :phase
    normalization=NaN # A number or :peak
)

    all(x -> x in get_theta(cut), extrema(theta)) || error("some requested theta are outside those of cut")
    quantity = Symbol(lowercase(string(quantity)))
    pol isa Symbol || pol isa AbstractString && (pol = Symbol(lowercase(string(pol))))
    if normalization isa Symbol
        normalization == :peak || error("Illegal value for Normalization")
    elseif normalization isa Number
        if isnan(normalization)
            # Set appropriate default depending on quantity to be plotted
            if quantity in (:db, :phase)
                normalization = 0.0
            elseif quantity in (:linear, :power)
                normalization = 1.0
            end
        end
    else
        error("Illegal type for normalization")
    end

    # set a default value for an attribute with `-->`.  Force it with `:=`.
    xguide --> "θ (deg)"
    yguide --> (quantity == :db ? "Amplitude (dB)" :
                quantity == :linear ? "Field Amplitude" :
                quantity == :power ? "Power Amplitude" :
                quantity == :phase ? "Phase (deg)" :
                error("Illegal value for quantity: $quantity"))
    icomp = get_icomp(cut)
    pol_labels = (("E_θ", "E_ϕ"), ("E_R", "E_L"), ("E_h", "E_v"))[icomp]
    evec = get_evec(cut)
    _, imaxnorm = findmax(_norm², evec)
    Evecmax = evec[imaxnorm]
    icopol = abs2(Evecmax[1]) > abs2(Evecmax[2]) ? 1 : 2
    if isequal(normalization, :peak)
        isequal(quantity, :phase) && error(":peak normalization may not be requested for :phase plot")
        Esqmax = maximum(abs2, Evecmax)
        Emax = sqrt(Esqmax)
        normdb = 10 * log10(Esqmax)
    else
        Emax = Esqmax = normdb = normalization
    end
    # Select polarizations to plot
    if isequal(pol, :both)
        ipols = (icopol, 3 - icopol)
    elseif isequal(pol, :copol)
        ipols = (icopol,)
    elseif isequal(pol, :xpol)
        ipols = (3 - icopol,)
    elseif isequal(pol, 1)
        ipols = (1,)
    elseif isequal(pol, 2)
        ipols = (2,)
    else
        error("illegal value for pol: $(pol)")
    end

    # Add a series for each phi cut and each selected polarization
    for ϕ in phi
        iϕ = findfirst(≈(ϕ), get_phi(cut))
        isnothing(iϕ) && error("ϕ = $(ϕ)° is not present in the cut to be plotted")
        spl = CubicSpline((@view evec[:, iϕ]), get_theta(cut))
        efield = spl(theta)
        for (ii, ipol) in enumerate(ipols)
            if quantity == :db
                y = [10 * log10(abs2(e[ipol])) - normdb for e in efield]
            elseif quantity == :power
                y = [abs2(e[ipol]) / Esqmax for e in efield]
            elseif quantity == :linear
                y = [abs(e[ipol]) / Emax for e in efield]
            else
                y = [rad2deg(angle(e[ipol])) - normalization for e in efield]
            end
            @series begin
                linestyle --> (:solid, :dash)[ii]
                label --> "$(pol_labels[ipol]), ϕ = $ϕ"
                theta, y
            end
        end
    end
end
