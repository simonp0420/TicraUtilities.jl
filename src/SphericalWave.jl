using OffsetArrays: OffsetArray, OffsetVector
using Printf: @printf
using AssociatedLegendrePolynomials: Plm!, LegendreNormCoeff, LegendreOrthoNorm, legendre, legendre!
using ForwardDiff: Dual, value, partials
derivative(x::Dual) = partials(x,1)
using FFTW: ifft!, fftshift
using OffsetArrays: OffsetArray, Origin
using LinearAlgebra: norm
using QuadGK: quadgk, kronrod
using FastGaussQuadrature: gausslegendre
using Interpolations: cubic_spline_interpolation
using Dates: now


"""
    SWEQPartition
Struct containing all the data read from a Ticra-compatible Q-type SWE file for one frequency.
The `qsmns` field contains the Q coefficients. It is indexed as `qsmns[s,m,n]` where `s ∈ {1,2}`, `m ∈ -mmax:mmax`, 
`n ∈ {1, ..., mmax}`.  But not all entries are used.  The only possible nonzero entries are where 
`n ∈ {max(1, abs(m)), ..., nmax}`.  Note that if the coefficients stored in the SWE files are `Q'`, and
the cofficients stored in a `SWEQPartition` instance are `Q`, then `Q = sqrt(8π) * conj(Q')`, as discussed
in the File Formats section of the Ticra official documentation.

### Fields
* `prgtag::String`: Program tag and time stamp, from the program that created the file.
* `idstrg::String`: Identification text.
* `nthe::Int`: Number of θ-samples over 360°. Must be even and ≥ 4. 
* `nphi::Int`: Number of ϕ-samples over 360°. Must be ≥ 3. 
* `nmax::Int`: Maximum value for polar index `n` in the SW expansion. `1 ≤ nmax ≤ nthe÷2`.
* `mmax::Int`: Maximum value for azimuthal index `|m|` in the SW expansion. 0 ≤ mmax ≤ min((nphi-1)÷2, nmax).
* `t4::String`: Dummy string.
* `t5::String`: Dummy string.  If parsed, contains 5 real numbers.
* `t6::String`: Dummy string.  If parsed, contains 5 real numbers.
* `t7::String`: Dummy string.
* `t8::String`: Dummy string.
* `qsmns::OffsetArray{ComplexF64, 3}`: Contains the Q coefficients, assuming exp(jωt) time variation
  and using the Ticra normalization convention.  The array has axes `(1:2, -mmax:mmax, 1:nmax)`, meaning
  that some of the entries (those with `n < abs(m)`) are always zero.
* `powerms::OffsetVector{Float64, Vector{Float64}}`: The array has axes `(0:mmax)` and the `m`th element
  contains the total power (one-half the sum of the magnitude squared) of all modes with `±m` as the m index.
"""
@kwdef struct SWEQPartition
    prgtag::String
    idstrg::String
    nthe::Int
    nphi::Int
    nmax::Int
    mmax::Int
    t4::String = "dummy t4"
    t5::String = "1.0 2.0 3.0 4.0 5.0"
    t6::String = "1.0 2.0 3.0 4.0 5.0"
    t7::String = "dummy t7"
    t8::String = "dummy t8"
    qsmns::OffsetArray{ComplexF64, 3, Array{ComplexF64, 3}} = 
        OffsetArray(zeros(ComplexF64, (2, 2mmax+1, nmax)), 1:2, -mmax:mmax, 1:nmax)
    powerms::OffsetVector{Float64, Vector{Float64}}  = OffsetArray(zeros(mmax+1), 0:mmax)
end

import Base.show
function show(io::IO, swe::SWEQPartition)
    println(io, "SWEQPartition")
    println(io, "  prgtag   $(swe.prgtag)")
    println(io, "  idstrg   $(swe.idstrg)")
    println(io, "  nthe     $(swe.nthe)")
    println(io, "  nphi     $(swe.nphi)")
    println(io, "  nmax     $(swe.nmax)")
    println(io, "  mmax     $(swe.mmax)")
    println(io, string("  qsmns    OffsetArray{ComplexF64}(",
    first(axes(swe.qsmns,1)), ":", last(axes(swe.qsmns,1)), ",",
    first(axes(swe.qsmns,2)), ":", last(axes(swe.qsmns,2)), ",",
    first(axes(swe.qsmns,3)), ":", last(axes(swe.qsmns,3)), ")"))
    println(io, string("  powerms  OffsetArray{Float64}(",
    first(axes(swe.powerms,1)), ":", last(axes(swe.powerms,1)), ")"))
    return nothing
end




"""
    read_sphfile(fname) -> Vector{SWEQPartition}

Read the SWE coefficients from a Q-type spherical wave expansion file. 

In the process of reading the data, the coefficients in the file (Q′) are conjugated and then
multiplied by the factor √(8π) to achieve Ticra-standard normalization.  
Each element of the returned vector corresponds to a particular operating frequency partition
in the file.  If there is only a single partition in the file, then instead of returning a 1-element
vector, the single element of type `SWEQPartition` is returned as a scalar.
"""
function read_sphfile(fname::AbstractString)
    normfactor = sqrt(8π)
    open(fname, "r") do io
        sps = SWEQPartition[]
        while !eof(io)
            # Read header info of next partition:
            prgtag, idstrg = (readline(io) for _ in 1:2)
            nthe, nphi, nmax, mmax = parse.(Int, split(readline(io)))
            t4, t5, t6, t7, t8 = (readline(io) for _ in 4:8)

            qsmns = OffsetArray(zeros(ComplexF64, 2, 2mmax+1, nmax), 1:2, -mmax:mmax, 1:nmax)
            powerms = OffsetArray(zeros(mmax+1), 0:mmax)
            for absm in 0:mmax
                # Read Q coefficients for current value of absm
                words = split(readline(io))
                absm2 = parse(Int, words[1])
                absm2 == absm || error("absm ≠ m")
                miszero = iszero(absm)
                powerm = parse(Float64, words[2])
                powerms[absm] = powerm * 8π
                nmin = max(1, absm)
                for n in nmin:nmax
                    q1r, q1i, q2r, q2i = parse.(Float64, split(readline(io)))
                    qsmns[1, -absm, n] = normfactor * complex(q1r, q1i) |> conj # s = 1
                    qsmns[2, -absm, n] = normfactor * complex(q2r, q2i) |> conj # s = 2
                    if !miszero
                        q1r, q1i, q2r, q2i = parse.(Float64, split(readline(io)))
                        qsmns[1, absm, n] = normfactor * complex(q1r, q1i) |> conj # s = 1
                        qsmns[2, absm, n] = normfactor * complex(q2r, q2i) |> conj # s = 2
                    end
                end
            end
            push!(sps, SWEQPartition(; prgtag, idstrg, nthe, nphi, nmax, mmax, t4, t5, t6, t7, t8, qsmns, powerms))
        end
        if length(sps) > 1
            return sps
        else
            return sps[1]
        end
    end
end



"""
    write_sphfile(fname, qs::Vector{SWEQPartition})
    write_sphfile(fname, qs::SWEQPartition)

Write SWE coefficients to a Q-type spherical wave expansion file.

In the process of writing the data, the coefficients in the file (Q) are conjugated and then
multiplied by the factor 1/sqrt(8π) to become Q′ and achieve consistency with Ticra-standard normalization.  
"""
function write_sphfile(fname::AbstractString, sps::Vector{SWEQPartition})
    ipwrnormfactor = inv(8π)
    inormfactor = sqrt(ipwrnormfactor)
    open(fname, "w") do io
        for sp in sps # Loop over partitions
            (; prgtag, idstrg, nthe, nphi, nmax, mmax) = sp
            (; t4, t5, t6, t7, t8, qsmns, powerms) = sp
        
            # Write header info of next partition:
            println(io, prgtag)
            println(io, idstrg)
            @printf(io, "%6i %5i %5i %5i\n", nthe, nphi, nmax, mmax)
            foreach(t -> println(io, t), (t4, t5, t6, t7, t8))

            for absm in 0:mmax
                miszero = iszero(absm)
                nmin = max(1, absm)
                # Write Q coefficients for current value of absm
                @printf(io, "%6i %23.17G\n", absm, ipwrnormfactor*powerms[absm])
                for n in nmin:nmax
                    q1r, q1i = inormfactor * qsmns[1, -absm, n] |> conj |> reim # s = 1
                    q2r, q2i = inormfactor * qsmns[2, -absm, n] |> conj |> reim # s = 2
                    @printf(io, " %23.16E %23.16E %23.16E %23.16E\n", q1r, q1i, q2r, q2i)
                    if !miszero
                        q1r, q1i = inormfactor * qsmns[1, absm, n] |> conj |> reim # s = 1
                        q2r, q2i = inormfactor * qsmns[2, absm, n] |> conj |> reim # s = 2
                        @printf(io, " %23.16E %23.16E %23.16E %23.16E\n", q1r, q1i, q2r, q2i)
                    end
                end
            end
        end
    end
    return
end

write_sphfile(fname::AbstractString, qs::SWEQPartition) = write_sphfile(fname, [qs])



# Functions for associated Legendre functions, both normalized and unnormalized.
# See documentation of AssociatedLegendrePolynomials.jl for calling details.
const NMMAX = 500
const Qcoef = LegendreNormCoeff{LegendreOrthoNorm,Float64}(NMMAX)

"""
    P̄nm(n, m, x)

Normalized associated Legendre function following Ticra definition.
See documentation for `AssociatedLegendrePolynomials.Plm` for calling conventions.
"""
@inline function P̄nm(n::Int, m::Int, x::Number)
    m < 0 && throw(ArgumentError("m must be nonnegative"))
    p = legendre(Qcoef, n, m, x)
    # Remove Condon–Shortley phase factor
    return isodd(m) ? -p : p
end

@inline function P̄nm(n::UnitRange, m::UnitRange, x::Number)
    any(<(0),m) && throw(ArgumentError("all m values must be nonnegative"))
    p = legendre(Qcoef, n, m, x)
    # Remove Condon–Shortley phase factor
    for mp1 in axes(p,2) 
        isodd(mp1) && continue # 1-based indexing!
        p[:,mp1] .*= -1
    end
    return p
end

#P̄nm(n, m, x) = legendre(Qcoef, n, m, x)

"""
    P̄nm!(Λ, n, m, x)

Normalized associated Legendre function following Ticra definition.
See documentation for `AssociatedLegendrePolynomials.Plm!` for calling conventions.
"""
function P̄nm!(Λ::Matrix, n::Int, m::Int, x::Number)
    legendre!(Qcoef, Λ, n, m, x)
    # Remove Condon–Shortley phase factor
    for mp1 in axes(Λ, 2)
        isodd(mp1) && continue # 1-based indexing!
        Λ[:,mp1] .*= -1
    end
    return nothing
end

#P̄nm!(Λ, n, m, x) = legendre!(Qcoef, Λ, n, m, x)
#Pnm(n, m, x) = Plm(n, m, x)
#Pnm!(Λ, n, m, x) = Plm!(Λ, n, m, x)



"""
    cut2sph_adaptive(cut::TicraCut; keywords...) -> s::SWEQPartition

Convert a `TicraCut` object to a `SWEQPartition` using a very accurate
but slower adaptive Gauss-Kronrod quadrature for the θ integrals.

## Keyword Arguments (and their default values)
* `pwrtol=1e-10`: The power tolerance.  Spherical modes are included until the excluded
  modes' power is less than `pwrtol` times the total modal power. A zero or negative value
  precludes removal of any modes.
* `mmax=$(NMMAX)`: An upper limit for the `m` (azimuthal) mode index to be included.
  The actual limit will be set to `min(Nϕ÷2, mmax)` for odd `Nϕ`, and `min(Nϕ÷2-1, mmax)`
  for even `Nϕ`, where `Nϕ` is the number of ϕ = constant polar cuts in the cut object.
* `nmax=$(NMMAX)`: An upper limit for the `n` (polar) mode index to be included.
  The actual limit will be the lesser of `nmax` and `Nθ-1` where `Nθ` is the number of 
  θ values included in each ϕ = constant polar cut.
"""
function cut2sph_adaptive(cut::TicraCut; mmax=NMMAX, nmax=NMMAX, pwrtol=1.e-10, gkorder=25)
    get_ncomp(cut) == 2 || error("Must have only 2 polarization components")
    cutpwr = power(cut)
    cutθϕ = deepcopy(cut); convert_cut!(cutθϕ, 1) # Convert to Eθ and Eϕ
    θs = get_theta(cutθϕ); Nθ = length(θs); Δθ = deg2rad(θs[2] - θs[1])
    ϕs = get_phi(cutθϕ);   Nϕ = length(ϕs); Δϕ = deg2rad(ϕs[2] - ϕs[1])
    eθ = first.(get_evec(cutθϕ))
    eϕ = last.(get_evec(cutθϕ))

    # Perform ϕ integration:
    ifft!(eθ, 2);  eθ .*= 2Δϕ * Nϕ # undo scaling and include leading "2"
    ifft!(eϕ, 2);  eϕ .*= 2Δϕ * Nϕ # undo scaling and include leading "2"
    eθ = fftshift(eθ, 2)
    eϕ = fftshift(eϕ, 2)
    Nϕo2 = Nϕ ÷ 2
    nmax = min(Nθ - 1, nmax)
    nrange = 1:nmax
    if isodd(Nϕ)
        mabsmax = min(Nϕo2, mmax)
        Eθmat = OffsetArray(eθ, 1:Nθ, -Nϕo2:Nϕo2)
        Eϕmat = OffsetArray(eϕ, 1:Nθ, -Nϕo2:Nϕo2)
    else
        mabsmax = min(Nϕo2 - 1, mmax)
        Eθmat = OffsetArray(eθ, 1:Nθ, -Nϕo2:(Nϕo2-1))
        Eϕmat = OffsetArray(eϕ, 1:Nθ, -Nϕo2:(Nϕo2-1))
    end
    mrange = -mabsmax:mabsmax
    qsmns = OffsetArray(zeros(ComplexF64, (2, 2mabsmax+1, nmax)), 1:2, -mabsmax:mabsmax, 1:nmax)
    cfactor = sqrt(π)/360 # Includes 1/sqrt(2) needed for f1 and f2, and π/180 for dθ in degrees

    for m in mrange
        Eθfun = cubic_spline_interpolation(θs, view(Eθmat, :, m))
        Eϕfun = cubic_spline_interpolation(θs, view(Eϕmat, :, m))
        mabs = abs(m)
        mfactor = 1
        (m > 0 && isodd(m)) && (mfactor = -1)
        negjⁿ⁺¹ = negjⁿ = negj = complex(0,-1)
        for n in nrange
            negjⁿ = negjⁿ⁺¹
            negjⁿ⁺¹ = negjⁿ * negj
            n < mabs && continue
            nfactor = inv(sqrt(n*(n+1)))
            cmn = cfactor * mfactor * nfactor
            f1factor = negjⁿ⁺¹ * cmn 
            f2factor = negjⁿ * cmn 
            int, errest = quadgk(0.0, 180.0; atol=1e-10, order=gkorder) do θ
                (θ < θs[1] || θ > θs[end]) && return @SVector[0.0,0.0]
                sin²θ = sind(θ)^2
                result = P̄nm(n, mabs, Dual(cosd(θ), 1.0))
                pnm = value(result)
                pnm′ = derivative(result)
                f1θc = f1factor * (-negj) * (m * pnm)
                f1ϕc = f1factor * (sin²θ * pnm′)
                f2θc = f2factor * ((-sin²θ) * pnm′)
                f2ϕc = f2factor * (-negj) * (m * pnm)
                Eθ = Eθfun(θ)
                Eϕ = Eϕfun(θ)
                integrand = @SVector[f1θc * Eθ + f1ϕc * Eϕ, f2θc * Eθ + f2ϕc * Eϕ]
                return integrand
            end
            qsmns[1,m,n], qsmns[2,m,n] = int
        end
    end

    # The Q coefficients have now been calculated.
    (qpwr, powerm, qsmn) = _filter_qmodes_by_power(qsmns, pwrtol)

    # Create output SWEQPartition
    date, clock = split(string(now()), 'T')
    funcname = nameof(var"#self#")
    prgtag = string(funcname, " ", date, " ", clock)
    idstrg = "Spherical Wave Q-Coefficients"
    mmax = last(axes(qsmn, 2))
    nmax = last(axes(qsmn, 3))
    #nthe = (Nθ - 1) * 2
    nthe = max(2*nmax, 6)
    #nphi = Nϕ
    nphi = max(2(mmax + 1), 4)
    t4 = t5 = t6 = t7 = "Dummy Text"
    return SWEQPartition(; prgtag, idstrg, nthe, nphi, nmax, mmax, t4, t5, t6, t7, qsmns=qsmn, powerms=powerm)
end

const GKORDER = 500  # Default value for nonadaptive GK quadrature
const gkdefaults = kronrod(GKORDER, 0, 180)


"""
    cut2sph(cut::TicraCut; keywords...) -> s::SWEQPartition

Convert a `TicraCut` object to a `SWEQPartition` using a fixed-order Gauss-Kronrod
quadrature scheme for the θ integrals.

## Keyword Arguments (and their default values)
* `mmax=$(NMMAX)`: An upper limit for the `m` (azimuthal) mode index to be included.
  The actual limit will be set to `min(Nϕ÷2, mmax)` for odd `Nϕ`, and `min(Nϕ÷2-1, mmax)`
  for even `Nϕ`, where `Nϕ` is the number of ϕ = constant polar cuts in the cut object.
* `nmax=$(NMMAX)`: An upper limit for the `n` (polar) mode index to be included.
  The actual limit will be the lesser of `nmax` and `Nθ-1` where `Nθ` is the number of 
  θ values included in each ϕ = constant polar cut.
* `pwrtol=1e-10`: The power tolerance.  Spherical modes are included until the excluded
  modes' power is less than `pwrtol` times the total modal power.  A zero or negative value
  precludes removal of any modes.
* `gkorder=$(GKORDER)`: The order of Gauss-Kronrod quadrature to use for the θ integrals.  a
  warning will be printed if the order is not large enough for accurate calculation of
  all requested modal Q coefficients.  Note: Decreasing `gkorder` from its default value
  will not speed up the calculations--it will slow them down.  This is because the GK
  nodes and weights are precomputed for the default value, but must be computed on the 
  fly for any other value.
"""
function cut2sph(cut::TicraCut; pwrtol=1e-10, mmax=NMMAX, nmax=NMMAX, gkorder=GKORDER)
    get_ncomp(cut) == 2 || error("Must have only 2 polarization components")
    cutpwr = power(cut)
    cutθϕ = deepcopy(cut); convert_cut!(cutθϕ, 1) # Convert to Eθ and Eϕ
    θs = get_theta(cutθϕ); Nθ = length(θs); Δθ = deg2rad(θs[2] - θs[1])
    ϕs = get_phi(cutθϕ);   Nϕ = length(ϕs); Δϕ = deg2rad(ϕs[2] - ϕs[1])
    eθ = first.(get_evec(cutθϕ))
    eϕ = last.(get_evec(cutθϕ))
    # Perform ϕ integration:
    ifft!(eθ, 2);  eθ .*= 2Δϕ * Nϕ # undo scaling
    ifft!(eϕ, 2);  eϕ .*= 2Δϕ * Nϕ # undo scaling
    eθ = fftshift(eθ, 2)
    eϕ = fftshift(eϕ, 2)
    Nϕo2 = Nϕ ÷ 2
    nmax = min(Nθ - 1, nmax)
    if nmax > NMMAX
        @warn "Reducing nmax to NMMAX = $(NMMAX)."
        nmax = NMMAX
    end
    if gkorder == GKORDER
        θnodes, wts, gwts = gkdefaults
    else
        θnodes, wts, gwts = kronrod(gkorder, 0, 180)
    end
    lnodes = length(θnodes)

    erange = isodd(Nϕ) ? (-Nϕo2:Nϕo2) : (-Nϕo2:(Nϕo2-1))
    mabsmax = min(last(erange), mmax)
    mabsmax > NMMAX && error("mabsmax is $(mabsmax) but may not exceed NMMAX = $(NMMAX).")
    eθoffset = OffsetArray(eθ, 1:Nθ, erange)
    eϕoffset = OffsetArray(eϕ, 1:Nθ, erange)
    Eθ, Eϕ = (OffsetArray(zeros(ComplexF64, lnodes, Nϕ), 1:lnodes, erange) for _ in 1:2)

    mrange = -mabsmax:mabsmax
    for m in mrange
        Eθfun = cubic_spline_interpolation(θs, view(eθoffset, :, m))
        Eϕfun = cubic_spline_interpolation(θs, view(eϕoffset, :, m))
        @inbounds for i in 1:lnodes
            Eθ[i,m] = Eθfun(θnodes[i])
            Eϕ[i,m] = Eϕfun(θnodes[i])
        end
    end
    result_parent = Matrix{typeof(Dual(1.0,1.0))}(undef, nmax+1, mabsmax+1) # storage for legendre functions
    result = Origin(0)(result_parent)
    qsmns = OffsetArray(zeros(ComplexF64, (2, 2mabsmax+1, nmax)), 1:2, -mabsmax:mabsmax, 1:nmax)
    qsmns_low = OffsetArray(zeros(ComplexF64, (2, 2mabsmax+1, nmax)), 1:2, -mabsmax:mabsmax, 1:nmax)
    cfactor = sqrt(π)/360 # Includes 1/sqrt(2) needed for f1 and f2, and π/180 for dθ in degrees
    @inbounds for i in 1:lnodes
        θ = θnodes[i]
        (θ < 1.e-5 || θ > 180 - 1.e-5) && continue # No contribution from endpoints
        sinθ, cosθ = sincosd(θ)
        sin²θ = sinθ * sinθ
        P̄nm!(result_parent, nmax, mabsmax, Dual(cosθ, 1.0))
        for m in mrange
            mabs = abs(m)
            mfactor =  (m > 0 && isodd(m)) ? -1 : 1
            Eθᵢₘ = Eθ[i, m]
            Eϕᵢₘ = Eϕ[i, m]
            negj = negjⁿ = negjⁿ⁺¹ = complex(0,-1)
            for n in 1:nmax
                negjⁿ = negjⁿ⁺¹
                negjⁿ⁺¹ = negjⁿ * negj
                n < mabs && continue
                nfactor = inv(sqrt(n*(n+1)))
                cmn = cfactor * mfactor * nfactor
                pnm = value(result[n,mabs])
                mpnm = m*pnm
                pnm′ = derivative(result[n,mabs])
                f1factor = negjⁿ⁺¹ * cmn 
                jf1factor = complex(-imag(f1factor), real(f1factor))
                f1θconj = jf1factor * mpnm
                sin²θpnm′ = sin²θ * pnm′
                f1ϕconj = f1factor * sin²θpnm′
                f2factor = negjⁿ * cmn 
                jf2factor = complex(-imag(f2factor), real(f2factor))
                f2θconj = f2factor * (-sin²θpnm′)
                f2ϕconj = jf2factor * mpnm
                q1 = f1θconj * Eθᵢₘ + f1ϕconj * Eϕᵢₘ
                q2 = f2θconj * Eθᵢₘ + f2ϕconj * Eϕᵢₘ
                # Compute full-order Gaussian-Kronrod quadrature:
                qsmns[1,m,n] += wts[i] * q1
                qsmns[2,m,n] += wts[i] * q2
                if iseven(i)
                    # Compute lower-order Gaussian quadrature:
                    ii = i÷2
                    qsmns_low[1,m,n] += gwts[ii] * q1
                    qsmns_low[2,m,n] += gwts[ii] * q2
                end
            end
        end
    end

    # The Q coefficients have now been calculated. Check GK convergence:
    gkabserrtol = 1e-8
    gkabserr = maximum(abs(q - qlow) for (q,qlow) in zip(qsmns,qsmns_low))
    if gkabserr > gkabserrtol
        @warn """Quadrature Nonconvergence
                 Max Gauss-Kronrod error = $(round(gkabserr; sigdigits=4)) > tolerance = $(gkabserrtol). 
                 Recommendation: Increase gkorder or reduce nmax
        """
    end

    # Eliminate modes with negligible power
    (qpwr, powerm, qsmn) = _filter_qmodes_by_power(qsmns, pwrtol)

    # Create output SWEQPartition
    date, clock = split(string(now()), 'T')
    funcname = nameof(var"#self#")
    prgtag = string(funcname, " ", date, " ", clock)
    idstrg = "Spherical Wave Q-Coefficients"
    mmax = last(axes(qsmn, 2))
    nmax = last(axes(qsmn, 3))
    #nthe = (Nθ - 1) * 2
    nthe = max(2*nmax, 6)
    #nphi = Nϕ
    nphi = max(2(mmax + 1), 4)
    t4 = t5 = t6 = t7 = "Dummy Text"
    return SWEQPartition(; prgtag, idstrg, nthe, nphi, nmax, mmax, t4, t5, t6, t7, qsmns=qsmn, powerms=powerm)
end


"""
    cut2sph_gauss(cut::TicraCut; keywords...) -> s::SWEQPartition

Convert a `TicraCut` object to a `SWEQPartition` using a fixed-order Gauss-Kronrod
quadrature scheme for the θ integrals.

## Keyword Arguments (and their default values)
* `mmax=$(NMMAX)`: An upper limit for the `m` (azimuthal) mode index to be included.
  The actual limit will be set to `min(Nϕ÷2, mmax)` for odd `Nϕ`, and `min(Nϕ÷2-1, mmax)`
  for even `Nϕ`, where `Nϕ` is the number of ϕ = constant polar cuts in the cut object.
* `nmax=$(NMMAX)`: An upper limit for the `n` (polar) mode index to be included.
  The actual limit will be the lesser of `nmax` and `Nθ-1` where `Nθ` is the number of 
  θ values included in each ϕ = constant polar cut.
* `pwrtol=1e-10`: The power tolerance.  Spherical modes are included until the excluded
  modes' power is less than `pwrtol` times the total modal power.  A zero or negative value
  precludes removal of any modes.
* `gaussorder=1200: The order of Gauss-Legendre quadrature to use for the θ integrals.
  There is no test performed to determine if this value is large enough for accurate calculation of
  all requested modal Q coefficients.  However the default value is sufficient for accurate results
  in most circumstances.
"""
function cut2sph_gauss(cut::TicraCut; pwrtol=1e-10, mmax=NMMAX, nmax=NMMAX, gaussorder=1200)
    get_ncomp(cut) == 2 || error("Must have only 2 polarization components")
    cutpwr = power(cut)
    cutθϕ = deepcopy(cut); convert_cut!(cutθϕ, 1) # Convert to Eθ and Eϕ
    θs = get_theta(cutθϕ); Nθ = length(θs); Δθ = deg2rad(θs[2] - θs[1])
    ϕs = get_phi(cutθϕ);   Nϕ = length(ϕs); Δϕ = deg2rad(ϕs[2] - ϕs[1])
    eθ = first.(get_evec(cutθϕ))
    eϕ = last.(get_evec(cutθϕ))
    # Perform ϕ integration:
    ifft!(eθ, 2);  eθ .*= 2Δϕ * Nϕ # undo scaling
    ifft!(eϕ, 2);  eϕ .*= 2Δϕ * Nϕ # undo scaling
    eθ = fftshift(eθ, 2)
    eϕ = fftshift(eϕ, 2)
    Nϕo2 = Nϕ ÷ 2
    nmax = min(Nθ - 1, nmax)
    if nmax > NMMAX
        @warn "Reducing nmax to NMMAX = $(NMMAX)."
        nmax = NMMAX
    end

    θnodes, wts  = gausslegendre(gaussorder)
    lnodes = length(θnodes)

    erange = isodd(Nϕ) ? (-Nϕo2:Nϕo2) : (-Nϕo2:(Nϕo2-1))
    mabsmax = min(last(erange), mmax)
    mabsmax > NMMAX && error("mabsmax is $(mabsmax) but may not exceed NMMAX = $(NMMAX).")
    eθoffset = OffsetArray(eθ, 1:Nθ, erange)
    eϕoffset = OffsetArray(eϕ, 1:Nθ, erange)
    Eθ, Eϕ = (OffsetArray(zeros(ComplexF64, lnodes, Nϕ), 1:lnodes, erange) for _ in 1:2)

    mrange = -mabsmax:mabsmax
    for m in mrange
        Eθfun = cubic_spline_interpolation(θs, view(eθoffset, :, m))
        Eϕfun = cubic_spline_interpolation(θs, view(eϕoffset, :, m))
        @inbounds for i in 1:lnodes
            θ = 90 * (θnodes[i] + 1)
            Eθ[i,m] = Eθfun(θ)
            Eϕ[i,m] = Eϕfun(θ)
        end
    end
    result_parent = Matrix{typeof(Dual(1.0,1.0))}(undef, nmax+1, mabsmax+1) # storage for legendre functions
    result = Origin(0)(result_parent)
    qsmns = OffsetArray(zeros(ComplexF64, (2, 2mabsmax+1, nmax)), 1:2, -mabsmax:mabsmax, 1:nmax)
    cfactor = sqrt(π)/360 # Includes 1/sqrt(2) needed for f1 and f2, and π/180 for dθ in degrees
    @inbounds for i in 1:lnodes
        θ = 90 * (θnodes[i] + 1)
        #(θ < 1.e-5 || θ > 180 - 1.e-5) && continue # No contribution from endpoints
        sinθ, cosθ = sincosd(θ)
        sin²θ = sinθ * sinθ
        P̄nm!(result_parent, nmax, mabsmax, Dual(cosθ, 1.0))
        for m in mrange
            mabs = abs(m)
            mfactor = 1
            if m > 0 && isodd(m)
                mfactor = -1
            end
            Eθᵢₘ = Eθ[i, m]
            Eϕᵢₘ = Eϕ[i, m]
            negj = negjⁿ = negjⁿ⁺¹ = complex(0,-1)
            for n in 1:nmax
                negjⁿ = negjⁿ⁺¹
                negjⁿ⁺¹ = negjⁿ * negj
                n < mabs && continue
                nfactor = inv(sqrt(n*(n+1)))
                cmn = cfactor * mfactor * nfactor
                pnm = value(result[n,mabs])
                mpnm = m*pnm
                pnm′ = derivative(result[n,mabs])
                f1factor = negjⁿ⁺¹ * cmn 
                jf1factor = complex(-imag(f1factor), real(f1factor))
                f1θconj = jf1factor * mpnm
                f1ϕconj = f1factor * (sin²θ * pnm′)
                f2factor = negjⁿ * cmn 
                jf2factor = complex(-imag(f2factor), real(f2factor))
                f2θconj = f2factor * ((-sin²θ) * pnm′)
                f2ϕconj = jf2factor * mpnm
                q1 = f1θconj * Eθᵢₘ + f1ϕconj * Eϕᵢₘ
                q2 = f2θconj * Eθᵢₘ + f2ϕconj * Eϕᵢₘ
                # Compute Gaussian quadrature:
                wt = 90 * wts[i]
                qsmns[1,m,n] += wt * q1
                qsmns[2,m,n] += wt * q2
            end
        end
    end

    # Eliminate modes with negligible power
    (qpwr, powerm, qsmn) = _filter_qmodes_by_power(qsmns, pwrtol)

    # Create output SWEQPartition
    date, clock = split(string(now()), 'T')
    funcname = nameof(var"#self#")
    prgtag = string(funcname, " ", date, " ", clock)
    idstrg = "Spherical Wave Q-Coefficients"
    nthe = (Nθ - 1) * 2
    mmax = last(axes(qsmn, 2))
    nmax = last(axes(qsmn, 3))
    nphi = Nϕ
    t4 = t5 = t6 = t7 = "Dummy Text"
    return SWEQPartition(; prgtag, idstrg, nthe, nphi, nmax, mmax, t4, t5, t6, t7, qsmns=qsmn, powerms=powerm)
end


"""
    cut2sph_gauss_generic(cut::TicraCut; keywords...) -> s::SWEQPartition

Convert a `TicraCut` object to a `SWEQPartition` using a fixed-order Gauss-Kronrod
quadrature scheme for the θ integrals.

## Keyword Arguments (and their default values)
* `mmax=$(NMMAX)`: An upper limit for the `m` (azimuthal) mode index to be included.
  The actual limit will be set to `min(Nϕ÷2, mmax)` for odd `Nϕ`, and `min(Nϕ÷2-1, mmax)`
  for even `Nϕ`, where `Nϕ` is the number of ϕ = constant polar cuts in the cut object.
* `nmax=$(NMMAX)`: An upper limit for the `n` (polar) mode index to be included.
  The actual limit will be the lesser of `nmax` and `Nθ-1` where `Nθ` is the number of 
  θ values included in each ϕ = constant polar cut.
* `pwrtol=1e-10`: The power tolerance.  Spherical modes are included until the excluded
  modes' power is less than `pwrtol` times the total modal power.  A zero or negative value
  precludes removal of any modes.
* `gaussorder=1200: The order of Gauss-Legendre quadrature to use for the θ integrals.
  There is no test performed to determine if this value is large enough for accurate calculation of
  all requested modal Q coefficients.  However the default value is sufficient for accurate results
  in most circumstances.
"""
function cut2sph_gauss_generic(cut::TicraCut; pwrtol=1e-10, mmax=NMMAX, nmax=NMMAX, gaussorder=1200)
    get_ncomp(cut) == 2 || error("Must have only 2 polarization components")
    cutpwr = power(cut)
    cutθϕ = deepcopy(cut); convert_cut!(cutθϕ, 1) # Convert to Eθ and Eϕ
    θs = get_theta(cutθϕ); Nθ = length(θs); Δθ = deg2rad(θs[2] - θs[1])
    ϕs = get_phi(cutθϕ);   Nϕ = length(ϕs); Δϕ = deg2rad(ϕs[2] - ϕs[1])
    eθ = first.(get_evec(cutθϕ))
    eϕ = last.(get_evec(cutθϕ))
    # Perform ϕ integration:
    ifft!(eθ, 2);  eθ .*= 2Δϕ * Nϕ # undo scaling
    ifft!(eϕ, 2);  eϕ .*= 2Δϕ * Nϕ # undo scaling
    eθ = fftshift(eθ, 2)
    eϕ = fftshift(eϕ, 2)
    Nϕo2 = Nϕ ÷ 2
    nmax = min(Nθ - 1, nmax)
    if nmax > NMMAX
        @warn "Reducing nmax to NMMAX = $(NMMAX)."
        nmax = NMMAX
    end

    T = typeof(float(pwrtol))
    CT = complex(T)

    θnodes_orig, wts_orig  = gausslegendre(gaussorder)
    θnodes, wts  = (θnodes_orig), (wts_orig)

    lnodes = length(θnodes)

    erange = isodd(Nϕ) ? (-Nϕo2:Nϕo2) : (-Nϕo2:(Nϕo2-1))
    mabsmax = min(last(erange), mmax)
    mabsmax > NMMAX && error("mabsmax is $(mabsmax) but may not exceed NMMAX = $(NMMAX).")
    eθoffset = OffsetArray(eθ, 1:Nθ, erange)
    eϕoffset = OffsetArray(eϕ, 1:Nθ, erange)
    Eθ, Eϕ = (OffsetArray(zeros(CT, lnodes, Nϕ), 1:lnodes, erange) for _ in 1:2)

    mrange = -mabsmax:mabsmax
    for m in mrange
        Eθfun = cubic_spline_interpolation(θs, view(eθoffset, :, m))
        Eϕfun = cubic_spline_interpolation(θs, view(eϕoffset, :, m))
        #Eθfun = scale(interpolate(view(eθoffset, :, m), BSpline(Quadratic())), θs)
        #Eϕfun = scale(interpolate(view(eϕoffset, :, m), BSpline(Quadratic())), θs)
        #Eθfun = pchip(θs, view(eθoffset, :, m))
        #Eϕfun = pchip(θs, view(eϕoffset, :, m))
        @inbounds for i in 1:lnodes
            θ = 90 * (θnodes[i] + 1)
            Eθ[i,m] = Eθfun(θ)
            Eϕ[i,m] = Eϕfun(θ)
        end
    end
    result_parent = Matrix{typeof(Dual(T(1.0),T(1.0)))}(undef, nmax+1, mabsmax+1) # storage for legendre functions
    result = Origin(0)(result_parent)
    qsmns = OffsetArray(zeros(CT, (2, 2mabsmax+1, nmax)), 1:2, -mabsmax:mabsmax, 1:nmax)
    cfactor = sqrt(T(BigFloat(π)))/360 # Includes 1/sqrt(2) needed for f1 and f2, and π/180 for dθ in degrees
    @inbounds for i in 1:lnodes
        θ = T(90 * (θnodes[i] + 1))
        #(θ < 1.e-5 || θ > 180 - 1.e-5) && continue # No contribution from endpoints
        sinθ, cosθ = sincosd(θ)
        sin²θ = sinθ * sinθ
        P̄nm!(result_parent, nmax, mabsmax, Dual(cosθ, one(T)))
        for m in mrange
            mabs = abs(m)
            mfactor = 1
            if m > 0 && isodd(m)
                mfactor = -1
            end
            Eθᵢₘ = Eθ[i, m]
            Eϕᵢₘ = Eϕ[i, m]
            negj = negjⁿ = negjⁿ⁺¹ = complex(0,-1)
            for n in 1:nmax
                negjⁿ = negjⁿ⁺¹
                negjⁿ⁺¹ = negjⁿ * negj
                n < mabs && continue
                nfactor = inv(sqrt(T(n)*(n+1)))
                cmn = cfactor * mfactor * nfactor
                pnm = value(result[n,mabs])
                mpnm = m*pnm
                pnm′ = derivative(result[n,mabs])
                f1factor = negjⁿ⁺¹ * cmn 
                jf1factor = complex(-imag(f1factor), real(f1factor))
                f1θconj = jf1factor * mpnm
                f1ϕconj = f1factor * (sin²θ * pnm′)
                f2factor = negjⁿ * cmn 
                jf2factor = complex(-imag(f2factor), real(f2factor))
                f2θconj = f2factor * ((-sin²θ) * pnm′)
                f2ϕconj = jf2factor * mpnm
                q1 = f1θconj * Eθᵢₘ + f1ϕconj * Eϕᵢₘ
                q2 = f2θconj * Eθᵢₘ + f2ϕconj * Eϕᵢₘ
                # Compute Gaussian quadrature:
                wt = 90 * wts[i]
                qsmns[1,m,n] += wt * q1
                qsmns[2,m,n] += wt * q2
            end
        end
    end

    # Eliminate modes with negligible power
    (qpwr, powerm, qsmn) = _filter_qmodes_by_power(qsmns, pwrtol)

    # Create output SWEQPartition
    date, clock = split(string(now()), 'T')
    funcname = nameof(var"#self#")
    prgtag = string(funcname, " ", date, " ", clock)
    idstrg = "Spherical Wave Q-Coefficients"
    nthe = (Nθ - 1) * 2
    mmax = last(axes(qsmn, 2))
    nmax = last(axes(qsmn, 3))
    nphi = Nϕ
    t4 = t5 = t6 = t7 = "Dummy Text"
    return SWEQPartition(; prgtag, idstrg, nthe, nphi, nmax, mmax, t4, t5, t6, t7, qsmns=qsmn, powerms=powerm)
end


function _filter_qmodes_by_power(qsmns, pwrtol)
    srange, mrange, nrange = axes(qsmns)
    mabsmax = last(mrange)
    nmax = last(nrange)
    powerms = Origin(0)(zeros(mabsmax+1))
    powerms[0] = 0.5 * sum(abs2, (qsmns[s,0,n] for s in srange for n in nrange))
    for mabs in 1:mabsmax
        powerms[mabs] = 0.5 * sum(abs2, (qsmns[s,m,n] for s in 1:2 for m in (-mabs,mabs) for n in nrange))
    end
    qpwr = sum(powerms)
    pwrtol ≤ 0 && return (qpwr, powerms, qsmns)

    # Retain only m modes with nonneglible power
    mretain = mabsmax
    if pwrtol > 0
        for mtest in 0:mabsmax
            mretain = mtest
            qpwr - sum(powerms[m] for m in 0:mretain) < pwrtol && break
        end
        qpwr = sum(powerms[m] for m in 1:mretain)
    end

    # Retain only n modes with nonnegligible power
    nretain = nmax
    if pwrtol > 0
        for ntest in 1:nmax
            nretain = ntest
            qsretained = (qsmns[s,m,n] for s in 1:2 for m in -mretain:mretain for n in 1:nretain)
            qpwr - 0.5*sum(abs2, qsretained) < pwrtol && break
        end
    end

    mretain = min(mretain, nretain)

    if (mabsmax ≠ mretain) || (nmax ≠ nretain)
        mabsmax = mretain
        nmax = nretain
        qsmn = OffsetArray(zeros(ComplexF64, 2, 2mabsmax + 1, nmax), 1:2, -mabsmax:mabsmax, 1:nmax)
        @inbounds for n in 1:nmax, m in -mabsmax:mabsmax, s in 1:2
            qsmn[s,m,n] = qsmns[s,m,n]
        end
        powerm = Origin(0)(zeros(mabsmax+1))
        powerm[0] = 0.5 * sum(abs2, qsmn[s,0,n] for s in 1:2, n in 1:nmax)
        for mabs in 1:mabsmax
            powerm[mabs] = 0.5 * sum(abs2, qsmn[s,m,n] for s in 1:2, m in (-mabs,mabs), n in 1:nmax)
        end
        qpwr = sum(powerm)
    else
        powerm = powerms
        qsmn = qsmns
    end                

    return (qpwr, powerm, qsmn)
end


"""
    sph2cut(swefile; theta, phi, ipol) -> cut::TicraCut

    sph2cut(swe; theta, phi, ipol) -> cut::TicraCut

Convert a set of Q-type spherical wave modal coefficients to far-field electric field 
values, returned as a `TicraCut` object. 

The single positional input argument can be either a string containing the name 
of a Ticra-compatible Q-type spherical wave file, or the returned value from reading 
such a file with `read_sphfile`.

## Optional Keyword Arguments:
- `theta`: An abstract range denoting the desired polar angles (colattitudes) 
  in degrees at which the field should be evaluated. If an empty range is provided (the default), then
  the values will be determined automatically by examining the modal content in `swe`.
- `phi`: An abstract range denoting the desired azimuthal angles in degrees at 
  which the field should be evaluated.  If an empty range is provided (the default), then
  the values will be determined automatically by examining the modal content in `swe`.
- `ipol`:  An integer in the range 0:3 denoting the desired polarization decomposition 
  for the computed field values. The meanings of these values are:
    - 0: (Default value) Choose the basis among choices 1, 2, or 3 that produces the largest peak copol magnitude.
    - 1: Use a (θ̂, ϕ̂) basis.
    - 2: Use a (R̂, L̂) (i.e., circular polarization) basis.
    - 3: Use a (ĥ, v̂) (Ludwig 3) basis.

## Usage Example
    cut = sph2cut("testfile.sph"; phi=0:5:355, theta=0:1:180, ipol=2)
"""
function sph2cut(sphfile::AbstractString; kwargs...)
    swe = read_sphfile(sphfile)
    sph2cut(swe; kwargs...)
end

function sph2cut(swe::SWEQPartition; 
    theta::AbstractRange=0.0:-1.0:1.0, 
    phi::AbstractRange=0.0:-1.0:1.0,
    ipol::Int = 0)

    (0 ≤ ipol ≤ 3) || Throw(ArgumentError("polarization must be either 0,1,2, or 3"))

    if isempty(theta)
        Nθ = swe.nthe ÷ 2 + 1
        θs = range(0.0, 180.0, Nθ)
    else
        Nθ = length(theta)
        θs = range(first(theta), last(theta), Nθ)
    end
    if isempty(phi)
        Nϕ = swe.nphi
        ϕs = range(0.0, 360-360/Nϕ, Nϕ)
    else
        Nϕ = length(phi)
        ϕs = range(first(phi), last(phi), Nϕ)
    end

    Es = _q2evec(swe.qsmns, θs, ϕs)
    
    if iszero(ipol)
        # Find the maximum norm E-field
        _, imaxnorm = findmax(_norm², Es)
        Eθϕ_maxnorm = Es[imaxnorm]
        # Check which polarization decomposition produces largest copol:
        Eθϕ = Eθϕ_maxnorm
        (θ̂,ϕ̂), (R̂,L̂), (ĥ,v̂) = _pol_basis_vectors(ϕ_maxnorm)
        ERL = @SMatrix[R̂ ⋅ θ̂   R̂ ⋅ ϕ̂; L̂ ⋅ θ̂   L̂ ⋅ ϕ̂] * Eθϕ
        Ehv = @SMatrix[ĥ ⋅ θ̂   ĥ ⋅ ϕ̂; v̂ ⋅ θ̂   v̂ ⋅ ϕ̂] * Eθϕ
        icomp = argmax(maximum(abs2.(x)) for x in (Eθϕ, ERL, Ehv))
    else
        icomp = ipol
    end

    # Convert to desired polarization basis
    if icomp > 1 # Other than θ̂, ϕ̂
        θ̂, ϕ̂ = _pol_basis_vectors(0.0)[1]
        # Re-express field vectors in the selected polarization components
        @inbounds for iϕ in 1:Nϕ
            b̂₁, b̂₂ = _pol_basis_vectors(ϕs[iϕ])[icomp]
            polmat = @SMatrix[b̂₁ ⋅ θ̂   b̂₁ ⋅ ϕ̂; b̂₂ ⋅ θ̂   b̂₂ ⋅ ϕ̂]
            for iθ in 1:Nθ
                Es[iθ, iϕ] = polmat * Es[iθ, iϕ]
            end
        end
    end

    date, clock = split(string(now()), 'T')
    funcname = nameof(var"#self#")
    # Create the TicraCut object:
    ncomp = 2
    icut = 1
    icomp = icomp
    text = [string("Phi = ", ϕ, ", created by ", funcname, " on ", date, " at ", clock) for ϕ in ϕs]
    theta = θs
    phi = ϕs
    evec = Es
    cut = TicraCut(;ncomp, icut, icomp, text, theta, phi, evec)

    return cut
end # function


function _q2evec(qsmns, θs, ϕs)
    mabsmax = last(axes(qsmns, 2))
    nmax = last(axes(qsmns, 3))

    # Prepare Channel for multithreading
    TD = typeof(Dual(1.0,1.0))
    chnl = Channel{Matrix{TD}}(Threads.nthreads())
    foreach(1:Threads.nthreads()) do _
        put!(chnl, Matrix{TD}(undef, nmax+1, mabsmax+1))
    end

    Nθ = length(θs)
    Nϕ = length(ϕs)

    # Precompute complex exponentials for fast lookup
    expmjmϕs_parent = ones(ComplexF64, Nϕ, mabsmax+1)
    expmjmϕs = OffsetArray(expmjmϕs_parent, 1:Nϕ, 0:mabsmax)
    expmjmϕs[:,1] .= (cis(deg2rad(-ϕ)) for ϕ in ϕs)
    for m in 2:mabsmax
        for iϕ in axes(expmjmϕs,1)
            expmjmϕs[iϕ, m] = expmjmϕs[iϕ, m-1] * expmjmϕs[iϕ, 1]
        end
    end

    (θ̂, ϕ̂) = _pol_basis_vectors(0)[1] # θ̂ and ϕ̂ are independent of ϕ in (θ, ϕ) basis
    Es = zeros(SVector{2, ComplexF64}, Nθ, Nϕ)
    cfactor = 1/(2*sqrt(π)) # Includes 1/sqrt(2) needed for f1 and f2

    Threads.@threads for iθ in eachindex(θs) # 57.6 msec
    #@inbounds for iθ in eachindex(θs) # 366.23 msec
        result_parent = take!(chnl)
        result = Origin(0)(result_parent)
        θ = θs[iθ]
        θis0 = iszero(θ)
        θis180 = θ == 180
        sinθ, cosθ = sincosd(θ)
        (θis0 || θis180) || P̄nm!(result_parent, nmax, mabsmax, Dual(cosθ, 1.0))

        @inbounds for iϕ in eachindex(ϕs)
            Eθϕ = @SVector[complex(0.0,0.0), complex(0.0,0.0)] # Initialization
            nsign = 1
            jⁿ = complex(0,1)
            jⁿ⁺¹ = jⁿ * complex(0,1)
            @inbounds for n in 1:nmax

                if θis0 || θis180
                    mPfactor = sqrt(n * (n+1) * (2n+1) / 8)
                    mrange = -1:2:1
                else
                    # General, non-endpoint θ
                    mlim = min(mabsmax, n)
                    mrange = -mlim:1:mlim
                end # Limiting cases

                nfactor = inv(sqrt(n*(n+1)))

                for m in mrange
                    mabs = abs(m)
                    if (θis0 || θis180)
                        mP = sign(m) * mPfactor
                        dP = mPfactor
                        if θis180
                            mP *= nsign
                            dP *= -nsign
                        end
                    else
                        res = result[n,mabs]
                        mP = m * value(res) / sinθ
                        dP = -sinθ * derivative(res)
                    end
                    mfactor = (m > 0 && isodd(m)) ? -1 : 1
                    cisfact = m ≥ 0 ? expmjmϕs[iϕ, m] : conj(expmjmϕs[iϕ, abs(m)])
                    #cisfact = cis(-m * deg2rad(ϕ))
                    cmn = (cfactor * mfactor * nfactor) * cisfact
                    f1factor = -jⁿ⁺¹ * cmn 
                    f2factor = jⁿ * cmn 
                    f⃗₁ = f1factor * (im * mP * θ̂  +       dP * ϕ̂)
                    f⃗₂ = f2factor * (dP * θ̂       -  im * mP * ϕ̂)
                    Eθϕ += qsmns[1,m,n] * f⃗₁ + qsmns[2,m,n] * f⃗₂
                end # m loop

                # Update n loop variables
                nsign = -nsign
                jⁿ = jⁿ⁺¹
                jⁿ⁺¹ *= 1im

            end # n loop

            Es[iθ, iϕ] = Eθϕ # Store result
 
        end # ϕ loop

        put!(chnl, result_parent)

    end # θ loop

    return Es
end # function




function compute_single_mode(cut::TicraCut, s::Int, m::Int, n::Int)
    get_ncomp(cut) == 2 || error("Must have only 2 polarization components")
    cutpwr = power(cut)
    cutθϕ = deepcopy(cut); convert_cut!(cutθϕ, 1) # Convert to Eθ and Eϕ
    θs = get_theta(cutθϕ); Nθ = length(θs); Δθ = deg2rad(θs[2] - θs[1])
    (first(θs) == 0 && last(θs) == 180) || error("bad theta limits")
    ϕs = get_phi(cutθϕ);   Nϕ = length(ϕs); Δϕ = deg2rad(ϕs[2] - ϕs[1])
    eθ = first.(get_evec(cutθϕ))
    eϕ = last.(get_evec(cutθϕ))
    # Perform ϕ integration:
    ifft!(eθ, 2);  eθ .*= 2Δϕ * Nϕ # undo scaling
    ifft!(eϕ, 2);  eϕ .*= 2Δϕ * Nϕ # undo scaling
    eθ = fftshift(eθ, 2)
    eϕ = fftshift(eϕ, 2)
    Nϕo2 = Nϕ ÷ 2

    T = Float64
    CT = complex(T)

    erange = isodd(Nϕ) ? (-Nϕo2:Nϕo2) : (-Nϕo2:(Nϕo2-1))
    eθoffset = OffsetArray(eθ, 1:Nθ, erange)
    eϕoffset = OffsetArray(eϕ, 1:Nθ, erange)

    Eθfun = cubic_spline_interpolation(θs, view(eθoffset, :, m))
    Eϕfun = cubic_spline_interpolation(θs, view(eϕoffset, :, m))

    cfactor = sqrt(T(BigFloat(π)))/360 # Includes 1/sqrt(2) needed for f1 and f2, and π/180 for dθ in degrees
    mfactor = 1
    if m > 0 && isodd(m)
        mfactor = -1
    end
    negj = complex(0,-1)
    negjⁿ = (negj)^n
    negjⁿ⁺¹ = negjⁿ * negj
    nfactor = inv(sqrt(T(n)*(n+1)))
    cmn = cfactor * mfactor * nfactor
    coefficient = quadgk(0.0, 180.0; atol=1e-12) do θ
        sinθ, cosθ = sincosd(θ)
        sin²θ = sinθ * sinθ
        res = P̄nm(n, abs(m), Dual(cosθ, one(T)))
        Eθᵢₘ = Eθfun(θ)
        Eϕᵢₘ = Eϕfun(θ)
        pnm = value(res)
        mpnm = m*pnm
        pnm′ = derivative(res)
        if s == 1
            f1factor = negjⁿ⁺¹ * cmn 
            jf1factor = complex(-imag(f1factor), real(f1factor))
            f1θconj = jf1factor * mpnm
            f1ϕconj = f1factor * (sin²θ * pnm′)
            q1 = f1θconj * Eθᵢₘ + f1ϕconj * Eϕᵢₘ
            return q1
        elseif s == 2
            f2factor = negjⁿ * cmn 
            jf2factor = complex(-imag(f2factor), real(f2factor))
            f2θconj = f2factor * ((-sin²θ) * pnm′)
            f2ϕconj = jf2factor * mpnm
            q2 = f2θconj * Eθᵢₘ + f2ϕconj * Eϕᵢₘ
            return q2
        else
            error("illegal value of s")
        end
    end
    return coefficient
end
