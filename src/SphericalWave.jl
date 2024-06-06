using OffsetArrays: OffsetArray, OffsetVector
using Printf: @printf
using AssociatedLegendrePolynomials: Plm!, LegendreNormCoeff, LegendreOrthoNorm, legendre, legendre!
using ForwardDiff: Dual, value, partials
derivative(x::Dual) = partials(x,1)
using FFTW: ifft!, fftshift
using OffsetArrays: OffsetArray, Origin
using LinearAlgebra: norm
using QuadGK: quadgk, kronrod
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

"""
    read_sphfile(fname) -> Vector{SWEQPartition}

Read the SWE coefficients from a Q-type spherical wave expansion file.

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
                powerms[absm] = powerm
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
"""
function write_sphfile(fname::AbstractString, sps = Vector{SWEQPartition})
    inormfactor = inv(sqrt(8π))
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
                @printf(io, "%6i %23.17G\n", absm, powerms[absm])
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

write_sphfile(fname, qs::SWEQPartition) = write_sphfile(fname, [qs])



# Functions for associated Legendre functions, both normalized and unnormalized.
# See documentation of AssociatedLegendrePolynomials.jl for calling details.
const NMMAX = 500
const Qcoef = LegendreNormCoeff{LegendreOrthoNorm,Float64}(NMMAX)
P̄nm(n, m, x) = legendre(Qcoef, n, m, x)
P̄nm!(Λ, n, m, x) = legendre!(Qcoef, Λ, n, m, x)
#Pnm(n, m, x) = Plm(n, m, x)
#Pnm!(Λ, n, m, x) = Plm!(Λ, n, m, x)

"""
    cut2sph_trap(cut::TicraCut; keywords...) -> s::SWEQPartition

Convert a `TicraCut` object to a `SWEQPartition` using a fast but less accurate
trapezoidal rule for the θ quadrature.

## Keyword Arguments (and their default values)
* `mmax=$(NMMAX)`: An upper limit for the `m` (azimuthal) mode index to be included.
  The actual limit will be set to `min(Nϕ÷2, mmax)` for odd `Nϕ`, and `min(Nϕ÷2-1, mmax)`
  for even `Nϕ`, where `Nϕ` is the number of ϕ = constant polar cuts in the cut object.
* `nmax=$(NMMAX)`: An upper limit for the `n` (polar) mode index to be included.
  The actual limit will be the lesser of `nmax` and `Nθ-1` where `Nθ` is the number of 
  θ values included in each ϕ = constant polar cut.
"""
function cut2sph_trap(cut::TicraCut; mmax=NMMAX, nmax=NMMAX)
    get_ncomp(cut) == 2 || error("Must have only 2 polarization components")
    cutθϕ = deepcopy(cut); convert_cut!(cutθϕ, 1) # Convert to Eθ and Eϕ
    θs = get_theta(cutθϕ); Nθ = length(θs); Δθ = deg2rad(θs[2] - θs[1])
    ϕs = get_phi(cutθϕ);   Nϕ = length(ϕs); Δϕ = deg2rad(ϕs[2] - ϕs[1])
    eθ = first.(get_evec(cutθϕ))
    eϕ = last.(get_evec(cutθϕ))
    # Perform ϕ integration:
    ifft!(eθ, 2);  eθ .*= Δϕ * Nϕ # undo scaling
    ifft!(eϕ, 2);  eϕ .*= Δϕ * Nϕ # undo scaling
    eθ = fftshift(eθ, 2)
    eϕ = fftshift(eϕ, 2)
    Nϕo2 = Nϕ ÷ 2
    nmax = min(Nθ - 1, nmax)
    nrange = 1:nmax
    if isodd(Nϕ)
        mabsmax = min(Nϕo2, mmax)
        Eθ = OffsetArray(eθ, 1:Nθ, -Nϕo2:Nϕo2)
        Eϕ = OffsetArray(eϕ, 1:Nθ, -Nϕo2:Nϕo2)
    else
        mabsmax = min(Nϕo2 - 1, mmax)
        Eθ = OffsetArray(eθ, 1:Nθ, -Nϕo2:(Nϕo2-1))
        Eϕ = OffsetArray(eϕ, 1:Nθ, -Nϕo2:(Nϕo2-1))
    end
    mrange = -mabsmax:mabsmax
    result_parent = Matrix{typeof(Dual(1.0,1.0))}(undef, nmax+1, mabsmax+1) # storage for legendre functions
    result = Origin(0)(result_parent)
    qsmns = OffsetArray(zeros(ComplexF64, (2, 2mabsmax+1, nmax)), 1:2, -mabsmax:mabsmax, 1:nmax)
    cfactor = Δθ * inv(sqrt(π)) # Includes extra factor 2Δθ/sqrt(2) needed for f1 and f2 and the θ quadrature
    for iθ in 1:Nθ
        intfactor = 1.0 # Trapezoid rule factor
        (iθ == 1 || iθ == Nθ) && (intfactor = 0.5)
        θ = θs[iθ]
        iszero(θ) && (θ = 1.e-5) # Avoid derivative singularity
        θ == 180 && (θ = 180.0 - 1.e-5) # Avoid derivative singularity
        sin²θ = sind(θ)^2
        P̄nm!(result_parent, nmax, mabsmax, Dual(cosd(θ), 1.0))
        for m in mrange
            mabs = abs(m)
            mfactor = 1
            if m > 0 && isodd(m)
                mfactor = -1
            end

            negjⁿ = negj = complex(0,-1)
            negjⁿ⁺¹ = complex(-1,0)
            for n in nrange
                n < mabs && continue
                negjⁿ⁺¹ = negjⁿ * negj
                nfactor = inv(sqrt(n*(n+1)))
                cmn = cfactor * mfactor * nfactor
                pnm = value(result[n,mabs])
                pnm′ = derivative(result[n,mabs])
                f1factor = negjⁿ⁺¹ * cmn 
                f1θconj = f1factor * negj * m * pnm 
                f1ϕconj = f1factor * (-sin²θ) * pnm′ 
                f2factor = negjⁿ * cmn 
                f2θconj = f2factor * (-sin²θ) * pnm′ 
                f2ϕconj = f2factor * (-negj) * m * pnm 
                qsmns[1,m,n] += intfactor * (f1θconj * Eθ[iθ, m] + f1ϕconj * Eϕ[iθ, m])
                qsmns[2,m,n] += intfactor * (f2θconj * Eθ[iθ, m] + f2ϕconj * Eϕ[iθ, m])
                negjⁿ = negjⁿ⁺¹
            end
        end
    end

    # The Q coefficients have now been calculated.
    powerms = Origin(0)(zeros(mabsmax+1))
    powerms[0] = 0.5 * sum(abs2, (qsmns[s,0,n] for s in 1:2 for n in nrange))
    for mabs in 1:mabsmax
        powerms[mabs] = 0.5 * sum(abs2, (qsmns[s,m,n] for s in 1:2 for m in (-mabs,mabs) for n in nrange))
    end
    qpwr = sum(powerms)
    return (qpwr, powerms, qsmns)
end




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
        negjⁿ = negj = complex(0,-1)
        negjⁿ⁺¹ = complex(-1,0)
        for n in nrange
            n < mabs && continue
            negjⁿ⁺¹ = negjⁿ * negj
            nfactor = inv(sqrt(n*(n+1)))
            cmn = cfactor * mfactor * nfactor
            f1factor = negjⁿ⁺¹ * cmn 
            f2factor = negjⁿ * cmn 
            int, errest = quadgk(0.0, 180.0; rtol=1e-10, atol=1e-10, order=gkorder) do θ
                (θ < θs[1] || θ > θs[end]) && return @SVector[0.0,0.0]
                sin²θ = sind(θ)^2
                result = P̄nm(n, mabs, Dual(cosd(θ), 1.0))
                pnm = value(result)
                pnm′ = derivative(result)
                f1θc = f1factor * negj * (m * pnm)
                f1ϕc = f1factor * ((-sin²θ) * pnm′)
                f2θc = f2factor * ((-sin²θ) * pnm′)
                f2ϕc = f2factor * (-negj) * (m * pnm)
                Eθ = Eθfun(θ)
                Eϕ = Eϕfun(θ)
                integrand = @SVector[f1θc * Eθ + f1ϕc * Eϕ, f2θc * Eθ + f2ϕc * Eϕ]
                return integrand
            end
            qsmns[1,m,n], qsmns[2,m,n] = int
            negjⁿ = negjⁿ⁺¹
        end
    end

    # The Q coefficients have now been calculated.
    (qpwr, powerm, qsmn) = _filter_qmodes_by_power(qsmns, pwrtol)

    # Create output SWEQPartition
    date, clock = split(string(now()), 'T')
    fname = _caller_name()
    prgtag = string(fname, " ", date, " ", clock)
    idstrg = "Spherical Wave Q-Coefficients"
    nthe = 2nmax
    mmax = mabsmax
    nphi = 2mmax + 1
    t4 = t5 = t6 = t7 = "Dummy Text"
    return SWEQPartition(; prgtag, idstrg, nthe, nphi, nmax, mmax, t4, t5, t6, t7, qsmns=qsmn, powerms=powerm)
end

const GKORDER = 400  # Default value for nonadaptive GK quadrature
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
    nrange = 1:nmax
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
        sin²θ = sind(θ)^2
        P̄nm!(result_parent, nmax, mabsmax, Dual(cosd(θ), 1.0))
        for m in mrange
            mabs = abs(m)
            mfactor = 1
            if m > 0 && isodd(m)
                mfactor = -1
            end

            negjⁿ = negj = complex(0,-1)
            negjⁿ⁺¹ = complex(-1,0)
            for n in nrange
                n < mabs && continue
                negjⁿ⁺¹ = negjⁿ * negj
                nfactor = inv(sqrt(n*(n+1)))
                cmn = cfactor * mfactor * nfactor
                pnm = value(result[n,mabs])
                pnm′ = derivative(result[n,mabs])
                f1factor = negjⁿ⁺¹ * cmn 
                f1θconj = f1factor * negj * (m * pnm)
                f1ϕconj = f1factor * ((-sin²θ) * pnm′)
                f2factor = negjⁿ * cmn 
                f2θconj = f2factor * ((-sin²θ) * pnm′)
                f2ϕconj = f2factor * (-negj) * (m * pnm)
                q1 = f1θconj * Eθ[i, m] + f1ϕconj * Eϕ[i, m]
                q2 = f2θconj * Eθ[i, m] + f2ϕconj * Eϕ[i, m]
                # Compute full-order Gaussian-Kronrod quadrature:
                qsmns[1,m,n] += wts[i] * q1
                qsmns[2,m,n] += wts[i] * q2
                if iseven(i)
                    # Compute lower-order Gaussian quadrature:
                    ii = i÷2
                    qsmns_low[1,m,n] += gwts[ii] * q1
                    qsmns_low[2,m,n] += gwts[ii] * q2
                end
                negjⁿ = negjⁿ⁺¹
            end
        end
    end

    # The Q coefficients have now been calculated. Check GK convergence:
    gkabserrtol = 1e-6
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
    fname = _caller_name()
    prgtag = string(fname, " ", date, " ", clock)
    idstrg = "Spherical Wave Q-Coefficients"
    nthe = 2nmax
    mmax = mabsmax
    nphi = 2mmax + 1
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


function _caller_name()
    s = string(stacktrace()[2])
    i = findfirst('(', s)
    isnothing(i) && return ""
    return s[1:i-1]
end