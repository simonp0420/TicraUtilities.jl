using OffsetArrays: OffsetArray, OffsetVector, Origin
using Printf: @printf
using AssociatedLegendrePolynomials: Plm!, LegendreNormCoeff, LegendreOrthoNorm, legendre, legendre!, LegendreOrthoCoeff
using ForwardDiff: Dual, value, partials
derivative(x::Dual) = partials(x, 1)
using FFTW: fft!, bfft!, ifft!, fftshift
using LinearAlgebra: norm
using Dates: now


"""
    SPHQPartition
Struct containing all the data read from a Ticra-compatible Q-type SPH file for one frequency.
The `qsmns` field contains the Q coefficients. It is indexed as `qsmns[s,m,n]` where `s ∈ {1,2}`, `m ∈ -mmax:mmax`, 
`n ∈ {1, ..., mmax}`.  But not all entries are used.  The only possible nonzero entries are where 
`n ∈ {max(1, abs(m)), ..., nmax}`.  Note that if the coefficients stored in the SPH files are `Q'`, and
the cofficients stored in a `SPHQPartition` instance are `Q`, then `Q = sqrt(8π) * conj(Q')`, as discussed
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
* `powerms::OffsetArray{Float64, 1}`: The vector has axes `(0:mmax)` and the `m`th element
  contains the total power (one-half the sum of the magnitude squared) of all modes with `±m` as the m index.
"""
@kwdef struct SPHQPartition
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
    qsmns::OffsetArray{ComplexF64,3,Array{ComplexF64,3}} =
        OffsetArray(zeros(ComplexF64, (2, 2mmax + 1, nmax)), 1:2, -mmax:mmax, 1:nmax)
    powerms::OffsetVector{Float64,Vector{Float64}} = OffsetArray(zeros(mmax + 1), 0:mmax)
end

"""
    get_prgtag(s::SPHQPartition)

Return the `prgtag` string of `s`.
"""
get_prgtag(s::SPHQPartition) = s.prgtag

"""
    get_idstrg(s::SPHQPartition)

Return the `idstrg` string of `s`.
"""
get_idstrg(s::SPHQPartition) = s.idstrg

"""
    get_nthe(s::SPHQPartition)

Return the integer `nthe` (number of theta points) associated with `s`.
"""
get_nthe(s::SPHQPartition) = s.nthe

"""
    get_nphi(s::SPHQPartition)

Return the integer `nphi` (number of phi points) associated with `s`.
"""
get_nphi(s::SPHQPartition) = s.nphi

"""
    get_nmax(s::SPHQPartition)

Return the integer `nmax` (maximum `n` index) associated with `s`.
"""
get_nmax(s::SPHQPartition) = s.nmax

"""
    get_mmax(s::SPHQPartition)

Return the integer `mmax` (maximum `m` index) associated with `s`.
"""
get_mmax(s::SPHQPartition) = s.mmax

"""
    get_t4(s::SPHQPartition)

Return the string `t4` associated with `s`.
"""
get_t4(s::SPHQPartition) = s.t4

"""
    get_t5(s::SPHQPartition)

Return the string `t5` associated with `s`.
"""
get_t5(s::SPHQPartition) = s.t5

"""
    get_t6(s::SPHQPartition)

Return the string `t6` associated with `s`.
"""
get_t6(s::SPHQPartition) = s.t6

"""
    get_t7(s::SPHQPartition)

Return the string `t7` associated with `s`.
"""
get_t7(s::SPHQPartition) = s.t7

"""
    get_t8(s::SPHQPartition)

Return the string `t8` associated with `s`.
"""
get_t8(s::SPHQPartition) = s.t8

"""
    get_qsmns(s::SPHQPartition)

Return the offset array `qsmns::OffsetArray{ComplexF64, 3}` that contains the 
Q coefficients, assuming exp(jωt) time variation and using the Ticra normalization 
convention.  The array has axes `(1:2, -mmax:mmax, 1:nmax)`, meaning  that some of 
the entries (those with `n < abs(m)`) are always zero.
"""
get_qsmns(s::SPHQPartition) = s.qsmns

"""
    get_powerms(s::SPHQPartition)

Return the offset vector `powerms`. The vector has axes `(0:mmax)` (i.e. it can be indexed 
with integers `m` ranging from `0` to `mmax`) and its `m`th element contains the total power
(one-half the sum of the magnitude squared) of all modes with `±m` as the `m` index.
"""
get_powerms(s::SPHQPartition) = s.powerms

import Base.show
function show(io::IO, mime::MIME"text/plain", sph::SPHQPartition)
    println(io, "SPHQPartition")
    println(io, "  prgtag   $(sph.prgtag)")
    println(io, "  idstrg   $(sph.idstrg)")
    println(io, "  nthe     $(sph.nthe)")
    println(io, "  nphi     $(sph.nphi)")
    println(io, "  nmax     $(sph.nmax)")
    println(io, "  mmax     $(sph.mmax)")
    println(io, string("  qsmns    OffsetArray{ComplexF64}(",
        first(axes(sph.qsmns, 1)), ":", last(axes(sph.qsmns, 1)), ",",
        first(axes(sph.qsmns, 2)), ":", last(axes(sph.qsmns, 2)), ",",
        first(axes(sph.qsmns, 3)), ":", last(axes(sph.qsmns, 3)), ")"))
    println(io, string("  powerms  OffsetArray{Float64}(",
        first(axes(sph.powerms, 1)), ":", last(axes(sph.powerms, 1)), ")"))
    return nothing
end

function show(io::IO, sph::SPHQPartition)
    print(io, "SPHQPartition with nmax = $(sph.nmax), mmax = $(sph.mmax)")
    return nothing
end




"""
    read_sphfile(fname) -> Vector{SPHQPartition}

Read the SPH coefficients from a Q-type spherical wave expansion file. 

In the process of reading the data, the coefficients in the file (Q′) are conjugated and then
multiplied by the factor √(8π) to achieve Ticra-standard normalization.  
Each element of the returned vector corresponds to a particular operating frequency partition
in the file.  If there is only a single partition in the file, then instead of returning a 1-element
vector, the single element of type `SPHQPartition` is returned as a scalar.
"""
function read_sphfile(fname::AbstractString)
    normfactor = sqrt(8π)
    open(fname, "r") do io
        sps = SPHQPartition[]
        while !eof(io)
            # Read header info of next partition:
            prgtag, idstrg = (readline(io) for _ in 1:2)
            nthe, nphi, nmax, mmax = parse.(Int, split(readline(io)))
            t4, t5, t6, t7, t8 = (readline(io) for _ in 4:8)

            qsmns = OffsetArray(zeros(ComplexF64, 2, 2mmax + 1, nmax), 1:2, -mmax:mmax, 1:nmax)
            powerms = OffsetArray(zeros(mmax + 1), 0:mmax)
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
            push!(sps, SPHQPartition(; prgtag, idstrg, nthe, nphi, nmax, mmax, t4, t5, t6, t7, t8, qsmns, powerms))
        end
        if length(sps) > 1
            return sps
        else
            return sps[1]
        end
    end
end



"""
    write_sphfile(fname, qs::Vector{SPHQPartition})
    write_sphfile(fname, qs::SPHQPartition)

Write SPH coefficients to a Q-type spherical wave expansion file.

In the process of writing the data, the coefficients in the file (Q) are conjugated and then
multiplied by the factor 1/sqrt(8π) to become Q′ and achieve consistency with Ticra-standard normalization.  
"""
function write_sphfile(fname::AbstractString, sps::Vector{SPHQPartition})
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
                @printf(io, "%6i %23.17G\n", absm, ipwrnormfactor * powerms[absm])
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

write_sphfile(fname::AbstractString, qs::SPHQPartition) = write_sphfile(fname, [qs])



# Functions for associated Legendre functions, both normalized and unnormalized.
# See documentation of AssociatedLegendrePolynomials.jl for calling details.
NMMAX::Int = 500
Qcoef::LegendreOrthoCoeff{Float64} = LegendreNormCoeff{LegendreOrthoNorm,Float64}(NMMAX)

function _checkmn(m, n)
    global NMMAX, Qcoef
    nmmax = max(maximum(abs.(m)), maximum(n))
    if nmmax > NMMAX
        NMMAX = nmmax
        Qcoef = LegendreNormCoeff{LegendreOrthoNorm,Float64}(NMMAX)
    end
    return
end

"""
    P̄nm(n, m, x)

Normalized associated Legendre function following Ticra definition.
See documentation for `AssociatedLegendrePolynomials.Plm` for calling conventions.
"""
@inline function P̄nm(n::Int, m::Int, x::Number)
    _checkmn(m, n)
    m < 0 && throw(ArgumentError("m must be nonnegative"))
    p = legendre(Qcoef, n, m, x)
    # Remove Condon–Shortley phase factor
    return isodd(m) ? -p : p
end

@inline function P̄nm(n::UnitRange, m::UnitRange, x::Number)
    any(<(0), m) && throw(ArgumentError("all m values must be nonnegative"))
    _checkmn(m, n)
    p = legendre(Qcoef, n, m, x)
    # Remove Condon–Shortley phase factor
    for mp1 in axes(p, 2)
        isodd(mp1) && continue # 1-based indexing!
        p[:, mp1] .*= -1
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
    _checkmn(m, n)
    legendre!(Qcoef, Λ, n, m, x)
    # Remove Condon–Shortley phase factor
    for mp1 in axes(Λ, 2)
        isodd(mp1) && continue # 1-based indexing!
        Λ[:, mp1] .*= -1
    end
    return nothing
end

#P̄nm!(Λ, n, m, x) = legendre!(Qcoef, Λ, n, m, x)
#Pnm(n, m, x) = Plm(n, m, x)
#Pnm!(Λ, n, m, x) = Plm!(Λ, n, m, x)



"""
    sph2cut(sphfile::AbstractString; theta, phi, ipol) -> cut::Cut

    sph2cut(sph:SPHQPartition; theta, phi, ipol) -> cut::Cut

    sph2cut(sphs:Vector{SPHQPartition}; theta, phi, ipol) -> cuts::Vector{Cut}

Convert a set of Q-type spherical wave modal coefficients to far-field electric field 
values, returned as a `Cut` object. 

The single positional input argument can be either a string containing the name 
of a Ticra-compatible Q-type spherical wave file, or the returned value from reading 
such a file with `read_sphfile`.

## Optional Keyword Arguments:
- `theta`: An abstract range denoting the desired polar angles (colattitudes) 
  in degrees at which the field should be evaluated. If an empty range is provided (the default), then
  the values will be determined automatically by examining the modal content in `sph`.
- `phi`: An abstract range denoting the desired azimuthal angles in degrees at 
  which the field should be evaluated.  If an empty range is provided (the default), then
  the values will be determined automatically by examining the modal content in `sph`.
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
    sph = read_sphfile(sphfile)
    sph2cut(sph; kwargs...)
end

sph2cut(sphs::AbstractVector{SPHQPartition}; kwargs...) = [sph2cut(sph; kwargs...) for sph in sphs]

function sph2cut(sph::SPHQPartition;
    theta::AbstractRange=0.0:-1.0:1.0,
    phi::AbstractRange=0.0:-1.0:1.0,
    ipol::Int=0)

    (0 ≤ ipol ≤ 3) || Throw(ArgumentError("polarization must be either 0,1,2, or 3"))

    if isempty(theta)
        Nθ = sph.nthe ÷ 2 + 1
        θs = range(0.0, 180.0, Nθ)
    else
        Nθ = length(theta)
        θs = range(first(theta), last(theta), Nθ)
    end
    if isempty(phi)
        Nϕ = sph.nphi
        ϕs = range(0.0, 360 - 360 / Nϕ, Nϕ)
    else
        Nϕ = length(phi)
        ϕs = range(first(phi), last(phi), Nϕ)
    end

    Es = _q2evec(sph.qsmns, θs, ϕs)

    if iszero(ipol)
        # Find the maximum norm E-field
        _, imaxnorm = findmax(_norm², Es)
        Eθϕ_maxnorm = Es[imaxnorm]
        ϕ_maxnorm = ϕs[last(Tuple(imaxnorm))]
        # Check which polarization decomposition produces largest copol:
        Eθϕ = Eθϕ_maxnorm
        (θ̂, ϕ̂), (R̂, L̂), (ĥ, v̂) = _pol_basis_vectors(ϕ_maxnorm)
        ERL = @SMatrix[R̂⋅θ̂ R̂⋅ϕ̂; L̂⋅θ̂ L̂⋅ϕ̂] * Eθϕ
        Ehv = @SMatrix[ĥ⋅θ̂ ĥ⋅ϕ̂; v̂⋅θ̂ v̂⋅ϕ̂] * Eθϕ
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
            polmat = @SMatrix[b̂₁⋅θ̂ b̂₁⋅ϕ̂; b̂₂⋅θ̂ b̂₂⋅ϕ̂]
            for iθ in 1:Nθ
                Es[iθ, iϕ] = polmat * Es[iθ, iϕ]
            end
        end
    end

    date, clock = split(string(now()), 'T')
    funcname = nameof(var"#self#")
    # Create the Cut object:
    ncomp = 2
    icut = 1
    icomp = icomp
    text = [string("Phi = ", ϕ, ", created by ", funcname, " on ", date, " at ", clock) for ϕ in ϕs]
    theta = θs
    phi = ϕs
    evec = Es
    cut = Cut(; ncomp, icut, icomp, text, theta, phi, evec)

    return cut
end # function


function _q2evec(qsmns, θs, ϕs)
    mabsmax = last(axes(qsmns, 2))
    nmax = last(axes(qsmns, 3))
    _checkmn(mabsmax, nmax)

    # Prepare Channel for multithreading
    TD = typeof(Dual(1.0, 1.0))
    chnl = Channel{Matrix{TD}}(Threads.nthreads())
    foreach(1:Threads.nthreads()) do _
        put!(chnl, Matrix{TD}(undef, nmax + 1, mabsmax + 1))
    end

    Nθ = length(θs)
    Nϕ = length(ϕs)

    # Precompute complex exponentials for fast lookup
    expmjmϕs_parent = ones(ComplexF64, Nϕ, mabsmax + 1)
    expmjmϕs = OffsetArray(expmjmϕs_parent, 1:Nϕ, 0:mabsmax)
    expmjmϕs[:, 1] .= (cis(deg2rad(-ϕ)) for ϕ in ϕs)
    for m in 2:mabsmax
        for iϕ in axes(expmjmϕs, 1)
            expmjmϕs[iϕ, m] = expmjmϕs[iϕ, m-1] * expmjmϕs[iϕ, 1]
        end
    end

    (θ̂, ϕ̂) = _pol_basis_vectors(0)[1] # θ̂ and ϕ̂ are independent of ϕ in (θ, ϕ) basis
    Es = zeros(SVector{2,ComplexF64}, Nθ, Nϕ)
    cfactor = 1 / (2 * sqrt(π)) # Includes 1/sqrt(2) needed for f1 and f2

    Threads.@threads for iθ in eachindex(θs)
        result_parent = take!(chnl)
        result = Origin(0)(result_parent)
        θ = θs[iθ]
        θis0 = iszero(θ)
        θis180 = θ == 180
        sinθ, cosθ = sincosd(θ)
        (θis0 || θis180) || P̄nm!(result_parent, nmax, mabsmax, Dual(cosθ, 1.0))

        @inbounds for iϕ in eachindex(ϕs)
            Eθϕ = @SVector[complex(0.0, 0.0), complex(0.0, 0.0)] # Initialization
            nsign = 1
            jⁿ = complex(0, 1)
            jⁿ⁺¹ = jⁿ * complex(0, 1)
            @inbounds for n in 1:nmax

                if θis0 || θis180
                    mPfactor = sqrt(n * (n + 1) * (2n + 1) / 8)
                    mrange = -1:2:1
                else
                    # General, non-endpoint θ
                    mlim = min(mabsmax, n)
                    mrange = -mlim:1:mlim
                end # Limiting cases

                nfactor = inv(sqrt(n * (n + 1)))

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
                        res = result[n, mabs]
                        mP = m * value(res) / sinθ
                        dP = -sinθ * derivative(res)
                    end
                    mfactor = (m > 0 && isodd(m)) ? -1 : 1
                    cisfact = m ≥ 0 ? expmjmϕs[iϕ, m] : conj(expmjmϕs[iϕ, abs(m)])
                    #cisfact = cis(-m * deg2rad(ϕ))
                    cmn = (cfactor * mfactor * nfactor) * cisfact
                    f1factor = -jⁿ⁺¹ * cmn
                    f2factor = jⁿ * cmn
                    f⃗₁ = f1factor * (im * mP * θ̂ + dP * ϕ̂)
                    f⃗₂ = f2factor * (dP * θ̂ - im * mP * ϕ̂)
                    Eθϕ += qsmns[1, m, n] * f⃗₁ + qsmns[2, m, n] * f⃗₂
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



# 
# Beginning of implementation of algorithms from Hansen, 1988: "Spherical Near-Field Antenna Measurements".
# Equation numbers referenced below refer to those of this book.
#

"""
    _Δⁿₙₘ(m::Integer,n::Integer)

Starting value of Δ iterations from Hansen Eq. (4.95), after correcting sign error in denominator
"""
function _Δⁿₙₘ(m::Integer, n::Integer)
    n > 0 || throw(ArgumentError("n must be positive"))
    mabs = abs(m)
    mabs ≤ n || throw(ArgumentError("|m| must be less than or equal to n"))
    return sqrt(prod((n - mabs + i) / 2i for i in 1:(n+mabs)) / 2.0^(n - mabs)) # Hansen (4.95) after correction
end

"""
    _Δⁿₘₚₘ(m′::Integer, m::Integer,n::Integer)

Δ calculation from Hansen Eq. (4.94). This is not efficient and is only for test purposes.
"""
function _Δⁿₘₚₘ(m′::Integer, m::Integer, n::Integer)
    n > 0 || throw(ArgumentError("n must be positive"))
    -n ≤ m ≤ n || throw(ArgumentError("|m| must be less than or equal to n"))
    0 ≤ m′ ≤ n || throw(ArgumentError("m′ must be nonnegative and less than or equal to n"))
    Δip1 = 0.0
    Δi = _Δⁿₙₘ(m, n)
    for i in n:-1:m′+1
        Δim1 = -(sqrt((n + i + 1) * (n - i)) * Δip1 + 2m * Δi) / sqrt((n + i) * (n - i + 1))
        Δip1, Δi = Δi, Δim1
    end
    return Δi
end

Π(j) = isodd(j) ? 0.0 : 2.0 / (1 - j^2) # Eq. (4.84)

"""
    cut2sph(cut::Cut; keywords...) -> sph::SPHQPartition

    cut2sph(cuts::AbstractVector{Cut}; keywords...) -> sphs::Vector{SPHQPartition}

    cut2sph(cutfile::AbstractString; kwargs...) -> s::SPHQPartition

Convert a `Cut` object to a `SPHQPartition` using recursive FFT/IFFT methods from
the Hansen 1988 book "Spherical Near-Field Antenna Measurements.

The single positional input argument can be either a string containing the name 
of a Ticra-compatible, spherical polar cut file, or the returned value of type `Cut` 
that results from reading such a file with `read_cutfile`.  The output of this function
can be passed to `write_sphfile` to create a Ticra-compatible file of Q-type 
spherical wave coefficients.

If the input cuts extend in θ only to θ₀ < 180°, then it will be assumed that
the fields are identically zero for θ₀ < θ ≤ 180°.

## Keyword Arguments (and their default values)
* `mmax=1000`: An upper limit for the `m` (azimuthal) mode index to be included.
  The actual limit will be set to `min(Nϕ÷2, mmax)` for odd `Nϕ`, and `min(Nϕ÷2-1, mmax)`
  for even `Nϕ`, where `Nϕ` is the number of ϕ = constant polar cuts in the cut object.
* `nmax=1000`: An upper limit for the `n` (polar) mode index to be included.
  The actual limit will be the lesser of `nmax` and `Nθ-1` where `Nθ` is the number of 
  θ values included in each ϕ = constant polar cut.
* `pwrtol=0.0`: The power tolerance.  Spherical modes are included until the excluded
  modes' power is less than `pwrtol` times the total modal power.  A zero or negative value
  precludes removal of any modes.
"""
function cut2sph(cutfile::AbstractString; kwargs...)
    cut = read_cutfile(cutfile)
    cut2sph(cut; kwargs...)
end

cut2sph(cuts::AbstractVector{<:Cut}; kwargs...) = [cut2sph(cut; kwargs...) for cut in cuts]

function cut2sph(cut::Cut; pwrtol=0.0, mmax=1000, nmax=1000)
    get_ncomp(cut) == 2 || error("Cut must have only 2 polarization components")
    cutθϕ = deepcopy(cut)
    convert_cut!(cutθϕ, 1) # Convert to Eθ and Eϕ
    θs = get_theta(cutθϕ)
    Nθ = length(θs)
    Δθ = θs[2] - θs[1]
    iszero(first(θs)) || error("First θ value in cut must be zero")
    180 / Δθ ≈ round(Int, 180 / Δθ) || error("Δθ must divide evenly into 180 in cut")
    ϕs = get_phi(cutθϕ)
    Nϕ = length(ϕs)
    if 180 ≠ last(θs)
        # Extend cut to θ = 180 with zeros
        newθs = (last(θs)+Δθ):Δθ:180
        newzeros = zeros(eltype(cut.evec), length(newθs), Nϕ)
        cutθϕ.evec = vcat(cutθϕ.evec, newzeros)
        cutθϕ.theta = first(θs):Δθ:180
        θs = get_theta(cutθϕ)
        Nθ = length(θs)
    end

    # Step 1: CP components using Eθ and Eϕ components
    W₊₁ₘ = [first(e) + im * last(e) for e in get_evec(cutθϕ)]
    W₋₁ₘ = [first(e) - im * last(e) for e in get_evec(cutθϕ)]

    # Step 2: ϕ integration. Eq. (4.127):
    ifft!(W₊₁ₘ, 2)
    ifft!(W₋₁ₘ, 2)

    # Create storage for extended samples in theta
    Nθe = round(Int, 360 / Δθ)
    Nθe == 2(Nθ - 1) || error("Nθe = $Nθe is not equal to 2(Nθ-1) = $(2(Nθ-1))")
    vp1 = zeros(ComplexF64, Nθe)  # For one column of w
    vm1 = zeros(ComplexF64, Nθe)  # For one column of w

    N = min(Nθ - 1, nmax) # Maximum modal n value
    Nϕo2 = Nϕ ÷ 2
    M = isodd(Nϕ) ? min(Nϕo2, mmax) : min(Nϕo2 - 1, mmax) # Max modal m value

    b̃p1 = zeros(ComplexF64, 4N) # Extended Fourier coefficients
    b̃m1 = zeros(ComplexF64, 4N)

    # Step 12: Periodic extension of the Pi function
    Π̃ = ComplexF64[Π(j > 2N ? j - 4N : j) for j in 0:(4N-1)]
    ftΠ̃ = fft!(Π̃)

    qsmns = OffsetArray(zeros(ComplexF64, (2, 2M + 1, N)), 1:2, -M:M, 1:N)
    sqroots = Origin(0)(Float64[sqrt(p) for p in 0:2N+1])
    cfactors = [-sqrt((n + 0.5) * π) * (im)^n for n in 1:N]
    for m in -M:M
        mfactorp1 = (1.0im)^(1 + m)
        # Fill extended vectors
        mμsign = iseven(1 - m) ? 1 : -1
        column = m ≥ 0 ? m + 1 : Nϕ + m + 1
        vp1 .= zero(ComplexF64)
        vm1 .= zero(ComplexF64)
        vp1[1:Nθ] .= @view W₊₁ₘ[:, column]
        vm1[1:Nθ] .= @view W₋₁ₘ[:, column]
        for (i, ie) in enumerate(Nθ+1:Nθe)
            vp1[ie] = mμsign * W₊₁ₘ[Nθ-i, column]
            vm1[ie] = mμsign * W₋₁ₘ[Nθ-i, column]
        end

        # Step 3, Eq. 4.128
        bp1 = ifft!(vp1)
        bm1 = ifft!(vm1)
        # Eq. 4.87:
        b̃p1 .= zero(ComplexF64)
        b̃m1 .= zero(ComplexF64)
        Np1 = N + 1
        b̃p1[1:Np1] .= @view bp1[1:Np1]
        b̃p1[Np1] *= 0.5
        b̃m1[1:Np1] .= @view bm1[1:Np1]
        b̃m1[Np1] *= 0.5
        Nm1 = N - 1
        b̃p1[(end-Nm1):end] .= @view bp1[(end-Nm1):end]
        b̃p1[end-Nm1] *= 0.5
        b̃m1[(end-Nm1):end] .= @view bm1[(end-Nm1):end]
        b̃m1[end-Nm1] *= 0.5
        # Step 4: Compute K sequences using convolution
        fft!(b̃p1)
        b̃p1 .*= ftΠ̃
        fft!(b̃m1)
        b̃m1 .*= ftΠ̃
        Kp1 = ifft!(b̃p1)
        Km1 = ifft!(b̃m1)
        for n in max(1, abs(m)):N
            Δiₘ, Δip1ₘ = _Δⁿₙₘ(m, n), 0.0 # For (m', m) recursion
            Δi₁, Δip1₁ = _Δⁿₙₘ(1, n), 0.0  # For (m',μ) recursion (1 subscript means μ=+1)
            sp1 = sm1 = zero(ComplexF64) # Sums for μ = ±1
            mprimefact = 1 # See Eq. (A2.32)
            for i in n:-1:1 # i plays role of m′
                Δi₋₁ = mprimefact * Δi₁ # μ=-1 via Eq. (A2.32)
                sp1 += Δiₘ * Δi₁ * Kp1[i+1]
                sm1 += Δiₘ * Δi₋₁ * Km1[i+1]
                # Recursion:
                root1 = sqroots[n+i+1] * sqroots[n-i]
                root2 = sqroots[n+i] * sqroots[n-i+1]
                Δim1ₘ = -(root1 * Δip1ₘ + 2m * Δiₘ) / root2
                Δip1ₘ, Δiₘ = Δiₘ, Δim1ₘ
                Δim1₁ = -(root1 * Δip1₁ + 2 * Δi₁) / root2
                Δip1₁, Δi₁ = Δi₁, Δim1₁
                mprimefact = -mprimefact
            end
            # Do m′ = 0:
            Δi₋₁ = mprimefact * Δi₁
            sp1 = Δiₘ * Δi₁ * Kp1[1] + 2 * sp1
            sm1 = Δiₘ * Δi₋₁ * Km1[1] + 2 * sm1

            nsign = iseven(n) ? -1 : 1
            cfactor = nsign * cfactors[n] * mfactorp1
            qsmns[1, m, n] = cfactor * (sp1 - sm1)
            qsmns[2, m, n] = cfactor * (sp1 + sm1)
        end
    end


    # Eliminate modes with negligible power
    (qpwr, powerms, qsmn) = _filter_qmodes_by_power(qsmns, pwrtol)


    # Create output SPHQPartition
    date, clock = split(string(now()), 'T')
    funcname = nameof(var"#self#")
    prgtag = string(funcname, " ", date, " ", clock)
    idstrg = "Spherical Wave Q-Coefficients"
    nthe = Nθe
    M = last(axes(qsmn, 2))
    N = last(axes(qsmn, 3))
    nphi = Nϕ
    t4 = t5 = t6 = t7 = "Dummy Text"
    return SPHQPartition(; prgtag, idstrg, nthe, nphi, nmax=N, mmax=M, t4, t5, t6, t7, qsmns=qsmn, powerms)
end

function _filter_qmodes_by_power(qsmns, pwrtol)
    srange, mrange, nrange = axes(qsmns)
    mabsmax = last(mrange)
    nmax = last(nrange)
    powerms = Origin(0)(zeros(mabsmax + 1))
    powerms[0] = 0.5 * sum(abs2, (qsmns[s, 0, n] for s in srange for n in nrange))
    for mabs in 1:mabsmax
        powerms[mabs] = 0.5 * sum(abs2, (qsmns[s, m, n] for s in 1:2 for m in (-mabs, mabs) for n in nrange))
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
            qsretained = (qsmns[s, m, n] for s in 1:2 for m in -mretain:mretain for n in 1:nretain)
            qpwr - 0.5 * sum(abs2, qsretained) < pwrtol && break
        end
    end

    mretain = min(mretain, nretain)

    if (mabsmax ≠ mretain) || (nmax ≠ nretain)
        mabsmax = mretain
        nmax = nretain
        qsmn = OffsetArray(zeros(ComplexF64, 2, 2mabsmax + 1, nmax), 1:2, -mabsmax:mabsmax, 1:nmax)
        @inbounds for n in 1:nmax, m in -mabsmax:mabsmax, s in 1:2
            qsmn[s, m, n] = qsmns[s, m, n]
        end
        powerm = Origin(0)(zeros(mabsmax + 1))
        powerm[0] = 0.5 * sum(abs2, qsmn[s, 0, n] for s in 1:2, n in 1:nmax)
        for mabs in 1:mabsmax
            powerm[mabs] = 0.5 * sum(abs2, qsmn[s, m, n] for s in 1:2, m in (-mabs, mabs), n in 1:nmax)
        end
        qpwr = sum(powerm)
    else
        powerm = powerms
        qsmn = qsmns
    end

    return (qpwr, powerm, qsmn)
end

