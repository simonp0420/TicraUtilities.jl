"""
    SphericalWave

Module for spherical wave coefficients
"""
module SphericalWave

export read_sphfile, write_sphfile

using OffsetArrays: OffsetArray, OffsetVector
using Printf: @printf


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



end # module