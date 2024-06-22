"""
  TicraExi

Data from a Ticra excitation file.
"""
struct TicraExi
  header::Vector{String}
  ampdb::Vector{Float64}
  phsdeg:: Vector{Float64}
  id::Vector{String}
end

"""
  get_header(t::TicraExi)

Return a vector of header strings.
"""
get_header(t::TicraExi) = t.header

"""
  amplitude_db(t::TicraExi)

Return a vector of amplitudes in dB
"""
amplitude_db(t::TicraExi) = t.ampdb

"""
  phase_deg(t::TicraExi)

Return a vector of phases in degrees.
"""
phase_deg(t::TicraExi) = t.phsdeg

"""
  get_ids(t::TicraExi)

Return a vector of id strings.
"""
get_ids(t::TicraExi) = t.id


import Base.complex

"""
  complex(t::TicraExi)

Return a vector complex excitation amplitudes.
"""
complex(exi::TicraExi) = [10^(ampdb/20) * cis(deg2rad(phs))  for (ampdb,phs) in zip(exi.ampdb,exi.phsdeg)]


"""
  read_exi(filename::AbstractString)
  
Reads contents of a Ticra-compatible excitation file and returns a `TicraExi` object.
"""
function read_exi(filename::AbstractString)
    open(filename,"r") do fid
        header = String[]; ampdb = Float64[]; phsdeg = Float64[]; id = String[]
        inheader = true
        for line in eachline(fid)
            if inheader
                if length(line) ≥ 4 && isequal(view(line, 1:4), "++++")
                    inheader = false
                    continue
                end
                push!(header, rstrip(line))
            else
                (i,a,p) = split(line)
                push!(id,i)
                push!(ampdb, parse(Float64,a))
                push!(phsdeg, parse(Float64,p))
            end
        end
        return TicraExi(header,ampdb,phsdeg,id)
    end
end

using Printf: @printf
"""
  write_exi(filename::AbstractString, t::TicraExi)
  
Create a Ticra-compatible excitation file from a `TicraExi` object.
"""
function write_exi(filename::AbstractString, t::TicraExi)
    open(filename, "w") do fid
        for line in t.header
            println(fid, rstrip(line))
        end
        println(fid, "++++")
        foreach(t.id, t.ampdb, t.phsdeg) do id, ampdb, phsdeg
            @printf(fid, "%s %13.6f %13.6f\n", id, ampdb, phsdeg)
        end
    end
    return nothing
end
