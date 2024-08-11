"""
    Ticra  A module to work with Ticra files 
"""
module Ticra

export
    amplitude_db, 
    CoorSys,
    EllipticalRim,
    eval_cut,
    header,
    ids,
    parse_tor_file,
    phase_deg,  
    phscen,
    power,
    read_exi,
    read_station,
    read_tabulatedrimxyold,
    read_tfile,
    read_cutfile,
    read_cutfiles,
    SuperEllipticalRim,
    TabulatedRimXYold,
    Station,
    TorObj,
    write_exi,
    write_station,
    write_ticra_cut,
    write_tor_object
    #==
       CorrList,
       curvature,
       field_rl_scan,
       insiderim, 
       make_cut_from_eh,
       newfeedfrequency, 
       phi, 
       readfeeddiam, 
       readfieldcut,
       read_ticra_sfc,
       Rim,
       search_name_index,
       smooth_sfc_to_given_curvature,
       swefile,
       TicraSfc, 
       write_feed_file,
       write_ticra_sfc,

include("eval_horn_primary_champ.jl")
=##

include("Cut.jl")
include("Station.jl")
include("Tfile.jl")
include("Exi.jl")
include("TorFile.jl")
include("Geom.jl")

#==
include("TicraSfc.jl")
include("readfieldcut.jl")
include("readfeeddiam.jl")
include("field_rl_scan.jl")
include("newfeedfrequency.jl")
include("plotcut.jl")
include("Champ.jl")
include("make_cut_from_eh.jl")
==#
end # module
