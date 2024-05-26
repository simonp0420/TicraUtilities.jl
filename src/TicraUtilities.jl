module TicraUtilities

export
    amplitude_db, 
    convert_cut!,
    CoorSys,
    EllipticalRim,
    eval_cut,
    get_theta, 
    get_phi,
    get_p1,
    get_p2,
    get_p3,
    get_ncomp,
    get_icut,
    get_icomp,
    get_text,
    header,
    ids,
    parse_tor_file,
    phase_deg,  
    phscen,
    power,
    read_exi,
    read_sphfile,
    read_station,
    read_tabulatedrimxyold,
    read_tfile,
    read_ticra_cut,
    read_ticra_cuts,
    SuperEllipticalRim,
    TabulatedRimXYold,
    TicraStation,
    TorObj,
    write_exi,
    write_sphfile,
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

include("TicraCut.jl")
include("TicraStation.jl")
include("Tfile.jl")
include("TicraExi.jl")
include("TorFile.jl")
include("TicraGeom.jl")
include("SphericalWave.jl")

using .SphericalWave: read_sphfile, write_sphfile

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
