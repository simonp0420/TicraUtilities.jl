module TicraUtilities

export
    amplitude_db, 
    asym2sym,
    convert_cut!,
    cut2sph, 
    cut2sph_adaptive, 
    cut2sph_gauss, 
    eval_cut,
    get_evec,
    get_header,
    get_icomp,
    get_icut,
    get_ids,
    get_ncomp,
    get_phi,
    get_text,
    get_theta, 
    maxdb,
    normalize2dir!,
    parse_tor_file,
    phase_deg,  
    phscen,
    power,
    read_cut,
    read_cuts,
    read_exi,
    read_sphfile,
    read_station,
    read_surface,
    read_tabulatedrimxyold,
    read_tfile,
    sor_efficiency,
    sph2cut,
    sym2asym,
    write_exi,
    write_sphfile,
    write_station,
    write_surface,
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
include("SphericalWave.jl")
include("Surface.jl")

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
