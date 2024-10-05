module TicraUtilities

export
    amplitude_db, 
    asym2sym,
    convert_cut,
    convert_cut!,
    Cut,
    cut2sph, 
    eh2bor1cut,
    Exi,
    get_ampdb,
    get_att,
    get_evec,
    get_goal,
    get_header,
    get_icomp,
    get_icut,
    get_id,
    get_ids,
    get_idstrg,
    get_ipol,
    get_mmax,
    get_name,
    get_ncomp,
    get_nmax,
    get_nphi,
    get_npoint,
    get_nthe,
    get_objtype,
    get_phi,
    get_phsdeg,
    get_powerms,
    get_prgtag,
    get_propname,
    get_propval,
    get_qsmns,
    get_rot,
    get_t4,
    get_t5,
    get_t6,
    get_t7,
    get_t8,
    get_text,
    get_theta, 
    get_u,
    get_v,
    get_weight,
    get_x,
    get_y,
    get_z,
    maximum_db,
    normalize!,
    phase_deg,  
    phscen,
    power,
    read_cutfile,
    read_exifile,
    read_sphfile,
    read_stationfile,
    read_surface,
    read_tabulatedrimxyold,
    #read_tfile,
    read_torfile,
    sor_efficiency,
    sph2cut,
    SPHQPartition,
    Station,
    Surface,
    sym2asym,
    TorObj,
    write_cutfile,
    write_exifile,
    write_sphfile,
    write_stationfile,
    write_surface,
    write_cutfile,
    write_torfile

include("Cut.jl")
include("Station.jl")
#include("Tfile.jl")
include("Exi.jl")
include("TorFile.jl")
include("Geom.jl")
include("SphericalWave.jl")
include("Surface.jl")


using PrecompileTools: @setup_workload, @compile_workload
@compile_workload begin
    cutfile = joinpath(@__DIR__, "..", "test", "test.cut")
    cut = read_cutfile(cutfile)
    adb = amplitude_db(cut)
    scut = asym2sym(cut)
    cut1 = convert_cut!(cut,1)    
    cut2 = convert_cut!(cut,2)
    sph_julia = cut2sph(cutfile)
    evec = get_evec(cut)
    icomp = get_icomp(cut)
    icut = get_icut(cut)
    ncomp = get_ncomp(cut)
    phi = get_phi(cut)
    text = get_text(cut)
    theta = get_theta(cut)
    mdb = maximum_db(cut)
    normalize!(cut)
    torfile = joinpath(@__DIR__, "..", "test", "tabulated_rim_tor_file.tor")
    torparsed = read_torfile(torfile)
    pdeg = phase_deg(cut, 1)
    pcen = phscen(cutfile)
    p = power(scut)
    cuts = read_cutfile(cutfile)
    exifile = joinpath(@__DIR__, "..", "test", "beam_A14R.exi")
    exi = read_exifile(exifile)
    sphfile = joinpath(@__DIR__, "..", "test", "center_element_rhcp_excited_q.sph")
    sph = read_sphfile(sphfile)
    stafile = joinpath(@__DIR__, "..", "test", "scenario2_coverage.sta")
    stations = read_stationfile(stafile)
    sfcfile = joinpath(@__DIR__, "..", "test", "parent_parabola.sfc")
    sfc = read_surface(sfcfile)

    t1 = sor_efficiency(cutfile; F=40.0, D=18.0, Oc=0.4, pol=:l3h, dz=0.0)
end
# include("make_cut_from_eh.jl")


end # module
