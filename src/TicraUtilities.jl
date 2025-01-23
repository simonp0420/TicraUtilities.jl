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
    get_class,
    get_evec,
    get_freqs,
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
    get_sff,
    get_sfr,
    get_srf,
    get_srr,
    get_t4,
    get_t5,
    get_t6,
    get_t7,
    get_t8,
    get_text,
    get_theta, 
    get_title,
    get_u,
    get_v,
    get_weight,
    get_x,
    get_y,
    get_z,
    issym,
    maximum_db,
    normalize!,
    parse_tor_struct,
    phase_deg,  
    phscen,
    power,
    read_cutfile,
    read_exifile,
    read_sphfile,
    read_stationfile,
    read_surface,
    read_tabulatedrimxyold,
    read_tepfile,
    read_torfile,
    sor_efficiency,
    sph2cut,
    SPHQPartition,
    Station,
    Surface,
    sym2asym,
    symsqueeze,
    TEPperiodic,
    TEPscatter,
    tepp2s,
    teps2p,
    TorObj,
    write_cutfile,
    write_exifile,
    write_sphfile,
    write_stationfile,
    write_surface,
    write_tepfile,
    write_torfile

include("Cut.jl")
include("Station.jl")
#include("Tfile.jl")
include("Exi.jl")
include("TorFile.jl")
#include("Geom.jl")
include("SphericalWave.jl")
include("Surface.jl")
include("TEP.jl")


using PrecompileTools: @setup_workload, @compile_workload
using Logging: with_logger, NullLogger
@compile_workload begin
    cutfile = joinpath(pkgdir(@__MODULE__), "test", "test.cut")
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
    t5 = read_cutfile(joinpath(pkgdir(@__MODULE__), "test", "issue20_cut.cut"))
    at5 = with_logger(() -> sym2asym(t5), NullLogger()) 
    st5 = asym2sym(at5)
    at6 = sym2asym(st5)

    #=
    normalize!(cut)
    torfile = joinpath(pkgdir(@__MODULE__), "test", "tabulated_rim_tor_file.tor")
    torparsed = read_torfile(torfile)
    pdeg = phase_deg(cut, 1)
    pcen = phscen(cutfile)
    p = power(scut)
    cuts = read_cutfile(cutfile)
    exifile = joinpath(pkgdir(@__MODULE__), "test", "beam_A14R.exi")
    exi = read_exifile(exifile)
    sphfile = joinpath(pkgdir(@__MODULE__), "test", "center_element_rhcp_excited_q.sph")
    sph = read_sphfile(sphfile)
    stafile = joinpath(pkgdir(@__MODULE__), "test", "scenario2_coverage.sta")
    stations = read_stationfile(stafile)
    sfcfile = joinpath(pkgdir(@__MODULE__), "test", "parent_parabola.sfc")
    sfc = read_surface(sfcfile)

    t1 = sor_efficiency(cutfile; F=40.0, D=18.0, Oc=0.4, pol=:l3h, dz=0.0)

    tep_p1 = read_tepfile(joinpath(pkgdir(@__MODULE__), "test", "ticra_tools_twister.tep"))
    tfile = tempname()
    write_tepfile(tfile, tep_p1)
    d = 15u"mm"    
    tep_s2 = tepp2s(tep_p1, d)
    freqs = get_freqs(tep_p1)
    tep_p3 = teps2p(tep_s2, d, freqs)

    tep_scatter = read_tepfile(joinpath(pkgdir(@__MODULE__), "..", "test", "tepscatter1freq.tep"))
    tfile = tempname()
    write_tepfile(tfile, tep_scatter)
=#    
end

# include("make_cut_from_eh.jl")


end # module
