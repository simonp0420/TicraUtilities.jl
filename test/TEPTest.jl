using SafeTestsets

@safetestset "TEPperiodic" begin 
    using TicraUtilities

    tep_periodic = read_tepfile("ticra_tools_twister.tep")
    tfile = tempname()
    write_tepfile(tfile, tep_periodic)
    tep_periodic2 = read_tepfile(tfile)
    @test get_sff(tep_periodic) == get_sff(tep_periodic2)
    @test get_sfr(tep_periodic) == get_sfr(tep_periodic2)
    @test get_srr(tep_periodic) == get_srr(tep_periodic2)
    @test get_srf(tep_periodic) == get_srf(tep_periodic2)
    @test get_theta(tep_periodic) == get_theta(tep_periodic2)
    @test get_phi(tep_periodic) == get_phi(tep_periodic2)
    @test get_freqs(tep_periodic) == get_freqs(tep_periodic2)
    @test get_name(tep_periodic) == get_name(tep_periodic2)
    @test get_class(tep_periodic) == get_class(tep_periodic2)

end

