using SafeTestsets
using Test

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

@safetestset "TEPscatter" begin 
    using TicraUtilities

    tep_scatter = read_tepfile("tepscatter1freq.tep")
    tfile = tempname()
    write_tepfile(tfile, tep_scatter)
    tep_scatter2 = read_tepfile(tfile)
    @test get_sff(tep_scatter) == get_sff(tep_scatter2)
    @test get_sfr(tep_scatter) == get_sfr(tep_scatter2)
    @test get_srr(tep_scatter) == get_srr(tep_scatter2)
    @test get_srf(tep_scatter) == get_srf(tep_scatter2)
    @test get_theta(tep_scatter) == get_theta(tep_scatter2)
    @test get_phi(tep_scatter) == get_phi(tep_scatter2)
    @test get_title(tep_scatter) == get_title(tep_scatter2)
end

using Unitful
@testset "tepp2s and teps2p" begin 
    using TicraUtilities
    tep_p1 = read_tepfile("ticra_tools_twister.tep")
    d = 15u"mm"    
    tep_s2 = tepp2s(tep_p1, d)
    freqs = get_freqs(tep_p1)
    tep_p3 = teps2p(tep_s2, d, freqs)
    @test get_sff(tep_p1) ≈ get_sff(tep_p3)
    @test get_srf(tep_p1) ≈ get_srf(tep_p3)
    @test get_srr(tep_p1) ≈ get_srr(tep_p3)
    @test get_sfr(tep_p1) ≈ get_sfr(tep_p3)
    tep_s4 = tepp2s(tep_p3, d)
    for i in eachindex(tep_s2, tep_s4)
        @test get_sff(tep_s2[i]) ≈ get_sff(tep_s4[i])
        @test get_srf(tep_s2[i]) ≈ get_srf(tep_s4[i])
        @test get_srr(tep_s2[i]) ≈ get_srr(tep_s4[i])
        @test get_sfr(tep_s2[i]) ≈ get_sfr(tep_s4[i])
    end
end
