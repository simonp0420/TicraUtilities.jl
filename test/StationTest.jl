using SafeTestsets

@safetestset "Station" begin 
    using TicraUtilities

    stadat = only(read_stationfile("scenario2_coverage.sta"))
    @test get_npoint(stadat) == 397
    npoint = get_npoint(stadat)
    @test get_u(stadat)[npoint] == -0.00750
    @test get_v(stadat)[npoint] == -0.02000
    @test get_goal(stadat)[npoint] == 0.441500015E+02
    @test get_weight(stadat)[npoint] == -0.100000000E+01
    @test get_ipol(stadat)[npoint] == 4
    @test get_rot(stadat)[npoint] == 0.0
    @test get_att(stadat)[npoint] == 0.0
    @test get_id(stadat)[npoint] == ""

    tfile = tempname()
    write_stationfile(tfile, stadat)
    junk = only(read_stationfile(tfile))
    @test get_npoint(stadat) == get_npoint(junk)
    @test get_u(stadat) == get_u(junk)
    @test get_v(stadat) == get_v(junk)
    @test get_goal(stadat) == get_goal(junk)
    @test get_weight(stadat) == get_weight(junk)
    @test get_ipol(stadat) == get_ipol(junk)
    @test get_rot(stadat) == get_rot(junk)
    @test get_att(stadat) == get_att(junk)
    @test get_id(stadat) == get_id(junk)
end

