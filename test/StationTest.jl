using SafeTestsets

@safetestset "Station" begin 
    using TicraUtilities

    stadat = read_stationfile("scenario2_coverage.sta")
    @test  length(stadat) == 1
    @test stadat[1].npoint == 397
    npoint = stadat[1].npoint
    @test stadat[1].u[npoint] == -0.00750
    @test stadat[1].v[npoint] == -0.02000
    @test stadat[1].goal[npoint] == 0.441500015E+02
    @test stadat[1].weight[npoint] == -0.100000000E+01
    @test stadat[1].ipol[npoint] == 4
    @test stadat[1].rot[npoint] == 0.0
    @test stadat[1].att[npoint] == 0.0
    @test stadat[1].ID[npoint] == ""

    tfile = tempname()
    write_stationfile(tfile, stadat)
    junk = read_stationfile(tfile)
    @test length(junk) == 1
    @test stadat[1].npoint == junk[1].npoint
    @test stadat[1].u == junk[1].u
    @test stadat[1].v == junk[1].v
    @test stadat[1].goal == junk[1].goal
    @test stadat[1].weight == junk[1].weight
    @test stadat[1].ipol == junk[1].ipol
    @test stadat[1].rot == junk[1].rot
    @test stadat[1].att == junk[1].att
    @test stadat[1].ID == junk[1].ID
end

