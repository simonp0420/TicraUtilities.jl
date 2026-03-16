using SafeTestsets

@safetestset "read_ffdfile" begin 
    using TicraUtilities
    t1fi = read_ffdfile(joinpath(@__DIR__, "ffd", "dipole_findep.ffd"))
    t1fd = read_ffdfile(joinpath(@__DIR__, "ffd", "dipole_fdep_1freq.ffd"))
    t3 = read_ffdfile(joinpath(@__DIR__, "ffd", "dipole_fdep_3freqs.ffd"))
    @test length(t3) == 3

    @test t3[1].theta == t3[2].theta == t3[3].theta == t1fi.theta == t1fd.theta == 0:2:180
    @test t3[1].phi == t3[2].phi == t3[3].phi == t1fi.phi == t1fd.phi == -180:2:180

    @test t3[1].evec == t3[2].evec == t3[3].evec == t1fi.evec == t1fd.evec
    @test t1fi.evec[begin][1] == 0.000132235762166748 + 0.0006963056517931414im
    @test t1fi.evec[begin][2] == 0.002024889625826101 - 0.0004704676594853991im
    @test t1fd.evec[end][1] == -0.0006904994350891429 - 0.0007043391827337916im
    @test t1fd.evec[end][2] == -0.003426592661168071 - 7.536174210720495e-5im
    @test t1fd.evec[25, 60][1] == -0.09103464550148568 - 0.6764057884083713im
    @test t1fd.evec[25, 60][2] == 0.0001648293857863699 - 0.0008244348013471126im

    @test t1fi.frequency == 0
    @test t1fd.frequency == 1e10
    @test (t3[1].frequency, t3[2].frequency, t3[3].frequency) == (1e10, 2e10, 3e10)

end

