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

@safetestset "write_ffdfile" begin
    using TicraUtilities

    # Single-frequency file
    ffd1 = read_ffdfile(joinpath(@__DIR__, "ffd", "dipole_findep.ffd"))
    tfile = joinpath(tempdir(), "temp.ffd")
    write_ffdfile(tfile, ffd1)
    ffd1b = read_ffdfile(tfile)
    @test ffd1b.frequency == ffd1.frequency
    @test ffd1b.theta == ffd1.theta
    @test ffd1b.phi == ffd1.phi
    @test ffd1b.evec ≈ ffd1.evec

    # Multi-frequency file (vector input)
    ffds = read_ffdfile(joinpath(@__DIR__, "ffd", "dipole_fdep_3freqs.ffd"))
    tfile = joinpath(tempdir(), "temp2.ffd")
    write_ffdfile(tfile, ffds)
    ffds_b = read_ffdfile(tfile)
    @test length(ffds_b) == length(ffds)
    for (ffd, ffd_b) in zip(ffds, ffds_b)
        @test ffd_b.frequency == ffd.frequency
        @test ffd_b.theta == ffd.theta
        @test ffd_b.phi == ffd.phi
        @test ffd_b.evec ≈ ffd.evec
    end

    # Writing a vector of length 1 should round-trip as a single Ffd
    tfile = joinpath(tempdir(), "temp3.ffd")
    write_ffdfile(tfile, [ffd1])
    ffd1c = read_ffdfile(tfile)
    @test ffd1c isa Ffd
    @test ffd1c.frequency == ffd1.frequency
    @test ffd1c.theta == ffd1.theta
    @test ffd1c.phi == ffd1.phi
    @test ffd1c.evec ≈ ffd1.evec
end

@safetestset "ffd2cut" begin
    using TicraUtilities

    # Convert from file and write cut file
    ffdfile = joinpath(@__DIR__, "ffd", "dipole_findep.ffd")
    cutfile = joinpath(tempdir(), "temp.cut")
    cut = ffd2cut(ffdfile, cutfile)
    cut2 = read_cutfile(cutfile)
    @test isapprox(cut, cut2)

    # Default behaviour should reorder phi cuts to start at 0
    @test cut.phi[1] == 0.0
    @test cut.phi[end] == 358.0
    @test length(cut.phi) == 180

    # phi0=false should preserve original phi ordering
    ffd = read_ffdfile(ffdfile)
    cut_no_phi0 = ffd2cut(ffd; phi0=false)
    @test cut_no_phi0.phi == ffd.phi

    # Converting a vector of Ffd objects should return a vector of Cut objects
    ffds = read_ffdfile(joinpath(@__DIR__, "ffd", "dipole_fdep_3freqs.ffd"))
    cuts = ffd2cut(ffds; phi0=false)
    @test length(cuts) == length(ffds)
    for (ffd, cut) in zip(ffds, cuts)
        @test occursin(string(ffd.frequency), cut.text[1])
    end

    # Calling with same input/output file should throw
    @test_throws ArgumentError ffd2cut(ffdfile, ffdfile)
end

