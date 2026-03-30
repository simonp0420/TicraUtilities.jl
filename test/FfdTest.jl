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

@safetestset "cut2ffd" begin
    using TicraUtilities

    # Single Cut to frequency-independent Ffd
    cut = read_cutfile(joinpath(@__DIR__, "test.cut"))
    ffd = cut2ffd(cut)
    @test ffd.frequency == 0.0
    @test ffd.theta == cut.theta
    @test ffd.phi == cut.phi
    @test size(ffd.evec) == size(cut.evec)

    # Single Cut to frequency-dependent Ffd
    ffd_freq = cut2ffd(cut; frequency=1e9)
    @test ffd_freq.frequency == 1e9
    @test ffd_freq.theta == cut.theta
    @test ffd_freq.phi == cut.phi
    @test size(ffd_freq.evec) == size(cut.evec)

    # Vector of Cuts to vector of Ffds (multifrequency)
    cuts = [cut, cut, cut]  # Use same cut for simplicity
    frequencies = [1e9, 2e9, 3e9]
    ffds = cut2ffd(cuts; frequency=frequencies)
    @test length(ffds) == length(cuts)
    for (ffd, f) in zip(ffds, frequencies)
        @test ffd.frequency == f
        @test ffd.theta == cut.theta
        @test ffd.phi == cut.phi
        @test size(ffd.evec) == size(cut.evec)
    end

    # File-based tests

    # Create freq-dependent Ffd
    cutfile = joinpath(@__DIR__, "test.cut")
    ffdfile = joinpath(tempdir(), "temp.ffd")
    ffds_file = cut2ffd(cutfile, ffdfile; frequency=1e9)
    ffd_read = read_ffdfile(ffdfile)
    @test ffd_read.frequency == 1e9
    @test ffd_read.theta == cut.theta
    @test ffd_read.phi == cut.phi
    @test size(ffd_read.evec) == size(cut.evec)

    # Create freq-independent Ffd
    cutfile = joinpath(@__DIR__, "test.cut")
    ffdfile = joinpath(tempdir(), "temp.ffd")
    ffds_file = cut2ffd(cutfile, ffdfile)
    for (lcount, line) in enumerate(eachline(ffdfile))
        if lcount == 9
            @test !occursin("Frequency", line)
            break
        end
    end
    ffd_read = read_ffdfile(ffdfile)
    @test ffd_read.frequency == 0.0
    @test ffd_read.theta == cut.theta
    @test ffd_read.phi == cut.phi
    @test size(ffd_read.evec) == size(cut.evec)

    # Multifrequency file - create a temp cut file with multiple cuts
    cuts_temp = [cut, cut, cut]
    cutfile_multi = joinpath(tempdir(), "temp_multi.cut")
    write_cutfile(cutfile_multi, cuts_temp)
    ffdfile_multi = joinpath(tempdir(), "temp_multi.ffd")
    @test_throws TypeError cut2ffd(cutfile_multi, ffdfile_multi)
    ffds_file_multi = cut2ffd(cutfile_multi, ffdfile_multi; frequency=frequencies)
    ffds_read_multi = read_ffdfile(ffdfile_multi)
    @test length(ffds_read_multi) == length(ffds_file_multi)
    for (ffd, f) in zip(ffds_read_multi, frequencies)
        @test ffd.frequency == f
    end
end

@safetestset "ffd2sph" begin
    using TicraUtilities

    # Single Ffd to SPH
    ffd = read_ffdfile(joinpath(@__DIR__, "ffd", "dipole_findep.ffd"))
    sph = ffd2sph(ffd)
    @test sph isa SPHQPartition
    @test sph.frequency == 0.0

    # Multifrequency Ffds to vector of SPHs
    ffds = read_ffdfile(joinpath(@__DIR__, "ffd", "dipole_fdep_3freqs.ffd"))
    sphs = ffd2sph(ffds)
    @test length(sphs) == length(ffds)
    for (sph, ffd) in zip(sphs, ffds)
        @test sph.frequency == ffd.frequency
    end

    # File-based tests with style :ticra
    ffdfile = joinpath(@__DIR__, "ffd", "dipole_findep.ffd")
    sphfile_ticra = joinpath(tempdir(), "temp_ticra.sph")
    sph_file = ffd2sph(ffdfile, sphfile_ticra; style=:ticra)
    open(sphfile_ticra, "r") do fid
        foreach((i) -> readline(fid), 1:8)
        line9 = readline(fid)
        @test !occursin("Frequency", line9)
    end
    sph_read = read_sphfile(sphfile_ticra)
    @test sph_read.frequency == 0.0

    # File-based tests with style :hfss (frequency-dependent)
    ffd_freq = read_ffdfile(joinpath(@__DIR__, "ffd", "dipole_fdep_1freq.ffd"))
    ffd_temp_file = joinpath(tempdir(), "temp_freq.ffd")
    write_ffdfile(ffd_temp_file, ffd_freq)
    sphfile_hfss = joinpath(tempdir(), "temp_hfss.swef")
    sph_file_hfss = ffd2sph(ffd_temp_file, sphfile_hfss; style=:hfss)
    open(sphfile_hfss, "r") do fid
        foreach((i) -> readline(fid), 1:8)
        line9 = readline(fid)
        @test occursin("Frequency", line9)
    end
    sph_read_hfss = read_sphfile(sphfile_hfss)
    @test sph_read_hfss.frequency == ffd_freq.frequency
end

@testset "sph2ffd" begin
    using TicraUtilities

    # Single SPH to Ffd
    sph = read_sphfile(joinpath(@__DIR__, "tc4p506_champ3.sph")) |> first
    ffd = sph2ffd(sph)
    @test ffd isa Ffd
    @test ffd.frequency == sph.frequency

    # SPH with custom theta/phi
    ffd_custom = sph2ffd(sph; theta=0:2:180, phi=0:5:355)
    @test ffd_custom.theta == 0:2:180
    @test ffd_custom.phi == 0:5:355

    # Vector of SPHs to vector of Ffds
    # For multifrequency, we can create from the 3freq ffd, convert to sph, then back
    ffds_orig = read_ffdfile(joinpath(@__DIR__, "ffd", "dipole_fdep_3freqs.ffd"))
    sphs_multi = ffd2sph(ffds_orig)
    ffds_back = sph2ffd(sphs_multi)
    @test length(ffds_back) == length(ffds_orig)
    for (ffd_back, ffd_orig) in zip(ffds_back, ffds_orig)
        @test ffd_back.frequency == ffd_orig.frequency
        @test ffd_back.theta == ffd_orig.theta
        @test length(ffd_back.phi) + 1 == length(ffd_orig.phi)
    end

    # File-based tests with Ticra style
    sphfile_ticra = joinpath(@__DIR__, "tc4p506_champ3.sph")
    ffds_from_file = sph2ffd(sphfile_ticra; frequency = (1:10) * 1e9)
    @test [ffd.frequency for ffd in ffds_from_file] == (1:10) * 1e9

    # For HFSS style, create a temp file
    sph_temp_hfss = joinpath(tempdir(), "temp_hfss.swef")
    @reset sph.frequency = 1e9
    write_sphfile(sph_temp_hfss, sph; style=:hfss)
    ffd_hfss = sph2ffd(sph_temp_hfss)
    @test ffd_hfss.frequency == 1e9

    # Test with frequency override
    ffd_override = sph2ffd(sph; frequency=2e9)
    @test ffd_override.frequency == 2e9
end

@safetestset "additional_hfss_interop" begin
    using TicraUtilities
    using LinearAlgebra: norm

    # Test get_frequency on vector of Ffds
    ffds = read_ffdfile(joinpath(@__DIR__, "ffd", "dipole_fdep_3freqs.ffd"))
    ffd_freqs = get_frequency.(ffds)
    @test ffd_freqs == [1e10, 2e10, 3e10]

    # Test get_frequency on single Ffd
    ffd_single = read_ffdfile(joinpath(@__DIR__, "ffd", "dipole_findep.ffd"))
    @test get_frequency(ffd_single) == 0.0

    # Test reading SWEF file
    swef_file = joinpath(@__DIR__, "ffd", "center_element_rhcp_excited_q.swef")
    swef = read_sphfile(swef_file)
    @test swef isa SPHQPartition
    @test get_frequency(swef) > 0.0  # Should have frequency

    # Test writing SPH file with style=:hfss and frequency
    sph = read_sphfile(joinpath(@__DIR__, "tc4p506_champ3.sph")) |> first
    new_swef_file = joinpath(tempdir(), "temp_hfss_freq.swef")
    write_sphfile(new_swef_file, sph, style=:hfss, frequency=1e9)
    # Check that line 9 contains frequency
    open(new_swef_file, "r") do fid
        foreach((i) -> readline(fid), 1:8)
        line9 = readline(fid)
        @test occursin("Frequency", line9)
    end
    swef_read = read_sphfile(new_swef_file)
    @test get_frequency(swef_read) == 1e9

    # Test ffd2sph file-to-file with style=:hfss
    ffds_file = joinpath(@__DIR__, "ffd", "dipole_fdep_3freqs.ffd")
    swef_file_out = joinpath(tempdir(), "temp_multi.swef")
    ffd2sph(ffds_file, swef_file_out, style=:hfss)
    # Check that it has frequency lines
    open(swef_file_out, "r") do fid
        lines = readlines(fid)
        @test any(occursin("Frequency", line) for line in lines)
    end

    # Test sph2ffd file-to-file
    new_ffds_file = joinpath(tempdir(), "temp_from_swef.ffd")
    sph2ffd(swef_file_out, new_ffds_file)
    ffds_from_swef = read_ffdfile(new_ffds_file)
    @test length(ffds_from_swef) == 3
    @test all(ffd -> get_frequency(ffd) > 0, ffds_from_swef)

    # Test maximum norm difference between orig and reconstructed Ffd objects using get_evec
    ffd_orig = read_ffdfile(joinpath(@__DIR__, "ffd", "dipole_findep.ffd"))
    sph_from_ffd = ffd2sph(ffd_orig)
    ffd_reconstructed = sph2ffd(sph_from_ffd, phi=get_phi(ffd_orig))
    max_diff = maximum(norm, get_evec(ffd_reconstructed) - get_evec(ffd_orig))
    @test max_diff < 1e-3  # Should be small due to default pwrtol
    sph_from_ffd = ffd2sph(ffd_orig, pwrtol=0)
    ffd_reconstructed = sph2ffd(sph_from_ffd, phi=get_phi(ffd_orig))
    max_diff = maximum(norm, get_evec(ffd_reconstructed) - get_evec(ffd_orig))
    @test max_diff < 1e-14  # Should be small due to default pwrtol

end

