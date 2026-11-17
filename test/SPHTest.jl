using SafeTestsets

@safetestset "read_write_sphfile" begin 
    using TicraUtilities
    q1s = read_sphfile(joinpath(@__DIR__, "tc4p506_champ3.sph"))
    tfile = joinpath(tempdir(), "temp.sph")
    write_sphfile(tfile, q1s)
    q2s = read_sphfile(tfile)

    @test length(q1s) == length(q2s)

    equal_test_fields = (:prgtag, :idstrg, :t4, :t5, :t6, :t7, :t8, :nthe, :nphi, :nmax, :mmax)
    array_fields =  (:qsmns, :powerms)
    for (q1, q2) in zip(q1s, q2s)
        for f in equal_test_fields
            @test getproperty(q1, f) == getproperty(q2, f)
        end
        for f in array_fields
            @test axes(getproperty(q1, f)) == axes(getproperty(q2, f))
            @test getproperty(q1, f) ≈ getproperty(q2, f)
        end
    end

end

@safetestset "sph2cut" begin 
    using TicraUtilities
    using LinearAlgebra: norm
    sphfile = joinpath(@__DIR__, "center_element_rhcp_excited_q.sph")
    cut_test = sph2cut(sphfile; phi=0:5:355, theta=0:1:180, ipol=2)
    cut_good = read_cutfile(joinpath(@__DIR__, "center_element_rhcp_excited_q.cut"))
    Eθ_err =  maximum(norm, first.(cut_test.evec) - first.(cut_good.evec))
    Eϕ_err =  maximum(norm, last.(cut_test.evec) - last.(cut_good.evec))

    @test Eθ_err < 1e-9
    @test Eϕ_err < 1e-9
end

@safetestset "cut2sph" begin 
    using TicraUtilities
    using Logging: with_logger, NullLogger
    cutfile = joinpath(@__DIR__, "issue20_cut.cut")
    
    @test with_logger(() -> cut2sph(cutfile), NullLogger()) isa SPHQPartition
end

@safetestset "Deltas" begin 
    using TicraUtilities
    @test TicraUtilities._Δⁿₙₘ(4,5) ≈ sqrt(10) / 32
    @test TicraUtilities._Δⁿₙₘ(1,1) ≈ 0.5
    @test TicraUtilities._Δⁿₙₘ(0,1) ≈ 1/√2
    @test TicraUtilities._Δⁿₙₘ(0,5) ≈ 3√7/16

    @test abs(TicraUtilities._Δⁿₘₚₘ(0,0,1)) < 1e-15
    @test TicraUtilities._Δⁿₘₚₘ(1,0,1) ≈ 1/√2
    @test TicraUtilities._Δⁿₘₚₘ(0,1,1) ≈ -1/√2
    @test TicraUtilities._Δⁿₘₚₘ(1,1,1) ≈ 0.5
    @test TicraUtilities._Δⁿₘₚₘ(0,0,4) ≈ 3/8
    @test abs(TicraUtilities._Δⁿₘₚₘ(0,3,4))  < 1e-15
    @test TicraUtilities._Δⁿₘₚₘ(0,4,4) ≈ √70/16
    @test abs(TicraUtilities._Δⁿₘₚₘ(1,0,4)) < 1e-15
    @test TicraUtilities._Δⁿₘₚₘ(2,0,4) ≈ -√10/8
    @test TicraUtilities._Δⁿₘₚₘ(2,1,4) ≈ -√2/8
    @test TicraUtilities._Δⁿₘₚₘ(2,2,4) ≈ 0.25
    @test TicraUtilities._Δⁿₘₚₘ(2,2,4) ≈ 0.25

    @test abs(TicraUtilities._Δⁿₘₚₘ(0,0,5)) < 1e-15
    @test TicraUtilities._Δⁿₘₚₘ(1,0,5) ≈ √30/16
    @test abs(TicraUtilities._Δⁿₘₚₘ(2,0,5)) < 1e-15
    @test TicraUtilities._Δⁿₘₚₘ(3,0,5) ≈ -√35/16
    @test abs(TicraUtilities._Δⁿₘₚₘ(0,4,5)) < 1e-15
    @test TicraUtilities._Δⁿₘₚₘ(1,4,5) ≈ √21/16
    @test TicraUtilities._Δⁿₘₚₘ(0,5,5) ≈ -3*√7/16
    @test TicraUtilities._Δⁿₘₚₘ(0,5,5) ≈ -3*√7/16
end

@safetestset "Cut Reconstruction" begin 
    using TicraUtilities: read_cutfile, cut2sph, sph2cut
    using LinearAlgebra: norm
    cutfile = joinpath(@__DIR__, "center_element_rhcp_excited.cut")
    cut1 = read_cutfile(cutfile)
    sph1 = cut2sph(cutfile)
    cut2 = sph2cut(sph1)
    maxerr = maximum(norm, cut1.evec - cut2.evec)
    @test maxerr < 0.0008

    cut1.evec = cut1.evec[1:91,:]
    cut1.theta = 0:90
    sph2 = cut2sph(cut1)
    maxerr2 = maximum(abs, sph2.qsmns - sph1.qsmns)
    @test maxerr2 < 1e-15
end


@safetestset "relative modal power budget" begin
    using TicraUtilities
    filter_modes = TicraUtilities._filter_qmodes_by_power
    q = TicraUtilities.OffsetArray(zeros(ComplexF64, 2, 3, 2), 1:2, -1:1, 1:2)
    q[1,0,1], q[1,1,1], q[1,0,2] = sqrt(2), sqrt(0.12), sqrt(0.12)
    reference_axes = nothing
    for scale in (1e-6, 1.0, 1e6)
        input = scale * q
        pwr, powers, kept = filter_modes(input, 0.1)
        total = 0.5sum(abs2, input)
        @test total - pwr <= 0.1total + 1e-14total
        @test pwr ≈ 0.5sum(abs2, kept)
        @test pwr ≈ sum(powers)
        if isnothing(reference_axes)
            reference_axes = axes(kept)
        else
            @test axes(kept) == reference_axes
        end
    end
    @test filter_modes(q, 0.0)[3] === q
    @test filter_modes(q, -1.0)[3] === q
end

@safetestset "axisymmetric spherical modes and truncated three-component cuts" begin
    using TicraUtilities
    q = TicraUtilities.OffsetArray(zeros(ComplexF64, 2, 1, 1), 1:2, 0:0, 1:1)
    q[2,0,1] = 1
    sph = SPHQPartition(; prgtag="test", idstrg="m=0", nthe=12, nphi=8, nmax=1, mmax=0, qsmns=q)
    cut = sph2cut(sph; theta=0:30:180, phi=0:45:315, ipol=1)
    @test all(iszero, cut.evec[[1,end],:])
    @test all(cut.evec[:,j] ≈ cut.evec[:,1] for j in 2:8)
    @test abs(cut.evec[4,1][1]) ≈ sqrt(3/(16π))
    for ipol in (-1,4)
        @test_throws ArgumentError sph2cut(sph; ipol)
    end
    ev = [TicraUtilities.SVector{3,ComplexF64}(e[1], e[2], 123) for e in cut.evec[1:4,:]]
    partial = Cut(; ncomp=3, icut=1, icomp=1, text=cut.text, theta=0:30:90, phi=cut.phi, evec=ev)
    result = @test_logs (:warn, r"Discarding 3rd slot") cut2sph(partial; nmax=2, mmax=0)
    @test result.nmax == 2
    @test result.mmax == 0
end

@safetestset "nmax limits output modes without changing lower coefficients" begin
    using TicraUtilities
    theta, phi = 0:15:180, 0:45:315
    ev = [TicraUtilities.SVector{2,ComplexF64}(sind(t)*cosd(t)^2,0) for t in theta,p in phi]
    cut = Cut(; ncomp=2, icut=1, icomp=1, text=fill("",8), theta, phi, evec=ev)
    full = cut2sph(cut; mmax=0)
    for nmax in (1,2,3,4,6)
        limited = cut2sph(cut; nmax, mmax=0)
        @test limited.nmax == nmax
        for n in 1:nmax, s in 1:2
            @test limited.qsmns[s,0,n] ≈ full.qsmns[s,0,n] atol=1e-14
        end
    end
    @test cut2sph(cut; nmax=1).mmax <= 1
end
