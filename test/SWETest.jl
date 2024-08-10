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
    sphfile = joinpath(@__DIR__, "center_element_rhcp_excited_q.swe")
    cut_test = sph2cut(sphfile; phi=0:5:355, theta=0:1:180, ipol=2)
    cut_good = read_ticra_cut(joinpath(@__DIR__, "center_element_rhcp_excited_q.cut"))
    Eθ_err =  maximum(norm, first.(cut_test.evec) - first.(cut_good.evec))
    Eϕ_err =  maximum(norm, last.(cut_test.evec) - last.(cut_good.evec))

    @test Eθ_err < 1e-9
    @test Eϕ_err < 1e-9
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

