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

