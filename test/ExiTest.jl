using SafeTestsets

@safetestset "Exi" begin 
    using TicraUtilities
    t = read_exifile(joinpath(@__DIR__, "beam_A14R.exi"))
    @test length(get_header(t)) == 1
    @test get_header(t)[1] == "DirecTV-15S A14R"
    @test length(get_ids(t)) == length(amplitude_db(t)) == length(phase_deg(t)) == 16
    @test get_ampdb(t) == [-22.079632, -18.455905, -13.189809, -12.546182, -20.599867,
                      -13.271168, -8.728724, -7.783794, -11.498586, -10.453275,
                      -8.245438, -9.028067, -17.107155, -22.580741, -17.524788, -13.861297]
    @test get_phsdeg(t) == [-21.860822, 17.101443, -4.950916, -18.354877, -17.122747, 17.979759,
                        9.269317, -6.876699, -13.262062, 17.23366, 2.516957, -7.826971,
                        0.33419, 25.276105, 10.616932, -7.269771]
    @test get_ids(t) == ["elem0102", "elem0124", "elem0125", "elem0126", "elem0127", "elem0147",
                   "elem0148", "elem0149", "elem0150", "elem0173", "elem0174", "elem0175",
                   "elem0176", "elem0198", "elem0199", "elem0200"]
    @test complex(t)[[begin,end]] == [0.07304811253768777 - 0.02930714779532224im, 0.20110825515817327 - 0.025654723899872im]

    testfile = tempname()
    write_exifile(testfile, t)
    t2 = read_exifile(testfile)
    for field in fieldnames(typeof(t2))
        @test getfield(t, field) == getfield(t2, field)
    end
end

