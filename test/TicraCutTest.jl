using SafeTestsets

@safetestset "read_ticra_cuts" begin 
    using TicraUtilities
    t2 = read_ticra_cuts(joinpath(@__DIR__, "test2.cut"))
    @test length(t2) == 3
    @test t2[1].theta == t2[2].theta == t2[3].theta == 0.0:0.5:180.0
    @test t2[1].phi == t2[2].phi == t2[3].phi == 0.0:5.0:355.0
    @test t2[1].p1[1,1] ≈ 0.16835983292 + 0.0092799907903im
    @test t2[1].p2[1,1] ≈ 1.5392084724 + 3.7955862332im
    @test t2[1].p1[end,1] ≈ -0.12828987268 + 0.0010199989877im
    @test t2[1].p2[end,1] ≈ -0.0011099988984 + 0.0064499935988im
    @test t2[3].p1[1,1] ≈ 0.16835983292 + 0.0092799907903im
    @test t2[3].p2[1,1] ≈ 1.5392084724 + 3.7955862332im
    @test t2[3].p1[end,end] ≈ -0.1261598748 + 0.023279976896im
    @test t2[3].p2[end,end] ≈ -0.0022099978067 + 0.0061599938867im
    @test t2[1].icomp == 2
    @test t2[2].icomp == 3
    @test t2[3].icomp == 4
    @test t2[1].ncomp == t2[2].ncomp == t2[3].ncomp == 2
    @test t2[1].icut == t2[2].icut == t2[3].icut == 1

end

@safetestset "read_ticra_cut" begin 
    using TicraUtilities
    t = read_ticra_cut(joinpath(@__DIR__, "test.cut"))
    @test t.theta == 0.0:0.5:180.0
    @test t.phi == 0.0:5.0:355.0
    @test t.p1[1,1] ≈ 0.16835983292 + 0.0092799907903im
    @test t.p2[1,1] ≈ 1.5392084724 + 3.7955862332im
    @test t.p1[end,1] ≈ -0.12828987268 + 0.0010199989877im
    @test t.p2[end,1] ≈ -0.0011099988984 + 0.0064499935988im
    @test t.p1[end,end] ≈ -0.1261598748 + 0.023279976896im
    @test t.p2[end,end] ≈ -0.0022099978067 + 0.0061599938867im
end

@safetestset "TicraCut amplitude_db" begin 
    using TicraUtilities
    t = read_ticra_cut(joinpath(@__DIR__, "test.cut"))
    @test 10*log10(abs2(t.p1[10,10])) == amplitude_db(t,1)[10,10]
    @test 10*log10(abs2(t.p2[10,10])) == amplitude_db(t,2)[10,10]
    @test 10*log10(abs2(t.p2[10,10])) == amplitude_db(t,"copol")[10,10]
    @test 10*log10(abs2(t.p1[10,10])) == amplitude_db(t,"crosspol")[10,10]
end

@safetestset "TicraCut phase_deg" begin 
    using TicraUtilities
    t = read_ticra_cut(joinpath(@__DIR__, "test.cut"))
    @test rad2deg(angle(t.p1[10,10])) == phase_deg(t,1)[10,10]
    @test rad2deg(angle(t.p2[10,10])) == phase_deg(t,2)[10,10]
    @test rad2deg(angle(t.p2[10,10])) == phase_deg(t,"copol")[10,10]
    @test rad2deg(angle(t.p1[10,10])) == phase_deg(t,"crosspol")[10,10]
end

@safetestset "TicraCut power" begin 
    using TicraUtilities
    t = read_ticra_cut(joinpath(@__DIR__, "test.cut"))
    @test abs(power(t) - 4π) ≤ 1.e-5
    ts = read_ticra_cuts(joinpath(@__DIR__, "test2.cut"))
    for t1 in ts
        @test abs(power(t1) - 4π) ≤ 1.e-5
    end
end

@safetestset "write_ticra_cut" begin 
    using TicraUtilities
    cut1 = read_ticra_cut(joinpath(@__DIR__, "test.cut"))
    tfile = joinpath(tempdir(), "temp.cut")
    write_ticra_cut(tfile, cut1, "Cut file normalized to directivity")
    cut2 = read_ticra_cut(tfile)
    for fn in fieldnames(TicraUtilities.TicraCut)
        @test getfield(cut1, fn) == getfield(cut2, fn)
    end

    cuts1 = read_ticra_cuts(joinpath(@__DIR__, "test2.cut"))
    tfile = joinpath(tempdir(),"temp.cut")
    write_ticra_cut(tfile, cuts1, "Cut file normalized to directivity")
    cuts2 = read_ticra_cuts(tfile)
    for (i,(cut1, cut2)) in enumerate(zip(cuts1, cuts2))
        for fn in fieldnames(TicraUtilities.TicraCut)
            if getfield(cut1, fn) ≠ getfield(cut2, fn)
                @show i
                @show fn
                @show (getfield(cut1, fn) , getfield(cut2, fn))
            end
            @test getfield(cut1, fn) == getfield(cut2, fn)
        end
    end
end

@safetestset "TicraCut phscen" begin 
    using TicraUtilities
    (x,y,z0,z90) = phscen(joinpath(@__DIR__, "test.cut"))
    @test x ≈ -6.740495718988637e-5
    @test y ≈ 0.0002328819050067705
    @test z0 ≈ 0.07914423679589953
    @test z90 ≈ 0.06397213096541265
end


@safetestset "eval_cut" begin 
    using TicraUtilities
    (c, xn, sp, slh, et, pc, xpd) = eval_cut(joinpath(@__DIR__, "test.cut"), 11.80285, 15.0)
    @test c ≈ 12.246790349581193
    @test length(xn) == 2232
    @test all(xn[1:10] .≈ [27.70884584330901, 27.719167304847517, 27.77196806190906, 27.867909509973337, 28.004986266580865,
                       28.181480124837947, 28.395836074712978, 28.644739904741375, 28.92386281246011, 29.227489414900464])
    @test sp ≈ 6.122080463458496
    @test slh == 0.0
    @test length(et) == 72
    @test all((et[begin],et[end]) .≈ (1.4331621843858997, 1.4238487987626876))
    @test pc ≈ -0.010691508997658557
    @test length(xpd) == 2232
    @test all((xpd[begin], xpd[end], sum(xpd)/length(xpd)) .≈ 
              (27.70884584330901, 26.590458769432246, 30.348100641028992))
end
