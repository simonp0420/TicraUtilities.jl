using SafeTestsets

@safetestset "read_cutfile" begin 
    using TicraUtilities
    t2 = read_cutfile(joinpath(@__DIR__, "test2.cut"))
    @test length(t2) == 3
    @test t2[1].theta == t2[2].theta == t2[3].theta == 0.0:0.5:180.0
    @test t2[1].phi == t2[2].phi == t2[3].phi == 0.0:5.0:355.0
    @test t2[1].evec[1,1][1] ≈ 0.16835983292 + 0.0092799907903im
    @test t2[1].evec[1,1][2] ≈ 1.5392084724 + 3.7955862332im
    @test t2[1].evec[end,1][1] ≈ -0.12828987268 + 0.0010199989877im
    @test t2[1].evec[end,1][2] ≈ -0.0011099988984 + 0.0064499935988im
    @test t2[3].evec[1,1][1] ≈ 0.16835983292 + 0.0092799907903im
    @test t2[3].evec[1,1][2] ≈ 1.5392084724 + 3.7955862332im
    @test t2[3].evec[end,end][1] ≈ -0.1261598748 + 0.023279976896im
    @test t2[3].evec[end,end][2] ≈ -0.0022099978067 + 0.0061599938867im
    @test t2[1].icomp == 2
    @test t2[2].icomp == 3
    @test t2[3].icomp == 4
    @test t2[1].ncomp == t2[2].ncomp == t2[3].ncomp == 2
    @test t2[1].icut == t2[2].icut == t2[3].icut == 1

    @test read_cutfile(joinpath(@__DIR__, "single_cut.cut")) isa Cut # Issue #17
    @test read_cutfile(joinpath(@__DIR__, "issue20_cut.cut")) isa Cut # Issue #20
end

@safetestset "read_cutfile_2" begin 
    using TicraUtilities
    t = read_cutfile(joinpath(@__DIR__, "test.cut"))
    @test t.theta == 0.0:0.5:180.0
    @test t.phi == 0.0:5.0:355.0
    @test t.evec[1,1][1] ≈ 0.16835983292 + 0.0092799907903im
    @test t.evec[1,1][2] ≈ 1.5392084724 + 3.7955862332im
    @test t.evec[end,1][1] ≈ -0.12828987268 + 0.0010199989877im
    @test t.evec[end,1][2] ≈ -0.0011099988984 + 0.0064499935988im
    @test t.evec[end,end][1] ≈ -0.1261598748 + 0.023279976896im
    @test t.evec[end,end][2] ≈ -0.0022099978067 + 0.0061599938867im
end


@safetestset "sym_asym_sym" begin 
    using TicraUtilities
    using Logging: with_logger, NullLogger
    t1 = read_cutfile(joinpath(@__DIR__, "test.cut"))
    t2 = asym2sym(t1)
    t3 = sym2asym(t2)
    @test t3 ≈ t1

    t4 = read_cutfile(joinpath(@__DIR__, "single_cut.cut"))
    @test sym2asym(t4) isa Cut

    t5 = read_cutfile(joinpath(@__DIR__, "issue20_cut.cut"))
    at5 = with_logger(() -> sym2asym(t5), NullLogger()) 
    st5 = asym2sym(at5)
    at6 = sym2asym(st5)
    st6 = asym2sym(at6)
    @test at6 ≈ at5
    @test st6 ≈ st5

end

@safetestset "Cut amplitude_db" begin 
    using TicraUtilities
    t = read_cutfile(joinpath(@__DIR__, "test.cut"))
    @test 10*log10(abs2(t.evec[10,10][1])) == amplitude_db(t,1)[10,10]
    @test 10*log10(abs2(t.evec[10,10][2])) == amplitude_db(t,2)[10,10]
    @test 10*log10(abs2(t.evec[10,10][2])) == amplitude_db(t,:copol)[10,10]
    @test 10*log10(abs2(t.evec[10,10][1])) == amplitude_db(t,:xpol)[10,10]
    @test 10*log10(abs2(t.evec[10,10][2])) == amplitude_db(t,"copol")[10,10]
    @test 10*log10(abs2(t.evec[10,10][1])) == amplitude_db(t,"xpol")[10,10]
end

@safetestset "Cut phase_deg" begin 
    using TicraUtilities
    t = read_cutfile(joinpath(@__DIR__, "test.cut"))
    @test rad2deg(angle(t.evec[10,10][1])) == phase_deg(t,1)[10,10]
    @test rad2deg(angle(t.evec[10,10][2])) == phase_deg(t,2)[10,10]
    @test rad2deg(angle(t.evec[10,10][2])) == phase_deg(t,:copol)[10,10]
    @test rad2deg(angle(t.evec[10,10][1])) == phase_deg(t,:xpol)[10,10]
    @test rad2deg(angle(t.evec[10,10][2])) == phase_deg(t,"copol")[10,10]
    @test rad2deg(angle(t.evec[10,10][1])) == phase_deg(t,"xpol")[10,10]
end

@safetestset "Cut power" begin 
    using TicraUtilities
    t = read_cutfile(joinpath(@__DIR__, "test.cut"))
    @test abs(power(t) - 4π) ≤ 1.e-5
    @test maximum_db(t) ≈ 12.246790349581193
    ts = read_cutfile(joinpath(@__DIR__, "test2.cut"))
    for t1 in ts
        @test abs(power(t1) - 4π) ≤ 1.e-5
    end
end

@safetestset "write_cutfile" begin 
    using TicraUtilities

    cut1 = read_cutfile(joinpath(@__DIR__, "test.cut"))
    tfile = joinpath(tempdir(), "temp.cut")
    write_cutfile(tfile, cut1, "Cut file normalized to directivity")
    cut2 = read_cutfile(tfile)
    for fn in fieldnames(TicraUtilities.Cut)
        @test getfield(cut1, fn) == getfield(cut2, fn)
    end

    cut1_3components = TicraUtilities._add_3rd_component(cut1)
    tfile = joinpath(tempdir(), "temp.cut")
    write_cutfile(tfile, cut1_3components, "Cut file normalized to directivity")
    cut2_3components = read_cutfile(tfile)
    for fn in fieldnames(TicraUtilities.Cut)
        if fn == :evec
            @test getfield(cut1_3components, fn) ≈ getfield(cut2_3components, fn)
        else
            @test getfield(cut1_3components, fn) == getfield(cut2_3components, fn)
        end
    end

    cuts1 = read_cutfile(joinpath(@__DIR__, "test2.cut"))
    tfile = joinpath(tempdir(),"temp.cut")
    write_cutfile(tfile, cuts1, "Cut file normalized to directivity")
    cuts2 = read_cutfile(tfile)
    for (i,(cut1, cut2)) in enumerate(zip(cuts1, cuts2))
        for fn in fieldnames(TicraUtilities.Cut)
            if getfield(cut1, fn) ≠ getfield(cut2, fn)
                @show i
                @show fn
                @show (getfield(cut1, fn) , getfield(cut2, fn))
            end
            @test getfield(cut1, fn) == getfield(cut2, fn)
        end
    end

    cuts1_3components = [TicraUtilities._add_3rd_component(cut) for cut in cuts1]
    tfile = joinpath(tempdir(),"temp.cut")
    write_cutfile(tfile, cuts1_3components, "Cut file normalized to directivity")
    cuts2_3components = read_cutfile(tfile)
    for (i,(cut1, cut2)) in enumerate(zip(cuts1_3components, cuts2_3components))
        for fn in fieldnames(TicraUtilities.Cut)
            if fn == :evec
                @test getfield(cut1, fn) ≈ getfield(cut2, fn)
            else
                if getfield(cut1, fn) ≠ getfield(cut2, fn)
                    @show i
                    @show fn
                    @show (getfield(cut1, fn) , getfield(cut2, fn))
                end
                @test getfield(cut1, fn) == getfield(cut2, fn)
            end            
        end
    end


end

@safetestset "Cut phscen" begin 
    using TicraUtilities
    (x,y,z0,z90) = phscen(joinpath(@__DIR__, "test.cut"))
    @test x ≈ -6.740495718988637e-5
    @test y ≈ 0.0002328819050067705
    @test z0 ≈ 0.07932537034841126
    @test z90 ≈ 0.06417366223154171
end


@safetestset "eval_cut" begin 
    using TicraUtilities
    (c, xn, sp, slh, et, pc, xpd) = TicraUtilities.eval_cut(joinpath(@__DIR__, "test.cut"), 11.80285, 15.0)
    @test c ≈ 12.246790349581193
    @test length(xn) == 2232
    @test all(xn[1:10] .≈ [27.70884584330901, 27.719167304847517, 27.77196806190906, 27.867909509973337, 28.004986266580865,
                       28.181480124837947, 28.395836074712978, 28.644739904741375, 28.92386281246011, 29.227489414900464])
    @test sp ≈ 6.122090227848759
    @test slh == 0.0
    @test length(et) == 72
    @test all((et[begin],et[end]) .≈ (1.4331621843858997, 1.4238487987626876))
    @test pc ≈ -0.011224035005254658
    @test length(xpd) == 2232
    @test all((xpd[begin], xpd[end], sum(xpd)/length(xpd)) .≈ 
              (27.70884584330901, 26.590458769432246, 30.348100641028992))
end

@safetestset "convert_cut" begin 
    using TicraUtilities
    cutfile = joinpath(@__DIR__, "test.cut")
    cut2 = read_cutfile(cutfile)
    cut12 = convert_cut(cut2, 1)
    cut21 = convert_cut(cut12, 2)
    cut32 = convert_cut(cut2, 3)
    cut31 = deepcopy(cut12); convert_cut!(cut31, 3)
    cut13 = convert_cut(cut32, 1)
    cut23 = convert_cut(cut31, 2)
    @test cut32.evec ≈ cut31.evec
    @test cut13.evec ≈ cut12.evec
    @test cut23.evec ≈ cut2.evec
    @test cut21.evec ≈ cut2.evec
end


@safetestset "eh2bor1cut" begin 
    using TicraUtilities
    cutfile = joinpath(@__DIR__, "ticra_hpol_horn.cut")
    cuteh = read_cutfile(cutfile)
    fe = get_evec(cuteh, 1)[:, 1] # ϕ = 0° cut is E-plane for horizontal pol
    fh = get_evec(cuteh, 1)[:, end] # ϕ = 90° cut is H-plane for horizontal pol
    cutrhcp = eh2bor1cut(get_theta(cuteh), fe, fh; pol=:rhcp, xpd=-30)
    copol_amp_db = amplitude_db(cutrhcp, 1)[1,1]
    xpol_amp_db = amplitude_db(cutrhcp, 2)[1,1]
    @test copol_amp_db - xpol_amp_db ≈ 30.0

    fe_func(θ) = fe[findfirst(x -> abs(x-θ) ≤ 1e-8, get_theta(cuteh))]
    fh_func(θ) = fh[findfirst(x -> abs(x-θ) ≤ 1e-8, get_theta(cuteh))]
    cutrhcp2 = eh2bor1cut(get_theta(cuteh), fe_func, fh_func; pol=:rhcp, xpd=-30)
    @test cutrhcp ≈ cutrhcp2
end

@safetestset "showcut" begin 
    using TicraUtilities
    cutfile = joinpath(@__DIR__, "ticra_hpol_horn.cut")
    cuteh = read_cutfile(cutfile)
    iobuf = IOBuffer()
    show(iobuf, MIME"text/plain"(), cuteh)
    str = String(take!(iobuf))
    str2 = """Cut
             ncomp\t2
             icut \t1
             icomp\t3
             text
                 \tField data in cuts
                 \tField data in cuts
                 \tField data in cuts
             theta\t0.0:0.5:180.0
             phi  \t0.0:45.0:90.0
             evec \t361×3 Matrix{StaticArraysCore.SVector{2, ComplexF64}}
           """
    @test str == str2
end