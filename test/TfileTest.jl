using SafeTestsets

@safetestset "read_tfile" begin 
    using TicraUtilities
    t = read_tfile(joinpath(@__DIR__, "Tfile_scenario1_f13.fmt"))
    @test t.nel == 60
    @test t.npart == 1
    @test t.nsta == [191]
    @test length(t.p1) == length(t.p2) == t.npart
    @test t.p1[1][1,1] == 18.01576189 - 8.525185888im
    @test t.p1[1][end,1] == 0.127669761 + 0.1206916781im
    @test t.p1[1][1,end] == 0.008447299661 + 0.01567232165im
    @test t.p1[1][end,end] == 0.001329830687 - 0.0003094218168im
    @test t.p2[1][1,1] == 142.511723 + 12.07070806im
    @test t.p2[1][end,1] == 0.02885447223 + 0.01976923336im
    @test t.p2[1][1,end] == 0.06710393568 - 0.09723884395im
    @test t.p2[1][end,end] == 0.001444810561 + 5.980552217e-5im
end

