using SafeTestsets

@safetestset "TorFile" begin 
    using TicraUtilities
    t = read_torfile(joinpath(@__DIR__, "tabulated_rim_tor_file.tor"))
    @test length(t) == 2
    @test t[1].name == "east_9m_rim"
    @test t[1].objtype == "tabulated_rim_xy"
    @test t[1].propname == 
        ["file_name", "unit", "number_of_points", "translation", "polar_origin"]
    @test t[1].propval ==
        ["h_9m_scalloped_rim.xyz", "in", "112", "struct(x: 0.0 in, y: 0.0 in)",
        "struct(status: automatic, x: 0.0 in, y: 0.0 in)"]
    @test t[2].name == "extsub_sw.rim"
    @test t[2].objtype == "tabulated_rim_xy"
    @test t[2].propname == ["file_name", "unit", "number_of_points", "translation", "polar_origin"]
    @test t[2].propval ==
    ["extsub_sw.rim", "in", "360", "struct(x: 0.0 in, y: 0.0 in)",
    "struct(status: automatic, x: 0.0 in, y: 0.0 in)"]
    @test TicraUtilities._parse_tor_xy_struct(t[2].propval[4]) == (0.0, 0.0, "in", "in")
end

