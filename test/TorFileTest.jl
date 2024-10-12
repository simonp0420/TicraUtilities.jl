#using SafeTestsets
using Test
using Unitful: @u_str
@testset "TorFile" begin 
    using TicraUtilities
    t = read_torfile(joinpath(@__DIR__, "tabulated_rim_tor_file.tor"))
    @test length(t) == 2
    @test get_name(t[1]) == "east_9m_rim"
    @test get_objtype(t[1]) == "tabulated_rim_xy"
    @test get_propname(t[1]) == 
        ["file_name", "unit", "number_of_points", "translation", "polar_origin"]
    @test get_propval(t[1]) ==
        ["h_9m_scalloped_rim.xyz", "in", "112", "struct(x: 0.0 in, y: 0.0 in)",
        "struct(status: automatic, x: 0.0 in, y: 0.0 in)"]
    @test get_name(t[2]) == "extsub_sw.rim"
    @test get_objtype(t[2]) == "tabulated_rim_xy"
    @test get_propname(t[2]) == ["file_name", "unit", "number_of_points", "translation", "polar_origin"]
    @test get_propval(t[2]) ==
    ["extsub_sw.rim", "in", "360", "struct(x: 0.0 in, y: 0.0 in)",
    "struct(status: automatic, x: 0.0 in, y: 0.0 in)"]

    tfile = joinpath(tempdir(), "temp.tor")
    write_torfile(tfile, t)
    t = read_torfile(tfile)
    @test length(t) == 2
    @test get_name(t[1]) == "east_9m_rim"
    @test get_objtype(t[1]) == "tabulated_rim_xy"
    @test get_propname(t[1]) == 
        ["file_name", "unit", "number_of_points", "translation", "polar_origin"]
    @test get_propval(t[1]) ==
        ["h_9m_scalloped_rim.xyz", "in", "112", "struct(x: 0.0 in, y: 0.0 in)",
        "struct(status: automatic, x: 0.0 in, y: 0.0 in)"]
    @test get_name(t[2]) == "extsub_sw.rim"
    @test get_objtype(t[2]) == "tabulated_rim_xy"
    @test get_propname(t[2]) == ["file_name", "unit", "number_of_points", "translation", "polar_origin"]
    @test get_propval(t[2]) ==
    ["extsub_sw.rim", "in", "360", "struct(x: 0.0 in, y: 0.0 in)",
    "struct(status: automatic, x: 0.0 in, y: 0.0 in)"]
    

    @test TicraUtilities.parse_tor_struct(t[1].propval[end]) == (; status = "automatic", x = 0.0u"inch", y = 0.0u"inch")
end

