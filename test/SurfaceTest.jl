using SafeTestsets

@safetestset "Surface" begin 
    using TicraUtilities
    sfc = read_surface(joinpath(@__DIR__, "parent_parabola.sfc"))
    @test get_x(sfc) == range(0.1981760000E+02, 0.1232176000E+03, length=301)
    @test get_y(sfc) == range(-0.5170000000E+02, 0.5170000000E+02, length=301)
    zs = get_z(sfc)
    @test zs[1,1] == 0.8230755705E+01
    @test zs[end] == 0.4793928731E+02
    @test size(zs) == (301,301)
    @test get_text(sfc)  == "Surface z-values for reflector: main_reflector, with unit: in."
    fname = tempname()
    write_surface(fname, sfc)
    sfc2 = read_surface(fname)
    @test get_x(sfc) == get_x(sfc2)
    @test get_y(sfc) == get_y(sfc2)
    @test get_z(sfc) == get_z(sfc2)
    @test get_text(sfc) == get_text(sfc2)
    sfc_plus = sfc + sfc2
    @test get_x(sfc_plus) == get_x(sfc2)
    @test get_y(sfc_plus) == get_y(sfc2)
    @test get_z(sfc_plus) ≈ 2 * get_z(sfc2)
    sfc_minus = sfc - sfc2
    @test get_x(sfc_minus) == get_x(sfc2)
    @test get_y(sfc_minus) == get_y(sfc2)
    @test get_z(sfc_minus) ≈ 0 * get_z(sfc2)
end

