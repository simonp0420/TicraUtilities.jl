using SafeTestsets

@safetestset "sor_efficiency" begin 
    using TicraUtilities
    cd(@__DIR__)
    cutfile = joinpath(@__DIR__, "center_element_rhcp_excited.cut")
    t1 = sor_efficiency(cutfile, F=10, D=18, Oc=0.4, pol=:l3h, dz=-0.1)
    @test t1.ηₛ ≈ 0.8727422758602933
    @test t1.ηₓ ≈ 0.4985615310506499
    @test t1.ηᵢₗ ≈ 0.8519300087313848
    @test t1.ηₚ ≈ 0.9635018253504125

    # Higher resolution for the ϕ integrals
    cut = read_cutfile(cutfile)
    nphi = 2 * length(get_phi(cut))
    Δphi = 360 / nphi
    phi = range(0,360-Δphi, nphi)
    swe_julia = cut2sph(cutfile)
    newcut = sph2cut(swe_julia; phi)
    t2 = sor_efficiency(newcut, F=10, D=18, Oc=0.4, pol=:l3h, dz=-0.1)
    @test maximum(abs(t1[i]-t2[i]) for i in 1:4) < 1e-6
    @test abs(prod(t1) - prod(t2)) < 1e-7
end
