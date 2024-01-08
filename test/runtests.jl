using TwoBandHamiltonians
using Test
using Unitful

# Hand-coded reference formulas for GappedDirac
gappeddirac_vx_cc(m,kx,ky) = kx/sqrt(kx^2+ky^2+m^2)
gappeddirac_vx_vv(m,kx,ky) = -kx/sqrt(kx^2+ky^2+m^2)
gappeddirac_vx_cv(m,kx,ky) = (m*kx/sqrt(kx^2+ky^2+m^2) + im*ky) / (kx + im*ky)
gappeddirac_vx_vc(m,kx,ky) = conj(gappeddirac_vx_cv(m,kx,ky))

gappeddirac_vy_cc(m,kx,ky) = ky/sqrt(kx^2+ky^2+m^2)
gappeddirac_vy_vv(m,kx,ky) = -ky/sqrt(kx^2+ky^2+m^2)
gappeddirac_vy_cv(m,kx,ky) = (m*ky/sqrt(kx^2+ky^2+m^2) - im*kx) / (kx + im*ky)
gappeddirac_vy_vc(m,kx,ky) = conj(gappeddirac_vy_cv(m,kx,ky))

gappeddirac_dx_cc(m,kx,ky) = ky * (1 -m/sqrt(kx^2+ky^2+m^2)) / (2*(kx^2 + ky^2))
gappeddirac_dx_vv(m,kx,ky) = -gappeddirac_dx_cc(m,kx,ky)
gappeddirac_dx_cv(m,kx,ky) = (ky/sqrt(kx^2+ky^2+m^2) - im*kx*m / (kx^2+ky^2+m^2)) / (2*(kx + im*ky))
gappeddirac_dx_vc(m,kx,ky) = conj(gappeddirac_dx_cv(m,kx,ky))

gappeddirac_dy_cc(m,kx,ky) = kx * (-1 + m/sqrt(kx^2+ky^2+m^2)) / (2*(kx^2 + ky^2))
gappeddirac_dy_vv(m,kx,ky) = - gappeddirac_dy_cc(m,kx,ky)
gappeddirac_dy_cv(m,kx,ky) = (-kx/sqrt(kx^2+ky^2+m^2) - im*m*ky/(kx^2+ky^2+m^2)) / (2*(kx + im*ky))
gappeddirac_dy_vc(m,kx,ky) = conj(gappeddirac_dy_cv(m,kx,ky))


const vels = [
    vx_cc,
    vx_cv,
    vx_vc,
    vx_vv,
    vy_cc,
    vy_cv,
    vy_vc,
    vy_vv]

const vels_ref = [
    gappeddirac_vx_cc,
    gappeddirac_vx_cv,
    gappeddirac_vx_vc,
    gappeddirac_vx_vv,
    gappeddirac_vy_cc,
    gappeddirac_vy_cv,
    gappeddirac_vy_vc,
    gappeddirac_vy_vv]

const dipoles = [
    dx_cc,
    dx_cv,
    dx_vc,
    dx_vv,
    dy_cc,
    dy_cv,
    dy_vc,
    dy_vv]

const dipoles_ref = [
    gappeddirac_dx_cc,
    gappeddirac_dx_cv,
    gappeddirac_dx_vc,
    gappeddirac_dx_vv,
    gappeddirac_dy_cc,
    gappeddirac_dy_cv,
    gappeddirac_dy_vc,
    gappeddirac_dy_vv]

@testset "TwoBandHamiltonians" begin
    krange  = -1:0.01:1
    ks      = krange[krange .!= 0.0]
    @testset "GappedDirac" begin 
        m = 0.1
        h = GappedDirac(m)

        @testset "$v" for (v,vref) in zip(vels,vels_ref)
            @test [v(h,kx,ky) ≈ vref(m,kx,ky) for kx=ks,ky=ks] |> all
        end

        @testset "$d" for (d,dref) in zip(dipoles,dipoles_ref)
            @test [d(h,kx,ky) ≈ dref(m,kx,ky) for kx=ks,ky=ks] |> all
        end
    end
    
end
