using TwoBandHamiltonians
using Test
using Unitful

import TwoBandHamiltonians.getvels_x
import TwoBandHamiltonians.getvels_y

const vels = [vx_cc,vx_cv,vx_vc,vx_vv,vy_cc,vy_cv,vy_vc,vy_vv]


@testset "TwoBandHamiltonians" begin
    krange  = -1:0.01:1
    ks      = krange[krange .!= 0.0]
    @testset "GappedDirac" begin 
        m = 0.1
        h = GappedDirac(m)
        href = GappedDiracOld(m)
        vels_ref = [getvels_x(href)...,getvels_y(href)...]

        @testset "$v" for (v,vref) in zip(vels,vels_ref)
            @test [v(h,kx,ky) ≈ vref(kx,ky) for kx=ks,ky=ks] |> all
        end
    end
    
end
