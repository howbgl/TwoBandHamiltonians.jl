export GeneralTwoBand
"""
    GeneralTwoBand{T} <: Hamiltonian{T}

Supertype of all 2x2 Hamiltonians with all matrixelements via dispatch.

# Idea
The idea is that all Hamiltonians of the form

```math
\\hat{H} = \\vec{h}(\\vec{k})\\cdot\\vec{\\sigma}
```

can be diagonalized analytically and hence most desired matrixelements such as velocities or
dipoles can be expressed solely through 
    
```math
\\vec{h}(\\vec{k})=[h_x(\\vec{k}),h_y(\\vec{k}),h_z(\\vec{k})]
```

and its derivatives with respect to ``k_\\mu``. Any particular Hamiltonian deriving form
GeneralTwoBand{T} must then only implement ``\\vec{h}(\\vec{k})`` and its derivatives.

# See also
[`GappedDirac`](@ref TwoBandHamiltonians.GappedDirac)
"""
abstract type GeneralTwoBand{T} <: Hamiltonian{T} end

export hvec
"""
    hvec(h,kx,ky)

Returns ``\\vec{h}(\\vec{k})`` defining the Hamiltonian at ``\vec{k}=[k_x,k_y]``.
"""
hvec(h::GeneralTwoBand,kx,ky) = SA[hx(h,kx,ky),hy(h,kx,ky),hz(h,kx,ky)]

export Δϵ,ϵ
Δϵ(h::GeneralTwoBand,kx,ky) = 2ϵ(h,kx,ky)
ϵ(h::GeneralTwoBand,kx,ky)  = sqrt(hx(h,kx,ky)^2 + hy(h,kx,ky)^2 + hz(h,kx,ky)^2)


Δϵ(hx::Number,hy::Number,hz::Number)    = 2sqrt(hx^2 + hy^2 + hz^2)
ϵ(hx::Number,hy::Number,hz::Number)     = sqrt(hx^2 + hy^2 + hz^2)

Δϵ(h::SVector{3,<:Number})  = 2sqrt(h[1]^2 + h[2]^2 + h[3]^2)
ϵ(h::SVector{3,<:Number})   = sqrt(h[1]^2 + h[2]^2 + h[3]^2)


export σx_cv,σy_cv,σz_cv,σx_vc,σy_vc,σz_vc
σx_cv(h::GeneralTwoBand,kx,ky) = σx_cv(hx(h,kx,ky),hy(h,kx,ky),hz(h,kx,ky))
σy_cv(h::GeneralTwoBand,kx,ky) = σy_cv(hx(h,kx,ky),hy(h,kx,ky),hz(h,kx,ky))
σz_cv(h::GeneralTwoBand,kx,ky) = σz_cv(hx(h,kx,ky),hy(h,kx,ky),hz(h,kx,ky))
σx_vc(h::GeneralTwoBand,kx,ky) = conj(σx_cv(h,kx,ky))
σy_vc(h::GeneralTwoBand,kx,ky) = conj(σy_cv(h,kx,ky))
σz_vc(h::GeneralTwoBand,kx,ky) = conj(σz_cv(h,kx,ky))

σx_cv(hx::Number,hy::Number,hz::Number) = (im * hy + hx*hz/ϵ(hx,hy,hz)) / (hx + im*hy)
σy_cv(hx::Number,hy::Number,hz::Number) = (-im * hx + hy*hz/ϵ(hx,hy,hz)) / (hx + im*hy)
σz_cv(hx::Number,hy::Number,hz::Number) = (-hx + im*hy) / ϵ(hx,hy,hz)
σx_vc(hx::Number,hy::Number,hz::Number) = (-im * hy + hx*hz/ϵ(hx,hy,hz)) / (hx - im*hy)
σy_vc(hx::Number,hy::Number,hz::Number) = (im * hx + hy*hz/ϵ(hx,hy,hz)) / (hx - im*hy)
σz_vc(hx::Number,hy::Number,hz::Number) = (-hx - im*hy) / ϵ(hx,hy,hz)

σx_cv(h::SVector{3,<:Number}) = (im * h[2] + h[1]*h[3]/ϵ(h)) / (h[1] + im*h[2])
σy_cv(h::SVector{3,<:Number}) = (-im * h[1] + h[2]*h[3]/ϵ(h)) / (h[1] + im*h[2])
σz_cv(h::SVector{3,<:Number}) = (-h[1] + im*h[2]) / ϵ(h)
σx_vc(h::SVector{3,<:Number}) = (-im * h[2] + h[1]*h[3]/ϵ(h)) / (h[1] - im*h[2])
σy_vc(h::SVector{3,<:Number}) = (im * h[1] + h[2]*h[3]/ϵ(h)) / (h[1] - im*h[2])
σz_vc(h::SVector{3,<:Number}) = (-h[1] - im*h[2]) / ϵ(h)


export vμ_cv,vμ_vc
function vμ_cv(h::SVector{3,<:Number},dh::SVector{3,<:Number})
    return dh[1]*σx_cv(h) + dh[2]*σy_cv(h) + dh[3]*σz_cv(h) 
end
function vμ_vc(h::SVector{3,<:Number},dh::SVector{3,<:Number})
    return dh[1]*σx_vc(h) + dh[2]*σy_vc(h) + dh[3]*σz_vc(h) 
end

export vμ_cc,vμ_vv
function vμ_cc(h::SVector{3,<:Number},dh::SVector{3,<:Number})
    return (h[1]*dh[1] + h[2]*dh[2] + h[3]*dh[3]) / ϵ(h)
end
function vμ_vv(h::SVector{3,<:Number},dh::SVector{3,<:Number})
    return -(h[1]*dh[1] + h[2]*dh[2] + h[3]*dh[3]) / ϵ(h)
end

export vx_cv,vx_vc,vx_cc,vx_vv
vx_cv(h::GeneralTwoBand,kx,ky) = vμ_cv(hvec(h,kx,ky),dhdkx(h,kx,ky))
vx_vc(h::GeneralTwoBand,kx,ky) = vμ_vc(hvec(h,kx,ky),dhdkx(h,kx,ky))
vx_cc(h::GeneralTwoBand,kx,ky) = vμ_cc(hvec(h,kx,ky),dhdkx(h,kx,ky))
vx_vv(h::GeneralTwoBand,kx,ky) = vμ_vv(hvec(h,kx,ky),dhdkx(h,kx,ky))

export vy_cv,vy_vc,vy_cc,vy_vv
vy_cv(h::GeneralTwoBand,kx,ky) = vμ_cv(hvec(h,kx,ky),dhdky(h,kx,ky))
vy_vc(h::GeneralTwoBand,kx,ky) = vμ_vc(hvec(h,kx,ky),dhdky(h,kx,ky))
vy_cc(h::GeneralTwoBand,kx,ky) = vμ_cc(hvec(h,kx,ky),dhdky(h,kx,ky))
vy_vv(h::GeneralTwoBand,kx,ky) = vμ_vv(hvec(h,kx,ky),dhdky(h,kx,ky))



n2(h::SVector{3,<:Number}) = 2ϵ(h)*(ϵ(h) + h[3])

function dμ_cv(h::SVector{3,<:Number},dh::SVector{3,<:Number})
    -im*vμ_cv(h,dh) / Δϵ(h)
end
function dμ_vc(h::SVector{3,<:Number},dh::SVector{3,<:Number})
    im*vμ_vc(h,dh) / Δϵ(h)
end

export dx_cc,dy_cc,dx_vv,dy_vv,dx_cv,dx_vc,dy_cv,dy_vc
function dx_cc(h::SVector{3,<:Number},jac::SMatrix{3,3,<:Number})
    return (jac[1,1]*h[2] - h[1]*jac[2,1]) / n2(h)
end
function dy_cc(h::SVector{3,<:Number},jac::SMatrix{3,3,<:Number})
    return (jac[1,2]*h[2] - h[1]*jac[2,2]) / n2(h)
end
function dx_vv(h::SVector{3,<:Number},jac::SMatrix{3,3,<:Number})
    return -(jac[1,1]*h[2] - h[1]*jac[2,1]) / n2(h)
end
function dy_vv(h::SVector{3,<:Number},jac::SMatrix{3,3,<:Number})
    return -(jac[1,2]*h[2] - h[1]*jac[2,2]) / n2(h)
end

dx_cv(h::GeneralTwoBand,kx,ky) = dμ_cv(hvec(h,kx,ky),dhdkx(h,kx,ky))
dx_vc(h::GeneralTwoBand,kx,ky) = dμ_vc(hvec(h,kx,ky),dhdkx(h,kx,ky))
dx_cc(h::GeneralTwoBand,kx,ky) = dx_cc(hvec(h,kx,ky),jac(h,kx,ky))
dx_vv(h::GeneralTwoBand,kx,ky) = dx_vv(hvec(h,kx,ky),jac(h,kx,ky))

dy_cv(h::GeneralTwoBand,kx,ky) = dμ_cv(hvec(h,kx,ky),dhdky(h,kx,ky))
dy_vc(h::GeneralTwoBand,kx,ky) = dμ_vc(hvec(h,kx,ky),dhdky(h,kx,ky))
dy_cc(h::GeneralTwoBand,kx,ky) = dy_cc(hvec(h,kx,ky),jac(h,kx,ky))
dy_vv(h::GeneralTwoBand,kx,ky) = dy_vv(hvec(h,kx,ky),jac(h,kx,ky))


