export GappedDirac
struct GappedDirac{T<:Real} <: Hamiltonian{T}
    Δ::T
end


function scalegapped_dirac(mass::Unitful.Energy,fermivelocity::Unitful.Velocity)
    tc = uconvert(u"fs",Unitful.ħ/mass)
    lc = uconvert(u"nm",fermivelocity*tc)
    us = UnitScaling(tc,lc)
    return us,GappedDirac(mass)
end

getparams(h::GappedDirac) = (Δ=h.Δ,vF=one(h.Δ))

getϵ(h::GappedDirac)     = (kx,ky) -> sqrt(kx^2+ky^2+h.Δ^2)
getdx_cc(h::GappedDirac) = (kx,ky) -> ky * (one(h.t1) -h.Δ/sqrt(kx^2+ky^2+h.Δ^2)) / (2*(kx^2 + ky^2))
getdx_cv(h::GappedDirac) = (kx,ky) -> (ky/sqrt(kx^2+ky^2+h.Δ^2) - im*kx*h.Δ / (kx^2+ky^2+h.Δ^2)) / (2*(kx + im*ky))
getdx_vc(h::GappedDirac) = (kx,ky) -> (ky/sqrt(kx^2+ky^2+h.Δ^2) + im*kx*h.Δ / (kx^2+ky^2+h.Δ^2)) / (2*(kx - im*ky))
getdx_vv(h::GappedDirac) = (kx,ky) -> -ky * (one(h.t1) -h.Δ/sqrt(kx^2+ky^2+h.Δ^2)) / (2*(kx^2 + ky^2))

getvx_cc(h::GappedDirac) = (kx,ky) -> kx/sqrt(kx^2+ky^2+h.Δ^2)
getvx_cv(h::GappedDirac) = (kx,ky) -> (h.Δ*kx/sqrt(kx^2+ky^2+h.Δ^2) + im*ky) / (kx + im*ky)
getvx_vc(h::GappedDirac) = (kx,ky) -> (h.Δ*kx/sqrt(kx^2+ky^2+h.Δ^2) - im*ky) / (kx - im*ky)
getvx_vv(h::GappedDirac) = (kx,ky) -> -kx/sqrt(kx^2+ky^2+h.Δ^2)
getvy_cc(h::GappedDirac) = (kx,ky) -> ky/sqrt(kx^2+ky^2+h.Δ^2)
getvy_cv(h::GappedDirac) = (kx,ky) -> (h.Δ*ky/sqrt(kx^2+ky^2+h.Δ^2) - im*kx) / (kx + im*ky)
getvy_vc(h::GappedDirac) = (kx,ky) -> (h.Δ*ky/sqrt(kx^2+ky^2+h.Δ^2) + im*kx) / (kx - im*ky)
getvy_vv(h::GappedDirac) = (kx,ky) -> -ky/sqrt(kx^2+ky^2+h.Δ^2)

getΔϵ(h::GappedDirac)   = (kx,ky) -> 2sqrt(kx^2+ky^2+h.Δ^2)
getΔvx(h::GappedDirac)  = (kx,ky) -> 2kx/sqrt(kx^2+ky^2+h.Δ^2)
getΔvy(h::GappedDirac)  = (kx,ky) -> 2ky/sqrt(kx^2+ky^2+h.Δ^2)
getΔv(h::GappedDirac)   = (getΔvx(h),getΔvy(h))

getgxx(h::GappedDirac)  = (kx,ky) -> 2*(h.Δ^2+ky^2) / (h.Δ^2+kx^2+ky^2)^(3/2)
getgxy(h::GappedDirac)  = (kx,ky) -> 2*(-kx*ky) / (h.Δ^2+kx^2+ky^2)^(3/2)
getgyx(h::GappedDirac)  = getgxy(h)
getgyy(h::GappedDirac)  = (kx,ky) -> 2*(h.Δ^2+kx^2) / (h.Δ^2+kx^2+ky^2)^(3/2)


getdipoles_x(h::GappedDirac)    = (getdx_cc(h),getdx_cv(h),getdx_vc(h),getdx_vv(h))
getvels_x(h::GappedDirac)       = (getvx_cc(h),getvx_cv(h),getvx_vc(h),getvx_vv(h))
getvels_y(h::GappedDirac)       = (getvy_cc(h),getvy_cv(h),getvy_vc(h),getvy_vv(h))

function getcurvature(h::GappedDirac)

    gxx = getgxx(h)
    gxy = getgxy(h)
    gyx = getgyx(h)
    gyy = getgyy(h)

    return (kx,ky) -> [[gxx(kx,ky), gyx(kx,ky)], [gxy(kx,ky), gyy(kx,ky)]] # column-major!
end

function printparamsSI(h::GappedDirac,us::UnitScaling;digits=3)

    p   = getparams(h)
    Δ   = round(energySI(p.Δ,us),sigditis=digits)
    vF  = round(velocitySI(p.vF,us),sigdigits=digits)
    
    return "m = $Δ ($(h.Δ))\nvF = $vF ($(1.0))"
end

