export GappedDiracOld
"""
    GappedDiracOld{T<:Real} <: Hamiltonian{T}

Legacy hard-coded. Holds the parameter of a dimensionless massive Dirac Hamiltonian.

The Hamiltonian reads ``\\hat{H} =  k_x\\sigma_x + k_y\\sigma_y + \\Delta\\sigma_z`` such 
that ``2\\Delta`` is the bandgap at ``\\vec{k}=0``

# Examples
```jldoctest
julia> h = GappedDiracOld(0.1)
GappedDiracOld{Float64}(0.1)
```

# See also
[`UnitScaling(timescale,lengthscale)`](@ref)
"""
struct GappedDiracOld{T<:Real} <: Hamiltonian{T}
    Δ::T
end


function scalegapped_dirac(mass::Unitful.Energy,fermivelocity::Unitful.Velocity)
    tc = uconvert(u"fs",Unitful.ħ/mass)
    lc = uconvert(u"nm",fermivelocity*tc)
    us = UnitScaling(tc,lc)
    return us,GappedDiracOld(mass)
end

getparams(h::GappedDiracOld) = (Δ=h.Δ,vF=one(h.Δ))

getϵ(h::GappedDiracOld)     = (kx,ky) -> sqrt(kx^2+ky^2+h.Δ^2)
getdx_cc(h::GappedDiracOld) = (kx,ky) -> ky * (one(h.t1) -h.Δ/sqrt(kx^2+ky^2+h.Δ^2)) / (2*(kx^2 + ky^2))
getdx_cv(h::GappedDiracOld) = (kx,ky) -> (ky/sqrt(kx^2+ky^2+h.Δ^2) - im*kx*h.Δ / (kx^2+ky^2+h.Δ^2)) / (2*(kx + im*ky))
getdx_vc(h::GappedDiracOld) = (kx,ky) -> (ky/sqrt(kx^2+ky^2+h.Δ^2) + im*kx*h.Δ / (kx^2+ky^2+h.Δ^2)) / (2*(kx - im*ky))
getdx_vv(h::GappedDiracOld) = (kx,ky) -> -ky * (one(h.t1) -h.Δ/sqrt(kx^2+ky^2+h.Δ^2)) / (2*(kx^2 + ky^2))

getvx_cc(h::GappedDiracOld) = (kx,ky) -> kx/sqrt(kx^2+ky^2+h.Δ^2)
getvx_cv(h::GappedDiracOld) = (kx,ky) -> (h.Δ*kx/sqrt(kx^2+ky^2+h.Δ^2) + im*ky) / (kx + im*ky)
getvx_vc(h::GappedDiracOld) = (kx,ky) -> (h.Δ*kx/sqrt(kx^2+ky^2+h.Δ^2) - im*ky) / (kx - im*ky)
getvx_vv(h::GappedDiracOld) = (kx,ky) -> -kx/sqrt(kx^2+ky^2+h.Δ^2)
getvy_cc(h::GappedDiracOld) = (kx,ky) -> ky/sqrt(kx^2+ky^2+h.Δ^2)
getvy_cv(h::GappedDiracOld) = (kx,ky) -> (h.Δ*ky/sqrt(kx^2+ky^2+h.Δ^2) - im*kx) / (kx + im*ky)
getvy_vc(h::GappedDiracOld) = (kx,ky) -> (h.Δ*ky/sqrt(kx^2+ky^2+h.Δ^2) + im*kx) / (kx - im*ky)
getvy_vv(h::GappedDiracOld) = (kx,ky) -> -ky/sqrt(kx^2+ky^2+h.Δ^2)

getΔϵ(h::GappedDiracOld)   = (kx,ky) -> 2sqrt(kx^2+ky^2+h.Δ^2)
getΔvx(h::GappedDiracOld)  = (kx,ky) -> 2kx/sqrt(kx^2+ky^2+h.Δ^2)
getΔvy(h::GappedDiracOld)  = (kx,ky) -> 2ky/sqrt(kx^2+ky^2+h.Δ^2)
getΔv(h::GappedDiracOld)   = (getΔvx(h),getΔvy(h))

getgxx(h::GappedDiracOld)  = (kx,ky) -> 2*(h.Δ^2+ky^2) / (h.Δ^2+kx^2+ky^2)^(3/2)
getgxy(h::GappedDiracOld)  = (kx,ky) -> 2*(-kx*ky) / (h.Δ^2+kx^2+ky^2)^(3/2)
getgyx(h::GappedDiracOld)  = getgxy(h)
getgyy(h::GappedDiracOld)  = (kx,ky) -> 2*(h.Δ^2+kx^2) / (h.Δ^2+kx^2+ky^2)^(3/2)


getdipoles_x(h::GappedDiracOld)    = (getdx_cc(h),getdx_cv(h),getdx_vc(h),getdx_vv(h))
getvels_x(h::GappedDiracOld)       = (getvx_cc(h),getvx_cv(h),getvx_vc(h),getvx_vv(h))
getvels_y(h::GappedDiracOld)       = (getvy_cc(h),getvy_cv(h),getvy_vc(h),getvy_vv(h))

function getcurvature(h::GappedDiracOld)

    gxx = getgxx(h)
    gxy = getgxy(h)
    gyx = getgyx(h)
    gyy = getgyy(h)

    return (kx,ky) -> [[gxx(kx,ky), gyx(kx,ky)], [gxy(kx,ky), gyy(kx,ky)]] # column-major!
end

function printparamsSI(h::GappedDiracOld,us::UnitScaling;digits=3)

    p   = getparams(h)
    Δ   = round(energySI(p.Δ,us),sigditis=digits)
    vF  = round(velocitySI(p.vF,us),sigdigits=digits)
    
    return "m = $Δ ($(h.Δ))\nvF = $vF ($(1.0))"
end

