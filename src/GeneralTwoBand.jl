export GeneralTwoBand

abstract type GeneralTwoBand{T} <: Hamiltonian{T} end

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

export σx_cc,σy_cc,σz_cc
σx_cc(h::GeneralTwoBand,kx,ky) = zero(eltype(promote(kx,ky)))
σy_cc(h::GeneralTwoBand,kx,ky) = zero(eltype(promote(kx,ky)))
σz_cc(h::GeneralTwoBand,kx,ky) = zero(eltype(promote(kx,ky)))

export σx_vv,σy_vv,σz_vv
σx_vv(h::GeneralTwoBand,kx,ky) = zero(eltype(promote(kx,ky)))
σy_vv(h::GeneralTwoBand,kx,ky) = zero(eltype(promote(kx,ky)))
σz_vv(h::GeneralTwoBand,kx,ky) = zero(eltype(promote(kx,ky)))

σx_cv(hx::Number,hy::Number,hz::Number) = (im * hy + hx*hz/ϵ(hx,hy,hz)) / (hx + im*hy)
σy_cv(hx::Number,hy::Number,hz::Number) = (-im * hx + hy*hz/ϵ(hx,hy,hz)) / (hx + im*hy)
σz_cv(hx::Number,hy::Number,hz::Number) = (-hx + im*hy) / ϵ(hx,hy,hz)
σx_vc(hx::Number,hy::Number,hz::Number) = conj(σx_cv(hx,hy,hz))
σy_vc(hx::Number,hy::Number,hz::Number) = conj(σy_cv(hx,hy,hz))
σz_vc(hx::Number,hy::Number,hz::Number) = conj(σz_cv(hx,hy,hz))

σx_cc(hx::Number,hy::Number,hz::Number) = zero(eltype(promote(hx,hy,hz)))
σy_cc(hx::Number,hy::Number,hz::Number) = zero(eltype(promote(hx,hy,hz)))
σz_cc(hx::Number,hy::Number,hz::Number) = zero(eltype(promote(hx,hy,hz)))

σx_vv(hx::Number,hy::Number,hz::Number) = zero(eltype(promote(hx,hy,hz)))
σy_vv(hx::Number,hy::Number,hz::Number) = zero(eltype(promote(hx,hy,hz)))
σz_vv(hx::Number,hy::Number,hz::Number) = zero(eltype(promote(hx,hy,hz)))

σx_cv(h::SVector{3,<:Number}) = (im * h[2] + h[1]*h[3]/ϵ(h)) / (h[1] + im*h[2])
σy_cv(h::SVector{3,<:Number}) = (-im * h[1] + h[2]*h[3]/ϵ(h)) / (h[1] + im*h[2])
σz_cv(h::SVector{3,<:Number}) = (-h[1] + im*h[2]) / ϵ(h)
σx_vc(h::SVector{3,<:Number}) = conj(σx_cv(h))
σy_vc(h::SVector{3,<:Number}) = conj(σy_cv(h))
σz_vc(h::SVector{3,<:Number}) = conj(σz_cv(h))


export vμ_cv,vμ_vc
function vμ_cv(h::SVector{3,<:Number},dh::SVector{3,<:Number})
    return dh[1]*σx_cv(h) + dh[2]*σy_cv(h) + dh[3]*σz_cv(h) 
end
vμ_vc(h::SVector{3,<:Number},dh::SVector{3,<:Number}) = conj(vμ_cv(h,dh))

export vμ_cc,vμ_vv
function vμ_cc(h::SVector{3,<:Number},dh::SVector{3,<:Number})
    return (h[1]*dh[1] + h[2]*dh[2] +h[3]*dh[3]) / ϵ(h)
end
function vμ_vv(h::SVector{3,<:Number},dh::SVector{3,<:Number})
    return -(h[1]*dh[1] + h[2]*dh[2] +h[3]*dh[3]) / ϵ(h)
end
