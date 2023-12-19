export GeneralTwoBand

abstract type GeneralTwoBand{T} <: Hamiltonian{T} end

gethvec(h::GeneralTwoBand) = [gethx(h),gethy(h),gethz(h)]

export Δϵ
export ϵ
Δϵ(h::GeneralTwoBand,kx,ky) = 2ϵ(h,kx,ky)
ϵ(h::GeneralTwoBand,kx,ky)  = sqrt(hx(h,kx,ky)^2 + hy(h,kx,ky)^2 + hz(h,kx,ky)^2)


export σx_cv,σy_cv,σz_cv
σx_cv(h::GeneralTwoBand,kx,ky) = σx_cv(hx(h,kx,ky),hy(h,kx,ky),hz(h,kx,ky))
σy_cv(h::GeneralTwoBand,kx,ky) = σy_cv(hx(h,kx,ky),hy(h,kx,ky),hz(h,kx,ky))
σz_cv(h::GeneralTwoBand,kx,ky) = σz_cv(hx(h,kx,ky),hy(h,kx,ky),hz(h,kx,ky))
σx_vc(h::GeneralTwoBand,kx,ky) = conj(σx_cv(h,kx,ky))
σy_vc(h::GeneralTwoBand,kx,ky) = conj(σy_cv(h,kx,ky))
σz_vc(h::GeneralTwoBand,kx,ky) = conj(σz_cv(h,kx,ky))

σx_cc(h::GeneralTwoBand,kx,ky) = zero(eltype(promote(kx,ky)))
σy_cc(h::GeneralTwoBand,kx,ky) = zero(eltype(promote(kx,ky)))
σz_cc(h::GeneralTwoBand,kx,ky) = zero(eltype(promote(kx,ky)))

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



function getΔϵ(h::GeneralTwoBand)
    hx,hy,hz = gethvec(h)
    return (kx,ky) -> 2sqrt(hx(kx,ky)^2 + hy(kx,ky)^2 + hz(kx,ky)^2)
end

export DiracTry
struct DiracTry{T<:Real} <: GeneralTwoBand{T}
    Δ::T
end
gethx(h::DiracTry) = (kx,ky) -> kx
gethy(h::DiracTry) = (kx,ky) -> ky
gethz(h::DiracTry) = (kx,ky) -> h.Δ

export hx,hy,hz
hx(h::DiracTry,kx,ky) = kx
hy(h::DiracTry,kx,ky) = ky
hz(h::DiracTry,kx,ky) = h.Δ