
export GappedDirac
struct GappedDirac{T<:Real} <: GeneralTwoBand{T} 
    m::T
end

hx(h::GappedDirac,kx,ky) = kx
hy(h::GappedDirac,kx,ky) = ky
hz(h::GappedDirac,kx,ky) = h.m