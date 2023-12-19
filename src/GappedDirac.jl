
export GappedDirac
struct GappedDirac{T<:Real} <: GeneralTwoBand{T} 
    m::T
end

hx(h::GappedDirac,kx,ky) = kx
hy(h::GappedDirac,kx,ky) = ky
hz(h::GappedDirac,kx,ky) = h.m

dhdkx(h::GappedDirac,kx,ky) = SA[one(typeof(h.m)),zero(typeof(h.m)),zero(typeof(h.m))]
dhdky(h::GappedDirac,kx,ky) = SA[zero(typeof(h.m)),one(typeof(h.m)),zero(typeof(h.m))]

# Jacobian ∂h_i/∂k_j
jac(h::GappedDirac,kx,ky) = SA[
    one(typeof(h.m)) zero(typeof(h.m)) zero(typeof(h.m))
    zero(typeof(h.m)) one(typeof(h.m)) zero(typeof(h.m))
    zero(typeof(h.m)) zero(typeof(h.m)) zero(typeof(h.m))]