
export GappedDirac
"""
    GappedDirac{T<:Real} <: GeneralTwoBand{T}

Holds the parameter of a dimensionless massive Dirac Hamiltonian.

The Hamiltonian reads 
```julia
\\hat{H} = k_x\\sigma_x + k_y\\sigma_y + m\\sigma_z
``` 
such that ``\vec{h}=[k_x,k_y,m]``.

# See also
[`GeneralTwoBand`](@ref TwoBandHamiltonians.GeneralTwoBand)
"""
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