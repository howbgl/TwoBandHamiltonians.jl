
export HexWarpDirac
struct HexWarpDirac{T<:Real} <: GeneralTwoBand{T} 
    A::T
    R::T
end

hx(h::HexWarpDirac,kx,ky) = h.A*ky
hy(h::HexWarpDirac,kx,ky) = -h.A*kx
hz(h::HexWarpDirac,kx,ky) = 2h.R*(kx^3 - 3kx*ky^2)

dhdkx(h::HexWarpDirac,kx,ky) = SA[zero(h.A),-h.A,6h.R*(kx^2 - ky^2)]
dhdky(h::HexWarpDirac,kx,ky) = SA[h.A,zero(h.A),-12h.R*kx*ky]

# Jacobian ∂h_i/∂k_j
jac(h::HexWarpDirac,kx,ky) = SA[
    zero(h.A)               h.A
    -h.A                    zero(h.A)
    6h.R*(kx^2 - ky^2)      -12h.R*kx*ky]

export gethvec
gethvec(h::HexWarpDirac) = let A=h.A,R=h.R
    (kx,ky) -> SA[A*ky,-A*kx,2R*(kx^3 - 3kx*ky^2)]
end

export getdhdx,getdhdy
getdhdx(h::HexWarpDirac) = let A=h.A,R=h.R
    (kx,ky) -> SA[zero(A),-A,6R*(kx^2 - ky^2)]
end
getdhdy(h::HexWarpDirac) = let A=h.A,R=h.R
    (kx,ky) -> SA[A,zero(A),-12R*kx*ky]
end

export getjac
getjac(h::HexWarpDirac) = let A=h.A,R=h.R
    (kx,ky) -> SA[
        zero(A)             A
        -A                  zero(A)
        6R*(kx^2 - ky^2)    -12R*kx*ky]
end



function printparamsSI(h::HexWarpDirac,us::UnitScaling;digits=3)
    
    return """
        A  = $(round(h.A,sigdigits=digits))
        R = $(round(h.R,sigdigits=digits))"""
end