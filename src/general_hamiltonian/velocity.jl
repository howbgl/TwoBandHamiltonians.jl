
export vμ_cv,vμ_vc
function vμ_cv(h::SVector{3,<:Number},dh::SVector{3,<:Number})
    return dh[1]*σx_cv(h) + dh[2]*σy_cv(h) + dh[3]*σz_cv(h) 
end
function vμ_vc(h::SVector{3,<:Number},dh::SVector{3,<:Number})
    return dh[1]*σx_vc(h) + dh[2]*σy_vc(h) + dh[3]*σz_vc(h) 
end

export getvμ_cv,getvμ_vc
getvμ_cv(hvec::Function,dhdμ::Function) = (kx,ky) -> sum(σvec_cv(hvec(kx,ky)) .* dhdμ(kx,ky))
getvμ_vc(hvec::Function,dhdμ::Function) = (kx,ky) -> sum(σvec_vc(hvec(kx,ky)) .* dhdμ(kx,ky))


export vμ_cc,vμ_vv
function vμ_cc(h::SVector{3,<:Number},dh::SVector{3,<:Number})
    return (h[1]*dh[1] + h[2]*dh[2] + h[3]*dh[3]) / ϵ(h)
end
function vμ_vv(h::SVector{3,<:Number},dh::SVector{3,<:Number})
    return -(h[1]*dh[1] + h[2]*dh[2] + h[3]*dh[3]) / ϵ(h)
end

export getvμ_cc,getvμ_vv
function getvμ_cc(hvec::Function,dhdμ::Function)
    return let Δϵ_loc = getΔϵ(hvec)
        (kx,ky) -> sum(hvec(kx,ky) .* dhdμ(kx,ky)) / Δϵ_loc(kx,ky)
    end
end
function getvμ_vv(hvec::Function,dhdμ::Function)
    return let Δϵ_loc = getΔϵ(hvec)
        (kx,ky) -> sum(hvec(kx,ky) .* dhdμ(kx,ky)) / Δϵ_loc(kx,ky)
    end
end


export vx_cv,vx_vc,vx_cc,vx_vv
"""
    vx_cv(h,kx,ky)

Returns the velocity operator matrix element ⟨ψ+|vx|ψ-⟩ at ``\\vec{k}=[k_x,k_y]``.

# Example
```jldoctest
julia> h = GappedDirac(1.0); vx_cv(h,1.0,-1.0)
0.7886751345948129 - 0.21132486540518708im
```
"""
vx_cv(h::GeneralTwoBand,kx,ky) = vμ_cv(hvec(h,kx,ky),dhdkx(h,kx,ky))
"""
    vx_vc(h,kx,ky)

Returns the velocity operator matrix element ⟨ψ-|vx|ψ+⟩ at ``\\vec{k}=[k_x,k_y]``.

# Example
```jldoctest
julia> h = GappedDirac(1.0); vx_vc(h,1.0,-1.0)
0.7886751345948129 + 0.21132486540518708im
```
"""
vx_vc(h::GeneralTwoBand,kx,ky) = vμ_vc(hvec(h,kx,ky),dhdkx(h,kx,ky))
"""
    vx_cc(h,kx,ky)

Returns the velocity operator matrix element ⟨ψ+|vx|ψ+⟩ at ``\\vec{k}=[k_x,k_y]``.

# Example
```jldoctest
julia> h = GappedDirac(1.0); vx_cc(h,1.0,-1.0)
0.5773502691896258
```
"""
vx_cc(h::GeneralTwoBand,kx,ky) = vμ_cc(hvec(h,kx,ky),dhdkx(h,kx,ky))
"""
    vx_vv(h,kx,ky)

Returns the velocity operator matrix element ⟨ψ-|vx|ψ-⟩ at ``\\vec{k}=[k_x,k_y]``.

# Example
```jldoctest
julia> h = GappedDirac(1.0); vx_vv(h,1.0,-1.0)
-0.5773502691896258
```
"""
vx_vv(h::GeneralTwoBand,kx,ky) = vμ_vv(hvec(h,kx,ky),dhdkx(h,kx,ky))

export vy_cv,vy_vc,vy_cc,vy_vv
"""
    vy_cv(h,kx,ky)

Returns the velocity operator matrix element ⟨ψ+|vy|ψ-⟩ at ``\\vec{k}=[k_x,k_y]``.

# Example
```jldoctest
julia> h = GappedDirac(1.0); vy_cv(h,1.0,-1.0)
0.21132486540518708 - 0.7886751345948129im
```
"""
vy_cv(h::GeneralTwoBand,kx,ky) = vμ_cv(hvec(h,kx,ky),dhdky(h,kx,ky))
"""
    vy_vc(h,kx,ky)

Returns the velocity operator matrix element ⟨ψ-|vy|ψ+⟩ at ``\\vec{k}=[k_x,k_y]``.

# Example
```jldoctest
julia> h = GappedDirac(1.0); vy_vc(h,1.0,-1.0)
0.21132486540518708 + 0.7886751345948129im
```
"""
vy_vc(h::GeneralTwoBand,kx,ky) = vμ_vc(hvec(h,kx,ky),dhdky(h,kx,ky))
"""
    vy_cc(h,kx,ky)

Returns the velocity operator matrix element ⟨ψ+|vy|ψ+⟩ at ``\\vec{k}=[k_x,k_y]``.

# Example
```jldoctest
julia> h = GappedDirac(1.0); vy_cc(h,1.0,-1.0)
-0.5773502691896258
```
"""
vy_cc(h::GeneralTwoBand,kx,ky) = vμ_cc(hvec(h,kx,ky),dhdky(h,kx,ky))
"""
    vy_vv(h,kx,ky)

Returns the velocity operator matrix element ⟨ψ-|vy|ψ-⟩ at ``\\vec{k}=[k_x,k_y]``.

# Example
```jldoctest
julia> h = GappedDirac(1.0); vy_vv(h,1.0,-1.0)
0.5773502691896258
```
"""
vy_vv(h::GeneralTwoBand,kx,ky) = vμ_vv(hvec(h,kx,ky),dhdky(h,kx,ky))

export getvx_cc,getvx_cv,getvx_vc,getvx_vv,getvy_cc,getvy_cv,getvy_vc,getvy_vv

getvx_cc(h::GeneralTwoBand) = getvμ_cc(gethvec(h),getdhdx(h))
getvx_cv(h::GeneralTwoBand) = getvμ_cv(gethvec(h),getdhdx(h))
getvx_vc(h::GeneralTwoBand) = getvμ_vc(gethvec(h),getdhdx(h))
getvx_vv(h::GeneralTwoBand) = getvμ_vv(gethvec(h),getdhdx(h))
getvy_cc(h::GeneralTwoBand) = getvμ_cc(gethvec(h),getdhdy(h))
getvy_cv(h::GeneralTwoBand) = getvμ_cv(gethvec(h),getdhdy(h))
getvy_vc(h::GeneralTwoBand) = getvμ_vc(gethvec(h),getdhdy(h))
getvy_vv(h::GeneralTwoBand) = getvμ_vv(gethvec(h),getdhdy(h))
