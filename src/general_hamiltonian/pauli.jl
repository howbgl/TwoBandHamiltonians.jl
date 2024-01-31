
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

export σvec_cv,σvec_vc
σvec_cv(h::SVector{3,<:Number}) = SA[σx_cv(h),σy_cv(h),σz_cv(h)]
σvec_vc(h::SVector{3,<:Number}) = SA[σx_vc(h),σy_vc(h),σz_vc(h)]
