
export UnitScaling
"""
    UnitScaling(timescale,lengthscale)

Represents a physical length- and time-scale used for non-dimensionalization of a system.

# Examples
```jldoctest
julia> using Unitful; us = UnitScaling(u"1.0s",u"1.0m")
UnitScaling{Float64}(1.0e15, 1.0e9)
```
"""
struct UnitScaling{T<:Real}
    timescale::T
    lengthscale::T
end
function UnitScaling(timescale::Unitful.Time,lengthscale::Unitful.Length) 
    return UnitScaling(ustrip(u"fs",timescale),ustrip(u"nm",lengthscale))
end


function getparams(us::UnitScaling)
    return (timescale=Quantity(us.timescale,u"fs"),
            lengthscale=Quantity(us.lengthscale,u"nm"))
end
gettimescale(us::UnitScaling)       = Quantity(us.timescale,u"fs")
getlengthscale(us::UnitScaling)     = Quantity(us.lengthscale,u"nm")

function energySI(en,us::UnitScaling)
    tc,lc = getparams(us)
    return uconvert(u"meV",en*Unitful.ħ/tc)
end
function electricfieldSI(field,us::UnitScaling)
    tc,lc   = getparams(us)
    e   = uconvert(u"C",1u"eV"/1u"V")
    return uconvert(u"MV/cm",field*Unitful.ħ/(e*tc*lc))
end
function timeSI(time,us::UnitScaling)
    tc,lc = getparams(us)
    return uconvert(u"fs",time*tc)
end
function lengthSI(length,us::UnitScaling)
    tc,lc = getparams(us)
    return uconvert(u"Å",length*lc)
end
function frequencySI(ν,us::UnitScaling)
    tc,lc = getparams(us)
    return uconvert(u"THz",ν/tc)
end
function velocitySI(v,us::UnitScaling)
    tc,lc = getparams(us)
    return uconvert(u"m/s",v*lc/tc)
end
function wavenumberSI(k,us::UnitScaling)
    tc,lc = getparams(us)
    return uconvert(u"Å^-1",k/lc)
end