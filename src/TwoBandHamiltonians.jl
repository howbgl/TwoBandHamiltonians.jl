module TwoBandHamiltonians

using LinearAlgebra
using Unitful
using StaticArrays

abstract type Hamiltonian{T} end

include("UnitScaling.jl")
include("GeneralTwoBand.jl")
include("GappedDirac.jl")

end
