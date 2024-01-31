module TwoBandHamiltonians

using LinearAlgebra
using Unitful
using StaticArrays

abstract type Hamiltonian{T} end

include("UnitScaling.jl")
include("general_hamiltonian/GeneralTwoBand.jl")
include("GappedDirac.jl")
include("HexWarpDirac.jl")

end
