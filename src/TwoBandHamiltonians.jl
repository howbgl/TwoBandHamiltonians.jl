module TwoBandHamiltonians

using Unitful

abstract type Hamiltonian{T} end

include("UnitScaling.jl")
include("GappedDirac.jl")

end
