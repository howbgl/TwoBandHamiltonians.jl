var documenterSearchIndex = {"docs":
[{"location":"reference/","page":"-","title":"-","text":"CurrentModule = TwoBandHamiltonians","category":"page"},{"location":"reference/#Reference","page":"-","title":"Reference","text":"","category":"section"},{"location":"reference/","page":"-","title":"-","text":"Modules = [TwoBandHamiltonians]","category":"page"},{"location":"reference/#TwoBandHamiltonians.GappedDirac","page":"-","title":"TwoBandHamiltonians.GappedDirac","text":"GappedDirac{T<:Real} <: GeneralTwoBand{T}\n\nHolds the parameter of a dimensionless massive Dirac Hamiltonian.\n\nThe Hamiltonian reads \n\n\\hat{H} = k_x\\sigma_x + k_y\\sigma_y + m\\sigma_z\n\nsuch that ech=k_xk_ym.\n\nSee also\n\nGeneralTwoBand\n\n\n\n\n\n","category":"type"},{"location":"reference/#TwoBandHamiltonians.GappedDiracOld","page":"-","title":"TwoBandHamiltonians.GappedDiracOld","text":"GappedDiracOld{T<:Real} <: Hamiltonian{T}\n\nLegacy hard-coded. Holds the parameter of a dimensionless massive Dirac Hamiltonian.\n\nThe Hamiltonian reads hatH =  k_xsigma_x + k_ysigma_y + Deltasigma_z such  that 2Delta is the bandgap at veck=0\n\nExamples\n\njulia> h = GappedDiracOld(0.1)\nGappedDiracOld{Float64}(0.1)\n\nSee also\n\nUnitScaling(timescale,lengthscale)\n\n\n\n\n\n","category":"type"},{"location":"reference/#TwoBandHamiltonians.GeneralTwoBand","page":"-","title":"TwoBandHamiltonians.GeneralTwoBand","text":"GeneralTwoBand{T} <: Hamiltonian{T}\n\nSupertype of all 2x2 Hamiltonians with all matrixelements via dispatch.\n\nThe idea is that all Hamiltonians of the form\n\nhatH = vech(veck)cdotvecsigma\n\ncan be diagonalized analytically and hence most desired matrixelements such as velocities or dipoles can be expressed solely through \n\nvech(veck)=h_x(veck)h_y(veck)h_z(veck)\n\nand its derivatives with respect to k_mu. Any particular Hamiltonian deriving form GeneralTwoBand{T} must then only implement vech(veck) and its derivatives.\n\nSee also\n\nGappedDirac\n\n\n\n\n\n","category":"type"},{"location":"reference/#TwoBandHamiltonians.UnitScaling","page":"-","title":"TwoBandHamiltonians.UnitScaling","text":"UnitScaling(timescale,lengthscale)\n\nRepresents a physical length- and time-scale used for non-dimensionalization of a system.\n\nExamples\n\njulia> using Unitful; us = UnitScaling(u\"1.0s\",u\"1.0m\")\nUnitScaling{Float64}(1.0e15, 1.0e9)\n\nFurther information\n\nSee here\n\n\n\n\n\n","category":"type"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = TwoBandHamiltonians","category":"page"},{"location":"#TwoBandHamiltonians","page":"Home","title":"TwoBandHamiltonians","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for TwoBandHamiltonians, a package that contains various operator matrix elements of 2x2 Hamiltonians.","category":"page"},{"location":"#Contents","page":"Home","title":"Contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"}]
}
