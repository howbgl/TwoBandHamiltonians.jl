```@meta
CurrentModule = TwoBandHamiltonians
```

# TwoBandHamiltonians

Documentation for [TwoBandHamiltonians](https://github.com/howbgl/TwoBandHamiltonians.jl),
a package that defines various operator matrix elements of 2x2 Hamiltonians via
multiple-dispatch.

## Idea

Consider Hamiltonians defined in k-space which are of the general form

```math
\hat{H} = \vec{h}(\vec{k})\cdot\vec{\sigma}
```

Since this is a simple 2x2 matrix, one can express all desired operator matrix elements
in terms of ``\\vec{h}(\\vec{k})`` and its derivatives.

Below we list the functions implemented so far.

## Eigenstate basis

At the moment all matrix elements are given in the basis of the Hamiltonian's eigenstates
describing valence- and conduction band states.
**Some** matrix elements depend on the eigenstates' gauge, so we give the explicit form of
the eigenstates used here in the following:

```math
\hat{H}\ket{\psi_\mu} = \epsilon_\mu\ket{\psi_\mu}
```
with ``\\epsilon_{\\mu}(\\vec{k})=\\epsilon_{\\pm}(\\vec{k})=\\pm\\norm{\\vec{h}(\\vec{k})}.

## Contents

```@contents
```

## Index

```@index
```
