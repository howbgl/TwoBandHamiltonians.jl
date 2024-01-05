# TwoBandHamiltonians

```@meta
CurrentModule = TwoBandHamiltonians
```

Documentation for [TwoBandHamiltonians](https://github.com/howbgl/TwoBandHamiltonians.jl),
a package that defines various operator matrix elements of 2x2 Hamiltonians via
multiple-dispatch.

## Idea

Consider Hamiltonians defined in k-space which are of the general form

$$
\hat{H} = \vec{h}(\vec{k})\cdot\vec{\sigma}
$$

Since this is a simple 2x2 matrix, one can express all desired operator matrix elements
in terms of $\vec{h}(\vec{k})$ and its derivatives.

Below we list the functions implemented so far.

## Eigenstate basis

At the moment all matrix elements are given in the basis of the Hamiltonian's eigenstates
describing valence- and conduction band states. **Some** matrix elements depend on the eigenstates' gauge, so we give the explicit form of the eigenstates used here in the following:

$$
\hat{H}\ket{\psi_\mu} = \epsilon_\mu\ket{\psi_\mu}
$$
with
$$
\epsilon_{\mu}(\vec{k})\equiv\epsilon_{\pm}(\vec{k})=\pm\|\vec{h}(\vec{k})\|\equiv\pm h(\vec{k})
$$

In $\sigma_z$-standardbasis ($\mu\in\{\uparrow,\downarrow\}$) the eigenstates read

$$
\braket{\mu|\psi_+} = \frac{1}{\sqrt{h(h+h_z)}}\begin{pmatrix}
    h + h_z\\h_x + i\,h
    _y
\end{pmatrix}
$$
$$
\braket{\mu|\psi_-} = \frac{1}{\sqrt{h(h+h_z)}}\begin{pmatrix}
    -h_x + i\, h_y\\h + h_z
\end{pmatrix}
$$

## Contents

```@contents
```

## Index

```@index
```
