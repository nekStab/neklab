<img src="imgs/logo_neklab.png" style="align:center; width:512px" />

### Status

[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)
[![GitHub release](https://img.shields.io/github/release/nekStab/neklab.svg)](https://github.com/nekStab/neklab/releases/latest)
[![last-commit](https://img.shields.io/github/last-commit/nekStab/neklab)](https://github.com/nekStab/neklab/commits/main)


| **Documentation** | [link](https://nekstab.github.io/neklab/index.html) |
|:-------------:|:----:|
| **Contact**   | [jean-christophe.loiseau@ensam.eu](mailto:jean-christophe.loiseau@ensam.eu)    |


**Scope :** `neklab` is a toolbox for the massively spectral element solver `Nek5000` intended to extend its capabilities for performing large-scale linear stability and bifurcation analysis.
It makes use of [`LightKrylov`](https://github.com/nekStab/LightKrylov), our abstract linear algebra package as backend.
`neklab` can be used to perform the following analyses:

- **Fixed points and periodic orbits computations -** Uses a *time-stepper* based Newton-Krylov solver to compute particular solutions of the incompressible Navier-Stokes equations.

- **Modal stability analysis and Floquet multipliers computation -** Uses a *time-stepper* approach and a Krylov-Schur algorithm to compute the leading eigenpairs of the linearized Navier-Stokes operator (for fixed points) and of the monodromy matrix (for periodic orbits).

- **Non-modal stability analysis -** Uses an iterative Golub-Kahan factorization to compute the leading singular triplets of the exponential propagator $\mathrm{exp}(\tau A)$.

- **Optimal Time Dependent modes -** Computes the OTD modes for an arbitrary time-evolution using the algorithm from **ADD REFERENCE**.

## Description

## Capabilities

### Examples

A collection of classical examples are present in the [`examples`](https://github.com/nekStab/neklab/tree/main/examples) folder.

#### [Backward facing step](https://github.com/nekStab/neklab/tree/main/examples/back_fstep)

- [`transient_growth`](https://github.com/nekStab/neklab/tree/main/examples/back_fstep/transient_growth): Computation of the four leading singular triplets of the exponential propagator $\mathrm{exp}(\tau A)$, where $A$ is the linearized Navier-Stokes operator and $\tau$ the desired time horizon.

- [`gramian`](https://github.com/nekStab/neklab/tree/main/examples/back_fstep/gramian) (experimental): Computation of the leading eigenpairs of the observability gramian. It uses [`LightROM`](https://github.com/nekStab/LightROM) to time-march the corresponding differential Lyapunov equation until a steady state is reached using dynamical low-rank approximation.

#### [Cylinder flow](https://github.com/nekStab/neklab/tree/main/examples/cylinder)

- [`dns`](https://github.com/nekStab/neklab/tree/README/examples/cylinder/dns): Runs a standard direct numerical simulation of the two-dimensional past a circular cylinder at Reynolds 40 and 180 using vanilla `Nek5000`.

- [`newton`](https://github.com/nekStab/neklab/tree/README/examples/cylinder/newton): Illustrate how a stationary solution ($Re = 40$) and a periodic orbit ($Re = 180$) of the incompressible Navier-Stokes equations can be computed using a *time-stepper* based Newton-Krylov solver.

- [`stability`](https://github.com/nekStab/neklab/tree/main/examples/cylinder/stability): Computation of the leading left (direct) and right (adjoint) eigenvectors of the linearized Navier-Stokes operator at $Re = 50$ using `LightKrylov`'s Krylov-Schur eigenvalue algorithm.

- [`resolvent`](https://github.com/nekStab/neklab/tree/README/examples/cylinder/resolvent) (experimental): Computation of the leading singular triplets of the resolvent operator at $Re = 50$ using a *time-stepper* formulation of the problem.

#### [Poiseuille flow](https://github.com/nekStab/neklab/tree/main/examples/poiseuille)

#### [Rayleigh-Bénard convection](https://github.com/nekStab/neklab/tree/main/examples/rayBen/)

#### [Annular thermosyphon](https://github.com/nekStab/neklab/tree/main/examples/thermosyphon/)

## Installation

### Dependencies

## Contributing

## Acknowledgment

### Related projects
