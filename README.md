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

#### Backward facing step

- [`transient_growth`](https://github.com/nekStab/neklab/tree/main/examples/back_fstep/transient_growth): Computation of the four leading singular triplets of the exponential propagator $\mathrm{exp}(\tau A)$, where $A$ is the linearized Navier-Stokes operator and $\tau$ the desired time horizon.

- [`gramian`](https://github.com/nekStab/neklab/tree/main/examples/back_fstep/gramian) (experimental): Computation of the leading eigenpairs of the observability/controlability gramian. It computes a low-rank approximation in the frequency domain using a *time-stepper* formulation of the resolvent problem.

#### Cylinder flow

- [`dns`](https://github.com/nekStab/neklab/tree/README/examples/cylinder/dns): Runs a standard direct numerical simulation of the two-dimensional past a circular cylinder at Reynolds 40 and 180 using vanilla `Nek5000`.

- [`newton`](https://github.com/nekStab/neklab/tree/README/examples/cylinder/newton): Illustrate how a stationary solution ($Re = 40$) and a periodic orbit ($Re = 180$) of the incompressible Navier-Stokes equations can be computed using a *time-stepper* based Newton-Krylov solver.

- [`stability`](https://github.com/nekStab/neklab/tree/main/examples/cylinder/stability): Computation of the leading left (direct) and right (adjoint) eigenvectors of the linearized Navier-Stokes operator at $Re = 50$ using `LightKrylov`'s Krylov-Schur eigenvalue algorithm.

- [`resolvent`](https://github.com/nekStab/neklab/tree/README/examples/cylinder/resolvent) (experimental): Computation of the leading singular triplets of the resolvent operator at $Re = 50$ using a *time-stepper* formulation of the problem.

#### Poiseuille flow

- [`stability`](https://github.com/nekStab/neklab/tree/main/examples/poiseuille/stability): Illustrate how to use `neklab` to compute the leading eigenpairs of the linearized Navier-Stokes operator for a 2D-periodic (streamwise or spanwise) flow.

- [`OTD_steady`](https://github.com/nekStab/neklab/tree/main/examples/poiseuille/OTD_steady): Computation of the leading Optimal Time Dependent (OTD) modes.

#### Rayleigh-Bénard convection

- [`baseflow`](https://github.com/nekStab/neklab/tree/README/examples/rayBen/baseflow): Computation of the unstable steady solution for the Rayleigh-Bénard convection, including coupling between the momentum and temprature equations.

#### Annular thermosyphon

- [`baseflow`](https://github.com/nekStab/neklab/tree/main/examples/thermosyphon/baseflow): Computation of the unstable steady solution for the Rayleigh-Bénard convection in an annular thermosyphon, including coupling between the momentum and temprature equations.


## Installation

## Contributing

`neklab` is currently developped and maintained by a team of two:

- [Jean-Christophe Loiseau](https://loiseaujc.github.io/): Assistant Professor of Applied maths and Fluid dynamics at DynFluid, Arts et Métiers Institute of Technology, Paris, France.
- [Simon Kern](https://github.com/Simkern/): PhD in Fluid dynamics (KTH, Sweden, 2023) and currently postdoctoral researcher at DynFluid.

[Ricardo Frantz](https://github.com/ricardofrantz), PhD in Fluid Dynamics (Arts et Métiers, France, 2022) and now working in Switzerland, has also immensely contributed to the initial development of `neklab`.

Contributions are more than welcomed! More information can be found in the following pages:

- [Guidelines](https://github.com/nekStab/neklab/blob/main/CONTRIBUTING.md)
- [Issues](https://github.com/nekStab/neklab/issues)
- [Workflow](https://github.com/nekStab/neklab/blob/main/WORKFLOW.md)
- [Style guide](https://github.com/nekStab/neklab/blob/main/STYLE_GUIDE.md)
- [Code of conduct](https://github.com/nekStab/neklab/blob/main/CODE_OF_CONDUCT.md)
- [Licence](https://github.com/nekStab/neklab/blob/main/LICENSE)

## Acknowledgment

The development of `neklab` is part of an on-going research project funded by [Agence Nationale pour la Recherche](https://anr.fr/en/) (ANR) under the grant agreement ANR-22-CE46-0008. The project started in January 2023 and will run until December 2026.
We are also very grateful to the [fortran-lang](https://fortran-lang.org/) community and the maintainers of [`stdlib`](https://github.com/fortran-lang/stdlib), in particular to @perazz, @jalvesz and @jvdp1 for their awesome work on the `stdlib_linalg` module which greatly simplified the developlement of `LightKrylov`, our linear algebra backend.

### Related projects

- [`KTH Framework`](https://github.com/KTH-Nek5000/KTH_Framework): Another toolbox for `Nek5000` developed at KTH (Sweden) with a scope similar to `neklab`.
- [`NekLab Docker`](https://github.com/eduardomartini/neklab_docker): A docker container for `neklab` by @eduardomartini.
