# JOSS LightKrylov - neklab example

This branch of the larger [neklab][https://github.com/nekStab/neklab] repository contains the source code for the numerical examples presented in the reference below showcasing the integration of the [LightKrylov][https://github.com/nekStab/LightKrylov] abstract linear algebra toolbox with the massively parallel open-source CFD code [Nek5000][https://github.com/Nek5000/Nek5000].

## Usage
All examples contain necessary files for code compilation and running the case except those, that can be easily recreated using standard Nek tools. Moreover, to reduce a number of binary files in the repository, we do not include multiple copies of the mesh file `1cyl.re2`. This mesh file must be linked or copied into the run directory prior to execution.

Files provided with each example.
* setup source file `###.usr`.
* runtime parameters file `###.par`.
* required `SIZE` file containing definitions of static arrays dimensions.
* Solution fields for different Reynolds numbers `BF_1cyl*.fld` to be used either as initial conditions for the Newton-GMRES fixed-point iteration or baseflow field for the eigenvalue computation. 

To compile the code:
* ensure that `Nek5000` has been successfully cloned and the `LightKrylov`-specific changes have been executed. This is most easily acheived by running the `Nek5000_setup.sh` script in the neklab root directory (Note: the script must be executable).
* ensure that `LightKrylov` has been successfully cloned and installed on the system. This is most easily acheived by running the `LightKrylov_setup.sh` script in the neklab root directory (Note: the script must be executable). We recommend running the test suite to check that everything works correctly.
* navigate to the folder of the case of interest and build the case using script `makeneklab` found in the `app` directory of the neklab root directory.

To run the case:
* navigate to the folder of the case of interest.
* copy or link the mesh file `1cyl.re2` from the example root directory.
* generate the processor map file `1cyl.map` using the tool `genmap` distributed together with `Nek5000`.
* run the code using the parallel executable provided in `Nek5000/bin/`.

# Example list
* `Newton-GMRES`: Newton-Krylov iteration to find the steady fixed-point of the nonlinear Navier-Stokes equations for the supercritical cylinder flow at Re=100.
* `eigs`: Krylov-Schur eigenvalue iteration to find the leading eigenpair of the exponential propagator of the linearized Navier-Stokes operator around the steady fixed-point of the supercritical cylinder flow at Re=100. 

# References
Kern et al. (2025). LightKrylov: Lightweight implementation of Krylov subspace techniques in modern Fortran. Journal of Open Source Software, 2
¿VOL? (¿ISSUE?), ¿PAGE? https://doi.org/10.xxxxxx/draft.