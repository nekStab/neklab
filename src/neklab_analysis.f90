      module neklab_analysis
         use stdlib_stats_distribution_normal, only: normal => rvs_normal
         use stdlib_optval, only: optval
         use stdlib_linalg, only: diag, eye
         use stdlib_logger, only: information_level, warning_level, debug_level, error_level, all_level, success
         use LightKrylov, only: atol_dp, dp, eigs, svds, save_eigenspectrum
         use LightKrylov, only: kexpm, gmres_rdp
         use LightKrylov, only: initialize_krylov_subspace, zero_basis
         use LightKrylov, only: linear_combination, innerprod, rand_basis
         use LightKrylov, only: newton, newton_dp_opts
         use LightKrylov_Logger
         use LightKrylov_NewtonKrylov, only: dynamic_tol_dp, constant_atol_dp
         use neklab_vectors
         use neklab_linops
         use neklab_utils
         use neklab_nek_setup
         use neklab_systems
      
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
      
         private
         character(len=*), parameter, private :: this_module = 'neklab_analysis'
      
         public :: linear_stability_analysis_fixed_point
         public :: transient_growth_analysis_fixed_point
         public :: newton_fixed_point_iteration
         public :: newton_periodic_orbit
      
      contains
      
         subroutine linear_stability_analysis_fixed_point(exptA, kdim, nev, adjoint)
            type(exptA_linop), intent(inout) :: exptA
      !! Operator whose stability properties are to be investigated.
            integer, intent(in) :: kdim
      !! Maximum dimension of the Krylov subspace.
            integer, intent(in) :: nev
      !! Desired number of eigenpairs to converge.
            logical, intent(in), optional :: adjoint
      !! Whether direct or adjoint analysis should be conducted.
      
      ! Eigenvalue computation related variables.
            type(nek_dvector), allocatable :: eigvecs(:)
            complex(kind=dp), allocatable :: eigvals(:)
            real(kind=dp), allocatable :: residuals(:)
            integer :: info
      
      ! Miscellaneous.
            real(kind=dp) :: alpha
            integer :: i
            logical :: adjoint_
            character(len=3) :: file_prefix
      
      ! Optional parameters.
            if (present(adjoint)) then
               adjoint_ = adjoint
            else
               adjoint_ = .false.
            end if
      
      ! Allocate eigenvectors and initialize Krylov basis.
            allocate (eigvecs(nev)); call zero_basis(eigvecs)
      
      ! Run the eigenvalue analysis.
            call eigs(exptA, eigvecs, eigvals, residuals, info, kdim=kdim, transpose=adjoint_)
      
      ! Transform eigenspectrum to continuous-time representation.
            eigvals = log(eigvals)/exptA%tau
      
      ! Determine the file prefix.
            file_prefix = merge("adj", "dir", adjoint_)
      
      ! Save eigenspectrum to disk.
            call save_eigenspectrum(eigvals, residuals, trim(file_prefix)//"_eigenspectrum.npy")
      
      ! Export eigenfunctions to disk.
            call outpost_dnek(eigvecs(:nev), file_prefix)
      
            return
         end subroutine linear_stability_analysis_fixed_point
      
         subroutine transient_growth_analysis_fixed_point(exptA, nsv, kdim)
            type(exptA_linop), intent(inout) :: exptA
      !! Operator whose singular value decomposition needs to be computed.
            integer, intent(in) :: nsv
      !! Desired number of singular triplets.
            integer, intent(in) :: kdim
      !! Maximum dimension of the Krylov subspace in LightKrylov.
      
      ! Singular value decomposition.
            type(nek_dvector), allocatable :: U(:), V(:)
            real(kind=dp), allocatable :: S(:), residuals(:)
            integer :: info
      
      ! Miscellaneous.
            integer :: i, j
            character(len=3) :: file_prefix
      
      ! Allocate singular vectors.
            allocate (U(nsv)); call zero_basis(U)
            allocate (V(nsv)); call zero_basis(V)
      
      ! Call to LightKrylov.
            call svds(exptA, U, S, V, residuals, info, kdim=kdim)
      
      ! Save singular spectrum to disk.
            if (nid == 0) then
               open (unit=1234, file="singular_spectrum.dat")
               write (1234, *) S
               close (1234)
            end if
      
      ! Export optimal perturbations and optimal responses.
            file_prefix = "prt"; call outpost_dnek(V(:nsv), file_prefix)
            file_prefix = "rsp"; call outpost_dnek(U(:nsv), file_prefix)
      
            return
         end subroutine transient_growth_analysis_fixed_point
      
         subroutine newton_fixed_point_iteration(sys, bf, tol)
            type(nek_system), intent(inout) :: sys
      !! System for which a fixed point is sought
            type(nek_dvector), intent(inout) :: bf
      !! Initial guess for the fixed point
            real(dp), intent(inout) :: tol
      !! Absolute tolerance for the Newton solver
      
      ! Misc
            integer :: info
            type(newton_dp_opts) :: opts
      !type(gmres_dp_opts)  :: gmres_opts
            character(len=3) :: file_prefix
      
      ! Set up logging
            call logger_setup(nio=0, log_level=information_level, log_stdout=.false., log_timestamp=.true.)
      
      ! Define options for the Newton solver
            opts = newton_dp_opts(maxiter=30, ifbisect=.true.)
      
      ! Call to LightKrylov.
            call newton(sys, bf, gmres_rdp, info, tolerance=tol, options=opts, scheduler=constant_atol_dp)
      
      ! Outpost initial condition.
            file_prefix = 'nwt'
            call outpost_dnek(bf, file_prefix)
      
            return
         end subroutine newton_fixed_point_iteration
      
         subroutine newton_periodic_orbit(sys, bf, tol)
            type(nek_system_upo), intent(inout) :: sys
      !! System for which a fixed point is sought
            type(nek_ext_dvector), intent(inout) :: bf
      !! Initial guess for the fixed point
            real(dp), intent(inout) :: tol
      !! Absolute tolerance for the Newton solver
      
      ! Misc
            integer :: info
            type(newton_dp_opts) :: opts
      !type(gmres_dp_opts)  :: gmres_opts
            character(len=3) :: file_prefix
      
      ! Set up logging
            call logger_setup(nio=0, log_level=information_level, log_stdout=.false., log_timestamp=.true.)
      
      ! Define options for the Newton solver
            opts = newton_dp_opts(maxiter=30, ifbisect=.true.)
      
      ! Call to LightKrylov.
            call newton(sys, bf, gmres_rdp, info, tolerance=tol, options=opts, scheduler=constant_atol_dp)
      
      ! Outpost initial condition.
            file_prefix = 'nwt'
            call outpost_ext_dnek(bf, file_prefix)
      
            return
         end subroutine newton_periodic_orbit
      
      end module neklab_analysis
