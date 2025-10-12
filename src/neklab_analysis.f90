      module neklab_analysis
         use stdlib_stats_distribution_normal, only: normal => rvs_normal
         use stdlib_optval, only: optval
         use stdlib_logger, only: information_level, warning_level, debug_level, error_level, all_level, success
         use LightKrylov, only: dp, eigs, svds, save_eigenspectrum
         use LightKrylov, only: gmres_rdp
         use LightKrylov, only: initialize_krylov_subspace, orthonormalize_basis, zero_basis, rand_basis
         use LightKrylov, only: linear_combination, innerprod
         use LightKrylov, only: newton, newton_dp_opts, GMRES_dp_opts
         use LightKrylov_Logger
         use LightKrylov_Timing, only: timer => global_lightkrylov_timer
         use LightKrylov_AbstractVectors, only: abstract_vector_rdp
         use LightKrylov_AbstractLinops, only: abstract_exptA_linop_rdp
         use LightKrylov_AbstractSystems, only: abstract_system_rdp
         use LightKrylov_NewtonKrylov, only: newton_dp_metadata
         use neklab_vectors
         use neklab_linops
         use neklab_utils
         use neklab_nek_setup
         use neklab_otd
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
         public :: otd_analysis
      
      contains
      
         subroutine linear_stability_analysis_fixed_point(exptA, kdim, nev, adjoint, X0)
            class(abstract_exptA_linop_rdp), intent(inout) :: exptA
      !! Operator whose stability properties are to be investigated.
            integer, intent(in) :: kdim
      !! Maximum dimension of the Krylov subspace.
            integer, intent(in) :: nev
      !! Desired number of eigenpairs to converge.
            logical, intent(in), optional :: adjoint
      !! Whether direct or adjoint analysis should be conducted.
		type(nek_dvector), optional, intent(in) :: X0
      !! Initial guess for the eigenvectors
      
      ! Eigenvalue computation related variables.
            character(len=*), parameter :: this_procedure = 'stability_main'
            type(nek_dvector), allocatable :: eigvecs(:)
            complex(kind=dp), allocatable :: eigvals(:)
            real(kind=dp), allocatable :: residuals(:)
            integer :: info
      
      ! Miscellaneous.
            real(kind=dp) :: alpha
            integer :: i
            logical :: adjoint_
            character(len=3) :: file_prefix
      
      ! Set up logging
            call logger_setup(nio=0, log_level=information_level, log_stdout=.false., log_timestamp=.true.)
      ! Initialize timers
            call timer%initialize()
            call timer%add_timer('Linear Stability Fixed Point', start=.true.)
      
      ! Optional parameters.
            if (present(adjoint)) then
               adjoint_ = adjoint
            else
               adjoint_ = .false.
            end if
      
      ! Allocate eigenvectors and initialize Krylov basis.
            allocate (eigvecs(nev)); call zero_basis(eigvecs)
      
      ! Run the eigenvalue analysis.
            call eigs(exptA, eigvecs, eigvals, residuals, info, x0=X0, kdim=kdim,
     &                  transpose=adjoint_, write_intermediate=.true.)
      
      ! Transform eigenspectrum to continuous-time representation.
            eigvals = log(eigvals)/exptA%tau
      
      ! Determine the file prefix.
            file_prefix = merge("adj", "dir", adjoint_)
      
      ! Save eigenspectrum to disk.
            call save_eigenspectrum(eigvals, residuals, trim(file_prefix)//"_eigenspectrum.npy")
      
      ! Export eigenfunctions to disk.
            call outpost_dnek(eigvecs(:nev), file_prefix)

            call nek_log_message('Exiting eigenvalue computation.', this_module, this_procedure)

      ! Finalize exptA timings
            call exptA%finalize_timer()
      ! Finalize timing
            call logger_setup(logfile='lightkrylov_tmr.log', nio=0, log_level=warning_level, log_stdout=.false., log_timestamp=.true.)
            call timer%finalize()
      
		call nek_log_message('Exiting eigenvalue computation.', this_module, this_procedure)
      
         end subroutine linear_stability_analysis_fixed_point
      
         subroutine transient_growth_analysis_fixed_point(exptA, nsv, kdim)
            class(abstract_exptA_linop_rdp), intent(inout) :: exptA
      !! Operator whose singular value decomposition needs to be computed.
            integer, intent(in) :: nsv
      !! Desired number of singular triplets.
            integer, intent(in) :: kdim
      !! Maximum dimension of the Krylov subspace in LightKrylov.
      
      ! Singular value decomposition.
            character(len=*), parameter :: this_procedure = 'transient_growth_main'
            type(nek_dvector), allocatable :: U(:), V(:)
            real(kind=dp), allocatable :: S(:), residuals(:)
            integer :: info
      
      ! Miscellaneous.
            integer :: i, j
            character(len=3) :: file_prefix
      
      ! Set up logging
            call logger_setup(nio=0, log_level=information_level, log_stdout=.false., log_timestamp=.true.)
      ! Initialize timers
            call timer%initialize()
            call timer%add_timer('Transient Growth Fixed Point', start=.true.)

      ! Allocate singular vectors.
            allocate (U(nsv)); call initialize_krylov_subspace(U)
            allocate (V(nsv)); call initialize_krylov_subspace(V)
      
      ! Call to LightKrylov.
            call svds(exptA, U, S, V, residuals, info, kdim=kdim, write_intermediate=.true.)
      
      ! Save singular spectrum to disk.
            if (nid == 0) then
               open (unit=1234, file="singular_spectrum.dat")
               write (1234, *) S
               close (1234)
            end if
      
      ! Export optimal perturbations and optimal responses.
            file_prefix = "prt"; call outpost_dnek(V(:nsv), file_prefix)
            file_prefix = "rsp"; call outpost_dnek(U(:nsv), file_prefix)

            call nek_log_message('Exiting transient growth computation.', this_module, this_procedure)
      
         end subroutine transient_growth_analysis_fixed_point
      
         subroutine newton_fixed_point_iteration(sys, bf, tol, tol_mode, input_is_fixed_point)
            class(abstract_system_rdp), intent(inout) :: sys
      !! System for which a fixed point is sought
            class(abstract_vector_rdp), intent(inout) :: bf
      !! Initial guess for the fixed point
            real(dp), intent(inout) :: tol
      !! Absolute tolerance for the Newton solver
            integer, optional, intent(in) :: tol_mode
      !! constant or dynamic tolerances?
            logical, optional, intent(out) :: input_is_fixed_point
      !! optional flag to return whether the intial condition is a fixed point (and no new solution is computed)
      
      ! Misc
            character(len=*), parameter :: this_procedure = 'newton_main'
            integer :: info, tol_mode_
            type(newton_dp_opts) :: newton_opts
            type(gmres_dp_opts) :: GMRES_opts
            character(len=3) :: file_prefix
            type(newton_dp_metadata) :: meta
      
		tol_mode_ = optval(tol_mode, 1)

      	call nek_log_message('Starting newton iteration.', this_module, this_procedure)
      
      ! Define options for the Newton solver
         newton_opts = newton_dp_opts(maxiter=40, ifbisect=.false.)
      ! Define options for the GMRES solver
         GMRES_opts = GMRES_dp_opts(kdim = 30, maxiter = 10)
      
      ! Call to LightKrylov.
            if (tol_mode_ == 1) then
               call newton(sys, bf, GMRES_rdp, info, atol=tol, options=newton_opts, linear_solver_options=GMRES_opts, scheduler=nek_constant_tol, meta=meta)
            else
		   call newton(sys, bf, GMRES_rdp, info, atol=tol, options=newton_opts, linear_solver_options=GMRES_opts, scheduler=nek_dynamic_tol, meta=meta)
		end if
      
      ! Outpost initial condition.
            file_prefix = 'nwt'
            call set_fldindex(file_prefix, 1)
            call outpost_nek(bf, file_prefix)
            if (present(input_is_fixed_point)) then
               input_is_fixed_point = meta%input_is_fixed_point
            end if

		call nek_log_message('Exiting newton iteration.', this_module, this_procedure)
      
         end subroutine newton_fixed_point_iteration
      
         subroutine otd_analysis(OTD, opts_)
            type(nek_otd), intent(inout) :: OTD
            type(otd_opts), optional, intent(in) :: opts_
            type(otd_opts) :: opts
      ! internal
            character(len=*), parameter :: this_procedure = 'OTD_main'
            real(dp), dimension(:), allocatable :: sigma
            real(dp), dimension(:, :), allocatable :: Lr, Phi, svec, G
            complex(dp), dimension(:), allocatable :: lambda
            complex(dp), dimension(:, :), allocatable :: eigvec
            type(nek_dvector), allocatable :: Lu(:)
      ! Misc
            integer :: i, j, r, log_level
            character(len=3) :: file_prefix
            character(len=128) :: msg
      
            if (present(opts_)) then
               opts = opts_
            else
               opts = otd_opts()
            end if
      
      ! Set up logging
            call logger_setup(nio=0, log_level=information_level, log_stdout=.false., log_timestamp=.true.)
      
      ! initialize OTD structure
            call OTD%init(opts)
      
      ! Allocate memory
            r = OTD%r
            allocate (sigma(r), svec(r, r)); sigma = 0.0_dp; svec = 0.0_dp
            allocate (lambda(r), eigvec(r, r)); lambda = 0.0_dp; eigvec = 0.0_dp
            allocate (Lr(r, r), Phi(r, r)); Lr = 0.0_dp; Phi = 0.0_dp
            allocate (Lu(r), source=OTD%baseflow); call zero_basis(Lu)
      
      ! Intgrate the nonlinear equations forward
            time = 0.0_dp
            do istep = 1, nsteps
               call nek_advance()
               if (istep >= opts%startstep) then

      ! load perturbations
                  do i = 1, r
                     call nek2vec(OTD%basis(i), vxp(:, i:i), vyp(:, i:i), vzp(:, i:i), prp(:, i:i), tp(:, :, i:i))
                  end do

      ! orthonormalize
                  if ((istep <= opts%startstep + 10) .or.
     &               mod(istep, opts%orthostep) == 0 .or.
     &               mod(istep, opts%printstep) == 0 .or.
     &               mod(istep, opts%iostep) == 0) then
      
                     if (log_level <= debug_level) then
                        allocate (G(r, r)); G = 0.0_dp
                        G = innerprod(OTD%basis, OTD%basis)
                        write (msg, '(A,I5,A,*(1X,E10.3))') 'Step ', istep, ': norm.  err pre: ',  (G(i,i) - 1.0_dp, i=1, r)
                        call nek_log_information(msg, this_module, this_procedure)
                        write (msg, '(A,I5,A,*(1X,E10.3))') 'Step ', istep, ': ortho. err pre: ', ((G(i,j), j=i+1, r), i=1, r)
                        call nek_log_information(msg, this_module, this_procedure)
                     end if
         
                     call orthonormalize_basis(OTD%basis)

                     if (log_level <= debug_level) then
                        write (msg, '(A,I5,A,*(1X,E10.3))') 'Step ', istep, ': norm.  err post:',  (G(i,i) - 1.0_dp, i=1, r)
                        call nek_log_debug(msg, this_module, this_procedure)
                        write (msg, '(A,I5,A,*(1X,E10.3))') 'Step ', istep, ': ortho. err post:', ((G(i,j), j=i+1, r), i=1, r)
                        call nek_log_debug(msg, this_module, this_procedure)
                     end if
                  end if

      ! compute Lu
                  do i = 1, r
                     if (opts%trans) then
                        call OTD%apply_rmatvec(OTD%basis(i), Lu(i))
                     else
                        call OTD%apply_matvec(OTD%basis(i), Lu(i))
                     end if
                  end do

      ! compute reduced operator
                  Lr = innerprod(OTD%basis, Lu)
      
                  Phi = 0.0_dp
                  do i = 1, r
                     do j = i + 1, r
                        Phi(i, j) = Lr(i, j)
                        Phi(j, 1) = -Lr(i, j)
                     end do
                  end do

      ! output projected modes
                  if (mod(istep, opts%printstep) == 0) then
                     call OTD%spectral_analysis(Lr, sigma, svec, lambda, eigvec, ifprint=.true.)
                  end if

      ! at the end of the step we copy data back to nek2vec
                  do i = 1, r
                     call vec2nek(vxp(:, i:i), vyp(:, i:i), vzp(:, i:i), prp(:, i:i), tp(:, :, i:i), OTD%basis(i))
                  end do

      ! project basis vectors and output modes
                  if (mod(istep, opts%iostep) == 0) then
                     if (mod(istep, opts%printstep) /= 0) then
                        call OTD%spectral_analysis(Lr, sigma, svec, lambda, eigvec, ifprint=.false.)
                     end if
                     call OTD%outpost_OTDmodes(eigvec)
                  end if

      ! output basis vectors
                  if (mod(istep, opts%iorststep) == 0) then
                     write (file_prefix, '(A)') 'rst'
                     call outpost_dnek(OTD%basis, file_prefix)
                  end if

      ! set the forcing
                  call OTD%generate_forcing(Lr, Phi)

               end if ! istep >= otd_startstep
            end do ! istep ... nsteps

            call nek_log_message('Exiting OTD computation.', this_module, this_procedure)
         end subroutine otd_analysis
      
         end module neklab_analysis
