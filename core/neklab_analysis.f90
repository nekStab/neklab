      module neklab_analysis
      use LightKrylov, only: wp, eigs, save_eigenspectrum, initialize_krylov_subspace
      use neklab_vectors
      use neklab_linops
      implicit none
      include "SIZE"
      include "TOTAL"
      include "ADJOINT"
      
      private
      
      public :: linear_stability_analysis_fixed_point
      
      contains
      
      subroutine linear_stability_analysis_fixed_point(exptA, kdim, nev, adjoint)
      type(exponential_propagator), intent(in) :: exptA
      !! Operator whose stability properties are to be investigated.
      integer, intent(in) :: kdim
      !! Maximum dimension of the Krylov subspace.
      integer, intent(in) :: nev
      !! Desired number of eigenpairs to converge.
      logical, intent(in), optional :: adjoint
      !! Whether direct or adjoint analysis should be conducted.
      
      ! Eigenvalue computation related variables.
      type(nek_dvector), allocatable :: eigvecs(:)
      complex(kind=wp) :: eigvals(kdim)
      real(kind=wp) :: residuals(kdim)
      integer :: info
      
      ! Miscellaneous.
      real(kind=wp) :: alpha
      integer :: i
      logical :: transpose
      
      ! Optional parameters.
      if (present(adjoint)) then
      transpose = adjoint
      else
      transpose = .false.
      end if
      
      ! Allocate eigenvectors and initialize Krylov basis.
      allocate (eigvecs(kdim+1)); call initialize_krylov_subspace(eigvecs)
      call eigvecs(1)%rand(); alpha = eigvecs(1)%norm()
      call eigvecs(1)%scal(1.0_wp/alpha)
      
      ! Run the eigenvalue analysis.
      call eigs(exptA, eigvecs, eigvals, residuals, info, nev=nev)
      
      ! Transform eigenspectrum to continuous-time representation.
      eigvals = log(eigvals)/exptA%tau
      
      ! Save eigenspectrum to disk.
      call save_eigenspectrum(eigvals%re, eigvals%im, residuals, "eigenspectrum.npy")
      
      ! Export eigenfunctions to disk.
      do i = 1, nev
      call outpost_vec(eigvecs(i), "eig")
      end do
      
      return
      end subroutine linear_stability_analysis_fixed_point
      
      end module neklab_analysis
