      module neklab_analysis
      use LightKrylov, only: dp, eigs, save_eigenspectrum, initialize_krylov_subspace
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
      complex(kind=dp) , allocatable :: eigvals(:)
      real(kind=dp)    , allocatable :: residuals(:)
      integer :: info
      
      ! Miscellaneous.
      real(kind=dp) :: alpha
      integer :: i
      logical :: transpose
      character(len=3) :: file_prefix
      
      ! Optional parameters.
      if (present(adjoint)) then
      transpose = adjoint
      else
      transpose = .false.
      end if
      
      ! Allocate eigenvectors and initialize Krylov basis.
      allocate (eigvecs(nev)); call initialize_krylov_subspace(eigvecs)
      
      ! Run the eigenvalue analysis.
      call eigs(exptA, eigvecs, eigvals, residuals, info, kdim=kdim, transpose=adjoint)
      
      ! Transform eigenspectrum to continuous-time representation.
      eigvals = log(eigvals)/exptA%tau
      
      ! Determine the file prefix.
      file_prefix = merge("adj", "dir", transpose)

      ! Save eigenspectrum to disk.
      call save_eigenspectrum(eigvals, residuals, trim(file_prefix) // "_eigenspectrum.npy")
      
      ! Export eigenfunctions to disk.
      do i = 1, nev
      call outpost_vec(eigvecs(i), file_prefix)
      end do
      
      return
      end subroutine linear_stability_analysis_fixed_point
      
      end module neklab_analysis
