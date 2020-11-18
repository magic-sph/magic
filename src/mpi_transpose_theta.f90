#include "perflib_preproc.cpp"
!----------------------------------------------------------------------------------
   !
module mpi_transpose_theta
   ! 
   ! This module contains the implementation of theta-parallel transpositions
   !

   use constants, only: zero
   use precision_mod
   use parallel_mod
   use mem_alloc
   use truncation
   use blocking, only: lm_balance, lo_map, st_map, llm, ulm
   use mpi_transp, only: type_mpitransp
   use fft, only: fft_many, ifft_many
   use LMmapping

   implicit none

   private

   ! --------------------------------------------------------------------------------------
   type, abstract, public :: type_mpitransp_theta
      !
      !-- Abstract interface for the θ<->m transposer
      !   
      !   
      !   Author: Rafael Lago (MPCDF) November 2020
      !
      character(:), ALLOCATABLE :: name
      logical :: isnonblocking
      integer :: n_fields
   contains
      procedure(create_if), deferred :: create
      procedure(destroy_if), deferred :: destroy
      procedure(m2th_if), deferred :: m2th
      procedure(th2m_if), deferred :: th2m
      procedure :: waitall => waitall_dummy
   end type type_mpitransp_theta
   
   interface 

      subroutine create_if(this, n_fields)
         import
         class(type_mpitransp_theta) :: this
         integer, intent(in) :: n_fields
      end subroutine create_if

      subroutine destroy_if(this)
         import
         class(type_mpitransp_theta) :: this
      end subroutine destroy_if

      subroutine m2th_if(this, f_mloc, f_thetaloc)
         !-- Transposition from (m_loc,θ_glb) to (θ_loc,m_glb).
         import
         class(type_mpitransp_theta) :: this
         complex(cp), intent(in)  :: f_mloc(n_theta_max, n_m_loc)
         complex(cp), intent(out) :: f_thetaloc(n_theta_loc, n_m_max)
      end subroutine m2th_if

      subroutine th2m_if(this, f_thetaloc, f_mloc)
         !-- Transposition from (θ_loc,m_glb) to (m_loc,θ_glb)
         import
         class(type_mpitransp_theta) :: this
         complex(cp), intent(in)  :: f_thetaloc(n_theta_loc, n_m_max)
         complex(cp), intent(out) :: f_mloc(n_theta_max, n_m_loc)
      end subroutine th2m_if
      
   end interface
   
   ! --------------------------------------------------------------------------------------
   type, extends(type_mpitransp_theta) :: type_mpitransp_theta_a2av
      !
      !-- Class for a single (θ,m) plane
      !   
      !   Author: Rafael Lago (MPCDF) November 2020
      !
      complex(cp), allocatable  :: sendbuf(:)
      complex(cp), allocatable  :: recvbuf(:)
      
      integer, allocatable :: m2th_sendcount(:),m2th_recvcount(:)
      integer, allocatable :: m2th_senddispl(:),m2th_recvdispl(:)
      integer, allocatable :: th2m_sendcount(:),th2m_recvcount(:)
      integer, allocatable :: th2m_senddispl(:),th2m_recvdispl(:)
   contains
      procedure :: create  => create_a2av
      procedure :: destroy => destroy_a2av
      procedure :: m2th    => m2th_a2av
      procedure :: th2m    => th2m_a2av
   end type type_mpitransp_theta_a2av

  public :: transform_m2phi_new, transform_phi2m_new, type_mpitransp_theta_a2av

contains

   !----------------------------------------------------------------------------------
   subroutine waitall_dummy(this)
      class(type_mpitransp_theta) :: this
      ! Nothing to wait...
      return
   end subroutine

   !----------------------------------------------------------------------------------
   subroutine create_a2av(this, n_fields)
      class(type_mpitransp_theta_a2av) :: this
      integer, intent(in) :: n_fields
      integer :: irank, pos
      
      this%n_fields = n_fields
      
      allocate(this%sendbuf(n_m_loc*n_theta_max))
      allocate(this%recvbuf(n_theta_loc*n_m_max))
      
      allocate(this%m2th_sendcount(0:n_ranks_theta-1))
      allocate(this%m2th_recvcount(0:n_ranks_theta-1))
      allocate(this%m2th_senddispl(0:n_ranks_theta-1))
      allocate(this%m2th_recvdispl(0:n_ranks_theta-1))
      
      this%m2th_sendcount = 0
      this%m2th_recvcount = 0
      this%m2th_senddispl = 0
      this%m2th_recvdispl = 0
      
      pos = 1
      do irank=0,n_ranks_theta-1
         this%m2th_senddispl(irank) = pos-1
         pos = pos + dist_theta(irank,0) * n_m_loc
         
         this%m2th_sendcount(irank) = pos - this%m2th_senddispl(irank) - 1
         this%m2th_recvdispl(irank) = sum(this%m2th_recvcount)
         this%m2th_recvcount(irank) = dist_m(irank,0) * n_theta_loc
      end do
      
      allocate(this%th2m_sendcount(0:n_ranks_m-1))
      allocate(this%th2m_recvcount(0:n_ranks_m-1))
      allocate(this%th2m_senddispl(0:n_ranks_m-1))
      allocate(this%th2m_recvdispl(0:n_ranks_m-1))
      
      this%th2m_sendcount = 0
      this%th2m_recvcount = 0
      this%th2m_senddispl = 0
      this%th2m_recvdispl = 0
      
      pos = 1
      do irank=0,n_ranks_m-1
         this%th2m_senddispl(irank) = pos-1
         pos = pos + n_theta_loc*dist_m(irank,0)
         
         this%th2m_sendcount(irank) = pos - this%th2m_senddispl(irank) - 1
         this%th2m_recvdispl(irank) = irank*n_m_loc*dist_theta(irank,0)
         this%th2m_recvcount(irank) = n_m_loc*dist_theta(irank,0)
      end do
      
      this%name = "a2av"
      this%isnonblocking = .FALSE.
   end subroutine
   
   !----------------------------------------------------------------------------------
   subroutine destroy_a2av(this)
      class(type_mpitransp_theta_a2av) :: this
      
      deallocate(this%sendbuf)
      deallocate(this%recvbuf)
      
      deallocate(this%m2th_sendcount)
      deallocate(this%m2th_recvcount)
      deallocate(this%m2th_senddispl)
      deallocate(this%m2th_recvdispl)
      
      deallocate(this%th2m_sendcount)
      deallocate(this%th2m_recvcount)
      deallocate(this%th2m_senddispl)
      deallocate(this%th2m_recvdispl)
      
      this%name = "[destroyed]"
      this%isnonblocking = .FALSE.
   end subroutine

   !----------------------------------------------------------------------------------
   subroutine m2th_a2av(this, f_mloc, f_thetaloc)
      !   
      !   Author: Rafael Lago (MPCDF) December 2020
      !
      !-- TODO this with mpi_type to stride the data and check if performance 
      !--      improves
      !
      class(type_mpitransp_theta_a2av) :: this
      !-- Input variable
      complex(cp), intent(in) :: f_mloc(n_theta_max, n_m_loc)

      !-- Output variable
      complex(cp), intent(out) :: f_thetaloc(n_theta_loc, n_m_max)
      
      integer :: irank, j, pos, n_t, l_t, u_t, m_idx, n_m
      
      PERFON('m2thS')
      pos = 1
      do irank=0,n_ranks_theta-1
         !-- Copy each theta chunk so that the send buffer is contiguous
         !-- TODO check performance of this; implementing this with mpi_type
         !   striding the data will probably be faster
         n_t = dist_theta(irank,0)
         l_t = dist_theta(irank,1)
         u_t = dist_theta(irank,2)
         do j=1, n_m_loc
            this%sendbuf(pos:pos + n_t - 1) = f_mloc(l_t:u_t,j)
            pos = pos + n_t
         end do
      end do
      PERFOFF
      
#ifdef WITH_MPI
      PERFON('m2thW')
      call MPI_Alltoallv(this%sendbuf, this%m2th_sendcount, this%m2th_senddispl, MPI_DEF_COMPLEX, &
           &             this%recvbuf, this%m2th_recvcount, this%m2th_recvdispl, MPI_DEF_COMPLEX, &
           &             comm_theta, ierr)
      PERFOFF
#endif
      
      !-- Now we reorder the receiver buffer. If the m distribution looks like:
      !   coord_r 0: 0, 4,  8, 12, 16
      !   coord_r 1: 1, 5,  9, 13
      !   coord_r 2: 2, 6, 10, 14
      !   coord_r 3: 3, 7, 11, 15
      !   then the columns of recvbuf are ordered as 0,4,8,12,16,1,5,9,13(...)
      !   and so forth. m_arr will contain this ordering (+1):
      PERFON('m2thS')
      do irank=0,n_ranks_theta-1
         pos = this%m2th_recvdispl(irank)+1
         do n_m=1,dist_m(irank,0)
            m_idx = dist_m(irank,n_m)/minc+1
            f_thetaloc(:,m_idx)=this%recvbuf(pos:pos+n_theta_loc-1)
            pos = pos+n_theta_loc
         end do
      end do
      PERFOFF
      
   end subroutine m2th_a2av

   !----------------------------------------------------------------------------------
   subroutine th2m_a2av(this, f_thetaloc, f_mloc)
      !   
      !   Author: Rafael Lago (MPCDF) December 2020
      !
      !-- TODO this with mpi_type to stride the data
      !
      class(type_mpitransp_theta_a2av) :: this
      !-- Input variable
      complex(cp), intent(in) :: f_thetaloc(n_theta_loc, n_m_max)

      !-- Output variable
      complex(cp), intent(out) :: f_mloc(n_theta_max, n_m_loc)
      
      !-- Local variables
      integer :: irank, j, itheta, m, pos, l_t, u_t, n_t, n_m
      
      PERFON('th2mS')
      pos = 1
      do irank=0,n_ranks_m-1
         !-- Copy each m which belongs to the irank-th coord_r into the send buffer
         !   column-wise. That will simplify a lot things later
         !
         !@>TODO check performance of this; implementing this with mpi_type
         !  striding the data could be faster
         do j=1,dist_m(irank,0)
            m = dist_m(irank,j)/minc
            do itheta=1,n_theta_loc
               this%recvbuf(pos) = f_thetaloc(itheta,m+1)
               pos = pos + 1
            end do
         end do
      end do
      
      PERFOFF
      
#ifdef WITH_MPI
      PERFON('th2mW')
      print *, "yes, claling!!!"
      call MPI_Alltoallv(this%recvbuf, this%th2m_sendcount, this%th2m_senddispl, MPI_DEF_COMPLEX,   &
           &             this%recvbuf, this%th2m_recvcount, this%th2m_recvdispl, MPI_DEF_COMPLEX, &
           &             comm_m, ierr)
      PERFOFF
#endif

      PERFON('th2mS')
      do irank=0,n_ranks_theta-1
         pos = this%th2m_recvdispl(irank)+1
         l_t=dist_theta(irank,1)
         u_t=dist_theta(irank,2)
         n_t=dist_theta(irank,0)
         do n_m=1,n_m_loc
            f_mloc(l_t:u_t,n_m)=this%sendbuf(pos:pos+n_t-1)
            pos = pos+n_t
         end do
      end do
      PERFOFF
      
   end subroutine th2m_a2av
   
   !----------------------------------------------------------------------------------
   subroutine transform_m2phi_new(transposer, fL, f)
      !-- Transforms from (θ,m) space into (φ,θ) space including transpositions 
      !   and FFT. 
      !   
      !   Author: Rafael Lago (MPCDF) April 2020
      !   TODO: some functions might requires this transformation for multiple
      !      fields at once (e.g. vector transform). In that case, a non-blocking 
      !      transpose_theta_m would make more sense. This would make this 
      !      function obsolete.
      !   TODO: there is a lot of room for immprovement here (e.g. in-place, 
      !     use just one intermediate buffer, vector transform, etc)
      !
      
      !-- Input variables
      class(type_mpitransp_theta), intent(inout) :: transposer
      complex(cp), intent(inout) :: fL(n_theta_max,n_m_loc)
      
      !-- Output variables
      real(cp),    intent(out)   :: f(n_theta_loc,n_phi_max)
      
      !-- Local variables
      complex(cp) :: lF(n_theta_loc,n_m_max)
      complex(cp) :: Ff(n_theta_loc,n_phi_max/2+1)
      
      call transposer%m2th(fL, lF)
      !-- TODO: The FFT must be performed for an array with the dimensions of 
      !   F_loc which may end up paded with zeroes.
      !   Is there any way to tell MKL to perform a "truncated" FFT?
      Ff = zero
      Ff(1:n_theta_loc,1:n_m_max) = lF
      
      call ifft_many(Ff, f)

   end subroutine transform_m2phi_new
!----------------------------------------------------------------------------------
   subroutine transform_phi2m_new(transposer, f, fL)
      !-- Transforms from (φ,θ) space into (θ,m) space including transpositions 
      !   and FFT. 
      !   
      !   Author: Rafael Lago (MPCDF) April 2020
      !   TODO: some functions might requires this transformation for multiple
      !      fields at once (e.g. vector transform). In that case, a non-blocking 
      !      transpose_theta_m would make more sense. This would make this 
      !      function obsolete.
      !   TODO: there is a lot of room for immprovement here (e.g. in-place, 
      !     use just one intermediate buffer, vector transform, etc)
      !
      
      !-- Input variables
      class(type_mpitransp_theta), intent(inout) :: transposer
      real(cp),    intent(inout) :: f(n_theta_loc,n_phi_max)
      
      !-- Output variables
      complex(cp), intent(out) :: fL(n_theta_max,n_m_loc)
      
      !-- Local variables
      complex(cp) :: lF(n_theta_loc,n_m_max)
      complex(cp) :: Ff(n_theta_loc,n_phi_max/2+1)

      call fft_many(f, Ff)
      lF(1:n_theta_loc,1:n_m_max) = Ff(1:n_theta_loc,1:n_m_max)
      call transposer%th2m(lF, fL)

   end subroutine transform_phi2m_new
   
   
end module mpi_transpose_theta
