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
      ! TODO: make FFT plans belong to this class and allow them to have their dedicated
      ! plan for computing FFTs
      !
      character(:), ALLOCATABLE :: name
      logical :: isnonblocking
      integer :: max_buff
   contains
      procedure(create_if), deferred :: create
      procedure(destroy_if), deferred :: destroy
      procedure(m2th_if),  deferred :: m2th
      procedure(th2m_if),  deferred :: th2m
      procedure(m2phi_if), deferred :: m2phi
      procedure(phi2m_if), deferred :: phi2m
      procedure :: waitall => waitall_dummy
   end type type_mpitransp_theta
   
   interface 

      subroutine create_if(this, max_buff)
         import
         class(type_mpitransp_theta) :: this
         integer, intent(in) :: max_buff
      end subroutine create_if

      subroutine destroy_if(this)
         import
         class(type_mpitransp_theta) :: this
      end subroutine destroy_if

      subroutine m2th_if(this, f_mloc, f_thetaloc)
         !-- Transposition from (m_loc,θ_glb) to (θ_loc,m_glb).
         import
         class(type_mpitransp_theta) :: this
         complex(cp), target, intent(in)  :: f_mloc(n_theta_max, n_m_loc)
         complex(cp), target, intent(out) :: f_thetaloc(n_theta_loc, n_m_max)
      end subroutine m2th_if

      subroutine th2m_if(this, f_thetaloc, f_mloc)
         !-- Transposition from (θ_loc,m_glb) to (m_loc,θ_glb)
         import
         class(type_mpitransp_theta) :: this
         complex(cp), target, intent(in)  :: f_thetaloc(n_theta_loc, n_m_max)
         complex(cp), target, intent(out) :: f_mloc(n_theta_max, n_m_loc)
      end subroutine th2m_if
      
      subroutine m2phi_if(this, f_m, f_phi)
         !-- Transform from (θ_glb,m_loc) to (θ_loc,φ_glb)
         import
         class(type_mpitransp_theta) :: this
         complex(cp), target, intent(inout) :: f_m(n_theta_max,n_m_loc)
         real(cp),    target, intent(out)   :: f_phi(n_theta_loc,n_phi_max)
      end subroutine m2phi_if

      subroutine phi2m_if(this, f_phi, f_m)
         !-- Transform from (θ_loc,φ_glb) to (θ_glb,m_loc)
         import
         class(type_mpitransp_theta) :: this
         real(cp),   target, intent(inout) :: f_phi(n_theta_loc,n_phi_max)
         complex(cp),target, intent(out)   :: f_m(n_theta_max,n_m_loc)
      end subroutine phi2m_if
      
   end interface
   
   !   Workarounds the lack of array of pointers in fortran
   ! --------------------------------------------------------------------------------------
   type cmplx_pointer_wrapper
      complex(cp), contiguous, pointer :: p(:,:)
   end type cmplx_pointer_wrapper
   type real_pointer_wrapper
      real(cp), contiguous, pointer :: p(:,:)
   end type real_pointer_wrapper
   
   ! --------------------------------------------------------------------------------------
   type, extends(type_mpitransp_theta) :: type_mpitransp_theta_a2av
      !
      !-- Class for a single (θ,m) plane
      !   
      !   Author: Rafael Lago (MPCDF) November 2020
      !
      integer :: fftlen
      
      integer, allocatable :: m2th_sendcount(:),m2th_recvcount(:)
      integer, allocatable :: m2th_senddispl(:),m2th_recvdispl(:)
      integer, allocatable :: th2m_sendcount(:),th2m_recvcount(:)
      integer, allocatable :: th2m_senddispl(:),th2m_recvdispl(:)
   contains
      procedure :: create  => create_a2av
      procedure :: destroy => destroy_a2av
      procedure :: m2th    => m2th_a2av
      procedure :: th2m    => th2m_a2av
      procedure :: m2phi   => m2phi_a2av
      procedure :: phi2m   => phi2m_a2av
   end type type_mpitransp_theta_a2av
   
   ! --------------------------------------------------------------------------------------
   type, extends(type_mpitransp_theta) :: type_mpitransp_theta_a2ab
      !
      !-- Class for a multiple (θ,m) planes
      !   
      !   Author: Rafael Lago (MPCDF) December 2020
      !
      type(cmplx_pointer_wrapper), allocatable :: m_loc_ptr(:)
      type(cmplx_pointer_wrapper), allocatable :: th_loc_ptr(:)
      type(real_pointer_wrapper),  allocatable :: phi_ptr(:)
      
      integer :: n_buffered, direction, fftlen
      
      integer, allocatable :: m2th_sendcount(:),m2th_recvcount(:)
      integer, allocatable :: m2th_senddispl(:),m2th_recvdispl(:)
      integer, allocatable :: th2m_sendcount(:),th2m_recvcount(:)
      integer, allocatable :: th2m_senddispl(:),th2m_recvdispl(:)
   contains
      procedure :: create  => create_a2ab
      procedure :: destroy => destroy_a2ab
      procedure :: m2th    => m2th_start_a2ab
      procedure :: th2m    => th2m_start_a2ab
      procedure :: m2phi   => m2phi_start_a2ab
      procedure :: phi2m   => phi2m_start_a2ab
      procedure :: waitall => waitall_a2ab
   end type type_mpitransp_theta_a2ab
   
   integer, parameter :: DIRECTION_NONE = 0
   integer, parameter :: DIRECTION_M2TH = -1
   integer, parameter :: DIRECTION_TH2M =  1
   integer, parameter :: DIRECTION_M2PHI = -2
   integer, parameter :: DIRECTION_PHI2M =  2
   
   public :: type_mpitransp_theta_a2av, type_mpitransp_theta_a2ab

contains

   !----------------------------------------------------------------------------------
   subroutine waitall_dummy(this)
      class(type_mpitransp_theta) :: this
      ! Nothing to wait...
      return
   end subroutine

   !
   !  TYPE_MPITRANSP_THETA_A2AV
   !  
   !  Simplest implementation (A2AV, 1 double buffer only)
   !   
   !   
   !   Author: Rafael Lago (MPCDF) December 2020
   !
   !==================================================================================
   subroutine create_a2av(this, max_buff)
      class(type_mpitransp_theta_a2av) :: this
      integer, intent(in) :: max_buff
      integer :: irank, pos
      
      if (max_buff>1) then
         print *, "type_mpitransp_theta_a2av error: this class is only suitable for &
           & max_buff=1. Aborting..."
          stop
      end if
      this%max_buff = 1
      
      this%fftlen = max(n_m_max, n_phi_max/2+1)
      
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
      complex(cp), target, intent(in) :: f_mloc(n_theta_max, n_m_loc)

      !-- Output variable
      complex(cp), target, intent(out) :: f_thetaloc(n_theta_loc, n_m_max)
      
      integer :: irank, j, pos, n_t, l_t, u_t, m_idx, n_m
      complex(cp)  :: m_loc_buffer(n_m_loc*n_theta_max)
      complex(cp)  :: th_loc_buffer(n_theta_loc*n_m_max)
      
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
            m_loc_buffer(pos:pos + n_t - 1) = f_mloc(l_t:u_t,j)
            pos = pos + n_t
         end do
      end do
      PERFOFF
      
      call mpi_barrier(comm_theta, ierr)
#ifdef WITH_MPI
      PERFON('m2thW')
      call MPI_Alltoallv(m_loc_buffer, this%m2th_sendcount, this%m2th_senddispl, MPI_DEF_COMPLEX, &
           &             th_loc_buffer, this%m2th_recvcount, this%m2th_recvdispl, MPI_DEF_COMPLEX, &
           &             comm_theta, ierr)
      PERFOFF
#endif
      
      !-- Now we reorder the receiver buffer. If the m distribution looks like:
      !   coord_r 0: 0, 4,  8, 12, 16
      !   coord_r 1: 1, 5,  9, 13
      !   coord_r 2: 2, 6, 10, 14
      !   coord_r 3: 3, 7, 11, 15
      !   then the columns of th_loc_buffer are ordered as 0,4,8,12,16,1,5,9,13(...)
      !   and so forth. m_arr will contain this ordering (+1):
      PERFON('m2thS')
      do irank=0,n_ranks_theta-1
         pos = this%m2th_recvdispl(irank)+1
         do n_m=1,dist_m(irank,0)
            m_idx = dist_m(irank,n_m)/minc+1
            f_thetaloc(:,m_idx)=th_loc_buffer(pos:pos+n_theta_loc-1)
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
      complex(cp), target, intent(in) :: f_thetaloc(n_theta_loc, n_m_max)

      !-- Output variable
      complex(cp), target, intent(out) :: f_mloc(n_theta_max, n_m_loc)
      
      !-- Local variables
      integer :: irank, j, itheta, m, pos, l_t, u_t, n_t, n_m
      complex(cp)  :: m_loc_buffer(n_m_loc*n_theta_max)
      complex(cp)  :: th_loc_buffer(n_theta_loc*n_m_max)
      
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
            th_loc_buffer(pos:pos+n_theta_loc-1) = f_thetaloc(1:n_theta_loc,m+1)
            pos = pos + n_theta_loc
         end do
      end do
      PERFOFF
      
      call mpi_barrier(comm_theta, ierr)
#ifdef WITH_MPI
      PERFON('th2mW')
      call MPI_Alltoallv(th_loc_buffer, this%th2m_sendcount, this%th2m_senddispl, MPI_DEF_COMPLEX,   &
           &             m_loc_buffer, this%th2m_recvcount, this%th2m_recvdispl, MPI_DEF_COMPLEX, &
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
            f_mloc(l_t:u_t,n_m)=m_loc_buffer(pos:pos+n_t-1)
            pos = pos+n_t
         end do
      end do
      PERFOFF
      
   end subroutine th2m_a2av
   
   !----------------------------------------------------------------------------------
   subroutine m2phi_a2av(this, f_m, f_phi)
      !-- Transform from (θ_glb,m_loc) to (θ_loc,φ_glb) space including transpositions 
      !   and FFT. 
      !   
      !   Author: Rafael Lago (MPCDF) December 2020
      !-- Input variables
      class(type_mpitransp_theta_a2av) :: this
      complex(cp), target, intent(inout) :: f_m(n_theta_max,n_m_loc)
      
      !-- Output variables
      real(cp),    target, intent(out)   :: f_phi(n_theta_loc,n_phi_max)
      complex(cp), allocatable :: fft_buffer(:,:)
      
      allocate(fft_buffer(n_theta_loc, this%fftlen))
      call this%m2th(f_m, fft_buffer)
      fft_buffer(:,n_m_max+1:) = zero
      call ifft_many(fft_buffer, f_phi)
      deallocate(fft_buffer)

   end subroutine m2phi_a2av
   
   !----------------------------------------------------------------------------------
   subroutine phi2m_a2av(this, f_phi, f_m)
      !-- Transform from (θ_loc,φ_glb) to (θ_glb,m_loc) space including transpositions 
      !   and FFT. 
      !   
      !   Author: Rafael Lago (MPCDF) December 2020
      !-- Input variables
      class(type_mpitransp_theta_a2av) :: this
      real(cp),    target, intent(inout) :: f_phi(n_theta_loc,n_phi_max)
      
      !-- Output variables
      complex(cp), target, intent(out)   :: f_m(n_theta_max,n_m_loc)
      complex(cp), allocatable :: fft_buffer(:,:)

      allocate(fft_buffer(n_theta_loc, this%fftlen))
      call fft_many(f_phi, fft_buffer)
      call this%th2m(fft_buffer(1:n_theta_loc,1:n_m_max), f_m)
      deallocate(fft_buffer)

   end subroutine phi2m_a2av
   
   !
   !  TYPE_MPITRANSP_THETA_BUFFFLEX
   !  
   !  Buffered flexible implementation. Uses A2AV for at most max_buff, but handles 
   !  fewer fields as well
   !   
   !   Author: Rafael Lago (MPCDF) December 2020
   !
   !==================================================================================
   subroutine create_a2ab(this, max_buff)
      class(type_mpitransp_theta_a2ab) :: this
      integer, intent(in) :: max_buff
      integer :: irank, pos
      
      this%max_buff = max_buff
      
      allocate(this%m_loc_ptr(max_buff))
      allocate(this%th_loc_ptr(max_buff))
      allocate(this%phi_ptr(max_buff))
      
      this%fftlen = max(n_m_max, n_phi_max/2+1)
      this%n_buffered = 0
      
      allocate(this%m2th_sendcount(0:n_ranks_theta-1))
      allocate(this%m2th_recvcount(0:n_ranks_theta-1))
      allocate(this%m2th_senddispl(0:n_ranks_theta-1))
      allocate(this%m2th_recvdispl(0:n_ranks_theta-1))
      
      this%m2th_sendcount = 0
      this%m2th_recvcount = 0
      this%m2th_senddispl = 0
      this%m2th_recvdispl = 0
      
      ! Counts and displacements are computed for a single field
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
      
      if (n_phi_max<n_m_max) then
         print *, "Error creating type_mpitransp_theta_a2ab: n_phi_max<n_m_max &
         & is not allowed for this transposition mode! n_phi_max=",n_phi_max, " and&
         & n_m_max=", n_m_max
         stop
      end if
      
      this%name = "a2ab"
      this%direction = DIRECTION_NONE
      this%isnonblocking = .FALSE.
      
   end subroutine
   
   !----------------------------------------------------------------------------------
   subroutine destroy_a2ab(this)
      class(type_mpitransp_theta_a2ab) :: this
      integer :: i
      
      do i=1, this%max_buff
         if (associated(this%m_loc_ptr(i)%p)) nullify(this%m_loc_ptr(i)%p)
         if (associated(this%th_loc_ptr(i)%p)) nullify(this%th_loc_ptr(i)%p)
         if (associated(this%phi_ptr(i)%p)) nullify(this%phi_ptr(i)%p)
      end do
      deallocate(this%m_loc_ptr)
      deallocate(this%th_loc_ptr)
      deallocate(this%phi_ptr)
      
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
   subroutine m2th_start_a2ab(this, f_mloc, f_thetaloc)
      !   
      !   Author: Rafael Lago (MPCDF) December 2020
      !
      class(type_mpitransp_theta_a2ab) :: this
      !-- Input variable
      complex(cp), target, intent(in) :: f_mloc(n_theta_max, n_m_loc)

      !-- Output variable
      complex(cp), target, intent(out) :: f_thetaloc(n_theta_loc, n_m_max)
      
      
      if (this%n_buffered == this%max_buff) then
         call this%waitall()
      end if
      
      if (this%n_buffered == 0) then
         this%direction = DIRECTION_M2TH
      else if (this%direction /= DIRECTION_M2TH) then 
         print *, "Error in type_mpitransp_theta_a2ab%m2th: do not switch directions before calling waitall! "//&
         & " Wrong direction: ", this%direction
         ! Try to force a segfault, to get a traceback...
         deallocate(this%m_loc_ptr(1)%p)
         this%m_loc_ptr(1)%p(-4, -1) = -5
         stop
      end if
      
      this%n_buffered = this%n_buffered + 1
      this%m_loc_ptr(this%n_buffered)%p  => f_mloc
      this%th_loc_ptr(this%n_buffered)%p => f_thetaloc
      
   end subroutine m2th_start_a2ab
   
   !----------------------------------------------------------------------------------
   subroutine th2m_start_a2ab(this, f_thetaloc, f_mloc)
      !   
      !   Author: Rafael Lago (MPCDF) December 2020
      !
      class(type_mpitransp_theta_a2ab) :: this
      !-- Input variable
      complex(cp), target, intent(in) :: f_thetaloc(n_theta_loc, n_m_max)

      !-- Output variable
      complex(cp), target, intent(out) :: f_mloc(n_theta_max, n_m_loc)
      
      
      if (this%n_buffered == this%max_buff) then
         call this%waitall()
      end if
      
      if (this%n_buffered == 0) then
         this%direction = DIRECTION_TH2M
      else if (this%direction /= DIRECTION_TH2M) then 
         print *, "Error in type_mpitransp_theta_a2ab%th2m: do not switch directions before calling waitall! "//&
         & " Wrong direction: ", this%direction
         ! Try to force a segfault, to get a traceback...
         deallocate(this%m_loc_ptr(1)%p)
         this%m_loc_ptr(1)%p(-4, -1) = -5
         stop
      end if
      
      this%n_buffered = this%n_buffered + 1
      this%m_loc_ptr(this%n_buffered)%p  => f_mloc
      this%th_loc_ptr(this%n_buffered)%p => f_thetaloc
      
   end subroutine th2m_start_a2ab
   
   
   !----------------------------------------------------------------------------------
   subroutine m2phi_start_a2ab(this, f_m, f_phi)
      !-- Transform from (θ_glb,m_loc) to (θ_loc,φ_glb) space including transpositions 
      !   and FFT. 
      !   
      !   Author: Rafael Lago (MPCDF) December 2020
      !-- Input variables
      class(type_mpitransp_theta_a2ab) :: this
      complex(cp), target, intent(inout) :: f_m(n_theta_max,n_m_loc)
      
      !-- Output variables
      real(cp), target, intent(out)   :: f_phi(n_theta_loc,n_phi_max)
      
      if (this%n_buffered == this%max_buff) then
         call this%waitall()
      end if
      
      if (this%n_buffered == 0) then
         this%direction = DIRECTION_M2PHI
      else if (this%direction /= DIRECTION_M2PHI) then 
         print *, "Error in type_mpitransp_theta_a2ab%m2phi: do not switch directions before calling waitall! "//&
         & " Wrong direction: ", this%direction
         ! Try to force a segfault, to get a traceback...
         deallocate(this%m_loc_ptr(1)%p)
         this%m_loc_ptr(1)%p(-4, -1) = -5
         stop
      end if
      
      this%n_buffered = this%n_buffered + 1
      this%m_loc_ptr(this%n_buffered)%p => f_m
      this%phi_ptr(this%n_buffered)%p   => f_phi

   end subroutine m2phi_start_a2ab
   
   !----------------------------------------------------------------------------------
   subroutine phi2m_start_a2ab(this, f_phi, f_m)
      !-- Transform from (θ_loc,φ_glb) to (θ_glb,m_loc) space including transpositions 
      !   and FFT. 
      !   
      !   Author: Rafael Lago (MPCDF) December 2020
      !-- Input variables
      class(type_mpitransp_theta_a2ab) :: this
      real(cp), target, intent(inout) :: f_phi(n_theta_loc,n_phi_max)
      
      !-- Output variables
      complex(cp), target, intent(out) :: f_m(n_theta_max,n_m_loc)
      
      if (this%n_buffered == this%max_buff) then
         call this%waitall()
      end if
      
      if (this%n_buffered == 0) then
         this%direction = DIRECTION_PHI2M
      else if (this%direction /= DIRECTION_PHI2M) then 
         print *, "Error in type_mpitransp_theta_a2ab%phi2m: do not switch directions before calling waitall! "//&
         & " Wrong direction: ", this%direction
         
         ! Try to force a segfault, to get a traceback...
         deallocate(this%m_loc_ptr(1)%p)
         this%m_loc_ptr(1)%p(-4, -1) = -5
         stop
      end if
      
      this%n_buffered = this%n_buffered + 1
      this%m_loc_ptr(this%n_buffered)%p => f_m
      this%phi_ptr(this%n_buffered)%p   => f_phi
      
   end subroutine phi2m_start_a2ab
      
   !----------------------------------------------------------------------------------
   subroutine m2th_finish_a2ab(this, n_b)
      !   
      !   Author: Rafael Lago (MPCDF) December 2020
      !
      !-- TODO this with mpi_type to stride the data and check if performance 
      !--      improves
      !
      class(type_mpitransp_theta_a2ab) :: this
      integer, intent(in) :: n_b
      
      !-- Local variables
      integer :: irank, pos, n_t, l_t, u_t, m_idx, n_m, n_f, ierr
      complex(cp)  :: m_loc_buffer(n_m_loc*n_theta_max*n_b)
      complex(cp)  :: th_loc_buffer(n_theta_loc*n_m_max*n_b)
      
      PERFON('m2thS')
      pos = 1
      do irank=0,n_ranks_theta-1
         !-- Copy each theta chunk so that the send buffer is contiguous
         !@>TODO check performance of this; implementing this with mpi_type striding the data could be faster
         !@>TODO copy this data over during the _start call; it may prevent cache repagination
         n_t = dist_theta(irank,0)
         l_t = dist_theta(irank,1)
         u_t = dist_theta(irank,2)
         do n_f = 1, n_b 
            do n_m=1, n_m_loc
               m_loc_buffer(pos:pos + n_t - 1) = this%m_loc_ptr(n_f)%p(l_t:u_t,n_m)
               pos = pos + n_t
            end do
         end do
      end do
      PERFOFF
      
      
      call mpi_barrier(comm_theta, ierr)
#ifdef WITH_MPI
      PERFON('m2thW')
      call MPI_Alltoallv(m_loc_buffer, this%m2th_sendcount*n_b, this%m2th_senddispl*n_b, MPI_DEF_COMPLEX, &
           &             th_loc_buffer, this%m2th_recvcount*n_b, this%m2th_recvdispl*n_b, MPI_DEF_COMPLEX, &
           &             comm_theta, ierr)
      PERFOFF
#endif
      
      !-- Now we reorder the receiver buffer
      PERFON('m2thS')
      do irank=0,n_ranks_theta-1
         pos = this%m2th_recvdispl(irank)*n_b+1
         do n_f=1,n_b
            do n_m=1,dist_m(irank,0)
               m_idx = dist_m(irank,n_m)/minc+1
               this%th_loc_ptr(n_f)%p(1:n_theta_loc,m_idx)=th_loc_buffer(pos:pos+n_theta_loc-1)
               pos = pos+n_theta_loc
            end do
         end do
      end do
      PERFOFF
      
   end subroutine m2th_finish_a2ab

   !----------------------------------------------------------------------------------
   subroutine th2m_finish_a2ab(this, n_b)
      !   
      !   Author: Rafael Lago (MPCDF) December 2020
      !
      !-- TODO this with mpi_type to stride the data
      !
      class(type_mpitransp_theta_a2ab) :: this
      integer, intent(in) :: n_b
      
      !-- Local variables
      integer :: irank, itheta, m, pos, l_t, u_t, n_t, n_m, n_f, ierr
      complex(cp)  :: m_loc_buffer(n_m_loc*n_theta_max*n_b)
      complex(cp)  :: th_loc_buffer(n_theta_loc*n_m_max*n_b)
      
      PERFON('th2mS')
      pos = 1
      do irank=0,n_ranks_m-1
         !-- Copy each m which belongs to the irank-th coord_r into the send buffer
         !   column-wise. That will simplify a lot things later
         !
         !@>TODO check performance of this; implementing this with mpi_type striding the data could be faster
         !@>TODO copy this data over during the _start call; it may prevent cache repagination
         do n_f=1, n_b
            do n_m=1,dist_m(irank,0)
               m = dist_m(irank,n_m)/minc
               th_loc_buffer(pos:pos+n_theta_loc-1) = this%th_loc_ptr(n_f)%p(1:n_theta_loc,m+1)
               pos = pos + n_theta_loc
            end do
         end do
      end do
      PERFOFF

      call mpi_barrier(comm_theta, ierr)
#ifdef WITH_MPI
      PERFON('th2mW')
      call MPI_Alltoallv(th_loc_buffer, this%th2m_sendcount*n_b, this%th2m_senddispl*n_b, MPI_DEF_COMPLEX,   &
           &             m_loc_buffer, this%th2m_recvcount*n_b, this%th2m_recvdispl*n_b, MPI_DEF_COMPLEX, &
           &             comm_m, ierr)
      PERFOFF
#endif

      PERFON('th2mS')
      do irank=0,n_ranks_theta-1
         pos = this%th2m_recvdispl(irank)*n_b+1
         l_t=dist_theta(irank,1)
         u_t=dist_theta(irank,2)
         n_t=dist_theta(irank,0)
         do n_f=1, n_b
            do n_m=1,n_m_loc
               this%m_loc_ptr(n_f)%p(l_t:u_t,n_m)=m_loc_buffer(pos:pos+n_t-1)
               pos = pos + n_t
            end do
         end do
      end do
      PERFOFF
      
   end subroutine th2m_finish_a2ab
   
   !----------------------------------------------------------------------------------
   subroutine waitall_a2ab(this)
      class(type_mpitransp_theta_a2ab) :: this
      
      complex(cp), target, allocatable :: theta_buffer(:,:,:)
      integer :: n_f
      
      if (this%n_buffered == 0) return
      
      if (this%direction == DIRECTION_M2TH) then
         call m2th_finish_a2ab(this, this%n_buffered)
         
      else if (this%direction == DIRECTION_M2PHI) then
         allocate(theta_buffer(n_theta_loc,this%fftlen,this%n_buffered))
         do n_f=1, this%n_buffered
            this%th_loc_ptr(n_f)%p => theta_buffer(:,:,n_f)
         end do
         
         call m2th_finish_a2ab(this, this%n_buffered)
         
         !@>TODO have a flexible fft interface which can handle from 1 to n_buffered transforms
         do n_f=1, this%n_buffered
            theta_buffer(:,n_m_max+1:,n_f) = zero
            call ifft_many(theta_buffer(:,:,n_f), this%phi_ptr(n_f)%p)
            nullify( this%th_loc_ptr(n_f)%p )
         end do
      
         deallocate(theta_buffer)
         
      else if (this%direction == DIRECTION_TH2M) then
         call th2m_finish_a2ab(this, this%n_buffered)
      
      else if (this%direction == DIRECTION_PHI2M) then
         allocate(theta_buffer(n_theta_loc,this%fftlen,this%n_buffered))
         !@>TODO have a flexible fft interface which can handle from 1 to n_buffered transforms
         do n_f=1, this%n_buffered
            call fft_many(this%phi_ptr(n_f)%p, theta_buffer(:,:,n_f))
            this%th_loc_ptr(n_f)%p => theta_buffer(:,:,n_f)
         end do
         
         call th2m_finish_a2ab(this, this%n_buffered)
         do n_f=1, this%n_buffered
            nullify( this%th_loc_ptr(n_f)%p )
         end do
         deallocate(theta_buffer)
      else
         print *, "Unknow direction in type_mpitransp_theta_a2ab%waitall: ",  \
            this%direction, ". Aborting..."
         stop
      end if
      
      this%n_buffered = 0
      this%direction = DIRECTION_NONE
   end subroutine
   
end module mpi_transpose_theta
