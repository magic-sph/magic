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
      integer :: n_fields
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
      complex(cp), pointer :: p(:,:)
   end type cmplx_pointer_wrapper
   type real_pointer_wrapper
      real(cp), pointer :: p(:,:)
   end type real_pointer_wrapper
   
   ! --------------------------------------------------------------------------------------
   type, extends(type_mpitransp_theta) :: type_mpitransp_theta_a2av
      !
      !-- Class for a single (θ,m) plane
      !   
      !   Author: Rafael Lago (MPCDF) November 2020
      !
      complex(cp), allocatable  :: m_loc_buffer(:)
      complex(cp), allocatable  :: th_loc_buffer(:)
      complex(cp), allocatable  :: fft_buffer(:,:)
      
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
      complex(cp), allocatable  :: m_loc_buffer(:)
      complex(cp), allocatable  :: th_loc_buffer(:)
      complex(cp), allocatable  :: fft_buffer(:,:)
      
      type(cmplx_pointer_wrapper), allocatable :: m_loc_ptr(:)
      type(cmplx_pointer_wrapper), allocatable :: th_loc_ptr(:)
      type(real_pointer_wrapper),  allocatable :: phi_ptr(:)
      
      logical, allocatable :: m_loc_status(:)
      logical, allocatable :: th_loc_status(:)
      
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
   
   public :: transform_m2phi_new, transform_phi2m_new, type_mpitransp_theta_a2av, type_mpitransp_theta_a2ab

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
   subroutine create_a2av(this, n_fields)
      class(type_mpitransp_theta_a2av) :: this
      integer, intent(in) :: n_fields
      integer :: irank, pos
      
      if (n_fields>1) then
         print *, "type_mpitransp_theta_a2av error: this class is only suitable for &
           & n_fields=1. Aborting..."
          stop
      end if
      this%n_fields = 1
      
      allocate(this%m_loc_buffer(n_m_loc*n_theta_max))
      allocate(this%th_loc_buffer(n_theta_loc*n_m_max))
      
      this%fftlen = max(n_m_max, n_phi_max/2+1)
      allocate(this%fft_buffer(n_theta_loc, this%fftlen))
      
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
      
      deallocate(this%m_loc_buffer)
      deallocate(this%th_loc_buffer)
      deallocate(this%fft_buffer)
      
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
            this%m_loc_buffer(pos:pos + n_t - 1) = f_mloc(l_t:u_t,j)
            pos = pos + n_t
         end do
      end do
      PERFOFF
      
#ifdef WITH_MPI
      PERFON('m2thW')
      call MPI_Alltoallv(this%m_loc_buffer, this%m2th_sendcount, this%m2th_senddispl, MPI_DEF_COMPLEX, &
           &             this%th_loc_buffer, this%m2th_recvcount, this%m2th_recvdispl, MPI_DEF_COMPLEX, &
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
            f_thetaloc(:,m_idx)=this%th_loc_buffer(pos:pos+n_theta_loc-1)
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
               this%th_loc_buffer(pos) = f_thetaloc(itheta,m+1)
               pos = pos + 1
            end do
         end do
      end do
      
      PERFOFF
      
#ifdef WITH_MPI
      PERFON('th2mW')
      call MPI_Alltoallv(this%th_loc_buffer, this%th2m_sendcount, this%th2m_senddispl, MPI_DEF_COMPLEX,   &
           &             this%m_loc_buffer, this%th2m_recvcount, this%th2m_recvdispl, MPI_DEF_COMPLEX, &
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
            f_mloc(l_t:u_t,n_m)=this%m_loc_buffer(pos:pos+n_t-1)
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
      
      call this%m2th(f_m, this%fft_buffer)
      if (n_m_max<this%fftlen) then
         this%fft_buffer(:,n_m_max+1:) = zero
      end if
      call ifft_many(this%fft_buffer, f_phi)

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

      call fft_many(f_phi, this%fft_buffer)
      call this%th2m(this%fft_buffer(1:n_theta_loc,1:n_m_max), f_m)

   end subroutine phi2m_a2av
   
   
   !
   !  TYPE_MPITRANSP_THETA_BUFFFLEX
   !  
   !  Buffered flexible implementation. Uses A2AV for at most n_fields, but handles 
   !  fewer fields as well
   !   
   !   Author: Rafael Lago (MPCDF) December 2020
   !
   !==================================================================================
   subroutine create_a2ab(this, n_fields)
      class(type_mpitransp_theta_a2ab) :: this
      integer, intent(in) :: n_fields
      integer :: irank, pos
      
      this%n_fields = n_fields
      
      allocate(this%m_loc_ptr(n_fields))
      allocate(this%th_loc_ptr(n_fields))
      allocate(this%phi_ptr(n_fields))
      
      allocate(this%m_loc_status(n_fields))
      allocate(this%th_loc_status(n_fields))
      
      allocate(this%m_loc_buffer(n_m_loc*n_theta_max*n_fields))
      allocate(this%th_loc_buffer(n_theta_loc*n_m_max*n_fields))
      
      this%fftlen = max(n_m_max, n_phi_max/2+1)
      allocate(this%fft_buffer(n_theta_loc, this%fftlen))
      
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
      
      this%name = "bufflex"
      this%direction = DIRECTION_NONE
      this%isnonblocking = .FALSE.
      
   end subroutine
   
   !----------------------------------------------------------------------------------
   subroutine destroy_a2ab(this)
      class(type_mpitransp_theta_a2ab) :: this
      
      deallocate(this%m_loc_ptr)
      deallocate(this%th_loc_ptr)
      deallocate(this%phi_ptr)
      
      deallocate(this%m_loc_status)
      deallocate(this%th_loc_status)
      
      deallocate(this%m_loc_buffer)
      deallocate(this%th_loc_buffer)
      deallocate(this%fft_buffer)
      
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
      
      this%n_buffered = this%n_buffered + 1
      
      if (this%n_buffered > this%n_fields) then
         print *, "Error in type_mpitransp_theta_a2ab%m2th: number of fields &
          &exceeded the initially allocated resources. Please, call waitall &
          &before m2th or initialize type_mpitransp_theta_a2ab with a larger&
          & buffer size"
         stop
      else if (this%n_buffered == 1) then
         this%direction = DIRECTION_M2TH
      else if (this%direction /= DIRECTION_M2TH) then 
         print *, "Error in type_mpitransp_theta_a2ab%m2th: wrong direction;&
            & do not switch directions before calling waitall"
         stop
      end if
      
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
      
      this%n_buffered = this%n_buffered + 1
      
      if (this%n_buffered > this%n_fields) then
         print *, "Error in type_mpitransp_theta_a2ab%th2m: number of fields&
            &exceeded the initially allocated resources. Please, call waitall &
            &before th2m or initialize type_mpitransp_theta_a2ab with a larger&
            & buffer size"
         stop
      else if (this%n_buffered == 1) then
         this%direction = DIRECTION_TH2M
      else if (this%direction /= DIRECTION_TH2M) then 
         print *, "Error in type_mpitransp_theta_a2ab%th2m: wrong direction;&
            & do not switch directions before calling waitall"
         stop
      end if
      
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
      
      if (this%n_buffered == this%n_fields) then
         call this%waitall()
!          print *, "Error in type_mpitransp_theta_a2ab%m2phi: number of fields &
!           &exceeded the initially allocated resources. Please, call waitall &
!           &before m2phi or initialize type_mpitransp_theta_a2ab with a larger&
!           & buffer size"
!          stop
      end if
      
      if (this%n_buffered == 0) then
         this%direction = DIRECTION_M2PHI
      else if (this%direction /= DIRECTION_M2PHI) then 
         print *, "Error in type_mpitransp_theta_a2ab%m2phi: wrong direction;&
            & do not switch directions before calling waitall"
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
      
      if (this%n_buffered == this%n_fields) then
         call this%waitall()
!          print *, "Error in type_mpitransp_theta_a2ab%phi2m: number of fields&
!             &exceeded the initially allocated resources. Please, call waitall &
!             &before phi2m or initialize type_mpitransp_theta_a2ab with a larger&
!             & buffer size"
!          stop
      end if
      
      if (this%n_buffered == 0) then
         this%direction = DIRECTION_PHI2M
      else if (this%direction /= DIRECTION_PHI2M) then 
         print *, "Error in type_mpitransp_theta_a2ab%phi2m: wrong direction;&
            & do not switch directions before calling waitall"
         stop
      end if
      
      this%n_buffered = this%n_buffered + 1
      this%m_loc_ptr(this%n_buffered)%p => f_m
      this%phi_ptr(this%n_buffered)%p   => f_phi
      
   end subroutine phi2m_start_a2ab
      
   !----------------------------------------------------------------------------------
   subroutine m2th_finish_a2ab(this)
      !   
      !   Author: Rafael Lago (MPCDF) December 2020
      !
      !-- TODO this with mpi_type to stride the data and check if performance 
      !--      improves
      !
      class(type_mpitransp_theta_a2ab) :: this
      
      !-- Local variables
      integer :: irank, pos, n_t, l_t, u_t, m_idx, n_m, n_f, n_b
      integer, allocatable :: m2th_senddispl(:), m2th_recvdispl(:)
      
      PERFON('m2thS')
      n_b = this%n_buffered
      pos = 1
      do irank=0,n_ranks_theta-1
         !-- Copy each theta chunk so that the send buffer is contiguous
         !-- TODO check performance of this; implementing this with mpi_type
         !   striding the data will probably be faster
         n_t = dist_theta(irank,0)
         l_t = dist_theta(irank,1)
         u_t = dist_theta(irank,2)
         do n_f = 1, n_b 
            do n_m=1, n_m_loc
               this%m_loc_buffer(pos:pos + n_t - 1) = this%m_loc_ptr(n_f)%p(l_t:u_t,n_m)
               pos = pos + n_t
            end do
         end do
      end do
      PERFOFF
      
      
#ifdef WITH_MPI
      PERFON('m2thW')
      call MPI_Alltoallv(this%m_loc_buffer, this%m2th_sendcount*n_b, this%m2th_senddispl*n_b, MPI_DEF_COMPLEX, &
           &             this%th_loc_buffer, this%m2th_recvcount*n_b, this%m2th_recvdispl*n_b, MPI_DEF_COMPLEX, &
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
               this%th_loc_ptr(n_f)%p(1:n_theta_loc,m_idx)=this%th_loc_buffer(pos:pos+n_theta_loc-1)
               pos = pos+n_theta_loc
            end do
         end do
      end do
      PERFOFF
      
   end subroutine m2th_finish_a2ab

   !----------------------------------------------------------------------------------
   subroutine th2m_finish_a2ab(this)
      !   
      !   Author: Rafael Lago (MPCDF) December 2020
      !
      !-- TODO this with mpi_type to stride the data
      !
      class(type_mpitransp_theta_a2ab) :: this
      
      !-- Local variables
      integer :: irank, itheta, m, pos, l_t, u_t, n_t, n_m, n_f, n_b
      
      PERFON('th2mS')
      n_b = this%n_buffered
      pos = 1
      do irank=0,n_ranks_m-1
         !-- Copy each m which belongs to the irank-th coord_r into the send buffer
         !   column-wise. That will simplify a lot things later
         !
         !@>TODO check performance of this; implementing this with mpi_type
         !  striding the data could be faster
         do n_f=1, n_b
            do n_m=1,dist_m(irank,0)
               m = dist_m(irank,n_m)/minc
               do itheta=1,n_theta_loc
                  this%th_loc_buffer(pos) = this%th_loc_ptr(n_f)%p(itheta,m+1)
                  pos = pos + 1
               end do
            end do
         end do
      end do
      PERFOFF
      
      
#ifdef WITH_MPI
      PERFON('th2mW')
      call MPI_Alltoallv(this%th_loc_buffer, this%th2m_sendcount*n_b, this%th2m_senddispl*n_b, MPI_DEF_COMPLEX,   &
           &             this%m_loc_buffer, this%th2m_recvcount*n_b, this%th2m_recvdispl*n_b, MPI_DEF_COMPLEX, &
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
               this%m_loc_ptr(n_f)%p(l_t:u_t,n_m)=this%m_loc_buffer(pos:pos+n_t-1)
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
         call m2th_finish_a2ab(this)
         
      else if (this%direction == DIRECTION_M2PHI) then
         allocate(theta_buffer(n_theta_loc,this%fftlen,this%n_buffered))
         do n_f=1, this%n_buffered
            this%th_loc_ptr(n_f)%p => theta_buffer(:,:,n_f)
         end do
         
         call m2th_finish_a2ab(this)
         
         !@>TODO have a flexible fft interface which can handle from 1 to n_buffered transforms
         do n_f=1, this%n_buffered
            theta_buffer(:,n_m_max+1:,n_f) = zero
            call ifft_many(theta_buffer(:,:,n_f), this%phi_ptr(n_f)%p)
            nullify( this%th_loc_ptr(n_f)%p )
         end do
      
         deallocate(theta_buffer)
         
      else if (this%direction == DIRECTION_TH2M) then
         call th2m_finish_a2ab(this)
      
      else if (this%direction == DIRECTION_PHI2M) then
         allocate(theta_buffer(n_theta_loc,this%fftlen,this%n_buffered))
         !@>TODO have a flexible fft interface which can handle from 1 to n_buffered transforms
         do n_f=1, this%n_buffered
            call fft_many(this%phi_ptr(n_f)%p, theta_buffer(:,:,n_f))
            this%th_loc_ptr(n_f)%p => theta_buffer(:,:,n_f)
         end do
         
         call th2m_finish_a2ab(this)
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

   !
   !   Auxiliary functions
   !  
   !   Author: Rafael Lago (MPCDF)
   !
   !==================================================================================
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
!       -- TODO: The FFT must be performed for an array with the dimensions of 
!         F_loc which may end up paded with zeroes.
!         Is there any way to tell MKL to perform a "truncated" FFT?
      call transposer%waitall()
      
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
      call transposer%waitall()

   end subroutine transform_phi2m_new
   
   
end module mpi_transpose_theta
