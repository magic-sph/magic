module mod_mpiatoav_new
   ! 
   !-- First alltoallv implementation
   !>@TODO skip the reordering of the local data into the buffer; we can copy it
   ! into the the output array directly
   !>@TODO use MPI_TYPES to create a strided array! We will still need the 
   !  reordering into the buffer, but we can force n_fields to remain the last
   !  dimension thus reducing considerably the memory reshuffling
   !  

   use precision_mod
   use parallel_mod
   use mem_alloc
   use truncation
   use blocking, only: lm_balance, lo_map, st_map, llm, ulm
   use mpi_transp, only: type_mpitransp
   use fft
   use LMmapping

   implicit none

   private

   type, extends(type_mpitransp) :: type_mpiatoav_new
      integer :: ncount
      complex(cp), pointer :: send_buf(:,:,:)
      complex(cp), pointer :: recv_buf(:)
      integer, allocatable :: send_count(:)
      integer, allocatable :: recv_count(:)
      integer, allocatable :: send_disp(:)
      integer, allocatable :: recv_disp(:)
      integer :: send_sum, recv_sum
      logical :: initialized = .false.

   contains
      procedure :: create_comm => create_mlo_atoav_new
      procedure :: destroy_comm => destroy_mlo_atoav_new
      procedure :: transp_lm2r => transp_lm2r_atoav_new_dummy
      procedure :: transp_r2lm => transp_r2lm_atoav_new_dummy
      procedure :: transp_lm2r_dist => transp_lm2r_atoav_new_dist
      procedure :: transp_r2lm_dist => transp_r2lm_atoav_new_dist
      procedure :: reorder_rloc2buffer
      procedure :: reorder_buffer2lmloc
      procedure :: reorder_buffer2rloc
      procedure :: reorder_lmloc2buffer
   end type type_mpiatoav_new
   
   logical, private :: mlo_atoav_new_initialized = .false.
   integer, private :: n_atoav_new_obj = 0
   
   integer, private, allocatable :: lmr2buf(:)
   integer, private, allocatable :: buf2mlo(:)
   integer, private, allocatable :: send_lm_count(:)
   integer, private, allocatable :: recv_mlo_count(:)
   
   public :: type_mpiatoav_new


contains

   !-- Transposition from (n_lm_loc,nRStart) to (n_mlo_loc,n_r_max).
   !   
   !   Author: Rafael Lago (MPCDF) April 2020
   subroutine initialize_mlo_atoav_new
      !-- Initialize the MPI types for the transposition from ML to Radial
      !   
      !   This transposition assumes that all m's are kept within the comm_r
      !   Thus, comm_mlo is not used, just comm_r instead.
      !   
      !   Author: Rafael Lago (MPCDF) May 2020
      !    
      !-- Local variables
      integer :: i, j, l, m, mlo, icoord_mo, icoord_molo, icoord_lo, ierr
      integer :: disp
      
      !-- Help variables for computing lmr2buf
      integer :: dest(n_lm_loc) ! holds rank_R which will receive the corresponding lm point
      integer :: counter(0:n_ranks_r-1)
      
      if (mlo_atoav_new_initialized) return
      
      allocate(send_lm_count(0:n_ranks_lo-1))
      allocate(recv_mlo_count(0:n_ranks_lo-1))
      
      ! First we count only (l,m) points here
      send_lm_count = 0
      do i=1,n_lm_loc
         l = map_dist_st%lm2l(i)
         m = map_dist_st%lm2m(i)
         
         icoord_lo   = map_mlo%ml2coord_lo(m,l)
         dest(i) = icoord_lo
         send_lm_count(icoord_lo) = send_lm_count(icoord_lo) + 1
      end do
      
      ! Builds the reordering into the send_buf
      counter = 1
      allocate(lmr2buf(n_lm_loc))
      allocate(buf2mlo(n_mlo_loc))
      do i=1,n_lm_loc
         icoord_lo = dest(i)
         if (icoord_lo>0) then
            lmr2buf(i) = SUM(send_lm_count(0:icoord_lo-1)) + counter(icoord_lo)
         else
            lmr2buf(i) = counter(icoord_lo)
         end if
         if (icoord_lo == coord_r) then
            l = map_dist_st%lm2l(i)
            m = map_dist_st%lm2m(i)
            mlo = map_mlo%ml2i(m,l)
            buf2mlo(counter(icoord_lo)) = mlo
         end if
         counter(icoord_lo) = counter(icoord_lo) + 1
      end do
      
      recv_mlo_count = dist_r(:,0)
      
      ! Some few sanity checks
      ! Once iimplement the bypassing of reodering local data, most of these
      ! checks will be invalid!
      if (minval(buf2mlo)/=1)  then
         write (*,'(A,I0,A)') ' * Wrong buf2mlo! min=', minval(buf2mlo), ', should be 1'
         call abortRun(' * initialize_mlo_atoav_new: failed sanity check test')
      end if
      if (maxval(buf2mlo)/=n_mlo_loc)  then
         write (*,'(A,I0,A,I0)') ' * Wrong buf2mlo! max=', maxval(buf2mlo), ', should be ',n_mlo_loc
         call abortRun(' * initialize_mlo_atoav_new: failed sanity check test')
      end if
      if (minval(lmr2buf)/=1)  then
         write (*,'(A,I0,A)') ' * Wrong lmr2buf! min=', minval(lmr2buf), ', should be 1'
         call abortRun(' * initialize_mlo_atoav_new: failed sanity check test')
      end if
      if (maxval(lmr2buf)/=n_lm_loc)  then
         write (*,'(A,I0,A,I0)') ' * Wrong lmr2buf! max=', maxval(lmr2buf), ', should be ',n_lm_loc
         call abortRun(' * initialize_mlo_atoav_new: failed sanity check test')
      end if
      if (sum(recv_mlo_count)/=n_r_max)  then
         write (*,'(A,I0,A,I0)') ' * Wrong recv_mlo_count! sum=', sum(recv_mlo_count), ', should be ',n_r_max
         call abortRun(' * initialize_mlo_atoav_new: failed sanity check test')
      end if
      if (sum(send_lm_count)/=n_lm_loc)  then
         write (*,'(A,I0,A,I0)') ' * Wrong send_lm_count! sum=', sum(send_lm_count), ', should be ',n_lm_loc
         call abortRun(' * initialize_mlo_atoav_new: failed sanity check test')
      end if
      
      mlo_atoav_new_initialized = .true.
   end subroutine initialize_mlo_atoav_new
   
   !----------------------------------------------------------------------------
   subroutine reorder_rloc2buffer(this, arr_Rloc)
   !-- Reorder input variable into send buffer
   !
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoav_new) :: this
      complex(cp), intent(in)  :: arr_Rloc(n_lm_loc, nRstart:nRstop, *)
      integer :: i
      
      !>@TODO optimize the order of this loop
      do i=1,n_lm_loc
         this%send_buf(1:n_r_loc, 1:this%n_fields, lmr2buf(i)) = &
            arr_Rloc(i,nRstart:nRstop, 1:this%n_fields)
      end do
   end subroutine
   
   !----------------------------------------------------------------------------
   subroutine reorder_buffer2rloc(this, arr_Rloc)
   !-- Reorder receive buffer into output variable
   !   OBS: in the lm->r transposition, the buffer names are flipped.
   !     recv_buffer is actually send_buffer and vice versa.
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoav_new) :: this
      complex(cp), intent(out)  :: arr_Rloc(n_lm_loc, nRstart:nRstop, *)
      integer :: i
      
      !>@TODO optimize the order of this loop
      do i=1,n_lm_loc
         arr_Rloc(i,nRstart:nRstop, 1:this%n_fields) = &
         &  this%send_buf(1:n_r_loc, 1:this%n_fields, lmr2buf(i))
      end do
   end subroutine
   
   !----------------------------------------------------------------------------
   subroutine reorder_buffer2lmloc(this, arr_LMloc)
   !-- Reorder receive buffer onto output variable
   !
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoav_new) :: this
      complex(cp), intent(out)  :: arr_LMloc(n_mlo_loc, n_r_max, *)
      integer :: icoord_lo, lr, ur, lb, ub, i, n

      !>@TODO optimize the order of this loop
      do icoord_lo=0,n_ranks_r-1
         lr = dist_r(icoord_lo,1)
         ur = dist_r(icoord_lo,2)
         n  = dist_r(icoord_lo,0)
         lb = this%recv_disp(icoord_lo)+1
         ub = lb + n*this%n_fields-1
         do i=1, n_mlo_loc
            arr_LMloc(buf2mlo(i), lr:ur, 1:this%n_fields) = reshape(this%recv_buf(lb:ub), (/n, this%n_fields/))
            lb = ub + 1
            ub = lb + n*this%n_fields-1
         end do
      end do
   end subroutine
   
   !----------------------------------------------------------------------------
   subroutine reorder_lmloc2buffer(this, arr_LMloc)
   !-- Reorder input variable into send buffer
   !   OBS: in the lm->r transposition, the buffer names are flipped.
   !     recv_buffer is actually send_buffer and vice versa.
   !
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoav_new) :: this
      complex(cp), intent(in)  :: arr_LMloc(n_mlo_loc, n_r_max, *)
      integer :: icoord_lo, lr, ur, lb, ub, i, n
      
      !>@TODO optimize the order of this loop
      do icoord_lo=0,n_ranks_r-1
         lr = dist_r(icoord_lo,1)
         ur = dist_r(icoord_lo,2)
         n  = dist_r(icoord_lo,0)
         lb = this%recv_disp(icoord_lo)+1
         ub = lb + n*this%n_fields-1
         do i=1, n_mlo_loc
            this%recv_buf(lb:ub) = reshape(arr_LMloc(buf2mlo(i), lr:ur, 1:this%n_fields), (/n*this%n_fields/))
            lb = ub + 1
            ub = lb + n*this%n_fields-1
         end do
      end do
   end subroutine
   
   !----------------------------------------------------------------------------
   subroutine finalize_mlo_atoav_new
   !
   !   Author: Rafael Lago (MPCDF) January 2018
   !   
      integer :: i, ierr
      
      deallocate(lmr2buf)
      deallocate(buf2mlo)
      mlo_atoav_new_initialized = .false.
   end subroutine finalize_mlo_atoav_new
   

   !----------------------------------------------------------------------------
   subroutine create_mlo_atoav_new(this, n_fields)
   !   
   !   Author: Rafael Lago (MPCDF) April 2020

      class(type_mpiatoav_new) :: this
      integer, intent(in) :: n_fields
      integer :: i, buffer_size
      
      if (this%initialized) return 
      
      if (.not. mlo_atoav_new_initialized) call initialize_mlo_atoav_new
      this%n_fields = n_fields
      allocate(this%send_count(0:n_ranks_lo-1))
      allocate(this%recv_count(0:n_ranks_lo-1))
      allocate(this%send_disp(0:n_ranks_lo-1))
      allocate(this%recv_disp(0:n_ranks_lo-1))
      
      this%send_count = send_lm_count*n_r_loc*n_fields
      this%recv_count = recv_mlo_count*n_mlo_loc*n_fields
      
      this%send_disp = 0
      this%recv_disp = 0
      do i=1,n_ranks_lo-1
         this%send_disp(i) = this%send_disp(i-1)+this%send_count(i-1)
         this%recv_disp(i) = this%recv_disp(i-1)+this%recv_count(i-1)
      end do
      this%send_sum = sum(this%send_count)
      this%recv_sum = sum(this%recv_count)
      buffer_size = max(this%send_sum, this%recv_sum)
      
      allocate(this%send_buf(n_r_loc, n_fields, n_lm_loc))
      allocate(this%recv_buf(n_r_max*n_fields*n_mlo_loc))
      this%send_buf = cmplx(-1.0, -1.0)
      this%recv_buf = cmplx(-1.0, -1.0)
      
      n_atoav_new_obj = n_atoav_new_obj + 1
      
   end subroutine create_mlo_atoav_new
   
   !----------------------------------------------------------------------------
   subroutine destroy_mlo_atoav_new(this)
   !   
   !   Author: Rafael Lago (MPCDF) April 2020
   !   
      class(type_mpiatoav_new) :: this
      integer :: i, ierr
      
      if (.not. this%initialized) return
      
      nullify(this%send_buf)
      nullify(this%recv_buf)
!       nullify(this%buffer)
      
      deallocate(this%send_count)
      deallocate(this%recv_count)
      deallocate(this%send_disp)
      deallocate(this%recv_disp)
      
      n_atoav_new_obj = n_atoav_new_obj - 1
      this%initialized = .false.
      
      if (n_atoav_new_obj==0) call finalize_mlo_atoav_new
      
   end subroutine destroy_mlo_atoav_new
   
   !----------------------------------------------------------------------------
   subroutine transp_lm2r_atoav_new_dist(this, arr_LMloc, arr_Rloc)
   !   
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoav_new)     :: this
      complex(cp), intent(in)  :: arr_LMloc(n_mlo_loc, n_r_max, *)
      complex(cp), intent(out) :: arr_Rloc(n_lm_loc, nRstart:nRstop, *)
      
      !-- Local variables:
      integer :: ierr

      call reorder_lmloc2buffer(this, arr_LMloc)
      call mpi_alltoallv(this%recv_buf, this%recv_count, this%recv_disp, MPI_DEF_COMPLEX, &
                         this%send_buf, this%send_count, this%send_disp, MPI_DEF_COMPLEX, &
                         comm_r, ierr)
      call this%reorder_buffer2rloc(arr_Rloc)
      
   end subroutine transp_lm2r_atoav_new_dist
   
   !----------------------------------------------------------------------------
   subroutine transp_r2lm_atoav_new_dist(this, arr_Rloc, arr_LMloc)
   !   
   !   Author: Rafael Lago (MPCDF) May 2020
   !   
      class(type_mpiatoav_new)     :: this
      complex(cp), intent(in)  :: arr_Rloc(n_lm_loc, nRstart:nRstop, *)
      complex(cp), intent(out) :: arr_LMloc(n_mlo_loc, n_r_max, *)
      
      !-- Local variables:
      integer :: ierr
      
      call this%reorder_rloc2buffer(arr_Rloc)
      call mpi_alltoallv(this%send_buf, this%send_count, this%send_disp, MPI_DEF_COMPLEX, &
                         this%recv_buf, this%recv_count, this%recv_disp, MPI_DEF_COMPLEX, &
                         comm_r, ierr)
      call this%reorder_buffer2lmloc(arr_LMloc)

   end subroutine transp_r2lm_atoav_new_dist
   
   !----------------------------------------------------------------------------
   subroutine transp_lm2r_atoav_new_dummy(this, arr_LMloc, arr_Rloc)
      class(type_mpiatoav_new) :: this
      complex(cp), intent(in) :: arr_LMloc(llm:ulm,1:n_r_max,*)
      complex(cp), intent(out) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)
      print*, "Dummy transp_lm2r_atoav_new_dummy, not yet implemented!"
   end subroutine transp_lm2r_atoav_new_dummy

   subroutine transp_r2lm_atoav_new_dummy(this, arr_Rloc, arr_LMloc)
      class(type_mpiatoav_new) :: this
      complex(cp), intent(in) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)
      complex(cp), intent(out) :: arr_LMloc(llm:ulm,1:n_r_max,*)
      print*, "Dummy transp_r2lm_atoav_new_dummy, not yet implemented!"
   end subroutine transp_r2lm_atoav_new_dummy
   
end module mod_mpiatoav_new
