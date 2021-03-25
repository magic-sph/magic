#include "perflib_preproc.cpp"
module mod_mpiatoaw_2d

   use precision_mod
   use parallel_mod
   use mem_alloc
   use truncation
   use blocking, only: lm_balance, lo_map, st_map, llm, ulm
   use mpi_transp, only: type_mpitransp
   use fft
   use LMmapping
!    use mod_mpiatoav_2d, only: initialize_mlo_atoav_2d, lmr2buf, buf2mlo, send_lm_count, recv_mlo_count


   type, extends(type_mpitransp) :: type_mpiatoaw_2d
      logical :: initialized = .false.

   contains
      procedure :: create_comm => create_mlo_atoaw_2d
      procedure :: destroy_comm => destroy_mlo_atoaw_2d
!       procedure :: transp_lm2r => transp_lm2r_atoaw_2d_dummy
!       procedure :: transp_r2lm => transp_r2lm_atoaw_2d_dummy
      procedure :: transp_lm2r => transp_lm2r_atoaw_2d_dist
      procedure :: transp_r2lm => transp_r2lm_atoaw_2d_dist
      procedure :: transp_lm2r_dist => transp_lm2r_atoaw_2d_dist
      procedure :: transp_r2lm_dist => transp_r2lm_atoaw_2d_dist
   end type type_mpiatoaw_2d
   
   !
   ! These are the same for all objects of the type type_mpiatoaw_2d
   integer, private, allocatable :: lmloc_type(:,:),  rloc_type(:,:)
   logical, private :: mlo_atoaw_2d_initialized = .false.
   integer, private :: n_atoaw_2d_obj = 0
   
   public :: type_mpiatoaw_2d
   
contains
   
   subroutine initialize_mlo_atoaw_2d
      !-- Initialize the MPI types for the transposition from ML to Radial
      !   
      !   This is taylored for the lfirst distribution of lm-points.
      !   It assumes that all ranks in comm_r have the same n_mlo_loc
      !   
      !   Author: Rafael Lago (MPCDF) March 2021
      !    
      !-- TODO
      !
      !-- Local variables
      integer :: icoord_lo
      integer :: lr, nr, m, l, k
      
      integer (kind=mpi_address_kind) :: lb, extend, bytesCMPLX
      integer :: blocklenghts(max(n_mlo_array, n_lm_loc))
      integer :: indexes(n_mlo_array), col_type, ext_type, mtx_type
      
      if (mlo_atoaw_2d_initialized) return

      call MPI_Type_Get_Extent(MPI_DEF_COMPLEX, lb, bytesCMPLX, ierr) 

      allocate(lmloc_type(0:n_ranks_r-1,3))
      allocate( rloc_type(0:n_ranks_r-1,3))
      
      blocklenghts = 1                  ! Fixed...
      lmloc_type(:,1) = 0               ! Type
      lmloc_type(:,2) = 1               ! Count
      lmloc_type(:,3) = 0               ! Displacement
      rloc_type(:,1) = 0                ! Type
      rloc_type(:,2) = 1                ! Count
      rloc_type(:,3) = 0                ! Displacement
      
      
      do icoord_lo=0,n_ranks_r-1
         nr = dist_r(icoord_lo,0)
         lr = dist_r(icoord_lo,1)
         call mpi_type_create_subarray(2, (/n_mlo_loc,n_r_max/), (/n_mlo_loc,nr/), &
               (/0,lr-1/), MPI_ORDER_FORTRAN, MPI_DEF_COMPLEX, lmloc_type(icoord_lo,1), ierr)
         call MPI_Type_commit(lmloc_type(icoord_lo,1),ierr)

         k = 0
         indexes = -1
         do i = 1, n_mlo_array
            m = dist_mlo(coord_mo,icoord_lo,i,1)
            l = dist_mlo(coord_mo,icoord_lo,i,2)
            if ((m>=0).and.(l>=0)) then
               k = k + 1
               indexes(k) = map_dist_st%lm2(l,m)-1
            end if
         end do
         
         if (k>0) then
            call MPI_Type_indexed(k,blocklenghts(1:k), indexes(1:k), MPI_DEF_COMPLEX, col_type, ierr)
            extend = int(n_lm_loc*bytesCMPLX,kind=mpi_address_kind)
            call MPI_Type_create_resized(col_type, lb, extend, ext_type, ierr)
            call MPI_Type_contiguous(n_r_loc, ext_type, rloc_type(icoord_lo, 1), ierr)
            
            call MPI_Type_commit(rloc_type(icoord_lo, 1),ierr)
            call MPI_Type_free(col_type,ierr)
            call MPI_Type_free(ext_type,ierr)
         else
            print *, " - Strange error, I didn't take that into account! - Lago"
            stop
         end if
      end do
      
      mlo_atoaw_2d_initialized = .true.

   end subroutine initialize_mlo_atoaw_2d
!----------------------------------------------------------------------------
   subroutine finalize_mlo_atoaw_2d
      !   
      !   Author: Rafael Lago (MPCDF) March 2021
      !   
      integer :: i, ierr
      
      do i=0,n_ranks_r-1
         call MPI_Type_free(lmloc_type(i,1),ierr)
         call MPI_Type_free(rloc_type(i,1),ierr)
      end do
      deallocate(lmloc_type)
      deallocate(rloc_type)
      
      mlo_atoaw_2d_initialized = .false.
      
   end subroutine finalize_mlo_atoaw_2d
!----------------------------------------------------------------------------
   subroutine create_mlo_atoaw_2d(this, n_fields)
      !   
      !   Author: Rafael Lago (MPCDF) March 2021
      !
      class(type_mpiatoaw_2d) :: this
      integer, intent(in) :: n_fields
      
      if (this%initialized) return
      
      if (.not. mlo_atoaw_2d_initialized) call initialize_mlo_atoaw_2d
      this%n_fields = n_fields
      n_atoaw_2d_obj = n_atoaw_2d_obj + 1
      
   end subroutine create_mlo_atoaw_2d
!----------------------------------------------------------------------------
   subroutine destroy_mlo_atoaw_2d(this)
      !   
      !   Author: Rafael Lago (MPCDF) March 2021
      !
      class(type_mpiatoaw_2d) :: this
      integer :: i, ierr
      
      if (.not. this%initialized) return
      
      n_atoaw_2d_obj = n_atoaw_2d_obj - 1
      this%initialized = .false.
      if (n_atoaw_2d_obj==0) call finalize_mlo_atoaw_2d
      
   end subroutine destroy_mlo_atoaw_2d
!----------------------------------------------------------------------------
   subroutine transp_lm2r_atoaw_2d_dist(this, arr_LMloc, arr_Rloc)
      !   
      !   Author: Rafael Lago (MPCDF) March 2021
      !
      class(type_mpiatoaw_2d)  :: this
      complex(cp), intent(in)  :: arr_LMloc(n_mlo_loc, n_r_max, *)
      complex(cp), intent(out) :: arr_Rloc(n_lm_loc, nRstart:nRstop, *)
      
      integer :: ierr
      
      PERFON('lm2rW')
      call mpi_alltoallw(arr_LMloc, lmloc_type(:,2)*this%n_fields, lmloc_type(:,3), lmloc_type(:,1), &
      &     arr_Rloc, rloc_type(:,2)*this%n_fields, rloc_type(:,3), rloc_type(:,1), comm_r, ierr)
      PERFOFF
      
   end subroutine transp_lm2r_atoaw_2d_dist
!----------------------------------------------------------------------------
   subroutine transp_r2lm_atoaw_2d_dist(this, arr_Rloc, arr_LMloc)
      !   
      !   Author: Rafael Lago (MPCDF) March 2021
      !   
      class(type_mpiatoaw_2d)  :: this
      complex(cp), intent(in)  :: arr_Rloc(n_lm_loc, nRstart:nRstop, *)
      complex(cp), intent(out) :: arr_LMloc(n_mlo_loc, n_r_max, *)
      
      integer :: ierr
      
      PERFON('r2lmW')
      call mpi_alltoallw(arr_Rloc, rloc_type(:,2)*this%n_fields, rloc_type(:,3), rloc_type(:,1), &
      &     arr_LMloc, lmloc_type(:,2)*this%n_fields, lmloc_type(:,3), lmloc_type(:,1), comm_r, ierr)
      PERFOFF

   end subroutine transp_r2lm_atoaw_2d_dist
!----------------------------------------------------------------------------
   subroutine transp_lm2r_atoaw_2d_dummy(this, arr_LMloc, arr_Rloc)
      class(type_mpiatoaw_2d) :: this
      complex(cp), intent(in) :: arr_LMloc(llm:ulm,1:n_r_max,*)
      complex(cp), intent(out) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)
      print*, "Dummy transp_lm2r_atoaw_2d_dummy, not yet implemented!"
   end subroutine transp_lm2r_atoaw_2d_dummy
!----------------------------------------------------------------------------
   subroutine transp_r2lm_atoaw_2d_dummy(this, arr_Rloc, arr_LMloc)
      class(type_mpiatoaw_2d) :: this
      complex(cp), intent(in) :: arr_Rloc(1:lm_max,nRstart:nRstop,*)
      complex(cp), intent(out) :: arr_LMloc(llm:ulm,1:n_r_max,*)
      print*, "Dummy transp_r2lm_atoaw_2d_dummy, not yet implemented!"
   end subroutine transp_r2lm_atoaw_2d_dummy
!----------------------------------------------------------------------------
end module mod_mpiatoaw_2d
