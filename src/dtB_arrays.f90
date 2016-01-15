module dtB_arrays_mod

   use truncation, only: lmP_max_dtB
   use mem_alloc, only: bytes_allocated
   use precision_mod
   use constants, only: zero

   implicit none

   private

   type, public :: dtB_arrays_t
      !----- Local dtB output stuff:
      complex(cp), allocatable :: BtVrLM(:),BpVrLM(:),BrVtLM(:),BrVpLM(:), &
                                  &               BtVpLM(:), BpVtLM(:)
      complex(cp), allocatable :: BtVpCotLM(:),BpVtCotLM(:),BtVpSn2LM(:), &
                                  &               BpVtSn2LM(:)
      complex(cp), allocatable :: BrVZLM(:),BtVZLM(:),BtVZcotLM(:),       &
                                  &               BtVZsn2LM(:)
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: set_zero
   end type dtB_arrays_t


contains

   subroutine initialize(this)

      class(dtB_arrays_t) :: this

      allocate( this%BtVrLM(lmP_max_dtB) )
      allocate( this%BpVrLM(lmP_max_dtB) )
      allocate( this%BrVtLM(lmP_max_dtB) )
      allocate( this%BrVpLM(lmP_max_dtB) )
      allocate( this%BtVpLM(lmP_max_dtB) )
      allocate( this%BpVtLM(lmP_max_dtB) )
      allocate( this%BtVpCotLM(lmP_max_dtB) )
      allocate( this%BpVtCotLM(lmP_max_dtB) )
      allocate( this%BtVpSn2LM(lmP_max_dtB) )
      allocate( this%BpVtSn2LM(lmP_max_dtB) )
      allocate( this%BrVZLM(lmP_max_dtB) )
      allocate( this%BtVZLM(lmP_max_dtB) )
      allocate( this%BtVZcotLM(lmP_max_dtB) )
      allocate( this%BtVZsn2LM(lmP_max_dtB) )
      bytes_allocated = bytes_allocated+ 14*lmP_max_dtB*SIZEOF_DEF_COMPLEX

   end subroutine initialize
!----------------------------------------------------------------------------
   subroutine finalize(this)

      class(dtB_arrays_t) :: this

      deallocate( this%BtVrLM )
      deallocate( this%BpVrLM )
      deallocate( this%BrVtLM )
      deallocate( this%BrVpLM )
      deallocate( this%BtVpLM )
      deallocate( this%BpVtLM )
      deallocate( this%BtVpCotLM )
      deallocate( this%BpVtCotLM )
      deallocate( this%BtVpSn2LM )
      deallocate( this%BpVtSn2LM )
      deallocate( this%BrVZLM )
      deallocate( this%BtVZLM )
      deallocate( this%BtVZcotLM )
      deallocate( this%BtVZsn2LM )

   end subroutine finalize
!----------------------------------------------------------------------------
   subroutine set_zero(this)

      class(dtB_arrays_t) :: this

      integer :: l

      do l=1,lmP_max_dtB
         this%BtVrLM = zero
         this%BpVrLM = zero
         this%BrVtLM = zero
         this%BrVpLM = zero
         this%BtVpLM = zero
         this%BpVtLM = zero
         this%BrVZLM = zero
         this%BtVZLM = zero
         this%BtVpCotLM = zero
         this%BpVtCotLM = zero
         this%BtVZcotLM = zero
         this%BtVpSn2LM = zero
         this%BpVtSn2LM = zero
         this%BtVZsn2LM = zero
      end do

   end subroutine set_zero
!----------------------------------------------------------------------------
end module dtB_arrays_mod
