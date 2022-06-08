module dtB_arrays_mod

   use truncation, only: lmP_max_dtB
   use mem_alloc, only: bytes_allocated
   use precision_mod
   use constants, only: zero

   implicit none

   private

   type, public :: dtB_arrays_t
      !----- Local dtB output stuff:
      complex(cp), allocatable :: BtVrLM(:), BpVrLM(:), BrVtLM(:)
      complex(cp), allocatable :: BtVpLM(:), BpVtLM(:), BrVpLM(:)
      complex(cp), allocatable :: BpVtBtVpCotLM(:), BpVtBtVpSn2LM(:)
      complex(cp), allocatable :: BrVZLM(:), BtVZLM(:), BtVZsn2LM(:)
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: set_zero
   end type dtB_arrays_t

contains

   subroutine initialize(this)

      class(dtB_arrays_t) :: this

      allocate( this%BtVrLM(lmP_max_dtB), this%BpVrLM(lmP_max_dtB) )
      allocate( this%BrVtLM(lmP_max_dtB), this%BrVpLM(lmP_max_dtB) )
      allocate( this%BtVpLM(lmP_max_dtB), this%BpVtLM(lmP_max_dtB) )
      allocate( this%BpVtBtVpCotLM(lmP_max_dtB), this%BpVtBtVpSn2LM(lmP_max_dtB) )
      allocate( this%BrVZLM(lmP_max_dtB), this%BtVZLM(lmP_max_dtB) )
      allocate( this%BtVZsn2LM(lmP_max_dtB) )
      bytes_allocated = bytes_allocated+ 11*lmP_max_dtB*SIZEOF_DEF_COMPLEX

   end subroutine initialize
!----------------------------------------------------------------------------
   subroutine finalize(this)

      class(dtB_arrays_t) :: this

      deallocate( this%BtVrLM, this%BpVrLM, this%BrVtLM )
      deallocate( this%BrVpLM, this%BtVpLM, this%BpVtLM )
      deallocate( this%BpVtBtVpCotLM, this%BpVtBtVpSn2LM )
      deallocate( this%BrVZLM, this%BtVZLM, this%BtVZsn2LM )

   end subroutine finalize
!----------------------------------------------------------------------------
   subroutine set_zero(this)

      class(dtB_arrays_t) :: this

      this%BtVrLM(:) = zero
      this%BpVrLM(:) = zero
      this%BrVtLM(:) = zero
      this%BrVpLM(:) = zero
      this%BtVpLM(:) = zero
      this%BpVtLM(:) = zero
      this%BrVZLM(:) = zero
      this%BtVZLM(:) = zero
      this%BpVtBtVpCotLM(:) = zero
      this%BpVtBtVpSn2LM(:) = zero
      this%BtVZsn2LM(:) = zero

   end subroutine set_zero
!----------------------------------------------------------------------------
end module dtB_arrays_mod
