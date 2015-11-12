module dtB_arrays_mod

   use truncation, only: lmP_max_dtB
   use precision_mod, only: cp

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
end module dtB_arrays_mod
