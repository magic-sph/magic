module TO_arrays_mod

   use truncation, only: l_max
   use precision_mod, only: cp

   implicit none

   private

   type, public :: TO_arrays_t
      !----- Local TO output stuff:
      real(cp), allocatable :: dzRstrLM(:),dzAstrLM(:)
      real(cp), allocatable :: dzCorLM(:),dzLFLM(:)

   contains
      procedure :: initialize
      procedure :: finalize
   end type TO_arrays_t

contains

   subroutine initialize(this)

      class(TO_arrays_t) :: this

      allocate( this%dzRstrLM(l_max+2),this%dzAstrLM(l_max+2) )
      allocate( this%dzCorLM(l_max+2),this%dzLFLM(l_max+2) )

   end subroutine initialize
!----------------------------------------------------------------------------
   subroutine finalize(this)

      class(TO_arrays_t) :: this

      deallocate( this%dzRstrLM,this%dzAstrLM )
      deallocate( this%dzCorLM,this%dzLFLM )

   end subroutine finalize
!----------------------------------------------------------------------------
end module TO_arrays_mod
