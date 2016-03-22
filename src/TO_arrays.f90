module TO_arrays_mod

   use truncation, only: l_max
   use mem_alloc, only: bytes_allocated
   use precision_mod

   implicit none

   private

   type, public :: TO_arrays_t
      !----- Local TO output stuff:
      real(cp), allocatable :: dzRstrLM(:),dzAstrLM(:)
      real(cp), allocatable :: dzCorLM(:),dzLFLM(:)

   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: set_zero
   end type TO_arrays_t

contains

   subroutine initialize(this)

      class(TO_arrays_t) :: this

      allocate( this%dzRstrLM(l_max+2),this%dzAstrLM(l_max+2) )
      allocate( this%dzCorLM(l_max+2),this%dzLFLM(l_max+2) )
      bytes_allocated = bytes_allocated+4*(l_max+2)*SIZEOF_DEF_REAL

   end subroutine initialize
!----------------------------------------------------------------------------
   subroutine finalize(this)

      class(TO_arrays_t) :: this

      deallocate( this%dzRstrLM,this%dzAstrLM )
      deallocate( this%dzCorLM,this%dzLFLM )

   end subroutine finalize
!----------------------------------------------------------------------------
   subroutine set_zero(this)

      class(TO_arrays_t) :: this

      !-- Local variable:
      integer :: l

      do l=0,l_max
         this%dzRstrLM(l+1)=0.0_cp
         this%dzAstrLM(l+1)=0.0_cp
         this%dzCorLM(l+1) =0.0_cp
         this%dzLFLM(l+1)  =0.0_cp
      end do

   end subroutine set_zero
!----------------------------------------------------------------------------
end module TO_arrays_mod
