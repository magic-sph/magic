module time_array

   use precision_mod
   use constants, only: zero
   use mem_alloc, only: bytes_allocated

   implicit none

   private

   type, public :: type_tarray

      complex(cp), allocatable :: impl(:,:,:)
      complex(cp), pointer :: expl(:,:,:)
      complex(cp), allocatable :: old(:,:,:)

   contains

      procedure :: initialize
      procedure :: set_initial_values
      procedure :: finalize

   end type type_tarray

contains

   subroutine initialize(this, nMstart, nMstop, n_r_max, norder_imp, norder_exp, &
              &          norder_imp_lin)

      class(type_tarray) :: this

      !-- Input variables
      integer, intent(in) :: nMstart
      integer, intent(in) :: nMstop
      integer, intent(in) :: n_r_max
      integer, intent(in) :: norder_imp
      integer, intent(in) :: norder_exp
      integer, intent(in) :: norder_imp_lin

      allocate( this%impl(nMstart:nMstop,n_r_max,norder_imp_lin-1) )
      !allocate( this%expl(nMstart:nMstop,n_r_max,norder_exp) )
      allocate( this%old(nMstart:nMstop,n_r_max,norder_imp-1) )

      bytes_allocated = bytes_allocated + (nMstop-nMstart+1)*n_r_max*(norder_exp+&
      &                 norder_imp+norder_imp_lin-2)

      call this%set_initial_values()

   end subroutine initialize
!----------------------------------------------------------------------------------
   subroutine set_initial_values(this)

      class(type_tarray) :: this

      this%old(:,:,:) =zero
      !this%expl(:,:,:)=zero
      this%impl(:,:,:)=zero

   end subroutine set_initial_values
!----------------------------------------------------------------------------------
   subroutine finalize(this)

      class(type_tarray) :: this

      deallocate( this%old, this%impl )

   end subroutine finalize
!----------------------------------------------------------------------------------
end module time_array
