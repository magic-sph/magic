module time_array
   !
   ! This module defines two types that are defined to store the implicit/explicit
   ! terms at the different sub-stage/step.
   !

   use precision_mod
   use constants, only: zero
   use mem_alloc, only: bytes_allocated

   implicit none

   private

   type, public :: type_tarray
      integer :: llm ! Start index of first dimension
      integer :: ulm  ! Stop index of first dimension
      integer :: nRStart ! Start index of second dimension
      integer :: nRStop  ! Stop index of second dimension
      complex(cp), allocatable :: impl(:,:,:) ! Array that contains the implicit states
      complex(cp), pointer :: expl(:,:,:) ! Array that contains the explicit states
      complex(cp), allocatable :: old(:,:,:) ! Array that contains the old states
      logical :: l_exp
   contains
      procedure :: initialize
      procedure :: finalize
   end type type_tarray

   type, public :: type_tscalar
      real(cp), allocatable :: impl(:) ! Implicit states
      real(cp), allocatable :: old(:)  ! Old states
      real(cp), allocatable :: expl(:) ! Explicit states
   contains
      procedure :: initialize => initialize_scalar
      procedure :: finalize => finalize_scalar
   end type type_tscalar

contains

   subroutine initialize(this, llm, ulm, nRstart, nRstop, nold, nexp, nimp, &
              &          l_allocate_exp)
      !
      ! Memory allocation of the arrays and initial values set to zero
      !

      class(type_tarray) :: this

      !-- Input variables
      integer,           intent(in) :: llm ! Lower boundary of first dimension
      integer,           intent(in) :: ulm ! Upper boundary of first dimension
      integer,           intent(in) :: nRstart ! Lower boundary of second dimension
      integer,           intent(in) :: nRstop  ! Upper boundary of second dimension
      integer,           intent(in) :: nold ! Number of old states
      integer,           intent(in) :: nexp ! Number of explicit states
      integer,           intent(in) :: nimp ! Number of implicit states
      logical, optional, intent(in) :: l_allocate_exp ! A boolean to specify whether the explicit state has to be allocated

      !-- Local variable
      logical :: l_allocate

      if ( present(l_allocate_exp) ) then
         l_allocate = l_allocate_exp
      else
         l_allocate = .false.
      end if

      this%llm = llm
      this%ulm = ulm
      this%nRstart = nRstart
      this%nRstop = nRstop
      this%l_exp = l_allocate

      allocate( this%impl(llm:ulm,nRstart:nRstop,nimp) )
      allocate( this%old(llm:ulm,nRstart:nRstop,nold) )

      bytes_allocated = bytes_allocated + (ulm-llm+1)*(nRstop-nRstart+1)*(&
      &                 nold+nimp)*SIZEOF_DEF_COMPLEX

      this%old(:,:,:) =zero
      this%impl(:,:,:)=zero

      if ( l_allocate ) then
         allocate( this%expl(llm:ulm,nRstart:nRstop,nexp) )
         bytes_allocated = bytes_allocated + (ulm-llm+1)*(nRstop-nRstart+1)*nexp* &
         &                 SIZEOF_DEF_COMPLEX
         this%expl(:,:,:)=zero
      end if

   end subroutine initialize
!----------------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Memory deallocation
      !

      class(type_tarray) :: this

      deallocate( this%old, this%impl )
      if ( this%l_exp ) deallocate( this%expl )

   end subroutine finalize
!----------------------------------------------------------------------------------
   subroutine initialize_scalar(this, nold, nexp, nimp )
      !
      ! This routine initializes time stepper help arrays in case of a scalar
      !

      class(type_tscalar) :: this

      !-- Input variables
      integer, intent(in) :: nold ! Number of old states
      integer, intent(in) :: nexp ! Number of explicit states
      integer, intent(in) :: nimp ! Number of implicit states

      allocate( this%expl(nexp), this%old(nold), this%impl(nimp) )
      bytes_allocated = bytes_allocated + (nexp+nold+nimp)*SIZEOF_DEF_REAL

      this%expl(:)=0.0_cp
      this%old(:) =0.0_cp
      this%impl(:)=0.0_cp

   end subroutine initialize_scalar
!----------------------------------------------------------------------------------
   subroutine finalize_scalar(this)
      !
      ! Memory deallocation
      !

      class(type_tscalar) :: this

      deallocate( this%old, this%expl, this%impl )

   end subroutine finalize_scalar
!----------------------------------------------------------------------------------
end module time_array
