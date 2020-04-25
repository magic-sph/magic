module time_array
   !
   ! This module defines two types that are defined to store the implicit/explicit
   ! terms at the different sub-stage/step.
   !

   use precision_mod
   use constants, only: zero
   use mem_alloc, only: bytes_allocated
   use truncation     ! DELETEMEEE
   use mpi_thetap_mod ! DELETEMEEE

   implicit none

   private

   type, public :: type_tarray
      complex(cp), allocatable :: impl(:,:,:) ! Array that contains the implicit states ! DEPRECATED
      complex(cp), pointer :: expl(:,:,:) ! Array that contains the explicit states     ! DEPRECATED
      complex(cp), allocatable :: old(:,:,:) ! Array that contains the old states       ! DEPRECATED
      complex(cp), allocatable :: impl_dist(:,:,:) ! Array that contains the implicit states
      complex(cp), pointer :: expl_dist(:,:,:) ! Array that contains the explicit states
      complex(cp), allocatable :: old_dist(:,:,:) ! Array that contains the old states
      logical :: l_exp
      integer :: nold, nexp, nimp, n_r, n_mlo
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: slice_all
      procedure :: gather_all
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

   subroutine initialize(this, llm, ulm, n_r_max, nold, nexp, nimp, &
              &          l_allocate_exp)
      !
      ! Memory allocation of the arrays and initial values set to zero
      !@>TODO: if you're passing the dimensions by argument, you have to save them 
      !>  in the object for later use; it makes no sense to use scope-restricted
      !>  variabes here and later rely on global variables instead later
      !@>TODO add n_mlo_loc as input argument here

      class(type_tarray) :: this

      !-- Input variables
      integer,           intent(in) :: llm ! Lower boundary of first dimension DEPRECATED
      integer,           intent(in) :: ulm ! Upper boundary of first dimension DEPRECATED
      integer,           intent(in) :: n_r_max ! Second dimension
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
      
      this%nold = nold
      this%nexp = nexp
      this%nimp = nimp
      this%n_r = n_r_max
      this%n_mlo = n_mlo_loc !>@TODO: receive as input argument

      this%l_exp = l_allocate

      allocate( this%impl(llm:ulm,n_r_max,nimp) )
      allocate( this%old(llm:ulm,n_r_max,nold) )
      allocate( this%impl_dist(n_mlo_loc,n_r_max,nimp) )
      allocate( this%old_dist(n_mlo_loc,n_r_max,nold) )

      bytes_allocated = bytes_allocated + (ulm-llm+1)*n_r_max*(&
      &                 nold+nimp)*SIZEOF_DEF_COMPLEX

      this%old(:,:,:) =zero
      this%impl(:,:,:)=zero

      if ( l_allocate ) then
         allocate( this%expl(llm:ulm,n_r_max,nexp) )
         allocate( this%expl_dist(n_mlo_loc,n_r_max,nexp) )
         bytes_allocated = bytes_allocated + (ulm-llm+1)*n_r_max*nexp*SIZEOF_DEF_COMPLEX
         this%expl(:,:,:)=zero
      end if

   end subroutine initialize
!----------------------------------------------------------------------------------
   subroutine gather_all(this)
      class(type_tarray) :: this
      integer :: i, j
      
      do i=1,this%nimp
         call transform_new2old(this%impl_dist(:,:,i), this%impl(:,:,i))
      end do

      do i=1,this%nold
         call transform_new2old(this%old_dist(:,:,i), this%old(:,:,i))
      end do
         
      if (associated(this%expl)) then
         do i=1,this%nexp
            call transform_new2old(this%expl_dist(:,:,i), this%expl(:,:,i))
         end do
      end if
      
   end subroutine gather_all
!----------------------------------------------------------------------------------
   subroutine slice_all(this)
      class(type_tarray) :: this
      integer :: i
      
      do i=1,this%nimp
         call transform_old2new(this%impl(:,:,i), this%impl_dist(:,:,i))
      end do

      do i=1,this%nold
         call transform_old2new(this%old(:,:,i), this%old_dist(:,:,i))
      end do
         
      if (associated(this%expl)) then
         do i=1,this%nexp
            call transform_old2new(this%expl(:,:,i), this%expl_dist(:,:,i))
         end do
      end if
      
   end subroutine slice_all
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
