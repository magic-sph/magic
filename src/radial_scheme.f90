module radial_scheme
   !
   ! This is an abstract type that defines the radial scheme used in MagIC
   !

   use precision_mod

   implicit none

   private

   type, abstract, public :: type_rscheme

      integer :: n_max
      integer :: order
      integer :: order_boundary
      integer :: nRmax
      real(cp) :: rnorm
      real(cp) :: boundary_fac
      real(cp), allocatable :: rMat(:,:)
      real(cp), allocatable :: drMat(:,:)
      real(cp), allocatable :: d2rMat(:,:)
      real(cp), allocatable :: d3rMat(:,:)
      real(cp), allocatable :: d4rMat(:,:)
      real(cp), allocatable :: dr(:,:)
      real(cp), allocatable :: ddr(:,:)
      real(cp), allocatable :: dddr(:,:)
      real(cp), allocatable :: ddddr(:,:)
      real(cp), allocatable :: dr_top(:,:)
      real(cp), allocatable :: dr_bot(:,:)
      real(cp), allocatable :: ddr_top(:,:)
      real(cp), allocatable :: ddr_bot(:,:)
      real(cp), allocatable :: drx(:) ! First derivative of non-linear mapping (see Bayliss and Turkel, 1990)

      character(len=72) :: version

   contains

      procedure(initialize_if),  deferred :: initialize
      procedure(empty_if),       deferred :: finalize
      procedure(get_der_mat_if), deferred :: get_der_mat
      procedure(get_grid_if),    deferred :: get_grid
      procedure(robin_bc_if),    deferred :: robin_bc
      procedure :: costf1_complex
      procedure :: costf1_complex_1d
      procedure :: costf1_real_1d
      generic :: costf1 => costf1_complex, costf1_complex_1d, costf1_real_1d

   end type type_rscheme

   interface

      subroutine empty_if(this, gpu_dct)
         import
         class(type_rscheme) :: this
         logical, optional, intent(in) :: gpu_dct   ! compute DCT-I on GPU with hipfort_hipfft
      end subroutine empty_if

      subroutine get_grid_if(this,n_r_max,ricb,rcmb,ratio1,ratio2,r)
         import
         class(type_rscheme) :: this
         !-- Input quantities:
         integer,  intent(in) :: n_r_max    ! Number of grid points
         real(cp), intent(inout) :: ratio1  ! Nboudary/Nbulk
         real(cp), intent(in) :: ratio2     ! drMin/drMax
         real(cp), intent(in) :: ricb       ! inner boundary
         real(cp), intent(in) :: rcmb       ! outer boundary
         !-- Output quantities:
         real(cp), intent(out) :: r(n_r_max) ! radius
      end subroutine get_grid_if

      subroutine initialize_if(this,n_r_max,order,order_boundary, gpu_dct)
         import
         class(type_rscheme) :: this
         integer, intent(in) :: n_r_max
         integer, intent(in) :: order
         integer, intent(in) :: order_boundary
         logical, optional, intent(in) :: gpu_dct   ! compute DCT-I on GPU with hipfort_hipfft
      end subroutine initialize_if

      subroutine get_der_mat_if(this,n_r_max)
         import
         class(type_rscheme) :: this
         integer, intent(in) :: n_r_max
      end subroutine get_der_mat_if

      subroutine robin_bc_if(this,atop,btop,rhs_top,abot,bbot,rhs_bot,f)
         import
         class(type_rscheme) :: this
         real(cp),    intent(in) :: atop, btop, abot, bbot
         complex(cp), intent(in) :: rhs_bot, rhs_top
         complex(cp), intent(inout) :: f(:)
      end subroutine robin_bc_if

   end interface

contains

   subroutine costf1_complex(this,f,n_f_max,n_f_start,n_f_stop,gpu_dct,work_array)

      class(type_rscheme) :: this

      !-- Input variables:
      integer,  intent(in) :: n_f_max            ! number of columns in f,f2
      integer,  intent(in) :: n_f_start,n_f_stop ! columns to be transformed
      logical, optional, intent(in) :: gpu_dct   ! compute DCT-I on GPU with hipfort_hipfft

      !-- Output variables:
      complex(cp), intent(inout) :: f(n_f_max,this%nRmax) ! data/coeff input
      complex(cp), optional, target, intent(inout) :: work_array(n_f_max,this%nRmax)

   end subroutine costf1_complex
!------------------------------------------------------------------------------
   subroutine costf1_real_1d(this,f,gpu_dct)

      class(type_rscheme) :: this

      real(cp), intent(inout) :: f(this%nRmax)   ! data/coeff input
      logical, optional, intent(in) :: gpu_dct   ! compute DCT-I on GPU with hipfort_hipfft

   end subroutine costf1_real_1d
!------------------------------------------------------------------------------
   subroutine costf1_complex_1d(this,f,gpu_dct)

      class(type_rscheme) :: this

      complex(cp), intent(inout) :: f(this%nRmax)   ! data/coeff input
      logical, optional, intent(in) :: gpu_dct   ! compute DCT-I on GPU with hipfort_hipfft

   end subroutine costf1_complex_1d
!------------------------------------------------------------------------------
end module radial_scheme
