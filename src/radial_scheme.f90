module radial_scheme

   use precision_mod

   implicit none

   private

   type, abstract, public :: type_rscheme

      integer :: n_max
      integer :: order
      real(cp) :: rnorm
      real(cp) :: boundary_fac
      real(cp), allocatable :: rMat(:,:)
      real(cp), allocatable :: drMat(:,:)
      real(cp), allocatable :: d2rMat(:,:)
      real(cp), allocatable :: d3rMat(:,:)
      real(cp), allocatable :: dr(:,:)
      real(cp), allocatable :: ddr(:,:)
      real(cp), allocatable :: dddr(:,:)
      real(cp), allocatable :: dr_top(:,:)
      real(cp), allocatable :: dr_bot(:,:)
      real(cp), allocatable :: ddr_top(:,:)
      real(cp), allocatable :: ddr_bot(:,:)
      real(cp), allocatable :: dddr_top(:,:)
      real(cp), allocatable :: dddr_bot(:,:)
      real(cp), allocatable :: drx(:)   ! First derivative of non-linear mapping (see Bayliss and Turkel, 1990)
      real(cp), allocatable :: ddrx(:)  ! Second derivative of non-linear mapping
      real(cp), allocatable :: dddrx(:) ! Third derivative of non-linear mapping

      character(len=72) :: version

   contains

      procedure(initialize_if), deferred :: initialize
      procedure(empty_if), deferred :: finalize
      procedure(get_der_mat_if), deferred :: get_der_mat
      procedure(get_grid_if), deferred :: get_grid
      !procedure(get_dr_real_1d_if), deferred :: get_dr_complex
      procedure :: costf1_complex
      procedure :: costf1_real
      procedure :: costf1_complex_1d
      procedure :: costf1_real_1d
      generic :: costf1 => costf1_complex, costf1_real, costf1_complex_1d, costf1_real_1d

   end type type_rscheme

   interface 

      subroutine empty_if(this)
         import
         class(type_rscheme) :: this
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

      subroutine initialize_if(this,n_r_max,order)

         import
         class(type_rscheme) :: this
         integer, intent(in) :: n_r_max
         integer, intent(in) :: order

      end subroutine initialize_if

      subroutine get_der_mat_if(this,n_r_max)

         import
         class(type_rscheme) :: this
         integer, intent(in) :: n_r_max

      end subroutine get_der_mat_if

   end interface

contains

   subroutine costf1_complex(this,f,n_f_max,n_f_start,n_f_stop)

      class(type_rscheme) :: this

      !-- Input variables:
      integer,  intent(in) :: n_f_max            ! number of columns in f,f2
      integer,  intent(in) :: n_f_start,n_f_stop ! columns to be transformed

      !-- Output variables:
      complex(cp), intent(inout) :: f(n_f_max,*) ! data/coeff input

   end subroutine

   subroutine costf1_real(this,f,n_f_max,n_f_start,n_f_stop,f2)

      class(type_rscheme) :: this

      !-- Input variables:
      integer,  intent(in) :: n_f_max            ! number of columns in f,f2
      integer,  intent(in) :: n_f_start,n_f_stop ! columns to be transformed

      !-- Output variables:
      real(cp), intent(inout) :: f(n_f_max,*)   ! data/coeff input
      real(cp), intent(out) :: f2(n_f_max,*)    ! work array of the same size as f

   end subroutine

   subroutine costf1_real_1d(this,f)

      class(type_rscheme) :: this

      real(cp), intent(inout) :: f(*)   ! data/coeff input

   end subroutine

   subroutine costf1_complex_1d(this,f)

      class(type_rscheme) :: this

      complex(cp), intent(inout) :: f(*)   ! data/coeff input

   end subroutine

end module radial_scheme
