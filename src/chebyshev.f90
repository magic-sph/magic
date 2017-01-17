module chebyshev

   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: lm_max, lm_max_real
   use fft_fac_mod, only: fft_fac_complex, fft_fac_real
   use constants, only: half, one, two, three, four, pi, sin36, &
       &                cos36, sin60, sin72, cos72
   use LMLoop_data, only: llm,ulm
   use radial_scheme, only: type_rscheme
   use logic, only: l_newmap, l_PV, l_TO
   use useful, only: factorise
   use chebyshev_polynoms_mod, only: cheb_grid
   use cosine_transform_odd, only: costf_odd_t

   implicit none

   private

   type, public, extends(type_rscheme) :: type_cheb_odd

      real(cp) :: alpha1 !Input parameter for non-linear map to define degree of spacing (0.0:2.0)
      real(cp) :: alpha2 !Input parameter for non-linear map to define central point of different spacing (-1.0:1.0)
      type(costf_odd_t) :: chebt_oc
      real(cp), allocatable :: r_cheb(:)
      complex(cp), pointer :: work_costf(:,:)

   contains

      procedure :: initialize
      procedure :: finalize
      procedure :: get_der_mat
      procedure :: get_grid => initialize_mapping
      procedure :: costf1_complex_1d
      procedure :: costf1_complex
      procedure :: costf1_real
      procedure :: costf1_real_1d

   end type type_cheb_odd

contains

   subroutine initialize(this, n_r_max, order)
      !
      !  Purpose of this subroutine is to calculate and store several     
      !  values that will be needed for a fast cosine transform of the    
      !  first kind. The actual transform is performed by the             
      !  subroutine costf1.                                               
      !

      class(type_cheb_odd) :: this
      
      integer, intent(in) :: n_r_max
      integer, intent(in) :: order ! This is going to be n_cheb_max

      !-- Local variables:
      integer :: ni,nd

      this%rnorm = sqrt(two/real(n_r_max-1,kind=cp))
      this%n_max = order  ! n_cheb_max
      this%boundary_fac = half
      this%version = 'cheb'
      this%nRmax = n_r_max

      allocate( this%rMat(n_r_max,n_r_max) )
      allocate( this%drMat(n_r_max,n_r_max) )
      allocate( this%d2rMat(n_r_max,n_r_max) )
      allocate( this%d3rMat(n_r_max,n_r_max) )
      allocate( this%r_cheb(n_r_max) )
      bytes_allocated=bytes_allocated+(4*n_r_max*n_r_max+n_r_max)*SIZEOF_DEF_REAL

      allocate( this%work_costf(1:ulm-llm+1,n_r_max) )
      bytes_allocated=bytes_allocated+n_r_max*(ulm-llm+1)*SIZEOF_DEF_COMPLEX

      ni = 2*n_r_max+2
      nd = 2*n_r_max+5

      call this%chebt_oc%initialize(n_r_max, ni, nd)

      allocate ( this%drx(n_r_max), this%ddrx(n_r_max), this%dddrx(n_r_max) )
      bytes_allocated=bytes_allocated+3*n_r_max*SIZEOF_DEF_REAL

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine initialize_mapping(this, n_r_max, ricb, rcmb, ratio1, ratio2, r)

      class(type_cheb_odd) :: this

      !-- Input variables:
      integer,  intent(in) :: n_r_max
      real(cp), intent(in) :: ricb
      real(cp), intent(in) :: rcmb
      real(cp), intent(inout) :: ratio1
      real(cp), intent(in) :: ratio2

      !-- Output variable:
      real(cp), intent(out) :: r(n_r_max)

      !-- Local variables:
      integer :: n_r
      real(cp) :: lambd,paraK,paraX0 !parameters of the nonlinear mapping

      if ( l_newmap ) then
         this%alpha1=ratio1
         this%alpha2=ratio2
         paraK=atan(this%alpha1*(one+this%alpha2))/atan(this%alpha1*(one-this%alpha2))
         paraX0=(paraK-one)/(paraK+one)
         lambd=atan(this%alpha1*(one-this%alpha2))/(one-paraX0)
      else
         this%alpha1=0.0_cp
         this%alpha2=0.0_cp
      end if

      call cheb_grid(ricb,rcmb,n_r_max-1,r,this%r_cheb,this%alpha1,this%alpha2, &
           &         paraX0,lambd)

      if ( l_newmap ) then

         do n_r=1,n_r_max
            this%drx(n_r) =                          (two*this%alpha1) /   &
            &    ((one+this%alpha1**2*(two*r(n_r)-ricb-rcmb-this%alpha2)**2)* &
            &    lambd)
            this%ddrx(n_r) = -(8.0_cp*this%alpha1**3*(two*r(n_r)-ricb-rcmb-this%alpha2)) / &
            &    ((one+this%alpha1**2*(-two*r(n_r)+ricb+rcmb+this%alpha2)**2)**2*     &
            &    lambd)
            this%dddrx(n_r) =        (16.0_cp*this%alpha1**3*(-one+three*this%alpha1**2* &
            &                     (-two*r(n_r)+ricb+rcmb+this%alpha2)**2)) / &
            &    ((one+this%alpha1**2*(-two*r(n_r)+ricb+rcmb+this%alpha2)**2)**3* &
            &    lambd)
         end do

      else

         do n_r=1,n_r_max
            this%drx(n_r)  =two/(rcmb-ricb)
            this%ddrx(n_r) =0.0_cp
            this%dddrx(n_r)=0.0_cp
         end do

      end if

   end subroutine initialize_mapping
!------------------------------------------------------------------------------
   subroutine finalize(this)

      class(type_cheb_odd) :: this


      deallocate( this%rMat, this%drMat, this%d2rMat, this%d3rMat )
      deallocate( this%r_cheb )
      deallocate( this%drx )
      deallocate( this%ddrx, this%dddrx )
      deallocate( this%work_costf )

      call this%chebt_oc%finalize()

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine get_der_mat(this, n_r_max)
      !
      !  Construct Chebychev polynomials and their first, second,
      !  and third derivative up to degree n_r at n_r points x
      !  in the interval [a,b]. Since the Chebs are only defined
      !  in [-1,1] we have to use a map, mapping the points x
      !  points y in the interval [-1,1]. This map is executed
      !  by the subroutine cheb_grid and has to be done
      !  before calling this program.
      !

      class(type_cheb_odd) :: this

      !-- Input variables:
      integer, intent(in) :: n_r_max

      !-- Local variables:
      integer :: n,k   ! counter
      ! in [-1,1] to x-derivatives in [a,b]

      !-- definition of map_fac:
      !   d Cheb(y) / d x = d y / d x * d Cheb(y) / d y
      !                   = map_fac * d Cheb(y) / d y

      !-- construction of chebs and derivatives with recursion:
      do k=1,n_r_max  ! do loop over the n_r grid points !

         !----- set first two chebs:
         this%rMat(1,k)=one
         this%rMat(2,k)=this%r_cheb(k)
         this%drMat(1,k)=0.0_cp
         this%drMat(2,k)=this%drx(k)
         this%d2rMat(1,k)=0.0_cp
         this%d2rMat(2,k)=this%ddrx(k)
         this%d3rMat(1,k)=0.0_cp
         this%d3rMat(2,k)=this%dddrx(k)

         !----- now construct the rest with a recursion:
         do n=3,n_r_max ! do loop over the (n-1) order of the chebs

            this%rMat(n,k)=    two*this%r_cheb(k)*this%rMat(n-1,k)-this%rMat(n-2,k)
            this%drMat(n,k)=   two*this%drx(k)*this%rMat(n-1,k) + &
            &                        two*this%r_cheb(k)*this%drMat(n-1,k) - &
            &                                 this%drMat(n-2,k)
            this%d2rMat(n,k)=  two*this%ddrx(k)*this%rMat(n-1,k) + &
            &                four*this%drx(k)*this%drMat(n-1,k) + &
            &                       two*this%r_cheb(k)*this%d2rMat(n-1,k) - &
            &                                this%d2rMat(n-2,k)
            this%d3rMat(n,k)=  two*this%dddrx(k)*this%rMat(n-1,k) + &
            &              6.0_cp*this%ddrx(k)*this%drMat(n-1,k) + &
            &             6.0_cp*this%drx(k)*this%d2rMat(n-1,k) + &
            &                       two*this%r_cheb(k)*this%d3rMat(n-1,k) - &
            &                                this%d3rMat(n-2,k)

         end do

      end do

      !-- This transposition is needed to bring those matrices in alignement
      !-- with the fortran column-major storage (see update routines)
      this%rMat  =transpose(this%rMat)
      this%drMat =transpose(this%drMat)
      this%d2rMat=transpose(this%d2rMat)
      this%d3rMat=transpose(this%d3rMat)

   end subroutine get_der_mat
!------------------------------------------------------------------------------
   subroutine costf1_complex(this,f,n_f_max,n_f_start,n_f_stop,work_array)
      !
      !  Purpose of this subroutine is to perform a multiple
      !  cosine transforms for n+1 datapoints
      !  on the columns numbered n_f_start to n_f_stop in the array
      !  f(n_f_max,n+1)
      !  Depending whether the input f contains data or coeff arrays
      !  coeffs or data are returned in f.
      !
      class(type_cheb_odd) :: this

      !-- Input variables:
      integer,  intent(in) :: n_f_max            ! number of columns in f,f2
      integer,  intent(in) :: n_f_start,n_f_stop ! columns to be transformed
    
      !-- Output variables:
      complex(cp), intent(inout) :: f(n_f_max,this%nRmax) ! data/coeff input
      complex(cp), optional, target, intent(inout) :: work_array(n_f_max,this%nRmax)
    
      !-- Local variables:
      complex(cp), pointer :: work(:,:)


      if ( present(work_array) ) then
         work(1:,1:) => work_array(1:n_f_max,1:this%nRmax)
      else
         work(1:,1:) => this%work_costf(1:n_f_max,1:)
      end if

      call this%chebt_oc%costf1(f,n_f_max,n_f_start,n_f_stop,work(:,1:this%nRmax))
    
   end subroutine costf1_complex
!------------------------------------------------------------------------------
   subroutine costf1_complex_1d(this,f)
    
      class(type_cheb_odd) :: this

      !-- Output variables:
      complex(cp), intent(inout) :: f(this%nRmax)   ! data/coeff input
    
      !-- Local variables:
      complex(cp) :: work1d(this%nRmax)
    
      call this%chebt_oc%costf1(f, work1d)

   end subroutine costf1_complex_1d
!------------------------------------------------------------------------------
   subroutine costf1_real(this,f,n_f_max,n_f_start,n_f_stop,f2)

      class(type_cheb_odd) :: this
    
      !-- Input variables:
      integer,  intent(in) :: n_f_max            ! number of columns in f,f2
      integer,  intent(in) :: n_f_start,n_f_stop ! columns to be transformed
    
      !-- Output variables:
      real(cp), intent(inout) :: f(n_f_max,this%nRmax) ! data/coeff input
      real(cp), intent(out) :: f2(n_f_max,this%nRmax)  ! work array of the same size as f
      call this%chebt_oc%costf1(f,n_f_max,n_f_start,n_f_stop,f2)

   end subroutine costf1_real
!------------------------------------------------------------------------------
   subroutine costf1_real_1d(this,f)
    
      class(type_cheb_odd) :: this

      !-- Output variables:
      real(cp), intent(inout) :: f(this%nRmax)   ! data/coeff input
    
      !-- Local variables:
      real(cp) :: work1d_real(this%nrMax)

      call this%chebt_oc%costf1(f,work1d_real)

   end subroutine costf1_real_1d
!------------------------------------------------------------------------------
end module chebyshev
