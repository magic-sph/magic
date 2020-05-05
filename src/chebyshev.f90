module chebyshev

   use precision_mod
   use mem_alloc, only: bytes_allocated
   use constants, only: half, one, two, three, four, pi
   use radial_scheme, only: type_rscheme
   use useful, only: factorise
   use chebyshev_polynoms_mod, only: cheb_grid
   use cosine_transform_odd, only: costf_odd_t
   use num_param, only: map_function
   use truncation, only: n_mlo_loc

   implicit none

   private

   type, public, extends(type_rscheme) :: type_cheb_odd
      real(cp) :: alpha1 !Input parameter for non-linear map to define degree of spacing (0.0:2.0)
      real(cp) :: alpha2 !Input parameter for non-linear map to define central point of different spacing (-1.0:1.0)
      logical :: l_map
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
      procedure :: robin_bc
   end type type_cheb_odd

contains

   subroutine initialize(this, n_r_max, order, order_boundary)
      !
      !  Purpose of this subroutine is to calculate and store several     
      !  values that will be needed for a fast cosine transform of the    
      !  first kind. The actual transform is performed by the             
      !  subroutine costf1.                                               
      !

      class(type_cheb_odd) :: this
      
      integer, intent(in) :: n_r_max
      integer, intent(in) :: order ! This is going to be n_cheb_max
      integer, intent(in) :: order_boundary ! this is used to determine whether mappings are used
            
      !-- Local variables:
      integer :: ni,nd

      this%rnorm = sqrt(two/real(n_r_max-1,kind=cp))
      this%n_max = order  ! n_cheb_max
      this%boundary_fac = half
      this%version = 'cheb'
      this%nRmax = n_r_max
      this%order_boundary=order_boundary

      if ( order_boundary == 1 ) then
         this%l_map=.true.
      else
         this%l_map=.false.
      end if

      allocate( this%rMat(n_r_max,n_r_max) )
      allocate( this%drMat(n_r_max,n_r_max) )
      allocate( this%d2rMat(n_r_max,n_r_max) )
      allocate( this%d3rMat(n_r_max,n_r_max) )
      allocate( this%r_cheb(n_r_max) )
      bytes_allocated=bytes_allocated+(4*n_r_max*n_r_max+n_r_max)*SIZEOF_DEF_REAL

      allocate( this%work_costf(n_mlo_loc,n_r_max) )
      bytes_allocated=bytes_allocated+n_r_max*n_mlo_loc*SIZEOF_DEF_COMPLEX

      allocate( this%dr_top(n_r_max,1), this%dr_bot(n_r_max,1) )
      bytes_allocated=bytes_allocated+2*n_r_max*SIZEOF_DEF_REAL
      this%dr_top(:,:)=0.0_cp
      this%dr_bot(:,:)=0.0_cp

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
      real(cp) :: lambd,paraK,paraX0 !parameters of the nonlinear mapping

      !--
      !-- There's possibly an issue when the Chebyshev mapping was used in
      !-- the old grid and a different mapping is used on the new one
      !--

      if ( this%l_map ) then
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
           &         paraX0,lambd,this%l_map)

      if ( this%l_map ) then

         !-- Tangent mapping (see Bayliss et al. 1992)
         if ( index(map_function, 'TAN') /= 0 .or.      &
         &    index(map_function, 'BAY') /= 0 ) then
            this%drx(:)  =two*this%alpha1*(rcmb - ricb)/(lambd*(this%alpha1**2* &
            &             (this%alpha2*(rcmb-ricb)-two*r(:)+rcmb+ricb)**2+      &
            &             (rcmb-ricb)**2))
            this%ddrx(:) =8.0_cp*this%alpha1**3*(rcmb-ricb)*(this%alpha2*       &
            &             (rcmb-ricb)-two*r(:)+rcmb+ricb)/(lambd*(              &
            &             this%alpha1**2*(this%alpha2*(rcmb-ricb)-two*r(:)      &
            &             +rcmb+ricb)**2+(rcmb-ricb)**2)**2)
            this%dddrx(:)=16.0_cp*this%alpha1**3*(rcmb-ricb)*(three*            &
            &             this%alpha1**2*(this%alpha2*(rcmb-ricb)-two*r(:)      &
            &             +rcmb+ricb)**2-(rcmb-ricb)**2)/(lambd*(this%alpha1**2*&
            &             (this%alpha2*(rcmb-ricb)-two*r(:)+rcmb+ricb)**2+      &
            &             (rcmb-ricb)**2)**3)

         !-- Arcsin mapping (see Kosloff and Tal-Ezer, 1993)
         else if ( index(map_function, 'ARCSIN') /= 0 .or. &
         &         index(map_function, 'KTL') /= 0 ) then
            this%drx(:)  =two*asin(this%alpha1)/this%alpha1*sqrt(one-           &
            &             this%alpha1**2*this%r_cheb(:)**2)/(rcmb-ricb)
            this%ddrx(:) =-four*asin(this%alpha1)**2*this%r_cheb(:)/            &
            &              (rcmb-ricb)**2
            this%dddrx(:)=-8.0_cp*asin(this%alpha1)**3*sqrt(one-this%alpha1**2* &
            &             this%r_cheb(:)**2)/this%alpha1/(rcmb-ricb)**3
         end if

      else !-- Regular affine mapping between ricb and rcmb

         this%drx(:)  =two/(rcmb-ricb)
         this%ddrx(:) =0.0_cp
         this%dddrx(:)=0.0_cp

      end if

   end subroutine initialize_mapping
!------------------------------------------------------------------------------
   subroutine finalize(this)

      class(type_cheb_odd) :: this

      deallocate( this%rMat, this%drMat, this%d2rMat, this%d3rMat )
      deallocate( this%r_cheb, this%drx, this%ddrx, this%dddrx )
      deallocate( this%work_costf, this%dr_top, this%dr_bot )

      call this%chebt_oc%finalize()

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine robin_bc(this,atop,btop,rhs_top,abot,bbot,rhs_bot,f)
      !
      ! This subroutine is used to determine the two boundary points of a field
      ! f subject to two Robin boundary conditions of the form:
      !
      ! atop*df/dr+btop*f = rhs_top;  abot*df/dr+bbot*f = rhs_bot
      !
      ! The method follows Canuto, SIAM, 1986 (p. 818)

      class(type_cheb_odd) :: this

      !-- Input variables
      real(cp), intent(in) :: atop, btop, abot, bbot
      complex(cp), intent(in) :: rhs_top, rhs_bot

      !-- In/out variables: only the boundary points are changed
      complex(cp), intent(inout) :: f(:)

      !-- Local variables
      integer :: n_r, nRmax
      real(cp) :: sum_top_r, sum_top_i, sum_bot_r, sum_bot_i
      real(cp) :: val_top_r, val_top_i, val_bot_r, val_bot_i

      nRmax = size(f)

      !-- First construct the sums that will be needed afterwards
      sum_top_r = 0.0_cp
      sum_top_i = 0.0_cp
      sum_bot_r = 0.0_cp
      sum_bot_i = 0.0_cp
      do n_r=2,nRmax-1
         sum_top_r = sum_top_r+this%dr_top(n_r,1)*real(f(n_r))
         sum_top_i = sum_top_i+this%dr_top(n_r,1)*aimag(f(n_r))
         sum_bot_r = sum_bot_r+this%dr_bot(n_r,1)*real(f(n_r))
         sum_bot_i = sum_bot_i+this%dr_bot(n_r,1)*aimag(f(n_r))
      end do

      !-- First get the values at the bottom boundary
      val_bot_r=( abot*this%dr_bot(1,1)*(real(rhs_top)-atop*sum_top_r)-          &
      &          (atop*this%dr_top(1,1)+btop)*(real(rhs_bot)-abot*sum_bot_r) ) / &
      &         (abot*atop*this%dr_bot(1,1)*this%dr_top(nRmax,1)-                &
      &         (atop*this%dr_top(1,1)+btop)*(abot*this%dr_bot(nRmax,1)+bbot))
      val_bot_i=( abot*this%dr_bot(1,1)*(aimag(rhs_top)-atop*sum_top_i)-         &
      &          (atop*this%dr_top(1,1)+btop)*(aimag(rhs_bot)-abot*sum_bot_i) ) /&
      &         (abot*atop*this%dr_bot(1,1)*this%dr_top(nRmax,1)-                &
      &         (atop*this%dr_top(1,1)+btop)*(abot*this%dr_bot(nRmax,1)+bbot))

      !-- Then get the values at the top boundary
      val_top_r=(real(rhs_top)-atop*(sum_top_r+this%dr_top(nRmax,1)*val_bot_r))/ &
      &         (atop*this%dr_top(1,1)+btop)
      val_top_i=(aimag(rhs_top)-atop*(sum_top_i+this%dr_top(nRmax,1)*val_bot_i))/&
      &         (atop*this%dr_top(1,1)+btop)

      !-- Finally assemble the complex numbers
      f(1)     = cmplx(val_top_r, val_top_i, kind=cp)
      f(nRmax) = cmplx(val_bot_r, val_bot_i, kind=cp)

   end subroutine robin_bc
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
      real(cp) :: coeff, diff
      integer :: n,k   ! counter
      ! in [-1,1] to x-derivatives in [a,b]

      !-- definition of map_fac:
      !   d Cheb(y) / d x = d y / d x * d Cheb(y) / d y
      !                   = map_fac * d Cheb(y) / d y

      !-- construction of chebs and derivatives with recursion:
      do k=1,n_r_max  ! do loop over the n_r grid points !

         !----- set first two chebs:
         this%rMat(1,k)  =one
         this%rMat(2,k)  =this%r_cheb(k)
         this%drMat(1,k) =0.0_cp
         this%drMat(2,k) =this%drx(k)
         this%d2rMat(1,k)=0.0_cp
         this%d2rMat(2,k)=this%ddrx(k)
         this%d3rMat(1,k)=0.0_cp
         this%d3rMat(2,k)=this%dddrx(k)

         !----- now construct the rest with a recursion:
         do n=3,n_r_max ! do loop over the (n-1) order of the chebs

            this%rMat(n,k) =two*this%r_cheb(k)*this%rMat(n-1,k)-this%rMat(n-2,k)
            this%drMat(n,k)=    two*this%drx(k)*this%rMat(n-1,k) + &
            &               two*this%r_cheb(k)*this%drMat(n-1,k) - &
            &                                  this%drMat(n-2,k)
            this%d2rMat(n,k)=  two*this%ddrx(k)*this%rMat(n-1,k) + &
            &                 four*this%drx(k)*this%drMat(n-1,k) + &
            &              two*this%r_cheb(k)*this%d2rMat(n-1,k) - &
            &                                 this%d2rMat(n-2,k)
            this%d3rMat(n,k)=  two*this%dddrx(k)*this%rMat(n-1,k) + &
            &               6.0_cp*this%ddrx(k)*this%drMat(n-1,k) + &
            &               6.0_cp*this%drx(k)*this%d2rMat(n-1,k) + &
            &               two*this%r_cheb(k)*this%d3rMat(n-1,k) - &
            &                                  this%d3rMat(n-2,k)

         end do

      end do

      !-- This transposition is needed to bring those matrices in alignement
      !-- with the fortran column-major storage (see update routines)
      this%rMat  =transpose(this%rMat)
      this%drMat =transpose(this%drMat)
      this%d2rMat=transpose(this%d2rMat)
      this%d3rMat=transpose(this%d3rMat)

      !-- Compute a vector that allows the computation of the first derivative
      !-- on the boundary point
      this%dr_top(1,1)=(two*(n_r_max-1)*(n_r_max-1)+one)/6.0_cp
      do k=2,n_r_max
         diff = two*sin( half*(k-1)*pi/(n_r_max-1) ) * sin( half*(k-1)*pi/(n_r_max-1) ) 
         if ( mod(k,2) == 0 ) then
            coeff=-one
         else
            coeff=one
         end if
         this%dr_top(k,1)=two * coeff/diff
      end do
      !-- Factor half for the last one
      this%dr_top(n_r_max,1) = half*this%dr_top(n_r_max,1) 

      !-- dr bot is the reverted vector with a -1 multiplication
      !-- The "flipping-trick" is used to minimize the round-off error
      !-- See Baltensperger & Trummer, SIAM, 2003, p. 1478
      do k=1,n_r_max
         this%dr_bot(k,1)=-this%dr_top(n_r_max+1-k,1)
      end do

      !-- Finally multiply by the mapping to get the proper derivative of r
      this%dr_top(:,1) = this%dr_top(:,1)*this%drx(1)
      this%dr_bot(:,1) = this%dr_bot(:,1)*this%drx(n_r_max)

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
