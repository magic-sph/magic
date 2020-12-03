module integration
   !
   ! Radial integration functions
   !

   use precision_mod
   use constants, only: half, one, two
   use radial_scheme, only: type_rscheme
   use cosine_transform_odd

   implicit none

   private

   public :: rIntIC, rInt_R, simps

contains

   real(cp) function rIntIC(f,nRmax,drFac,chebt)
      !
      !  This function performs the radial integral over a
      !  function f that is given on the appropriate nRmax
      !  radial Chebychev grid points.
      !

      !-- Input variables:
      integer,           intent(in) :: nRmax
      real(cp),          intent(inout) :: f(nRmax)
      real(cp),          intent(in) :: drFac
      type(costf_odd_t), intent(in) :: chebt

      !-- Local variables:
      integer :: nCheb,nChebInt
      real(cp) :: weight
      real(cp) :: work(nRmax)

      call chebt%costf1(f,work)
      f(1)    =half*f(1)
      f(nRmax)=half*f(nRmax)

      !-- Sum contribution:
      rIntIC=f(1)           ! This is zero order contribution
      do nCheb=2,nRmax      ! Only even chebs for IC
         nChebInt=2*nCheb-1
         weight  =-one/real(nChebInt*(nChebInt-2),cp)
         rIntIC  =rIntIC+weight*f(nCheb)
      end do

      !-- Remaining renormalisation:
      rIntIC=sqrt(two/real(nRmax-1,cp))*rIntIC/drFac

   end function rIntIC
!------------------------------------------------------------------------------
   real(cp) function rInt_R(f,r,r_scheme) result(rInt)
      !
      !   Same as function rInt but for a radial dependent mapping function
      !   dr_fac2.
      !

      !-- Input variables:
      real(cp),            intent(in) :: f(:)    ! Input function
      real(cp),            intent(in) :: r(:)    ! Radius
      class(type_rscheme), intent(in) :: r_scheme! Radial scheme (FD or Cheb)

      !-- Local variables
      real(cp), allocatable :: f2(:)
      integer :: nCheb, nRmax


      !--- Integrals:
      if ( r_scheme%version == 'cheb' ) then

         nRmax=size(f)
         allocate( f2(nRmax) )
         f2(:)=f(:)/r_scheme%drx(:)

         !-- Transform to cheb space:
         call r_scheme%costf1(f2)
         f2(1)    =half*f2(1)
         f2(nRmax)=half*f2(nRmax)

         !-- Sum contribution:
         rInt=f2(1)            ! This is zero order contribution
         !do nCheb=3,r_scheme%n_max,2 ! Only even chebs contribute
         do nCheb=3,nRmax,2 ! Only even chebs contribute
            rInt=rInt-one/real(nCheb*(nCheb-2),cp)*f2(nCheb)
         end do

         !-- Remaining renormalisation:
         rInt=two*sqrt(two/real(nRmax-1,cp))*rInt

         deallocate( f2 )

      else

         rInt=simps(f,r)

      end if

   end function rInt_R
!------------------------------------------------------------------------------
   real(cp) function simps(f,r) result(rInt)
      !
      ! Simpson's method to integrate a function
      !

      !-- Input variables:
      real(cp), intent(in) :: f(:)    ! Input function
      real(cp), intent(in) :: r(:)    ! Radius

      !-- Local variables
      real(cp) :: h1, h2
      integer :: n_r, n_r_max

      n_r_max=size(f)

      if ( mod(n_r_max,2)==1 ) then ! Odd number (Simpson ok)

         rInt = 0.0_cp
         do n_r=2,n_r_max-1,2
            h2=r(n_r+1)-r(n_r)
            h1=r(n_r)-r(n_r-1)
            rInt=rInt+(h1+h2)/6.0_cp*( f(n_r-1)*(two*h1-h2)/h1     +&
            &                      f(n_r)  *(h1+h2)*(h1+h2)/(h1*h2)+&
            &                      f(n_r+1)*(two*h2-h1)/h2 )
         end do

         rInt = -rInt

      else ! Even number (twice simpson + trapz on the first and last points)

         rInt = half*(r(2)-r(1))*(f(2)+f(1))
         do n_r=3,n_r_max,2
            h2=r(n_r+1)-r(n_r)
            h1=r(n_r)-r(n_r-1)
            rInt=rInt+(h1+h2)/6.0_cp*( f(n_r-1)*(two*h1-h2)/h1     +&
            &                      f(n_r)  *(h1+h2)*(h1+h2)/(h1*h2)+&
            &                      f(n_r+1)*(two*h2-h1)/h2 )
         end do
         rInt = rInt+half*(r(n_r_max)-r(n_r_max-1))*(f(n_r_max)+f(n_r_max-1))
         do n_r=2,n_r_max-1,2
            h2=r(n_r+1)-r(n_r)
            h1=r(n_r)-r(n_r-1)
            rInt=rInt+(h1+h2)/6.0_cp*( f(n_r-1)*(two*h1-h2)/h1     +&
            &                      f(n_r)  *(h1+h2)*(h1+h2)/(h1*h2)+&
            &                      f(n_r+1)*(two*h2-h1)/h2 )
         end do
         rInt = -half*rInt

      end if

   end function simps
!------------------------------------------------------------------------------
end module integration
