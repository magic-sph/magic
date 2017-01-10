module integration
   !
   ! Radial integration functions
   !

   use precision_mod
   use constants, only: half, one, two
   use cosine_transform_odd

   implicit none

   private

   public :: rInt, rIntIC, rInt_R, rInt_fd

contains

   real(cp) function rInt(f,nRmax,drFac,chebt)
      !
      !  This function performs the radial integral over a
      !  function f that is given on the appropriate nRmax
      !  radial Chebychev grid points.
      !
      !  .. note:: drFac maps radius to cheb space [-1,1]
      !            drFac=two/(rMax-rMin)
      !

      !-- Input variables:
      integer,           intent(in) :: nRmax
      real(cp),          intent(in) :: f(nRmax)
      real(cp),          intent(in) :: drFac
      type(costf_odd_t), intent(in) :: chebt

      !-- Local variables:
      integer :: nCheb,nR
      real(cp) :: work(nRmax)
      real(cp) :: f2(nRmax)

      !-- Save input:
      do nR=1,nRmax
         f2(nR)=f(nR)
      end do

      !-- Transform to cheb space:
      call chebt%costf1(f2,work)
      f2(1)    =half*f2(1)
      f2(nRmax)=half*f2(nRmax)

      !-- Sum contribution:
      rInt=f2(1)           ! This is zero order contribution
      do nCheb=3,nRmax,2  ! Only even chebs contribute
         rInt=rInt-one/real(nCheb*(nCheb-2),cp)*f2(nCheb)
      end do

      !-- Remaining renormalisation:
      rInt=two/drFac*sqrt(two/real(nRmax-1,cp))*rInt

   end function rInt
!------------------------------------------------------------------------------
   real(cp) function rInt_fd(f,r,nRmax) result(rInt)

      !-- Input variables
      integer,  intent(in) :: nRmax
      real(cp), intent(in) :: f(nRmax)
      real(cp), intent(in) :: r(nRmax)

      !-- Local variables
      real(cp) :: h1, h2
      integer :: nR

      if ( mod(nRmax,2)==1 ) then ! Odd number (Simpson ok)

         rInt = 0.0_cp
         do nR=2,nRmax-1,2
            h2=r(nR+1)-r(nR)
            h1=r(nR)-r(nR-1)
            rInt=rInt+(h1+h2)/6.*( f(nR-1)*(two*h1-h2)/h1         +&
            &                      f(nR)  *(h1+h2)*(h1+h2)/(h1*h2)+&
            &                      f(nR+1)*(two*h2-h1)/h2 )
         end do

         print*, 'here'

      else ! Even number (twice simpson + trapz on the first and last points)

         rInt = half*(r(2)-r(1))*(f(2)+f(1))
         do nR=3,nRmax,2
            h2=r(nR+1)-r(nR)
            h1=r(nR)-r(nR-1)
            rInt=rInt+(h1+h2)/6.*( f(nR-1)*(two*h1-h2)/h1         +&
            &                      f(nR)  *(h1+h2)*(h1+h2)/(h1*h2)+&
            &                      f(nR+1)*(two*h2-h1)/h2 )
         end do
         rInt = rInt+half*(r(nRmax)-r(nRmax-1))*(f(nRmax)+f(nRmax-1))
         do nR=2,nRmax-1,2
            h2=r(nR+1)-r(nR)
            h1=r(nR)-r(nR-1)
            rInt=rInt+(h1+h2)/6.*( f(nR-1)*(two*h1-h2)/h1         +&
            &                      f(nR)  *(h1+h2)*(h1+h2)/(h1*h2)+&
            &                      f(nR+1)*(two*h2-h1)/h2 )
         end do
         rInt = half*rInt

      end if

   end function rInt_fd
!------------------------------------------------------------------------------
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
   real(cp) function rInt_R(f,n_r_max,n_cheb_max,dr_fac,chebt)
      !
      !   Same as function rInt but for a radial dependent mapping function
      !   dr_fac2.
      !

      !-- Input variables:
      integer,           intent(in) :: n_r_max,n_cheb_max
      real(cp),          intent(in) :: f(n_r_max)
      real(cp),          intent(in) :: dr_fac(n_r_max)
      type(costf_odd_t), intent(in) :: chebt
              
      !-- Local variables
      real(cp) :: f2(n_r_max)
      real(cp) :: work(n_r_max)
      integer :: nR,nCheb
                 
      !--- Integrals:
      do nR=1,n_r_max
         f2(nR)=f(nR)/dr_fac(nR)
      end do

      !-- Transform to cheb space:
      call chebt%costf1(f2,work)
      f2(1)      =half*f2(1)
      f2(n_r_max)=half*f2(n_r_max)

      !-- Sum contribution:
      rInt_R=f2(1)            ! This is zero order contribution
      do nCheb=3,n_cheb_max,2 ! Only even chebs contribute
         rInt_R=rInt_R-one/real(nCheb*(nCheb-2),cp)*f2(nCheb)
      end do

      !-- Remaining renormalisation:
      rInt_R=two*sqrt(two/real(n_r_max-1,cp))*rInt_R

   end function rInt_R
!------------------------------------------------------------------------------
end module integration
