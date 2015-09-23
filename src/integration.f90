module integration

   use precision_mod
   use const, only: half, one, two
   use cosine_transform, only: costf1

   implicit none

   private

   public :: rInt, rIntIC, rInt_R

contains

   real(cp) function rInt(f,nRmax,drFac,i_costf_init,d_costf_init)
      !
      !  This function performs the radial integral over a
      !  function f that is given on the appropriate nRmax
      !  radial Chebychev grid points.
      !  The arrays i_costf_init,d_costf_init are
      !  defined by calling init_costf1.
      !
      !  .. Note:: drFac maps radius to cheb space [-1,1]
      !            drFac=two/(rMax-rMin)
      !

      !-- Input variables:
      real(cp), intent(in) :: f(*)
      integer,  intent(in) :: nRmax
      real(cp), intent(in) :: drFac
      integer,  intent(in) :: i_costf_init(*)
      real(cp), intent(in) :: d_costf_init(*)

      !-- Local variables:
      integer :: nCheb,nR
      real(cp) :: WORK(nRmax)
      real(cp) :: f2(nRmax)

      !-- Save input:
      do nR=1,nRmax
         f2(nR)=f(nR)
      end do

      !-- Transform to cheb space:
      call costf1(f2,work,i_costf_init,d_costf_init)
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
   real(cp) function rIntIC(f,nRmax,drFac,i_costf_init,d_costf_init)
      !
      !  This function performs the radial integral over a
      !  function f that is given on the appropriate nRmax
      !  radial Chebychev grid points.
      !  The arrays i_costf_init,d_costf_init are
      !  defined by calling init_costf1.
      !

      !-- Input variables:
      real(cp), intent(inout) :: f(*)
      integer,  intent(in) :: nRmax
      real(cp), intent(in) :: drFac
      integer,  intent(in) :: i_costf_init(*)
      real(cp), intent(in) :: d_costf_init(*)

      !-- Local variables:
      integer :: nCheb,nChebInt
      real(cp) :: weight
      real(cp) :: work(nRmax)

      call costf1(f,work,i_costf_init,d_costf_init)
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
   real(cp) function rInt_R(f,n_r_max,n_cheb_max,dr_fac, &
                                    i_costf_init,d_costf_init)
      !
      !   Same as function rInt but for a radial dependent mapping function
      !   dr_fac2.
      !

      !-- Input variables:
      real(cp), intent(in) :: f(*)
      integer,  intent(in) :: n_r_max,n_cheb_max
      real(cp), intent(in) :: dr_fac(*)
      real(cp), intent(in) :: d_costf_init(*)
      integer,  intent(in) :: i_costf_init(*)
              
      !-- Local variables
      real(cp) :: f2(n_r_max)
      real(cp) :: work(n_r_max)
      integer :: nR,nCheb
                 
      !--- Integrals:
      do nR=1,n_r_max
          f2(nR)=f(nR)/dr_fac(nR)
      end do

      !-- Transform to cheb space:
      call costf1(f2,work,i_costf_init,d_costf_init)
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
