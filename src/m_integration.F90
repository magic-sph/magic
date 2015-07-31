!$Id$
module integration

   implicit none

   private

   public :: rInt, rIntIC, rInt_R

contains

   real(kind=8) function rInt(f,nRmax,drFac,i_costf_init,d_costf_init)
      !-------------------------------------------------------------
      !  This function performs the radial integral over a
      !  function f that is given on the appropriate nRmax
      !  radial Chebychev grid points.
      !  The arrays i_costf_init,d_costf_init are
      !  defined by calling init_costf1.
      !  Note: drFac maps radius to cheb space [-1,1]
      !        drFac=2.D0/(rMax-rMin)
      !-------------------------------------------------------------

      !-- Input variables:
      real(kind=8), intent(in) :: f(*)
      integer,      intent(in) :: nRmax
      real(kind=8), intent(in) :: drFac
      integer,      intent(in) :: i_costf_init(*)
      real(kind=8), intent(in) :: d_costf_init(*)

      !-- Local variables:
      integer :: nCheb,nR
      real(kind=8) :: WORK(nRmax)
      real(kind=8) :: f2(nRmax)

      !-- Save input:
      do nR=1,nRmax
         f2(nR)=f(nR)
      end do

      !-- Transform to cheb space:
      call costf1(f2,1,1,1,work,i_costf_init,d_costf_init)
      f2(1)    =0.5D0*f2(1)
      f2(nRmax)=0.5D0*f2(nRmax)

      !-- Sum contribution:
      rInt=f2(1)           ! This is zero order contribution
      do nCheb=3,nRmax,2  ! Only even chebs contribute
         rInt=rInt-1.D0/dble(nCheb*(nCheb-2))*f2(nCheb)
      end do

      !-- Remaining renormalisation:
      rInt=2.D0/drFac*dsqrt(2.D0/dble(nRmax-1))*rInt

   end function rInt
!------------------------------------------------------------------------------
   real(kind=8) function rIntIC(f,nRmax,drFac,i_costf_init,d_costf_init)
      !----------------------------------------------------------------
      !  This function performs the radial integral over a
      !  function f that is given on the appropriate nRmax
      !  radial Chebychev grid points.
      !  The arrays i_costf_init,d_costf_init are
      !  defined by calling init_costf1.
      !----------------------------------------------------------------

      !-- Input variables:
      real(kind=8), intent(inout) :: f(*)
      integer,      intent(in) :: nRmax
      real(kind=8), intent(in) :: drFac
      integer,      intent(in) :: i_costf_init(*)
      real(kind=8), intent(in) :: d_costf_init(*)

      !-- Local variables:
      integer :: nCheb,nChebInt
      real(kind=8) :: weight
      real(kind=8) :: work(nRmax)

      call costf1(f,1,1,1,work,i_costf_init,d_costf_init)
      f(1)    =0.5D0*f(1)
      f(nRmax)=0.5D0*f(nRmax)

      !-- Sum contribution:
      rIntIC=f(1)           ! This is zero order contribution
      do nCheb=2,nRmax      ! Only even chebs for IC
         nChebInt=2*nCheb-1
         weight  =-1.D0/dble(nChebInt*(nChebInt-2))
         rIntIC  =rIntIC+weight*f(nCheb)
      end do

      !-- Remaining renormalisation:
      rIntIC=dsqrt(2.D0/dble(nRmax-1))*rIntIC/drFac

   end function rIntIC
!------------------------------------------------------------------------------
   real(kind=8) function rInt_R(f,n_r_max,n_cheb_max,dr_fac, &
                                    i_costf_init,d_costf_init)
      !----------------------------------------------------------------------
      !   Same as function rInt but for a radial dependent mapping function
      !   dr_fac2.
      !----------------------------------------------------------------------

      !-- Input variables:
      real(kind=8), intent(in) :: f(*)
      integer,      intent(in) :: n_r_max,n_cheb_max
      real(kind=8), intent(in) :: dr_fac(*)
      real(kind=8), intent(in) :: d_costf_init(*)
      integer,      intent(in) :: i_costf_init(*)
              
      !-- Local variables
      real(kind=8) :: f2(n_r_max)
      real(kind=8) :: work(n_r_max)
      integer :: nR,nCheb
                 
      !--- Integrals:
      do nR=1,n_r_max
          f2(nR)=f(nR)/dr_fac(nR)
      end do

      !-- Transform to cheb space:
      call costf1(f2,1,1,1,work,i_costf_init,d_costf_init)
      f2(1)      =0.5D0*f2(1)
      f2(n_r_max)=0.5D0*f2(n_r_max)

      !-- Sum contribution:
      rInt_R=f2(1)            ! This is zero order contribution
      do nCheb=3,n_cheb_max,2 ! Only even chebs contribute
         rInt_R=rInt_R-1.D0/dble(nCheb*(nCheb-2))*f2(nCheb)
      end do

      !-- Remaining renormalisation:
      rInt_R=2.D0*dsqrt(2.D0/dble(n_r_max-1))*rInt_R

   end function rInt_R
!------------------------------------------------------------------------------
end module integration
