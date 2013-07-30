!$Id$
!********************************************************************
SUBROUTINE get_PolTorRms(Pol,drPol,Tor, &
     &                   PolRms,TorRms,PolAsRms,TorAsRms,map)
  !********************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/17/02  by JW. --------------!

  !--------------------------------------------------------------------

  !  calculates integral PolRms=sqrt( Integral (pol^2 dV) )
  !  calculates integral TorRms=sqrt( Integral (tor^2 dV) )
  !  plus axisymmetric parts.
  !  integration in theta,phi by summation of spherical harmonics
  !  integration in r by using Chebycheff integrals

  !  The mapping map gives the mapping lm to l,m for the input
  !  arrays Pol,drPol and Tor

  !  Output: PolRms,TorRms,PolAsRms,TorAsRms

  !--------------------------------------------------------------------

  USE truncation
  USE radial_functions
  USE blocking
  USE horizontal_data
  USE const
  USE usefull, ONLY: cc2real
  USE integration, ONLY: rInt_R
  USE LMmapping, only: mappings
  IMPLICIT NONE

  !-- Input:
  COMPLEX(kind=8),INTENT(IN) :: Pol(lm_max,n_r_max)   ! Poloidal field Potential
  COMPLEX(kind=8),INTENT(IN) :: drPol(lm_max,n_r_max) ! Radial derivative of Pol
  COMPLEX(kind=8),INTENT(IN) :: Tor(lm_max,n_r_max)   ! Toroidal field Potential
  TYPE(mappings),INTENT(IN) :: map

  !-- Output:
  REAL(kind=8),INTENT(OUT) :: PolRms,PolAsRms
  REAL(kind=8),INTENT(OUT) :: TorRms,TorAsRms

  !-- Local:
  REAL(kind=8) :: PolRmsTemp,TorRmsTemp
  REAL(kind=8) :: PolRms_r(n_r_max)
  REAL(kind=8) :: TorRms_r(n_r_max)
  REAL(kind=8) :: PolAsRms_r(n_r_max)
  REAL(kind=8) :: TorAsRms_r(n_r_max)

  INTEGER :: n_r,lm,l,m
  REAL(kind=8) :: fac

  !-- end of declaration
  !---------------------------------------------------------------------

  DO n_r=1,n_r_max

     PolRms_r(n_r)  =0.D0
     TorRms_r(n_r)  =0.D0
     PolAsRms_r(n_r)=0.D0
     TorAsRms_r(n_r)=0.D0

     DO lm=2,lm_max
        l=map%lm2l(lm)
        m=map%lm2m(lm)
        PolRmsTemp= dLh(st_map%lm2(l,m)) * (                       &
             dLh(st_map%lm2(l,m))*or2(n_r)*cc2real(Pol(lm,n_r),m) + &
             cc2real(drPol(lm,n_r),m) )
        TorRmsTemp=   dLh(st_map%lm2(l,m))*cc2real(Tor(lm,n_r),m)
        IF ( m == 0 ) THEN  ! axisymmetric part
           PolAsRms_r(n_r)=PolAsRms_r(n_r) + PolRmsTemp
           TorAsRms_r(n_r)=TorAsRms_r(n_r) + TorRmsTemp
        ELSE
           PolRms_r(n_r)  =PolRms_r(n_r)   + PolRmsTemp
           TorRms_r(n_r)  =TorRms_r(n_r)   + TorRmsTemp
        END IF
     END DO    ! do loop over lms in block
     PolRms_r(n_r)=PolRms_r(n_r) + PolAsRms_r(n_r)
     TorRms_r(n_r)=TorRms_r(n_r) + TorAsRms_r(n_r)
  END DO    ! radial grid points

  !-- Radial Integrals:
  PolRms  =rInt_R(PolRms_r  ,n_r_max,n_r_max,drx, &
       i_costf_init,d_costf_init)
  TorRms  =rInt_R(TorRms_r  ,n_r_max,n_r_max,drx, &
       i_costf_init,d_costf_init)
  PolAsRms=rInt_R(PolAsRms_r,n_r_max,n_r_max,drx, &
       i_costf_init,d_costf_init)
  TorAsRms=rInt_R(TorAsRms_r,n_r_max,n_r_max,drx, &
       i_costf_init,d_costf_init)
  fac=1.D0/vol_oc
  PolRms  =SQRT(fac*PolRms)
  TorRms  =SQRT(fac*TorRms)
  PolAsRms=SQRT(fac*PolAsRms)
  TorAsRms=SQRT(fac*TorAsRms)


  RETURN
end SUBROUTINE get_PolTorRms

!-----------------------------------------------------------------------------
