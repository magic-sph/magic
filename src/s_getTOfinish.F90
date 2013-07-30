!$Id$
!***********************************************************************
    SUBROUTINE getTOfinish(nR,dtLast,zAS,dzAS,ddzAS,dzRstrLM, &
                           dzAstrLM,dzCorLM,dzLFLM)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!-----------------------------------------------------------------------
!  This program was previously part of s_getTO.f.
!  It has now been seperated to get it out of the theta-block loop.
!-----------------------------------------------------------------------

    USE truncation
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE torsional_oscillations
    USE blocking
    USE horizontal_data
    USE logic
    USE output_data

    IMPLICIT NONE

!-- Input of variables:
    INTEGER :: nR
    REAL(kind=8) :: dtLast

!-- Input of toroidal flow scalar and second radial derivative:
    REAL(kind=8) :: zAS(l_max+1)
    REAL(kind=8) :: dzAS(l_max+1) ! anelastic
    REAL(kind=8) :: ddzAS(l_max+1)

!-- Input of arrays from s_getTO.f:
    REAL(kind=8) :: dzRstrLM(l_max+2),dzAstrLM(l_max+2)
    REAL(kind=8) :: dzCorLM(l_max+2),dzLFLM(l_max+2)

!-- Local stuff:
    INTEGER :: l,lS,lA,lm

!-- End of declaration
!---------------------------------------------------------------------------


!------ When all thetas are done calculate viscous stress in LM space:
    dzStrLMr(1,nR) =0.D0
    dzRstrLMr(1,nR)=0.D0
    dzAstrLMr(1,nR)=0.D0
    dzCorLMr(1,nR) =0.D0
    dzLFLMr(1,nR)  =0.D0
    dzdVpLMr(1,nR) =0.D0
    dzddVpLMr(1,nR)=0.D0
    DO l=1,l_max
        lS=(l-1)+1
        lA=(l+1)+1
        lm=lm2(l,0)
        dzStrLMr(l+1,nR)= hdif_V(lm) * (                      &
                                                 ddzAS(l+1) - &
                                        beta(nR)* dzAS(l+1) - &
           (dLh(lm)*or2(nR)+dbeta(nR)+2.d0*beta(nR)*or1(nR))* &
                               zAS(l+1) )
    !---- -r**2/(l(l+1)) 1/sin(theta) dtheta sin(theta)**2
    !     minus sign to bring stuff on the RHS of NS equation !
        dzRstrLMr(l+1,nR)=-r(nR)*r(nR)/dLh(lm) * ( &
                       dTheta1S(lm)*dzRstrLM(lS) - &
                       dTheta1A(lm)*dzRstrLM(lA) )
        dzAstrLMr(l+1,nR)=-r(nR)*r(nR)/dLh(lm) * ( &
                       dTheta1S(lm)*dzAstrLM(lS) - &
                       dTheta1A(lm)*dzAstrLM(lA) )
        dzCorLMr(l+1,nR) =-r(nR)*r(nR)/dLh(lm) * ( &
                        dTheta1S(lm)*dzCorLM(lS) - &
                        dTheta1A(lm)*dzCorLM(lA) )
        dzLFLMr(l+1,nR)  = r(nR)*r(nR)/dLh(lm) * ( &
                         dTheta1S(lm)*dzLFLM(lS) - &
                         dTheta1A(lm)*dzLFLM(lA) )
        dzdVpLMr(l+1,nR) =(zAS(l+1)-dzdVpLMr(l+1,nR))/dtLast
        dzddVpLMr(l+1,nR)=(zAS(l+1)/dtLast+dzddVpLMr(l+1,nR))/dtLast
    END DO


    RETURN
    end SUBROUTINE getTOfinish

!-----------------------------------------------------------------------------
