!$Id$
!**********************************************************************
    SUBROUTINE getTOnext(zAS,br,bt,bp,lTONext,lTONext2,    &
                 dt,dtLast,nR,nThetaStart,nThetaBlockSize, &
                                     BsLast,BpLast,BzLast)
!***********************************************************************

!-----------------------------------------------------------------------
!  Preparing TO calculation by storing flow and magnetic field
!  contribution to build time derivative.
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
    REAL(kind=8) :: dt,dtLast
    INTEGER :: nR
    INTEGER :: nThetaStart
    INTEGER :: nThetaBlockSize
    LOGICAL :: lTONext,lTONext2

!-- Input of toriodal flow scalar:
    REAL(kind=8) :: zAS(l_max+1)
!-- Input of field components:
    REAL(kind=8) :: br(nrp,nfs),bt(nrp,nfs),bp(nrp,nfs)

!-- Output:
    REAL(kind=8) :: BsLast(n_phi_maxStr,n_theta_maxStr,n_r_maxStr)
    REAL(kind=8) :: BpLast(n_phi_maxStr,n_theta_maxStr,n_r_maxStr)
    REAL(kind=8) :: BzLast(n_phi_maxStr,n_theta_maxStr,n_r_maxStr)

!-- Local:
    INTEGER :: l,lm
    INTEGER :: nTheta,nThetaBlock
    INTEGER :: nPhi

    REAL(kind=8) :: sinT,cosT
    REAL(kind=8) :: BsF1,BsF2,BpF1,BzF1,BzF2


!-- End of declaration
!---------------------------------------------------------------------------

    IF ( lVerbose ) WRITE(*,*) '! Starting getTOnext!',dtLast

    nTheta=nThetaStart-1

    IF ( lTONext2 .AND. nThetaStart == 1 ) THEN

        dzddVpLMr(1,nR)=0.D0
        DO l=1,l_max
            lm=lm2(l,0)
            dzddVpLMr(l+1,nR)=zAS(l+1)
        END DO

    ELSE IF ( lTOnext ) THEN

        DO nThetaBlock=1,nThetaBlockSize
            nTheta=nTheta+1
            sinT =sinTheta(nTheta)
            cosT =cosTheta(nTheta)
            BsF1 =sinT*or2(nR)
            BsF2 =cosT/sinT*or1(nR)
            BpF1 =or1(nR)/sinT
            BzF1 =cosT*or2(nR)
            BzF2 =or1(nR)

        !--- Get zonal means of velocity and derivatives:
            DO nPhi=1,n_phi_max
                BsLast(nPhi,nTheta,nR)=BsF1*br(nPhi,nThetaBlock) + &
                                       BsF2*bt(nPhi,nThetaBlock)
                BpLast(nPhi,nTheta,nR)=BpF1*bp(nPhi,nThetaBlock)
                BzLast(nPhi,nTheta,nR)=BzF1*br(nPhi,nThetaBlock) - &
                                       BzF2*bt(nPhi,nThetaBlock)
            END DO
                          
        END DO ! Loop over thetas in block !
                
        IF ( nThetaStart == 1 ) THEN
            dzdVpLMr(1,nR) =0.D0
            dzddVpLMr(1,nR)=0.D0
            DO l=1,l_max
                lm=lm2(l,0)
                dzdVpLMr(l+1,nR) = zAS(l+1)
                dzddVpLMr(l+1,nR)= ( dzddVpLMr(l+1,nR) - &
                                   ((dtLast+dt)/dt)*zAS(l+1) )/dtLast
            END DO
        END IF

    END IF

    IF ( lVerbose ) WRITE(*,*) '! End of getTOnext!'

    RETURN
    end SUBROUTINE getTOnext

!-----------------------------------------------------------------------------
