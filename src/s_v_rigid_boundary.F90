!$Id$
!***********************************************************************
    subroutine v_rigid_boundary(nR,omega,lDeriv,vrr,vtr,vpr, &
         &                      cvrr,dvrdtr,dvrdpr,dvtdpr,dvpdpr, &
         &                      nThetaStart)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------------------------------------------------------------+
!  |                                                                   |
!  |  Purpose of this subroutine is to set the velocities and their    |
!  |  derivatives at a fixed boundary.                                 |
!  |  While vt is zero, since we only allow for rotation about the     |
!  |  z-axix, vp= r sin(theta) v_phi = r**2 sin(theta)**2 omega        |
!  |  cvr= r**2 * radial component of (\curl v) =                      |
!  |       r**2  2 cos(theta) omega                                    |
!  |                                                                   |
!  +-------------------------------------------------------------------+

    USE truncation
    USE radial_functions
    USE blocking
    USE horizontal_data

    IMPLICIT NONE

!-- Input of variables:
    INTEGER,INTENT(in) :: nR            ! no of radial grid point
    LOGICAL,INTENT(in) :: lDeriv         ! derivatives required ?
    INTEGER,INTENT(in) :: nThetaStart    ! no of theta to start with
            
!-- Input of boundary rotation rate
    REAL(kind=8),INTENT(in) :: omega

!-- output:
    REAL(kind=8),INTENT(OUT) :: vrr(nrp,nfs)
    REAL(kind=8),INTENT(OUT) :: vpr(nrp,nfs)
    REAL(kind=8),INTENT(OUT) :: vtr(nrp,nfs)
    REAL(kind=8),INTENT(OUT) :: cvrr(nrp,nfs)
    REAL(kind=8),INTENT(OUT) :: dvrdtr(nrp,nfs)
    REAL(kind=8),INTENT(OUT) :: dvrdpr(nrp,nfs)
    REAL(kind=8),INTENT(OUT) :: dvtdpr(nrp,nfs)
    REAL(kind=8),INTENT(OUT) :: dvpdpr(nrp,nfs)

!-- local:
    REAL(kind=8) :: r2
    INTEGER :: nThetaCalc,nTheta,nThetaNHS
    INTEGER :: nPhi

!-- end of declaration
!----------------------------------------------------------------------------


    IF ( nR == n_r_cmb ) THEN
        r2=r_cmb*r_cmb
    ELSE IF ( nR == n_r_icb ) THEN
        r2=r_icb*r_icb
    ELSE
        WRITE(*,*)
        WRITE(*,*) '! v_rigid boundary called for a grid'
        WRITE(*,*) '! points which is not a boundary !  '
        RETURN
    END IF

    nThetaCalc=nThetaStart-1
    DO nTheta=1,sizeThetaB
        nThetaCalc=nThetaCalc+1
        nThetaNHS =(nThetaCalc+1)/2 ! northern hemisphere=odd n_theta
        DO nPhi=1,n_phi_max
            vrr(nPhi,nTheta)=0.D0
            vtr(nPhi,nTheta)=0.D0
            vpr(nPhi,nTheta)=r2*rho0(nR)*sn2(nThetaNHS)*omega
            IF ( lDeriv ) THEN
                cvrr(nPhi,nTheta)  = &
                    r2*rho0(nR)*2.D0*cosTheta(nThetaCalc)*omega
                dvrdtr(nPhi,nTheta)=0.D0
                dvrdpr(nPhi,nTheta)=0.D0
                dvtdpr(nPhi,nTheta)=0.D0
                dvpdpr(nPhi,nTheta)=0.D0
            END IF
        END DO
    END DO


    RETURN
    end subroutine v_rigid_boundary

!-------------------------------------------------------------------------
