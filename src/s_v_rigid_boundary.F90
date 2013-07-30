!$Id$
!***********************************************************************
    subroutine v_rigid_boundary(nR,omega,lDeriv,vrr,vtr,vpr, &
                           cvrr,dvrdtr,dvrdpr,dvtdpr,dvpdpr, &
                                                nThetaStart)
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
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+

    USE truncation
    USE radial_functions
    USE blocking
    USE horizontal_data

    IMPLICIT NONE

!-- Input of variables:
    INTEGER :: nR            ! no of radial grid point
    LOGICAL :: lDeriv         ! derivatives required ?
    INTEGER :: nThetaStart    ! no of theta to start with
            
!-- Input of boundary rotation rate
    REAL(kind=8) :: omega

!-- output:
    REAL(kind=8) :: vrr(nrp,nfs)
    REAL(kind=8) :: vpr(nrp,nfs)
    REAL(kind=8) :: vtr(nrp,nfs)
    REAL(kind=8) :: cvrr(nrp,nfs)
    REAL(kind=8) :: dvrdtr(nrp,nfs)
    REAL(kind=8) :: dvrdpr(nrp,nfs)
    REAL(kind=8) :: dvtdpr(nrp,nfs)
    REAL(kind=8) :: dvpdpr(nrp,nfs)

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
