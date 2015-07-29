!$Id$
!***********************************************************************
    SUBROUTINE get_RAS(Blm,Br,rT,nThetaStart,sizeThetaBlock)
!***********************************************************************

!  +-------------------------------------------------------------------+
!  |                                                                   |
!  |  Purpose of this subroutine is to calculate the axisymmetric      |
!  |  radial component Br of an axisymmetric ploidal field Blm         |
!  |  given in spherical harmonic space (l,m=0).                       |
!  |                                                                   |
!  +-------------------------------------------------------------------+

    USE truncation
    USE radial_functions
    USE num_param
    USE blocking
    USE horizontal_data
    USE logic

    IMPLICIT NONE

    INTEGER :: nThetaStart    ! first theta to be treated
    INTEGER :: sizeThetaBlock ! last theta
    REAL(kind=8) :: rT             ! radius
    COMPLEX(kind=8) :: Blm(lm_max_dtB)         ! field in (l,m)-space for rT

!-- Output:
    REAL(kind=8) :: Br(*)

!-- Local:
    INTEGER :: lm,l
    INTEGER :: nTheta,nThetaN
    REAL(kind=8) :: fac
    REAL(kind=8) :: sign
    REAL(kind=8) :: Br_1,Br_n,Br_s

!-- end of declaration
!----------------------------------------------------------------------

    fac=1.D0/(rT*rT)

    DO nTheta=1,sizeThetaBlock,2 ! loop over thetas in northers HS
        nThetaN=(nThetaStart+nTheta)/2

        sign=-1.d0
        Br_n=0.D0
        Br_s=0.D0
        DO l=0,l_max
            lm=lm2(l,0)
            sign=-sign
            Br_1=REAL(Blm(l+1))*DBLE(l*(l+1))*Plm(lm,nThetaN)
            Br_n=Br_n+Br_1
            Br_s=Br_s+sign*Br_1
        END DO  ! Loop over degree
        Br(nTheta)  =fac*Br_n
        Br(nTheta+1)=fac*Br_s

    END DO        ! Loop over colatitudes

    RETURN
    end SUBROUTINE get_RAS

!---------------------------------------------------------------------------
