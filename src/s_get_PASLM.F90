!$Id$
!***********************************************************************
    SUBROUTINE get_PASLM(Tlm,Bp,rT,nThetaStart,sizeThetaBlock)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------------------------------------------------------------+
!  |                                                                   |
!  |  Purpose of this subroutine is to calculated the axisymmetric     |
!  |  phi component Bp of an axisymmetric toroidal field Tlm           |
!  |  given in spherical harmonic space (l,m=0).                       |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+

    USE truncation
    USE radial_functions
    USE num_param
    USE blocking
    USE horizontal_data
    USE logic

    IMPLICIT NONE

    INTEGER :: nThetaStart    ! first theta to be treated
    INTEGER :: sizeThetaBlock ! size of theta block
    REAL(kind=8) :: rT             ! radius
    COMPLEX(kind=8) :: Tlm(lm_max_dtB) ! field in (l,m)-space for rT

!-- Output:
    REAL(kind=8) :: Bp(*)

!-- Local:
    INTEGER :: lm,l
    INTEGER :: nTheta,nThetaN
    REAL(kind=8) :: fac
    REAL(kind=8) :: sign
    REAL(kind=8) :: Bp_1,Bp_n,Bp_s

!-- end of declaration
!----------------------------------------------------------------------

    DO nTheta=1,sizeThetaBlock,2 ! loop over thetas in northers HS

        nThetaN=(nThetaStart+nTheta)/2
        fac=osn1(nThetaN)/rT

        sign=-1.d0
        Bp_n=0.D0
        Bp_s=0.D0
        DO l=0,l_max
            lm=lm2(l,0)
            sign=-sign
            Bp_1=-REAL(Tlm(l+1))*dPlm(lm,nThetaN)
            Bp_n=Bp_n+Bp_1
            Bp_s=Bp_s-sign*Bp_1
        END DO  ! Loop over degree
        Bp(nTheta)  =fac*Bp_n
        Bp(nTheta+1)=fac*Bp_s

    END DO        ! Loop over colatitudes

    RETURN
    end SUBROUTINE get_PASLM
!---------------------------------------------------------------------------
