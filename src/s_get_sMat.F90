!$Id$
!***********************************************************************
    SUBROUTINE get_Smat(dt,l,hdif,sMat,sPivot)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  |  Purpose of this subroutine is to contruct the time step matricies|
!  |  sMat(i,j) and s0mat for the entropy equation.                    |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+
      
    USE truncation
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE algebra, ONLY: sgefa

    IMPLICIT NONE

    REAL(kind=8) :: dt
    INTEGER :: l
    REAL(kind=8) :: hdif

!-- Output
    REAL(kind=8) :: sMat(n_r_max,n_r_max)
    INTEGER :: sPivot(n_r_max)

!-- Local variables:
    INTEGER :: info,nCheb,nR
    REAL(kind=8) :: O_dt,dLh

!-- end of declaration
!-----------------------------------------------------------------------

    O_dt=1.D0/dt

    dLh=DBLE(l*(l+1))

!----- Boundary coditions:
    DO nCheb=1,n_cheb_max
        IF ( ktops == 1 ) then
            sMat(1,nCheb)=cheb_norm
        ELSE
            sMat(1,nCheb)=cheb_norm*dcheb(nCheb,1)
        END IF
        IF ( kbots == 1 ) THEN
            sMat(n_r_max,nCheb)=cheb_norm*cheb(nCheb,n_r_max)
        ELSE
            sMat(n_r_max,nCheb)=cheb_norm*dcheb(nCheb,n_r_max)
        END IF
    END DO
    IF ( n_cheb_max < n_r_max ) THEN ! fill with zeros !
        DO nCheb=n_cheb_max+1,n_r_max
            sMat(1,nCheb)      =0.D0
            sMat(n_r_max,nCheb)=0.D0
        END DO
    END IF

!----- Other points:
    DO nCheb=1,n_r_max
        DO nR=2,n_r_max-1
            sMat(nR,nCheb)= cheb_norm * (                    &
                                       O_dt*cheb(nCheb,nR) - &
              alpha*opr*hdif*kappa(nR)*(  d2cheb(nCheb,nR) + &
              ( PolFac*beta(nR)+2.d0*or1(nR)+dLkappa(nR) ) * &
                                           dcheb(nCheb,nR) - &
                   dLh*or2(nR)*             cheb(nCheb,nR) ) )
        END DO
    END DO

!----- Factor for highest and lowest cheb:
    DO nR=1,n_r_max
        sMat(nR,1)      =0.5D0*sMat(nR,1)
        sMat(nR,n_r_max)=0.5D0*sMat(nR,n_r_max)
    END DO

!----- LU decomposition:
    CALL sgefa(sMat,n_r_max,n_r_max,sPivot,info)
    IF ( info /= 0 ) THEN
        WRITE(*,*) 'Singular matrix sMat!'
        STOP '31'
    END IF

            
    RETURN
    end SUBROUTINE get_Smat

!-----------------------------------------------------------------------------
