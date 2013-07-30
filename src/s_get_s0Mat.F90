!$Id$
!***********************************************************************
    SUBROUTINE get_s0Mat(dt,sMat,sPivot)
!***********************************************************************

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  |  Purpose of this subroutine is to contruct the time step matrix   |
!  |  sMat0                                                            |
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

!-- Output
    REAL(kind=8) :: sMat(n_r_max,n_r_max)
    INTEGER :: sPivot(n_r_max)

!-- Local variables:
    INTEGER :: info,nCheb,nR
    REAL(kind=8) :: O_dt

!-- end of declaration
!-----------------------------------------------------------------------

    O_dt=1.D0/dt

!----- Boundary condition:
    DO nCheb=1,n_cheb_max

        IF ( ktops == 1 ) THEN
        !--------- Constant entropy at CBM:
            sMat(1,nCheb)=cheb_norm
        ELSE
        !--------- Constant flux at CBM:
            sMat(1,nCheb)=cheb_norm*dcheb(nCheb,1)
        END IF
        IF ( kbots == 1 ) THEN
        !--------- Constant entropy at ICB:
            sMat(n_r_max,nCheb)=cheb_norm*cheb(nCheb,n_r_max)
        ELSE
        !--------- Constant flux at ICB:
            sMat(n_r_max,nCheb)=cheb_norm*dcheb(nCheb,n_r_max)
        END IF
    END DO
    IF ( n_cheb_max < n_r_max ) THEN ! fill with zeros !
        DO nCheb=n_cheb_max+1,n_r_max
            sMat(1,nCheb)      =0.D0
            sMat(n_r_max,nCheb)=0.D0
        END DO
    END IF

    DO nCheb=1,n_r_max
        DO nR=2,n_r_max-1
            sMat(nR,nCheb)= cheb_norm * (                  &
                                     O_dt*cheb(nCheb,nR) - & 
               alpha*opr*kappa(nR)*(    d2cheb(nCheb,nR) + &
               (PolFac*beta(nR)+2.D0*or1(nR)+dLkappa(nR))* &
                                         dcheb(nCheb,nR) ) )
        END DO
    END DO

!----- Factors for highest and lowest cheb mode:
    DO nR=1,n_r_max
        sMat(nR,1)      =0.5D0*sMat(nR,1)
        sMat(nR,n_r_max)=0.5D0*sMat(nR,n_r_max)
    END DO

!---- LU decomposition:
    CALL sgefa(sMat,n_r_max,n_r_max,sPivot,info)
    IF ( info /= 0 ) THEN
        WRITE(*,*) '! Singular matrix sMat0!'
        STOP '28'
    END IF


    RETURN
    end SUBROUTINE get_s0Mat

!-----------------------------------------------------------------------------
