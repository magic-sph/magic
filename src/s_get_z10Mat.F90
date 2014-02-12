!$Id$
!***********************************************************************
#ifdef WITH_PRECOND_Z10
    SUBROUTINE get_z10Mat(dt,l,hdif,zMat,zPivot,zMat_fac)
#else
    SUBROUTINE get_z10Mat(dt,l,hdif,zMat,zPivot)
#endif
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------------------------------------------------------------+
!  |                                                                   |
!  |  Purpose of this subroutine is to construct and LU-decompose the  |
!  |  inversion matrix z10mat for the implicit time step of the       |
!  |  toroidal velocity potential z of degree l=1 and order m=0.       |
!  |  This differs from the the normal zmat only if either the ICB or  |
!  |  CMB have no-slip boundary condition and inner core or mantle are |
!  |  chosen to rotate freely (either kbotv=1 and/or ktopv=1).         |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+
      
    USE truncation
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE logic
    USE const
    USE algebra, ONLY: sgefa

    IMPLICIT NONE

    REAL(kind=8) :: dt
    INTEGER :: l
    REAL(kind=8) :: hdif

!-- Output: z10Mat and pivot_z10
    REAL(kind=8) :: zMat(n_r_max,n_r_max)
    INTEGER :: zPivot(n_r_max)
#ifdef WITH_PRECOND_Z10
    REAL(kind=8),INTENT(out) :: zMat_fac(n_r_max)
#endif

!-- local variables:
    INTEGER :: nR,nCheb,info
    REAL(kind=8) :: O_dt,dLh

!-- end of declaration
!-----------------------------------------------------------------------

    O_dt=1.D0/dt
    dLh=DBLE(l*(l+1))

!-- Boundary conditions:
               
    DO nCheb=1,n_cheb_max

    !----- CMB condition:
    !        Note opposite sign of vicouse torques (-dz+2 z /r )
    !        for CMB and ICB!

        IF ( ktopv == 1 ) THEN  ! free slip
            zMat(1,nCheb)=                   cheb_norm * &
                 ( (2.D0*or1(1)+beta(1))*cheb(nCheb,1) - &
                                        dcheb(nCheb,1) )
        ELSE IF ( ktopv == 2 ) THEN ! no slip
            IF ( l_SRMA ) THEN
                zMat(1,nCheb)= cheb_norm * c_z10_omega_ma*cheb(nCheb,1)
            ELSE IF ( l_rot_ma ) then
                zMat(1,nCheb)= cheb_norm *               ( &
                          c_dt_z10_ma*O_dt*cheb(nCheb,1) - &
                       alpha*( 2.d0*or1(1)*cheb(nCheb,1) - &
                                          dcheb(nCheb,1) ) )
            ELSE
                zMat(1,nCheb)= cheb_norm*cheb(nCheb,1)
            END IF
        END IF

    !----- ICB condition:
        IF ( kbotv == 1 ) THEN  ! free slip
            zMat(n_r_max,nCheb)=                          cheb_norm * &
            ( (2.D0*or1(n_r_max)+beta(n_r_max))*cheb(nCheb,n_r_max) - &
                                               dcheb(nCheb,n_r_max) )
        ELSE IF ( kbotv == 2 ) THEN ! no slip
            IF ( l_SRIC ) THEN
                zMat(n_r_max,nCheb)= cheb_norm * &
                                    c_z10_omega_ic*cheb(nCheb,n_r_max)
            ELSE if ( l_rot_ic ) THEN     !  time integration of z10
                zMat(n_r_max,nCheb)= cheb_norm *             (     &
                            c_dt_z10_ic*O_dt*cheb(nCheb,n_r_max) + &
                   alpha*( 2.D0*or1(n_r_max)*cheb(nCheb,n_r_max) - &
                                            dcheb(nCheb,n_r_max) ) )
            ELSE
                zMat(n_r_max,nCheb)= cheb_norm * cheb(nCheb,n_r_max)
            END IF
        END IF

    END DO

!----- Other points: (same as zMat)
    DO nCheb=1,n_r_max
        DO nR=2,n_r_max-1
            zMat(nR,nCheb)=                   cheb_norm * ( &
                          O_dt*dLh*or2(nR)*cheb(nCheb,nR) - &
                        alpha*hdif*dLh*visc(nR)*or2(nR) * ( &
                                         d2cheb(nCheb,nR) + &
             (dLvisc(nR)- beta(nR))*      dcheb(nCheb,nR) - &
            (dLvisc(nR)*beta(nR)+2.d0*dLvisc(nR)*or1(nR)  + &
              dLh*or2(nR)+dbeta(nR)+2.D0*beta(nR)*or1(nR))* &
                                           cheb(nCheb,nR) ) )

        END DO
    END DO

!-- Normalisation
    DO nR=1,n_r_max
        zMat(nR,1)      =0.5D0*zMat(nR,1)
        zMat(nR,n_r_max)=0.5D0*zMat(nR,n_r_max)
    END DO

!-- Fill up with zeros:
    DO nCheb=n_cheb_max+1,n_r_max
        zMat(1,nCheb)      =0.D0
        zMat(n_r_max,nCheb)=0.D0
    END DO

#ifdef WITH_PRECOND_Z10
  ! compute the linesum of each line
  DO nR=1,n_r_max
     zMat_fac(nR)=1.0D0/MAXVAL(ABS(zMat(nR,:)))
     zMat(nR,:) = zMat(nR,:)*zMat_fac(nR)
  END DO
#endif

!-- LU-decomposition of z10mat:
    CALL sgefa(zMat,n_r_max,n_r_max,zPivot,info)
    IF ( info /= 0 ) THEN
        WRITE(*,*) 'ERROR MESSAGE FROM SUBROUTINE GET_z10MAT:'
        WRITE(*,*) 'singular matrix z10Mat!'
        STOP
    END IF

    RETURN
    end SUBROUTINE get_z10Mat

!-----------------------------------------------------------------------
