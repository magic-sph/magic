!$Id$
!***********************************************************************
    SUBROUTINE j_cond(lm0,aj0,aj0_ic)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  |  Purpose of this subroutine is to solve the dIFfusion equation    |
!  |  for an initial toroidal magnetic field.                          |
!  |                                                                   |
!  |  CALLed in s_init_b.f                                             |
!  |  CALLs
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+

    USE truncation
    USE radial_functions
    USE physical_parameters
    USE init_fields
    USE horizontal_data
    USE logic
    USE matrices
    USE algebra, ONLY: cgesl,sgefa

    IMPLICIT NONE

!-- input:
    INTEGER :: lm0

!-- output:
    COMPLEX(kind=8) :: aj0(*)
    COMPLEX(kind=8) :: aj0_ic(*)

!-- local:
    INTEGER :: n_cheb,n_r,info,n_r_real
    COMPLEX(kind=8) :: rhs(n_r_tot)
    COMPLEX(kind=8) :: work_l(n_r_max)
    COMPLEX(kind=8) :: work_l_ic(n_r_ic_max)

!-- END of declaration
!--------------------------------------------------------------------------

    n_r_real=n_r_max
    IF ( l_cond_ic ) n_r_real=n_r_real+n_r_ic_max

!-- Set matrix:

!----- Outer core:

    DO n_cheb=1,n_cheb_max
        DO n_r=2,n_r_max-1
                            jMat(n_r,n_cheb,1)= cheb_norm * &
            hdIF_B(lm0)*dLh(lm0)*opm*lambda(n_r)*or2(n_r) * &
                    (                  d2cheb(n_cheb,n_r) + &
                          dLlambda(n_r)*dcheb(n_cheb,n_r) - &
                       dLh(lm0)*or2(n_r)*cheb(n_cheb,n_r) )
        END DO
    END DO
     
!----- boundary conditions:

!----- CMB:
    DO n_cheb=1,n_cheb_max     ! should be bpeaktop at CMB
        jMat(1,n_cheb,1)=cheb_norm*cheb(n_cheb,1)
    END DO
    DO n_cheb=n_cheb_max+1,n_r_max
        jMat(1,n_cheb,1)=0.d0
    END DO

!----- ICB:
    IF ( l_cond_ic ) then  ! matching condition at inner core:

        DO n_cheb=1,n_cheb_max
            jMat(n_r_max,n_cheb,1)=cheb_norm*cheb(n_cheb,n_r_max)
            jMat(n_r_max+1,n_cheb,1)= cheb_norm * &
                                sigma_ratio*dcheb(n_cheb,n_r_max)
        END DO
        DO n_cheb=n_cheb_max+1,n_r_max
            jMat(n_r_max,n_cheb,1)=0.d0
            jMat(n_r_max+1,n_cheb,1)=0.d0
        END DO

    ELSE

        DO n_cheb=1,n_cheb_max
            jMat(n_r_max,n_cheb,1)=cheb(n_cheb,n_r_max)*cheb_norm
        END DO
        DO n_cheb=n_cheb_max+1,n_r_max
            jMat(n_r_max,n_cheb,1)=0.d0
        END DO

    END IF
     
    DO n_r=1,n_r_max
        jMat(n_r,1,1)      =0.5*jMat(n_r,1,1)
        jMat(n_r,n_r_max,1)=0.5*jMat(n_r,n_r_max,1)
    END DO

!----- Inner core:

    IF ( l_cond_ic ) then

        DO n_cheb=1,n_r_ic_max ! counts even IC cheb modes
            DO n_r=2,n_r_ic_max ! counts IC radial grid point
                           jMat(n_r_max+n_r,n_r_max+n_cheb,1) = &
                cheb_norm_ic*dLh(lm0)*or3(n_r_max)*opm*O_sr * ( &
                              r_ic(n_r)*d2cheb_ic(n_cheb,n_r) + &
                         2.d0*D_lP1(lm0)*dcheb_ic(n_cheb,n_r) )
            END DO
        END DO

    !-------- boundary conditions:
        DO n_cheb=1,n_cheb_ic_max
            jMat(n_r_max,n_r_max+n_cheb,1)= &
                         -cheb_norm_ic*cheb_ic(n_cheb,1)
            jMat(n_r_max+1,n_r_max+n_cheb,1)= -cheb_norm_ic * ( &
                                           dcheb_ic(n_cheb,1) + &
                    D_lP1(lm0)*or1(n_r_max)*cheb_ic(n_cheb,1) )
        END DO
        DO n_cheb=n_r_max+n_cheb_ic_max+1,n_r_tot
            jMat(n_r_max,n_cheb,1)=0.d0
            jMat(n_r_max+1,n_cheb,1)=0.d0
        END DO

    !-------- normalization for lowest Cheb mode:
        DO n_r=n_r_max+1,n_r_tot
            jMat(n_r,n_r_max+1,1)=0.5*jMat(n_r,n_r_max+1,1)
        END DO

    !-------- fill matrix up with zeros:
        DO n_cheb=n_r_max+1,n_r_tot
            DO n_r=1,n_r_max-1
                jMat(n_r,n_cheb,1)=0.d0
            END DO
        END DO
        DO n_cheb=1,n_r_max
            DO n_r=n_r_max+2,n_r_tot
                jMat(n_r,n_cheb,1)=0.d0
            END DO
        END DO

    END IF ! conducting inner core ?

!----- invert matrix:
    CALL sgefa(jMat(1,1,1),n_r_tot,n_r_real,jPivot(1,1),info)
    IF ( info /= 0 ) then
        write(*,*) 'Singular matrix jMat in j_cond.'
        stop
    END IF
     
!----- zero RHS, except BC's
    DO n_r=2,n_r_real-1
        rhs(n_r)=CMPLX(0.d0,0.d0,KIND=KIND(0d0))
    END DO
    rhs(1)= bpeaktop                             ! Outer boundary
    IF ( .NOT. l_cond_ic ) rhs(n_r_max)=bpeakbot  ! Inner boundary
     
!----- solve linear system:
    CALL cgesl(jMat(1,1,1),n_r_tot,n_r_real,jPivot(1,1),rhs)

!----- copy result for OC:
    DO n_cheb=1,n_cheb_max
        aj0(n_cheb)=rhs(n_cheb)
    END DO
    DO n_cheb=n_cheb_max+1,n_r_max
        aj0(n_cheb)=CMPLX(0.d0,0.d0,KIND=KIND(0d0))
    END DO

!----- transform to radial space:
    CALL costf1(aj0,1,1,1,work_l,i_costf_init,d_costf_init)
     
    IF ( l_cond_ic ) then
         
    !----- copy result for IC:
        DO n_cheb=1,n_cheb_ic_max
            aj0_ic(n_cheb)=rhs(n_r_max+n_cheb)
        END DO
        DO n_cheb=n_cheb_ic_max+1,n_r_ic_max
            aj0_ic(n_cheb)=CMPLX(0.d0,0.d0,KIND=KIND(0d0))
        END DO
         
    !----- transform to radial space:
    !  Note: this is assuming that aj0_ic is an even function !
        CALL costf1(aj0_ic,1,1,1,work_l_ic, &
                    i_costf1_ic_init,d_costf1_ic_init)

    END IF
     
    RETURN
    END SUBROUTINE j_cond

!--------------------------------------------------------------------------------
