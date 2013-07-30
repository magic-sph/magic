!$Id$
!***********************************************************************
    subroutine get_b_nl_bcs(bc,br_vt_lm,br_vp_lm, &
                            lm_min_b,lm_max_b,b_nl_bc,aj_nl_bc)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------------------------------------------------------------+
!  |                                                                   |
!  |  Purpose of this subroutine is to calculate the nonlinear term    |
!  |  of the magnetic boundary condition for a conducting mantle in    |
!  |  physical space (phi,theta), assuming that the conductance        |
!  |  of the mantle is much smaller than that of the core.             |
!  |  Calculation is performed for the theta block:                    |
!  |          n_theta_min<=n_theta<=n_theta_min+n_theta_block-1        |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+
      
    USE truncation
    USE physical_parameters
    USE blocking
    USE horizontal_data

    IMPLICIT NONE

!-- input:
    CHARACTER(len=3),INTENT(IN) :: bc             ! Distinguishes 'CMB' and 'ICB'
    INTEGER,INTENT(IN) :: lm_min_b,lm_max_b  ! limits of lm-block
    COMPLEX(kind=8),INTENT(IN) :: br_vt_lm(lmP_max)     ! [br*vt/(r**2*sin(theta)**2)]
    COMPLEX(kind=8),INTENT(IN) :: br_vp_lm(lmP_max)     ! [br*vp/(r**2*sin(theta)**2)

!-- output:
    COMPLEX(kind=8),INTENT(OUT) :: b_nl_bc(lm_min_b:lm_max_b)  ! nonlinear bc for b
    COMPLEX(kind=8),INTENT(OUT) :: aj_nl_bc(lm_min_b:lm_max_b) ! nonlinear bc for aj

!-- local variables:
    INTEGER :: l,m       ! degree and order
    INTEGER :: lm        ! position of degree and order
    INTEGER :: lmP       ! same as lm but for l running to l_max+1
    INTEGER :: lmPS,lmPA ! lmP for l-1 and l+1
    REAL(kind=8) :: fac

!-- end of declaration
!-----------------------------------------------------------------------
    write(*,"(2A)") "In get_b_nl_bcs with bc=",bc

    IF ( bc == 'CMB' ) THEN

        fac=conductance_ma*prmag

        DO lm=lm_min_b,lm_max_b
            l   =lm2l(lm)
            m   =lm2m(lm)
            lmP =lm2lmP(lm)
            lmPS=lmP2lmPS(lmP)
            lmPA=lmP2lmPA(lmP)
            IF ( l > m ) THEN
                b_nl_bc(lm)= fac/dLh(lm) * (    &
                  dTheta1S(lm)*br_vt_lm(lmPS) - &
                  dTheta1A(lm)*br_vt_lm(lmPA) + &
                  dPhi(lm)*br_vp_lm(lmP)   )
                aj_nl_bc(lm)=-fac/dLh(lm) * (   &
                  dTheta1S(lm)*br_vp_lm(lmPS) - &
                  dTheta1A(lm)*br_vp_lm(lmPA) - &
                  dPhi(lm)*br_vt_lm(lmP)    )
            ELSE IF ( l == m ) THEN
                b_nl_bc(lm)= fac/dLh(lm) * (    &
                - dTheta1A(lm)*br_vt_lm(lmPA) + &
                  dPhi(lm)*br_vp_lm(lmP)   )
                aj_nl_bc(lm)=-fac/dLh(lm) * (   &
                - dTheta1A(lm)*br_vp_lm(lmPA) - &
                  dPhi(lm)*br_vt_lm(lmP)   )
            END IF
        END DO

    ELSE IF ( bc == 'ICB' ) THEN

        fac=sigma_ratio*prmag

        Do lm=lm_min_b,lm_max_b
            l   =lm2l(lm)
            m   =lm2m(lm)
            lmP =lm2lmP(lm)
            lmPS=lmP2lmPS(lmP)
            lmPA=lmP2lmPA(lmP)
            IF ( l > m ) THEN
                aj_nl_bc(lm)=-fac/dLh(lm) * (   &
                  dTheta1S(lm)*br_vp_lm(lmPS) - &
                  dTheta1A(lm)*br_vp_lm(lmPA) - &
                  dPhi(lm)*br_vt_lm(lmP)    )
            ELSE IF ( l == m ) THEN
                aj_nl_bc(lm)=-fac/dLh(lm) * (   &
                - dTheta1A(lm)*br_vp_lm(lmPA) - &
                  dPhi(lm)*br_vt_lm(lmP)    )
            END IF
        END DO

    ELSE
        WRITE(*,*) 'Wrong input of bc into s_get_b_nl_bcs.f'
        STOP
    END IF
            
    RETURN
    end subroutine get_b_nl_bcs

!-----------------------------------------------------------------------
