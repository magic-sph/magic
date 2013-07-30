!$Id$
!***********************************************************************
    SUBROUTINE updateWP(w,dw,ddw,dwdt,dwdtLast, &
                          p,dp,dpdt,dpdtLast,s, &
      workA,workB,w1,coex,dt,nLMB,lRmsNext,nTh)
!***********************************************************************

!-----------------------------------------------------------------------

!  updates the poloidal velocity potential w, the pressure p,  and
!  their derivatives
!  adds explicit part to time derivatives of w and p

!-------------------------------------------------------------------------

    USE truncation
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE blocking
    USE horizontal_data
    USE logic
    USE matrices
    USE RMS
    USE algebra, ONLY: cgeslML

    IMPLICIT NONE

!-- Input/output of scalar fields:
    COMPLEX(kind=8) :: w(lm_max,n_r_max)
    COMPLEX(kind=8) :: dw(lm_max,n_r_max)
    COMPLEX(kind=8) :: ddw(lm_max,n_r_max)
    COMPLEX(kind=8) :: dwdt(lm_max,n_r_max)
    COMPLEX(kind=8) :: dwdtLast(lm_max,n_r_max)
    COMPLEX(kind=8) :: p(lm_max,n_r_max)
    COMPLEX(kind=8) :: dp(lm_max,n_r_max)
    COMPLEX(kind=8) :: dpdt(lm_max,n_r_max)
    COMPLEX(kind=8) :: dpdtLast(lm_max,n_r_max)
    COMPLEX(kind=8) :: s(lm_max,n_r_max)
!-- Output: updated w,dw,ddw,p,dp,dwdtLast,dpdtLast

!-- Input of other variables:
    REAL(kind=8) :: w1        ! weight for time step !
    REAL(kind=8) :: coex      ! factor depending on alpha
    REAL(kind=8) :: dt        ! time step
    INTEGER :: nLMB     ! block number
    LOGICAL :: lRmsNext
    INTEGER :: nTh      ! thread number

!-- Input of recycled work arrays:
    COMPLEX(kind=8) :: workA(lm_max,n_r_max)
    COMPLEX(kind=8) :: workB(lm_max,n_r_max)

!-- Local variables:
    REAL(kind=8) :: w2                  ! weight of second time step
    REAL(kind=8) :: O_dt
    INTEGER :: l1,m1              ! degree and order
    INTEGER :: lm1,lm,lmB         ! position of (l,m) in array
    INTEGER :: lmStart,lmStop ! max and min number of orders m
    INTEGER :: lmStart_real      ! range of lm for real array
    INTEGER :: lmStop_real       !
    INTEGER :: lmStart_00        ! excluding l=0,m=0
    INTEGER :: nLMB2
    INTEGER :: nR                ! counts radial grid points
    INTEGER :: n_cheb             ! counts cheb modes

    COMPLEX(kind=8) :: Dif(lm_max),Pre(lm_max),Buo(lm_max)

    COMPLEX(kind=8) :: rhs1(2*n_r_max,sizeLMB2max)

     
!-- end of declaration
!---------------------------------------------------------------------

    IF ( .NOT. l_update_v ) RETURN

    lmStart     =lmStartB(nLMB)
    lmStop      =lmStopB(nLMB)
    lmStart_00  =MAX(2,lmStart)
    lmStart_real=2*lmStart_00-1
    lmStop_real =2*lmStop

    w2  =1.D0-w1
    O_dt=1.D0/dt

    DO nLMB2=1,nLMBs2(nLMB)
        lmB=0
        DO lm=1,sizeLMB2(nLMB2,nLMB)
            lm1=lm22lm(lm,nLMB2,nLMB)
            l1 =lm22l(lm,nLMB2,nLMB)
            m1 =lm22m(lm,nLMB2,nLMB)
            IF ( l1 > 0 ) THEN
                IF ( .NOT. lWPmat(l1) ) THEN
                    CALL get_wpMat(dt,l1,hdif_V(lm1), &
                                   wpMat(1,1,l1),wpPivot(1,l1))
                    lWPmat(l1)=.TRUE.
                END IF
                lmB=lmB+1
                rhs1(1,lmB)        =0.D0
                rhs1(n_r_max,lmB)  =0.D0
                rhs1(n_r_max+1,lmB)=0.D0
                rhs1(2*n_r_max,lmB)=0.D0
                DO nR=2,n_r_max-1
                    rhs1(nR,lmB)=                         &
                        O_dt*dLh(lm1)*or2(nR)*w(lm1,nR) + &
                           rho0(nR)*agrav(nR)*s(lm1,nR) + &
                                        w1*dwdt(lm1,nR) + &
                                    w2*dwdtLast(lm1,nR)
                    rhs1(nR+n_r_max,lmB)=                 &
                      -O_dt*dLh(lm1)*or2(nR)*dw(lm1,nR) + &
                                        w1*dpdt(lm1,nR) + &
                                    w2*dpdtLast(lm1,nR)
                END DO
            END IF
        END DO
        IF ( lmB > 0 )                                     &
        CALL cgeslML(wpMat(1,1,l1),2*n_r_max,2*n_r_max,    &
                     wpPivot(1,l1),rhs1,2*n_r_max,lmB)

        IF ( lRmsNext ) THEN ! Store old w
            DO nR=1,n_r_max
                DO lm=1,sizeLMB2(nLMB2,nLMB)
                    lm1=lm22lm(lm,nLMB2,nLMB)
                    workB(lm1,nR)=w(lm1,nR)
                END DO
            END DO
        END IF

        lmB=0
        DO lm=1,sizeLMB2(nLMB2,nLMB)
            lm1=lm22lm(lm,nLMB2,nLMB)
            l1 =lm22l(lm,nLMB2,nLMB)
            m1 =lm22m(lm,nLMB2,nLMB)
            IF ( l1 > 0 ) THEN
                lmB=lmB+1
                IF ( m1 > 0 ) then
                    DO n_cheb=1,n_cheb_max
                        w(lm1,n_cheb)=rhs1(n_cheb,lmB)
                        p(lm1,n_cheb)=rhs1(n_r_max+n_cheb,lmB)
                    END DO
                ELSE
                    DO n_cheb=1,n_cheb_max
                        w(lm1,n_cheb)= &
                            CMPLX(REAL(rhs1(n_cheb,lmB)),0.D0,KIND=KIND(0d0))
                        p(lm1,n_cheb)= &
                            CMPLX(REAL(rhs1(n_r_max+n_cheb,lmB)),0.D0,KIND=KIND(0d0))
                    END DO
                END IF
            END IF
        END DO

    END DO   ! end of loop over lm1 blocks

!-- set cheb modes > n_cheb_max to zero (dealiazing)
    DO n_cheb=n_cheb_max+1,n_r_max
        DO lm1=lmStart_00,lmStop
            w(lm1,n_cheb)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            p(lm1,n_cheb)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        END DO
    END DO

!-- Transform to radial space and get radial derivatives
!   using dwdtLast, dpdtLast as work arrays:
    CALL costf1(   w,lm_max_real,lmStart_real,lmStop_real, &
                       dwdtLast,i_costf_init,d_costf_init)
    CALL get_dddr(                         w,dw,ddw,workA, &
                     lm_max_real,lmStart_real,lmStop_real, &
                     n_r_max,n_cheb_max,dwdtLast,dpdtLast, &
                 i_costf_init,d_costf_init,drx,ddrx,dddrx)
    CALL costf1(   p,lm_max_real,lmStart_real,lmStop_real, &
                       dwdtLast,i_costf_init,d_costf_init)
    CALL get_dr(p,dp,lm_max_real,lmStart_real,lmStop_real, &
                     n_r_max,n_cheb_max,dwdtLast,dpdtLast, &
                            i_costf_init,d_costf_init,drx)

!-- Calculate explicit time step part:
    IF ( ra /= 0.D0 ) THEN
        DO nR=n_r_cmb+1,n_r_icb-1
            DO lm1=lmStart_00,lmStop
                Dif(lm1)=   hdif_V(lm1)*dLh(lm1)*or2(nR)*visc(nR) * &
                           (                   ddw(lm1,nR)        + &
                       (2.d0*dLvisc(nR)-beta(nR)/3.d0)*dw(lm1,nR) - &
                 (dLh(lm1)*or2(nR)+4.d0/3.d0*(dbeta(nR)           + &
                                               dLvisc(nR)*beta(nR)+ &
                (3.d0*dLvisc(nR)+beta(nR))*or1(nR)))              * &
                                                        w(lm1,nR) )
                Pre(lm1)=-dp(lm1,nR)+beta(nR)*p(lm1,nR)
                Buo(lm1)=rho0(nR)*rgrav(nR)*s(lm1,nR)
                dwdtLast(lm1,nR)=dwdt(lm1,nR) - &
                                 coex*(Pre(lm1)+Buo(lm1)+Dif(lm1))
                dpdtLast(lm1,nR)=     dpdt(lm1,nR)  - coex*( &
                                dLh(lm1)*or2(nR)*p(lm1,nR) + &
                   hdif_V(lm1)*visc(nR)*dLh(lm1)*or2(nR) * ( &
                                   -workA(lm1,nR)          + &
                (beta(nR)-dLvisc(nR))*ddw(lm1,nR)          + &
                 ( dLh(lm1)*or2(nR)+dLvisc(nR)*beta(nR)    + &
                                                   dbeta(nR) &
                   +2.d0*(dLvisc(nR)+beta(nR))*or1(nR))    * &
                                       dw(lm1,nR)          - &
                                       dLh(lm1)*or2(nR)    * &
               (2.d0*or1(nR)+2.d0/3.d0*beta(nR)+dLvisc(nR))* &
                                        w(lm1,nR)          )  )
                IF ( lRmsNext ) THEN
                    workB(lm1,nR)=O_dt*dLh(lm1)*or2(nR) * &
                                  ( w(lm1,nR)-workB(lm1,nR) )
                    IF ( l_RMStest ) &
                        workB(lm1,nR)=workB(lm1,nR)-Dif(lm1)
                END IF
            END DO
            IF ( lRmsNext ) THEN
                CALL hInt2Pol(Dif,nR,lmStart_00,lmStop,DifPolLMr, &
                              DifPol2hInt(nR,nTh),DifPolAs2hInt(nR,nTh))
                CALL hInt2Pol(workB(1,nR),nR,lmStart_00,lmStop, &
                                                     dtVPolLMr, &
                      dtVPol2hInt(nR,nTh),dtVPolAs2hInt(nR,nTh))
            END IF
        END DO

    ELSE  ! no s-contribution !

        DO nR=n_r_cmb+1,n_r_icb-1
            DO lm1=lmStart_00,lmStop
                Dif(lm1)=     hdif_V(lm1)*dLh(lm1)*or2(nR)*visc(nR)* &
                            (                   ddw(lm1,nR)        + &
                        (2.d0*dLvisc(nR)-beta(nR)/3.d0)*dw(lm1,nR) - &
                  (dLh(lm1)*or2(nR)+4.d0/3.d0*(dbeta(nR)           + &
                                                dLvisc(nR)*beta(nR)+ &
                 (3.d0*dLvisc(nR)+beta(nR))*or1(nR)))              * &
                                                         w(lm1,nR) )
                Pre(lm1)=-dp(lm1,nR)+beta(nR)*p(lm1,nR)
                dwdtLast(lm1,nR)=dwdt(lm1,nR)  - &
                                 coex*(Pre(lm1)+Dif(lm1))
                dpdtLast(lm1,nR)=        dpdt(lm1,nR) - coex*( &
                                  dLh(lm1)*or2(nR)*p(lm1,nR) + &
                     hdif_V(lm1)*visc(nR)*dLh(lm1)*or2(nR) * ( &
                                     -workA(lm1,nR)          + &
                  (beta(nR)-dLvisc(nR))*ddw(lm1,nR)          + &
                   ( dLh(lm1)*or2(nR)+dLvisc(nR)*beta(nR)    + &
                                                     dbeta(nR) &
                     +2.d0*(dLvisc(nR)+beta(nR))*or1(nR))    * &
                                         dw(lm1,nR)          - &
                                         dLh(lm1)*or2(nR)    * &
                 (2.d0*or1(nR)+2.d0/3.d0*beta(nR)+dLvisc(nR))* &
                                          w(lm1,nR)          )  )
                IF ( lRmsNext ) THEN
                    workB(lm1,nR)=O_dt*dLh(lm1)*or2(nR) * &
                        ( w(lm1,nR)-workB(lm1,nR) )
                    IF ( l_RMStest ) &
                        workB(lm1,nR)=workB(lm1,nR)-Dif(lm1)
                END IF
            END DO
            IF ( lRmsNext ) THEN
                CALL hInt2Pol(Dif,nR,lmStart_00,lmStop,DifPolLMr, &
                              DifPol2hInt(nR,nTh),DifPolAs2hInt(nR,nTh))
                CALL hInt2Pol(workB(1,nR),nR,lmStart_00,lmStop, &
                                                     dtVPolLMr, &
                     dtVPol2hInt(nR,nTh),dtVpolAs2hInt(nR,nTh))
            END IF
        END DO


    END IF


!  Note: workA=dddw not needed beyond this point!
            
     
    RETURN
    end SUBROUTINE updateWP

!------------------------------------------------------------------------------
