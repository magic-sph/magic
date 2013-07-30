!$Id$
!*************************************************************************
    SUBROUTINE updateS(s,ds,dVSrLM,dsdt,dsdtLast, &
                       workA,workB,w1,coex,dt,nLMB)
!*************************************************************************

!-------------------------------------------------------------------------

!  updates the entropy field s and its radial derivatives
!  adds explicit part to time derivatives of s

!-------------------------------------------------------------------------

    USE truncation
    USE radial_functions
    USE physical_parameters
    USE init_fields
    USE blocking
    USE horizontal_data
    USE logic
    USE matrices
    USE output_data
    USE const
    USE algebra, ONLY: cgeslML,sgesl

    IMPLICIT NONE

!-- Input of variables:
    REAL(kind=8) :: w1        ! weight for time step !
    REAL(kind=8) :: coex      ! factor depending on alpha
    REAL(kind=8) :: dt        ! time step
    INTEGER :: nLMB

!-- Input/output of scalar fields:
    COMPLEX(kind=8) :: s(lm_max,n_r_max)
    COMPLEX(kind=8) :: ds(lm_max,n_r_max)
    COMPLEX(kind=8) :: dVSrLM(lm_max,n_r_max)
    COMPLEX(kind=8) :: dsdt(lm_max,n_r_max)
    COMPLEX(kind=8) :: dsdtLast(lm_max,n_r_max)
!-- Output: udpated s,ds,dsdtLast

!-- Input of recycled work arrays:
    COMPLEX(kind=8) :: workA(lm_max,n_r_max)
    COMPLEX(kind=8) :: workB(lm_max,n_r_max)

!-- Input/output of time stepping matricies:
! include 'c_mat.f'
            
!-- Local variables:
    REAL(kind=8) :: w2                  ! weight of second time step
    REAL(kind=8) :: O_dt
    INTEGER :: l1,m1              ! degree and order
    INTEGER :: lm1,lmB,lm         ! position of (l,m) in array
    INTEGER :: lmStart_real      ! range of lm for real array
    INTEGER :: lmStop_real       !
    INTEGER :: lmStart,lmStop
    INTEGER :: nLMB2
    INTEGER :: nR                 ! counts radial grid points
    INTEGER :: n_cheb             ! counts cheb modes
    REAL(kind=8) ::  rhs(n_r_max)       ! real RHS for l=m=0
    COMPLEX(kind=8) :: rhs1(2*n_r_max,sizeLMB2max) ! comples RHS for l>0

!-- end of declaration
!---------------------------------------------------------------------
     
    IF ( .NOT. l_update_s ) RETURN

    lmStart     =lmStartB(nLMB)
    lmStop      =lmStopB(nLMB)
    lmStart_real=2*lmStart-1
    lmStop_real =2*lmStop
    w2  =1.-w1
    O_dt=1.D0/dt


!--- Finish calculation of dsdt:
    CALL get_drNS(                       dVSrLM,workA, &
                 lm_max_real,lmStart_real,lmStop_real, &
                             n_r_max,n_cheb_max,workB, &
                        i_costf_init,d_costf_init,drx)

    DO nR=1,n_r_max
        DO lm=lmStart,lmStop
            dsdt(lm,nR)=orho1(nR)*(dsdt(lm,nR)-or2(nR)*workA(lm,nR))
        END DO
    END DO

    DO nLMB2=1,nLMBs2(nLMB)
        lmB=0
        DO lm=1,sizeLMB2(nLMB2,nLMB)
            lm1=lm22lm(lm,nLMB2,nLMB)
            l1 =lm22l(lm,nLMB2,nLMB)
            m1 =lm22m(lm,nLMB2,nLMB)
            IF ( l1 == 0 ) THEN
                IF ( .NOT. lSmat(l1) ) THEN
                    CALL get_s0Mat(dt,s0Mat,s0Pivot)
                    lSmat(l1)=.TRUE.
                END IF
                rhs(1)=      REAL(tops(0,0))
                rhs(n_r_max)=REAL(bots(0,0))
                DO nR=2,n_r_max-1
                    rhs(nR)=REAL(s(1,nR))*O_dt+ &
                          w1*REAL(dsdt(1,nR)) + &
                          w2*REAL(dsdtLast(1,nR))
                END DO
                CALL sgesl(s0Mat,n_r_max,n_r_max,s0Pivot,rhs)
            ELSE
                IF ( .NOT. lSmat(l1) ) THEN
                    CALL get_sMat(dt,l1,hdif_S(lm1), &
                                  sMat(1,1,l1),sPivot(1,l1))
                    lSmat(l1)=.TRUE.
                END IF
                lmB=lmB+1
                rhs1(1,lmB)=      tops(l1,m1)
                rhs1(n_r_max,lmB)=bots(l1,m1)
                DO nR=2,n_r_max-1
                    rhs1(nR,lmB)=s(lm1,nR)*O_dt + &
                                w1*dsdt(lm1,nR) + &
                                w2*dsdtLast(lm1,nR)
                END DO
            END IF
        END DO
        IF ( lmB > 0 ) &
        CALL cgeslML(sMat(1,1,l1),n_r_max,n_r_max, &
                     sPivot(1,l1),rhs1,2*n_r_max,lmB)
        lmB=0
        DO lm=1,sizeLMB2(nLMB2,nLMB)
            lm1=lm22lm(lm,nLMB2,nLMB)
            l1 =lm22l(lm,nLMB2,nLMB)
            m1 =lm22m(lm,nLMB2,nLMB)
            IF ( l1 == 0 ) THEN
                DO n_cheb=1,n_cheb_max
                    s(1,n_cheb)=rhs(n_cheb)
                END DO
            ELSE
                lmB=lmB+1
                IF ( m1 > 0 ) THEN
                    DO n_cheb=1,n_cheb_max
                        s(lm1,n_cheb)=rhs1(n_cheb,lmB)
                    END DO
                ELSE
                    DO n_cheb=1,n_cheb_max
                        s(lm1,n_cheb)= &
                             CMPLX(REAL(rhs1(n_cheb,lmB)),0.D0,KIND=KIND(0d0))
                    END DO
                END IF
            END IF
        END DO

    END DO     ! loop over lm blocks
     
!-- set cheb modes > n_cheb_max to zero (dealiazing)
    DO n_cheb=n_cheb_max+1,n_r_max
        DO lm1=lmStart,lmStop
            s(lm1,n_cheb)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        END DO
    END DO

!-- Get radial derivatives of s: workA,dsdtLast used as work arrays
    CALL costf1(s,lm_max_real,lmStart_real,lmStop_real, &
                dsdtLast,i_costf_init,d_costf_init)
    CALL get_ddr(s,ds,workA,lm_max_real,lmStart_real,lmStop_real, &
                               n_r_max,n_cheb_max,workB,dsdtLast, &
                              i_costf_init,d_costf_init,drx,ddrx)

!-- Calculate explicit time step part:
    DO nR=n_r_cmb+1,n_r_icb-1
        DO lm1=lmStart,lmStop
            dsdtLast(lm1,nR)=dsdt(lm1,nR) -                            &
                                    coex*opr*hdif_S(lm1) * kappa(nR) * &
                              (                        workA(lm1,nR) + &
               (PolFac*beta(nR)+2.D0*or1(nR)+dLkappa(nR))*ds(lm1,nR) - &
                            dLh(lm1)*or2(nR)*              s(lm1,nR) )
        END DO
    END DO

!-- workA=dds not needed further after this point, used as work array later


    RETURN
    end SUBROUTINE updateS

!-------------------------------------------------------------------------------
