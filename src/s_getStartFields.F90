!$Id$
!***********************************************************************
    SUBROUTINE getStartFields(time,dt,dtNew,n_time_step)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/15/02  by JW. --------------!

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  |  Purpose of this subroutine is to initialize the fields and       |
!  |  other auxiliary parameters.                                      |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+

    USE truncation
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE init_fields
    USE Grenoble
    USE blocking
    USE horizontal_data
    USE logic
    USE fields
    USE fieldsLast
    USE output_data
    USE const
    USE usefull, ONLY: cc2real

    IMPLICIT NONE

!---- Output variables:
    REAL(kind=8) :: time,dt,dtNew
    INTEGER :: n_time_step

!-- Local variables:
    INTEGER :: nR,lm,l1m0,nLMB,l,m
    INTEGER :: lmStart,lmStop,lmStartReal,lmStopReal
    REAL(kind=8) :: coex
    REAL(kind=8) :: d_omega_ma_dt,d_omega_ic_dt
    CHARACTER(len=76) :: message

    REAL(kind=8) :: sEA,sES,sAA

    REAL(kind=8) :: s0(n_r_max),ds0(n_r_max)
    REAL(kind=8) :: w1(n_r_max),w2(n_r_max)

    COMPLEX(kind=8) :: workA(lm_max,n_r_max)
    COMPLEX(kind=8) :: workB(lm_max,n_r_max)


!-- end of declaration
!---------------------------------------------------------------
     

!-- Start with setting fields to zero:
!   Touching the fields with the appropriate processor
!   for the LM-distribute parallel region (LMLoop) makes
!   sure that they are located close the individual
!   processors in memory:
     
    !$OMP PARALLEL DO SCHEDULE(STATIC,1)   &
    !$OMP  PRIVATE(lmStart,lmStop,lm,nR)
    DO nLMB=1,nLMBs ! Blocking of loop over all (l,m)

        lmStart=lmStartB(nLMB)
        lmStop =lmStopB(nLMB)
        DO lm=lmStart,lmStop

            IF ( l_conv ) THEN
                DO nR=1,n_r_max
                    w(lm,nR)       =zero
                    dw(lm,nR)      =zero
                    ddw(lm,nR)     =zero
                    dwdtLast(lm,nR)=zero
                    z(lm,nR)       =zero
                    dz(lm,nR)      =zero
                    dzdtLast(lm,nR)=zero
                    p(lm,nR)       =zero
                    dp(lm,nR)      =zero
                    dpdtLast(lm,nR)=zero
                END DO
            END IF
            IF ( l_heat ) THEN
                DO nR=1,n_r_max
                    s(lm,nR)       =zero
                    ds(lm,nR)      =zero
                    dsdtLast(lm,nR)=zero
                END DO
            END IF
            IF ( l_mag ) THEN
                DO nR=1,n_r_max
                    b(lm,nR)       =zero
                    db(lm,nR)      =zero
                    ddb(lm,nR)     =zero
                    dbdtLast(lm,nR)=zero
                    aj(lm,nR)      =zero
                    dj(lm,nR)      =zero
                    ddj(lm,nR)     =zero
                    djdtLast(lm,nR)=zero
                END DO
            END IF
            IF ( l_cond_ic ) THEN
                DO nR=1,n_r_ic_max
                    b_ic(lm,nR)       =zero
                    db_ic(lm,nR)      =zero
                    ddb_ic(lm,nR)     =zero
                    dbdt_icLast(lm,nR)=zero
                    aj_ic(lm,nR)      =zero
                    dj_ic(lm,nR)      =zero
                    ddj_ic(lm,nR)     =zero
                    djdt_icLast(lm,nR)=zero
                END DO
            END IF

        END DO

    END DO

    !$OMP   END PARALLEL DO   ! END OF SMP PARALLEL LOOP OVER LM blocks !


    IF ( l_start_file ) THEN
        CALL readStartFields(   w,dwdtLast,z,dzdtLast, &
                                p,dpdtLast,s,dsdtLast, &
                               b,dbdtLast,aj,djdtLast, &
                   b_ic,dbdt_icLast,aj_ic,djdt_icLast, &
                                    omega_ic,omega_ma, &
          lorentz_torque_icLast,lorentz_torque_maLast, &
                            time,dt,dtNew,n_time_step)
        IF ( dt > 0.D0 ) THEN
            WRITE(message,'(''! Using old time step:'',D16.6)') dt
        ELSE
            dt=dtMax
            WRITE(message,'(''! Using dtMax time step:'',D16.6)') dtMax
        END IF
    ELSE
        time =0.D0
        dt   =dtMax
        dtNew=dtMax
        n_time_step=0
        WRITE(message,'(''! Using dtMax time step:'',D16.6)') dtMax
    END IF
    CALL logWrite(message)

!---- Computations for the Nusselt number if we are anelastic
    IF (l_heat) THEN
        CALL s_cond(s0)
        CALL get_dr(s0,ds0,1,1,1,n_r_max,n_cheb_max, &
                    w1,w2,i_costf_init,d_costf_init,drx)
    END IF

! For computing the Nusselt number
    topcond=-1.D0/DSQRT(4.D0*pi)*ds0(1)
    botcond=-1.D0/DSQRT(4.D0*pi)*ds0(n_r_max)
!----- Get radial derivatives and initialize:

    !$OMP PARALLEL DO SCHEDULE(STATIC,1)  &
    !$OMP  PRIVATE(lmStart,lmStop,lmStartReal,lmStopReal)

    DO nLMB=1,nLMBs ! Blocking of loop over all (l,m)

        lmStart=lmStartB(nLMB)
        lmStop =lmStopB(nLMB)
        lmStartReal=2*lmStart-1
        lmStopReal =2*lmStop

    !----- Initialize/add magnetic field:
        IF ( ( imagcon /= 0 .OR. init_b1 /= 0 .OR. lGrenoble ) &
              .AND. ( l_mag .OR. l_mag_LF ) ) THEN
            CALL initB(b,aj,b_ic,aj_ic, &
                 lorentz_torque_icLast, &
                 lorentz_torque_maLast, &
                        lmStart,lmStop)
        END IF

    !----- Initialize/add velocity, set IC and ma rotation:
        IF ( l_conv .OR. l_mag_kin .OR. l_SRIC .OR. l_SRMA ) &
            CALL initV(w,z,omega_ic,omega_ma,lmStart,lmStop)

                
    !----- Initialize/add entropy:
        IF ( ( init_s1 /= 0 .OR. impS /= 0 ) .AND. l_heat ) &
            CALL initS(s,lmStart,lmStop)


        IF ( l_conv .OR. l_mag_kin ) THEN
            CALL get_ddr(      w,dw,ddw,lm_max_real, &
                             lmStartReal,lmStopReal, &
                     n_r_max,n_cheb_max,workA,workB, &
                 i_costf_init,d_costf_init,drx,ddrx)
            CALL get_dr(           z,dz,lm_max_real, &
                             lmStartReal,lmStopReal, &
                     n_r_max,n_cheb_max,workA,workB, &
                      i_costf_init,d_costf_init,drx)
        END IF

        IF ( l_mag .OR. l_mag_kin  ) THEN
            CALL get_ddr(      b,db,ddb,lm_max_real, &
                             lmStartReal,lmStopReal, &
                     n_r_max,n_cheb_max,workA,workB, &
                 i_costf_init,d_costf_init,drx,ddrx)
            CALL get_ddr(     aj,dj,ddj,lm_max_real, &
                             lmStartReal,lmStopReal, &
                     n_r_max,n_cheb_max,workA,workB, &
                 i_costf_init,d_costf_init,drx,ddrx)
        END IF
        IF ( l_cond_ic ) then
            CALL get_ddr_even(b_ic,db_ic,ddb_ic,lm_max_real, &
                                     lmStartReal,lmStopReal, &
             n_r_ic_max,n_cheb_ic_max,dr_fac_ic,workA,workB, &
                          i_costf1_ic_init,d_costf1_ic_init, &
                          i_costf2_ic_init,d_costf2_ic_init)
            CALL get_ddr_even(aj_ic,dj_ic,ddj_ic,lm_max_real, &
                                      lmStartReal,lmStopReal, &
              n_r_ic_max,n_cheb_ic_max,dr_fac_ic,workA,workB, &
                           i_costf1_ic_init,d_costf1_ic_init, &
                           i_costf2_ic_init,d_costf2_ic_init)
        END IF

        IF ( l_heat ) THEN
        !-- Get radial derivatives of entropy:
            CALL get_dr(              s,ds,lm_max_real, &
                                lmStartReal,lmStopReal, &
                        n_r_max,n_cheb_max,workA,workB, &
                         i_costf_init,d_costf_init,drx)
        END IF

    END DO ! Loop over LM blocks

    !$OMP   END PARALLEL DO   ! END OF SMP PARALLEL LOOP OVER LM blocks !

!--- Get symmetry properties of tops excluding l=m=0:
    sES=0.D0
    sEA=0.D0
    sAA=0.D0
    DO m=0,l_max,minc
        DO l=m,l_max
            IF ( l > 0 ) THEN
                IF ( MOD(l+m,2) == 0 ) THEN
                    sES=sES+cc2real(tops(l,m),m)
                ELSE
                    sEA=sEA+cc2real(tops(l,m),m)
                END IF
                IF ( m /= 0 ) sAA=sAA+cc2real(tops(l,m),m)
            END IF
        END DO
    END DO
    IF ( sEA+sES == 0 ) THEN
        WRITE(message,'(''! Only l=m=0 comp. in tops:'')')
        CALL logWrite(message)
    ELSE
        sEA=DSQRT(sEA/(sEA+sES))
        sAA=DSQRT(sAA/(sEA+sES))
        WRITE(message,'(''! Rel. RMS equ. asym. tops:'',D16.6)') sEA
        CALL logWrite(message)
        WRITE(message,'(''! Rel. RMS axi. asym. tops:'',D16.6)') sAA
        CALL logWrite(message)
    END IF

!----- Get changes in mantle and ic rotation rate:
    IF ( .NOT. l_mag_LF ) THEN
        lorentz_torque_icLast=0.D0
        lorentz_torque_maLast=0.D0
    END IF
    IF ( l_z10mat ) THEN
        l1m0=lm2(1,0)
        coex=-2.D0*(alpha-1.D0)
        IF ( .NOT. l_SRMA .AND. ktopv == 2 .AND. l_rot_ma ) THEN
            d_omega_ma_dt=LFfac*c_lorentz_ma*lorentz_torque_maLast
            d_omega_ma_dtLast=d_omega_ma_dt -              &
                   coex * ( 2.d0*or1(1)*REAL(z(l1m0,1)) - &
                                        REAL(dz(l1m0,1)) )
        END IF
        IF ( .NOT. l_SRIC .AND. kbotv == 2 .AND. l_rot_ic ) THEN
            d_omega_ic_dt=LFfac*c_lorentz_ic*lorentz_torque_icLast
            d_omega_ic_dtLast=                    d_omega_ic_dt + &
              coex * ( 2.D0*or1(n_r_max)*REAL(z(l1m0,n_r_max)) - &
                                         REAL(dz(l1m0,n_r_max)) )
        END IF
    ELSE
        d_omega_ma_dtLast=0.D0
        d_omega_ic_dtLast=0.D0
    END IF


    RETURN
    end SUBROUTINE getStartFields


!---------------------------------------------------------------------------
!-- end of subroutine getStartFields
!---------------------------------------------------------------------------
