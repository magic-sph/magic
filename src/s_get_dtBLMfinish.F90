!$Id$
!***********************************************************************
    SUBROUTINE get_dtBLMfinish(time,n_time_step, &
               TstrRLM,TadvRLM,TomeRLM,omega_ic, &
              b,ddb,aj,dj,ddj,b_ic,db_ic,ddb_ic, &
                             aj_ic,dj_ic,ddj_ic)
!***********************************************************************

!  Parallelization note: this routine accesses the r-distributed
!  fields TstrRLM,TadvRLM,TomeRLM, and the fields in c_dtB.v
!  and also the LM-distributed scalar fields  b,ddb,aj,dj,ddj,b_ic,
!  db_ic,ddb_ic,aj_ic,dj_ic,ddj_ic

!  +-------------------------------------------------------------------+
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+
      
    USE truncation
    USE radial_functions
    USE physical_parameters
    USE blocking
    USE horizontal_data
    USE logic
    USE dtB_mod

    IMPLICIT NONE

!-- Input of constant parameters:
! include 'truncation.f'    ! contains mins,nmaf,nlaf
! include 'c_horizontal.f'  ! contains
! include 'c_phys_param.f'  ! includes conductance_ma,prmag,sigma_ratio
! include 'c_radial.f'      ! includes radial functions
! include 'c_logic.f'
! include 'c_blocking.f'

!-- Input of variables:
    REAL(kind=8) :: time
    INTEGER :: n_time_step
!-- Input from stuff calculated in s_get_dtB.f in s_radialLoopG.f:
    COMPLEX(kind=8) :: TstrRLM(lm_max_dtB,n_r_max_dtB)
    COMPLEX(kind=8) :: TadvRLM(lm_max_dtB,n_r_max_dtB)
    COMPLEX(kind=8) :: TomeRLM(lm_max_dtB,n_r_max_dtB)
    REAL(kind=8) :: omega_ic

!-- Input of scalar fields:
    COMPLEX(kind=8) :: b(lm_maxMag,n_r_maxMag)
    COMPLEX(kind=8) :: ddb(lm_maxMag,n_r_maxMag)
    COMPLEX(kind=8) :: aj(lm_maxMag,n_r_maxMag)
    COMPLEX(kind=8) :: dj(lm_maxMag,n_r_maxMag)
    COMPLEX(kind=8) :: ddj(lm_maxMag,n_r_maxMag)
    COMPLEX(kind=8) :: b_ic(lm_maxMag,n_r_ic_maxMag)
    COMPLEX(kind=8) :: db_ic(lm_maxMag,n_r_ic_maxMag)
    COMPLEX(kind=8) :: ddb_ic(lm_maxMag,n_r_ic_maxMag)
    COMPLEX(kind=8) :: aj_ic(lm_maxMag,n_r_ic_maxMag)
    COMPLEX(kind=8) :: dj_ic(lm_maxMag,n_r_ic_maxMag)
    COMPLEX(kind=8) :: ddj_ic(lm_maxMag,n_r_ic_maxMag)

!-- Input/Output:
!   Parallelization note:
!   On input the following fields in c_dtB.f are R-distributed
!   and have to be collected onto the processor
!   executing this routine:
!            PstrLM,PadvLM,TstrLM,TadvLM,TomeLM
! include 'c_dtB.f'

!-- Local variables:
    INTEGER :: nLMB,lm,nR    ! position of degree and order
    INTEGER :: lmStart,lmStop! limits of lm-block
    INTEGER :: lmStart_real,lmStop_real

    COMPLEX(kind=8) :: workA(lm_max,n_r_max),workB(lm_max,n_r_max)

!-- end of declaration
!-----------------------------------------------------------------------


    DO nLMB=1,nLMBs ! Blocking of loop over all (l,m)

        lmStart=lmStartB(nLMB)
        lmStop =lmStopB(nLMB)
        lmStart_real=2*lmStart-1
        lmStop_real =2*lmStop

        IF ( l_cond_ic ) THEN
            DO nR=1,n_r_ic_max
                DO lm=lmStart,lmStop
                    PadvLMIC(lm,nR)=-omega_ic*dPhi(lm)*b_ic(lm,nR)
                    TadvLMIC(lm,nR)=-omega_ic*dPhi(lm)*aj_ic(lm,nR)
                    PdifLMIC(lm,nR)=opm*O_sr * ( ddb_ic(lm,nR) + &
                        2.D0*D_lP1(lm)*O_r_ic(nR)*db_ic(lm,nR) )
                    TdifLMIC(lm,nR)=opm*O_sr * ( ddj_ic(lm,nR) + &
                        2.D0*D_lP1(lm)*O_r_ic(nR)*dj_ic(lm,nR) )
                END DO
            END DO
        END IF

        DO nR=1,n_r_max
            DO lm=lmStart,lmStop
                PdifLM(lm,nR)= opm*lambda(nR)*hdif_B(lm) * &
                    (ddb(lm,nR)-dLh(lm)*or2(nR)*b(lm,nR))
                TdifLM(lm,nR)= opm*lambda(nR)*hdif_B(lm) * &
                   ( ddj(lm,nR) + dLlambda(nR)*dj(lm,nR) - &
                               dLh(lm)*or2(nR)*aj(lm,nR) )
            END DO
        END DO

        CALL get_drNS(TstrRLM,workA,lm_max_real, &
                       lmStart_real,lmStop_real, &
                       n_r_max,n_cheb_max,workB, &
                  i_costf_init,d_costf_init,drx)

        DO nR=1,n_r_max
            DO lm=lmStart,lmStop
                TstrLM(lm,nR)=TstrLM(lm,nR)+or1(nR)*workA(lm,nR)
            END DO
        END DO

        CALL get_drNS(TomeRLM,workA,lm_max_real, &
                       lmStart_real,lmStop_real, &
                      n_r_max,n_cheb_max,workB, &
                  i_costf_init,d_costf_init,drx)

        DO nR=1,n_r_max
            DO lm=lmStart,lmStop
                TomeLM(lm,nR)=TomeLM(lm,nR)+or1(nR)*workA(lm,nR)
            END DO
        END DO

        CALL get_drNS(TadvRLM,workA,lm_max_real, &
                       lmStart_real,lmStop_real, &
                       n_r_max,n_cheb_max,workB, &
                  i_costf_init,d_costf_init,drx)

        DO nR=1,n_r_max
            DO lm=lmStart,lmStop
                TadvLM(lm,nR)=TadvLM(lm,nR)+or1(nR)*workA(lm,nR)
            END DO
        END DO

    END DO


    IF ( l_DTrMagSpec .AND. n_time_step > 1 ) THEN
        CALL rBrSpec(time,PstrLM,PadvLMIC,'rBrProSpec',.FALSE.)
        CALL rBrSpec(time,PadvLM,PadvLMIC,'rBrAdvSpec',.TRUE.)
        CALL rBrSpec(time,PdifLM,PdifLMIC,'rBrDifSpec',.TRUE.)
        DO nR=1,n_r_max
            DO lm=1,lm_max
                PstrLM(lm,nR)=PstrLM(lm,nR)-PadvLM(lm,nR)
            END DO
        END DO
        CALL rBrSpec(time,PstrLM,PadvLMIC,'rBrDynSpec',.FALSE.)
        CALL rBpSpec(time,TstrLM,TadvLMIC,'rBpProSpec',.FALSE.)
        CALL rBpSpec(time,TadvLM,TadvLMIC,'rBpAdvSpec',.TRUE.)
        CALL rBpSpec(time,TdifLM,TdifLMIC,'rBpDifSpec',.TRUE.)
        DO nR=1,n_r_max
            DO lm=1,lm_max
                TstrLM(lm,nR)=TstrLM(lm,nR)-TadvLM(lm,nR)
            END DO
        END DO
        CALL rBpSpec(time,TstrLM,TadvLMIC,'rBpDynSpec',.FALSE.)
    END IF


    RETURN
    end SUBROUTINE get_dtBLMfinish

!-----------------------------------------------------------------------
