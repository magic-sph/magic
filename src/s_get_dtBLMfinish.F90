!$Id$
!***********************************************************************
SUBROUTINE get_dtBLMfinish(time,n_time_step, &
     &                     omega_ic, &
     &                     b,ddb,aj,dj,ddj,b_ic,db_ic,ddb_ic, &
     &                     aj_ic,dj_ic,ddj_ic)
  !***********************************************************************

  USE truncation
  USE radial_functions,ONLY: O_r_ic,lambda,or2,dLlambda,i_costf_init,d_costf_init,&
       &drx,or1
  USE physical_parameters,ONLY: opm,O_sr
  USE blocking, ONLY:lo_map,st_map
  USE horizontal_data, ONLY: dPhi,D_lP1,dLh,hdif_B
  USE logic,ONLY: l_cond_ic,l_DTrMagSpec
  USE dtB_mod,ONLY: PdifLM,TdifLM,PstrLM,TstrLM,TomeLM,PadvLM,TadvLM,&
       & PadvLMIC,TadvLMIC,PdifLMIC,TdifLMIC,TstrRLM,TadvRLM,TomeRLM,&
       & PdifLM_LMloc,TdifLM_LMloc,PadvLMIC_LMloc,TadvLMIC_LMloc,PdifLMIC_LMloc,TdifLMIC_LMloc,&
       &dtB_gather_Rloc_on_rank0
  USE LMLoop_data, ONLY: llmMag,ulmMag,llm,ulm,llm_real,ulm_real
  USE communications,ONLY: gather_all_from_lo_to_rank0,gt_OC,gt_IC
  USE parallel_mod,only: rank
  IMPLICIT NONE

  !-- Input of variables:
  REAL(kind=8),intent(IN) :: time
  INTEGER,intent(IN) :: n_time_step
  REAL(kind=8),intent(IN) :: omega_ic

  !-- Input of scalar fields:
  COMPLEX(kind=8),intent(IN) :: b(llmMag:ulmMag,n_r_maxMag)
  COMPLEX(kind=8),intent(IN) :: ddb(llmMag:ulmMag,n_r_maxMag)
  COMPLEX(kind=8),intent(IN) :: aj(llmMag:ulmMag,n_r_maxMag)
  COMPLEX(kind=8),intent(IN) :: dj(llmMag:ulmMag,n_r_maxMag)
  COMPLEX(kind=8),intent(IN) :: ddj(llmMag:ulmMag,n_r_maxMag)
  COMPLEX(kind=8),intent(IN) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
  COMPLEX(kind=8),intent(IN) :: db_ic(llmMag:ulmMag,n_r_ic_maxMag)
  COMPLEX(kind=8),intent(IN) :: ddb_ic(llmMag:ulmMag,n_r_ic_maxMag)
  COMPLEX(kind=8),intent(IN) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
  COMPLEX(kind=8),intent(IN) :: dj_ic(llmMag:ulmMag,n_r_ic_maxMag)
  COMPLEX(kind=8),intent(IN) :: ddj_ic(llmMag:ulmMag,n_r_ic_maxMag)

  !-- Local variables:
  INTEGER :: nLMB,nR    ! position of degree and order
  INTEGER :: lmStart,lmStop! limits of lm-block
  INTEGER :: lmStart_real,lmStop_real

  COMPLEX(kind=8) :: workA(lm_max,n_r_max),workB(lm_max,n_r_max)
  COMPLEX(kind=8),DIMENSION(lm_max) :: temp_global

  INTEGER :: l,m,lm
  !-- end of declaration
  !-----------------------------------------------------------------------
  
  ! gathering TstrRLM,TadvRLM and TomeRLM on rank0,
  ! they are then in st_map order
  call dtB_gather_Rloc_on_rank0

  !DO nLMB=1,nLMBs ! Blocking of loop over all (l,m)

  !   lmStart=lmStartB(nLMB)
  !   lmStop =lmStopB(nLMB)
  !   lmStart_real=2*lmStart-1
  !   lmStop_real =2*lmStop

  IF ( l_cond_ic ) THEN
     DO nR=1,n_r_ic_max
        DO lm=llm,ulm
           l=lo_map%lm2l(lm)
           m=lo_map%lm2m(lm)
           PadvLMIC_LMloc(lm,nR)=-omega_ic*dPhi(st_map%lm2(l,m))*b_ic(lm,nR)
           TadvLMIC_LMloc(lm,nR)=-omega_ic*dPhi(st_map%lm2(l,m))*aj_ic(lm,nR)
           PdifLMIC_LMloc(lm,nR)=opm*O_sr * ( ddb_ic(lm,nR) + &
                2.D0*D_lP1(st_map%lm2(l,m))*O_r_ic(nR)*db_ic(lm,nR) )
           TdifLMIC_LMloc(lm,nR)=opm*O_sr * ( ddj_ic(lm,nR) + &
                2.D0*D_lP1(st_map%lm2(l,m))*O_r_ic(nR)*dj_ic(lm,nR) )
        END DO
     END DO
  END IF

  DO nR=1,n_r_max
     DO lm=llm,ulm
        l=lo_map%lm2l(lm)
        m=lo_map%lm2m(lm)
        PdifLM_LMloc(lm,nR)= opm*lambda(nR)*hdif_B(st_map%lm2(l,m)) * &
             (ddb(lm,nR)-dLh(st_map%lm2(l,m))*or2(nR)*b(lm,nR))
        TdifLM_LMloc(lm,nR)= opm*lambda(nR)*hdif_B(st_map%lm2(l,m)) * &
             ( ddj(lm,nR) + dLlambda(nR)*dj(lm,nR) - &
             dLh(st_map%lm2(l,m))*or2(nR)*aj(lm,nR) )
     END DO
  END DO

  IF (rank.EQ.0) THEN
     CALL get_drNS(TstrRLM,workA,lm_max_real, &
          1,lm_max_real, &
          n_r_max,n_cheb_max,workB, &
          i_costf_init,d_costf_init,drx)

     DO nR=1,n_r_max
        DO lm=1,lm_max
           TstrLM(lm,nR)=TstrLM(lm,nR)+or1(nR)*workA(lm,nR)
        END DO
     END DO
     
     CALL get_drNS(TomeRLM,workA,lm_max_real, &
          1,lm_max_real, &
          n_r_max,n_cheb_max,workB, &
          i_costf_init,d_costf_init,drx)
     
     DO nR=1,n_r_max
        DO lm=1,lm_max
           TomeLM(lm,nR)=TomeLM(lm,nR)+or1(nR)*workA(lm,nR)
        END DO
     END DO
     
     CALL get_drNS(TadvRLM,workA,lm_max_real, &
          1,lm_max_real, &
          n_r_max,n_cheb_max,workB, &
          i_costf_init,d_costf_init,drx)
     
     DO nR=1,n_r_max
        DO lm=1,lm_max
           TadvLM(lm,nR)=TadvLM(lm,nR)+or1(nR)*workA(lm,nR)
        END DO
     END DO
  END IF
  !END DO

  ! PdifLM and TdifLM need to be gathered over lm
  CALL gather_all_from_lo_to_rank0(gt_OC,PdifLM_LMloc,PdifLM)
  CALL gather_all_from_lo_to_rank0(gt_OC,TdifLM_LMloc,TdifLM)
     
  IF ( l_DTrMagSpec .AND. n_time_step > 1 ) THEN

     ! also gather PadvLMIC,TadvLMIC,PdifLMIC and TdifLMIC
     CALL gather_all_from_lo_to_rank0(gt_IC,PadvLMIC_LMloc,PadvLMIC)
     CALL gather_all_from_lo_to_rank0(gt_IC,TadvLMIC_LMloc,TadvLMIC)
     CALL gather_all_from_lo_to_rank0(gt_IC,PdifLMIC_LMloc,PdifLMIC)
     CALL gather_all_from_lo_to_rank0(gt_IC,TdifLMIC_LMloc,TdifLMIC)

     IF (rank.EQ.0) THEN
        CALL rBrSpec(time,PstrLM,PadvLMIC,'rBrProSpec',.FALSE.,st_map)
        CALL rBrSpec(time,PadvLM,PadvLMIC,'rBrAdvSpec',.TRUE.,st_map)
        CALL rBrSpec(time,PdifLM,PdifLMIC,'rBrDifSpec',.TRUE.,st_map)
        DO nR=1,n_r_max
           DO lm=1,lm_max
              PstrLM(lm,nR)=PstrLM(lm,nR)-PadvLM(lm,nR)
           END DO
        END DO
        CALL rBrSpec(time,PstrLM,PadvLMIC,'rBrDynSpec',.FALSE.,st_map)

        CALL rBpSpec(time,TstrLM,TadvLMIC,'rBpProSpec',.FALSE.,st_map)
        CALL rBpSpec(time,TadvLM,TadvLMIC,'rBpAdvSpec',.TRUE.,st_map)
        CALL rBpSpec(time,TdifLM,TdifLMIC,'rBpDifSpec',.TRUE.,st_map)
        DO nR=1,n_r_max
           DO lm=1,lm_max
              TstrLM(lm,nR)=TstrLM(lm,nR)-TadvLM(lm,nR)
           END DO
        END DO
        CALL rBpSpec(time,TstrLM,TadvLMIC,'rBpDynSpec',.FALSE.,st_map)
     END IF
  END IF


  RETURN
end SUBROUTINE get_dtBLMfinish

!-----------------------------------------------------------------------
