!$Id$
MODULE fields_average_mod
  USE truncation
  USE radial_functions
  USE num_param
  USE blocking,ONLY: lmStartB,lmStopB
  USE horizontal_data
  USE logic
  use kinetic_energy
  use magnetic_energy
  USE output_data
  use parallel_mod
#if (FFTLIB==JW)
  USE fft_JW
#elif (FFTLIB==MKL)
  USE fft_MKL
#endif
  USE const
  USE LMLoop_data,ONLY: llm,ulm,llm_real,ulm_real,llmMag,ulmMag
  USE communications, ONLY: get_global_sum, gather_from_lo_to_rank0,&
       &gather_all_from_lo_to_rank0,gt_OC,gt_IC
  USE write_special,only: write_Bcmb
  IMPLICIT NONE

  COMPLEX(kind=8),ALLOCATABLE :: w_ave(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: z_ave(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: s_ave(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: b_ave(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: aj_ave(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: b_ic_ave(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: aj_ic_ave(:,:)

  contains
    SUBROUTINE initialize_fields_average_mod

      ALLOCATE( w_ave(llm:ulm,n_r_max) )
      ALLOCATE( z_ave(llm:ulm,n_r_max) )
      ALLOCATE( s_ave(llm:ulm,n_r_max) )
      ALLOCATE( b_ave(llm:ulm,n_r_max) )
      ALLOCATE( aj_ave(llm:ulm,n_r_max) )
      ALLOCATE( b_ic_ave(llm:ulm,n_r_ic_max) )
      ALLOCATE( aj_ic_ave(llm:ulm,n_r_ic_max) )

    END SUBROUTINE initialize_fields_average_mod
  !**********************************************************************
    SUBROUTINE fields_average(nAve,l_stop_time,                        &
       &                      time_passed,time_norm,omega_ic,omega_ma, &
       &                      w,z,s,b,aj,b_ic,aj_ic)
    !***********************************************************************

    !   !------------ This is release 2 level 2  --------------!
    !   !------------ Created on 2/6/02  by JW. --------------!

    !  +-------------+----------------+------------------------------------+
    !  |                                                                   |
    !  |  This subroutine averages fields b and v over time.               |
    !  |                                                                   |
    !  +-------------------------------------------------------------------+

    !-- Input of variables:
    INTEGER,intent(IN) :: nAve    ! number for averaged time steps
    LOGICAL,intent(IN) :: l_stop_time     ! true if this is the last time step
    REAL(kind=8),intent(IN) :: time_passed     ! time passed since last log
    REAL(kind=8),INTENT(IN) :: time_norm       ! time passed since start of time loop
    REAL(kind=8),intent(IN) :: omega_ic,omega_ma

    !-- Input of scalar fields:
    COMPLEX(kind=8),INTENT(IN) :: w(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(IN) :: z(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(IN) :: s(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(IN) :: b(llmMag:ulmMag,n_r_maxMag)
    COMPLEX(kind=8),intent(IN) :: aj(llmMag:ulmMag,n_r_maxMag)
    COMPLEX(kind=8),intent(IN) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
    COMPLEX(kind=8),intent(IN) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)


    !-- Local stuff:
    ! fields for the gathering
    COMPLEX(kind=8),DIMENSION(1:lm_maxMag,n_r_ic_maxMag) :: b_ic_ave_global,&
         & db_ic_ave_global,ddb_ic_ave_global,aj_ic_ave_global,dj_ic_ave_global
    COMPLEX(kind=8),DIMENSION(1:lm_maxMag,n_r_maxMag) :: b_ave_global
    COMPLEX(kind=8),DIMENSION(1:lm_max) :: db_ave_global,aj_ave_global,&
         & w_ave_global,dw_ave_global,z_ave_global, s_ave_global

    !----- Time averaged fields:
    COMPLEX(kind=8) :: dw_ave(llm:ulm,n_r_max)
    COMPLEX(kind=8) :: ds_ave(llm:ulm,n_r_max)
    COMPLEX(kind=8) :: db_ave(llm:ulm,n_r_max)
    COMPLEX(kind=8) :: db_ic_ave(llm:ulm,n_r_ic_max)
    COMPLEX(kind=8) :: ddb_ic_ave(llm:ulm,n_r_ic_max)
    COMPLEX(kind=8) :: dj_ic_ave(llm:ulm,n_r_ic_max)

    !----- Work array:
    COMPLEX(kind=8) :: workA(lm_max,n_r_max)
    COMPLEX(kind=8) :: workA_LMloc(llm:ulm,n_r_max)

    !----- Fields in grid space:
    COMPLEX(kind=8) :: Br(ncp,nfs),Bt(ncp,nfs),Bp(ncp,nfs) ! B field comp.
    COMPLEX(kind=8) :: Vr(ncp,nfs),Vt(ncp,nfs),Vp(ncp,nfs) ! B field comp.
    COMPLEX(kind=8) :: Sr(ncp,nfs)                         ! entropy

    !----- Help arrays for fields:
    COMPLEX(kind=8) :: dLhb(lm_max),bhG(lm_max),bhC(lm_max)
    COMPLEX(kind=8) :: dLhw(lm_max),vhG(lm_max),vhC(lm_max)

    !----- Energies of time average field:
    REAL(kind=8) :: ekinR(n_r_max)     ! kinetic energy w radius
    REAL(kind=8) :: e_kin_p_ave,e_kin_t_ave
    REAL(kind=8) :: e_kin_p_as_ave,e_kin_t_as_ave
    REAL(kind=8) :: e_mag_p_ave,e_mag_t_ave
    REAL(kind=8) :: e_mag_p_as_ave,e_mag_t_as_ave
    REAL(kind=8) :: e_mag_p_ic_ave,e_mag_t_ic_ave
    REAL(kind=8) :: e_mag_p_as_ic_ave,e_mag_t_as_ic_ave
    REAL(kind=8) :: e_mag_os_ave,e_mag_as_os_ave
    REAL(kind=8) :: Dip,DipCMB,e_cmb,elsAnel

    INTEGER :: lm,nR,nThetaB,nThetaStart
    INTEGER :: n_e_sets,n_spec

    CHARACTER(len=80) :: outFile
    INTEGER :: nOut,n_cmb_sets

    LOGICAL :: lGraphHeader

    REAL(kind=8) :: time
    REAL(kind=8) :: dt_norm

    INTEGER :: nBpotSets,nVpotSets,nTpotSets
    INTEGER :: lmStart,lmStop,lmStart_real, lmStop_real
#ifdef WITH_MPI
    !CHARACTER(len=MPI_MAX_ERROR_STRING) :: error_string
    !integer :: length_of_error
#endif
    !-- end of declaration
    !---------------------------------------------------------------


    !-- Initialise average for first time step:

    IF ( nAve.EQ.1 ) THEN  

       !zero=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
       IF ( n_graphs.GT.0 ) THEN
          IF ( l_conv ) THEN
             w_ave=zero
             z_ave=zero
          END IF
          IF ( l_heat ) THEN
             s_ave=zero
          END IF
          IF ( l_mag ) THEN
             b_ave=zero
             aj_ave=zero
             IF ( l_cond_ic ) THEN
                b_ic_ave=zero
                aj_ic_ave=zero
             END IF
          END IF
       END IF

    END IF  ! First step

    !-- Add new time step:

    IF ( l_conv ) THEN
       DO nR=1,n_r_max
          DO lm=llm,ulm
             w_ave(lm,nR)=w_ave(lm,nR) + time_passed*w(lm,nR)
             z_ave(lm,nR)=z_ave(lm,nR) + time_passed*z(lm,nR)
          END DO
       END DO
    END IF
    IF ( l_heat ) THEN
       DO nR=1,n_r_max
          DO lm=llm,ulm
             s_ave(lm,nR)=s_ave(lm,nR) + time_passed*s(lm,nR)
          END DO
       END DO
    END IF
    IF ( l_mag ) THEN
       DO nR=1,n_r_max
          DO lm=llm,ulm
             b_ave(lm,nR) =b_ave(lm,nR)  + time_passed*b(lm,nR)
             aj_ave(lm,nR)=aj_ave(lm,nR) + time_passed*aj(lm,nR)
          END DO
       END DO
       IF ( l_cond_ic ) THEN
          DO nR=1,n_r_ic_max
             DO lm=llm,ulm
                b_ic_ave(lm,nR) =b_ic_ave(lm,nR) + time_passed*b_ic(lm,nR)
                aj_ic_ave(lm,nR)=aj_ic_ave(lm,nR)+ time_passed*aj_ic(lm,nR)
             END DO
          END DO
       END IF
    END IF

    !--- Output, intermediate output every 10th averaging to save result
    !    will be overwritten.
    IF ( l_stop_time .OR. MOD(nAve,10).EQ.0 ) THEN

       !WRITE(*,"(A,2ES22.15)") "w_ave = ",get_global_sum( w_ave )
       time   =-1.D0  ! This signifies averaging in output files!
       dt_norm=1.D0/time_norm

       IF ( l_conv ) THEN
          DO nR=1,n_r_max
             DO lm=llm,ulm
                w_ave(lm,nR)=dt_norm*w_ave(lm,nR)
                z_ave(lm,nR)=dt_norm*z_ave(lm,nR)
             END DO
          END DO
       END IF
       IF ( l_heat ) THEN
          DO nR=1,n_r_max
             DO lm=llm,ulm
                s_ave(lm,nR)=dt_norm*s_ave(lm,nR)
             END DO
          END DO
       END IF
       IF ( l_mag ) THEN
          DO nR=1,n_r_max
             DO lm=llm,ulm
                b_ave(lm,nR) =dt_norm*b_ave(lm,nR)
                aj_ave(lm,nR)=dt_norm*aj_ave(lm,nR)
             END DO
          END DO
       END IF
       IF ( l_cond_ic ) THEN
          DO nR=1,n_r_ic_max
             DO lm=llm,ulm
                b_ic_ave(lm,nR) =dt_norm*b_ic_ave(lm,nR)
                aj_ic_ave(lm,nR)=dt_norm*aj_ic_ave(lm,nR)
             END DO
          END DO
       END IF

       !----- Get the radial derivatives:
       lmStart=lmStartB(rank+1)
       lmStop = lmStopB(rank+1)
       lmStart_real=2*lmStart-1
       lmStop_real =2*lmStop

       CALL get_drNS(w_ave,dw_ave,ulm_real-llm_real+1,lmStart_real-llm_real+1,lmStop_real-llm_real+1,&
            &        n_r_max,n_cheb_max,workA_LMloc,        &
            &        i_costf_init,d_costf_init,drx)
       CALL get_drNS(b_ave,db_ave,ulm_real-llm_real+1,lmStart_real-llm_real+1,lmStop_real-llm_real+1,&
            &        n_r_max,n_cheb_max,workA_LMloc,           &
            &        i_costf_init,d_costf_init,drx)
       IF ( l_heat ) THEN
          CALL get_drNS(s_ave,ds_ave,ulm_real-llm_real+1,lmStart_real-llm_real+1,lmStop_real-llm_real+1,&
               &        n_r_max,n_cheb_max,workA_LMloc,     &
               &        i_costf_init,d_costf_init,drx)
       END IF
       IF ( l_cond_ic ) THEN
          CALL get_ddrNS_even(b_ic_ave,db_ic_ave,ddb_ic_ave,        &
               &              ulm_real-llm_real+1,lmStart_real-llm_real+1,lmStop_real-llm_real+1,        &
               &              n_r_ic_max,n_cheb_ic_max,dr_fac_ic,workA_LMloc,        &
               &              i_costf1_ic_init,d_costf1_ic_init,        &
               &              i_costf2_ic_init,d_costf2_ic_init)
          CALL get_drNS_even(aj_ic_ave,dj_ic_ave,                   &
               &             ulm_real-llm_real+1,lmStart_real-llm_real+1,lmStop_real-llm_real+1,          &
               &             n_r_ic_max,n_cheb_ic_max,dr_fac_ic,workA_LMloc,          &
               &             i_costf1_ic_init,d_costf1_ic_init,          &
               &             i_costf2_ic_init,d_costf2_ic_init)
       END IF

       !----- Get averaged spectra:
       !      Note: average spectra will be in file no 0
       n_spec=0
       CALL spectrum(time,n_spec,w_ave,dw_ave,z_ave, &
            &        b_ave,db_ave,aj_ave,            &
            &        b_ic_ave,db_ic_ave,aj_ic_ave)  

       IF ( l_heat ) THEN
          CALL spectrumC(time,n_spec,s_ave,ds_ave)
       END IF
       IF ( l_save_out ) THEN
          OPEN(nLF,FILE=log_file,STATUS='UNKNOWN',                  &
               &             POSITION='APPEND')
       END IF

       !----- Write averaged energies into log-file at end of run:
       IF ( l_stop_time ) THEN 

          !----- Calculate energies of averaged field:
          n_e_sets=1
          CALL get_e_kin(time,.FALSE.,.TRUE.,n_e_sets, &
               &         w_ave,dw_ave,z_ave,           &
               &         e_kin_p_ave,e_kin_t_ave,      &
               &         e_kin_p_as_ave,e_kin_t_as_ave,&
               &         eKinR)

          CALL get_e_mag(time,.FALSE.,.TRUE.,n_e_sets,                  &
               &         b_ave,db_ave,aj_ave,                           &
               &         b_ic_ave,db_ic_ave,aj_ic_ave,                  &
               &         e_mag_p_ave,e_mag_t_ave,                       &
               &         e_mag_p_as_ave,e_mag_t_as_ave,                 &
               &         e_mag_p_ic_ave,e_mag_t_ic_ave,                 &
               &         e_mag_p_as_ic_ave,e_mag_t_as_ic_ave,           &
               &         e_mag_os_ave,e_mag_as_os_ave,e_cmb,Dip,DipCMB, &
               &         elsAnel)

          IF (rank.EQ.0) THEN
             !----- Output of energies of averaged field:
             WRITE(nLF,'(/,                                            &
                  &           '' ! ENERGIES OF TIME AVERAGED FIELD'')')
             WRITE(nLF,                                                &
                  &        '('' !  (total,poloidal,toroidal,total density)'')')
             WRITE(nLF,'(1P,'' !  Kinetic energies:'',4D16.6)')        &
                  &           (e_kin_p_ave+e_kin_t_ave),e_kin_p_ave,e_kin_t_ave,     &
                  &           (e_kin_p_ave+e_kin_t_ave)/vol_oc
             WRITE(nLF,'(1P,'' !  OC Mag  energies:'',4D16.6)')        &
                  &           (e_mag_p_ave+e_mag_t_ave),e_mag_p_ave,e_mag_t_ave,     &
                  &           (e_mag_p_ave+e_mag_t_ave)/vol_oc
             WRITE(nLF,'(1P,'' !  IC Mag  energies:'',4D16.6)')        &
                  &           (e_mag_p_ic_ave+e_mag_t_ic_ave),                       &
                  &           e_mag_p_ic_ave,e_mag_t_ic_ave,                         &
                  &           (e_mag_p_ic_ave+e_mag_t_ic_ave)/vol_ic
             WRITE(nLF,'(1P,'' !  OS Mag  energies:'',D16.6)')         &
                  &           e_mag_os_ave
             WRITE(nLF,'(/,'' !  AXISYMMETRIC PARTS:'')')
             WRITE(nLF,                                                &
                  &        '('' !  (total,poloidal,toroidal,total density)'')')
             WRITE(nLF,'(1P,'' !  Kinetic AS energies:'',4D16.6)')     &
                  &           (e_kin_p_as_ave+e_kin_t_as_ave),                       &
                  &           e_kin_p_as_ave,e_kin_t_as_ave,                         &
                  &           (e_kin_p_as_ave+e_kin_t_as_ave)/vol_oc
             WRITE(nLF,'(1P,'' !  OC Mag  AS energies:'',4D16.6)')     &
                  &           (e_mag_p_as_ave+e_mag_t_as_ave),                       &
                  &           e_mag_p_as_ave,e_mag_t_as_ave,                         &
                  &           (e_mag_p_as_ave+e_mag_t_as_ave)/vol_oc
             WRITE(nLF,'(1P,'' !  IC Mag  AS energies:'',4D16.6)')     &
                  &           (e_mag_p_as_ic_ave+e_mag_t_as_ic_ave),                 &
                  &           e_mag_p_as_ic_ave,e_mag_t_as_ic_ave,                   &
                  &           (e_mag_p_as_ic_ave+e_mag_t_as_ic_ave)/vol_ic
             WRITE(nLF,'(1P,'' !  OS Mag  AS energies:'',D16.6)')      &
                  &           e_mag_os_ave
             WRITE(nLF,'(1P,'' !  Relative ax. dip. E:'',D16.6)')      &
                  &           Dip           
          END IF
       END IF ! End of run ?
          
       !----- Construct name of graphic file and open it:
       ! For the graphic file of the average fields, we gather them
       ! on rank 0 and use the old serial output routine.

       IF (rank.EQ.0) THEN
          IF ( ngform.EQ.0 ) THEN
             graph_file='G_ave.'//tag
             OPEN(n_graph_file,FILE=graph_file,STATUS='UNKNOWN',FORM='UNFORMATTED')
          ELSE
             graph_file='g_ave.'//tag
             OPEN(n_graph_file,FILE=graph_file,STATUS='UNKNOWN',FORM='FORMATTED')
          END IF

          !----- Write header into graphic file:
          lGraphHeader=.TRUE.
          CALL graphOut(time,0,ngform,Vr,Vt,Vp,   &
               &        Br,Bt,Bp,Sr,              &
               &        0,sizeThetaB,lGraphHeader)
       END IF

       !----- Transform and output of data:
       CALL gather_all_from_lo_to_rank0(gt_OC,b_ave,b_ave_global)
       !----- Outer core:
       DO nR=1,n_r_max
          ! b_ave is different as it is again used later for graphOut_IC
          CALL gather_from_lo_to_rank0(db_ave(llm,nR),db_ave_global)
          CALL gather_from_lo_to_rank0(aj_ave(llm,nR),aj_ave_global)
          CALL gather_from_lo_to_rank0(w_ave(llm,nR),w_ave_global)
          CALL gather_from_lo_to_rank0(dw_ave(llm,nR),dw_ave_global)
          CALL gather_from_lo_to_rank0(z_ave(llm,nR),z_ave_global)
          CALL gather_from_lo_to_rank0(s_ave(llm,nR),s_ave_global)

          IF (rank.EQ.0) THEN
             CALL legPrep(b_ave_global(1,nR),db_ave_global,db_ave_global, &
                  &       aj_ave_global,aj_ave_global,dLh,lm_max,  &
                  &       l_max,minc,r(nR),.FALSE.,.TRUE.,       &
                  &       dLhb,bhG,bhC,dLhb,bhG,bhC)
             CALL legPrep(w_ave_global,dw_ave_global,dw_ave_global, &
                  &       z_ave_global,z_ave_global,dLh,lm_max,    &
                  &       l_max,minc,r(nR),.FALSE.,.TRUE.,       &
                  &       dLhw,vhG,vhC,dLhb,bhG,bhC)

             DO nThetaB=1,nThetaBs  
                nThetaStart=(nThetaB-1)*sizeThetaB+1

                !-------- Transform to grid space:
                CALL legTF(dLhb,bhG,bhC,dLhw,vhG,vhC,                  &
                     &     l_max,minc,nThetaStart,sizeThetaB,                  &
                     &     Plm,dPlm,lm_max,ncp,.TRUE.,.FALSE.,                  &
                     &     Br,Bt,Bp,Br,Br,Br)
                CALL legTF(dLhw,vhG,vhC,dLhw,vhG,vhC,                  &
                     &     l_max,minc,nThetaStart,sizeThetaB,                  &
                     &     Plm,dPlm,lm_max,ncp,.TRUE.,.FALSE.,                  &
                     &     Vr,Vt,Vp,Br,Br,Br)
                CALL legTF(s_ave_global,vhG,vhC,dLhw,vhG,vhC,           &
                     &     l_max,minc,nThetaStart,sizeThetaB,           &
                     &     Plm,dPlm,lm_max,ncp,.FALSE.,.FALSE.,           &
                     &     Sr,Vt,Vp,Br,Br,Br)
                CALL fft_thetab(Br,1)
                CALL fft_thetab(Bp,1)
                CALL fft_thetab(Bt,1)
                CALL fft_thetab(Vr,1)
                CALL fft_thetab(Vt,1)
                CALL fft_thetab(Vp,1)
                CALL fft_thetab(Sr,1)

                !-------- Graphic output:
                CALL graphOut(time,nR,ngform,Vr,Vt,Vp, &
                     &        Br,Bt,Bp,Sr,&
                     &        nThetaStart,sizeThetaB,lGraphHeader)
             END DO
          END IF
       END DO

       !----- Inner core: Transform is included in graphOut_IC!
       IF ( l_mag .AND. n_r_ic_max.GT.0 ) THEN
          CALL gather_all_from_lo_to_rank0(gt_IC,b_ic_ave,b_ic_ave_global)
          CALL gather_all_from_lo_to_rank0(gt_IC,db_ic_ave,db_ic_ave_global)
          CALL gather_all_from_lo_to_rank0(gt_IC,ddb_ic_ave,ddb_ic_ave_global)
          CALL gather_all_from_lo_to_rank0(gt_IC,aj_ic_ave,aj_ic_ave_global)
          CALL gather_all_from_lo_to_rank0(gt_IC,dj_ic_ave,dj_ic_ave_global)

          IF (rank.EQ.0) THEN
             CALL graphOut_IC(ngform,b_ic_ave_global,db_ic_ave_global,&
                  &           ddb_ic_ave_global,aj_ic_ave_global,&
                  &           dj_ic_ave_global,b_ave_global)
          END IF
       END IF

       if (rank.eq.0) CLOSE(n_graph_file)  ! close graphic output file !

       !----- Write info about graph-file into STDOUT and log-file:
       IF ( l_stop_time ) THEN
          IF (rank.EQ.0) WRITE(nLF,'(/,'' ! WRITING AVERAGED GRAPHIC FILE !'')')
       END IF

       !--- Store time averaged poloidal magnetic coeffs at cmb
       IF (rank.EQ.0) THEN
          IF ( l_mag) THEN
             outFile='B_coeff_cmb_ave.'//tag
             nOut   =93
             n_cmb_sets=0
             !CALL write_Bcmb(time,b(1,n_r_cmb),lm_max,l_max,           &
             !     &           l_max_cmb,minc,lm2,n_cmb_sets,outFile,nOut)
             CALL write_Bcmb(time,b_ave_global(1,n_r_cmb),1,lm_max,l_max,           &
                  &          l_max_cmb,minc,lm2,n_cmb_sets,outFile,nOut)
          ENDIF
       END IF

       !--- Store potentials of averaged field:
       !    dw_ave and db_ave used as work arrays here.
       nBpotSets=-1
       nVpotSets=-1
       nTpotSets=-1
       IF ( l_mag) THEN
          CALL storePotW(time,b_ave,aj_ave,b_ic_ave,aj_ic_ave,      &
               &         workA_LMloc,dw_ave,db_ave,nBpotSets,'Bpot_ave.', &
               &         omega_ma,omega_ic)
       ENDIF
       CALL storePotW(time,w_ave,z_ave,b_ic_ave,aj_ic_ave,          &
            &         workA_LMloc,dw_ave,db_ave,nVpotSets,'Vpot_ave.', &
            &                                      omega_ma,omega_ic)
       CALL storePotW(time,s_ave,z_ave,b_ic_ave,aj_ic_ave,          &
            &         workA_LMloc,dw_ave,db_ave,nVpotSets,'Tpot_ave.', &
            &                                      omega_ma,omega_ic)

       IF ( l_save_out ) CLOSE(nLF)

       ! now correct the stored average fields by the factor which has been
       ! applied before
       IF ( l_conv ) THEN
          DO nR=1,n_r_max
             DO lm=llm,ulm
                w_ave(lm,nR)=w_ave(lm,nR)*time_norm
                z_ave(lm,nR)=z_ave(lm,nR)*time_norm
             END DO
          END DO
       END IF
       IF ( l_heat ) THEN
          DO nR=1,n_r_max
             DO lm=llm,ulm
                s_ave(lm,nR)=s_ave(lm,nR)*time_norm
             END DO
          END DO
       END IF
       IF ( l_mag ) THEN
          DO nR=1,n_r_max
             DO lm=llm,ulm
                b_ave(lm,nR) =b_ave(lm,nR)*time_norm
                aj_ave(lm,nR)=aj_ave(lm,nR)*time_norm
             END DO
          END DO
       END IF
       IF ( l_cond_ic ) THEN
          DO nR=1,n_r_ic_max
             DO lm=llm,ulm
                b_ic_ave(lm,nR) =b_ic_ave(lm,nR)*time_norm
                aj_ic_ave(lm,nR)=aj_ic_ave(lm,nR)*time_norm
             END DO
          END DO
       END IF


    END IF ! last time step ?


    RETURN
  END SUBROUTINE fields_average

  !------------------------------------------------------------------------------
END MODULE fields_average_mod
