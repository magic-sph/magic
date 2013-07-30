!$Id$
MODULE fields_average_mod
  USE truncation
  USE radial_functions
  USE num_param
  USE blocking
  USE horizontal_data
  USE logic
  use kinetic_energy
  use magnetic_energy
  USE output_data
#if (FFTLIB==JW)
  USE fft_JW
#elif (FFTLIB==MKL)
  USE fft_MKL
#endif
  USE const

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

    ALLOCATE( w_ave(lm_max_ave,n_r_max_ave) )
    ALLOCATE( z_ave(lm_max_ave,n_r_max_ave) )
    ALLOCATE( s_ave(lm_max_ave,n_r_max_ave) )
    ALLOCATE( b_ave(lm_max_ave,n_r_max_ave) )
    ALLOCATE( aj_ave(lm_max_ave,n_r_max_ave) )
    ALLOCATE( b_ic_ave(lm_max_ave,n_r_ic_max_ave) )
    ALLOCATE( aj_ic_ave(lm_max_ave,n_r_ic_max_ave) )

    END SUBROUTINE initialize_fields_average_mod
  !**********************************************************************
    SUBROUTINE fields_average(nAve,l_stop_time,                     &
       &     time_passed,time_norm,omega_ic,omega_ma,                     &
       &                       w,z,s,b,aj,b_ic,aj_ic)
    !***********************************************************************

    !   !------------ This is release 2 level 2  --------------!
    !   !------------ Created on 2/6/02  by JW. --------------!


    !  +-------------+----------------+------------------------------------+
    !  |                                                                   |
    !  |  This subroutine averages fields b and v over time.               |
    !  |                                                                   |
    !  +-------------------------------------------------------------------+
    !  |  ruler                                                            |
    !  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
    !--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+


    !-- Input of variables:
    INTEGER :: nAve    ! number for averaged time steps
    LOGICAL :: l_stop_time     ! true if this is the last time step
    REAL(kind=8) :: time_passed     ! time step
    REAL(kind=8) :: time_norm
    REAL(kind=8) :: omega_ic,omega_ma

    !-- Input of scalar fields:
    COMPLEX(kind=8) :: w(lm_max,n_r_max)
    COMPLEX(kind=8) :: z(lm_max,n_r_max)
    COMPLEX(kind=8) :: s(lm_max,n_r_max)
    COMPLEX(kind=8) :: b(lm_maxMag,n_r_maxMag)
    COMPLEX(kind=8) :: aj(lm_maxMag,n_r_maxMag)
    COMPLEX(kind=8) :: b_ic(lm_maxMag,n_r_ic_maxMag)
    COMPLEX(kind=8) :: aj_ic(lm_maxMag,n_r_ic_maxMag)


    !-- Local stuff:

    !----- Time averaged fields:
    COMPLEX(kind=8) :: dw_ave(lm_max_ave,n_r_max_ave)
    COMPLEX(kind=8) :: ds_ave(lm_max_ave,n_r_max_ave)
    COMPLEX(kind=8) :: db_ave(lm_max_ave,n_r_max_ave)
    COMPLEX(kind=8) :: db_ic_ave(lm_max_ave,n_r_ic_max_ave)
    COMPLEX(kind=8) :: ddb_ic_ave(lm_max_ave,n_r_ic_max_ave)
    COMPLEX(kind=8) :: dj_ic_ave(lm_max_ave,n_r_ic_max_ave)

    !----- Work array:
    COMPLEX(kind=8) :: workA(lm_max,n_r_max)

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

    !-- end of declaration
    !---------------------------------------------------------------


    !-- Initialise average for first time step:

    IF ( nAve.EQ.1 ) THEN  

       !zero=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
       IF ( n_graphs.GT.0 ) THEN
          IF ( l_conv ) THEN
             DO nR=1,n_r_max
                DO lm=1,lm_max
                   w_ave(lm,nR)=zero
                   z_ave(lm,nR)=zero
                END DO
             END DO
          END IF
          IF ( l_heat ) THEN
             DO nR=1,n_r_max
                DO lm=1,lm_max
                   s_ave(lm,nR)=zero
                END DO
             END DO
          END IF
          IF ( l_mag ) THEN
             DO nR=1,n_r_max
                DO lm=1,lm_max
                   b_ave(lm,nR) =zero
                   aj_ave(lm,nR)=zero
                END DO
             END DO
             IF ( l_cond_ic ) THEN
                DO nR=1,n_r_ic_max
                   DO lm=1,lm_max
                      b_ic_ave(lm,nR) =zero
                      aj_ic_ave(lm,nR)=zero
                   END DO
                END DO
             END IF
          END IF
       END IF

    END IF  ! First step

    !-- Add new time step:

    IF ( l_conv ) THEN
       DO nR=1,n_r_max
          DO lm=1,lm_max
             w_ave(lm,nR)=w_ave(lm,nR) + time_passed*w(lm,nR)
             z_ave(lm,nR)=z_ave(lm,nR) + time_passed*z(lm,nR)
          END DO
       END DO
    END IF
    IF ( l_heat ) THEN
       DO nR=1,n_r_max
          DO lm=1,lm_max
             s_ave(lm,nR)=s_ave(lm,nR) + time_passed*s(lm,nR)
          END DO
       END DO
    END IF
    IF ( l_mag ) THEN
       DO nR=1,n_r_max
          DO lm=1,lm_max
             b_ave(lm,nR) =b_ave(lm,nR)  + time_passed*b(lm,nR)
             aj_ave(lm,nR)=aj_ave(lm,nR) + time_passed*aj(lm,nR)
          END DO
       END DO
       IF ( l_cond_ic ) THEN
          DO nR=1,n_r_ic_max
             DO lm=1,lm_max
                b_ic_ave(lm,nR) =b_ic_ave(lm,nR) +                     &
                     &                               time_passed*b_ic(lm,nR)
                aj_ic_ave(lm,nR)=aj_ic_ave(lm,nR)+                     &
                     &                               time_passed*aj_ic(lm,nR)
             END DO
          END DO
       END IF
    END IF


    !--- Output, intermediate output every 10th averaging to save result
    !    will be overwritten.
    IF ( l_stop_time .OR. MOD(nAve,10).EQ.0 ) THEN

       time   =-1.D0  ! This signifies averaging in output files!
       dt_norm=1.D0/time_norm

       IF ( l_conv ) THEN
          DO nR=1,n_r_max
             DO lm=1,lm_max
                w_ave(lm,nR)=dt_norm*w_ave(lm,nR)
                z_ave(lm,nR)=dt_norm*z_ave(lm,nR)
             END DO
          END DO
       END IF
       IF ( l_heat ) THEN
          DO nR=1,n_r_max
             DO lm=1,lm_max
                s_ave(lm,nR)=dt_norm*s_ave(lm,nR)
             END DO
          END DO
       END IF
       IF ( l_mag ) THEN
          DO nR=1,n_r_max
             DO lm=1,lm_max
                b_ave(lm,nR) =dt_norm*b_ave(lm,nR)
                aj_ave(lm,nR)=dt_norm*aj_ave(lm,nR)
             END DO
          END DO
       END IF
       IF ( l_cond_ic ) THEN
          DO nR=1,n_r_ic_max
             DO lm=1,lm_max
                b_ic_ave(lm,nR) =dt_norm*b_ic_ave(lm,nR)
                aj_ic_ave(lm,nR)=dt_norm*aj_ic_ave(lm,nR)
             END DO
          END DO
       END IF

       !----- Get the radial derivatives:
       CALL get_drNS(w_ave,dw_ave,lm_max_real,1,lm_max_real,        &
            &                      n_r_max,n_cheb_max,workA,        &
            &                 i_costf_init,d_costf_init,drx)
       CALL get_drNS(b_ave,db_ave,lm_max_real,1,2*lm_max,           &
            &                   n_r_max,n_cheb_max,workA,           &
            &              i_costf_init,d_costf_init,drx)
       IF ( l_heat ) THEN
          CALL get_drNS(s_ave,ds_ave,lm_max_real,1,lm_max_real,     &
               &                      n_r_max,n_cheb_max,workA,     &
               &                 i_costf_init,d_costf_init,drx)
       END IF
       IF ( l_cond_ic ) THEN
          CALL get_ddrNS_even(b_ic_ave,db_ic_ave,ddb_ic_ave,        &
               &                                lm_max_real,1,lm_max_real,        &
               &                 n_r_ic_max,n_cheb_ic_max,dr_fac_ic,workA,        &
               &                        i_costf1_ic_init,d_costf1_ic_init,        &
               &                        i_costf2_ic_init,d_costf2_ic_init)
          CALL get_drNS_even(aj_ic_ave,dj_ic_ave,                   &
               &                              lm_max_real,1,lm_max_real,          &
               &               n_r_ic_max,n_cheb_ic_max,dr_fac_ic,workA,          &
               &                      i_costf1_ic_init,d_costf1_ic_init,          &
               &                      i_costf2_ic_init,d_costf2_ic_init)
       END IF

       !----- Get averaged spectra:
       !      Note: average spectra will be in file no 0
       n_spec=0
       CALL spectrum(time,n_spec,w_ave,dw_ave,z_ave,                &
            &                              b_ave,db_ave,aj_ave,                &
            &                     b_ic_ave,db_ic_ave,aj_ic_ave)  

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
          CALL get_e_kin(time,.FALSE.,.TRUE.,n_e_sets,              &
               &                                 w_ave,dw_ave,z_ave, &
               &                            e_kin_p_ave,e_kin_t_ave, &
               &                      e_kin_p_as_ave,e_kin_t_as_ave, &
               &                                              eKinR)

          CALL get_e_mag(time,.FALSE.,.TRUE.,n_e_sets,              &
               &                                b_ave,db_ave,aj_ave,              &
               &                       b_ic_ave,db_ic_ave,aj_ic_ave,              &
               &                            e_mag_p_ave,e_mag_t_ave,              &
               &                      e_mag_p_as_ave,e_mag_t_as_ave,              &
               &                      e_mag_p_ic_ave,e_mag_t_ic_ave,              &
               &                e_mag_p_as_ic_ave,e_mag_t_as_ic_ave,              &
               &      e_mag_os_ave,e_mag_as_os_ave,e_cmb,Dip,DipCMB,              &
               &                                            elsAnel)

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

       END IF ! End of run ?

       !----- Construct name of graphic file and open it:
       IF ( ngform.EQ.0 ) THEN
          graph_file='G_ave.'//tag
          OPEN(n_graph_file,FILE=graph_file,STATUS='UNKNOWN',FORM='UNFORMATTED')
       ELSE
          graph_file='g_ave.'//tag
          OPEN(n_graph_file,FILE=graph_file,STATUS='UNKNOWN',FORM='FORMATTED')
       END IF

       !----- Write header into graphic file:
       lGraphHeader=.TRUE.
       CALL graphOut(time,0,ngform,Vr,Vt,Vp,                        &
            &                              Br,Bt,Bp,Sr,                        &
            &                0,sizeThetaB,lGraphHeader)

       !----- Transform and output of data:

       !----- Outer core:

       !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1)                                   &
       !$OMP  PRIVATE(nThetaB,nThetaStart,Vr,Vt,Vp,Br,Bt,Bp,Sr,                &
       !$OMP                         dLhb,bhG,bhC,dLhw,vhG,vhC)

       DO nR=1,n_r_max

          CALL legPrep(b_ave(1,nR),db_ave(1,nR),db_ave(1,nR),       &
               &                      aj_ave(1,nR),aj_ave(1,nR),dLh,lm_max,       &
               &                           l_max,minc,r(nR),.FALSE.,.TRUE.,       &
               &                                 dLhb,bhG,bhC,dLhb,bhG,bhC)
          CALL legPrep(w_ave(1,nR),dw_ave(1,nR),dw_ave(1,nR),       &
               &                        z_ave(1,nR),z_ave(1,nR),dLh,lm_max,       &
               &                           l_max,minc,r(nR),.FALSE.,.TRUE.,       &
               &                                 dLhw,vhG,vhC,dLhb,bhG,bhC)

          DO nThetaB=1,nThetaBs  
             nThetaStart=(nThetaB-1)*sizeThetaB+1

             !-------- Transform to grid space:
             CALL legTF(dLhb,bhG,bhC,dLhw,vhG,vhC,                  &
                  &              l_max,minc,nThetaStart,sizeThetaB,                  &
                  &             Plm,dPlm,lm_max,ncp,.TRUE.,.FALSE.,                  &
                  &                              Br,Bt,Bp,Br,Br,Br)
             CALL legTF(dLhw,vhG,vhC,dLhw,vhG,vhC,                  &
                  &              l_max,minc,nThetaStart,sizeThetaB,                  &
                  &             Plm,dPlm,lm_max,ncp,.TRUE.,.FALSE.,                  &
                  &                              Vr,Vt,Vp,Br,Br,Br)
             CALL legTF(s_ave(1,nR),vhG,vhC,dLhw,vhG,vhC,           &
                  &                     l_max,minc,nThetaStart,sizeThetaB,           &
                  &                   Plm,dPlm,lm_max,ncp,.FALSE.,.FALSE.,           &
                  &                                     Sr,Vt,Vp,Br,Br,Br)
             CALL fft_thetab(Br,1)
             CALL fft_thetab(Bp,1)
             CALL fft_thetab(Bt,1)
             CALL fft_thetab(Vr,1)
             CALL fft_thetab(Vt,1)
             CALL fft_thetab(Vp,1)
             CALL fft_thetab(Sr,1)

             !-------- Grafik output:
             CALL graphOut(time,nR,ngform,Vr,Vt,Vp,                    &
                  &                                  Br,Bt,Bp,Sr,                    &
                  &          nThetaStart,sizeThetaB,lGraphHeader)

          END DO
       END DO

       !$OMP END PARALLEL DO   !

       !----- Inner core: Transform is included in graphOut_IC!
       IF ( l_mag .AND. n_r_ic_max.GT.0 )                           &
            &     CALL graphOut_IC(ngform,b_ic_ave,db_ic_ave,ddb_ic_ave,       &
            &                                 aj_ic_ave,dj_ic_ave,b_ave)

       CLOSE(n_graph_file)


       !----- Write info about graph-file into STDOUT and log-file:
       IF ( l_stop_time ) THEN
          WRITE(*,'(/,'' ! WRITING AVERAGED GRAPHIC FILE !'')')
          WRITE(nLF,'(/,'' ! WRITING AVERAGED GRAPHIC FILE !'')')
       END IF

       !--- Store time averaged poloidal magnetic coeffs at cmb
       IF ( l_mag) THEN
          outFile='B_coeff_cmb_ave.'//tag
          nOut   =93
          n_cmb_sets=0
          CALL write_Bcmb(time,b_ave(1,n_r_cmb),lm_max,l_max,           &
               &           l_max_cmb,minc,lm2,n_cmb_sets,outFile,nOut)
       ENDIF


       !--- Store potentials of averaged field:
       !    dw_ave and db_ave used as work arrays here.
       nBpotSets=-1
       nVpotSets=-1
       nTpotSets=-1
       IF ( l_mag) THEN
          CALL storePotW(time,b_ave,aj_ave,b_ic_ave,aj_ic_ave,      &
               &                  workA,dw_ave,db_ave,nBpotSets,'Bpot_ave.',      &
               &                                          omega_ma,omega_ic)
       ENDIF
       CALL storePotW(time,w_ave,z_ave,b_ic_ave,aj_ic_ave,          &
            &              workA,dw_ave,db_ave,nVpotSets,'Vpot_ave.',          &
            &                                      omega_ma,omega_ic)
       CALL storePotW(time,s_ave,z_ave,b_ic_ave,aj_ic_ave,          &
            &              workA,dw_ave,db_ave,nVpotSets,'Tpot_ave.',          &
            &                                      omega_ma,omega_ic)

       IF ( l_save_out ) CLOSE(nLF)

       IF ( l_conv ) THEN
          DO nR=1,n_r_max
             DO lm=1,lm_max
                w_ave(lm,nR)=w_ave(lm,nR)/dt_norm
                z_ave(lm,nR)=z_ave(lm,nR)/dt_norm
             END DO
          END DO
       END IF
       IF ( l_heat ) THEN
          DO nR=1,n_r_max
             DO lm=1,lm_max
                s_ave(lm,nR)=s_ave(lm,nR)/dt_norm
             END DO
          END DO
       END IF
       IF ( l_mag ) THEN
          DO nR=1,n_r_max
             DO lm=1,lm_max
                b_ave(lm,nR) =b_ave(lm,nR)/dt_norm
                aj_ave(lm,nR)=aj_ave(lm,nR)/dt_norm
             END DO
          END DO
       END IF
       IF ( l_cond_ic ) THEN
          DO nR=1,n_r_ic_max
             DO lm=1,lm_max
                b_ic_ave(lm,nR) =b_ic_ave(lm,nR)/dt_norm
                aj_ic_ave(lm,nR)=aj_ic_ave(lm,nR)/dt_norm
             END DO
          END DO
       END IF


    END IF ! last time step ?


    RETURN
  END SUBROUTINE fields_average

  !------------------------------------------------------------------------------
END MODULE fields_average_mod
