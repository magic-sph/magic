!$Id$
MODULE outPar_mod
  USE truncation, ONLY:n_r_max, l_max, lm_max
  USE blocking, ONLY: nfs,nThetaBs,sizeThetaB
  USE logic,ONLY: l_viscBcCalc,l_anel
  USE horizontal_data, ONLY: gauss
  USE fields, ONLY: s_Rloc
  USE physical_parameters, ONLY: ek,prmag
  USE const, ONLY: pi,mass
  USE radial_functions, ONLY:r,or2,sigma,rho0
  USE radial_data, ONLY: n_r_icb,nRstart,nRstop
  use parallel_mod
  USE output_data,only: tag
  IMPLICIT NONE

  REAL(kind=8),ALLOCATABLE :: dlVMeanR(:),dlVcMeanR(:)
  REAL(kind=8),ALLOCATABLE :: dlVu2MeanR(:),dlVu2cMeanR(:)
  REAL(kind=8),ALLOCATABLE :: RolMeanR(:),RolMeanRu2(:),RmMeanR(:)
  REAL(kind=8),ALLOCATABLE :: sMeanR(:),Svar(:),Mvar(:)
  REAL(kind=8),ALLOCATABLE :: uhMeanR(:),duhMeanR(:)
  REAL(kind=8),ALLOCATABLE :: gradT2MeanR(:)


contains
  SUBROUTINE initialize_outPar_mod

    ALLOCATE( dlVMeanR(n_r_max),dlVcMeanR(n_r_max) )
    ALLOCATE( sMeanR(n_r_max),Svar(nRstart:nRstop),Mvar(nRstart:nRstop) )
    ALLOCATE( dlVu2MeanR(n_r_max),dlVu2cMeanR(n_r_max) )
    ALLOCATE( RolMeanR(n_r_max),RolMeanRu2(n_r_max),RmMeanR(n_r_max) )
    ALLOCATE( uhMeanR(n_r_max),duhMeanR(n_r_max),gradT2MeanR(n_r_max) )

    dlVMeanR     =0.D0
    dlVcMeanR    =0.D0
    dlVu2MeanR   =0.D0
    dlVu2cMeanR  =0.D0
    RolMeanR     =0.d0
    RolMeanRu2   =0.d0
    RmMeanR      =0.d0
    sMeanR       =0.d0
    uhMeanR      =0.d0
    duhMeanR     =0.d0
    gradT2MeanR  =0.d0


  END SUBROUTINE initialize_outPar_mod

  !***********************************************************************
  SUBROUTINE outPar(timePassed,timeNorm,n_time_step,l_stop_time, &
                    ekinR,RolRu2,dlVR,dlVRc,dlVRu2,dlVRu2c, &
                    uhLMr,duhLMr,gradsLMr,RmR)
    !***********************************************************************

    !--- Input of variables
    IMPLICIT NONE

    REAL(kind=8),INTENT(IN) :: timePassed,timeNorm
    LOGICAL,INTENT(IN) :: l_stop_time
    INTEGER,INTENT(IN) :: n_time_step
    REAL(kind=8),INTENT(IN) :: RolRu2(n_r_max),dlVRu2(n_r_max),dlVRu2c(n_r_max)
    REAL(kind=8),INTENT(IN) :: dlVR(n_r_max),dlVRc(n_r_max)
    REAL(kind=8),INTENT(IN) :: ekinR(n_r_max)     ! kinetic energy w radius
    REAL(kind=8),INTENT(IN) :: uhLMr(l_max+1,nRstart:nRstop)
    REAL(kind=8),INTENT(IN) :: duhLMr(l_max+1,nRstart:nRstop)
    REAL(kind=8),INTENT(IN) :: gradsLMr(l_max+1,nRstart:nRstop)
    REAL(kind=8),INTENT(OUT):: RmR(n_r_max)

    !-- Further counter
    INTEGER :: nR,lm,n

    !--- Property parameters:
    REAL(kind=8) :: Mtmp
    REAL(kind=8),DIMENSION(n_r_max) :: ReR, RoR, RolR

    CHARACTER(len=76) :: filename

    ! For horizontal velocity
    INTEGER :: nTheta,nThetaStart,nThetaBlock,nThetaNHS
    REAL(kind=8),DIMENSION(nRstart:nRstop) :: duhR,uhR,gradT2R,sR
    REAL(kind=8),DIMENSION(n_r_max) :: duhR_global,uhR_global,gradT2R_global,sR_global
    REAL(kind=8),DIMENSION(n_r_max) :: Svar_global
    REAL(kind=8),DIMENSION(nfs) :: duh,uh,gradT2

    INTEGER :: i,sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1)


    !--- end of declaration
    !-----------------------------------------------------------------
    IF ( l_viscBcCalc ) THEN
       !IF (rank.EQ.0) THEN
       !DO nR=1,n_r_max
       DO nR=nRstart,nRstop
          !sR(nR)=0.D0
          ! Mean entropy profile
          sR(nR) = SUM(REAL(s_Rloc(:,nR)))
          !DO lm=1,lm_max
          !   sR(nR) = sR(nR)+REAL(s_Rloc(lm,nR))
          !END DO
          ! calculate entropy/temperature variance:
          IF (n_time_step .LE. 1) THEN
             Mvar(nR)       =sR(nR)
             Svar(nR)       =0.d0
          ELSE
             Mtmp      =Mvar(nR)
             Mvar(nR)  =Mvar(nR) + (sR(nR)-Mvar(nR))/n_time_step
             Svar(nR)  =Svar(nR) + (sR(nR)-Mtmp)*(sR(nR)-Mvar(nR))
          END IF
          !WRITE(*,"(A,I3,A,3ES20.12)") "sR,Svar,Mvar (",nR,") = ",sR(nR),Svar(nR),Mvar(nR)
       END DO
       !END IF

        DO nR=nRstart,nRstop
            uhR(nR) =0.d0
            gradT2R(nR)=0.d0
            duhR(nR)=0.d0
            DO n=1,nThetaBs ! Loop over theta blocks
                nTheta=(n-1)*sizeThetaB
                nThetaStart=nTheta+1
                CALL lmAS2pt(duhLMr(1,nR),duh,nThetaStart,sizeThetaB)
                CALL lmAS2pt(uhLMr(1,nR),uh,nThetaStart,sizeThetaB)
                CALL lmAS2pt(gradsLMr(1,nR),gradT2,nThetaStart,sizeThetaB)
                DO nThetaBlock=1,sizeThetaB
                    nTheta=nTheta+1
                    nThetaNHS=(nTheta+1)/2
                    duhR(nR)=duhR(nR)+gauss(nThetaNHS)*duh(nThetaBlock)
                    uhR(nR) =uhR(nR) +gauss(nThetaNHS)* uh(nThetaBlock)
                    gradT2R(nR)=gradT2R(nR)+gauss(nThetaNHS)*gradT2(nThetaBlock)
                END DO
            END DO
        END DO
        duhR=0.5d0*duhR ! Normalisation for the theta integration
        uhR =0.5d0* uhR ! Normalisation for the theta integration
        gradT2R =0.5d0*gradT2R ! Normalisation for the theta integration

        sendcount  = (nRstop-nRstart+1)
        recvcounts = nr_per_rank
        recvcounts(n_procs-1) = (nr_per_rank+1)
        DO i=0,n_procs-1
            displs(i) = i*nr_per_rank
        END DO
        CALL MPI_GatherV(duhR,sendcount,MPI_DOUBLE_PRECISION,&
            &           duhR_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
            &           0,MPI_COMM_WORLD,ierr)
        CALL MPI_GatherV(uhR,sendcount,MPI_DOUBLE_PRECISION,&
            &           uhR_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
            &           0,MPI_COMM_WORLD,ierr)
        CALL MPI_GatherV(gradT2R,sendcount,MPI_DOUBLE_PRECISION,&
            &           gradT2R_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
            &           0,MPI_COMM_WORLD,ierr)
        CALL MPI_GatherV(sR,sendcount,MPI_DOUBLE_PRECISION,&
            &           sR_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
            &           0,MPI_COMM_WORLD,ierr)
        CALL MPI_GatherV(Svar,sendcount,MPI_DOUBLE_PRECISION,&
            &           Svar_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
            &           0,MPI_COMM_WORLD,ierr)
     END IF

    IF (rank.EQ.0) THEN
       DO nR=1,n_r_max
          ReR(nR)=SQRT(2.D0*ekinR(nR)*or2(nR)/(4*pi*mass))
          RoR(nR)=ReR(nR)*ek
          IF ( dlVR(nR) /= 0d0 ) THEN
              RolR(nR)=RoR(nR)/dlVR(nR)
          ELSE
              RolR(nR)=RoR(nR)
          END IF
          !WRITE(*,"(A,I3,4ES20.12)") "outPar: ",nR,dlVR(nR),dlVRc(nR),dlVRu2(nR),dlVRu2c(nR)
          RmR(nR)=ReR(nR)*prmag*sigma(nR)*r(nR)*r(nR)
       END DO

       dlVMeanR   =dlVMeanR   +timePassed*dlVR
       dlVcMeanR  =dlVcMeanR  +timePassed*dlVRc
       dlVu2MeanR =dlVu2MeanR +timePassed*dlVRu2
       dlVu2cMeanR=dlVu2cMeanR+timePassed*dlVRu2c
       RolMeanR   =RolMeanR   +timePassed*RolR
       RolMeanRu2 =RolMeanRu2 +timePassed*RolRu2
       RmMeanR    =RmMeanR    +timePassed*RmR*dsqrt(mass/rho0)*or2
       !WRITE(*,"(A,ES20.12)") "dlVcMeanR(n_r_icb) = ",dlVcMeanR(n_r_icb)
       ! this is to get u2 value for RmR(r) to plot in parrad.tag
       ! and also remove r**2, so it has to be volume-averaged 
       ! like RolR
       IF ( l_viscBcCalc ) THEN
          sMeanR     =sMeanR    +timePassed*sR_global
          uhMeanR    =uhMeanR   +timePassed*uhR_global
          !WRITE(*,"(A,3ES20.12)") "duhR_global,timePassed,duhMeanR = ",SUM(duhR_global),timePassed,SUM(duhMeanR)
          duhMeanR   =duhMeanR  +timePassed*duhR_global
          gradT2MeanR=gradT2MeanR+timePassed*gradT2R_global
       END IF

       IF ( l_stop_time ) THEN 

          dlVMeanR   =dlVMeanR/timeNorm
          dlVcMeanR  =dlVcMeanR/timeNorm
          RolMeanR   =RolMeanR/timeNorm
          IF ( l_anel ) THEN
             dlVu2MeanR =dlVu2MeanR/timeNorm
             dlVu2cMeanR=dlVu2cMeanR/timeNorm
             RolMeanRu2 =RolMeanRu2/timeNorm
          ELSE
             dlVu2MeanR =dlVMeanR
             dlVu2cMeanR=dlVcMeanR
             RolMeanRu2 =RolMeanR
          END IF
          RmMeanR    =RmMeanR/timeNorm

          IF ( l_viscBcCalc ) THEN
             sMeanR     =sMeanR/timeNorm
             Svar_global=Svar_global/(n_time_step-1)
             duhMeanR   =duhMeanR/timeNorm
             uhMeanR    =uhMeanR/timeNorm
             gradT2MeanR=gradT2MeanR/timeNorm
          END IF

          !----- Output into parrad file:
          filename='parR.'//tag
          OPEN(99,FILE=filename,STATUS='UNKNOWN')
          DO nR=1,n_r_max
             WRITE(99,'(D20.10,8D12.4)')       &
                        &   r(nR),             &! 1) radius
                        &   RmMeanR(nR),       &! 2) magnetic Reynolds number
                        &   RolMeanR(nR),      &! 3) local Rossby number
                        &   RolMeanRu2(nR),    &! 4) u squared local Rossby number
                        &   dlVMeanR(nR),      &! 5) local length scale
                        &   dlVcMeanR(nR),     &! 6) conv. local length scale
                        &   dlVu2MeanR(nR),    &! 7) u squared local length scale 
                        &   dlVu2cMeanR(nR)     ! 8) u squared conv. local length scale
          END DO
          CLOSE(99)

          IF ( l_viscBcCalc ) THEN
             filename='bLayersR.'//tag
             OPEN(99,FILE=filename,STATUS='UNKNOWN')
             DO nR=1,n_r_max
                WRITE(99,'(D20.10,6D12.4)')            &
                        &   r(nR),                     &! 1) radius
                        &   sMeanR(nR)/SQRT(4.D0*pi),  &! 2) entropy
                        &   Svar_global(nR),           &! 3) entropy variance
                        &   uhMeanR(nR),               &! 4) uh
                        &   duhMeanR(nR),              &! 5) duh/dr
                        &   gradT2MeanR(nR)             ! 6) (grad T)**2
             END DO
             CLOSE(99)
          END IF

       END IF ! l_stop_time ?

    END IF ! rank0

    RETURN 
  END SUBROUTINE outPar
  !-----------------------------------------------------------------------
END MODULE outPar_mod
