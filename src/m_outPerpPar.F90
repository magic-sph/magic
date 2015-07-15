!$Id$
MODULE outPerpPar_mod
  USE truncation, ONLY: n_r_max,l_max
  USE blocking, ONLY: nThetaBs,sizeThetaB,nfs
  USE horizontal_data, ONLY: gauss
  USE const, ONLY: pi
  USE output_data, ONLY: tag, n_perpPar_file, perpPar_file
  USE radial_data, ONLY: nRstart,nRstop
  USE radial_functions, ONLY: r,dr_fac,i_costf_init,d_costf_init
  USE parallel_mod
  USE logic, ONLY: l_save_out
  USE num_param, ONLY: tScale
 
  IMPLICIT NONE

  REAL(kind=8),ALLOCATABLE :: EperpMeanR(:),EparMeanR(:)
  REAL(kind=8),ALLOCATABLE :: EperpaxiMeanR(:),EparaxiMeanR(:)


CONTAINS
 !********************************************
 SUBROUTINE initialize_outPerpPar_mod

    ALLOCATE( EperpMeanR(n_r_max),EparMeanR(n_r_max) )
    ALLOCATE( EperpaxiMeanR(n_r_max),EparaxiMeanR(n_r_max) )

    EperpMeanR   =0.D0
    EparMeanR    =0.D0
    EperpaxiMeanR=0.D0
    EparaxiMeanR =0.D0

 END SUBROUTINE initialize_outPerpPar_mod
 !********************************************
 SUBROUTINE outPerpPar(time,timePassed,timeNorm,l_stop_time, &
                 &     EperpLMr,EparLMr,EperpaxiLMr,EparaxiLMr)

    USE integration, ONLY: rInt

    IMPLICIT NONE

    !--- Input of variables
    REAL(kind=8),INTENT(IN) :: time,timePassed,timeNorm
    LOGICAL,INTENT(IN) :: l_stop_time
    REAL(kind=8),INTENT(IN) :: EparLMr(l_max+1,nRstart:nRstop)
    REAL(kind=8),INTENT(IN) :: EperpLMr(l_max+1,nRstart:nRstop)
    REAL(kind=8),INTENT(IN) :: EparaxiLMr(l_max+1,nRstart:nRstop)
    REAL(kind=8),INTENT(IN) :: EperpaxiLMr(l_max+1,nRstart:nRstop)

    !--- Local variables
    INTEGER :: nR,n,nTheta,nThetaStart,nThetaBlock,nThetaNHS
    CHARACTER(len=76) :: filename

    REAL(kind=8),DIMENSION(nRstart:nRstop) :: EperpR,EparR,EperpaxiR,EparaxiR
    REAL(kind=8),DIMENSION(n_r_max) :: EperpR_global,EparR_global
    REAL(kind=8),DIMENSION(n_r_max) :: EperpaxiR_global,EparaxiR_global
    REAL(kind=8),DIMENSION(nfs) :: Eperp,Epar,Eperpaxi,Eparaxi
    REAL(kind=8) :: EperpT,EparT,EperpaxT,EparaxT

    INTEGER :: i,sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1)
    !--- end of declaration

    DO nR=nRstart,nRstop
       EperpR(nR)   =0.d0
       EparR(nR)    =0.d0
       EparaxiR(nR) =0.d0
       EperpaxiR(nR)=0.d0
       DO n=1,nThetaBs ! Loop over theta blocks
          nTheta=(n-1)*sizeThetaB
          nThetaStart=nTheta+1
          CALL lmAS2pt(EperpLMr(1,nR),Eperp,nThetaStart,sizeThetaB)
          CALL lmAS2pt(EparLMr(1,nR),Epar,nThetaStart,sizeThetaB)
          CALL lmAS2pt(EperpaxiLMr(1,nR),Eperpaxi,nThetaStart,sizeThetaB)
          CALL lmAS2pt(EparaxiLMr(1,nR),Eparaxi,nThetaStart,sizeThetaB)
          DO nThetaBlock=1,sizeThetaB
             nTheta=nTheta+1
             nThetaNHS=(nTheta+1)/2
             EperpR(nR)=EperpR(nR)+gauss(nThetaNHS)*Eperp(nThetaBlock)
             EparR(nR) =EparR(nR) +gauss(nThetaNHS)* Epar(nThetaBlock)
             EperpaxiR(nR)=EperpaxiR(nR)+gauss(nThetaNHS)*Eperpaxi(nThetaBlock)
             EparaxiR(nR)=EparaxiR(nR)+gauss(nThetaNHS)*Eparaxi(nThetaBlock)
          END DO
       END DO
    END DO
    EperpR   =0.5d0*EperpR    ! Normalisation for the theta integration
    EparR    =0.5d0*EparR     ! Normalisation for the theta integration
    EperpaxiR=0.5d0*EperpaxiR ! Normalisation for the theta integration
    EparaxiR =0.5d0*EparaxiR  ! Normalisation for the theta integration

    sendcount  = (nRstop-nRstart+1)
    recvcounts = nr_per_rank
    recvcounts(n_procs-1) = (nr_per_rank+1)
    DO i=0,n_procs-1
       displs(i) = i*nr_per_rank
    END DO

    CALL MPI_GatherV(EperpR,sendcount,MPI_DOUBLE_PRECISION,&
        &           EperpR_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
        &           0,MPI_COMM_WORLD,ierr)
    CALL MPI_GatherV(EparR,sendcount,MPI_DOUBLE_PRECISION,&
        &           EparR_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
        &           0,MPI_COMM_WORLD,ierr)
    CALL MPI_GatherV(EperpaxiR,sendcount,MPI_DOUBLE_PRECISION,&
        &           EperpaxiR_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
        &           0,MPI_COMM_WORLD,ierr)
    CALL MPI_GatherV(EparaxiR,sendcount,MPI_DOUBLE_PRECISION,&
        &           EparaxiR_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
        &           0,MPI_COMM_WORLD,ierr)


    IF (rank == 0) THEN
       EperpT  =4.D0*pi*rInt(EperpR_global*r**2,n_r_max,dr_fac,i_costf_init,d_costf_init)
       EparT   =4.D0*pi*rInt(EparR_global*r**2,n_r_max,dr_fac,i_costf_init,d_costf_init)
       EperpaxT=4.D0*pi*rInt(EperpaxiR_global*r**2,n_r_max,dr_fac,i_costf_init,d_costf_init)
       EparaxT =4.D0*pi*rInt(EparaxiR_global*r**2,n_r_max,dr_fac,i_costf_init,d_costf_init)

       !-- Output
       IF ( l_save_out ) THEN
          OPEN(n_perpPar_file,FILE=perpPar_file,STATUS='UNKNOWN',POSITION='APPEND')
       END IF
       WRITE(n_perpPar_file,'(1P,D20.12,4D16.8)') &
            &  time*tScale,     & ! 1
            &  EperpT,EparT,    & ! 2,3
            &  EperpaxT,EparaxT   ! 4,5
       IF ( l_save_out ) CLOSE(n_perpPar_file)


       EperpMeanR    =EperpMeanR     +timePassed*EperpR_global
       EparMeanR     =EparMeanR      +timePassed*EparR_global
       EperpaxiMeanR =EperpaxiMeanR  +timePassed*EperpaxiR_global
       EparaxiMeanR  =EparaxiMeanR   +timePassed*EparaxiR_global
       IF ( l_stop_time ) THEN
           EperpMeanR     =EperpMeanR/timeNorm
           EparMeanR      =EparMeanR/timeNorm
           EperpaxiMeanR  =EperpaxiMeanR/timeNorm
           EparaxiMeanR   =EparaxiMeanR/timeNorm
           filename='perpParR.'//tag
           OPEN(99,FILE=filename,STATUS='UNKNOWN')
           DO nR=1,n_r_max
              WRITE(99,'(D20.10,4D20.12)')      &
                         &   r(nR),             &! 1) radius
                         &   EperpMeanR(nR),    &! 2) E perpendicular
                         &   EparMeanR(nR),     &! 3) E parallel
                         &   EperpaxiMeanR(nR), &! 4) E perp (axisymetric)
                         &   EparaxiMeanR(nR)    ! 5) E parallel (axisymetric)
           END DO
           CLOSE(99)
       END IF
    END IF


 END SUBROUTINE outPerpPar
 !********************************************
END MODULE outPerpPar_mod
