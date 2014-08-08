!$Id$
!********************************************************************
SUBROUTINE dtVrms(time,nRMS_sets)
  !********************************************************************

  !  +-------------+----------------+------------------------------------+
  !  |  For testing RMS balance and determining the necessary values of  |
  !  |  rDea and rCut one can use the l_RMStest=true options.            |
  !  |  In this case the poloidal and toroidal RMS dtV which also        |
  !  |  include the diffusion effects should be identical (as close as   |
  !  |  desired) to the RMS sum of forces stored in the Geo value and    |
  !  |  the Mag value, respectively.                                     |
  !  |  An additional tests is the Arc value which should be identical   |
  !  |  to the poloidal kinetic energy.                                  |
  !  |                                                                   |
  !  |  N.B. The second test with the Arc value cannot work so easily in |
  !  |       the anelastic version as the density enters now in the      |
  !  |       force balance. The first test should work, though.          |
  !  +-------------------------------------------------------------------+

  USE truncation
  USE radial_functions
  USE physical_parameters
  USE blocking
  USE logic
  USE RMS,ONLY: CorPolLMr,CorPol2hInt,CorPolAs2hInt,&
       & CorTor2hInt,CorTorAs2hInt,&
       & AdvPolLMr,AdvPol2hInt,AdvPolAs2hInt,&
       & AdvTor2hInt,AdvTorAs2hInt,&
       & LFPolLMr,LFPol2hInt,LFPolAs2hInt,&
       & LFTor2hInt,LFTorAs2hInt,&
       & BuoLMr,Buo2hInt,BuoAs2hInt,&
       & PreLMr,Pre2hInt,PreAs2hInt,&
       & GeoLMr,Geo2hInt,GeoAs2hInt,&
       & MagLMr,Mag2hInt,MagAs2hInt,&
       & ArcLMr,Arc2hInt, ArcAs2hInt,&
       & DifPolLMr,DifPol2hInt,DifPolAs2hInt,&
       & DifTor2hInt,DifTorAs2hInt,&
       & dtVPolLMr,dtVPol2hInt,dtVPolAs2hInt,&
       & dtVTor2hInt,dtVTorAs2hInt
  USE output_data
  USE const
  USE parallel_mod
  USE integration, ONLY: rInt_R
  USE communications, only:myallgather

  IMPLICIT NONE

  !-- Input variable:
  REAL(kind=8),INTENT(IN) :: time
  INTEGER,INTENT(INOUT) :: nRMS_sets

  !-- Output:
  REAL(kind=8) :: dtVPolRms,dtVPolAsRms
  REAL(kind=8) :: dtVTorRms,dtVTorAsRms
  REAL(kind=8) :: CorPolRms,CorPolAsRms
  REAL(kind=8) :: CorTorRms,CorTorAsRms
  REAL(kind=8) :: AdvPolRms,AdvPolAsRms
  REAL(kind=8) :: AdvTorRms,AdvTorAsRms
  REAL(kind=8) :: LFPolRms, LFPolAsRms
  REAL(kind=8) :: LFTorRms, LFTorAsRms
  REAL(kind=8) :: DifPolRms,DifPolAsRms
  REAL(kind=8) :: DifTorRms,DifTorAsRms
  REAL(kind=8) :: BuoRms,   BuoAsRms
  REAL(kind=8) :: PreRms,   PreAsRms
  REAL(kind=8) :: GeoRms,   GeoAsRms
  REAL(kind=8) :: MagRms,   MagAsRms
  REAL(kind=8) :: ArcRms,   ArcAsRms

  !-- Local:
  INTEGER :: nR,nRC
  !INTEGER :: n
  REAL(kind=8) :: volC
  REAL(kind=8) :: bal1,bal2,bal3

  COMPLEX(kind=8) :: workA(lm_max,n_r_max),workB(lm_max,n_r_max)
  INTEGER,DIMENSION(0:n_procs-1) :: recvcounts,displs
  REAL(kind=8) :: global_sum(n_r_max)
  INTEGER :: irank,sendcount
  !-- end of declaration
  !---------------------------------------------------------------------

  ! First gather all needed arrays on rank 0
  ! some more arrays to gather for the dtVrms routine
  ! we need some more fields for the dtBrms routine
  sendcount  = (nRstop-nRstart+1)*lm_max
  recvcounts = nr_per_rank*lm_max
  recvcounts(n_procs-1) = (nr_per_rank+1)*lm_max
  DO irank=0,n_procs-1
     displs(irank) = irank*nr_per_rank*lm_max
  END DO
  ! CorPolLMr,dtVPolLMr,AdvPolLMr,LFPolLMr,DifPolLMr,BuoLMr
  ! PreLMr,GeoLMr,MagLMr,ArcLMr
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_COMPLEX,&
       & CorPolLMr,recvcounts,displs,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_COMPLEX,&
       & AdvPolLMr,recvcounts,displs,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_COMPLEX,&
       & LFPolLMr,recvcounts,displs,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_COMPLEX,&
       & BuoLMr,recvcounts,displs,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_COMPLEX,&
       & PreLMr,recvcounts,displs,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_COMPLEX,&
       & GeoLMr,recvcounts,displs,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_COMPLEX,&
       & MagLMr,recvcounts,displs,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_COMPLEX,&
       & ArcLMr,recvcounts,displs,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)

  ! The following fields are only 1D and R distributed.
  sendcount  = (nRstop-nRstart+1)
  recvcounts = nr_per_rank
  recvcounts(n_procs-1) = (nr_per_rank+1)
  DO irank=0,n_procs-1
     displs(irank) = irank*nr_per_rank
  END DO
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
       & CorPol2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
       & CorPolAs2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
       & CorTor2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
       & CorTorAs2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
       & AdvPol2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
       & AdvPolAs2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
       & AdvTor2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
       & AdvTorAs2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
       & LFPol2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
       & LFPolAs2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
       & LFTor2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
       & LFTorAs2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
       & Buo2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
       & BuoAs2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
       & Pre2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
       & PreAs2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
       & Geo2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
       & GeoAs2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
       & Mag2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
       & MagAs2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
       & Arc2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  CALL MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
       & ArcAs2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  ! The following fields are LM distributed and have to be gathered:
  ! dtVPolLMr, DifPolLMr

  CALL myAllGather(dtVPolLMr,lm_max,n_r_max)
  CALL myAllGather(DifPolLMr,lm_max,n_r_max)

  CALL mpi_reduce(dtVPol2hInt(1,1),global_sum,n_r_max,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &0,MPI_COMM_WORLD,ierr)
  IF (rank.EQ.0) dtVPol2hInt(:,1)=global_sum
  CALL mpi_reduce(dtVPolAs2hInt(1,1),global_sum,n_r_max,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &0,MPI_COMM_WORLD,ierr)
  IF (rank.EQ.0) dtVPolAs2hInt(:,1)=global_sum
  CALL mpi_reduce(dtVTor2hInt(1,1),global_sum,n_r_max,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &0,MPI_COMM_WORLD,ierr)
  IF (rank.EQ.0) dtVTor2hInt(:,1)=global_sum
  CALL mpi_reduce(dtVTorAs2hInt(1,1),global_sum,n_r_max,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &0,MPI_COMM_WORLD,ierr)
  IF (rank.EQ.0) dtVTorAs2hInt(:,1)=global_sum
  CALL mpi_reduce(DifPol2hInt(1,1),global_sum,n_r_max,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &0,MPI_COMM_WORLD,ierr)
  IF (rank.EQ.0) DifPol2hInt(:,1)=global_sum
  CALL mpi_reduce(DifPolAs2hInt(1,1),global_sum,n_r_max,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &0,MPI_COMM_WORLD,ierr)
  IF (rank.EQ.0) DifPolAs2hInt(:,1)=global_sum
  CALL mpi_reduce(DifTor2hInt(1,1),global_sum,n_r_max,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &0,MPI_COMM_WORLD,ierr)
  IF (rank.EQ.0) DifTor2hInt(:,1)=global_sum
  CALL mpi_reduce(DifTorAs2hInt(1,1),global_sum,n_r_max,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &0,MPI_COMM_WORLD,ierr)
  IF (rank.EQ.0) DifTorAs2hInt(:,1)=global_sum

  IF (rank.EQ.0) THEN

     !WRITE(*,"(A,ES22.14)") "dtVPol2hInt = ",SUM(dtVPol2hInt)
     IF ( nRMS_sets == 0 ) THEN
        nRMS_sets=nRMS_sets+1
        !--- Initialize new cut-back grid:
        CALL init_rNB(r,n_r_max,n_cheb_max,rCut,rDea, &
             rC,n_r_maxC,n_cheb_maxC,nCut, &
             dr_facC,i_costf_initC,nDi_costf1, &
             d_costf_initC,nDd_costf1)
     END IF
     nRC=nCut+1
     volC=4.D0/3.D0*pi*(r(1+nCut)**3-r(n_r_max-nCut)**3)

     IF ( l_conv ) THEN

        !------ Coriolis force
        IF ( l_corr ) THEN
           CALL get_drNS_R( CorPolLMr(1,nRC),workA(1,nRC), &
                &           lm_max_real,1,lm_max_real,n_r_maxC,n_cheb_maxC, &
                &           dr_facC,workB,i_costf_initC,d_costf_initC)
           DO nR=1,n_r_maxC
              CALL hInt2dPol(workA(1,nR+nCut),2,lm_max, &
                   CorPol2hInt(nR+nCut),CorPolAs2hInt(nR+nCut),st_map)
           END DO
           CorPolRms=rInt_R(CorPol2hInt(nRC),n_r_maxC, &
                n_cheb_maxC,dr_facC, &
                i_costf_initC,d_costf_initC)
           CorPolAsRms=rInt_R(CorPolAs2hInt(nRC),n_r_maxC, &
                n_cheb_maxC,dr_facC, &
                i_costf_initC,d_costf_initC)
           CorTorRms  =rInt_R(CorTor2hInt(nRC),n_r_maxC, &
                n_cheb_maxC,dr_facC, &
                i_costf_initC,d_costf_initC)
           CorTorAsRms=rInt_R(CorTorAs2hInt(nRC),n_r_maxC, &
                n_cheb_maxC,dr_facC, &
                i_costf_initC,d_costf_initC)
           CorPolRms  =DSQRT(CorPolRms  /volC)
           CorPolAsRms=DSQRT(CorPolAsRms/volC)
           CorTorRms  =DSQRT(CorTorRms  /volC)
           CorTorAsRms=DSQRT(CorTorAsRms/volC)
        END IF

        !------ Advection:
        IF ( l_conv_nl ) THEN
           CALL get_drNS_R(AdvPolLMr(1,nRC),workA(1,nRC), &
                &          lm_max_real,1,lm_max_real,n_r_maxC,n_cheb_maxC, &
                &          dr_facC,workB,i_costf_initC,d_costf_initC)
           DO nR=1,n_r_maxC
              CALL hInt2dPol(workA(1,nR+nCut),2,lm_max, &
                   &         AdvPol2hInt(nR+nCut),AdvPolAs2hInt(nR+nCut),st_map)
           END DO
           AdvPolRms  =rInt_R(AdvPol2hInt(nRC),n_r_maxC, &
                n_cheb_maxC,dr_facC, &
                i_costf_initC,d_costf_initC)
           AdvPolAsRms=rInt_R(AdvPolAs2hInt(nRC),n_r_maxC, &
                n_cheb_maxC,dr_facC, &
                i_costf_initC,d_costf_initC)
           AdvTorRms  =rInt_R(AdvTor2hInt(nRC),n_r_maxC, &
                n_cheb_maxC,dr_facC, &
                i_costf_initC,d_costf_initC)
           AdvTorAsRms=rInt_R(AdvTorAs2hInt(nRC),n_r_maxC, &
                n_cheb_maxC,dr_facC, &
                i_costf_initC,d_costf_initC)
           AdvPolRms  =DSQRT(AdvPolRms  /volC)
           AdvPolAsRms=DSQRT(AdvPolAsRms/volC)
           AdvTorRms  =DSQRT(AdvTorRms  /volC)
           AdvTorAsRms=DSQRT(AdvTorAsRms/volC)
        END IF

        !------ Lorentz force:
        IF ( l_mag_LF ) THEN
           CALL get_drNS_R(LFPolLMr(1,nRC),workA(1,nRC), &
                &          lm_max_real,1,lm_max_real,n_r_maxC,n_cheb_maxC, &
                &          dr_facC,workB,i_costf_initC,d_costf_initC)
           DO nR=1,n_r_maxC
              CALL hInt2dPol( workA(1,nR+nCut),2,lm_max, &
                   &          LFPol2hInt(nR+nCut),LFPolAs2hInt(nR+nCut),st_map)
           END DO
           LFPolRms  =rInt_R(LFPol2hInt(nRC),n_r_maxC, &
                n_cheb_maxC,dr_facC, &
                i_costf_initC,d_costf_initC)
           LFPolAsRms=rInt_R(LFPolAs2hInt(nRC),n_r_maxC, &
                n_cheb_maxC,dr_facC, &
                i_costf_initC,d_costf_initC)
           LFTorRms  =rInt_R(LFTor2hInt(nRC),n_r_maxC, &
                n_cheb_maxC,dr_facC, &
                i_costf_initC,d_costf_initC)
           LFTorAsRms=rInt_R(LFTorAs2hInt(nRC),n_r_maxC, &
                n_cheb_maxC,dr_facC, &
                i_costf_initC,d_costf_initC)
           LFPolRms  =DSQRT(LFPolRms  /volC)
           LFPolAsRms=DSQRT(LFPolAsRms/volC)
           LFTorRms  =DSQRT(LFTorRms  /volC)
           LFTorAsRms=DSQRT(LFTorAsRms/volC)
        ELSE
           LFPolRms  =0.D0
           LFPolAsRms=0.D0
           LFTorRms  =0.D0
           LFTorAsRms=0.D0
        END IF

        !------ Buoyancy:
        IF ( l_heat ) THEN
           CALL get_drNS_R(BuoLMr(1,nRC),workA(1,nRC), &
                &          lm_max_real,1,lm_max_real,n_r_maxC,n_cheb_maxC, &
                &          dr_facC,workB,i_costf_initC,d_costf_initC)
           DO nR=1,n_r_maxC
              CALL hInt2dPol(workA(1,nR+nCut),2,lm_max, &
                   &         Buo2hInt(nR+nCut),BuoAs2hInt(nR+nCut),st_map)
              Buo2hInt(nR)  =Buo2hInt(nR)*rgrav(nR)**2
              BuoAs2hInt(nR)=BuoAs2hInt(nR)*rgrav(nR)**2
           END DO
           BuoRms  =rInt_R(Buo2hInt(nRC),n_r_maxC, &
                n_cheb_maxC,dr_facC, &
                i_costf_initC,d_costf_initC)
           BuoAsRms=rInt_R(BuoAs2hInt(nRC),n_r_maxC, &
                n_cheb_maxC,dr_facC, &
                i_costf_initC,d_costf_initC)
           BuoRms  =DSQRT(BuoRms  /volC)
           BuoAsRms=DSQRT(BuoAsRms/volC)
        ELSE
           BuoRms=0.D0
           BuoAsRms=0.D0
        END IF

        !------ Pressure gradient:
        CALL get_drNS_R(       PreLMr(1,nRC),workA(1,nRC), &
             lm_max_real,1,lm_max_real,n_r_maxC,n_cheb_maxC, &
             dr_facC,workB,i_costf_initC,d_costf_initC)
        DO nR=1,n_r_maxC
           CALL hInt2dPol(     workA(1,nR+nCut),2,lm_max, &
                Pre2hInt(nR+nCut),PreAs2hInt(nR+nCut),st_map)
        END DO
        PreRms  =rInt_R(Pre2hInt(nRC),n_r_maxC, &
             n_cheb_maxC,dr_facC, &
             i_costf_initC,d_costf_initC)
        PreAsRms=rInt_R(PreAs2hInt(nRC),n_r_maxC, &
             n_cheb_maxC,dr_facC, &
             i_costf_initC,d_costf_initC)
        PreRms  =DSQRT(PreRms  /volC)
        PreAsRms=DSQRT(PreAsRms/volC)

        !------ Geostrophic balance:
        CALL get_drNS_R(       GeoLMr(1,nRC),workA(1,nRC), &
             lm_max_real,1,lm_max_real,n_r_maxC,n_cheb_maxC, &
             dr_facC,workB,i_costf_initC,d_costf_initC)
        DO nR=1,n_r_maxC
           CALL hInt2dPol(     workA(1,nR+nCut),2,lm_max, &
                Geo2hInt(nR+nCut),GeoAs2hInt(nR+nCut),st_map)
        END DO
        GeoRms  =rInt_R(Geo2hInt(nRC),n_r_maxC, &
             n_cheb_maxC,dr_facC, &
             i_costf_initC,d_costf_initC)
        GeoAsRms=rInt_R(GeoAs2hInt(nRC),n_r_maxC, &
             n_cheb_maxC,dr_facC, &
             i_costf_initC,d_costf_initC)
        GeoRms  =DSQRT(GeoRms  /volC)
        GeoAsRms=DSQRT(GeoAsRms/volC)

        !------ Magnetostrophic balance:
        IF ( .NOT. l_RMStest ) THEN
           CALL get_drNS_R(       MagLMr(1,nRC),workA(1,nRC), &
                lm_max_real,1,lm_max_real,n_r_maxC,n_cheb_maxC, &
                dr_facC,workB,i_costf_initC,d_costf_initC)
           DO nR=1,n_r_maxC
              CALL hInt2dPol(     workA(1,nR+nCut),2,lm_max, &
                   Mag2hInt(nR+nCut),MagAs2hInt(nR+nCut),st_map)
           END DO
        END IF
        MagRms  =rInt_R(Mag2hInt(nRC),n_r_maxC, &
             n_cheb_maxC,dr_facC, &
             i_costf_initC,d_costf_initC)
        MagAsRms=rInt_R(MagAs2hInt(nRC),n_r_maxC, &
             n_cheb_maxC,dr_facC, &
             i_costf_initC,d_costf_initC)
        MagRms  =DSQRT(MagRms  /volC)
        MagAsRms=DSQRT(MagAsRms/volC)

        !------ Archemidian balance:
        CALL get_drNS_R( ArcLMr(1,nRC),workA(1,nRC), &
             &           lm_max_real,1,lm_max_real,n_r_maxC,n_cheb_maxC, &
             &           dr_facC,workB,i_costf_initC,d_costf_initC)
        DO nR=1,n_r_maxC
           CALL hInt2dPol( workA(1,nR+nCut),2,lm_max, &
                &          Arc2hInt(nR+nCut),ArcAs2hInt(nR+nCut),st_map)
        END DO
        ArcRms  =rInt_R(Arc2hInt(nRC),n_r_maxC, &
             n_cheb_maxC,dr_facC, &
             i_costf_initC,d_costf_initC)
        ArcAsRms=rInt_R(ArcAs2hInt(nRC),n_r_maxC, &
             n_cheb_maxC,dr_facC, &
             i_costf_initC,d_costf_initC)
        IF ( l_RMStest ) THEN
           ArcRms  =ArcRms/2.D0
           ArcAsRms=ArcAsRms/2.D0
        ELSE
           ArcRms  =DSQRT(ArcRms  /volC)
           ArcAsRms=DSQRT(ArcAsRms/volC)
        END IF

        !------ Diffusion:
        CALL get_drNS_R(DifPolLMr(1,nRC),workA(1,nRC), &
             &          lm_max_real,1,lm_max_real,n_r_maxC,n_cheb_maxC, &
             &          dr_facC,workB,i_costf_initC,d_costf_initC)
        DO nR=1,n_r_maxC
           CALL hInt2dPol( workA(1,nR+nCut),2,lm_max, &
                &          DifPol2hInt(nR+nCut,1),DifPolAs2hInt(nR+nCut,1),lo_map)
           !DO n=2,nThreadsLMmax
           !   DifPol2hInt(nR+nCut,1)  =DifPol2hInt(nR+nCut,1) + &
           !        DifPol2hInt(nR+nCut,n)
           !   DifPolAs2hInt(nR+nCut,1)=DifPolAs2hInt(nR+nCut,1) + &
           !        DifPolAs2hInt(nR+nCut,n)
           !   DifTor2hInt(nR+nCut,1)  =DifTor2hInt(nR+nCut,1) + &
           !        DifTor2hInt(nR+nCut,n)
           !   DifTorAs2hInt(nR+nCut,1)=DifTorAs2hInt(nR+nCut,1) + &
           !        DifTorAs2hInt(nR+nCut,n)
           !END DO
        END DO
        DifPolRms  =rInt_R(DifPol2hInt(nRC,1),n_r_maxC, &
             n_cheb_maxC,dr_facC, &
             i_costf_initC,d_costf_initC)
        DifPolAsRms=rInt_R(DifPolAs2hInt(nRC,1),n_r_maxC, &
             n_cheb_maxC,dr_facC, &
             i_costf_initC,d_costf_initC)
        DifTorRms  =rInt_R(DifTor2hInt(nRC,1),n_r_maxC, &
             n_cheb_maxC,dr_facC, &
             i_costf_initC,d_costf_initC)
        DifTorAsRms=rInt_R(DifTorAs2hInt(nRC,1),n_r_maxC, &
             n_cheb_maxC,dr_facC, &
             i_costf_initC,d_costf_initC)
        DifPolRms  =SQRT(DifPolRms  /volC)
        DifPolAsRms=SQRT(DifPolAsRms/volC)
        DifTorRms  =SQRT(DifTorRms  /volC)
        DifTorAsRms=SQRT(DifTorAsRms/volC)

        !------ Flow changes: Inertia - Advection
        CALL get_drNS_R( dtVPolLMr(1,nRC),workA(1,nRC), &
             &           lm_max_real,1,lm_max_real,n_r_maxC,n_cheb_maxC, &
             &           dr_facC,workB,i_costf_initC,d_costf_initC)
        DO nR=1,n_r_maxC
           !WRITE(*,"(A,I3,3ES22.14)") "workA = ",nR+nCut,SUM(workA(2:lm_max,nR+nCut)),&
           !     & dtVPol2hInt(nR+nCut,1)
           CALL hInt2dPol( workA(1,nR+nCut),2,lm_max, &
                &          dtVPol2hInt(nR+nCut,1),dtVPolAs2hInt(nR+nCut,1),lo_map)
           !DO n=2,nThreadsLMmax
           !   dtVPol2hInt(nR+nCut,1)  =dtVPol2hInt(nR+nCut,1) + &
           !        &                   dtVPol2hInt(nR+nCut,n)
           !   dtVPolAs2hInt(nR+nCut,1)=dtVPolAs2hInt(nR+nCut,1) + &
           !        &                   dtVPolAs2hInt(nR+nCut,n)
           !   dtVTor2hInt(nR+nCut,1)  =dtVTor2hInt(nR+nCut,1) + &
           !        &                   dtVTor2hInt(nR+nCut,n)
           !   dtVTorAs2hInt(nR+nCut,1)=dtVTorAs2hInt(nR+nCut,1) + &
           !        &                   dtVTorAs2hInt(nR+nCut,n)
           !END DO
           !WRITE(*,"(A,I3,ES22.14)") "after = ",nR+nCut,dtVPol2hInt(nR+nCut,1)
        END DO
        !WRITE(*,"(A,I3,ES22.14)") "dtVPol2hInt(nRC,1) = ",nRC,dtVPol2hInt(nRC,1)
        dtVPolRms  =rInt_R(dtVPol2hInt(nRC,1),n_r_maxC, &
             &             n_cheb_maxC,dr_facC, &
             &             i_costf_initC,d_costf_initC)
        !WRITE(*,"(A,ES22.14)") "dtVPolRms = ",dtVPolRms
        dtVPolAsRms=rInt_R(dtVPolAs2hInt(nRC,1),n_r_maxC, &
             &             n_cheb_maxC,dr_facC, &
             &             i_costf_initC,d_costf_initC)
        dtVTorRms  =rInt_R(dtVTor2hInt(nRC,1),n_r_maxC, &
             &             n_cheb_maxC,dr_facC, &
             &             i_costf_initC,d_costf_initC)
        dtVTorAsRms=rInt_R(dtVTorAs2hInt(nRC,1),n_r_maxC, &
             &             n_cheb_maxC,dr_facC, &
             &             i_costf_initC,d_costf_initC)
        dtVPolRms  =SQRT(dtVPolRms  /volC)
        dtVPolAsRms=SQRT(dtVPolAsRms/volC)
        dtVTorRms  =SQRT(dtVTorRms  /volC)
        dtVTorAsRms=SQRT(dtVTorAsRms/volC)

        !----- Output:
        IF ( l_save_out) THEN
           OPEN(n_dtvrms_file,FILE=dtvrms_file,FORM='FORMATTED', &
                STATUS='UNKNOWN',POSITION='APPEND')
        END IF
        IF ( l_RMStest ) THEN
           WRITE(n_dtvrms_file,'(1P,D20.12,15D16.8)') &
                time, &
                dtVPolRms,dtVTorRms, &
                CorPolRms,CorTorRms, &
                LFPolRms,LFTorRms, &
                AdvPolRms,AdvTorRms, &
                DifPolRms,DifTorRms, &
                BuoRms,PreRms, &
                GeoRms, &
                MagRms, &
                ArcRms
        ELSE
           WRITE(n_dtvrms_file,'(1P,D20.12,15D16.8)') &
                time, &
                dtVPolRms,dtVTorRms, & !2,3
                CorPolRms,CorTorRms, & !4,5
                LFPolRms,LFTorRms, &   !6,7
                AdvPolRms,AdvTorRms, & !8,9
                DifPolRms,DifTorRms, & !10,11
                BuoRms,PreRms, &       !12,13
                GeoRms/(CorPolRms+PreRms), & !14
                MagRms/(CorPolRms+PreRms+LFPolRms), &  !15
                ArcRms/(CorPolRms+PreRms+LFPolRms+BuoRms) !16
        END IF
        IF ( l_save_out) THEN
           CLOSE(n_dtvrms_file)
        END IF
        IF ( l_save_out) THEN
           OPEN(n_dtvasrms_file,FILE=dtvasrms_file,FORM='FORMATTED', &
                STATUS='UNKNOWN',POSITION='APPEND')
        END IF
        IF ( l_RMStest ) THEN
           WRITE(n_dtvasrms_file,'(1P,D20.12,15D16.8)') &
                time, &
                dtVPolAsRms,dtVTorAsRms, &
                CorPolAsRms,CorTorAsRms, &
                LFPolAsRms,LFTorAsRms, &
                AdvPolAsRms,AdvTorAsRms, &
                DifPolAsRms,DifTorAsRms, &
                BuoAsRms,PreAsRms, &
                GeoAsRms, &
                MagAsRms, &
                ArcAsRms
        ELSE
           IF ( PreAsRms/= 0d0 ) THEN
              bal1=GeoAsRms/(CorPolAsRms+PreAsRms)
              bal2=MagAsRms/(CorPolAsRms+PreAsRms+LFPolAsRms)
              bal3=ArcAsRms/(CorPolAsRms+PreAsRms+LFPolAsRms+BuoAsRms)
           ELSE
              bal1=0d0
              bal2=0d0
              bal3=0d0
           END IF
           WRITE(n_dtvasrms_file,'(1P,D20.12,15D16.8)') &
                time, &
                dtVPolAsRms,dtVTorAsRms, &
                CorPolAsRms,CorTorAsRms, &
                LFPolAsRms,LFTorAsRms, &
                AdvPolAsRms,AdvTorAsRms, &
                DifPolAsRms,DifTorAsRms, &
                BuoAsRms,PreAsRms, &
                bal1,bal2,bal3
        END IF
        IF ( l_save_out) THEN
           CLOSE(n_dtvasrms_file)
        END IF

     END IF

  END IF
  RETURN
end SUBROUTINE dtVrms

!-----------------------------------------------------------------------------
