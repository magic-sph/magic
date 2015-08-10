!$Id$
module out_RMS

   use truncation, only: lm_max, n_r_max, lm_max_dtB, n_r_max_dtB, &
                         n_cheb_max, lm_maxMag, n_theta_max, minc, &
                         n_r_maxMag, n_phi_max
   use parallel_mod, only: MPI_IN_PLACE, n_procs, MPI_SUM,     & 
                           MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, &
                           MPI_DOUBLE_PRECISION, ierr,         &
                           nr_per_rank, rank
   use radial_data, only: nRstop, nRstart
   use radial_functions, only: n_r_maxC, i_costf_init, d_costf_init,  &
                               drx, r, r_CMB, i_costf_initC, rgrav,   &
                               rC, n_cheb_maxC, d_costf_initC,        &
                               nDi_costf1, nDd_costf1, dr_facC, nCut
   use physical_parameters, only: ra, ek, pr, prmag, radratio
   use blocking, only: st_map, nThetaBs, nfs, sizeThetaB, lo_map, lm2
   use logic, only: l_save_out, l_RMStest, l_heat, l_conv_nl, l_mag_LF, &
                    l_conv, l_corr
   use RMS, only: CorPolLMr, CorPol2hInt, CorPolAs2hInt, CorTor2hInt,   &
                  CorTorAs2hInt, AdvPolLMr, AdvPol2hInt, AdvPolAs2hInt, &
                  AdvTor2hInt, AdvTorAs2hInt, LFPolLMr, LFPol2hInt,     &
                  LFPolAs2hInt, LFTor2hInt, LFTorAs2hInt, BuoLMr,       &
                  Buo2hInt, BuoAs2hInt, PreLMr, Pre2hInt, PreAs2hInt,   &
                  GeoLMr, Geo2hInt, GeoAs2hInt, MagLMr, Mag2hInt,       &
                  MagAs2hInt, ArcLMr, Arc2hInt, ArcAs2hInt, DifPolLMr,  &
                  DifPol2hInt, DifPolAs2hInt, DifTor2hInt,              &
                  DifTorAs2hInt, dtVPolLMr, dtVPol2hInt, dtVPolAs2hInt, &
                  dtVTor2hInt, dtVTorAs2hInt, dtBPolLMr, dtBPol2hInt,   &
                  dtBPolAs2hInt, dtBTor2hInt, dtBTorAs2hInt

   use dtB_mod, only: PstrLM, TstrLM, PadvLM, TadvLM, TomeLM, PdifLM,  &
                      TdifLM, PstrRms,TstrRms, PstrAsRms, TstrAsRms,   &
                      PadvRms, TadvRms, PadvAsRms, TadvAsRms, PdifRms, &
                      TdifRms, PdifAsRms, TdifAsRms, TomeRms,          &
                      TomeAsRms, dtB_gather_Rloc_on_rank0
                                                                  
   use num_param, only: tScale
   use horizontal_data, only: phi, theta_ord
   use output_data, only: TAG, runid, n_dtdrms_file, dtdrms_file,        &
                          n_dtvasrms_file, dtvasrms_file, n_dtvrms_file, &
                          dtvrms_file, rCut, rDea, n_dtbrms_file,        &
                          dtbrms_file
   use const, only: pi, vol_oc
   use integration, only: rInt_R
   use communications, only: myallgather
   use RMS_helpers, only: init_rNB, hInt2dPol, get_PolTorRms, get_PASLM, &
                          get_RAS
   use radial_der, only: get_drNS

   implicit none

   private

   public :: dtVrms, dtBrms

contains
 
   subroutine dtVrms(time,nRMS_sets)
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

      !-- Input variable:
      real(kind=8), intent(in) :: time
      integer,      intent(inout) :: nRMS_sets
    
      !-- Output:
      real(kind=8) :: dtVPolRms,dtVPolAsRms
      real(kind=8) :: dtVTorRms,dtVTorAsRms
      real(kind=8) :: CorPolRms,CorPolAsRms
      real(kind=8) :: CorTorRms,CorTorAsRms
      real(kind=8) :: AdvPolRms,AdvPolAsRms
      real(kind=8) :: AdvTorRms,AdvTorAsRms
      real(kind=8) :: LFPolRms, LFPolAsRms
      real(kind=8) :: LFTorRms, LFTorAsRms
      real(kind=8) :: DifPolRms,DifPolAsRms
      real(kind=8) :: DifTorRms,DifTorAsRms
      real(kind=8) :: BuoRms,   BuoAsRms
      real(kind=8) :: PreRms,   PreAsRms
      real(kind=8) :: GeoRms,   GeoAsRms
      real(kind=8) :: MagRms,   MagAsRms
      real(kind=8) :: ArcRms,   ArcAsRms
    
      !-- Local:
      integer :: nR,nRC
      !integer :: n
      real(kind=8) :: volC
      real(kind=8) :: bal1,bal2,bal3
    
      complex(kind=8) :: workA(lm_max,n_r_max),workB(lm_max,n_r_max)
      integer :: recvcounts(0:n_procs-1),displs(0:n_procs-1)
      real(kind=8) :: global_sum(n_r_max)
      integer :: irank,sendcount
    
      ! First gather all needed arrays on rank 0
      ! some more arrays to gather for the dtVrms routine
      ! we need some more fields for the dtBrms routine
      sendcount  = (nRstop-nRstart+1)*lm_max
      recvcounts = nr_per_rank*lm_max
      recvcounts(n_procs-1) = (nr_per_rank+1)*lm_max
      do irank=0,n_procs-1
         displs(irank) = irank*nr_per_rank*lm_max
      end do
      ! CorPolLMr,dtVPolLMr,AdvPolLMr,LFPolLMr,DifPolLMr,BuoLMr
      ! PreLMr,GeoLMr,MagLMr,ArcLMr
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_COMPLEX,&
           & CorPolLMr,recvcounts,displs,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_COMPLEX,&
           & AdvPolLMr,recvcounts,displs,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_COMPLEX,&
           & LFPolLMr,recvcounts,displs,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_COMPLEX,&
           & BuoLMr,recvcounts,displs,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_COMPLEX,&
           & PreLMr,recvcounts,displs,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_COMPLEX,&
           & GeoLMr,recvcounts,displs,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_COMPLEX,&
           & MagLMr,recvcounts,displs,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_COMPLEX,&
           & ArcLMr,recvcounts,displs,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
    
      ! The following fields are only 1D and R distributed.
      sendcount  = (nRstop-nRstart+1)
      recvcounts = nr_per_rank
      recvcounts(n_procs-1) = (nr_per_rank+1)
      do irank=0,n_procs-1
         displs(irank) = irank*nr_per_rank
      end do
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
           & CorPol2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
           & CorPolAs2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
           & CorTor2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
           & CorTorAs2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
           & AdvPol2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
           & AdvPolAs2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
           & AdvTor2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
           & AdvTorAs2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
           & LFPol2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
           & LFPolAs2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
           & LFTor2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
           & LFTorAs2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
           & Buo2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
           & BuoAs2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
           & Pre2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
           & PreAs2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
           & Geo2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
           & GeoAs2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
           & Mag2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
           & MagAs2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
           & Arc2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,&
           & ArcAs2hInt,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    
      ! The following fields are LM distributed and have to be gathered:
      ! dtVPolLMr, DifPolLMr
    
      call myAllGather(dtVPolLMr,lm_max,n_r_max)
      call myAllGather(DifPolLMr,lm_max,n_r_max)
    
      call mpi_reduce(dtVPol2hInt(1,1),global_sum,n_r_max,MPI_DOUBLE_PRECISION, &
           &          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) dtVPol2hInt(:,1)=global_sum
      call mpi_reduce(dtVPolAs2hInt(1,1),global_sum,n_r_max,MPI_DOUBLE_PRECISION,&
           &          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) dtVPolAs2hInt(:,1)=global_sum
      call mpi_reduce(dtVTor2hInt(1,1),global_sum,n_r_max,MPI_DOUBLE_PRECISION, &
           &          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) dtVTor2hInt(:,1)=global_sum
      call mpi_reduce(dtVTorAs2hInt(1,1),global_sum,n_r_max,MPI_DOUBLE_PRECISION, &
           &          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) dtVTorAs2hInt(:,1)=global_sum
      call mpi_reduce(DifPol2hInt(1,1),global_sum,n_r_max,MPI_DOUBLE_PRECISION, &
           &          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) DifPol2hInt(:,1)=global_sum
      call mpi_reduce(DifPolAs2hInt(1,1),global_sum,n_r_max,MPI_DOUBLE_PRECISION, &
           &          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) DifPolAs2hInt(:,1)=global_sum
      call mpi_reduce(DifTor2hInt(1,1),global_sum,n_r_max,MPI_DOUBLE_PRECISION, &
           &          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) DifTor2hInt(:,1)=global_sum
      call mpi_reduce(DifTorAs2hInt(1,1),global_sum,n_r_max,MPI_DOUBLE_PRECISION,&
           &          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) DifTorAs2hInt(:,1)=global_sum
    
      if ( rank == 0 ) then
    
         !write(*,"(A,ES22.14)") "dtVPol2hInt = ",SUM(dtVPol2hInt)
         if ( nRMS_sets == 0 ) then
            nRMS_sets=nRMS_sets+1
            !--- Initialize new cut-back grid:
            call init_rNB(r,n_r_max,n_cheb_max,rCut,rDea, &
                            rC,n_r_maxC,n_cheb_maxC,nCut, &
                        dr_facC,i_costf_initC,nDi_costf1, &
                                d_costf_initC,nDd_costf1)
         end if
         nRC=nCut+1
         volC=4.D0/3.D0*pi*(r(1+nCut)**3-r(n_r_max-nCut)**3)
    
         if ( l_conv ) then
    
            !------ Coriolis force
            if ( l_corr ) then
               call get_drNS( CorPolLMr(1,nRC),workA(1,nRC),        &
                    &           lm_max,1,lm_max,n_r_maxC,n_cheb_maxC, &
                    &           workB,i_costf_initC,d_costf_initC,dr_facC)
               do nR=1,n_r_maxC
                  call hInt2dPol(workA(1,nR+nCut),2,lm_max,CorPol2hInt(nR+nCut), &
                                 CorPolAs2hInt(nR+nCut),st_map)
               end do
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
               CorPolRms  =dsqrt(CorPolRms  /volC)
               CorPolAsRms=dsqrt(CorPolAsRms/volC)
               CorTorRms  =dsqrt(CorTorRms  /volC)
               CorTorAsRms=dsqrt(CorTorAsRms/volC)
            end if
    
            !------ Advection:
            if ( l_conv_nl ) then
               call get_drNS(AdvPolLMr(1,nRC),workA(1,nRC),          &
                    &          lm_max,1,lm_max,n_r_maxC,n_cheb_maxC,   &
                    &          workB,i_costf_initC,d_costf_initC,dr_facC)
               do nR=1,n_r_maxC
                  call hInt2dPol(workA(1,nR+nCut),2,lm_max,AdvPol2hInt(nR+nCut), &
                                 AdvPolAs2hInt(nR+nCut),st_map)
               end do
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
               AdvPolRms  =dsqrt(AdvPolRms  /volC)
               AdvPolAsRms=dsqrt(AdvPolAsRms/volC)
               AdvTorRms  =dsqrt(AdvTorRms  /volC)
               AdvTorAsRms=dsqrt(AdvTorAsRms/volC)
            end if
    
            !------ Lorentz force:
            if ( l_mag_LF ) then
               call get_drNS(LFPolLMr(1,nRC),workA(1,nRC),lm_max,1,lm_max,&
                    &        n_r_maxC,n_cheb_maxC,workB,i_costf_initC,    &
                    &        d_costf_initC,dr_facC)
               do nR=1,n_r_maxC
                  call hInt2dPol( workA(1,nR+nCut),2,lm_max,LFPol2hInt(nR+nCut), &
                                  LFPolAs2hInt(nR+nCut),st_map)
               end do
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
               LFPolRms  =dsqrt(LFPolRms  /volC)
               LFPolAsRms=dsqrt(LFPolAsRms/volC)
               LFTorRms  =dsqrt(LFTorRms  /volC)
               LFTorAsRms=dsqrt(LFTorAsRms/volC)
            else
               LFPolRms  =0.D0
               LFPolAsRms=0.D0
               LFTorRms  =0.D0
               LFTorAsRms=0.D0
            end if
    
            !------ Buoyancy:
            if ( l_heat ) then
               call get_drNS(BuoLMr(1,nRC),workA(1,nRC),lm_max,1,lm_max, &
                    &        n_r_maxC,n_cheb_maxC,workB,i_costf_initC,   &
                    &        d_costf_initC,dr_facC)
               do nR=1,n_r_maxC
                  call hInt2dPol(workA(1,nR+nCut),2,lm_max, Buo2hInt(nR+nCut), &
                                 BuoAs2hInt(nR+nCut),st_map)
                  Buo2hInt(nR)  =Buo2hInt(nR)*rgrav(nR)**2
                  BuoAs2hInt(nR)=BuoAs2hInt(nR)*rgrav(nR)**2
               end do
               BuoRms  =rInt_R(Buo2hInt(nRC),n_r_maxC, &
                                  n_cheb_maxC,dr_facC, &
                          i_costf_initC,d_costf_initC)
               BuoAsRms=rInt_R(BuoAs2hInt(nRC),n_r_maxC, &
                                    n_cheb_maxC,dr_facC, &
                            i_costf_initC,d_costf_initC)
               BuoRms  =dsqrt(BuoRms  /volC)
               BuoAsRms=dsqrt(BuoAsRms/volC)
            else
               BuoRms=0.D0
               BuoAsRms=0.D0
            end if
    
            !------ Pressure gradient:
            call get_drNS(PreLMr(1,nRC),workA(1,nRC),lm_max,1,  &
                 &        lm_max,n_r_maxC,n_cheb_maxC,workB,    &
                 &        i_costf_initC,d_costf_initC,dr_facC)
            do nR=1,n_r_maxC
               call hInt2dPol(workA(1,nR+nCut),2,lm_max,Pre2hInt(nR+nCut), &
                             PreAs2hInt(nR+nCut),st_map)
            end do
            PreRms  =rInt_R(Pre2hInt(nRC),n_r_maxC, &
                               n_cheb_maxC,dr_facC, &
                       i_costf_initC,d_costf_initC)
            PreAsRms=rInt_R(PreAs2hInt(nRC),n_r_maxC, &
                                 n_cheb_maxC,dr_facC, &
                         i_costf_initC,d_costf_initC)
            PreRms  =dsqrt(PreRms  /volC)
            PreAsRms=dsqrt(PreAsRms/volC)
    
            !------ Geostrophic balance:
            call get_drNS(GeoLMr(1,nRC),workA(1,nRC),lm_max,1, &
                 &        lm_max,n_r_maxC,n_cheb_maxC,workB,   &
                 &        i_costf_initC,d_costf_initC,dr_facC)
            do nR=1,n_r_maxC
               call hInt2dPol(workA(1,nR+nCut),2,lm_max,Geo2hInt(nR+nCut), &
                              GeoAs2hInt(nR+nCut),st_map)
            end do
            GeoRms  =rInt_R(Geo2hInt(nRC),n_r_maxC, &
                               n_cheb_maxC,dr_facC, &
                       i_costf_initC,d_costf_initC)
            GeoAsRms=rInt_R(GeoAs2hInt(nRC),n_r_maxC, &
                                 n_cheb_maxC,dr_facC, &
                         i_costf_initC,d_costf_initC)
            GeoRms  =dsqrt(GeoRms  /volC)
            GeoAsRms=dsqrt(GeoAsRms/volC)
    
            !------ Magnetostrophic balance:
            if ( .not. l_RMStest ) then
               call get_drNS(MagLMr(1,nRC),workA(1,nRC),lm_max,1, &
                    &        lm_max,n_r_maxC,n_cheb_maxC,workB,   &
                    &        i_costf_initC,d_costf_initC,dr_facC)
               do nR=1,n_r_maxC
                  call hInt2dPol(workA(1,nR+nCut),2,lm_max,Mag2hInt(nR+nCut), &
                                 MagAs2hInt(nR+nCut),st_map)
               end do
            end if
            MagRms  =rInt_R(Mag2hInt(nRC),n_r_maxC, &
                               n_cheb_maxC,dr_facC, &
                       i_costf_initC,d_costf_initC)
            MagAsRms=rInt_R(MagAs2hInt(nRC),n_r_maxC, &
                                 n_cheb_maxC,dr_facC, &
                         i_costf_initC,d_costf_initC)
            MagRms  =dsqrt(MagRms  /volC)
            MagAsRms=dsqrt(MagAsRms/volC)
    
            !------ Archemidian balance:
            call get_drNS(ArcLMr(1,nRC),workA(1,nRC),           &
                 &        lm_max,1,lm_max,n_r_maxC,n_cheb_maxC, &
                 &        workB,i_costf_initC,d_costf_initC,dr_facC)
            do nR=1,n_r_maxC
               call hInt2dPol( workA(1,nR+nCut),2,lm_max,Arc2hInt(nR+nCut), &
                               ArcAs2hInt(nR+nCut),st_map)
            end do
            ArcRms  =rInt_R(Arc2hInt(nRC),n_r_maxC, &
                               n_cheb_maxC,dr_facC, &
                       i_costf_initC,d_costf_initC)
            ArcAsRms=rInt_R(ArcAs2hInt(nRC),n_r_maxC, &
                                 n_cheb_maxC,dr_facC, &
                         i_costf_initC,d_costf_initC)
            if ( l_RMStest ) then
               ArcRms  =ArcRms/2.D0
               ArcAsRms=ArcAsRms/2.D0
            else
               ArcRms  =dsqrt(ArcRms  /volC)
               ArcAsRms=dsqrt(ArcAsRms/volC)
            end if
    
            !------ Diffusion:
            call get_drNS(DifPolLMr(1,nRC),workA(1,nRC),        &
                 &        lm_max,1,lm_max,n_r_maxC,n_cheb_maxC, &
                 &        workB,i_costf_initC,d_costf_initC,dr_facC)
            do nR=1,n_r_maxC
               call hInt2dPol( workA(1,nR+nCut),2,lm_max,DifPol2hInt(nR+nCut,1), &
                               DifPolAs2hInt(nR+nCut,1),lo_map)
               !do n=2,nThreadsLMmax
               !   DifPol2hInt(nR+nCut,1)  =DifPol2hInt(nR+nCut,1) + &
               !        DifPol2hInt(nR+nCut,n)
               !   DifPolAs2hInt(nR+nCut,1)=DifPolAs2hInt(nR+nCut,1) + &
               !        DifPolAs2hInt(nR+nCut,n)
               !   DifTor2hInt(nR+nCut,1)  =DifTor2hInt(nR+nCut,1) + &
               !        DifTor2hInt(nR+nCut,n)
               !   DifTorAs2hInt(nR+nCut,1)=DifTorAs2hInt(nR+nCut,1) + &
               !        DifTorAs2hInt(nR+nCut,n)
               !end do
            end do
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
            call get_drNS( dtVPolLMr(1,nRC),workA(1,nRC),lm_max,1,lm_max,&
                 &         n_r_maxC,n_cheb_maxC,workB,i_costf_initC,     &
                 &         d_costf_initC,dr_facC)
            do nR=1,n_r_maxC
               call hInt2dPol( workA(1,nR+nCut),2,lm_max,dtVPol2hInt(nR+nCut,1), &
                               dtVPolAs2hInt(nR+nCut,1),lo_map)
               !do n=2,nThreadsLMmax
               !   dtVPol2hInt(nR+nCut,1)  =dtVPol2hInt(nR+nCut,1) + &
               !        &                   dtVPol2hInt(nR+nCut,n)
               !   dtVPolAs2hInt(nR+nCut,1)=dtVPolAs2hInt(nR+nCut,1) + &
               !        &                   dtVPolAs2hInt(nR+nCut,n)
               !   dtVTor2hInt(nR+nCut,1)  =dtVTor2hInt(nR+nCut,1) + &
               !        &                   dtVTor2hInt(nR+nCut,n)
               !   dtVTorAs2hInt(nR+nCut,1)=dtVTorAs2hInt(nR+nCut,1) + &
               !        &                   dtVTorAs2hInt(nR+nCut,n)
               !end do
               !write(*,"(A,I3,ES22.14)") "after = ",nR+nCut,dtVPol2hInt(nR+nCut,1)
            end do
            !write(*,"(A,I3,ES22.14)") "dtVPol2hInt(nRC,1) = ",nRC,dtVPol2hInt(nRC,1)
            dtVPolRms  =rInt_R(dtVPol2hInt(nRC,1),n_r_maxC, &
                                       n_cheb_maxC,dr_facC, &
                               i_costf_initC,d_costf_initC)
            !write(*,"(A,ES22.14)") "dtVPolRms = ",dtVPolRms
            dtVPolAsRms=rInt_R(dtVPolAs2hInt(nRC,1),n_r_maxC, &
                                         n_cheb_maxC,dr_facC, &
                                 i_costf_initC,d_costf_initC)
            dtVTorRms  =rInt_R(dtVTor2hInt(nRC,1),n_r_maxC, &
                                       n_cheb_maxC,dr_facC, &
                               i_costf_initC,d_costf_initC)
            dtVTorAsRms=rInt_R(dtVTorAs2hInt(nRC,1),n_r_maxC, &
                                         n_cheb_maxC,dr_facC, &
                                i_costf_initC,d_costf_initC)
            dtVPolRms  =SQRT(dtVPolRms  /volC)
            dtVPolAsRms=SQRT(dtVPolAsRms/volC)
            dtVTorRms  =SQRT(dtVTorRms  /volC)
            dtVTorAsRms=SQRT(dtVTorAsRms/volC)
    
            !----- Output:
            if ( l_save_out) then
               open(n_dtvrms_file, file=dtvrms_file, form='formatted', &
                    status='unknown', position='append')
            end if
            if ( l_RMStest ) then
               write(n_dtvrms_file,'(1P,D20.12,15D16.8)')  &
                    time, dtVPolRms, dtVTorRms, CorPolRms, &
                    CorTorRms, LFPolRms, LFTorRms,         &
                    AdvPolRms, AdvTorRms, DifPolRms,       &
                    DifTorRms, BuoRms,PreRms, GeoRms,      &
                    MagRms, ArcRms
            else
               write(n_dtvrms_file,'(1P,D20.12,15D16.8)')    &
                    time, dtVPolRms, dtVTorRms, CorPolRms,   &
                    CorTorRms, LFPolRms,LFTorRms, AdvPolRms, &
                    AdvTorRms, DifPolRms, DifTorRms, BuoRms, &
                    PreRms, GeoRms/(CorPolRms+PreRms),       & 
                    MagRms/(CorPolRms+PreRms+LFPolRms),      &  
                    ArcRms/(CorPolRms+PreRms+LFPolRms+BuoRms) 
            end if
            if ( l_save_out) then
               close(n_dtvrms_file)
            end if
            if ( l_save_out) then
               open(n_dtvasrms_file, file=dtvasrms_file, form='formatted', &
                    status='unknown', position='append')
            end if
            if ( l_RMStest ) then
               write(n_dtvasrms_file,'(1P,D20.12,15D16.8)')      &
                    time, dtVPolAsRms, dtVTorAsRms, CorPolAsRms, &
                    CorTorAsRms, LFPolAsRms, LFTorAsRms,         &
                    AdvPolAsRms, AdvTorAsRms, DifPolAsRms,       &
                    DifTorAsRms, BuoAsRms, PreAsRms, GeoAsRms,   &
                    MagAsRms, ArcAsRms
            else
               if ( PreAsRms/= 0d0 ) then
                  bal1=GeoAsRms/(CorPolAsRms+PreAsRms)
                  bal2=MagAsRms/(CorPolAsRms+PreAsRms+LFPolAsRms)
                  bal3=ArcAsRms/(CorPolAsRms+PreAsRms+LFPolAsRms+BuoAsRms)
               else
                  bal1=0d0
                  bal2=0d0
                  bal3=0d0
               end if
               write(n_dtvasrms_file,'(1P,D20.12,15D16.8)')      &
                    time, dtVPolAsRms, dtVTorAsRms, CorPolAsRms, &
                    CorTorAsRms, LFPolAsRms, LFTorAsRms,         &
                    AdvPolAsRms,AdvTorAsRms, DifPolAsRms,        &
                    DifTorAsRms, BuoAsRms,PreAsRms, bal1, bal2,  &
                    bal3
            end if
            if ( l_save_out) then
               close(n_dtvasrms_file)
            end if
    
         end if
    
      end if

   end subroutine dtVrms
!----------------------------------------------------------------------------
   subroutine dtBrms(time)

      !-- Input of variables:
      real(kind=8), intent(in) :: time
    
      !-- Local
      integer :: nR,n,l1m0,lm
      character(len=80) :: fileName
    
      real(kind=8) :: dtBPolRms,dtBPolAsRms
      real(kind=8) :: dtBTorRms,dtBTorAsRms
      real(kind=8) :: DstrRms
      real(kind=8) :: DadvRms
      real(kind=8) :: DdifRms
      real(kind=8) :: DdynRms
      real(kind=8) :: PdynRms,PdynAsRms
      real(kind=8) :: TdynRms,TdynAsRms
      real(kind=8) :: dummy1,dummy2,dummy3
    
      !----- Note: five full additional fields needed here!
      !      This may cause memory problems.
      complex(kind=8) :: PdynLM(lm_max_dtB,n_r_max_dtB)
      complex(kind=8) :: drPdynLM(lm_max_dtB,n_r_max_dtB)
      complex(kind=8) :: TdynLM(lm_max_dtB,n_r_max_dtB)
      complex(kind=8) :: workA(lm_max_dtB,n_r_max_dtB)
      complex(kind=8) :: workB(lm_max_dtB,n_r_max_dtB)
    
      !-- For new movie output
      integer :: nField,nFields,nFieldSize
      integer :: nTheta,nThetaN,nThetaS,nThetaStart
      integer :: nPos
      real(kind=8) :: dumm(12),rS
      real(kind=8) :: fOut(n_theta_max*n_r_max)
      real(kind=8) :: outBlock(nfs)
      character(len=80) :: version
      logical :: lRmsMov
    
      real(kind=8) :: global_sum(n_r_max)
    
      call myAllGather(dtBPolLMr,lm_maxMag,n_r_maxMag)
      call mpi_reduce(dtBPol2hInt(1,1),global_sum,n_r_max,MPI_DOUBLE_PRECISION, &
                      MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) dtBPol2hInt(:,1)=global_sum
      call mpi_reduce(dtBPolAs2hInt(1,1),global_sum,n_r_max,MPI_DOUBLE_PRECISION, &
                      MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) dtBPolAs2hInt(:,1)=global_sum
    
      call mpi_reduce(dtBTor2hInt(1,1),global_sum,n_r_max,MPI_DOUBLE_PRECISION, &
                      MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) then
         dtBTor2hInt(:,1)=global_sum
      end if
    
    
      if ( rank == 0 ) then
         l1m0=lm2(1,0)
    
         !--- Stretching
         call get_drNS(PstrLM,workA,lm_max,1,lm_max,n_r_max, &
              &        n_cheb_max,workB,i_costf_init,d_costf_init,drx)
         !--- Finalize rms poloidal and toroidal stretching:
         call get_PolTorRms(PstrLM,workA,TstrLM,PstrRms,TstrRms,PstrAsRms, &
              &             TstrAsRms,st_map)
    
         !--- Calculate dipole stretching and copy for total dynamo term:
         do nR=1,n_r_max
            do lm=1,lm_max
               PdynLM(lm,nR)  =PstrLM(lm,nR)
               drPdynLM(lm,nR)=workA(lm,nR)
               if ( lm /= l1m0 ) then
                  PstrLM(lm,nR)=cmplx(0.D0,0.D0,kind=kind(0d0))
                  workA(lm,nR) =cmplx(0.D0,0.D0,kind=kind(0d0))
               end if
            end do
         end do
         !--- Get dipole stretching
         call get_PolTorRms(PstrLM,workA,TstrLM,DstrRms,dummy1,dummy2,dummy3,st_map)
    
         !--- Finalize advection
         call get_drNS(PadvLM,workA,lm_max,1,lm_max,n_r_max, &
              &        n_cheb_max,workB,i_costf_init,d_costf_init,drx)
         call get_PolTorRms(PadvLM,workA,TadvLM,PadvRms,TadvRms, &
              &             PadvAsRms,TadvAsRms,st_map)
         do nR=1,n_r_max
            do lm=1,lm_max
               PstrLM(lm,nR)  =PdynLM(lm,nR)
               PdynLM(lm,nR)  =PdynLM(lm,nR)-PadvLM(lm,nR)
               drPdynLM(lm,nR)=drPdynLM(lm,nR)+workA(lm,nR)
               TdynLM(lm,nR)  =TstrLM(lm,nR)-TadvLM(lm,nR)
               if ( lm /= l1m0 ) then
                  PadvLM(lm,nR)=cmplx(0.D0,0.D0,kind=kind(0d0))
                  workA(lm,nR) =cmplx(0.D0,0.D0,kind=kind(0d0))
               end if
            end do
         end do
         !--- Get dipole advection and total dynamo terms:
         call get_PolTorRms(PadvLM,workA,TadvLM,DadvRms, &
                            dummy1,dummy2,dummy3,st_map)
         call get_PolTorRms(PdynLM,drPdynLM,TdynLM,PdynRms, &
                            TdynRms,PdynAsRms,TdynAsRms,st_map)
         do nR=1,n_r_max
            do lm=1,lm_max
               PadvLM(lm,nR)=PdynLM(lm,nR)
               if ( lm /= l1m0 ) then
                  PdynLM(lm,nR)  =cmplx(0.D0,0.D0,kind=kind(0d0))
                  drPdynLM(lm,nR)=cmplx(0.D0,0.D0,kind=kind(0d0))
               end if
            end do
         end do
         !--- Get dipole dynamo terms:
         call get_PolTorRms(PdynLM,drPdynLM,TdynLM,DdynRms, &
                            dummy1,dummy2,dummy3,st_map)
    
         !--- Diffusion:
         call get_drNS(PdifLM,workA,lm_max,1,lm_max,n_r_max, &
              &        n_cheb_max,workB,i_costf_init,d_costf_init,drx)
         call get_PolTorRms(PdifLM,workA,TdifLM,PdifRms,TdifRms,&
                            PdifAsRms,TdifAsRms,st_map)
    
         do nR=1,n_r_max
            do lm=1,lm_max
               PdynLM(lm,nR)=PdifLM(lm,nR)
               if ( lm /= st_map%lm2(1,0) ) then
                  PdifLM(lm,nR)=cmplx(0.D0,0.D0,kind=kind(0d0))
                  workA(lm,nR) =cmplx(0.D0,0.D0,kind=kind(0d0))
               end if
            end do
         end do
         call get_PolTorRms(PdifLM,workA,TdifLM,DdifRms,dummy1,dummy2,dummy3,st_map)
    
    
         !--- Omega effect: (different mapping for PdifLM,workA and TomeLM)
         call get_PolTorRms(PdifLM,workA,TomeLM,dummy1,TomeRms,dummy2,TomeAsRms,st_map)
    
    
         !--- B changes:
         call get_drNS(dtBPolLMr,workA,lm_max,1,lm_max,n_r_max, &
              &        n_cheb_max,workB,i_costf_init,d_costf_init,drx)
         do nR=1,n_r_max
            call hInt2dPol(workA(1,nR),2,lm_max,dtBPol2hInt(nR,1), &
                           dtBPolAs2hInt(nR,1),lo_map)
         end do
         dtBPolRms  =rInt_R(dtBPol2hInt(1,1),n_r_max,   &
              &             n_r_max,drx,i_costf_init,d_costf_init)
         dtBPolAsRms=rInt_R(dtBPolAs2hInt(1,1),n_r_max, &
              &             n_r_max,drx,i_costf_init,d_costf_init)
         dtBTorRms  =rInt_R(dtBTor2hInt(1,1),n_r_max,   &
              &             n_r_max,drx,i_costf_init,d_costf_init)
         dtBTorAsRms=rInt_R(dtBTorAs2hInt(1,1),n_r_max, &
              &             n_r_max,drx,i_costf_init,d_costf_init)
         dtBPolRms  =sqrt(dtBPolRms  /vol_oc)
         dtBPolAsRms=sqrt(dtBPolAsRms/vol_oc)
         dtBTorRms  =sqrt(dtBTorRms  /vol_oc)
         dtBTorAsRms=sqrt(dtBTorAsRms/vol_oc)
    
         !-- Output:
         if ( l_save_out) then
            open(n_dtbrms_file, file=dtbrms_file, form='formatted', &
                 status='unknown', position='append')
         end if
         write(n_dtbrms_file,'(1P,D20.12,12D16.8)')              &
              time, dtBPolRms, dtBTorRms, PstrRms, TstrRms,      &     
              PadvRms,TadvRms, PdifRms,TdifRms, TomeRms/TstrRms, &
              TomeRms,  PdynRms,TdynRms 
         if ( l_save_out) then
            close(n_dtbrms_file)
         end if
    
         if ( l_save_out) then
            open(n_dtdrms_file, file=dtdrms_file, form='formatted', &
                 status='unknown', position='append')
         end if
         write(n_dtdrms_file,'(1P,D20.12,3D16.8)') &
              time, DstrRms, DadvRms, DdifRms
         if ( l_save_out) then
            close(n_dtdrms_file)
         end if
    
         !-- Output of movie files for axisymmetric toroidal field changes:
         !   Tstr,Tome,Tdyn=Tstr+Tadv,
         lRmsMov=.false.
         if ( lRmsMov ) then
    
            nFieldSize=n_theta_max*n_r_max
            nFields=7
            fileName='dtTas_mov.'//TAG
            open(90, file=fileName, status='unknown', form='unformatted')
    
            !------ Write header
            version='JW_Movie_Version_2'
            write(90) version
            dumm(1)=112           ! type of input
            dumm(2)=3             ! marker for constant phi plane
            dumm(3)=0.D0          ! surface constant
            dumm(4)=nFields       ! no of fields
            write(90) (real(dumm(n),4),n=1,4)
    
            !------ Define marker for output fields stored in movie field
            dumm(1)=101           ! Field marker for AS Br stretching
            dumm(2)=102           ! Field marker for AS Br dynamo term
            dumm(3)=103           ! Field marker for AS Br diffusion
            dumm(4)=104           ! Field marker for AS Bp stretching
            dumm(5)=105           ! Field marker for AS Bp dynamo term
            dumm(6)=106           ! Field marker for AS Bp omega effect
            dumm(7)=107           ! Field marker for AS Bp diffusion
            write(90) (sngl(dumm(n)),n=1,nFields)
    
            !------ Now other info about grid and parameters:
            write(90) runid        ! run identifier
            dumm( 1)=n_r_max          ! total number of radial points
            dumm( 2)=n_r_max          ! no of radial point in outer core
            dumm( 3)=n_theta_max      ! no. of theta points
            dumm( 4)=n_phi_max        ! no. of phi points
            dumm( 5)=minc             ! imposed symmetry
            dumm( 6)=ra               ! control parameters
            dumm( 7)=ek               ! (for information only)
            dumm( 8)=pr               !      -"-
            dumm( 9)=prmag            !      -"-
            dumm(10)=radratio         ! ratio of inner / outer core
            dumm(11)=tScale           ! timescale
            write(90) (sngl(dumm(n)),     n=1,11)
            write(90) (sngl(r(n)/r_CMB),  n=1,n_r_max)
            write(90) (sngl(theta_ord(n)),n=1,n_theta_max)
            write(90) (sngl(phi(n)),      n=1,n_phi_max)
    
            dumm(1)=1    ! time frame number for movie
            dumm(2)=0.D0 ! time
            dumm(3)=0.D0
            dumm(4)=0.D0
            dumm(5)=0.D0
            dumm(6)=0.D0
            dumm(7)=0.D0
            dumm(8)=0.D0
            write(90) (sngl(dumm(n)),n=1,8)
    
            !------ Loop over different output field:
            do nField=1,nFields
    
               !------ Loop over r and theta:
               do nR=1,n_r_max ! Loop over radial points
                  rS=r(nR)
                  do n=1,nThetaBs ! Loop over theta blocks
                     nThetaStart=(n-1)*sizeThetaB+1
    
                     !------ Convert from lm to theta block and store in outBlock:
                     if ( nField == 1 ) then
                        call get_RAS(PstrLM(1,nR),outBlock,rS,nThetaStart,sizeThetaB)
                     else if ( nField == 2 ) then
                        ! Note that PadvLM stores PdynLM=PstrLM+PadvLM at this point!
                        call get_RAS(PadvLM(1,nR),outBlock,rS,nThetaStart,sizeThetaB)
                     else if ( nField == 3 ) then
                        ! Note that PdynLM stores PdifLM at this point!
                        call get_RAS(PdynLM(1,nR),outBlock,rS,nThetaStart,sizeThetaB)
                     else if ( nField == 4 ) then
                        call get_PASLM(TstrLM(1,nR),outBlock,rS,nThetaStart,sizeThetaB)
                     else if ( nField == 5 ) then
                        call get_PASLM(TdynLM(1,nR),outBlock,rS,nThetaStart,sizeThetaB)
                     else if ( nField == 6 ) then
                        call get_PASLM(TomeLM(1,nR),outBlock,rS,nThetaStart,sizeThetaB)
                     else if ( nField == 7 ) then
                        call get_PASLM(TdifLM(1,nR),outBlock,rS,nThetaStart,sizeThetaB)
                     end if
    
                     !------ Storage of field in fout for theta block
                     do nTheta=1,sizeThetaB,2
                        !-- Convert to correct order in theta grid points 
                        !-- and store of fOut:
                        nThetaN=(nThetaStart+nTheta)/2
                        nPos=(nR-1)*n_theta_max+nThetaN
                        fOut(nPos)=outBlock(nTheta)
                        nThetaS=n_theta_max-nThetaN+1
                        nPos=(nR-1)*n_theta_max+nThetaS
                        fOut(nPos)=outBlock(nTheta+1)
                     end do ! Loop over thetas in block
    
                  end do ! Loop over theta blocks
    
               end do ! Loop over R
    
               !------ Output of field:
               write(90) (sngl(fOut(nPos)),nPos=1,nFieldSize)
    
            end do ! Loop over different fields
    
         end if ! output of mov fields ?

      end if

   end subroutine dtBrms
!----------------------------------------------------------------------------
end module out_RMS
