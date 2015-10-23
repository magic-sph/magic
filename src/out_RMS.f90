module out_RMS

   use parallel_mod
   use precision_mod
   use truncation, only: lm_max, n_r_max, lm_max_dtB, n_r_max_dtB, &
                         n_cheb_max, lm_maxMag, n_theta_max, minc, &
                         n_r_maxMag, n_phi_max
   use radial_data, only: nRstop, nRstart
   use radial_functions, only: chebt_oc, drx, r, r_CMB, rgrav
                               
   use physical_parameters, only: ra, ek, pr, prmag, radratio
   use blocking, only: st_map, nThetaBs, nfs, sizeThetaB, lo_map, lm2
   use logic, only: l_save_out, l_RMStest, l_heat, l_conv_nl, l_mag_LF, &
                    l_conv, l_corr
   use RMS

   use dtB_mod, only: PstrLM, TstrLM, PadvLM, TadvLM, TomeLM, PdifLM,  &
                      TdifLM, PstrRms,TstrRms, PstrAsRms, TstrAsRms,   &
                      PadvRms, TadvRms, PadvAsRms, TadvAsRms, PdifRms, &
                      TdifRms, PdifAsRms, TdifAsRms, TomeRms,          &
                      TomeAsRms
                                                                  
   use num_param, only: tScale
   use horizontal_data, only: phi, theta_ord
   use output_data, only: TAG, runid, n_dtdrms_file, dtdrms_file,        &
                          n_dtvasrms_file, dtvasrms_file, n_dtvrms_file, &
                          dtvrms_file, n_dtbrms_file, dtbrms_file
   use constants, only: pi, vol_oc, zero, half, four, third
   use integration, only: rInt_R
#ifdef WITH_MPI
   use communications, only: myAllGather
#endif
   use RMS_helpers, only: hInt2dPol, get_PolTorRms, get_PASLM, get_RAS
   use radial_der, only: get_drNS
   use cosine_transform_odd, only: costf_odd_t

   implicit none

   private

   public :: dtVrms, dtBrms

contains
 
   subroutine dtVrms(time,nRMS_sets)
      !
      !  For testing RMS balance and determining the necessary values of  
      !  rDea and rCut one can use the l_RMStest=true options.            
      !  In this case the poloidal and toroidal RMS dtV which also        
      !  include the diffusion effects should be identical (as close as   
      !  desired) to the RMS sum of forces stored in the Geo value and    
      !  the Mag value, respectively.                                     
      !  An additional tests is the Arc value which should be identical   
      !  to the poloidal kinetic energy.                                  
      !                                                                   
      !  .. note:: The second test with the Arc value cannot work so easily in 
      !            the anelastic version as the density enters now in the      
      !            force balance. The first test should work, though.         
      !

      !-- Input variable:
      real(cp), intent(in) :: time
      integer, intent(inout) :: nRMS_sets
    
      !-- Output:
      real(cp) :: dtVPolRms,dtVPolAsRms
      real(cp) :: dtVTorRms,dtVTorAsRms
      real(cp) :: CorPolRms,CorPolAsRms
      real(cp) :: CorTorRms,CorTorAsRms
      real(cp) :: AdvPolRms,AdvPolAsRms
      real(cp) :: AdvTorRms,AdvTorAsRms
      real(cp) :: LFPolRms, LFPolAsRms
      real(cp) :: LFTorRms, LFTorAsRms
      real(cp) :: DifPolRms,DifPolAsRms
      real(cp) :: DifTorRms,DifTorAsRms
      real(cp) :: BuoRms,   BuoAsRms
      real(cp) :: PreRms,   PreAsRms
      real(cp) :: GeoRms,   GeoAsRms
      real(cp) :: MagRms,   MagAsRms
      real(cp) :: ArcRms,   ArcAsRms
    
      !-- Local:
      integer :: nR,nRC
      !integer :: n
      real(cp) :: volC
      real(cp) :: bal1,bal2,bal3
    
      complex(cp) :: workA(lm_max,n_r_max),workB(lm_max,n_r_max)
      integer :: recvcounts(0:n_procs-1),displs(0:n_procs-1)
      real(cp) :: global_sum(n_r_max)
      integer :: irank,sendcount
    
      ! First gather all needed arrays on rank 0
      ! some more arrays to gather for the dtVrms routine
      ! we need some more fields for the dtBrms routine
#ifdef WITH_MPI
      sendcount  = (nRstop-nRstart+1)*lm_max
      recvcounts = nr_per_rank*lm_max
      recvcounts(n_procs-1) = (nr_per_rank+1)*lm_max
      do irank=0,n_procs-1
         displs(irank) = irank*nr_per_rank*lm_max
      end do
      ! CorPolLMr,dtVPolLMr,AdvPolLMr,LFPolLMr,DifPolLMr,BuoLMr
      ! PreLMr,GeoLMr,MagLMr,ArcLMr
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_COMPLEX,&
           & CorPolLMr,recvcounts,displs,MPI_DEF_COMPLEX,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_COMPLEX,&
           & AdvPolLMr,recvcounts,displs,MPI_DEF_COMPLEX,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_COMPLEX,&
           & LFPolLMr,recvcounts,displs,MPI_DEF_COMPLEX,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_COMPLEX,&
           & BuoLMr,recvcounts,displs,MPI_DEF_COMPLEX,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_COMPLEX,&
           & PreLMr,recvcounts,displs,MPI_DEF_COMPLEX,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_COMPLEX,&
           & GeoLMr,recvcounts,displs,MPI_DEF_COMPLEX,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_COMPLEX,&
           & MagLMr,recvcounts,displs,MPI_DEF_COMPLEX,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_COMPLEX,&
           & ArcLMr,recvcounts,displs,MPI_DEF_COMPLEX,MPI_COMM_WORLD,ierr)
    
      ! The following fields are only 1D and R distributed.
      sendcount  = (nRstop-nRstart+1)
      recvcounts = nr_per_rank
      recvcounts(n_procs-1) = (nr_per_rank+1)
      do irank=0,n_procs-1
         displs(irank) = irank*nr_per_rank
      end do
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & CorPol2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & CorPolAs2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & CorTor2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & CorTorAs2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & AdvPol2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & AdvPolAs2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & AdvTor2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & AdvTorAs2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & LFPol2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & LFPolAs2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & LFTor2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & LFTorAs2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & Buo2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & BuoAs2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & Pre2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & PreAs2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & Geo2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & GeoAs2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & Mag2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & MagAs2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & Arc2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & ArcAs2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
    
      ! The following fields are LM distributed and have to be gathered:
      ! dtVPolLMr, DifPolLMr
    
      call myAllGather(dtVPolLMr,lm_max,n_r_max)
      call myAllGather(DifPolLMr,lm_max,n_r_max)
    
      call MPI_Reduce(dtVPol2hInt(1,1),global_sum,n_r_max,MPI_DEF_REAL, &
           &          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) dtVPol2hInt(:,1)=global_sum
      call MPI_Reduce(dtVPolAs2hInt(1,1),global_sum,n_r_max,MPI_DEF_REAL,&
           &          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) dtVPolAs2hInt(:,1)=global_sum
      call MPI_Reduce(dtVTor2hInt(1,1),global_sum,n_r_max,MPI_DEF_REAL, &
           &          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) dtVTor2hInt(:,1)=global_sum
      call MPI_Reduce(dtVTorAs2hInt(1,1),global_sum,n_r_max,MPI_DEF_REAL, &
           &          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) dtVTorAs2hInt(:,1)=global_sum
      call MPI_Reduce(DifPol2hInt(1,1),global_sum,n_r_max,MPI_DEF_REAL, &
           &          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) DifPol2hInt(:,1)=global_sum
      call MPI_Reduce(DifPolAs2hInt(1,1),global_sum,n_r_max,MPI_DEF_REAL, &
           &          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) DifPolAs2hInt(:,1)=global_sum
      call MPI_Reduce(DifTor2hInt(1,1),global_sum,n_r_max,MPI_DEF_REAL, &
           &          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) DifTor2hInt(:,1)=global_sum
      call MPI_Reduce(DifTorAs2hInt(1,1),global_sum,n_r_max,MPI_DEF_REAL,&
           &          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) DifTorAs2hInt(:,1)=global_sum
#endif
    
      if ( rank == 0 ) then
    
         !write(*,"(A,ES22.14)") "dtVPol2hInt = ",SUM(dtVPol2hInt)
         if ( nRMS_sets == 0 ) then
            nRMS_sets=nRMS_sets+1
         end if
         nRC=nCut+1
         volC=four*third*pi*(r(1+nCut)**3-r(n_r_max-nCut)**3)
    
         if ( l_conv ) then
    
            !------ Coriolis force
            if ( l_corr ) then
               call get_drNS( CorPolLMr(1,nRC),workA(1,nRC),        &
                    &           lm_max,1,lm_max,n_r_maxC,n_cheb_maxC, &
                    &           workB,chebt_RMS,dr_facC)
               do nR=1,n_r_maxC
                  call hInt2dPol(workA(1,nR+nCut),2,lm_max,CorPol2hInt(nR+nCut), &
                                 CorPolAs2hInt(nR+nCut),st_map)
               end do
               CorPolRms=rInt_R(CorPol2hInt(nRC),n_r_maxC, &
                                      n_cheb_maxC,dr_facC, &
                              chebt_RMS)
               CorPolAsRms=rInt_R(CorPolAs2hInt(nRC),n_r_maxC, &
                                          n_cheb_maxC,dr_facC, &
                                  chebt_RMS)
               CorTorRms  =rInt_R(CorTor2hInt(nRC),n_r_maxC, &
                                        n_cheb_maxC,dr_facC, &
                                chebt_RMS)
               CorTorAsRms=rInt_R(CorTorAs2hInt(nRC),n_r_maxC, &
                                          n_cheb_maxC,dr_facC, &
                                  chebt_RMS)
               CorPolRms  =sqrt(CorPolRms  /volC)
               CorPolAsRms=sqrt(CorPolAsRms/volC)
               CorTorRms  =sqrt(CorTorRms  /volC)
               CorTorAsRms=sqrt(CorTorAsRms/volC)
            end if
    
            !------ Advection:
            if ( l_conv_nl ) then
               call get_drNS(AdvPolLMr(1,nRC),workA(1,nRC),          &
                    &          lm_max,1,lm_max,n_r_maxC,n_cheb_maxC,   &
                    &          workB,chebt_RMS,dr_facC)
               do nR=1,n_r_maxC
                  call hInt2dPol(workA(1,nR+nCut),2,lm_max,AdvPol2hInt(nR+nCut), &
                                 AdvPolAs2hInt(nR+nCut),st_map)
               end do
               AdvPolRms  =rInt_R(AdvPol2hInt(nRC),n_r_maxC, &
                                        n_cheb_maxC,dr_facC, &
                                chebt_RMS)
               AdvPolAsRms=rInt_R(AdvPolAs2hInt(nRC),n_r_maxC, &
                                          n_cheb_maxC,dr_facC, &
                                  chebt_RMS)
               AdvTorRms  =rInt_R(AdvTor2hInt(nRC),n_r_maxC, &
                                        n_cheb_maxC,dr_facC, &
                                chebt_RMS)
               AdvTorAsRms=rInt_R(AdvTorAs2hInt(nRC),n_r_maxC, &
                                          n_cheb_maxC,dr_facC, &
                                  chebt_RMS)
               AdvPolRms  =sqrt(AdvPolRms  /volC)
               AdvPolAsRms=sqrt(AdvPolAsRms/volC)
               AdvTorRms  =sqrt(AdvTorRms  /volC)
               AdvTorAsRms=sqrt(AdvTorAsRms/volC)
            end if
    
            !------ Lorentz force:
            if ( l_mag_LF ) then
               call get_drNS(LFPolLMr(1,nRC),workA(1,nRC),lm_max,1,lm_max,&
                    &        n_r_maxC,n_cheb_maxC,workB,chebt_RMS,dr_facC)
               do nR=1,n_r_maxC
                  call hInt2dPol( workA(1,nR+nCut),2,lm_max,LFPol2hInt(nR+nCut), &
                                  LFPolAs2hInt(nR+nCut),st_map)
               end do
               LFPolRms  =rInt_R(LFPol2hInt(nRC),n_r_maxC, &
                                      n_cheb_maxC,dr_facC, &
                              chebt_RMS)
               LFPolAsRms=rInt_R(LFPolAs2hInt(nRC),n_r_maxC, &
                                        n_cheb_maxC,dr_facC, &
                                chebt_RMS)
               LFTorRms  =rInt_R(LFTor2hInt(nRC),n_r_maxC, &
                                      n_cheb_maxC,dr_facC, &
                              chebt_RMS)
               LFTorAsRms=rInt_R(LFTorAs2hInt(nRC),n_r_maxC, &
                                        n_cheb_maxC,dr_facC, &
                                chebt_RMS)
               LFPolRms  =sqrt(LFPolRms  /volC)
               LFPolAsRms=sqrt(LFPolAsRms/volC)
               LFTorRms  =sqrt(LFTorRms  /volC)
               LFTorAsRms=sqrt(LFTorAsRms/volC)
            else
               LFPolRms  =0.0_cp
               LFPolAsRms=0.0_cp
               LFTorRms  =0.0_cp
               LFTorAsRms=0.0_cp
            end if
    
            !------ Buoyancy:
            if ( l_heat ) then
               call get_drNS(BuoLMr(1,nRC),workA(1,nRC),lm_max,1,lm_max, &
                    &        n_r_maxC,n_cheb_maxC,workB,chebt_RMS,dr_facC)
               do nR=1,n_r_maxC
                  call hInt2dPol(workA(1,nR+nCut),2,lm_max, Buo2hInt(nR+nCut), &
                                 BuoAs2hInt(nR+nCut),st_map)
                  Buo2hInt(nR)  =Buo2hInt(nR)*rgrav(nR)**2
                  BuoAs2hInt(nR)=BuoAs2hInt(nR)*rgrav(nR)**2
               end do
               BuoRms  =rInt_R(Buo2hInt(nRC),n_r_maxC, &
                                  n_cheb_maxC,dr_facC, &
                          chebt_RMS)
               BuoAsRms=rInt_R(BuoAs2hInt(nRC),n_r_maxC, &
                                    n_cheb_maxC,dr_facC, &
                            chebt_RMS)
               BuoRms  =sqrt(BuoRms  /volC)
               BuoAsRms=sqrt(BuoAsRms/volC)
            else
               BuoRms=0.0_cp
               BuoAsRms=0.0_cp
            end if
    
            !------ Pressure gradient:
            call get_drNS(PreLMr(1,nRC),workA(1,nRC),lm_max,1,  &
                 &        lm_max,n_r_maxC,n_cheb_maxC,workB,    &
                 &        chebt_RMS,dr_facC)
            do nR=1,n_r_maxC
               call hInt2dPol(workA(1,nR+nCut),2,lm_max,Pre2hInt(nR+nCut), &
                             PreAs2hInt(nR+nCut),st_map)
            end do
            PreRms  =rInt_R(Pre2hInt(nRC),n_r_maxC, &
                               n_cheb_maxC,dr_facC, &
                       chebt_RMS)
            PreAsRms=rInt_R(PreAs2hInt(nRC),n_r_maxC, &
                                 n_cheb_maxC,dr_facC, &
                         chebt_RMS)
            PreRms  =sqrt(PreRms  /volC)
            PreAsRms=sqrt(PreAsRms/volC)
    
            !------ Geostrophic balance:
            call get_drNS(GeoLMr(1,nRC),workA(1,nRC),lm_max,1, &
                 &        lm_max,n_r_maxC,n_cheb_maxC,workB,   &
                 &        chebt_RMS,dr_facC)
            do nR=1,n_r_maxC
               call hInt2dPol(workA(1,nR+nCut),2,lm_max,Geo2hInt(nR+nCut), &
                              GeoAs2hInt(nR+nCut),st_map)
            end do
            GeoRms  =rInt_R(Geo2hInt(nRC),n_r_maxC, &
                               n_cheb_maxC,dr_facC, &
                       chebt_RMS)
            GeoAsRms=rInt_R(GeoAs2hInt(nRC),n_r_maxC, &
                                 n_cheb_maxC,dr_facC, &
                         chebt_RMS)
            GeoRms  =sqrt(GeoRms  /volC)
            GeoAsRms=sqrt(GeoAsRms/volC)
    
            !------ Magnetostrophic balance:
            if ( .not. l_RMStest ) then
               call get_drNS(MagLMr(1,nRC),workA(1,nRC),lm_max,1, &
                    &        lm_max,n_r_maxC,n_cheb_maxC,workB,   &
                    &        chebt_RMS,dr_facC)
               do nR=1,n_r_maxC
                  call hInt2dPol(workA(1,nR+nCut),2,lm_max,Mag2hInt(nR+nCut), &
                                 MagAs2hInt(nR+nCut),st_map)
               end do
            end if
            MagRms  =rInt_R(Mag2hInt(nRC),n_r_maxC, &
                               n_cheb_maxC,dr_facC, &
                       chebt_RMS)
            MagAsRms=rInt_R(MagAs2hInt(nRC),n_r_maxC, &
                                 n_cheb_maxC,dr_facC, &
                         chebt_RMS)
            MagRms  =sqrt(MagRms  /volC)
            MagAsRms=sqrt(MagAsRms/volC)
    
            !------ Archemidian balance:
            call get_drNS(ArcLMr(1,nRC),workA(1,nRC),           &
                 &        lm_max,1,lm_max,n_r_maxC,n_cheb_maxC, &
                 &        workB,chebt_RMS,dr_facC)
            do nR=1,n_r_maxC
               call hInt2dPol( workA(1,nR+nCut),2,lm_max,Arc2hInt(nR+nCut), &
                               ArcAs2hInt(nR+nCut),st_map)
            end do
            ArcRms  =rInt_R(Arc2hInt(nRC),n_r_maxC, &
                               n_cheb_maxC,dr_facC, &
                       chebt_RMS)
            ArcAsRms=rInt_R(ArcAs2hInt(nRC),n_r_maxC, &
                                 n_cheb_maxC,dr_facC, &
                         chebt_RMS)
            if ( l_RMStest ) then
               ArcRms  =half*ArcRms
               ArcAsRms=half*ArcAsRms
            else
               ArcRms  =sqrt(ArcRms  /volC)
               ArcAsRms=sqrt(ArcAsRms/volC)
            end if
    
            !------ Diffusion:
            call get_drNS(DifPolLMr(1,nRC),workA(1,nRC),        &
                 &        lm_max,1,lm_max,n_r_maxC,n_cheb_maxC, &
                 &        workB,chebt_RMS,dr_facC)
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
                               chebt_RMS)
            DifPolAsRms=rInt_R(DifPolAs2hInt(nRC,1),n_r_maxC, &
                                         n_cheb_maxC,dr_facC, &
                                 chebt_RMS)
            DifTorRms  =rInt_R(DifTor2hInt(nRC,1),n_r_maxC, &
                                       n_cheb_maxC,dr_facC, &
                               chebt_RMS)
            DifTorAsRms=rInt_R(DifTorAs2hInt(nRC,1),n_r_maxC, &
                                         n_cheb_maxC,dr_facC, &
                                 chebt_RMS)
            DifPolRms  =sqrt(DifPolRms  /volC)
            DifPolAsRms=sqrt(DifPolAsRms/volC)
            DifTorRms  =sqrt(DifTorRms  /volC)
            DifTorAsRms=sqrt(DifTorAsRms/volC)
    
            !------ Flow changes: Inertia - Advection
            call get_drNS( dtVPolLMr(1,nRC),workA(1,nRC),lm_max,1,lm_max,&
                 &         n_r_maxC,n_cheb_maxC,workB,chebt_RMS,dr_facC)
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
                               chebt_RMS)
            !write(*,"(A,ES22.14)") "dtVPolRms = ",dtVPolRms
            dtVPolAsRms=rInt_R(dtVPolAs2hInt(nRC,1),n_r_maxC, &
                                         n_cheb_maxC,dr_facC, &
                                 chebt_RMS)
            dtVTorRms  =rInt_R(dtVTor2hInt(nRC,1),n_r_maxC, &
                                       n_cheb_maxC,dr_facC, &
                               chebt_RMS)
            dtVTorAsRms=rInt_R(dtVTorAs2hInt(nRC,1),n_r_maxC, &
                                         n_cheb_maxC,dr_facC, &
                                chebt_RMS)
            dtVPolRms  =sqrt(dtVPolRms  /volC)
            dtVPolAsRms=sqrt(dtVPolAsRms/volC)
            dtVTorRms  =sqrt(dtVTorRms  /volC)
            dtVTorAsRms=sqrt(dtVTorAsRms/volC)
    
            !----- Output:
            if ( l_save_out) then
               open(n_dtvrms_file, file=dtvrms_file, form='formatted', &
                    status='unknown', position='append')
            end if
            if ( l_RMStest ) then
               write(n_dtvrms_file,'(1P,ES20.12,15ES16.8)')  &
                    time, dtVPolRms, dtVTorRms, CorPolRms,   &
                    CorTorRms, LFPolRms, LFTorRms,           &
                    AdvPolRms, AdvTorRms, DifPolRms,         &
                    DifTorRms, BuoRms,PreRms, GeoRms,        &
                    MagRms, ArcRms
            else
               write(n_dtvrms_file,'(1P,ES20.12,15ES16.8)')  &
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
               write(n_dtvasrms_file,'(1P,ES20.12,15ES16.8)')    &
                    time, dtVPolAsRms, dtVTorAsRms, CorPolAsRms, &
                    CorTorAsRms, LFPolAsRms, LFTorAsRms,         &
                    AdvPolAsRms, AdvTorAsRms, DifPolAsRms,       &
                    DifTorAsRms, BuoAsRms, PreAsRms, GeoAsRms,   &
                    MagAsRms, ArcAsRms
            else
               if ( PreAsRms/= 0.0_cp ) then
                  bal1=GeoAsRms/(CorPolAsRms+PreAsRms)
                  bal2=MagAsRms/(CorPolAsRms+PreAsRms+LFPolAsRms)
                  bal3=ArcAsRms/(CorPolAsRms+PreAsRms+LFPolAsRms+BuoAsRms)
               else
                  bal1=0.0_cp
                  bal2=0.0_cp
                  bal3=0.0_cp
               end if
               write(n_dtvasrms_file,'(1P,ES20.12,15ES16.8)')    &
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
      real(cp), intent(in) :: time
    
      !-- Local
      integer :: nR,n,l1m0,lm
      character(len=80) :: fileName
    
      real(cp) :: dtBPolRms,dtBPolAsRms
      real(cp) :: dtBTorRms,dtBTorAsRms
      real(cp) :: DstrRms
      real(cp) :: DadvRms
      real(cp) :: DdifRms
      real(cp) :: DdynRms
      real(cp) :: PdynRms,PdynAsRms
      real(cp) :: TdynRms,TdynAsRms
      real(cp) :: dummy1,dummy2,dummy3
    
      !----- Note: five full additional fields needed here!
      !      This may cause memory problems.
      complex(cp) :: PdynLM(lm_max_dtB,n_r_max_dtB)
      complex(cp) :: drPdynLM(lm_max_dtB,n_r_max_dtB)
      complex(cp) :: TdynLM(lm_max_dtB,n_r_max_dtB)
      complex(cp) :: workA(lm_max_dtB,n_r_max_dtB)
      complex(cp) :: workB(lm_max_dtB,n_r_max_dtB)
    
      !-- For new movie output
      integer :: nField,nFields,nFieldSize
      integer :: nTheta,nThetaN,nThetaS,nThetaStart
      integer :: nPos
      real(cp) :: dumm(12),rS
      real(cp) :: fOut(n_theta_max*n_r_max)
      real(cp) :: outBlock(nfs)
      character(len=80) :: version
      logical :: lRmsMov
    
      real(cp) :: global_sum(n_r_max)
    
#ifdef WITH_MPI
      call myAllGather(dtBPolLMr,lm_maxMag,n_r_maxMag)
      call MPI_Reduce(dtBPol2hInt(1,1),global_sum,n_r_max,MPI_DEF_REAL, &
                      MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) dtBPol2hInt(:,1)=global_sum
      call MPI_Reduce(dtBPolAs2hInt(1,1),global_sum,n_r_max,MPI_DEF_REAL, &
                      MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) dtBPolAs2hInt(:,1)=global_sum
    
      call MPI_Reduce(dtBTor2hInt(1,1),global_sum,n_r_max,MPI_DEF_REAL, &
                      MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) then
         dtBTor2hInt(:,1)=global_sum
      end if
#endif
    
    
      if ( rank == 0 ) then
         l1m0=lm2(1,0)
    
         !--- Stretching
         call get_drNS(PstrLM,workA,lm_max,1,lm_max,n_r_max, &
              &        n_cheb_max,workB,chebt_oc,drx)
         !--- Finalize rms poloidal and toroidal stretching:
         call get_PolTorRms(PstrLM,workA,TstrLM,PstrRms,TstrRms,PstrAsRms, &
              &             TstrAsRms,st_map)
    
         !--- Calculate dipole stretching and copy for total dynamo term:
         do nR=1,n_r_max
            do lm=1,lm_max
               PdynLM(lm,nR)  =PstrLM(lm,nR)
               drPdynLM(lm,nR)=workA(lm,nR)
               if ( lm /= l1m0 ) then
                  PstrLM(lm,nR)=zero
                  workA(lm,nR) =zero
               end if
            end do
         end do
         !--- Get dipole stretching
         call get_PolTorRms(PstrLM,workA,TstrLM,DstrRms,dummy1,dummy2,dummy3,st_map)
    
         !--- Finalize advection
         call get_drNS(PadvLM,workA,lm_max,1,lm_max,n_r_max, &
              &        n_cheb_max,workB,chebt_oc,drx)
         call get_PolTorRms(PadvLM,workA,TadvLM,PadvRms,TadvRms, &
              &             PadvAsRms,TadvAsRms,st_map)
         do nR=1,n_r_max
            do lm=1,lm_max
               PstrLM(lm,nR)  =PdynLM(lm,nR)
               PdynLM(lm,nR)  =PdynLM(lm,nR)-PadvLM(lm,nR)
               drPdynLM(lm,nR)=drPdynLM(lm,nR)+workA(lm,nR)
               TdynLM(lm,nR)  =TstrLM(lm,nR)-TadvLM(lm,nR)
               if ( lm /= l1m0 ) then
                  PadvLM(lm,nR)=zero
                  workA(lm,nR) =zero
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
                  PdynLM(lm,nR)  =zero
                  drPdynLM(lm,nR)=zero
               end if
            end do
         end do
         !--- Get dipole dynamo terms:
         call get_PolTorRms(PdynLM,drPdynLM,TdynLM,DdynRms, &
                            dummy1,dummy2,dummy3,st_map)
    
         !--- Diffusion:
         call get_drNS(PdifLM,workA,lm_max,1,lm_max,n_r_max, &
              &        n_cheb_max,workB,chebt_oc,drx)
         call get_PolTorRms(PdifLM,workA,TdifLM,PdifRms,TdifRms,&
                            PdifAsRms,TdifAsRms,st_map)
    
         do nR=1,n_r_max
            do lm=1,lm_max
               PdynLM(lm,nR)=PdifLM(lm,nR)
               if ( lm /= st_map%lm2(1,0) ) then
                  PdifLM(lm,nR)=zero
                  workA(lm,nR) =zero
               end if
            end do
         end do
         call get_PolTorRms(PdifLM,workA,TdifLM,DdifRms,dummy1,dummy2,dummy3,st_map)
    
    
         !--- Omega effect: (different mapping for PdifLM,workA and TomeLM)
         call get_PolTorRms(PdifLM,workA,TomeLM,dummy1,TomeRms,dummy2,TomeAsRms,st_map)
    
    
         !--- B changes:
         call get_drNS(dtBPolLMr,workA,lm_max,1,lm_max,n_r_max, &
              &        n_cheb_max,workB,chebt_oc,drx)
         do nR=1,n_r_max
            call hInt2dPol(workA(1,nR),2,lm_max,dtBPol2hInt(nR,1), &
                           dtBPolAs2hInt(nR,1),lo_map)
         end do
         dtBPolRms  =rInt_R(dtBPol2hInt(1,1),n_r_max,   &
              &             n_r_max,drx,chebt_oc)
         dtBPolAsRms=rInt_R(dtBPolAs2hInt(1,1),n_r_max, &
              &             n_r_max,drx,chebt_oc)
         dtBTorRms  =rInt_R(dtBTor2hInt(1,1),n_r_max,   &
              &             n_r_max,drx,chebt_oc)
         dtBTorAsRms=rInt_R(dtBTorAs2hInt(1,1),n_r_max, &
              &             n_r_max,drx,chebt_oc)
         dtBPolRms  =sqrt(dtBPolRms  /vol_oc)
         dtBPolAsRms=sqrt(dtBPolAsRms/vol_oc)
         dtBTorRms  =sqrt(dtBTorRms  /vol_oc)
         dtBTorAsRms=sqrt(dtBTorAsRms/vol_oc)
    
         !-- Output:
         if ( l_save_out) then
            open(n_dtbrms_file, file=dtbrms_file, form='formatted', &
                 status='unknown', position='append')
         end if
         write(n_dtbrms_file,'(1P,ES20.12,12ES16.8)')            &
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
         write(n_dtdrms_file,'(1P,ES20.12,3ES16.8)') &
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
            dumm(3)=0.0_cp          ! surface constant
            dumm(4)=nFields       ! no of fields
            write(90) (real(dumm(n),kind=outp),n=1,4)
    
            !------ Define marker for output fields stored in movie field
            dumm(1)=101           ! Field marker for AS Br stretching
            dumm(2)=102           ! Field marker for AS Br dynamo term
            dumm(3)=103           ! Field marker for AS Br diffusion
            dumm(4)=104           ! Field marker for AS Bp stretching
            dumm(5)=105           ! Field marker for AS Bp dynamo term
            dumm(6)=106           ! Field marker for AS Bp omega effect
            dumm(7)=107           ! Field marker for AS Bp diffusion
            write(90) (real(dumm(n),kind=outp),n=1,nFields)
    
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
            write(90) (real(dumm(n),kind=outp),     n=1,11)
            write(90) (real(r(n)/r_CMB,kind=outp),  n=1,n_r_max)
            write(90) (real(theta_ord(n),kind=outp),n=1,n_theta_max)
            write(90) (real(phi(n),kind=outp),      n=1,n_phi_max)
    
            dumm(1)=1    ! time frame number for movie
            dumm(2)=0.0_cp ! time
            dumm(3)=0.0_cp
            dumm(4)=0.0_cp
            dumm(5)=0.0_cp
            dumm(6)=0.0_cp
            dumm(7)=0.0_cp
            dumm(8)=0.0_cp
            write(90) (real(dumm(n),kind=outp),n=1,8)
    
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
               write(90) (real(fOut(nPos),kind=outp),nPos=1,nFieldSize)
    
            end do ! Loop over different fields
    
         end if ! output of mov fields ?

      end if

   end subroutine dtBrms
!----------------------------------------------------------------------------
end module out_RMS
