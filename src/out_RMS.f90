module out_RMS

   use parallel_mod
   use precision_mod
   use truncation, only: lm_max, n_r_max, lm_max_dtB, n_r_max_dtB, &
                         n_cheb_max, lm_maxMag, n_theta_max, minc, &
                         n_r_maxMag, n_phi_max, l_max
   use radial_data, only: nRstop, nRstart
   use radial_functions, only: chebt_oc, drx, r, r_CMB, rgrav, dr_fac
   use physical_parameters, only: ra, ek, pr, prmag, radratio
   use blocking, only: st_map, nThetaBs, nfs, sizeThetaB, lo_map, lm2, &
                       lm2m
   use logic, only: l_save_out, l_heat, l_conv_nl, l_mag_LF, l_conv, l_corr
   use RMS

   use dtB_mod, only: PstrLM, TstrLM, PadvLM, TadvLM, TomeLM, PdifLM,  &
                      TdifLM, PstrRms,TstrRms, PstrAsRms, TstrAsRms,   &
                      PadvRms, TadvRms, PadvAsRms, TadvAsRms, PdifRms, &
                      TdifRms, PdifAsRms, TdifAsRms, TomeRms,          &
                      TomeAsRms
                                                                  
   use num_param, only: tScale
   use horizontal_data, only: phi, theta_ord
   use output_data, only: TAG, runid, n_dtvrms_file, dtvrms_file, &
                          n_dtbrms_file, dtbrms_file
   use constants, only: pi, vol_oc, zero, half, four, third
   use integration, only: rInt_R, rInt
#ifdef WITH_MPI
   use communications, only: myAllGather
#endif
   use RMS_helpers, only: hInt2dPol, get_PolTorRms, get_PASLM, get_RAS, &
                          hInt2dPolLM
   use radial_der, only: get_drNS
   use cosine_transform_odd, only: costf_odd_t

   implicit none

   private

   public :: dtVrms, dtBrms

contains
 
   subroutine dtVrms(time,nRMS_sets)
      !
      ! This routine calculates and stores the different contributions
      ! of the forces entering the Navier-Stokes equation.
      !

      !-- Input variable:
      real(cp), intent(in) :: time
      integer, intent(inout) :: nRMS_sets
    
      !-- Output:
      real(cp) :: dtV_Rms,dtVRmsL
      real(cp) :: CorRms,CorRmsL
      real(cp) :: AdvRms,AdvRmsL
      real(cp) :: LFRms,LFRmsL
      real(cp) :: DifRms,DifRmsL
      real(cp) :: BuoRms,BuoRmsL
      real(cp) :: PreRms,PreRmsL
      real(cp) :: GeoRms,GeoRmsL
      real(cp) :: MagRms,MagRmsL
      real(cp) :: ArcRms,ArcRmsL
      real(cp) :: CLFRms,CLFRmsL
      real(cp) :: PLFRms,PLFRmsL
    
      !-- Local:
      integer :: nR,nRC,l,n
      real(cp) :: volC
      real(cp) :: Rms(n_r_max),Dif2hInt(n_r_max),dtV2hInt(n_r_max)
    
      complex(cp) :: workA(lm_max,n_r_max),workB(lm_max,n_r_max)
      integer :: recvcounts(0:n_procs-1),displs(0:n_procs-1)
      real(cp) :: global_sum(l_max+1,n_r_max)
      integer :: irank,sendcount
    
      ! First gather all needed arrays on rank 0
      ! some more arrays to gather for the dtVrms routine
      ! we need some more fields for the dtBrms routine
#ifdef WITH_MPI
    
      ! The following fields are only 1D and R distributed.
      sendcount  = (nRstop-nRstart+1)*(l_max+1)
      recvcounts = nr_per_rank*(l_max+1)
      recvcounts(n_procs-1) = (nr_per_rank+1)*(l_max+1)
      do irank=0,n_procs-1
         displs(irank) = irank*nr_per_rank*(l_max+1)
      end do
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & Cor2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & Adv2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & LF2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & Buo2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & Pre2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & Geo2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & Mag2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & Arc2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & CLF2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(MPI_IN_PLACE,sendcount,MPI_DEF_REAL,&
           & PLF2hInt,recvcounts,displs,MPI_DEF_REAL,MPI_COMM_WORLD,ierr)
    
      ! The following fields are LM distributed and have to be gathered:
      ! dtVPolLMr, DifPolLMr
    
      call myAllGather(dtVPolLMr,lm_max,n_r_max)
      call myAllGather(DifPolLMr,lm_max,n_r_max)
    
      call MPI_Reduce(dtVPol2hInt(:,:,1),global_sum,n_r_max*(l_max+1), &
           &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) dtVPol2hInt(:,:,1)=global_sum
      call MPI_Reduce(dtVTor2hInt(:,:,1),global_sum,n_r_max*(l_max+1), &
           &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) dtVTor2hInt(:,:,1)=global_sum
      call MPI_Reduce(DifPol2hInt(:,:,1),global_sum,n_r_max*(l_max+1), &
           &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) DifPol2hInt(:,:,1)=global_sum
      call MPI_Reduce(DifTor2hInt(:,:,1),global_sum,n_r_max*(l_max+1), &
           &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) DifTor2hInt(:,:,1)=global_sum
#endif
    
      if ( rank == 0 ) then
    
         !write(*,"(A,ES22.14)") "dtVPol2hInt = ",SUM(dtVPol2hInt)
         if ( nRMS_sets == 0 ) then
            nRMS_sets=nRMS_sets+1
         end if
         nRC=nCut+1
         volC=four*third*pi*(r(1+nCut)**3-r(n_r_max-nCut)**3)
    
         CorRms=0.0_cp
         if ( l_corr ) then
            do l=0,l_max
               !-- Copy each mode on radial array
               do nR=1,n_r_max
                  Rms(nR)=Cor2hInt(l,nR)
               end do
               !-- Integrate in radius
               CorRmsL=rInt_R(Rms(nRC),n_r_maxC,n_cheb_maxC,dr_facC,chebt_RMS)
               !-- Add for total Rms
               CorRms = CorRms+CorRmsL
               !-- Finish Rms for mode l
               CorRmsL=sqrt(CorRmsL/volC)
            end do
         end if
         CorRms=sqrt(CorRms/volC)

         !-- Advection
         AdvRms=0.0_cp
         if ( l_conv_nl ) then
            do l=0,l_max
               do nR=1,n_r_max
                  Rms(nR)=Adv2hInt(l,nR)
               end do
               AdvRmsL=rInt_R(Rms(nRC),n_r_maxC,n_cheb_maxC,dr_facC,chebt_RMS)
               AdvRms =AdvRms+AdvRmsL
               AdvRmsL=sqrt(AdvRmsL/volC)
            end do
         end if
         AdvRms=sqrt(AdvRms/volC)

         !-- Lorentz force
         LFRms=0.0_cp
         if ( l_mag_LF ) then
            do l=0,l_max
               do nR=1,n_r_max
                  Rms(nR)=LF2hInt(l,nR)
               end do
               LFRmsL=rInt_R(Rms(nRC),n_r_maxC,n_cheb_maxC,dr_facC,chebt_RMS)
               LFRms =LFRms+LFRmsL
               LFRmsL=sqrt(LFRmsL/volC)
            end do
         end if
         LFRms=sqrt(LFRms/volC)
    
         !-- Buoyancy
         BuoRms=0.0_cp
         if ( l_conv_nl ) then
            do l=0,l_max
               do nR=1,n_r_max
                  Rms(nR)=Buo2hInt(l,nR)
               end do
               BuoRmsL=rInt_R(Rms(nRC),n_r_maxC,n_cheb_maxC,dr_facC,chebt_RMS)
               BuoRms =BuoRms+BuoRmsL
               BuoRmsL=sqrt(BuoRmsL/volC)
            end do
         end if
         BuoRms=sqrt(BuoRms/volC)
    
         !-- Pressure gradient
         PreRms=0.0_cp
         do l=0,l_max
            do nR=1,n_r_max
               Rms(nR)=Pre2hInt(l,nR)
            end do
            PreRmsL=rInt_R(Rms(nRC),n_r_maxC,n_cheb_maxC,dr_facC,chebt_RMS)
            PreRms =PreRms+PreRmsL
            PreRmsL=sqrt(PreRmsL/volC)
         end do
         PreRms=sqrt(PreRms/volC)

         !-- Geostrophic balance
         GeoRms=0.0_cp
         do l=0,l_max
            do nR=1,n_r_max
               Rms(nR)=Geo2hInt(l,nR)
            end do
            GeoRmsL=rInt_R(Rms(nRC),n_r_maxC,n_cheb_maxC,dr_facC,chebt_RMS)
            GeoRms =GeoRms+GeoRmsL
            GeoRmsL=sqrt(GeoRmsL/volC)
         end do
         GeoRms=sqrt(GeoRms/volC)

         !-- Magnetostrophic balance
         MagRms=0.0_cp
         if ( l_mag_LF ) then
            do l=0,l_max
               do nR=1,n_r_max
                  Rms(nR)=Mag2hInt(l,nR)
               end do
               MagRmsL=rInt_R(Rms(nRC),n_r_maxC,n_cheb_maxC,dr_facC,chebt_RMS)
               MagRms =MagRms+MagRmsL
               MagRmsL=sqrt(MagRmsL/volC)
            end do
         end if
         MagRms=sqrt(MagRms/volC)

         !-- Coriolis/Lorentz balance:
         CLFRms=0.0_cp
         if ( l_mag_LF ) then
            do l=0,l_max
               do nR=1,n_r_max
                  Rms(nR)=CLF2hInt(l,nR)
               end do
               CLFRmsL=rInt_R(Rms(nRC),n_r_maxC,n_cheb_maxC,dr_facC,chebt_RMS)
               CLFRms =CLFRms+CLFRmsL
               CLFRmsL=sqrt(CLFRmsL/volC)
            end do
         end if
         CLFRms=sqrt(CLFRms/volC)

         !-- Pressure/Lorentz balance:
         PLFRms=0.0_cp
         if ( l_mag_LF ) then
            do l=0,l_max
               do nR=1,n_r_max
                  Rms(nR)=PLF2hInt(l,nR)
               end do
               PLFRmsL=rInt_R(Rms(nRC),n_r_maxC,n_cheb_maxC,dr_facC,chebt_RMS)
               PLFRms =PLFRms+PLFRmsL
               PLFRmsL=sqrt(PLFRmsL/volC)
            end do
         end if
         PLFRms=sqrt(PLFRms/volC)

         !-- Archimedian balance:
         ArcRms=0.0_cp
         if ( l_conv_nl ) then
            do l=0,l_max
               do nR=1,n_r_max
                  Rms(nR)=Arc2hInt(l,nR)
               end do
               ArcRmsL=rInt_R(Rms(nRC),n_r_maxC,n_cheb_maxC,dr_facC,chebt_RMS)
               ArcRms =ArcRms+ArcRmsL
               ArcRmsL=sqrt(ArcRmsL/volC)
            end do
         end if
         ArcRms=sqrt(ArcRms/volC)

         !-- Diffusion
         DifRms=0.0_cp
         call get_drNS(DifPolLMr(1,nRC),workA(1,nRC),        &
              &        lm_max,1,lm_max,n_r_maxC,n_cheb_maxC, &
              &        workB,chebt_RMS,dr_facC)
         do nR=1,n_r_maxC
            call hInt2dPol( workA(1,nR+nCut),2,lm_max,DifPol2hInt(1,nR+nCut,1), &
                             lo_map )
         end do
         do l=0,l_max
            do nR=1,n_r_maxC
               Dif2hInt(nR+nCut)=0.0_cp
               do n=1,1
                  Dif2hInt(nR+nCut)=Dif2hInt(nR+nCut) +        &
                                    DifPol2hInt(l,nR+nCut,n) + &
                                    DifTor2hInt(l,nR+nCut,n)
               end do
           end do
           DifRmsL=rInt_R(Dif2hInt(nRC),n_r_maxC,n_cheb_maxC,dr_facC,chebt_RMS)
           DifRms =DifRms+DifRmsL
           DifRmsL=sqrt(DifRmsL/volC)
         end do
         DifRms=sqrt(DifRms/volC)

         !-- Flow changes
         dtV_Rms=0.0_cp
         call get_drNS( dtVPolLMr(1,nRC),workA(1,nRC),lm_max,1,lm_max,&
              &         n_r_maxC,n_cheb_maxC,workB,chebt_RMS,dr_facC)
         do nR=1,n_r_maxC
            call hInt2dPol( workA(1,nR+nCut),2,lm_max,dtVPol2hInt(1,nR+nCut,1), &
                            lo_map)
         end do
         do l=0,l_max
            do nR=1,n_r_maxC
               dtV2hInt(nR+nCut)=0.0_cp
               do n=1,1
                  dtV2hInt(nR+nCut)=dtV2hInt(nR+nCut) +         &
                                    dtVPol2hInt(l,nR+nCut,n) +  &
                                    dtVTor2hInt(l,nR+nCut,n)
               end do
            end do
            dtVRmsL=rInt_R(dtV2hInt(nRC),n_r_maxC,n_cheb_maxC,dr_facC,chebt_RMS)
            dtV_Rms =dtV_Rms+dtVRmsL
            dtVRmsL=sqrt(dtVRmsL/volC)
         end do
         dtV_Rms=sqrt(dtV_Rms/volC)
    
         !----- Output:
         if ( l_save_out) then
            open(n_dtvrms_file, file=dtvrms_file, form='formatted', &
                 status='unknown', position='append')
         end if
         write(n_dtvrms_file,'(1P,ES20.12,7ES16.8,5ES14.6)')&
              time, dtV_Rms, CorRms, LFRms, AdvRms, DifRms, &
              BuoRms, PreRms, GeoRms/(CorRms+PreRms),       &
              MagRms/(CorRms+PreRms+LFRms),                 &
              ArcRms/(CorRms+PreRms+LFRms+BuoRms),          &
              CLFRms/(CorRms+LFRms), PLFRms/(PreRms+LFRms)
         if ( l_save_out) then
            close(n_dtvrms_file)
         end if
    
      end if

   end subroutine dtVrms
!----------------------------------------------------------------------------
   subroutine dtBrms(time)

      !-- Input of variables:
      real(cp), intent(in) :: time
    
      !-- Local
      integer :: nR,n,l1m0,l1m1,lm,m
      character(len=80) :: fileName
    
      real(cp) :: dtBPolRms,dtBPolAsRms
      real(cp) :: dtBTorRms,dtBTorAsRms
      real(cp) :: DdynRms,DdynAsRms
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

      real(cp) :: dtBP(n_r_max),dtBPAs(n_r_max)
      real(cp) :: dtBT(n_r_max),dtBTAs(n_r_max)
    
      !-- For new movie output
      integer :: nField,nFields,nFieldSize
      integer :: nTheta,nThetaN,nThetaS,nThetaStart
      integer :: nPos
      real(cp) :: dumm(12),rS
      real(cp) :: fOut(n_theta_max*n_r_max)
      real(cp) :: outBlock(nfs)
      character(len=80) :: version
      logical :: lRmsMov
    
      real(cp) :: global_sum(lm_max,n_r_max)
    
#ifdef WITH_MPI
      call myAllGather(dtBPolLMr,lm_maxMag,n_r_maxMag)

      call MPI_Reduce(dtBPol2hInt(:,:,1),global_sum,n_r_max*lm_max, &
           &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) dtBPol2hInt(:,:,1)=global_sum
      call MPI_Reduce(dtBTor2hInt(:,:,1),global_sum,n_r_max*lm_max, &
           &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) dtBTor2hInt(:,:,1)=global_sum
#endif
    
    
      if ( rank == 0 ) then
         l1m0=lm2(1,0)
         l1m1=lm2(1,1)
    
         !--- Stretching
         call get_drNS(PstrLM,workA,lm_max,1,lm_max,n_r_max, &
              &        n_cheb_max,workB,chebt_oc,drx)
    
         !--- Add to the total dynamo term
         do nR=1,n_r_max
            do lm=1,lm_max
               PdynLM(lm,nR)  =PstrLM(lm,nR)
               drPdynLM(lm,nR)=workA(lm,nR)
            end do
         end do

         !-- Finalize advection
         call get_drNS(PadvLM,workA,lm_max,1,lm_max,n_r_max, &
                        n_cheb_max,workB,chebt_oc,drx)

         !-- Add to total dynamo term:
         do nR=1,n_r_max
            do lm=1,lm_max
               PdynLM(lm,nR)  =PdynLM(lm,nR)-PadvLM(lm,nR)
               drPdynLM(lm,nR)=drPdynLM(lm,nR)-workA(lm,nR)
               TdynLM(lm,nR)  =TstrLM(lm,nR)-TadvLM(lm,nR)
            end do
         end do

         !--- Get RMS values of the total dynamo term:
         call get_PolTorRms(PdynLM,drPdynLM,TdynLM,PdynRms,TdynRms,PdynAsRms, &
                            TdynAsRms,st_map)

         !--- Finalize diffusion:
         call get_drNS(PdifLM,workA,lm_max,1,lm_max,n_r_max, &
              &        n_cheb_max,workB,chebt_oc,drx)

         !-- Get RMS values for diffusion
         call get_PolTorRms(PdifLM,workA,TdifLM,PdifRms,TdifRms,&
                            PdifAsRms,TdifAsRms,st_map)

         !--- Get Omega effect rms: total toroidal field changes due to zonal flow
         !    (this is now stretching plus advection, changed May 23 2013):
         !    TomeAsRms is the rms of the more classical Omega effect which
         !    decribes the creation of axisymmetric azimuthal field by zonal flow.
         call get_PolTorRms(PdifLM,workA,TomeLM,dummy1,TomeRms,dummy2, &
                            TomeAsRms,st_map)

         !--- B changes:
         call get_drNS(dtBPolLMr,workA,lm_max,1,lm_max,n_r_max, &
              &        n_cheb_max,workB,chebt_oc,drx)

         do nR=1,n_r_max
            call hInt2dPolLM(workA(1,nR),2,lm_max,dtBPol2hInt(1,nR,1),lo_map)
            dtBP(nR)  =0.0_cp
            dtBT(nR)  =0.0_cp
            dtBPAs(nR)=0.0_cp
            dtBTAs(nR)=0.0_cp
            do n=1,1
               do lm=1,lm_max
                  m=lm2m(lm)
                  dtBP(nR)=dtBP(nR)+dtBPol2hInt(lm,nR,n)
                  dtBT(nR)=dtBT(nR)+dtBTor2hInt(lm,nR,n)
                  if ( m == 0 ) then
                     dtBPAs(nR)=dtBPAs(nR)+dtBPol2hInt(lm,nR,n)
                     dtBTAs(nR)=dtBTAs(nR)+dtBTor2hInt(lm,nR,n)
                  end if
               end do
            end do
         end do

         dtBPolRms  =rInt(dtBP,n_r_max,dr_fac,chebt_oc)
         dtBPolAsRms=rInt(dtBPAs,n_r_max,dr_fac,chebt_oc)
         dtBTorRms  =rInt(dtBT,n_r_max,dr_fac,chebt_oc)
         dtBTorAsRms=rInt(dtBTAs,n_r_max,dr_fac,chebt_oc)

         dtBPolRms  =sqrt(dtBPolRms  /vol_oc)
         dtBPolAsRms=sqrt(dtBPolAsRms/vol_oc)
         dtBTorRms  =sqrt(dtBTorRms  /vol_oc)
         dtBTorAsRms=sqrt(dtBTorAsRms/vol_oc)


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
                        call get_RAS(PdynLM(1,nR),outBlock,rS,nThetaStart,sizeThetaB)
                     else if ( nField == 3 ) then
                        ! Note that PdynLM stores PdifLM at this point!
                        call get_RAS(PdifLM(1,nR),outBlock,rS,nThetaStart,sizeThetaB)
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

         !-- Get dipole dynamo contribution:
         do nR=1,n_r_max
            do lm=1,lm_max
               if ( lm/=l1m0 .and. lm/=l1m1 ) THEN
                  PdynLM(lm,nR)  =zero
                  drPdynLM(lm,nR)=zero
               end if
            end do
         end do
         !-- Get dipole dynamo terms:
         call get_PolTorRms(PdynLM,drPdynLM,TdynLM,DdynRms,dummy1, &
                            DdynAsRms,dummy3,st_map)

         !-- Output:
         if ( l_save_out) then
            open(n_dtbrms_file, file=dtbrms_file, form='formatted', &
                 status='unknown', position='append')
         end if
         write(n_dtbrms_file,'(1P,ES20.12,10ES16.8)')            &
              time, dtBPolRms, dtBTorRms, PdynRms, TdynRms,      &
              PdifRms, TdifRms, TomeRms/TdynRms,                 &
              TomeAsRms/TdynRms,  DdynRms,DdynAsRms
         if ( l_save_out) then
            close(n_dtbrms_file)
         end if


      end if

   end subroutine dtBrms
!----------------------------------------------------------------------------
end module out_RMS
