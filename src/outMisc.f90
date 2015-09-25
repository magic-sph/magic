module outMisc_mod

   use parallel_mod
   use precision_mod
   use truncation, only: l_max, n_r_max, lm_max
   use radial_data, only: n_r_icb, n_r_cmb, nRstart, nRstop
   use radial_functions, only: botcond, r_icb, dr_fac, i_costf_init, &
                               d_costf_init, topcond, kappa, r_cmb,  &
                               temp0, r, rho0, dtemp0
   use physical_parameters, only: epsS
   use num_param, only: lScale
   use blocking, only: lo_map, nThetaBs, nfs, sizeThetaB
   use horizontal_data, only: gauss
   use logic, only: l_save_out, l_anelastic_liquid, l_par, &
                    l_hel, l_heat
   use output_data, only: tag, misc_file, n_misc_file
   use Egeos_mod, only: getEgeos
   use constants, only: pi, vol_oc, osq4pi, sq4pi, one, two, four
   use useful, only: cc2real
   use integration, only: rInt,rInt_R
   use LMLoop_data,only: llm,ulm
   use legendre_spec_to_grid, only: lmAS2pt

   implicit none

   private

   public :: outMisc

contains

   subroutine outMisc(timeScaled,HelLMr,Hel2LMr,HelnaLMr,Helna2LMr, &
     &             nLogs,w,dw,ddw,z,dz,s,ds,Geos,dpFlow,dzFlow)

      !-- Input of variables:
      real(cp),    intent(in) :: timeScaled
      real(cp),    intent(in) :: HelLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(in) :: Hel2LMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(in) :: HelnaLMr(l_max+1,nRstart:nRstop)
      real(cp),    intent(in) :: Helna2LMr(l_max+1,nRstart:nRstop)
      integer,     intent(in) :: nLogs
    
      !-- Input of scalar fields:
      complex(cp), intent(in) :: s(llm:ulm,n_r_max)
      complex(cp), intent(in) :: ds(llm:ulm,n_r_max)
      !---- Fields transfered to getEgeos:
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dw(llm:ulm,n_r_max)
      complex(cp), intent(in) :: ddw(llm:ulm,n_r_max)
      complex(cp), intent(in) :: z(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dz(llm:ulm,n_r_max)
    
      !-- Output: (and stuff written in misc.TAG files)
      real(cp),    intent(out) :: Geos
      real(cp),    intent(out) :: dpFlow,dzFlow
    
      !-- Local stuff:
      integer :: nTheta,nThetaStart,nThetaBlock,nThetaNHS,n
      real(cp) :: HelNr(nRstart:nRstop), HelSr(nRstart:nRstop)
      real(cp) :: HelnaNr(nRstart:nRstop), HelnaSr(nRstart:nRstop)
      real(cp) :: Hel2Nr(nRstart:nRstop), Hel2Sr(nRstart:nRstop)
      real(cp) :: Helna2Nr(nRstart:nRstop), Helna2Sr(nRstart:nRstop)
      real(cp) :: HelEAr(nRstart:nRstop)
      real(cp) :: HelNr_global(n_r_max), HelSr_global(n_r_max)
      real(cp) :: HelnaNr_global(n_r_max), HelnaSr_global(n_r_max)
      real(cp) :: Helna2Nr_global(n_r_max), Helna2Sr_global(n_r_max)
      real(cp) :: Hel2Nr_global(n_r_max), Hel2Sr_global(n_r_max)
      real(cp) :: HelEAr_global(n_r_max)
      real(cp) :: Hel(nfs), Hel2(nfs), Helna(nfs), Helna2(nfs), r2
      real(cp) :: HelN,HelS
      real(cp) :: HelnaN,HelnaS
      real(cp) :: HelnaRMSN,HelnaRMSS
      real(cp) :: HelRMSN,HelRMSS,HelEA,HelRMS,HelnaRMS
      real(cp) :: Egeos,EkNTC,EkSTC,Ekin
      real(cp) :: CVzOTC,CVorOTC,CHelOTC
      real(cp) :: topnuss,botnuss
      real(cp) :: topflux,botflux
    
      integer :: n_r
      integer :: i,sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1),ierr
    
    
      if ( l_hel )  then
    
         !------ Integration of Helicity, on input the Helicity is
         !       already axisymmetric !
         do n_r=nRstart,nRstop
            r2=r(n_r)*r(n_r)
            HelNr(n_r) =0.0_cp
            HelSr(n_r) =0.0_cp
            HelnaNr(n_r) =0.0_cp
            HelnaSr(n_r) =0.0_cp
            HelEAr(n_r)=0.0_cp
            Hel2Nr(n_r) =0.0_cp
            Hel2Sr(n_r) =0.0_cp
            Helna2Nr(n_r) =0.0_cp
            Helna2Sr(n_r) =0.0_cp
    
            do n=1,nThetaBs ! Loop over theta blocks
               nTheta=(n-1)*sizeThetaB
               nThetaStart=nTheta+1
               call lmAS2pt(HelLMr(1,n_r),Hel,nThetaStart,sizeThetaB)
               call lmAS2pt(Hel2LMr(1,n_r),Hel2,nThetaStart,sizeThetaB)
               call lmAS2pt(HelnaLMr(1,n_r),Helna,nThetaStart,sizeThetaB)
               call lmAS2pt(Helna2LMr(1,n_r),Helna2,nThetaStart,sizeThetaB)
               do nThetaBlock=1,sizeThetaB
                  nTheta=nTheta+1
                  nThetaNHS=(nTheta+1)/2
    
                  !------ Integration over theta:
                  if ( mod(nTheta,2) == 1 ) then ! NHS
                     Hel2Nr(n_r)=Hel2Nr(n_r)+gauss(nThetaNHS)*r2*Hel2(nThetaBlock)
                     Helna2Nr(n_r)=Helna2Nr(n_r)+gauss(nThetaNHS)*r2*Helna2(nThetaBlock)
                     HelEAr(n_r)=HelEAr(n_r)+gauss(nThetaNHS)*r2*Hel(nThetaBlock)
                     HelNr(n_r) =HelNr(n_r)+gauss(nThetaNHS)*r2*Hel(nThetaBlock)
                     HelnaNr(n_r) =HelnaNr(n_r)+gauss(nThetaNHS)*r2*Helna(nThetaBlock)
                  else
                     Hel2Sr(n_r)=Hel2Sr(n_r)+gauss(nThetaNHS)*r2*Hel2(nThetaBlock)
                     Helna2Sr(n_r)=Helna2Sr(n_r)+gauss(nThetaNHS)*r2*Helna2(nThetaBlock)
                     HelEAr(n_r)=HelEAr(n_r)-gauss(nThetaNHS)*r2*Hel(nThetaBlock)
                     HelSr(n_r) =HelSr(n_r)+gauss(nThetaNHS)*r2*Hel(nThetaBlock)
                     HelnaSr(n_r)=HelnaSr(n_r)+gauss(nThetaNHS)*r2*Helna(nThetaBlock)
                  end if
               end do
            end do
    
         end do
    
         ! Now we have to gather the results on rank 0 for
         ! the arrays: Hel2Nr,Helna2Nr,HelEAr,HelNr,HelnaNr
         ! Hel2Sr,Helna2Sr,HelSr,HelnaSr
    
         sendcount  = (nRstop-nRstart+1)
         recvcounts = nr_per_rank
         recvcounts(n_procs-1) = (nr_per_rank+1)
         do i=0,n_procs-1
            displs(i) = i*nr_per_rank
         end do
#ifdef WITH_MPI
         call MPI_GatherV(Hel2Nr,sendcount,MPI_DEF_REAL,&
              &           Hel2Nr_global,recvcounts,displs,MPI_DEF_REAL,&
              &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(Helna2Nr,sendcount,MPI_DEF_REAL,&
              &           Helna2Nr_global,recvcounts,displs,MPI_DEF_REAL,&
              &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(HelEAr,sendcount,MPI_DEF_REAL,&
              &           HelEAr_global,recvcounts,displs,MPI_DEF_REAL,&
              &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(HelNr,sendcount,MPI_DEF_REAL,&
              &           HelNr_global,recvcounts,displs,MPI_DEF_REAL,&
              &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(HelnaNr,sendcount,MPI_DEF_REAL,&
              &           HelnaNr_global,recvcounts,displs,MPI_DEF_REAL,&
              &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(HelSr,sendcount,MPI_DEF_REAL,&
              &           HelSr_global,recvcounts,displs,MPI_DEF_REAL,&
              &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(Helna2Sr,sendcount,MPI_DEF_REAL,&
              &           Helna2Sr_global,recvcounts,displs,MPI_DEF_REAL,&
              &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(Hel2Sr,sendcount,MPI_DEF_REAL,&
              &           Hel2Sr_global,recvcounts,displs,MPI_DEF_REAL,&
              &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(HelnaSr,sendcount,MPI_DEF_REAL,&
              &           HelnaSr_global,recvcounts,displs,MPI_DEF_REAL,&
              &           0,MPI_COMM_WORLD,ierr)
#else
         Hel2Nr_global=Hel2Nr
         Helna2Nr_global=Helna2Nr
         HelEAr_global=HelEAr
         HelNr_global=HelNr
         HelnaNr_global=HelnaNr
         HelSr_global=HelSr
         Helna2Sr_global=Helna2Sr
         Hel2Sr_global=Hel2Sr
         HelnaSr_global=HelnaSr
#endif
    
         if ( rank == 0 ) then
            !------ Integration over r without the boundaries and normalization:
            HelN  =rInt(HelNr_global,n_r_max,dr_fac,i_costf_init,d_costf_init)
            HelS  =rInt(HelSr_global,n_r_max,dr_fac,i_costf_init,d_costf_init)
            HelnaN=rInt(HelnaNr_global,n_r_max,dr_fac,i_costf_init,d_costf_init)
            HelnaS=rInt(HelnaSr_global,n_r_max,dr_fac,i_costf_init,d_costf_init)
            HelEA =rInt(HelEAr_global,n_r_max,dr_fac,i_costf_init,d_costf_init)
            HelRMSN=rInt(Hel2Nr_global,n_r_max,dr_fac,i_costf_init,d_costf_init)
            HelRMSS=rInt(Hel2Sr_global,n_r_max,dr_fac,i_costf_init,d_costf_init)
            HelnaRMSN=rInt(Helna2Nr_global,n_r_max,dr_fac,i_costf_init,d_costf_init)
            HelnaRMSS=rInt(Helna2Sr_global,n_r_max,dr_fac,i_costf_init,d_costf_init)
    
            HelN  =two*pi*HelN/(vol_oc/2) ! Note integrated over half spheres only !
            HelS  =two*pi*HelS/(vol_oc/2) ! Factor 2*pi is from phi integration
            HelnaN=two*pi*HelnaN/(vol_oc/2) ! Note integrated over half spheres only !
            HelnaS=two*pi*HelnaS/(vol_oc/2) ! Factor 2*pi is from phi integration
            HelEA =two*pi*HelEA/vol_oc
            HelRMSN=sqrt(two*pi*HelRMSN/(vol_oc/2))
            HelRMSS=sqrt(two*pi*HelRMSS/(vol_oc/2))
            HelnaRMSN=sqrt(two*pi*HelnaRMSN/(vol_oc/2))
            HelnaRMSS=sqrt(two*pi*HelnaRMSS/(vol_oc/2))
            HelRMS=HelRMSN+HelRMSS
            HelnaRMS=HelnaRMSN+HelnaRMSS
    
            if ( HelnaRMS /= 0 ) then
               HelnaN =HelnaN/HelnaRMSN
               HelnaS =HelnaS/HelnaRMSS
            else
               HelnaN =0.0_cp
               HelnaS =0.0_cp
            end if
            if ( HelRMS /= 0 ) then
               HelN =HelN/HelRMSN
               HelS =HelS/HelRMSS
               HelEA=HelEA/HelRMS
            else
               HelN =0.0_cp
               HelS =0.0_cp
               HelEA=0.0_cp
            end if
         end if
      else
         HelN     =0.0_cp
         HelS     =0.0_cp
         HelEA    =0.0_cp
         HelRMSN  =0.0_cp
         HelRMSS  =0.0_cp
         HelnaN   =0.0_cp
         HelnaS   =0.0_cp
         HelnaRMSN=0.0_cp
         HelnaRMSS=0.0_cp
      end if
    
      if ( l_par ) then
         call getEgeos(timeScaled,nLogs,w,dw,ddw,z,dz, &
              &        Egeos,EkNTC,EkSTC,Ekin,         &
              &        dpFlow,dzFlow,CVzOTC,CVorOTC,CHelOTC)
         if ( Ekin > 0.0_cp ) then
            Geos=Egeos/Ekin ! Output, relative geostrophic kinetic Energy
         else
            Geos=0.0_cp
            Ekin=-one ! Only used for ratio, must thus be non-zero
         end if
      else
         Egeos  =0.0_cp
         EkNTC  =0.0_cp
         EkSTC  =0.0_cp
         Ekin   =-one ! Only used for ratio, must thus be non-zero
         dpFlow =0.0_cp
         dzFlow =0.0_cp
         Geos   =0.0_cp
         CVzOTC =0.0_cp
         CVorOTC=0.0_cp
         CHelOTC=0.0_cp
      end if
    
      if ( rank == 0 ) then
         !-- Evaluate nusselt numbers (boundary heat flux density):
         if ( topcond/=0.0_cp .and. l_heat ) then
            if ( l_anelastic_liquid ) then
               botnuss=-osq4pi/botcond*real(ds(1,n_r_icb))/lScale+one
               topnuss=-osq4pi/topcond*real(ds(1,n_r_cmb))/lScale+one
               botflux=-rho0(n_r_max)*(real(ds(1,n_r_max))*osq4pi+ &
                        one/epsS*dtemp0(n_r_max))*r_icb**2*four*pi*kappa(n_r_max)
               topflux=-rho0(1)*(real(ds(1,1))*osq4pi+one/epsS*dtemp0(1))*r_cmb**2* &
                        four*pi*kappa(1)
            else
               botnuss=-osq4pi/botcond*real(ds(1,n_r_icb))/lScale
               topnuss=-osq4pi/topcond*real(ds(1,n_r_cmb))/lScale
               botflux=-rho0(n_r_max)*temp0(n_r_max)*real(ds(1,n_r_max))/lScale* &
                        r_icb**2*sq4pi*kappa(n_r_max)
               topflux=-rho0(1)*temp0(1)*real(ds(1,1))/lScale*r_cmb**2* &
                        sq4pi*kappa(1)
            end if
         else
            botnuss=0.0_cp
            topnuss=0.0_cp
            botflux=0.0_cp
            topflux=0.0_cp
         end if
    
         if ( l_save_out ) then
            open(n_misc_file, file=misc_file, status='unknown', position='append')
         end if

         write(n_misc_file,'(1P,ES20.12,21ES16.8)')       &
              & timeScaled, botnuss, topnuss,             &
              & real(s(1,n_r_icb)), real(s(1,n_r_cmb)),   &
              & HelN, HelS, HelRMSN, HelRMSS,             &
              & Geos, EkNTC/Ekin, EkSTC/Ekin, Ekin,       & !10-13
              & CVzOTC, CVorOTC, CHelOTC,                 & !14-16   
              & HelnaN, HelnaS, HelnaRMSN, HelnaRMSS,     &
              & botflux, topflux

         if ( l_save_out ) close(n_misc_file)
         !--- NOTE: Ekin can be compared with energy in e_kin.TAG to
         !    get an idea of the precision of cylindrical integration in getEgeos.
      end if
    
   end subroutine outMisc
!---------------------------------------------------------------------------
end module outMisc_mod
