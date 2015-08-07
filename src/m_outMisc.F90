!$Id$
module outMisc_mod

   use parallel_mod
   use truncation, only: l_max, n_r_max, lm_max
   use radial_data, only: n_r_icb, n_r_cmb, nRstart, nRstop
   use radial_functions, only: botcond, r_icb, dr_fac, i_costf_init, &
                               d_costf_init, topcond, kappa, r_cmb,  &
                               temp0, r, rho0, dtemp0
   use physical_parameters, only: epsS
   use num_param, only: lScale
   use blocking, only: lo_map, nThetaBs, nfs, sizeThetaB
   use horizontal_data, only: gauss
   use logic, only: l_save_out, l_anelastic_liquid, l_prms, l_par, &
                    l_hel, l_heat
   use output_data, only: tag, misc_file, n_misc_file
   use Egeos_mod, only: getEgeos
   use const, only: pi, vol_oc
   use useful, only: cc2real
   use integration, only: rInt,rInt_R
   use LMLoop_data,only: llm,ulm
   use legendre_spec_to_grid, only: lmAS2pt

   implicit none

   private

   public :: outMisc

contains

   subroutine outMisc(timeScaled,HelLMr,Hel2LMr,HelnaLMr,Helna2LMr, &
     &             nLogs,w,dw,ddw,z,dz,s,ds,p,Geos,dpFlow,dzFlow)

      !-- Input of variables:
      real(kind=8),    intent(in) :: timeScaled
      real(kind=8),    intent(in) :: HelLMr(l_max+1,nRstart:nRstop)
      real(kind=8),    intent(in) :: Hel2LMr(l_max+1,nRstart:nRstop)
      real(kind=8),    intent(in) :: HelnaLMr(l_max+1,nRstart:nRstop)
      real(kind=8),    intent(in) :: Helna2LMr(l_max+1,nRstart:nRstop)
      integer,         intent(in) :: nLogs
    
      !-- Input of scalar fields:
      complex(kind=8), intent(in) :: s(llm:ulm,n_r_max)
      complex(kind=8), intent(in) :: ds(llm:ulm,n_r_max)
      !---- Fields transfered to getEgeos:
      complex(kind=8), intent(in) :: w(llm:ulm,n_r_max)
      complex(kind=8), intent(in) :: dw(llm:ulm,n_r_max)
      complex(kind=8), intent(in) :: ddw(llm:ulm,n_r_max)
      complex(kind=8), intent(in) :: z(llm:ulm,n_r_max)
      complex(kind=8), intent(in) :: dz(llm:ulm,n_r_max)
      complex(kind=8), intent(in) :: p(llm:ulm,n_r_max)
    
      !-- Output: (and stuff written in misc.TAG files)
      real(kind=8),    intent(out) :: Geos
      real(kind=8),    intent(out) :: dpFlow,dzFlow
    
      !-- Local stuff:
      integer :: nTheta,nThetaStart,nThetaBlock,nThetaNHS,n,lm44
      logical :: lm44_is_local
      real(kind=8) :: pplot_global(n_r_max), pplot(n_r_max)
      real(kind=8) :: HelNr(nRstart:nRstop), HelSr(nRstart:nRstop)
      real(kind=8) :: HelnaNr(nRstart:nRstop), HelnaSr(nRstart:nRstop)
      real(kind=8) :: Hel2Nr(nRstart:nRstop), Hel2Sr(nRstart:nRstop)
      real(kind=8) :: Helna2Nr(nRstart:nRstop), Helna2Sr(nRstart:nRstop)
      real(kind=8) :: HelEAr(nRstart:nRstop)
      real(kind=8) :: HelNr_global(n_r_max), HelSr_global(n_r_max)
      real(kind=8) :: HelnaNr_global(n_r_max), HelnaSr_global(n_r_max)
      real(kind=8) :: Helna2Nr_global(n_r_max), Helna2Sr_global(n_r_max)
      real(kind=8) :: Hel2Nr_global(n_r_max), Hel2Sr_global(n_r_max)
      real(kind=8) :: HelEAr_global(n_r_max)
      complex(kind=8) :: p44_local(n_r_max)
      real(kind=8) :: Hel(nfs), Hel2(nfs), Helna(nfs), Helna2(nfs), r2
      real(kind=8) :: HelN,HelS
      real(kind=8) :: HelnaN,HelnaS
      real(kind=8) :: HelnaRMSN,HelnaRMSS
      real(kind=8) :: HelRMSN,HelRMSS,HelEA,HelRMS,HelnaRMS
      real(kind=8) :: Egeos,EkNTC,EkSTC,Ekin
      real(kind=8) :: CVzOTC,CVorOTC,CHelOTC
      real(kind=8) :: topnuss,botnuss
      real(kind=8) :: topflux,botflux
    
      integer :: n_r,m,lm,mytag,status(MPI_STATUS_SIZE)
      real(kind=8) :: osq4pi
      integer :: i,sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1),ierr
    
      character(len=76) :: filename2
    
    
      if ( l_hel )  then
    
         !------ Integration of Helicity, on input the Helicity is
         !       already axisymmetric !
         do n_r=nRstart,nRstop
            r2=r(n_r)*r(n_r)
            HelNr(n_r) =0.D0
            HelSr(n_r) =0.D0
            HelnaNr(n_r) =0.D0
            HelnaSr(n_r) =0.D0
            HelEAr(n_r)=0.D0
            Hel2Nr(n_r) =0.D0
            Hel2Sr(n_r) =0.D0
            Helna2Nr(n_r) =0.D0
            Helna2Sr(n_r) =0.D0
    
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
         call MPI_GatherV(Hel2Nr,sendcount,MPI_DOUBLE_PRECISION,&
              &           Hel2Nr_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
              &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(Helna2Nr,sendcount,MPI_DOUBLE_PRECISION,&
              &           Helna2Nr_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
              &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(HelEAr,sendcount,MPI_DOUBLE_PRECISION,&
              &           HelEAr_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
              &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(HelNr,sendcount,MPI_DOUBLE_PRECISION,&
              &           HelNr_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
              &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(HelnaNr,sendcount,MPI_DOUBLE_PRECISION,&
              &           HelnaNr_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
              &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(HelSr,sendcount,MPI_DOUBLE_PRECISION,&
              &           HelSr_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
              &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(Helna2Sr,sendcount,MPI_DOUBLE_PRECISION,&
              &           Helna2Sr_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
              &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(Hel2Sr,sendcount,MPI_DOUBLE_PRECISION,&
              &           Hel2Sr_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
              &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(HelnaSr,sendcount,MPI_DOUBLE_PRECISION,&
              &           HelnaSr_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
              &           0,MPI_COMM_WORLD,ierr)
    
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
    
            HelN  =2.D0*pi*HelN/(vol_oc/2) ! Note integrated over half spheres only !
            HelS  =2.D0*pi*HelS/(vol_oc/2) ! Factor 2*pi is from phi integration
            HelnaN  =2.D0*pi*HelnaN/(vol_oc/2) ! Note integrated over half spheres only !
            HelnaS  =2.D0*pi*HelnaS/(vol_oc/2) ! Factor 2*pi is from phi integration
            HelEA =2.D0*pi*HelEA/vol_oc
            HelRMSN=dsqrt(2.D0*pi*HelRMSN/(vol_oc/2))
            HelRMSS=dsqrt(2.D0*pi*HelRMSS/(vol_oc/2))
            HelnaRMSN=dsqrt(2.D0*pi*HelnaRMSN/(vol_oc/2))
            HelnaRMSS=dsqrt(2.D0*pi*HelnaRMSS/(vol_oc/2))
            HelRMS=HelRMSN+HelRMSS
            HelnaRMS=HelnaRMSN+HelnaRMSS
    
            if ( HelnaRMS /= 0 ) then
               HelnaN =HelnaN/HelnaRMSN
               HelnaS =HelnaS/HelnaRMSS
            else
               HelnaN =0.D0
               HelnaS =0.D0
            end if
            if ( HelRMS /= 0 ) then
               HelN =HelN/HelRMSN
               HelS =HelS/HelRMSS
               HelEA=HelEA/HelRMS
            else
               HelN =0.D0
               HelS =0.D0
               HelEA=0.D0
            end if
         end if
      else
         HelN     =0.D0
         HelS     =0.D0
         HelEA    =0.D0
         HelRMSN  =0.D0
         HelRMSS  =0.D0
         HelnaN   =0.D0
         HelnaS   =0.D0
         HelnaRMSN=0.D0
         HelnaRMSS=0.D0
      end if
    
      if ( l_par ) then
         call getEgeos(timeScaled,nLogs,w,dw,ddw,z,dz, &
              &        Egeos,EkNTC,EkSTC,Ekin,         &
              &        dpFlow,dzFlow,CVzOTC,CVorOTC,CHelOTC)
         if ( Ekin > 0.d0 ) then
            Geos=Egeos/Ekin ! Output, relative geostrophic kinetic Energy
         else
            Geos=0.d0
            Ekin=-1.D0 ! Only used for ratio, must thus be non-zero
         end if
      else
         Egeos  =0.D0
         EkNTC  =0.D0
         EkSTC  =0.D0
         Ekin   =-1.D0 ! Only used for ratio, must thus be non-zero
         dpFlow =0.D0
         dzFlow =0.D0
         Geos   =0.D0
         CVzOTC =0.D0
         CVorOTC=0.D0
         CHelOTC=0.D0
      end if
    
      if ( rank == 0 ) then
         !-- Evaluate nusselt numbers (boundary heat flux density):
         osq4pi =1.D0/dsqrt(4.D0*pi)
         if ( topcond/=0.D0 .and. l_heat ) then
            if ( l_anelastic_liquid ) then
               botnuss=-osq4pi/botcond*real(ds(1,n_r_icb))/lScale+1.D0
               topnuss=-osq4pi/topcond*real(ds(1,n_r_cmb))/lScale+1.D0
               botflux=-rho0(n_r_max)*(real(ds(1,n_r_max))*osq4pi+ &
                        1.D0/epsS*dtemp0(n_r_max))*r_icb**2*4.D0*pi*kappa(n_r_max)
               topflux=-rho0(1)*(real(ds(1,1))*osq4pi+1.D0/epsS*dtemp0(1))*r_cmb**2* &
                        4.D0*pi*kappa(1)
            else
               botnuss=-osq4pi/botcond*real(ds(1,n_r_icb))/lScale
               topnuss=-osq4pi/topcond*real(ds(1,n_r_cmb))/lScale
               botflux=-rho0(n_r_max)*temp0(n_r_max)*real(ds(1,n_r_max))/lScale* &
                        r_icb**2*dsqrt(4.D0*pi)*kappa(n_r_max)
               topflux=-rho0(1)*temp0(1)*real(ds(1,1))/lScale*r_cmb**2* &
                        dsqrt(4.D0*pi)*kappa(1)
            end if
         else
            botnuss=0.D0
            topnuss=0.D0
            botflux=0.D0
            topflux=0.D0
         end if
    
         if ( l_save_out ) then
            open(n_misc_file, file=misc_file, status='unknown', position='APPEND')
         end if

         write(n_misc_file,'(1P,D20.12,21D16.8)')         &
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
    
      if ( l_prms ) then
         do n_r=1,n_r_max
            pplot(n_r)=0.D0
            do lm=llm,ulm
               m=lo_map%lm2m(lm)
               pplot(n_r)=pplot(n_r)+cc2real(p(lm,n_r),m)
            end do
         end do
         call MPI_Reduce(pplot,pplot_global,n_r_max,MPI_DOUBLE_PRECISION,&
              & MPI_SUM,0,MPI_COMM_WORLD,ierr)
         ! Send the p(4,4) value to rank 0
         lm44=lo_map%lm2(4,4)
         lm44_is_local=(llm <= lm44).and.(lm44 <= ulm)
         mytag=120
         if ( lm44_is_local .and. ( rank /= 0 )) then
            ! copy one row of p into a vector to send
            ! it to rank 0
            p44_local=p(lm44,:)
            call MPI_Send(p44_local,n_r_max,MPI_DOUBLE_complex,0,mytag, &
                       &  MPI_COMM_WORLD,ierr)
         end if
         if ( rank == 0 ) then
            if ( .not. lm44_is_local ) then
               call MPI_Recv(p44_local,n_r_max,MPI_DOUBLE_complex, &
                    & MPI_ANY_SOURCE,mytag,MPI_COMM_WORLD,status,ierr)
            else
               p44_local=p(lm44,:)
            end if
            filename2='p.'//TAG
            open(94, file=filename2, status='UNKNOWN')
            do n_r=1,n_r_max
               pplot_global(n_r)=dsqrt(pplot_global(n_r)/lm_max)
               write(94,*) r(n_r),pplot_global(n_r),real(p44_local(n_r))
            end do
            close(94)
         end if
      end if
    
   end subroutine outMisc
!---------------------------------------------------------------------------
end module outMisc_mod
