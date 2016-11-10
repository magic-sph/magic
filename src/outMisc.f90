module outMisc_mod
   !
   ! This module contains several subroutines that can compute and store
   ! various informations: helicity, heat transfer.
   !

   use parallel_mod
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: l_max, n_r_max, lm_max
   use radial_data, only: n_r_icb, n_r_cmb, nRstart, nRstop
   use radial_functions, only: r_icb, dr_fac, chebt_oc, kappa,   &
       &                       r_cmb,temp0, r, rho0, dLtemp0,    &
       &                       dLalpha0, beta, orho1, alpha0,    &
       &                       otemp1, drx
   use physical_parameters, only: ViscHeatFac, ThExpNb, ogrun
   use num_param, only: lScale
   use blocking, only: nThetaBs, nfs, sizeThetaB
   use horizontal_data, only: gauss
   use logic, only: l_save_out, l_anelastic_liquid, l_heat, &
       &            l_temperature_diff, l_chemical_conv
   use output_data, only: tag, heat_file, n_heat_file, helicity_file, &
       &                  n_helicity_file
   use constants, only: pi, vol_oc, osq4pi, sq4pi, one, two, four
   use start_fields, only: topcond, botcond, deltacond, topxicond, botxicond, &
       &                   deltaxicond
   use useful, only: cc2real
   use integration, only: rInt, rInt_R
   use LMLoop_data,only: llm,ulm
   use legendre_spec_to_grid, only: lmAS2pt

   implicit none

   private

   real(cp), allocatable :: TMeanR(:), SMeanR(:), PMeanR(:), XiMeanR(:)

   public :: outHelicity, outHeat, initialize_outMisc_mod

contains

   subroutine initialize_outMisc_mod

      if ( l_heat .or. l_chemical_conv ) then
         allocate( TMeanR(n_r_max) )
         allocate( SMeanR(n_r_max) )
         allocate( PMeanR(n_r_max) )
         allocate( XiMeanR(n_r_max) )
         TMeanR(:)  = 0.0_cp
         SMeanR(:)  = 0.0_cp
         PMeanR(:)  = 0.0_cp
         XiMeanR(:) = 0.0_cp
      end if
      bytes_allocated=bytes_allocated+5*n_r_max*SIZEOF_DEF_REAL

   end subroutine initialize_outMisc_mod
!---------------------------------------------------------------------------
   subroutine outHelicity(timeScaled,HelLMr,Hel2LMr,HelnaLMr,Helna2LMr)
      !
      ! This subroutine is used to store informations about kinetic 
      ! helicity
      !

      !-- Input of variables:
      real(cp), intent(in) :: timeScaled
      real(cp), intent(in) :: HelLMr(l_max+1,nRstart:nRstop)
      real(cp), intent(in) :: Hel2LMr(l_max+1,nRstart:nRstop)
      real(cp), intent(in) :: HelnaLMr(l_max+1,nRstart:nRstop)
      real(cp), intent(in) :: Helna2LMr(l_max+1,nRstart:nRstop)
    
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
    
      integer :: n_r
      integer :: i,sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1),ierr
    
    
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
         HelN  =rInt(HelNr_global,n_r_max,dr_fac,chebt_oc)
         HelS  =rInt(HelSr_global,n_r_max,dr_fac,chebt_oc)
         HelnaN=rInt(HelnaNr_global,n_r_max,dr_fac,chebt_oc)
         HelnaS=rInt(HelnaSr_global,n_r_max,dr_fac,chebt_oc)
         HelEA =rInt(HelEAr_global,n_r_max,dr_fac,chebt_oc)
         HelRMSN=rInt(Hel2Nr_global,n_r_max,dr_fac,chebt_oc)
         HelRMSS=rInt(Hel2Sr_global,n_r_max,dr_fac,chebt_oc)
         HelnaRMSN=rInt(Helna2Nr_global,n_r_max,dr_fac,chebt_oc)
         HelnaRMSS=rInt(Helna2Sr_global,n_r_max,dr_fac,chebt_oc)
 
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
    
         if ( l_save_out ) then
            open(n_helicity_file, file=helicity_file, status='unknown',  &
               & position='append')
         end if

         write(n_helicity_file,'(1P,ES20.12,8ES16.8)')    &
              & timeScaled,HelN, HelS, HelRMSN, HelRMSS,  &
              & HelnaN, HelnaS, HelnaRMSN, HelnaRMSS

         if ( l_save_out ) close(n_helicity_file)
            
      end if
    
   end subroutine outHelicity
!---------------------------------------------------------------------------
   subroutine outHeat(time,timePassed,timeNorm,l_stop_time,s,ds,p,dp,xi,dxi)
      !
      ! This subroutine is used to store informations about heat transfer
      ! (i.e. Nusselt number, temperature, entropy, ...)
      !

      !-- Input of variables:
      real(cp),    intent(in) :: time
      real(cp),    intent(in) :: timePassed
      real(cp),    intent(in) :: timeNorm
      logical,     intent(in) :: l_stop_time
    
      !-- Input of scalar fields:
      complex(cp), intent(in) :: s(llm:ulm,n_r_max)
      complex(cp), intent(in) :: ds(llm:ulm,n_r_max)
      complex(cp), intent(in) :: p(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dp(llm:ulm,n_r_max)
      complex(cp), intent(in) :: xi(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dxi(llm:ulm,n_r_max)
    
      !-- Local stuff:
      real(cp) :: rhoprime(n_r_max)
      real(cp) :: topnuss,botnuss,deltanuss
      real(cp) :: topsherwood,botsherwood,deltasherwood
      real(cp) :: toptemp,bottemp
      real(cp) :: topxi,botxi
      real(cp) :: toppres,botpres,mass
      real(cp) :: topentropy,botentropy
      real(cp) :: topflux,botflux
      character(len=76) :: filename
      integer :: n_r, filehandle

      if ( rank == 0 ) then

         if ( l_anelastic_liquid ) then
            do n_r=1,n_r_max
               TMeanR(n_r)   = TMeanR(n_r)+timePassed*osq4pi*real(s(1,n_r))
               PMeanR(n_r)   = PMeanR(n_r)+timePassed*osq4pi*real(p(1,n_r))
               SMeanR(n_r)   = otemp1(n_r)*TMeanR(n_r)-ViscHeatFac*ThExpNb* &
               &               alpha0(n_r)*orho1(n_r)*PMeanR(n_r)
               rhoprime(n_r) = osq4pi*ThExpNb*alpha0(n_r)*( -rho0(n_r)* &
                  &            real(s(1,n_r))+ViscHeatFac*   &
                  &            ogrun*real(p(1,n_r)) )
            end do
         else
            do n_r=1,n_r_max
               SMeanR(n_r)   = SMeanR(n_r)+timePassed*osq4pi*real(s(1,n_r))
               PMeanR(n_r)   = PMeanR(n_r)+timePassed*osq4pi*real(p(1,n_r))
               TMeanR(n_r)   = temp0(n_r)*SMeanR(n_r)+ViscHeatFac*ThExpNb* &
                  &            alpha0(n_r)*temp0(n_r)*orho1(n_r)*PMeanR(n_r)
               rhoprime(n_r) = osq4pi*ThExpNb*alpha0(n_r)*( -rho0(n_r)* &
                  &            temp0(n_r)*real(s(1,n_r))+ViscHeatFac*   &
                  &            ogrun*real(p(1,n_r)) )
            end do
         end if

         !-- Evaluate nusselt numbers (boundary heat flux density):
         toppres=osq4pi*real(p(1,n_r_cmb))
         botpres=osq4pi*real(p(1,n_r_icb))
         if ( topcond/=0.0_cp ) then

            if ( l_temperature_diff ) then

               botnuss=-osq4pi/botcond*temp0(n_r_icb)*( dLtemp0(n_r_icb)*   &
                 &      real(s(1,n_r_icb)) + real(ds(1,n_r_icb)) +          &
                 &      ViscHeatFac*ThExpNb*alpha0(n_r_icb)*orho1(n_r_icb)*(&
                 &   ( dLalpha0(n_r_icb)+dLtemp0(n_r_icb)-beta(n_r_icb) )*  &  
                 &      real(p(1,n_r_icb)) + real(dp(1,n_r_icb)) ) ) / lScale
               topnuss=-osq4pi/topcond*temp0(n_r_cmb)*( dLtemp0(n_r_cmb)*   &
                 &      real(s(1,n_r_cmb)) + real(ds(1,n_r_cmb)) +          &
                 &      ViscHeatFac*ThExpNb*alpha0(n_r_cmb)*orho1(n_r_cmb)*(&  
                 &   ( dLalpha0(n_r_cmb)+dLtemp0(n_r_cmb)-beta(n_r_cmb) )*  &  
                 &      real(p(1,n_r_cmb)) + real(dp(1,n_r_cmb)) ) ) / lScale

               botflux=four*pi*r_icb**2*kappa(n_r_icb)*rho0(n_r_icb) *      &
                 &     botnuss*botcond*lScale
               topflux=four*pi*r_cmb**2*kappa(n_r_cmb)*rho0(n_r_cmb) *      &
                 &     topnuss*topcond*lScale

               botentropy=osq4pi*real(s(1,n_r_icb))
               topentropy=osq4pi*real(s(1,n_r_cmb))

               bottemp   =temp0(n_r_icb)*botentropy+ViscHeatFac*ThExpNb*   &
                 &        orho1(n_r_icb)*temp0(n_r_icb)*alpha0(n_r_icb)*   &
                 &        botpres
               toptemp   =temp0(n_r_cmb)*topentropy+ViscHeatFac*ThExpNb*   &
                 &        orho1(n_r_cmb)*temp0(n_r_cmb)*alpha0(n_r_cmb)*   &
                 &        toppres
               deltanuss = deltacond/(bottemp-toptemp)
            else ! entropy diffusion

               if ( l_anelastic_liquid ) then
                  botnuss=-osq4pi/botcond*real(ds(1,n_r_icb))/lScale
                  topnuss=-osq4pi/topcond*real(ds(1,n_r_cmb))/lScale
                  botflux=-rho0(n_r_max)*real(ds(1,n_r_max))*osq4pi &
                           *r_icb**2*four*pi*kappa(n_r_max)
                  topflux=-rho0(1)*real(ds(1,1))*osq4pi &
                           *r_cmb**2*four*pi*kappa(1)

                  bottemp=osq4pi*real(s(1,n_r_icb))
                  toptemp=osq4pi*real(s(1,n_r_cmb))

                  botentropy=otemp1(n_r_icb)*bottemp-ViscHeatFac*ThExpNb*   &
                 &        orho1(n_r_icb)*alpha0(n_r_icb)*botpres
                  topentropy=otemp1(n_r_cmb)*toptemp-ViscHeatFac*ThExpNb*   &
                 &        orho1(n_r_cmb)**alpha0(n_r_cmb)*toppres
                  deltanuss = deltacond/(bottemp-toptemp)
               else
                  botnuss=-osq4pi/botcond*real(ds(1,n_r_icb))/lScale
                  topnuss=-osq4pi/topcond*real(ds(1,n_r_cmb))/lScale
                  botflux=-rho0(n_r_max)*temp0(n_r_max)*real(ds(1,n_r_max))/lScale* &
                           r_icb**2*sq4pi*kappa(n_r_max)
                  topflux=-rho0(1)*temp0(1)*real(ds(1,1))/lScale*r_cmb**2* &
                           sq4pi*kappa(1)

                  botentropy=osq4pi*real(s(1,n_r_icb))
                  topentropy=osq4pi*real(s(1,n_r_cmb))

                  bottemp   =temp0(n_r_icb)*botentropy+ViscHeatFac*ThExpNb*   &
                    &        orho1(n_r_icb)*temp0(n_r_icb)*alpha0(n_r_icb)*   &
                    &        botpres
                  toptemp   =temp0(n_r_cmb)*topentropy+ViscHeatFac*ThExpNb*   &
                    &        orho1(n_r_cmb)*temp0(n_r_cmb)*alpha0(n_r_cmb)*   &
                    &        toppres
                  if ( botentropy /= topentropy ) then
                     deltanuss = deltacond/(botentropy-topentropy)
                  else
                     deltanuss = one
                  end if
               end if

            end if
         else
            botnuss   =one
            topnuss   =one
            botflux   =0.0_cp
            topflux   =0.0_cp
            bottemp   =0.0_cp
            toptemp   =0.0_cp
            botentropy=0.0_cp
            topentropy=0.0_cp
            deltanuss =one
         end if

         if ( l_chemical_conv ) then
            if ( topxicond/=0.0_cp ) then
               do n_r=1,n_r_max
                  XiMeanR(n_r)  = XiMeanR(n_r)+timePassed*osq4pi*real(xi(1,n_r))
               end do
               topxi=osq4pi*real(xi(1,n_r_cmb))
               botxi=osq4pi*real(xi(1,n_r_icb))
               botsherwood=-osq4pi/botxicond*real(dxi(1,n_r_icb))/lScale
               topsherwood=-osq4pi/topxicond*real(dxi(1,n_r_cmb))/lScale
               deltasherwood = deltaxicond/(botxi-topxi)
            else
               topxi=0.0_cp
               botxi=0.0_cp
               botsherwood=one
               topsherwood=one
               deltasherwood=one
            end if
         else
            topxi=0.0_cp
            botxi=0.0_cp
            botsherwood=one
            topsherwood=one
            deltasherwood=one
         end if

         mass=four*pi*rInt_R(rhoprime*r**2,n_r_max,n_r_max,drx,chebt_oc)
    
         if ( l_save_out ) then
            open(n_heat_file, file=heat_file, status='unknown', position='append')
         end if

         !-- avoid too small number in output
         if ( abs(toppres) <= 1e-11_cp ) toppres=0.0_cp

         write(n_heat_file,'(1P,ES20.12,16ES16.8)')           &
              & time, botnuss, topnuss, deltanuss,            &
              & bottemp, toptemp, botentropy, topentropy,     &
              & botflux, topflux, toppres, mass, topsherwood, &
              & botsherwood, deltasherwood, botxi, topxi

         if ( l_save_out ) close(n_heat_file)

         if ( l_stop_time ) then
            SMeanR(:)=SMeanR(:)/timeNorm
            TMeanR(:)=TMeanR(:)/timeNorm
            PMeanR(:)=PMeanR(:)/timeNorm
            XiMeanR(:)=XiMeanR(:)/timeNorm

            rhoPrime(:)=ThExpNb*alpha0(:)*(-rho0(:)*temp0(:)*SMeanR(:)+ &
               &         ViscHeatFac*ogrun*PMeanR(:) )

            filename='heatR.'//tag
            open(newunit=filehandle, file=filename, status='unknown')
            do n_r=1,n_r_max
               write(filehandle, '(ES20.10,5ES15.7)' ) &
               &      r(n_r),SMeanR(n_r),TMeanR(n_r),  &
               &      PMeanR(n_r),rhoprime(n_r),XiMeanR(n_r)
            end do

            close(filehandle)
         end if

      end if ! rank == 0
    
   end subroutine outHeat
!---------------------------------------------------------------------------
end module outMisc_mod
