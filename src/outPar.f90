!$Id$
module outPar_mod

   use parallel_mod
   use truncation, only: n_r_max, n_r_maxMag, l_max, lm_max, &
                         l_maxMag
   use blocking, only: nfs, nThetaBs, sizeThetaB, lm2m
   use logic, only: l_viscBcCalc, l_anel, l_fluxProfs, l_mag_nl, &
                    l_perpPar, l_save_out
   use horizontal_data, only: gauss
   use fields, only: s_Rloc, ds_Rloc
   use physical_parameters, only: ek,prmag,OhmLossFac,ViscHeatFac,opr
   use const, only: pi,mass
   use radial_functions, only: r, or2, sigma, rho0, kappa, temp0, &
                               dr_fac, i_costf_init, d_costf_init
   use radial_data, only: n_r_icb, nRstart, nRstop, nRstartMag, &
                          nRstopMag

   use num_param, only: tScale
   use output_data, only: tag, n_perpPar_file, perpPar_file
   use useful, only: cc2real
   use integration, only: rInt
   use legendre_spec_to_grid, only: lmAS2pt

   implicit none

   private

   real(kind=8), allocatable :: dlVMeanR(:),dlVcMeanR(:)
   real(kind=8), allocatable :: dlVu2MeanR(:),dlVu2cMeanR(:)
   real(kind=8), allocatable :: RolMeanR(:),RolMeanRu2(:),RmMeanR(:)
   real(kind=8), allocatable :: sMeanR(:),Svar(:),Mvar(:)
   real(kind=8), allocatable :: uhMeanR(:),duhMeanR(:)
   real(kind=8), allocatable :: gradT2MeanR(:)
   real(kind=8), allocatable :: fcondMeanR(:),fconvMeanR(:),fkinMeanR(:)
   real(kind=8), allocatable :: fviscMeanR(:)
   real(kind=8), allocatable :: fresMeanR(:), fpoynMeanR(:)
   real(kind=8), allocatable :: EperpMeanR(:),EparMeanR(:)
   real(kind=8), allocatable :: EperpaxiMeanR(:),EparaxiMeanR(:)

   public initialize_outPar_mod, outPar, outPerpPar

contains
   subroutine initialize_outPar_mod

      allocate( dlVMeanR(n_r_max),dlVcMeanR(n_r_max) )
      allocate( dlVu2MeanR(n_r_max),dlVu2cMeanR(n_r_max) )
      allocate( RolMeanR(n_r_max),RolMeanRu2(n_r_max),RmMeanR(n_r_max) )

      dlVMeanR(:)     =0.D0
      dlVcMeanR(:)    =0.D0
      dlVu2MeanR(:)   =0.D0
      dlVu2cMeanR(:)  =0.D0
      RolMeanR(:)     =0.d0
      RolMeanRu2(:)   =0.d0
      RmMeanR(:)      =0.d0

      if ( l_viscBcCalc ) then
         allocate( sMeanR(n_r_max),Svar(nRstart:nRstop),Mvar(nRstart:nRstop) )
         allocate( uhMeanR(n_r_max),duhMeanR(n_r_max),gradT2MeanR(n_r_max) )
         sMeanR(:)       =0.d0
         uhMeanR(:)      =0.d0
         duhMeanR(:)     =0.d0
         gradT2MeanR(:)  =0.d0
      end if

      if ( l_fluxProfs ) then
         allocate( fcondMeanR(n_r_max),fconvMeanR(n_r_max),fkinMeanR(n_r_max) )
         allocate( fviscMeanR(n_r_max) )
         fcondMeanR(:)   =0.d0
         fconvMeanR(:)   =0.d0
         fkinMeanR(:)    =0.d0
         fviscMeanR(:)   =0.d0
         if ( l_mag_nl ) then
            allocate( fresMeanR(n_r_max),fpoynMeanR(n_r_max) )
            fresMeanR(:)    =0.D0
            fpoynMeanR(:)   =0.D0
         end if
      end if

      if ( l_perpPar ) then
         allocate( EperpMeanR(n_r_max),EparMeanR(n_r_max) )
         allocate( EperpaxiMeanR(n_r_max),EparaxiMeanR(n_r_max) )

         EperpMeanR(:)   =0.D0
         EparMeanR(:)    =0.D0
         EperpaxiMeanR(:)=0.D0
         EparaxiMeanR(:) =0.D0
      end if

   end subroutine initialize_outPar_mod
!-----------------------------------------------------------------------
   subroutine outPar(timePassed,timeNorm,nLogs,l_stop_time,       &
                     ekinR,RolRu2,dlVR,dlVRc,dlVRu2,dlVRu2c,      &
                     uhLMr,duhLMr,gradsLMr,fconvLMr,fkinLMr,      &
                     fviscLMr,fpoynLMr,fresLMr,RmR)

      !--- Input of variables
      real(kind=8), intent(in) :: timePassed,timeNorm
      logical,      intent(in) :: l_stop_time
      integer,      intent(in) :: nLogs
      real(kind=8), intent(in) :: RolRu2(n_r_max),dlVRu2(n_r_max),dlVRu2c(n_r_max)
      real(kind=8), intent(in) :: dlVR(n_r_max),dlVRc(n_r_max)
      real(kind=8), intent(in) :: ekinR(n_r_max)     ! kinetic energy w radius
      real(kind=8), intent(in) :: uhLMr(l_max+1,nRstart:nRstop)
      real(kind=8), intent(in) :: duhLMr(l_max+1,nRstart:nRstop)
      real(kind=8), intent(in) :: gradsLMr(l_max+1,nRstart:nRstop)
      real(kind=8), intent(in) :: fkinLMr(l_max+1,nRstart:nRstop)
      real(kind=8), intent(in) :: fconvLMr(l_max+1,nRstart:nRstop)
      real(kind=8), intent(in) :: fviscLMr(l_max+1,nRstart:nRstop)
      real(kind=8), intent(in) :: fpoynLMr(l_maxMag+1,nRstartMag:nRstopMag)
      real(kind=8), intent(in) :: fresLMr(l_maxMag+1,nRstartMag:nRstopMag)

      !--- Output of variables
      real(kind=8), intent(out):: RmR(n_r_max)

      !-- Local variables
      integer :: nR,n,m,lm
      real(kind=8) :: ReR(n_r_max), RoR(n_r_max), RolR(n_r_max)
      character(len=76) :: filename
      integer :: nTheta,nThetaStart,nThetaBlock,nThetaNHS
      real(kind=8) :: duhR(nRstart:nRstop), uhR(nRstart:nRstop)
      real(kind=8) :: gradT2R(nRstart:nRstop), sR(nRstart:nRstop), sR2(nRstart:nRstop)
      real(kind=8) :: fkinR(nRstart:nRstop), fcR(nRstart:nRstop)
      real(kind=8) :: fconvR(nRstart:nRstop), fviscR(nRstart:nRstop)
      real(kind=8) :: fresR(nRstartMag:nRstopMag),fpoynR(nRstartMag:nRstopMag)
      real(kind=8) :: duhR_global(n_r_max), uhR_global(n_r_max)
      real(kind=8) :: gradT2R_global(n_r_max), sR_global(n_r_max)
      real(kind=8) :: Svar_global(n_r_max)
      real(kind=8) :: fkinR_global(n_r_max), fcR_global(n_r_max)
      real(kind=8) :: fconvR_global(n_r_max), fviscR_global(n_r_max)
      real(kind=8) :: fresR_global(n_r_maxMag), fpoynR_global(n_r_maxMag)
      real(kind=8) :: duh(nfs), uh(nfs), gradT2(nfs)
      real(kind=8) :: fkin(nfs), fconv(nfs), fvisc(nfs), fres(nfs), fpoyn(nfs)

      integer :: i,sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1)


      if ( l_viscBcCalc ) then
         do nR=nRstart,nRstop
            sR(nR) = real(s_Rloc(1,nR))
            ! calculate entropy/temperature variance:
            sR2(nR)=0.D0
            do lm=1,lm_max
              m=lm2m(lm)
              sR2(nR)=sR2(nR)+cc2real(s_Rloc(lm,nR),m)
            end do
            if ( nLogs  <=  1) then
               Mvar(nR)=sR(nR)
               Svar(nR)=sR2(nR)-sR(nR)**2
            else
               Mvar(nR)=Mvar(nR)+(sR(nR)-Mvar(nR))/nLogs
               Svar(nR)=Svar(nR)+(sR2(nR)-Mvar(nR)**2)
            end if
         end do

         do nR=nRstart,nRstop
            uhR(nR) =0.d0
            gradT2R(nR)=0.d0
            duhR(nR)=0.d0
            do n=1,nThetaBs ! Loop over theta blocks
               nTheta=(n-1)*sizeThetaB
               nThetaStart=nTheta+1
               call lmAS2pt(duhLMr(1,nR),duh,nThetaStart,sizeThetaB)
               call lmAS2pt(uhLMr(1,nR),uh,nThetaStart,sizeThetaB)
               call lmAS2pt(gradsLMr(1,nR),gradT2,nThetaStart,sizeThetaB)
               do nThetaBlock=1,sizeThetaB
                  nTheta=nTheta+1
                  nThetaNHS=(nTheta+1)/2
                  duhR(nR)=duhR(nR)+gauss(nThetaNHS)*duh(nThetaBlock)
                  uhR(nR) =uhR(nR) +gauss(nThetaNHS)* uh(nThetaBlock)
                  gradT2R(nR)=gradT2R(nR)+gauss(nThetaNHS)*gradT2(nThetaBlock)
               end do
            end do
         end do
         duhR=0.5d0*duhR ! Normalisation for the theta integration
         uhR =0.5d0* uhR ! Normalisation for the theta integration
         gradT2R =0.5d0*gradT2R ! Normalisation for the theta integration

         sendcount  = (nRstop-nRstart+1)
         recvcounts = nr_per_rank
         recvcounts(n_procs-1) = (nr_per_rank+1)
         do i=0,n_procs-1
            displs(i) = i*nr_per_rank
         end do
         call MPI_GatherV(duhR,sendcount,MPI_DOUBLE_PRECISION,              &
             &           duhR_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
             &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(uhR,sendcount,MPI_DOUBLE_PRECISION,              &
             &           uhR_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
             &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(gradT2R,sendcount,MPI_DOUBLE_PRECISION,              &
             &           gradT2R_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
             &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(sR,sendcount,MPI_DOUBLE_PRECISION,              &
             &           sR_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
             &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(Svar,sendcount,MPI_DOUBLE_PRECISION,              &
             &           Svar_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
             &           0,MPI_COMM_WORLD,ierr)
      end if

      if ( l_fluxProfs ) then
         do nR=nRstart,nRstop
            fcR(nR)=-real(ds_Rloc(1,nR))*kappa(nR)*rho0(nR)* &
                     temp0(nR)*r(nR)*r(nR)*dsqrt(4.D0*pi)
         end do
         do nR=nRstart,nRstop
            fkinR(nR) =0.d0
            fconvR(nR)=0.d0
            fviscR(nR)=0.d0
            do n=1,nThetaBs ! Loop over theta blocks
               nTheta=(n-1)*sizeThetaB
               nThetaStart=nTheta+1
               call lmAS2pt(fkinLMr(1,nR),fkin,nThetaStart,sizeThetaB)
               call lmAS2pt(fconvLMr(1,nR),fconv,nThetaStart,sizeThetaB)
               call lmAS2pt(fviscLMr(1,nR),fvisc,nThetaStart,sizeThetaB)
               do nThetaBlock=1,sizeThetaB
                  nTheta=nTheta+1
                  nThetaNHS=(nTheta+1)/2
                  fkinR(nR) =fkinR(nR) +gauss(nThetaNHS)* fkin(nThetaBlock)
                  fconvR(nR)=fconvR(nR)+gauss(nThetaNHS)*fconv(nThetaBlock)
                  fviscR(nR)=fviscR(nR)+gauss(nThetaNHS)*fvisc(nThetaBlock)
               end do
            end do
         end do

         if ( l_mag_nl ) then
            do nR=nRstart,nRstop
               fresR(nR) =0.d0
               fpoynR(nR)=0.d0
               do n=1,nThetaBs ! Loop over theta blocks
                  nTheta=(n-1)*sizeThetaB
                  nThetaStart=nTheta+1
                  call lmAS2pt(fpoynLMr(1,nR),fpoyn,nThetaStart,sizeThetaB)
                  call lmAS2pt(fresLMr(1,nR),fres,nThetaStart,sizeThetaB)
                  do nThetaBlock=1,sizeThetaB
                     nTheta=nTheta+1
                     nThetaNHS=(nTheta+1)/2
                     fpoynR(nR)=fpoynR(nR)+gauss(nThetaNHS)*fpoyn(nThetaBlock)
                     fresR(nR) =fresR(nR) +gauss(nThetaNHS)*fres(nThetaBlock)
                  end do
               end do
            end do
         end if

         sendcount  = (nRstop-nRstart+1)
         recvcounts = nr_per_rank
         recvcounts(n_procs-1) = (nr_per_rank+1)
         do i=0,n_procs-1
            displs(i) = i*nr_per_rank
         end do
         call MPI_GatherV(fkinR,sendcount,MPI_DOUBLE_PRECISION,&
             &           fkinR_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
             &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(fconvR,sendcount,MPI_DOUBLE_PRECISION,&
             &           fconvR_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
             &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(fviscR,sendcount,MPI_DOUBLE_PRECISION,&
             &           fviscR_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
             &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(fcR,sendcount,MPI_DOUBLE_PRECISION,&
             &           fcR_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
             &           0,MPI_COMM_WORLD,ierr)

         if ( l_mag_nl ) then
            call MPI_GatherV(fpoynR,sendcount,MPI_DOUBLE_PRECISION,&
                &           fpoynR_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
                &           0,MPI_COMM_WORLD,ierr)
            call MPI_GatherV(fresR,sendcount,MPI_DOUBLE_PRECISION,&
                &           fresR_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
                &           0,MPI_COMM_WORLD,ierr)
         end if
      end if


      if ( rank == 0 ) then
         do nR=1,n_r_max
            ReR(nR)=SQRT(2.D0*ekinR(nR)*or2(nR)/(4*pi*mass))
            RoR(nR)=ReR(nR)*ek
            if ( dlVR(nR) /= 0d0 ) then
               RolR(nR)=RoR(nR)/dlVR(nR)
            else
               RolR(nR)=RoR(nR)
            end if
            RmR(nR)=ReR(nR)*prmag*sigma(nR)*r(nR)*r(nR)
         end do

         dlVMeanR   =dlVMeanR   +timePassed*dlVR
         dlVcMeanR  =dlVcMeanR  +timePassed*dlVRc
         dlVu2MeanR =dlVu2MeanR +timePassed*dlVRu2
         dlVu2cMeanR=dlVu2cMeanR+timePassed*dlVRu2c
         RolMeanR   =RolMeanR   +timePassed*RolR
         RolMeanRu2 =RolMeanRu2 +timePassed*RolRu2
         RmMeanR    =RmMeanR    +timePassed*RmR*dsqrt(mass/rho0)*or2
         !write(*,"(A,ES20.12)") "dlVcMeanR(n_r_icb) = ",dlVcMeanR(n_r_icb)
         ! this is to get u2 value for RmR(r) to plot in parrad.tag
         ! and also remove r**2, so it has to be volume-averaged 
         ! like RolR
         if ( l_viscBcCalc ) then
            sMeanR     =sMeanR    +timePassed*sR_global
            uhMeanR    =uhMeanR   +timePassed*uhR_global
            duhMeanR   =duhMeanR  +timePassed*duhR_global
            gradT2MeanR=gradT2MeanR+timePassed*gradT2R_global
         end if

         if ( l_fluxProfs ) then
            fkinMeanR =fkinMeanR  +timePassed*fkinR_global
            fcondMeanR=fcondMeanR +timePassed*fcR_global
            fconvMeanR=fconvMeanR +timePassed*fconvR_global
            fviscMeanR=fviscMeanR +timePassed*fviscR_global
            if ( l_mag_nl ) then
               fresMeanR =fresMeanR +timePassed*fresR_global
               fpoynMeanR=fpoynMeanR+timePassed*fpoynR_global
            end if
         end if

         if ( l_stop_time ) then 

            dlVMeanR   =dlVMeanR/timeNorm
            dlVcMeanR  =dlVcMeanR/timeNorm
            RolMeanR   =RolMeanR/timeNorm
            if ( l_anel ) then
               dlVu2MeanR =dlVu2MeanR/timeNorm
               dlVu2cMeanR=dlVu2cMeanR/timeNorm
               RolMeanRu2 =RolMeanRu2/timeNorm
            else
               dlVu2MeanR =dlVMeanR
               dlVu2cMeanR=dlVcMeanR
               RolMeanRu2 =RolMeanR
            end if
            RmMeanR    =RmMeanR/timeNorm

            if ( l_viscBcCalc ) then
               sMeanR     =sMeanR/timeNorm
               Svar_global=Svar_global/(nLogs)
               duhMeanR   =duhMeanR/timeNorm
               uhMeanR    =uhMeanR/timeNorm
               gradT2MeanR=gradT2MeanR/timeNorm
            end if

            if ( l_fluxProfs ) then
               fkinMeanR =ViscHeatFac*fkinMeanR/timeNorm
               fcondMeanR=opr*fcondMeanR/timeNorm
               fconvMeanR=fconvMeanR/timeNorm
               fviscMeanR=ViscHeatFac*fviscMeanR/timeNorm
               if ( l_mag_nl ) then
                  fresMeanR =OhmLossFac*fresMeanR/timeNorm
                  fpoynMeanR=prmag*OhmLossFac*fpoynMeanR/timeNorm
               end if
            end if

            !----- Output into parrad file:
            filename='parR.'//tag
            open(99, file=filename, status='UNKNOWN')
            do nR=1,n_r_max
               write(99,'(D20.10,8D12.4)')       &
                          &   r(nR),             &! 1) radius
                          &   RmMeanR(nR),       &! 2) magnetic Reynolds number
                          &   RolMeanR(nR),      &! 3) local Rossby number
                          &   RolMeanRu2(nR),    &! 4) u squared local Rossby number
                          &   dlVMeanR(nR),      &! 5) local length scale
                          &   dlVcMeanR(nR),     &! 6) conv. local length scale
                          &   dlVu2MeanR(nR),    &! 7) u squared local length scale 
                          &   dlVu2cMeanR(nR)     ! 8) u squared conv. local length scale
            end do
            close(99)

            if ( l_viscBcCalc ) then
               filename='bLayersR.'//tag
               open(99, file=filename, status='UNKNOWN')
               do nR=1,n_r_max
                  write(99,'(D20.10,6D20.12)')           &
                          &   r(nR),                     &! 1) radius
                          &   sMeanR(nR)/SQRT(4.D0*pi),  &! 2) entropy
                          &   Svar_global(nR)/(4.D0*pi), &! 3) entropy variance
                          &   uhMeanR(nR),               &! 4) uh
                          &   duhMeanR(nR),              &! 5) duh/dr
                          &   gradT2MeanR(nR)             ! 6) (grad T)**2
               end do
               close(99)
            end if

            if ( l_fluxProfs ) then
               filename='fluxesR.'//tag
               open(99, file=filename, status='UNKNOWN')
               do nR=1,n_r_max
                  write(99,'(D20.10,7D20.12)')           &
                          &   r(nR),                     &! 1) radius
                          &   fcondMeanR(nR),            &! 2) Fcond
                          &   fconvMeanR(nR),            &! 3) Fconv
                          &   fkinMeanR(nR),             &! 4) Fkin
                          &   fviscMeanR(nR),            &! 5) Fvisc
                          &   fpoynMeanR(nR),            &! 6) Fpoyn
                          &   fresMeanR(nR)               ! 7) Fres
               end do
               close(99)
            end if

         end if ! l_stop_time ?

      end if ! rank0

   end subroutine outPar
!----------------------------------------------------------------------------
   subroutine outPerpPar(time,timePassed,timeNorm,l_stop_time, &
                 &     EperpLMr,EparLMr,EperpaxiLMr,EparaxiLMr)


      !--- Input of variables
      real(kind=8), intent(in) :: time,timePassed,timeNorm
      logical,      intent(in) :: l_stop_time
      real(kind=8), intent(in) :: EparLMr(l_max+1,nRstart:nRstop)
      real(kind=8), intent(in) :: EperpLMr(l_max+1,nRstart:nRstop)
      real(kind=8), intent(in) :: EparaxiLMr(l_max+1,nRstart:nRstop)
      real(kind=8), intent(in) :: EperpaxiLMr(l_max+1,nRstart:nRstop)

      !--- Local variables
      integer :: nR,n,nTheta,nThetaStart,nThetaBlock,nThetaNHS
      character(len=76) :: filename

      real(kind=8) ::EperpaxiR(nRstart:nRstop), EparaxiR(nRstart:nRstop)
      real(kind=8) :: EperpR(nRstart:nRstop), EparR(nRstart:nRstop)
      real(kind=8) :: EperpR_global(n_r_max), EparR_global(n_r_max)
      real(kind=8) :: EperpaxiR_global(n_r_max), EparaxiR_global(n_r_max)
      real(kind=8) :: Eperp(nfs), Epar(nfs), Eperpaxi(nfs), Eparaxi(nfs)
      real(kind=8) :: EperpT,EparT,EperpaxT,EparaxT

      integer :: i,sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1)

      do nR=nRstart,nRstop
         EperpR(nR)   =0.d0
         EparR(nR)    =0.d0
         EparaxiR(nR) =0.d0
         EperpaxiR(nR)=0.d0
         do n=1,nThetaBs ! Loop over theta blocks
            nTheta=(n-1)*sizeThetaB
            nThetaStart=nTheta+1
            call lmAS2pt(EperpLMr(1,nR),Eperp,nThetaStart,sizeThetaB)
            call lmAS2pt(EparLMr(1,nR),Epar,nThetaStart,sizeThetaB)
            call lmAS2pt(EperpaxiLMr(1,nR),Eperpaxi,nThetaStart,sizeThetaB)
            call lmAS2pt(EparaxiLMr(1,nR),Eparaxi,nThetaStart,sizeThetaB)
            do nThetaBlock=1,sizeThetaB
               nTheta=nTheta+1
               nThetaNHS=(nTheta+1)/2
               EperpR(nR)=EperpR(nR)+gauss(nThetaNHS)*Eperp(nThetaBlock)
               EparR(nR) =EparR(nR) +gauss(nThetaNHS)* Epar(nThetaBlock)
               EperpaxiR(nR)=EperpaxiR(nR)+gauss(nThetaNHS)*Eperpaxi(nThetaBlock)
               EparaxiR(nR)=EparaxiR(nR)+gauss(nThetaNHS)*Eparaxi(nThetaBlock)
            end do
         end do
      end do
      EperpR   =0.5d0*EperpR    ! Normalisation for the theta integration
      EparR    =0.5d0*EparR     ! Normalisation for the theta integration
      EperpaxiR=0.5d0*EperpaxiR ! Normalisation for the theta integration
      EparaxiR =0.5d0*EparaxiR  ! Normalisation for the theta integration

      sendcount  = (nRstop-nRstart+1)
      recvcounts = nr_per_rank
      recvcounts(n_procs-1) = (nr_per_rank+1)
      do i=0,n_procs-1
         displs(i) = i*nr_per_rank
      end do

      call MPI_GatherV(EperpR,sendcount,MPI_DOUBLE_PRECISION,&
          &           EperpR_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
          &           0,MPI_COMM_WORLD,ierr)
      call MPI_GatherV(EparR,sendcount,MPI_DOUBLE_PRECISION,&
          &           EparR_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
          &           0,MPI_COMM_WORLD,ierr)
      call MPI_GatherV(EperpaxiR,sendcount,MPI_DOUBLE_PRECISION,&
          &           EperpaxiR_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
          &           0,MPI_COMM_WORLD,ierr)
      call MPI_GatherV(EparaxiR,sendcount,MPI_DOUBLE_PRECISION,&
          &           EparaxiR_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
          &           0,MPI_COMM_WORLD,ierr)


      if ( rank == 0 ) then
         EperpT  =4.D0*pi*rInt(EperpR_global*r**2,n_r_max,dr_fac, &
                               i_costf_init,d_costf_init)
         EparT   =4.D0*pi*rInt(EparR_global*r**2,n_r_max,dr_fac, &
                               i_costf_init,d_costf_init)
         EperpaxT=4.D0*pi*rInt(EperpaxiR_global*r**2,n_r_max,dr_fac, &
                               i_costf_init,d_costf_init)
         EparaxT =4.D0*pi*rInt(EparaxiR_global*r**2,n_r_max,dr_fac, &
                               i_costf_init,d_costf_init)

         !-- Output
         if ( l_save_out ) then
            open(n_perpPar_file, file=perpPar_file, status='UNKNOWN', position='APPEND')
         end if
         write(n_perpPar_file,'(1P,D20.12,4D16.8)') &
              &  time*tScale,     & ! 1
              &  EperpT,EparT,    & ! 2,3
              &  EperpaxT,EparaxT   ! 4,5
         if ( l_save_out ) close(n_perpPar_file)

         EperpMeanR    =EperpMeanR     +timePassed*EperpR_global
         EparMeanR     =EparMeanR      +timePassed*EparR_global
         EperpaxiMeanR =EperpaxiMeanR  +timePassed*EperpaxiR_global
         EparaxiMeanR  =EparaxiMeanR   +timePassed*EparaxiR_global
         if ( l_stop_time ) then
             EperpMeanR     =EperpMeanR/timeNorm
             EparMeanR      =EparMeanR/timeNorm
             EperpaxiMeanR  =EperpaxiMeanR/timeNorm
             EparaxiMeanR   =EparaxiMeanR/timeNorm
             filename='perpParR.'//tag
             open(99, file=filename, status='UNKNOWN')
             do nR=1,n_r_max
                write(99,'(D20.10,4D20.12)')      &
                           &   r(nR),             &! 1) radius
                           &   EperpMeanR(nR),    &! 2) E perpendicular
                           &   EparMeanR(nR),     &! 3) E parallel
                           &   EperpaxiMeanR(nR), &! 4) E perp (axisymetric)
                           &   EparaxiMeanR(nR)    ! 5) E parallel (axisymetric)
             end do
             close(99)
         end if
      end if

   end subroutine outPerpPar
!----------------------------------------------------------------------------
end module outPar_mod
