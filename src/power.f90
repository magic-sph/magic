module power

   use parallel_mod
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_ic_maxMag, n_r_max, n_r_ic_max, l_max, &
       &                 n_r_maxMag
   use radial_data, only: n_r_icb, n_r_cmb, nRstart, nRstop
   use radial_functions, only: r_cmb, r_icb, r, chebt_oc, chebt_ic, &
       &                       or2, O_r_ic2, drx, lambda, temp0,    &
       &                       O_r_ic, rgrav, r_ic, dr_fac_ic,      &
       &                       alpha0, orho1, otemp1
   use physical_parameters, only: kbotv, ktopv, opm, LFfac, BuoFac, &
       &                          ChemFac, ThExpNb, ViscHeatFac
   use num_param, only: tScale, eScale
   use blocking, only: lo_map, st_map, lmStartB, lmStopB, nfs,      &
       &               nThetaBs, sizeThetaB
   use horizontal_data, only: dLh, gauss
   use logic, only: l_rot_ic, l_SRIC, l_rot_ma, l_SRMA, l_save_out, &
       &            l_conv, l_cond_ic, l_heat, l_mag, l_TP_form,    &
       &            l_chemical_conv, l_anelastic_liquid
   use output_data, only: tag
   use useful, only: cc2real,cc22real
   use LMLoop_data,only: llm,ulm,llmMag,ulmMag
   use integration, only: rInt_R,rIntIC
   use outRot, only: get_viscous_torque
   use constants, only: one, two, half
   use legendre_spec_to_grid, only: lmAS2pt

   implicit none

   private

   real(cp), allocatable :: buoMeanR(:)       ! Buoyancy power (thermal)
   real(cp), allocatable :: buo_chem_MeanR(:) ! Buoyancy power (chemical)
   real(cp), allocatable :: viscHeatMeanR(:)  ! Viscous dissipation
   real(cp), allocatable :: ohmDissR(:)       ! Ohmic dissipation
   integer :: n_power_file
   character(len=72) :: power_file

   public :: initialize_output_power, get_power, finalize_output_power

contains
  
   subroutine initialize_output_power
      !
      ! Memory allocation
      !

      allocate( buoMeanR(n_r_max) )
      allocate( ohmDissR(n_r_max) )
      allocate( viscHeatMeanR(n_r_max) )
      allocate( buo_chem_MeanR(n_r_max) )
      bytes_allocated = bytes_allocated+4*n_r_max*SIZEOF_DEF_REAL

      buoMeanR(:)       = 0.0_cp
      ohmDissR(:)       = 0.0_cp
      viscHeatMeanR(:)  = 0.0_cp
      buo_chem_MeanR(:) = 0.0_cp

      power_file='power.'//tag
      if ( rank == 0 .and. (.not. l_save_out) ) then
         open(newunit=n_power_file, file=power_file, status='new')
      end if

   end subroutine initialize_output_power
!----------------------------------------------------------------------------
   subroutine finalize_output_power

      deallocate( buoMeanR, ohmDissR, viscHeatMeanR, buo_chem_MeanR )

      if ( rank == 0 .and. (.not. l_save_out) ) then
         close(n_power_file)
      end if

   end subroutine finalize_output_power
!----------------------------------------------------------------------------
   subroutine get_power(time,timePassed,timeNorm,l_stop_time, &
              &         omega_IC,omega_MA,lorentz_torque_IC,  &
              &         lorentz_torque_MA,w,ddw,z,dz,s,p,     &
              &         xi,b,ddb,aj,dj,db_ic,ddb_ic,aj_ic,    &
              &         dj_ic,viscLMr,viscDiss,ohmDiss)
      !
      !  This subroutine calculates power and dissipation of
      !  the core/mantle system.
      !  Energy input into the outer core is by buoyancy and
      !  possibly viscous accelarations at the boundaries if
      !  the rotation rates of inner core or mantle are prescribed
      !  and kept fixed.
      !  The losses are due to Ohmic and viscous dissipation.
      !  If inner core and mantel are allowed to change their
      !  rotation rates due to viscous forces this power is not
      !  lost from the system and has to be respected.
      !
      !  The output is written into a file  power.TAG.
      !

      !-- Input of variables:
      logical,     intent(in) :: l_stop_time
      real(cp),    intent(in) :: time,timePassed,timeNorm
      real(cp),    intent(in) :: omega_IC,omega_MA
      real(cp),    intent(in) :: lorentz_torque_IC,lorentz_torque_MA
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(in) :: ddw(llm:ulm,n_r_max)
      complex(cp), intent(in) :: z(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dz(llm:ulm,n_r_max)
      complex(cp), intent(in) :: s(llm:ulm,n_r_max)
      complex(cp), intent(in) :: p(llm:ulm,n_r_max)
      complex(cp), intent(in) :: xi(llm:ulm,n_r_max)
      complex(cp), intent(in) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: ddb(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: aj(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: dj(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: db_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: ddb_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: dj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      real(cp),    intent(in) :: viscLMr(l_max+1,nRstart:nRstop)

      !-- Output:
      real(cp),    intent(out) :: viscDiss,ohmDiss

      !-- local:
      integer :: n_r,lm,l,m,l1m0,n
      integer :: nTheta,nThetaStart,nThetaBlock,nThetaNHS

      real(cp) :: r_ratio
      real(cp) :: viscHeatR(nRstart:nRstop)
      real(cp) :: viscHeatR_global(n_r_max)
      real(cp) :: visc(nfs)
      real(cp) :: curlB2,buoy,curlB2_IC,buoy_chem,viscHeat
      real(cp) :: curlB2_r(n_r_max),curlB2_r_global(n_r_max)
      real(cp) :: buoy_r(n_r_max),buoy_r_global(n_r_max)
      real(cp) :: buoy_chem_r(n_r_max),buoy_chem_r_global(n_r_max)
      real(cp) :: curlB2_rIC(n_r_ic_max),curlB2_rIC_global(n_r_ic_max)
      real(cp) :: viscous_torque_ic,viscous_torque_ma

      complex(cp) :: laplace,Bh

      character(len=76) :: fileName
      character(len=7), save :: marker
      real(cp) :: z10ICB,z10CMB,drz10ICB,drz10CMB
      real(cp) :: powerIC,powerMA
      real(cp) :: powerDiffOld,powerDiffT
      real(cp), save :: powerDiff, eDiffInt,tStart

      logical :: rank_has_l1m0
      integer :: sr_tag, fileHandle
#ifdef WITH_MPI
      integer :: i,sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1)
      integer :: status(MPI_STATUS_SIZE)
#endif

      if ( marker /= 'started' ) then
         tStart   =time
         powerDiff=0.0_cp
         eDiffInt =0.0_cp
      end if

      do n_r=nRstart,nRstop
         viscHeatR(n_r)=0.0_cp
         do n=1,nThetaBs ! Loop over theta blocks
            nTheta=(n-1)*sizeThetaB
            nThetaStart=nTheta+1
            call lmAS2pt(viscLMr(1,n_r),visc,nThetaStart,sizeThetaB)
            do nThetaBlock=1,sizeThetaB
               nTheta=nTheta+1
               nThetaNHS=(nTheta+1)/2
               viscHeatR(n_r)=viscHeatR(n_r)+gauss(nThetaNHS)*visc(nThetaBlock)
            end do
         end do
      end do

      do n_r=1,n_r_max

         !if ( l_conv ) then
         !   curlU2_r(n_r)=0.0_cp
         !   !do lm=2,lm_max
         !   do lm=max(2,llm),ulm
         !      l=lo_map%lm2l(lm)
         !      m=lo_map%lm2m(lm)
         !      laplace=dLh(st_map%lm2(l,m))*or2(n_r)*w(lm,n_r)-ddw(lm,n_r)
         !      curlU2_r(n_r)=     curlU2_r(n_r) + dLh(st_map%lm2(l,m)) * ( &
         !           dLh(st_map%lm2(l,m))*or2(n_r)*cc2real(z(lm,n_r),m)  + &
         !           cc2real(dz(lm,n_r),m) + &
         !           cc2real(laplace,m)    )
         !   end do
         !end if

         if ( l_mag ) then
            curlB2_r(n_r)=0.0_cp
            do lm=max(2,llm),ulm
               l=lo_map%lm2l(lm)
               m=lo_map%lm2m(lm)
               laplace=dLh(st_map%lm2(l,m))*or2(n_r)*b(lm,n_r)-ddb(lm,n_r)
               curlB2_r(n_r)=curlB2_r(n_r) +  &
               &             dLh(st_map%lm2(l,m))*lambda(n_r)*(              &
               &             dLh(st_map%lm2(l,m))*or2(n_r)*                  &
               &             cc2real(aj(lm,n_r),m) + cc2real(dj(lm,n_r),m) + &
               &             cc2real(laplace,m)    )
            end do
         end if

         if ( l_heat ) then
            buoy_r(n_r)=0.0_cp
            if ( l_TP_form ) then
               do lm=max(2,llm),ulm
                  l=lo_map%lm2l(lm)
                  m=lo_map%lm2m(lm)
                  buoy_r(n_r)=buoy_r(n_r) +                                 &
                  &           dLh(st_map%lm2(l,m))*BuoFac*rgrav(n_r)*       &
                  &           ( otemp1(n_r)*cc22real(w(lm,n_r),s(lm,n_r),m) &
                  &           -ViscHeatFac*ThExpNb*alpha0(n_r)*orho1(n_r)*  &
                  &            cc22real(w(lm,n_r),p(lm,n_r),m) )
               end do
            else
               do lm=max(2,llm),ulm
                  l=lo_map%lm2l(lm)
                  m=lo_map%lm2m(lm)
                  buoy_r(n_r)=buoy_r(n_r) +                         &
                  &           dLh(st_map%lm2(l,m))*BuoFac*          &
                  &           rgrav(n_r)*cc22real(w(lm,n_r),s(lm,n_r),m)
               end do
            end if
         end if
         if ( l_chemical_conv ) then
            buoy_chem_r(n_r)=0.0_cp
            do lm=max(2,llm),ulm
               l=lo_map%lm2l(lm)
               m=lo_map%lm2m(lm)
               buoy_chem_r(n_r)=buoy_chem_r(n_r) +                      &
               &                dLh(st_map%lm2(l,m))*ChemFac*           &
               &                rgrav(n_r)*cc22real(w(lm,n_r),xi(lm,n_r),m)
            end do
         end if

      end do    ! radial grid points

#ifdef WITH_MPI
      if ( l_mag ) call MPI_Reduce(curlB2_r, curlB2_r_global, n_r_max,&
                   &               MPI_DEF_REAL, MPI_SUM, 0,          &
                   &               MPI_COMM_WORLD, ierr)
      if ( l_heat ) call MPI_Reduce(buoy_r, buoy_r_global, n_r_max,   &
                    &               MPI_DEF_REAL, MPI_SUM, 0,         &
                    &               MPI_COMM_WORLD, ierr)
      if ( l_chemical_conv ) call MPI_Reduce(buoy_chem_r, buoy_chem_r_global, &
                             &               n_r_max, MPI_DEF_REAL, MPI_SUM,  &
                             &               0, MPI_COMM_WORLD, ierr)
      !if ( l_conv ) call MPI_Reduce(curlU2_r,curlU2_r_global,n_r_max,&
      !     & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      if ( l_conv ) then
         sendcount  = (nRstop-nRstart+1)
         recvcounts = nr_per_rank
         recvcounts(n_procs-1) = (nr_per_rank+1)
         do i=0,n_procs-1
            displs(i) = i*nr_per_rank
         end do
         call MPI_GatherV(viscHeatR, sendcount, MPI_DEF_REAL,                &
              &           viscHeatR_global, recvcounts,displs, MPI_DEF_REAL, &
              &           0, MPI_COMM_WORLD, ierr)
      end if

#else
      !if ( l_conv ) curlU2_r_global=curlU2_r
      if ( l_mag )  curlB2_r_global=curlB2_r
      if ( l_heat ) buoy_r_global  =buoy_r
      if ( l_chemical_conv ) buoy_chem_r_global = buoy_chem_r

      if ( l_conv ) viscHeatR_global(:)=viscHeatR(:)
#endif

      if ( rank == 0 ) then
         !-- Transform to cheb space:
         if ( l_conv ) then
            !curlU2MeanR=curlU2MeanR+timePassed*curlU2_r_global*eScale
            !curlU2=rInt_R(curlU2_r_global,n_r_max,n_r_max,drx,chebt_oc)
            viscHeatMeanR=viscHeatMeanR+timePassed*viscHeatR_global*eScale
            viscHeat=rInt_R(viscHeatR_global,n_r_max,n_r_max,drx,chebt_oc)
            viscHeat=eScale*viscHeat
         else
            viscHeat=0.0_cp
         end if
         if ( l_mag )  then
            ohmDissR=ohmDissR+timePassed*curlB2_r_global*LFfac*opm*eScale
            curlB2=rInt_R(curlB2_r_global,n_r_max,n_r_max,drx,chebt_oc)
            curlB2=LFfac*opm*eScale*curlB2
         else
            curlB2=0.0_cp
         end if
         if ( l_heat ) then
            buoMeanR=buoMeanR+timePassed*buoy_r_global*eScale
            buoy=rInt_R(buoy_r_global,n_r_max,n_r_max,drx,chebt_oc)
            buoy=eScale*buoy
         else
            buoy=0.0_cp
         end if
         if ( l_chemical_conv ) then
            buo_chem_MeanR=buo_chem_MeanR+timePassed*buoy_chem_r_global*eScale
            buoy_chem=rInt_R(buoy_chem_r_global,n_r_max,n_r_max,drx,chebt_oc)
            buoy_chem=eScale*buoy_chem
         else
            buoy_chem=0.0_cp
         end if
      end if

      !-- Inner core:

      if ( l_cond_ic .and. l_mag ) then

         do n_r=1,n_r_ic_max
            r_ratio=r_ic(n_r)/r_ic(1)
            curlB2_rIC(n_r)=0.0_cp
            do lm=max(2,llm),ulm
               l=lo_map%lm2l(lm)
               m=lo_map%lm2m(lm)
               Bh=(l+one)*O_r_ic(n_r)*aj_ic(lm,n_r)+dj_ic(lm,n_r)
               laplace=-ddb_ic(lm,n_r) - two*(l+one)*O_r_ic(n_r)*db_ic(lm,n_r)
               curlB2_rIC(n_r)=curlB2_rIC(n_r) +                          &
               &               dLh(st_map%lm2(l,m))*r_ratio**(2*l+2) *  ( &
               &               dLh(st_map%lm2(l,m))*O_r_ic2(n_r)*         &
               &               cc2real(aj_ic(lm,n_r),m) +                 &
               &               cc2real(Bh,m) + cc2real(laplace,m)       )
            end do
         end do    ! radial grid points

#ifdef WITH_MPI
         call MPI_Reduce(curlB2_rIC, curlB2_rIC_global, n_r_ic_max, &
              &          MPI_DEF_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
#else
         curlB2_rIC_global=curlB2_rIC
#endif

         if ( rank == 0 ) then
            curlB2_IC=rIntIC(curlB2_rIC_global,n_r_ic_max,dr_fac_ic,chebt_ic)
            curlB2_IC=LFfac*opm*eScale*curlB2_IC
         end if
      else
         if ( rank == 0 ) then
            curlB2_IC=0.0_cp
         end if
      end if  ! conducting inner core ?

      !-- Add up and correct:
      !  A correction is neccesaary when using the integral over
      !  (curl U)**2 to calculate the dippisated energy and inner core
      !  or mantle are allowed to rotate.
      !  The only flow component here is l=1,m=0, i.e. lm=2
      !  If the rotation rates of inner core of mantle are kept
      !  fixed, the viscous dissipation may actually be a power source if
      !  the radial derivatives drz10 are negative (positive) at the
      !  ICB (CMB)! The energy transfere is described by the very
      !  correction terms.
      l1m0=lo_map%lm2(1,0)
      rank_has_l1m0=.false. ! set default
      sr_tag=46378 !arbitray send-recv tag
      if ( lmStartB(rank+1) <= l1m0 .and. lmStopB(rank+1) >= l1m0 ) then
         !if (rank /= 0) then
         !   PRINT*,"in get_power, l1m0 is not on rank 0!"
         !   stop
         !end if
         if ( l_rot_IC ) then
            z10ICB  =real(z(l1m0,n_r_ICB))
            drz10ICB=real(dz(l1m0,n_r_ICB))
         else
            z10ICB  =0.0_cp
            drz10ICB=0.0_cp
         end if
         if ( l_rot_MA ) then
            z10CMB  =real(z(l1m0,n_r_CMB))
            drz10CMB=real(dz(l1m0,n_r_CMB))
         else
            z10CMB  =0.0_cp
            drz10CMB=0.0_cp
         end if

#ifdef WITH_MPI
         if ( rank /= 0 ) then
            ! send data to rank 0
            call MPI_Send(z10ICB, 1, MPI_DEF_COMPLEX, 0, sr_tag, &
                 &        MPI_COMM_WORLD, ierr)
            call MPI_Send(drz10ICB, 1, MPI_DEF_COMPLEX, 0, sr_tag+1, &
                 &        MPI_COMM_WORLD, ierr)
            call MPI_Send(z10CMB, 1, MPI_DEF_COMPLEX, 0, sr_tag+2, &
                 &        MPI_COMM_WORLD, ierr)
            call MPI_Send(drz10CMB, 1, MPI_DEF_COMPLEX, 0, sr_tag+3, &
                 &        MPI_COMM_WORLD, ierr)
         end if
#endif
         rank_has_l1m0=.true.
      end if

      if ( rank == 0 ) then
#ifdef WITH_MPI
         if ( .not. rank_has_l1m0 ) then
            ! receive data from the source ranks
            call MPI_Recv(z10ICB, 1, MPI_DEF_COMPLEX, MPI_ANY_SOURCE, &
                 &        sr_tag, MPI_COMM_WORLD, status, ierr)
            call MPI_Recv(drz10ICB, 1, MPI_DEF_COMPLEX, MPI_ANY_SOURCE, &
                 &        sr_tag+1, MPI_COMM_WORLD, status, ierr)
            call MPI_Recv(z10CMB, 1, MPI_DEF_COMPLEX, MPI_ANY_SOURCE, &
                 &        sr_tag+2, MPI_COMM_WORLD, status, ierr)
            call MPI_Recv(drz10CMB, 1, MPI_DEF_COMPLEX, MPI_ANY_SOURCE, &
                 &        sr_tag+3, MPI_COMM_WORLD, status, ierr)
         end if
#endif

         if ( l_conv ) then
            viscDiss= -viscHeat
            if ( l_rot_IC ) viscDiss=viscDiss - two*z10ICB*drz10ICB
            if ( l_rot_MA ) viscDiss=viscDiss + two*z10CMB*drz10CMB
         else
            viscDiss=0.0_cp
         end if

         !--- If the inner core or mantle rotation rate is allowed to change due
         !    to viscous drag, the power transfer has to be taken into account:

         !-- Calculating viscous torques:
         if ( l_rot_ic .and. kbotv == 2 ) then
            call get_viscous_torque(viscous_torque_ic,z10ICB,drz10ICB,r_icb)
         else
            viscous_torque_ic=0.0_cp
         end if
         if ( l_rot_ma .and. ktopv == 2 ) then
            call get_viscous_torque(viscous_torque_ma,z10CMB,drz10CMB,r_cmb)
         else
            viscous_torque_ma=0.0_cp
         end if

         if ( l_rot_IC .and. .not. l_SRIC ) then
            if ( kbotv == 2 ) then
               call get_viscous_torque(viscous_torque_ic,z10ICB,drz10ICB,r_icb)
            else
               viscous_torque_ic=0.0_cp
            end if
            powerIC=omega_IC*(viscous_torque_ic+lorentz_torque_ic)
         else
            powerIC=0.0_cp
         end if
         if ( l_rot_MA ) then
            if ( ktopv == 2 ) then
               call get_viscous_torque(viscous_torque_ma,z10CMB,drz10CMB,r_cmb)
            else
               viscous_torque_ma=0.0_cp
            end if
            powerMA=omega_MA*(-viscous_torque_ma+lorentz_torque_ma)
         else
            powerMA=0.0_cp
         end if

         !--- Because the two systems are coupled only the total 
         !--- ohmic dissipation in useful:
         ohmDiss=-curlB2-curlB2_IC

         powerDiffOld=powerDiff
         powerDiff   =(buoy+buoy_chem+powerIC+powerMA+viscDiss+ohmDiss)

         if ( marker == 'started' ) then
            powerDiffT  =1.5_cp*powerDiff-half*powerDiffOld
            eDiffInt=eDiffInt+timePassed*timePassed*powerDiffT
            if ( l_save_out ) then
               open(newunit=n_power_file, file=power_file, status='unknown', &
               &    position='append')
            end if
            write(n_power_file,'(1P,ES20.12,10ES16.8)')  &
            &    time*tScale, buoy, buoy_chem,           &
            &     -two*z10ICB*drz10ICB,                  &
            &    two*z10CMB*drz10CMB, viscDiss,          &
            &    ohmDiss, powerMA, powerIC, powerDiff,   &
            &    eDiffInt/timeNorm
            if ( l_save_out ) close(n_power_file)
         else
            marker='started'
         end if

         if ( l_stop_time ) then
            buoMeanR(:)=buoMeanR(:)/timeNorm
            buoMeanR(1)      =0.0_cp ! Ensure this is really zero on the boundaries
            buoMeanR(n_r_max)=0.0_cp
            if ( l_chemical_conv ) then
               buo_chem_MeanR(:)=buo_chem_MeanR(:)/timeNorm
               buo_chem_MeanR(1)=0.0_cp
               buo_chem_MeanR(n_r_max)=0.0_cp
            end if
            ohmDissR(:)     =ohmDissR(:)/timeNorm
            viscHeatMeanR(:)=viscHeatMeanR(:)/timeNorm
            fileName='powerR.'//tag
            open(newunit=fileHandle, file=fileName, status='unknown')
            do n_r=1,n_r_max
               write(fileHandle,'(ES20.10,4ES15.7)')  &
                    &   r(n_r),               & ! 1) radius
                    &   buoMeanR(n_r),        & ! 2) Buo power
                    &   buo_chem_MeanR(n_r),  & ! 3) Chem power
                    &   viscHeatMeanR(n_r),   & ! 4) Viscous heating
                    &   ohmDissR(n_r)           ! 5) Ohmic dissipation
            end do
            close(fileHandle)
         end if
      end if

  end subroutine get_power
!----------------------------------------------------------------------------
end module power
