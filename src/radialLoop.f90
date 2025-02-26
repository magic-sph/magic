module radialLoop

   use precision_mod
   use mem_alloc, only: memWrite, bytes_allocated
   use truncation, only: lm_max, lm_maxMag, l_max, l_maxMag
   use radial_data,only: nRstart, nRstop, nRstartMag, nRstopMag
   use time_schemes, only: type_tscheme
   use rIteration, only: rIter_t
   use rIter_mod, only: rIter_single_t

   implicit none

   private

   public :: initialize_radialLoop, finalize_radialLoop, radialLoopG

   class(rIter_t), pointer :: rIter

contains

   subroutine initialize_radialLoop

      integer(lip) :: local_bytes_used

      local_bytes_used = bytes_allocated
      allocate( rIter_single_t :: rIter )
      call rIter%initialize()
      local_bytes_used = bytes_allocated-local_bytes_used

      call memWrite('radialLoop.f90', local_bytes_used)

   end subroutine initialize_radialLoop
!----------------------------------------------------------------------------
   subroutine finalize_radialLoop

      call rIter%finalize()
      deallocate(rIter)

   end subroutine finalize_radialLoop
!----------------------------------------------------------------------------
   subroutine radialLoopG(l_graph,l_frame,time,timeStage,tscheme,dtLast,   &
              &          lTOCalc,lTONext,lTONext2,lHelCalc,lPowerCalc,     & !lMagHelCalc
              &          lRmsCalc,lPressCalc,lPressNext,lViscBcCalc,       &
              &          lFluxProfCalc,lPerpParCalc,lGeosCalc,lHemiCalc,   &
              &          l_probe_out,dsdt,dwdt,dzdt,dpdt,dxidt,dphidt,dbdt,&
              &          djdt,dVxVhLM,dVxBhLM,dVSrLM,dVXirLM,              &
              &          lorentz_torque_ic,lorentz_torque_ma,br_vt_lm_cmb, &
              &          br_vp_lm_cmb,br_vt_lm_icb,br_vp_lm_icb,dtrkc,dthkc)
!======= !PNS
!    subroutine radialLoopG(l_graph,l_frame,time,timeStage,tscheme,dtLast, &
!               &          lTOCalc,lTONext,lTONext2,lHelCalc,              &
!               &          lMagHelCalc,lPowerCalc,                         &
!               &          lRmsCalc,lPressCalc,lPressNext,lViscBcCalc,     &
!               &          lFluxProfCalc,lPerpParCalc,l_probe_out,dsdt,    &
!               &          dwdt,dzdt,dpdt,dxidt,dbdt,djdt,dVxVhLM,dVxBhLM, &
!               &          dVSrLM,dVXirLM,lorentz_torque_ic,               &
!               &          lorentz_torque_ma,br_vt_lm_cmb,br_vp_lm_cmb,    &
!               &          br_vt_lm_icb,br_vp_lm_icb,                      &
!               &          HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,viscLMr,      &
!               &          magHelLMr, uhLMr,                               &
!               &          duhLMr,gradsLMr,fconvLMr,fkinLMr,fviscLMr,      &
!               &          fpoynLMr,fresLMr,EperpLMr,EparLMr,              &
!               &          EperpaxiLMr,EparaxiLMr,dtrkc,dthkc)
! >>>>>>> pns
      !
      !  This subroutine performs the actual time-stepping.
      !

      !--- Input of variables:
      logical,             intent(in) :: l_graph,l_frame
      logical,             intent(in) :: lTOcalc,lTONext,lTONext2,lHelCalc
      logical,             intent(in) :: lPowerCalc,lGeosCalc
!      logical,             intent(in) :: lMagHelCalc !PNS
      logical,             intent(in) :: lViscBcCalc,lFluxProfCalc,lPerpParCalc
      logical,             intent(in) :: lRmsCalc,lHemiCalc
      logical,             intent(in) :: l_probe_out
      logical,             intent(in) :: lPressCalc
      logical,             intent(in) :: lPressNext
      real(cp),            intent(in) :: time,timeStage,dtLast
      class(type_tscheme), intent(in) :: tscheme

      !---- Output of explicit time step:
      !---- dVSrLM and dVxBhLM are output of contributions to explicit time step that
      !     need a further treatment (radial derivatives required):
      complex(cp), intent(out) :: dwdt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dzdt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dpdt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dsdt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dxidt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dphidt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dVSrLM(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dVXirLM(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dbdt(lm_maxMag,nRstartMag:nRstopMag)
      complex(cp), intent(out) :: djdt(lm_maxMag,nRstartMag:nRstopMag)
      complex(cp), intent(out) :: dVxVhLM(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dVxBhLM(lm_maxMag,nRstartMag:nRstopMag)

      ! PNS
      ! !---- Output for axisymmetric helicity:
      ! real(cp),    intent(out) :: HelLMr(l_max+1,nRstart:nRstop)
      ! real(cp),    intent(out) :: Hel2LMr(l_max+1,nRstart:nRstop)
      ! real(cp),    intent(out) :: HelnaLMr(l_max+1,nRstart:nRstop)
      ! real(cp),    intent(out) :: Helna2LMr(l_max+1,nRstart:nRstop)
      ! real(cp),    intent(out) :: uhLMr(l_max+1,nRstart:nRstop)
      ! real(cp),    intent(out) :: duhLMr(l_max+1,nRstart:nRstop)
      ! real(cp),    intent(out) :: viscLMr(l_max+1,nRstart:nRstop)
      ! real(cp),    intent(out) :: magHelLMr(l_max+1,nRstart:nRstop)
      ! real(cp),    intent(out) :: gradsLMr(l_max+1,nRstart:nRstop)
      ! real(cp),    intent(out) :: fkinLMr(l_max+1,nRstart:nRstop)
      ! real(cp),    intent(out) :: fconvLMr(l_max+1,nRstart:nRstop)
      ! real(cp),    intent(out) :: fviscLMr(l_max+1,nRstart:nRstop)
      ! real(cp),    intent(out) :: fresLMr(l_maxMag+1,nRstartMag:nRstopMag)
      ! real(cp),    intent(out) :: fpoynLMr(l_maxMag+1,nRstartMag:nRstopMag)
      ! real(cp),    intent(out) :: EperpLMr(l_max+1,nRstart:nRstop)
      ! real(cp),    intent(out) :: EparLMr(l_max+1,nRstart:nRstop)
      ! real(cp),    intent(out) :: EperpaxiLMr(l_max+1,nRstart:nRstop)
      ! real(cp),    intent(out) :: EparaxiLMr(l_max+1,nRstart:nRstop)

      !---- Output of nonlinear products for nonlinear
      !     magnetic boundary conditions (needed in s_updateB.f):
      complex(cp), intent(out) :: br_vt_lm_cmb(lm_max) ! product br*vt at CMB
      complex(cp), intent(out) :: br_vp_lm_cmb(lm_max) ! product br*vp at CMB
      complex(cp), intent(out) :: br_vt_lm_icb(lm_max) ! product br*vt at ICB
      complex(cp), intent(out) :: br_vp_lm_icb(lm_max) ! product br*vp at ICB
      real(cp),    intent(out) :: lorentz_torque_ma,lorentz_torque_ic

      !---- Output for Courant criteria:
      real(cp),    intent(out) :: dtrkc(nRstart:nRstop),dthkc(nRstart:nRstop)

      call rIter%radialLoop(l_graph,l_frame,time,timeStage,tscheme,dtLast,      &
           &             lTOCalc,lTONext,lTONext2,lHelCalc,lPowerCalc,          &
           &             lRmsCalc,lPressCalc,lPressNext,lViscBcCalc,            &
           &             lFluxProfCalc,lPerpParCalc,lGeosCalc,lHemiCalc,        &
           &             l_probe_out,dsdt,dwdt,dzdt,dpdt,dxidt,dphidt,dbdt,djdt,&
           &             dVxVhLM,dVxBhLM,dVSrLM,dVXirLM,lorentz_torque_ic,      &
           &             lorentz_torque_ma, br_vt_lm_cmb,br_vp_lm_cmb,          &
           &             br_vt_lm_icb,br_vp_lm_icb,dtrkc,dthkc)

! ======= PNS
!       !--- Local variables:
!       integer :: nR,nr_Mag,nBc,lm
!       !integer :: nTheta,nThetaB,nThetaLast
!       integer :: nThetaStart!,nThetaStop
!       logical :: lDeriv,lMagNlBc
!       logical :: lGraphHeader    ! Write header into graph file

!       PERFON('rloop')
!       !LIKWID_ON('rloop')
!       lGraphHeader=l_graph
!       if ( lGraphHeader ) then
! #ifdef WITH_MPI
!          call graphOut_mpi_header(time,nR,nThetaStart,sizeThetaB)
! #else
!          call graphOut_header(time)
! #endif
!       end if

!       if ( rank == 0 ) then
!          dtrkc(n_r_cmb)=1.e10_cp
!          dthkc(n_r_cmb)=1.e10_cp
!       elseif (rank == n_procs-1) then
!          dtrkc(n_r_icb)=1.e10_cp
!          dthkc(n_r_icb)=1.e10_cp
!       end if

!       !------ Set nonlinear terms that are possibly needed at the boundaries.
!       !       They may be overwritten by get_td later.
!       do lm=1,lm_max
!          if ( rank == 0 ) then
!             dVSrLM(lm,n_r_cmb) =zero
!             if ( l_chemical_conv ) dVXirLM(lm,n_r_cmb)=zero
!             if ( l_mag ) dVxBhLM(lm,n_r_cmb)=zero
!             if ( l_double_curl ) dVxVhLM(lm,n_r_cmb)=zero
!          elseif (rank == n_procs-1) then
!             dVSrLM(lm,n_r_icb) =zero
!             if ( l_chemical_conv ) dVXirLM(lm,n_r_icb)=zero
!             if ( l_mag ) dVxBhLM(lm,n_r_icb)=zero
!             if ( l_double_curl ) dVxVhLM(lm,n_r_icb)=zero
!          end if
!       end do

!       !------ Having to calculate non-linear boundary terms?
!       lMagNlBc=.false.
!       if ( ( l_mag_nl .or. l_mag_kin ) .and.                          &
!            &       ( ktopv == 1 .or. l_cond_ma .or.                   &
!            &          ( ktopv == 2 .and. l_rot_ma ) ) .or.            &
!            &       ( kbotv == 1 .or. l_cond_ic .or.                   &
!            &          ( kbotv == 2 .and. l_rot_ic ) ) )               &
!            &     lMagNlBc=.true.

!       !--- Start the big do loop over the radial threads:

!       !nThreadsRmax=1
!       do nR=nRstart,nRstop
!          !IF( nTh > nThreadsRmax ) nThreadsRmax=nTh
!          if ( lVerbose ) then
!             write(*,'(/," ! Starting radial level ",i4)') nR
!          end if

!          !nR = nRC
!          nBc = 0
!          lDeriv = .true.

!          if ( nR == n_r_cmb ) then
!             !nR  = n_r_cmb
!             nBc = ktopv
!             lDeriv= lTOCalc .or. lHelCalc .or. l_frame .or. lPerpParCalc   &
!             &       .or. lViscBcCalc .or. lFluxProfCalc .or. lRmsCalc .or. &
!             &       lPowerCalc
!          else if ( nR == n_r_icb ) then
!             nBc = kbotv
!             lDeriv= lTOCalc .or. lHelCalc .or. l_frame  .or. lPerpParCalc  &
!             &       .or. lViscBcCalc .or. lFluxProfCalc .or. lRmsCalc .or. &
!             &       lPowerCalc
!          end if

!          if ( l_mag .or. l_mag_LF ) then
!             nR_Mag=nR
!          else
!             nR_Mag=1
!          end if

!          call this_rIteration%set_steering_variables(lTOCalc,lTOnext,     &
!               & lTOnext2,lDeriv,lRmsCalc,lHelCalc, lMagHelCalc,           &
!               & lPowerCalc,l_frame,                                       &
!               & lMagNlBc,l_graph,lViscBcCalc,lFluxProfCalc,lPerpParCalc,  &
!               & lPressCalc, lPressNext, l_probe_out)

!          call this_rIteration%do_iteration(nR,nBc,time,timeStage,tscheme,dtLast,&
!               & dsdt(:,nR),dwdt(:,nR),dzdt(:,nR),dpdt(:,nR),dxidt(:,nR),        &
!               & dbdt(:,nR_Mag),djdt(:,nR_Mag),dVxVhLM(:,nR),dVxBhLM(:,nR_Mag),  &
!               & dVSrLM(:,nR),dVXirLM(:,nR),br_vt_lm_cmb,                        &
!               & br_vp_lm_cmb,br_vt_lm_icb,br_vp_lm_icb,lorentz_torque_ic,       &
!               & lorentz_torque_ma,HelLMr(:,nR),Hel2LMr(:,nR),HelnaLMr(:,nR),    &
!               & Helna2LMr(:,nR),viscLMr(:,nR),magHelLMr(:,nR),                  &
!               & uhLMr(:,nR),duhLMr(:,nR),                                       &
!               & gradsLMr(:,nR),fconvLMr(:,nR),fkinLMr(:,nR),fviscLMr(:,nR),     &
!               & fpoynLMr(:,nR_Mag),fresLMr(:,nR_Mag),EperpLMr(:,nR),            &
!               & EparLMr(:,nR),EperpaxiLMr(:,nR),EparaxiLMr(:,nR))

!          dtrkc(nR)=this_rIteration%dtrkc
!          dthkc(nR)=this_rIteration%dthkc

!       end do    ! Loop over radial levels

!       !----- Correct sign of mantel Lorentz torque (see above):
!       lorentz_torque_ma=-lorentz_torque_ma

!       !LIKWID_OFF('rloop')
!       PERFOFF
! >>>>>>> pns
   end subroutine radialLoopG
!----------------------------------------------------------------------------
end module radialLoop
