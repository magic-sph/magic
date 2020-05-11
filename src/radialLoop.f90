#include "perflib_preproc.cpp"
module radialLoop

   use precision_mod
   use mem_alloc, only: memWrite, bytes_allocated
   use truncation, only: n_lm_loc, n_lmMag_loc, &
       &                 nRstart,nRstop,n_r_cmb, nRstartMag, nRstopMag,   &
       &                 n_r_icb, n_lmP_loc
   use physical_parameters, only: ktopv, kbotv
   use blocking, only: nThetaBs, sizeThetaB
   use logic, only: l_dtB, l_mag, l_mag_LF, lVerbose, l_rot_ma, l_rot_ic, &
       &            l_cond_ic, l_mag_kin, l_cond_ma, l_mag_nl,            &
       &            l_single_matrix, l_double_curl, l_chemical_conv
   use constants, only: zero
   use parallel_mod, only: n_ranks_r, coord_r, l_master_rank
   use time_schemes, only: type_tscheme
#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif
   use rIteration_mod, only: rIteration_t
   use rIterThetaBlocking_mod, only: rIterThetaBlocking_t
#ifdef WITH_SHTNS
   use rIterThetaBlocking_shtns_mod, only: rIterThetaBlocking_shtns_t
#else
#endif
#ifdef WITH_MPI
   use graphOut_mod, only: graphOut_mpi_header
#else
   use graphOut_mod, only: graphOut_header
#endif

   implicit none

   private

   public :: initialize_radialLoop, finalize_radialLoop, radialLoopG

   class(rIteration_t), pointer :: this_rIteration

contains

   subroutine initialize_radialLoop

      character(len=100) :: this_type
      integer(lip) :: local_bytes_used

      local_bytes_used = bytes_allocated

#ifdef WITH_SHTNS
      allocate( rIterThetaBlocking_shtns_t :: this_rIteration )
#endif
      this_type = this_rIteration%getType()
      if ( l_master_rank ) write(*,"(2A)") " ! Using rIteration type: ",trim(this_type)
      call this_rIteration%initialize()
      select type (this_rIteration)
         class is (rIterThetaBlocking_t)
            !call this_rIteration%set_ThetaBlocking(nThetaBs,sizeThetaB)
         class default
            write(*,*) "This_rIteration has no matching type in radialLoop.f90"
      end select

      local_bytes_used = bytes_allocated-local_bytes_used

      call memWrite('radialLoop.f90', local_bytes_used)

   end subroutine initialize_radialLoop
!----------------------------------------------------------------------------
   subroutine finalize_radialLoop

      call this_rIteration%finalize()
      deallocate(this_rIteration)

   end subroutine finalize_radialLoop
!----------------------------------------------------------------------------
   subroutine radialLoopG(l_graph,l_frame,time,timeStage,tscheme,dtLast, &
              &          lTOCalc,lTONext,lTONext2,lHelCalc,lPowerCalc,   &
              &          lRmsCalc,lPressCalc,lPressNext,lViscBcCalc,     &
              &          lFluxProfCalc,lPerpParCalc,l_probe_out,dsdt,    &
              &          dwdt,dzdt,dpdt,dxidt,dbdt,djdt,dVxVhLM,dVxBhLM, &
              &          dVSrLM,dVXirLM,lorentz_torque_ic,               &
              &          lorentz_torque_ma,br_vt_lm_cmb,br_vp_lm_cmb,    &
              &          br_vt_lm_icb,br_vp_lm_icb,HelASr,Hel2ASr,       &
              &          HelnaASr,Helna2ASr,HelEAASr,viscAS,uhASr,       &
              &          duhASr,gradsASr,fconvASr,fkinASr,fviscASr,      &
              &          fpoynASr,fresASr,EperpASr,EparASr,              &
              &          EperpaxiASr,EparaxiASr,dtrkc,dthkc)
      !
      !  This subroutine performs the actual time-stepping.
      !

      !--- Input of variables:
      logical,             intent(in) :: l_graph,l_frame
      logical,             intent(in) :: lTOcalc,lTONext,lTONext2,lHelCalc
      logical,             intent(in) :: lPowerCalc
      logical,             intent(in) :: lViscBcCalc,lFluxProfCalc,lPerpParCalc
      logical,             intent(in) :: lRmsCalc
      logical,             intent(in) :: l_probe_out
      logical,             intent(in) :: lPressCalc
      logical,             intent(in) :: lPressNext
      real(cp),            intent(in) :: time,timeStage,dtLast
      class(type_tscheme), intent(in) :: tscheme

      !---- Output of explicit time step:
      !---- dVSrLM and dVxBhLM are output of contributions to explicit time step that
      !     need a further treatment (radial derivatives required):
      complex(cp), intent(out) :: dwdt(n_lm_loc,nRstart:nRstop)
      complex(cp), intent(out) :: dzdt(n_lm_loc,nRstart:nRstop)
      complex(cp), intent(out) :: dpdt(n_lm_loc,nRstart:nRstop)
      complex(cp), intent(out) :: dsdt(n_lm_loc,nRstart:nRstop)
      complex(cp), intent(out) :: dxidt(n_lm_loc,nRstart:nRstop)
      complex(cp), intent(out) :: dVSrLM(n_lm_loc,nRstart:nRstop)
      complex(cp), intent(out) :: dVXirLM(n_lm_loc,nRstart:nRstop)
      complex(cp), intent(out) :: dbdt(n_lmMag_loc,nRstartMag:nRstopMag)
      complex(cp), intent(out) :: djdt(n_lmMag_loc,nRstartMag:nRstopMag)
      complex(cp), intent(out) :: dVxVhLM(n_lm_loc,nRstart:nRstop)
      complex(cp), intent(out) :: dVxBhLM(n_lmMag_loc,nRstartMag:nRstopMag)
      real(cp),    intent(out) :: lorentz_torque_ma,lorentz_torque_ic

      !---- Output for axisymmetric helicity:
      real(cp),    intent(out) :: HelASr(2,nRstart:nRstop)
      real(cp),    intent(out) :: Hel2ASr(2,nRstart:nRstop)
      real(cp),    intent(out) :: HelnaASr(2,nRstart:nRstop)
      real(cp),    intent(out) :: Helna2ASr(2,nRstart:nRstop)
      real(cp),    intent(out) :: HelEAASr(nRstart:nRstop)
      real(cp),    intent(out) :: uhASr(nRstart:nRstop)
      real(cp),    intent(out) :: duhASr(nRstart:nRstop)
      real(cp),    intent(out) :: viscAS(nRstart:nRstop)
      real(cp),    intent(out) :: gradsASr(nRstart:nRstop)
      real(cp),    intent(out) :: fkinASr(nRstart:nRstop)
      real(cp),    intent(out) :: fconvASr(nRstart:nRstop)
      real(cp),    intent(out) :: fviscASr(nRstart:nRstop)
      real(cp),    intent(out) :: fresASr(nRstartMag:nRstopMag)
      real(cp),    intent(out) :: fpoynASr(nRstartMag:nRstopMag)
      real(cp),    intent(out) :: EperpASr(nRstart:nRstop)
      real(cp),    intent(out) :: EparASr(nRstart:nRstop)
      real(cp),    intent(out) :: EperpaxiASr(nRstart:nRstop)
      real(cp),    intent(out) :: EparaxiASr(nRstart:nRstop)

      !---- Output of nonlinear products for nonlinear
      !     magnetic boundary conditions (needed in s_updateB.f):
      complex(cp), intent(out) :: br_vt_lm_cmb(n_lmP_loc) ! product br*vt at CMB
      complex(cp), intent(out) :: br_vp_lm_cmb(n_lmP_loc) ! product br*vp at CMB
      complex(cp), intent(out) :: br_vt_lm_icb(n_lmP_loc) ! product br*vt at ICB
      complex(cp), intent(out) :: br_vp_lm_icb(n_lmP_loc) ! product br*vp at ICB

      !---- Output for Courant criteria:
      real(cp),intent(out) :: dtrkc(nRstart:nRstop),dthkc(nRstart:nRstop)


      !--- Local variables:
      integer :: nR,nr_Mag,nBc,lm
      !integer :: nTheta,nThetaB,nThetaLast
      integer :: nThetaStart!,nThetaStop
      logical :: lDeriv,lMagNlBc
      logical :: lGraphHeader    ! Write header into graph file

      PERFON('rloop')
      !LIKWID_ON('rloop')
      lGraphHeader=l_graph
      if ( lGraphHeader ) then
#ifdef WITH_MPI
         call graphOut_mpi_header(time,nR,nThetaStart,sizeThetaB)
#else
         call graphOut_header(time)
#endif
      end if

      if ( coord_r == 0 ) then
         dtrkc(n_r_cmb)=1.e10_cp
         dthkc(n_r_cmb)=1.e10_cp
      elseif (coord_r == n_ranks_r-1) then
         dtrkc(n_r_icb)=1.e10_cp
         dthkc(n_r_icb)=1.e10_cp
      end if

      !------ Set nonlinear terms that are possibly needed at the boundaries.
      !       They may be overwritten by get_td later.
      do lm=1,n_lm_loc
         if ( coord_r == 0 ) then
            dVSrLM(lm,n_r_cmb) =zero
            if ( l_chemical_conv ) dVXirLM(lm,n_r_cmb)=zero
            if ( l_mag ) dVxBhLM(lm,n_r_cmb)=zero
            if ( l_double_curl ) dVxVhLM(lm,n_r_cmb)=zero
         elseif (coord_r == n_ranks_r-1) then
            dVSrLM(lm,n_r_icb) =zero
            if ( l_chemical_conv ) dVXirLM(lm,n_r_icb)=zero
            if ( l_mag ) dVxBhLM(lm,n_r_icb)=zero
            if ( l_double_curl ) dVxVhLM(lm,n_r_icb)=zero
         end if
      end do

      !------ Having to calculate non-linear boundary terms?
      lMagNlBc=.false.
      if ( ( l_mag_nl .or. l_mag_kin ) .and.                          &
           &       ( ktopv == 1 .or. l_cond_ma .or.                   &
           &          ( ktopv == 2 .and. l_rot_ma ) ) .or.            &
           &       ( kbotv == 1 .or. l_cond_ic .or.                   &
           &          ( kbotv == 2 .and. l_rot_ic ) ) )               &
           &     lMagNlBc=.true.

      !--- Start the big do loop over the radial threads:

      !nThreadsRmax=1
      do nR=nRstart,nRstop
         !IF( nTh > nThreadsRmax ) nThreadsRmax=nTh
         if ( lVerbose ) then
            write(*,'(/," ! Starting radial level ",i4)') nR
         end if

         !nR = nRC
         nBc = 0
         lDeriv = .true.

         if ( nR == n_r_cmb ) then
            !nR  = n_r_cmb
            nBc = ktopv
            lDeriv= lTOCalc .or. lHelCalc .or. l_frame .or. lPerpParCalc   &
            &       .or. lViscBcCalc .or. lFluxProfCalc .or. lRmsCalc .or. &
            &       lPowerCalc
         else if ( nR == n_r_icb ) then
            nBc = kbotv
            lDeriv= lTOCalc .or. lHelCalc .or. l_frame  .or. lPerpParCalc  &
            &       .or. lViscBcCalc .or. lFluxProfCalc .or. lRmsCalc .or. &
            &       lPowerCalc
         end if

         if ( l_mag .or. l_mag_LF ) then
            nR_Mag=nR
         else
            nR_Mag=1
         end if

         call this_rIteration%set_steering_variables(lTOCalc,lTOnext,     &
              & lTOnext2,lDeriv,lRmsCalc,lHelCalc,lPowerCalc,l_frame,     &
              & lMagNlBc,l_graph,lViscBcCalc,lFluxProfCalc,lPerpParCalc,  &
              & lPressCalc, lPressNext, l_probe_out)

         call this_rIteration%do_iteration(nR,nBc,time,timeStage,tscheme,dtLast,&
              & dsdt(:,nR),dwdt(:,nR),dzdt(:,nR),dpdt(:,nR),dxidt(:,nR),        &
              & dbdt(:,nR_Mag),djdt(:,nR_Mag),dVxVhLM(:,nR),dVxBhLM(:,nR_Mag),  &
              & dVSrLM(:,nR),dVXirLM(:,nR),br_vt_lm_cmb,                        &
              & br_vp_lm_cmb,br_vt_lm_icb,br_vp_lm_icb,lorentz_torque_ic,       &
              & lorentz_torque_ma,HelASr(:,nR),Hel2ASr(:,nR),HelnaASr(:,nR),    &
              & Helna2ASr(:,nR),HelEAASr(nR),viscAS(nR),uhASr(nR),duhASr(nR),   &
              & gradsASr(nR),fconvASr(nR),fkinASr(nR),fviscASr(nR),             &
              & fpoynASr(nR_Mag),fresASr(nR_Mag),EperpASr(nR),                  &
              & EparASr(nR),EperpaxiASr(nR),EparaxiASr(nR))

         dtrkc(nR)=this_rIteration%dtrkc
         dthkc(nR)=this_rIteration%dthkc
         
      end do    ! Loop over radial levels

      !----- Correct sign of mantel Lorentz torque (see above):
      lorentz_torque_ma=-lorentz_torque_ma
      
      !LIKWID_OFF('rloop')
      PERFOFF
   end subroutine radialLoopG
!----------------------------------------------------------------------------
end module radialLoop
