!$Id$
#include "perflib_preproc.cpp"
module radialLoop

   use truncation, only: lm_max, lm_maxMag, l_max, l_maxMag, lmP_max
   use physical_parameters, only: ktopv, kbotv
   use blocking, only: nThetaBs, sizeThetaB
   use logic, only: l_dtB, l_mag, l_mag_LF, lVerbose, l_rot_ma, l_rot_ic, &
                    l_cond_ic, l_mag_kin, l_cond_ma, l_mag_nl
   use output_data, only: ngform
   use const, only: zero
   use parallel_mod, only: rank, n_procs
   use radial_data,only: nRstart,nRstop,n_r_cmb, nRstartMag, nRstopMag, &
                         n_r_icb
#if (FFTLIB==JW)
   use fft_JW
#elif (FFTLIB==MKL)
   use fft_MKL
#endif
#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif
   use rIteration_mod, only: rIteration_t
   use rIterThetaBlocking_mod, only: rIterThetaBlocking_t
   use rIterThetaBlocking_seq_mod, only: rIterThetaBlocking_seq_t
   use rIterThetaBlocking_OpenMP_mod, only: rIterThetaBlocking_OpenMP_t
   use graphOut_mod, only: graphOut_mpi_header

   implicit none

   private
   !---- Nonlinear terms, field components and help arrays
   !     for Legendre transform, field components, and dtB and 
   !     TO output. These are all stored in COMMON BLOCKS
   !     that are thread privat rather than using explicit 
   !     PRIVATE statements in the PARALLEL do clause.


   public :: initialize_radialLoop,finalize_radialLoop,radialLoopG

   CLASS(rIteration_t), pointer :: this_rIteration

contains

   subroutine initialize_radialLoop

      character(len=100) :: this_type

#ifdef WITHOMP
      allocate( rIterThetaBlocking_OpenMP_t :: this_rIteration )
#else
      allocate( rIterThetaBlocking_seq_t :: this_rIteration )
#endif
      this_type = this_rIteration%getType()
      write(*,"(2A)") "Using rIteration type: ",trim(this_type)
      call this_rIteration%initialize()
      select type (this_rIteration)
      CLASS is (rIterThetaBlocking_t)
         call this_rIteration%set_ThetaBlocking(nThetaBs,sizeThetaB)
      CLASS default
         print*,"this_rIteration has no matching type in m_radialLoop.F90"
      end select

   end subroutine initialize_radialLoop
!----------------------------------------------------------------------------
   subroutine finalize_radialLoop

      call this_rIteration%finalize()
      deallocate(this_rIteration)

   end subroutine finalize_radialLoop
!----------------------------------------------------------------------------
   subroutine radialLoopG(l_graph,l_cour,l_frame,time,dt,dtLast,        &
       &                 lTOCalc,lTONext,lTONext2,lHelCalc,lRmsCalc,    &
       &                 lViscBcCalc,lFluxProfCalc,lPerpParCalc,        &
       &                 dsdt,dwdt,dzdt,dpdt,dbdt,djdt,dVxBhLM,dVSrLM,  &
       &                 lorentz_torque_ic,lorentz_torque_ma,           &
       &                 br_vt_lm_cmb,br_vp_lm_cmb,                     &
       &                 br_vt_lm_icb,br_vp_lm_icb,                     &
       &                 HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,uhLMr,       &
       &                 duhLMr,gradsLMr,fconvLMr,fkinLMr,fviscLMr,     &
       &                 fpoynLMr,fresLMr,EperpLMr,EparLMr,             &
       &                 EperpaxiLMr,EparaxiLMr,dtrkc,dthkc)
      !  +-------------+----------------+------------------------------------+
      !  |                                                                   |
      !  |  This subroutine performs the actual time-stepping.               |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+

      !--- Input of variables:
      logical,      intent(in) :: l_graph,l_cour,l_frame
      logical,      intent(in) :: lTOcalc,lTONext,lTONext2,lHelCalc
      logical,      intent(in) :: lViscBcCalc,lFluxProfCalc,lPerpParCalc
      logical,      intent(in) :: lRmsCalc
      real(kind=8), intent(in) :: time,dt,dtLast

      !---- Output of explicit time step:
      !---- dVSrLM and dVxBhLM are output of contributions to explicit time step that
      !     need a further treatment (radial derivatives required):
      complex(kind=8), intent(out) :: dwdt(lm_max,nRstart:nRstop)
      complex(kind=8), intent(out) :: dzdt(lm_max,nRstart:nRstop)
      complex(kind=8), intent(out) :: dpdt(lm_max,nRstart:nRstop)
      complex(kind=8), intent(out) :: dsdt(lm_max,nRstart:nRstop)
      complex(kind=8), intent(out) :: dVSrLM(lm_max,nRstart:nRstop)
      complex(kind=8), intent(out) :: dbdt(lm_maxMag,nRstartMag:nRstopMag)
      complex(kind=8), intent(out) :: djdt(lm_maxMag,nRstartMag:nRstopMag)
      complex(kind=8), intent(out) :: dVxBhLM(lm_maxMag,nRstartMag:nRstopMag)
      real(kind=8),    intent(out) :: lorentz_torque_ma,lorentz_torque_ic

      !---- Output for axisymmetric helicity:
      real(kind=8),    intent(out) :: HelLMr(l_max+1,nRstart:nRstop)
      real(kind=8),    intent(out) :: Hel2LMr(l_max+1,nRstart:nRstop)
      real(kind=8),    intent(out) :: HelnaLMr(l_max+1,nRstart:nRstop)
      real(kind=8),    intent(out) :: Helna2LMr(l_max+1,nRstart:nRstop)
      real(kind=8),    intent(out) :: uhLMr(l_max+1,nRstart:nRstop)
      real(kind=8),    intent(out) :: duhLMr(l_max+1,nRstart:nRstop)
      real(kind=8),    intent(out) :: gradsLMr(l_max+1,nRstart:nRstop)
      real(kind=8),    intent(out) :: fkinLMr(l_max+1,nRstart:nRstop)
      real(kind=8),    intent(out) :: fconvLMr(l_max+1,nRstart:nRstop)
      real(kind=8),    intent(out) :: fviscLMr(l_max+1,nRstart:nRstop)
      real(kind=8),    intent(out) :: fresLMr(l_maxMag+1,nRstartMag:nRstopMag)
      real(kind=8),    intent(out) :: fpoynLMr(l_maxMag+1,nRstartMag:nRstopMag)
      real(kind=8),    intent(out) :: EperpLMr(l_max+1,nRstart:nRstop)
      real(kind=8),    intent(out) :: EparLMr(l_max+1,nRstart:nRstop)
      real(kind=8),    intent(out) :: EperpaxiLMr(l_max+1,nRstart:nRstop)
      real(kind=8),    intent(out) :: EparaxiLMr(l_max+1,nRstart:nRstop)

      !---- Output of nonlinear products for nonlinear
      !     magnetic boundary conditions (needed in s_updateB.f):
      complex(kind=8), intent(out) :: br_vt_lm_cmb(lmP_max) ! product br*vt at CMB
      complex(kind=8), intent(out) :: br_vp_lm_cmb(lmP_max) ! product br*vp at CMB
      complex(kind=8), intent(out) :: br_vt_lm_icb(lmP_max) ! product br*vt at ICB
      complex(kind=8), intent(out) :: br_vp_lm_icb(lmP_max) ! product br*vp at ICB

      !---- Output for Courant criteria:
      real(kind=8),intent(out) :: dtrkc(nRstart:nRstop),dthkc(nRstart:nRstop)


      !--- Local variables:
      integer :: nR,nr_Mag,nBc,lm
      !integer :: nTheta,nThetaB,nThetaLast
      integer :: nThetaStart!,nThetaStop
      logical :: lDeriv,lOutBc,lMagNlBc
      logical :: lGraphHeader    ! Write header into graph file
      logical :: isRadialBoundaryPoint


      PERFON('rloop')
      !LIKWID_ON('rloop')
      lGraphHeader=l_graph
      if ( lGraphHeader ) then
         call graphOut_mpi_header(time,nR,ngform,nThetaStart,sizeThetaB)
      end if

      if ( l_cour ) then
         if ( rank == 0 ) then
            dtrkc(n_r_cmb)=1.D10
            dthkc(n_r_cmb)=1.D10
         elseif (rank == n_procs-1) then
            dtrkc(n_r_icb)=1.D10
            dthkc(n_r_icb)=1.D10
         end if
      end if

      !------ Set nonlinear terms that are possibly needed at the boundaries.
      !       They may be overwritten by get_td later.
      do lm=1,lm_max
         if ( rank == 0 ) then
            dVSrLM(lm,n_r_cmb) =zero
            if ( l_mag ) then
               dVxBhLM(lm,n_r_cmb)=zero
            end if
         elseif (rank == n_procs-1) then
            dVSrLM(lm,n_r_icb) =zero
            if ( l_mag ) then
               dVxBhLM(lm,n_r_icb)=zero
            end if
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

      !------ When boundary output, Courant criterion, or non-magnetic 
      !       boundary conditions are required I have to calculate 
      !       the fields at the boundaries. This is done in one thread and 
      !       is triggered by lOutBc=.true.
      lOutBc=.false.
      if ( lTOCalc .or. lHelCalc .or. l_frame .or.         &
           & l_cour .or. l_dtB .or. lMagNlBc .or. l_graph  &
           & .or. lPerpParCalc .or. lViscBcCalc .or.       &
           & lFluxProfCalc) lOutBc=.true.

      !nRstart=n_r_cmb
      !nRstop =n_r_icb-1

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
         isRadialBoundaryPoint=(nR == n_r_cmb).or.(nR == n_r_icb)

         if ( nR == n_r_cmb ) then 
            if ( lOutBc ) then
               !nR  = n_r_cmb
               nBc = ktopv
               lDeriv= lTOCalc .or. lHelCalc .or. l_frame .or. lPerpParCalc &
            &          .or. lViscBcCalc .or. lFluxProfCalc 
            else
               cycle   ! Nothing needs to be done by thread one !
            end if
         elseif ( nR == n_r_icb ) then
            if ( lOutBc ) then
               !nR = n_r_icb
               nBc = kbotv
               lDeriv= lTOCalc .or. lHelCalc .or. l_frame  .or. lPerpParCalc &
            &          .or. lViscBcCalc .or. lFluxProfCalc
            else
               cycle
            end if
         end if

         if ( l_mag .or. l_mag_LF ) then
            nR_Mag=nR
         else
            nR_Mag=1
         end if

         call this_rIteration%set_steering_variables(l_cour,lTOCalc,lTOnext,lTOnext2,&
              & lDeriv,lRmsCalc,lHelCalc,l_frame,lMagNlBc,l_graph,lViscBcCalc,       &
              & lFluxProfCalc,lPerpParCalc)

         call this_rIteration%do_iteration(nR,nBc,time,dt,dtLast,                  &
              & dsdt(:,nR),dwdt(:,nR),dzdt(:,nR),dpdt(:,nR),dbdt(:,nR_Mag),        &
              & djdt(:,nR_Mag),dVxBhLM(:,nR_Mag),dVSrLM(:,nR),br_vt_lm_cmb,        &
              & br_vp_lm_cmb,br_vt_lm_icb,br_vp_lm_icb,lorentz_torque_ic,          &
              & lorentz_torque_ma,HelLMr(:,nR),Hel2LMr(:,nR),HelnaLMr(:,nR),       &
              & Helna2LMr(:,nR),uhLMr(:,nR),duhLMr(:,nR),gradsLMr(:,nR),           &
              & fconvLMr(:,nR),fkinLMr(:,nR),fviscLMr(:,nR),fpoynLMr(:,nR),        &
              & fresLMr(:,nR),EperpLMr(:,nR),EparLMr(:,nR),EperpaxiLMr(:,nR),      &
              & EparaxiLMr(:,nR))

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
