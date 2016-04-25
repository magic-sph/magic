#include "perflib_preproc.cpp"
module step_time_mod

#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif

   use fields
   use fieldsLast
   use parallel_mod
   use precision_mod
   use constants, only: zero, one, half
   use mem_alloc, only: bytes_allocated, memWrite
   use truncation, only: n_r_max, l_max, l_maxMag, n_r_maxMag, &
       &                 lm_max, lmP_max, lm_maxMag
   use num_param, only: n_time_steps, runTimeLimit, tEnd, dtMax, &
       &                dtMin, tScale, alpha, runTime
   use radial_data, only: nRstart, nRstop, nRstartMag, nRstopMag, &
       &                  n_r_icb, n_r_cmb
   use blocking, only: nLMBs, lmStartB, lmStopB
   use logic, only: l_mag, l_mag_LF, l_dtB, l_RMS, l_hel, l_TO,        &
       &            l_TOmovie, l_r_field, l_cmb_field, l_storeTpot,    &
       &            l_storeVpot, l_storeBpot, l_HTmovie, l_DTrMagSpec, &
       &            lVerbose, l_time_hits, l_b_nl_icb, l_b_nl_cmb,     &
       &            l_FluxProfs, l_ViscBcCalc, l_perpPar, l_HT, l_dtB, &
       &            l_dtBmovie, l_heat, l_conv, l_movie,l_true_time,   &
       &            l_runTimeLimit, l_save_out, l_dt_cmb_field,        &
       &            l_chemical_conv
   use movie_data, only: t_movieS
   use radialLoop, only: radialLoopG
   use LMLoop_data, only: llm, ulm, llmMag, ulmMag, lm_per_rank, &
       &                  lm_on_last_rank
   use LMLoop_mod, only: LMLoop
   use output_data, only: tag, n_graph_step, n_graphs, n_t_graph, t_graph, &
       &                  n_spec_step, n_specs, n_t_spec, t_spec,          &
       &                  n_movie_step, n_movie_frames, n_t_movie, t_movie,&
       &                  n_TOmovie_step, n_TOmovie_frames, n_t_TOmovie,   &
       &                  t_TOmovie, n_Bpot_step, n_Bpots, n_t_Bpot,       &
       &                  t_Bpot, n_Vpot_step, n_Vpots, n_t_Vpot, t_Vpot,  &
       &                  n_Tpot_step, n_Tpots, n_t_Tpot, t_Tpot,          &
       &                  n_rst_step, n_rsts, n_t_rst, t_rst, n_stores,    &
       &                  n_log_step, n_logs, n_t_log, t_log, n_cmb_step,  &
       &                  n_cmbs, n_t_cmb, t_cmb, n_r_field_step,          &
       &                  n_r_fields, n_t_r_field, t_r_field, n_TO_step,   &
       &                  n_TOs, n_t_TO, t_TO, n_TOZ_step, n_TOZs,         &
       &                  n_t_TOZ, t_TOZ, l_graph_time, graph_file,        &
#ifdef WITH_MPI
       &                  nLF, log_file, graph_mpi_fh, n_log_file,         &
#else
       &                  nLF, log_file, n_graph_file, n_log_file,         &
#endif
       &                  n_time_hits
   use output_mod, only: output
   use charmanip, only: capitalize, dble2str
   use useful, only: l_correct_step, safeOpen, safeClose, logWrite
   use communications, only: get_global_sum, r2lo_redist_start, lm2r_type,  &
       &                     lo2r_redist_wait, r2lm_type, lo2r_field,       &
       &                     lo2r_flow, scatter_from_rank0_to_lo, lo2r_xi,  &
       &                     r2lo_redist_wait, r2lo_flow, r2lo_s, r2lo_xi,  &
       &                     r2lo_b, lo2r_s
   use courant_mod, only: dt_courant
   use nonlinear_bcs, only: get_b_nl_bcs
   use timing ! Everything is needed

   implicit none 

   private

   !DIR$ ATTRIBUTES ALIGN:64 :: dwdt_Rloc,dzdt_Rloc,dpdt_Rloc,dsdt_Rloc,dVSrLM_Rloc
   complex(cp), allocatable, target  :: dflowdt_Rloc_container(:,:,:)
   complex(cp), allocatable, target  :: dsdt_Rloc_container(:,:,:)
   complex(cp), allocatable, target  :: dxidt_Rloc_container(:,:,:)
   complex(cp), allocatable, target  :: dbdt_Rloc_container(:,:,:)
   complex(cp), pointer :: dwdt_Rloc(:,:),dzdt_Rloc(:,:)
   complex(cp), pointer :: dpdt_Rloc(:,:), dsdt_Rloc(:,:), dVSrLM_Rloc(:,:)
   complex(cp), pointer :: dxidt_Rloc(:,:), dVXirLM_Rloc(:,:)

   !DIR$ ATTRIBUTES ALIGN:64 :: djdt_Rloc,dbdt_Rloc,dVxBhLM_Rloc
   complex(cp), pointer :: djdt_Rloc(:,:), dVxBhLM_Rloc(:,:)
   complex(cp), pointer :: dbdt_Rloc(:,:)

   ! The same arrays, but now the LM local part
   complex(cp), allocatable, target  :: dflowdt_LMloc_container(:,:,:)
   complex(cp), allocatable, target  :: dsdt_LMloc_container(:,:,:)
   complex(cp), allocatable, target  :: dxidt_LMloc_container(:,:,:)
   complex(cp), allocatable, target  :: dbdt_LMloc_container(:,:,:)
   complex(cp), pointer :: dwdt_LMloc(:,:), dzdt_LMloc(:,:)
   complex(cp), pointer :: dpdt_LMloc(:,:), dsdt_LMloc(:,:), dVSrLM_LMloc(:,:)
   complex(cp), pointer :: dxidt_LMloc(:,:), dVXirLM_LMloc(:,:)
   complex(cp), pointer :: dbdt_LMloc(:,:), djdt_LMloc(:,:), dVxBhLM_LMloc(:,:)

   complex(cp), allocatable :: dbdt_CMB_LMloc(:)

   public :: initialize_step_time,step_time

contains

   subroutine initialize_step_time

      !-- Local variables
      integer :: nR,lm
      integer(lip) :: local_bytes_used

      local_bytes_used = bytes_allocated

      allocate( dflowdt_Rloc_container(lm_max,nRstart:nRstop,1:3) )
      dwdt_Rloc(1:lm_max,nRstart:nRstop) => dflowdt_Rloc_container(:,:,1)
      dzdt_Rloc(1:lm_max,nRstart:nRstop) => dflowdt_Rloc_container(:,:,2)
      dpdt_Rloc(1:lm_max,nRstart:nRstop) => dflowdt_Rloc_container(:,:,3)

      allocate( dsdt_Rloc_container(lm_max,nRstart:nRstop,1:2) )
      dsdt_Rloc(1:lm_max,nRstart:nRstop)   => dsdt_Rloc_container(:,:,1)
      dVSrLM_Rloc(1:lm_max,nRstart:nRstop) => dsdt_Rloc_container(:,:,2)
      bytes_allocated = bytes_allocated+ &
                        5*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX

      if ( l_chemical_conv ) then
         allocate( dxidt_Rloc_container(lm_max,nRstart:nRstop,1:2) )
         dxidt_Rloc(1:lm_max,nRstart:nRstop)   => dxidt_Rloc_container(:,:,1)
         dVXirLM_Rloc(1:lm_max,nRstart:nRstop) => dxidt_Rloc_container(:,:,2)
         bytes_allocated = bytes_allocated+ &
                           2*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
      else
         allocate( dxidt_Rloc_container(1,1,1:2) )
         dxidt_Rloc(1:1,1:1)   => xi_Rloc_container(:,:,1)
         dVXirLM_Rloc(1:1,1:1) => xi_Rloc_container(:,:,2)
      end if

      ! the magnetic part
      allocate( dbdt_Rloc_container(lm_maxMag,nRstartMag:nRstopMag,1:3) )
      dbdt_Rloc(1:lm_maxMag,nRstartMag:nRstopMag)   => dbdt_Rloc_container(:,:,1)
      djdt_Rloc(1:lm_maxMag,nRstartMag:nRstopMag)   => dbdt_Rloc_container(:,:,2)
      dVxBhLM_Rloc(1:lm_maxMag,nRstartMag:nRstopMag)=> dbdt_Rloc_container(:,:,3)
      bytes_allocated = bytes_allocated+ &
                        3*lm_maxMag*(nRstopMag-nRstartMag+1)*SIZEOF_DEF_COMPLEX

      ! first touch
      do nR=nRstart,nRstop
         !$OMP PARALLEL do 
         do lm=1,lm_max
            if ( l_mag ) then
               dbdt_Rloc(lm,nR)=zero
               djdt_Rloc(lm,nR)=zero
               dVxBhLM_Rloc(lm,nR)=zero
            end if
            dwdt_Rloc(lm,nR)=zero
            dzdt_Rloc(lm,nR)=zero
            dsdt_Rloc(lm,nR)=zero
            dpdt_Rloc(lm,nR)=zero
            dVSrLM_Rloc(lm,nR)=zero
            if ( l_chemical_conv ) then
               dxidt_Rloc(lm,nR)  =zero
               dVXirLM_Rloc(lm,nR)=zero
            end if
         end do
         !$OMP END PARALLEL DO
      end do
      !call print_address("dbdt_Rloc"//C_NULL_CHAR,dbdt_Rloc)
      !call print_address("djdt_Rloc"//C_NULL_CHAR,djdt_Rloc)
      !call print_address("dsdt_Rloc"//C_NULL_CHAR,dsdt_Rloc)
      !call print_address("dVSrLM_Rloc"//C_NULL_CHAR,dVSrLM_Rloc)
      !call print_address("dVxBhLM"//C_NULL_CHAR,dVxBhLM_Rloc)

      !do lm=1,lm_maxMag,4
      !   write(str,"(A,I3,A)") "djdt_Rloc(",lm,")"//C_NULL_CHAR
      !   call print_address(str,djdt_Rloc(lm,nRStartMag))
      !end do


      ! The same arrays, but now the LM local part
      allocate(dflowdt_LMloc_container(llm:ulm,n_r_max,1:3))
      dwdt_LMloc(llm:ulm,1:n_r_max) => dflowdt_LMloc_container(:,:,1)
      dzdt_LMloc(llm:ulm,1:n_r_max) => dflowdt_LMloc_container(:,:,2)
      dpdt_LMloc(llm:ulm,1:n_r_max) => dflowdt_LMloc_container(:,:,3)

      allocate(dsdt_LMloc_container(llm:ulm,n_r_max,1:2))
      dsdt_LMloc(llm:ulm,1:n_r_max)   => dsdt_LMloc_container(:,:,1)
      dVSrLM_LMloc(llm:ulm,1:n_r_max) => dsdt_LMloc_container(:,:,2)
      bytes_allocated = bytes_allocated+ 5*(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX

      if ( l_chemical_conv ) then
         allocate(dxidt_LMloc_container(llm:ulm,n_r_max,1:2))
         dxidt_LMloc(llm:ulm,1:n_r_max)   => dxidt_LMloc_container(:,:,1)
         dVXirLM_LMloc(llm:ulm,1:n_r_max) => dxidt_LMloc_container(:,:,2)
         bytes_allocated = bytes_allocated+2*(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
      end if

      allocate(dbdt_LMloc_container(llmMag:ulmMag,n_r_maxMag,1:3))
      dbdt_LMloc(llmMag:ulmMag,1:n_r_maxMag)    => dbdt_LMloc_container(:,:,1)
      djdt_LMloc(llmMag:ulmMag,1:n_r_maxMag)    => dbdt_LMloc_container(:,:,2)
      dVxBhLM_LMloc(llmMag:ulmMag,1:n_r_maxMag) => dbdt_LMloc_container(:,:,3)
      bytes_allocated = bytes_allocated+ &
                        3*(ulmMag-llmMag+1)*n_r_maxMag*SIZEOF_DEF_COMPLEX

      ! Only when l_dt_cmb_field is requested
      ! There might be a way to allocate only when needed
      allocate ( dbdt_CMB_LMloc(llmMag:ulmMag) )
      bytes_allocated = bytes_allocated+(ulmMag-llmMag+1)*SIZEOF_DEF_COMPLEX

      local_bytes_used = bytes_allocated-local_bytes_used
      call memWrite('step_time.f90', local_bytes_used)

   end subroutine initialize_step_time
!-------------------------------------------------------------------------------
   subroutine step_time(time,dt,dtNew,n_time_step)
      !
      !  This subroutine performs the actual time-stepping.
      !

      !-- Input from initialization:
      !   time and n_time_step updated and returned to magic.f
      real(cp), intent(inout) :: time
      real(cp), intent(inout) :: dt,dtNew
      integer,  intent(inout) :: n_time_step

      !--- Local variables:

      !--- Logicals controlling output/calculation:
      logical :: l_graph          ! 
      logical :: l_spectrum
      logical :: l_cour           ! Check Courant criteria
      logical :: lCourChecking    ! Ongoing Courant criteria check
      logical :: l_store          ! Store output in restart file 
      logical :: l_new_rst_file   ! Use new rst file
      logical :: l_log            ! Log output
      logical :: l_stop_time      ! Stop time stepping
      logical :: l_frame          ! Movie frame output
      logical :: lTOframe         ! TO movie frame output
      logical :: l_cmb            ! Store set of b at CMB
      logical :: l_r              ! Store coeff at various depths
      logical :: lHelCalc         ! Calculate helicity for output
      logical :: lviscBcCalc      ! Calculate horizontal velocity and (grad T)**2
      logical :: lFluxProfCalc    ! Calculate radial flux components
      logical :: lperpParCalc     ! Calculate perpendicular and parallel Ekin
      logical :: lTOCalc          ! Calculate TO stuff
      logical :: lTONext,lTONext2 ! TO stuff for next steps
      logical :: lTOframeNext,lTOframeNext2
      logical :: lTOZhelp,lTOZwrite
      logical :: l_logNext,l_logNext2
      logical :: l_Bpot,l_Vpot,l_Tpot
      logical :: lRmsCalc,lRmsNext
      logical :: lMat             ! update matricies

      !--- Counter:
      integer :: n                ! Counter
      integer :: n_graph          ! No. of graphic file
      integer :: n_frame          ! No. of movie frames
      integer :: n_cmb_sets       ! No. of stored sets of b at CMB

      !--- Stuff needed to construct output files:
      character(len=20) :: string

      character(len=255) :: message
      character(len=76) :: SIG    

      !--- Courant criteria/diagnosis:
      real(cp) :: dtr,dth
      !-- Saves values for time step
      real(cp) :: dtrkc_Rloc(nRstart:nRstop), dthkc_Rloc(nRstart:nRstop) 

      !--- Explicit part of time stepping, calculated in s_radialLoopG.f and
      !    passed to LMLoop.f where the time step is preformed.
      !    Note that the respective arrays for the changes in inner-core
      !    magnetic field are calculated in s_updateB.f and are only
      !    needed there.

      !--- Lorentz torques:
      real(cp) :: lorentz_torque_ma,lorentz_torque_ic

      !-- Arrays for m_outMisc.F90 and m_outPar.F90
      real(cp) :: HelLMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: Hel2LMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: HelnaLMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: Helna2LMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: uhLMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: duhLMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: gradsLMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: fconvLMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: fkinLMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: fviscLMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: fpoynLMr_Rloc(l_maxMag+1,nRstartMag:nRstopMag)
      real(cp) :: fresLMr_Rloc(l_maxMag+1,nRstartMag:nRstopMag)
      real(cp) :: EperpLMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: EparLMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: EperpaxiLMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: EparaxiLMr_Rloc(l_max+1,nRstart:nRstop)

      !--- Nonlinear magnetic boundary conditions needed in s_updateB.f :
      complex(cp) :: br_vt_lm_cmb(lmP_max)    ! product br*vt at CMB
      complex(cp) :: br_vp_lm_cmb(lmP_max)    ! product br*vp at CMB
      complex(cp) :: br_vt_lm_icb(lmP_max)    ! product br*vt at ICB
      complex(cp) :: br_vp_lm_icb(lmP_max)    ! product br*vp at ICB
      complex(cp) :: b_nl_cmb(lm_max)         ! nonlinear bc for b at CMB
      complex(cp) :: aj_nl_cmb(lm_max)        ! nonlinear bc for aj at CMB
      complex(cp) :: aj_nl_icb(lm_max)        ! nonlinear bc for dr aj at ICB

      !--- Various stuff for time control:
      real(cp) :: timeLast
      real(cp) :: dtLast
      real(cp) :: w1,w2,w2New,coex
      integer :: n_time_steps_go,n_time_cour
      logical :: l_new_dt        ! causes call of matbuild !
      logical :: l_new_dtNext    ! causes call of matbuild !
      logical :: l_new_dtHit     ! causes call of matbuild !
      integer :: n_dt_changed    
      integer :: n_dt_check        
      real(cp) :: timeScaled        ! Scaled time for output.
      integer :: nPercent         ! percentage of finished time stepping
      real(cp) :: tenth_n_time_steps

      !-- Interupt procedure:
      integer :: signals(4)
      integer :: n_stop_signal     ! =1 causes run to stop
      integer :: n_graph_signal    ! =1 causes output of graphic file
      integer :: n_rst_signal      ! =1 causes output of rst file
      integer :: n_spec_signal     ! =1 causes output of a spec file
      integer :: old_stop_signal,old_graph_signal,old_rst_signal,old_spec_signal

      !--- Timing
      integer :: runTimePassed(4)
      integer :: runTimeRstart(4),runTimeRstop(4)
      integer :: runTimeTstart(4),runTimeTstop(4)
      integer :: runTimeR(4),runTimeLM(4),runTimeT(4)
      integer :: runTimeTL(4),runTimeTM(4)
      integer :: nTimeT,nTimeTL,nTimeTM,nTimeR,nTimeLM

      logical, parameter :: DEBUG_OUTPUT=.false.
      integer :: lmStart,lmStop,lmStart_00
      integer :: nR_i1,nR_i2
      ! MPI related variables
      integer :: info
      !integer :: lmStart_on_rank,lmStop_on_rank,nR
      integer, allocatable :: recvcounts(:),displs(:)
      !character(len=MPI_max_ERROR_STRING) :: error_string
      integer :: nLMB
      !integer :: nR,length_of_error

      complex(cp), pointer :: ptr_dbdt_CMB(:)
      !real(cp) :: start_time, end_time

      !integer :: signal_window
      integer(lip) :: time_in_ms


      if ( lVerbose ) write(*,'(/,'' ! STARTING STEP_TIME !'')')

#ifdef WITH_MPI
      ! allocate the buffers for MPI gathering
      allocate(recvcounts(0:n_procs-1),displs(0:n_procs-1))
      call MPI_INFO_CREATE(info,ierr)
#endif

      l_log       =.false.
      l_stop_time =.false. 
      l_new_dt    =.true.   ! Invokes calculation of t-step matricies
      l_new_dtNext=.true.  
      w2New       =-half*dtNew/dt
      n_dt_changed=0        ! No of time steps since dt changed
      n_dt_check  =4        ! No of courant checks after dt changed

      tenth_n_time_steps=real(n_time_steps,kind=cp)/10.0_cp
      nPercent=9

      !---- Set Lorentz torques to zero:
      lorentz_torque_ic=0.0_cp
      lorentz_torque_ma=0.0_cp

      !---- Counter for output files/sets:
      n_graph   =0    ! No. of graphic file
      n_frame   =0    ! No. of movie frames
      n_cmb_sets=0    ! No. of store dt_b sets at CMB

      !---- Prepare signalling via file signal
      signals=0
      n_stop_signal =0     ! Stop signal returned to calling program
      n_graph_signal=0     ! Graph signal returned to calling program
      n_spec_signal=0      ! Spec signal
      n_rst_signal=0
      if ( rank == 0 ) then
         message='signal'//'.'//tag
         open(19, file=trim(message), status='unknown')
         write(19,'(A3)') 'NOT'
         close(19)
      end if
      !call MPI_Win_create(signals,4*SIZEOF_integer,SIZEOF_integer,info,&
      !     & MPI_COMM_WORLD,signal_window,ierr)

      !-- STARTING THE TIME STEPPING LOOP:
      if ( rank == 0 ) then
         write(*,*)
         write(*,*) '! Starting time integration!'
      end if
      nTimeT =0
      nTimeTL=0
      nTimeTM=0
      nTimeR =0
      nTimeLM=0
      do n=1,4
         runTime(n)  =0
         runTimeT(n) =0
         runTimeTM(n)=0
         runTimeTL(n)=0
         runTimeR(n) =0
         runTimeLM(n)=0
      end do

      !!!!! Time loop starts !!!!!!
      n_time_cour=-2 ! Causes a Courant check after first update
      if ( n_time_steps == 1 ) then
         n_time_steps_go=1 ! Output only, for example G-file/movie etc.
      else if ( n_time_steps == 2 ) then
         n_time_steps_go=2 ! 
      else
         n_time_steps_go=n_time_steps+1  ! Last time step for output only !
      end if

#ifdef WITH_MPI
      call mpi_barrier(MPI_COMM_WORLD,ierr)
#endif

      PERFON('tloop')
      !LIKWID_ON('tloop')
      outer: do n_time_step=1,n_time_steps_go 
         n_time_cour=n_time_cour+1
         
         if ( lVerbose .or. DEBUG_OUTPUT ) then 
            write(*,*)
            write(*,*) '! Starting time step ',n_time_step
         end if

         call wallTime(runTimeTstart)

         ! =================================== BARRIER ======================
         !PERFON('barr_0')
         !call MPI_Barrier(MPI_COMM_WORLD,ierr)
         !PERFOFF
         ! ==================================================================

         ! Here now comes the block where the LM distributed fields
         ! are redistributed to Rloc distribution which is needed for the radialLoop.
         ! s,ds
         ! z,dz
         ! w,dw,ddw,p,dp
         ! b,db,ddb,aj,dj,ddj
         ! b_ic,db_ic, ddb_ic,aj_ic,dj_ic,ddj_ic

         ! Waiting for the completion before we continue to the radialLoop
         ! put the waits before signals to avoid cross communication
         PERFON('lo2r_wt')
         if (l_heat) then
            call lo2r_redist_wait(lo2r_s)
         end if
         if (l_chemical_conv) then
            call lo2r_redist_wait(lo2r_xi)
         end if
         if (l_conv) then
            call lo2r_redist_wait(lo2r_flow)
         end if

         if ( l_mag ) then
            call lo2r_redist_wait(lo2r_field)
         end if

#ifdef WITH_MPI
         ! Broadcast omega_ic and omega_ma
         call MPI_Bcast(omega_ic,1,MPI_DEF_REAL,rank_with_l1m0, &
                        MPI_COMM_WORLD,ierr)
         call MPI_Bcast(omega_ma,1,MPI_DEF_REAL,rank_with_l1m0, &
                        MPI_COMM_WORLD,ierr)
#endif
         PERFOFF

#ifdef WITH_MPI
         ! =================================== BARRIER ======================
         PERFON('barr_1')
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
         PERFOFF
         ! ==================================================================
#endif

         PERFON('signals')
         !This dealing with a signal file is quite expensive
         ! as the file can be read only on one rank and the result
         ! must be distributed to all other ranks.
         if ( rank == 0 ) then
            !----- Signalling via file signal:
            message='signal'//'.'//tag
            open(19, file=trim(message), status='old')
            read(19,*) SIG
            close(19)
            if ( len(trim(SIG)) > 0 ) then ! Non blank string ?
               call capitalize(SIG)

               old_stop_signal=n_stop_signal
               if ( index(SIG,'END')/=0 ) signals(1)=1  !n_stop_signal=1
               old_graph_signal=n_graph_signal
               if ( index(SIG,'GRA')/=0 ) then 
                  !n_graph_signal=1
                  signals(2)=1
                  open(19, file=trim(message), status='unknown')
                  write(19,'(A3)') 'NOT'
                  close(19)
               else
                  !n_graph_signal=0
                  signals(2)=0
               end if
               old_rst_signal=n_rst_signal
               if ( index(SIG,'RST')/=0 ) then
                  signals(3)=1
                  !n_rst_signal=1
                  open(19, file=trim(message), status='unknown')
                  write(19,'(A3)') 'NOT'
                  close(19)
               else
                  signals(3)=0
                  !n_rst_signal=0
               end if
               old_spec_signal=n_spec_signal
               if ( index(SIG,'SPE')/=0 ) then
                  signals(4)=1
                  !n_spec_signal=1
                  open(19, file=trim(message), status='unknown')
                  write(19,'(A3)') 'NOT'
                  close(19)
               else
                  signals(4)=0
                  !n_spec_signal=0
               end if
            end if
         end if
         ! Only broadcast the results from the signal file if
         ! something changed. For this we need one-sided communication
         ! because only process 0 knows if the communication is needed.
         !write(*,"(A)") "Win_fence 1 start"
         !PERFON('fence1')
         !call MPI_Win_fence(0,signal_window,ierr)
         !PERFOFF
         !write(*,"(A)") "Win_fence 1 end"

         ! Broadcast the results from the signal file to all processes
         ! =======> THIS IS A GLOBAL SYNCHRONIZATION POINT <==========
#if 0
         if ( rank == 0 ) then
            if ((old_stop_signal /= n_stop_signal) .or.      &
                 & (old_graph_signal /= n_graph_signal) .or. &
                 & (old_rst_signal /= n_rst_signal) .or.     &
                 & (old_spec_signal /= n_spec_signal)) then
               do iRank=1,n_procs-1
                  write(*,"(A,I4)") "MPI_putting from rank 0 to rank ",iRank
                  call MPI_Put(signals,4,MPI_integer,&
                       &       iRank,0,4,MPI_integer,signal_window,ierr)
               end do
            end if
         end if
#endif
#ifdef WITH_MPI
         call MPI_Bcast(signals,4,MPI_integer,0,MPI_COMM_WORLD,ierr)
         !call MPI_Bcast(n_stop_signal,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
         !call MPI_Bcast(n_graph_signal,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
         !call MPI_Bcast(n_spec_signal,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
         !call MPI_Bcast(n_rst_signal,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
#endif
         !write(*,"(A)") "Win_fence 2 start"
         !PERFON('fence2')
         !call MPI_Win_fence(0,signal_window,ierr)
         !PERFOFF
         !write(*,"(A)") "Win_fence 2 end"
         n_stop_signal=signals(1)
         n_graph_signal=signals(2)
         n_rst_signal=signals(3)
         n_spec_signal=signals(4)
         PERFOFF

#ifdef WITH_MPI
         PERFON('barr_2')
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
         PERFOFF
#endif

         PERFON('chk_stop')
         !--- Various reasons to stop the time integration:
         if ( l_runTimeLimit ) then
#ifdef WITH_MPI
            time_in_ms=time2ms(runTime)
            call MPI_Allreduce(MPI_IN_PLACE,time_in_ms,1,MPI_integer8, &
                 &             MPI_MAX,MPI_COMM_WORLD,ierr)
            call ms2time(time_in_ms,runTime)
#endif
            if ( lTimeLimit(runTime,runTimeLimit) ) then
               write(message,'("! Run time limit exeeded !")')
               call logWrite(message)
               l_stop_time=.true.
            end if
         end if
         if ( (n_stop_signal > 0) .or. (n_time_step == n_time_steps_go) ) then
            l_stop_time=.true.   ! last time step !
         end if

         !--- Another reasons to stop the time integration:
         if ( time >= tEND .and. tEND /= 0.0_cp ) l_stop_time=.true.
         PERFOFF
         !PERFON('logics')
         !-- Checking logic for output: 
         l_graph= l_correct_step(n_time_step-1,time,timeLast,n_time_steps,       &
              &                  n_graph_step,n_graphs,n_t_graph,t_graph,0) .or. &
              &             n_graph_signal == 1
         !l_graph=.false.
         n_graph_signal=0   ! reset interrupt signal !
         l_spectrum=                                                             &
              &         l_correct_step(n_time_step-1,time,timeLast,n_time_steps, &
              &           n_spec_step,n_specs,n_t_spec,t_spec,0) .or.            &
              &           n_spec_signal == 1
         l_frame= l_movie .and. (                                                &
              &        l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
              &        n_movie_step,n_movie_frames,n_t_movie,t_movie,0) .or.     &
              &              n_time_steps_go == 1 )
         if ( l_mag .or. l_mag_LF ) then
            l_dtB=( l_frame .and. l_dtBmovie ) .or.         &
                 &              ( l_log .and. l_DTrMagSpec ) 
         end if
         l_HT  = l_frame .and. l_HTmovie

         lTOframe=l_TOmovie .and.                                                &
              &     l_correct_step(n_time_step-1,time,timeLast,n_time_steps,     &
              &     n_TOmovie_step,n_TOmovie_frames,n_t_TOmovie,t_TOmovie,0)

         l_Bpot=l_storeBpot .and. (                                              &
              &        l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
              &                       n_Bpot_step,n_Bpots,n_t_Bpot,t_Bpot,0).or. &
              &            n_time_steps == 1 )
         l_Vpot=l_storeVpot .and. (                                              &
              &        l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
              &                       n_Vpot_step,n_Vpots,n_t_Vpot,t_Vpot,0).or. &
              &            n_time_steps == 1 )
         l_Tpot=l_storeTpot .and. (                                              &
              &        l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
              &                       n_Tpot_step,n_Tpots,n_t_Tpot,t_Tpot,0).or. &
              &            n_time_steps == 1 )

         l_cour=.true.

         l_new_rst_file=                                                         & 
              &        l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
              &                       n_rst_step,n_rsts,n_t_rst,t_rst,0) .or.    &
              &        n_rst_signal == 1
         n_rst_signal=0
         l_store= l_new_rst_file .or.                                            &
              &        l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
              &                       0,n_stores,0,t_rst,0)

         l_log= l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
              &                       n_log_step,n_logs,n_t_log,t_log,0)
         l_cmb= l_cmb_field .and.                                                &
              &        l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
              &                       n_cmb_step,n_cmbs,n_t_cmb,t_cmb,0)
         l_r= l_r_field .and.                                                    &
              &        l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
              &                       n_r_field_step,n_r_fields,n_t_r_field,     &
              &                       t_r_field,0)
         l_logNext=.false.
         if ( n_time_step+1 <= n_time_steps+1 )                       &
              &        l_logNext=                                     &
              &        l_correct_step(n_time_step,time+dt,timeLast,   &
              &              n_time_steps,n_log_step,n_logs,n_t_log,t_log,0)
         l_logNext2=.false.
         if ( n_time_step+2 <= n_time_steps+1 )                         &
              &        l_logNext2=                                      &
              &        l_correct_step(n_time_step+1,time+2*dt,timeLast, &
              &         n_time_steps,n_log_step,n_logs,n_t_log,t_log,0)
         lTOCalc= n_time_step > 2 .and. l_TO .and.                   &
              &          l_correct_step(n_time_step-1,time,timeLast, &
              &          n_time_steps,n_TO_step,n_TOs,n_t_TO,t_TO,0)
         lTOnext     =.false.
         lTOframeNext=.false.
         if ( n_time_step+1 <= n_time_steps+1 ) then
            lTONext= l_TO .and.                                            &
                 &           l_correct_step(n_time_step,time+dt,timeLast,  &
                 &            n_time_steps,n_TO_step,n_TOs,n_t_TO,t_TO,0)
            lTOframeNext= l_TOmovie .and.                                   &
                 &           l_correct_step(n_time_step,time+dt,timeLast,   &
                 &          n_time_steps,n_TOmovie_step,n_TOmovie_frames,   &
                 &                               n_t_TOmovie,t_TOmovie,0)
         end if
         lTONext      =lTOnext.or.lTOframeNext 
         lTONext2     =.false.
         lTOframeNext2=.false.
         if ( n_time_step+2 <= n_time_steps+1 ) then
            lTONext2= l_TO .and.                                              &
                 &           l_correct_step(n_time_step+1,time+2*dt,timeLast, &
                 &                                    n_time_steps,n_TO_step, &
                 &                                       n_TOs,n_t_TO,t_TO,0)
            lTOframeNext2= l_TOmovie .and.                                    &
                 &           l_correct_step(n_time_step+1,time+2*dt,timeLast, &
                 &                               n_time_steps,n_TOmovie_step, &
                 &                  n_TOmovie_frames,n_t_TOmovie,t_TOmovie,0)
         end if
         lTONext2=lTOnext2.or.lTOframeNext2 
         lTOZhelp= n_time_step > 2 .and. l_TO .and.                         &
              &                l_correct_step(n_time_step-1,time,timeLast,  &
              &            n_time_steps,n_TOZ_step,n_TOZs,n_t_TOZ,t_TOZ,0)
         if ( lTOZhelp ) lTOZwrite=.true.

         lRmsCalc=l_RMS.and.l_log.and.(n_time_step > 1)
         if ( l_mag .or. l_mag_LF ) l_dtB = l_dtB .or. lRmsCalc
         lRmsNext=l_RMS.and.l_logNext ! Used for storing in update routines !

         if ( n_time_step == 1 ) l_log=.true.

         if ( l_stop_time ) then                  ! Programm stopped by kill -30
            l_new_rst_file=.true.                 ! Write rst-file and some
            if ( n_stores > 0 ) l_store=.true.    ! diagnostics before dying ! 
            l_log=.true.
            lRmsNext=.false.
         end if

         lHelCalc=l_hel.and.l_log

         lperpParCalc=l_perpPar.and.l_log

         lFluxProfCalc =l_FluxProfs.and.l_log

         lViscBcCalc =l_ViscBcCalc.and.l_log

         if ( l_graph ) then  ! write graphic output !
            PERFON('graph')
            n_graph=n_graph+1     ! increase counter for graphic file
            if ( l_graph_time ) then 
               call dble2str(time,string)
               graph_file='G_t='//trim(string)//'.'//tag
            else
               write(string, *) n_graph
               graph_file='G_'//trim(adjustl(string))//'.'//tag
            end if
            if ( rank == 0 ) then
               write(*,'(1p,/,A,/,A,ES20.10,/,A,i15,/,A,A)')&
                    &" ! Storing graphic file:",            &
                    &"             at time=",timeScaled,    &
                    &"            step no.=",n_time_step,   &
                    &"           into file=",graph_file
               call safeOpen(nLF,log_file)
               write(nLF,'(1p,/,A,/,A,ES20.10,/,A,i15,/,A,A)') &
                    &" ! Storing graphic file:",               &
                    &"             at time=",timeScaled,       &
                    &"            step no.=",n_time_step,      &
                    &"           into file=",graph_file
               call safeClose(nLF)
            end if
#ifdef WITH_MPI
            call MPI_File_open(MPI_COMM_WORLD,graph_file,             &
                 &             IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),  &
                 &             MPI_INFO_NULL,graph_mpi_fh,ierr)
#else
            open(n_graph_file,file=graph_file,status='new',form='unformatted')
#endif
            !call MPI_ERROR_STRING(ierr,error_string,length_of_error,ierr)
            !PRINT*,"MPI_FILE_OPEN returned: ",trim(error_string)
            PERFOFF
         end if

         if ( DEBUG_OUTPUT ) then
            do nLMB=1+rank*nLMBs_per_rank,min((rank+1)*nLMBs_per_rank,nLMBs)
               lmStart=lmStartB(nLMB)
               lmStop=lmStopB(nLMB)
               lmStart_00  =max(2,lmStart)
               
               !do nR=1,n_r_max
               !   write(*,"(A,I2,A,2ES20.12)") "dw_LMloc for nR=",nR," is ", &
               !        SUM( dw_LMloc(lmStart:lmStop,nR) )
               !end do
               write(*,"(A,I3,6ES20.12)") "start w: ",nLMB,        &
                    & GET_GLOBAL_SUM( w_LMloc(lmStart:lmStop,:) ), &
                    & GET_GLOBAL_SUM( dw_LMloc(lmStart:lmStop,:) ),&
                    & GET_GLOBAL_SUM( ddw_LMloc(lmStart:lmStop,:) )
               write(*,"(A,I3,4ES20.12)") "start z: ",nLMB,       &
                    & GET_GLOBAL_SUM( z_LMloc(lmStart:lmStop,:) ),&
                    & GET_GLOBAL_SUM( dz_LMloc(lmStart:lmStop,:) )
               write(*,"(A,I3,4ES20.12)") "start s: ",nLMB,       &
                    & GET_GLOBAL_SUM( s_LMloc(lmStart:lmStop,:) ),&
                    & GET_GLOBAL_SUM( ds_LMloc(lmStart:lmStop,:) )
               write(*,"(A,I3,4ES20.12)") "start p: ",nLMB,       &
                    & GET_GLOBAL_SUM( p_LMloc(lmStart:lmStop,:) ),&
                    & GET_GLOBAL_SUM( dp_LMloc(lmStart_00:lmStop,:) )
               if ( l_mag ) then
                  write(*,"(A,I3,8ES20.12)") "start b: ",nLMB,        &
                       & GET_GLOBAL_SUM( b_LMloc(lmStart:lmStop,:) ), &
                       & GET_GLOBAL_SUM( db_LMloc(lmStart:lmStop,:) ),&
                       & GET_GLOBAL_SUM( ddb_LMloc(lmStart:lmStop,:) )
                  write(*,"(A,I3,8ES20.12)") "start aj: ",nLMB,        &
                       & GET_GLOBAL_SUM( aj_LMloc(lmStart:lmStop,:) ), &
                       & GET_GLOBAL_SUM( dj_LMloc(lmStart:lmStop,:) ), &
                       & GET_GLOBAL_SUM( ddj_LMloc(lmStart:lmStop,:) )
                  write(*,"(A,I3,8ES20.12)") "start b_ic: ",nLMB,        &
                       & GET_GLOBAL_SUM( b_ic_LMloc(lmStart:lmStop,:) ), &
                       & GET_GLOBAL_SUM( db_ic_LMloc(lmStart:lmStop,:) ),&
                       & GET_GLOBAL_SUM( ddb_ic_LMloc(lmStart:lmStop,:) )
                  write(*,"(A,I3,8ES20.12)") "start aj_ic: ",nLMB,       &
                       & GET_GLOBAL_SUM( aj_ic_LMloc(lmStart:lmStop,:) ),&
                       & GET_GLOBAL_SUM( dj_ic_LMloc(lmStart:lmStop,:) ),&
                       & GET_GLOBAL_SUM( ddj_ic_LMloc(lmStart:lmStop,:) )
               end if
            end do
         end if
         
         !--- Now the real work starts with the radial loop that calculates
         !    the nonlinear terms:
         if ( lVerbose ) then 
            write(*,*)
            write(*,*) '! Starting radial loop!'
         end if

         !PERFOFF
         ! =============================== BARRIER ===========================
         !PERFON('barr_2')
         !call MPI_Barrier(MPI_COMM_WORLD,ierr)
         !PERFOFF
         ! ===================================================================

         call wallTime(runTimeRstart)
         call radialLoopG(l_graph,l_cour,l_frame,time,dt,dtLast,               &
              &           lTOCalc,lTONext,lTONext2,lHelCalc,lRmsCalc,          &
              &           lViscBcCalc,lFluxProfCalc,lperpParCalc,              &
              &           dsdt_Rloc,dwdt_Rloc,dzdt_Rloc,dpdt_Rloc,dxidt_Rloc,  &
              &           dbdt_Rloc,djdt_Rloc,dVxBhLM_Rloc,dVSrLM_Rloc,        &
              &           dVXirLM_Rloc,lorentz_torque_ic,                      &
              &           lorentz_torque_ma,br_vt_lm_cmb,br_vp_lm_cmb,         &
              &           br_vt_lm_icb,br_vp_lm_icb,HelLMr_Rloc,Hel2LMr_Rloc,  &
              &           HelnaLMr_Rloc,Helna2LMr_Rloc,uhLMr_Rloc,duhLMr_Rloc, &
              &           gradsLMr_Rloc,fconvLMr_Rloc,fkinLMr_Rloc,            &
              &           fviscLMr_Rloc,fpoynLMr_Rloc,fresLMr_Rloc,            &
              &           EperpLMr_Rloc,EparLMr_Rloc,EperpaxiLMr_Rloc,         &
              &           EparaxiLMr_Rloc,dtrkc_Rloc,dthkc_Rloc)

         if ( lVerbose ) write(*,*) '! r-loop finished!'
         if ( .not.l_log ) then
            call wallTime(runTimeRstop)
            if ( .not.lNegTime(runTimeRstart,runTimeRstop) ) then
               nTimeR=nTimeR+1
               call subTime(runTimeRstart,runTimeRstop,runTimePassed)
               call addTime(runTimeR,runTimePassed)
            end if
         end if

         !---------------------------------------
         !--- Gather all r-distributed arrays ---
         !---------------------------------------
         ! we have here the following arrays from radialLoopG:
         !  dsdt,dwdt,dzdt,dpdt,dbdt,djdt,dVxBhLM,dVSrLM
         !  TstrRLM,TadvRLM,TomeRLM, 
         !  HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,dtrkc,dthkc
         !
         ! gather in-place, we need allgatherV because auf the unequal
         ! number of points on the processes (last block is one larger)

         if ( DEBUG_OUTPUT ) then
            nR_i1=max(1,nRstart)
            nR_i2=min(n_r_max,nRstop)
            write(*,"(A,10ES20.12)") "middl: dwdt,dsdt,dzdt,dpdt = ",&
                 & GET_GLOBAL_SUM( dwdt_Rloc(:,nR_i1:nR_i2) ),       &
                 & GET_GLOBAL_SUM( dsdt_Rloc(:,nR_i1:nR_i2) ),       &
                 & GET_GLOBAL_SUM( dzdt_Rloc(:,nR_i1:nR_i2) ),       &
                 & GET_GLOBAL_SUM( dpdt_Rloc(:,nR_i1:nR_i2) ),       &
                 & GET_GLOBAL_SUM( dVSrLM_Rloc )
            if ( l_mag ) then
               write(*,"(A,6ES20.12)") "middl: dbdt,djdt,dVxBhLM = ",&
                    & GET_GLOBAL_SUM( dbdt_Rloc(:,nR_i1:nR_i2) ),    &
                    & GET_GLOBAL_SUM( djdt_Rloc(:,nR_i1:nR_i2) ),    &
                    & GET_GLOBAL_SUM( dVxBhLM_Rloc )
            end if
         end if

         ! ===================================== BARRIER =======================
         !PERFON('barr_rad')
         !call MPI_Barrier(MPI_COMM_WORLD,ierr)
         !PERFOFF
         ! =====================================================================
         if ( lVerbose ) write(*,*) "! start r2lo redistribution"

         PERFON('r2lo_dst')
         call r2lo_redist_start(r2lo_flow,dflowdt_Rloc_container,dflowdt_LMloc_container)
         call r2lo_redist_wait(r2lo_flow)

         if ( l_heat ) then
            call r2lo_redist_start(r2lo_s,dsdt_Rloc_container,dsdt_LMloc_container)
            call r2lo_redist_wait(r2lo_s)
         end if

         if ( l_chemical_conv ) then
            call r2lo_redist_start(r2lo_xi,dxidt_Rloc_container,dxidt_LMloc_container)
            call r2lo_redist_wait(r2lo_xi)
         end if

         if ( l_mag ) then
            call r2lo_redist_start(r2lo_b,dbdt_Rloc_container,dbdt_LMloc_container)
            call r2lo_redist_wait(r2lo_b)
         end if

#ifdef WITH_MPI
         ! ------------------
         ! also exchange the lorentz_torques which are only set at the boundary points
         ! but are needed on all processes.
         call MPI_Bcast(lorentz_torque_ic,1,MPI_DEF_REAL, &
              &         n_procs-1,MPI_COMM_WORLD,ierr)
         call MPI_Bcast(lorentz_torque_ma,1,MPI_DEF_REAL, &
              &         0,MPI_COMM_WORLD,ierr)
#endif
         PERFOFF
         if ( lVerbose ) write(*,*) "! r2lo redistribution finished"

         if ( DEBUG_OUTPUT ) then
            write(*,"(A,8ES20.12)")                                           &
                 & "lo_arr middl: dzdt_LMloc,z_LMloc,dz_LMloc,dzdtLast_lo = ",&
                 & GET_GLOBAL_SUM( dzdt_LMloc(:,2:n_r_max-1) ),               &
                 & GET_GLOBAL_SUM( z_LMloc ),GET_GLOBAL_SUM( dz_LMloc ),      &
                 & GET_GLOBAL_SUM( dzdtLast_lo )
            write(*,"(A,8ES20.12)") "lo_arr middl: dsdt,s,ds,dsdtLast = ", &
                 & GET_GLOBAL_SUM( dsdt_LMloc(:,2:n_r_max-1) ),            &
                 & GET_GLOBAL_SUM( s_LMloc ),GET_GLOBAL_SUM( ds_LMloc ),   &
                 & GET_GLOBAL_SUM( dsdtLast_LMloc )
            write(*,"(A,10ES20.12)")                            &
                 & "lo_arr middl: dwdt,w,dw,ddw,dwdtLast = ",   &
                 & GET_GLOBAL_SUM( dwdt_LMloc(:,2:n_r_max-1) ), &
                 & GET_GLOBAL_SUM( w_LMloc ),                   &
                 & GET_GLOBAL_SUM( dw_LMloc ),                  &
                 & GET_GLOBAL_SUM( ddw_LMloc ),                 &
                 & GET_GLOBAL_SUM( dwdtLast_LMloc )
            if ( l_mag ) then
               write(*,"(A,10ES20.12)")                            &
                    &"lo_arr middl: dbdt,b,db,ddb,dbdtLast = ",    &
                    & GET_GLOBAL_SUM( dbdt_LMloc(:,2:n_r_max-1) ), &              
                    & GET_GLOBAL_SUM( b_LMloc ),                   &
                    & GET_GLOBAL_SUM( db_LMloc ),                  &
                    & GET_GLOBAL_SUM( ddb_LMloc ),                 &
                    & GET_GLOBAL_SUM( dbdtLast_LMloc )
               write(*,"(A,12ES20.12)")                                     &
                    & "lo_arr middl: djdt,aj,dj,ddj,djdtLast,dVxBhLM = ",   &
                    & GET_GLOBAL_SUM( djdt_LMloc(:,2:n_r_max-1) ),          &
                    & GET_GLOBAL_SUM( aj_LMloc ),GET_GLOBAL_SUM( dj_LMloc ),&
                    & GET_GLOBAL_SUM( ddj_LMloc ),                          &
                    & GET_GLOBAL_SUM( djdtLast_LMloc ),                     &
                    & GET_GLOBAL_SUM( dVxBhLM_LMloc )
            end if
            write(*,"(A,2ES20.12)") "middl: dtrkc,dthkc = ", &
                  SUM(dtrkc_Rloc),SUM(dthkc_Rloc)
         end if

         !PERFOFF

         !--- Output before update of fields in LMLoop:
         ! =================================== BARRIER ======================
         !PERFON('barr_4')
         !call MPI_Barrier(MPI_COMM_WORLD,ierr)
         !PERFOFF
         ! ==================================================================
         if ( lVerbose ) write(*,*) "! start output"
         PERFON('output')
         if ( nRstart <= n_r_cmb .and. l_cmb .and. l_dt_cmb_field ) then
            ptr_dbdt_CMB => dbdt_Rloc(:,n_r_cmb)
            call scatter_from_rank0_to_lo(ptr_dbdt_CMB, dbdt_CMB_LMloc)
         end if
         if ( lVerbose ) write(*,*) "! start real output"
         call output(time,dt,dtNew,n_time_step,l_stop_time,                &
              &      l_Bpot,l_Vpot,l_Tpot,l_log,l_graph,lRmsCalc,          &
              &      l_store,l_new_rst_file,                               &
              &      l_spectrum,lTOCalc,lTOframe,lTOZwrite,                &
              &      l_frame,n_frame,l_cmb,n_cmb_sets,l_r,                 &
              &      lorentz_torque_ic,lorentz_torque_ma,dbdt_CMB_LMloc,   &
              &      HelLMr_Rloc,Hel2LMr_Rloc,HelnaLMr_Rloc,Helna2LMr_Rloc,&
              &      uhLMr_Rloc,duhLMr_Rloc,gradsLMr_Rloc,fconvLMr_Rloc,   &
              &      fkinLMr_Rloc,fviscLMr_Rloc,fpoynLMr_Rloc,             & 
              &      fresLMr_Rloc,EperpLMr_Rloc,EparLMr_Rloc,              & 
              &      EperpaxiLMr_Rloc,EparaxiLMr_Rloc)
         PERFOFF
         if ( lVerbose ) write(*,*) "! output finished"

         if ( l_graph ) then
#ifdef WITH_MPI
            PERFON('graph')
            call MPI_File_close(graph_mpi_fh,ierr)
            !call MPI_ERROR_STRING(ierr,error_string,length_of_error,ierr)
            !PRINT*,"MPI_FILE_CLOSE returned: ",trim(error_string)
            PERFOFF
#else
            close(n_graph_file)
#endif
         end if
         ! =================================== BARRIER ======================
         !PERFON('barr_5')
         !call MPI_Barrier(MPI_COMM_WORLD,ierr)
         !PERFOFF
         ! ==================================================================

         !----- Finish time stepping, the last step is only for output!
         if ( l_stop_time ) exit outer  ! END OF TIME INTEGRATION

         !------ Nonlinear magnetic boundary conditions:
         !       For stressfree conducting boundaries
         PERFON('nl_m_bnd')
         if ( l_b_nl_cmb ) then
            b_nl_cmb(1) =(1.0_cp,1.0_cp)
            aj_nl_cmb(1)=(1.0_cp,1.0_cp)
            call get_b_nl_bcs('CMB', br_vt_lm_cmb,br_vp_lm_cmb,              &
                 &            2,lm_max,b_nl_cmb(2:lm_max),aj_nl_cmb(2:lm_max))
         end if
         if ( l_b_nl_icb ) then
            aj_nl_icb(1)=(1.0_cp,1.0_cp)
            call get_b_nl_bcs('ICB', br_vt_lm_icb,br_vp_lm_icb,              &
                 &            2,lm_max,b_nl_cmb(2:lm_max),aj_nl_icb(2:lm_max))
         end if
         PERFOFF

         PERFON('t_check')
         !---- Time-step check and change if needed (l_new_dtNext=.true.)
         !     I anticipate the dt change here that is only used in 
         !     the next time step cause its needed for coex=(alpha-1)/w2 !
         dtLast=dt
         dt=dtNew        ! Update to new time step
         w2=w2New        ! Weight for time-derivatives of last time step
         w1=one-w2      ! Weight for currect time step
         l_new_dt    =l_new_dtNext
         l_new_dtNext=.false.

         !------ Checking Courant criteria, l_new_dt and dtNew are output
         if ( l_cour ) then
            !PRINT*,"dtrkc_Rloc = ",dtrkc_Rloc
            call dt_courant(dtr,dth,l_new_dtNext,dt,dtNew,dtMax, &
                 &          dtrkc_Rloc,dthkc_Rloc)
         end if

         !------ Checking whether we have to hit a specific time for output,
         !       dtNew is changed accordingly and l_new_dt is set .true.
         if ( l_true_time .and. l_time_hits ) then
            call check_time_hits(l_new_dtHit,time,dt,dtNew)
            l_new_dtNext=l_new_dtNext .or. l_new_dtHit
         end if

         !----- Stop if time step has become too small:
         if ( dtNew < dtMin ) then
            if ( rank == 0 ) then
               write(*,'(1p,/,A,ES14.4,/,A)')            &
                    &" ! TIME STEP TOO SMALL, dt=",dtNew,&
                    &" ! I THUS stop THE RUN !"
               call safeOpen(nLF,log_file)
               write(nLF,'(1p,/,A,ES14.4,/,A)')          &
                    &" ! TIME STEP TOO SMALL, dt=",dtNew,&
                    &" ! I THUS stop THE RUN !"
               call safeClose(nLF)
            end if
            stop
         end if
         if ( l_new_dtNext ) then
            !------ Writing info and getting new weights:
            w2New=-half*dtNew/dt ! Weight I will be using for next update !
            n_dt_changed=0
            lCourChecking=.true.
            if ( rank == 0 ) then
               write(*,'(1p,/,A,ES18.10,/,A,i9,/,A,ES15.8,/,A,ES15.8)')  &
                    &" ! Changing time step at time=",(time+dt)*tScale,  &
                    &"                 time step no=",n_time_step+1,     &
                    &"                      last dt=",dt*tScale,         &
                    &"                       new dt=",dtNew*tScale
               call safeOpen(nLF,log_file)
               write(nLF,'(1p,/,A,ES18.10,/,A,i9,/,A,ES15.8,/,A,ES15.8)') &
                    &" ! Changing time step at time=",(time+dt)*tScale,   &
                    &"                 time step no=",n_time_step+1,      &
                    &"                      last dt=",dt*tScale,          &
                    &"                       new dt=",dtNew*tScale
               call safeClose(nLF)
            end if
         else
            w2New=-half ! Normal weight if timestep is not changed !
            n_dt_changed=n_dt_changed+1
            if ( n_dt_changed <= n_dt_check  ) then
               lCourChecking=.true.
            else
               lCourChecking=.false.
            end if
         end if
         PERFOFF


         !----- UPDATING THE FIELDS:
         !      This is the second parallel part. Here we parallize over lm.

         !----- Advancing time:
         coex  =(alpha-one)/w2New
         timeLast        =time               ! Time of last time step
         time            =time+dt            ! Update time
         timeScaled      =time*tScale

         lMat=.false.
         if ( l_new_dt ) then
            !----- Calculate matricies for new time step if dt /= dtLast
            lMat=.true.
            if ( rank == 0 ) then
               write(*,'(1p,/,'' ! BUILDING MATRICIES AT STEP/TIME:'',   &
                    &              i8,ES16.6)') n_time_step,timeScaled
            end if
         end if

         ! =================================== BARRIER ======================
         !PERFON('barr_6')
         !call MPI_Barrier(MPI_COMM_WORLD,ierr)
         !PERFOFF
         ! ==================================================================

         call wallTime(runTimeRstart)
         if ( lVerbose ) write(*,*) '! starting lm-loop!'
         if ( lVerbose ) write(*,*) '! No of lm-blocks:',nLMBs

         call LMLoop(w1,coex,time,dt,lMat,lRmsNext,dVxBhLM_LMloc,   &
              &      dVSrLM_LMloc,dVXirLM_LMloc,dsdt_LMloc,         &
              &      dwdt_LMloc,dzdt_LMloc,dpdt_LMloc,dxidt_LMloc,  &
              &      dbdt_LMloc,djdt_LMloc,                         &
              &      lorentz_torque_ma,lorentz_torque_ic,           &
              &      b_nl_cmb,aj_nl_cmb,aj_nl_icb)

         if ( lVerbose ) write(*,*) '! lm-loop finished!'
         call wallTime(runTimeRstop)

         if ( .not.lNegTime(runTimeRstart,runTimeRstop) ) then
            nTimeLM=nTimeLM+1
            call subTime(runTimeRstart,runTimeRstop,runTimePassed)
            call addTime(runTimeLM,runTimePassed)
         end if

         if (DEBUG_OUTPUT) then
            write(*,"(A,6ES20.12)") "lo_arr end: z_LMloc,dz_LMloc,dzdtLast_lo = ",&
                 & GET_GLOBAL_SUM( z_LMloc ),               &
                 & GET_GLOBAL_SUM( dz_LMloc ),              &
                 & GET_GLOBAL_SUM( dzdtLast_lo )
            write(*,"(A,6ES20.12)") "lo_arr end: s,ds,dsdtLast = ",&
                 & GET_GLOBAL_SUM( s_LMloc ),                      &
                 & GET_GLOBAL_SUM( ds_LMloc ),                     &
                 & GET_GLOBAL_SUM( dsdtLast_LMloc )
            write(*,"(A,6ES20.12)") "lo_arr end: w,dw,dwdtLast = ",&
                 & GET_GLOBAL_SUM( w_LMloc ),                      &
                 & GET_GLOBAL_SUM( dw_LMloc ),                     &
                 & GET_GLOBAL_SUM( dwdtLast_LMloc )
            write(*,"(A,4ES20.12)") "w(bnd_r) = ",            &
                 & GET_GLOBAL_SUM( w_LMloc(:,n_r_icb) ),      &
                 & GET_GLOBAL_SUM( w_LMloc(:,n_r_cmb) )
            write(*,"(A,6ES20.12)") "lo_arr end: p,dpdtLast = ",&
                 & GET_GLOBAL_SUM( p_LMloc ),                   &
                 & GET_GLOBAL_SUM( dpdtLast_LMloc )
            write(*,"(A,6ES20.12)") "lo_arr end: b,db,dbdtLast = ",&
                 & GET_GLOBAL_SUM( b_LMloc ),                      &
                 & GET_GLOBAL_SUM( db_LMloc ),                     &
                 & GET_GLOBAL_SUM( dbdtLast_LMloc )
            write(*,"(A,6ES20.12)") "lo_arr end: aj,dj,djdtLast = ",&
                 & GET_GLOBAL_SUM( aj_LMloc ),                      &
                 & GET_GLOBAL_SUM( dj_LMloc ),                      &
                 & GET_GLOBAL_SUM( djdtLast_LMloc )
         end if

         !----- Timing and info of advancement:
         ! =================================== BARRIER ======================
         !start_time=MPI_Wtime()
         !PERFON('barr_lm')
         !call MPI_Barrier(MPI_COMM_WORLD,ierr)
         !PERFOFF
         !end_time=MPI_Wtime()
         !write(*,"(A,I4,A,F10.6,A)") " lm barrier on rank ",rank,&
         !     & " takes ",end_time-start_time," s."
         ! ==================================================================
         call wallTime(runTimeTstop)
         if ( .not.lNegTime(runTimeTstart,runTimeTstop) ) then
            call subTime(runTimeTstart,runTimeTstop,runTimePassed)
            if ( lMat .and. .not.l_log ) then
               nTimeTM=nTimeTM+1
               call addTime(runTimeTM,runTimePassed)
            else if ( .not.lMat .and. l_log ) then
               nTimeTL=nTimeTL+1
               call addTime(runTimeTL,runTimePassed)
            else if ( .not.lMat .and. .not.l_log ) then
               nTimeT=nTimeT+1
               call addTime(runTimeT,runTimePassed)
            end if
            call addTime(runTime,runTimePassed)
         end if
         if ( real(n_time_step,cp)+tenth_n_time_steps*real(nPercent,cp) >=  &
            & real(n_time_steps,cp)  .or. n_time_steps < 31 ) then
            write(message,'(" ! Time step finished:",i6)') n_time_step
            call logWrite(message)
            if ( real(n_time_step,cp)+tenth_n_time_steps*real(nPercent,cp) >= &
               & real(n_time_steps,cp) .and. n_time_steps >= 10 ) then
               write(message,'(" ! This is           :",i3,"%")') (10-nPercent)*10
               call logWrite(message)
               nPercent=nPercent-1
            end if
            do n=1,4
               runTimePassed(n)=runTimeT(n)
            end do
            if ( nTimeT > 0 ) then
               call meanTime(runTimePassed,nTimeT)
               if ( rank == 0 ) then
                  call writeTime(6,'! Mean wall time for time step:', runTimePassed)
                  call safeOpen(nLF,log_file)
                  call writeTime(nLF,'! Mean wall time for time step:', runTimePassed)
                  call safeClose(nLF)
               end if
            end if
         end if

      end do outer ! end of time stepping !

      !LIKWID_OFF('tloop')
      PERFOFF

      if ( l_movie ) then
         if ( rank == 0 ) then
            if (n_frame > 0) then
               write(*,'(1p,/,/,A,i10,3(/,A,ES16.6))')                    &
                    &" !  No of stored movie frames: ",n_frame,           &
                    &" !     starting at time: ",t_movieS(1)*tScale,      &
                    &" !       ending at time: ",t_movieS(n_frame)*tScale,&
                    &" !      with step width: ",(t_movieS(2)-t_movieS(1))*tScale
               call safeOpen(nLF,log_file)
               write(nLF,'(1p,/,/,A,i10,3(/,A,ES16.6))')                  &
                    &" !  No of stored movie frames: ",n_frame,           &
                    &" !     starting at time: ",t_movieS(1)*tScale,      &
                    &" !       ending at time: ",t_movieS(n_frame)*tScale,&
                    &" !      with step width: ",(t_movieS(2)-t_movieS(1))*tScale
               call safeClose(nLF)
            else
               write(*,'(1p,/,/,A,i10,3(/,A,ES16.6))')          &
                    &" !  No of stored movie frames: ",n_frame, &
                    &" !     starting at time: ",0.0_cp,        &
                    &" !       ending at time: ",0.0_cp,        &
                    &" !      with step width: ",0.0_cp
               call safeOpen(nLF,log_file)
               write(nLF,'(1p,/,/,A,i10,3(/,A,ES16.6))')        &
                    &" !  No of stored movie frames: ",n_frame, &
                    &" !     starting at time: ",0.0_cp,        &
                    &" !       ending at time: ",0.0_cp,        &
                    &" !      with step width: ",0.0_cp
               call safeClose(nLF)
            end if
         end if
      end if

      if ( l_cmb_field ) then
         write(message,'(A,i9)') " !  No of stored sets of b at CMB: ",n_cmb_sets
         call logWrite(message)
      end if

      call meanTime(runTimeR, nTimeR)
      call meanTime(runTimeLM,nTimeLM)
      call meanTime(runTimeTM,nTimeTM)
      call meanTime(runTimeTL,nTimeTL)
      call meanTime(runTimeT,nTimeT)
      if ( rank == 0 ) then
         call writeTime(6,'! Mean wall time for r Loop                 :',runTimeR)
         call writeTime(6,'! Mean wall time for LM Loop                :',runTimeLM)
         call writeTime(6,'! Mean wall time for t-step with matrix calc:',runTimeTM)
         call writeTime(6,'! Mean wall time for t-step with log output :',runTimeTL)
         call writeTime(6,'! Mean wall time for pure t-step            :',runTimeT)
         call safeOpen(nLF,log_file)
         call writeTime(nLF,'! Mean wall time for r Loop                 :',runTimeR)
         call writeTime(nLF,'! Mean wall time for LM Loop                :',runTimeLM)
         call writeTime(nLF,'! Mean wall time for t-step with matrix calc:',runTimeTM)
         call writeTime(nLF,'! Mean wall time for t-step with log output :',runTimeTL)
         call writeTime(nLF,'! Mean wall time for pure t-step            :',runTimeT)
         call safeClose(nLF)
      end if
      !-- Write output for variable conductivity test:
      !       if ( imagcon == -10 ) then
      !          message='testVarCond.'//tag
      !          open(99,file=message,status='unknown')
      !           do nR=n_r_max,1,-1             ! Diffusive toroidal field
      !             write(99,*) r(nR),real(aj(4,nR)),jVarCon(nR)
      !           end do
      !          close(99)
      !       end if


      !-- WORK IS DONE !
#ifdef WITH_MPI
      call MPI_INFO_FREE(info,ierr)
      deallocate(recvcounts,displs)
#endif

   end subroutine step_time
!------------------------------------------------------------------------------
   subroutine check_time_hits(l_new_dt,time,dt,dt_new)
      !
      !  Checks whether a certain dt is required to hit a                 
      !  specific output-time.                                            
      !

      !-- Output: ev. modified dt
      logical,  intent(out) :: l_new_dt ! signfies change of dt !
      real(cp), intent(inout) :: time,dt,dt_new
       
      !-- Local variables:
      integer :: n_dt_hit
      integer, parameter :: n_dt_hit_max=10
      real(cp) ::  dt_hit(n_dt_hit_max) ! dt for different hit times
      integer :: n                          ! counter
      real(cp) ::  time_new             ! Next time step

      time_new=time+dt
      l_new_dt=.false.

      n_dt_hit=7

      do n=1,n_dt_hit
         dt_hit(n)=0.0_cp
      end do

      do n=1,n_time_hits
         if ( t_rst(n) > time .and. t_rst(n) < time_new ) &
              dt_hit(1)=t_rst(n)-time
         if ( t_graph(n) > time .and. t_graph(n) < time_new ) &
              dt_hit(2)=t_graph(n)-time
         if ( t_log(n) > time .and. t_log(n) < time_new ) &
              dt_hit(3)=t_log(n)-time
         if ( t_spec(n) > time .and. t_spec(n) < time_new ) &
              dt_hit(4)=t_spec(n)-time
         if ( t_cmb(n) > time .and. t_cmb(n) < time_new ) &
              dt_hit(5)=t_cmb(n)-time
         if ( t_movie(n) > time .and. t_movie(n) < time_new ) &
              dt_hit(6)=t_movie(n)-time
         if ( t_TO(n) > time .and. t_TO(n) < time_new ) &
              dt_hit(7)=t_TO(n)-time
         if ( t_TOmovie(n) > time .and. t_TOmovie(n) < time_new ) &
              dt_hit(7)=t_TOmovie(n)-time
      end do

      do n=1,n_dt_hit
         if ( dt_hit(n) /= 0.0_cp .and. dt_hit(n) < dt_new ) then
            l_new_dt=.true.
            dt_new=dt_hit(n)
         end if
      end do

      if ( l_new_dt ) then
         if ( dt_new < dtMin ) dt_new=dtMin
         time_new=time+dt_new
         write(*, '(/," ! TIME STEP CHANGED TO HIT TIME:",1p,2ES16.6)') &
         &     time_new*tScale,time*tScale
         if ( rank == 0 ) then
            if ( l_save_out ) then
               open(n_log_file, file=log_file, status='unknown', position='append')
               write(n_log_file,                                                &
                    &     '(/," ! TIME STEP CHANGED TO HIT TIME:",1p,2ES16.6)') &
                    &     time_new*tScale,time*tScale
               close(n_log_file)
            else
               write(n_log_file,                                               &
                    &    '(/," ! TIME STEP CHANGED TO HIT TIME:",1p,2ES16.6)') &
                    &    time_new*tScale,time*tScale
            end if
         end if
      end if

   end subroutine check_time_hits
!------------------------------------------------------------------------------
end module step_time_mod
