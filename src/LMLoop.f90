#include "perflib_preproc.cpp"
module LMLoop_mod

#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif

   use fields
   use fieldsLast
   use omp_lib
   use precision_mod
   use parallel_mod, only: rank
   use truncation, only: l_max, lm_max, n_r_max, n_r_maxMag
   use radial_data, only: n_r_icb, n_r_cmb
   use blocking, only: lmStartB, lmStopB
   use logic, only: l_mag, l_conv, l_anelastic_liquid, lVerbose, l_heat
   use matrices, only: lZ10mat, lSmat, lZmat, lWPmat, lBmat
   use output_data, only: nLF, log_file
   use timing, only: wallTime,subTime,writeTime
   use LMLoop_data, only: llm, ulm, llmMag, ulmMag
   use debugging,  only: debug_write
   use communications, only: GET_GLOBAL_SUM, lo2r_redist_start, &
                            lo2r_s, lo2r_z, lo2r_p, lo2r_b,     &
                            lo2r_aj, lo2r_w
   use updateS_mod, only: initialize_updateS,updateS,updateS_ala
   use updateZ_mod, only: initialize_updateZ,updateZ
   use updateWP_mod, only: initialize_updateWP,updateWP
   use updateB_mod, only: initialize_updateB,updateB
   use useful, only: safeOpen, safeClose

   implicit none

   private

   public :: LMLoop,initialize_LMLoop

contains

   subroutine initialize_LMLoop

      call initialize_updateS
      call initialize_updateZ
      call initialize_updateWP
      call initialize_updateB

   end subroutine initialize_LMLoop
!----------------------------------------------------------------------------
   subroutine LMLoop(w1,coex,time,dt,lMat,lRmsNext,               &
       &            dVxBhLM,dVSrLM,dsdt,dwdt,dzdt,dpdt,dbdt,djdt, &
       &            lorentz_torque_ma,lorentz_torque_ic,          &
       &            b_nl_cmb,aj_nl_cmb,aj_nl_icb,n_time_step)
      !
      !  This subroutine performs the actual time-stepping.
      !
      !

      !-- Input of variables:
      real(cp),    intent(in) :: w1,coex
      real(cp),    intent(in) :: dt,time
      logical,     intent(in) :: lMat
      logical,     intent(in) :: lRmsNext

      !--- Input from radialLoop:
      !    These fields are provided in the R-distributed space!
      ! for djdt in update_b
      complex(cp), intent(inout) :: dVxBhLM(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(inout) :: dVSrLM(llm:ulm,n_r_max)   ! for dsdt in update_s
      integer,     intent(in) :: n_time_step

      !--- Input from radialLoop and then redistributed:
      complex(cp), intent(inout) :: dsdt(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dwdt(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dzdt(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dpdt(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dbdt(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(inout) :: djdt(llmMag:ulmMag,n_r_maxMag)
      real(cp),    intent(in) :: lorentz_torque_ma,lorentz_torque_ic
      complex(cp), intent(in) :: b_nl_cmb(lm_max)   ! nonlinear bc for b at CMB
      complex(cp), intent(in)  :: aj_nl_cmb(lm_max)  ! nonlinear bc for aj at CMB
      complex(cp), intent(in)  :: aj_nl_icb(lm_max)  ! nonlinear bc for dr aj at ICB

      !--- Local counter
      integer :: nLMB
      integer :: lmStart,lmStop
      integer :: l
      integer :: tStart(4),tStop(4),tPassed(4)

      logical,parameter :: DEBUG_OUTPUT=.false.
      !--- Inner core rotation from last time step
      real(cp), save :: omega_icLast
      complex(cp) :: sum_dwdt

      PERFON('LMloop')
      !LIKWID_ON('LMloop')
      if ( lVerbose ) call safeOpen(nLF,log_file)

      omega_icLast=omega_ic

      if ( lMat ) then ! update matrices:
      !---- The following logicals tell whether the respective inversion
      !     matrices have been updated. lMat=.true. when a general
      !     update is necessary. These logicals are THREADPRIVATE and
      !     stored in the module matrices in m_mat.F90:
         lZ10mat=.false.
         do l=0,l_max
            lSmat(l) =.false.
            lZmat(l) =.false.
            lWPmat(l)=.false.
            lBmat(l) =.false.
         end do
      end if

      !nThreadsLMmax = 1
      nLMB=1+rank
      !nTh=1
      if ( lVerbose ) then
         write(*,'(/," ! lm block no:",i3)') nLMB
         call wallTime(tStart)
      end if

      if ( DEBUG_OUTPUT ) then
         lmStart=lmStartB(nLMB)
         lmStop =lmStopB(nLMB)
      end if
      !call debug_write(dwdt(lmStart:lmStop,:),lmStop-lmStart+1,n_r_max, &
      !                                        "dwdt",n_time_step*1000+nLMB*100,"E")

      if ( l_heat ) then ! dp,workA usead as work arrays
         if ( DEBUG_OUTPUT ) then
            write(*,"(A,I2,8ES20.12)") "s_before: ",nLMB,   &
                 & GET_GLOBAL_SUM( s_LMloc(:,:) ),          &
                 & GET_GLOBAL_SUM( ds_LMloc(:,:) ),         &
                 & GET_GLOBAL_SUM( dsdt(:,:) ),             &
                 & GET_GLOBAL_SUM( dsdtLast_LMloc(:,:) )
         end if
         !call debug_write(dsdt,ulm-llm+1,n_r_max,"dsdt_LMloc", &
         !                        n_time_step*1000+nLMB*100,"E")
         PERFON('up_S')
         if ( l_anelastic_liquid ) then
            call updateS_ala(s_LMloc,ds_LMloc,w_LMloc,dVSrLM,dsdt,    & 
                 &       dsdtLast_LMloc,w1,coex,dt,nLMB)
         else
            call updateS(s_LMloc,ds_LMloc,dVSrLM,dsdt,dsdtLast_LMloc, &
                 &       w1,coex,dt,nLMB)

         end if
         PERFOFF
         ! Here one could start the redistribution of s_LMloc,ds_LMloc etc. with a 
         ! nonblocking send
         !call MPI_Barrier(MPI_COMM_WORLD,ierr)
         !PERFON('rdstSst')
         call lo2r_redist_start(lo2r_s,s_LMloc_container,s_Rloc_container)
         !PERFOFF

         if ( DEBUG_OUTPUT ) then
            write(*,"(A,I2,4ES20.12)") "s_after : ",nLMB,  &
                 & GET_GLOBAL_SUM( s_LMloc(:,:) ),         &
                 & GET_GLOBAL_SUM( ds_LMloc(:,:) )
            write(*,"(A,I2,8ES22.14)") "s_after(bnd_r): ",nLMB, &
                 & GET_GLOBAL_SUM( s_LMloc(:,n_r_icb) ),        &
                 & GET_GLOBAL_SUM( s_LMloc(:,n_r_cmb) ),        &
                 & GET_GLOBAL_SUM( ds_LMloc(:,n_r_icb) ),       &
                 & GET_GLOBAL_SUM( ds_LMloc(:,n_r_cmb) )
         end if
      end if
      if ( l_conv ) then
         if ( DEBUG_OUTPUT ) then
            write(*,"(A,I2,6ES20.12)") "z_before: ",nLMB,   &
                 & GET_GLOBAL_SUM( z_LMloc(:,:) ),          &
                 & GET_GLOBAL_SUM( dz_LMloc(:,:) ),         &
                 & GET_GLOBAL_SUM( dzdtLast_lo(:,:) )
         end if
         PERFON('up_Z')
         ! dp, dVSrLM, workA used as work arrays
         !call updateZ( z_LMloc, dz_LMloc, dzdt, dzdtLast_LMloc, time, &
         call updateZ( z_LMloc, dz_LMloc, dzdt, dzdtLast_lo, time, &
              &        omega_ma,d_omega_ma_dtLast,                 &
              &        omega_ic,d_omega_ic_dtLast,                 &
              &        lorentz_torque_ma,lorentz_torque_maLast,    &
              &        lorentz_torque_ic,lorentz_torque_icLast,    &
              &        w1,coex,dt,lRmsNext)
         PERFOFF

         !call MPI_Barrier(MPI_COMM_WORLD,ierr)
         !PERFON('rdstZst')
         call lo2r_redist_start(lo2r_z,z_LMloc_container,z_Rloc_container)
         !PERFOFF

         if ( DEBUG_OUTPUT ) then
            !do lm=lmStart,lmStop
            !   write(*,"(A,I4,6ES20.12)") "z_after : ",lm,SUM( z_LMloc(lm,:) ),&
            !        & SUM( dz_LMloc(lm,:) ),SUM( dzdtLast_lo(lm,:) )
            !end do
            write(*,"(A,I2,6ES20.12)") "z_after: ",nLMB,  &
                 & GET_GLOBAL_SUM( z_LMloc(:,:) ),        &
                 & GET_GLOBAL_SUM( dz_LMloc(:,:) ),       &
                 & GET_GLOBAL_SUM( dzdtLast_lo(:,:) )
         end if
         ! dVSrLM, workA used as work arrays
         !call debug_write(dwdt(:,:),lmStop-lmStart+1,n_r_max, &
         !                "dwdt",n_time_step*1000+nLMB*100,"E")
         if ( DEBUG_OUTPUT ) then
            sum_dwdt=GET_GLOBAL_SUM( dwdt(:,:) )
            write(*,"(A,I2,8ES22.14,4(I3,F19.16))") "wp_before: ",nLMB,&
                 & GET_GLOBAL_SUM( w_LMloc(:,:) ),                     &
                 & GET_GLOBAL_SUM( p_LMloc(:,:) ),                     &
                 & GET_GLOBAL_SUM( dwdtLast_LMloc(:,:) ),              &
                 & GET_GLOBAL_SUM( dpdtLast_LMloc(:,:) ),              &
                 & exponent(real(sum_dwdt)),fraction(real(sum_dwdt)),  &
                 & exponent(aimag(sum_dwdt)),fraction(aimag(sum_dwdt))
         end if
         PERFON('up_WP')
         call updateWP( w_LMloc, dw_LMloc, ddw_LMloc, dwdt, dwdtLast_LMloc, &
              &         p_LMloc, dp_LMloc, dpdt, dpdtLast_LMloc, s_LMloc,   &
              &         w1,coex,dt,nLMB,lRmsNext)
         PERFOFF

         !call MPI_Barrier(MPI_COMM_WORLD,ierr)
         !PERFON('rdstWPst')
         call lo2r_redist_start(lo2r_w,w_LMloc_container,w_Rloc_container)
         call lo2r_redist_start(lo2r_p,p_LMloc_container,p_Rloc_container)
         !PERFOFF

         if ( DEBUG_OUTPUT ) then
            write(*,"(A,I2,12ES22.14)") "wp_after: ",nLMB,  &
                 & GET_GLOBAL_SUM( w_LMloc(:,:) ),          &
                 & GET_GLOBAL_SUM( p_LMloc(:,:) ),          &
                 & GET_GLOBAL_SUM( dwdtLast_LMloc(:,:) ),   &
                 & GET_GLOBAL_SUM( dpdtLast_LMloc(:,:) ),   &
                 &GET_GLOBAL_SUM( dw_LMloc(:,:) )
            write(*,"(A,I2,4ES22.14)") "wp_after(bnd_r): ",nLMB, &
                 & GET_GLOBAL_SUM( w_LMloc(:,n_r_icb) ),         &
                 & GET_GLOBAL_SUM( w_LMloc(:,n_r_cmb) )
         end if
      end if
      if ( l_mag ) then ! dwdt,dpdt used as work arrays
         if ( DEBUG_OUTPUT ) then
            write(*,"(A,I2,14ES20.12)") "b_before: ",nLMB,   &
                 & GET_GLOBAL_SUM(  b_LMloc(:,:) ),          & 
                 & GET_GLOBAL_SUM( aj_LMloc(:,:) ),          &
                 & GET_GLOBAL_SUM( b_ic_LMloc(:,:) ),        &
                 & GET_GLOBAL_SUM( aj_ic_LMloc(:,:) ),       &
                 & GET_GLOBAL_SUM( dbdt_icLast_LMloc(:,:) ), &
                 & GET_GLOBAL_SUM( djdt_icLast_LMloc(:,:) ), &
                 & GET_GLOBAL_SUM( dVxBhLM(:,:) )
         end if
         !LIKWID_ON('up_B')
         PERFON('up_B')
         call updateB( b_LMloc,db_LMloc,ddb_LMloc,aj_LMloc,dj_LMloc,ddj_LMloc, &
              &        dVxBhLM, dbdt, dbdtLast_LMloc, djdt, djdtLast_LMloc,    &
              &        b_ic_LMloc, db_ic_LMloc, ddb_ic_LMloc, aj_ic_LMloc,     &
              &        dj_ic_LMloc, ddj_ic_LMloc, dbdt_icLast_LMloc,           &
              &        djdt_icLast_LMloc, b_nl_cmb, aj_nl_cmb, aj_nl_icb,      &
              &        omega_icLast, w1, coex, dt, time, nLMB, lRmsNext )
         PERFOFF
         !LIKWID_OFF('up_B')
         call lo2r_redist_start(lo2r_b, b_LMloc_container, b_Rloc_container)
         call lo2r_redist_start(lo2r_aj, aj_LMloc_container, aj_Rloc_container)

         if ( DEBUG_OUTPUT ) then
            write(*,"(A,I2,8ES20.12)") "b_after: ",nLMB, &
                 & GET_GLOBAL_SUM(  b_LMloc(:,:) ),      & 
                 & GET_GLOBAL_SUM( aj_LMloc(:,:) ),      &
                 & GET_GLOBAL_SUM( dbdtLast_LMloc(:,:) ),&
                 & GET_GLOBAL_SUM( djdtLast_LMloc(:,:) )
            write(*,"(A,I2,8ES20.12)") "b_ic_after: ",nLMB, &
                 & GET_GLOBAL_SUM( b_ic_LMloc(:,:) ),       &
                 & GET_GLOBAL_SUM( aj_ic_LMloc(:,:) ),      &
                 & GET_GLOBAL_SUM( dbdt_icLast_LMloc(:,:) ),&
                 & GET_GLOBAL_SUM( djdt_icLast_LMloc(:,:) )
         end if
      end if

      if ( lVerbose ) then
         call wallTime(tStop)
         call subTime(tStart,tStop,tPassed)
         !write(nLF,*) '! Thread no:',nTh
         write(nLF,*) 'lmStart,lmStop:',lmStartB(nLMB),lmStopB(nLMB)
         call writeTime(nLF,'! Time for thread:',tPassed)
      end if


      lorentz_torque_maLast=lorentz_torque_ma
      lorentz_torque_icLast=lorentz_torque_ic

      if ( lVerbose ) call safeClose(nLF)

      !LIKWID_OFF('LMloop')
      PERFOFF
   end subroutine LMLoop
!--------------------------------------------------------------------------------
end module LMLoop_mod
