module storeCheckPoints
   !
   ! This module contains several subroutines that can be used to store the
   ! rst_#.tag files
   !

   use precision_mod
   use parallel_mod, only: rank
   use communications, only: gt_OC, gt_IC, gather_all_from_lo_to_rank0
   use truncation, only: n_r_max,n_r_ic_max,minc,nalias,n_theta_max,n_phi_tot, &
       &                 lm_max,lm_maxMag,n_r_maxMag,n_r_ic_maxMag,l_max,      &
       &                 fd_stretch, fd_ratio
   use radial_functions, only: rscheme_oc, alph1, alph2
   use physical_parameters, only: ra, pr, prmag, radratio, ek, sigma_ratio, &
       &                          raxi, sc
   use LMLoop_data, only: llm,ulm, llmMag, ulmMag
   use num_param, only: tScale
   use fieldsLast, only: d_omega_ma_dtLast,d_omega_ic_dtLast, &
       &                 lorentz_torque_maLast,lorentz_torque_icLast
   use init_fields, only: inform,omega_ic1,omegaOsz_ic1,tOmega_ic1, &
       &                  omega_ic2,omegaOsz_ic2,tOmega_ic2,        &
       &                  omega_ma1,omegaOsz_ma1,tOmega_ma1,        &
       &                  omega_ma2,omegaOsz_ma2,tOmega_ma2
   use logic, only: l_heat, l_mag, l_cond_ic, l_chemical_conv, l_save_out
   use output_data, only: tag, log_file, n_log_file
   use charmanip, only: dble2str

   implicit none

   private
 
   public :: store

contains

   subroutine store(time,dt,dtNew,n_time_step,l_stop_time,l_new_rst_file,  &
              &     w,z,p,s,xi,b,aj,b_ic,aj_ic,dwdtLast,dzdtLast,dpdtLast, &
              &     dsdtLast,dxidtLast,dbdtLast,djdtLast,dbdt_icLast,      &
              &     djdt_icLast)
      !
      ! store results on disc file (restart file)
      ! In addition to the magnetic field and velocity potentials
      ! we store the time derivative terms
      ! djdt(lm,nR),dbdt(lm,nR), ......
      !

      !-- Input of variables:
      real(cp),    intent(in) :: time,dt,dtNew
      integer,     intent(in) :: n_time_step
      logical,     intent(in) :: l_stop_time
      logical,     intent(in) :: l_new_rst_file

      !-- Input of scalar fields to be stored:
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(in) :: z(llm:ulm,n_r_max)
      complex(cp), intent(in) :: p(llm:ulm,n_r_max)
      complex(cp), intent(in) :: s(llm:ulm,n_r_max)
      complex(cp), intent(in) :: xi(llm:ulm,n_r_max)
      complex(cp), intent(in) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: aj(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: dwdtLast(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dzdtLast(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dpdtLast(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dsdtLast(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dxidtLast(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dbdtLast(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: djdtLast(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: dbdt_icLast(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: djdt_icLast(llmMag:ulmMag,n_r_ic_maxMag)

      !-- Local variables
      complex(cp), allocatable :: workA(:,:), workB(:,:), workC(:,:)
      complex(cp), allocatable :: workD(:,:), workE(:,:)

      integer :: n_rst_file
      character(len=72) :: string,rst_file

      if ( l_stop_time .or. .not.l_new_rst_file ) then
         rst_file="rst_end."//tag
      else if ( l_new_rst_file ) then
         call dble2str(time,string)
         rst_file='rst_t='//trim(string)//'.'//tag
      end if

      !-- Write parameters:
      if ( .not. l_chemical_conv ) then
         if ( .not. l_heat ) then
            inform=21
         else
            inform=22
         end if
      else
         if ( .not. l_heat ) then
            inform=23
         else
            inform=24
         end if
      end if

      if ( rank == 0 ) then
         open(newunit=n_rst_file, file=rst_file, status='unknown', &
         &    form='unformatted')

         !-- Write the header of the file
         write(n_rst_file) time*tScale,dt*tScale,ra,pr,prmag,ek,radratio, &
         &              inform,n_r_max,n_theta_max,n_phi_tot,minc,nalias, &
         &                                        n_r_ic_max,sigma_ratio

         !-- Store radius and scheme version (FD or CHEB)
         if ( rscheme_oc%version == 'cheb' ) then
            write(n_rst_file) rscheme_oc%version, rscheme_oc%n_max, &
            &                 rscheme_oc%order_boundary, alph1, alph2
         else
            write(n_rst_file) rscheme_oc%version, rscheme_oc%order, &
            &                 rscheme_oc%order_boundary, fd_stretch, fd_ratio
         end if
      end if

      !-- Memory allocation of global arrays to write outputs
      if ( rank == 0 ) then
         allocate( workA(lm_max,n_r_max), workB(lm_max,n_r_max) )
         allocate( workC(lm_max,n_r_max) )
         if ( l_chemical_conv .or. l_heat .or. l_mag ) allocate( workD(lm_max,n_r_max) )
         if ( l_chemical_conv .and. l_heat ) allocate( workE(lm_max,n_r_max) )
      else
         allocate( workA(1,1), workB(1,1), workC(1,1) )
         if ( l_chemical_conv .or. l_heat .or. l_mag ) allocate( workD(1,1) )
         if ( l_chemical_conv .and. l_heat ) allocate( workE(1,1) )
      end if

      !-- Gather fields on rank 0
      call gather_all_from_lo_to_rank0(gt_OC,w,workA)
      call gather_all_from_lo_to_rank0(gt_OC,z,workB)
      call gather_all_from_lo_to_rank0(gt_OC,p,workC)
      if ( .not. l_chemical_conv  ) then
         if ( l_heat ) call gather_all_from_lo_to_rank0(gt_OC,s,workD)
      else
         if ( l_heat ) then
            call gather_all_from_lo_to_rank0(gt_OC,s,workD)
            call gather_all_from_lo_to_rank0(gt_OC,xi,workE)
         else
            call gather_all_from_lo_to_rank0(gt_OC,xi,workD)
         end if
      end if

      !-- Write output
      if ( rank == 0 ) then
         if ( .not. l_chemical_conv ) then
            if ( l_heat ) then
               write(n_rst_file) workA,workB,workC,workD
            else
               write(n_rst_file) workA,workB,workC
            end if
         else
            if ( l_heat ) then
               write(n_rst_file) workA,workB,workC,workD,workE
            else
               write(n_rst_file) workA,workB,workC,workD
            end if
         end if
      end if

      !-- Gather d?/dt fields on rank 0
      call gather_all_from_lo_to_rank0(gt_OC,dwdtLast,workA)
      call gather_all_from_lo_to_rank0(gt_OC,dzdtLast,workB)
      call gather_all_from_lo_to_rank0(gt_OC,dpdtLast,workC)
      if ( .not. l_chemical_conv  ) then
         if ( l_heat ) call gather_all_from_lo_to_rank0(gt_OC,dsdtLast,workD)
      else
         if ( l_heat ) then
            call gather_all_from_lo_to_rank0(gt_OC,dsdtLast,workD)
            call gather_all_from_lo_to_rank0(gt_OC,dxidtLast,workE)
         else
            call gather_all_from_lo_to_rank0(gt_OC,dxidtLast,workD)
         end if
      end if

      !-- Write output
      if ( rank == 0 ) then
         if ( .not. l_chemical_conv ) then
            if ( l_heat ) then
               write(n_rst_file) workD,workA,workB,workC
            else
               write(n_rst_file) workA,workB,workC
            end if
         else
            if ( l_heat ) then
               write(n_rst_file) workD,workA,workB,workC,workE
            else
               write(n_rst_file) workA,workB,workC,workD
            end if
         end if

         if ( l_chemical_conv ) then
            write(n_rst_file) raxi,sc
         end if
      end if

      if ( l_chemical_conv .and. l_heat ) deallocate( workE )

      if ( l_mag ) then
         !-- Gather magnetic field
         call gather_all_from_lo_to_rank0(gt_OC,b,workA)
         call gather_all_from_lo_to_rank0(gt_OC,aj,workB)
         call gather_all_from_lo_to_rank0(gt_OC,dbdtLast,workC)
         call gather_all_from_lo_to_rank0(gt_OC,djdtLast,workD)

         !-- Write magnetic field:
         if ( rank == 0 ) write(n_rst_file) workA,workB,workC,workD

         !-- Inner core
         if ( l_cond_ic ) then
            deallocate( workA, workB, workC, workD )
            allocate( workA(lm_max,n_r_ic_max), workB(lm_max,n_r_ic_max) )
            allocate( workC(lm_max,n_r_ic_max), workD(lm_max,n_r_ic_max) )

            !-- Gather inner core magnetic field
            call gather_all_from_lo_to_rank0(gt_IC,b_ic,workA)
            call gather_all_from_lo_to_rank0(gt_IC,aj_ic,workB)
            call gather_all_from_lo_to_rank0(gt_IC,dbdt_icLast,workC)
            call gather_all_from_lo_to_rank0(gt_IC,djdt_icLast,workD)

            !-- Write IC magnetic field:
            if ( rank == 0 ) write(n_rst_file) workA,workB,workC,workD
         end if

      end if

      if ( rank == 0 ) then
         !-- Store Lorentz-torques and rotation rates:
         write(n_rst_file) lorentz_torque_icLast,lorentz_torque_maLast, &
         &                 omega_ic1,omegaOsz_ic1,tOmega_ic1,           &
         &                 omega_ic2,omegaOsz_ic2,tOmega_ic2,           &
         &                 omega_ma1,omegaOsz_ma1,tOmega_ma1,           &
         &                 omega_ma2,omegaOsz_ma2,tOmega_ma2,           &
         &                 dtNew

         close(n_rst_file)

         write(*,'(/,1P,A,/,A,ES20.10,/,A,I15,/,A,A)')&
         &    " ! Storing restart file:",             &
         &    "             at time=",time,           &
         &    "            step no.=",n_time_step,    &
         &    "           into file=",rst_file

         if ( l_save_out ) then
            open(newunit=n_log_file, file=log_file, status='unknown', &
            &    position='append')
         end if
         write(n_log_file,'(/,1P,A,/,A,ES20.10,/,A,I15,/,A,A)') &
         &    " ! Storing restart file:",                       &
         &    "             at time=",time,                     &
         &    "            step no.=",n_time_step,              &
         &    "           into file=",rst_file
         if ( l_save_out ) close(n_log_file)

      end if

      deallocate( workA,workB,workC )
      if ( l_chemical_conv .or. l_heat .or. l_mag ) deallocate( workD )



   end subroutine store
!----------------------------------------------------------------------
end module storeCheckPoints
