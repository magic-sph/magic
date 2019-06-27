module storeCheckPoints
   !
   ! This module contains several subroutines that can be used to store the
   ! checkpoint_#.tag files
   !

   use precision_mod
   use parallel_mod, only: rank
   use communications, only: gt_OC, gt_IC, gather_all_from_lo_to_rank0
   use truncation, only: n_r_max,n_r_ic_max,minc,nalias,n_theta_max,n_phi_tot, &
       &                 lm_max,lm_maxMag,n_r_maxMag,n_r_ic_maxMag,l_max,      &
       &                 fd_stretch, fd_ratio
   use radial_functions, only: rscheme_oc
   use physical_parameters, only: ra, pr, prmag, radratio, ek, sigma_ratio, &
       &                          raxi, sc
   use LMLoop_data, only: llm,ulm, llmMag, ulmMag
   use num_param, only: tScale, alph1, alph2
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

   subroutine store(time,dt,dtNew,n_time_step,l_stop_time,l_new_rst_file,     &
              &     l_ave_file,w,z,p,s,xi,b,aj,b_ic,aj_ic,dwdtLast,dzdtLast,  &
              &     dpdtLast,dsdtLast,dxidtLast,dbdtLast,djdtLast,dbdt_icLast,&
              &     djdt_icLast)
      !
      ! This subroutine stores the results in a checkpoint file.
      ! In addition to the magnetic field and velocity potentials
      ! we also store the time derivative terms djdt(lm,nR),dbdt(lm,nR), ...
      ! to allow to restart with 2nd order Adams-Bashforth scheme.
      ! To minimize the memory imprint, a gather/write strategy has been adopted
      ! here. This implies that only one global array dimension(lm_max,n_r_max)
      ! is required.
      !

      !-- Input of variables:
      real(cp),    intent(in) :: time,dt,dtNew
      integer,     intent(in) :: n_time_step
      logical,     intent(in) :: l_stop_time
      logical,     intent(in) :: l_ave_file
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
      complex(cp), allocatable :: work(:,:)

      integer :: n_rst_file, version
      character(len=72) :: string,rst_file

      version = 1

      if ( l_ave_file ) then
         rst_file="checkpoint_ave."//tag
      else
         if ( l_stop_time .or. .not.l_new_rst_file ) then
            rst_file="checkpoint_end."//tag
         else if ( l_new_rst_file ) then
            call dble2str(time,string)
            rst_file='checkpoint_t='//trim(string)//'.'//tag
         end if
      end if

      if ( rank == 0 ) then
         open(newunit=n_rst_file, file=rst_file, status='unknown', &
         &    form='unformatted', access='stream')

         !-- Write the header of the file
         write(n_rst_file) version
         write(n_rst_file) time*tScale,dt*tScale,n_time_step
         write(n_rst_file) ra,pr,raxi,sc,prmag,ek,radratio,sigma_ratio
         write(n_rst_file) n_r_max,n_theta_max,n_phi_tot,minc,nalias, &
         &                 n_r_ic_max

         !-- Store radius and scheme version (FD or CHEB)
         if ( rscheme_oc%version == 'cheb' ) then
            write(n_rst_file) rscheme_oc%version, rscheme_oc%n_max, &
            &                 rscheme_oc%order_boundary, alph1, alph2
         else
            write(n_rst_file) rscheme_oc%version, rscheme_oc%order, &
            &                 rscheme_oc%order_boundary, fd_stretch, fd_ratio
         end if

         !-- Store Lorentz-torques and rotation rates:
         write(n_rst_file) lorentz_torque_icLast,lorentz_torque_maLast, &
         &                 omega_ic1,omegaOsz_ic1,tOmega_ic1,           &
         &                 omega_ic2,omegaOsz_ic2,tOmega_ic2,           &
         &                 omega_ma1,omegaOsz_ma1,tOmega_ma1,           &
         &                 omega_ma2,omegaOsz_ma2,tOmega_ma2,           &
         &                 dtNew

         !-- Write logical to know how many fields are stored
         write(n_rst_file) l_heat, l_chemical_conv, l_mag, l_cond_ic
      end if

      !-- Memory allocation of global arrays to write outputs
      if ( rank == 0 ) then
         allocate( work(lm_max,n_r_max) )
      else
         allocate( work(1,1) )
      end if

      !-- Gather fields on rank 0 and write

      !-- Poloidal flow
      call gather_all_from_lo_to_rank0(gt_OC,w,work)
      if ( rank == 0 ) write(n_rst_file) work
      call gather_all_from_lo_to_rank0(gt_OC,dwdtLast,work)
      if ( rank == 0 ) write(n_rst_file) work
      !-- Toroidal flow
      call gather_all_from_lo_to_rank0(gt_OC,z,work)
      if ( rank == 0 ) write(n_rst_file) work
      call gather_all_from_lo_to_rank0(gt_OC,dzdtLast,work)
      if ( rank == 0 ) write(n_rst_file) work
      !-- Pressure
      call gather_all_from_lo_to_rank0(gt_OC,p,work)
      if ( rank == 0 ) write(n_rst_file) work
      call gather_all_from_lo_to_rank0(gt_OC,dpdtLast,work)
      if ( rank == 0 ) write(n_rst_file) work
      !-- Entropy
      if ( l_heat) then
         call gather_all_from_lo_to_rank0(gt_OC,s,work)
         if ( rank == 0 ) write(n_rst_file) work
         call gather_all_from_lo_to_rank0(gt_OC,dsdtLast,work)
         if ( rank == 0 ) write(n_rst_file) work
      end if
      !-- Chemical composition
      if ( l_chemical_conv  ) then
         call gather_all_from_lo_to_rank0(gt_OC,xi,work)
         if ( rank == 0 ) write(n_rst_file) work
         call gather_all_from_lo_to_rank0(gt_OC,dxidtLast,work)
         if ( rank == 0 ) write(n_rst_file) work
      end if

      !-- Outer core magnetic field
      if ( l_mag ) then
         call gather_all_from_lo_to_rank0(gt_OC,b,work)
         if ( rank == 0 ) write(n_rst_file) work
         call gather_all_from_lo_to_rank0(gt_OC,dbdtLast,work)
         if ( rank == 0 ) write(n_rst_file) work
         call gather_all_from_lo_to_rank0(gt_OC,aj,work)
         if ( rank == 0 ) write(n_rst_file) work
         call gather_all_from_lo_to_rank0(gt_OC,djdtLast,work)
         if ( rank == 0 ) write(n_rst_file) work
      end if

      !-- Inner core magnetic field
      if ( l_mag .and. l_cond_ic ) then
         deallocate( work )
         if ( rank == 0 ) then
            allocate ( work(lm_max, n_r_ic_max) )
         else
            allocate ( work(1,1) )
         end if

         call gather_all_from_lo_to_rank0(gt_IC,b_ic,work)
         if ( rank == 0 ) write(n_rst_file) work
         call gather_all_from_lo_to_rank0(gt_IC,dbdt_icLast,work)
         if ( rank == 0 ) write(n_rst_file) work
         call gather_all_from_lo_to_rank0(gt_IC,aj_ic,work)
         if ( rank == 0 ) write(n_rst_file) work
         call gather_all_from_lo_to_rank0(gt_IC,djdt_icLast,work)
         if ( rank == 0 ) write(n_rst_file) work
      end if

      !-- Deallocate work array
      deallocate( work)


      !-- Close checkpoint file and display a message in the log file
      if ( rank == 0 ) then

         close(n_rst_file)

         write(*,'(/,1P,A,/,A,ES20.10,/,A,I15,/,A,A)')&
         &    " ! Storing checkpoint file:",          &
         &    "             at time=",time,           &
         &    "            step no.=",n_time_step,    &
         &    "           into file=",rst_file

         if ( l_save_out ) then
            open(newunit=n_log_file, file=log_file, status='unknown', &
            &    position='append')
         end if
         write(n_log_file,'(/,1P,A,/,A,ES20.10,/,A,I15,/,A,A)') &
         &    " ! Storing checkpoint file:",                    &
         &    "             at time=",time,                     &
         &    "            step no.=",n_time_step,              &
         &    "           into file=",rst_file
         if ( l_save_out ) close(n_log_file)

      end if


   end subroutine store
!----------------------------------------------------------------------
end module storeCheckPoints
