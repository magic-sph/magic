module storeCheckPoints
   !
   ! This module contains several subroutines that can be used to store the
   ! checkpoint_#.tag files
   !

   use precision_mod
   use parallel_mod
   use mpi_alltoall_mod, only: type_mpiatoav
   use communications, only: gt_OC, gt_IC, gather_all_from_lo_to_rank0
   use truncation, only: n_r_max,n_r_ic_max,minc,nalias,n_theta_max,n_phi_tot, &
       &                 lm_max,lm_maxMag,n_r_maxMag,n_r_ic_maxMag,l_max,      &
       &                 fd_stretch, fd_ratio, nRstart, nRstop, nRstartMag,    &
       &                 nRstopMag, n_r_loc
   use radial_functions, only: rscheme_oc, r
   use physical_parameters, only: ra, pr, prmag, radratio, ek, sigma_ratio, &
       &                          raxi, sc
   use blocking, only: llm, ulm, llmMag, ulmMag
   use num_param, only: tScale, alph1, alph2
   use init_fields, only: inform,omega_ic1,omegaOsz_ic1,tOmega_ic1, &
       &                  omega_ic2,omegaOsz_ic2,tOmega_ic2,        &
       &                  omega_ma1,omegaOsz_ma1,tOmega_ma1,        &
       &                  omega_ma2,omegaOsz_ma2,tOmega_ma2
   use logic, only: l_heat, l_mag, l_cond_ic, l_chemical_conv, l_save_out, &
       &            l_double_curl
   use output_data, only: tag, log_file, n_log_file
   use charmanip, only: dble2str
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray, type_tscalar

   implicit none

   private

#ifdef WITH_MPI
   public :: store, store_mpi
#else
   public :: store
#endif

contains

   subroutine store(time,tscheme,n_time_step,l_stop_time,l_new_rst_file,    &
              &     l_ave_file,w,z,p,s,xi,b,aj,b_ic,aj_ic,dwdt,dzdt,        &
              &     dpdt,dsdt,dxidt,dbdt,djdt,dbdt_ic,djdt_ic,domega_ma_dt, &
              &     domega_ic_dt,lorentz_torque_ma_dt,lorentz_torque_ic_dt)
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
      real(cp),            intent(in) :: time
      class(type_tscheme), intent(in) :: tscheme
      integer,             intent(in) :: n_time_step
      logical,             intent(in) :: l_stop_time
      logical,             intent(in) :: l_ave_file
      logical,             intent(in) :: l_new_rst_file

      !-- Input of scalar fields to be stored:
      complex(cp),         intent(in) :: w(llm:ulm,n_r_max)
      complex(cp),         intent(in) :: z(llm:ulm,n_r_max)
      complex(cp),         intent(in) :: p(llm:ulm,n_r_max)
      complex(cp),         intent(in) :: s(llm:ulm,n_r_max)
      complex(cp),         intent(in) :: xi(llm:ulm,n_r_max)
      complex(cp),         intent(in) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(cp),         intent(in) :: aj(llmMag:ulmMag,n_r_maxMag)
      complex(cp),         intent(in) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp),         intent(in) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      type(type_tarray),   intent(in) :: dwdt, dzdt, dpdt, dsdt, dxidt, dbdt
      type(type_tarray),   intent(in) :: djdt, dbdt_ic, djdt_ic
      type(type_tscalar),  intent(in) :: domega_ic_dt, domega_ma_dt
      type(type_tscalar),  intent(in) :: lorentz_torque_ic_dt, lorentz_torque_ma_dt

      !-- Local variables
      complex(cp), allocatable :: work(:,:)
      logical :: l_press_store

      integer :: n_rst_file, version, n_o
      character(len=72) :: string,rst_file

      version = 2
      l_press_store = ( .not. l_double_curl ) 

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
         write(n_rst_file) time*tScale
         write(n_rst_file) tscheme%family, tscheme%nexp, tscheme%nimp, tscheme%nold
         write(n_rst_file) tscheme%dt(:)*tScale
         write(n_rst_file) n_time_step
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

         write(n_rst_file) r(:)

         !-- Store Lorentz-torques and rotation rates:
         if ( tscheme%family == 'MULTISTEP' ) then
            do n_o=2,tscheme%nexp
               write(n_rst_file) domega_ic_dt%expl(n_o)
            end do
            do n_o=2,tscheme%nimp
               write(n_rst_file) domega_ic_dt%impl(n_o)
            end do
            do n_o=2,tscheme%nold
               write(n_rst_file) domega_ic_dt%old(n_o)
            end do
            do n_o=2,tscheme%nexp
               write(n_rst_file) domega_ma_dt%expl(n_o)
            end do
            do n_o=2,tscheme%nimp
               write(n_rst_file) domega_ma_dt%impl(n_o)
            end do
            do n_o=2,tscheme%nold
               write(n_rst_file) domega_ma_dt%old(n_o)
            end do

            do n_o=2,tscheme%nexp
               write(n_rst_file) lorentz_torque_ic_dt%expl(n_o)
            end do
            do n_o=2,tscheme%nimp
               write(n_rst_file) lorentz_torque_ic_dt%impl(n_o)
            end do
            do n_o=2,tscheme%nold
               write(n_rst_file) lorentz_torque_ic_dt%old(n_o)
            end do
            do n_o=2,tscheme%nexp
               write(n_rst_file) lorentz_torque_ma_dt%expl(n_o)
            end do
            do n_o=2,tscheme%nimp
               write(n_rst_file) lorentz_torque_ma_dt%impl(n_o)
            end do
            do n_o=2,tscheme%nold
               write(n_rst_file) lorentz_torque_ma_dt%old(n_o)
            end do
         end if
         write(n_rst_file) omega_ic1,omegaOsz_ic1,tOmega_ic1,     &
         &                 omega_ic2,omegaOsz_ic2,tOmega_ic2,     &
         &                 omega_ma1,omegaOsz_ma1,tOmega_ma1,     &
         &                 omega_ma2,omegaOsz_ma2,tOmega_ma2

         !-- Write logical to know how many fields are stored
         write(n_rst_file) l_heat, l_chemical_conv, l_mag, l_press_store, &
         &                 l_cond_ic
      end if

      !-- Memory allocation of global arrays to write outputs
      if ( rank == 0 ) then
         allocate( work(lm_max,n_r_max) )
      else
         allocate( work(1,1) )
      end if

      !-- Gather fields on rank 0 and write

      !-- Poloidal flow
      call write_one_field(n_rst_file, tscheme, w, dwdt, work)

      !-- Toroidal flow
      call write_one_field(n_rst_file, tscheme, z, dzdt, work)

      !-- Pressure
      if ( l_press_store ) call write_one_field(n_rst_file, tscheme, p, dpdt, work)

      !-- Entropy
      if ( l_heat) call write_one_field(n_rst_file, tscheme, s, dsdt, work)

      !-- Chemical composition
      if ( l_chemical_conv) call write_one_field(n_rst_file, tscheme, xi, dxidt, work)

      !-- Outer core magnetic field
      if ( l_mag ) then
         call write_one_field(n_rst_file, tscheme, b, dbdt, work)
         call write_one_field(n_rst_file, tscheme, aj, djdt, work)
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
         if ( tscheme%family == 'MULTISTEP' ) then
            do n_o=2,tscheme%nexp
               call gather_all_from_lo_to_rank0(gt_IC,dbdt_ic%expl(:,:,n_o),work)
               if ( rank == 0 ) write(n_rst_file) work
            end do
            do n_o=2,tscheme%nimp
               call gather_all_from_lo_to_rank0(gt_IC,dbdt_ic%impl(:,:,n_o),work)
               if ( rank == 0 ) write(n_rst_file) work
            end do
            do n_o=2,tscheme%nold
               call gather_all_from_lo_to_rank0(gt_IC,dbdt_ic%old(:,:,n_o),work)
               if ( rank == 0 ) write(n_rst_file) work
            end do
         end if
         call gather_all_from_lo_to_rank0(gt_IC,aj_ic,work)
         if ( rank == 0 ) write(n_rst_file) work
         if ( tscheme%family == 'MULTISTEP' ) then
            do n_o=2,tscheme%nexp
               call gather_all_from_lo_to_rank0(gt_IC,djdt_ic%expl(:,:,n_o),work)
               if ( rank == 0 ) write(n_rst_file) work
            end do
            do n_o=2,tscheme%nimp
               call gather_all_from_lo_to_rank0(gt_IC,djdt_ic%impl(:,:,n_o),work)
               if ( rank == 0 ) write(n_rst_file) work
            end do
            do n_o=2,tscheme%nold
               call gather_all_from_lo_to_rank0(gt_IC,djdt_ic%old(:,:,n_o),work)
               if ( rank == 0 ) write(n_rst_file) work
            end do
         end if

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
!-----------------------------------------------------------------------------------
   subroutine write_one_field(fh, tscheme, w, dwdt, work)

      !-- Input variables
      integer,             intent(in) :: fh ! file unit
      class(type_tscheme), intent(in) :: tscheme ! Time scheme
      complex(cp),         intent(in) :: w(llm:ulm,n_r_max) ! field
      type(type_tarray),   intent(in) :: dwdt

      !-- Output variables
      complex(cp),      intent(inout) :: work(:,:)

      !-- Local variables
      integer :: n_o

      call gather_all_from_lo_to_rank0(gt_OC, w, work)
      if ( rank == 0 ) write(fh) work

      if ( tscheme%family == 'MULTISTEP' ) then
         do n_o=2,tscheme%nexp
            call gather_all_from_lo_to_rank0(gt_OC, dwdt%expl(:,:,n_o), work)
            if ( rank == 0 ) write(fh) work
         end do

         do n_o=2,tscheme%nimp
            call gather_all_from_lo_to_rank0(gt_OC, dwdt%impl(:,:,n_o), work)
            if ( rank == 0 ) write(fh) work
         end do

         do n_o=2,tscheme%nold
            call gather_all_from_lo_to_rank0(gt_OC, dwdt%old(:,:,n_o), work)
            if ( rank == 0 ) write(fh) work
         end do
      end if

   end subroutine write_one_field
!-----------------------------------------------------------------------------------
#ifdef WITH_MPI
   subroutine store_mpi(time,tscheme,n_time_step,l_stop_time,l_new_rst_file,  &
              &         l_ave_file,w,z,p,s,xi,b,aj,b_ic,aj_ic,dwdt,dzdt,dpdt, &
              &         dsdt,dxidt,dbdt,djdt,dbdt_ic,djdt_ic,domega_ma_dt,    &
              &         domega_ic_dt,lorentz_torque_ma_dt,lorentz_torque_ic_dt)
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
      real(cp),            intent(in) :: time
      class(type_tscheme), intent(in) :: tscheme
      integer,             intent(in) :: n_time_step
      logical,             intent(in) :: l_stop_time
      logical,             intent(in) :: l_ave_file
      logical,             intent(in) :: l_new_rst_file

      !-- Input of scalar fields to be stored:
      complex(cp),         intent(in) :: w(lm_max,nRstart:nRstop)
      complex(cp),         intent(in) :: z(lm_max,nRstart:nRstop)
      complex(cp),         intent(in) :: p(lm_max,nRstart:nRstop)
      complex(cp),         intent(in) :: s(lm_max,nRstart:nRstop)
      complex(cp),         intent(in) :: xi(lm_max,nRstart:nRstop)
      complex(cp),         intent(in) :: b(lm_maxMag,nRstartMag:nRstopMag)
      complex(cp),         intent(in) :: aj(lm_maxMag,nRstartMag:nRstopMag)
      complex(cp),         intent(in) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp),         intent(in) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      type(type_tarray),   intent(in) :: dwdt, dzdt, dpdt, dsdt, dxidt, dbdt
      type(type_tarray),   intent(in) :: djdt, dbdt_ic, djdt_ic
      type(type_tscalar),  intent(in) :: domega_ic_dt, domega_ma_dt
      type(type_tscalar),  intent(in) :: lorentz_torque_ic_dt, lorentz_torque_ma_dt

      !-- Local variables
      complex(cp), allocatable :: work(:,:)

      type(type_mpiatoav) :: lo2r
      logical :: l_press_store
      integer :: version, info, fh, datatype
      character(len=72) :: string, rst_file
      integer :: istat(MPI_STATUS_SIZE)
      integer :: arr_size(2), arr_loc_size(2), arr_start(2), n_o
      integer(lip) :: disp, offset, size_tmp

      version = 2
      l_press_store = (.not. l_double_curl)

      call lo2r%create_comm(1)
      allocate( work(lm_max,nRstart:nRstop) )

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

      !-- MPI-IO setup
      call mpiio_setup(info)

      !-- Open file
      call MPI_File_Open(MPI_COMM_WORLD, rst_file, ior(MPI_MODE_WRONLY, &
           &             MPI_MODE_CREATE), info, fh, ierr)

      disp = 0
      call MPI_File_Set_View(fh, disp, MPI_BYTE, MPI_BYTE, "native", &
           &                 info, ierr)

      !-- Only rank=0 writes the header of the file
      if ( rank == 0 ) then
         !-- Write the header of the file
         call MPI_File_Write(fh, version, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, time*tScale, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, tscheme%family, len(tscheme%family), MPI_CHARACTER,&
              &              istat, ierr)
         call MPI_File_Write(fh, tscheme%nexp, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, tscheme%nimp, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, tscheme%nold, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, tscheme%dt*tScale, size(tscheme%dt), MPI_DEF_REAL, &
              &              istat, ierr)
         call MPI_File_Write(fh, n_time_step, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, ra, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, pr, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, raxi, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, sc, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, prmag, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, ek, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, radratio, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, sigma_ratio, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, n_r_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, n_theta_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, n_phi_tot, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, minc, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, nalias, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, n_r_ic_max, 1, MPI_INTEGER, istat, ierr)

         !-- Store radius and scheme version (FD or CHEB)
         call MPI_File_Write(fh, rscheme_oc%version, len(rscheme_oc%version), &
              &              MPI_CHARACTER, istat, ierr)
         if ( rscheme_oc%version == 'cheb' ) then
            call MPI_File_Write(fh, rscheme_oc%n_max, 1, MPI_INTEGER, istat, ierr)
            call MPI_File_Write(fh, rscheme_oc%order_boundary, 1, MPI_INTEGER, &
                 &              istat, ierr)
            call MPI_File_Write(fh, alph1, 1, MPI_DEF_REAL, istat, ierr)
            call MPI_File_Write(fh, alph2, 1, MPI_DEF_REAL, istat, ierr)
         else
            call MPI_File_Write(fh, rscheme_oc%order, 1, MPI_INTEGER, istat, ierr)
            call MPI_File_Write(fh, rscheme_oc%order_boundary, 1, MPI_INTEGER, &
                 &              istat, ierr)
            call MPI_File_Write(fh, fd_stretch, 1, MPI_DEF_REAL, istat, ierr)
            call MPI_File_Write(fh, fd_ratio, 1, MPI_DEF_REAL, istat, ierr)
         end if

         call MPI_File_Write(fh, r, n_r_max, MPI_DEF_REAL, istat, ierr)

         !-- Store Lorentz-torques and rotation rates:
         if ( tscheme%family == 'MULTISTEP' ) then
            do n_o=2,tscheme%nexp
               call MPI_File_Write(fh, domega_ic_dt%expl(n_o), 1, MPI_DEF_REAL, &
                    &              istat, ierr)
            end do
            do n_o=2,tscheme%nimp
               call MPI_File_Write(fh, domega_ic_dt%impl(n_o), 1, MPI_DEF_REAL, &
                    &              istat, ierr)
            end do
            do n_o=2,tscheme%nold
               call MPI_File_Write(fh, domega_ic_dt%old(n_o), 1, MPI_DEF_REAL, &
                    &              istat, ierr)
            end do
            do n_o=2,tscheme%nexp
               call MPI_File_Write(fh, domega_ma_dt%expl(n_o), 1, MPI_DEF_REAL, &
                    &              istat, ierr)
            end do
            do n_o=2,tscheme%nimp
               call MPI_File_Write(fh, domega_ma_dt%impl(n_o), 1, MPI_DEF_REAL, &
                    &              istat, ierr)
            end do
            do n_o=2,tscheme%nold
               call MPI_File_Write(fh, domega_ma_dt%old(n_o), 1, MPI_DEF_REAL, &
                    &              istat, ierr)
            end do

            do n_o=2,tscheme%nexp
               call MPI_File_Write(fh, lorentz_torque_ic_dt%expl(n_o), 1, &
                    &              MPI_DEF_REAL, istat, ierr)
            end do
            do n_o=2,tscheme%nimp
               call MPI_File_Write(fh, lorentz_torque_ic_dt%impl(n_o), 1, &
                    &              MPI_DEF_REAL, istat, ierr)
            end do
            do n_o=2,tscheme%nold
               call MPI_File_Write(fh, lorentz_torque_ic_dt%old(n_o), 1, &
                    &              MPI_DEF_REAL, istat, ierr)
            end do
            do n_o=2,tscheme%nexp
               call MPI_File_Write(fh, lorentz_torque_ma_dt%expl(n_o), 1, &
                    &              MPI_DEF_REAL, istat, ierr)
            end do
            do n_o=2,tscheme%nimp
               call MPI_File_Write(fh, lorentz_torque_ma_dt%impl(n_o), 1, &
                    &              MPI_DEF_REAL, istat, ierr)
            end do
            do n_o=2,tscheme%nold
               call MPI_File_Write(fh, lorentz_torque_ma_dt%old(n_o), 1, &
                    &              MPI_DEF_REAL, istat, ierr)
            end do
         end if
         call MPI_File_Write(fh, omega_ic1, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, omegaOsz_ic1, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, tOmega_ic1, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, omega_ic2, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, omegaOsz_ic2, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, tOmega_ic2, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, omega_ma1, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, omegaOsz_ma1, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, tOmega_ma1, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, omega_ma2, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, omegaOsz_ma2, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, tOmega_ma2, 1, MPI_DEF_REAL, istat, ierr)

         !-- Write logical to know how many fields are stored
         call MPI_File_Write(fh, l_heat, 1, MPI_LOGICAL, istat, ierr)
         call MPI_File_Write(fh, l_chemical_conv, 1, MPI_LOGICAL, istat, ierr)
         call MPI_File_Write(fh, l_mag, 1, MPI_LOGICAL, istat, ierr)
         call MPI_File_Write(fh, l_press_store, 1, MPI_LOGICAL, istat, ierr)
         call MPI_File_Write(fh, l_cond_ic, 1, MPI_LOGICAL, istat, ierr)

         !-- Rank 0 gets the displacement
         call MPI_File_get_position(fh, offset, ierr)
         call MPI_File_get_byte_offset(fh, offset, disp, ierr)
      end if

      !-- Broadcast the displacement
      call MPI_Bcast(disp, 1, MPI_OFFSET, 0, MPI_COMM_WORLD, ierr)

      arr_size(1) = lm_max
      arr_size(2) = n_r_max
      arr_loc_size(1) = lm_max
      arr_loc_size(2) = n_r_loc
      arr_start(1) = 0
      arr_start(2) = nRstart-1
      call MPI_Type_Create_Subarray(2, arr_size, arr_loc_size, arr_start, &
           &                        MPI_ORDER_FORTRAN, MPI_DEF_COMPLEX,   &
           &                        datatype, ierr)
      call MPI_Type_Commit(datatype, ierr)

      !-- Set the view after the header
      call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
           &                 info, ierr)

      size_tmp=int(lm_max,kind=lip)*int(n_r_max,kind=lip)* &
      &        int(SIZEOF_DEF_COMPLEX,kind=lip)

      !--------------------
      !-- Now finally write the fields
      !--------------------

      !-- Poloidal potential: w
      call write_one_field_mpi(fh, info, datatype, tscheme, w, dwdt, &
           &                   work, lo2r, size_tmp,  disp)

      !-- Toroidal potential: z
      call write_one_field_mpi(fh, info, datatype, tscheme, z, dzdt, &
           &                   work, lo2r, size_tmp,  disp)

      !-- Pressure: p
      if ( l_press_store ) then
         call write_one_field_mpi(fh, info, datatype, tscheme, p, dpdt, &
              &                   work, lo2r, size_tmp,  disp)
      end if

      !-- Entropy: s
      if ( l_heat ) then
         call write_one_field_mpi(fh, info, datatype, tscheme, s, dsdt, &
              &                   work, lo2r, size_tmp,  disp)
      end if

      !-- Chemical composition: xi
      if ( l_chemical_conv ) then
         call write_one_field_mpi(fh, info, datatype, tscheme, xi, dxidt, &
              &                   work, lo2r, size_tmp,  disp)
      end if

      !-- Outer core magnetic field:
      if ( l_mag ) then
         call write_one_field_mpi(fh, info, datatype, tscheme, b, dbdt, &
              &                   work, lo2r, size_tmp,  disp)
         call write_one_field_mpi(fh, info, datatype, tscheme, aj, djdt, &
              &                   work, lo2r, size_tmp,  disp)
      end if

      !-- Displacement at the end of the file
      offset = 0
      call MPI_File_Seek(fh, offset, MPI_SEEK_END, ierr)
      call MPI_File_get_byte_offset(fh, offset, disp, ierr)
      call MPI_File_Set_View(fh, disp, MPI_BYTE, MPI_BYTE, "native", &
           &                 info, ierr)

      call MPI_Type_Free(datatype, ierr)
      deallocate( work )

      !-- Inner core magnetic field (only written by rank 0 for now)
      if ( l_mag .and. l_cond_ic ) then

         if ( rank == 0 ) then
            allocate ( work(lm_max, n_r_ic_max) )
         else
            allocate ( work(1,1) )
         end if

         call gather_all_from_lo_to_rank0(gt_IC,b_ic,work)
         if ( rank == 0 ) then
            call MPI_File_Write(fh, work, lm_max*n_r_ic_max, MPI_DEF_COMPLEX, &
                 &              istat, ierr)
         end if
         if ( tscheme%family == 'MULTISTEP' ) then
            do n_o=2,tscheme%nexp
               call gather_all_from_lo_to_rank0(gt_IC,dbdt_ic%expl(:,:,n_o),work)
               if ( rank == 0 ) then
                  call MPI_File_Write(fh, work, lm_max*n_r_ic_max, MPI_DEF_COMPLEX, &
                       &              istat, ierr)
               end if
            end do
            do n_o=2,tscheme%nimp
               call gather_all_from_lo_to_rank0(gt_IC,dbdt_ic%impl(:,:,n_o),work)
               if ( rank == 0 ) then
                  call MPI_File_Write(fh, work, lm_max*n_r_ic_max, MPI_DEF_COMPLEX, &
                       &              istat, ierr)
               end if
            end do
            do n_o=2,tscheme%nold
               call gather_all_from_lo_to_rank0(gt_IC,dbdt_ic%old(:,:,n_o),work)
               if ( rank == 0 ) then
                  call MPI_File_Write(fh, work, lm_max*n_r_ic_max, MPI_DEF_COMPLEX, &
                       &              istat, ierr)
               end if
            end do
         end if

         call gather_all_from_lo_to_rank0(gt_IC,aj_ic,work)
         if ( rank == 0 ) then
            call MPI_File_Write(fh, work, lm_max*n_r_ic_max, MPI_DEF_COMPLEX, &
                 &              istat, ierr)
         end if
         if ( tscheme%family == 'MULTISTEP' ) then
            do n_o=2,tscheme%nexp
               call gather_all_from_lo_to_rank0(gt_IC,djdt_ic%expl(:,:,n_o),work)
               if ( rank == 0 ) then
                  call MPI_File_Write(fh, work, lm_max*n_r_ic_max, MPI_DEF_COMPLEX, &
                       &              istat, ierr)
               end if
            end do
            do n_o=2,tscheme%nimp
               call gather_all_from_lo_to_rank0(gt_IC,djdt_ic%impl(:,:,n_o),work)
               if ( rank == 0 ) then
                  call MPI_File_Write(fh, work, lm_max*n_r_ic_max, MPI_DEF_COMPLEX, &
                       &              istat, ierr)
               end if
            end do
            do n_o=2,tscheme%nold
               call gather_all_from_lo_to_rank0(gt_IC,djdt_ic%old(:,:,n_o),work)
               if ( rank == 0 ) then
                  call MPI_File_Write(fh, work, lm_max*n_r_ic_max, MPI_DEF_COMPLEX, &
                       &              istat, ierr)
               end if
            end do
         end if

         deallocate( work ) 

      end if

      !-- Destroy the MPI communicator
      call lo2r%destroy_comm()

      !-- Close file 
      call MPI_Info_free(info, ierr)
      call MPI_File_close(fh, ierr)

      !-- Close checkpoint file and display a message in the log file
      if ( rank == 0 ) then

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

   end subroutine store_mpi
!-----------------------------------------------------------------------------------
   subroutine write_one_field_mpi(fh, info, datatype, tscheme, w, dwdt, &
              &                   work, lo2r, size_tmp,  disp)

      !-- Input variables
      integer,             intent(in) :: fh ! file unit
      integer,             intent(in) :: info ! MPI file info
      integer,             intent(in) :: datatype ! MPI file info
      class(type_tscheme), intent(in) :: tscheme
      integer(lip),        intent(in) :: size_tmp
      complex(cp),         intent(in) :: w(llm:ulm, n_r_max) ! field
      type(type_tarray),   intent(in) :: dwdt
      type(type_mpiatoav), intent(in) :: lo2r

      !-- Output variables
      integer(lip),        intent(inout) :: disp
      complex(cp),         intent(inout) :: work(llm:ulm, n_r_max) 

      !-- Local variables
      integer :: n_o
      integer :: istat(MPI_STATUS_SIZE)

      call MPI_File_Write_all(fh, w, lm_max*n_r_loc, &
           &                  MPI_DEF_COMPLEX, istat, ierr)
      disp = disp+size_tmp
      call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
           &                 info, ierr)

      if ( tscheme%family == 'MULTISTEP' ) then

         do n_o=2,tscheme%nexp
            call lo2r%transp_lm2r(dwdt%expl(:,:,n_o), work)
            call MPI_File_Write_all(fh, work, lm_max*n_r_loc, &
                 &                  MPI_DEF_COMPLEX, istat, ierr)
            disp = disp+size_tmp
            call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
                 &                 info, ierr)
         end do

         do n_o=2,tscheme%nimp
            call lo2r%transp_lm2r(dwdt%impl(:,:,n_o), work)
            call MPI_File_Write_all(fh, work, lm_max*n_r_loc, &
                 &                  MPI_DEF_COMPLEX, istat, ierr)
            disp = disp+size_tmp
            call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
                 &                 info, ierr)
         end do

         do n_o=2,tscheme%nold
            call lo2r%transp_lm2r(dwdt%old(:,:,n_o), work)
            call MPI_File_Write_all(fh, work, lm_max*n_r_loc, &
                 &                  MPI_DEF_COMPLEX, istat, ierr)
            disp = disp+size_tmp
            call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
                 &                 info, ierr)
         end do

      end if

   end subroutine write_one_field_mpi
!-----------------------------------------------------------------------------------
#endif
end module storeCheckPoints
