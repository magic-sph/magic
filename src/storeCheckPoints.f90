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
       &                 fd_stretch, fd_ratio
   use radial_functions, only: rscheme_oc
   use physical_parameters, only: ra, pr, prmag, radratio, ek, sigma_ratio, &
       &                          raxi, sc
   use blocking, only: llm, ulm, llmMag, ulmMag
   use radial_data, only: nRstart, nRstop, nRstartMag, nRstopMag
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

#ifdef WITH_MPI
   public :: store, store_mpi
#else
   public :: store
#endif

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
#ifdef WITH_MPI
   subroutine store_mpi(time,dt,dtNew,n_time_step,l_stop_time,l_new_rst_file, &
              &         l_ave_file,w,z,p,s,xi,b,aj,b_ic,aj_ic,dwdtLast,       &
              &         dzdtLast,dpdtLast,dsdtLast,dxidtLast,dbdtLast,        &
              &         djdtLast,dbdt_icLast,djdt_icLast)
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
      complex(cp), intent(in) :: w(lm_max,nRstart:nRstop)
      complex(cp), intent(in) :: z(lm_max,nRstart:nRstop)
      complex(cp), intent(in) :: p(lm_max,nRstart:nRstop)
      complex(cp), intent(in) :: s(lm_max,nRstart:nRstop)
      complex(cp), intent(in) :: xi(lm_max,nRstart:nRstop)
      complex(cp), intent(in) :: b(lm_maxMag,nRstartMag:nRstopMag)
      complex(cp), intent(in) :: aj(lm_maxMag,nRstartMag:nRstopMag)
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

      type(type_mpiatoav) :: lo2r
      integer :: version, info, fh, datatype
      character(len=72) :: string, rst_file
      integer :: istat(MPI_STATUS_SIZE)
      integer :: arr_size(2), arr_loc_size(2), arr_start(2)
      integer(kind=MPI_OFFSET_KIND) :: disp, offset

      version = 1

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
         call MPI_File_Write(fh, dt*tScale, 1, MPI_DEF_REAL, istat, ierr)
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

         !-- Store Lorentz-torques and rotation rates:
         call MPI_File_Write(fh, lorentz_torque_icLast, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, lorentz_torque_maLast, 1, MPI_DEF_REAL, istat, ierr)
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
         call MPI_File_Write(fh, dtNew, 1, MPI_DEF_REAL, istat, ierr)

         !-- Write logical to know how many fields are stored
         call MPI_File_Write(fh, l_heat, 1, MPI_LOGICAL, istat, ierr)
         call MPI_File_Write(fh, l_chemical_conv, 1, MPI_LOGICAL, istat, ierr)
         call MPI_File_Write(fh, l_mag, 1, MPI_LOGICAL, istat, ierr)
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
      arr_loc_size(2) = nR_per_rank
      arr_start(1) = 0
      arr_start(2) = nRstart-1
      call MPI_Type_Create_Subarray(2, arr_size, arr_loc_size, arr_start, &
           &                        MPI_ORDER_FORTRAN, MPI_DEF_COMPLEX,   &
           &                        datatype, ierr)
      call MPI_Type_Commit(datatype, ierr)

      !-- Set the view after the header
      call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
           &                 info, ierr)

      !-- Now finally write the fields
      !-- Poloidal potential: w
      call MPI_File_Write_all(fh, w, lm_max*nR_per_rank, &
           &                  MPI_DEF_COMPLEX, istat, ierr)
      disp = disp+lm_max*n_r_max*SIZEOF_DEF_COMPLEX
      call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
           &                 info, ierr)
      call lo2r%transp_lm2r(dwdtLast, work)
      call MPI_File_Write_all(fh, work, lm_max*nR_per_rank, &
           &                  MPI_DEF_COMPLEX, istat, ierr)
      disp = disp+lm_max*n_r_max*SIZEOF_DEF_COMPLEX
      call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
           &                 info, ierr)

      !-- Toroidal potential: z
      call MPI_File_Write_all(fh, z, lm_max*nR_per_rank, &
           &                  MPI_DEF_COMPLEX, istat, ierr)
      disp = disp+lm_max*n_r_max*SIZEOF_DEF_COMPLEX
      call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
           &                 info, ierr)
      call lo2r%transp_lm2r(dzdtLast, work)
      call MPI_File_Write_all(fh, work, lm_max*nR_per_rank, &
           &                  MPI_DEF_COMPLEX, istat, ierr)
      disp = disp+lm_max*n_r_max*SIZEOF_DEF_COMPLEX
      call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
           &                 info, ierr)

      !-- Pressure: p
      call MPI_File_Write_all(fh, p, lm_max*nR_per_rank, &
           &                  MPI_DEF_COMPLEX, istat, ierr)
      disp = disp+lm_max*n_r_max*SIZEOF_DEF_COMPLEX
      call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
           &                 info, ierr)
      call lo2r%transp_lm2r(dpdtLast, work)
      call MPI_File_Write_all(fh, work, lm_max*nR_per_rank, &
           &                  MPI_DEF_COMPLEX, istat, ierr)
      disp = disp+lm_max*n_r_max*SIZEOF_DEF_COMPLEX
      call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
           &                 info, ierr)

      !-- Entropy: s
      if ( l_heat ) then
         call MPI_File_Write_all(fh, s, lm_max*nR_per_rank, &
              &                  MPI_DEF_COMPLEX, istat, ierr)
         disp = disp+lm_max*n_r_max*SIZEOF_DEF_COMPLEX
         call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
              &                 info, ierr)
         call lo2r%transp_lm2r(dsdtLast, work)
         call MPI_File_Write_all(fh, work, lm_max*nR_per_rank, &
              &                  MPI_DEF_COMPLEX, istat, ierr)
         disp = disp+lm_max*n_r_max*SIZEOF_DEF_COMPLEX
         call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
              &                 info, ierr)
      end if

      !-- Chemical composition: xi
      if ( l_chemical_conv ) then
         call MPI_File_Write_all(fh, xi, lm_max*nR_per_rank, &
              &                  MPI_DEF_COMPLEX, istat, ierr)
         disp = disp+lm_max*n_r_max*SIZEOF_DEF_COMPLEX
         call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
              &                 info, ierr)
         call lo2r%transp_lm2r(dxidtLast, work)
         call MPI_File_Write_all(fh, work, lm_max*nR_per_rank, &
              &                  MPI_DEF_COMPLEX, istat, ierr)
         disp = disp+lm_max*n_r_max*SIZEOF_DEF_COMPLEX
         call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
              &                 info, ierr)
      end if

      !-- Outer core magnetic field:
      if ( l_mag ) then
         call MPI_File_Write_all(fh, b, lm_max*nR_per_rank, &
              &                  MPI_DEF_COMPLEX, istat, ierr)
         disp = disp+lm_max*n_r_max*SIZEOF_DEF_COMPLEX
         call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
              &                 info, ierr)
         call lo2r%transp_lm2r(dbdtLast, work)
         call MPI_File_Write_all(fh, work, lm_max*nR_per_rank, &
              &                  MPI_DEF_COMPLEX, istat, ierr)
         disp = disp+lm_max*n_r_max*SIZEOF_DEF_COMPLEX
         call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
              &                 info, ierr)
         call MPI_File_Write_all(fh, aj, lm_max*nR_per_rank, &
              &                  MPI_DEF_COMPLEX, istat, ierr)
         disp = disp+lm_max*n_r_max*SIZEOF_DEF_COMPLEX
         call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
              &                 info, ierr)
         call lo2r%transp_lm2r(djdtLast, work)
         call MPI_File_Write_all(fh, work, lm_max*nR_per_rank, &
              &                  MPI_DEF_COMPLEX, istat, ierr)
         disp = disp+lm_max*n_r_max*SIZEOF_DEF_COMPLEX
         call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, datatype, "native", &
              &                 info, ierr)
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
         call gather_all_from_lo_to_rank0(gt_IC,dbdt_icLast,work)
         if ( rank == 0 ) then
            call MPI_File_Write(fh, work, lm_max*n_r_ic_max, MPI_DEF_COMPLEX, &
                 &              istat, ierr)
         end if
         call gather_all_from_lo_to_rank0(gt_IC,aj_ic,work)
         if ( rank == 0 ) then
            call MPI_File_Write(fh, work, lm_max*n_r_ic_max, MPI_DEF_COMPLEX, istat, &
                 &              ierr)
         end if
         call gather_all_from_lo_to_rank0(gt_IC,djdt_icLast,work)
         if ( rank == 0 ) then
            call MPI_File_Write(fh, work, lm_max*n_r_ic_max, MPI_DEF_COMPLEX, istat, &
                 &              ierr)
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
#endif
!----------------------------------------------------------------------
end module storeCheckPoints
