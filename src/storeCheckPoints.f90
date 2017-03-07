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
 
#ifdef WITH_HDF5
   public :: store, storeHdf5_serial, storeHdf5_parallel
#else
   public :: store
#endif

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
#ifdef WITH_HDF5
   subroutine storeHdf5_serial(time,dt,dtNew,w,z,p,s,b,aj,b_ic,aj_ic, &
                                 dwdtLast,dzdtLast,dpdtLast,dsdtLast, &
                                 dbdtLast,djdtLast,dbdt_icLast,djdt_icLast)

      use hdf5
      use hdf5Helpers, only: writeHdf5_attribute

      !--- Input variables
      real(cp),    intent(in) :: time,dt,dtNew
      complex(cp), intent(in) :: w(lm_max,n_r_max)
      complex(cp), intent(in) :: z(lm_max,n_r_max)
      complex(cp), intent(in) :: p(lm_max,n_r_max)
      complex(cp), intent(in) :: s(lm_max,n_r_max)
      complex(cp), intent(in) :: b(lm_maxMag,n_r_maxMag)
      complex(cp), intent(in) :: aj(lm_maxMag,n_r_maxMag)
      complex(cp), intent(in) :: b_ic(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: aj_ic(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: dwdtLast(lm_max,n_r_max)
      complex(cp), intent(in) :: dzdtLast(lm_max,n_r_max)
      complex(cp), intent(in) :: dpdtLast(lm_max,n_r_max)
      complex(cp), intent(in) :: dsdtLast(lm_max,n_r_max)
      complex(cp), intent(in) :: dbdtLast(lm_maxMag,n_r_maxMag)
      complex(cp), intent(in) :: djdtLast(lm_maxMag,n_r_maxMag)
      complex(cp), intent(in) :: dbdt_icLast(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: djdt_icLast(lm_maxMag,n_r_ic_maxMag)

      !--- HDF5 file identifier
      integer(HID_T) :: file_id  

      !--- HDF5 dataset
      integer(HID_T) :: dset1_id      ! Dataset identifier
      integer(HID_T) :: dspace1_id    ! Dataspace identifier
      integer(HSIZE_T) :: dims(2) ! Dataset dimensions

      !--- HDF5 Groups
      integer(HID_T) :: groupParams_id,groupFields_id   ! Groups identifiers
      integer(HID_T) :: groupDtFields_id,groupTorque_id 

      !--- HDF5 Attributes
      integer(HID_T) :: aspace_id,attr_id

      !--- HDF5 Type
      integer(HID_T) :: type_id
      integer(HSIZE_T) :: re_size,im_size,complex_t_size,offset
      integer     ::   error

      ! Initialize FORTRAN interface.
      call h5open_f(error)

      ! Create a new file using default properties.
      call h5fcreate_f(rst_file, H5F_ACC_TRUNC_F, file_id, error)

      ! Create a group for physical parameters
      call h5gcreate_f(file_id, '/Params', groupParams_id, error)

      ! Create/write/close attributes
      call h5screate_f(H5S_SCALAR_F,aspace_id,error)
      call writeHdf5_attribute(groupParams_id,aspace_id,'Time',time*tScale)
      call writeHdf5_attribute(groupParams_id,aspace_id,'dt',dt)
      call writeHdf5_attribute(groupParams_id,aspace_id,'Ra',ra)
      call writeHdf5_attribute(groupParams_id,aspace_id,'Pr',pr)
      call writeHdf5_attribute(groupParams_id,aspace_id,'Ek',ek)
      call writeHdf5_attribute(groupParams_id,aspace_id,'Radratio',radratio)
      call writeHdf5_attribute(groupParams_id,aspace_id,'Prmag',prmag)
      call writeHdf5_attribute(groupParams_id,aspace_id,'sigma_ratio',sigma_ratio)
      call h5sclose_f(aspace_id, error)

      ! Close group
      call h5gclose_f(groupParams_id, error)

      ! Create a group for torque and rotation rates
      call h5gcreate_f(file_id, '/Torque', groupTorque_id, error)
      call h5screate_f(H5S_SCALAR_F,aspace_id,error)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'lorentz_torque_ic',  &
                              lorentz_torque_icLast)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'lorentz_torque_ma',  &
                              lorentz_torque_maLast)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'omega_ic1',omega_ic1)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'omegaOsz_ic1',omegaOsz_ic1)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'tOmega_ic1',tOmega_ic1)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'omega_ic2',omega_ic2)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'omegaOsz_ic2',omegaOsz_ic2)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'tOmega_ic2',tOmega_ic2)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'omega_ma1',omega_ma1)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'omegaOsz_ma1',omegaOsz_ma1)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'tOmega_ma1',tOmega_ma1)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'omega_ma2',omega_ma2)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'omegaOsz_ma2',omegaOsz_ma2)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'tOmega_ma2',tOmega_ma2)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'dtNew',dtNew)
      call h5sclose_f(aspace_id, error)

      ! Close group
      call h5gclose_f(groupTorque_id, error)

      ! Create a group for the fields
      call h5gcreate_f(file_id, '/Fields', groupFields_id, error)

      ! Create a group for the time-derivative of the fields
      call h5gcreate_f(file_id, '/dtFields', groupDtFields_id, error)

      ! Create/write/close attributes
      call h5screate_f(H5S_SCALAR_F,aspace_id,error)
      call writeHdf5_attribute(groupFields_id,aspace_id,'n_r_max',n_r_max)
      call writeHdf5_attribute(groupFields_id,aspace_id,'n_theta_max',n_theta_max)
      call writeHdf5_attribute(groupFields_id,aspace_id,'n_phi_tot',n_phi_tot)
      call writeHdf5_attribute(groupFields_id,aspace_id,'minc',minc)
      call writeHdf5_attribute(groupFields_id,aspace_id,'nalias',nalias)
      call writeHdf5_attribute(groupFields_id,aspace_id,'n_r_ic_max',n_r_ic_max)
      call h5sclose_f(aspace_id, error)

      ! Create a HF compound type  to store Fortran complex
      call h5tget_size_f(H5T_NATIVE_DOUBLE,re_size,error)
      call h5tget_size_f(H5T_NATIVE_DOUBLE,im_size,error)
      complex_t_size = re_size+im_size
      call h5tcreate_f(H5T_COMPOUND_F,complex_t_size,type_id,error)
      offset = 0
      call h5tinsert_f(type_id, 'real', offset, H5T_NATIVE_DOUBLE, error)
      offset = offset + re_size
      call h5tinsert_f(type_id, 'imag', offset, H5T_NATIVE_DOUBLE, error)

      ! Create the dataspace for physical arrays
      dims = [lm_max, n_r_max]
      call h5screate_simple_f(2, dims, dspace1_id, error)

      ! Create and write the dataset with default properties.
      call h5dcreate_f(groupFields_id,'w_pol',type_id, dspace1_id,dset1_id, error)
      call h5dwrite_f(dset1_id, type_id, C_LOC(w), error)
      call h5dclose_f(dset1_id, error)

      call h5dcreate_f(groupDtFields_id,'dwdtLast',type_id,dspace1_id,dset1_id, error)
      call h5dwrite_f(dset1_id, type_id, C_LOC(dwdtLast), error)
      call h5dclose_f(dset1_id, error)

      call h5dcreate_f(groupFields_id,'z_tor',type_id,dspace1_id,dset1_id,error)
      call h5dwrite_f(dset1_id, type_id, C_LOC(z), error)
      call h5dclose_f(dset1_id, error)

      call h5dcreate_f(groupDtFields_id,'dzdtLast',type_id,dspace1_id,dset1_id,error)
      call h5dwrite_f(dset1_id, type_id, C_LOC(dzdtLast), error)
      call h5dclose_f(dset1_id, error)

      call h5dcreate_f(groupFields_id,'pressure',type_id,dspace1_id,dset1_id,error)
      call h5dwrite_f(dset1_id, type_id, C_LOC(p), error)
      call h5dclose_f(dset1_id, error)

      call h5dcreate_f(groupDtFields_id,'dpdtLast',type_id,dspace1_id,dset1_id,error)
      call h5dwrite_f(dset1_id, type_id, C_LOC(dpdtLast), error)
      call h5dclose_f(dset1_id, error)


      if ( l_heat ) then
         call h5dcreate_f(groupFields_id,'entropy',type_id,dspace1_id,dset1_id,error)
         call h5dwrite_f(dset1_id, type_id, C_LOC(s), error)
         call h5dclose_f(dset1_id, error)

         call h5dcreate_f(groupDtFields_id,'dsdtLast',type_id,dspace1_id,dset1_id,error)
         call h5dwrite_f(dset1_id, type_id, C_LOC(dsdtLast), error)
         call h5dclose_f(dset1_id, error)
      end if

      ! Terminate access to the data space.
      call h5sclose_f(dspace1_id, error)

      if ( l_mag ) then
         dims = [lm_maxMag, n_r_maxMag]
         call h5screate_simple_f(2, dims, dspace1_id, error)

         call h5dcreate_f(groupFields_id,'b_pol',type_id,dspace1_id,dset1_id,error)
         call h5dwrite_f(dset1_id, type_id, C_LOC(b), error)
         call h5dclose_f(dset1_id, error)

         call h5dcreate_f(groupDtFields_id,'dbdtLast',type_id,dspace1_id,dset1_id, error)
         call h5dwrite_f(dset1_id, type_id, C_LOC(dbdtLast), error)
         call h5dclose_f(dset1_id, error)

         call h5dcreate_f(groupFields_id,'aj_tor',type_id,dspace1_id,dset1_id,error)
         call h5dwrite_f(dset1_id, type_id, C_LOC(aj), error)
         call h5dclose_f(dset1_id, error)

         call h5dcreate_f(groupDtFields_id,'djdtLast',type_id,dspace1_id,dset1_id,error)
         call h5dwrite_f(dset1_id, type_id, C_LOC(djdtLast), error)
         call h5dclose_f(dset1_id, error)

         call h5sclose_f(dspace1_id, error)

         if ( l_cond_ic ) then
            dims = [lm_maxMag, n_r_ic_maxMag]
            call h5screate_simple_f(2, dims, dspace1_id, error)

            call h5dcreate_f(groupFields_id,'b_ic_pol',type_id,dspace1_id,dset1_id,error)
            call h5dwrite_f(dset1_id, type_id, C_LOC(b_ic), error)
            call h5dclose_f(dset1_id, error)

            call h5dcreate_f(groupDtFields_id,'dbdt_icLast',type_id,dspace1_id, &
                             dset1_id,error)
            call h5dwrite_f(dset1_id, type_id, C_LOC(dbdt_icLast), error)
            call h5dclose_f(dset1_id, error)

            call h5dcreate_f(groupFields_id,'aj_ic_tor',type_id,dspace1_id, &
                             dset1_id,error)
            call h5dwrite_f(dset1_id, type_id, C_LOC(aj_ic), error)
            call h5dclose_f(dset1_id, error)

            call h5dcreate_f(groupDtFields_id,'djdt_icLast',type_id,dspace1_id, &
                             dset1_id,error)
            call h5dwrite_f(dset1_id, type_id, C_LOC(djdt_icLast), error)
            call h5dclose_f(dset1_id, error)

            call h5sclose_f(dspace1_id, error)
         end if
      end if


      ! End access to the dataset and release resources used by it.


      ! close groups
      call h5gclose_f(groupFields_id, error)
      call h5gclose_f(groupDtFields_id, error)

      ! Close the file.
      call h5fclose_f(file_id, error)

      ! Close FORTRAN interface.
      call h5close_f(error)

   end subroutine storeHdf5_serial
!----------------------------------------------------------------------
   subroutine storeHdf5_parallel(time,dt,dtNew,w,z,p,s,b,aj,b_ic,aj_ic, &
                                 dwdtLast,dzdtLast,dpdtLast,dsdtLast,   &
                                 dbdtLast,djdtLast,dbdt_icLast,djdt_icLast)

      use parallel_mod
      use hdf5
      use hdf5Helpers, only: writeHdf5_attribute, write_dataset
      use LMLoop_data, only: llm, ulm, llmMag, ulmMag

      !--- Input variables
      real(cp),    intent(in) :: time,dt,dtNew
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(in) :: z(llm:ulm,n_r_max)
      complex(cp), intent(in) :: p(llm:ulm,n_r_max)
      complex(cp), intent(in) :: s(llm:ulm,n_r_max)
      complex(cp), intent(in) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: aj(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: dwdtLast(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dzdtLast(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dpdtLast(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dsdtLast(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dbdtLast(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: djdtLast(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: dbdt_icLast(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: djdt_icLast(llmMag:ulmMag,n_r_ic_maxMag)

      integer :: l,m,lm
      !--- HDF5 file identifier
      integer(HID_T) :: file_id  
      integer(HID_T) :: plist_id

      !--- HDF5 dataset
      integer(HSIZE_T) :: dims_full(2) ! Dataset dimensions

      !--- HDF5 Groups
      integer(HID_T) :: groupParams_id,groupFields_id   ! Groups identifiers
      integer(HID_T) :: groupDtFields_id,groupTorque_id 

      !--- HDF5 Attributes
      integer(HID_T) :: aspace_id,attr_id

      !--- HDF5 Type
      integer(HID_T) :: type_id
      integer(HSIZE_T) :: re_size,im_size,complex_t_size,offset
      integer :: error

      ! Attributes must be written collectively in HDF5: broadcasting
      ! some scalars is required
      call MPI_Bcast(omega_ic1,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(omega_ic2,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(omegaOsz_ic1,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(omegaOsz_ic2,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(tOmega_ic1,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(tOmega_ic2,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(omega_ma1,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(omega_ma2,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(omegaOsz_ma1,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(omegaOsz_ma2,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(tOmega_ma1,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(tOmega_ma2,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)

      ! Initialize FORTRAN interface.
      call h5open_f(error)

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

      ! Create a new file collectively
      call h5fcreate_f(rst_file, H5F_ACC_TRUNC_F, file_id, error, access_prp=plist_id)
      call h5pclose_f(plist_id, error)

      ! Create a group for physical parameters
      call h5gcreate_f(file_id, '/Params', groupParams_id, error)
 
      ! Create/write/close attributes
      call h5screate_f(H5S_SCALAR_F,aspace_id,error)
      call writeHdf5_attribute(groupParams_id,aspace_id,'Time',time*tScale)
      call writeHdf5_attribute(groupParams_id,aspace_id,'dt',dt)
      call writeHdf5_attribute(groupParams_id,aspace_id,'Ra',ra)
      call writeHdf5_attribute(groupParams_id,aspace_id,'Pr',pr)
      call writeHdf5_attribute(groupParams_id,aspace_id,'Ek',ek)
      call writeHdf5_attribute(groupParams_id,aspace_id,'Radratio',radratio)
      call writeHdf5_attribute(groupParams_id,aspace_id,'Prmag',prmag)
      call writeHdf5_attribute(groupParams_id,aspace_id,'sigma_ratio',sigma_ratio)
      call h5sclose_f(aspace_id, error)

      ! Close group
      call h5gclose_f(groupParams_id, error)

      ! Create a group for torque and rotation rates
      call h5gcreate_f(file_id, '/Torque', groupTorque_id, error)
      call h5screate_f(H5S_SCALAR_F,aspace_id,error)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'lorentz_torque_ic',  &
                              lorentz_torque_icLast)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'lorentz_torque_ma',  &
                              lorentz_torque_maLast)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'omega_ic1',omega_ic1)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'omegaOsz_ic1',omegaOsz_ic1)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'tOmega_ic1',tOmega_ic1)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'omega_ic2',omega_ic2)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'omegaOsz_ic2',omegaOsz_ic2)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'tOmega_ic2',tOmega_ic2)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'omega_ma1',omega_ma1)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'omegaOsz_ma1',omegaOsz_ma1)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'tOmega_ma1',tOmega_ma1)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'omega_ma2',omega_ma2)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'omegaOsz_ma2',omegaOsz_ma2)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'tOmega_ma2',tOmega_ma2)
      call writeHdf5_attribute(groupTorque_id,aspace_id,'dtNew',dtNew)

      call h5sclose_f(aspace_id, error)

      ! Close group
      call h5gclose_f(groupTorque_id, error)

      ! Create a group for the fields
      call h5gcreate_f(file_id, '/Fields', groupFields_id, error)

      ! Create a group for the time-derivative of the fields
      call h5gcreate_f(file_id, '/dtFields', groupDtFields_id, error)

      ! Create/write/close attributes
      call h5screate_f(H5S_SCALAR_F,aspace_id,error)
      call writeHdf5_attribute(groupFields_id,aspace_id,'n_r_max',n_r_max)
      call writeHdf5_attribute(groupFields_id,aspace_id,'n_theta_max',n_theta_max)
      call writeHdf5_attribute(groupFields_id,aspace_id,'n_phi_tot',n_phi_tot)
      call writeHdf5_attribute(groupFields_id,aspace_id,'minc',minc)
      call writeHdf5_attribute(groupFields_id,aspace_id,'nalias',nalias)
      call writeHdf5_attribute(groupFields_id,aspace_id,'n_r_ic_max',n_r_ic_max)
      call h5sclose_f(aspace_id, error)

      ! Create a HF compound type  to store Fortran complex
      call h5tget_size_f(H5T_NATIVE_DOUBLE,re_size,error)
      call h5tget_size_f(H5T_NATIVE_DOUBLE,im_size,error)
      complex_t_size = re_size+im_size
      call h5tcreate_f(H5T_COMPOUND_F,complex_t_size,type_id,error)
      offset = 0
      call h5tinsert_f(type_id, 'real', offset, H5T_NATIVE_DOUBLE, error)
      offset = offset + re_size
      call h5tinsert_f(type_id, 'imag', offset, H5T_NATIVE_DOUBLE, error)

      ! Create the dataspace for physical arrays
      dims_full = [lm_max, n_r_max]

      call write_dataset(groupFields_id, 'w_pol', type_id, w, n_r_max, dims_full)
      call write_dataset(groupFields_id, 'z_tor', type_id, z, n_r_max, dims_full)
      call write_dataset(groupFields_id, 'pressure', type_id, p, n_r_max, dims_full)

      call write_dataset(groupDtFields_id, 'dwdtLast', type_id, dwdtLast, n_r_max, &
                         dims_full)
      call write_dataset(groupDtFields_id, 'dzdtLast', type_id, dzdtLast, n_r_max, &
                         dims_full)
      call write_dataset(groupDtFields_id, 'dpdtLast', type_id, dpdtLast, n_r_max, &
                         dims_full)

      if ( l_heat ) then
         call write_dataset(groupFields_id, 'entropy', type_id, s, n_r_max, dims_full)
         call write_dataset(groupDtFields_id, 'dsdtLast', type_id, dsdtLast, &
                            n_r_max, dims_full)
      end if

      if ( l_mag ) then
         dims_full = [lm_maxMag, n_r_maxMag]
         call write_dataset(groupFields_id, 'b_pol', type_id, b, n_r_maxMag, dims_full)
         call write_dataset(groupFields_id, 'aj_tor', type_id, aj, n_r_maxMag, dims_full)

         call write_dataset(groupDtFields_id, 'dbdtLast', type_id, dbdtLast, &
                            n_r_maxMag, dims_full)
         call write_dataset(groupDtFields_id, 'djdtLast', type_id, djdtLast, &
                            n_r_maxMag, dims_full)

         if ( l_cond_ic ) then
            dims_full = [lm_maxMag, n_r_ic_maxMag]
            call write_dataset(groupFields_id, 'b_ic_pol', type_id, b_ic, &
                               n_r_ic_maxMag, dims_full)
            call write_dataset(groupFields_id, 'aj_ic_tor', type_id, aj_ic, & 
                               n_r_ic_maxMag, dims_full)

            call write_dataset(groupDtFields_id, 'dbdt_icLast', type_id, dbdt_icLast, &
                               n_r_ic_maxMag, dims_full)
            call write_dataset(groupDtFields_id, 'djdt_icLast', type_id, djdt_icLast, &
                               n_r_ic_maxMag, dims_full)
         end if
      end if

      ! close groups
      call h5gclose_f(groupFields_id, error)
      call h5gclose_f(groupDtFields_id, error)

      ! Close the file.
      call h5fclose_f(file_id, error)

      ! Close FORTRAN interface.
      call h5close_f(error)

   end subroutine storeHdf5_parallel
#endif
end module storeCheckPoints
