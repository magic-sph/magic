module storeCheckPoints
   !
   ! This module contains several subroutines that can be used to store the
   ! rst_#.TAG files
   !

   use precision_mod
   use truncation, only: n_r_max,n_r_ic_max,minc,nalias,n_theta_max,n_phi_tot, &
                         lm_max,lm_maxMag,n_r_maxMag,n_r_ic_maxMag,l_max
   use physical_parameters, only: ra,pr,prmag,radratio,ek,sigma_ratio
   use num_param, only: tScale
   use fieldsLast, only: d_omega_ma_dtLast,d_omega_ic_dtLast, &
                         lorentz_torque_maLast,lorentz_torque_icLast
   use init_fields, only: inform,omega_ic1,omegaOsz_ic1,tOmega_ic1, &
                          omega_ic2,omegaOsz_ic2,tOmega_ic2,        &
                          omega_ma1,omegaOsz_ma1,tOmega_ma1,        &
                          omega_ma2,omegaOsz_ma2,tOmega_ma2
   use logic, only: l_heat,l_mag,l_cond_ic
   use output_data, only: n_rst_file,rst_file

   implicit none

   private
 
#ifdef WITH_HDF5
   public :: store, storeHdf5_serial, storeHdf5_parallel
#else
   public :: store
#endif

contains

   subroutine store(time,dt,dtNew,w,z,p,s,b,aj,b_ic,aj_ic, &
                    dwdtLast,dzdtLast,dpdtLast,dsdtLast,   &
                    dbdtLast,djdtLast,dbdt_icLast,djdt_icLast)
      !
      ! store results on disc file (restart file)
      ! In addition to the magnetic field and velocity potentials
      ! we store the time derivative terms
      ! djdt(lm,nR),dbdt(lm,nR), ......
      !

      !-- Input of variables:
      real(cp),    intent(in) :: time,dt,dtNew

      !-- Input of scalar fields to be stored:
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

      !-- Write parameters:
      if ( .not. l_heat ) then
         inform=11
      else
         inform=12
      end if

      write(n_rst_file) time*tScale,dt*tScale,ra,pr,prmag,ek,radratio, &
                     inform,n_r_max,n_theta_max,n_phi_tot,minc,nalias, &
                                               n_r_ic_max,sigma_ratio

      if ( l_heat ) then
         write(n_rst_file) w,z,p,s
         write(n_rst_file) dsdtLast,dwdtLast,dzdtLast,dpdtLast
      else
         write(n_rst_file) w,z,p
         write(n_rst_file) dwdtLast,dzdtLast,dpdtLast
      end if

      !-- Write magnetic field:
      if ( l_mag ) write(n_rst_file) b,aj,dbdtLast,djdtLast

      !-- Write IC magnetic field:
      if ( l_mag .and. l_cond_ic ) &
          write(n_rst_file) b_ic,aj_ic,dbdt_icLast,djdt_icLast

      !-- Store Lorentz-torques and rotation rates:
      write(n_rst_file) lorentz_torque_icLast, &
                        lorentz_torque_maLast, &
            omega_ic1,omegaOsz_ic1,tOmega_ic1, &
            omega_ic2,omegaOsz_ic2,tOmega_ic2, &
            omega_ma1,omegaOsz_ma1,tOmega_ma1, &
            omega_ma2,omegaOsz_ma2,tOmega_ma2, &
                                        dtNew

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
