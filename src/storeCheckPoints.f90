#include "assertions.cpp"
#include "perflib_preproc.cpp"
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
   use assertions

   implicit none

   private

   public :: store
#ifdef WITH_HDF5
   public :: storeHdf5
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
   subroutine storeHdf5(time,dt,dtNew,w,z,p,s,b,aj,b_ic,aj_ic, &
                        dwdtLast,dzdtLast,dpdtLast,dsdtLast,   &
                        dbdtLast,djdtLast,dbdt_icLast,djdt_icLast)

      use parallel_mod
      use hdf5
      use hdf5Helpers
      use blocking, only: lo_map
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

      !--- HDF5 file identifier
      integer(HID_T) :: file_id
      integer(HID_T) :: plist_id

      !--- HDF5 Groups
      integer(HID_T) :: groupParams_id,groupFields_id   ! Groups identifiers
      integer(HID_T) :: groupDtFields_id,groupTorque_id

      !--- HDF5 Type
      integer(HID_T) :: complex_t
      integer :: error

      ! Attributes must be written collectively in HDF5: broadcasting
      ! some scalars is required
#ifdef WITH_MPI
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
#endif

      PERFON('HDF5-IO')

      ! Initialize FORTRAN interface.
      call h5open_f(error)
      assert_equal(error, 0, "h5open_f")

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      assert_equal(error, 0, "h5pcreate_f")

#ifdef WITH_MPI
      call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
      assert_equal(error, 0, "h5pset_fapl_mpio_f")
#endif

      ! Create a new file collectively
      call h5fcreate_f(rst_file, H5F_ACC_TRUNC_F, file_id, error, access_prp=plist_id)
      assert_equal(error, 0, "h5fcreate_f")

      call h5pclose_f(plist_id, error)
      assert_equal(error, 0, "h5pclose_f")

      ! Create a group for physical parameters
      call h5gcreate_f(file_id, '/Params', groupParams_id, error)
      assert_equal(error, 0, "h5gcreate_f")
 
      ! Create/write/close attributes
      call writeHdf5_attribute(groupParams_id,'Time',time*tScale)
      call writeHdf5_attribute(groupParams_id,'dt',dt)
      call writeHdf5_attribute(groupParams_id,'Ra',ra)
      call writeHdf5_attribute(groupParams_id,'Pr',pr)
      call writeHdf5_attribute(groupParams_id,'Ek',ek)
      call writeHdf5_attribute(groupParams_id,'Radratio',radratio)
      call writeHdf5_attribute(groupParams_id,'Prmag',prmag)
      call writeHdf5_attribute(groupParams_id,'sigma_ratio',sigma_ratio)

      ! Close group
      call h5gclose_f(groupParams_id, error)
      assert_equal(error, 0, "h5gclose_f")

      ! Create a group for torque and rotation rates
      call h5gcreate_f(file_id, '/Torque', groupTorque_id, error)
      assert_equal(error, 0, "h5gcreate_f")
      call writeHdf5_attribute(groupTorque_id,'lorentz_torque_ic',  &
                              lorentz_torque_icLast)
      call writeHdf5_attribute(groupTorque_id,'lorentz_torque_ma',  &
                              lorentz_torque_maLast)
      call writeHdf5_attribute(groupTorque_id,'omega_ic1',omega_ic1)
      call writeHdf5_attribute(groupTorque_id,'omegaOsz_ic1',omegaOsz_ic1)
      call writeHdf5_attribute(groupTorque_id,'tOmega_ic1',tOmega_ic1)
      call writeHdf5_attribute(groupTorque_id,'omega_ic2',omega_ic2)
      call writeHdf5_attribute(groupTorque_id,'omegaOsz_ic2',omegaOsz_ic2)
      call writeHdf5_attribute(groupTorque_id,'tOmega_ic2',tOmega_ic2)
      call writeHdf5_attribute(groupTorque_id,'omega_ma1',omega_ma1)
      call writeHdf5_attribute(groupTorque_id,'omegaOsz_ma1',omegaOsz_ma1)
      call writeHdf5_attribute(groupTorque_id,'tOmega_ma1',tOmega_ma1)
      call writeHdf5_attribute(groupTorque_id,'omega_ma2',omega_ma2)
      call writeHdf5_attribute(groupTorque_id,'omegaOsz_ma2',omegaOsz_ma2)
      call writeHdf5_attribute(groupTorque_id,'tOmega_ma2',tOmega_ma2)
      call writeHdf5_attribute(groupTorque_id,'dtNew',dtNew)


      ! Close group
      call h5gclose_f(groupTorque_id, error)
      assert_equal(error, 0, "h5gclose_f")

      ! Create a group for the fields
      call h5gcreate_f(file_id, '/Fields', groupFields_id, error)
      assert_equal(error, 0, "h5gcreate_f")

      ! Create a group for the time-derivative of the fields
      call h5gcreate_f(file_id, '/dtFields', groupDtFields_id, error)
      assert_equal(error, 0, "h5gcreate_f")

      ! Create/write/close attributes
      call writeHdf5_attribute(groupFields_id,'n_r_max',n_r_max)
      call writeHdf5_attribute(groupFields_id,'n_theta_max',n_theta_max)
      call writeHdf5_attribute(groupFields_id,'n_phi_tot',n_phi_tot)
      call writeHdf5_attribute(groupFields_id,'minc',minc)
      call writeHdf5_attribute(groupFields_id,'nalias',nalias)
      call writeHdf5_attribute(groupFields_id,'n_r_ic_max',n_r_ic_max)

      ! Create a HDF5 compound type  to store Fortran complex variables
      complex_t = hdf5_fortran_complex_type()

      ! save this datatype in the file
      call h5tcommit_f(file_id, "complex_t", complex_t, error)
      assert_equal(error, 0, "h5tcommit_f")

      ! Save the current LM-mapping
      call write_dataset(file_id, 'lm_map', H5T_NATIVE_INTEGER, lo_map%lm2)

      call write_dataset(groupFields_id, 'w_pol',    complex_t, w, [llm, 1], [lm_max, n_r_max])
      call write_dataset(groupFields_id, 'z_tor',    complex_t, z, [llm, 1], [lm_max, n_r_max])
      call write_dataset(groupFields_id, 'pressure', complex_t, p, [llm, 1], [lm_max, n_r_max])

      call write_dataset(groupDtFields_id, 'dwdtLast', complex_t, dwdtLast, [llm, 1], [lm_max, n_r_max])
      call write_dataset(groupDtFields_id, 'dzdtLast', complex_t, dzdtLast, [llm, 1], [lm_max, n_r_max])
      call write_dataset(groupDtFields_id, 'dpdtLast', complex_t, dpdtLast, [llm, 1], [lm_max, n_r_max])

      if ( l_heat ) then
         call write_dataset(groupFields_id,   'entropy',  complex_t, s,        [llm, 1], [lm_max, n_r_max])
         call write_dataset(groupDtFields_id, 'dsdtLast', complex_t, dsdtLast, [llm, 1], [lm_max, n_r_max])
      end if

      if ( l_mag ) then
         call write_dataset(groupFields_id, 'b_pol',  complex_t, b,  [llmMag, 1], [lm_maxMag, n_r_maxMag])
         call write_dataset(groupFields_id, 'aj_tor', complex_t, aj, [llmMag, 1], [lm_maxMag, n_r_maxMag])
         call write_dataset(groupDtFields_id, 'dbdtLast', complex_t, dbdtLast, [llmMag, 1], [lm_maxMag, n_r_maxMag])
         call write_dataset(groupDtFields_id, 'djdtLast', complex_t, djdtLast, [llmMag, 1], [lm_maxMag, n_r_maxMag])

         if ( l_cond_ic ) then
            call write_dataset(groupFields_id,   'aj_ic_tor',   complex_t, aj_ic,       [llmMag, 1], [lm_maxMag, n_r_ic_maxMag])
            call write_dataset(groupDtFields_id, 'dbdt_icLast', complex_t, dbdt_icLast, [llmMag, 1], [lm_maxMag, n_r_ic_maxMag])
            call write_dataset(groupDtFields_id, 'djdt_icLast', complex_t, djdt_icLast, [llmMag, 1], [lm_maxMag, n_r_ic_maxMag])
            call write_dataset(groupFields_id,   'b_ic_pol',    complex_t, b_ic,        [llmMag, 1], [lm_maxMag, n_r_ic_maxMag])
         end if
      end if

      ! close groups
      call h5gclose_f(groupFields_id, error)
      assert_equal(error, 0, "h5gclose_f")
      call h5gclose_f(groupDtFields_id, error)
      assert_equal(error, 0, "h5gclose_f")

      ! release datatype
      call h5tclose_f(complex_t, error)
      assert_equal(error, 0, "h5tclose_f")

      ! Close the file.
      call h5fclose_f(file_id, error)
      assert_equal(error, 0, "h5fclose_f")

      ! Close FORTRAN interface.
      call h5close_f(error)
      assert_equal(error, 0, "h5close_f")

      PERFOFF

   end subroutine storeHdf5
#endif
end module storeCheckPoints
