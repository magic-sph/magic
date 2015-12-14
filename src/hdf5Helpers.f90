#include "perflib_preproc.cpp"
module hdf5Helpers
   !
   ! This module contains several useful tools to manipulate HDF5 files
   !

   use precision_mod
   use blocking, only: st_map, lo_map
   use LMLoop_data, only: llm, ulm
   use hdf5

   implicit none

   private

   interface readHdf5_attribute
      module procedure readHdf5_attr_int
      module procedure readHdf5_attr_dble
   end interfacE readHdf5_attribute

   interface writeHdf5_attribute
      module procedure writeHdf5_attr_int
      module procedure writeHdf5_attr_dble
   end interface writeHdf5_attribute

   public :: readHdf5_attribute, writeHdf5_attribute, write_dataset

contains

#define check_error(msg) call x_check_error(msg, error, __FILE__, __LINE__)
   subroutine x_check_error(msg, error, file, line)
      use parallel_mod, only : rank
      use, intrinsic :: iso_fortran_env, only : error_unit
      character(len=*), intent(in) :: msg
      integer, intent(in) :: error
      character(len=*), intent(in) :: file
      integer, intent(in) :: line

      if (error /= 0) then
        write(error_unit, '(a,i0,a,i0)') "Task #", rank, ": " // msg // ": error /= 0 at " // file // ":", line
        stop 1
      endif
   end subroutine

   subroutine write_dataset(loc_id, dataset_name, dataset_type, dat, dim1, dims_full)
      use, intrinsic :: iso_fortran_env, only : error_unit
      use parallel_mod, only : rank

      !--- Input variables
      integer,          intent(in) :: dim1
      integer(HID_T),   intent(in) :: loc_id
      integer(HID_T),   intent(in) :: dataset_type
      character(len=*), intent(in) :: dataset_name
      complex(cp), intent(in) :: dat(llm:ulm,dim1)
      complex(cp), target :: dat_transposed(dim1, llm:ulm)
      integer(HSIZE_T), intent(in) :: dims_full(2)

      !--- Local variables
      integer(HSIZE_T) :: dims(2)
      integer(HSIZE_T) :: off(2)
      integer(HID_T) :: dspace_id, dset_id, memspace, plist_id
      integer :: error

      dat_transposed(:,:) = transpose(dat)

      PERFON('write_dataset')

      ! global layout of the dataset in the file
      dims = [dims_full(2), dims_full(1)]
      call h5screate_simple_f(2, dims, dspace_id, error)
      check_error("h5screate_simple_f")

      ! create the dataset in the file
      call h5dcreate_f(loc_id, dataset_name, dataset_type, dspace_id, dset_id, error)
      check_error("h5dcreate_f()")

      call h5sclose_f(dspace_id, error)
      check_error("h5sclose_f()")


      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      check_error("h5pcreate_f")
#ifdef WITH_MPI
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
      check_error("h5pset_dxpl_mpio_f")
#endif

      ! describe the hyperslab of the local data in the file
      off  = [0, llm - 1]
      dims = [dim1, ulm - llm + 1]
      call h5dget_space_f(dset_id, dspace_id, error)
      check_error("h5dget_space_f")
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, off, dims, error)
      check_error("h5sselect_hyperslab_f")

      ! describe the layout of 'dat' in memory
      call h5screate_simple_f(2, dims, memspace, error)
      check_error("h5screate_simple_f()")

      ! collectively write
      call h5dwrite_f(dset_id, dataset_type, C_LOC(dat_transposed(:,:)), error, &
                      file_space_id=dspace_id, mem_space_id=memspace, &
                      xfer_prp=plist_id)
      check_error("h5dwrite_f")

      ! cleanup
      call h5sclose_f(memspace, error)
      check_error("h5sclose_f")
      call h5sclose_f(dspace_id, error)
      check_error("h5sclose_f()")
      call h5pclose_f(plist_id, error)
      check_error("h5pclose_f")
      call h5dclose_f(dset_id, error)
      check_error("h5dclose_f")

      PERFOFF
   end subroutine write_dataset
!------------------------------------------------------------------------------
   subroutine readHdf5_attr_dble(loc_id,attr_name,attr_value)

      !--- Input variables
      character(len=*), intent(in) :: attr_name
      integer(HID_T),   intent(in) :: loc_id

      !--- Output variables
      real(cp), intent(out) :: attr_value

      !--- Local variables
      integer(HSIZE_T),dimension(1) :: adims = [1]  ! Attribute dimensions
      integer(HID_T) :: attr_id
      integer :: error
      logical :: attr_exists

      call h5aexists_f(loc_id, attr_name, attr_exists, error)
      if ( attr_exists ) then
         call h5aopen_f(loc_id, attr_name, attr_id, error)
         call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, attr_value, adims, error)
         call h5aclose_f(attr_id, error)
      else
         attr_value=0.0_cp
      end if

   end subroutine readHdf5_attr_dble
!------------------------------------------------------------------------------
   subroutine readHdf5_attr_int(loc_id,attr_name,attr_value)

      !--- Input variables
      character(len=*), intent(in) :: attr_name
      integer(HID_T),   intent(in) :: loc_id

      !--- Output variables
      integer, intent(out) :: attr_value

      !--- Local variables
      integer(HSIZE_T), dimension(1) :: adims = [1]  ! Attribute dimensions
      integer(HID_T) :: attr_id
      integer :: error
      logical :: attr_exists

      call h5aexists_f(loc_id, attr_name, attr_exists, error)
      if ( attr_exists ) then
         call h5aopen_f(loc_id, attr_name, attr_id, error)
         call h5aread_f(attr_id, H5T_NATIVE_integer, attr_value, adims, error)
         call h5aclose_f(attr_id, error)
      else
         attr_value=0
      end if

   end subroutine readHdf5_attr_int
!------------------------------------------------------------------------------
   subroutine writeHdf5_attr_dble(loc_id,aspace_id,attr_name,attr_value)

      !--- Input variables
      character(len=*),intent(in) :: attr_name
      integer(HID_T),  intent(in) :: aspace_id
      integer(HID_T),  intent(in) :: loc_id
      real(cp),    intent(in) :: attr_value

      !--- Local variables
      integer(HSIZE_T), dimension(1) :: adims = [1]  ! Attribute dimensions
      integer(HID_T) :: attr_id
      integer :: error
      logical :: attr_exists

      call h5acreate_f(loc_id, attr_name, H5T_NATIVE_DOUBLE, aspace_id, attr_id, error)
      call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_value, adims, error)
      call h5aclose_f(attr_id, error)

   end subroutine writeHdf5_attr_dble
!------------------------------------------------------------------------------
   subroutine writeHdf5_attr_int(loc_id,aspace_id,attr_name,attr_value)

      !--- Input variables
      character(len=*),intent(in) :: attr_name
      integer(HID_T),  intent(in) :: aspace_id
      integer(HID_T),  intent(in) :: loc_id
      integer,         intent(in) :: attr_value

      !--- Local variables
      integer(HSIZE_T), dimension(1) :: adims = [1]  ! Attribute dimensions
      integer(HID_T) :: attr_id
      integer :: error
      logical :: attr_exists

      call h5acreate_f(loc_id, attr_name, H5T_NATIVE_integer, aspace_id, attr_id, error)
      call h5awrite_f(attr_id, H5T_NATIVE_integer, attr_value, adims, error)
      call h5aclose_f(attr_id, error)

   end subroutine writeHdf5_attr_int
!------------------------------------------------------------------------------
end module hdf5Helpers
