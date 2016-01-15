#include "perflib_preproc.cpp"
module hdf5Helpers
   !
   ! This module contains several useful tools to manipulate HDF5 files
   !

   use, intrinsic :: iso_fortran_env, only : error_unit
   use precision_mod
   use blocking, only: st_map, lo_map
   use LMLoop_data, only: llm, ulm
   use hdf5

   implicit none

   private

   interface readHdf5_attribute
      module procedure readHdf5_attr_int
      module procedure readHdf5_attr_dble
      !module procedure readHdf5_attr_string
   end interfacE readHdf5_attribute

   interface writeHdf5_attribute
      module procedure writeHdf5_attr_int
      module procedure writeHdf5_attr_dble
      module procedure writeHdf5_attr_string
   end interface writeHdf5_attribute

   public :: readHdf5_attribute, writeHdf5_attribute, write_dataset

contains

#define check(error, msg) call x_check(error, msg, __FILE__, __LINE__)
   subroutine x_check(error, msg, file, line)
      use parallel_mod, only : rank, mpi_abort, MPI_COMM_WORLD
      logical, intent(in) :: error
      character(len=*), intent(in) :: msg
      character(len=*), intent(in) :: file
      integer, intent(in) :: line

      integer ierr

      if (error) then
        write(error_unit, '(a,i0,a,i0)') "Task #", rank, ": " // msg // ": error condition at " // file // ":", line
        call mpi_abort(MPI_COMM_WORLD, 1, ierr)
        stop 1
      endif
   end subroutine

   subroutine write_dataset(loc_id, dataset_name, dataset_type, dat, dim1, dims_full)
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
      check(error /= 0, "h5screate_simple_f")

      ! create the dataset in the file
      call h5dcreate_f(loc_id, dataset_name, dataset_type, dspace_id, dset_id, error)
      check(error /= 0, "h5dcreate_f()")

      call h5sclose_f(dspace_id, error)
      check(error /= 0, "h5sclose_f()")


      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      check(error /= 0, "h5pcreate_f")
#ifdef WITH_MPI
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
      check(error /= 0, "h5pset_dxpl_mpio_f")
#endif

      ! describe the hyperslab of the local data in the file
      off  = [0, llm - 1]
      dims = [dim1, ulm - llm + 1]
      call h5dget_space_f(dset_id, dspace_id, error)
      check(error /= 0, "h5dget_space_f")
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, off, dims, error)
      check(error /= 0, "h5sselect_hyperslab_f")

      ! describe the layout of 'dat' in memory
      call h5screate_simple_f(2, dims, memspace, error)
      check(error /= 0, "h5screate_simple_f()")

      ! collectively write
      call h5dwrite_f(dset_id, dataset_type, C_LOC(dat_transposed(:,:)), error, &
                      file_space_id=dspace_id, mem_space_id=memspace, &
                      xfer_prp=plist_id)
      check(error /= 0, "h5dwrite_f")

      ! cleanup
      call h5sclose_f(memspace, error)
      check(error /= 0, "h5sclose_f")
      call h5sclose_f(dspace_id, error)
      check(error /= 0, "h5sclose_f()")
      call h5pclose_f(plist_id, error)
      check(error /= 0, "h5pclose_f")
      call h5dclose_f(dset_id, error)
      check(error /= 0, "h5dclose_f")

      PERFOFF
   end subroutine write_dataset
!------------------------------------------------------------------------------
   subroutine readHdf5_attr_dble(loc_id,attr_name,attr_value)
      use parallel_mod, only : rank, mpi_abort, MPI_COMM_WORLD

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
      check(error /= 0, "h5aexists_f")
      if ( attr_exists ) then
         call h5aopen_f(loc_id, attr_name, attr_id, error)
         check(error /= 0, "h5aopen_f")
         call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, attr_value, adims, error)
         check(error /= 0, "h5aread_f")
         call h5aclose_f(attr_id, error)
         check(error /= 0, "h5aclose_f")
      else
         write(error_unit, '(a,i0,a,i0)') "Task #", rank, ": attribute '" // attr_name // "' does not exist, aborting!"
         call mpi_abort(MPI_COMM_WORLD, 1, error)
         stop 1
      end if

   end subroutine readHdf5_attr_dble
!------------------------------------------------------------------------------
   subroutine readHdf5_attr_int(loc_id,attr_name,attr_value)
      use parallel_mod, only : rank, mpi_abort, MPI_COMM_WORLD

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
      check(error /= 0, "h5aexists_f")
      if ( attr_exists ) then
         call h5aopen_f(loc_id, attr_name, attr_id, error)
         check(error /= 0, "h5aopen_f")
         call h5aread_f(attr_id, H5T_NATIVE_integer, attr_value, adims, error)
         check(error /= 0, "h5aread_f")
         call h5aclose_f(attr_id, error)
         check(error /= 0, "h5aclose_f")
      else
         write(error_unit, '(a,i0,a,i0)') "Task #", rank, ": attribute '" // attr_name // "' does not exist, aborting!"
         call mpi_abort(MPI_COMM_WORLD, 1, error)
         stop 1
      end if

   end subroutine readHdf5_attr_int
!------------------------------------------------------------------------------
   subroutine readHdf5_attr_string(loc_id,attr_name,attr_value)
      use parallel_mod, only : rank, mpi_abort, MPI_COMM_WORLD

      !--- Input variables
      character(len=*), intent(in) :: attr_name
      integer(HID_T),   intent(in) :: loc_id

      !--- Output variables
      character(len=*), intent(out) :: attr_value

      !--- Local variables
      integer(HSIZE_T), dimension(1) :: adims = [1]  ! Attribute dimensions
      integer(HID_T) :: attr_id, dtype_id
      integer(SIZE_T) :: attribute_string_length
      integer :: error
      logical :: attr_exists

      call h5aexists_f(loc_id, attr_name, attr_exists, error)
      check(error /= 0, "h5aexists_f")
      if ( attr_exists ) then
         call h5aopen_f(loc_id, attr_name, attr_id, error)
         check(error /= 0, "h5aopen_f")

         ! get length of string attribute
         call h5aget_type_f(attr_id, dtype_id, error)
         check(error /= 0, "h5aget_type_f")
         call h5tget_size_f(dtype_id, attribute_string_length, error)
         check(error /= 0, "h5tget_size_f")
         check(attribute_string_length > len(attr_value), "attribute longer than supplied string")

         attr_value = ""

         call h5aread_f(attr_id, dtype_id, attr_value(1:attribute_string_length), adims, error)
         check(error /= 0, "h5aread_f")
         call h5tclose_f(dtype_id, error)
         check(error /= 0, "h5tclose_f")
         call h5aclose_f(attr_id, error)
         check(error /= 0, "h5aclose_f")
      else
         write(error_unit, '(a,i0,a,i0)') "Task #", rank, ": attribute '" // attr_name // "' does not exist, aborting!"
         call mpi_abort(MPI_COMM_WORLD, 1, error)
         stop 1
      end if

   end subroutine readHdf5_attr_string
!------------------------------------------------------------------------------
   subroutine writeHdf5_attr_dble(loc_id,attr_name,attr_value)

      !--- Input variables
      character(len=*),intent(in) :: attr_name
      integer(HID_T),  intent(in) :: loc_id
      real(cp),    intent(in) :: attr_value

      !--- Local variables
      integer(HSIZE_T), dimension(1) :: adims = [1]  ! Attribute dimensions
      integer(HID_T) :: attr_id, dspace_id
      integer :: error
      logical :: attr_exists

      call h5screate_f(H5S_SCALAR_F,dspace_id,error)
      check(error /= 0, "h5screate_f")
      call h5acreate_f(loc_id, attr_name, H5T_NATIVE_DOUBLE, dspace_id, attr_id, error)
      check(error /= 0, "h5acreate_f")

      call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_value, adims, error)
      check(error /= 0, "h5awrite_f")

      call h5aclose_f(attr_id, error)
      check(error /= 0, "h5aclose_f")
      call h5sclose_f(dspace_id, error)
      check(error /= 0, "h5sclose_f")

   end subroutine writeHdf5_attr_dble
!------------------------------------------------------------------------------
   subroutine writeHdf5_attr_int(loc_id,attr_name,attr_value)

      !--- Input variables
      character(len=*),intent(in) :: attr_name
      integer(HID_T),  intent(in) :: loc_id
      integer,         intent(in) :: attr_value

      !--- Local variables
      integer(HSIZE_T), dimension(1) :: adims = [1]  ! Attribute dimensions
      integer(HID_T) :: attr_id, dspace_id
      integer :: error
      logical :: attr_exists

      call h5screate_f(H5S_SCALAR_F,dspace_id,error)
      check(error /= 0, "h5screate_f")
      call h5acreate_f(loc_id, attr_name, H5T_NATIVE_integer, dspace_id, attr_id, error)
      check(error /= 0, "h5acreate_f")

      call h5awrite_f(attr_id, H5T_NATIVE_integer, attr_value, adims, error)
      check(error /= 0, "h5awrite_f")

      call h5aclose_f(attr_id, error)
      check(error /= 0, "h5aclose_f")
      call h5sclose_f(dspace_id, error)
      check(error /= 0, "h5sclose_f")

   end subroutine writeHdf5_attr_int
!------------------------------------------------------------------------------
   subroutine writeHdf5_attr_string(loc_id,attr_name,attr_value)

      !--- Input variables
      character(len=*),intent(in) :: attr_name
      integer(HID_T),  intent(in) :: loc_id
      character(len=*),intent(in) :: attr_value

      !--- Local variables
      integer(HSIZE_T), dimension(1) :: adims = [1]  ! Attribute dimensions
      integer(HID_T) :: attr_id, dspace_id, dtype_id
      integer :: error
      logical :: attr_exists

      call h5screate_f(H5S_SCALAR_F,dspace_id,error)
      check(error /= 0, "h5screate_f")

      ! create datatype for this specific string length
      call h5tcopy_f(H5T_NATIVE_CHARACTER, dtype_id, error)
      check(error /= 0, "h5tcopy_f")
      call h5tset_size_f(dtype_id, int(len(trim(attr_value)), kind=size_t), error)
      check(error /= 0, "h5tset_size_f")

      call h5acreate_f(loc_id, attr_name, dtype_id, dspace_id, attr_id, error)
      check(error /= 0, "h5acreate_f")

      call h5awrite_f(attr_id, dtype_id, attr_value, adims, error)
      check(error /= 0, "h5awrite_f")

      call h5aclose_f(attr_id, error)
      check(error /= 0, "h5aclose_f")
      call h5tclose_f(dtype_id, error)
      check(error /= 0, "h5tclose_f")
      call h5sclose_f(dspace_id, error)
      check(error /= 0, "h5sclose_f")

   end subroutine writeHdf5_attr_string
!------------------------------------------------------------------------------
end module hdf5Helpers
