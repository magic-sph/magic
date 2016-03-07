#include "assertions.cpp"
#include "perflib_preproc.cpp"
module hdf5Helpers
   !
   ! This module contains several useful tools to manipulate HDF5 files
   !
   use assertions
   use precision_mod
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
      module procedure writeHdf5_attr_string
   end interface writeHdf5_attribute

   interface write_dataset
      module procedure write_distributed_dataset
      module procedure write_replicated_dataset
   end interface

   public :: readHdf5_attribute, writeHdf5_attribute
   public :: read_dataset, interpolate_dataset, write_dataset
   public :: dataset_shape
   public :: hdf5_fortran_complex_type
   public :: HSIZE_T, HID_T

contains

!------------------------------------------------------------------------------

   ! Store a replicated dataset (where all task have the same array with identical values)
   subroutine write_replicated_dataset(loc_id, dataset_name, dataset_type, dat)
      use parallel_mod, only : rank, MPI_COMM_WORLD
      integer(HID_T),   intent(in) :: loc_id
      character(len=*), intent(in) :: dataset_name
      integer(HID_T),   intent(in) :: dataset_type
      integer, intent(in), target  :: dat(:,:)

      !--- Local variables
      integer(HID_T) :: dspace_id, dset_id, memspace_id, plist_id
      integer :: error

      PERFON('write_replicated_dataset')

      ! layout of the dataset in the file
      call h5screate_simple_f(2, int(shape(dat), kind=HSIZE_T), dspace_id, error)
      assert_equal(error, 0, "h5screate_simple_f")

      ! create the dataset in the file
      call h5dcreate_f(loc_id, dataset_name, dataset_type, dspace_id, dset_id, error)
      assert_equal(error, 0, "h5dcreate_f()")

      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      assert_equal(error, 0, "h5pcreate_f")
#ifdef WITH_MPI
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
      assert_equal(error, 0, "h5pset_dxpl_mpio_f")
#endif

      ! describe the layout of 'dat' in memory
      call h5screate_simple_f(2, int(shape(dat), kind=HSIZE_T), memspace_id, error)
      assert_equal(error, 0, "h5screate_simple_f()")

      if (rank /= 0) then
        ! only rank 0 actually writes this array
        call h5sselect_none_f(dspace_id, error)
        call h5sselect_none_f(memspace_id, error)
      endif

      ! collectively write
      call h5dwrite_f(dset_id, dataset_type, C_LOC(dat), error, &
                      mem_space_id=memspace_id, file_space_id=dspace_id, &
                      xfer_prp=plist_id)
      assert_equal(error, 0, "h5dwrite_f")

      ! cleanup
      call h5sclose_f(memspace_id, error)
      assert_equal(error, 0, "h5sclose_f")
      call h5sclose_f(dspace_id, error)
      assert_equal(error, 0, "h5sclose_f()")
      call h5pclose_f(plist_id, error)
      assert_equal(error, 0, "h5pclose_f")
      call h5dclose_f(dset_id, error)
      assert_equal(error, 0, "h5dclose_f")

      PERFOFF

   end subroutine write_replicated_dataset

!------------------------------------------------------------------------------

   subroutine write_distributed_dataset(loc_id, dataset_name, dataset_type, dat, lbound_local, dims_global)
      !--- Input variables
      integer(HID_T),   intent(in) :: loc_id
      character(len=*), intent(in) :: dataset_name
      integer(HID_T),   intent(in) :: dataset_type
      integer,          intent(in) :: lbound_local(2), dims_global(2)
      complex(cp), intent(in) :: dat(lbound_local(1):,lbound_local(2):)
      complex(cp), target :: dat_transposed(size(dat, dim=2), size(dat, dim=1))

      !--- Local variables
      integer(HSIZE_T) :: dims(2)
      integer(HSIZE_T) :: off(2)
      integer(HID_T) :: dspace_id, dset_id, memspace_id, plist_id
      integer :: error

      PERFON('write_distributed_dataset')

      ! global layout of the dataset in the file
      dims = [dims_global(2), dims_global(1)]
      call h5screate_simple_f(2, dims, dspace_id, error)
      assert_equal(error, 0, "h5screate_simple_f")

      ! create the dataset in the file
      call h5dcreate_f(loc_id, dataset_name, dataset_type, dspace_id, dset_id, error)
      assert_equal(error, 0, "h5dcreate_f()")

      call h5sclose_f(dspace_id, error)
      assert_equal(error, 0, "h5sclose_f()")


      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      assert_equal(error, 0, "h5pcreate_f")
#ifdef WITH_MPI
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
      assert_equal(error, 0, "h5pset_dxpl_mpio_f")
#endif

      ! describe the hyperslab of the local data in the file
      off  = [lbound_local(2) - 1, lbound_local(1) - 1]
      dims = [size(dat, dim=2), size(dat, dim=1)]
      call h5dget_space_f(dset_id, dspace_id, error)
      assert_equal(error, 0, "h5dget_space_f")
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, off, dims, error)
      assert_equal(error, 0, "h5sselect_hyperslab_f")

      ! describe the layout of 'dat' in memory
      call h5screate_simple_f(2, dims, memspace_id, error)
      assert_equal(error, 0, "h5screate_simple_f()")

      ! collectively write
      dat_transposed(:,:) = transpose(dat)
      call h5dwrite_f(dset_id, dataset_type, C_LOC(dat_transposed(:,:)), error, &
                      mem_space_id=memspace_id, file_space_id=dspace_id, &
                      xfer_prp=plist_id)
      assert_equal(error, 0, "h5dwrite_f")

      ! cleanup
      call h5sclose_f(memspace_id, error)
      assert_equal(error, 0, "h5sclose_f")
      call h5sclose_f(dspace_id, error)
      assert_equal(error, 0, "h5sclose_f()")
      call h5pclose_f(plist_id, error)
      assert_equal(error, 0, "h5pclose_f")
      call h5dclose_f(dset_id, error)
      assert_equal(error, 0, "h5dclose_f")

      PERFOFF

   end subroutine write_distributed_dataset

!------------------------------------------------------------------------------

   subroutine read_dataset(loc_id, dataset_name, dataset_type, dat)
      integer(HID_T),   intent(in) :: loc_id
      character(len=*), intent(in) :: dataset_name
      integer(HID_T),   intent(in) :: dataset_type
      integer, intent(out), target  :: dat(:,:)

      !--- Local variables
      logical :: dataspace_is_simple
      integer(HID_T) :: dspace_id, dset_id, memspace_id, plist_id
      integer(HSIZE_T) :: dat_shape(2), dataset_shape(2), max_dims(2), offset(2)
      integer :: ndim
      integer :: error

      PERFON('read_dataset')

      ! open dataset
      call h5dopen_f(loc_id, dataset_name, dset_id, error)
      assert_equal(error, 0, "h5dopen_f")

      ! Get its dataspace
      call h5dget_space_f(dset_id, dspace_id, error)
      assert_equal(error, 0, "h5dget_space_f")

      ! Some sanity checks if dat(:,:) and the dataset are compatible
      call h5sis_simple_f(dspace_id, dataspace_is_simple, error)
      assert_equal(error, 0, "h5sis_simple_f(), dataset: " // dataset_name)
      assert_true(dataspace_is_simple, "dataset " // dataset_name // " has complex dataspace")

      call h5sget_simple_extent_ndims_f(dspace_id, ndim, error)
      assert_equal(error, 0, "h5sget_simple_extent_ndims_f(), dataset " // dataset_name)
      assert_equal(ndim, 2, "dataset " // dataset_name // " has unsupported rank")

      dat_shape = shape(dat)
      call h5sget_simple_extent_dims_f(dspace_id, dataset_shape, max_dims, error)
      assert_false(error == -1, "h5sget_simple_extent_dims_f(), dataset " // dataset_name)
      assert_true(all(dataset_shape .eq. dat_shape), "wrong shape for dataset " // dataset_name)

      ! Release the dataspace
      call h5sclose_f(dspace_id, error)
      assert_equal(error, 0, "h5sclose_f(), dataset " // dataset_name)

      ! Select hyperslab in the file.
      call h5screate_simple_f(2, dat_shape, memspace_id, error)
      assert_equal(error, 0, "h5screate_simple_f(), dataset " // dataset_name)

      call h5dget_space_f(dset_id, dspace_id, error)
      assert_equal(error, 0, "h5dget_space_f(), dataset " // dataset_name)

      offset(:) = 0
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, dat_shape, error)
      assert_equal(error, 0, "h5sselect_hyperslab_f(), dataset " // dataset_name)

      ! Create property list for collective dataset write.
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      assert_equal(error, 0, "h5pcreate_f(), dataset " // dataset_name)

      ! Do the read
      call h5dread_f(dset_id, dataset_type, dat, dat_shape, error, &
                     mem_space_id=memspace_id, file_space_id=dspace_id, &
                     xfer_prp=plist_id)
      assert_equal(error, 0, "h5dread_f(), dataset " // dataset_name)

      ! Release allocated objects
      call h5dclose_f(dset_id, error)
      assert_equal(error, 0, "h5dclose_f(), dataset " // dataset_name)

      call h5sclose_f(dspace_id, error)
      assert_equal(error, 0, "h5sclose_f(), dataset " // dataset_name)

      call h5sclose_f(memspace_id, error)
      assert_equal(error, 0, "h5sclose_f(), dataset " // dataset_name)

      call h5pclose_f(plist_id, error)
      assert_equal(error, 0, "h5pclose_f(), dataset " // dataset_name)

      PERFOFF
   end subroutine read_dataset

!------------------------------------------------------------------------------

   subroutine interpolate_dataset(loc_id, dataset_name, dataset_type, dat, lbound_local, file_lm_map, lBc, l_IC)
      use, intrinsic :: ieee_arithmetic
      use parallel_mod, only : rank, MPI_COMM_WORLD
      use blocking, only: lo_map
      use readCheckPoints_helper, only : mapDataR_copy
      integer,          intent(in) :: lbound_local(2)
      integer(HID_T),   intent(in) :: loc_id
      integer(HID_T),   intent(in) :: dataset_type
      character(len=*), intent(in) :: dataset_name
      complex(cp), intent(out) :: dat(lbound_local(1):,lbound_local(2):)
      integer, intent(in) :: file_lm_map(0:, 0:)
      logical, intent(in) :: lBc,l_IC

      !--- Local variables
      complex(cp), target :: dat_transposed(lbound_local(2):lbound_local(2) + size(dat, dim=2) - 1, &
                                            lbound_local(1):lbound_local(1) + size(dat, dim=1) - 1)
      logical :: dataspace_is_simple, lm_missing
      integer(HID_T) :: dspace_id, dset_id, memspace_id, plist_id
      integer(HSIZE_T) :: blockshape(2), dataset_shape(2), max_dims(2), file_offset(2)
      type(c_ptr) :: ptr
      integer :: ndim
      integer :: error
      integer :: l, m
      integer :: lm, lm_file_start
      integer(HSIZE_T) :: nblock

      PERFON('interpolate_dataset')

      !write(*,'(a,i0,a,i0,a,i0,a,i0,a)') "Reading " // &
      !           dataset_name // "(", lbound_local(1), ":", lbound_local(1) + size(dat, dim=1) - 1, ", ", &
      !                                lbound_local(2), ":", lbound_local(2) + size(dat, dim=2) - 1, ")"

      ! initialize output array with sNaN
      ! to abort any computiation with unitilized data, in
      ! case there is an error in this function
      dat(:, :) = cmplx(ieee_value(1.0_cp,ieee_signaling_nan), &
                        ieee_value(1.0_cp,ieee_signaling_nan))

      ! Open the dataset
      call h5dopen_f(loc_id, dataset_name, dset_id, error)
      assert_equal(error, 0, "h5dopen_f(), dataset " // dataset_name)

      ! Get its dataspace
      call h5dget_space_f(dset_id, dspace_id, error)
      assert_equal(error, 0, "h5dget_space_f(), dataset " // dataset_name)

      ! Make sanity checks about the shape
      call h5sis_simple_f(dspace_id, dataspace_is_simple, error)
      assert_equal(error, 0, "h5sis_simple_f(), dataset " // dataset_name)
      assert_true(dataspace_is_simple, "data field: " // dataset_name)

      call h5sget_simple_extent_ndims_f(dspace_id, ndim, error)
      assert_equal(error, 0, "h5sget_simple_extent_ndims_f(), dataset " // dataset_name)
      assert_equal(ndim, 2, "dataset " // dataset_name // " has unsupported rank")

      call h5sget_simple_extent_dims_f(dspace_id, dataset_shape, max_dims, error)
      assert_false(error == -1, "h5sget_simple_extent_dims_f(), dataset " // dataset_name)

      ! Release the dataspace
      call h5sclose_f(dspace_id, error)
      assert_equal(error, 0, "h5sclose_f(), dataset " // dataset_name)

      ! Create property list for collective dataset write.
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      assert_equal(error, 0, "h5pcreate_f(), dataset " // dataset_name)

      ! read the dataset in contigous LM blocks
      lm = lbound_local(1)
      do while (lm < lbound_local(1) + size(dat, dim=1))

         ! physical (l,m) of the current simulation's mapping
         l = lo_map%lm2l(lm)
         m = lo_map%lm2m(lm)

         lm_missing = .false.
         if (l < 0 .or. l > ubound(file_lm_map, dim=1)) then
           lm_missing = .true.
         else if (m < 0 .or. m > ubound(file_lm_map, dim=2)) then
           lm_missing = .true.
         else
           lm_file_start = file_lm_map(l, m)
           if (lm_file_start == 0) then
             lm_missing = .true.
           endif
         endif

         if (lm_missing) then
            write(error_unit, '(a,i0,a,i0,a,i0,a)') &
                "Task #", rank, ": restart file contains no entry for (l,m) = (", l, ",", m, ")"
            call mpi_abort(MPI_COMM_WORLD, 1, error)
            stop 1
         endif

         ! find the end of a contigous block to read from the restart file
         ! i.e. increase nblock until there is a gap or until we reach the end of our (l,m) region
         nblock = 1
         do while (lm + nblock <= lbound_local(1) + size(dat, dim=1) - 1)
             if (file_lm_map(lo_map%lm2l(lm + nblock), lo_map%lm2m(lm + nblock)) == lm_file_start + nblock) then
               nblock = nblock + 1
             else
               exit
             endif
         end do
         ! contigous block in dataset(:, lm : lm + nblock - 1)

         blockshape = [dataset_shape(1), nblock]
         file_offset  = [0, lm_file_start - 1]

         ! describe local memory layout
         call h5screate_simple_f(2, blockshape, memspace_id, error)
         assert_equal(error, 0, "h5screate_simple_f(), dataset " // dataset_name)

         ! describe layout in file
         call h5dget_space_f(dset_id, dspace_id, error)
         assert_equal(error, 0, "h5dget_space_f(), dataset " // dataset_name)
         call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, file_offset, blockshape, error)
         assert_equal(error, 0, "h5sselect_hyperslab_f(), dataset " // dataset_name)

         ! Do the read
         !write(*,'(a,i0,a,i0)') " reading LM block ", lm, ":", lm + nblock - 1
         ptr = c_loc(dat_transposed(:, lm : lm + nblock))
         call h5dread_f(dset_id, dataset_type, ptr, error, &
                        mem_space_id=memspace_id, file_space_id=dspace_id, &
                        xfer_prp=plist_id)
         assert_equal(error, 0, "h5dread_f(), dataset " // dataset_name)

         call h5sclose_f(dspace_id, error)
         assert_equal(error, 0, "h5sclose_f(), dataset " // dataset_name)

         call h5sclose_f(memspace_id, error)
         assert_equal(error, 0, "h5sclose_f(), dataset " // dataset_name)

         lm = lm + nblock
      end do

      ! Release allocated objects
      call h5dclose_f(dset_id, error)
      assert_equal(error, 0, "h5dclose_f(), dataset " // dataset_name)

      call h5pclose_f(plist_id, error)
      assert_equal(error, 0, "h5pclose_f(), dataset " // dataset_name)

      ! interpolate and transpose 'dat_transposed' into result array 'dat'
      do lm = lbound_local(1), lbound_local(1) + size(dat, dim=1) - 1
        dat(lm,:) = mapDataR_copy(dat_transposed(:,lm), size(dat, dim=2), int(dataset_shape(1)), lBc, l_IC)
      end do

      PERFOFF
   end subroutine interpolate_dataset

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
      assert_equal(error, 0, "h5aexists_f")
      if ( attr_exists ) then
         call h5aopen_f(loc_id, attr_name, attr_id, error)
         assert_equal(error, 0, "h5aopen_f")
         call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, attr_value, adims, error)
         assert_equal(error, 0, "h5aread_f")
         call h5aclose_f(attr_id, error)
         assert_equal(error, 0, "h5aclose_f")
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
      assert_equal(error, 0, "h5aexists_f")
      if ( attr_exists ) then
         call h5aopen_f(loc_id, attr_name, attr_id, error)
         assert_equal(error, 0, "h5aopen_f")
         call h5aread_f(attr_id, H5T_NATIVE_integer, attr_value, adims, error)
         assert_equal(error, 0, "h5aread_f")
         call h5aclose_f(attr_id, error)
         assert_equal(error, 0, "h5aclose_f")
      else
         write(error_unit, '(a,i0,a,i0)') "Task #", rank, ": attribute '" // attr_name // "' does not exist, aborting!"
         call mpi_abort(MPI_COMM_WORLD, 1, error)
         stop 1
      end if

   end subroutine readHdf5_attr_int

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

      call h5screate_f(H5S_SCALAR_F,dspace_id,error)
      assert_equal(error, 0, "h5screate_f")
      call h5acreate_f(loc_id, attr_name, H5T_NATIVE_DOUBLE, dspace_id, attr_id, error)
      assert_equal(error, 0, "h5acreate_f")

      call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_value, adims, error)
      assert_equal(error, 0, "h5awrite_f")

      call h5aclose_f(attr_id, error)
      assert_equal(error, 0, "h5aclose_f")
      call h5sclose_f(dspace_id, error)
      assert_equal(error, 0, "h5sclose_f")

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

      call h5screate_f(H5S_SCALAR_F,dspace_id,error)
      assert_equal(error, 0, "h5screate_f")
      call h5acreate_f(loc_id, attr_name, H5T_NATIVE_integer, dspace_id, attr_id, error)
      assert_equal(error, 0, "h5acreate_f")

      call h5awrite_f(attr_id, H5T_NATIVE_integer, attr_value, adims, error)
      assert_equal(error, 0, "h5awrite_f")

      call h5aclose_f(attr_id, error)
      assert_equal(error, 0, "h5aclose_f")
      call h5sclose_f(dspace_id, error)
      assert_equal(error, 0, "h5sclose_f")

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

      call h5screate_f(H5S_SCALAR_F,dspace_id,error)
      assert_equal(error, 0, "h5screate_f")

      ! create datatype for this specific string length
      call h5tcopy_f(H5T_NATIVE_CHARACTER, dtype_id, error)
      assert_equal(error, 0, "h5tcopy_f")
      call h5tset_size_f(dtype_id, int(len(trim(attr_value)), kind=size_t), error)
      assert_equal(error, 0, "h5tset_size_f")

      call h5acreate_f(loc_id, attr_name, dtype_id, dspace_id, attr_id, error)
      assert_equal(error, 0, "h5acreate_f")

      call h5awrite_f(attr_id, dtype_id, attr_value, adims, error)
      assert_equal(error, 0, "h5awrite_f")

      call h5aclose_f(attr_id, error)
      assert_equal(error, 0, "h5aclose_f")
      call h5tclose_f(dtype_id, error)
      assert_equal(error, 0, "h5tclose_f")
      call h5sclose_f(dspace_id, error)
      assert_equal(error, 0, "h5sclose_f")

   end subroutine writeHdf5_attr_string

!------------------------------------------------------------------------------

   subroutine dataset_shape(loc_id, dataset_name, dset_shape)
      integer(HID_T),   intent(in)  :: loc_id
      character(len=*), intent(in)  :: dataset_name
      integer(HSIZE_T), intent(out) :: dset_shape(:)

      !--- Local variables
      logical :: dataspace_is_simple
      integer(HID_T) :: dspace_id, dset_id
      integer(HSIZE_T) :: max_dims(size(dset_shape))
      integer :: ndim
      integer :: error

      ! open dataset
      call h5dopen_f(loc_id, dataset_name, dset_id, error)
      assert_equal(error, 0, "h5dopen_f")

      ! Get its dataspace
      call h5dget_space_f(dset_id, dspace_id, error)
      assert_equal(error, 0, "h5dget_space_f")

      ! Some sanity checks if dat(:,:) and the dataset are compatible
      call h5sis_simple_f(dspace_id, dataspace_is_simple, error)
      assert_equal(error, 0, "h5sis_simple_f(), dataset: " // dataset_name)
      assert_true(dataspace_is_simple, "dataset " // dataset_name // " has complex dataspace")

      call h5sget_simple_extent_ndims_f(dspace_id, ndim, error)
      assert_equal(error, 0, "h5sget_simple_extent_ndims_f(), dataset " // dataset_name)
      assert_equal(ndim, size(dset_shape), "dataset "//dataset_name//" has wrong rank for 'dset_shape'")

      call h5sget_simple_extent_dims_f(dspace_id, dset_shape, max_dims, error)
      assert_false(error == -1, "h5sget_simple_extent_dims_f(), dataset " // dataset_name)

      call h5sclose_f(dspace_id, error)
      assert_equal(error, 0, "h5sclose_f(), dataset " // dataset_name)

      call h5dclose_f(dset_id, error)
      assert_equal(error, 0, "h5dclose_f(), dataset " // dataset_name)

   end subroutine dataset_shape

!------------------------------------------------------------------------------

   ! Create a HF compound type that represents Fortran complex data
   function hdf5_fortran_complex_type() result(complex_t)
      integer(HID_T) :: complex_t
      integer(HSIZE_T) :: re_size,im_size,complex_t_size,offset
      integer :: error

      call h5tget_size_f(H5T_NATIVE_DOUBLE,re_size,error)
      assert_equal(error, 0, "h5tget_size_f")

      call h5tget_size_f(H5T_NATIVE_DOUBLE,im_size,error)
      assert_equal(error, 0, "h5tget_size_f")

      complex_t_size = re_size+im_size
      call h5tcreate_f(H5T_COMPOUND_F,complex_t_size,complex_t,error)
      assert_equal(error, 0, "h5tcreate_f")

      offset = 0
      call h5tinsert_f(complex_t, 'real', offset, H5T_NATIVE_DOUBLE, error)
      assert_equal(error, 0, "h5tinsert_f")

      offset = offset + re_size
      call h5tinsert_f(complex_t, 'imag', offset, H5T_NATIVE_DOUBLE, error)
      assert_equal(error, 0, "h5tinsert_f")

   end function

end module hdf5Helpers
