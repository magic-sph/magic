!$Id$
module hdf5Helpers

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
!---------------
   subroutine write_dataset(loc_id, dataset_name, dataset_type, dat, dim1, dims_full)

      use blocking, only: st_map, lo_map
      use LMLoop_data, only: llm, ulm

      !--- Input variables
      integer,          intent(in) :: dim1
      integer(HID_T),   intent(in) :: loc_id
      integer(HID_T),   intent(in) :: dataset_type
      character(len=*), intent(in) :: dataset_name
      complex(kind=8),  intent(in) :: dat(llm:ulm,dim1)
      integer(HSIZE_T), intent(in) :: dims_full(2)

      !--- Local variables
      integer(HSIZE_T) :: dims_loc(2)
      integer(HSSIZE_T) :: off(2)
      integer(HID_T) :: dspace_id, dset_id, memspace, plist_id
      integer :: error
      integer :: l,m,lm

      call h5screate_simple_f(2, dims_full, dspace_id, error)
      call h5dcreate_f(loc_id, dataset_name, dataset_type, dspace_id, dset_id, error)
      call h5sclose_f(dspace_id, error)

      dims_loc = [1,dim1]

      call h5screate_simple_f(2, dims_loc, memspace, error)
      call h5dget_space_f(dset_id, dspace_id, error)
      do lm=llm,ulm
         l = lo_map%lm2l(lm)
         m = lo_map%lm2m(lm)
         off(1)=st_map%lm2(l,m)-1
         off(2)=0
         call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, off, dims_loc, error)
         call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
         call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
         call h5dwrite_f(dset_id, dataset_type, C_LOC(dat(lm,:)), error, &
                         file_space_id=dspace_id, mem_space_id=memspace, &
                         xfer_prp=plist_id)
      enddo
      call h5sclose_f(memspace, error)
      call h5sclose_f(dspace_id, error)
      call h5dclose_f(dset_id, error)
      call h5pclose_f(plist_id, error)

   end subroutine write_dataset
!---------------
   subroutine readHdf5_attr_dble(loc_id,attr_name,attr_value)

      !--- Input variables
      character(len=*),intent(in) :: attr_name
      integer(HID_T),intent(in) :: loc_id
      real(kind=8),intent(out) :: attr_value

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
         attr_value=0.D0
      end if

   end subroutine readHdf5_attr_dble
!---------------
   subroutine readHdf5_attr_int(loc_id,attr_name,attr_value)

      !--- Input variables
      character(len=*),intent(in) :: attr_name
      integer(HID_T),intent(in) :: loc_id
      integer,intent(out) :: attr_value

      !--- Local variables
      integer(HSIZE_T),dimension(1) :: adims = [1]  ! Attribute dimensions
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
!---------------
   subroutine writeHdf5_attr_dble(loc_id,aspace_id,attr_name,attr_value)

      !--- Input variables
      character(len=*),intent(in) :: attr_name
      integer(HID_T),  intent(in) :: aspace_id
      integer(HID_T),  intent(in) :: loc_id
      real(kind=8),    intent(in) :: attr_value

      !--- Local variables
      integer(HSIZE_T),DIMENSION(1) :: adims = [1]  ! Attribute dimensions
      integer(HID_T) :: attr_id
      integer :: error
      logical :: attr_exists

      call h5acreate_f(loc_id, attr_name, H5T_NATIVE_DOUBLE, aspace_id, attr_id, error)
      call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_value, adims, error)
      call h5aclose_f(attr_id, error)

   end subroutine writeHdf5_attr_dble
!---------------
   subroutine writeHdf5_attr_int(loc_id,aspace_id,attr_name,attr_value)

      !--- Input variables
      character(len=*),intent(in) :: attr_name
      integer(HID_T),  intent(in) :: aspace_id
      integer(HID_T),  intent(in) :: loc_id
      integer,         intent(in) :: attr_value

      !--- Local variables
      integer(HSIZE_T),DIMENSION(1) :: adims = [1]  ! Attribute dimensions
      integer(HID_T) :: attr_id
      integer :: error
      logical :: attr_exists

      call h5acreate_f(loc_id, attr_name, H5T_NATIVE_integer, aspace_id, attr_id, error)
      call h5awrite_f(attr_id, H5T_NATIVE_integer, attr_value, adims, error)
      call h5aclose_f(attr_id, error)

   end subroutine writeHdf5_attr_int
!---------------
end module hdf5Helpers
