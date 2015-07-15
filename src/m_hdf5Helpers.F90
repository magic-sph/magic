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

   public :: readHdf5_attribute, writeHdf5_attribute

contains
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
