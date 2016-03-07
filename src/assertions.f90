module assertions
   use, intrinsic :: iso_fortran_env, only : error_unit
   use parallel_mod, only : rank
   use mpi, only : mpi_abort, MPI_COMM_WORLD
   implicit none

   interface x_assert_equal
      module procedure x_assert_equal_integer
   end interface

   public

   contains

   subroutine x_assert_equal_integer(value, expected_value, msg, file, line)
      integer, intent(in) :: value, expected_value
      character(len=*), intent(in) :: msg
      character(len=*), intent(in) :: file
      integer, intent(in) :: line

      integer ierr

      if (value /= expected_value) then
         write(error_unit, '(a,i0,a,i0,a)')    "Task #", rank, ": " //file // ":", line, ", " // msg
         write(error_unit, '(a,i0,a,i0,a,i0)') "Task #", rank, "  assertion failed, expected ", expected_value, ", but got ", value
         call mpi_abort(MPI_COMM_WORLD, 1, ierr)
         stop 1
      endif
#ifdef DEBUG
      write(error_unit, '(a,i0,a,i0,a,i0)') "Task #", rank, "  assertion at " //file// ":", line, " valid, got expected value ", &
        expected_value
#endif
   end subroutine

   subroutine x_assert_true(value, msg, file, line)
      logical, intent(in) :: value
      character(len=*), intent(in) :: msg
      character(len=*), intent(in) :: file
      integer, intent(in) :: line

      integer ierr

      if (.not. value) then
         write(error_unit, '(a,i0,a,i0,a)') "Task #", rank, ": " //file // ":", line, ", " // msg
         write(error_unit, '(a,i0,a,i0)')   "Task #", rank, "  assertion failed"
         call mpi_abort(MPI_COMM_WORLD, 1, ierr)
         stop 1
      endif
#ifdef DEBUG
      write(error_unit, '(a,i0,a,i0,a)') "Task #", rank, "  assertion at " //file// ":", line, " valid"
#endif
   end subroutine

   subroutine x_assert_false(value, msg, file, line)
      logical, intent(in) :: value
      character(len=*), intent(in) :: msg
      character(len=*), intent(in) :: file
      integer, intent(in) :: line

      call x_assert_true(.not. value, msg, file, line)
   end subroutine

end module
