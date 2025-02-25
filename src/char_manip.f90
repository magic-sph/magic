module charmanip
   !
   ! This module contains several useful routines to manipule character strings
   !

   use precision_mod

   implicit none

contains

   subroutine capitalize(string)
      !
      !   Convert lower-case letters into capital letters
      !

      !-- Input variables
      character(len=*), intent(inout) :: string

      !-- Local variables
      character :: a
      integer ::  i

      do i = 1, len(string)
         a = string(i:i)
         if ( a >= 'a' .and. a <= 'z' ) string(i:i) = char(ichar(a)-32)
      end do

   end subroutine capitalize
!------------------------------------------------------------------------------
   subroutine delete_string(string,string_del,length)
      !
      !   Deletes string_del from string and returns new length of string.
      !

      !-- Input variables
      character(len=*), intent(inout) :: string
      character(len=*), intent(in) :: string_del

      !-- Output variables
      integer,          intent(out) :: length

      !-- Local variables
      integer :: length_del
      integer :: n,i
      integer :: pos

      length=len_trim(string)
      length_del=len(string_del)

      do n=1,length

         pos=index(trim(string),string_del)

         if ( pos == 0 ) return

         string(pos:length-length_del) = string(pos+length_del:length)
         length=length-length_del

         do i=length+1,length+length_del
            string(i:i)=' '
         end do

      end do

   end subroutine delete_string
!------------------------------------------------------------------------------
   subroutine dble2str(num, str)
      !
      !  converts a dble number num into a character str
      !

      !-- Input variable
      real(cp), intent(in) :: num

      !-- Output variable
      character(len=*), intent(out) :: str

      !-- Local variables
      character(len=72) :: work
      integer :: i

      write(work, '(F20.12)') num
      write(str, '(A8)') trim(adjustl(work))
      i = index(str,'.')
      str(i:i)='_'

   end subroutine dble2str
!------------------------------------------------------------------------
   subroutine write_long_string(prefix, long_string, out_unit)
      !
      ! This subroutine is used to split a long string (with len(str)>85) into
      ! a multi-lines string for a cleaner printout.
      !

      !-- Input variables
      character(len=*), intent(in) :: prefix
      character(len=*), intent(in) :: long_string
      integer,          intent(in) :: out_unit

      !-- Local variables
      integer :: len_long_str, len_prefix, istart, iend, len_per_line
      integer, parameter :: line_width=85
      character(len=8) :: fmtstr

      len_long_str = len(long_string)
      len_prefix = len(prefix)
      len_per_line = line_width-len_prefix
      write(fmtstr,'(A,i2)') 'A', len_prefix-2

      if ( len_prefix+len_long_str <= 85 ) then
         write(out_unit, '(A,A)') prefix, long_string
      else
         istart = 1; iend = len_per_line
         write(out_unit, '(A,A)') prefix, long_string(istart:iend)
         istart = iend+1
         do while (iend < len_long_str-len_per_line)
            iend = istart+len_per_line
            write(out_unit, '(A,'//trim(fmtstr)//',A)') ' !','',long_string(istart:iend)
            istart = iend+1
         end do
         write(out_unit, '(A,'//trim(fmtstr)//',A)') ' !','',long_string(istart:)
      end if

   end subroutine write_long_string
!------------------------------------------------------------------------
end module charmanip
