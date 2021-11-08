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
      character(len=26), parameter :: LOWER_CASE='abcdefghijklmnopqrstuvwxyz'
      character(len=26), parameter :: UPPER_CASE='ABCDEFGHIJKLMNOPQRSTUVWXYZ'

      integer :: i, n

      do i = 1, len(string)
      ! -- Find location of letter in lower case constant string
         n = index( LOWER_CASE, string( i:i ) )
      ! -- If current substring is a lower case letter,
      !    make it upper case
         if ( n /= 0 ) string( i:i ) = UPPER_CASE( n:n )
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

      length=len(trim(string))
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
   subroutine str2dble(string,num)
      !
      !   interprets next word in string as an dble real number
      !   deletes leading blanks and next_word from string
      !

      !-- Input variable:
      character(len=*), intent(in) :: string  ! input

      !-- Output variable:
      real(cp), intent(out) :: num     ! output

      !-- Local variable:
      integer :: fileHandle

      read(string,*) num

      open(newunit=fileHandle, file='.helpfile', status='unknown')
      write(fileHandle,*) num
      close(fileHandle)

   end subroutine str2dble
!------------------------------------------------------------------------------
   integer function length_to_blank(string)
      !
      !   determines number of characters before first blank in string
      !

      !-- Input variable
      character(len=*), intent(in) :: string

      !-- Local variable
      integer :: i

      length_to_blank=0
      do i=1,len(string)
         if( string(i:i) == ' ' ) then
            length_to_blank=i-1
            exit
         end if
      end do

   end function length_to_blank
!------------------------------------------------------------------------------
   integer function length_to_char(string,char)

      !-- Input variables:
      character(len=*), intent(in) :: string
      character(len=1), intent(in) :: char

      !-- Local variables:
      logical :: isDetected
      integer :: i

      isDetected=.false.
      length_to_char=0
      do i=1,len(string)
         if( string(i:i) == char ) then
            length_to_char=i-1
            isDetected=.true.
            exit
         end if
      end do

      if ( .not. isDetected ) length_to_char=-1 ! char not found !

   end function length_to_char
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
         do while (iend <= len_long_str-len_per_line)
            iend = istart+len_per_line
            write(out_unit, '(A,'//trim(fmtstr)//',A)') ' !','',long_string(istart:iend)
            istart = iend+1
         end do
         write(out_unit, '(A,'//trim(fmtstr)//',A)') ' !','',long_string(istart:)
      end if

   end subroutine write_long_string
!------------------------------------------------------------------------
end module charmanip
