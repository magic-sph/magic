!$Id$
module charmanip

  implicit none

  contains
!------------------------------------------------------------------------
    subroutine capitalize(string)
!   Convert lower-case letters into capital letters
!------------------------------------------------------------------------

    implicit none

    character(len=*) :: string

    character(len=26),parameter :: LOWER_CASE='abcdefghijklmnopqrstuvwxyz'
    character(len=26),parameter :: UPPER_CASE='ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    integer :: i, n

    DO i = 1, len(string)
    ! -- Find location of letter in lower case constant string
        n = index( LOWER_CASE, string( i:i ) )
    ! -- If current substring is a lower case letter,
    !    make it upper case
        IF ( n /= 0 ) string( i:i ) = UPPER_CASE( n:n )
    END DO

    return
    end subroutine capitalize
!------------------------------------------------------------------------
    subroutine delete_string(string,string_del,length)
!   Deletes string_del from string and returns new length of string.
!------------------------------------------------------------------------

    implicit none

    character(len=*) :: string
    character(len=*) :: string_del
    integer :: length,length_del

    integer :: n,i

    integer :: pos

    length=len(trim(string))
    length_del=len(string_del)

    do n=1,length
    
        pos=index(trim(string),string_del)
    
        if ( pos == 0 ) return
    
        string(pos:length-length_del) = &
        string(pos+length_del:length)
        length=length-length_del
                   
        do i=length+1,length+length_del
            string(i:i)=' '
        end do
    
    end do

    return
    end subroutine delete_string
!------------------------------------------------------------------------
    subroutine str2dble(string,num)
!   interprets next word in string as an dble real number
!   deletes leading blanks and next_word from string
!------------------------------------------------------------------------

    implicit none

    character(len=*) :: string  ! input
    real(kind=8) :: num            ! output

    read(string,*) num

    open(99,file='.helpfile',status='unknown')
    write(99,*) num
    close(99)

    return
    end subroutine str2dble
!------------------------------------------------------------------------
    integer function length_to_blank(string)
!   determines number of characters before first blank in string
!------------------------------------------------------------------------

    implicit none

    character(len=*) :: string

    integer :: i

    length_to_blank=0
    do i=1,len(string)
        if( string(i:i) == ' ' ) goto 1
    end do

    1 length_to_blank=i-1

    return
    end function length_to_blank
!------------------------------------------------------------------------
    integer function length_to_char(string,char)
!------------------------------------------------------------------------

    implicit none

    character(len=*) :: string
    character(len=1) :: char

    integer :: i

    length_to_char=0
    do i=1,len(string)
        if( string(i:i) == char ) goto 1
    end do

    length_to_char=-1 ! char not found !
    return

    1 length_to_char=i-1
    return

    end function length_to_char
!------------------------------------------------------------------------
    subroutine dble2str(num, str)
!  converts a dble number num into a character str

    implicit none

    real(kind=8) :: num
    character(len=72) :: work
    character(len=*) :: str
    integer :: i

    write(work, '(F20.12)') num
    write(str, '(A8)') trim(adjustl(work))
    i = index(str,'.')
    str(i:i)='_'

    return

    end subroutine dble2str
!------------------------------------------------------------------------
end module charmanip
