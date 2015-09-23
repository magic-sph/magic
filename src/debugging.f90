module debugging

   use precision_mod

   implicit none
 
   private
 
   interface debug_write
      module procedure debug_write_2D,debug_write_1D
   end interface
 
   public :: debug_write
 
contains

   subroutine debug_write_2D(arr,dim1,dim2,label,timestep,form)

      !-- Input variables
      integer,          intent(in) :: dim1, dim2
      complex(cp),      intent(in) :: arr(dim1,dim2)
      character(len=*), intent(in) :: label
      integer,          intent(in) :: timestep
      character(len=1), optional, intent(in) :: form
 
      !-- Local variables
      character(len=50) :: filename
      logical :: write_unformatted,write_exponent
      integer :: i,j
 
      if (present(form)) then
         if ((form == 'U').or.(form == 'u')) then
            write_unformatted=.true.
         elseif (form == 'E') then
            write_unformatted=.false.
            write_exponent=.true.
         else
            write_unformatted=.false.
            write_exponent=.false.
         end if
      else
         write_unformatted=.true.
      end if
 
      if (write_unformatted) then
         write(filename,"(A,I4.4,A)") trim(label),timestep,".dat"
         open(732,file=trim(filename),form="unformatted")
         write(732) arr
         close(732)
      else
         write(filename,"(A,I4.4,A)") trim(label),timestep,".txt"
         open(732,file=trim(filename))
         do j=1,dim2
            do i=1,dim1
               if (write_exponent) then
                  write(732,"(2(I4,F21.18))") exponent(real(arr(i,j))),       &
                       & fraction(real(arr(i,j))), exponent(aimag(arr(i,j))), &
                       & fraction(aimag(arr(i,j)))
               else
                  write(732,"(2ES22.15)") arr(i,j)
               end if
            end do
         end do
         close(732)
      end if

   end subroutine debug_write_2D
!------------------------------------------------------------------------------
   subroutine debug_write_1D(arr,dim1,label,timestep,form)

      !-- Input variables:
      integer,          intent(in) :: dim1
      complex(cp),      intent(in) :: arr(dim1)
      character(len=*), intent(in) :: label
      integer,          intent(in) :: timestep
      character(len=1), optional, intent(in) :: form
 
      !-- Local variables
      character(len=50) :: filename
      logical :: write_unformatted,write_exponent
      integer :: i
 
      if (present(form)) then
         if ((form == 'U').or.(form == 'u')) then
            write_unformatted=.true.
         elseif (form == 'E') then
            write_unformatted=.false.
            write_exponent=.true.
         else
            write_unformatted=.false.
            write_exponent=.false.
         end if
      else
         write_unformatted=.true.
      end if
 
      if (write_unformatted) then
         write(filename,"(A,I4.4,A)") trim(label),timestep,".dat"
         open(732,file=trim(filename),form="unformatted")
         write(732) arr
         close(732)
      else
         write(filename,"(A,I4.4,A)") trim(label),timestep,".txt"
         open(732,file=trim(filename))
         do i=1,dim1
            if (write_exponent) then
               write(732,"(2(I4,F21.18))") exponent(real(arr(i))),      &
                    & fraction(real(arr(i))), exponent(aimag(arr(i))),  &
                    & fraction(aimag(arr(i)))
            else
               write(732,"(2ES22.15)") arr(i)
            end if
         end do
         close(732)
      end if

   end subroutine debug_write_1D
!------------------------------------------------------------------------------
end module debugging
