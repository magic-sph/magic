MODULE debugging
  implicit none

  INTERFACE debug_write
     MODULE PROCEDURE debug_write_2D,debug_write_1D
  END INTERFACE

  PRIVATE
  public :: debug_write

CONTAINS
SUBROUTINE debug_write_2D(arr,dim1,dim2,label,timestep,form)
  implicit none
  INTEGER :: dim1, dim2
  COMPLEX(kind=8) :: arr(dim1,dim2)
  CHARACTER(len=*) :: label
  integer :: timestep
  CHARACTER(len=1),optional :: form

  character(len=50) :: filename
  LOGICAL :: write_unformatted,write_exponent
  INTEGER :: i,j

  IF (PRESENT(form)) THEN
     IF ((form.EQ.'U').OR.(form.EQ.'u')) THEN
        write_unformatted=.true.
     ELSEIF (form.EQ.'E') THEN
        write_unformatted=.false.
        write_exponent=.true.
     ELSE
        write_unformatted=.false.
        write_exponent=.false.
     END IF
  ELSE
     write_unformatted=.true.
  END IF

  IF (write_unformatted) THEN
     WRITE(filename,"(A,I4.4,A)") TRIM(label),timestep,".dat"
     OPEN(732,file=TRIM(filename),form="unformatted")
     WRITE(732) arr
     CLOSE(732)
  ELSE
     WRITE(filename,"(A,I4.4,A)") TRIM(label),timestep,".txt"
     OPEN(732,file=TRIM(filename))
     DO j=1,dim2
        DO i=1,dim1
           IF (write_exponent) THEN
              WRITE(732,"(2(I4,F21.18))") EXPONENT(REAL(arr(i,j))),FRACTION(REAL(arr(i,j))),&
                   &EXPONENT(aimag(arr(i,j))),FRACTION(aimag(arr(i,j)))
           ELSE
              WRITE(732,"(2ES22.15)") arr(i,j)
           END IF
        END DO
     END DO
     CLOSE(732)
  END IF
END SUBROUTINE debug_write_2D

SUBROUTINE debug_write_1D(arr,dim1,label,timestep,form)
  implicit none
  INTEGER :: dim1
  COMPLEX(kind=8) :: arr(dim1)
  CHARACTER(len=*) :: label
  integer :: timestep
  CHARACTER(len=1),optional :: form

  character(len=50) :: filename
  LOGICAL :: write_unformatted,write_exponent
  INTEGER :: i

  IF (PRESENT(form)) THEN
     IF ((form.EQ.'U').OR.(form.EQ.'u')) THEN
        write_unformatted=.true.
     ELSEIF (form.EQ.'E') THEN
        write_unformatted=.false.
        write_exponent=.true.
     ELSE
        write_unformatted=.false.
        write_exponent=.false.
     END IF
  ELSE
     write_unformatted=.true.
  END IF

  IF (write_unformatted) THEN
     WRITE(filename,"(A,I4.4,A)") TRIM(label),timestep,".dat"
     OPEN(732,file=TRIM(filename),form="unformatted")
     WRITE(732) arr
     CLOSE(732)
  ELSE
     WRITE(filename,"(A,I4.4,A)") TRIM(label),timestep,".txt"
     OPEN(732,file=TRIM(filename))
     DO i=1,dim1
        IF (write_exponent) THEN
           WRITE(732,"(2(I4,F21.18))") EXPONENT(REAL(arr(i))),FRACTION(REAL(arr(i))),&
                &EXPONENT(AIMAG(arr(i))),FRACTION(AIMAG(arr(i)))
        ELSE
           WRITE(732,"(2ES22.15)") arr(i)
        END IF
     END DO
     CLOSE(732)
  END IF
END SUBROUTINE debug_write_1D
END MODULE debugging
