!$Id: s_write_Bcmb.F90 430 2013-02-19 14:29:48Z gastine $
!********************************************************************
MODULE write_special
  IMPLICIT NONE

CONTAINS
  SUBROUTINE write_Bcmb(time,b,llm,ulm,l_max,l_max_cmb,minc, &
       &                lm2,n_cmb_sets,cmb_file,n_cmb_file)
    !********************************************************************

    !    !------------ This is release 2 level 2  --------------!
    !    !------------ Created on 1/18/02  by JW. -----------

    !  +-------------+----------------+------------------------------------+
    !  |                                                                   |
    !  | Each call of this subroutine writes time and the poloidal magnetic|
    !  | potential coefficients b at the CMB up to degree and order        |
    !  | l_max_cmb at the end of output file $cmb_file.                    |
    !  | The parameters l_max_cmb, minc and the number of stored coeffs.   |
    !  | are written into the first line of $cmb_file.                     |
    !  | Each further set contains:                                        |
    !  |          time,real(b(l=0,m=0)),imag(b(l=0,m=0)),                  |
    !  |               real(b(l=1,m=0)),imag(b(l=1,m=0)),                  |
    !  |                    ...................                            |
    !  |                                                                   |
    !  | Real and imaginary part of b(*) for all orders m<=l are written   |
    !  | for a specific degree l, then for the degrees l+1, l+2, l_max_cmb.|
    !  | Logical l_save_cmb can be used to control save output:            |
    !  |     l_save_cmb=.true. : Outputfile $cmb_file is opened and closed |
    !  |                         at each call in this subroutine.          |
    !  |     l_save_cmb=.false.: Outputfile must be open before this       |
    !  |                         subroutine is called.                     |
    !  |                         Normally, the $cmb_file is opened at the  |
    !  |                         beginning of the run and closed after the |
    !  |                         run is finished.                          |
    !  |                                                                   |
    !  +-------------------------------------------------------------------+

    IMPLICIT NONE

    !-- Input variables:
    !INTEGER,INTENT(IN) :: lm_max         ! Leading dimension of b
    INTEGER,INTENT(IN) :: llm,ulm
    INTEGER,intent(IN) :: l_max          ! Max degree of b(*,*)
    INTEGER,intent(INOUT) :: l_max_cmb      ! Max degree of output
    INTEGER,intent(IN) :: minc           ! Basic wave-number
    INTEGER,intent(IN) :: lm2(0:l_max,0:l_max) ! Gives position of (l,m) coeff
    REAL(kind=8),intent(IN) ::  time           ! Time
    COMPLEX(kind=8),intent(IN) :: b(llm:ulm) ! Poloidal field potential
    CHARACTER(len=*),intent(IN):: cmb_file ! Name of output file
    INTEGER,intent(IN) :: n_cmb_file     ! Output unit for $cmb_file

    !-- Output:
    INTEGER,intent(INOUT) :: n_cmb_sets     ! Total no. of cmb sets,
    ! should be set to zero before first call

    !-- Local variables:
    INTEGER :: n,n_out ! counter
    INTEGER :: l,m            ! degree and order
    INTEGER :: lm             ! position of (l,m) in b(*,n_r_cmb)
    INTEGER :: m_max_cmb      ! Max order of output
    INTEGER :: lm_max_cmb     ! Max no of combinations l,m for output
    INTEGER :: n_data         ! No of output data
    INTEGER :: n_r_cmb        ! Position of cmb-radius on grid

    !INTEGER, PARAMETER :: n_data_max=40000     ! Maximum no of output data
    REAL(kind=8),DIMENSION(:),allocatable ::  out ! Output array


    !-- end of declaration
    !---------------------------------------------------------------------
    !WRITE(*,"(2A)") "write_Bcmb: ",TRIM(cmb_file)


    !--- Definition of max degree for output
    IF ( l_max < l_max_cmb ) l_max_cmb=l_max

    !--- Define postition of CMB on radial grid:
    n_r_cmb=1

    !--- Calculate no of data for l_max_cmb:
    m_max_cmb=(l_max_cmb/minc)*minc
    lm_max_cmb= m_max_cmb*(l_max_cmb+1)/minc &
         &     -m_max_cmb*(m_max_cmb-minc)/(2*minc) &
         &     +l_max_cmb-m_max_cmb+1
    n_data=2*lm_max_cmb-l_max_cmb-2
    ALLOCATE(out(n_data))
    !IF ( n_data > n_data_max ) THEN
    !   WRITE(*,*)
    !   WRITE(*,*) ' Dimension n_data_max too small'
    !   WRITE(*,*) ' in subroutine write_b_cmb!'
    !   WRITE(*,*) ' Should be at least:',n_data
    !   STOP
    !END IF

    !--- Increase no. of cmb_sets:
    n_cmb_sets=n_cmb_sets+1

    !--- Open output file name:
    OPEN(n_cmb_file,FILE=cmb_file,POSITION='APPEND', &
         FORM='UNFORMATTED')

    !--- If this is the first set write l_max_cmb and minc into file:
    !WRITE(*,"(A,I3)") "write_Bcmb: n_cmb_sets = ",n_cmb_sets
    IF ( n_cmb_sets == 1 ) WRITE(n_cmb_file) l_max_cmb,minc,n_data

    !--- Write b(*) into output array out(*):
    n_out=0

    !--- Axisymmetric part: (m=0) only real part stored
    DO l=1,l_max_cmb
       lm=lm2(l,0)
       n_out=n_out+1
       out(n_out)=REAL(b(lm))
    END DO

    !--- Non-axisymmetric part: store real and imag part
    DO m=minc,l_max_cmb,minc
       DO l=m,l_max_cmb
          lm=lm2(l,m)
          n_out=n_out+1
          out(n_out)=REAL(b(lm))
          n_out=n_out+1
          out(n_out)=AIMAG(b(lm))
          IF ( n_out > n_data ) THEN
             WRITE(*,*)
             WRITE(*,*) ' n_out larger than n_data'
             WRITE(*,*) ' in subroutine write_b_cmb!'
             WRITE(*,*) ' Should not happen!'
             STOP
          END IF
       END DO
    END DO

    !--- Finally write output array out(*) into cmb_file:
    WRITE(n_cmb_file) time,(out(n),n=1,n_out)

    !--- Close cmb_file
    CLOSE(n_cmb_file)
    deallocate(out)

    RETURN
  end SUBROUTINE write_Bcmb
END MODULE write_special
!----------------------------------------------------------------------
