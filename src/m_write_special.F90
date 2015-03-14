!$Id: s_write_Bcmb.F90 430 2013-02-19 14:29:48Z gastine $
!********************************************************************
MODULE write_special
  
  USE logic, ONLY: l_save_out

  IMPLICIT NONE

CONTAINS
  SUBROUTINE write_Bcmb(time,b,llm,ulm,l_max,l_max_cmb,minc, &
       &                lm2,n_cmb_sets,cmb_file,n_cmb_file)
    !********************************************************************

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

    !--- Increase no. of cmb_sets:
    n_cmb_sets=n_cmb_sets+1

    !--- Open output file name:
    IF ( l_save_out ) OPEN(n_cmb_file,FILE=cmb_file,POSITION='APPEND', &
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
    IF ( l_save_out ) CLOSE(n_cmb_file)

    deallocate(out)

    RETURN
  end SUBROUTINE write_Bcmb
!----------------------------------------------------------------------
  SUBROUTINE write_coeff_r(time,w,dw,ddw,z,r,          &
     &                     llm,ulm,l_max,l_max_r,minc, &
     &                     lm2,n_sets,file,n_file,nVBS)
  !********************************************************************

  !  +-------------+----------------+------------------------------------+
  !  |                                                                   |
  !  | Each call of this subroutine writes time and the poloidal and     |
  !  | toroidal coeffitients w,dw,z at a specific radius up to degree    |
  !  | and order l_max_r at the end of output file $file.                |
  !  | If the input is magnetic field (nVBS=2) ddw is stored as well.    |
  !  | If the input is entropy field (nVBS=3) only ddw=s is stored.      |
  !  | The parameters l_max_r, minc, the number of stored coeffs and     |
  !  | radius in the outer core are written into the first line of $file.|
  !  | Each further set contains:                                        |
  !  |          time,real(w(l=0,m=0)),imag(w(l=0,m=0)),                  |
  !  |               real(w(l=1,m=0)),imag(w(l=1,m=0)),                  |
  !  |                    ...................                            |
  !  |               real(dw(l=0,m=0)),imag(dw(l=0,m=0)),                |
  !  |               real(dw(l=1,m=0)),imag(dw(l=1,m=0)),                |
  !  |                    ...................                            |
  !  |               real(z(l=0,m=0)),imag(z(l=0,m=0)),                  |
  !  |               real(z(l=1,m=0)),imag(z(l=1,m=0)),                  |
  !  |                    ...................                            |
  !  |               real(ddw(l=0,m=0)),imag(ddw(l=0,m=0)),              |
  !  |               real(ddw(l=1,m=0)),imag(ddw(l=1,m=0)),              |
  !  |                                                                   |
  !  | Real and imaginary part of w(*) for all orders m<=l are written   |
  !  | for a specific degree l, then for the degrees l+1, l+2, l_max_r.  |
  !  | Logical l_save_r can be used to control save output:              |
  !  |     l_save_out=.true. : Outputfile $file is opened and closed     |
  !  |                         at each call in this subroutine.          |
  !  |     l_save_out=.false.: Outputfile must be open before this       |
  !  |                         subroutine is called.                     |
  !  |                         Normally, the $file is opened at the      |
  !  |                         beginning of the run and closed after the |
  !  |                         run is finished.                          |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+

    IMPLICIT NONE

    !-- Input variables:
    INTEGER,INTENT(in) :: llm,ulm    
    INTEGER,INTENT(in) :: l_max          ! Max degree of b(*,*)
    INTEGER,INTENT(inout) :: l_max_r     ! Max degree of output
    INTEGER,INTENT(in) :: minc           ! Basic wave-number
    INTEGER,INTENT(in) :: lm2(0:l_max,0:l_max)
    REAL(KIND=8),INTENT(in) ::  r              ! radius of coeffs
    REAL(KIND=8),INTENT(in) ::  time           ! Time
    COMPLEX(KIND=8),INTENT(in) :: w(llm:ulm)   ! Poloidal field potential
    COMPLEX(KIND=8),INTENT(in) :: dw(llm:ulm)  ! dr of Poloidal field potential
    COMPLEX(KIND=8),INTENT(in) :: ddw(llm:ulm) ! dr^2 of Poloidal field potential
    COMPLEX(KIND=8),INTENT(in) :: z(llm:ulm)   ! Toroidal field potential
    CHARACTER(LEN=*),INTENT(in) :: file     ! Name of output file
    INTEGER,INTENT(in) :: n_file         ! Output unit for $file
    INTEGER,INTENT(in) :: nVBS           ! True if output is flow

    !-- Output:
    INTEGER,INTENT(inout) :: n_sets         ! Total no. of cmb sets,
    ! should be set to zero before first call

    !-- Local variables:
    INTEGER :: n,n_out        ! counter
    INTEGER :: l,m            ! degree and order
    INTEGER :: lm             ! position of (l,m) in v(*),..
    INTEGER :: m_max_r        ! Max order of output
    INTEGER :: lm_max_r       ! Max no of combinations l,m for output
    INTEGER :: n_data         ! No of output data

    REAL(KIND=8),ALLOCATABLE,DIMENSION(:) ::  out! Output array

    !-- end of declaration
    !---------------------------------------------------------------------


    !--- Definition of max degree for output
    IF ( l_max < l_max_r ) l_max_r=l_max

    !--- Calculate no of data for l_max_r:
    m_max_r=(l_max_r/minc)*minc
    lm_max_r=m_max_r*(l_max_r+1)/minc- &
         m_max_r*(m_max_r-minc)/(2*minc) + &
         l_max_r-m_max_r+1
    n_data=2*lm_max_r-l_max_r-2
    !--- JW 10.Apr.2014: corrected dimension check for different output:
    IF ( nVBS.EQ.1 ) THEN
       ALLOCATE(out(3*n_data))
    ELSE IF ( nVBS.EQ.2 ) THEN
       ALLOCATE(out(4*n_data))
    ELSE IF ( nVBS.EQ.3 ) THEN
       ALLOCATE(out(n_data+1))
    END IF

    !--- Increase no. of sets:
    n_sets=n_sets+1

    !--- Open output file with name $file:
    IF ( l_save_out ) THEN
        OPEN(n_file,FILE=file,FORM='unformatted',STATUS='unknown',POSITION='APPEND')
    END IF

    !--- If this is the first set write, l_max_r and minc into first line:
    IF ( n_sets == 1 ) THEN
       WRITE(n_file) l_max_r,minc,n_data,r
    END IF


    !--- Write b(*) into output array out(*):

    n_out=0

    IF ( nVBS.EQ.3 ) THEN

       !--- Axisymmetric part of s: (m=0) only real part stored
       DO l=0,l_max ! start with l=0
          lm=lm2(l,0)
          IF ( l <= l_max_r ) THEN
             n_out=n_out+1
             out(n_out)=REAL(w(lm))
          END IF
       END DO

    ELSE

       !--- Axisymmetric part of w: (m=0) only real part stored
       DO l=1,l_max ! start with l=1
          lm=lm2(l,0)
          IF ( l <= l_max_r ) THEN
             n_out=n_out+1
             out(n_out)=REAL(w(lm))
          END IF
       END DO

    END IF

    !--- Non-axisymmetric part of w: store real and imag part
    DO m=minc,l_max_r,minc
       DO l=m,l_max
          lm=lm2(l,m)
          IF ( l <= l_max_r ) THEN
             n_out=n_out+1
             out(n_out)=REAL(w(lm))
             n_out=n_out+1
             out(n_out)=AIMAG(w(lm))
          END IF
       END DO
    END DO


    IF ( nVBS.NE.3 ) THEN
    !-- Now output for flow or magnetic field only:

       !--- Axisymmetric part of dw: (m=0) only real part stored
       DO l=1,l_max
          lm=lm2(l,0)
          IF ( l <= l_max_r ) THEN
             n_out=n_out+1
             out(n_out)=REAL(dw(lm))
          END IF
       END DO

       !--- Non-axisymmetric part of dv: store real and imag part
       DO m=minc,l_max_r,minc
          DO l=m,l_max
             lm=lm2(l,m)
             IF ( l <= l_max_r ) THEN
                n_out=n_out+1
                out(n_out)=REAL(dw(lm))
                n_out=n_out+1
                out(n_out)=AIMAG(dw(lm))
             END IF
          END DO
       END DO

       !--- Axisymmetric part of z: (m=0) only real part stored
       DO l=1,l_max
          lm=lm2(l,0)
          IF ( l <= l_max_r ) THEN
             n_out=n_out+1
             out(n_out)=REAL(z(lm))
          END IF
       END DO

       !--- Non-axisymmetric part of z: store real and imag part
       DO m=minc,l_max_r,minc
          DO l=m,l_max
             lm=lm2(l,m)
             IF ( l <= l_max_r ) THEN
                n_out=n_out+1
                out(n_out)=REAL(z(lm))
                n_out=n_out+1
                out(n_out)=AIMAG(z(lm))
             END IF
          END DO
       END DO

    END IF

    !--- If this is a magnetic field I also store the second radial derivative
    !    of the poloidal potential to caluclate diffusion:

    IF ( nVBS.EQ.2 ) THEN

       !--- Axisymmetric part of ddw: (m=0) only real part stored
       DO l=1,l_max
          lm=lm2(l,0)
          IF ( l <= l_max_r ) THEN
             n_out=n_out+1
             out(n_out)=REAL(ddw(lm))
          END IF
       END DO

       !--- Non-axisymmetric part of ddw: store real and imag part
       DO m=minc,l_max_r,minc
          DO l=m,l_max
             lm=lm2(l,m)
             IF ( l <= l_max_r ) then
                n_out=n_out+1
                out(n_out)=REAL(ddw(lm))
                n_out=n_out+1
                out(n_out)=AIMAG(ddw(lm))
             END IF
          END DO
       END DO

    END IF

    !--- Finally write output array out(*) into file:
    WRITE(n_file) time,(out(n),n=1,n_out)

    !--- Close file
    IF ( l_save_out ) CLOSE(n_file)

    DEALLOCATE(out)

    RETURN
  end SUBROUTINE write_coeff_r

END MODULE write_special
!----------------------------------------------------------------------
