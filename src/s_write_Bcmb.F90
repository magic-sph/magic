!$Id$
!********************************************************************
    SUBROUTINE write_Bcmb(time,b,lm_max,l_max,l_max_cmb,minc, &
                          lm2,n_cmb_sets,cmb_file,n_cmb_file)
!********************************************************************

!    !------------ This is release 2 level 2  --------------!
!    !------------ Created on 1/18/02  by JW. -----------

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  | Each call of this subroutine writes time and the poloidal magnetic|
!  | potential coeffitients b at the CMB up to degree and order        |
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
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+


    IMPLICIT NONE


!-- Input variables:
    INTEGER :: lm_max         ! Leading dimension of b
    INTEGER :: l_max          ! Max degree of b(*,*)
    INTEGER :: l_max_cmb      ! Max degree of output
    INTEGER :: minc           ! Basic wave-number
    INTEGER :: lm2(0:l_max,0:l_max) ! Gives position of (l,m) coeff
    REAL(kind=8) ::  time           ! Time
    COMPLEX(kind=8) :: b(lm_max) ! Poloidal field potential
    CHARACTER(len=*) :: cmb_file ! Name of output file
    INTEGER :: n_cmb_file     ! Output unit for $cmb_file

!-- Output:
    INTEGER :: n_cmb_sets     ! Total no. of cmb sets,
! should be set to zero before first call

!-- Local variables:
    INTEGER :: n,n_out ! counter
    INTEGER :: l,m            ! degree and order
    INTEGER :: lm             ! position of (l,m) in b(*,n_r_cmb)
    INTEGER :: m_max_cmb      ! Max order of output
    INTEGER :: lm_max_cmb     ! Max no of combinations l,m for output
    INTEGER :: n_data         ! No of output data
    INTEGER :: n_r_cmb        ! Position of cmb-radius on grid

    INTEGER, PARAMETER :: n_data_max=40000     ! Maximum no of output data
    REAL(kind=8) ::  out(n_data_max)! Output array


!-- end of declaration
!---------------------------------------------------------------------

            
!--- Definition of max degree for output
    IF ( l_max < l_max_cmb ) l_max_cmb=l_max

!--- Define postition of CMB on radial grid:
    n_r_cmb=1

!--- Calculate no of data for l_max_cmb:
    m_max_cmb=(l_max_cmb/minc)*minc
    lm_max_cmb=m_max_cmb*(l_max_cmb+1)/minc-         &
               m_max_cmb*(m_max_cmb-minc)/(2*minc) + &
               l_max_cmb-m_max_cmb+1
    n_data=2*lm_max_cmb-l_max_cmb-2
    IF ( n_data > n_data_max ) THEN
        WRITE(*,*)
        WRITE(*,*) ' Dimension n_data_max too small'
        WRITE(*,*) ' in subroutine write_b_cmb!'
        WRITE(*,*) ' Should be at least:',n_data
        STOP
    END IF

!--- Increase no. of cmb_sets:
    n_cmb_sets=n_cmb_sets+1

!--- Open output file name:
    OPEN(n_cmb_file,FILE=cmb_file,POSITION='APPEND', &
         FORM='UNFORMATTED')

!--- If this is the first set write l_max_cmb and minc into first line:
    IF ( n_cmb_sets == 1 ) WRITE(n_cmb_file) l_max_cmb,minc,n_data

!--- Write b(*) into output array out(*):
    n_out=0

!--- Axisymmetric part: (m=0) only real part stored
    DO l=1,l_max
        lm=lm2(l,0)
        IF ( l <= l_max_cmb ) THEN
            n_out=n_out+1
            out(n_out)=REAL(b(lm))
        END IF
    END DO

!--- Non-axisymmetric part: store real and imag part
    DO m=minc,l_max_cmb,minc
        DO l=m,l_max
            lm=lm2(l,m)
            IF ( l <= l_max_cmb ) THEN
                n_out=n_out+1
                out(n_out)=REAL(b(lm))
                n_out=n_out+1
                out(n_out)=AIMAG(b(lm))
                IF ( n_out > n_data_max ) THEN
                    WRITE(*,*)
                    WRITE(*,*) ' Dimension n_data_max too small'
                    WRITE(*,*) ' in subroutine write_b_cmb!'
                    WRITE(*,*) ' Should be at least:',n_data
                    STOP
                END IF
            END IF
        END DO
    END DO


!--- Finally write output array out(*) into cmb_file:
    WRITE(n_cmb_file) time,(out(n),n=1,n_out)

!--- Close cmb_file
    CLOSE(n_cmb_file)


    RETURN
    end SUBROUTINE write_Bcmb

!----------------------------------------------------------------------
