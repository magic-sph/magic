!$Id$
!********************************************************************
SUBROUTINE write_coeff_r(time,w,dw,ddw,z,r, &
     &                   lm_max,l_max,l_max_r,minc, &
     &                   lm2,n_sets,file,n_file,l_save_out,lV)
  !********************************************************************

  !  +-------------+----------------+------------------------------------+
  !  |                                                                   |
  !  | Each call of this subroutine writes time and the poloidal and     |
  !  | toroidal coeffitients w,dw,z at a specific radius up to degree    |
  !  | and order l_max_r at the end of output file $file.                |
  !  | If the input is magnetic field (lV=.FALSE.) ddw is stored as well.|
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
  integer :: lm_max         ! Leading dimension of b
  integer :: l_max          ! Max degree of b(*,*)
  integer :: l_max_r      ! Max degree of output
  integer :: minc           ! Basic wave-number
  INTEGER :: lm2(0:l_max,0:l_max)
  real(kind=8) ::  r              ! radius of coeffs
  real(kind=8) ::  time           ! Time
  complex(kind=8) :: w(lm_max)   ! Poloidal field potential
  complex(kind=8) :: dw(lm_max)  ! dr of Poloidal field potential
  complex(kind=8) :: ddw(lm_max) ! dr^2 of Poloidal field potential
  complex(kind=8) :: z(lm_max)   ! Toroidal field potential
  character(len=*) :: file     ! Name of output file
  integer :: n_file         ! Output unit for $file
  logical :: l_save_out     ! Controlls output
  LOGICAL :: lV             ! True if output is flow

  !-- Output:
  integer :: n_sets         ! Total no. of cmb sets,
  ! should be set to zero before first call

  !-- Local variables:
  integer :: n,n_out        ! counter
  integer :: l,m            ! degree and order
  integer :: lm             ! position of (l,m) in v(*),..
  integer :: m_max_r        ! Max order of output
  integer :: lm_max_r       ! Max no of combinations l,m for output
  integer :: n_data         ! No of output data

  integer, parameter :: n_data_max=50000     ! Maximum no of output data
  real(kind=8) ::  out(n_data_max)! Output array


  !-- end of declaration
  !---------------------------------------------------------------------


  !--- Definition of max degree for output
  if ( l_max < l_max_r ) l_max_r=l_max

  !--- Calculate no of data for l_max_r:
  m_max_r=(l_max_r/minc)*minc
  lm_max_r=m_max_r*(l_max_r+1)/minc- &
       m_max_r*(m_max_r-minc)/(2*minc) + &
       l_max_r-m_max_r+1
  n_data=2*lm_max_r-l_max_r-2
  if ( 4*n_data > n_data_max ) then
     write(*,*)
     write(*,*) ' Dimension n_data_max too small'
     write(*,*) ' in subroutine write_coeff_r !'
     write(*,*) ' Should be at least:',3*n_data
     stop
  end if

  !--- Increase no. of sets:
  n_sets=n_sets+1

  !--- Open output file with name $file:
  if ( l_save_out ) then
     open(n_file,file=file,form='unformatted',status='unknown',position='APPEND')
     !------ Position file after last set:
  end if


  !--- If this is the first set write, l_max_r and minc into first line:
  if ( n_sets == 1 ) then
     write(n_file) l_max_r,minc,n_data,r
  end if


  !--- Write b(*) into output array out(*):

  n_out=0

  !--- Axisymmetric part of w: (m=0) only real part stored
  do l=1,l_max
     lm=lm2(l,0)
     if ( l <= l_max_r ) then
        n_out=n_out+1
        out(n_out)=REAL(w(lm))
     end if
  end do

  !--- Non-axisymmetric part of w: store real and imag part
  do m=minc,l_max_r,minc
     do l=m,l_max
        lm=lm2(l,m)
        if ( l <= l_max_r ) then
           n_out=n_out+1
           out(n_out)=REAL(w(lm))
           n_out=n_out+1
           out(n_out)=AIMAG(w(lm))
        end if
     end do
  end do

  !--- Axisymmetric part of dw: (m=0) only real part stored
  do l=1,l_max
     lm=lm2(l,0)
     if ( l <= l_max_r ) then
        n_out=n_out+1
        out(n_out)=REAL(dw(lm))
     end if
  end do

  !--- Non-axisymmetric part of dv: store real and imag part
  do m=minc,l_max_r,minc
     do l=m,l_max
        lm=lm2(l,m)
        if ( l <= l_max_r ) then
           n_out=n_out+1
           out(n_out)=REAL(dw(lm))
           n_out=n_out+1
           out(n_out)=AIMAG(dw(lm))
        end if
     end do
  end do

  !--- Axisymmetric part of z: (m=0) only real part stored
  do l=1,l_max
     lm=lm2(l,0)
     if ( l <= l_max_r ) then
        n_out=n_out+1
        out(n_out)=REAL(z(lm))
     end if
  end do

  !--- Non-axisymmetric part of z: store real and imag part
  do m=minc,l_max_r,minc
     do l=m,l_max
        lm=lm2(l,m)
        if ( l <= l_max_r ) then
           n_out=n_out+1
           out(n_out)=REAL(z(lm))
           n_out=n_out+1
           out(n_out)=AIMAG(z(lm))
        end if
     end do
  end do

  !--- If this is a magnetic field I also store the second radial derivative
  !    of the poloidal potential to caluclate diffusion:

  IF ( .NOT. lV ) THEN

     !--- Axisymmetric part of ddw: (m=0) only real part stored
     do l=1,l_max
        lm=lm2(l,0)
        if ( l <= l_max_r ) then
           n_out=n_out+1
           out(n_out)=REAL(ddw(lm))
        end if
     end do

     !--- Non-axisymmetric part of ddw: store real and imag part
     do m=minc,l_max_r,minc
        do l=m,l_max
           lm=lm2(l,m)
           if ( l <= l_max_r ) then
              n_out=n_out+1
              out(n_out)=REAL(ddw(lm))
              n_out=n_out+1
              out(n_out)=AIMAG(ddw(lm))
           end if
        end do
     end do

  END IF

  !--- Finally write output array out(*) into file:
  write(n_file) time,(out(n),n=1,n_out)

  !--- Close file
  if ( l_save_out ) close(n_file)


  RETURN
end SUBROUTINE write_coeff_r


!----------------------------------------------------------------------
