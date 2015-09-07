!$Id$
module out_coeff
  
   use precision_mod
   use logic, only: l_save_out

   implicit none

   private

   public :: write_Bcmb, write_coeff_r

contains
!----------------------------------------------------------------------
   subroutine write_Bcmb(time,b,llm,ulm,l_max,l_max_cmb,minc, &
        &                lm2,n_cmb_sets,cmb_file,n_cmb_file)
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

      !-- Input variables:
      integer,          intent(in) :: llm,ulm
      integer,          intent(in) :: l_max               ! Max degree of b(*,*)
      integer,          intent(in) :: minc                ! Basic wave-number
      integer,          intent(in) :: lm2(0:l_max,0:l_max)! Gives position of (l,m) coeff
      real(cp),         intent(in) ::  time               ! Time
      complex(cp),      intent(in) :: b(llm:ulm)          ! Poloidal field potential
      character(len=*), intent(in):: cmb_file             ! Name of output file
      integer,          intent(in) :: n_cmb_file          ! Output unit for $cmb_file

      !-- Output variables:
      integer, intent(inout) :: l_max_cmb      ! Max degree of output
      integer, intent(inout) :: n_cmb_sets     ! Total no. of cmb sets,
      ! should be set to zero before first call

      !-- Local variables:
      integer :: n,n_out ! counter
      integer :: l,m            ! degree and order
      integer :: lm             ! position of (l,m) in b(*,n_r_cmb)
      integer :: m_max_cmb      ! Max order of output
      integer :: lm_max_cmb     ! Max no of combinations l,m for output
      integer :: n_data         ! No of output data
      integer :: n_r_cmb        ! Position of cmb-radius on grid

      real(cp), allocatable ::  out(:) ! Output array

      !--- Definition of max degree for output
      if ( l_max < l_max_cmb ) l_max_cmb=l_max

      !--- Define postition of CMB on radial grid:
      n_r_cmb=1

      !--- Calculate no of data for l_max_cmb:
      m_max_cmb=(l_max_cmb/minc)*minc
      lm_max_cmb= m_max_cmb*(l_max_cmb+1)/minc &
           &     -m_max_cmb*(m_max_cmb-minc)/(2*minc) &
           &     +l_max_cmb-m_max_cmb+1
      n_data=2*lm_max_cmb-l_max_cmb-2

      allocate(out(n_data))

      !--- Increase no. of cmb_sets:
      n_cmb_sets=n_cmb_sets+1

      !--- Open output file name:
      if ( l_save_out .or. n_cmb_sets == 0 ) then
         open(n_cmb_file, file=cmb_file, position='append', form='unformatted')
      end if

      !--- If this is the first set write l_max_cmb and minc into file:
      if ( n_cmb_sets <= 1 ) write(n_cmb_file) l_max_cmb,minc,n_data

      !--- Write b(*) into output array out(*):
      n_out=0

      !--- Axisymmetric part: (m=0) only real part stored
      do l=1,l_max_cmb
         lm=lm2(l,0)
         n_out=n_out+1
         out(n_out)=real(b(lm))
      end do

      !--- Non-axisymmetric part: store real and imag part
      do m=minc,l_max_cmb,minc
         do l=m,l_max_cmb
            lm=lm2(l,m)
            n_out=n_out+1
            out(n_out)=real(b(lm))
            n_out=n_out+1
            out(n_out)=aimag(b(lm))
            if ( n_out > n_data ) then
               write(*,*)
               write(*,*) ' n_out larger than n_data'
               write(*,*) ' in subroutine write_b_cmb!'
               write(*,*) ' Should not happen!'
               stop
            end if
         end do
      end do

      !--- Finally write output array out(*) into cmb_file:
      write(n_cmb_file) time,(out(n),n=1,n_out)

      !--- Close cmb_file
      if ( l_save_out .or. n_cmb_sets == 0 ) then
         close(n_cmb_file)
      end if

      deallocate(out)

   end subroutine write_Bcmb
!----------------------------------------------------------------------
   subroutine write_coeff_r(time,w,dw,ddw,z,r,          &
      &                     llm,ulm,l_max,l_max_r,minc, &
      &                     lm2,n_sets,file,n_file,nVBS)
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

      !-- Input variables:
      integer,          intent(in) :: llm,ulm    
      integer,          intent(in) :: l_max        ! Max degree of b(*,*)
      integer,          intent(in) :: minc         ! Basic wave-number
      integer,          intent(in) :: lm2(0:l_max,0:l_max)
      real(cp),         intent(in) ::  r           ! radius of coeffs
      real(cp),         intent(in) ::  time        ! Time
      complex(cp),      intent(in) :: w(llm:ulm)   ! Poloidal field potential
      complex(cp),      intent(in) :: dw(llm:ulm)  ! dr of Poloidal field potential
      complex(cp),      intent(in) :: ddw(llm:ulm) ! dr^2 of Poloidal field potential
      complex(cp),      intent(in) :: z(llm:ulm)   ! Toroidal field potential
      character(len=*), intent(in) :: file         ! Name of output file
      integer,          intent(in) :: n_file       ! Output unit for $file
      integer,          intent(in) :: nVBS         ! True if output is flow

      !-- Output:
      integer, intent(inout) :: l_max_r     ! Max degree of output
      integer, intent(inout) :: n_sets      ! Total no. of cmb sets,
      ! should be set to zero before first call

      !-- Local variables:
      integer :: n,n_out        ! counter
      integer :: l,m            ! degree and order
      integer :: lm             ! position of (l,m) in v(*),..
      integer :: m_max_r        ! Max order of output
      integer :: lm_max_r       ! Max no of combinations l,m for output
      integer :: n_data         ! No of output data

      real(cp), allocatable ::  out(:)! Output array

      !--- Definition of max degree for output
      if ( l_max < l_max_r ) l_max_r=l_max

      !--- Calculate no of data for l_max_r:
      m_max_r=(l_max_r/minc)*minc
      lm_max_r=m_max_r*(l_max_r+1)/minc- &
           m_max_r*(m_max_r-minc)/(2*minc) + &
           l_max_r-m_max_r+1
      n_data=2*lm_max_r-l_max_r-2
      !--- JW 10.Apr.2014: corrected dimension check for different output:
      if ( nVBS == 1 ) then
         allocate(out(3*n_data))
      else if ( nVBS == 2 ) then
         allocate(out(4*n_data))
      else if ( nVBS == 3 ) then
         allocate(out(n_data+1))
      end if

      !--- Increase no. of sets:
      n_sets=n_sets+1

      !--- Open output file with name $file:
      if ( l_save_out ) then
         open(n_file, file=file, form='unformatted', status='unknown', &
              position='append')
      end if

      !--- If this is the first set write, l_max_r and minc into first line:
      if ( n_sets == 1 ) then
         write(n_file) l_max_r,minc,n_data,r
      end if

      !--- Write b(*) into output array out(*):
      n_out=0

      if ( nVBS == 3 ) then
         !--- Axisymmetric part of s: (m=0) only real part stored
         do l=0,l_max ! start with l=0
            lm=lm2(l,0)
            if ( l <= l_max_r ) then
               n_out=n_out+1
               out(n_out)=real(w(lm))
            end if
         end do
      else
         !--- Axisymmetric part of w: (m=0) only real part stored
         do l=1,l_max ! start with l=1
            lm=lm2(l,0)
            if ( l <= l_max_r ) then
               n_out=n_out+1
               out(n_out)=real(w(lm))
            end if
         end do
      end if

      !--- Non-axisymmetric part of w: store real and imag part
      do m=minc,l_max_r,minc
         do l=m,l_max
            lm=lm2(l,m)
            if ( l <= l_max_r ) then
               n_out=n_out+1
               out(n_out)=real(w(lm))
               n_out=n_out+1
               out(n_out)=aimag(w(lm))
            end if
         end do
      end do

      if ( nVBS /= 3 ) then
      !-- Now output for flow or magnetic field only:
         !--- Axisymmetric part of dw: (m=0) only real part stored
         do l=1,l_max
            lm=lm2(l,0)
            if ( l <= l_max_r ) then
               n_out=n_out+1
               out(n_out)=real(dw(lm))
            end if
         end do
         !--- Non-axisymmetric part of dv: store real and imag part
         do m=minc,l_max_r,minc
            do l=m,l_max
               lm=lm2(l,m)
               if ( l <= l_max_r ) then
                  n_out=n_out+1
                  out(n_out)=real(dw(lm))
                  n_out=n_out+1
                  out(n_out)=aimag(dw(lm))
               end if
            end do
         end do
         !--- Axisymmetric part of z: (m=0) only real part stored
         do l=1,l_max
            lm=lm2(l,0)
            if ( l <= l_max_r ) then
               n_out=n_out+1
               out(n_out)=real(z(lm))
            end if
         end do
         !--- Non-axisymmetric part of z: store real and imag part
         do m=minc,l_max_r,minc
            do l=m,l_max
               lm=lm2(l,m)
               if ( l <= l_max_r ) then
                  n_out=n_out+1
                  out(n_out)=real(z(lm))
                  n_out=n_out+1
                  out(n_out)=aimag(z(lm))
               end if
            end do
         end do
      end if

      !--- If this is a magnetic field I also store the second radial derivative
      !    of the poloidal potential to caluclate diffusion:
      if ( nVBS == 2 ) then
         !--- Axisymmetric part of ddw: (m=0) only real part stored
         do l=1,l_max
            lm=lm2(l,0)
            if ( l <= l_max_r ) then
               n_out=n_out+1
               out(n_out)=real(ddw(lm))
            end if
         end do
         !--- Non-axisymmetric part of ddw: store real and imag part
         do m=minc,l_max_r,minc
            do l=m,l_max
               lm=lm2(l,m)
               if ( l <= l_max_r ) then
                  n_out=n_out+1
                  out(n_out)=real(ddw(lm))
                  n_out=n_out+1
                  out(n_out)=aimag(ddw(lm))
               end if
            end do
         end do
      end if

      !--- Finally write output array out(*) into file:
      write(n_file) time,(out(n),n=1,n_out)

      !--- Close file
      if ( l_save_out ) then
         close(n_file)
      end if

      deallocate(out)

   end subroutine write_coeff_r
!----------------------------------------------------------------------
end module out_coeff
