module potreader_single

   implicit none

   real(kind=4) :: ra, ek, pr, prmag, radratio, sigma_ratio, omega_ma, omega_ic
   real(kind=4) :: raxi, sc
   real(kind=4) :: time
   integer :: n_r_max, l_max, n_r_ic_max, lm_max, minc, version
   real(kind=4), allocatable :: radius(:), rho0(:)
   complex(kind=4), allocatable :: pol(:,:), tor(:,:)
   complex(kind=4), allocatable :: pol_ic(:,:), tor_ic(:,:)

contains

   subroutine readPot(filename, endian, l_read_tor, l_read_ic, ver)

      !-- Input variables
      character(len=*), intent(in) :: filename
      character(len=1), intent(in) :: endian
      logical,          intent(in) :: l_read_tor
      logical,          intent(in) :: l_read_ic
      integer,          intent(in) :: ver

      if ( ver == 0 ) then

         if ( endian == 'B' ) then
            open(unit=10, file=filename, form='unformatted', convert='big_endian')
         else
            open(unit=10, file=filename, form='unformatted', convert='little_endian')
         end if

         version = 0
         raxi = 0
         sc = 0
         !-- Header
         read(10) l_max, n_r_max, n_r_ic_max, minc, lm_max
         read(10) ra, ek, pr, prmag, radratio, sigma_ratio, omega_ma, omega_ic
         read(10) time

      else ! out

         if ( endian == 'B' ) then
            open(unit=10, file=filename, form='unformatted', convert='big_endian', &
            &    access='stream')
         else
            open(unit=10, file=filename, form='unformatted', convert='little_endian',&
            &    access='stream')
         end if

         !-- Header
         read(10) version
         read(10) time
         read(10) ra, pr, raxi, sc, prmag, ek, radratio, sigma_ratio
         read(10) n_r_max, n_r_ic_max, l_max, minc, lm_max
         read(10) omega_ic, omega_ma

      end if

      !-- Radius and density
      if (allocated(radius)) deallocate (radius, rho0)
      allocate( radius(n_r_max), rho0(n_r_max) )
      read(10) radius, rho0

      !-- Poloidal potential
      if ( allocated(pol) ) deallocate(pol)
      allocate( pol(lm_max, n_r_max) )
      read(10) pol

      !-- Toroidal potential
      if ( l_read_tor ) then
         if ( allocated(tor) ) deallocate(tor)
         allocate( tor(lm_max, n_r_max) )
         read(10) tor
      end if

      !-- Inner core
      if ( l_read_ic ) then
         if ( allocated(pol_ic) ) deallocate(pol_ic)
         allocate( pol_ic(lm_max, n_r_ic_max) )
         read(10) pol_ic

         if ( allocated(tor_ic) ) deallocate(tor_ic)
         allocate( tor_ic(lm_max, n_r_ic_max) )
         read(10) tor_ic
      end if

      close(10)

   end subroutine readPot

end module potreader_single
