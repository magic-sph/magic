module potreader_single

   implicit none

   real(kind=4) :: ra, ek, pr, prmag, radratio, sigma_ratio, omega_ma, omega_ic
   real(kind=4) :: time
   integer :: n_r_max, l_max, n_r_ic_max, lm_max, minc
   real(kind=4), allocatable :: radius(:), rho0(:)
   complex(kind=4), allocatable :: pol(:,:), tor(:,:)

contains

   subroutine readPot(filename, endian, l_read_tor)

      !-- Input variables
      character(len=*), intent(in) :: filename
      character(len=1), intent(in) :: endian
      logical,          intent(in) :: l_read_tor

      if ( endian == 'B' ) then
         open(unit=10, file=filename, form='unformatted', convert='big_endian')
      else
         open(unit=10, file=filename, form='unformatted', convert='little_endian')
      end if

      !-- Header
      read(10) l_max, n_r_max, n_r_ic_max, minc, lm_max
      read(10) ra, ek, pr, prmag, radratio, sigma_ratio, omega_ma, omega_ic
      read(10) time

      !-- Radius and density
      if (allocated(radius)) then
         deallocate (radius)
         deallocate (rho0)
      endif
      allocate( radius(n_r_max), rho0(n_r_max) )
      read(10) radius, rho0

      !-- Poloidal potential
      if ( allocated(pol) ) then
         deallocate(pol)
      end if
      allocate( pol(lm_max, n_r_max) )
      read(10) pol

      !-- Toroidal potential
      if ( l_read_tor ) then
         if ( allocated(tor) ) then
            deallocate(tor)
         end if
         allocate( tor(lm_max, n_r_max) )
         read(10) tor
      end if
      close(10)
   end subroutine readPot

end module potreader_single
