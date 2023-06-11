module truncation
   !
   ! This module defines the grid points and the truncation 
   !

   use precision_mod, only: cp
   use logic, only: l_finite_diff, l_cond_ic
   use useful, only: abortRun

   implicit none

   integer :: n_r_max       ! number of radial grid points
   integer :: n_cheb_max    ! max degree-1 of cheb polynomia
   integer :: n_phi_tot     ! number of longitude grid points
   integer :: n_r_ic_max    ! number of grid points in inner core
   integer :: n_cheb_ic_max ! number of chebs in inner core
   integer :: minc          ! basic wavenumber, longitude symmetry  
   integer :: nalias        ! controls dealiasing in latitude
   logical :: l_axi         ! logical for axisymmetric calculations
   character(len=72) :: radial_scheme ! radial scheme (either Cheybev of FD)
   real(cp) :: fd_stretch    ! regular intervals over irregular intervals
   real(cp) :: fd_ratio      ! drMin over drMax (only when FD are used)
   real(cp) :: rcut_l        ! Cut-off radius when degree varies with radius
   integer :: fd_order       ! Finite difference order (for now only 2 and 4 are safe)
   integer :: fd_order_bound ! Finite difference order on the  boundaries
 
   !-- Derived quantities:
   integer :: n_phi_max   ! absolute number of phi grid-points
   integer :: n_theta_max ! number of theta grid-points
   integer :: nlat_padded ! number of theta grid-points with padding included
   integer :: n_theta_axi ! number of theta grid-points (axisymmetric models)
   integer :: l_max       ! max degree of Plms
   integer :: m_min       ! min order of Plms
   integer :: m_max       ! max order of Plms
   integer :: n_m_max     ! max number of ms (different oders)
   integer :: lm_max      ! number of l/m combinations
   integer :: n_r_tot     ! total number of radial grid points
 
   !--- Now quantities for magnetic fields:
   !    Set lMag=0 if you want to save this memory (see c_fields)!
   integer :: lMagMem       ! Memory for magnetic field calculation
   integer :: n_r_maxMag    ! Number of radial points to calculate magnetic field
   integer :: n_r_ic_maxMag ! Number of radial points to calculate IC magnetic field
   integer :: n_r_totMag    ! n_r_maxMag + n_r_ic_maxMag
   integer :: l_maxMag      ! Max. degree for magnetic field calculation
   integer :: lm_maxMag     ! Max. number of l/m combinations for magnetic field calculation
 
   !-- Movie memory control:
   integer :: ldtBMem        ! Memory for movie output
   integer :: lm_max_dtB     ! Number of l/m combinations for movie output
   integer :: n_r_max_dtB    ! Number of radial points for movie output
   integer :: n_r_ic_max_dtB ! Number of IC radial points for movie output
 
   !--- Memory control for stress output:
   integer :: lStressMem     ! Memory for stress output
   integer :: n_r_maxStr     ! Number of radial points for stress output
   integer :: n_theta_maxStr ! Number of theta points for stress output
   integer :: n_phi_maxStr   ! Number of phi points for stress output
 
contains

   subroutine initialize_truncation

      integer :: n_r_maxML,n_r_ic_maxML,n_r_totML,l_maxML,lm_maxML
      integer :: lm_max_dL,n_r_max_dL,n_r_ic_max_dL
      integer :: n_r_maxSL,n_theta_maxSL,n_phi_maxSL, l, m

      if ( .not. l_axi ) then
         if ( l_max == 0 ) then
            ! absolute number of phi grid-points
            n_phi_max=n_phi_tot/minc

            ! number of theta grid-points
            n_theta_max=n_phi_tot/2

            ! max degree and order of Plms
            l_max=(nalias*n_theta_max)/30
         else
            n_theta_max=(30*l_max)/nalias
            n_phi_tot  =2*n_theta_max

            !-- Try to get a better prime decomposition
            call prime_decomposition(n_phi_tot)
            n_phi_max  =n_phi_tot/minc
            n_theta_max=n_phi_tot/2
         end if

         if ( m_max == 0 ) then
            m_max=(l_max/minc)*minc
         end if
      else
         minc       =1
         n_phi_max  =1
         n_phi_tot  =1
         if ( l_max == 0 ) then
            n_theta_max=n_theta_axi
            l_max      =(nalias*n_theta_max)/30
         else
            n_theta_max=(30*l_max)/nalias
            n_theta_axi=n_theta_max
         end if
         m_min=0
         m_max=0
      end if

      if ( m_max > l_max ) m_max=l_max

      ! this will be possibly overwritten when SHTns is used
      nlat_padded=n_theta_max

      ! max number of ms (different oders)
      n_m_max=m_max/minc+1

      ! number of l/m combinations
      lm_max=0
      do m=m_min,m_max,minc
         do l=m,l_max
            lm_max = lm_max+1
         end do
      end do

      ! total number of radial grid points
      n_r_tot = n_r_max
      if ( l_cond_ic ) n_r_tot=n_r_max+n_r_ic_max

      !--- Now quantities for magnetic fields:
      !    Set lMag=0 if you want to save this memory (see c_fields)!
      n_r_maxML     = lMagMem*n_r_max
      n_r_ic_maxML  = lMagMem*n_r_ic_max
      n_r_totML     = lMagMem*n_r_tot
      l_maxML       = lMagMem*l_max
      lm_maxML      = lMagMem*lm_max
      n_r_maxMag    = max(1,n_r_maxML)
      n_r_ic_maxMag = max(1,n_r_ic_maxML)
      n_r_totMag    = max(1,n_r_totML)
      l_maxMag      = max(1,l_maxML)
      lm_maxMag     = max(1,lm_maxML)

      !-- Movie memory control:
      lm_max_dL    =ldtBMem*lm_max
      n_r_max_dL   =ldtBMem*n_r_max
      n_r_ic_max_dL=ldtBMem*n_r_ic_max
      lm_max_dtB    =max(lm_max_DL,1) 
      n_r_max_dtB   =max(n_r_max_DL,1)
      n_r_ic_max_dtB=max(n_r_ic_max_DL,1)

      !--- Memory control for stress output:
      n_r_maxSL     =lStressMem*n_r_max
      n_theta_maxSL =lStressMem*n_theta_max
      n_phi_maxSL   =lStressMem*n_phi_max
      n_r_maxStr    =max(n_r_maxSL,1)
      n_theta_maxStr=max(n_theta_maxSL,1)
      n_phi_maxStr  =max(n_phi_maxSL,1)

   end subroutine initialize_truncation
!--------------------------------------------------------------------------------
   subroutine checkTruncation
      !  This function checks truncations and writes it
      !  into STDOUT and the log-file.                                 
      !  MPI: called only by the processor responsible for output !  

      if ( minc < 1 ) then
         call abortRun('! Wave number minc should be > 0!')
      end if
      if ( mod(n_phi_tot,minc) /= 0 .and. (.not. l_axi) ) then
         call abortRun('! Number of longitude grid points n_phi_tot must be a multiple of minc')
      end if
      if ( mod(n_phi_max,4) /= 0 .and. (.not. l_axi) ) then
         call abortRun('! Number of longitude grid points n_phi_tot/minc must be a multiple of 4')
      end if
      if ( mod(n_phi_tot,16) /= 0 .and. (.not. l_axi) ) then
         call abortRun('! Number of longitude grid points n_phi_tot must be a multiple of 16')
      end if
      if ( n_phi_max/2 <= 2 .and. (.not. l_axi) ) then
         call abortRun('! Number of longitude grid points n_phi_tot must be larger than 2*minc')
      end if

      if ( n_theta_max <= 2 ) then
         call abortRun('! Number of latitude grid points n_theta_max must be larger than 2')
      end if
      if ( mod(n_theta_max,4) /= 0 ) then
         call abortRun('! Number n_theta_max of latitude grid points be a multiple must be a multiple of 4')
      end if
      if ( n_cheb_max > n_r_max ) then
         call abortRun('! n_cheb_max should be <= n_r_max!')
      end if
      if ( n_cheb_max < 1 ) then
         call abortRun('! n_cheb_max should be > 1!')
      end if

   end subroutine checkTruncation
!-----------------------------------------------------------------------------
   subroutine prime_decomposition(nlon)
      !
      ! This subroutine tries to find the closest nlon which allow a
      ! prime decomposition of the form 2**i * 3**j * 4**k
      !

      integer, intent(inout) :: nlon

      !-- Local variables
      real(cp) :: dist, res, dist_min
      integer :: i, j, k, i0, j0, k0

      dist_min = 100.0_cp
      do i=0,12
         do j=0,6
            do k=0,6
               res = real(2**i * 3**j * 5**k, cp)
               dist = real(res-nlon, cp)
               if ( dist >= 0 .and. dist < dist_min ) then
                  i0 = i
                  j0 = j
                  k0 = k
                  dist_min=dist
               end if
            end do
         end do
      end do

      nlon = 2**i0 * 3**j0 * 5**k0

   end subroutine prime_decomposition
!-----------------------------------------------------------------------------
end module truncation
