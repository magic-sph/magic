!***************************************************************************
!  defines grid points and truncation 
!***************************************************************************

module truncation
  implicit none
  INTEGER :: n_r_max       ! number of radial grid points, 
  INTEGER :: n_cheb_max    ! max degree-1 of cheb polynomia
  INTEGER :: n_phi_tot     ! number of longitude grid points
  INTEGER :: n_r_ic_max    ! number of grid points in inner core
  INTEGER :: n_cheb_ic_max ! number of chebs in inner core
  INTEGER :: minc          ! basic wavenumber, longitude symmetry  
  INTEGER :: nalias        ! controlls dealiasing in latitude and 



  !-- Derived quantities:
  INTEGER :: n_phi_max     ! absolute number of phi grid-points
  INTEGER :: n_theta_max   ! number of theta grid-points
  INTEGER :: l_max,m_max   ! max degree and order of Plms
  INTEGER :: n_m_max     ! max number of ms (different oders)
  INTEGER :: lm_max         ! number of l/m combinations
  INTEGER :: lmP_max        ! number of l/m combination if l runs to l_max+1
  INTEGER :: lm_max_real    ! number of l/m combination for real representation (cos/sin)
  INTEGER :: nrp,ncp       ! dimension of phi points in for real/complex arrays
  INTEGER :: n_r_tot        ! total number of radial grid points

  !--- Now quantities for magnetic fields:
  !    Set lMag=0 if you want to save this memory (see c_fields)!
  INTEGER :: lMagMem
  INTEGER :: n_r_maxMag,n_r_ic_maxMag,n_r_totMag,l_maxMag,lm_maxMag

  !-- Values for averaged fields, only necessary 
  !   if average is desirev (lAveMem=0).
  INTEGER :: lAveMem
  INTEGER :: n_r_max_ave,n_r_ic_max_ave,lm_max_ave

  !-- Movie memory control:
  INTEGER :: lMovieMem
  INTEGER, PARAMETER :: n_movie_work=1
  INTEGER :: ldtBMem
  INTEGER :: lm_max_dtB,n_r_max_dtB,n_r_ic_max_dtB
  INTEGER :: lmP_max_dtB

  !--- Memory control for stress output:
  INTEGER :: lStressMem
  INTEGER :: n_r_maxStr,n_theta_maxStr,n_phi_maxStr

  !--- Memory control for Geotrophys output:
  INTEGER, PARAMETER :: lGeos=0
  INTEGER :: n_r_maxGeos,lm_maxGeos,nrpGeos,ncpGeos

  !--------------------------------------------------------------------------
CONTAINS
  SUBROUTINE initialize_truncation

    INTEGER :: n_r_maxML,n_r_ic_maxML,n_r_totML,l_maxML,lm_maxML
    INTEGER :: n_r_max_AL,n_r_ic_max_AL,lm_max_AL
    INTEGER :: lm_max_dL,lmP_max_dL,n_r_max_dL,n_r_ic_max_dL
    INTEGER :: n_r_maxSL,n_theta_maxSL,n_phi_maxSL
    INTEGER :: n_r_maxGL,lm_maxGL,nrpGL

    ! absolute number of phi grid-points
    n_phi_max=n_phi_tot/minc

    ! number of theta grid-points
    n_theta_max=n_phi_tot/2

    ! max degree and order of Plms
    l_max=(nalias*n_theta_max)/30 
    m_max=(l_max/minc)*minc

    ! max number of ms (different oders)
    n_m_max=m_max/minc+1

    ! number of l/m combinations
    lm_max=m_max*(l_max+1)/minc - m_max*(m_max-minc)/(2*minc)+(l_max+1-m_max)
    ! number of l/m combination if l runs to l_max+1
    lmP_max=lm_max+n_m_max

    ! number of l/m combination 
    ! for real representation (cos/sin)
    lm_max_real=2*lm_max

    nrp=n_phi_max+2
    ncp=nrp/2

    ! total number of radial grid points
    n_r_tot=n_r_max+n_r_ic_max


    !--- Now quantities for magnetic fields:
    !    Set lMag=0 if you want to save this memory (see c_fields)!
    n_r_maxML     = lMagMem*n_r_max
    n_r_ic_maxML  = lMagMem*n_r_ic_max
    n_r_totML     = lMagMem*n_r_tot
    l_maxML       = lMagMem*l_max
    lm_maxML      = lMagMem*lm_max
    n_r_maxMag    = MAX(1,n_r_maxML)
    n_r_ic_maxMag = MAX(1,n_r_ic_maxML)
    n_r_totMag    = MAX(1,n_r_totML)
    l_maxMag      = MAX(1,l_maxML)
    lm_maxMag     = MAX(1,lm_maxML)

    !-- Values for averaged fields, only necessary 
    !   if average is desired (lAveMem=0).
    n_r_max_AL   =lAveMem*n_r_max
    n_r_ic_max_AL=lAveMem*n_r_ic_max
    lm_max_AL    =lAveMem*lm_max
    n_r_max_ave   =MAX(n_r_max_AL,1)
    n_r_ic_max_ave=MAX(n_r_ic_max_AL,1)
    lm_max_ave    =MAX(lm_max_AL,1)

    !-- Movie memory control:
    lm_max_dL    =ldtBMem*lm_max
    lmP_max_dL   =ldtBMem*lmP_max
    n_r_max_dL   =ldtBMem*n_r_max
    n_r_ic_max_dL=ldtBMem*n_r_ic_max
    lm_max_dtB    =MAX(lm_max_DL,1) 
    lmP_max_dtB   =MAX(lmP_max_DL,1)
    n_r_max_dtB   =MAX(n_r_max_DL,1)
    n_r_ic_max_dtB=MAX(n_r_ic_max_DL,1)

    !--- Memory control for stress output:
    n_r_maxSL     =lStressMem*n_r_max
    n_theta_maxSL =lStressMem*n_theta_max
    n_phi_maxSL   =lStressMem*n_phi_max
    n_r_maxStr    =MAX(n_r_maxSL,1)
    n_theta_maxStr=MAX(n_theta_maxSL,1)
    n_phi_maxStr  =MAX(n_phi_maxSL,1)

    !--- Memory control for Geotrophys output:
    n_r_maxGL  =lGeos*n_r_max
    lm_maxGL   =lGeos*lm_max
    nrpGL      =lGeos*nrp
    n_r_maxGeos=MAX(n_r_maxGL,2)
    lm_maxGeos =MAX(lm_maxGL,2)
    nrpGeos    =MAX(nrpGL,2)
    ncpGeos    =nrpGeos/2 
  END SUBROUTINE initialize_truncation

!***********************************************************************
  SUBROUTINE checkTruncation
!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  | Checks truncations defined in truncation.f and writes then        |
!  | into STDOUT and the log-file.                                     |
!  | MPI: called only by the processor responsible for output !        |
!  |                                                                   |
!  +-------------------------------------------------------------------+

    IMPLICIT NONE

!-----------------------------------------------------------------------

    if ( minc < 1 ) then
        write(*,*)
        write(*,*) '! Wave number minc should be > 0!'
        stop
    end if
    if ( mod(n_phi_tot,minc) /= 0 ) then
        write(*,*)
        write(*,*) '! Number of longitude grid points n_phi_tot'
        write(*,*) '! must be a multiple of minc !!'
        stop
    end if
    if ( mod(n_phi_max,4) /= 0 ) then
        write(*,*)
        write(*,*) '! Number of longitude grid points n_phi_tot/minc'
        write(*,*) '! must be a multiple of 4!!'
        stop
    end if
    if ( mod(n_phi_tot,16) /= 0 ) then
        write(*,*)
        write(*,*) '! Number of longitude grid points n_phi_tot'
        write(*,*) '! must be a multiple of 16!!'
        stop
    end if
    if ( n_phi_max/2 <= 2 ) then
        write(*,*)
        write(*,*) '! Number of longitude grid points n_phi_tot'
        write(*,*) '! must be larger than 2*minc!!',n_phi_max
        stop
    end if

!-- Checking radial grid:
    if ( mod(n_r_max-1,4) /= 0 ) then
        write(*,*)
        write(*,*) '! Number n_r_max-1 should be a multiple '
        write(*,*) '! of 4 !!'
        stop
    end if

    if ( n_theta_max <= 2 ) then
        write(*,*)
        write(*,*) '! Number of latitude grid points n_theta_max'
        write(*,*) '! must be larger than 2!!'
        stop
    end if
    if ( mod(n_theta_max,4) /= 0 ) then
        write(*,*)
        write(*,*) &
        '! Number n_theta_max of latitude grid points be a multiple '
        write(*,*) '! must be a multiple of 4 !!'
        stop
    end if
    if ( n_cheb_max > n_r_max ) then
        write(*,*)
        write(*,*) &
        '! n_cheb_max should be <= n_r_max!'
        stop
    end if
    if ( n_cheb_max < 1 ) then
        write(*,*)
        write(*,*) '! n_cheb_max should be > 1!'
        stop
    end if
    if ( (n_phi_max+1)*n_theta_max > lm_max_real*(n_r_max+2) ) then
        write(*,*)
        write(*,*) &
        '! (n_phi_max+1)*n_theta_max > lm_max_real*(n_r_max+2) !'
        stop '12'
    end if
    IF ( n_r_tot > 2*n_r_max ) THEN
        WRITE(*,*) 'Increase dimension of rhs1 and rhs2 to n_r_tot=', &
        n_r_tot
        STOP
    END IF


    RETURN
    end SUBROUTINE checkTruncation
!---------------------------------------------------------------------
end module truncation
