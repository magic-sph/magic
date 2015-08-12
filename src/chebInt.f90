!$Id$
module chebInt_mod

   use precision_mod, only: cp
   use chebyshev_polynoms_mod, only: cheb_grid
   use init_costf, only: init_costf1
   use cosine_transform, only: costf1
   use radial_der, only: get_dcheb
   use const, only: two, four, half

   implicit none

   private

   public :: chebIntInit, chebInt, chebIntD

contains 

   subroutine chebIntInit(zMin,zMax,zNorm,nNorm,nGridPointsMax, &
                          z,nGridPoints,i_costf_init,d_costf_init)

      !-- Input variables:
      real(cp), intent(in) ::  zMin,zMax  ! integration interval !
      real(cp), intent(in) ::  zNorm      ! norm interval length
      integer,  intent(in) :: nNorm       ! suggested number of grid points
                                              ! for norm length
                                              ! will be adjusted to nGridPoints
      integer :: nGridPointsMax ! dimension of z on input

      !-- Output variables:
      ! grid points, dimension at >= nGridPointsMax
      real(cp), intent(out) ::  z(nGridPointsMax) 
      !  number of used grid points
      integer,  intent(out) :: nGridPoints
      ! array needed for cheb transforms, dimension >= 2*nGridPointsMax+2
      integer,  intent(out) :: i_costf_init(2*nGridPointsMax+2) 
      ! array needed for cheb transforms, dimension >= 2*nGridPointsMax
      real(cp), intent(out) ::  d_costf_init(2*nGridPointsMax+5)

      !-- Local variables:
      integer, parameter :: nChebMax=722
      integer :: n
      real(cp) :: zCheb(nChebMax)
      integer, parameter :: nGridMin=4*6+1

      !-- Adjust number of z points:
      n=int(real(nNorm-1,kind=cp)/four*((zMax-zMin)/zNorm))
      if ( n < 2 )  then
         n=2
      else if ( n == 7 )  then
         n=8
      else if ( n == 11 ) then
         n=12
      else if ( n == 13 ) then
         n=12
      else if ( n == 14 ) then
         n=15
      else if ( n == 17 ) then
         n=16
      else if ( n == 19 ) then
         n=20
      else if ( n == 21 .or. n == 22 .or. n == 23 ) then
         n=24
      !---    After this we use more course steps:
      else if ( n > 24  .and. n <= 30 )  then
         n=30
      else if ( n > 30  .and. n <= 40 )  then
         n=40
      else if ( n > 40  .and. n <= 50 )  then
         n=50
      else if ( n > 50  .and. n <= 60 )  then
         n=60
      else if ( n > 60  .and. n <= 75 )  then
         n=75
      else if ( n > 75  .and. n <= 90 )  then
         n=90
      else if ( n > 90  .and. n <= 100 ) then
         n=100
      else if ( n > 100 .and. n <= 120 ) then
         n=120
      else if ( n > 120 .and. n <= 135 ) then
         n=135
      else if ( n > 135 .and. n <= 150 ) then
         n=150
      else if ( n > 135 .and. n <= 150 ) then
         n=150
      else if ( n > 150 .and. n <= 180 ) then
         n=180
      else if ( n > 180 .and. n <= 200 ) then
         n=200
      else if ( n > 200 .and. n <= 225 ) then
         n=225
      else if ( n > 225 .and. n <= 250 ) then
         n=250
      else if ( n > 250 ) then
         ! Maximum number of grid points set to 1001:
         write(*,*) 'Sorry, no provision for more than 1001 points!'
         write(*,*) 'Sorry, no provision for more than 1001 points!'
         write(*,*) 'Sorry, no provision for more than 1001 points!'
         n=200
      end if
      nGridPoints=4*n+1
      ! New minimum number of grid points (see above):
      nGridPoints=max(nGridPoints,nGridMin)

      if ( nGridPointsMax < nGridPoints ) then
         write(*,*) '! nGridPointsMax too small in chebIntInit!'
         write(*,*) '! Should be at least:',nGridPoints
         stop
      end if

      if ( nChebMax < nGridPoints ) then
         write(*,*) '! Increase nChebMax in s_chebIntInit.f!'
         write(*,*) '! Should be at least:',nGridPoints
         stop
      end if

      !-- Calculate nGridPoints grid points z(*) in interval [zMin,zMax]
      !   zCheb are the grid points in the cheb interval [-1,1]
      !   These are not really needed for the integration.
      call cheb_grid(zMin,zMax,nGridPoints-1,z,zCheb,0.0_cp,0.0_cp,0.0_cp,0.0_cp)

      !-- Initialize fast cos transform for chebs:
      call init_costf1(nGridPoints,i_costf_init,2*nGridPointsMax+2, &
                       d_costf_init,2*nGridPointsMax+5)

   end subroutine chebIntInit
!------------------------------------------------------------------------------
   real(cp) function chebInt(f,zMin,zMax,nGridPoints,nGridPointsMax, &
                                 i_costf_init,d_costf_init)

      !-- Input variables:
      real(cp), intent(in) ::  f(*)            ! function on grid points
      real(cp), intent(in) ::  zMin,zMax       ! integration boundaries
      integer,  intent(in) :: nGridPoints      ! No of grid points
      integer,  intent(in) :: nGridPointsMax   ! No of max grid points
      integer,  intent(in) :: i_costf_init(2*nGridPointsMax+2) ! help array
      real(cp), intent(in) :: d_costf_init(2*nGridPointsMax+5) ! help array

      !-- Local variables:
      integer, parameter :: nWorkMax=722  ! dimension for work array
      real(cp) :: work(nWorkMax)      ! work array
      real(cp) :: fr(nWorkMax)        ! function in cheb space
      real(cp) :: chebNorm
      integer :: nCheb              ! counter for chebs
      integer :: nGrid              ! counter for grid points

      chebNorm=sqrt(two/real(nGridPoints-1,cp))

      if ( nWorkMax < nGridPoints ) then
         write(*,*) '! Increase nWorkMax in chebInt!'
         write(*,*) '! Should be at least:',nGridPoints!'
         stop
      end if

      !-- Copy function:
      do nGrid=1,nGridPoints
         fr(nGrid)=f(nGrid)
      end do

      !-- Transform to cheb space:
      call costf1(fr,work,i_costf_init,d_costf_init)
      fr(1)          =half*fr(1)
      fr(nGridPoints)=half*fr(nGridPoints)

      !-- Integration:
      chebInt=0.0_cp
      do nCheb=1,nGridPoints,2  ! only even chebs contribute
         chebInt=chebInt - (zMax-zMin)/real(nCheb*(nCheb-2),cp)*fr(nCheb)
      end do

      !-- Normalize with interval:
      chebInt=chebNorm*chebInt/(zMax-zMin)

   end function chebInt
!------------------------------------------------------------------------------
   real(cp) function chebIntD(f,lDeriv,zMin,zMax,nGridPoints, &
                                  nGridPointsMax,i_costf_init,    &
                                  d_costf_init)

      !-- Input variables:
      real(cp), intent(inout) ::  f(*)           ! function on grid points
      logical,  intent(in) :: lDeriv            
      real(cp), intent(in) ::  zMin,zMax         ! integration boundaries
      integer,  intent(in) :: nGridPoints        ! No of grid points
      integer,  intent(in) :: nGridPointsMax     ! No of max grid points

      integer,  intent(in) :: i_costf_init(2*nGridPointsMax+2) ! help array
      real(cp), intent(in) :: d_costf_init(2*nGridPointsMax+5) ! help array

      !-- Local variables:
      integer, parameter :: nWorkMax=722  ! dimension for work array
      real(cp) :: work(nWorkMax)      ! work array
      real(cp) :: chebNorm
      real(cp) :: drFac               ! transform fac from cheb space
      integer :: nCheb              ! counter for chebs

      chebNorm=sqrt(two/real(nGridPoints-1,cp))

      if ( nWorkMax < nGridPoints ) then
         write(*,*) '! Increase nWorkMax in chebIntD!'
         write(*,*) '! Should be at least:',nGridPoints!'
         stop
      end if

      !-- Transform to cheb space:
      call costf1(f,work,i_costf_init,d_costf_init)

      !----- Copy:
      if ( lDeriv ) then
         do nCheb=1,nGridPoints
            work(nCheb)=f(nCheb)
         end do
      end if

      !-- Integration:
      f(1)          =half*f(1)
      f(nGridPoints)=half*f(nGridPoints)
      chebIntD=0.0_cp
      do nCheb=1,nGridPoints,2  ! only even chebs contribute
         chebIntD=chebIntD - (zMax-zMin)/real(nCheb*(nCheb-2),cp)*f(nCheb)
      end do
      !-- Normalize with intervall:
      chebIntD=chebNorm*chebIntD/(zMax-zMin)

      !-- Get derivatives:
      if ( lDeriv ) then
         drFac=two/(zMax-zMin)
         call get_dcheb(work,f,nGridPointsMax,nGridPoints,drFac)
         !-- Transform back to grid space:
         call costf1(f,work,i_costf_init,d_costf_init)
      end if

   end function chebIntD
!------------------------------------------------------------------------------
end module chebInt_mod
