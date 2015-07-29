!$Id$
!*************************************************************************
SUBROUTINE store_movie_frame(n_r,vr,vt,vp,br,bt,bp,sr,drSr, &
     &                       dvrdp,dvpdr,dvtdr,dvrdt,cvr,cbr,cbt, &
     &                       n_theta_start,n_theta_block,bCMB)
  !*************************************************************************

  !-------------------------------------------------------------------------
  !  Controls output of movie frames.
  !  Usually called from radialLoop.
  !-------------------------------------------------------------------------

  USE truncation
  USE movie_data, ONLY: frames,n_movie_fields,n_movies,n_movie_surface,n_movie_const,&
       &n_movie_field_type,n_movie_field_start,n_movie_field_stop
  IMPLICIT NONE

  !-- Input:
  INTEGER,INTENT(IN) :: n_r                ! radial grid point no.
  INTEGER,INTENT(IN) :: n_theta_start      ! start theta no.
  INTEGER,INTENT(IN) :: n_theta_block      ! size of theta block

  !----- Field components (block):
  real(kind=8),INTENT(IN) :: vr(nrp,*),vt(nrp,*),vp(nrp,*)
  real(kind=8),INTENT(IN) :: br(nrp,*),bt(nrp,*),bp(nrp,*)
  real(kind=8),INTENT(IN) :: sr(nrp,*),drSr(nrp,*)
  real(kind=8),INTENT(IN) :: dvrdp(nrp,*),dvpdr(nrp,*)
  real(kind=8),INTENT(IN) :: dvtdr(nrp,*),dvrdt(nrp,*)
  real(kind=8),INTENT(IN) :: cvr(nrp,*)
  real(kind=8),INTENT(IN) :: cbr(nrp,*),cbt(nrp,*)

  !----- Poloidal field potential:
  complex(kind=8),INTENT(IN) :: bCMB(lm_max)

  !-- Output: frames stored in array frames(*)


  !-- Local:
  integer :: n_movie        ! No. of movie
  integer :: n_field        ! No. of field
  integer :: n_surface      ! Surface (1=r,2=theta,3=phi)
  integer :: n_const        ! Gives surface
  integer :: n_field_type   ! Numbers field types
  integer :: n_store_last   ! Position i in frame(i) were field starts-1
  integer :: n_theta
  integer :: n_theta_cal
  integer :: n_theta_const
  integer :: n_field_size
  integer :: n_fields
  logical :: lThetaFound


  !-- end of declaration
  !----------------------------------------------------------------------

  do n_movie=1,n_movies

     n_fields =n_movie_fields(n_movie)
     n_surface=n_movie_surface(n_movie)
     n_const  =n_movie_const(n_movie)

     IF ( n_surface == -1 ) THEN ! Earth Surface

        if ( n_r /= 1 ) cycle  ! not CMB radius

        do n_field=1,n_fields
           n_field_type=n_movie_field_type(n_field,n_movie)
           n_store_last=n_movie_field_start(n_field,n_movie)-1
           IF ( n_store_last >= 0 ) THEN
              CALL store_fields_sur(n_store_last,n_field_type, &
                   &                n_theta_start,n_theta_block,bCMB)
           END IF
        end do

     else if ( n_surface == 0 ) then ! 3d

        do n_field=1,n_fields
           n_field_type=n_movie_field_type(n_field,n_movie)
           n_store_last=n_movie_field_start(n_field,n_movie)-1
           IF ( n_store_last >= 0 ) THEN
              CALL store_fields_3d(vr,vt,vp,br,bt,bp,sr,drSr, &
                   &               dvrdp,dvpdr,dvtdr,dvrdt,cvr,cbr,cbt, &
                   &               n_r,n_store_last,n_field_type, &
                   &               n_theta_start,n_theta_block)
           END IF
        end do

     else if ( n_surface == 1 ) then ! Surface r=constant

        if ( n_r /= n_const ) cycle  ! not desired radius

        DO n_field=1,n_fields
           n_field_type=n_movie_field_type(n_field,n_movie)
           n_store_last=n_movie_field_start(n_field,n_movie)-1
           IF ( n_store_last >= 0 ) THEN
              CALL store_fields_r(vr,vt,vp,br,bt,bp,sr,drSr, &
                   &              dvrdp,dvpdr,dvtdr,dvrdt,cvr, &
                   &              n_r,n_store_last,n_field_type, &
                   &              n_theta_start,n_theta_block)
           END IF
        END DO

     else if ( n_surface == 2 ) then ! Surface theta=constant

        !------ Test whether n_theta_movie is in the current theta block
        !       and find its position n_theta_movie_c:
        lThetaFound=.FALSE.
        DO n_theta=1,n_theta_block
           n_theta_cal = n_theta_start + n_theta - 1
           IF ( n_theta_cal == n_const ) THEN
              n_theta_const=n_theta
              lThetaFound=.TRUE.
              exit
           END IF
        END DO
        IF ( .NOT. lThetaFound) CYCLE        ! Theta not found !

        do n_field=1,n_fields
           n_field_type=n_movie_field_type(n_field,n_movie)
           n_store_last=n_movie_field_start(n_field,n_movie)-1
           IF ( n_store_last >= 0 ) THEN
              CALL store_fields_t(vr,vt,vp,br,bt,bp,sr,drSr, &
                   &              dvrdp,dvpdr,dvtdr,dvrdt,cvr,cbt, &
                   &              n_r,n_store_last,n_field_type, &
                   &              n_const,n_theta_const)
             END IF
        end do

     else if ( IABS(n_surface) == 3 ) then  ! Surface phi=const.

        do n_field=1,n_fields
           n_field_type=n_movie_field_type(n_field,n_movie)
           n_store_last=n_movie_field_start(n_field,n_movie)-1
           n_field_size=(n_movie_field_stop(n_field,n_movie) - &
                n_store_last)/2
           IF ( n_store_last >= 0 ) THEN
              CALL store_fields_p(vr,vt,vp,br,bp,bt,sr,drSr, &
                   &              dvrdp,dvpdr,dvtdr,dvrdt,cvr,cbr,cbt, &
                   &              n_r,n_store_last,n_field_type, &
                   &              n_const,n_field_size, &
                   &              n_theta_start,n_theta_block)
           END IF
        end do  ! Do loop over field for one movie


     END IF  ! Surface ?

  end do  ! Do loop over movies !


  return
end subroutine store_movie_frame


!--- END OF SUBROUTINE store_movie_frame
!-----------------------------------------------------------------------


!***********************************************************************
subroutine store_fields_sur(n_store_last,n_field_type, &
     &                      n_theta_start,n_theta_block,bCMB)
  !***********************************************************************

  !  +-------------+----------------+------------------------------------+
  !  |                                                                   |
  !  |  Purpose of this subroutine is to store movie frames for          |
  !  |  surfaces r=const. into array frame(*,*)                          |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+

  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE blocking
  USE horizontal_data
  USE movie_data
  USE output_data

  IMPLICIT NONE

  !--- Input:
  INTEGER,INTENT(IN) :: n_store_last     ! Start position for storing -1
  integer,INTENT(IN) :: n_field_type     ! Defines field type
  integer,INTENT(IN) :: n_theta_start    ! Beginning of theta block
  integer,INTENT(IN) :: n_theta_block    ! Size of theta block
  complex(kind=8),INTENT(IN) :: bCMB(lm_max)


  !--- Local:
  integer :: n_theta
  integer :: n_theta_b
  integer :: n_theta_cal
  integer :: n_phi
  integer :: n_o


  !----- Magnetic field at surface (theta-blocks):
  real(kind=8) :: br_sur(nrp,nfs) ! Radial magnetic field in (phi,theta)-space
  real(kind=8) :: bt_sur(nrp,nfs) ! Latitudinal magnetic field
  real(kind=8) :: bp_sur(nrp,nfs) ! Azimuthal magnetic field.


  !-- End of declaration
  !----------------------------------------------------------------------


  CALL get_B_surface(br_sur,bt_sur,bp_sur,bCMB,n_theta_start,n_theta_block)

  if ( n_field_type == 1 ) then

     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_store_last+(n_theta-1)*n_phi_max
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=br_sur(n_phi,n_theta_b)
        end do
     end do

  else if ( n_field_type == 2 ) then

     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_store_last+(n_theta-1)*n_phi_max
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=bt_sur(n_phi,n_theta_b)
        end do
     end do

  else if ( n_field_type == 3 ) then

     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_store_last+(n_theta-1)*n_phi_max
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=bp_sur(n_phi,n_theta_b)
        end do
     end do

  end if


  return
end subroutine store_fields_sur

!--- End of subroutine store_field_sur
!-----------------------------------------------------------------------


!***********************************************************************
subroutine store_fields_r(vr,vt,vp,br,bt,bp,sr,drSr, &
     dvrdp,dvpdr,dvtdr,dvrdt,cvr, &
     n_r,n_store_last,n_field_type, &
     n_theta_start,n_theta_block)
  !***********************************************************************

  !  +-------------+----------------+------------------------------------+
  !  |                                                                   |
  !  |  Purpose of this subroutine is to store movie frames for          |
  !  |  surfaces r=const. into array frame(*,*)                          |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+
  !  |  ruler                                                            |
  !  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
  !--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+


  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE blocking
  USE horizontal_data
  USE movie_data
  USE output_data

  IMPLICIT NONE

  !--- Input:

  !----- Fields (block):
  real(kind=8),INTENT(IN) :: vr(nrp,*),vt(nrp,*),vp(nrp,*)
  real(kind=8),INTENT(IN) :: br(nrp,*),bt(nrp,*),bp(nrp,*)
  real(kind=8),INTENT(IN) :: sr(nrp,*),drSr(nrp,*)
  real(kind=8),INTENT(IN) :: dvrdp(nrp,*),dvpdr(nrp,*)
  real(kind=8),INTENT(IN) :: dvtdr(nrp,*),dvrdt(nrp,*)
  real(kind=8),INTENT(IN) :: cvr(nrp,*)

  integer,INTENT(IN) :: n_r
  integer,INTENT(IN) :: n_store_last     ! Start position in frame(*)-1
  integer,INTENT(IN) :: n_field_type     ! Defines field type
  integer,INTENT(IN) :: n_theta_start    ! Beginning of theta block
  integer,INTENT(IN) :: n_theta_block    ! Size of theta block

  !--- Local:
  integer :: n_theta
  integer :: n_theta_b
  integer :: n_theta_cal
  integer :: n_phi
  integer :: n_o
  real(kind=8) ::  fac,fac_r


  !-- End of declaration
  !----------------------------------------------------------------------


  !--- Store data for all output thetas in the current block
  !    and all output phis:

  if ( n_field_type == 1 ) then

     fac=or2(n_r)
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_store_last+(n_theta-1)*n_phi_max
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=fac*br(n_phi,n_theta_b)
        end do
     end do

  else if ( n_field_type == 2 ) then

     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_store_last+(n_theta-1)*n_phi_max
        fac=or1(n_r)*O_sin_theta(n_theta_cal)
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=fac*bt(n_phi,n_theta_b)
        end do
     end do

  else if ( n_field_type == 3 ) then

     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_store_last+(n_theta-1)*n_phi_max
        fac=or1(n_r)*O_sin_theta(n_theta_cal)
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=fac*bp(n_phi,n_theta_b)
        end do
     end do

  else if ( n_field_type == 4 ) then

     fac=or2(n_r)*orho1(n_r)*vScale
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_store_last+(n_theta-1)*n_phi_max
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=fac*vr(n_phi,n_theta_b)
        end do
     end do

  else if ( n_field_type == 5 ) then

     fac_r=or1(n_r)*orho1(n_r)*vScale
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_store_last+(n_theta-1)*n_phi_max
        fac=fac_r*O_sin_theta(n_theta_cal)
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=fac*vt(n_phi,n_theta_b)
        end do
     end do

  else if ( n_field_type == 6 ) then

     fac_r=or1(n_r)*orho1(n_r)*vScale
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_store_last+(n_theta-1)*n_phi_max
        fac=fac_r*O_sin_theta(n_theta_cal)
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=fac*vp(n_phi,n_theta_b)
        end do
     end do

  else if ( n_field_type == 7 ) then

     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_store_last+(n_theta-1)*n_phi_max
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=sr(n_phi,n_theta_b)
        end do
     end do

  else if ( n_field_type == 16 ) then

     fac=or1(n_r)*orho1(n_r)*vScale
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_store_last+(n_theta-1)*n_phi_max
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=fac * ( &
                cosTheta(n_theta_cal)*or1(n_r)*cvr(n_phi,n_theta_b) - &
                or2(n_r)*dvrdp(n_phi,n_theta_b) + &
                dvpdr(n_phi,n_theta_b) - &
                beta(n_r)*vp(n_phi,n_theta_b)   )
        end do
     end do

  else if ( n_field_type == 17 ) then

     fac=-or2(n_r)*orho1(n_r)*vScale
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_store_last+(n_theta-1)*n_phi_max
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=fac * &
                vr(n_phi,n_theta_b)*drSr(n_phi,n_theta_b)
        end do
     end do

  else if ( n_field_type == 91 ) then

     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_store_last+(n_theta-1)*n_phi_max
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=drSr(n_phi,n_theta_b)
        end do
     end do

  else if ( n_field_type == 18 ) then

     !--- Helicity:
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_store_last+(n_theta-1)*n_phi_max
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)= &
                or4(n_r)*orho2(n_r)*vr(n_phi,n_theta_b) * &
                cvr(n_phi,n_theta_b) + &
                or2(n_r)*orho2(n_r)*O_sin_theta_E2(n_theta_cal)* ( &
                vt(n_phi,n_theta_b) * &
                ( or2(n_r)*dvrdp(n_phi,n_theta_b) - &
                dvpdr(n_phi,n_theta_b) + &
                beta(n_r)*vp(n_phi,n_theta_b)  ) + &
                vp(n_phi,n_theta_b) * &
                (          dvtdr(n_phi,n_theta_b) - &
                beta(n_r)*vt(n_phi,n_theta_b)    - &
                or2(n_r)*dvrdt(n_phi,n_theta_b) ) )
        end do
     end do

  else if ( n_field_type == 48 ) then

     !--- Radial component of vorticity:
     fac=or2(n_r)*orho1(n_r)*vScale
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_store_last+(n_theta-1)*n_phi_max
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=fac*cvr(n_phi,n_theta_b)
        end do
     end do

  end if


  return
end subroutine store_fields_r

!--- End of subroutine store_fields_r
!-----------------------------------------------------------------------




!***********************************************************************
subroutine store_fields_p(vr,vt,vp,br,bp,bt,sr,drSr, &
     &                    dvrdp,dvpdr,dvtdr,dvrdt,cvr,cbr,cbt, &
     &                    n_r,n_store_last,n_field_type, &
     &                    n_phi_const,n_field_size, &
     &                    n_theta_start,n_theta_block)
  !***********************************************************************

  !  +-------------+----------------+------------------------------------+
  !  |                                                                   |
  !  |  Purpose of this subroutine is to store movie frames for          |
  !  |  surfaces phi=const. into array frames(*,*)                          |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+

  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE blocking
  USE horizontal_data
  USE movie_data
  USE output_data

  IMPLICIT NONE

  !--- Input:
  
  !----- Fields (block):
  REAL(kind=8),INTENT(IN) :: vr(nrp,*),vt(nrp,*),vp(nrp,*)
  REAL(kind=8),INTENT(IN) :: br(nrp,*),bt(nrp,*),bp(nrp,*)
  REAL(kind=8),INTENT(IN) :: sr(nrp,*),drSr(nrp,*)
  REAL(kind=8),INTENT(IN) :: dvrdp(nrp,*),dvpdr(nrp,*)
  REAL(kind=8),INTENT(IN) :: dvtdr(nrp,*),dvrdt(nrp,*)
  REAL(kind=8),INTENT(IN) :: cvr(nrp,*)
  REAL(kind=8),INTENT(IN) :: cbr(nrp,*),cbt(nrp,*)

  INTEGER,INTENT(IN) :: n_r              ! No. of radial point
  INTEGER,INTENT(IN) :: n_store_last     ! Start position in frame(*)-1
  INTEGER,INTENT(IN) :: n_field_type     ! Defines field type
  INTEGER,INTENT(IN) :: n_phi_const      ! No. of surface phi
  INTEGER,INTENT(IN) :: n_field_size     ! Size of field
  INTEGER,INTENT(IN) :: n_theta_start    ! Beginning of theta block
  INTEGER,INTENT(IN) :: n_theta_block    ! Size of theta block

  !--- Output: writes fields into frames(*) in c_output.f

  !--- Local:
  integer :: n_phi_0
  integer :: n_phi_180
  integer :: n_theta,n_theta2
  integer :: n_theta_b
  integer :: n_theta_cal
  integer :: n_phi
  integer :: n_0,n_180
  real(kind=8) ::  phi_norm
  real(kind=8) ::  fac,fac_r

  real(kind=8) ::  fl(2)        ! Field for poloidal field lines


  !-- End of declaration
  !----------------------------------------------------------------------

  !--- Get phi no. for left and right halfspheres:
  n_phi_0=n_phi_const
  if ( mod(minc,2) == 1 ) then
     n_phi_180=n_phi_max/2+n_phi_0
  else
     n_phi_180=n_phi_0
  end if
  n_0=n_store_last+(n_r-1)*n_theta_max
  n_180=n_0+n_field_size
  !WRITE(*,"(A,4I5)") "store_movie: n_0,n_180,n_phi_0,n_phi_180=",n_0,n_180,n_phi_0,n_phi_180
  !WRITE(*,"(I3,A,I3,3(A,I5))") n_field_type,": Write on rank ",rank," frames indices ",&
  !     & n_store_last+(nRstart-1)*n_theta_max+1,&
  !     &" to ",n_store_last+nRstop*n_theta_max,", field_size = ",n_field_size
  phi_norm=1.d0/n_phi_max

  if ( n_field_type == 1 ) then

     fac=or2(n_r)
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_b+n_theta_start-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        frames(n_0+n_theta)  =fac*br(n_phi_0,n_theta_b)
        frames(n_180+n_theta)=fac*br(n_phi_180,n_theta_b)
     end do

  else if ( n_field_type == 2 ) then

     fac=or1(n_r)
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_b+n_theta_start-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        frames(n_0+n_theta)=fac*bt(n_phi_0,n_theta_b) * &
             O_sin_theta(n_theta_cal)
        frames(n_180+n_theta)=fac*bt(n_phi_180,n_theta_b) * &
             O_sin_theta(n_theta_cal)
     end do

  else if ( n_field_type == 3 ) then

     fac=or1(n_r)
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_b+n_theta_start-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        frames(n_0+n_theta)=fac*bp(n_phi_0,n_theta_b) * &
             O_sin_theta(n_theta_cal)
        frames(n_180+n_theta)=fac*bp(n_phi_180,n_theta_b) * &
             O_sin_theta(n_theta_cal)
     end do

  else if ( n_field_type == 4 ) then

     fac=or2(n_r)*orho1(n_r)*vScale
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_b+n_theta_start-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        frames(n_0+n_theta)=fac*vr(n_phi_0,n_theta_b)
        frames(n_180+n_theta)=fac*vr(n_phi_180,n_theta_b)
     end do

  else if ( n_field_type == 5 ) then

     fac=or1(n_r)*orho1(n_r)*vScale
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_b+n_theta_start-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        frames(n_0+n_theta)=fac*vt(n_phi_0,n_theta_b) * &
             O_sin_theta(n_theta_cal)
        frames(n_180+n_theta)=fac*vt(n_phi_180,n_theta_b) * &
             O_sin_theta(n_theta_cal)
     end do

  else if ( n_field_type == 6 ) then

     fac=or1(n_r)*orho1(n_r)*vScale
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_b+n_theta_start-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        frames(n_0+n_theta)=fac*vp(n_phi_0,n_theta_b) * &
             O_sin_theta(n_theta_cal)
        frames(n_180+n_theta)=fac*vp(n_phi_180,n_theta_b) * &
             O_sin_theta(n_theta_cal)
     end do

  else if ( n_field_type == 7 ) then

     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_b+n_theta_start-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        frames(n_0+n_theta)=sr(n_phi_0,n_theta_b)
        frames(n_180+n_theta)=sr(n_phi_180,n_theta_b)
     end do

  else if ( n_field_type == 8 ) then

     !--- Field for axisymmetric poloidal field lines:
     do n_theta_b=1,n_theta_block,2
        n_theta_cal=n_theta_b+n_theta_start-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_theta2=n_theta_cal2ord(n_theta_cal+1)
        CALL get_fl(fl,n_r,n_theta_cal,2,.FALSE.)
        !WRITE(*,"(2(A,I5),A)") "FL: Writing into frames(",n_0+n_theta,") and frames(",&
        !     &n_0+n_theta2,")"
        frames(n_0+n_theta) =fl(1)
        frames(n_0+n_theta2)=fl(2)
     end do

  else if ( n_field_type == 9 ) then

     !--- Axisymmetric B_phi:
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_b+n_theta_start-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        fl(1)=0.d0
        do n_phi=1,n_phi_max   ! Average over phis
           fl(1)=fl(1)+bp(n_phi,n_theta_b)
        end do
        frames(n_0+n_theta) =phi_norm*fl(1) * &
             or1(n_r)*O_sin_theta(n_theta_cal)
     end do

  else if ( n_field_type == 10 ) then

     !--- Field for axisymmetric velocity stream lines:
     do n_theta_b=1,n_theta_block,2
        n_theta_cal=n_theta_b+n_theta_start-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_theta2=n_theta_cal2ord(n_theta_cal+1)
        CALL get_sl(fl,n_r,n_theta_cal,2)
        frames(n_0+n_theta) =fl(1)
        frames(n_0+n_theta2)=fl(2)
     end do

  else if ( n_field_type == 11 ) then

     !--- Axisymmetric v_phi:
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_b+n_theta_start-1
        n_theta =n_theta_cal2ord(n_theta_cal)
        fl(1)=0.d0
        do n_phi=1,n_phi_max   ! Average over phis
           fl(1)=fl(1)+orho1(n_r)*vp(n_phi,n_theta_b)
        end do
        frames(n_0+n_theta) =phi_norm*fl(1) * &
             or1(n_r)*O_sin_theta(n_theta_cal)
     end do

  else if ( n_field_type == 12 ) then

     !--- Axisymmetric T:
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_b+n_theta_start-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        fl(1)=0.d0
        do n_phi=1,n_phi_max   ! Average over phis
           fl(1)=fl(1)+sr(n_phi,n_theta_b)
        end do
        frames(n_0+n_theta) =phi_norm*fl(1)
     end do

  else if ( n_field_type == 92 ) then

     !--- Axisymmetric dsdr:
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_b+n_theta_start-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        fl(1)=0.d0
        do n_phi=1,n_phi_max   ! Average over phis
           fl(1)=fl(1)+drSr(n_phi,n_theta_b)
        end do
        frames(n_0+n_theta) =phi_norm*fl(1)
     end do

  else if ( n_field_type == 16 ) then

     fac=or1(n_r)*orho1(n_r)*vScale
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_b+n_theta_start-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        frames(n_0+n_theta)=fac * ( &
             cosTheta(n_theta_cal)*or1(n_r)*cvr(n_phi_0,n_theta_b) - &
             or2(n_r)*dvrdp(n_phi_0,n_theta_b) + &
             dvpdr(n_phi_0,n_theta_b) - &
             beta(n_r)*vp(n_phi_0,n_theta_b)  )
        frames(n_180+n_theta)=fac * ( &
             cosTheta(n_theta_cal)*or1(n_r)*cvr(n_phi_180,n_theta_b) - &
             or2(n_r)*dvrdp(n_phi_180,n_theta_b) + &
             dvpdr(n_phi_180,n_theta_b) - &
             beta(n_r)*vp(n_phi_180,n_theta_b)  )
     end do

  else if ( n_field_type == 17 ) then

     fac=-or2(n_r)*orho1(n_r)*vScale
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_b+n_theta_start-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        frames(n_0+n_theta)=fac * &
             vr(n_phi_0,n_theta)*drSr(n_phi_0,n_theta_b)
        frames(n_180+n_theta)=fac * &
             vr(n_phi_180,n_theta)*drSr(n_phi_180,n_theta_b)
     end do

  else if ( n_field_type == 91 ) then

     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_b+n_theta_start-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        frames(n_0+n_theta)=drSr(n_phi_0,n_theta_b)
        frames(n_180+n_theta)=drSr(n_phi_180,n_theta_b)
     end do

  else if ( n_field_type == 18 ) then

     !--- Helicity:
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_b+n_theta_start-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        frames(n_0+n_theta)= &
             or4(n_r)*orho2(n_r)*vr(n_phi_0,n_theta_b) * &
             cvr(n_phi_0,n_theta_b) + &
             or2(n_r)*orho2(n_r)*O_sin_theta_E2(n_theta_cal)* ( &
             vt(n_phi_0,n_theta_b) * &
             ( or2(n_r)*dvrdp(n_phi_0,n_theta_b) - &
             dvpdr(n_phi_0,n_theta_b) + &
             beta(n_r)*   vp(n_phi_0,n_theta_b) ) + &
             vp(n_phi_0,n_theta_b) * &
             (          dvtdr(n_phi_0,n_theta_b) - &
             beta(n_r)*   vt(n_phi_0,n_theta_b) - &
             or2(n_r)*dvrdt(n_phi_0,n_theta_b) ) )
        frames(n_180+n_theta)= &
             or4(n_r)*orho2(n_r)*vr(n_phi_180,n_theta_b) * &
             cvr(n_phi_180,n_theta_b) + &
             or2(n_r)*orho2(n_r)*O_sin_theta_E2(n_theta_cal)* ( &
             vt(n_phi_180,n_theta_b) * &
             ( or2(n_r)*dvrdp(n_phi_180,n_theta_b) - &
             dvpdr(n_phi_180,n_theta_b) + &
             beta(n_r)*   vp(n_phi_180,n_theta_b) ) + &
             vp(n_phi_180,n_theta_b) * &
             (          dvtdr(n_phi_180,n_theta_b) - &
             beta(n_r)*  vt(n_phi_180,n_theta_b) - &
             or2(n_r)*dvrdt(n_phi_180,n_theta_b) ) )
     end do

  else if ( n_field_type == 19 ) then

     !--- Axisymmetric helicity:
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_b+n_theta_start-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        fl(1)=0.d0
        do n_phi=1,n_phi_max
           fl(1)=fl(1) + &
                or4(n_r)*orho2(n_r)*vr(n_phi,n_theta_b) * &
                cvr(n_phi,n_theta_b) + &
                or2(n_r)*orho2(n_r)*O_sin_theta_E2(n_theta_cal)*( &
                vt(n_phi,n_theta_b) * &
                ( or2(n_r)*dvrdp(n_phi,n_theta_b) - &
                dvpdr(n_phi,n_theta_b) + &
                beta(n_r)*   vp(n_phi,n_theta_b) ) + &
                vp(n_phi,n_theta_b) * &
                (        dvtdr(n_phi,n_theta_b) - &
                beta(n_r)*   vt(n_phi,n_theta_b) - &
                or2(n_r)*dvrdt(n_phi,n_theta_b) ) )
        end do
        frames(n_0+n_theta)=phi_norm*fl(1)
     end do


  else if ( n_field_type == 47 ) then

     !--- Phi component of vorticity:
     fac=vScale*orho1(n_r)*or1(n_r)
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_b+n_theta_start-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        frames(n_0+n_theta)= &
             fac*O_sin_theta(n_theta_cal)* &
             (          dvtdr(n_phi_0,n_theta_b) - &
             beta(n_r)*   vt(n_phi_0,n_theta_b) - &
             or2(n_r)*dvrdt(n_phi_0,n_theta_b) )
        frames(n_180+n_theta)= &
             fac*O_sin_theta(n_theta_cal)* &
             (          dvtdr(n_phi_180,n_theta_b) - &
             beta(n_r)*   vt(n_phi_180,n_theta_b) - &
             or2(n_r)*dvrdt(n_phi_180,n_theta_b) )
     end do

     !--- Phi component of Lorentz-Force:

  else if ( n_field_type == 54 ) then

     fac_r=LFfac*or3(n_r)
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_b+n_theta_start-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        fac=fac_r*O_sin_theta(n_theta_cal)
        frames(n_0+n_theta)= fac * &
             ( cbr(n_phi_0,n_theta_b)*bt(n_phi_0,n_theta_b) - &
             cbt(n_phi_0,n_theta_b)*br(n_phi_0,n_theta_b) )
        frames(n_180+n_theta)= fac * &
             ( cbr(n_phi_180,n_theta_b)*bt(n_phi_180,n_theta_b) - &
             cbt(n_phi_180,n_theta_b)*br(n_phi_180,n_theta_b) )
     end do

  end if


  return
end subroutine store_fields_p

!--- End of subroutine  store_fields_p
!-----------------------------------------------------------------------


!***********************************************************************
subroutine store_fields_t(vr,vt,vp,br,bt,bp,sr,drSr, &
     dvrdp,dvpdr,dvtdr,dvrdt,cvr,cbt, &
     n_r,n_store_last,n_field_type, &
     n_theta_const,n_theta)
  !***********************************************************************

  !  +-------------+----------------+------------------------------------+
  !  |                                                                   |
  !  |  Purpose of this subroutine is to store movie frames for          |
  !  |  surfaces r=const. into array frame(*,*)                          |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+
  !  |  ruler                                                            |
  !  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
  !--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+


  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE blocking
  USE horizontal_data
  USE movie_data
  USE output_data

  IMPLICIT NONE

  !--- Input:

  !----- Fields (block):
  real(kind=8),INTENT(IN) :: vr(nrp,*),vt(nrp,*),vp(nrp,*)
  real(kind=8),INTENT(IN) :: br(nrp,*),bt(nrp,*),bp(nrp,*)
  real(kind=8),INTENT(IN) :: sr(nrp,*),drSr(nrp,*)
  real(kind=8),INTENT(IN) :: dvrdp(nrp,*),dvpdr(nrp,*)
  real(kind=8),INTENT(IN) :: dvtdr(nrp,*),dvrdt(nrp,*)
  real(kind=8),INTENT(IN) :: cvr(nrp,*),cbt(nrp,*)

  integer,INTENT(IN) :: n_r              ! No. of radial grid point
  integer,INTENT(IN) :: n_store_last     ! Position in frame(*)-1
  integer,INTENT(IN) :: n_field_type     ! Defines field
  integer,INTENT(IN) :: n_theta_const    ! No. of theta to be stored
  integer,INTENT(IN) :: n_theta          ! No. of theta in block


  !--- Output: writes fields into frame(*,*) in c_output.f

  !--- Local:
  integer :: n_phi
  integer :: n_o

  real(kind=8) ::  fac


  !-- End of declaration
  !----------------------------------------------------------------------


  n_o=n_store_last+(n_r-1)*n_phi_max

  if ( n_field_type == 1 ) then

     fac=or2(n_r)
     do n_phi=1,n_phi_max
        frames(n_o+n_phi)=fac*br(n_phi,n_theta)
     end do

  else if ( n_field_type == 2 ) then

     fac=or1(n_r)*O_sin_theta(n_theta_const)
     do n_phi=1,n_phi_max
        frames(n_o+n_phi)=fac*bt(n_phi,n_theta)
     end do

  else if ( n_field_type == 3 ) then

     fac=or1(n_r)*O_sin_theta(n_theta_const)
     do n_phi=1,n_phi_max
        frames(n_o+n_phi)=fac*bp(n_phi,n_theta)
     end do

  else if ( n_field_type == 4 ) then

     fac=or2(n_r)*orho1(n_r)*vScale
     do n_phi=1,n_phi_max
        frames(n_o+n_phi)=fac*vr(n_phi,n_theta)
     end do

  else if ( n_field_type == 5 ) then

     fac=or1(n_r)*orho1(n_r)*O_sin_theta(n_theta_const)*vScale
     do n_phi=1,n_phi_max
        frames(n_o+n_phi)=fac*vt(n_phi,n_theta)
     end do

  else if ( n_field_type == 6 ) then

     fac=or1(n_r)*orho1(n_r)*O_sin_theta(n_theta_const)*vScale
     do n_phi=1,n_phi_max
        frames(n_o+n_phi)=fac*vp(n_phi,n_theta)
     end do

  else if ( n_field_type == 7 ) then

     do n_phi=1,n_phi_max
        frames(n_o+n_phi)=sr(n_phi,n_theta)
     end do

  else if ( n_field_type == 13 ) then

     fac=-or1(n_r)*O_sin_theta(n_theta_const)
     do n_phi=1,n_phi_max
        frames(n_o+n_phi)=fac*bt(n_phi,n_theta)
     end do

  else if ( n_field_type == 14 ) then

     fac=-or1(n_r)*O_sin_theta(n_theta_const)
     do n_phi=1,n_phi_max
        frames(n_o+n_phi)=fac*cbt(n_phi,n_theta)
     end do

  else if ( n_field_type == 15 ) then

     fac=or1(n_r)*orho1(n_r)*vScale
     do n_phi=1,n_phi_max
        frames(n_o+n_phi)=fac* ( &
             cosTheta(n_theta_const)*or1(n_r)*vr(n_phi,n_theta) - &
             vt(n_phi,n_theta) )
     end do

  else if ( n_field_type == 16 ) then

     fac=or1(n_r)*orho1(n_r)*vScale
     do n_phi=1,n_phi_max
        frames(n_phi+n_o)=fac * ( &
             cosTheta(n_theta_const)*or1(n_r)*cvr(n_phi,n_theta) - &
             or2(n_r)*dvrdp(n_phi,n_theta) + &
             dvpdr(n_phi,n_theta) - &
             beta(n_r)*   vp(n_phi,n_theta) )
     end do

  else if ( n_field_type == 17 ) then

     fac=-or2(n_r)*orho1(n_r)*vScale
     do n_phi=1,n_phi_max
        frames(n_o+n_phi)=fac * &
             vr(n_phi,n_theta)*drSr(n_phi,n_theta)
     end do

  else if ( n_field_type == 91 ) then

     fac=vScale
     do n_phi=1,n_phi_max
        frames(n_o+n_phi)=drSr(n_phi,n_theta)
     end do

  else if ( n_field_type == 18 ) then

     do n_phi=1,n_phi_max
        frames(n_o+n_phi)= &
             or4(n_r)*orho2(n_r)*vr(n_phi,n_theta) * &
             cvr(n_phi,n_theta) + &
             or2(n_r)*orho2(n_r)*O_sin_theta_E2(n_theta_const)*( &
             vt(n_phi,n_theta) * &
             ( or2(n_r)*dvrdp(n_phi,n_theta) - &
             dvpdr(n_phi,n_theta) + &
             beta(n_r)*   vp(n_phi,n_theta) ) + &
             vp(n_phi,n_theta) * &
             (          dvtdr(n_phi,n_theta) - &
             beta(n_r)*   vt(n_phi,n_theta) - &
             or2(n_r)*dvrdt(n_phi,n_theta) ) )
     end do

  end if


  return
end subroutine store_fields_t

!--- End of subroutine  store_fields_t
!-----------------------------------------------------------------------



!***********************************************************************
subroutine store_fields_3d(vr,vt,vp,br,bt,bp,sr,drSr, &
     dvrdp,dvpdr,dvtdr,dvrdt,cvr, &
     cbr,cbt,n_r,n_store_last,n_field_type, &
     n_theta_start,n_theta_block)
  !***********************************************************************

  !  +-------------+----------------+------------------------------------+
  !  |                                                                   |
  !  |  Purpose of this subroutine is to store movie frames for          |
  !  |  surfaces r=const. into array frame(*,*)                          |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+


  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE blocking
  USE horizontal_data
  USE movie_data
  USE output_data

  IMPLICIT NONE

  !--- Input:

  !----- Fields (block):
  real(kind=8),INTENT(IN) :: vr(nrp,*),vt(nrp,*),vp(nrp,*)
  real(kind=8),INTENT(IN) :: br(nrp,*),bt(nrp,*),bp(nrp,*)
  real(kind=8),INTENT(IN) :: sr(nrp,*),drSr(nrp,*)
  real(kind=8),INTENT(IN) :: dvrdp(nrp,*),dvpdr(nrp,*)
  real(kind=8),INTENT(IN) :: dvtdr(nrp,*),dvrdt(nrp,*)
  real(kind=8),INTENT(IN) :: cvr(nrp,*)
  real(kind=8),INTENT(IN) :: cbr(nrp,*),cbt(nrp,*)

  integer,INTENT(IN) :: n_r              ! No. of radial grid point
  integer,INTENT(IN) :: n_store_last     ! Position in frame(*)-1
  integer,INTENT(IN) :: n_field_type     ! Defines field
  integer,INTENT(IN) :: n_theta_start    ! No. of first theta to block
  integer,INTENT(IN) :: n_theta_block    ! Size of theta block


  !--- Output: writes fields into frame(*,*) in c_output.f

  !--- Local:
  integer :: n_phi
  integer :: n_theta,n_theta_b,n_theta_cal
  integer :: n_o,n_or

  real(kind=8) ::  fac,fac_r


  !-- End of declaration
  !----------------------------------------------------------------------

  n_or=n_store_last+(n_r-1)*n_theta_max*n_phi_max

  if ( n_field_type == 1 ) then

     fac=or2(n_r)
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_or+(n_theta-1)*n_phi_max
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=fac*br(n_phi,n_theta_b)
        end do
     end do

  else if ( n_field_type == 2 ) then

     fac_r=or1(n_r)
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_or+(n_theta-1)*n_phi_max
        fac=fac_r*O_sin_theta(n_theta_cal)
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=fac*bt(n_phi,n_theta_b)
        end do
     end do

  else if ( n_field_type == 3 ) then

     fac_r=or1(n_r)
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_or+(n_theta-1)*n_phi_max
        fac=fac_r*O_sin_theta(n_theta_cal)
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=fac*bp(n_phi,n_theta_b)
        end do
     end do

  else if ( n_field_type == 4 ) then

     fac=or2(n_r)*orho1(n_r)*vScale
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_or+(n_theta-1)*n_phi_max
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=fac*vr(n_phi,n_theta_b)
        end do
     end do

  else if ( n_field_type == 5 ) then

     fac_r=or1(n_r)*orho1(n_r)*vScale
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_or+(n_theta-1)*n_phi_max
        fac=fac_r*O_sin_theta(n_theta_cal)
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=fac*vt(n_phi,n_theta_b)
        end do
     end do

  else if ( n_field_type == 6 ) then

     fac_r=or1(n_r)*orho1(n_r)*vScale
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_or+(n_theta-1)*n_phi_max
        fac=fac_r*O_sin_theta(n_theta_cal)
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=fac*vp(n_phi,n_theta_b)
        end do
     end do

  else if ( n_field_type == 7 ) then

     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_or+(n_theta-1)*n_phi_max
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=sr(n_phi,n_theta_b)
        end do
     end do

  else if ( n_field_type == 15 ) then

     fac=or1(n_r)*orho1(n_r)*vScale
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_or+(n_theta-1)*n_phi_max
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=fac * ( &
                cosTheta(n_theta_cal)*or1(n_r)*vr(n_phi,n_theta_b) - &
                vt(n_phi,n_theta_b) )
        end do
     end do

  else if ( n_field_type == 16 ) then

     !--- Z-component of helicity
     fac=or1(n_r)*orho1(n_r)*vScale
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_or+(n_theta-1)*n_phi_max
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=fac * ( &
                cosTheta(n_theta_cal)*or1(n_r)*cvr(n_phi,n_theta_b) - &
                or2(n_r)*dvrdp(n_phi,n_theta_b) + &
                dvpdr(n_phi,n_theta_b) - &
                beta(n_r)*   vp(n_phi,n_theta_b) )
        end do
     end do

  else if ( n_field_type == 17 ) then

     fac=-or2(n_r)*orho1(n_r)*vScale
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_or+(n_theta-1)*n_phi_max
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=fac * &
                vr(n_phi,n_theta_b)*drSr(n_phi,n_theta_b)
        end do
     end do

  else if ( n_field_type == 91 ) then

     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_or+(n_theta-1)*n_phi_max
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=drSr(n_phi,n_theta_b)
        end do
     end do

  else if ( n_field_type == 18 ) then

     !--- Helicity
     fac=vScale*vScale*or2(n_r)*orho2(n_r)
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_or+(n_theta-1)*n_phi_max
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=fac * ( &
                or2(n_r)*vr(n_phi,n_theta_b) * &
                cvr(n_phi,n_theta_b) + &
                O_sin_theta_E2(n_theta_cal)* ( &
                vt(n_phi,n_theta_b) * &
                ( or2(n_r)*dvrdp(n_phi,n_theta_b) - &
                dvpdr(n_phi,n_theta_b) + &
                beta(n_r)*   vp(n_phi,n_theta_b) ) + &
                vp(n_phi,n_theta_b) * &
                (          dvtdr(n_phi,n_theta_b) - &
                beta(n_r)*   vt(n_phi,n_theta_b) - &
                or2(n_r)*dvrdt(n_phi,n_theta_b) ) ) )
        end do
     end do

  else if ( n_field_type == 47 ) then

     !--- Phi component of vorticity
     fac=or1(n_r)*orho1(n_r)*vScale
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_or+(n_theta-1)*n_phi_max
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=fac * &
                O_sin_theta(n_theta_cal)*vp(n_phi,n_theta_b) * &
                (          dvtdr(n_phi,n_theta_b) - &
                beta(n_r)*   vt(n_phi,n_theta_b) - &
                or2(n_r)*dvrdt(n_phi,n_theta_b) )
        end do
     end do

  else if ( n_field_type == 48 ) then

     !--- Radial component of vorticity:
     fac=or2(n_r)*orho1(n_r)*vScale
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_or+(n_theta-1)*n_phi_max
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=fac*cvr(n_phi,n_theta_b)
        end do
     end do

  else if ( n_field_type == 53 ) then

     !--- Omega effect: Br*dvp/dr
     fac_r=or3(n_r)*orho1(n_r)*vScale
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_or+(n_theta-1)*n_phi_max
        fac=fac_r*O_sin_theta(n_theta_cal)
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)=fac*br(n_phi,n_theta_b) * &
                ( dvpdr(n_phi,n_theta_b) - &
                (beta(n_r)+2.d0*or1(n_r))*vp(n_phi,n_theta_b) )
        end do
     end do

  else if ( n_field_type == 54 ) then

     !--- Phi component of Lorentz-Force:

     fac_r=LFfac*or3(n_r)
     do n_theta_b=1,n_theta_block
        n_theta_cal=n_theta_start+n_theta_b-1
        n_theta=n_theta_cal2ord(n_theta_cal)
        n_o=n_or+(n_theta-1)*n_phi_max
        fac=fac_r*O_sin_theta(n_theta_cal)
        do n_phi=1,n_phi_max
           frames(n_phi+n_o)= fac * &
                ( cbr(n_phi,n_theta_b)*bt(n_phi,n_theta_b) - &
                cbt(n_phi,n_theta_b)*br(n_phi,n_theta_b) )
        end do
     end do

  end if


  return
end subroutine store_fields_3d

!--- End of subroutine  store_fields_3d
!-----------------------------------------------------------------------
