!$Id$
!*************************************************************************
SUBROUTINE graphOut(time,n_r,FORMAT,vr,vt,vp,br,bt,bp,sr, &
     &              n_theta_start,n_theta_block_size,lGraphHeader)
  !*************************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/17/02  by JW. --------------!

  !-------------------------------------------------------------------------

  !  called in amhd

  !  output of components of velocity, magnetic field vector and
  !  entropy for graphics. Version May, 23, 2000.

  !  n_r: (input) for n_r = 0 a header is written
  !               for n_r > 0 values at radial level n_r are written

  !  format : (input) = 0: unformatted
  !                   = 1: formatted writing
  !                   =-1: comment lines are included into file for
  !                        easier reading (cannot be used for graphics
  !                        processing in this form)

  !  vr...sr: (input) arrays with grid-point values

  !  n_theta_start : (input) values are written for theta-points :
  !                  n_theta_start <= n_theta <= n_theta_start-1+n_theta_block

  !-------------------------------------------------------------------------

  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE blocking
  USE horizontal_data
  USE logic
  USE output_data

  IMPLICIT NONE

  REAL(kind=8) :: time

  INTEGER :: n_r                          ! radial grod point no.
  INTEGER :: format                       ! determins format
  INTEGER :: n_theta_start                ! start theta no.
  INTEGER :: n_theta_stop                 ! end theta no.
  INTEGER :: n_theta_block_size           ! size of theta block
  LOGICAL :: lGraphHeader

  REAL(kind=8) :: vr(nrp,*),vt(nrp,*),vp(nrp,*)
  REAL(kind=8) :: br(nrp,*),bt(nrp,*),bp(nrp,*)
  REAL(kind=8) :: sr(nrp,*)

  !-- Local:
  INTEGER :: n_phi   ! counter for longitude
  INTEGER :: n_theta ! counter for colatitude
  INTEGER :: n_theta_cal                  ! position of block colat in all colats

  INTEGER :: precision

  REAL(kind=8) :: fac,fac_r
  REAL(kind=4) :: dummy(n_phi_max,nfs)

  CHARACTER(len=20) :: version


  !-- End of Declaration
  !----------------------------------------------------------------------


  !-- Write header & colatitudes for n_r=0:

  IF ( lGraphHeader ) THEN

     IF ( format /= 0 ) THEN

        !----- Formatted output:

        version='Graphout_Version_6'

        !-------- Write parameters:
        WRITE(n_graph_file,'(a)')  version ! identifying header
        ! for IDL magsym5 routine
        WRITE(n_graph_file,'(A64)') runid
        IF ( format < 0 ) write(n_graph_file,'(/ &
             &  "Time,truncation,outputgrid:")')
        WRITE(n_graph_file,'(e12.5,1x,10i6)') &
             time,n_r_max,n_theta_max,n_phi_tot, &
             n_r_ic_max-1,minc,nThetaBs
        IF ( format < 0 ) write(n_graph_file, &
             & '(/"Ra, Ek, Pr, Pm, radratio, sigma_ratio:")')
        WRITE(n_graph_file,'(6e14.6)') &
             ra,ek,pr,prmag,radratio,sigma_ratio

        !-------- Write colatitudes:
        IF ( format < 0 ) THEN
           WRITE(n_graph_file, &
                &   '(/"Colatitudes:", /"point#, theta")')
           DO n_theta=1,n_theta_max/2   ! NHS
              WRITE(n_graph_file,'(i4, f9.4)') &
                   n_theta,theta_ord(n_theta)
           END DO
           WRITE(n_graph_file,'("--- equator ---")')
           DO n_theta=n_theta_max/2+1,n_theta_max   ! SHS
              WRITE(n_graph_file,'(i4, f9.4)') &
                   n_theta,theta_ord(n_theta)
           END DO
        ELSE IF ( format >= 0 ) then
           WRITE(n_graph_file,'(129(1x,f8.5))') &
                (theta_ord(n_theta),n_theta=1,n_theta_max)
        END IF

     ELSE

        !----- Unformatted output:

        version='Graphout_Version_7'

        !-------- Write parameters:
        WRITE(n_graph_file) version
        WRITE(n_graph_file) runid
        WRITE(n_graph_file)            SNGL(time), &
             FLOAT(n_r_max),FLOAT(n_theta_max), &
             FLOAT(n_phi_tot),FLOAT(n_r_ic_max-1), &
             FLOAT(minc),FLOAT(nThetaBs), &
             SNGL(ra),SNGL(ek),SNGL(pr),SNGL(prmag), &
             SNGL(radratio),SNGL(sigma_ratio)

        !-------- Write colatitudes:
        WRITE(n_graph_file) (SNGL(theta_ord(n_theta)), &
             n_theta=1,n_theta_max)

     END IF

     lGraphHeader=.FALSE.

  ELSE  ! Call not for writing header

     !*******************************************************************
     !  define CRITICAL section, so that the write statements do not
     !  get mixed up
     !*******************************************************************

     !-- Determine radius and thetas in this block:
     n_theta_stop=n_theta_start+n_theta_block_size-1
     IF ( format < 0 ) WRITE(n_graph_file, &
          &   '(/"--- Radial level, radius, i1, i2 ---")')

     IF ( format /= 0 ) THEN
        WRITE(n_graph_file,'(I4, F9.5, 2I4)') &
             n_r-1,r(n_r)/r(1), &
             n_theta_start,n_theta_stop
     ELSE
        WRITE(n_graph_file) &
             FLOAT(n_r-1),SNGL(r(n_r)/r(1)), &
             FLOAT(n_theta_start),FLOAT(n_theta_stop)
     END IF


     !-- Write entropy:
     precision=1
     IF ( format < 0 ) WRITE(n_graph_file,'(/" S: ")')
     DO n_theta=1,n_theta_block_size,2
        DO n_phi=1,n_phi_max ! do loop over phis
           dummy(n_phi,n_theta)  =SNGL(sr(n_phi,n_theta))   ! NHS
           dummy(n_phi,n_theta+1)=SNGL(sr(n_phi,n_theta+1)) ! SHS
        END DO
     END DO
     CALL graph_write(n_phi_max,n_theta_block_size,dummy, &
          precision,format,n_graph_file)

     precision=0

     !-- Calculate and write radial velocity:
     IF ( format < 0 ) WRITE(n_graph_file,'(/" Vr: ")')
     fac=or2(n_r)*vScale*orho1(n_r)
     DO n_theta=1,n_theta_block_size,2
        DO n_phi=1,n_phi_max
           dummy(n_phi,n_theta)  =SNGL(fac*vr(n_phi,n_theta))
           dummy(n_phi,n_theta+1)=SNGL(fac*vr(n_phi,n_theta+1))
        END DO
     END DO
     CALL graph_write(n_phi_max,n_theta_block_size,dummy, &
          precision,format,n_graph_file)


     !-- Calculate and write latitudinal velocity:
     IF ( format < 0 ) WRITE(n_graph_file,'(/" Vt: ")')
     fac_r=or1(n_r)*vScale*orho1(n_r)
     DO n_theta=1,n_theta_block_size,2
        n_theta_cal=n_theta_start+n_theta-1
        fac=fac_r*O_sin_theta(n_theta_cal)
        DO n_phi=1,n_phi_max
           dummy(n_phi,n_theta)  =SNGL(fac*vt(n_phi,n_theta))
           dummy(n_phi,n_theta+1)=SNGL(fac*vt(n_phi,n_theta+1))
        END DO
     END DO
     CALL graph_write(n_phi_max,n_theta_block_size,dummy, &
          precision,format,n_graph_file)

     !-- Calculate and write longitudinal velocity:
     IF ( format < 0 ) WRITE(n_graph_file,'(/" Vp: ")')
     fac_r=or1(n_r)*vScale*orho1(n_r)
     DO n_theta=1,n_theta_block_size,2
        n_theta_cal=n_theta_start+n_theta-1
        fac=fac_r*O_sin_theta(n_theta_cal)
        DO n_phi=1,n_phi_max
           dummy(n_phi,n_theta)  =SNGL(fac*vp(n_phi,n_theta))
           dummy(n_phi,n_theta+1)=SNGL(fac*vp(n_phi,n_theta+1))
        END DO
     END DO
     CALL graph_write(n_phi_max,n_theta_block_size,dummy, &
          precision,format,n_graph_file)

     IF ( l_mag ) THEN

        !-- Calculate and write radial magnetic field:
        IF ( format < 0 ) WRITE(n_graph_file,'(/'' Br: '')')
        fac=or2(n_r)
        DO n_theta=1,n_theta_block_size,2
           DO n_phi=1,n_phi_max
              dummy(n_phi,n_theta)  =SNGL(fac*br(n_phi,n_theta))
              dummy(n_phi,n_theta+1)=SNGL(fac*br(n_phi,n_theta+1))
           END DO
        END DO
        CALL graph_write(n_phi_max,n_theta_block_size,dummy, &
             precision,format,n_graph_file)

        !-- Calculate and write latitudinal magnetic field:
        IF ( format < 0 ) WRITE(n_graph_file,'(/'' Bt: '')')
        DO n_theta=1,n_theta_block_size,2
           n_theta_cal=n_theta_start+n_theta-1
           fac=or1(n_r)*O_sin_theta(n_theta_cal)
           DO n_phi=1,n_phi_max
              dummy(n_phi,n_theta)  =SNGL(fac*bt(n_phi,n_theta))
              dummy(n_phi,n_theta+1)=SNGL(fac*bt(n_phi,n_theta+1))
           END DO
        END DO
        CALL graph_write(n_phi_max,n_theta_block_size,dummy, &
             precision,format,n_graph_file)

        !-- Calculate and write longitudinal magnetic field:
        IF ( format < 0 ) WRITE(n_graph_file,'(/'' Bp: '')')
        DO n_theta=1,n_theta_block_size,2
           n_theta_cal=n_theta_start+n_theta-1
           fac=or1(n_r)*O_sin_theta(n_theta_cal)
           DO n_phi=1,n_phi_max
              dummy(n_phi,n_theta)  =SNGL(fac*bp(n_phi,n_theta))
              dummy(n_phi,n_theta+1)=SNGL(fac*bp(n_phi,n_theta+1))
           END DO
        END DO
        CALL graph_write(n_phi_max,n_theta_block_size,dummy, &
             precision,format,n_graph_file)

     END IF ! l_mag ?

     !-- End of CRITICAL section ********************************************

  END IF


  RETURN
end SUBROUTINE graphOut

!-----------------------------------------------------------------------
