!$Id$
!********************************************************************
SUBROUTINE rBpSpec(time,Tor,TorIC,fileRoot,lIC,map)
  !********************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/17/02  by JW. --------------!

  !--------------------------------------------------------------------
  !  Called from rank0, map gives the lm order of Tor and TorIC
  !--------------------------------------------------------------------

  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE blocking
  USE horizontal_data
  USE logic
  USE output_data
  USE usefull, ONLY: cc2real
  USE LMmapping,only:mappings
  IMPLICIT NONE

  REAL(kind=8) :: time

  COMPLEX(kind=8),INTENT(IN) :: Tor(lm_max,n_r_max)
  COMPLEX(kind=8),INTENT(IN) :: TorIC(lm_max,n_r_ic_max)
  CHARACTER(len=*),INTENT(IN) :: fileRoot
  LOGICAL :: lIC
  TYPE(mappings),intent(IN) :: map

  !-- Output:
  REAL(kind=8) :: e_t_AS(l_max,n_r_tot)
  REAL(kind=8) :: e_t(l_max,n_r_tot)

  !-- Local:
  CHARACTER(len=72) :: specFile
  INTEGER :: n_r,lm,l,m

  REAL(kind=8) :: fac,rRatio,amp
  REAL(kind=8) :: e_t_temp

  !-- Local function:
  LOGICAl :: lAS

  !-- end of declaration
  !---------------------------------------------------------------------

  fac=0.5D0*eScale/(16.D0*DATAN(1.D0))

  DO n_r=1,n_r_max
     DO l=1,6
        e_t(l,n_r)=0.D0
        e_t_AS(l,n_r) = 0.0D0
     END DO
     DO lm=2,lm_max
        l=map%lm2l(lm)
        IF ( l <= 6 ) THEN
           m=map%lm2m(lm)
           amp=REAL(Tor(lm,n_r))
           e_t_temp=dLh(st_map%lm2(l,m))*cc2real(Tor(lm,n_r),m)
           IF ( ABS(amp)/=0d0 ) THEN
              IF ( m == 0 ) e_t_AS(l,n_r)=fac*amp/ABS(amp)*e_t_temp
           END IF
           e_t(l,n_r)=e_t(l,n_r)+fac*e_t_temp
           !WRITE(*,"(A,4I4,6ES22.12)") "e_t_AS,e_t: ",n_r,lm,l,m,&
           !     & e_t_AS(l,n_r),e_t(l,n_r),e_t_temp,amp,Tor(lm,n_r)
        END IF
     END DO    ! do loop over lms in block
  END DO    ! radial grid points

  !-- Inner core:
  DO n_r=2,n_r_ic_max
     DO l=1,6
        e_t_AS(l,n_r_max-1+n_r)=0.D0
        e_t(l,n_r_max-1+n_r)   =0.D0
     END DO
  END DO
  IF ( lIC .AND. l_cond_ic ) then

     lAS=.true.
     IF ( trim(adjustl(fileRoot)) == 'rBrAdvSpec' ) lAS= .FALSE. 

     DO n_r=2,n_r_ic_max
        rRatio=r_ic(n_r)/r_ic(1)
        DO lm=2,lm_max
           l=map%lm2l(lm)
           IF ( l <= 6 ) THEN
              m=map%lm2m(lm)
              IF ( m /= 0 .OR. lAS ) THEN
                 e_t_temp= dLh(st_map%lm2(l,m))*rRatio**(2*l+2) &
                      &    * cc2real(TorIC(lm,n_r),m)
                 amp=REAL(TorIC(lm,n_r))
                 IF ( ABS(amp)/=0d0 ) THEN
                    IF ( m == 0 ) e_t_AS(l,n_r_max-1+n_r)= &
                         fac*amp/ABS(amp)*e_t_temp
                 END IF
                 e_t(l,n_r_max-1+n_r)=e_t(l,n_r_max-1+n_r) + &
                      fac*e_t_temp
              END IF
              !WRITE(*,"(A,4I4,2ES22.12)") "IC: e_t_AS,e_t: ",&
              !     & n_r_max-1+n_r,lm,l,m,e_t_AS(l,n_r_max-1+n_r),e_t(l,n_r_max-1+n_r)
           END IF
        END DO
     END DO

  END IF

  !-- Output into file:
  !     writing l=0/1/2 magnetic energy
  specFile=trim(adjustl(fileRoot))//'.'//tag
  OPEN(91,FILE=specFile,FORM='UNFORMATTED',STATUS='UNKNOWN', &
       POSITION='APPEND')

  !IF ( nLines == 0 ) WRITE(91) FLOAT(n_r_tot-1),SNGL(radratio)
  !IF ( nLines == 0 ) WRITE(91) &
  !     (SNGL(r(n_r)),n_r=1,n_r_max), (SNGL(r_ic(n_r)),n_r=2,n_r_ic_max)
  WRITE(91) SNGL(time),                          &
       (SNGL(e_t(1,n_r)),n_r=1,n_r_tot-1),    &
       (SNGL(e_t(2,n_r)),n_r=1,n_r_tot-1),    &
       (SNGL(e_t(3,n_r)),n_r=1,n_r_tot-1),    &
       (SNGL(e_t(4,n_r)),n_r=1,n_r_tot-1),    &
       (SNGL(e_t(5,n_r)),n_r=1,n_r_tot-1),    &
       (SNGL(e_t(6,n_r)),n_r=1,n_r_tot-1)
  WRITE(91) SNGL(time),                          &
       (SNGL(e_t_AS(1,n_r)),n_r=1,n_r_tot-1), &
       (SNGL(e_t_AS(2,n_r)),n_r=1,n_r_tot-1), &
       (SNGL(e_t_AS(3,n_r)),n_r=1,n_r_tot-1), &
       (SNGL(e_t_AS(4,n_r)),n_r=1,n_r_tot-1), &
       (SNGL(e_t_AS(5,n_r)),n_r=1,n_r_tot-1), &
       (SNGL(e_t_AS(6,n_r)),n_r=1,n_r_tot-1)

  CLOSE(91)


  RETURN
end subroutine rBpSpec

!----------------------------------------------------------------------------
