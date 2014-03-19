!$Id$
!********************************************************************
SUBROUTINE storePot(time,b,aj,b_ic,aj_ic, &
     &              nPotSets,root,omega_ma,omega_ic)
  !********************************************************************

  USE truncation
  USE radial_functions
  USE physical_parameters
  USE horizontal_data
  USE logic
  USE output_data
  USE charmanip, only: length_to_char

  IMPLICIT NONE

  !-- input: field via common blocks
  ! include 'truncation.f'
  ! include 'c_radial.f'
  ! include 'c_output.f'
  ! include 'c_phys_param.f'
  ! include 'c_horizontal.f'
  ! include 'c_logic.f'

  REAL(kind=8),INTENT(IN) :: time

  !-- Input of fields to be stored:
  COMPLEX(kind=8),INTENT(IN) :: b(lm_max,n_r_max)
  COMPLEX(kind=8),INTENT(IN) :: aj(lm_max,n_r_max)
  COMPLEX(kind=8),INTENT(IN) :: b_ic(lm_max,n_r_ic_max)
  COMPLEX(kind=8),INTENT(IN) :: aj_ic(lm_max,n_r_ic_max)
  INTEGER,INTENT(INOUT) :: nPotSets
  CHARACTER(80),INTENT(IN) :: root
  REAL(kind=8),INTENT(IN) :: omega_ma,omega_ic

  !-- Work arrays:
  COMPLEX(kind=8) :: workA(lm_max,n_r_max)
  COMPLEX(kind=8) :: workB(lm_max,n_r_max)
  COMPLEX(kind=8) :: workC(lm_max,n_r_max)

  CHARACTER(80) :: string
  INTEGER :: n_r,lm,n_cheb
  INTEGER :: lengthR
  CHARACTER(80) :: fileName
  LOGICAL :: lVB
  REAL(kind=8) :: chebNorm

  !-- end of declaration
  !---------------------------------------------------------------------

  nPotSets=nPotSets+1
  lengthR=length_to_char(root,'.')
  lVB=.FALSE.
  IF ( root(1:1) /= 'T' ) lVB= .TRUE. 

  !--- Copy:
  DO n_r=1,n_r_max
     DO lm =1,lm_max
        workA(lm,n_r)= b(lm,n_r)
        IF ( lVB ) workB(lm,n_r)=aj(lm,n_r)
     END DO
  END DO

  !--- Transform to Cheb-space:
  CALL costf1(workA,lm_max_real,1,lm_max_real, &
       workC,i_costf_init,d_costf_init)
  IF ( lVB ) &
       CALL costf1(workB,lm_max_real,1,lm_max_real, &
       workC,i_costf_init,d_costf_init)

  !--- Correct amplitude:
  chebNorm=DSQRT(2.D0/(n_r_max-1))
  DO n_cheb=1,n_cheb_max
     DO lm=1,lm_max
        IF ( n_cheb == 1 .OR. n_cheb == n_r_max ) THEN
           workA(lm,n_cheb)=chebNorm*0.5D0*workA(lm,n_cheb)
           IF ( lVB) workB(lm,n_cheb)=chebNorm*0.5D0*workB(lm,n_cheb)
        ELSE
           workA(lm,n_cheb)=chebNorm*workA(lm,n_cheb)
           IF ( lVB) workB(lm,n_cheb)=chebNorm*workB(lm,n_cheb)
        END IF
     END DO
  END DO

  !--- Write:
  IF ( nPotSets == 0 ) THEN ! nPotSets=-1 on call
     fileName=root(1:lengthR)//'.'//tag
  ELSE
     !------ Names including the time:
     !           IF ( l_graph_time ) THEN
     !              CALL dble2string(time,'_',6,string,length)
     !              fileName=root(1:lengthR)//'_t='//string(1:length)//'.'//tag
     !           ELSE
     !------ Numbered names:
     write(string, *) nPotSets
     fileName=root(1:lengthR)//'_'//trim(adjustl(string))//'.'//tag
     !           END IF
  END IF

  OPEN(99,FILE=fileName,FORM='UNFORMATTED',STATUS='UNKNOWN')

  WRITE(99) l_max,n_cheb_max,n_cheb_ic_max,minc,lm_max
  WRITE(99) SNGL(ra),SNGL(ek),SNGL(pr),SNGL(prmag), &
       SNGL(radratio),SNGL(sigma_ratio), &
       SNGL(omega_ma),SNGL(omega_ic)
  WRITE(99) SNGL(time), &
       ((CMPLX(REAL(workA(lm,n_cheb)),AIMAG(workA(lm,n_cheb)),KIND=KIND(0.D0)), &
       lm=1,lm_max),n_cheb=1,n_cheb_max)
  IF ( lVB ) &
       WRITE(99) SNGL(time), &
       ((CMPLX(REAL(workB(lm,n_cheb)),AIMAG(workB(lm,n_cheb)),KIND=KIND(0.D0)), &
       lm=1,lm_max),n_cheb=1,n_cheb_max)

  !-- Now inner core field
  IF ( root(1:1) == 'B' .AND. l_cond_ic ) THEN

     WRITE(*,*) 'WRITING IC DATA INTO FILE:',fileName

     DO n_r=1,n_r_ic_max
        DO lm =1,lm_max
           workA(lm,n_r)= b_ic(lm,n_r)
           workB(lm,n_r)=aj_ic(lm,n_r)
        END DO
     END DO

     CALL costf1(workA,lm_max_real,1,lm_max_real, &
          workC,i_costf1_ic_init,d_costf1_ic_init)
     CALL costf1(workB,lm_max_real,1,lm_max_real, &
          workC,i_costf1_ic_init,d_costf1_ic_init)

     chebNorm=DSQRT(2.D0/(n_r_ic_max-1))
     DO n_cheb=1,n_cheb_ic_max
        DO lm=1,lm_max
           IF ( n_cheb == 1 .OR. n_cheb == n_r_ic_max ) THEN
              workA(lm,n_cheb)=chebNorm*0.5D0*workA(lm,n_cheb)
              workB(lm,n_cheb)=chebNorm*0.5D0*workB(lm,n_cheb)
           ELSE
              workA(lm,n_cheb)=chebNorm*workA(lm,n_cheb)
              workB(lm,n_cheb)=chebNorm*workB(lm,n_cheb)
           END IF
        END DO
     END DO
     WRITE(99) SNGL(time), &
          ((CMPLX(REAL(workA(lm,n_cheb)),AIMAG(workA(lm,n_cheb)),KIND=KIND(0d0)), &
          lm=1,lm_max),n_cheb=1,n_cheb_ic_max)
     WRITE(99) SNGL(time), &
          ((CMPLX(REAL(workB(lm,n_cheb)),AIMAG(workB(lm,n_cheb)),KIND=KIND(0d0)), &
          lm=1,lm_max),n_cheb=1,n_cheb_ic_max)

  END IF

  CLOSE(99)

  RETURN
end SUBROUTINE storePot

!----------------------------------------------------------------------
