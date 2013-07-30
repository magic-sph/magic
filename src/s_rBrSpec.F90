!$Id$
!********************************************************************
    subroutine rBrSpec(time,Pol,PolIC,fileRoot,lIC)
!********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!--------------------------------------------------------------------

!--------------------------------------------------------------------

    USE truncation
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE blocking
    USE horizontal_data
    USE logic
    USE output_data
    USE charmanip, ONLY: length_to_blank
    USE usefull, ONLY: cc2real

    IMPLICIT NONE

    REAL(kind=8) :: time

    COMPLEX(kind=8) :: Pol(lm_max,n_r_max)
    COMPLEX(kind=8) :: PolIC(lm_max,n_r_ic_max)
    CHARACTER(len=*) fileRoot
    LOGICAL :: lIC

!-- Output:
    REAL(kind=8) :: e_p_AS(l_max,n_r_tot)
    REAL(kind=8) :: e_p(l_max,n_r_tot)

!-- Local:
    CHARACTER(len=72) :: specFile
    INTEGER :: n_r,lm,l,m

    REAL(kind=8) :: fac,O_r_icb_E_2,rRatio,amp
    REAL(kind=8) :: e_p_temp

    INTEGER :: length

!-- Local function:
    LOGICAL :: lAS

!-- end of declaration
!---------------------------------------------------------------------

! Factor energy scale, 1/2 (Energy)  and 1/(4 Pi) for surface to get density
!             1/r**2 applied below
    fac=0.5D0*eScale/(16.D0*DATAN(1.D0))
    length=length_to_blank(fileRoot)

    DO n_r=1,n_r_max
        DO l=1,6
            e_p(l,n_r)=0.D0
        END DO
        DO lm=2,lm_max
            l=lm2l(lm)
            IF ( l <= 6 ) THEN
                m=lm2m(lm)
                amp=REAL(Pol(lm,n_r))
                e_p_temp=dLh(lm)*dLh(lm)*or2(n_r)*cc2real(Pol(lm,n_r),m)
                IF ( m == 0 ) THEN
                    IF ( DABS(amp)/=0.d0 ) e_p_AS(l,n_r)=fac*amp/DABS(amp)*e_p_temp
                END IF
                e_p(l,n_r)=e_p(l,n_r)+fac*e_p_temp
            END IF
        END DO    ! do loop over lms in block
    END DO    ! radial grid points

!-- Inner core:

    IF ( lIC ) then

        lAS=.true.
        IF ( fileRoot(1:length) == 'rBrAdvSpec' ) lAS= .FALSE. 

        O_r_icb_E_2=1.d0/r_icb**2

        DO n_r=2,n_r_ic_max
            rRatio=r_ic(n_r)/r_ic(1)
            DO l=1,6
                        e_p(l,n_r_max-1+n_r)=0.D0
                e_p_AS(l,n_r_max-1+n_r)=0.D0
            END DO
            DO lm=2,lm_max
                l=lm2l(lm)
                IF ( l <= 6 ) THEN
                    m=lm2m(lm)
                    IF ( m /= 0 .OR. lAS ) THEN
                        IF( l_cond_ic ) THEN
                            e_p_temp=                               &
                                            dLh(lm)*rRatio**(2*l) * &
                               dLh(lm)*O_r_icb_E_2*cc2real(PolIC(lm,n_r),m)
                            amp=REAL(PolIC(lm,n_r))
                        ELSE
                            e_p_temp=                               &
                                dLh(lm)*O_r_icb_E_2*rRatio**(2*l) * &
                                dLh(lm)*cc2real(PolIC(lm,n_r_ICB),m)
                            amp=REAL(Pol(lm,n_r_ICB))
                        END IF
                        IF ( m == 0 ) THEN
                            IF ( DABS(amp)/=0d0) e_p_AS(l,n_r_max-1+n_r)= &
                                                 fac*amp/DABS(amp)*e_p_temp
                        END IF
                        e_p(l,n_r_max-1+n_r)=e_p(l,n_r_max-1+n_r) + &
                                             fac*e_p_temp
                    END IF
                END IF
            END DO
        END DO

    ELSE
        DO n_r=2,n_r_ic_max
            DO l=1,6
                e_p_AS(l,n_r_max-1+n_r)=0.d0
                e_p(l,n_r_max-1+n_r)   =0.d0
            END DO
        END DO
    END IF

!-- Output into file:
!     writing l=0/1/2 magnetic energy
    specFile=fileRoot(1:length)//'.'//tag
    OPEN(91,FILE=specFile,FORM='UNFORMATTED',STATUS='UNKNOWN', &
         POSITION='APPEND')

    !IF ( nLines == 0 ) WRITE(91) FLOAT(n_r_tot-1),SNGL(radratio)
    !IF ( nLines == 0 ) WRITE(91) &
    !    (SNGL(r(n_r)),n_r=1,n_r_max),(SNGL(r_ic(n_r)),n_r=2,n_r_ic_max)
    WRITE(91) SNGL(time),                          &
            (SNGL(e_p(1,n_r)),n_r=1,n_r_tot-1),    &
            (SNGL(e_p(2,n_r)),n_r=1,n_r_tot-1),    &
            (SNGL(e_p(3,n_r)),n_r=1,n_r_tot-1),    &
            (SNGL(e_p(4,n_r)),n_r=1,n_r_tot-1),    &
            (SNGL(e_p(5,n_r)),n_r=1,n_r_tot-1),    &
            (SNGL(e_p(6,n_r)),n_r=1,n_r_tot-1)
    WRITE(91) REAL(time),                          &
            (SNGL(e_p_AS(1,n_r)),n_r=1,n_r_tot-1), &
            (SNGL(e_p_AS(2,n_r)),n_r=1,n_r_tot-1), &
            (SNGL(e_p_AS(3,n_r)),n_r=1,n_r_tot-1), &
            (SNGL(e_p_AS(4,n_r)),n_r=1,n_r_tot-1), &
            (SNGL(e_p_AS(5,n_r)),n_r=1,n_r_tot-1), &
            (SNGL(e_p_AS(6,n_r)),n_r=1,n_r_tot-1)

    CLOSE(91)

    RETURN
    end subroutine rBrSpec

!----------------------------------------------------------------------------
