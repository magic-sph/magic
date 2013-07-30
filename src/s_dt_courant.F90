!$Id$
!***************************************************************************
    SUBROUTINE dt_courant(dt_r,dt_h,l_new_dt,dt,dt_new, &
                                     dtMax,dtrkc,dthkc)
!***************************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!---------------------------------------------------------------------------

! *** Check if Courant criterion based on combined
! *** fluid and Alfven velocity is satisfied
! *** Returns new value of time step dtnew

!     dtr,dth: (output) radial/horizontal Courant time step
!     n_time_step: (input) time step number
!     l_new_dt: (output) flag indicating that time step is changed (=1) or not (=0)
!     dt: (input) old time step
!     dtnew: (output) new time step
!     dtMin: (input) lower limit for time step (termination if dtnew < dtMin)
!     dtMax: (input) upper limit for time step
!     dtrkc: (input) radial Courant time step as function of radial level
!     dthkc: (input) horizontal Courant time step as function of radial level

!---------------------------------------------------------------------------

    USE truncation
    USE logic
    USE output_data

    IMPLICIT NONE

!-- Input:
    REAL(kind=8) :: dt
    REAL(kind=8) :: dtMax
    REAL(kind=8) :: dtrkc(*),dthkc(*)

!-- Output:
    LOGICAL :: l_new_dt
    REAL(kind=8) :: dt_new
    REAL(kind=8) :: dt_r,dt_h

!-- Local:
    INTEGER :: n_r
    REAL(kind=8) :: dt_rh,dt_2
    REAL(kind=8) :: dt_fac


!-- end of declaration
!---------------------------------------------------------------------------


    dt_fac=2.D0
    dt_r  =1000.D0*dtMax
    dt_h  =dt_r
    DO n_r=1,n_r_max
        dt_r=min(dtrkc(n_r),dt_r)
        dt_h=min(dthkc(n_r),dt_h)
    END DO

    dt_rh=MIN(dt_r,dt_h)
    dt_2 =MIN(0.5D0*(1.D0/dt_fac+1.D0)*dt_rh,dtMax)

    IF ( dt > dtMax ) THEN

        l_new_dt=.TRUE.
        dt_new=dtMax
        WRITE(6,'(1P,/," ! COURANT: dt=dtMax =",d12.4)') dtMax
        WRITE(6,'(1P,/," ! Think about changing dtMax !")')
        IF ( l_save_out ) THEN
            OPEN(n_log_file,FILE=log_file,STATUS='UNKNOWN', &
                 POSITION='APPEND')
            WRITE(n_log_file,'(1P,/," ! COURANT: dt=dtMax=",D12.4) &
          &       ') dtMax
            WRITE(n_log_file, &
          &       '(1P,/," ! Think about changing dtMax !")')
            CLOSE(n_log_file)
        ELSE
            WRITE(n_log_file,'(1P,/," ! COURANT: dt=dtMax=",D12.4) &
          &       ') dtMax
            WRITE(n_log_file, &
          &       '(1P,/," ! Think about changing dtMax !")')
        END IF

    ELSE IF ( dt > dt_rh ) THEN

        l_new_dt=.TRUE.
        dt_new  =dt_2
        WRITE(6,'(1P,/," ! COURANT: dt=",D11.4," > dt_r=",D12.4, &
          &       " and dt_h=",D12.4)') dt,dt_r,dt_h
        IF ( l_save_out ) THEN
            OPEN(n_log_file,FILE=log_file,STATUS='UNKNOWN', &
                 POSITION='APPEND')
            WRITE(n_log_file, &
           &      '(1P,/," ! COURANT: dt=",D11.4," > dt_r=",D12.4, &
           &      " and dt_h=",D12.4)') dt,dt_r,dt_h
            CLOSE(n_log_file)
        ELSE
            WRITE(n_log_file, &
           &      '(1P,/," ! COURANT: dt=",D11.4," > dt_r=",D12.4, &
           &      " and dt_h=",D12.4)') dt,dt_r,dt_h
        END IF

    ELSE IF ( dt_fac*dt < dt_rh .AND. dt < dtMax ) THEN
         
        l_new_dt=.TRUE.
        dt_new=dt_2
        WRITE(6,'(/," ! COURANT: ",F4.1,1P,"*dt=",D11.4, &
          &     " < dt_r=",D12.4," and dt_h=",D12.4)') &
          &     dt_fac,dt_fac*dt,dt_r,dt_h
        IF ( l_save_out ) THEN
            OPEN(n_log_file,FILE=log_file,STATUS='UNKNOWN', &
                 POSITION='APPEND')
            WRITE(n_log_file, &
              &   '(/," ! COURANT: ",F4.1,1P,"*dt=",D11.4, &
              &   " < dt_r=",D12.4," and dt_h=",D12.4)') &
              &   dt_fac,dt_fac*dt,dt_r,dt_h
            CLOSE(n_log_file)
        ELSE
            WRITE(n_log_file, &
              &   '(/," ! COURANT: ",F4.1,1P,"*dt=",D11.4, &
              &   " < dt_r=",D12.4," and dt_h=",D12.4)') &
              &   dt_fac,dt_fac*dt,dt_r,dt_h
        END IF
                   
    END IF

    IF ( dt == dt_new ) l_new_dt= .FALSE. 

    RETURN
    end SUBROUTINE dt_courant

!--------------------------------------------------------------------------
