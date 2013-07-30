!$Id$
module usefull
!************************************************************************
!  library with several usefull shit
!************************************************************************
  implicit none

  contains

!************************************************************************
    subroutine check_dim(dim,max_dim,dim_name,prog_name,stop_signal)
!************************************************************************
!  checks whether dim is not larger than max_dim !
!-----------------------------------------------------------------------

    implicit none

    character(len=*) :: prog_name,dim_name
    integer :: max_dim,dim
    integer :: stop_signal

    stop_signal=0  !? does this make sence ?
    if( dim > max_dim ) then
        write(*,'(2x,''MESSAGE FROM PROGRAM '',a)') prog_name
        write(*,'(2x,'' TOO SMALL DIMENSION '',a)') dim_name
        write(*,*) '   SHOULD BE AT LEAST',dim
        write(*,*) '   BUT IS ONLY       ',max_dim
        stop_signal=1
    endif

    return
    end subroutine check_dim
!-------------------------------------------------------------------------
    logical function l_correct_step(n,t,t_last,                &
                                    n_max,n_step,n_intervalls, &
                                    n_ts,times,n_eo)
!************************************************************************
! Suppose we have a (loop) maximum of n_max steps!
! If n_intervalls times in these steps a certain action should be carried out
! this can be invoked by l_correct_step=true if on input n_intervalls>0
! and n_step=0.
! Alternatively the action can be invoked every n_step steps if
! on input n_intervalls=0 and n_step>0.
! In both cases l_correct_step=true for n=n_max.

! The argument controlls whether in addition n should be
!       even: n_eo=2
!    or  odd: n_eo=1
!------------------------------------------------------------------------

    implicit none

!-- Input:
    integer      :: n            ! current step
    real(kind=8) ::  t            ! time at current step
    real(kind=8) ::  t_last       ! last time at current step
    integer      :: n_max        ! max number of steps
    integer      :: n_step       ! action intervall
    integer      :: n_intervalls ! number of actions
    integer      :: n_ts         ! number of times t
    real(kind=8) ::  times(*)     ! times where l_correct_step.eq.true
    integer      :: n_eo         ! even/odd controller

!-- Output: l_correct_step   !

!-- Local:
    integer :: n_delta      ! corrector for even/odd n
    integer :: n_offset     ! offset with no action
    integer :: n_t          ! counter for times
    integer :: n_steps      ! local step width


!-- End of declatation
!--------------------------------------------------------------------------


    if ( n_step /= 0 .AND. n_intervalls /= 0 ) then
        write(*,*) &
        '! ERROR MESSAGE FROM FUNCTION L_CORRECT_STEP:'
        write(*,*) &
        '! EITHER N_STEP OR N_INTERVALL HAVE TO BE ZERO!'
        stop
    end if

    l_correct_step=.false.

    if ( n_intervalls /= 0 ) then

        n_steps=n_max/n_intervalls
        if ( ( n_eo == 2 .AND. mod(n_step,2) /= 0 ) .OR. &
        ( n_eo == 1 .AND. mod(n_step,2) /= 1 ) ) then
            n_steps=n_steps+1
        end if

        n_offset=n_max-n_steps*n_intervalls

        if ( n > n_offset .AND. &
        mod(n-n_offset,n_steps) == 0 ) &
        l_correct_step= .TRUE. 

    else if ( n_step /= 0 ) then
                
        n_delta=0
        if ( ( n_eo == 1 .AND. mod(n,2) == 0 ) .OR. &
        ( n_eo == 2 .AND. mod(n,2) == 1 ) ) &
        n_delta=1

        if ( n == n_max .OR. &
        mod(n-n_delta,n_step) == 0 ) l_correct_step= .TRUE. 

    end if

    IF ( n_ts >= 1 ) THEN
        DO n_t=1,n_ts
            IF ( times(n_t) < t .AND. &
            times(n_t) >= t_last ) THEN
                l_correct_step=.TRUE.
                GOTO 100
            END IF
        END DO
    END IF
               
    100 continue

    return
    end function l_correct_step
!------------------------------------------------------------------------
    real(kind=8) function random(r)

!     random number generator

!     if ( r.eq.0 ) then
!        random(r) = next random number (between 0. and 1.)
!     if ( r.lt.0 ) then
!        random(r) = previous random number
!     if ( r.gt.0 ) then
!        random(r) = a new sequence of random numbers is started
!                  with seed r mod 1
!                  note: r must be a non-integer to get a different seq

!     called in sinit

!---------------------------------------------------------------------------

    implicit none

!-- input:
    real(kind=8) :: r

!-- local:
    integer :: ia1,ia0,ia1ma0,ic,ix1,ix0
    integer :: iy0,iy1

    save ix1, ix0


!-- end of declaration
!-------------------------------------------------------------------------
     
    ia1   =1536
    ia0   =1029
    ia1ma0=507
    ic    =1731

    if ( r == 0.d0 ) then

        iy0 = ia0*ix0
        iy1 = ia1*ix1 + ia1ma0*(ix0-ix1) + iy0
        iy0 = iy0 + ic

        ix0 = MOD(iy0,2048)
        iy1 = iy1 + (iy0-ix0)/2048
        ix1 = MOD(iy1,2048)

    else if ( r > 0.d0 ) then

        ix1 = IDINT(DMOD(r,1.D0)*4194304.D0 + 0.5D0)
        ix0 = MOD(ix1,2048)
        ix1 = (ix1-ix0)/2048

    end if
     
    random = ix1*2048 + ix0
    random = random / 4194304.d0

    return
    end function random
!-------------------------------------------------------------------------
    subroutine factorise(n,n_facs,fac,n_factors,factor)
!  +-------------+----------------+------------------------------------+
!  |  Purpose of this subroutine is factorize n into a number of       |
!  |  given factors fac(i).                                            |
!  +-------------------------------------------------------------------+
      
    implicit none

!-- Input:
    integer :: n         ! number to be factorised    !
    integer :: fac(*)    ! list of fators to be tried !
    integer :: n_facs    ! number of facs to be tried !

!-- Output:
    integer :: n_factors ! number of factors used
    integer :: factor(*) ! list of factors used

!-- Local variables:
    integer :: n_rest,n_fac
    integer :: factor_tot,factor_test


!-- end of declaration
!-----------------------------------------------------------------------


    if ( n < 1 )  then
        write(*,*) '! Error message from factorise:'
        write(*,*) '! n should be larger than 0.'
        stop
    else if ( n == 1 ) then
        n_factors=0
        factor(1)=0
        return
    end if

    n_factors=0    ! number of factors
    factor_tot=1   ! total factor
    n_rest=n       ! remaining

    do n_fac=1,n_facs  ! use new factor !

        factor_test=fac(n_fac)

        10 if ( mod(n_rest,factor_test) == 0 ) then

            n_rest=n_rest/factor_test
            n_factors=n_factors+1
            factor(n_factors)=factor_test
            factor_tot=factor_tot*factor_test

            if ( n_rest == 1 ) goto 20 ! finished !

            goto 10 ! try same factor again

        end if

    end do

    write(*,*) 'Sorry, no factorisation possible of:',n
    stop

    20 continue   ! finished, n_rest=1 !


    return
    end subroutine factorise
!-----------------------------------------------------------------------
    REAL(kind=8)  FUNCTION cc2real(c,m)
!---------------------------------------------------------------------

    IMPLICIT NONE

!-- Input:
    COMPLEX(kind=8) :: c
    INTEGER :: m

!---------------------------------------------------------------------

    IF ( m == 0 ) THEN
        cc2real=REAL(c)*REAL(c)
    ELSE
        cc2real=2.D0*(REAL(c)*REAL(c)+AIMAG(c)*AIMAG(c))
    END IF

    RETURN
    END FUNCTION cc2real
!---------------------------------------------------------------------
    REAL(kind=8) FUNCTION cc22real(c1,c2,m)
!---------------------------------------------------------------------

    IMPLICIT NONE

!-- Input:
    COMPLEX(kind=8) :: c1,c2
    INTEGER :: m

!---------------------------------------------------------------------

    IF ( m == 0 ) THEN
        cc22real=REAL(c1)*REAL(c2)
    ELSE
        cc22real=2.D0*(REAL(c1)*REAL(c2)+AIMAG(c1)*AIMAG(c2))
    END IF

    RETURN
    END FUNCTION cc22real
!---------------------------------------------------------------------
end module usefull
