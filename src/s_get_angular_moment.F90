!$Id$
!***********************************************************************
    SUBROUTINE get_angular_moment(z10,z11,omega_ic,omega_ma, &
      angular_moment_oc,angular_moment_ic,angular_moment_ma)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------------------------------------------------------------+
!  |                                                                   |
!  |    Calculates angular momentum of outer core, inner core and      |
!  |    mantle. For outer core we need z(l=1|m=0,1|r), for             |
!  |    inner core and mantle the respective rotation rates are needed.|
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+
      
    USE truncation
    USE radial_functions
    USE blocking
    USE const
    USE integration, ONLY: rInt_R

    IMPLICIT NONE

!-- Input of scalar fields:
    COMPLEX(kind=8) :: z10(n_r_max),z11(n_r_max)
    REAL(kind=8) :: omega_ic,omega_ma

!-- output:
    REAL(kind=8) :: angular_moment_oc(*)
    REAL(kind=8) :: angular_moment_ic(*)
    REAL(kind=8) :: angular_moment_ma(*)

!-- local variables:
    INTEGER :: n_r,n
    INTEGER :: l1m0,l1m1
    REAL(kind=8) :: f(n_r_max,3)
    REAL(kind=8) :: r_E_2             ! r**2
    REAL(kind=8) :: fac

!-- end of declaration
!-----------------------------------------------------------------------

     
!-- Calculating angular moments:

!-- Start with outer core:

!----- Construct radial function:
    l1m0=lm2(1,0)
    l1m1=lm2(1,1)
    do n_r=1,n_r_max
        r_E_2=r(n_r)*r(n_r)
        IF ( l1m1 > 0 ) THEN
            f(n_r,1)=r_E_2*REAL(z11(n_r))
            f(n_r,2)=r_E_2*AIMAG(z11(n_r))
        ELSE
            f(n_r,1)=0.D0
            f(n_r,2)=0.D0
        END IF
        f(n_r,3)=r_E_2*REAL(z10(n_r))
    end do

!----- Perform radial integral:
    DO n=1,3
        angular_moment_oc(n)=rInt_R(f(1,n),n_r_max,n_r_max,drx, &
                                    i_costf_init,d_costf_init)
    END DO

!----- Apply normalisation factors of chebs and other factors
!      plus the sign correction for y-component:
    fac=8.d0/3.d0*pi
    angular_moment_oc(1)= 2.d0*fac*y11_norm * &
                          angular_moment_oc(1)
    angular_moment_oc(2)=-2.d0*fac*y11_norm * &
                          angular_moment_oc(2)
    angular_moment_oc(3)=      fac*y10_norm * &
                          angular_moment_oc(3)

!----- Now inner core and mantle:
    angular_moment_ic(1)=0.d0
    angular_moment_ic(2)=0.d0
    angular_moment_ic(3)=c_moi_ic*omega_ic
    angular_moment_ma(1)=0.d0
    angular_moment_ma(2)=0.d0
    angular_moment_ma(3)=c_moi_ma*omega_ma

            
    return
    end SUBROUTINE get_angular_moment

!-----------------------------------------------------------------------
