!$Id$
!***********************************************************************
    SUBROUTINE get_wpMat(dt,l,hdif,wpMat,wpPivot)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  |  Purpose of this subroutine is to contruct the time step matricies|
!  |  zmat(i,j) and wpmat  for the NS equation.                        |
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
    USE logic
    USE algebra, ONLY: sgefa

    IMPLICIT NONE

    REAL(kind=8) :: dt
    INTEGER :: l
    REAL(kind=8) :: hdif

!-- Output:
    REAL(kind=8) :: wpMat(2*n_r_max,2*n_r_max)
    INTEGER :: wpPivot(2*n_r_max)

!-- local variables:
    INTEGER :: nR,nCheb,nR_p,nCheb_p
    INTEGER :: info
    REAL(kind=8) :: O_dt,dLh

!-- end of declaration
!-----------------------------------------------------------------------

    O_dt=1.D0/dt
    dLh =DBLE(l*(l+1))

!-- Now mode l>0

!----- Boundary conditions, see above:
    DO nCheb=1,n_cheb_max
        nCheb_p=nCheb+n_r_max

        wpMat(1,nCheb)        =cheb_norm*cheb(nCheb,1)
        wpMat(1,nCheb_p)      =0.D0
        wpMat(n_r_max,nCheb)  =cheb_norm*cheb(nCheb,n_r_max)
        wpMat(n_r_max,nCheb_p)=0.D0

        IF ( ktopv == 1 ) then  ! free slip !
            wpMat(n_r_max+1,nCheb)=   cheb_norm * ( &
                                  d2cheb(nCheb,1) - &
             (2.d0*or1(1)+beta(1))*dcheb(nCheb,1) )
        ELSE                    ! no slip, note exception for l=1,m=0
            wpMat(n_r_max+1,nCheb)=cheb_norm*dcheb(nCheb,1)
        END IF
        wpMat(n_r_max+1,nCheb_p)=0.D0

        IF ( kbotv == 1 ) then  ! free slip !
            wpMat(2*n_r_max,nCheb)=        cheb_norm * ( &
                                 d2cheb(nCheb,n_r_max) - &
                (2.d0*or1(n_r_max)+beta(n_r_max))*dcheb(nCheb,n_r_max))
        ELSE                 ! no slip, note exception for l=1,m=0
            wpMat(2*n_r_max,nCheb)=cheb_norm * dcheb(nCheb,n_r_max)
        END IF
        wpMat(2*n_r_max,nCheb_p)=0.D0

    END DO   !  loop over nCheb

    IF ( n_cheb_max < n_r_max ) then ! fill with zeros !
        DO nCheb=n_cheb_max+1,n_r_max
            nCheb_p=nCheb+n_r_max
            wpMat(1,nCheb)          =0.D0
            wpMat(n_r_max,nCheb)    =0.D0
            wpMat(n_r_max+1,nCheb)  =0.D0
            wpMat(2*n_r_max,nCheb)  =0.D0
            wpMat(1,nCheb_p)        =0.D0
            wpMat(n_r_max,nCheb_p)  =0.D0
            wpMat(n_r_max+1,nCheb_p)=0.D0
            wpMat(2*n_r_max,nCheb_p)=0.D0
        END DO
    END IF

!----- Other points:
    DO nCheb=1,n_r_max
        nCheb_p=nCheb+n_r_max
        DO nR=2,n_r_max-1
            nR_p=nR+n_r_max
            wpMat(nR,nCheb)=                cheb_norm * ( &
                        O_dt*dLh*or2(nR)*cheb(nCheb,nR) - &
                      alpha*hdif*visc(nR)*dLh*or2(nR) * ( &
                                       d2cheb(nCheb,nR) + &
            (2.*dLvisc(nR)-beta(nR)/3.D0)*   dcheb(nCheb,nR) - &
               ( dLh*or2(nR)+ 4.D0/3.D0*(dLvisc(nR)*beta(nR) + &
              (3.d0*dLvisc(nR)+beta(nR))*or1(nR)+dbeta(nR)) )* &
                                             cheb(nCheb,nR) ) )

            wpMat(nR,nCheb_p)= cheb_norm * &
                    alpha*(dcheb(nCheb,nR) &
               -  beta(nR)* cheb(nCheb,nR))
            wpMat(nR_p,nCheb)=               cheb_norm * ( &
                       -O_dt*dLh*or2(nR)*dcheb(nCheb,nR) - &
                       alpha*hdif*visc(nR)*dLh*or2(nR) * ( &
                                       -d3cheb(nCheb,nR) + &
            ( beta(nR)-dLvisc(nR) )*      d2cheb(nCheb,nR) + &
                ( dLh*or2(nR)+dbeta(nR)+dLvisc(nR)*beta(nR)+ &
                        2.D0*(dLvisc(nR)+beta(nR))*or1(nR))* &
                                           dcheb(nCheb,nR) - &
               dLh*or2(nR)*( 2.d0*or1(nR)                  + &
                            dLvisc(nR)+2.d0/3.d0*beta(nR) )* &
                                           cheb(nCheb,nR) ) )
            wpMat(nR_p,nCheb_p)= - cheb_norm * &
                                alpha*dLh*or2(nR)*cheb(nCheb,nR)
        END DO
    END DO


!----- Factor for highest and lowest cheb:
    DO nR=1,n_r_max
        nR_p=nR+n_r_max
        wpMat(nR,1)          =0.5D0*wpMat(nR,1)
        wpMat(nR,n_r_max)    =0.5D0*wpMat(nR,n_r_max)
        wpMat(nR,n_r_max+1)  =0.5D0*wpMat(nR,n_r_max+1)
        wpMat(nR,2*n_r_max)  =0.5D0*wpMat(nR,2*n_r_max)
        wpMat(nR_p,1)        =0.5D0*wpMat(nR_p,1)
        wpMat(nR_p,n_r_max)  =0.5D0*wpMat(nR_p,n_r_max)
        wpMat(nR_p,n_r_max+1)=0.5D0*wpMat(nR_p,n_r_max+1)
        wpMat(nR_p,2*n_r_max)=0.5D0*wpMat(nR_p,2*n_r_max)
    END DO

    CALL sgefa(wpMat,2*n_r_max,2*n_r_max,wpPivot,info)
    IF ( info /= 0 ) THEN
        WRITE(*,*) 'Singular matrix wpmat!'
        STOP '35'
    END IF

            
    RETURN
    end SUBROUTINE get_wpMat

!-----------------------------------------------------------------------------
