!$Id$
!***********************************************************************
    subroutine s_cond(s0)
!***********************************************************************

!     !------------ This is release 2 level 1  --------------!
!     !------------ Created on 1/17/02  by JW. --------------!

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  |  Purpose of this subroutine is to solve the entropy equation      |
!  |  for an the conductive (l=0,m=0)-mode.                            |
!  |  Output is the radial dependence of the solution in s0.           |
!  |                                                                   |
!  |  Called in s_init_s.f                                             |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+

    USE truncation
    USE radial_functions
    USE physical_parameters
    USE init_fields
    USE horizontal_data
    USE matrices
    USE algebra, ONLY: cgesl,sgefa

    IMPLICIT NONE

!-- output:
    REAL(kind=8),intent(OUT) :: s0(*)

!-- local:
    integer :: n_cheb,n_r,info
    complex(kind=8) :: rhs(n_r_max)
    real(kind=8) :: work(n_r_max)

!-- end of declaration
!--------------------------------------------------------------------

     
!-- Set Matrix:
    do n_cheb=1,n_r_max
        do n_r=2,n_r_max-1
            s0Mat(n_r,n_cheb)=cheb_norm*opr*kappa(n_r)* ( &
                                     d2cheb(n_cheb,n_r) + &
           ( 2.d0*or1(n_r)+beta(n_r)+otemp1(n_r)*dtemp0(n_r)+ &
                         dLkappa(n_r) )* dcheb(n_cheb,n_r)  )
        end do
    end do
     

!-- Set boundary conditions:
    do n_cheb=1,n_cheb_max
        if ( ktops == 1 .OR. kbots /= 1 ) then
            s0Mat(1,n_cheb)=cheb_norm
        else
            s0Mat(1,n_cheb)=dcheb(n_cheb,1)*cheb_norm
        end if
        if ( kbots == 1 ) then   ! Constant entropy at ICB
            s0Mat(n_r_max,n_cheb)=cheb(n_cheb,n_r_max)*cheb_norm
        else                     ! Constant heat flux at ICB
            s0Mat(n_r_max,n_cheb)=dcheb(n_cheb,n_r_max)*cheb_norm
        end if
    end do
     
!-- Fill with zeros:
    if ( n_cheb_max < n_r_max ) then
        do n_cheb=n_cheb_max+1,n_r_max
            s0Mat(1,n_cheb)=0.d0
            s0Mat(n_r_max,n_cheb)=0.d0
        end do
    end if
     
!-- Renormalize:
    do n_r=1,n_r_max
        s0Mat(n_r,1)=0.5*s0Mat(n_r,1)
        s0Mat(n_r,n_r_max)=0.5*s0Mat(n_r,n_r_max)
    end do
     
!-- Invert matrix:
    call sgefa(s0Mat,n_r_max,n_r_max,s0Pivot,info)
    if ( info /= 0 ) then
        write(*,*) '! Singular Matrix s0Mat in init_s!'
        stop '20'
    end if
     
!-- Set source terms in RHS:
    do n_r=2,n_r_max-1
        rhs(n_r)=-epsc*epscProf(n_r)*orho1(n_r)
    end do
     
!-- Set boundary values:
    if ( ktops /= 1 .AND. kbots /= 1 ) then
        rhs(1)=0.d0
    else
        rhs(1)=REAL(tops(0,0))
    end if
    rhs(n_r_max)=REAL(bots(0,0))
     
!-- Solve for s0:
    call cgesl(s0Mat,n_r_max,n_r_max,s0Pivot,rhs)
     
!-- Copy result to s0:
    do n_r=1,n_r_max
        s0(n_r)=REAL(rhs(n_r))
    end do

!-- Set cheb-modes > n_cheb_max to zero:
    if ( n_cheb_max < n_r_max ) then
        do n_cheb=n_cheb_max+1,n_r_max
            s0(n_cheb)=0.d0
        end do
    end if
     
!-- Transform to radial space:
    call costf1(s0,1,1,1,work,i_costf_init,d_costf_init)
     

    return
    end subroutine s_cond

!--------------------------------------------------------------------------------
