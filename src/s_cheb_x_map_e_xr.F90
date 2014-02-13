!$Id$
!*********************************************************
    subroutine cheb_x_map_e_xr(a,b,n,x,y, &
                               a1,a2,x0,lbd)
!*********************************************************

!    !------------ This is release 2 level 1  --------!
!    !------------ Created on 1/17/02  by JW. --------!


!---------------------------------------------------------
!   Given the intervall [a,b] the routine returns the
!   n+1 points that should be used to support a
!   Chebychev expansion. These are the n+1 extrema y(i) of
!   the Chebychev polynomial of degree n in the
!   intervall [-1,1].
!   The respective points mapped into the intervall of
!   queation [a,b] are the x(i).
!   NOTE: x(i) and y(i) are stored in the reversed order:
!    x(1)=b, x(n+1)=a, y(1)=1, y(n+1)=-1
!---------------------------------------------------------
     
    USE logic

    IMPLICIT NONE
     
!-- Input:
! include 'c_logic.f'

    REAL(kind=8),INTENT(IN) :: a,b  ! intervall boundaries
    INTEGER,INTENT(IN) :: n   ! degree of Cheb polynomial to
! be represented by the grid points
    REAL(kind=8),INTENT(out) :: x(*) ! grid points in intervall [a,b]
    REAL(kind=8),INTENT(OUT) :: y(*) ! grid points in intervall [-1,1]
    real(kind=8) :: a1,a2,x0,lbd
     
!-- Local variables:
    real(kind=8) :: bpa,bma,pi
    integer :: k
     
    pi=4.d0*datan(1.d0)
    bma=0.5d0*(b-a)
    bpa=0.5d0*(a+b)
     
    do k=1,n+1
        y(k)=dcos( pi*dble(k-1)/dble(n) )
        if (l_newmap) then
            x(k)=(a2+dtan(lbd*(y(k)-x0))/a1)/2 + bpa
        else
            x(k)=bma * y(k) + bpa
        end if
    end do
     
    return
    end subroutine cheb_x_map_e_xr
     
!---------------------------------------------------------
