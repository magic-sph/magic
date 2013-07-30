!$Id$
!*********************************************************
    subroutine cheb_x_map_e(a,b,n,x,y)
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
     
    implicit none
     
!-- Input:
    real(kind=8) :: a,b  ! intervall boundaries
    integer :: n   ! degree of Cheb polynomial to
! be represented by the grid points
    real(kind=8) :: x(*) ! grid points in intervall [a,b]
    real(kind=8) :: y(*) ! grid points in intervall [-1,1]
     
!-- Local variables:
    real(kind=8) :: bpa,bma,pi
    integer :: k
     
    pi=4.d0*datan(1.d0)
    bma=0.5d0*(b-a)
    bpa=0.5d0*(a+b)
     
    do k=1,n+1
        y(k)=dcos( pi*dble(k-1)/dble(n) )
        x(k)=bma * y(k) + bpa
    end do
     
    return
    end subroutine cheb_x_map_e
     
!---------------------------------------------------------
