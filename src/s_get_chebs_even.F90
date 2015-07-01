!$Id$
!******************************************************************
    subroutine get_chebs_even(n_r,a,b,y,n_r_max, &
                              cheb,dcheb,d2cheb,dim1,dim2)
!******************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!------------------------------------------------------------------

!  Construct even Chebychev polynomials and their first,
!  second and third derivative up to degree 2*(n_r/2) at
!  (n_r/2) points x in the intervall [a,b].
!  Since the Chebs are only defined in [-1,1] we have to
!  map the points x in [a,b] onto points y in the
!  intervall [-1,1]. This map is contructed by
!  the subroutine cheb_x_map_e.f which must be called
!  before entering this subroutine.
!  For even Chebs we need only half the point of the map,
!  these (n_r/2) points are in the intervall [1,0[ .
!  NOTE the reversed order in the points: y(1)=-1, y(n_r)=1.
!  y=0 is not used, which helps to avoid singularities.

!------------------------------------------------------------------
     
    implicit none
     
!-- INPUT:
    integer :: n_r ! number of grid points
! n_r grid points suffice for a cheb
! transform up to degree n_r-1
    integer :: n_r_max   ! max number of radial points, dims of y
    real(kind=8) :: a,b  ! intervall boundaries [a,b]
    real(kind=8) :: y(n_r_max) ! n_r grid points in intervall [a,b]
    integer :: dim1,dim2 ! dimensions of cheb,dcheb,......
     
!-- OUTPUT:
    real(kind=8) :: cheb(dim1,dim2)   ! cheb(i,j) is Chebychev pol.
! of degree i at grid point j
    real(kind=8) :: dcheb(dim1,dim2)  ! first derivative of cheb
    real(kind=8) :: d2cheb(dim1,dim2) ! second derivative o cheb
     
!-- INTERNAL VARIABLES:
    integer :: n,k   ! counter
    real(kind=8) :: map_fac ! maping factor to transfrom y-derivatives
! in [-1,1] to x-derivatives in [a,b]
    real(kind=8) :: last_cheb,last_dcheb,last_d2cheb

!-- End of declaration
!----------------------------------------------------------------
     
!-- d Cheb(y) / d x = d y / d x * d Cheb(y) / d y
!                   = map_fac * d Cheb(y) / d y
    map_fac=2.d0/(b-a)
     
!-- construction of chebs with recursion:
    do k=1,n_r
         
        cheb(1,k)=1.d0
        last_cheb=y(k)
        dcheb(1,k)=0.d0
        last_dcheb=map_fac
        d2cheb(1,k)=0.d0
        last_d2cheb=0.d0
         
        do n=2,n_r ! only even chebs stored !
             
        !-------- even chebs:
            cheb(n,k)=2.d0*y(k)*last_cheb-cheb(n-1,k)
            dcheb(n,k)=2.d0*map_fac*last_cheb + &
                         2.d0*y(k)*last_dcheb - &
                                   dcheb(n-1,k)
            d2cheb(n,k)=4.d0*map_fac*last_dcheb + &
                          2.d0*y(k)*last_d2cheb - &
                                    d2cheb(n-1,k)
             
        !-------- odd chebs: not stored but necessary for recursion
            last_cheb=2.d0*y(k)*cheb(n,k)-last_cheb
            last_dcheb=2.d0*map_fac*cheb(n,k) + &
                         2.d0*y(k)*dcheb(n,k) - &
                                     last_dcheb
            last_d2cheb=4.d0*map_fac*dcheb(n,k) + &
                          2.d0*y(k)*d2cheb(n,k) - &
                                    last_d2cheb

        end do
         
    end do


    return
    end subroutine get_chebs_even
! __________________________________________________________________
