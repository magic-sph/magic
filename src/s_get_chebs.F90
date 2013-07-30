!$Id$
!**********************************************************
    subroutine get_chebs(n_r,a,b,y,n_r_max,                  &
                         cheb,dcheb,d2cheb,d3cheb,dim1,dim2, &
                         map_fac1,map_fac2,map_fac3)
!**********************************************************

!    !------------ This is release 2 level 1  ----!
!    !------------ Created on 1/17/02  by JW. ----!

!----------------------------------------------------------
!  Construct Chebychev polynomials and their first, second,
!  and third derivative up to degree n_r at n_r points x
!  in the intervall [a,b]. Since the Chebs are only defined
!  in [-1,1] we have to use a map, mapping the points x
!  points y in the intervall [-1,1]. This map is executed
!  by the subroutine cheb_x_map_e.f and has to be done
!  before calling this program.
!----------------------------------------------------------
     
    use usefull, only: check_dim

    implicit none
     
!-- INPUT:
    integer :: n_r ! number of grid points
! n_r grid points suffice for a cheb
! transform up to degree n_r-1
    real(kind=8) :: a,b  ! intervall boundaries [a,b]
    integer :: n_r_max   ! leading dimension of
! cheb(i,j) and der. in calling routine
    real(kind=8) :: y(n_r_max) ! n_r grid points in intervall [a,b]
    integer :: dim1,dim2 ! dimensions of cheb,dcheb,....
    real(kind=8) :: map_fac1(n_r_max)
    real(kind=8) :: map_fac2(n_r_max)
    real(kind=8) :: map_fac3(n_r_max)
     
!-- OUTPUT:
    real(kind=8) :: cheb(dim1,dim2)   ! cheb(i,j) is Chebychev pol.
! of degree i at grid point j
    real(kind=8) :: dcheb(dim1,dim2)  ! first derivative of cheb
    real(kind=8) :: d2cheb(dim1,dim2) ! second derivative o cheb
    real(kind=8) :: d3cheb(dim1,dim2) ! third derivative of cheb
     
!-- LOCAL VARIABLES:
    integer :: n,k   ! counter
    integer :: stop_signal
    real(kind=8) :: map_fac ! maping factor to transfrom y-derivatives
! in [-1,1] to x-derivatives in [a,b]

!-- End of declaration
!--------------------------------------------------------------------


    call check_dim(n_r,n_r_max, &
                   'n_r_max','get_even_chebs',stop_signal)
    call check_dim(n_r,dim2,    &
                   'dim2','get_even_chebs',stop_signal)
    if( stop_signal == 1 ) stop
     
!-- definition of map_fac:
!   d Cheb(y) / d x = d y / d x * d Cheb(y) / d y
!                   = map_fac * d Cheb(y) / d y
    map_fac=2.d0/(b-a)
     
!-- construction of chebs and derivatives with recursion:
    do k=1,n_r  ! do loop over the n_r grid points !
         
    !----- set first two chebs:
        cheb(1,k)=1.d0
        cheb(2,k)=y(k)
        dcheb(1,k)=0.d0
        dcheb(2,k)=map_fac1(k)
        d2cheb(1,k)=0.d0
        d2cheb(2,k)=map_fac2(k)
        d3cheb(1,k)=0.d0
        d3cheb(2,k)=map_fac3(k)
         
    !----- now construct the rest with an recursion:
        do n=3,n_r ! do loop over the (n-1) order of the chebs

            cheb(n,k)=    2.d0*y(k)*cheb(n-1,k)-cheb(n-2,k)
            dcheb(n,k)=        2.d0*map_fac1(k)*cheb(n-1,k) + &
                                     2.d0*y(k)*dcheb(n-1,k) - &
                                               dcheb(n-2,k)
            d2cheb(n,k)=       2.d0*map_fac2(k)*cheb(n-1,k) + &
                              4.d0*map_fac1(k)*dcheb(n-1,k) + &
                                    2.d0*y(k)*d2cheb(n-1,k) - &
                                              d2cheb(n-2,k)
            d3cheb(n,k)=       2.d0*map_fac3(k)*cheb(n-1,k) + &
                              6.d0*map_fac2(k)*dcheb(n-1,k) + &
                             6.d0*map_fac1(k)*d2cheb(n-1,k) + &
                                    2.d0*y(k)*d3cheb(n-1,k) - &
                                              d3cheb(n-2,k)
             
        end do

    end do
     
    return
    end subroutine get_chebs
! __________________________________________________________________
