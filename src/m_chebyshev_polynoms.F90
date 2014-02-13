!$Id$
MODULE chebyshev_polynoms_mod
  USE usefull, ONLY: check_dim
  USE const,only: pi
  IMPLICIT NONE

  PRIVATE

  INTERFACE get_chebs
     module procedure get_chebs_recurr
  END INTERFACE get_chebs
  
  PUBLIC :: get_chebs
CONTAINS
!**********************************************************
  SUBROUTINE get_chebs_recurr(n_r,a,b,y,n_r_max,                  &
       &                      cheb,dcheb,d2cheb,d3cheb,dim1,dim2, &
       &                      map_fac1,map_fac2,map_fac3)
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

       !----- now construct the rest with a recursion:
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
  END SUBROUTINE get_chebs_recurr

  SUBROUTINE get_chebs_direct(n_r,a,b,y,n_r_max,  &
       &                      cheb,dcheb,d2cheb,d3cheb,dim1,dim2, &
       &                      map_fac1,map_fac2,map_fac3)
    !----------------------------------------------------------
    !  Construct Chebychev polynomials and their first, second,
    !  and third derivative up to degree n_r at n_r points x
    !  in the intervall [a,b]. Since the Chebs are only defined
    !  in [-1,1] we have to use a map, mapping the points x
    !  points y in the intervall [-1,1]. This map is executed
    !  by the subroutine cheb_x_map_e.f and has to be done
    !  before calling this program.
    !----------------------------------------------------------

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
    REAL(kind=8) :: local_cheb,local_dcheb,local_d2cheb,pos,spos,local_d3cheb
    real(kind=8) :: yk
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
    DO k=1,n_r  ! do loop over the n_r grid points !
    !   DO n=1,n_r
    !      cheb(n,k)=COS(pi*n*(k-1)/(n_r-1))
    !   END DO
    !END DO

       !----- set first two chebs:
       cheb(1,k)=1.d0
       ! cheb(2,k)=COS(pi*(k-1)/(n_r-1)) !y(k)
       cheb(2,k)=y(k)
       dcheb(1,k)=0.d0
       dcheb(2,k)=map_fac1(k)
       d2cheb(1,k)=0.d0
       d2cheb(2,k)=map_fac2(k)
       d3cheb(1,k)=0.d0
       d3cheb(2,k)=map_fac3(k)

       !----- now construct the rest with an recursion:
       do n=3,n_r ! do loop over the (n-1) order of the chebs

          local_cheb = COS(pi*(n-1)*(k-1)/(n_r-1))
          cheb(n,k)=local_cheb
          !cheb(n,k)=    2.d0*y(k)*cheb(n-1,k)-cheb(n-2,k)
          !IF (ABS(local_cheb).GT.0.0D0) THEN
          !   WRITE(*,"(A,2I3,3ES20.12,ES11.3)") "Error in cheb calculation: ",n,k,&
          !        & cheb(n,k),local_cheb,cheb(n,k)-local_cheb,(cheb(n,k)-local_cheb)/local_cheb
          !END IF
          IF ((k.gt.1) .and. MODULO((n-1)*(k-1),(n_r-1)).EQ.0) THEN
             local_dcheb=0.0D0
             !dcheb(n,k) = 0.0D0
          ELSEIF (k.EQ.1) THEN
             local_dcheb = map_fac1(k)*(n-1)**2
          ELSE
             local_dcheb = map_fac1(k)*(n-1)*SIN((n-1)*pi*(k-1)/(n_r-1))/SIN(pi*(k-1)/(n_r-1))
             !dcheb(n,k)=        2.d0*map_fac1(k)*cheb(n-1,k) + &
             !     2.d0*y(k)*dcheb(n-1,k) - &
             !     dcheb(n-2,k)
          END IF
          dcheb(n,k)=        2.d0*map_fac1(k)*cheb(n-1,k) + &
               2.d0*y(k)*dcheb(n-1,k) - &
               dcheb(n-2,k)

          !IF (ABS(local_dcheb).GT.0.0D0) THEN
          !WRITE(*,"(A,2I3,3ES20.12,ES11.3)") "Error in dcheb calculation: ",n,k,&
          !        & dcheb(n,k),local_dcheb,dcheb(n,k)-local_dcheb,&
          !        &(dcheb(n,k)-local_dcheb)/local_dcheb
          !END IF
          dcheb(n,k)=local_dcheb
          
          IF (2*(k-1).EQ.n_r-1) THEN
             IF (MODULO((n-1),4).EQ.0) THEN
                local_d2cheb = -(n-1)**2*map_fac1(k)**2
             ELSEIF (MODULO((n-1),2).eq.0) then
                local_d2cheb = (n-1)**2*map_fac1(k)**2
             ELSEIF (MODULO((n-1)+3,4).eq.0) then
                ! odd chebs and (n-1)=4r-3
                local_d2cheb = (n-1)*map_fac2(k)
             ELSE
                ! odd chebs and (n-1)=4r-1
                local_d2cheb = -(n-1)*map_fac2(k)
             END IF
          ELSEIF (k.EQ.n_r) THEN
             local_d2cheb=0.0D0
             d2cheb(n,k) = 0.0D0
          ELSE
             pos=pi*REAL(k-1,kind=8)/(n_r-1)
             local_d2cheb = map_fac1(k)**2*(n-1)*( COS(pos)*SIN((n-1)*pos) &
                  &                             -(n-1)*COS((n-1)*pos)*SIN(pos) &
                  &                           )/SIN(pos)**3&
                  &         +map_fac2(k)*(n-1)*SIN((n-1)*pos)/SIN(pos)
          END IF

          d2cheb(n,k)=       2.d0*map_fac2(k)*cheb(n-1,k) + &
               4.d0*map_fac1(k)*dcheb(n-1,k) + &
               2.d0*y(k)*d2cheb(n-1,k) - &
               d2cheb(n-2,k)
          !IF (ABS(local_d2cheb).GT.0.0D0) THEN
             WRITE(*,"(A,2I3,3ES20.12,ES11.3)") "Error in d2cheb calculation: ",n,k,&
                  & d2cheb(n,k),local_d2cheb,d2cheb(n,k)-local_d2cheb,&
                  &(d2cheb(n,k)-local_d2cheb)/local_d2cheb
          !END IF
          d2cheb(n,k) = local_d2cheb

          !pos = pi*REAL(k-1,kind=8)/(n_r-1)
          !spos = SIN((n-1)*pos)
          !yk= cos(pos)
          !WRITE(*,"(2I3,4ES11.3)") n,k,pos,spos,yk,SIN(pos)
          !IF (k.EQ.1) THEN
          !   local_d3cheb=0.0D0
          !ELSE
          !   IF (n.EQ.3) THEN
          !      local_d3cheb=0.0D0
          !   ELSE
          !      local_d3cheb =( ( 3*(n-1)*spos*yk**2* map_fac1(k)**3) &
          !           &          -( 3*(n-1)**2*COS((n-1)*pos)*SIN(pos)*yk* map_fac1(k)**3) &
          !           &          +(n-1)*spos*map_fac1(k)*(1.0-yk**2)&
          !           &            *( map_fac1(k)**2*(1-(n-1)**2) &
          !           &               +3*yk*map_fac2(k)&
          !           &             )&
          !           &          -(3*(n-1)**2*COS((n-1)*pos)*map_fac1(k)*map_fac2(k)*SIN(pos)**3) &
          !           &          +((n-1)*spos*map_fac3(k)*SIN(pos)**4)&
          !           &        )/SIN(pos)**5
          !   END IF
          !END IF

          d3cheb(n,k)=       2.d0*map_fac3(k)*cheb(n-1,k) + &
               6.d0*map_fac2(k)*dcheb(n-1,k) + &
               6.d0*map_fac1(k)*d2cheb(n-1,k) + &
               2.d0*y(k)*d3cheb(n-1,k) - &
               d3cheb(n-2,k)
          !IF (ABS(local_d3cheb).GT.0.0D0) THEN
          !WRITE(*,"(A,2I3,3ES20.12,ES11.3)") "Error in d3cheb calculation: ",n,k,&
          !        & d3cheb(n,k),local_d3cheb,d3cheb(n,k)-local_d3cheb,&
          !        &(d3cheb(n,k)-local_d3cheb)/local_d3cheb
          !END IF

       end do

    end do

    return
  END SUBROUTINE get_chebs_direct

END MODULE chebyshev_polynoms_mod
! __________________________________________________________________
