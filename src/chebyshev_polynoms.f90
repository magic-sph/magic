module chebyshev_polynoms_mod

   use precision_mod
   use logic, only: l_newmap
   use const, only: pi, half, one, two, four
 
   implicit none
 
   private
 
   interface get_chebs
      module procedure get_chebs_recurr
   end interface get_chebs
   
   public :: get_chebs, cheb_grid, get_chebs_even

contains

   subroutine get_chebs_recurr(n_r,a,b,y,n_r_max,                  &
        &                      cheb,dcheb,d2cheb,d3cheb,dim1,dim2, &
        &                      map_fac1,map_fac2,map_fac3)
      !
      !  Construct Chebychev polynomials and their first, second,
      !  and third derivative up to degree n_r at n_r points x
      !  in the interval [a,b]. Since the Chebs are only defined
      !  in [-1,1] we have to use a map, mapping the points x
      !  points y in the interval [-1,1]. This map is executed
      !  by the subroutine cheb_grid and has to be done
      !  before calling this program.
      !

      !-- Input variables:
      integer,  intent(in) :: n_r       ! number of grid points
      ! n_r grid points suffice for a cheb
      ! transform up to degree n_r-1
      real(cp), intent(in) :: a,b       ! interval boundaries [a,b]
      integer,  intent(in) :: n_r_max   ! leading dimension of
      ! cheb(i,j) and der. in calling routine
      real(cp), intent(in) :: y(n_r_max)! n_r grid points in interval [a,b]
      integer,  intent(in) :: dim1,dim2 ! dimensions of cheb,dcheb,....
      real(cp), intent(in) :: map_fac1(n_r_max)
      real(cp), intent(in) :: map_fac2(n_r_max)
      real(cp), intent(in) :: map_fac3(n_r_max)

      !-- Output variables:
      real(cp), intent(out) :: cheb(dim1,dim2)   ! cheb(i,j) is Chebychev pol.
      ! of degree i at grid point j
      real(cp), intent(out) :: dcheb(dim1,dim2)  ! first derivative of cheb
      real(cp), intent(out) :: d2cheb(dim1,dim2) ! second derivative o cheb
      real(cp), intent(out) :: d3cheb(dim1,dim2) ! third derivative of cheb

      !-- Local variables:
      integer :: n,k   ! counter
      real(cp) :: map_fac ! maping factor to transfrom y-derivatives
      ! in [-1,1] to x-derivatives in [a,b]

      !-- definition of map_fac:
      !   d Cheb(y) / d x = d y / d x * d Cheb(y) / d y
      !                   = map_fac * d Cheb(y) / d y
      map_fac=two/(b-a)

      !-- construction of chebs and derivatives with recursion:
      do k=1,n_r  ! do loop over the n_r grid points !

         !----- set first two chebs:
         cheb(1,k)=one
         cheb(2,k)=y(k)
         dcheb(1,k)=0.0_cp
         dcheb(2,k)=map_fac1(k)
         d2cheb(1,k)=0.0_cp
         d2cheb(2,k)=map_fac2(k)
         d3cheb(1,k)=0.0_cp
         d3cheb(2,k)=map_fac3(k)

         !----- now construct the rest with a recursion:
         do n=3,n_r ! do loop over the (n-1) order of the chebs

            cheb(n,k)=    two*y(k)*cheb(n-1,k)-cheb(n-2,k)
            dcheb(n,k)=        two*map_fac1(k)*cheb(n-1,k) + &
                                     two*y(k)*dcheb(n-1,k) - &
                                               dcheb(n-2,k)
            d2cheb(n,k)=       two*map_fac2(k)*cheb(n-1,k) + &
                             four*map_fac1(k)*dcheb(n-1,k) + &
                                    two*y(k)*d2cheb(n-1,k) - &
                                              d2cheb(n-2,k)
            d3cheb(n,k)=       two*map_fac3(k)*cheb(n-1,k) + &
                           6.0_cp*map_fac2(k)*dcheb(n-1,k) + &
                          6.0_cp*map_fac1(k)*d2cheb(n-1,k) + &
                                    two*y(k)*d3cheb(n-1,k) - &
                                              d3cheb(n-2,k)

         end do

      end do

   end subroutine get_chebs_recurr
!------------------------------------------------------------------------------
   subroutine get_chebs_even(n_r,a,b,y,n_r_max, &
                             cheb,dcheb,d2cheb,dim1,dim2)
      !
      !  Construct even Chebychev polynomials and their first,
      !  second and third derivative up to degree 2*(n_r/2) at
      !  (n_r/2) points x in the interval [a,b].
      !  Since the Chebs are only defined in [-1,1] we have to
      !  map the points x in [a,b] onto points y in the
      !  interval [-1,1]. This map is contructed by
      !  the subroutine cheb_x_map_e.f which must be called
      !  before entering this subroutine.
      !  For even Chebs we need only half the point of the map,
      !  these (n_r/2) points are in the interval [1,0[ .
      !  NOTE the reversed order in the points: y(1)=-1, y(n_r)=1.
      !  y=0 is not used, which helps to avoid singularities.
      !
       
      !-- Input variables:
      integer,  intent(in) :: n_r ! number of grid points
                                      ! n_r grid points suffice for a cheb
                                      ! transform up to degree n_r-1
      integer,  intent(in):: n_r_max     ! max number of radial points, dims of y
      real(cp), intent(in) :: a,b        ! interval boundaries [a,b]
      real(cp), intent(in) :: y(n_r_max) ! n_r grid points in interval [a,b]
      integer,  intent(in) :: dim1,dim2  ! dimensions of cheb,dcheb,......
       
      !-- Output variables:
      real(cp), intent(out) :: cheb(dim1,dim2)   ! cheb(i,j) is Chebychev pol.
                                                     ! of degree i at grid point j
      real(cp), intent(out) :: dcheb(dim1,dim2)  ! first derivative of cheb
      real(cp), intent(out) :: d2cheb(dim1,dim2) ! second derivative o cheb
       
      !-- Internal variables:
      integer :: n,k   ! counter
      real(cp) :: map_fac ! maping factor to transfrom y-derivatives
                              ! in [-1,1] to x-derivatives in [a,b]
      real(cp) :: last_cheb,last_dcheb,last_d2cheb

      !-- d Cheb(y) / d x = d y / d x * d Cheb(y) / d y
      !                   = map_fac * d Cheb(y) / d y
      map_fac=two/(b-a)
       
      !-- construction of chebs with recursion:
      do k=1,n_r
         cheb(1,k)=one
         last_cheb=y(k)
         dcheb(1,k)=0.0_cp
         last_dcheb=map_fac
         d2cheb(1,k)=0.0_cp
         last_d2cheb=0.0_cp
         do n=2,n_r ! only even chebs stored !
            !-- even chebs:
            cheb(n,k)=two*y(k)*last_cheb-cheb(n-1,k)
            dcheb(n,k)=two*map_fac*last_cheb + &
                         two*y(k)*last_dcheb - &
                                   dcheb(n-1,k)
            d2cheb(n,k)=four*map_fac*last_dcheb + &
                           two*y(k)*last_d2cheb - &
                                      d2cheb(n-1,k)
             
            !-- odd chebs: not stored but necessary for recursion
            last_cheb=two*y(k)*cheb(n,k)-last_cheb
            last_dcheb=two*map_fac*cheb(n,k) + &
                         two*y(k)*dcheb(n,k) - &
                                     last_dcheb
            last_d2cheb=four*map_fac*dcheb(n,k) + &
                           two*y(k)*d2cheb(n,k) - &
                                     last_d2cheb
         end do
      end do

   end subroutine get_chebs_even
!------------------------------------------------------------------------------
   subroutine get_chebs_direct(n_r,a,b,y,n_r_max,                 &
       &                      cheb,dcheb,d2cheb,d3cheb,dim1,dim2, &
       &                      map_fac1,map_fac2,map_fac3)
      !
      !  Construct Chebychev polynomials and their first, second,
      !  and third derivative up to degree n_r at n_r points x
      !  in the interval [a,b]. Since the Chebs are only defined
      !  in [-1,1] we have to use a map, mapping the points x
      !  points y in the interval [-1,1]. This map is executed
      !  by the subroutine cheb_grid and has to be done
      !  before calling this program.
      !

      !-- Input variables:
      integer,  intent(in) :: n_r       ! number of grid points
      ! n_r grid points suffice for a cheb
      ! transform up to degree n_r-1
      real(cp), intent(in) :: a,b       ! interval boundaries [a,b]
      integer,  intent(in) :: n_r_max   ! leading dimension of
      ! cheb(i,j) and der. in calling routine
      real(cp), intent(in) :: y(n_r_max)! n_r grid points in interval [a,b]
      integer,  intent(in):: dim1,dim2  ! dimensions of cheb,dcheb,....
      real(cp), intent(in) :: map_fac1(n_r_max)
      real(cp), intent(in) :: map_fac2(n_r_max)
      real(cp), intent(in) :: map_fac3(n_r_max)

      !-- Output variables:
      real(cp), intent(out) :: cheb(dim1,dim2)   ! cheb(i,j) is Chebychev pol.
      ! of degree i at grid point j
      real(cp), intent(out) :: dcheb(dim1,dim2)  ! first derivative of cheb
      real(cp), intent(out) :: d2cheb(dim1,dim2) ! second derivative o cheb
      real(cp), intent(out) :: d3cheb(dim1,dim2) ! third derivative of cheb

      !-- Local variables:
      integer :: n,k   ! counter
      real(cp) :: map_fac ! maping factor to transfrom y-derivatives
      real(cp) :: local_cheb,local_dcheb,local_d2cheb,pos
      !real(cp) :: spos,local_d3cheb,yk
      ! in [-1,1] to x-derivatives in [a,b]

      !-- definition of map_fac:
      !   d Cheb(y) / d x = d y / d x * d Cheb(y) / d y
      !                   = map_fac * d Cheb(y) / d y
      map_fac=two/(b-a)

      !-- construction of chebs and derivatives with recursion:
      do k=1,n_r  ! do loop over the n_r grid points !
      !   do n=1,n_r
      !      cheb(n,k)=cos(pi*n*(k-1)/(n_r-1))
      !   end do
      !end do

         !----- set first two chebs:
         cheb(1,k)=one
         ! cheb(2,k)=cos(pi*(k-1)/(n_r-1)) !y(k)
         cheb(2,k)=y(k)
         dcheb(1,k)=0.0_cp
         dcheb(2,k)=map_fac1(k)
         d2cheb(1,k)=0.0_cp
         d2cheb(2,k)=map_fac2(k)
         d3cheb(1,k)=0.0_cp
         d3cheb(2,k)=map_fac3(k)

         !----- now construct the rest with an recursion:
         do n=3,n_r ! do loop over the (n-1) order of the chebs

            local_cheb = cos(pi*(n-1)*(k-1)/(n_r-1))
            cheb(n,k)=local_cheb
            !cheb(n,k)=    two*y(k)*cheb(n-1,k)-cheb(n-2,k)
            !if (ABS(local_cheb) > 0.0_cp) then
            !   write(*,"(A,2I3,3ES20.12,ES11.3)") "Error in cheb calculation: ",n,k,&
            !        & cheb(n,k),local_cheb,cheb(n,k)-local_cheb,                    & 
            !        & (cheb(n,k)-local_cheb)/local_cheb
            !end if
            if ((k > 1) .and. modulo((n-1)*(k-1),(n_r-1)) == 0) then
               local_dcheb=0.0_cp
               !dcheb(n,k) = 0.0_cp
            elseif (k == 1) then
               local_dcheb = map_fac1(k)*(n-1)**2
            else
               local_dcheb = map_fac1(k)*(n-1)*sin((n-1)*pi*(k-1)/(n_r-1))/ &
                                               sin(pi*(k-1)/(n_r-1))
               !dcheb(n,k)=        two*map_fac1(k)*cheb(n-1,k) + &
               !     two*y(k)*dcheb(n-1,k) - &
               !     dcheb(n-2,k)
            end if
            dcheb(n,k)=  two*map_fac1(k)*cheb(n-1,k) + &
                               two*y(k)*dcheb(n-1,k) - &
                                         dcheb(n-2,k)

            !if (ABS(local_dcheb) > 0.0_cp) then
            !write(*,"(A,2I3,3ES20.12,ES11.3)") "Error in dcheb calculation: ",n,k,&
            !        & dcheb(n,k),local_dcheb,dcheb(n,k)-local_dcheb,&
            !        &(dcheb(n,k)-local_dcheb)/local_dcheb
            !end if
            dcheb(n,k)=local_dcheb
            
            if (2*(k-1) == n_r-1) then
               if (modulo((n-1),4) == 0) then
                  local_d2cheb = -(n-1)**2*map_fac1(k)**2
               elseif (modulo((n-1),2) == 0) then
                  local_d2cheb = (n-1)**2*map_fac1(k)**2
               elseif (modulo((n-1)+3,4) == 0) then
                  ! odd chebs and (n-1)=4r-3
                  local_d2cheb = (n-1)*map_fac2(k)
               else
                  ! odd chebs and (n-1)=4r-1
                  local_d2cheb = -(n-1)*map_fac2(k)
               end if
            elseif (k == n_r) then
               local_d2cheb=0.0_cp
               d2cheb(n,k) = 0.0_cp
            else
               pos=pi*real(k-1,cp)/(n_r-1)
               local_d2cheb = map_fac1(k)**2*(n-1)*( cos(pos)*sin((n-1)*pos)     &
                    &                             -(n-1)*cos((n-1)*pos)*sin(pos) &
                    &                           )/sin(pos)**3                    &
                    &         +map_fac2(k)*(n-1)*sin((n-1)*pos)/sin(pos)
            end if

            d2cheb(n,k)=       two*map_fac2(k)*cheb(n-1,k) + &
                             four*map_fac1(k)*dcheb(n-1,k) + &
                                    two*y(k)*d2cheb(n-1,k) - &
                                              d2cheb(n-2,k)
            !if (ABS(local_d2cheb) > 0.0_cp) then
               write(*,"(A,2I3,3ES20.12,ES11.3)") "Error in d2cheb calculation: ",n,k,&
                    & d2cheb(n,k),local_d2cheb,d2cheb(n,k)-local_d2cheb,&
                    &(d2cheb(n,k)-local_d2cheb)/local_d2cheb
            !end if
            d2cheb(n,k) = local_d2cheb

            !pos = pi*real(k-1,cp)/(n_r-1)
            !spos = sin((n-1)*pos)
            !yk= cos(pos)
            !write(*,"(2I3,4ES11.3)") n,k,pos,spos,yk,sin(pos)
            !if (k == 1) then
            !   local_d3cheb=0.0_cp
            !else
            !   if (n == 3) then
            !      local_d3cheb=0.0_cp
            !   else
            !      local_d3cheb =( ( 3*(n-1)*spos*yk**2* map_fac1(k)**3) &
            !           &          -( 3*(n-1)**2*cos((n-1)*pos)*sin(pos)*yk* map_fac1(k)**3) &
            !           &          +(n-1)*spos*map_fac1(k)*(1.0-yk**2)&
            !           &            *( map_fac1(k)**2*(1-(n-1)**2) &
            !           &               +3*yk*map_fac2(k)&
            !           &             )&
            !           &          -(3*(n-1)**2*cos((n-1)*pos)*map_fac1(k)*map_fac2(k)*sin(pos)**3) &
            !           &          +((n-1)*spos*map_fac3(k)*sin(pos)**4)&
            !           &        )/sin(pos)**5
            !   end if
            !end if

            d3cheb(n,k)=       two*map_fac3(k)*cheb(n-1,k) + &
                              6.0_cp*map_fac2(k)*dcheb(n-1,k) + &
                             6.0_cp*map_fac1(k)*d2cheb(n-1,k) + &
                                    two*y(k)*d3cheb(n-1,k) - &
                                              d3cheb(n-2,k)
            !if (ABS(local_d3cheb) > 0.0_cp) then
            !write(*,"(A,2I3,3ES20.12,ES11.3)") "Error in d3cheb calculation: ",n,k,&
            !        & d3cheb(n,k),local_d3cheb,d3cheb(n,k)-local_d3cheb,&
            !        &(d3cheb(n,k)-local_d3cheb)/local_d3cheb
            !end if

         end do

      end do

   end subroutine get_chebs_direct
!------------------------------------------------------------------------------
   subroutine cheb_grid(a,b,n,x,y,a1,a2,x0,lbd)
      !
      !   Given the interval [a,b] the routine returns the
      !   n+1 points that should be used to support a
      !   Chebychev expansion. These are the n+1 extrema y(i) of
      !   the Chebychev polynomial of degree n in the
      !   interval [-1,1].
      !   The respective points mapped into the interval of
      !   question [a,b] are the x(i).
      !   NOTE: x(i) and y(i) are stored in the reversed order:
      !    x(1)=b, x(n+1)=a, y(1)=1, y(n+1)=-1
      !
       
      !-- Input variables
      real(cp), intent(in) :: a,b   ! interval boundaries
      integer,  intent(in) :: n ! degree of Cheb polynomial to be represented by the grid points
      real(cp), intent(in) :: a1,a2,x0,lbd

      !-- Output variables
      real(cp), intent(out) :: x(*) ! grid points in interval [a,b]
      real(cp), intent(out) :: y(*) ! grid points in interval [-1,1]
       
      !-- Local variables:
      real(cp) :: bpa,bma
      integer :: k
       
      bma=half*(b-a)
      bpa=half*(a+b)
       
      do k=1,n+1
         y(k)=cos( pi*real(k-1,cp)/real(n,cp) )
         if ( l_newmap ) then
            x(k)=(a2+tan(lbd*(y(k)-x0))/a1)/2 + bpa
         else
            x(k)=bma * y(k) + bpa
         end if
      end do
        
   end subroutine cheb_grid
!------------------------------------------------------------------------------
end module chebyshev_polynoms_mod
