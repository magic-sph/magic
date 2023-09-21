module chebyshev_polynoms_mod

   use precision_mod
   use constants, only: pi, half, one, two, four
   use num_param, only: map_function

   implicit none

   private

   public :: cheb_grid, get_chebs_even

contains

   subroutine get_chebs_even(n_r,a,b,y,n_r_max,cheb,dcheb,d2cheb,dim1,dim2)
      !
      !  Construct even Chebyshev polynomials and their first, second and
      !  third derivative up to degree ``2*(n_r/2)`` at ``(n_r/2)`` points x
      !  in the interval :math:`[a,b]`. Since the polynoms are only defined
      !  in :math:`[-1,1]` we have to map the points x in [a,b] onto points y
      !  in the interval :math:`[-1,1]`. For even Chebs we need only half the
      !  point of the map, these ``(n_r/2)`` points are in the interval :math:`[1,0[`.
      !

      !-- Input variables:
      integer,  intent(in) :: n_r ! number of grid points
                                  ! n_r grid points suffice for a cheb
                                  ! transform up to degree n_r-1
      integer,  intent(in):: n_r_max     ! max number of radial points, dims of y
      real(cp), intent(in) :: a,b        ! interval boundaries :math:`[a,b]`
      real(cp), intent(in) :: y(n_r_max) ! n_r grid points in interval :math:`[a,b]`
      integer,  intent(in) :: dim1,dim2  ! dimensions of cheb,dcheb,......

      !-- Output variables:
      real(cp), intent(out) :: cheb(dim1,dim2)   ! ``cheb(i,j)`` is Chebyshev pol.
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
            dcheb(n,k)=two*map_fac*last_cheb + two*y(k)*last_dcheb - dcheb(n-1,k)
            d2cheb(n,k)=four*map_fac*last_dcheb + two*y(k)*last_d2cheb - d2cheb(n-1,k)

            !-- odd chebs: not stored but necessary for recursion
            last_cheb=two*y(k)*cheb(n,k)-last_cheb
            last_dcheb=two*map_fac*cheb(n,k) + two*y(k)*dcheb(n,k) - last_dcheb
            last_d2cheb=four*map_fac*dcheb(n,k) + two*y(k)*d2cheb(n,k) - last_d2cheb
         end do
      end do

   end subroutine get_chebs_even
!------------------------------------------------------------------------------
   subroutine cheb_grid(a,b,n,x,y,a1,a2,x0,lbd,l_map)
      !
      !   Given the interval :math:`[a,b]` the routine returns the n+1 points
      !   that should be used to support a Chebyshev expansion. These are
      !   the n+1 extrema ``y(i)`` of the Chebyshev polynomial of degree n in the
      !   interval :math:`[-1,1]`. The respective points mapped into the interval of
      !   question :math:`[a,b]` are the ``x(i)``.
      !
      !   .. note:: ``x(i)`` and ``y(i)`` are stored in the reversed order:
      !             ``x(1)=b``, ``x(n+1)=a``, ``y(1)=1``, ``y(n+1)=-1``
      !   ..
      !

      !-- Input variables
      real(cp), intent(in) :: a,b   ! interval boundaries
      integer,  intent(in) :: n ! degree of Cheb polynomial to be represented by the grid points
      real(cp), intent(in) :: a1,a2,x0,lbd
      logical,  intent(in) :: l_map ! Chebyshev mapping

      !-- Output variables
      real(cp), intent(out) :: x(*) ! grid points in interval :math:`[a,b]`
      real(cp), intent(out) :: y(*) ! grid points in interval :math:`[-1,1]`

      !-- Local variables:
      real(cp) :: bpa,bma,AJ,BJ
      integer :: k

      bma=half*(b-a)
      bpa=half*(a+b)
      if ( index(map_function, 'JAFARI') /= 0 .or. index(map_function, 'TT') /= 0 &
      &    .or. index(map_function, 'TEE') /= 0 ) then
         AJ=half*(asinh((one-a2)*a1)+asinh((one+a2)*a1))
         BJ=asinh((one-a2)*a1)
      end if

      do k=1,n+1
         y(k)=cos( pi*real(k-1,cp)/real(n,cp) )
         if ( l_map ) then
            if ( index(map_function, 'TAN') /= 0 .or.    &
            &    index(map_function, 'BAY') /= 0 ) then
               x(k)=bma*(a2+tan(lbd*(y(k)-x0))/a1) + bpa
            else if ( index(map_function, 'ARCSIN') /= 0 .or. &
            &         index(map_function, 'KTL') /= 0 ) then
               x(k)=bma*asin(a1*y(k))/asin(a1)+bpa
            elseif ( index(map_function, 'TT') /= 0 .or.    &
            &        index(map_function, 'TEE') /= 0 ) then
               x(k)=bma*(a2+sinh(AJ*(y(k)-one)+BJ)/a1)+bpa
            else if ( index(map_function, 'JAFARI') /= 0 ) then
               x(k)=bma*(a2+sinh(AJ*(tan(y(k)*atan(x0))/x0-one)+BJ)/a1)+bpa
            end if
         else
            x(k)=bma * y(k) + bpa
         end if
      end do

   end subroutine cheb_grid
!------------------------------------------------------------------------------
end module chebyshev_polynoms_mod
