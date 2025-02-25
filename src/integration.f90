module integration
   !
   ! Radial integration functions
   !

   use precision_mod
   use constants, only: half, one, two, pi
   use radial_scheme, only: type_rscheme
   use cosine_transform_odd

   implicit none

   private

   public :: rIntIC, rInt_R, simps, cylmean_otc, cylmean_itc

contains

   real(cp) function rIntIC(f,nRmax,drFac,chebt)
      !
      !  This function performs the radial integral over a
      !  function f that is given on the appropriate nRmax
      !  radial Chebychev grid points.
      !

      !-- Input variables:
      integer,           intent(in) :: nRmax
      real(cp),          intent(inout) :: f(nRmax)
      real(cp),          intent(in) :: drFac
      type(costf_odd_t), intent(in) :: chebt

      !-- Local variables:
      integer :: nCheb,nChebInt
      real(cp) :: weight
      real(cp) :: work(nRmax)

      call chebt%costf1(f,work)
      f(1)    =half*f(1)
      f(nRmax)=half*f(nRmax)

      !-- Sum contribution:
      rIntIC=f(1)           ! This is zero order contribution
      do nCheb=2,nRmax      ! Only even chebs for IC
         nChebInt=2*nCheb-1
         weight  =-one/real(nChebInt*(nChebInt-2),cp)
         rIntIC  =rIntIC+weight*f(nCheb)
      end do

      !-- Remaining renormalisation:
      rIntIC=sqrt(two/real(nRmax-1,cp))*rIntIC/drFac

   end function rIntIC
!------------------------------------------------------------------------------
   real(cp) function rInt_R(f,r,r_scheme) result(rInt)
      !
      !   Same as function rInt but for a radial dependent mapping function
      !   dr_fac2.
      !

      !-- Input variables:
      real(cp),            intent(in) :: f(:)    ! Input function
      real(cp),            intent(in) :: r(:)    ! Radius
      class(type_rscheme), intent(in) :: r_scheme! Radial scheme (FD or Cheb)

      !-- Local variables
      real(cp), allocatable :: f2(:)
      integer :: nCheb, nRmax


      !--- Integrals:
      if ( r_scheme%version == 'cheb' ) then

         nRmax=size(f)
         allocate( f2(nRmax) )
         f2(:)=f(:)/r_scheme%drx(:)

         !-- Transform to cheb space:
         call r_scheme%costf1(f2)
         f2(1)    =half*f2(1)
         f2(nRmax)=half*f2(nRmax)

         !-- Sum contribution:
         rInt=f2(1)            ! This is zero order contribution
         !do nCheb=3,r_scheme%n_max,2 ! Only even chebs contribute
         do nCheb=3,nRmax,2 ! Only even chebs contribute
            rInt=rInt-one/real(nCheb*(nCheb-2),cp)*f2(nCheb)
         end do

         !-- Remaining renormalisation:
         rInt=two*sqrt(two/real(nRmax-1,cp))*rInt

         deallocate( f2 )

      else

         rInt=simps(f,r)

      end if

   end function rInt_R
!------------------------------------------------------------------------------
   real(cp) function simps(f,r) result(rInt)
      !
      ! Simpson's method to integrate a function
      !

      !-- Input variables:
      real(cp), intent(in) :: f(:)    ! Input function
      real(cp), intent(in) :: r(:)    ! Radius

      !-- Local variables
      real(cp) :: h1, h2
      integer :: n_r, n_r_max

      n_r_max=size(f)

      if ( mod(n_r_max,2)==1 ) then ! Odd number (Simpson ok)

         rInt = 0.0_cp
         do n_r=2,n_r_max-1,2
            h2=r(n_r+1)-r(n_r)
            h1=r(n_r)-r(n_r-1)
            rInt=rInt+(h1+h2)/6.0_cp*( f(n_r-1)*(two*h1-h2)/h1     +&
            &                      f(n_r)  *(h1+h2)*(h1+h2)/(h1*h2)+&
            &                      f(n_r+1)*(two*h2-h1)/h2 )
         end do

         rInt = -rInt

      else ! Even number (twice simpson + trapz on the first and last points)

         rInt = half*(r(2)-r(1))*(f(2)+f(1))
         do n_r=3,n_r_max,2
            h2=r(n_r+1)-r(n_r)
            h1=r(n_r)-r(n_r-1)
            rInt=rInt+(h1+h2)/6.0_cp*( f(n_r-1)*(two*h1-h2)/h1     +&
            &                      f(n_r)  *(h1+h2)*(h1+h2)/(h1*h2)+&
            &                      f(n_r+1)*(two*h2-h1)/h2 )
         end do
         rInt = rInt+half*(r(n_r_max)-r(n_r_max-1))*(f(n_r_max)+f(n_r_max-1))
         do n_r=2,n_r_max-1,2
            h2=r(n_r+1)-r(n_r)
            h1=r(n_r)-r(n_r-1)
            rInt=rInt+(h1+h2)/6.0_cp*( f(n_r-1)*(two*h1-h2)/h1     +&
            &                      f(n_r)  *(h1+h2)*(h1+h2)/(h1*h2)+&
            &                      f(n_r+1)*(two*h2-h1)/h2 )
         end do
         rInt = -half*rInt

      end if

   end function simps
!------------------------------------------------------------------------------
   subroutine cylmean_otc(a,v,n_s_max,n_s_otc,r,s,theta,zDensIn)
      !
      ! This routine computes a z-averaging using Simpsons's rule outside T.C.
      !

      !-- Input variables
      integer,  intent(in) :: n_s_max
      integer,  intent(in) :: n_s_otc
      real(cp), intent(in) :: r(:)         ! Spherical radius
      real(cp), intent(in) :: s(n_s_max)   ! Cylindrical radius
      real(cp), intent(in) :: theta(:)     ! Colatitude
      real(cp), intent(in) :: a(:,:)
      real(cp), optional, intent(in) :: zDensIn

      !-- Output variable
      real(cp), intent(out) :: v(n_s_max)

      !-- Local variables
      integer :: n_z, nz, n_s, n_th, n_r, itr, n_r_max, n_theta_max
      integer :: n_r0, n_r1, n_r2, n_r3, n_th0, n_th1, n_th2, n_th3
      real(cp) :: zmin, zmax, dz, z, eps, r_cmb, rc, thet, r_icb
      real(cp) :: rr0, rr1, rr2, rr3, r10, r20, r30, r21, r31, r32
      real(cp) :: tt0, tt1, tt2, tt3, t10, t20, t30, t21, t31, t32
      real(cp) :: a01, a12, a23, a012, a123, tot, zDens
      real(cp) :: ait(0:3)
      real(cp), allocatable :: ac(:,:)

      n_r_max = size(r)
      n_theta_max = size(theta)
      eps=10.0_cp*epsilon(1.0_cp)
      r_cmb=r(1)
      r_icb=r(n_r_max)

      if ( present(zDensIn) ) then
         zDens=zDensIn
         zmax = sqrt(r_cmb*r_cmb-r_icb*r_icb) ! zmax
         zmin = 0.0_cp
         nz = 2*int(n_s_max*zDens*(zmax-zmin)/(2.0_cp*r_cmb)) ! Number of z points (one HS)
         allocate( ac(-nz:nz,n_s_max) )
      else
         zDens=1.0_cp
         allocate( ac(-n_s_max:n_s_max,n_s_max) )
      end if

      v(:) =0.0_cp

      !-- Loop over axial cylinders starts here
      !$omp parallel do default(shared) &
      !$omp private(n_s,zmax,zmin,nz,dz,n_z,z,rc,thet,n_r2,n_r,n_r3,n_r1,n_r0) &
      !$omp private(n_th,n_th1,n_th2,n_th3,n_th0,rr0,rr1,rr2,rr3,r10,r20,r30)  &
      !$omp private(r21,r31,r32,tt0,tt1,tt2,tt3,t10,t20,t30,t21,t31,t32,itr)   &
      !$omp private(a01,a12,a23,a012,a123,tot,ait)
      sLoop: do n_s=1,n_s_otc

         zmax = sqrt(r_cmb*r_cmb-s(n_s)*s(n_s)) ! zmax
         zmin = 0.0_cp
         nz = 2*int(n_s_max*(zmax-zmin)/(2.0_cp*r_cmb)) ! Number of z points (one HS)
         nz = max(int(zDens*nz), 4)  ! Minimum to 4 for Simpson integration
         dz = (zmax-zmin)/real(nz,cp)

         !-- Loop over z starts
         do n_z=-nz,nz
            z=zmin+dz*n_z
            rc=sqrt(s(n_s)*s(n_s)+z*z)   ! radius from center
            if (rc >= r_cmb) rc=r_cmb-eps
            if ( s(n_s)==0.0_cp ) then
               thet=0.0_cp
            else
               if (z > rc ) then
                  thet=0.0_cp
               else
                  thet=acos(z/rc)
               end if
            end if
            ac(n_z,n_s)=0.0_cp
            !
            !  **** Interpolate values from (theta,r)-grid onto equidistant
            !  **** (z,rax)-grid using a fourth-order Lagrangian scheme
            !
            !--  Find indices of radial grid levels that bracket rc
            n_r2=n_r_max-1
            rbracket: do n_r=n_r_max-1,1,-1
               if ( r(n_r) >= rc ) then
                  n_r2 = n_r
                  exit rbracket
               end if
            end do rbracket
            if(n_r2 == n_r_max-1) n_r2=n_r_max-2
            if(n_r2 == 1 ) n_r2=2
            n_r3=n_r2-1
            n_r1=n_r2+1
            n_r0=n_r2+2

            !-- Find indices of angular grid levels that bracket thet
            if ( thet < theta(1) ) then
               n_th1=1
            else
               n_th1=n_theta_max
               tbracket: do n_th=n_theta_max,1,-1
                  if( theta(n_th) <= thet) then
                     n_th1=n_th
                     exit tbracket
                  end if
               end do tbracket
            end if
            if ( n_th1 == n_theta_max ) n_th1=n_theta_max-2
            if ( n_th1 == n_theta_max-1 ) n_th1=n_theta_max-2
            if ( n_th1 == 1 ) n_th1=2
            n_th2=n_th1+1
            n_th3=n_th1+2
            n_th0=n_th1-1

            !--  Calculate differences in r for 4th-order interpolation
            rr0=rc-r(n_r0)
            rr1=rc-r(n_r1)
            rr2=rc-r(n_r2)
            rr3=rc-r(n_r3)
            r10= 1.0_cp/(r(n_r1)-r(n_r0))
            r20= 1.0_cp/(r(n_r2)-r(n_r0))
            r30= 1.0_cp/(r(n_r3)-r(n_r0))
            r21= 1.0_cp/(r(n_r2)-r(n_r1))
            r31= 1.0_cp/(r(n_r3)-r(n_r1))
            r32= 1.0_cp/(r(n_r3)-r(n_r2))

            !--  Calculate differences in theta for 4th-order interpolation
            tt0=thet-theta(n_th0)
            tt1=thet-theta(n_th1)
            tt2=thet-theta(n_th2)
            tt3=thet-theta(n_th3)
            t10=1.0_cp/(theta(n_th1)-theta(n_th0))
            t20=1.0_cp/(theta(n_th2)-theta(n_th0))
            t30=1.0_cp/(theta(n_th3)-theta(n_th0))
            t21=1.0_cp/(theta(n_th2)-theta(n_th1))
            t31=1.0_cp/(theta(n_th3)-theta(n_th1))
            t32=1.0_cp/(theta(n_th3)-theta(n_th2))

            !-- Loop over 4 neighboring grid angles
            do itr=0,3
               n_th=n_th0+itr

               !-- Interpolation in r-direction
               a01=(rr0*a(n_th,n_r1)-rr1*a(n_th,n_r0))*r10
               a12=(rr1*a(n_th,n_r2)-rr2*a(n_th,n_r1))*r21
               a23=(rr2*a(n_th,n_r3)-rr3*a(n_th,n_r2))*r32

               a012=(rr0*a12-rr2*a01)*r20
               a123=(rr1*a23-rr3*a12)*r31

               ait(itr)=(rr0*a123-rr3*a012)*r30
            end do

            !-- Interpolation in theta-direction
            a01=(tt0*ait(1)-tt1*ait(0))*t10
            a12=(tt1*ait(2)-tt2*ait(1))*t21
            a23=(tt2*ait(3)-tt3*ait(2))*t32

            a012=(tt0*a12-tt2*a01)*t20
            a123=(tt1*a23-tt3*a12)*t31
            ac(n_z,n_s)=ac(n_z,n_s)+(tt0*a123-tt3*a012)*t30
         end do

         !-- Simpson integration
         tot=ac(-nz,n_s)+ac(nz,n_s)
         do n_z=-nz+1,nz-1,2
            tot=tot+4.0_cp*ac(n_z,n_s)
         enddo
         do n_z=-nz+2,nz-2,2
            tot=tot+2.0_cp*ac(n_z,n_s)
         enddo
         v(n_s)=tot/(6.0_cp*nz)
      end do sLoop
      !$omp end parallel do

      !-- Equatorial point
      v(1)=0.5_cp*(a(n_theta_max/2,1)+a(n_theta_max/2+1,1))

      deallocate( ac )

   end subroutine cylmean_otc
!------------------------------------------------------------------------------------
   subroutine cylmean_itc(a,vn,vs,n_s_max,n_s_otc,r,s,theta, zDensIn)
      !
      ! This routine computes a z-averaging using Simpsons's rule inside T.C.
      !

      !-- Input variables
      integer,  intent(in) :: n_s_max
      integer,  intent(in) :: n_s_otc
      real(cp), intent(in) :: r(:)         ! Spherical radius
      real(cp), intent(in) :: s(n_s_max)         ! Cylindrical radius
      real(cp), intent(in) :: theta(:) ! Colatitude
      real(cp), intent(in) :: a(:,:)
      real(cp), optional, intent(in) :: zDensIn

      !-- Output variable
      real(cp), intent(out) :: vn(n_s_max)
      real(cp), intent(out) :: vs(n_s_max)

      !-- Local variables
      integer :: n_z, nz, n_s, n_th, n_r, itr, n_hs, n_r_max, n_theta_max
      integer :: n_r0, n_r1, n_r2, n_r3, n_th0, n_th1, n_th2, n_th3
      real(cp) :: zmin, zmax, dz, z, eps, r_cmb, r_icb, rc, thet, zDens
      real(cp) :: rr0, rr1, rr2, rr3, r10, r20, r30, r21, r31, r32
      real(cp) :: tt0, tt1, tt2, tt3, t10, t20, t30, t21, t31, t32
      real(cp) :: a01, a12, a23, a012, a123, tot1, tot2
      real(cp) ::  ait(0:3)
      real(cp), allocatable :: ac(:,:,:)

      n_r_max = size(r)
      n_theta_max = size(theta)
      eps=10.0_cp*epsilon(1.0_cp)
      r_cmb=r(1)
      r_icb=r(size(r))

      if ( present(zDensIn) ) then
         zDens=zDensIn
         zmax = sqrt(r_cmb*r_cmb-r_icb*r_icb) ! zmax
         zmin = 0.0_cp
         nz = 2*int(n_s_max*zDens*(zmax-zmin)/(2.0_cp*r_cmb)) ! Number of z points (one HS)
         allocate(ac(0:nz,n_s_max,2))
      else
         zDens=1.0_cp
         allocate(ac(0:n_s_max,n_s_max,2))
      end if

      vn(:)=0.0_cp
      vs(:)=0.0_cp
      ac(:,:,:)=0.0_cp

      !-- Loop over axial cylinders starts here
      !$omp parallel do default(shared) &
      !$omp private(n_s,zmax,zmin,nz,dz,n_z,z,rc,thet,n_r2,n_r,n_r3,n_r1,n_r0) &
      !$omp private(n_th,n_th1,n_th2,n_th3,n_th0,rr0,rr1,rr2,rr3,r10,r20,r30)  &
      !$omp private(r21,r31,r32,tt0,tt1,tt2,tt3,t10,t20,t30,t21,t31,t32,itr)   &
      !$omp private(a01,a12,a23,a012,a123,tot1,tot2,n_hs,ait)
      sLoop: do n_s=n_s_otc+1,n_s_max

         zmax = sqrt(r_cmb*r_cmb-s(n_s)*s(n_s)) ! zmax
         zmin = sqrt(r_icb*r_icb-s(n_s)*s(n_s))
         nz = 2*int(n_s_max*(zmax-zmin)/(2.0_cp*r_cmb)) ! Number of z points (one HS)
         nz = max(int(zDens*nz), 4)  ! Minimum to 4 for Simpson integration
         dz = (zmax-zmin)/real(nz,cp)

         !-- Loop over z starts
         do n_z=0,nz
            z=zmin+dz*n_z
            rc=sqrt(s(n_s)*s(n_s)+z*z)   ! radius from center
            if (rc >= r_cmb) rc=r_cmb-eps
            if (rc <= r_icb) rc=r_icb+eps
            if ( s(n_s)==0.0_cp ) then
               thet=0.0_cp
            else
               if (z > rc ) then
                  thet=0.0_cp
               else
                  thet=acos(z/rc)
               end if
            end if
            ac(n_z,n_s,1)=0.0_cp
            ac(n_z,n_s,2)=0.0_cp
            !
            !  **** Interpolate values from (theta,r)-grid onto equidistant
            !  **** (z,rax)-grid using a fourth-order Lagrangian scheme
            !
            !--  Find indices of radial grid levels that bracket rc
            n_r2=n_r_max-1
            rbracket: do n_r=n_r_max-1,1,-1
               if ( r(n_r) >= rc ) then
                  n_r2 = n_r
                  exit rbracket
               end if
            end do rbracket
            if (n_r2 == n_r_max-1) n_r2=n_r_max-2
            if (n_r2 == 1) n_r2=2
            n_r3=n_r2-1
            n_r1=n_r2+1
            n_r0=n_r2+2

            !-- Find indices of angular grid levels that bracket thet
            if ( thet < theta(1) ) then
               n_th1=1
            else
               n_th1=n_theta_max
               tbracket: do n_th=n_theta_max,1,-1
                  if( theta(n_th) <= thet) then
                     n_th1=n_th
                     exit tbracket
                  end if
               end do tbracket
            end if
            if ( n_th1 == n_theta_max ) n_th1=n_theta_max-2
            if ( n_th1 == n_theta_max-1 ) n_th1=n_theta_max-2
            if ( n_th1 == 1 ) n_th1=2
            n_th2=n_th1+1
            n_th3=n_th1+2
            n_th0=n_th1-1

            !--  Calculate differences in r for 4th-order interpolation
            rr0=rc-r(n_r0)
            rr1=rc-r(n_r1)
            rr2=rc-r(n_r2)
            rr3=rc-r(n_r3)
            r10= 1.0_cp/(r(n_r1)-r(n_r0))
            r20= 1.0_cp/(r(n_r2)-r(n_r0))
            r30= 1.0_cp/(r(n_r3)-r(n_r0))
            r21= 1.0_cp/(r(n_r2)-r(n_r1))
            r31= 1.0_cp/(r(n_r3)-r(n_r1))
            r32= 1.0_cp/(r(n_r3)-r(n_r2))

            !--  Calculate differences in theta for 4th-order interpolation
            tt0=thet-theta(n_th0)
            tt1=thet-theta(n_th1)
            tt2=thet-theta(n_th2)
            tt3=thet-theta(n_th3)
            t10=1.0_cp/(theta(n_th1)-theta(n_th0))
            t20=1.0_cp/(theta(n_th2)-theta(n_th0))
            t30=1.0_cp/(theta(n_th3)-theta(n_th0))
            t21=1.0_cp/(theta(n_th2)-theta(n_th1))
            t31=1.0_cp/(theta(n_th3)-theta(n_th1))
            t32=1.0_cp/(theta(n_th3)-theta(n_th2))

            !-- Loop over North/Sooth
            do n_hs=1,2

               !-- Loop over 4 neighboring grid angles
               do itr=0,3
                  n_th=n_th0+itr

                  if ( n_hs == 2 ) n_th = n_theta_max+1-n_th ! Southern hemisphere

                  !-- Interpolation in r-direction
                  a01=(rr0*a(n_th,n_r1)-rr1*a(n_th,n_r0))*r10
                  a12=(rr1*a(n_th,n_r2)-rr2*a(n_th,n_r1))*r21
                  a23=(rr2*a(n_th,n_r3)-rr3*a(n_th,n_r2))*r32

                  a012=(rr0*a12-rr2*a01)*r20
                  a123=(rr1*a23-rr3*a12)*r31

                  ait(itr)=(rr0*a123-rr3*a012)*r30
               end do

               !-- Interpolation in theta-direction
               a01=(tt0*ait(1)-tt1*ait(0))*t10
               a12=(tt1*ait(2)-tt2*ait(1))*t21
               a23=(tt2*ait(3)-tt3*ait(2))*t32

               a012=(tt0*a12-tt2*a01)*t20
               a123=(tt1*a23-tt3*a12)*t31
               ac(n_z,n_s,n_hs)=ac(n_z,n_s,n_hs)+(tt0*a123-tt3*a012)*t30
            end do

         end do

         !call interp_theta(a(:,1),ac(nz,n_s,1:2),r_cmb,s(n_s),theta)
         !if ( n_s == n_s_max ) then
         !   if ( abs(r_cmb*sin(theta(1))-s(n_s)) <= eps ) then
         !      ac(nz,n_s,1)=a(1,1)
         !      ac(nz,n_s,2)=a(n_theta_max,1)
         !   end if
         !end if

         !-- Simpson integration
         tot1=ac(0,n_s,1)+ac(nz,n_s,1)
         tot2=ac(0,n_s,2)+ac(nz,n_s,2)
         do n_z=1,nz-1,2
            tot1=tot1+4.0_cp*ac(n_z,n_s,1)
            tot2=tot2+4.0_cp*ac(n_z,n_s,2)
         enddo
         do n_z=2,nz-2,2
            tot1=tot1+2.0_cp*ac(n_z,n_s,1)
            tot2=tot2+2.0_cp*ac(n_z,n_s,2)
         enddo
         vn(n_s)=tot1/(3.0_cp*nz)
         vs(n_s)=tot2/(3.0_cp*nz)
      end do sLoop
      !$omp end parallel do

      deallocate(ac)

   end subroutine cylmean_itc
!------------------------------------------------------------------------------------
end module integration
