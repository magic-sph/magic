subroutine cylmean_otc(a,v,n_s_max,n_r_max,n_theta_max,r,s,theta)

   use iso_fortran_env, only: cp => real32

   implicit none

   !-- Input variables
   integer,  intent(in) :: n_r_max
   integer,  intent(in) :: n_theta_max
   integer,  intent(in) :: n_s_max
   real(cp), intent(in) :: r(n_r_max)         ! Spherical radius
   real(cp), intent(in) :: s(n_s_max)         ! Cylindrical radius
   real(cp), intent(in) :: theta(n_theta_max) ! Colatitude
   real(cp), intent(in) :: a(n_theta_max, n_r_max)

   !-- Output variable
   real(cp), intent(out) :: v(n_s_max)

   !-- Local variables
   integer :: n_z, nz, n_s, n_th, n_r, itr
   integer :: n_r0, n_r1, n_r2, n_r3, n_th0, n_th1, n_th2, n_th3
   real(cp) :: zmin, zmax, dz, z, eps, r_cmb, rc, thet, r_icb
   real(cp) :: rr0, rr1, rr2, rr3, r10, r20, r30, r21, r31, r32
   real(cp) :: tt0, tt1, tt2, tt3, t10, t20, t30, t21, t31, t32
   real(cp) :: a01, a12, a23, a012, a123, tot
   real(cp) :: ac(-n_s_max:n_s_max,n_s_max), ait(0:3)
   real(cp), parameter :: pi=acos(-1.0_cp)

   eps=10.0_cp*epsilon(1.0_cp)
   r_cmb=r(1)
   r_icb=r(n_r_max)

   v(:) =0.0_cp

   !-- Loop over axial cylinders starts here
   sLoop: do n_s=1,n_s_max

      if ( s(n_s) < r_icb ) exit sLoop

      zmax = sqrt(r_cmb*r_cmb-s(n_s)*s(n_s)) ! zmax
      zmin = 0.0_cp
      nz = 2*int(n_s_max*(zmax-zmin)/(2.0_cp*r_cmb)) ! Number of z points (one HS)
      nz = max(nz, 4)  ! Minimum to 4 for Simpson integration
      dz = (zmax-zmin)/real(nz,cp) 

      !-- Loop over z starts
      do n_z=-nz,nz
         z=zmin+dz*n_z                           
         rc=sqrt(s(n_s)*s(n_s)+z*z)   ! radius from center
         if (rc >= r_cmb) rc=r_cmb-eps
         thet=0.5_cp*pi-atan(z/s(n_s))  ! polar angle of point (rax,z)
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
         n_th1=n_theta_max
         tbracket: do n_th=n_theta_max,1,-1
            if( theta(n_th) <= thet) then
               n_th1=n_th
               exit tbracket
            end if
         end do tbracket
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

      !
      !  *** interpolation completed
      !
      !  *** simpson integration
      !
      tot=ac(-nz,n_s)+ac(nz,n_s)
      do n_z=-nz+1,nz-1,2
         tot=tot+4.0_cp*ac(n_z,n_s)
      enddo
      do n_z=-nz+2,nz-2,2
         tot=tot+2.0_cp*ac(n_z,n_s)
      enddo
      v(n_s)=tot/(6.0_cp*nz)
   end do sLoop

   !--  special case s=r_cmb
   v(1)=0.5_cp*(a(n_theta_max/2,1)+a(n_theta_max/2+1,1))

end subroutine cylmean_otc

subroutine cylmean_itc(a,vn,vs,n_s_max,n_r_max,n_theta_max,r,s,theta)

   use iso_fortran_env, only: cp => real32

   implicit none

   !-- Input variables
   integer,  intent(in) :: n_r_max
   integer,  intent(in) :: n_theta_max
   integer,  intent(in) :: n_s_max
   real(cp), intent(in) :: r(n_r_max)         ! Spherical radius
   real(cp), intent(in) :: s(n_s_max)         ! Cylindrical radius
   real(cp), intent(in) :: theta(n_theta_max) ! Colatitude
   real(cp), intent(in) :: a(n_theta_max, n_r_max)

   !-- Output variable
   real(cp), intent(out) :: vn(n_s_max)
   real(cp), intent(out) :: vs(n_s_max)

   !-- Local variables
   integer :: n_z, nz, n_s, n_th, n_r, itr, n_hs
   integer :: n_r0, n_r1, n_r2, n_r3, n_th0, n_th1, n_th2, n_th3
   real(cp) :: zmin, zmax, dz, z, eps, r_cmb, r_icb, rc, thet
   real(cp) :: rr0, rr1, rr2, rr3, r10, r20, r30, r21, r31, r32
   real(cp) :: tt0, tt1, tt2, tt3, t10, t20, t30, t21, t31, t32
   real(cp) :: a01, a12, a23, a012, a123, tot1, tot2
   real(cp) :: ac(0:n_s_max,n_s_max,2), ait(0:3)
   real(cp), parameter :: pi=acos(-1.0_cp)

   eps=10.0_cp*epsilon(1.0_cp)
   r_cmb=r(1)
   r_icb=r(n_r_max)

   vn(:)=0.0_cp
   vs(:)=0.0_cp

   !-- Loop over axial cylinders starts here
   sLoop: do n_s=1,n_s_max

      if ( s(n_s) >= r_icb ) cycle
      zmax = sqrt(r_cmb*r_cmb-s(n_s)*s(n_s)) ! zmax
      zmin = sqrt(r_icb*r_icb-s(n_s)*s(n_s))
      nz = 2*int(n_s_max*(zmax-zmin)/(2.0_cp*r_cmb)) ! Number of z points (one HS)
      nz = max(nz, 4)  ! Minimum to 4 for Simpson integration
      dz = (zmax-zmin)/real(nz,cp) 

      !-- Loop over z starts
      do n_z=0,nz
         z=zmin+dz*n_z                           
         rc=sqrt(s(n_s)*s(n_s)+z*z)   ! radius from center
         if (rc >= r_cmb) rc=r_cmb-eps
         if (rc <= r_icb) rc=r_icb+eps
         thet=0.5_cp*pi-atan(z/s(n_s))  ! polar angle of point (rax,z)
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
         if(n_r2 == n_r_max-1) n_r2=n_r_max-2
         if(n_r2 == 1) n_r2=2
         n_r3=n_r2-1
         n_r1=n_r2+1
         n_r0=n_r2+2

         !-- Find indices of angular grid levels that bracket thet
         n_th1=n_theta_max
         tbracket: do n_th=n_theta_max,1,-1
            if( theta(n_th) <= thet) then
               n_th1=n_th
               exit tbracket
            end if
         end do tbracket
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
      !
      !  *** interpolation completed
      !
      !  *** simpson integration
      !
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

end subroutine cylmean_itc

subroutine cylmean(a,v,n_s_max,n_r_max,n_theta_max,r,s,theta)

   use iso_fortran_env, only: cp => real32
   !$ use omp_lib

   implicit none

   !-- Input variables
   integer,  intent(in) :: n_r_max
   integer,  intent(in) :: n_theta_max
   integer,  intent(in) :: n_s_max
   real(cp), intent(in) :: r(n_r_max)         ! Spherical radius
   real(cp), intent(in) :: s(n_s_max)         ! Cylindrical radius
   real(cp), intent(in) :: theta(n_theta_max) ! Colatitude
   real(cp), intent(in) :: a(n_theta_max, n_r_max)

   !-- Output variable
   real(cp), intent(out) :: v(n_s_max)

   !-- Local variables
   integer :: n_z, nz, n_s, n_th, n_r, itr
   integer :: nZstart, nZstop
   integer :: n_r0, n_r1, n_r2, n_r3, n_th0, n_th1, n_th2, n_th3
   real(cp) :: zmin, zmax, dz, z, eps, r_cmb, r_icb, rc, thet
   real(cp) :: rr0, rr1, rr2, rr3, r10, r20, r30, r21, r31, r32
   real(cp) :: tt0, tt1, tt2, tt3, t10, t20, t30, t21, t31, t32
   real(cp) :: a01, a12, a23, a012, a123, tot
   real(cp) :: ac(-n_s_max:n_s_max,n_s_max), ait(0:3)
   real(cp), parameter :: pi=acos(-1.0_cp)

   eps=10.0_cp*epsilon(1.0_cp)
   r_cmb=r(1)
   r_icb=r(n_r_max)

   v(:)=0.0_cp

   !-- Loop over axial cylinders starts here
   !$omp parallel do default(shared) &
   !$omp private(n_s,zmax,zmin,nz,dz,nZstart,nZstop,n_z,z,rc,n_r) &
   !$omp private(thet,n_r0,n_r1,n_r2,n_r3,n_th0,n_th1,n_th2,n_th3)&
   !$omp private(rr0, rr1, rr2, rr3, r10, r20, r30, r21, r31, r32)&
   !$omp private(tt0, tt1, tt2, tt3, t10, t20, t30, t21, t31, t32)&
   !$omp private(a01, a12, a23, a012, a123, tot, itr, ait, n_th)
   sLoop: do n_s=1,n_s_max

      zmax = sqrt(r_cmb*r_cmb-s(n_s)*s(n_s)) ! zmax
      if ( s(n_s) >= r_icb ) then
         zmin = 0.0_cp
      else
         zmin = sqrt(r_icb*r_icb-s(n_s)*s(n_s))
      end if
      nz = 2*int(n_s_max*(zmax-zmin)/(2.0_cp*r_cmb)) ! Number of z points (one HS)
      nz = max(nz, 4)  ! Minimum to 4 for Simpson integration
      dz = (zmax-zmin)/real(nz,cp) 

      !-- Loop over z starts
      if ( s(n_s) >= r_icb ) then
         nZstart = -nz
         nZstop = nz
      else
         nZstart = 0
         nZstop = nz
      end if
      do n_z=nZstart,nZstop
         z=zmin+dz*n_z                           
         rc=sqrt(s(n_s)*s(n_s)+z*z)   ! radius from center
         if (rc >= r_cmb) rc=r_cmb-eps
         if (rc <= r_icb) rc=r_icb+eps
         thet=0.5_cp*pi-atan(z/s(n_s))  ! polar angle of point (rax,z)
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
         if(n_r2 == 1) n_r2=2
         n_r3=n_r2-1
         n_r1=n_r2+1
         n_r0=n_r2+2

         !-- Find indices of angular grid levels that bracket thet
         n_th1=n_theta_max
         tbracket: do n_th=n_theta_max,1,-1
            if( theta(n_th) <= thet) then
               n_th1=n_th
               exit tbracket
            end if
         end do tbracket
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

         !-- Only Northern hemisphere
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
      !
      !  *** interpolation completed
      !
      !  *** simpson integration
      !
      if ( s(n_s) >= r_icb ) then
         tot=ac(-nz,n_s)+ac(nz,n_s)
         do n_z=-nz+1,nz-1,2
            tot=tot+4.0_cp*ac(n_z,n_s)
         enddo
         do n_z=-nz+2,nz-2,2
            tot=tot+2.0_cp*ac(n_z,n_s)
         enddo
         v(n_s)=tot/(6.0_cp*nz)
      else
         tot=ac(0,n_s)+ac(nz,n_s)
         do n_z=1,nz-1,2
            tot=tot+4.0_cp*ac(n_z,n_s)
         enddo
         do n_z=2,nz-2,2
            tot=tot+2.0_cp*ac(n_z,n_s)
         enddo
         v(n_s)=tot/(3.0_cp*nz)
      end if

   end do sLoop
   !$omp end parallel do

   !--  special case s=r_cmb
   v(1)=0.5_cp*(a(n_theta_max/2,1)+a(n_theta_max/2+1,1))

end subroutine cylmean
