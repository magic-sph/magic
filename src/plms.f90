module plms_theta

   use precision_mod
   use constants, only: osq4pi, one, two
   use useful, only: abortRun

   implicit none

   private

   public :: plm_theta

contains

   subroutine plm_theta(theta,max_degree,min_order,max_order,m0, &
              &          plma,dtheta_plma,ndim_plma,norm)
      !
      !  This produces the :math:`P_{\ell}^m` for all degrees :math:`\ell` and
      !  orders :math:`m` for a given colatitude. :math:`\theta`, as well as
      !  :math:`\sin \theta d P_{\ell}^m /d\theta`. The :math:`P_{\ell}^m` are
      !  stored with a single lm index stored with m first.
      !
      !  Several normalisation are supported:
      !
      !     * n=0 -- surface normalised,
      !     * n=1 -- Schmidt normalised,
      !     * n=2 -- fully normalised.
      !

      !-- Input variables:
      real(cp), intent(in) :: theta      ! angle in radians
      integer,  intent(in) :: max_degree ! required max degree of plm
      integer,  intent(in) :: min_order  ! required min order of plm
      integer,  intent(in) :: max_order  ! required max order of plm
      integer,  intent(in) :: m0         ! basic wave number
      integer,  intent(in) :: ndim_plma  ! dimension of plma and dtheta_plma
      integer,  intent(in) :: norm       ! =0 fully normalised
                                         ! =1 Schmidt normalised

      !-- Output variables:
      real(cp), intent(out) :: plma(ndim_plma) ! associated Legendre polynomials at theta
      real(cp), intent(out) :: dtheta_plma(ndim_plma) ! their theta derivative

      !-- Local variables:
      real(cp) :: sq2,dnorm,fac,plm,plm1,plm2
      integer :: l,m,j,pos

      dnorm=one
      if ( norm == 2 ) dnorm=osq4pi
      sq2=sqrt(two)

      !-- calculate plms with recurrence relation, starting with
      !   the known plm(l=m):
      pos=0
      do m=min_order,max_order,m0

         fac=one
         do j=3,2*m+1,2
            fac=fac*real(j,cp)/real(j-1,cp)
         end do

         plm=sqrt(fac)
         if( sin(theta) /= 0.0_cp ) then
            plm=plm*sin(theta)**m
         else if( m /= 0 ) then
            plm=0.0_cp
         end if

         !-- plm for l=m:
         l=m
         if ( norm == 1 ) then
            dnorm=one/sqrt(real(2*l+1,cp))
            if ( m /= 0 ) dnorm=sq2*dnorm
         end if

         !-- Now store it:
         pos=pos+1
         plma(pos) = dnorm*plm

         plm1=0.0_cp

         !-- plm for l>m:
         do l=m+1,max_degree
            plm2=plm1
            plm1=plm
            plm= cos(theta)* sqrt( real( (2*l-1)*(2*l+1), cp ) /              &
            &                  real( (l-m)*(l+m), cp )  ) * plm1 -            &
            &         sqrt( (real(2*l+1,cp)*real(l+m-1,cp)*real(l-m-1,cp))  / &
            &            (real(2*l-3,cp)*real(l-m,cp)*real(l+m,cp))  ) * plm2
            if ( norm == 1 ) then
               dnorm=one/sqrt(real(2*l+1,cp))
               if ( m /= 0 ) dnorm=sq2*dnorm
            end if

            !----- Now store it:
            pos=pos+1
            if ( pos > ndim_plma ) then
               call abortRun('! Dimension ndim_plma too small in subroutine plm_theta')
            end if
            plma(pos) = dnorm*plm

         end do

         !-- additional plm(max_degree+1) necessary to calculate
         !   theta derivative, this is stored in dtheta_plma(pos)
         !   to avoid another local array:
         l=max_degree+1
         plm2=plm1
         plm1=plm
         plm= cos(theta)* sqrt( real( (2*l-1)*(2*l+1), cp ) /   &
         &                  real( (l-m)*(l+m), cp )  ) * plm1 - &
         &         sqrt( real( (2*l+1)*(l+m-1)*(l-m-1), cp ) /  &
         &           real( (2*l-3)*(l-m)*(l+m),cp ) ) * plm2
         if ( norm == 1 ) then
            dnorm=one/sqrt(real(2*l+1,cp))
            if ( m /= 0 ) dnorm=sq2*dnorm
         end if
         dtheta_plma(pos)=dnorm*plm

      end do    ! loop over order !

      !--  evaluate sin(theta) * theta derivative with recurrence relation
      !    using an l+1 and an l-1 plm:

      pos=0

      do m=min_order,max_order,m0

         !-------- l=m contribution:
         l=m
         pos=pos+1
         if ( pos > ndim_plma ) then
            call abortRun('! Dimension ndim_plma too small in subroutine plm_theta')
         end if
         if ( m < max_degree ) then
            if( norm == 0 .OR. norm == 2 ) then
               dtheta_plma(pos)= l/sqrt(real(2*l+3,cp)) * plma(pos+1)
            else if ( norm == 1 ) then
               dtheta_plma(pos)= l/sqrt(real(2*l+1,cp)) * plma(pos+1)
            end if
         else
            if( norm == 0 .OR. norm == 2 ) then
               dtheta_plma(pos)= l/sqrt(real(2*l+3,cp)) * dtheta_plma(pos)
            else if ( norm == 1 ) then
               dtheta_plma(pos)= l/sqrt(real(2*l+1,cp)) * dtheta_plma(pos)
            end if
         end if

         !-------- l>m contribution:
         do l=m+1,max_degree-1

            pos=pos+1
            if ( pos > ndim_plma ) then
               call abortRun('! Dimension ndim_plma too small in subroutine plm_theta')
            end if
            if( norm == 0 .OR. norm == 2 ) then
               dtheta_plma(pos)=                      &
               &  l*sqrt( real((l+m+1)*(l-m+1),cp) /  &
               &             real((2*l+1)*(2*l+3),cp) &
               &                ) * plma(pos+1)  -    &
               &  (l+1)*sqrt( real((l+m)*(l-m),cp) /  &
               &             real((2*l-1)*(2*l+1),cp) &
               &             ) * plma(pos-1)
            else if ( norm == 1 ) then
               dtheta_plma(pos)=                      &
               &    l*sqrt( real((l+m+1)*(l-m+1),cp)  &
               &                ) * plma(pos+1)  -    &
               &    (l+1)*sqrt( real((l+m)*(l-m),cp)  &
               &               ) * plma(pos-1)
               dtheta_plma(pos)=dtheta_plma(pos)/real(2*l+1,cp)
            end if

         end do ! loop over degree

         !-------- l=max_degree contribution, note
         !         usage of dtheta_plma(pos) instead of plma(pos+1) here:
         if ( m < max_degree ) then
            l=max_degree
            pos=pos+1
            if ( pos > ndim_plma ) then
               call abortRun('! Dimension ndim_plma too small in subroutine plm_theta')
            end if
            if( norm == 0 .OR. norm == 2 ) then
               dtheta_plma(pos)=                      &
               &  l*sqrt( real((l+m+1)*(l-m+1),cp) /  &
               &             real((2*l+1)*(2*l+3),cp) &
               &           ) * dtheta_plma(pos)  -    &
               &  (l+1)*sqrt( real((l+m)*(l-m),cp) /  &
               &             real((2*l-1)*(2*l+1),cp) &
               &             ) * plma(pos-1)
            else if ( norm == 1 ) then
               dtheta_plma(pos)=                      &
               &    l*sqrt( real((l+m+1)*(l-m+1),cp)  &
               &           ) * dtheta_plma(pos)  -    &
               &    (l+1)*sqrt( real((l+m)*(l-m),cp)  &
               &               ) * plma(pos-1)
               dtheta_plma(pos)=dtheta_plma(pos)/real(2*l+1,cp)
            end if
         end if

      end do ! loop over order

   end subroutine plm_theta
!------------------------------------------------------------------------
end module plms_theta
