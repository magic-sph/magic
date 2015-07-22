!$Id$
module plms_theta

   implicit none

   private

   public :: plm_theta, plm_thetaAS

contains
!*************************************************************************
   subroutine plm_theta(theta,max_degree,max_order,m0, &
                         plma,dtheta_plma,ndim_plma,norm)
      !------------------------------------------------------------------------
      !  The produces the plm for all degrees and orders for a given theta
      !  plus dtheta_plma=sin(theta)* (d plm)/(d theta)
      !  ndim_plma is the dimension of plma and dtheta_plma in the calling routine.
      !  max_degree is the required maximum degree of plm.
      !  The order is as follows:  plma(1)=plm(l=0,m=0),
      !                            plma(2)=plm(l=1,m=0),
      !                            plma(3)=plm(l=2,m=0),
      !                             ................
      !                            plma(max_degree+1)=plm(l=max_degree,m=0),
      !                            plma(max_degree+2)=plm(l=1,m=m0),
      !                            plma(max_degree+3)=plm(l=2,m=m0),
      !                             ................
      !  Same for dtheta_plma !
      !  Norm determins the normalisation: n=0 -- surface normalised,
      !                                    n=1 -- Schmidt normalised,
      !                                    n=2 -- fully normalised.
      !------------------------------------------------------------------------
        
      !-- input variables:
      real(kind=8), intent(in) :: theta ! angle in degrees
      integer,      intent(in) :: max_degree ! required max degree of plm
      integer,      intent(in) :: max_order  ! required max order of plm
      integer,      intent(in) :: m0         ! basic wave number
      integer,      intent(in) :: ndim_plma  ! dimension of plma and dtheta_plma
      integer,      intent(in) :: norm       ! =0 fully normalised
                                             ! =1 Schmidt normalised

      !-- Output variables:
      real(kind=8), intent(out) :: plma(ndim_plma) ! ass. legendres at theta
      real(kind=8), intent(out) :: dtheta_plma(ndim_plma) ! their theta derivative

      !-- Local variables:
      real(kind=8) :: sq2,dnorm,fac,plm,plm1,plm2
      integer :: l,m,j,pos
       
      dnorm=1.d0
      if ( norm == 2 ) dnorm=1.d0/dsqrt(16.d0*datan(1.d0))
      sq2=dsqrt(2.d0)

      !-- calculate plms with recurrence relation, starting with
      !   the known plm(l=m):
      pos=0
      do m=0,max_order,m0
           
         fac=1.d0
         do j=3,2*m+1,2
            fac=fac*dble(j)/dble(j-1)
         end do

         plm=dsqrt(fac)
         if( dsin(theta) /= 0.d0 ) then
            plm=plm*dsin(theta)**m
         else if( m /= 0 ) then
            plm=0.d0
         end if
           
      !-- plm for l=m:
         l=m
         if ( norm == 1 ) then
            dnorm=1.d0/dsqrt(dble(2*l+1))
            if ( m /= 0 ) dnorm=sq2*dnorm
         end if

      !-- Now store it:
         pos=pos+1
         plma(pos) = dnorm*plm
           
         plm1=0.d0
           
      !-- plm for l>m:
         do l=m+1,max_degree
            plm2=plm1
            plm1=plm
            plm= dcos(theta)* dsqrt( dble( (2*l-1)*(2*l+1) ) / &
                               dble( (l-m)*(l+m) )  ) * plm1 - &
                      dsqrt( (dble(2*l+1)*dble(l+m-1)*dble(l-m-1))  / &
                         (dble(2*l-3)*dble(l-m)*dble(l+m))  ) * plm2
            if ( norm == 1 ) then
               dnorm=1.d0/dsqrt(dble(2*l+1))
               if ( m /= 0 ) dnorm=sq2*dnorm
            end if
               
          !----- Now store it:
            pos=pos+1
            if ( pos > ndim_plma ) then
               write(*,*) '! Dimension ndim_plma too small'
               write(*,*) '! in subroutine plm_theta!'
               stop
            end if
            plma(pos) = dnorm*plm
               
         end do

      !-- additional plm(max_degree+1) necessary to calculate
      !   theta derivative, this is stored in dtheta_plma(pos)
      !   to avoid another local array:
         l=max_degree+1
         plm2=plm1
         plm1=plm
         plm= dcos(theta)* dsqrt( dble( (2*l-1)*(2*l+1) ) / &
                            dble( (l-m)*(l+m) )  ) * plm1 - &
                   dsqrt( dble( (2*l+1)*(l+m-1)*(l-m-1) ) / &
                     dble( (2*l-3)*(l-m)*(l+m) ) ) * plm2
         if ( norm == 1 ) then
            dnorm=1.d0/dsqrt(dble(2*l+1))
            if ( m /= 0 ) dnorm=sq2*dnorm
         end if
         dtheta_plma(pos)=dnorm*plm
           
      end do    ! loop over order !
       
      !--  evaluate sin(theta) * theta derivative with recurrence relation
      !    using an l+1 and an l-1 plm:
       
      pos=0

      do m=0,max_order,m0

      !-------- l=m contribution:
         l=m
         pos=pos+1
         if ( pos > ndim_plma ) then
            write(*,*) '! Dimension ndim_plma too small'
            write(*,*) '! in subroutine plm_theta!'
            stop
         end if
         if ( m < max_degree ) then
            if( norm == 0 .OR. norm == 2 ) then
               dtheta_plma(pos)= l/dsqrt(dble(2*l+3)) * plma(pos+1)
            else if ( norm == 1 ) then
               dtheta_plma(pos)= l/dsqrt(dble(2*l+1)) * plma(pos+1)
            end if
         else
            if( norm == 0 .OR. norm == 2 ) then
               dtheta_plma(pos)= l/dsqrt(dble(2*l+3)) * dtheta_plma(pos)
            else if ( norm == 1 ) then
               dtheta_plma(pos)= l/dsqrt(dble(2*l+1)) * dtheta_plma(pos)
            end if
         end if
               
      !-------- l>m contribution:
         do l=m+1,max_degree-1

            pos=pos+1
            if ( pos > ndim_plma ) then
               write(*,*) '! Dimension ndim_plma too small'
               write(*,*) '! in subroutine plm_theta!'
               stop
            end if
            if( norm == 0 .OR. norm == 2 ) then
               dtheta_plma(pos)=                   &
                  l*dsqrt( dble((l+m+1)*(l-m+1)) / &
                             dble((2*l+1)*(2*l+3)) &
                                ) * plma(pos+1)  - &
                  (l+1)*dsqrt( dble((l+m)*(l-m)) / &
                             dble((2*l-1)*(2*l+1)) &
                             ) * plma(pos-1)
            else if ( norm == 1 ) then
               dtheta_plma(pos)=                   &
                    l*dsqrt( dble((l+m+1)*(l-m+1)) &
                                ) * plma(pos+1)  - &
                    (l+1)*dsqrt( dble((l+m)*(l-m)) &
                               ) * plma(pos-1)
               dtheta_plma(pos)=dtheta_plma(pos)/dble(2*l+1)
            end if
                      
         end do ! loop over degree

      !-------- l=max_degree contribution, note
      !         usage of dtheta_plma(pos) instead of plma(pos+1) here:
         if ( m < max_degree ) then
            l=max_degree
            pos=pos+1
            if ( pos > ndim_plma ) then
               write(*,*) '! Dimension ndim_plma too small'
               write(*,*) '! in subroutine plm_theta!'
               stop
            end if
            if( norm == 0 .OR. norm == 2 ) then
               dtheta_plma(pos)=                   &
                  l*dsqrt( dble((l+m+1)*(l-m+1)) / &
                             dble((2*l+1)*(2*l+3)) &
                           ) * dtheta_plma(pos)  - &
                  (l+1)*dsqrt( dble((l+m)*(l-m)) / &
                             dble((2*l-1)*(2*l+1)) &
                             ) * plma(pos-1)
            else if ( norm == 1 ) then
               dtheta_plma(pos)=                   &
                    l*dsqrt( dble((l+m+1)*(l-m+1)) &
                           ) * dtheta_plma(pos)  - &
                    (l+1)*dsqrt( dble((l+m)*(l-m)) &
                               ) * plma(pos-1)
               dtheta_plma(pos)=dtheta_plma(pos)/dble(2*l+1)
            end if
         end if

      end do ! loop over order
     
   end subroutine plm_theta
!------------------------------------------------------------------------
   subroutine plm_thetaAS(theta,max_degree,plma,dtheta_plma,ndim_plma,norm)
      !------------------------------------------------------------------------
      !  The produces the plm for all degrees and order=0 for a given theta
      !  plus dtheta_plma=sin(theta)* (d plm)/(d theta)
      !
      !  Norm determins the normalisation: n=0 -- surface normalised,
      !                                    n=1 -- Schmidt normalised,
      !                                    n=2 -- fully normalised.
      !------------------------------------------------------------------------
       
      !-- Input variables:
      real(kind=8), intent(in) :: theta ! angle in degrees
      integer,      intent(in) :: max_degree ! required max degree of plm
      integer,      intent(in) :: ndim_plma  ! dimension of plma and dtheta_plma
      integer,      intent(in) :: norm       ! =0 fully normalised
                            ! =1 Schmidt normalised

      !-- Output variables:
      real(kind=8), intent(out) :: plma(ndim_plma) ! ass. legendres at theta
      real(kind=8), intent(out) :: dtheta_plma(ndim_plma) ! their theta derivative

      !-- Local variables
      real(kind=8) :: sq2,dnorm,fac,plm,plm1,plm2
      integer :: l,pos
       
      dnorm=1.d0
      if ( norm == 2 ) dnorm=1.d0/dsqrt(16.d0*datan(1.d0))
      sq2=dsqrt(2.d0)

      !-- calculate plms with recurrence relation, starting with
      !   the known plm(l=m):
      pos=0
      fac=1.d0
      plm=dsqrt(fac)
           
      !-- plm for l=m:
      l=0
      if ( norm == 1 ) dnorm=1.d0/dsqrt(dble(2*l+1))

      !-- Now store it:
      pos=pos+1
      plma(pos) = dnorm*plm
           
      plm1=0.d0
           
      !-- plm for l>m:
      do l=1,max_degree
         plm2=plm1
         plm1=plm
         plm= dcos(theta)* dsqrt( dble( (2*l-1)*(2*l+1) ) / &
                             dble( (l)*(l) )  ) * plm1 - &
                    dsqrt( dble( (2*l+1)*(l-1)*(l-1) ) / &
                    dble( (2*l-3)*(l)*(l) ) ) * plm2

         if( norm == 1 ) dnorm=1.d0/dsqrt(dble(2*l+1))
               
         !----- Now store it:
         pos=pos+1
         if ( pos > ndim_plma ) then
            write(*,*) '! Dimension ndim_plma too small'
            write(*,*) '! in subroutine plm_theta!'
            stop
         end if
         plma(pos) = dnorm*plm
           
      end do

      !-- additional plm(max_degree+1) necessary to calculate
      !   theta derivative, this is stored in dtheta_plma(pos)
      !   to avoid another local array:
      l=max_degree+1
      plm2=plm1
      plm1=plm
      plm= dcos(theta)* dsqrt( dble( (2*l-1)*(2*l+1) ) / &
                         dble( (l)*(l) )  ) * plm1 - &
                dsqrt( dble( (2*l+1)*(l-1)*(l-1) ) / &
                  dble( (2*l-3)*(l)*(l) ) ) * plm2

      if( norm == 1 ) dnorm=1.d0/dsqrt(dble(2*l+1))

      dtheta_plma(pos)=dnorm*plm
           
      !--  evaluate sin(theta) * theta derivative with recurrence relation
      !    using an l+1 and an l-1 plm:
      pos=0

      !-------- l=m contribution:
      l=0
      pos=pos+1
      if ( pos > ndim_plma ) then
         write(*,*) '! Dimension ndim_plma too small'
         write(*,*) '! in subroutine plm_theta!'
         stop
      end if
      if( norm == 0 .OR. norm == 2 ) then
         dtheta_plma(pos)= l/dsqrt(dble(2*l+3)) * &
                           plma(pos+1)
      else if ( norm == 1 ) then
         dtheta_plma(pos)= l/dsqrt(dble(2*l+1)) * &
                           plma(pos+1)
      end if
               
      !-------- l>m contribution:
      do l=1,max_degree-1

         pos=pos+1
         if ( pos > ndim_plma ) then
            write(*,*) '! Dimension ndim_plma too small'
            write(*,*) '! in subroutine plm_theta!'
            stop
         end if
         if( norm == 0 .OR. norm == 2 ) then
            dtheta_plma(pos)=  l*dsqrt( dble((l+1)*(l+1)) / &
                            dble((2*l+1)*(2*l+3)) &
                               ) * plma(pos+1)  - &
                 (l+1)*dsqrt( dble((l)*(l)) / &
                            dble((2*l-1)*(2*l+1)) &
                            ) * plma(pos-1)
         else if ( norm == 1 ) then
            dtheta_plma(pos)= l*dsqrt( dble((l+1)*(l+1)) &
                               ) * plma(pos+1)  - &
                   (l+1)*dsqrt( dble((l)*(l)) &
                              ) * plma(pos-1)
            dtheta_plma(pos)=dtheta_plma(pos)/dble(2*l+1)
         end if
              
      end do ! loop over degree

      !-------- l=max_degree contribution, note
      !         usage of dtheta_plma(pos) instead of plma(pos+1) here:
      l=max_degree
      pos=pos+1
      if ( pos > ndim_plma ) then
         write(*,*) '! Dimension ndim_plma too small'
         write(*,*) '! in subroutine plm_theta!'
         stop
      end if
      if( norm == 0 .or. norm == 2 ) then
         dtheta_plma(pos)= l*dsqrt( dble((l+1)*(l+1)) / &
                      dble((2*l+1)*(2*l+3)) &
                     ) * dtheta_plma(pos)  - &
          (l+1)*dsqrt( dble((l)*(l)) / &
                       dble((2*l-1)*(2*l+1)) &
                     ) * plma(pos-1)
      else if ( norm == 1 ) then
         dtheta_plma(pos)= l*dsqrt( dble((l+1)*(l+1)) ) * dtheta_plma(pos)  - &
                           (l+1)*dsqrt( dble((l)*(l)) ) * plma(pos-1)
         dtheta_plma(pos)=dtheta_plma(pos)/dble(2*l+1)
      end if

   end subroutine plm_thetaAS
!------------------------------------------------------------------------
end module plms_theta
