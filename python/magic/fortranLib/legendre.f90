module legendre

   implicit none

   integer :: n_m_max
   real(kind=8), allocatable :: Plm(:,:)
   real(kind=8), allocatable :: dPlm(:,:)
   real(kind=8), allocatable :: dPhi(:,:)
   real(kind=8), allocatable :: sinTh(:)
   integer, allocatable :: lStart(:)
   integer, allocatable :: lStop(:)
   logical, allocatable :: lmOdd(:)
  
contains

   subroutine initTransform(l_max,minc,lm_max, n_theta_max)

      !-- Input variables
      integer, intent(in) :: l_max
      integer, intent(in) :: minc
      integer, intent(in) :: n_theta_max
      integer, intent(in) :: lm_max

      !-- Local variables:
      integer ::  m_max, lm, lmP, l, m, lmP_max
      integer :: n_theta
      real(kind=8) :: colat, dpi
      real(kind=8) ::  gauss(n_theta_max)
      real(kind=8), allocatable :: plma(:),dtheta_plma(:)

      dpi=4.d0*atan(1.d0)

      m_max = (l_max/minc)*minc
      n_m_max = m_max/minc+1
      lmP_max = lm_max+n_m_max

      if ( allocated(plma) ) deallocate( plma, dtheta_plma)
      allocate(plma(1:lmP_max))
      allocate(dtheta_plma(1:lmP_max))

      if ( allocated(Plm) ) deallocate( Plm, dPlm, dPhi, sinTh)
      allocate(Plm(1:lm_max,1:n_theta_max/2))
      allocate(dPlm(1:lm_max,1:n_theta_max/2))
      allocate(dPhi(1:lm_max,1:n_theta_max/2))
      allocate(sinTh(n_theta_max))

      call gauleg(sinTh,gauss,n_theta_max)

      !-- Get Plm and dPlm
      do n_theta=1,n_theta_max/2  ! Loop over colat in NHS
         colat = sinTh(n_theta)
         call plm_theta(colat,l_max+1,m_max,minc,plma,dtheta_plma,lmP_max)
         lm = 0
         lmP = 0
         do m=0,m_max,minc
            do l=m,l_max
               lm = lm+1
               lmP = lmP+1
               Plm(lm,n_theta)  = (-1.d0)**(real(m,kind=8))*plma(lmP)
               ! True theta derivative !!!
               dPlm(lm,n_theta) = (-1.d0)**(real(m,kind=8))*dtheta_plma(lmP)/sin(colat)
               ! Add the theta dependence in dPhi to simplify the output
               dPhi(lm,n_theta) = real(m,kind=8)/sin(colat)
            end do
            lmP = lmP+1
         end do
      end do

      if ( allocated(lStart) ) deallocate( lStart, lStop, lmOdd )
      allocate(lStart(1:n_m_max))
      allocate(lStop(1:n_m_max))
      allocate(lmOdd(1:n_m_max))

      call getblocks(l_max,minc,n_m_max,lStart,lStop,lmOdd)

      deallocate( plma, dtheta_plma)

   end subroutine initTransform
!------------------------------------------------------------------------------
   subroutine gauleg(theta_ord,gauss,n_theta_max)

      !-- Input variable
      integer, intent(in) :: n_theta_max

      !-- Output variable
      real(kind=8), intent(out) :: theta_ord(n_theta_max),gauss(n_theta_max)

      !--Local variables
      integer ::M,I,J
      real(kind=8) :: dpi,XXM,XXL,EPS,ZZ,ZZ1
      real(kind=8) :: P1,P2,P3,PP

      dpi=4.d0*atan(1.d0)
      M=(n_theta_max+1)/2
      XXM=0.d0
      XXL=1.d0
      EPS=3.D-14

      do I=1,M
         ZZ=cos( dpi*( (real(I,kind=8)-0.25D0)/(real(n_theta_max,kind=8)+0.5D0)) )

         ZZ1=0
         do while (abs(ZZ-ZZ1)>EPS)
            P1=1.d0
            P2=0.d0
            do J=1,n_theta_max
               P3=P2
               P2=P1
               P1=( real(2*J-1,kind=8)*ZZ*P2-real(J-1,kind=8)*P3 )/real(J,kind=8)
            end do
            PP=real(n_theta_max,kind=8)*(ZZ*P1-P2)/(ZZ*ZZ-1.D0)
            ZZ1=ZZ
            ZZ=ZZ1-P1/PP
         end do

         theta_ord(I)=acos(XXM+XXL*ZZ)
         theta_ord(n_theta_max+1-I)=acos(XXM-XXL*ZZ)
         gauss(I)=2.D0*XXL/((1.D0-ZZ*ZZ)*PP*PP)
         gauss(n_theta_max+1-I)=gauss(I)
      end do

   end subroutine gauleg
!------------------------------------------------------------------------------
   subroutine plm_theta(theta,max_degree,max_order,m0,plma,dtheta_plma,ndim_plma)
        
      !-- Input variables
      real(kind=8), intent(in) :: theta
      integer,      intent(in) :: max_degree
      integer,      intent(in) :: max_order
      integer,      intent(in) :: m0
      integer,      intent(in) :: ndim_plma

      !-- Output variables
      real(kind=8), intent(out) :: plma(ndim_plma)
      real(kind=8), intent(out) :: dtheta_plma(ndim_plma)

      !-- Local variables
      real(kind=8) :: sq2,dnorm,fac,plm0,plm1,plm2
      integer :: l,m,j,pos
       
      sq2=sqrt(2.d0)
      dnorm = 1.d0/sqrt(16.d0*atan(1.0d0))

      pos=0
      do m=0,max_order,m0
         ! if (m == 0) then
            ! dnorm=1.d0
         ! else
            ! dnorm=sq2
         ! endif
          
         fac=1.d0
         do j=3,2*m+1,2
            fac=fac*real(j,kind=8)/real(j-1,kind=8)
         end do

         plm0=sqrt(fac)
         if( sin(theta) /= 0.d0 ) then
            plm0=plm0*(-sin(theta))**m
         elseif( m /= 0 ) then
            plm0=0.d0
         endif
          
         l=m
         pos=pos+1
         plma(pos) = dnorm*plm0
          
         plm1=0.d0
          
         do l=m+1,max_degree
            plm2=plm1
            plm1=plm0
            plm0= cos(theta)* sqrt( real( (2*l-1)*(2*l+1), kind=8 ) /        &
            &                      real( (l-m)*(l+m), kind=8 )  ) * plm1 -   &
            &                sqrt( real( (2*l+1)*(l+m-1)*(l-m-1), kind=8 ) / &
            &                      real( (2*l-3)*(l-m)*(l+m), kind=8 ) ) * plm2
             
            pos=pos+1
            plma(pos) = dnorm*plm0
              
         end do

         l=max_degree+1
         plm2=plm1
         plm1=plm0
         plm0= cos(theta)* sqrt( real( (2*l-1)*(2*l+1), kind=8 ) /        &
         &                      real( (l-m)*(l+m), kind=8 )  ) * plm1 -   &
         &                sqrt( real( (2*l+1)*(l+m-1)*(l-m-1), kind=8 ) / &
         &                      real( (2*l-3)*(l-m)*(l+m), kind=8 ) ) * plm2
         dtheta_plma(pos)=dnorm*plm0
      end do    ! loop over order !
       
      pos=0
      do m=0,max_order,m0
         l=m
         pos=pos+1
         if ( m < max_degree ) then
            dtheta_plma(pos)= l/sqrt(real(2*l+3,kind=8)) * plma(pos+1)
         else
            dtheta_plma(pos)= l/sqrt(real(2*l+3,kind=8)) * dtheta_plma(pos)
         end if
              
         do l=m+1,max_degree-1
            pos=pos+1
            dtheta_plma(pos)= l*sqrt( real((l+m+1)*(l-m+1),kind=8) / &
            &                         real((2*l+1)*(2*l+3),kind=8)   &
            &                                     ) * plma(pos+1)  - &
            &                 (l+1)*sqrt( real((l+m)*(l-m),kind=8) / &
            &                         real((2*l-1)*(2*l+1),kind=8)   &
            &                                    ) * plma(pos-1)
         end do ! loop over degree

         if ( m < max_degree ) then
            l=max_degree
            pos=pos+1
            dtheta_plma(pos)= l*sqrt( real((l+m+1)*(l-m+1),kind=8) / &
            &                         real((2*l+1)*(2*l+3),kind=8)   &
            &                       ) * dtheta_plma(pos)  -          &
            &                 (l+1)*sqrt( real((l+m)*(l-m),kind=8) / &
            &                         real((2*l-1)*(2*l+1),kind=8)   &
            &                           ) * plma(pos-1)
         end if
      end do ! loop over order
        
   end subroutine plm_theta
!------------------------------------------------------------------------------
   subroutine specspat_scal(br, n_phi_max, n_theta_max, n_m_max_loc, inputLM, &
              &             lm_max)

      !-- Input variables
      integer, intent(in) :: n_phi_max, n_theta_max
      integer, intent(in) :: lm_max, n_m_max_loc
      complex(kind=4) :: inputLM(lm_max)

      !-- Output variable
      complex(kind=8), intent(out) :: br(n_phi_max,n_theta_max)

      !-- Local variables
      integer :: nThetaNHS,nThetaN,nThetaS,n_m,lm,lms
      complex(kind=8) :: s12,z12

      nThetaNHS=0
      do nThetaN=1,n_theta_max/2
         nThetaS=n_theta_max-nThetaN+1
         nThetaNHS=nThetaNHS+1
         do n_m=1,n_m_max_loc
            lms=lStop(n_m)
            s12=(0.d0,0.d0)
            z12=(0.d0,0.d0)
            do lm=lStart(n_m),lms-1,2
               s12=s12+inputLM(lm)*Plm(lm,nThetaNHS)
               z12=z12+inputLM(lm+1)*Plm(lm+1,nThetaNHS)
            enddo
            if ( lmOdd(n_m) ) then
               s12=s12+inputLM(lms)*Plm(lms,nThetaNHS)
            endif
            br(n_m,nThetaN)=s12+z12
            br(n_m,nThetaS)=s12-z12
         enddo
      enddo

      !-- symmetrize
      if ( n_phi_max > 1 ) then
         br(n_m_max_loc+1:n_phi_max/2+1,:)=(0.d0,0.d0)
         br(n_phi_max/2+2:n_phi_max,:)    =conjg(br(n_phi_max/2:2:-1,:))
      end if

   end subroutine specspat_scal
!------------------------------------------------------------------------------
   subroutine specspat_equat_scal(br, n_phi_max, inputLM, lm_max)

      !-- Input variables
      integer, intent(in) :: n_phi_max
      integer, intent(in) :: lm_max
      complex(kind=4) :: inputLM(lm_max)

      !-- Output variable
      complex(kind=8), intent(out) :: br(n_phi_max)

      !-- Local variables
      integer :: nThetaNHS,n_m,lm,lms
      integer :: shapePlm(2)
      complex(kind=8) :: s12,z12

      shapePlm = shape(Plm)
      nThetaNHS = shapePlm(2)
      do n_m=1,n_m_max
         lms=lStop(n_m)
         s12=(0.d0,0.d0)
         z12=(0.d0,0.d0)
         do lm=lStart(n_m),lms-1,2
            s12=s12+inputLM(lm)*Plm(lm,nThetaNHS)
            z12=z12+inputLM(lm+1)*Plm(lm+1,nThetaNHS)
         enddo
         if ( lmOdd(n_m) ) then
            s12=s12+inputLM(lms)*Plm(lms,nThetaNHS)
         endif
         br(n_m)=s12+z12
      enddo

      !-- symmetrize
      br(n_m_max+1:n_phi_max/2+1)=(0.d0,0.d0)
      br(n_phi_max/2+2:n_phi_max)=conjg(br(n_phi_max/2:2:-1))

   end subroutine specspat_equat_scal
!------------------------------------------------------------------------------
   subroutine specspat_vec(bt, bp, n_phi_max, n_theta_max, n_m_max_loc, dpoldrLM, &
              &            torLM, lm_max)

      !-- Input variables
      integer, intent(in) :: n_phi_max, n_theta_max
      integer, intent(in) :: lm_max, n_m_max_loc
      complex(kind=4) :: dpoldrLM(lm_max)
      complex(kind=4) :: torLM(lm_max)

      !-- Output variables
      complex(kind=8), intent(out) :: bt(n_phi_max,n_theta_max)
      complex(kind=8), intent(out) :: bp(n_phi_max,n_theta_max)

      !-- Local variables
      integer :: nThetaNHS,nThetaN,nThetaS,n_m,lm,lms
      real(kind=8) :: PlmG(lm_max), PlmC(lm_max)
      complex(kind=8) :: vhG(lm_max), vhC(lm_max)
      complex(kind=8) :: vhN1,vhS1,vhN2,vhS2

      do lm=1,lm_max
         vhG(lm)=dpoldrLM(lm)-cmplx(-aimag(torLM(lm)),real(torLM(lm)),kind=8)
         vhC(lm)=dpoldrLM(lm)+cmplx(-aimag(torLM(lm)),real(torLM(lm)),kind=8)
      end do

      nThetaNHS=0
      do nThetaN=1,n_theta_max/2
         nThetaS=n_theta_max-nThetaN+1
         nThetaNHS=nThetaNHS+1
         do n_m=1,n_m_max_loc

            lms=lStop(n_m)
            do lm=lStart(n_m),lms-1,2
               PlmG(lm)=dPlm(lm,nThetaNHS)-dPhi(lm,nThetaNHS)*Plm(lm,nThetaNHS)
               PlmC(lm)=dPlm(lm,nThetaNHS)+dPhi(lm,nThetaNHS)*Plm(lm,nThetaNHS)
               PlmG(lm+1)=dPlm(lm+1,nThetaNHS)-dPhi(lm+1,nThetaNHS)*Plm(lm+1,nThetaNHS)
               PlmC(lm+1)=dPlm(lm+1,nThetaNHS)+dPhi(lm+1,nThetaNHS)*Plm(lm+1,nThetaNHS)
            end do
            if ( lmOdd(n_m) ) then
               PlmG(lms)=dPlm(lms,nThetaNHS)-dPhi(lms,nThetaNHS)*Plm(lms,nThetaNHS)
               PlmC(lms)=dPlm(lms,nThetaNHS)+dPhi(lms,nThetaNHS)*Plm(lms,nThetaNHS)
            endif
         enddo

         do n_m=1,n_m_max_loc
            lms=lStop(n_m)
            vhN1=(0.d0,0.d0)
            vhS1=(0.d0,0.d0)
            vhN2=(0.d0,0.d0)
            vhS2=(0.d0,0.d0)

            do lm=lStart(n_m),lms-1,2
               !-theta-component
               vhN1=vhN1+vhG(lm)*PlmG(lm)+vhG(lm+1)*PlmG(lm+1)
               vhS1=vhS1-vhG(lm)*PlmC(lm)+vhG(lm+1)*PlmC(lm+1)
               vhN2=vhN2+vhC(lm)*PlmC(lm)+vhC(lm+1)*PlmC(lm+1)
               vhS2=vhS2-vhC(lm)*PlmG(lm)+vhC(lm+1)*PlmG(lm+1)
            enddo
            if ( lmOdd(n_m) ) then
               vhN1=vhN1+vhG(lms)*PlmG(lms)
               vhS1=vhS1-vhG(lms)*PlmC(lms)
               vhN2=vhN2+vhC(lms)*PlmC(lms)
               vhS2=vhS2-vhC(lms)*PlmG(lms)
            endif
            bt(n_m,nThetaN)=0.5d0*(vhN1+vhN2)
            bt(n_m,nThetaS)=0.5d0*(vhS1+vhS2)

            bp(n_m,nThetaN)=0.5d0*cmplx(aimag(vhN1-vhN2),-real(vhN1-vhN2), kind=8)
            bp(n_m,nThetaS)=0.5d0*cmplx(aimag(vhS1-vhS2),-real(vhS1-vhS2), kind=8)
         enddo
      enddo

      !-- symmetrize
      if ( n_phi_max > 1 ) then
         bt(n_m_max_loc+1:n_phi_max/2+1,:)=(0.d0,0.d0)
         bt(n_phi_max/2+2:n_phi_max,:)    =conjg(bt(n_phi_max/2:2:-1,:))
         bp(n_m_max_loc+1:n_phi_max/2+1,:)=(0.d0,0.d0)
         bp(n_phi_max/2+2:n_phi_max,:)    =conjg(bp(n_phi_max/2:2:-1,:))
      end if

   end subroutine specspat_vec
!-------------------------------------------------------------------------------
   subroutine specspat_equat_vec(bt, bp, n_phi_max, dpoldrLM, torLM, lm_max)

      !-- Input variables
      integer, intent(in) :: n_phi_max
      integer, intent(in) :: lm_max
      complex(kind=4) :: dpoldrLM(lm_max)
      complex(kind=4) :: torLM(lm_max)

      !-- Output variables
      complex(kind=8), intent(out) :: bt(n_phi_max)
      complex(kind=8), intent(out) :: bp(n_phi_max)

      !-- Local variables
      integer :: nThetaNHS,n_m,lm,lms
      integer :: shapePlm(2)
      real(kind=8) :: PlmG(lm_max), PlmC(lm_max)
      complex(kind=8) :: vhG(lm_max), vhC(lm_max)
      complex(kind=8) :: vhN1,vhS1,vhN2,vhS2

      shapePlm = shape(Plm)
      nThetaNHS = shapePlm(2)

      do lm=1,lm_max
         vhG(lm)=dpoldrLM(lm)-cmplx(-aimag(torLM(lm)),real(torLM(lm)),kind=8)
         vhC(lm)=dpoldrLM(lm)+cmplx(-aimag(torLM(lm)),real(torLM(lm)),kind=8)
      end do

      do n_m=1,n_m_max

         lms=lStop(n_m)
         do lm=lStart(n_m),lms-1,2
            PlmG(lm)=dPlm(lm,nThetaNHS)-dPhi(lm,nThetaNHS)*Plm(lm,nThetaNHS)
            PlmC(lm)=dPlm(lm,nThetaNHS)+dPhi(lm,nThetaNHS)*Plm(lm,nThetaNHS)
            PlmG(lm+1)=dPlm(lm+1,nThetaNHS)-dPhi(lm+1,nThetaNHS)*Plm(lm+1,nThetaNHS)
            PlmC(lm+1)=dPlm(lm+1,nThetaNHS)+dPhi(lm+1,nThetaNHS)*Plm(lm+1,nThetaNHS)
         end do
         if ( lmOdd(n_m) ) then
            PlmG(lms)=dPlm(lms,nThetaNHS)-dPhi(lms,nThetaNHS)*Plm(lms,nThetaNHS)
            PlmC(lms)=dPlm(lms,nThetaNHS)+dPhi(lms,nThetaNHS)*Plm(lms,nThetaNHS)
         endif
      enddo

      do n_m=1,n_m_max
         lms=lStop(n_m)
         vhN1=(0.d0,0.d0)
         vhS1=(0.d0,0.d0)
         vhN2=(0.d0,0.d0)
         vhS2=(0.d0,0.d0)

         do lm=lStart(n_m),lms-1,2
            !-theta-component
            vhN1=vhN1+vhG(lm)*PlmG(lm)+vhG(lm+1)*PlmG(lm+1)
            vhS1=vhS1-vhG(lm)*PlmC(lm)+vhG(lm+1)*PlmC(lm+1)
            vhN2=vhN2+vhC(lm)*PlmC(lm)+vhC(lm+1)*PlmC(lm+1)
            vhS2=vhS2-vhC(lm)*PlmG(lm)+vhC(lm+1)*PlmG(lm+1)
         enddo
         if ( lmOdd(n_m) ) then
            vhN1=vhN1+vhG(lms)*PlmG(lms)
            vhS1=vhS1-vhG(lms)*PlmC(lms)
            vhN2=vhN2+vhC(lms)*PlmC(lms)
            vhS2=vhS2-vhC(lms)*PlmG(lms)
         endif
         bt(n_m)=0.5d0*(vhN1+vhN2)
         bt(n_m)=0.5d0*(vhS1+vhS2)

         bp(n_m)=0.5d0*cmplx(aimag(vhN1-vhN2),-real(vhN1-vhN2), kind=8)
         bp(n_m)=0.5d0*cmplx(aimag(vhS1-vhS2),-real(vhS1-vhS2), kind=8)
      enddo

      bt(n_m_max+1:n_phi_max/2+1)=(0.d0,0.d0)
      bt(n_phi_max/2+2:n_phi_max)=conjg(bt(n_phi_max/2:2:-1))
      bp(n_m_max+1:n_phi_max/2+1)=(0.d0,0.d0)
      bp(n_phi_max/2+2:n_phi_max)=conjg(bp(n_phi_max/2:2:-1))

   end subroutine specspat_equat_vec
!-------------------------------------------------------------------------------
   subroutine getblocks(l_max,minc,n_m_max,lStart,lStop,lmOdd)

      !-- Input variables
      integer, intent(in) :: l_max
      integer, intent(in) :: minc
      integer, intent(in) :: n_m_max

      !-- Output variables
      integer, intent(out) :: lStart(n_m_max),lStop(n_m_max)
      logical, intent(out) :: lmOdd(n_m_max)


      !-- Local variables
      integer :: mc, m

      lStart(1) =1
      lStop(1)  =l_max+1
      if ( mod(l_max,2) == 0 ) then
         lmOdd(1) =.true.
      else
         lmOdd(1) =.false.
      end if
      do mc=2,n_m_max
         m=(mc-1)*minc
         lStart(mc) =lStop(mc-1) +1
         lStop(mc)  =lStart(mc)  +l_max-m
         if ( mod(lStop(mc)-lStart(mc),2) == 0 ) then
            lmOdd(mc) =.true.
         else
            lmOdd(mc) =.false.
         end if
      end do

   end subroutine getblocks
!------------------------------------------------------------------------------
end module legendre
