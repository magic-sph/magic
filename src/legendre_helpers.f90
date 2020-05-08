module leg_helper_mod

   use precision_mod
   use truncation, only: l_axi
   use constants, only: zero, one, two
   
   implicit none

   private

   public :: legPrep, legPrep_IC

contains

   subroutine legPrep(w,dw,ddw,z,dz,dLh,lm_max,l_max,minc,r,lDeriv,lHor,dLhw, &
              &       vhG,vhC,dLhz,cvhG,cvhC)
      !
      !  Purpose of this subroutine is to prepare Legendre transforms     
      !  from (r,l,m) space to (r,theta,m) space by calculating           
      !  auxiliary arrays w, dw, ddw, ....... which contain               
      !  combinations of harmonic coeffs multiplied with (l,m)-dependend  
      !  factors as well as possible the radial dependencies.
      !
      !  lDeriv=.true. field derivatives required  for curl of field     
      !
    
      !-- Input variables:
      integer,     intent(in) :: lm_max
      complex(cp), intent(in) :: w(lm_max)
      complex(cp), intent(in) :: dw(lm_max)
      complex(cp), intent(in) :: ddw(lm_max)
      complex(cp), intent(in) :: z(lm_max)
      complex(cp), intent(in) :: dz(lm_max)
      real(cp),    intent(in) :: dLh(lm_max)
      integer,     intent(in) :: l_max,minc
      real(cp),    intent(in) :: r
      logical,     intent(in) :: lHor
      logical,     intent(in) :: lDeriv

      !-- Output variable:
      complex(cp), intent(out) :: dLhw(lm_max),dLhz(lm_max)
      complex(cp), intent(out) :: vhG(lm_max),vhC(lm_max)
      complex(cp), intent(out) :: cvhG(lm_max),cvhC(lm_max)

      !-- Local variables:
      integer :: lm,l,m
      real(cp) :: Or_e2
      complex(cp) :: help


      lm=0
      do m=0,l_max,minc
         do l=m,l_max
            lm=lm+1
            dLhw(lm)=dLh(lm)*w(lm)
            if ( lHor ) then
               vhG(lm) =dw(lm)-cmplx(-aimag(z(lm)),real(z(lm)),kind=cp)
               vhC(lm) =dw(lm)+cmplx(-aimag(z(lm)),real(z(lm)),kind=cp)
            end if
         end do
      end do
      dLhw(1)=zero
      if ( lHor ) then
         vhG(1) =zero
         vhC(1) =zero
      end if

      if ( lDeriv ) then
         Or_e2=one/r**2
         lm=0
         do m=0,l_max,minc
            do l=m,l_max
               lm=lm+1
               dLhz(lm)=dLh(lm)*z(lm)
               help=dLh(lm)*Or_e2*w(lm)-ddw(lm)
               cvhG(lm)=dz(lm)-cmplx(-aimag(help),real(help),kind=cp)
               cvhC(lm)=dz(lm)+cmplx(-aimag(help),real(help),kind=cp)
            end do
         end do
         dLhz(1)=zero
         cvhG(1)=zero
         cvhc(1)=zero
      end if

   end subroutine legPrep
!------------------------------------------------------------------------------
   subroutine legPrep_IC(w,dw,ddw,z,dz,dLh,lm_max,l_max,minc,r,r_ICB,lDeriv, &
              &          lHor,lCondIC,dLhw,vhG,vhC,dLhz,cvhG,cvhC)
      !
      !  Purpose of this subroutine is to prepare Legendre transforms     
      !  from (r,l,m) space to (r,theta,m) space by calculating           
      !  auxiliary arrays dLhw,vhG, ......... which contain               
      !  combinations of harmonic coeffs multiplied with (l,m)-dependend  
      !  factors as well as possible the radial dependencies.             
      !
      !  lDeriv=.true. field derivatives required  for curl of field     
      !
      !  Note: This routine is used for the inner core magnetic field     
      !  which has a special radial function ansatz. It can also be       
      !  used to prepare the calculation of a field in an insulating      
      !  inner core for lCondIC=.false.. For this the w has to be the     
      !  outer core poloidal field and nR is the grid point for the ICB.  
      !  In any case legTF can be used for the following Legendre         
      !  transform and fftJW for the Fourier transform.                   
      !

      !-- Input variables:
      integer,     intent(in) :: lm_max
      complex(cp), intent(in) :: w(lm_max)
      complex(cp), intent(in) :: dw(lm_max)
      complex(cp), intent(in) :: ddw(lm_max)
      complex(cp), intent(in) :: z(lm_max)
      complex(cp), intent(in) :: dz(lm_max)
      real(cp),    intent(in) :: dLh(lm_max)
      integer ,    intent(in) :: l_max,minc
      real(cp),    intent(in) :: r,r_ICB
      logical,     intent(in) :: lHor
      logical,     intent(in) :: lDeriv
      logical,     intent(in) :: lCondIC

      !-- Output variables:
      complex(cp), intent(out) :: dLhw(lm_max),dLhz(lm_max)
      complex(cp), intent(out) :: vhG(lm_max),vhC(lm_max)
      complex(cp), intent(out) :: cvhG(lm_max),cvhC(lm_max)

      !-- Local variables:
      integer :: lm,l,m
      complex(cp) :: help1,help2

      real(cp) :: rRatio,rDep(0:l_max),rDep2(0:l_max)


      rRatio  =r/r_ICB
      rDep(0) =rRatio
      rDep2(0)=one/r_ICB ! rDep2=rDep/r
      do l=1,l_max
         rDep(l) =rDep(l-1)*rRatio
         rDep2(l)=rDep2(l-1)*rRatio
      end do

      lm=0
      if ( .not. l_axi ) then
         do m=0,l_max,minc
            do l=m,l_max
               lm=lm+1
               dLhw(lm)=rDep(l)*dLh(lm)*w(lm)
               if ( lHor ) then
                  if ( lCondIC ) then
                     help1=rDep2(l)*((l+1)*w(lm)+r*dw(lm))
                     help2=rDep(l)*z(lm)
                     vhG(lm)=help1-cmplx(-aimag(help2),real(help2),kind=cp )
                     vhC(lm)=help1+cmplx(-aimag(help2),real(help2),kind=cp )
                  else
                     vhG(lm)=rDep2(l)*(l+1)*w(lm) ! Only poloidal
                     vhC(lm)=vhG(lm)
                  end if
               end if
            end do
         end do
      else
         do l=0,l_max
            lm=lm+1
            dLhw(lm)=rDep(l)*dLh(lm)*w(lm)
            if ( lHor ) then
               if ( lCondIC ) then
                  help1=rDep2(l)*((l+1)*w(lm)+r*dw(lm))
                  help2=rDep(l)*z(lm)
                  vhG(lm)=help1-cmplx(-aimag(help2),real(help2),kind=cp )
                  vhC(lm)=help1+cmplx(-aimag(help2),real(help2),kind=cp )
               else
                  vhG(lm)=rDep2(l)*(l+1)*w(lm) ! Only poloidal
                  vhC(lm)=vhG(lm)
               end if
            end if
         end do
      end if
      dLhw(1)=zero
      if ( lHor ) then
         vhG(1) =zero
         vhC(1) =zero
      end if

      if ( lDeriv ) then
         lm=0
         if ( .not. l_axi ) then
            do m=0,l_max,minc
               do l=m,l_max
                  lm=lm+1
                  if ( lCondIC ) then
                     dLhz(lm)=rDep(l)*dLh(lm)*z(lm)
                     help1=rDep2(l)*( (l+1)*z(lm)+r*dz(lm) )
                     help2=rDep2(l)*(-two*(l+1)*dw(lm)-r*ddw(lm))
                     cvhG(lm)=help1-cmplx(-aimag(help2),real(help2),kind=cp)
                     cvhC(lm)=help1+cmplx(-aimag(help2),real(help2),kind=cp)
                  else
                     dLhz(lm)=zero
                     cvhG(lm)=zero
                     cvhC(lm)=zero
                  end if
               end do
            end do
         else
            do l=0,l_max
               lm=lm+1
               if ( lCondIC ) then
                  dLhz(lm)=rDep(l)*dLh(lm)*z(lm)
                  help1=rDep(l)*( (l+1)*z(lm)/r+dz(lm) )
                  help2=rDep(l)*(-two*(l+1)/r*dw(lm)-ddw(lm))
                  cvhG(lm)=help1-cmplx(-aimag(help2),real(help2),kind=cp)
                  cvhC(lm)=help1+cmplx(-aimag(help2),real(help2),kind=cp)
               else
                  dLhz(lm)=zero
                  cvhG(lm)=zero
                  cvhC(lm)=zero
               end if
            end do
         end if
         dLhz(1)=zero
         cvhG(1)=zero
         cvhc(1)=zero
      end if

   end subroutine legPrep_IC
!------------------------------------------------------------------------------
end module leg_helper_mod
