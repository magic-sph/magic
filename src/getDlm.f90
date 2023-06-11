module getDlm_mod
   !
   ! This module is used to calculate the lengthscales. It computes both
   ! the integral lengthscale and the peak of the poloidal energy. It also
   ! stores the radial profiles of these lengthscales.
   !

   use parallel_mod
   use precision_mod
   use communications, only: reduce_radial
   use truncation, only: minc, m_max, l_max, n_r_max
   use radial_functions, only: or2, r, rscheme_oc, orho1
   use num_param, only: eScale
   use blocking, only: lo_map, llm, ulm
   use constants, only: pi, half
   use useful, only: cc2real, cc22real
   use integration, only: rInt_R
   use useful, only: abortRun

   implicit none

   private

   public :: getDlm

contains

   subroutine getDlm(w,dw,z,dl,dlR,dm,dlc,dlPolPeak,dlRc,dlPolPeakR,switch)
      !
      ! This routine is used to compute integral lengthscale using spectra
      !

      !-- Input variables:
      complex(cp),      intent(in) :: w(llm:ulm,n_r_max)
      complex(cp),      intent(in) :: dw(llm:ulm,n_r_max)
      complex(cp),      intent(in) :: z(llm:ulm,n_r_max)
      character(len=1), intent(in) :: switch

      !-- Output variables:
      real(cp), intent(out) :: dlR(n_r_max),dlRc(n_r_max),dlPolPeakR(n_r_max)
      real(cp), intent(out) :: dl,dlc,dm,dlPolPeak

      !-- Local variables:
      integer :: nR,lm,l,m,lFirst,lpeak(1)
      real(cp) :: e_p,e_t,e_m,e_l,fac
      real(cp) :: e_lr(n_r_max,0:l_max),e_lr_c(n_r_max,0:l_max),e_pol_lr(n_r_max,0:l_max)
      real(cp) :: e_lr_global(n_r_max,0:l_max),e_lr_c_global(n_r_max,0:l_max)
      real(cp) :: e_pol_lr_global(n_r_max,0:l_max), e_pol_l(l_max)
      real(cp) :: e_mr(n_r_max,0:l_max), e_mr_global(n_r_max,0:l_max)
      real(cp) :: ER(n_r_max),ELR(n_r_max),ERc(n_r_max),ELRc(n_r_max)
      real(cp) :: E,EL,EM,Ec,ELc
      real(cp) :: O_rho ! 1/rho (anelastic)

      e_lr(:,:)    =0.0_cp
      e_lr_c(:,:)  =0.0_cp
      e_pol_lr(:,:)=0.0_cp
      e_mr(:,:)    =0.0_cp
      if ( switch == 'B' ) then
         do nR=1,n_r_max
            do lm=llm,ulm
               l =lo_map%lm2l(lm)
               m =lo_map%lm2m(lm)

               e_p= real(l*(l+1),cp)*(real(l*(l+1),cp)*or2(nR)*cc2real(w(lm,nR),m) &
               &                      +cc2real(dw(lm,nR),m))
               e_t=real(l*(l+1),cp)*cc2real(z(lm,nR),m)

               e_lr(nR,l)    =e_lr(nR,l) + e_p+e_t
               e_lr_c(nR,l)  =0.0_cp
               e_pol_lr(nR,l)=e_pol_lr(nR,l)+e_p
               e_mr(nR,m)    =e_mr(nR,m) + e_p+e_t
            end do ! do loop over lms in block
            ! We have now a local sum over the local lm in
            ! e_lr(nR,l), e_mr(nR,m)
         end do    ! radial grid points

         lFirst=2
      else if ( switch == 'V' ) then
         do nR=1,n_r_max
            O_rho =orho1(nR)
            do lm=llm,ulm
               l =lo_map%lm2l(lm)
               m =lo_map%lm2m(lm)

               e_p=O_rho*real(l*(l+1),cp)*(real(l*(l+1),cp)*or2(nR)*cc2real(w(lm,nR),m) &
               &                           +cc2real(dw(lm,nR),m) )
               e_t=O_rho*real(l*(l+1),cp)*cc2real(z(lm,nR),m)
               if ( m /= 0 ) then
                  e_lr_c(nR,l)=e_lr_c(nR,l) + e_p+e_t
               end if
               e_lr(nR,l)    =e_lr(nR,l) + e_p+e_t
               e_pol_lr(nR,l)=e_pol_lr(nR,l)+e_p
               e_mr(nR,m)    =e_mr(nR,m) + e_p+e_t
            end do ! do loop over lms in block
         end do    ! radial grid points
         lFirst=1
      else
         call abortRun('Wrong switch in getDlm')
      end if

      ! reduce to rank 0
      call reduce_radial(e_lr, e_lr_global, 0)
      call reduce_radial(e_mr, e_mr_global, 0)
      call reduce_radial(e_lr_c, e_lr_c_global, 0)
      call reduce_radial(e_pol_lr, e_pol_lr_global, 0)

      if ( rank == 0 ) then
         !-- Radial Integrals:
         fac=half*eScale
         E  =0.0_cp
         EL =0.0_cp
         Ec =0.0_cp
         ELc=0.0_cp
         e_pol_l(:)=0.0_cp

         do l=lFirst,l_max
            e_l=0.0_cp
            e_l=fac*rInt_R(e_lr_global(:,l),r,rscheme_oc)
            E =E+e_l
            EL=EL+real(l,cp)*e_l
            e_l=0.0_cp
            e_l=fac*rInt_R(e_lr_c_global(:,l),r,rscheme_oc)
            Ec =Ec+e_l
            ELc=ELc+real(l,cp)*e_l
            e_pol_l(l)=fac*rInt_R(e_pol_lr_global(:,l),r,rscheme_oc)
         end do

         !-- Peak of poloidal energy
         lpeak(1) = maxloc(e_pol_l, dim=1)
         dlPolPeak = pi/real(lpeak(1),cp)
         if ( EL /= 0.0_cp ) then
            dl=pi*E/EL
         else
            dl=0.0_cp
         end if
         if ( switch == 'V' ) then
            if ( ELc /= 0.0_cp ) then
               dlc=pi*Ec/ELc
            else
               dlc=0.0_cp
            end if
         else if ( switch == 'B' ) then
            dlc=0.0_cp
         end if
         do nR=1,n_r_max
            if ( switch == 'V') then
               e_pol_l(:) = e_pol_lr_global(nR,1:l_max)
               lpeak(1) = maxloc(e_pol_l, dim=1)
               if ( maxval(e_pol_l) > 10.0_cp*epsilon(0.0_cp) ) then
                  dlPolPeakR(nR)=pi/real(lpeak(1),cp)
               else
                  dlPolPeakR(nR)=0.0_cp
               end if
            end if
            ER(nR)  =0.0_cp
            ELR(nR) =0.0_cp
            ERc(nR) =0.0_cp
            ELRc(nR)=0.0_cp
            do l=lFirst,l_max
               e_l=fac*e_lr_global(nR,l)
               ER(nR) =ER(nR)+e_l
               ELR(nR)=ELR(nR)+real(l,cp)*e_l
               if ( switch == 'V' ) then
                  e_l=fac*e_lr_c_global(nR,l)
                  ERc(nR) =ERc(nR)+e_l
                  ELRc(nR)=ELRc(nR)+real(l,cp)*e_l
               end if
            end do
            if ( switch == 'V' ) then
               if ( ELR(nR) > 10.0_cp*epsilon(0.0_cp) ) then
                  dlR(nR)=pi*ER(nR)/ELR(nR)
               else
                  dlR(nR)=0.0_cp
               end if
               if ( ELRc(nR) > 10.0_cp*epsilon(0.0_cp) ) then
                  dlRc(nR)=pi*ERc(nR)/ELRc(nR)
               else
                  dlRc(nR)=0.0_cp
               end if
            else if ( switch == 'B' ) then
               dlR(nR)=0.0_cp
               dlRc(nR)=0.0_cp
            end if
         end do
         E =0.0_cp
         EM=0.0_cp
         do m=minc,m_max,minc
            e_m=fac*rInt_R(e_mr_global(:,m),r,rscheme_oc)
            E =E +e_m
            EM=EM+real(m,cp)*e_m
         end do
         if ( EM /= 0.0_cp ) then
            dm=pi*E/EM
         else
            dm=0.0_cp
         end if
      end if

   end subroutine getDlm
!------------------------------------------------------------------------------
end module getDlm_mod
