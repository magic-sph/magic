!$Id$
module getDlm_mod

   use parallel_mod
   use truncation, only: minc, m_max, l_max, n_r_max
   use radial_functions, only: or2, drx, i_costf_init, d_costf_init, &
                               orho1
   use num_param, only: eScale
   use blocking, only: lo_map, st_map
   use horizontal_data, only: dLh
   use const, only: pi
   use LMLoop_data,only: llm, ulm
   use useful, only: cc2real, cc22real
   use integration, only: rInt_R
   
   implicit none
 
   private
 
   public :: getDlm

contains

   subroutine getDlm(w,dw,z,dl,dlR,dm,dlc,dlRc,switch)
      !--------------------------------------------------------------------
      !  calculates energy  = 1/2 Integral(B^2 dV)
      !  integration in theta,phi by summation over harmonic coeffs.
      !  integration in r by Chebycheff integrals
      !
      !  Output:
      !  enbp: Total poloidal        enbt: Total toroidal
      !  apome: Axisym. poloidal     atome: Axisym. toroidal
      !--------------------------------------------------------------------

      !-- Input variables:
      complex(kind=8),  intent(in) :: w(llm:ulm,n_r_max)
      complex(kind=8),  intent(in) :: dw(llm:ulm,n_r_max)
      complex(kind=8),  intent(in) :: z(llm:ulm,n_r_max)
      character(len=1), intent(in) :: switch

      !-- Output variables:
      real(kind=8),intent(out) :: dlR(n_r_max),dlRc(n_r_max)
      real(kind=8),intent(out) :: dl,dlc,dm

      !-- Local variables:
      integer :: nR,lm,l,m,lFirst
      real(kind=8) :: e_p,e_t,e_m,e_l
      real(kind=8) :: fac
      real(kind=8) :: e_lr(n_r_max,l_max),e_lr_c(n_r_max,l_max)
      real(kind=8) :: e_lr_global(n_r_max,l_max),e_lr_c_global(n_r_max,l_max)
      real(kind=8) :: e_mr(n_r_max,0:l_max)
      real(kind=8) :: e_mr_global(n_r_max,0:l_max)
      real(kind=8) :: ER(n_r_max),ELR(n_r_max)
      real(kind=8) :: E,EL,EM
      real(kind=8) :: ERc(n_r_max),ELRc(n_r_max)
      real(kind=8) :: Ec,ELc
      real(kind=8) :: O_rho ! 1/rho (anelastic)

      if ( switch == 'B' ) then
         do nR=1,n_r_max
            e_mr(nR,0) = 0.0D0
            do l=1,l_max
               e_lr(nR,l)=0.D0
               e_lr_c(nR,l)=0.D0
               e_mr(nR,l)=0.D0
            end do
            do lm=max(2,llm),ulm
               l =lo_map%lm2l(lm)
               m =lo_map%lm2m(lm)

               e_p= dLh(st_map%lm2(l,m)) *  ( &
                    dLh(st_map%lm2(l,m))*or2(nR)*cc2real(w(lm,nR),m) &
                    & + cc2real(dw(lm,nR),m) )
               e_t=dLh(st_map%lm2(l,m))*cc2real(z(lm,nR),m)

               e_lr(nR,l)=e_lr(nR,l) + e_p+e_t
               e_lr_c(nR,l)=0.D0
               e_mr(nR,m)=e_mr(nR,m) + e_p+e_t
            end do ! do loop over lms in block
            ! We have now a local sum over the local lm in
            ! e_lr(nR,l), e_mr(nR,m)
         end do    ! radial grid points

         lFirst=2
      else if ( switch == 'V' ) then
         do nR=1,n_r_max
            O_rho =orho1(nR)
            e_mr(nR,0) = 0.0D0
            do l=1,l_max
               e_lr(nR,l)=0.D0
               e_lr_c(nR,l)=0.D0
               e_mr(nR,l)=0.D0
            end do
            do lm=max(2,llm),ulm
               l =lo_map%lm2l(lm)
               m =lo_map%lm2m(lm)

               e_p= O_rho * dLh(st_map%lm2(l,m)) *  ( &
                    dLh(st_map%lm2(l,m))*or2(nR)*cc2real(w(lm,nR),m) &
                    & + cc2real(dw(lm,nR),m) )
               e_t=O_rho*dLh(st_map%lm2(l,m))*cc2real(z(lm,nR),m)
               if ( m /= 0 ) then
                  e_lr_c(nR,l)=e_lr_c(nR,l) + e_p+e_t
               end if
               e_lr(nR,l)=e_lr(nR,l) + e_p+e_t
               e_mr(nR,m)=e_mr(nR,m) + e_p+e_t
               !if (nR == n_r_icb) then
               !   write(*,"(A,3I4,10ES20.12)") "e_lr,e_mr,e_p,e_t = ",nR,l,m,&
               !        &e_lr(nR,l),e_mr(nR,m),&
               !        &e_p,e_t,w(lm,nR),dw(lm,nR),z(lm,nR)
               !end if
            end do ! do loop over lms in block
         end do    ! radial grid points
         lFirst=1
      else
         write(*,*) 'Wrong switch in s_getDlm.f'
         stop
      end if

      ! reduce to rank 0
      call MPI_Reduce(e_lr,e_lr_global,n_r_max*l_max,MPI_DOUBLE_PRECISION,&
           &          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_mr,e_mr_global,n_r_max*(l_max+1),MPI_DOUBLE_PRECISION,&
           &          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_lr_c,e_lr_c_global,n_r_max*l_max,MPI_DOUBLE_PRECISION,&
           &          MPI_SUM,0,MPI_COMM_WORLD,ierr)
         
      if ( rank == 0 ) then
         !-- Radial Integrals:
         fac=0.5D0*eScale
         E  =0.D0
         EL =0.D0
         Ec =0.D0
         ELc=0.D0

         do l=lFirst,l_max
            e_l=0.d0
            e_l=fac*rInt_R(e_lr_global(1,l),n_r_max,n_r_max,drx, &
                 &         i_costf_init,d_costf_init)
            !write(*,"(A,I5,ES20.12)") "getDlm: l,e_l = ",l,e_l
            E =E+e_l
            EL=EL+dble(l)*e_l
            e_l=0.d0
            e_l=fac*rInt_R(e_lr_c_global(1,l),n_r_max,n_r_max,drx, &
                           i_costf_init,d_costf_init)
            Ec =Ec+e_l
            ELc=ELc+dble(l)*e_l
         end do
         if ( EL /= 0d0 ) then
            !write(*,"(A,2ES20.12)") "getDlm: E,EL = ",E,EL
            dl=pi*E/EL
         else
            dl=0d0
         end if
         if ( switch == 'V' ) then
            if ( ELc /= 0d0 ) then
               dlc=pi*Ec/ELc
            else
               dlc=0d0
            end if
         else if ( switch == 'B' ) then
            dlc=0.d0
         end if
         do nR=1,n_r_max
            ER(nR)  =0.D0
            ELR(nR) =0.D0
            ERc(nR) =0.D0
            ELRc(nR)=0.D0
            do l=lFirst,l_max
               e_l=fac*e_lr_global(nR,l)
               ER(nR) =ER(nR)+e_l
               ELR(nR)=ELR(nR)+dble(l)*e_l
               if ( switch == 'V' ) then
                  e_l=fac*e_lr_c_global(nR,l)
                  ERc(nR) =ERc(nR)+e_l
                  ELRc(nR)=ELRc(nR)+dble(l)*e_l
               end if
            end do
            if ( switch == 'V' ) then
               if ( ELR(nR) /= 0d0 ) then
                  dlR(nR)=pi*ER(nR)/ELR(nR)
               else
                  dlR(nR)=0.d0
               end if
               if ( ELRc(nR) /= 0d0 ) then
                  dlRc(nR)=pi*ERc(nR)/ELRc(nR)
               else
                  dlRc(nR)=0.d0
               end if
               !write(*,"(I3,A,2ES20.12)") nR,": dlRc,dlR = ",dlRc(nR),dlR(nR)
            else if ( switch == 'B' ) then
               dlR(nR)=0.D0
               dlRc(nR)=0.D0
            end if
         end do
         E =0.D0
         EM=0.D0
         do m=minc,m_max,minc
            e_m=fac*rInt_R(e_mr_global(1,m),n_r_max,n_r_max,drx, &
                           i_costf_init,d_costf_init)
            E =E +e_m
            EM=EM+dble(m)*e_m
         end do
         if ( EM /= 0d0 ) then
            dm=pi*E/EM
         else
            dm=0d0
         end if
      end if

   end subroutine getDlm
!------------------------------------------------------------------------------
end module getDlm_mod
