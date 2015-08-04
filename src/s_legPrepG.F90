!$Id$
!***********************************************************************
subroutine legPrepG(nR,nBc,lDeriv,lRmsCalc,l_frame, &
     &              lTOnext,lTOnext2,lTOcalc, leg_helper)
   !  +-------------------------------------------------------------------+
   !  |                                                                   |
   !  |  Purpose of this subroutine is to prepare Legendre transforms     |
   !  |  from (r,l,m) space to (r,theta,m) space by calculating           |
   !  |  auxiliary arrays dpdw,dpddw, ....... dLhj which contain          |
   !  |  combinations of harmonic coeffs multiplied with (l,m)-dependend  |
   !  |  factors as well as the radial dependence.                        |
   !  |    nBc  =0 standard inner radial grid point                       |
   !  |    nBc  =1 radial velocity zero, spatial derivs not needed        |
   !  |    nBc  =2 all velocity comp. zero, spatial derivs not needed     |
   !  |   lDeriv=.true. field derivatives required                        |
   !  |                                                                   |
   !  +-------------------------------------------------------------------+

   use truncation, only: lm_max,l_max
   use radial_data, only: n_r_icb, n_r_cmb
   use radial_functions, only: or2
   use torsional_oscillations, only: ddzASL
   use Grenoble, only: lGrenoble, b0, db0, ddb0
   use blocking, only: lm2l, lm2m, lm2
   use horizontal_data, only: dLh
   use logic, only: l_conv, l_mag_kin, l_heat, l_mag, l_movie_oc,    &
       &            l_mag_LF, l_fluxProfs
   use fields, only: s_Rloc,ds_Rloc, z_Rloc,dz_Rloc, p_Rloc,dp_Rloc, &
       &             b_Rloc,db_Rloc,ddb_Rloc, aj_Rloc,dj_Rloc,       &
       &             w_Rloc,dw_Rloc,ddw_Rloc, omega_ic,omega_ma
   use const, only: zero
   use leg_helper_mod, only: leg_helper_t
 
   implicit none
 
   !-- Input of variables
   integer, intent(in) :: nR              ! radial level
   integer, intent(in) :: nBc             ! boundary condition
   logical, intent(in) :: lDeriv          ! get also field derivatives !
   logical, intent(in) :: lRmsCalc        ! Rms force balance ?
   logical, intent(in) :: l_frame         ! Calculate movie frame?
   logical, intent(in) :: lTOnext         ! for TO output
   logical, intent(in) :: lTOnext2
   logical, intent(in) :: lTOcalc
 
   !-- Input of scalar fields in LM-distributed space
   !   These have to be collected and stored in the
   !   R-distributed output variables
 
   !-- Output variables, R-distributed:
   type(leg_helper_t), intent(inout) :: leg_helper
 
   !-- Local variables:
   integer :: lm,l,m
   complex(kind=8) :: dbd
 
 
   if ( nR == n_r_icb ) leg_helper%omegaIC=omega_ic
   if ( nR == n_r_cmb ) leg_helper%omegaMA=omega_ma
   if ( l_conv .or. l_mag_kin ) then
 
      if ( l_heat ) then
         do lm=1,lm_max
            leg_helper%sR(lm) =s_Rloc(lm,nR)   ! used for plotting and Rms
            leg_helper%dsR(lm)=ds_Rloc(lm,nR)  ! used for plotting and Rms
         end do
      end if
      if ( lTOnext .or. lTOnext2 .or. lTOCalc ) then
         do lm=1,lm_max
            l=lm2l(lm)
            m=lm2m(lm)
            if ( l <= l_max .and. m == 0 ) then
               leg_helper%zAS(l+1)  =real(z_Rloc(lm,nR))   ! used in TO
               leg_helper%dzAS(l+1) =real(dz_Rloc(lm,nR))  ! used in TO (anelastic)
               leg_helper%ddzAS(l+1)=ddzASL(l+1,nR)    ! used in TO
            end if
         end do
      end if
      if ( lRmsCalc .or. l_fluxProfs ) then
         leg_helper%preR(1)=zero
         leg_helper%dpR(1)=zero
         do lm=2,lm_max
            leg_helper%preR(lm)=p_Rloc(lm,nR)    ! used for Rms in get_td (anelastic)
            leg_helper%dpR(lm)=dp_Rloc(lm,nR)  ! used for Rms in get_td
         end do
      end if
      if ( l_mag .and. l_frame .and. l_movie_oc .and. nR == n_r_cmb ) then
         leg_helper%bCMB(1)=zero ! used in s_store_movie_frame.f
         do lm=2,lm_max
            leg_helper%bCMB(lm)=b_Rloc(lm,nR)  ! used for movie output of surface field
         end do
      end if
 
      if ( nBc /= 2 ) then ! nBc=2 is flag for fixed boundary
         leg_helper%dLhw(1)=zero
         leg_helper%vhG(1) =zero
         leg_helper%vhC(1) =zero
         do lm=2,lm_max
            leg_helper%dLhw(lm)=dLh(lm)*w_Rloc(lm,nR)
            leg_helper%vhG(lm) =dw_Rloc(lm,nR) - &
                 cmplx(-aimag(z_Rloc(lm,nR)),real(z_Rloc(lm,nR)),kind=kind(0d0))
            leg_helper%vhC(lm) =dw_Rloc(lm,nR) + &
                 cmplx(-aimag(z_Rloc(lm,nR)),real(z_Rloc(lm,nR)),kind=kind(0d0))
         end do
      end if
 
      if ( lDeriv ) then
         leg_helper%dLhdw(1) =zero
         leg_helper%dLhz(1)  =zero
         leg_helper%dvhdrG(1)=zero
         leg_helper%dvhdrC(1)=zero
         do lm=2,lm_max
            leg_helper%dLhz(lm)  =dLh(lm)*z_Rloc(lm,nR)
            leg_helper%dLhdw(lm) =dLh(lm)*dw_Rloc(lm,nR)
            leg_helper%dvhdrG(lm)=ddw_Rloc(lm,nR) - &
                 cmplx(-aimag(dz_Rloc(lm,nR)),real(dz_Rloc(lm,nR)),kind=kind(0d0))
            leg_helper%dvhdrC(lm)=ddw_Rloc(lm,nR) + &
                 cmplx(-aimag(dz_Rloc(lm,nR)),real(dz_Rloc(lm,nR)),kind=kind(0d0))
         end do
      end if
 
   end if
 
   if ( l_mag .or. l_mag_LF ) then
 
      !PRINT*,"aj: ",SUM(ABS(aj(:,nR))),SUM(ABS(dLh))
      !PRINT*,"dj: ",SUM(ABS(dj(:,nR)))
      leg_helper%dLhb(1)=zero
      leg_helper%bhG(1) =zero
      leg_helper%bhC(1) =zero
      do lm=2,lm_max
         leg_helper%dLhb(lm)=dLh(lm)*b_Rloc(lm,nR)
         leg_helper%bhG(lm) =db_Rloc(lm,nR) - &
              cmplx(-aimag(aj_Rloc(lm,nR)),real(aj_Rloc(lm,nR)),kind=kind(0d0))
         leg_helper%bhC(lm) =db_Rloc(lm,nR) + &
              cmplx(-aimag(aj_Rloc(lm,nR)),real(aj_Rloc(lm,nR)),kind=kind(0d0))
      end do
      if ( lGrenoble ) then ! Add dipole imposed by inner core
         lm=lm2(1,0)
         leg_helper%dLhb(lm)=leg_helper%dLhb(lm)+dLh(lm)*b0(nR)
         leg_helper%bhG(lm) =leg_helper%bhG(lm)+db0(nR)
         leg_helper%bhC(lm) =leg_helper%bhC(lm)+db0(nR)
      end if
      if ( lDeriv ) then
         leg_helper%dLhj(1)=zero
         leg_helper%cbhG(1)=zero
         leg_helper%cbhC(1)=zero
         do lm=2,lm_max
            leg_helper%dLhj(lm)=dLh(lm)*aj_Rloc(lm,nR)
            dbd     =or2(nR)*leg_helper%dLhb(lm)-ddb_Rloc(lm,nR)
            leg_helper%cbhG(lm)=dj_Rloc(lm,nR)-cmplx(-aimag(dbd),real(dbd), &
                                                     kind=kind(0d0))
            leg_helper%cbhC(lm)=dj_Rloc(lm,nR)+cmplx(-aimag(dbd),real(dbd), &
                                                     kind=kind(0d0))
         end do
         if ( lGrenoble ) then ! Add dipole imposed by inner core
            lm=lm2(1,0)
            leg_helper%cbhG(lm)=leg_helper%cbhG(lm)+cmplx(0.D0,ddb0(nR),kind=kind(0d0))
            leg_helper%cbhC(lm)=leg_helper%cbhC(lm)-cmplx(0.D0,ddb0(nR),kind=kind(0d0))
         end if
      end if
 
   end if   ! magnetic terms required ?
 
 
end subroutine legPrepG
!------------------------------------------------------------------------
