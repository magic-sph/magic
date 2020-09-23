module nonlinear_lm_mod
   !
   ! This module is used to finish the assembly of the nonlinear terms in
   ! :math:`(\ell,m)` space. Derivatives along :math:`\theta` and :math:`\phi`
   ! are handled using recurrence relations.
   !

   use, intrinsic :: iso_c_binding
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: l_max, n_lm_loc, n_lmP_loc
   use logic, only : l_anel, l_conv_nl, l_corr, l_heat, l_anelastic_liquid, &
       &             l_mag_nl, l_mag_kin, l_mag_LF, l_conv, l_mag, l_RMS,   &
       &             l_chemical_conv, l_single_matrix, l_double_curl,       &
       &             l_adv_curl
   use radial_functions, only: r, or2, or1, beta, rho0, rgrav, epscProf, &
       &                       or4, temp0, alpha0, ogrun, orho1, dLvisc, &
       &                       visc, dbeta
   use physical_parameters, only: CorFac, ra, epsc, ViscHeatFac, &
       &                          OhmLossFac, n_r_LCR, epscXi,   &
       &                          BuoFac, ChemFac, ThExpNb
   use horizontal_data, only: dLh_loc, dTheta1S_loc, dTheta1A_loc, dPhi_loc,    &
       &                      dTheta2A_loc, dTheta3A_loc, dTheta4A_loc,         &
       &                      dTheta2S_loc, hdif_B, dTheta3S_loc, dTheta4S_loc, &
       &                      hdif_V
   use RMS, only: Adv2hInt, Pre2hInt, Cor2hInt, LF2hInt, Iner2hInt,    &
       &          Geo2hInt, Mag2hInt, ArcMag2hInt, CLF2hInt, PLF2hInt, &
       &          CIA2hInt, Arc2hInt, Buo_temp2hInt, Buo_xi2hint
   use constants, only: zero, two, four, third, three
   use fields, only: w_Rdist, dw_Rdist, ddw_Rdist, z_Rdist, dz_Rdist, s_Rdist, &
       &             p_Rdist, dp_Rdist, xi_Rdist
   use RMS_helpers, only: hIntRms
   use LMmapping, only: map_dist_st


   implicit none

   type :: nonlinear_lm_t
      !----- Nonlinear terms in lm-space:
      complex(cp), allocatable :: AdvrLM(:), AdvtLM(:), AdvpLM(:)
      complex(cp), allocatable :: LFrLM(:),  LFtLM(:),  LFpLM(:)
      complex(cp), allocatable :: VxBrLM(:), VxBtLM(:), VxBpLM(:)
      complex(cp), allocatable :: VSrLM(:),  VStLM(:),  VSpLM(:)
      complex(cp), allocatable :: VXirLM(:),  VXitLM(:),  VXipLM(:)
      complex(cp), allocatable :: ViscHeatLM(:), OhmLossLM(:)
      !----- RMS calculations
      complex(cp), allocatable :: Advt2LM(:), Advp2LM(:), PFt2LM(:), PFp2LM(:)
      complex(cp), allocatable :: LFt2LM(:), LFp2LM(:), CFt2LM(:), CFp2LM(:)
      complex(cp), allocatable :: dtVtLM(:), dtVpLM(:), dtVrLM(:)
      complex(cp), allocatable :: dpkindrLM(:)

   contains

      procedure :: initialize
      procedure :: finalize
      procedure :: set_zero
      procedure :: get_td

   end type nonlinear_lm_t

contains

   subroutine initialize(this,n_loc)
      !
      ! Memory allocation of ``get_td`` arrays
      !

      class(nonlinear_lm_t) :: this
      integer, intent(in) :: n_loc

      allocate( this%AdvrLM(n_loc), this%AdvtLM(n_loc) )
      allocate( this%AdvpLM(n_loc), this%LFrLM(n_loc) )
      allocate( this%LFtLM(n_loc), this%LFpLM(n_loc) )
      allocate( this%VxBrLM(n_loc), this%VxBtLM(n_loc) )
      allocate( this%VxBpLM(n_loc))
      bytes_allocated = bytes_allocated + 9*n_loc*SIZEOF_DEF_COMPLEX

      if ( l_anel) then
         allocate( this%ViscHeatLM(n_loc), this%OhmLossLM(n_loc) )
         bytes_allocated = bytes_allocated+2*n_loc*SIZEOF_DEF_COMPLEX
      end if

      if ( l_heat ) then
         allocate(this%VSrLM(n_loc),this%VStLM(n_loc),this%VSpLM(n_loc))
         bytes_allocated = bytes_allocated + 3*n_loc*SIZEOF_DEF_COMPLEX
      end if

      if ( l_chemical_conv ) then
         allocate(this%VXirLM(n_loc),this%VXitLM(n_loc),this%VXipLM(n_loc))
         bytes_allocated = bytes_allocated + 3*n_loc*SIZEOF_DEF_COMPLEX
      end if

      !-- RMS calculations
      if ( l_RMS ) then
         allocate( this%dtVrLM(n_loc), this%dtVtLM(n_loc) )
         allocate( this%dtVpLM(n_loc) )
         allocate( this%Advt2LM(n_loc), this%Advp2LM(n_loc) )
         allocate( this%LFt2LM(n_loc), this%LFp2LM(n_loc) )
         allocate( this%CFt2LM(n_loc), this%CFp2LM(n_loc) )
         allocate( this%PFt2LM(n_loc), this%PFp2LM(n_loc) )
         bytes_allocated = bytes_allocated + 11*n_loc*SIZEOF_DEF_COMPLEX
         if ( l_adv_curl ) then
            allocate( this%dpkindrLM(n_loc) )
            bytes_allocated = bytes_allocated + n_loc*SIZEOF_DEF_COMPLEX
         end if
      end if

   end subroutine initialize
!----------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Memory deallocation
      !

      class(nonlinear_lm_t) :: this

      deallocate( this%AdvrLM, this%AdvtLM, this%AdvpLM )
      deallocate( this%LFrLM, this%LFtLM, this%LFpLM )
      deallocate( this%VxBrLM, this%VxBtLM, this%VxBpLM )

      if ( l_anel ) deallocate( this%ViscHeatLM, this%OhmLossLM )

      if ( l_chemical_conv ) deallocate( this%VXirLM, this%VXitLM, this%VXipLM )

      if ( l_heat ) deallocate( this%VSrLM, this%VStLM, this%VSpLM )

      !-- RMS calculations
      if ( l_RMS ) then
         deallocate( this%Advt2LM, this%Advp2LM, this%LFt2LM, this%LFp2LM )
         deallocate( this%CFt2LM, this%CFp2LM, this%PFt2LM, this%PFp2LM )
         deallocate( this%dtVrLM, this%dtVtLM, this%dtVpLM )
         if ( l_adv_curl ) deallocate ( this%dpkindrLM )
      end if

   end subroutine finalize
!----------------------------------------------------------------------------
   subroutine set_zero(this, n_loc)
      !
      ! Set all the arrays to zero
      !

      class(nonlinear_lm_t) :: this
      integer, intent(in) :: n_loc

      !-- Local variable
      integer :: lm

      !$omp parallel do private(lm)
      do lm=1,n_loc
         this%AdvrLM(lm)=zero
         this%AdvtLM(lm)=zero
         this%AdvpLM(lm)=zero
         this%LFrLM(lm) =zero
         this%LFtLM(lm) =zero
         this%LFpLM(lm) =zero
         this%VxBrLM(lm)=zero
         this%VxBtLM(lm)=zero
         this%VxBpLM(lm)=zero
         if ( l_heat ) then
            this%VSrLM(lm) =zero
            this%VStLM(lm) =zero
            this%VSpLM(lm) =zero
         end if
         if ( l_anel ) then
            this%ViscHeatLM(lm)=zero
            this%OhmLossLM(lm) =zero
         end if

         if ( l_chemical_conv ) then
            this%VXirLM(lm)=zero
            this%VXitLM(lm)=zero
            this%VXipLM(lm)=zero
         end if

         if ( l_RMS ) then
            this%Advt2LM(lm)=zero
            this%Advp2LM(lm)=zero
            this%LFp2LM(lm) =zero
            this%LFt2LM(lm) =zero
            this%CFt2LM(lm) =zero
            this%CFp2LM(lm) =zero
            this%PFt2LM(lm) =zero
            this%PFp2LM(lm) =zero
            this%dtVtLM(lm) =zero
            this%dtVpLM(lm) =zero
            this%dtVrLM(lm) =zero
            if ( l_adv_curl ) this%dpkindrLM(lm)=zero
         end if
      end do
      !$omp end parallel do

   end subroutine set_zero
!----------------------------------------------------------------------------
   subroutine get_td(this,nR,nBc,lRmsCalc,lPressNext,dVSrLM,dVXirLM,    &
              &      dVxVhLM,dVxBhLM,dwdt,dzdt,dpdt,dsdt,dxidt,dbdt,djdt)
      !
      !  Purpose of this to calculate time derivatives
      !  ``dwdt``,``dzdt``,``dpdt``,``dsdt``,``dxidt``,``dbdt``,``djdt``
      !  and auxiliary arrays ``dVSrLM``, ``dVXirLM`` and ``dVxBhLM``, ``dVxVhLM``
      !  from non-linear terms in spectral form
      !

      !-- Input of variables:
      class(nonlinear_lm_t) :: this

      integer, intent(in) :: nR  ! Radial level
      integer, intent(in) :: nBc ! signifies boundary conditions
      logical, intent(in) :: lRmsCalc
      logical, intent(in) :: lPressNext

      !-- Output of variables:
      complex(cp), intent(out) :: dwdt(:),dzdt(:)
      complex(cp), intent(out) :: dpdt(:),dsdt(:)
      complex(cp), intent(out) :: dxidt(:)
      complex(cp), intent(out) :: dbdt(:),djdt(:)
      complex(cp), intent(out) :: dVxBhLM(:)
      complex(cp), intent(out) :: dVxVhLM(:)
      complex(cp), intent(out) :: dVSrLM(:)
      complex(cp), intent(out) :: dVXirLM(:)

      !-- Local variables:
      integer :: l,m,lm,lmS,lmA,lmP,lmPS,lmPA
      complex(cp) :: CorPol(n_lm_loc), dpdr(n_lm_loc)
      complex(cp) :: AdvPol(n_lm_loc),AdvTor(n_lm_loc)
      complex(cp) :: LFPol(n_lm_loc),LFTor(n_lm_loc)
      complex(cp) :: Geo(n_lm_loc),CLF(n_lm_loc),PLF(n_lm_loc)
      complex(cp) :: ArcMag(n_lm_loc),Mag(n_lm_loc),CIA(n_lm_loc),Arc(n_lm_loc)
      complex(cp) :: Buo_temp(n_lm_loc), Buo_xi(n_lm_loc)
      complex(cp) :: AdvPol_loc,CorPol_loc,AdvTor_loc,CorTor_loc
      complex(cp) :: dsdt_loc
      complex(cp) :: dxidt_loc
      integer :: lm_maybe_skip_first
      integer, pointer :: lm2l(:), lm2m(:), lm2lmP(:), lm2(:,:)
      integer, pointer :: lmP2lmPS(:), lmP2lmPA(:), lm2lmA(:), lm2lmS(:)

      lm2l(1:n_lm_loc) => map_dist_st%lm2l
      lm2m(1:n_lm_loc) => map_dist_st%lm2m
      lmP2lmPS(1:n_lmP_loc) => map_dist_st%lmP2lmPS
      lmP2lmPA(1:n_lmP_loc) => map_dist_st%lmP2lmPA
      lm2lmP(1:n_lm_loc) => map_dist_st%lm2lmP
      lm2lmS(1:n_lm_loc) => map_dist_st%lm2lmS
      lm2lmA(1:n_lm_loc) => map_dist_st%lm2lmA
      lm2(0:,0:) => map_dist_st%lm2

      lm_maybe_skip_first = 1
      if (map_dist_st%lm2(0,0) > 0) lm_maybe_skip_first = 2

      if (nBc == 0 .or. lRmsCalc ) then

         if ( l_conv ) then  ! Convection

            lm =lm2(0,0)   ! This is l=0,m=0
            if ( lm > 0 ) then
               lmA=lm2lmA(lm)
               lmP=lm2lmP(lm)
               lmPA=lmP2lmPA(lmP)
               if ( l_conv_nl ) then
                  AdvPol_loc=      or2(nR)*this%AdvrLM(lm)
                  AdvTor_loc=-dTheta1A_loc(lm)*this%AdvpLM(lmPA)
               else
                  AdvPol_loc=zero
                  AdvTor_loc=zero
               end if
               if ( l_corr ) then
                  CorPol_loc=two*CorFac*or1(nR) * dTheta2A_loc(lm)*z_Rdist(lmA,nR)
                  CorTor_loc= two*CorFac*or2(nR) * (                      &
                  &                dTheta3A_loc(lm)*dw_Rdist(lmA,nR) +    &
                  &        or1(nR)*dTheta4A_loc(lm)* w_Rdist(lmA,nR) )
               else
                  CorPol_loc=zero
                  CorTor_loc=zero
               end if

               if ( l_single_matrix ) then
                  dwdt(lm)=AdvPol_loc!+CorPol_loc
               else
                  dwdt(lm)=AdvPol_loc+CorPol_loc
               end if

               dzdt(lm)=AdvTor_loc+CorTor_loc

               if ( lRmsCalc ) then
                  if (l_heat) then
                     Buo_temp(lm) =BuoFac*rgrav(nR)*rho0(nR)*s_Rdist(lm,nR)
                  else
                     Buo_temp(lm) =0.0_cp
                  end if
                  if (l_chemical_conv) then
                     Buo_xi(lm) =ChemFac*rgrav(nR)*rho0(nR)*xi_Rdist(lm,nR)
                  else
                     Buo_xi(lm)=0.0_cp
                  end if
                  if ( l_mag_LF .and. nR>n_r_LCR ) then
                     LFPol(lm) =      or2(nR)*this%LFrLM(lm)
                     LFTor(lm) =-dTheta1A_loc(lm)*this%LFpLM(lmPA)
                     AdvPol(lm)=AdvPol_loc-LFPol(lm)
                     AdvTor(lm)=AdvTor_loc-LFTor(lm)
                  else
                     AdvPol(lm)=AdvPol_loc
                     AdvTor(lm)=AdvTor_loc
                  end if
                  CorPol(lm)=CorPol_loc

                  if ( l_double_curl ) then
                     !-- Recalculate the pressure gradient based on the poloidal
                     !-- equation equilibrium
                     dpdr(lm)=Buo_temp(lm)+Buo_xi(lm)+beta(nR)*p_Rdist(lm,nR)+ &
                     &        AdvPol_loc+CorPol_loc
                  else
                     dpdr(lm)=dp_Rdist(lm,nR)
                  end if

               end if ! lRmsCalc
            end if

            !
            !$omp parallel do default(shared) private(lm,l,m,lmS,lmA,lmP) &
            !$omp private(lmPS,lmPA,AdvPol_loc,CorPol_loc,AdvTor_loc,CorTor_loc)
            do lm=lm_maybe_skip_first,n_lm_loc
               l   =lm2l(lm)
               m   =lm2m(lm)
               lmS =lm2lmS(lm)
               lmA =lm2lmA(lm)
               lmP =lm2lmP(lm)
               lmPS=lmP2lmPS(lmP)
               lmPA=lmP2lmPA(lmP)

               if ( l_double_curl ) then ! Pressure is not needed

                  if ( l_corr ) then
                     if ( l < l_max .and. l > m ) then
                        CorPol_loc =two*CorFac*or2(nR)*orho1(nR)*(               &
                        &                dPhi_loc(lm)*(                          &
                        &        -ddw_Rdist(lm,nR)+beta(nR)*dw_Rdist(lm,nR)    + &
                        &             ( beta(nR)*or1(nR)+or2(nR))*               &
                        &                        dLh_loc(lm)*w_Rdist(lm,nR) )  + &
                        &         dTheta3A_loc(lm)*( dz_Rdist(lmA,nR)-           &
                        &                           beta(nR)*z_Rdist(lmA,nR) ) + &
                        &         dTheta3S_loc(lm)*( dz_Rdist(lmS,nR)-           &
                        &                           beta(nR)*z_Rdist(lmS,nR) ) + &
                        &          or1(nR)* (                                    &
                        &             dTheta4A_loc(lm)* z_Rdist(lmA,nR)          &
                        &            -dTheta4S_loc(lm)* z_Rdist(lmS,nR) ) )
                     else if ( l == l_max ) then
                        CorPol_loc =0.0_cp
                     else if ( l == m ) then
                        CorPol_loc =two*CorFac*or2(nR)*orho1(nR)*(               &
                        &                dPhi_loc(lm)*(                          &
                        &       -ddw_Rdist(lm,nR)+beta(nR)*dw_Rdist(lm,nR)     + &
                        &             ( beta(nR)*or1(nR)+or2(nR))*               &
                        &                       dLh_loc(lm)*w_Rdist(lm,nR) )   + &
                        &              dTheta3A_loc(lm)*( dz_Rdist(lmA,nR)-      &
                        &                         beta(nR)*z_Rdist(lmA,nR) )   + &
                        &          or1(nR)* (                                    &
                        &              dTheta4A_loc(lm)*   z_Rdist(lmA,nR) ) )
                     end if
                  else
                     CorPol_loc=zero
                  end if

                  if ( l_conv_nl ) then

                     if ( l > m ) then
                        dVxVhLM(lm)=      orho1(nR)*r(nR)*r(nR)* (     &
                        &        dTheta1S_loc(lm)*this%AdvtLM(lmPS) -  &
                        &        dTheta1A_loc(lm)*this%AdvtLM(lmPA) +  &
                        &             dPhi_loc(lm)*this%AdvpLM(lmP)  )
                     else if ( l == m ) then
                        dVxVhLM(lm)=      orho1(nR)*r(nR)*r(nR)* (     &
                        &      - dTheta1A_loc(lm)*this%AdvtLM(lmPA) +  &
                        &        dPhi_loc(lm)*this%AdvpLM(lmP)  )
                     end if

                     AdvPol_loc=dLh_loc(lm)*or4(nR)*orho1(nR)*this%AdvrLM(lmP)

                  else

                     AdvPol_loc =zero
                     dVxVhLM(lm)=zero

                  endif

               else ! We don't use the double curl

                  if ( l_corr .and. nBc /= 2 ) then
                     if ( l < l_max .and. l > m ) then
                        CorPol_loc =two*CorFac*or1(nR) * (      &
                        &       dPhi_loc(lm)*dw_Rdist(lm,nR) +  & ! phi-deriv of dw/dr
                        &   dTheta2A_loc(lm)*z_Rdist(lmA,nR) -  & ! sin(theta) dtheta z
                        &    dTheta2S_loc(lm)*z_Rdist(lmS,nR) )
                     else if ( l == l_max ) then
                        CorPol_loc=0.0_cp
                     else if ( l == m ) then
                        CorPol_loc = two*CorFac*or1(nR) * (      &
                        &        dPhi_loc(lm)*dw_Rdist(lm,nR)  + &
                        &    dTheta2A_loc(lm)*z_Rdist(lmA,nR) )
                     end if
                  else
                     CorPol_loc=zero
                  end if

                  if ( l_conv_nl ) then
                     AdvPol_loc=or2(nR)*this%AdvrLM(lmP)
                  else
                     AdvPol_loc=zero
                  endif

               end if ! Double curl or not for the poloidal equation

               dwdt(lm)=AdvPol_loc+CorPol_loc

               if ( lRmsCalc ) then ! RMS force balance

                  if ( l_anelastic_liquid ) then
                     if (l_heat) then
                        Buo_temp(lm) =BuoFac*alpha0(nR)*rgrav(nR)*(      &
                        &      rho0(nR)*s_Rdist(lm,nR)-ViscHeatFac*      &
                        &     (ThExpNb*alpha0(nR)*temp0(nR)+ogrun(nR))*  &
                        &     p_Rdist(lm,nR) )
                     else
                        Buo_temp(lm) =0.0_cp
                     end if
                     if (l_chemical_conv) then
                        Buo_xi(lm) =ChemFac*alpha0(nR)*rgrav(nR)*(       &
                        &    rho0(nR)*xi_Rdist(lm,nR)-ViscHeatFac*       &
                        &     (ThExpNb*alpha0(nR)*temp0(nR)+ogrun(nR))*  &
                        &     p_Rdist(lm,nR) )
                     else
                        Buo_xi(lm) =0.0_cp
                     end if
                  else
                     if (l_heat) then
                        Buo_temp(lm) =BuoFac*rho0(nR)*rgrav(nR)*s_Rdist(lm,nR)
                     else
                        Buo_temp(lm) =0.0_cp
                     end if
                     if (l_chemical_conv) then
                        Buo_xi(lm) =ChemFac*rho0(nR)*rgrav(nR)*xi_Rdist(lm,nR)
                     else
                        Buo_xi(lm)=0.0_cp
                     end if
                  end if

                  if ( l_double_curl ) then
                     ! In that case we have to recompute the Coriolis force
                     ! since we also want the pressure gradient
                     if ( l_corr .and. nBc /= 2 ) then
                        if ( l < l_max .and. l > m ) then
                           CorPol_loc = two*CorFac*or1(nR) * (      &
                           &        dPhi_loc(lm)*dw_Rdist(lm,nR) +  & ! phi-deriv of dw/dr
                           &    dTheta2A_loc(lm)*z_Rdist(lmA,nR) -  & ! sin(theta) dtheta z
                           &    dTheta2S_loc(lm)*z_Rdist(lmS,nR) )
                        else if ( l == l_max ) then
                           CorPol_loc= 0.0_cp
                        else if ( l == m ) then
                           CorPol_loc = two*CorFac*or1(nR) * (       &
                           &         dPhi_loc(lm)*dw_Rdist(lm,nR)  + &
                           &     dTheta2A_loc(lm)*z_Rdist(lmA,nR) )
                        end if
                     else
                        CorPol_loc=zero
                     end if

                     ! We also need to recompute AdvPol_loc here
                     if ( l_conv_nl ) then
                        AdvPol_loc=or2(nR)*this%AdvrLM(lmP)
                     else
                        AdvPol_loc=zero
                     endif

                     !-- Recalculate the pressure gradient based on the poloidal
                     !-- equation equilibrium
                     dpdr(lm)=Buo_temp(lm)+Buo_xi(lm)-this%dtVrLM(lmP)+          &
                     &      dLh_loc(lm)*or2(nR)*hdif_V(l)*visc(nR)* (            &
                     &                                        ddw_Rdist(lm,nR)+  &
                     &         (two*dLvisc(nR)-third*beta(nR))*dw_Rdist(lm,nR)-  &
                     &      ( dLh_loc(lm)*or2(nR)+four*third*( dbeta(nR)+        &
                     &          dLvisc(nR)*beta(nR)+(three*dLvisc(nR)+beta(nR))* &
                     &         or1(nR)) )*                      w_Rdist(lm,nR))+ &
                     &        beta(nR)*p_Rdist(lm,nR)+AdvPol_loc+CorPol_loc
                  else
                     dpdr(lm)=dp_Rdist(lm,nR)
                  end if

                  if ( l_mag_LF .and. nR>n_r_LCR ) then
                     LFPol(lm) =or2(nR)*this%LFrLM(lmP)
                     AdvPol(lm)=AdvPol_loc-LFPol(lm)
                  else
                     AdvPol(lm)=AdvPol_loc
                  end if
                  CorPol(lm)=CorPol_loc

               end if ! lRmsCalc

               if ( l_corr ) then
                  if ( l < l_max .and. l > m ) then
                     CorTor_loc=           two*CorFac*or2(nR) * (      &
                     &                 dPhi_loc(lm)*z_Rdist(lm,nR)   + &
                     &            dTheta3A_loc(lm)*dw_Rdist(lmA,nR)  + &
                     &    or1(nR)*dTheta4A_loc(lm)* w_Rdist(lmA,nR)  + &
                     &            dTheta3S_loc(lm)*dw_Rdist(lmS,nR)  - &
                     &    or1(nR)*dTheta4S_loc(lm)* w_Rdist(lmS,nR)  )
                  else if ( l == l_max ) then
                     CorTor_loc=0.0_cp
                  else if ( l == m ) then
                     CorTor_loc=  two*CorFac*or2(nR) * (       &
                     &         dPhi_loc(lm)*z_Rdist(lm,nR)   + &
                     &    dTheta3A_loc(lm)*dw_Rdist(lmA,nR)  + &
                     &    or1(nR)*dTheta4A_loc(lm)* w_Rdist(lmA,nR)  )
                  end if
               else
                  CorTor_loc=zero
               end if

               if ( l_conv_nl ) then
                  if ( l > m ) then
                     AdvTor_loc=   -dPhi_loc(lm)*this%AdvtLM(lmP)  + &
                     &          dTheta1S_loc(lm)*this%AdvpLM(lmPS) - &
                     &          dTheta1A_loc(lm)*this%AdvpLM(lmPA)
                  else if ( l == m ) then
                     AdvTor_loc=   -dPhi_loc(lm)*this%AdvtLM(lmP)  - &
                     &          dTheta1A_loc(lm)*this%AdvpLM(lmPA)
                  end if
               else
                  AdvTor_loc=zero
               end if

               dzdt(lm)=CorTor_loc+AdvTor_loc

               if ( lRmsCalc ) then
                  if ( l_mag_LF .and. nR>n_r_LCR ) then
                     !------ When RMS values are required, the Lorentz force is treated
                     !       separately:

                     if ( l > m ) then
                        !------- LFTor= 1/(E*Pm) * curl( curl(B) x B )_r
                        LFTor(lm) =   -dPhi_loc(lm)*this%LFtLM(lmP)  + &
                        &          dTheta1S_loc(lm)*this%LFpLM(lmPS) - &
                        &          dTheta1A_loc(lm)*this%LFpLM(lmPA)
                     else if ( l == m ) then
                        LFTor(lm) =   -dPhi_loc(lm)*this%LFtLM(lmP)  - &
                        &          dTheta1A_loc(lm)*this%LFpLM(lmPA)
                     end if
                     AdvTor(lm)=AdvTor_loc-LFTor(lm)
                  else
                     AdvTor(lm)=AdvTor_loc
                  end if
               end if

            end do
            !$omp end parallel do
            !

            if ( lRmsCalc ) then

               if ( l_conv_nl ) then
                  if ( l_adv_curl ) then
                     do lm=1,n_lm_loc
                        lmP =lm2lmP(lm)
                        AdvPol(lm)=AdvPol(lm)-this%dpkindrLM(lmP)
                     end do
                  end if
                  call hIntRms(AdvPol,nR,1,n_lm_loc,0,Adv2hInt(:,nR),map_dist_st, &
                       &       .false.)
                  call hIntRms(this%Advt2LM,nR,1,n_lmP_loc,1,Adv2hInt(:,nR), &
                       &       map_dist_st, .true.)
                  call hIntRms(this%Advp2LM,nR,1,n_lMP_loc,1,Adv2hInt(:,nR), &
                       &       map_dist_st, .true.)
                  do lm=1,n_lm_loc
                     lmP =lm2lmP(lm)
                     !-- Use Geo as work array
                     Geo(lm)=AdvPol(lm)-this%dtVrLM(lmP)
                  end do
                  call hIntRms(Geo,nR,1,n_lm_loc,0,Iner2hInt(:,nR),map_dist_st, &
                       &       .false.)
                  call hIntRms(this%Advt2LM-this%dtVtLM,nR,1,n_lmP_loc,1, &
                       &       Iner2hInt(:,nR),map_dist_st,.true.)
                  call hIntRms(this%Advp2LM-this%dtVpLM,nR,1,n_lmP_loc,1, &
                       &       Iner2hInt(:,nR),map_dist_st,.true.)
               end if

               if ( l_anelastic_liquid ) then
                  call hIntRms(dpdr,nR,1,n_lm_loc,0,Pre2hInt(:,nR),map_dist_st, &
                       &       .false.)
               else
                  ! rho* grad(p/rho) = grad(p) - beta*p
                  !-- Geo is used to store the pressure Gradient
                  if ( l_adv_curl ) then
                     do lm=1,n_lm_loc
                        lmP =lm2lmP(lm)
                        !-- Use Geo as work array
                        Geo(lm)=dpdr(lm)-beta(nR)*p_Rdist(lm,nR)- &
                        &       this%dpkindrLM(lmP)
                     end do
                  else
                     do lm=1,n_lm_loc
                        !-- Use Geo as work array
                        Geo(lm)=dpdr(lm)-beta(nR)*p_Rdist(lm,nR)
                     end do
                  end if
                  call hIntRms(Geo,nR,1,n_lm_loc,0,Pre2hInt(:,nR),map_dist_st, &
                       &       .false.)
               end if
               call hIntRms(this%PFt2LM,nR,1,n_lmP_loc,1,Pre2hInt(:,nR), &
                    &       map_dist_st,.true.)
               call hIntRms(this%PFp2LM,nR,1,n_lmP_loc,1,Pre2hInt(:,nR), &
                    &       map_dist_st,.true.)

               if ( l_heat ) then
                  call hIntRms(Buo_temp,nR,1,n_lm_loc,0,Buo_temp2hInt(:,nR),&
                       &       map_dist_st,.false.)
               end if
               if ( l_chemical_conv ) then
                  call hIntRms(Buo_xi,nR,1,n_lm_loc,0,Buo_xi2hInt(:,nR),&
                       &       map_dist_st,.false.)
               end if
               if ( l_corr ) then
                  call hIntRms(CorPol,nR,1,n_lm_loc,0,Cor2hInt(:,nR),       &
                       &       map_dist_st,.false.)
                  call hIntRms(this%CFt2LM,nR,1,n_lmP_loc,1,Cor2hInt(:,nR), &
                       &       map_dist_st,.true.)
                  calL hIntRms(this%CFp2LM,nR,1,n_lmP_loc,1,Cor2hInt(:,nR), &
                       &       map_dist_st,.true.)
               end if
               if ( l_mag_LF .and. nR>n_r_LCR ) then
                  call hIntRms(LFPol,nR,1,n_lm_loc,0,LF2hInt(:,nR),         &
                       &       map_dist_st,.false.)
                  call hIntRms(this%LFt2LM,nR,1,n_lmP_loc,1,LF2hInt(:,nR),  &
                       &       map_dist_st,.true.)
                  call hIntRms(this%LFp2LM,nR,1,n_lmP_loc,1,LF2hInt(:,nR),  &
                       &       map_dist_st,.true.)
               end if

               do lm=1,n_lm_loc
                  lmP =lm2lmP(lm)
                  Geo(lm)=CorPol(lm)-dpdr(lm)+beta(nR)*p_Rdist(lm,nR)
                  PLF(lm)=LFPol(lm)-dpdr(lm)+beta(nR)*p_Rdist(lm,nR)
                  if ( l_adv_curl ) then
                     Geo(lm)=Geo(lm)+this%dpkindrLM(lmP)
                     PLF(lm)=PLF(lm)+this%dpkindrLM(lmP)
                  end if
                  CLF(lm)=CorPol(lm)+LFPol(lm)
                  Mag(lm)=Geo(lm)+LFPol(lm)
                  Arc(lm)=Geo(lm)+Buo_temp(lm)+Buo_xi(lm)
                  ArcMag(lm)=Mag(lm)+Buo_temp(lm)+Buo_xi(lm)
                  CIA(lm)=ArcMag(lm)+AdvPol(lm)-this%dtVrLM(lmP)
                  !CIA(lm)=CorPol(lm)+Buo_temp(lm)+Buo_xi(lm)+AdvPol(lm)
               end do
               call hIntRms(Geo,nR,1,n_lm_loc,0,Geo2hInt(:,nR),map_dist_st,.false.)
               call hIntRms(CLF,nR,1,n_lm_loc,0,CLF2hInt(:,nR),map_dist_st,.false.)
               call hIntRms(PLF,nR,1,n_lm_loc,0,PLF2hInt(:,nR),map_dist_st,.false.)
               call hIntRms(Mag,nR,1,n_lm_loc,0,Mag2hInt(:,nR),map_dist_st,.false.)
               call hIntRms(Arc,nR,1,n_lm_loc,0,Arc2hInt(:,nR),map_dist_st,.false.)
               call hIntRms(ArcMag,nR,1,n_lm_loc,0,ArcMag2hInt(:,nR),map_dist_st, &
                    &       .false.)
               call hIntRms(CIA,nR,1,n_lm_loc,0,CIA2hInt(:,nR),map_dist_st,.false.)

               do lm=1,n_lm_loc
                  lmP =lm2lmP(lm)
                  Geo(lm)=-this%CFt2LM(lmP)-this%PFt2LM(lmP)
                  CLF(lm)=-this%CFt2LM(lmP)+this%LFt2LM(lmP)
                  PLF(lm)=this%LFt2LM(lmP)-this%PFt2LM(lmP)
                  Mag(lm)=Geo(lm)+this%LFt2LM(lmP)
                  Arc(lm)=Geo(lm)
                  ArcMag(lm)=Mag(lm)
                  CIA(lm)=ArcMag(lm)+this%Advt2LM(lmP)-this%dtVtLM(lmP)
                  !CIA(lm)=-this%CFt2LM(lmP)+this%Advt2LM(lmP)
               end do
               call hIntRms(Geo,nR,1,n_lm_loc,0,Geo2hInt(:,nR),map_dist_st,.true.)
               call hIntRms(CLF,nR,1,n_lm_loc,0,CLF2hInt(:,nR),map_dist_st,.true.)
               call hIntRms(PLF,nR,1,n_lm_loc,0,PLF2hInt(:,nR),map_dist_st,.true.)
               call hIntRms(Mag,nR,1,n_lm_loc,0,Mag2hInt(:,nR),map_dist_st,.true.)
               call hIntRms(Arc,nR,1,n_lm_loc,0,Arc2hInt(:,nR),map_dist_st,.true.)
               call hIntRms(ArcMag,nR,1,n_lm_loc,0,ArcMag2hInt(:,nR),map_dist_st, &
                    &       .true.)
               call hIntRms(CIA,nR,1,n_lm_loc,0,CIA2hInt(:,nR),map_dist_st,.true.)

               do lm=1,n_lm_loc
                  lmP =lm2lmP(lm)
                  Geo(lm)=-this%CFp2LM(lmP)-this%PFp2LM(lmP)
                  CLF(lm)=-this%CFp2LM(lmP)+this%LFp2LM(lmP)
                  PLF(lm)=this%LFp2LM(lmP)-this%PFp2LM(lmP)
                  Mag(lm)=Geo(lm)+this%LFp2LM(lmP)
                  Arc(lm)=Geo(lm)
                  ArcMag(lm)=Mag(lm)
                  CIA(lm)=ArcMag(lm)+this%Advp2LM(lmP)-this%dtVpLM(lmP)
                  !CIA(lm)=-this%CFp2LM(lmP)+this%Advp2LM(lmP)
               end do
               call hIntRms(Geo,nR,1,n_lm_loc,0,Geo2hInt(:,nR),map_dist_st,.true.)
               call hIntRms(CLF,nR,1,n_lm_loc,0,CLF2hInt(:,nR),map_dist_st,.true.)
               call hIntRms(PLF,nR,1,n_lm_loc,0,PLF2hInt(:,nR),map_dist_st,.true.)
               call hIntRms(Mag,nR,1,n_lm_loc,0,Mag2hInt(:,nR),map_dist_st,.true.)
               call hIntRms(Arc,nR,1,n_lm_loc,0,Arc2hInt(:,nR),map_dist_st,.true.)
               call hIntRms(ArcMag,nR,1,n_lm_loc,0,ArcMag2hInt(:,nR),map_dist_st, &
                    &       .true.)
               call hIntRms(CIA,nR,1,n_lm_loc,0,CIA2hInt(:,nR),map_dist_st,.true.)

            end if

            ! In case double curl is calculated dpdt is useless
            if ( (.not. l_double_curl) .or. lPressNext ) then
            !if ( .true. ) then
               !
               !$omp parallel do default(shared) private(lm,l,m,lmS,lmA,lmP) &
               !$omp private(lmPS,AdvPol_loc,CorPol_loc)
               !LIKWID_ON('td_cv2')
               do lm=lm_maybe_skip_first,n_lm_loc
                  l   =lm2l(lm)
                  m   =lm2m(lm)
                  lmS =lm2lmS(lm)
                  lmA =lm2lmA(lm)
                  lmP =lm2lmP(lm)
                  lmPS=lmP2lmPS(lmP)
                  lmPA=lmP2lmPA(lmP)

                  !------ Recycle CorPol and AdvPol:
                  if ( l_corr ) then
                     !
                     if ( l < l_max .and. l > m ) then
                        CorPol_loc=            two*CorFac*or2(nR) *      &
                        &           ( -dPhi_loc(lm)  * ( dw_Rdist(lm,nR) &
                        &            +or1(nR)*dLh_loc(lm)*w_Rdist(lm,nR) &
                        &                                         )      &
                        &              +dTheta3A_loc(lm)*z_Rdist(lmA,nR) &
                        &              +dTheta3S_loc(lm)*z_Rdist(lmS,nR) &
                        &           )

                     else if ( l == l_max ) then
                        CorPol_loc=0.0_cp
                     else if ( l == m ) then
                        CorPol_loc=              two*CorFac*or2(nR) *       &
                        &              ( -dPhi_loc(lm)  * ( dw_Rdist(lm,nR) &
                        &               +or1(nR)*dLh_loc(lm)*w_Rdist(lm,nR) &
                        &                                                 ) &
                        &                +dTheta3A_loc(lm)*z_Rdist(lmA,nR)  &
                        &              )

                     end if
                     !
                  else
                     CorPol_loc=zero
                  end if
                  if ( l_conv_nl ) then
                     !
                     if ( l > m ) then
                        AdvPol_loc= dTheta1S_loc(lm)*this%AdvtLM(lmPS) - &
                        &           dTheta1A_loc(lm)*this%AdvtLM(lmPA) + &
                        &               dPhi_loc(lm)*this%AdvpLM(lmP)
                     else if ( l == m ) then
                        AdvPol_loc=-dTheta1A_loc(lm)*this%AdvtLM(lmPA) + &
                        &               dPhi_loc(lm)*this%AdvpLM(lmP)
                     end if
                     !
                  else
                     AdvPol_loc=zero
                  end if
                  dpdt(lm)=AdvPol_loc+CorPol_loc

               end do ! lm loop
               !LIKWID_OFF('td_cv2')
               !$omp end parallel do
               !
            end if

         else
            do lm=lm_maybe_skip_first,n_lm_loc
               dwdt(lm) =0.0_cp
               dzdt(lm) =0.0_cp
               dpdt(lm) =0.0_cp
            end do
         end if ! l_conv ?

      end if

      if ( nBc == 0 ) then

         if ( l_heat ) then
            lm = lm2(0,0)
            if ( lm > 0 ) then
               dsdt_loc  =epsc*epscProf(nR)!+opr/epsS*divKtemp0(nR)
               dVSrLM(lm)=this%VSrLM(lm)
               if ( l_anel ) then
                  if ( l_anelastic_liquid ) then
                     if ( l_mag_nl ) then
                        dsdt_loc=dsdt_loc+ViscHeatFac*temp0(nR)*this%ViscHeatLM(lm)+ &
                        &        OhmLossFac*temp0(nR)*this%OhmLossLM(lm)
                     else
                        dsdt_loc=dsdt_loc+ViscHeatFac*temp0(nR)*this%ViscHeatLM(lm)
                     end if
                  else
                     if ( l_mag_nl ) then
                        dsdt_loc=dsdt_loc+ViscHeatFac*this%ViscHeatLM(lm)+ &
                        &        OhmLossFac*this%OhmLossLM(lm)
                     else
                        dsdt_loc=dsdt_loc+ViscHeatFac*this%ViscHeatLM(lm)
                     end if
                  end if
               end if
               dsdt(lm)=dsdt_loc
            end if

            !$omp parallel do default(shared) private(lm,lmP,dsdt_loc,l,m)
            do lm=lm_maybe_skip_first,n_lm_loc
               l   =lm2l(lm)
               m   =lm2m(lm)
               lmP =lm2lmP(lm)
               dVSrLM(lm)=this%VSrLM(lmP)
               dsdt_loc = dLh_loc(lm)*this%VStLM(lmP)
               if ( l_anel ) then
                  if ( l_anelastic_liquid ) then
                     dsdt_loc = dsdt_loc+ViscHeatFac*hdif_V(l)* &
                     &          temp0(nR)*this%ViscHeatLM(lmP)
                     if ( l_mag_nl ) then
                        dsdt_loc = dsdt_loc+OhmLossFac*hdif_B(l) * &
                        &          temp0(nR)*this%OhmLossLM(lmP)
                     end if
                  else
                     dsdt_loc = dsdt_loc+ViscHeatFac*hdif_V(l) * this%ViscHeatLM(lmP)
                     if ( l_mag_nl ) then
                        dsdt_loc = dsdt_loc+OhmLossFac*hdif_B(l) * this%OhmLossLM(lmP)
                     end if
                  end if
               end if
               dsdt(lm) = dsdt_loc
            end do
            !$omp end parallel do
         else
            do lm=lm_maybe_skip_first,n_lm_loc
               dsdt(lm)  =0.0_cp
               dVSrLM(lm)=0.0_cp
            end do
         end if

         if ( l_chemical_conv ) then

            lm = lm2(0,0)
            if ( lm > 0 ) then
               dVXirLM(lm)=this%VXirLM(lm)
               dxidt(lm)  =epscXi
            end if

            !$omp parallel do default(shared) private(lm,lmP)
            do lm=lm_maybe_skip_first,n_lm_loc
               lmP=lm2lmP(lm)
               dVXirLM(lm)=this%VXirLM(lmP)
               dxidt(lm)  =dLh_loc(lm)*this%VXitLM(lmP)
            end do
            !$omp end parallel do
         end if

         if ( l_mag_nl .or. l_mag_kin  ) then
            !

            !$omp parallel do default(shared) private(lm,lmP)
            do lm=lm_maybe_skip_first,n_lm_loc
               lmP =lm2lmP(lm)
               dbdt(lm)   = dLh_loc(lm)*this%VxBpLM(lmP)
               dVxBhLM(lm)=-dLh_loc(lm)*this%VxBtLM(lmP)*r(nR)*r(nR)
               djdt(lm)   = dLh_loc(lm)*or4(nR)*this%VxBrLM(lmP)
            end do
            !$omp end parallel do
            !
         else
            if ( l_mag ) then
               do lm=1,n_lm_loc
                  dbdt(lm)   =zero
                  djdt(lm)   =zero
                  dVxBhLM(lm)=zero
               end do
            end if
         end if

      else   ! boundary !
         !
         if ( l_mag_nl .or. l_mag_kin ) then

            lm = lm2(0,0)
            if ( lm > 0 ) then
               dVxBhLM(lm)=zero
               dVSrLM(lm) =zero
            end if
            !$omp parallel do default(shared) private(lm,lmP)
            do lm=lm_maybe_skip_first,n_lm_loc
               lmP =lm2lmP(lm)
               dVxBhLM(lm)=-dLh_loc(lm)*this%VxBtLM(lmP)*r(nR)*r(nR)
               dVSrLM(lm) =zero
            end do
            !$omp end parallel do
         else
            do lm=1,n_lm_loc
               if ( l_mag ) dVxBhLM(lm)=zero
               dVSrLM(lm) =zero
            end do
         end if
         if ( l_double_curl ) then
            do lm=1,n_lm_loc
               dVxVhLM(lm)=zero
            end do
         end if
         if ( l_chemical_conv ) then
            do lm=1,n_lm_loc
               dVXirLM(lm)=zero
            end do
         end if
      end if  ! boundary ? lvelo ?

   end subroutine get_td
!-----------------------------------------------------------------------------
end module nonlinear_lm_mod

module nonlinear_3D_lm_mod

   use, intrinsic :: iso_c_binding
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: l_max, n_lm_loc, n_lmP_loc, nRstart, nRstop, n_r_loc
   use logic, only : l_anel, l_conv_nl, l_corr, l_heat, l_anelastic_liquid, &
       &             l_mag_nl, l_mag_kin, l_mag_LF, l_conv, l_mag, l_RMS,   &
       &             l_chemical_conv, l_single_matrix, l_double_curl,       &
       &             l_adv_curl
   use radial_functions, only: r, or2, or1, beta, rho0, rgrav, epscProf, &
       &                       or4, temp0, alpha0, ogrun, orho1, dLvisc, &
       &                       visc, dbeta
   use physical_parameters, only: CorFac, ra, epsc, ViscHeatFac, &
       &                          OhmLossFac, n_r_LCR, epscXi,   &
       &                          BuoFac, ChemFac, ThExpNb
   use horizontal_data, only: dLh_loc, dTheta1S_loc, dTheta1A_loc, dPhi_loc,    &
       &                      dTheta2A_loc, dTheta3A_loc, dTheta4A_loc,         &
       &                      dTheta2S_loc, hdif_B, dTheta3S_loc, dTheta4S_loc, &
       &                      hdif_V
   use RMS, only: Adv2hInt, Pre2hInt, Cor2hInt, LF2hInt, Iner2hInt,    &
       &          Geo2hInt, Mag2hInt, ArcMag2hInt, CLF2hInt, PLF2hInt, &
       &          CIA2hInt, Arc2hInt, Buo_temp2hInt, Buo_xi2hint
   use constants, only: zero, two, four, third, three
   use fields, only: w_Rdist, dw_Rdist, ddw_Rdist, z_Rdist, dz_Rdist, s_Rdist, &
       &             p_Rdist, dp_Rdist, xi_Rdist
   use RMS_helpers, only: hIntRms
   use LMmapping, only: map_dist_st


   implicit none

   type :: nonlinear_3D_lm_t
      !----- Nonlinear terms in lm-space:
      complex(cp), allocatable :: AdvrLM(:,:), AdvtLM(:,:), AdvpLM(:,:)
      complex(cp), allocatable :: LFrLM(:,:),  LFtLM(:,:),  LFpLM(:,:)
      complex(cp), allocatable :: VxBrLM(:,:), VxBtLM(:,:), VxBpLM(:,:)
      complex(cp), allocatable :: VSrLM(:,:),  VStLM(:,:),  VSpLM(:,:)
      complex(cp), allocatable :: VXirLM(:,:),  VXitLM(:,:),  VXipLM(:,:)
      complex(cp), allocatable :: ViscHeatLM(:,:), OhmLossLM(:,:)
      !----- RMS calculations
      complex(cp), allocatable :: Advt2LM(:,:), Advp2LM(:,:), PFt2LM(:,:), PFp2LM(:,:)
      complex(cp), allocatable :: LFt2LM(:,:), LFp2LM(:,:), CFt2LM(:,:), CFp2LM(:,:)
      complex(cp), allocatable :: dtVtLM(:,:), dtVpLM(:,:), dtVrLM(:,:)
      complex(cp), allocatable :: dpkindrLM(:,:)

   contains

      procedure :: initialize
      procedure :: finalize
      procedure :: output
      procedure :: set_zero
      procedure :: get_td

   end type nonlinear_3D_lm_t

contains

   subroutine initialize(this, n_loc, nRstart, nRstop)

      class(nonlinear_3D_lm_t) :: this
      integer, intent(in) :: n_loc
      integer, intent(in) :: nRstart
      integer, intent(in) :: nRstop

      allocate( this%AdvrLM(n_loc,nRstart:nRstop), this%AdvtLM(n_loc,nRstart:nRstop) )
      allocate( this%AdvpLM(n_loc,nRstart:nRstop), this%LFrLM(n_loc,nRstart:nRstop) )
      allocate( this%LFtLM(n_loc,nRstart:nRstop), this%LFpLM(n_loc,nRstart:nRstop) )
      allocate( this%VxBrLM(n_loc,nRstart:nRstop), this%VxBtLM(n_loc,nRstart:nRstop) )
      allocate( this%VxBpLM(n_loc,nRstart:nRstop))
      bytes_allocated = bytes_allocated + 9*n_loc*n_r_loc*SIZEOF_DEF_COMPLEX

      if ( l_anel) then
         allocate( this%ViscHeatLM(n_loc,nRstart:nRstop) )
         allocate( this%OhmLossLM(n_loc,nRstart:nRstop) )
         bytes_allocated = bytes_allocated+2*n_loc*n_r_loc*SIZEOF_DEF_COMPLEX
      end if

      if ( l_heat ) then
         allocate(this%VSrLM(n_loc,nRstart:nRstop),this%VStLM(n_loc,nRstart:nRstop))
         allocate(this%VSpLM(n_loc,nRstart:nRstop))
         bytes_allocated = bytes_allocated + 3*n_loc*n_r_loc*SIZEOF_DEF_COMPLEX
      end if

      if ( l_chemical_conv ) then
         allocate(this%VXirLM(n_loc,nRstart:nRstop),this%VXitLM(n_loc,nRstart:nRstop))
         allocate(this%VXipLM(n_loc,nRstart:nRstop))
         bytes_allocated = bytes_allocated + 3*n_loc*n_r_loc*SIZEOF_DEF_COMPLEX
      end if

      !-- RMS calculations
      if ( l_RMS ) then
         allocate(this%dtVrLM(n_loc,nRstart:nRstop),this%dtVtLM(n_loc,nRstart:nRstop))
         allocate(this%dtVpLM(n_loc,nRstart:nRstop))
         allocate(this%Advt2LM(n_loc,nRstart:nRstop))
         allocate(this%Advp2LM(n_loc,nRstart:nRstop))
         allocate(this%LFt2LM(n_loc,nRstart:nRstop),this%LFp2LM(n_loc,nRstart:nRstop))
         allocate(this%CFt2LM(n_loc,nRstart:nRstop),this%CFp2LM(n_loc,nRstart:nRstop))
         allocate(this%PFt2LM(n_loc,nRstart:nRstop),this%PFp2LM(n_loc,nRstart:nRstop))
         bytes_allocated = bytes_allocated + 11*n_loc*n_r_loc*SIZEOF_DEF_COMPLEX
         if ( l_adv_curl ) then
            allocate( this%dpkindrLM(n_loc,nRstart:nRstop) )
            bytes_allocated = bytes_allocated + n_loc*n_r_loc*SIZEOF_DEF_COMPLEX
         end if
      end if

   end subroutine initialize
!----------------------------------------------------------------------------
   subroutine finalize(this)

      class(nonlinear_3D_lm_t) :: this

      deallocate( this%AdvrLM, this%AdvtLM, this%AdvpLM )
      deallocate( this%LFrLM, this%LFtLM, this%LFpLM )
      deallocate( this%VxBrLM, this%VxBtLM, this%VxBpLM )

      if ( l_anel ) deallocate( this%ViscHeatLM, this%OhmLossLM )

      if ( l_chemical_conv ) deallocate( this%VXirLM, this%VXitLM, this%VXipLM )

      if ( l_heat ) deallocate( this%VSrLM, this%VStLM, this%VSpLM )

      !-- RMS calculations
      if ( l_RMS ) then
         deallocate( this%Advt2LM, this%Advp2LM, this%LFt2LM, this%LFp2LM )
         deallocate( this%CFt2LM, this%CFp2LM, this%PFt2LM, this%PFp2LM )
         deallocate( this%dtVrLM, this%dtVtLM, this%dtVpLM )
         if ( l_adv_curl ) deallocate ( this%dpkindrLM )
      end if

   end subroutine finalize
!----------------------------------------------------------------------------
   subroutine set_zero(this)

      class(nonlinear_3D_lm_t) :: this

      this%AdvrLM(:,:)=zero
      this%AdvtLM(:,:)=zero
      this%AdvpLM(:,:)=zero
      this%LFrLM(:,:) =zero
      this%LFtLM(:,:) =zero
      this%LFpLM(:,:) =zero
      this%VxBrLM(:,:)=zero
      this%VxBtLM(:,:)=zero
      this%VxBpLM(:,:)=zero
      if ( l_heat ) then
         this%VSrLM(:,:) =zero
         this%VStLM(:,:) =zero
         this%VSpLM(:,:) =zero
      end if
      if ( l_anel ) then
         this%ViscHeatLM(:,:)=zero
         this%OhmLossLM(:,:) =zero
      end if

      if ( l_chemical_conv ) then
         this%VXirLM(:,:)=zero
         this%VXitLM(:,:)=zero
         this%VXipLM(:,:)=zero
      end if

      if ( l_RMS ) then
         this%Advt2LM(:,:)=zero
         this%Advp2LM(:,:)=zero
         this%LFp2LM(:,:) =zero
         this%LFt2LM(:,:) =zero
         this%CFt2LM(:,:) =zero
         this%CFp2LM(:,:) =zero
         this%PFt2LM(:,:) =zero
         this%PFp2LM(:,:) =zero
         this%dtVtLM(:,:) =zero
         this%dtVpLM(:,:) =zero
         this%dtVrLM(:,:) =zero
         if ( l_adv_curl ) this%dpkindrLM(:,:)=zero
      end if

   end subroutine set_zero
!----------------------------------------------------------------------------
   subroutine output(this)

      class(nonlinear_3d_lm_t) :: this

      write(*,"(A,6ES20.12)") "AdvrLM,AdvtLM,AdvpLM = ",&
           & sum(this%AdvrLM), sum(this%AdvtLM),        &
           & sum(this%AdvpLM)

   end subroutine output
!----------------------------------------------------------------------------
   subroutine get_td(this,nR,nBc,lRmsCalc,lPressNext,dVSrLM,dVXirLM,    &
              &      dVxVhLM,dVxBhLM,dwdt,dzdt,dpdt,dsdt,dxidt,dbdt,djdt)
      !
      !  Purpose of this to calculate time derivatives
      !  dwdt,dzdt,dpdt,dsdt,dxidt,dbdt,djdt
      !  and auxiliary arrays dVSrLM, dVXirLM and dVxBhLM, dVxVhLM
      !  from non-linear terms in spectral form,
      !  contained in flmw1-3,flms1-3, flmb1-3 (input)
      !

      !-- Input of variables:
      class(nonlinear_3d_lm_t) :: this

      integer, intent(in) :: nR
      integer, intent(in) :: nBc ! signifies boundary conditions
      logical, intent(in) :: lRmsCalc
      logical, intent(in) :: lPressNext

      !-- Output of variables:
      complex(cp), intent(out) :: dwdt(:),dzdt(:)
      complex(cp), intent(out) :: dpdt(:),dsdt(:)
      complex(cp), intent(out) :: dxidt(:)
      complex(cp), intent(out) :: dbdt(:),djdt(:)
      complex(cp), intent(out) :: dVxBhLM(:)
      complex(cp), intent(out) :: dVxVhLM(:)
      complex(cp), intent(out) :: dVSrLM(:)
      complex(cp), intent(out) :: dVXirLM(:)

      !-- Local variables:
      integer :: l,m,lm,lmS,lmA,lmP,lmPS,lmPA
      complex(cp) :: CorPol(n_lm_loc), dpdr(n_lm_loc)
      complex(cp) :: AdvPol(n_lm_loc),AdvTor(n_lm_loc)
      complex(cp) :: LFPol(n_lm_loc),LFTor(n_lm_loc)
      complex(cp) :: Geo(n_lm_loc),CLF(n_lm_loc),PLF(n_lm_loc)
      complex(cp) :: ArcMag(n_lm_loc),Mag(n_lm_loc),CIA(n_lm_loc),Arc(n_lm_loc)
      complex(cp) :: Buo_temp(n_lm_loc), Buo_xi(n_lm_loc)
      complex(cp) :: AdvPol_loc,CorPol_loc,AdvTor_loc,CorTor_loc
      complex(cp) :: dsdt_loc
      complex(cp) :: dxidt_loc
      integer :: lm_maybe_skip_first
      integer, pointer :: lm2l(:), lm2m(:), lm2lmP(:), lm2(:,:)
      integer, pointer :: lmP2lmPS(:), lmP2lmPA(:), lm2lmA(:), lm2lmS(:)

      lm2l(1:n_lm_loc) => map_dist_st%lm2l
      lm2m(1:n_lm_loc) => map_dist_st%lm2m
      lmP2lmPS(1:n_lmP_loc) => map_dist_st%lmP2lmPS
      lmP2lmPA(1:n_lmP_loc) => map_dist_st%lmP2lmPA
      lm2lmP(1:n_lm_loc) => map_dist_st%lm2lmP
      lm2lmS(1:n_lm_loc) => map_dist_st%lm2lmS
      lm2lmA(1:n_lm_loc) => map_dist_st%lm2lmA
      lm2(0:,0:) => map_dist_st%lm2

      lm_maybe_skip_first = 1
      if (map_dist_st%lm2(0,0) > 0) lm_maybe_skip_first = 2

      if (nBc == 0 .or. lRmsCalc ) then

         if ( l_conv ) then  ! Convection

            lm =lm2(0,0)   ! This is l=0,m=0
            if ( lm > 0 ) then
               lmA=lm2lmA(lm)
               lmP=lm2lmP(lm)
               lmPA=lmP2lmPA(lmP)
               if ( l_conv_nl ) then
                  AdvPol_loc=      or2(nR)*this%AdvrLM(lm,nR)
                  AdvTor_loc=-dTheta1A_loc(lm)*this%AdvpLM(lmPA,nR)
               else
                  AdvPol_loc=zero
                  AdvTor_loc=zero
               end if
               if ( l_corr ) then
                  CorPol_loc=two*CorFac*or1(nR) * dTheta2A_loc(lm)*z_Rdist(lmA,nR)
                  CorTor_loc= two*CorFac*or2(nR) * (                      &
                  &                dTheta3A_loc(lm)*dw_Rdist(lmA,nR) +    &
                  &        or1(nR)*dTheta4A_loc(lm)* w_Rdist(lmA,nR) )
               else
                  CorPol_loc=zero
                  CorTor_loc=zero
               end if

               if ( l_single_matrix ) then
                  dwdt(lm)=AdvPol_loc!+CorPol_loc
               else
                  dwdt(lm)=AdvPol_loc+CorPol_loc
               end if

               dzdt(lm)=AdvTor_loc+CorTor_loc

               if ( lRmsCalc ) then
                  if (l_heat) then
                     Buo_temp(lm) =BuoFac*rgrav(nR)*rho0(nR)*s_Rdist(lm,nR)
                  else
                     Buo_temp(lm) =0.0_cp
                  end if
                  if (l_chemical_conv) then
                     Buo_xi(lm) =ChemFac*rgrav(nR)*rho0(nR)*xi_Rdist(lm,nR)
                  else
                     Buo_xi(lm)=0.0_cp
                  end if
                  if ( l_mag_LF .and. nR>n_r_LCR ) then
                     LFPol(lm) =      or2(nR)*this%LFrLM(lm,nR)
                     LFTor(lm) =-dTheta1A_loc(lm)*this%LFpLM(lmPA,nR)
                     AdvPol(lm)=AdvPol_loc-LFPol(lm)
                     AdvTor(lm)=AdvTor_loc-LFTor(lm)
                  else
                     AdvPol(lm)=AdvPol_loc
                     AdvTor(lm)=AdvTor_loc
                  end if
                  CorPol(lm)=CorPol_loc

                  if ( l_double_curl ) then
                     !-- Recalculate the pressure gradient based on the poloidal
                     !-- equation equilibrium
                     dpdr(lm)=Buo_temp(lm)+Buo_xi(lm)+beta(nR)*p_Rdist(lm,nR)+ &
                     &        AdvPol_loc+CorPol_loc
                  else
                     dpdr(lm)=dp_Rdist(lm,nR)
                  end if

               end if ! lRmsCalc
            end if

            !
            !$omp parallel do default(shared) private(lm,l,m,lmS,lmA,lmP) &
            !$omp private(lmPS,lmPA,AdvPol_loc,CorPol_loc,AdvTor_loc,CorTor_loc)
            do lm=lm_maybe_skip_first,n_lm_loc
               l   =lm2l(lm)
               m   =lm2m(lm)
               lmS =lm2lmS(lm)
               lmA =lm2lmA(lm)
               lmP =lm2lmP(lm)
               lmPS=lmP2lmPS(lmP)
               lmPA=lmP2lmPA(lmP)

               if ( l_double_curl ) then ! Pressure is not needed

                  if ( l_corr ) then
                     if ( l < l_max .and. l > m ) then
                        CorPol_loc =two*CorFac*or2(nR)*orho1(nR)*(               &
                        &                dPhi_loc(lm)*(                          &
                        &        -ddw_Rdist(lm,nR)+beta(nR)*dw_Rdist(lm,nR)    + &
                        &             ( beta(nR)*or1(nR)+or2(nR))*               &
                        &                        dLh_loc(lm)*w_Rdist(lm,nR) )  + &
                        &         dTheta3A_loc(lm)*( dz_Rdist(lmA,nR)-           &
                        &                           beta(nR)*z_Rdist(lmA,nR) ) + &
                        &         dTheta3S_loc(lm)*( dz_Rdist(lmS,nR)-           &
                        &                           beta(nR)*z_Rdist(lmS,nR) ) + &
                        &          or1(nR)* (                                    &
                        &             dTheta4A_loc(lm)* z_Rdist(lmA,nR)          &
                        &            -dTheta4S_loc(lm)* z_Rdist(lmS,nR) ) )
                     else if ( l == l_max ) then
                        CorPol_loc =0.0_cp
                     else if ( l == m ) then
                        CorPol_loc =two*CorFac*or2(nR)*orho1(nR)*(               &
                        &                dPhi_loc(lm)*(                          &
                        &       -ddw_Rdist(lm,nR)+beta(nR)*dw_Rdist(lm,nR)     + &
                        &             ( beta(nR)*or1(nR)+or2(nR))*               &
                        &                       dLh_loc(lm)*w_Rdist(lm,nR) )   + &
                        &              dTheta3A_loc(lm)*( dz_Rdist(lmA,nR)-      &
                        &                         beta(nR)*z_Rdist(lmA,nR) )   + &
                        &          or1(nR)* (                                    &
                        &              dTheta4A_loc(lm)*   z_Rdist(lmA,nR) ) )
                     end if
                  else
                     CorPol_loc=zero
                  end if

                  if ( l_conv_nl ) then

                     if ( l > m ) then
                        dVxVhLM(lm)=      orho1(nR)*r(nR)*r(nR)* (        &
                        &        dTheta1S_loc(lm)*this%AdvtLM(lmPS,nR) -  &
                        &        dTheta1A_loc(lm)*this%AdvtLM(lmPA,nR) +  &
                        &             dPhi_loc(lm)*this%AdvpLM(lmP,nR)  )
                     else if ( l == m ) then
                        dVxVhLM(lm)=      orho1(nR)*r(nR)*r(nR)* (        &
                        &      - dTheta1A_loc(lm)*this%AdvtLM(lmPA,nR) +  &
                        &        dPhi_loc(lm)*this%AdvpLM(lmP,nR)  )
                     end if

                     AdvPol_loc=dLh_loc(lm)*or4(nR)*orho1(nR)*this%AdvrLM(lmP,nR)

                  else

                     AdvPol_loc =zero
                     dVxVhLM(lm)=zero

                  endif

               else ! We don't use the double curl

                  if ( l_corr .and. nBc /= 2 ) then
                     if ( l < l_max .and. l > m ) then
                        CorPol_loc =two*CorFac*or1(nR) * (      &
                        &       dPhi_loc(lm)*dw_Rdist(lm,nR) +  & ! phi-deriv of dw/dr
                        &   dTheta2A_loc(lm)*z_Rdist(lmA,nR) -  & ! sin(theta) dtheta z
                        &    dTheta2S_loc(lm)*z_Rdist(lmS,nR) )
                     else if ( l == l_max ) then
                        CorPol_loc=0.0_cp
                     else if ( l == m ) then
                        CorPol_loc = two*CorFac*or1(nR) * (      &
                        &        dPhi_loc(lm)*dw_Rdist(lm,nR)  + &
                        &    dTheta2A_loc(lm)*z_Rdist(lmA,nR) )
                     end if
                  else
                     CorPol_loc=zero
                  end if

                  if ( l_conv_nl ) then
                     AdvPol_loc=or2(nR)*this%AdvrLM(lmP,nR)
                  else
                     AdvPol_loc=zero
                  endif

               end if ! Double curl or not for the poloidal equation

               dwdt(lm)=AdvPol_loc+CorPol_loc

               if ( lRmsCalc ) then ! RMS force balance

                  if ( l_anelastic_liquid ) then
                     if (l_heat) then
                        Buo_temp(lm) =BuoFac*alpha0(nR)*rgrav(nR)*(      &
                        &      rho0(nR)*s_Rdist(lm,nR)-ViscHeatFac*      &
                        &     (ThExpNb*alpha0(nR)*temp0(nR)+ogrun(nR))*  &
                        &     p_Rdist(lm,nR) )
                     else
                        Buo_temp(lm) =0.0_cp
                     end if
                     if (l_chemical_conv) then
                        Buo_xi(lm) =ChemFac*alpha0(nR)*rgrav(nR)*(       &
                        &    rho0(nR)*xi_Rdist(lm,nR)-ViscHeatFac*       &
                        &     (ThExpNb*alpha0(nR)*temp0(nR)+ogrun(nR))*  &
                        &     p_Rdist(lm,nR) )
                     else
                        Buo_xi(lm) =0.0_cp
                     end if
                  else
                     if (l_heat) then
                        Buo_temp(lm) =BuoFac*rho0(nR)*rgrav(nR)*s_Rdist(lm,nR)
                     else
                        Buo_temp(lm) =0.0_cp
                     end if
                     if (l_chemical_conv) then
                        Buo_xi(lm) =ChemFac*rho0(nR)*rgrav(nR)*xi_Rdist(lm,nR)
                     else
                        Buo_xi(lm)=0.0_cp
                     end if
                  end if

                  if ( l_double_curl ) then
                     ! In that case we have to recompute the Coriolis force
                     ! since we also want the pressure gradient
                     if ( l_corr .and. nBc /= 2 ) then
                        if ( l < l_max .and. l > m ) then
                           CorPol_loc = two*CorFac*or1(nR) * (      &
                           &        dPhi_loc(lm)*dw_Rdist(lm,nR) +  & ! phi-deriv of dw/dr
                           &    dTheta2A_loc(lm)*z_Rdist(lmA,nR) -  & ! sin(theta) dtheta z
                           &    dTheta2S_loc(lm)*z_Rdist(lmS,nR) )
                        else if ( l == l_max ) then
                           CorPol_loc= 0.0_cp
                        else if ( l == m ) then
                           CorPol_loc = two*CorFac*or1(nR) * (       &
                           &         dPhi_loc(lm)*dw_Rdist(lm,nR)  + &
                           &     dTheta2A_loc(lm)*z_Rdist(lmA,nR) )
                        end if
                     else
                        CorPol_loc=zero
                     end if

                     ! We also need to recompute AdvPol_loc here
                     if ( l_conv_nl ) then
                        AdvPol_loc=or2(nR)*this%AdvrLM(lmP,nR)
                     else
                        AdvPol_loc=zero
                     endif

                     !-- Recalculate the pressure gradient based on the poloidal
                     !-- equation equilibrium
                     dpdr(lm)=Buo_temp(lm)+Buo_xi(lm)-this%dtVrLM(lmP,nR)+       &
                     &      dLh_loc(lm)*or2(nR)*hdif_V(l)*visc(nR)* (            &
                     &                                        ddw_Rdist(lm,nR)+  &
                     &         (two*dLvisc(nR)-third*beta(nR))*dw_Rdist(lm,nR)-  &
                     &      ( dLh_loc(lm)*or2(nR)+four*third*( dbeta(nR)+        &
                     &          dLvisc(nR)*beta(nR)+(three*dLvisc(nR)+beta(nR))* &
                     &         or1(nR)) )*                      w_Rdist(lm,nR))+ &
                     &        beta(nR)*p_Rdist(lm,nR)+AdvPol_loc+CorPol_loc
                  else
                     dpdr(lm)=dp_Rdist(lm,nR)
                  end if

                  if ( l_mag_LF .and. nR>n_r_LCR ) then
                     LFPol(lm) =or2(nR)*this%LFrLM(lmP,nR)
                     AdvPol(lm)=AdvPol_loc-LFPol(lm)
                  else
                     AdvPol(lm)=AdvPol_loc
                  end if
                  CorPol(lm)=CorPol_loc

               end if ! lRmsCalc

               if ( l_corr ) then
                  if ( l < l_max .and. l > m ) then
                     CorTor_loc=           two*CorFac*or2(nR) * (      &
                     &                 dPhi_loc(lm)*z_Rdist(lm,nR)   + &
                     &            dTheta3A_loc(lm)*dw_Rdist(lmA,nR)  + &
                     &    or1(nR)*dTheta4A_loc(lm)* w_Rdist(lmA,nR)  + &
                     &            dTheta3S_loc(lm)*dw_Rdist(lmS,nR)  - &
                     &    or1(nR)*dTheta4S_loc(lm)* w_Rdist(lmS,nR)  )
                  else if ( l == l_max ) then
                     CorTor_loc=0.0_cp
                  else if ( l == m ) then
                     CorTor_loc=  two*CorFac*or2(nR) * (       &
                     &         dPhi_loc(lm)*z_Rdist(lm,nR)   + &
                     &    dTheta3A_loc(lm)*dw_Rdist(lmA,nR)  + &
                     &    or1(nR)*dTheta4A_loc(lm)* w_Rdist(lmA,nR)  )
                  end if
               else
                  CorTor_loc=zero
               end if

               if ( l_conv_nl ) then
                  if ( l > m ) then
                     AdvTor_loc=   -dPhi_loc(lm)*this%AdvtLM(lmP,nR)  + &
                     &          dTheta1S_loc(lm)*this%AdvpLM(lmPS,nR) - &
                     &          dTheta1A_loc(lm)*this%AdvpLM(lmPA,nR)
                  else if ( l == m ) then
                     AdvTor_loc=   -dPhi_loc(lm)*this%AdvtLM(lmP,nR)  - &
                     &          dTheta1A_loc(lm)*this%AdvpLM(lmPA,nR)
                  end if
               else
                  AdvTor_loc=zero
               end if

               dzdt(lm)=CorTor_loc+AdvTor_loc

               if ( lRmsCalc ) then
                  if ( l_mag_LF .and. nR>n_r_LCR ) then
                     !------ When RMS values are required, the Lorentz force is treated
                     !       separately:

                     if ( l > m ) then
                        !------- LFTor= 1/(E*Pm) * curl( curl(B) x B )_r
                        LFTor(lm) =   -dPhi_loc(lm)*this%LFtLM(lmP,nR)  + &
                        &          dTheta1S_loc(lm)*this%LFpLM(lmPS,nR) - &
                        &          dTheta1A_loc(lm)*this%LFpLM(lmPA,nR)
                     else if ( l == m ) then
                        LFTor(lm) =   -dPhi_loc(lm)*this%LFtLM(lmP,nR)  - &
                        &          dTheta1A_loc(lm)*this%LFpLM(lmPA,nR)
                     end if
                     AdvTor(lm)=AdvTor_loc-LFTor(lm)
                  else
                     AdvTor(lm)=AdvTor_loc
                  end if
               end if

            end do
            !$omp end parallel do
            !

            if ( lRmsCalc ) then

               if ( l_conv_nl ) then
                  if ( l_adv_curl ) then
                     do lm=1,n_lm_loc
                        lmP =lm2lmP(lm)
                        AdvPol(lm)=AdvPol(lm)-this%dpkindrLM(lmP,nR)
                     end do
                  end if
                  call hIntRms(AdvPol,nR,1,n_lm_loc,0,Adv2hInt(:,nR),map_dist_st, &
                       &       .false.)
                  call hIntRms(this%Advt2LM(:,nR),nR,1,n_lmP_loc,1,Adv2hInt(:,nR), &
                       &       map_dist_st, .true.)
                  call hIntRms(this%Advp2LM(:,nR),nR,1,n_lMP_loc,1,Adv2hInt(:,nR), &
                       &       map_dist_st, .true.)
                  do lm=1,n_lm_loc
                     lmP =lm2lmP(lm)
                     !-- Use Geo as work array
                     Geo(lm)=AdvPol(lm)-this%dtVrLM(lmP,nR)
                  end do
                  call hIntRms(Geo,nR,1,n_lm_loc,0,Iner2hInt(:,nR),map_dist_st, &
                       &       .false.)
                  call hIntRms(this%Advt2LM(:,nR)-this%dtVtLM(:,nR),nR,1,n_lmP_loc,1,&
                       &       Iner2hInt(:,nR),map_dist_st,.true.)
                  call hIntRms(this%Advp2LM(:,nR)-this%dtVpLM(:,nR),nR,1,n_lmP_loc,1,&
                       &       Iner2hInt(:,nR),map_dist_st,.true.)
               end if

               if ( l_anelastic_liquid ) then
                  call hIntRms(dpdr,nR,1,n_lm_loc,0,Pre2hInt(:,nR),map_dist_st, &
                       &       .false.)
               else
                  ! rho* grad(p/rho) = grad(p) - beta*p
                  !-- Geo is used to store the pressure Gradient
                  if ( l_adv_curl ) then
                     do lm=1,n_lm_loc
                        lmP =lm2lmP(lm)
                        !-- Use Geo as work array
                        Geo(lm)=dpdr(lm)-beta(nR)*p_Rdist(lm,nR)- &
                        &       this%dpkindrLM(lmP,nR)
                     end do
                  else
                     do lm=1,n_lm_loc
                        !-- Use Geo as work array
                        Geo(lm)=dpdr(lm)-beta(nR)*p_Rdist(lm,nR)
                     end do
                  end if
                  call hIntRms(Geo,nR,1,n_lm_loc,0,Pre2hInt(:,nR),map_dist_st, &
                       &       .false.)
               end if
               call hIntRms(this%PFt2LM(:,nR),nR,1,n_lmP_loc,1,Pre2hInt(:,nR), &
                    &       map_dist_st,.true.)
               call hIntRms(this%PFp2LM(:,nR),nR,1,n_lmP_loc,1,Pre2hInt(:,nR), &
                    &       map_dist_st,.true.)

               if ( l_heat ) then
                  call hIntRms(Buo_temp,nR,1,n_lm_loc,0,Buo_temp2hInt(:,nR),&
                       &       map_dist_st,.false.)
               end if
               if ( l_chemical_conv ) then
                  call hIntRms(Buo_xi,nR,1,n_lm_loc,0,Buo_xi2hInt(:,nR),&
                       &       map_dist_st,.false.)
               end if
               if ( l_corr ) then
                  call hIntRms(CorPol,nR,1,n_lm_loc,0,Cor2hInt(:,nR),       &
                       &       map_dist_st,.false.)
                  call hIntRms(this%CFt2LM(:,nR),nR,1,n_lmP_loc,1,Cor2hInt(:,nR), &
                       &       map_dist_st,.true.)
                  calL hIntRms(this%CFp2LM(:,nR),nR,1,n_lmP_loc,1,Cor2hInt(:,nR), &
                       &       map_dist_st,.true.)
               end if
               if ( l_mag_LF .and. nR>n_r_LCR ) then
                  call hIntRms(LFPol,nR,1,n_lm_loc,0,LF2hInt(:,nR),         &
                       &       map_dist_st,.false.)
                  call hIntRms(this%LFt2LM(:,nR),nR,1,n_lmP_loc,1,LF2hInt(:,nR),  &
                       &       map_dist_st,.true.)
                  call hIntRms(this%LFp2LM(:,nR),nR,1,n_lmP_loc,1,LF2hInt(:,nR),  &
                       &       map_dist_st,.true.)
               end if

               do lm=1,n_lm_loc
                  lmP =lm2lmP(lm)
                  Geo(lm)=CorPol(lm)-dpdr(lm)+beta(nR)*p_Rdist(lm,nR)
                  PLF(lm)=LFPol(lm)-dpdr(lm)+beta(nR)*p_Rdist(lm,nR)
                  if ( l_adv_curl ) then
                     Geo(lm)=Geo(lm)+this%dpkindrLM(lmP,nR)
                     PLF(lm)=PLF(lm)+this%dpkindrLM(lmP,nR)
                  end if
                  CLF(lm)=CorPol(lm)+LFPol(lm)
                  Mag(lm)=Geo(lm)+LFPol(lm)
                  Arc(lm)=Geo(lm)+Buo_temp(lm)+Buo_xi(lm)
                  ArcMag(lm)=Mag(lm)+Buo_temp(lm)+Buo_xi(lm)
                  CIA(lm)=ArcMag(lm)+AdvPol(lm)-this%dtVrLM(lmP,nR)
                  !CIA(lm)=CorPol(lm)+Buo_temp(lm)+Buo_xi(lm)+AdvPol(lm)
               end do
               call hIntRms(Geo,nR,1,n_lm_loc,0,Geo2hInt(:,nR),map_dist_st,.false.)
               call hIntRms(CLF,nR,1,n_lm_loc,0,CLF2hInt(:,nR),map_dist_st,.false.)
               call hIntRms(PLF,nR,1,n_lm_loc,0,PLF2hInt(:,nR),map_dist_st,.false.)
               call hIntRms(Mag,nR,1,n_lm_loc,0,Mag2hInt(:,nR),map_dist_st,.false.)
               call hIntRms(Arc,nR,1,n_lm_loc,0,Arc2hInt(:,nR),map_dist_st,.false.)
               call hIntRms(ArcMag,nR,1,n_lm_loc,0,ArcMag2hInt(:,nR),map_dist_st, &
                    &       .false.)
               call hIntRms(CIA,nR,1,n_lm_loc,0,CIA2hInt(:,nR),map_dist_st,.false.)

               do lm=1,n_lm_loc
                  lmP =lm2lmP(lm)
                  Geo(lm)=-this%CFt2LM(lmP,nR)-this%PFt2LM(lmP,nR)
                  CLF(lm)=-this%CFt2LM(lmP,nR)+this%LFt2LM(lmP,nR)
                  PLF(lm)=this%LFt2LM(lmP,nR)-this%PFt2LM(lmP,nR)
                  Mag(lm)=Geo(lm)+this%LFt2LM(lmP,nR)
                  Arc(lm)=Geo(lm)
                  ArcMag(lm)=Mag(lm)
                  CIA(lm)=ArcMag(lm)+this%Advt2LM(lmP,nR)-this%dtVtLM(lmP,nR)
                  !CIA(lm)=-this%CFt2LM(lmP,nR)+this%Advt2LM(lmP,nR)
               end do
               call hIntRms(Geo,nR,1,n_lm_loc,0,Geo2hInt(:,nR),map_dist_st,.true.)
               call hIntRms(CLF,nR,1,n_lm_loc,0,CLF2hInt(:,nR),map_dist_st,.true.)
               call hIntRms(PLF,nR,1,n_lm_loc,0,PLF2hInt(:,nR),map_dist_st,.true.)
               call hIntRms(Mag,nR,1,n_lm_loc,0,Mag2hInt(:,nR),map_dist_st,.true.)
               call hIntRms(Arc,nR,1,n_lm_loc,0,Arc2hInt(:,nR),map_dist_st,.true.)
               call hIntRms(ArcMag,nR,1,n_lm_loc,0,ArcMag2hInt(:,nR),map_dist_st, &
                    &       .true.)
               call hIntRms(CIA,nR,1,n_lm_loc,0,CIA2hInt(:,nR),map_dist_st,.true.)

               do lm=1,n_lm_loc
                  lmP =lm2lmP(lm)
                  Geo(lm)=-this%CFp2LM(lmP,nR)-this%PFp2LM(lmP,nR)
                  CLF(lm)=-this%CFp2LM(lmP,nR)+this%LFp2LM(lmP,nR)
                  PLF(lm)=this%LFp2LM(lmP,nR)-this%PFp2LM(lmP,nR)
                  Mag(lm)=Geo(lm)+this%LFp2LM(lmP,nR)
                  Arc(lm)=Geo(lm)
                  ArcMag(lm)=Mag(lm)
                  CIA(lm)=ArcMag(lm)+this%Advp2LM(lmP,nR)-this%dtVpLM(lmP,nR)
                  !CIA(lm)=-this%CFp2LM(lmP,nR)+this%Advp2LM(lmP,nR)
               end do
               call hIntRms(Geo,nR,1,n_lm_loc,0,Geo2hInt(:,nR),map_dist_st,.true.)
               call hIntRms(CLF,nR,1,n_lm_loc,0,CLF2hInt(:,nR),map_dist_st,.true.)
               call hIntRms(PLF,nR,1,n_lm_loc,0,PLF2hInt(:,nR),map_dist_st,.true.)
               call hIntRms(Mag,nR,1,n_lm_loc,0,Mag2hInt(:,nR),map_dist_st,.true.)
               call hIntRms(Arc,nR,1,n_lm_loc,0,Arc2hInt(:,nR),map_dist_st,.true.)
               call hIntRms(ArcMag,nR,1,n_lm_loc,0,ArcMag2hInt(:,nR),map_dist_st, &
                    &       .true.)
               call hIntRms(CIA,nR,1,n_lm_loc,0,CIA2hInt(:,nR),map_dist_st,.true.)

            end if

            ! In case double curl is calculated dpdt is useless
            if ( (.not. l_double_curl) .or. lPressNext ) then
            !if ( .true. ) then
               !
               !$omp parallel do default(shared) private(lm,l,m,lmS,lmA,lmP) &
               !$omp private(lmPS,AdvPol_loc,CorPol_loc)
               !LIKWID_ON('td_cv2')
               do lm=lm_maybe_skip_first,n_lm_loc
                  l   =lm2l(lm)
                  m   =lm2m(lm)
                  lmS =lm2lmS(lm)
                  lmA =lm2lmA(lm)
                  lmP =lm2lmP(lm)
                  lmPS=lmP2lmPS(lmP)
                  lmPA=lmP2lmPA(lmP)

                  !------ Recycle CorPol and AdvPol:
                  if ( l_corr ) then
                     !
                     if ( l < l_max .and. l > m ) then
                        CorPol_loc=            two*CorFac*or2(nR) *      &
                        &           ( -dPhi_loc(lm)  * ( dw_Rdist(lm,nR) &
                        &            +or1(nR)*dLh_loc(lm)*w_Rdist(lm,nR) &
                        &                                         )      &
                        &              +dTheta3A_loc(lm)*z_Rdist(lmA,nR) &
                        &              +dTheta3S_loc(lm)*z_Rdist(lmS,nR) &
                        &           )

                     else if ( l == l_max ) then
                        CorPol_loc=0.0_cp
                     else if ( l == m ) then
                        CorPol_loc=              two*CorFac*or2(nR) *       &
                        &              ( -dPhi_loc(lm)  * ( dw_Rdist(lm,nR) &
                        &               +or1(nR)*dLh_loc(lm)*w_Rdist(lm,nR) &
                        &                                                 ) &
                        &                +dTheta3A_loc(lm)*z_Rdist(lmA,nR)  &
                        &              )

                     end if
                     !
                  else
                     CorPol_loc=zero
                  end if
                  if ( l_conv_nl ) then
                     !
                     if ( l > m ) then
                        AdvPol_loc= dTheta1S_loc(lm)*this%AdvtLM(lmPS,nR) - &
                        &           dTheta1A_loc(lm)*this%AdvtLM(lmPA,nR) + &
                        &               dPhi_loc(lm)*this%AdvpLM(lmP,nR)
                     else if ( l == m ) then
                        AdvPol_loc=-dTheta1A_loc(lm)*this%AdvtLM(lmPA,nR) + &
                        &               dPhi_loc(lm)*this%AdvpLM(lmP,nR)
                     end if
                     !
                  else
                     AdvPol_loc=zero
                  end if
                  dpdt(lm)=AdvPol_loc+CorPol_loc

               end do ! lm loop
               !LIKWID_OFF('td_cv2')
               !$omp end parallel do
               !
            end if

         else
            do lm=lm_maybe_skip_first,n_lm_loc
               dwdt(lm) =0.0_cp
               dzdt(lm) =0.0_cp
               dpdt(lm) =0.0_cp
            end do
         end if ! l_conv ?

      end if

      if ( nBc == 0 ) then

         if ( l_heat ) then
            lm = lm2(0,0)
            if ( lm > 0 ) then
               dsdt_loc  =epsc*epscProf(nR)!+opr/epsS*divKtemp0(nR)
               dVSrLM(lm)=this%VSrLM(lm,nR)
               if ( l_anel ) then
                  if ( l_anelastic_liquid ) then
                     if ( l_mag_nl ) then
                        dsdt_loc=dsdt_loc+ViscHeatFac*temp0(nR)*this%ViscHeatLM(lm,nR)+ &
                        &        OhmLossFac*temp0(nR)*this%OhmLossLM(lm,nR)
                     else
                        dsdt_loc=dsdt_loc+ViscHeatFac*temp0(nR)*this%ViscHeatLM(lm,nR)
                     end if
                  else
                     if ( l_mag_nl ) then
                        dsdt_loc=dsdt_loc+ViscHeatFac*this%ViscHeatLM(lm,nR)+ &
                        &        OhmLossFac*this%OhmLossLM(lm,nR)
                     else
                        dsdt_loc=dsdt_loc+ViscHeatFac*this%ViscHeatLM(lm,nR)
                     end if
                  end if
               end if
               dsdt(lm)=dsdt_loc
            end if

            !$omp parallel do default(shared) private(lm,lmP,dsdt_loc,l,m)
            do lm=lm_maybe_skip_first,n_lm_loc
               l   =lm2l(lm)
               m   =lm2m(lm)
               lmP =lm2lmP(lm)
               dVSrLM(lm)=this%VSrLM(lmP,nR)
               dsdt_loc = dLh_loc(lm)*this%VStLM(lmP,nR)
               if ( l_anel ) then
                  if ( l_anelastic_liquid ) then
                     dsdt_loc = dsdt_loc+ViscHeatFac*hdif_V(l)* &
                     &          temp0(nR)*this%ViscHeatLM(lmP,nR)
                     if ( l_mag_nl ) then
                        dsdt_loc = dsdt_loc+OhmLossFac*hdif_B(l) * &
                        &          temp0(nR)*this%OhmLossLM(lmP,nR)
                     end if
                  else
                     dsdt_loc = dsdt_loc+ViscHeatFac*hdif_V(l)*this%ViscHeatLM(lmP,nR)
                     if ( l_mag_nl ) then
                        dsdt_loc = dsdt_loc+OhmLossFac*hdif_B(l)*this%OhmLossLM(lmP,nR)
                     end if
                  end if
               end if
               dsdt(lm) = dsdt_loc
            end do
            !$omp end parallel do
         else
            do lm=lm_maybe_skip_first,n_lm_loc
               dsdt(lm)  =0.0_cp
               dVSrLM(lm)=0.0_cp
            end do
         end if

         if ( l_chemical_conv ) then

            lm = lm2(0,0)
            if ( lm > 0 ) then
               dVXirLM(lm)=this%VXirLM(lm,nR)
               dxidt(lm)  =epscXi
            end if

            !$omp parallel do default(shared) private(lm,lmP)
            do lm=lm_maybe_skip_first,n_lm_loc
               lmP=lm2lmP(lm)
               dVXirLM(lm)=this%VXirLM(lmP,nR)
               dxidt(lm)  =dLh_loc(lm)*this%VXitLM(lmP,nR)
            end do
            !$omp end parallel do
         end if

         if ( l_mag_nl .or. l_mag_kin  ) then
            !

            !$omp parallel do default(shared) private(lm,lmP)
            do lm=lm_maybe_skip_first,n_lm_loc
               lmP =lm2lmP(lm)
               dbdt(lm)   = dLh_loc(lm)*this%VxBpLM(lmP,nR)
               dVxBhLM(lm)=-dLh_loc(lm)*this%VxBtLM(lmP,nR)*r(nR)*r(nR)
               djdt(lm)   = dLh_loc(lm)*or4(nR)*this%VxBrLM(lmP,nR)
            end do
            !$omp end parallel do
            !
         else
            if ( l_mag ) then
               do lm=1,n_lm_loc
                  dbdt(lm)   =zero
                  djdt(lm)   =zero
                  dVxBhLM(lm)=zero
               end do
            end if
         end if

      else   ! boundary !
         !
         if ( l_mag_nl .or. l_mag_kin ) then

            lm = lm2(0,0)
            if ( lm > 0 ) then
               dVxBhLM(lm)=zero
               dVSrLM(lm) =zero
            end if
            !$omp parallel do default(shared) private(lm,lmP)
            do lm=lm_maybe_skip_first,n_lm_loc
               lmP =lm2lmP(lm)
               dVxBhLM(lm)=-dLh_loc(lm)*this%VxBtLM(lmP,nR)*r(nR)*r(nR)
               dVSrLM(lm) =zero
            end do
            !$omp end parallel do
         else
            do lm=1,n_lm_loc
               if ( l_mag ) dVxBhLM(lm)=zero
               dVSrLM(lm) =zero
            end do
         end if
         if ( l_double_curl ) then
            do lm=1,n_lm_loc
               dVxVhLM(lm)=zero
            end do
         end if
         if ( l_chemical_conv ) then
            do lm=1,n_lm_loc
               dVXirLM(lm)=zero
            end do
         end if
      end if  ! boundary ? lvelo ?

   end subroutine get_td
!-----------------------------------------------------------------------------
end module nonlinear_3D_lm_mod
