#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif

module nonlinear_lm_mod

   use, intrinsic :: iso_c_binding
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: lm_max, l_max, lm_maxMag, lmP_max
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
   use blocking, only: lm2l, lm2m, lm2lmP, lmP2lmPS, lmP2lmPA, lm2lmA, &
       &               lm2lmS, st_map, lo_map
   use horizontal_data, only: dLh, dTheta1S, dTheta1A, dPhi, dTheta2A, &
       &                      dTheta3A, dTheta4A, dTheta2S,  hdif_B,   &
       &                      dTheta3S, dTheta4S, hdif_V
   use RMS, only: Adv2hInt, Pre2hInt, Cor2hInt, LF2hInt, Iner2hInt,    &
       &          Geo2hInt, Mag2hInt, ArcMag2hInt, CLF2hInt, PLF2hInt, &
       &          CIA2hInt, Arc2hInt, Buo_temp2hInt, Buo_xi2hint
   use constants, only: zero, two, four, third, three
   use fields, only: w_Rloc, dw_Rloc, ddw_Rloc, z_Rloc, dz_Rloc, s_Rloc, &
       &             p_Rloc, dp_Rloc, xi_Rloc
   use RMS_helpers, only: hIntRms
   use communications, only: slice_FlmP_cmplx, gather_FlmP


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
      procedure :: output
      procedure :: set_zero
      procedure :: get_td
      procedure :: slice_all   !@>TODO temporary function to help transition, delete me!
      procedure :: gather_all  !@>TODO temporary function to help transition, delete me!

   end type nonlinear_lm_t

contains

   subroutine initialize(this,lmP_max)

      class(nonlinear_lm_t) :: this
      integer, intent(in) :: lmP_max

      allocate( this%AdvrLM(lmP_max), this%AdvtLM(lmP_max) )
      allocate( this%AdvpLM(lmP_max), this%LFrLM(lmP_max) )
      allocate( this%LFtLM(lmP_max), this%LFpLM(lmP_max) )
      allocate( this%VxBrLM(lmP_max), this%VxBtLM(lmP_max) )
      allocate( this%VxBpLM(lmP_max))
      bytes_allocated = bytes_allocated + 12*lmP_max*SIZEOF_DEF_COMPLEX

      if ( l_anel) then
         allocate( this%ViscHeatLM(lmP_max), this%OhmLossLM(lmP_max) )
         bytes_allocated = bytes_allocated+14*lmP_max*SIZEOF_DEF_COMPLEX
      end if

      if ( l_heat ) then
         allocate(this%VSrLM(lmP_max),this%VStLM(lmP_max),this%VSpLM(lmP_max))
         bytes_allocated = bytes_allocated + 3*lmP_max*SIZEOF_DEF_COMPLEX
      end if

      if ( l_chemical_conv ) then
         allocate(this%VXirLM(lmP_max),this%VXitLM(lmP_max),this%VXipLM(lmP_max))
         bytes_allocated = bytes_allocated + 3*lmP_max*SIZEOF_DEF_COMPLEX
      end if

      !-- RMS calculations
      if ( l_RMS ) then
         allocate( this%dtVrLM(lmP_max), this%dtVtLM(lmP_max) )
         allocate( this%dtVpLM(lmP_max) )
         allocate( this%Advt2LM(lmP_max), this%Advp2LM(lmP_max) )
         allocate( this%LFt2LM(lmP_max), this%LFp2LM(lmP_max) )
         allocate( this%CFt2LM(lmP_max), this%CFp2LM(lmP_max) )
         allocate( this%PFt2LM(lmP_max), this%PFp2LM(lmP_max) )
         bytes_allocated = bytes_allocated + 11*lmP_max*SIZEOF_DEF_COMPLEX
         if ( l_adv_curl ) then
            allocate( this%dpkindrLM(lmP_max) )
            bytes_allocated = bytes_allocated + lmP_max*SIZEOF_DEF_COMPLEX
         end if
      end if

   end subroutine initialize
   
!----------------------------------------------------------------------------
!@>TODO temporary function to help transition, delete me!
   subroutine slice_all(nl_lm_glb, nl_lm_loc)
      
      class(nonlinear_lm_t), intent(inout) :: nl_lm_loc
      class(nonlinear_lm_t), intent(inout) :: nl_lm_glb
   
      call slice_FlmP_cmplx(nl_lm_glb%AdvrLM, nl_lm_loc%AdvrLM)
      call slice_FlmP_cmplx(nl_lm_glb%AdvtLM, nl_lm_loc%AdvtLM)
      call slice_FlmP_cmplx(nl_lm_glb%AdvpLM, nl_lm_loc%AdvpLM)
      call slice_FlmP_cmplx(nl_lm_glb%LFrLM , nl_lm_loc%LFrLM )
      call slice_FlmP_cmplx(nl_lm_glb%LFtLM , nl_lm_loc%LFtLM )
      call slice_FlmP_cmplx(nl_lm_glb%LFpLM , nl_lm_loc%LFpLM )
      call slice_FlmP_cmplx(nl_lm_glb%VxBrLM, nl_lm_loc%VxBrLM)
      call slice_FlmP_cmplx(nl_lm_glb%VxBtLM, nl_lm_loc%VxBtLM)
      call slice_FlmP_cmplx(nl_lm_glb%VxBpLM, nl_lm_loc%VxBpLM)

      if ( l_anel) then
         call slice_FlmP_cmplx(nl_lm_glb%ViscHeatLM, nl_lm_loc%ViscHeatLM)
         call slice_FlmP_cmplx(nl_lm_glb%OhmLossLM,  nl_lm_loc%OhmLossLM )
      end if

      if ( l_heat ) then
         call slice_FlmP_cmplx(nl_lm_glb%VSrLM, nl_lm_loc%VSrLM)
         call slice_FlmP_cmplx(nl_lm_glb%VStLM, nl_lm_loc%VStLM)
         call slice_FlmP_cmplx(nl_lm_glb%VSpLM, nl_lm_loc%VSpLM)
      end if

      if ( l_chemical_conv ) then
         call slice_FlmP_cmplx(nl_lm_glb%VXirLM, nl_lm_loc%VXirLM)
         call slice_FlmP_cmplx(nl_lm_glb%VXitLM, nl_lm_loc%VXitLM)
         call slice_FlmP_cmplx(nl_lm_glb%VXipLM, nl_lm_loc%VXipLM)
      end if

      !-- RMS calculations
      if ( l_RMS ) then
         call slice_FlmP_cmplx(nl_lm_glb%dtVrLM , nl_lm_loc%dtVrLM )
         call slice_FlmP_cmplx(nl_lm_glb%dtVtLM , nl_lm_loc%dtVtLM )
         call slice_FlmP_cmplx(nl_lm_glb%dtVpLM , nl_lm_loc%dtVpLM )
         call slice_FlmP_cmplx(nl_lm_glb%Advt2LM, nl_lm_loc%Advt2LM)
         call slice_FlmP_cmplx(nl_lm_glb%Advp2LM, nl_lm_loc%Advp2LM)
         call slice_FlmP_cmplx(nl_lm_glb%LFt2LM , nl_lm_loc%LFt2LM )
         call slice_FlmP_cmplx(nl_lm_glb%LFp2LM , nl_lm_loc%LFp2LM )
         call slice_FlmP_cmplx(nl_lm_glb%CFt2LM , nl_lm_loc%CFt2LM )
         call slice_FlmP_cmplx(nl_lm_glb%CFp2LM , nl_lm_loc%CFp2LM )
         call slice_FlmP_cmplx(nl_lm_glb%PFt2LM , nl_lm_loc%PFt2LM )
         call slice_FlmP_cmplx(nl_lm_glb%PFp2LM , nl_lm_loc%PFp2LM )
         if ( l_adv_curl ) then
            call slice_FlmP_cmplx(nl_lm_glb%dpkindrLM, nl_lm_loc%dpkindrLM)
         end if
      end if
   
   end subroutine slice_all
   
!----------------------------------------------------------------------------
!@>TODO temporary function to help transition, delete me!
   subroutine gather_all(nl_lm_loc, nl_lm_glb)
   
      class(nonlinear_lm_t), intent(inout) :: nl_lm_loc
      class(nonlinear_lm_t), intent(inout) :: nl_lm_glb
   
      call gather_FlmP(nl_lm_loc%AdvrLM, nl_lm_glb%AdvrLM)
      call gather_FlmP(nl_lm_loc%AdvtLM, nl_lm_glb%AdvtLM)
      call gather_FlmP(nl_lm_loc%AdvpLM, nl_lm_glb%AdvpLM)
! !       call gather_FlmP(nl_lm_loc%LFrLM , nl_lm_glb%LFrLM )
! !       call gather_FlmP(nl_lm_loc%LFtLM , nl_lm_glb%LFtLM )
! !       call gather_FlmP(nl_lm_loc%LFpLM , nl_lm_glb%LFpLM )
! !       call gather_FlmP(nl_lm_loc%VxBrLM, nl_lm_glb%VxBrLM)
! !       call gather_FlmP(nl_lm_loc%VxBtLM, nl_lm_glb%VxBtLM)
! !       call gather_FlmP(nl_lm_loc%VxBpLM, nl_lm_glb%VxBpLM)
! ! 
! !       if ( l_anel) then
! !          call gather_FlmP(nl_lm_loc%ViscHeatLM, nl_lm_glb%ViscHeatLM)
! !          call gather_FlmP(nl_lm_loc%OhmLossLM,  nl_lm_glb%OhmLossLM )
! !       end if
! ! 
! !       if ( l_heat ) then
! !          call gather_FlmP(nl_lm_loc%VSrLM, nl_lm_glb%VSrLM)
! !          call gather_FlmP(nl_lm_loc%VStLM, nl_lm_glb%VStLM)
! !          call gather_FlmP(nl_lm_loc%VSpLM, nl_lm_glb%VSpLM)
! !       end if
! ! 
! !       if ( l_chemical_conv ) then
! !          call gather_FlmP(nl_lm_loc%VXirLM, nl_lm_glb%VXirLM)
! !          call gather_FlmP(nl_lm_loc%VXitLM, nl_lm_glb%VXitLM)
! !          call gather_FlmP(nl_lm_loc%VXipLM, nl_lm_glb%VXipLM)
! !       end if
! ! 
! !       !-- RMS calculations
! !       if ( l_RMS ) then
! !          call gather_FlmP(nl_lm_loc%dtVrLM , nl_lm_glb%dtVrLM )
! !          call gather_FlmP(nl_lm_loc%dtVtLM , nl_lm_glb%dtVtLM )
! !          call gather_FlmP(nl_lm_loc%dtVpLM , nl_lm_glb%dtVpLM )
! !          call gather_FlmP(nl_lm_loc%Advt2LM, nl_lm_glb%Advt2LM)
! !          call gather_FlmP(nl_lm_loc%Advp2LM, nl_lm_glb%Advp2LM)
! !          call gather_FlmP(nl_lm_loc%LFt2LM , nl_lm_glb%LFt2LM )
! !          call gather_FlmP(nl_lm_loc%LFp2LM , nl_lm_glb%LFp2LM )
! !          call gather_FlmP(nl_lm_loc%CFt2LM , nl_lm_glb%CFt2LM )
! !          call gather_FlmP(nl_lm_loc%CFp2LM , nl_lm_glb%CFp2LM )
! !          call gather_FlmP(nl_lm_loc%PFt2LM , nl_lm_glb%PFt2LM )
! !          call gather_FlmP(nl_lm_loc%PFp2LM , nl_lm_glb%PFp2LM )
! !          if ( l_adv_curl ) then
! !             call gather_FlmP(nl_lm_loc%dpkindrLM, nl_lm_glb%dpkindrLM)
! !          end if
! !       end if
   
   end subroutine gather_all
   
!----------------------------------------------------------------------------
   subroutine finalize(this)

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
   subroutine set_zero(this)

      class(nonlinear_lm_t) :: this

      !-- Local variable
      integer :: lm

#ifdef WITH_SHTNS
      !$omp parallel do private(lm)
#endif
      do lm=1,lmP_max
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
#ifdef WITH_SHTNS
      !$omp end parallel do
#endif

   end subroutine set_zero
!----------------------------------------------------------------------------
   subroutine output(this)

      class(nonlinear_lm_t) :: this

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
      class(nonlinear_lm_t) :: this

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
      complex(cp) :: CorPol(lm_max), dpdr(lm_max)
      complex(cp) :: AdvPol(lm_max),AdvTor(lm_max)
      complex(cp) :: LFPol(lm_max),LFTor(lm_max)
      complex(cp) :: Geo(lm_max),CLF(lm_max),PLF(lm_max)
      complex(cp) :: ArcMag(lm_max),Mag(lm_max),CIA(lm_max),Arc(lm_max)
      complex(cp) :: Buo_temp(lm_max)
      complex(cp) :: Buo_xi(lm_max)
      complex(cp) :: AdvPol_loc,CorPol_loc,AdvTor_loc,CorTor_loc
      complex(cp) :: dsdt_loc
#ifndef WITH_SHTNS
      complex(cp) :: dxidt_loc
#endif

      integer, parameter :: DOUBLE_COMPLEX_PER_CACHELINE=4


      !write(*,"(I3,A,4ES20.12)") nR,": get_td start: ",SUM(this%AdvrLM)

      !lm_chunksize=(((lm_max)/nThreads)/DOUBLE_COMPLEX_PER_CACHELINE) * &
      !             & DOUBLE_COMPLEX_PER_CACHELINE
      !lm_chunksize=4
      !write(*,"(A,I4)") "Using a chunksize of ",lm_chunksize

      if (nBc == 0 .or. lRmsCalc ) then

         if ( l_conv ) then  ! Convection

            lm =1   ! This is l=0,m=0
            lmA=lm2lmA(lm)
            lmP=1
            lmPA=lmP2lmPA(lmP)
            if ( l_conv_nl ) then
               AdvPol_loc=      or2(nR)*this%AdvrLM(lm)
               AdvTor_loc=-dTheta1A(lm)*this%AdvpLM(lmPA)
            else
               AdvPol_loc=zero
               AdvTor_loc=zero
            end if
            if ( l_corr ) then
               CorPol_loc=two*CorFac*or1(nR) * dTheta2A(lm)* z_Rloc(lmA,nR)
               CorTor_loc= two*CorFac*or2(nR) * (                 &
               &                dTheta3A(lm)*dw_Rloc(lmA,nR) +    &
               &        or1(nR)*dTheta4A(lm)* w_Rloc(lmA,nR) )
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
                  Buo_temp(lm) =BuoFac*rgrav(nR)*rho0(nR)*s_Rloc(lm,nR)
               else
                  Buo_temp(lm) =0.0_cp
               end if
               if (l_chemical_conv) then
                  Buo_xi(lm) =ChemFac*rgrav(nR)*rho0(nR)*xi_Rloc(lm,nR)
               else
                  Buo_xi(lm)=0.0_cp
               end if
               if ( l_mag_LF .and. nR>n_r_LCR ) then
                  LFPol(lm) =      or2(nR)*this%LFrLM(lm)
                  LFTor(lm) =-dTheta1A(lm)*this%LFpLM(lmPA)
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
                  dpdr(lm)=Buo_temp(lm)+Buo_xi(lm)+beta(nR)*p_Rloc(lm,nR)+ &
                  &        AdvPol_loc+CorPol_loc
               else
                  dpdr(lm)=dp_Rloc(lm,nR)
               end if

            end if ! lRmsCalc

            !PERFON('td_cv1')
            !$omp parallel do default(shared) private(lm,l,m,lmS,lmA,lmP) &
            !$omp private(lmPS,lmPA,AdvPol_loc,CorPol_loc,AdvTor_loc,CorTor_loc)
            do lm=2,lm_max
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
                        &                    dPhi(lm)*(                          &
                        &         -ddw_Rloc(lm,nR)+beta(nR)*dw_Rloc(lm,nR)     + &
                        &             ( beta(nR)*or1(nR)+or2(nR))*               &
                        &                            dLh(lm)*w_Rloc(lm,nR) )   + &
                        &             dTheta3A(lm)*( dz_Rloc(lmA,nR)-            &
                        &                            beta(nR)*z_Rloc(lmA,nR) ) + &
                        &             dTheta3S(lm)*( dz_Rloc(lmS,nR)-            &
                        &                            beta(nR)*z_Rloc(lmS,nR) ) + &
                        &          or1(nR)* (                                    &
                        &             dTheta4A(lm)* z_Rloc(lmA,nR)               &
                        &            -dTheta4S(lm)* z_Rloc(lmS,nR) ) )
                     else if ( l == l_max ) then
                        CorPol_loc =0.0_cp
                     else if ( l == m ) then
                        CorPol_loc =two*CorFac*or2(nR)*orho1(nR)*(               &
                        &                    dPhi(lm)*(                          &
                        &         -ddw_Rloc(lm,nR)+beta(nR)*dw_Rloc(lm,nR)     + &
                        &             ( beta(nR)*or1(nR)+or2(nR))*               &
                        &                            dLh(lm)*w_Rloc(lm,nR) )   + &
                        &             dTheta3A(lm)*( dz_Rloc(lmA,nR)-            &
                        &                            beta(nR)*z_Rloc(lmA,nR) ) + &
                        &          or1(nR)* (                                    &
                        &             dTheta4A(lm)* z_Rloc(lmA,nR) ) )
                     end if
                  else
                     CorPol_loc=zero
                  end if

                  if ( l_conv_nl ) then

                     if ( l > m ) then
                        dVxVhLM(lm)=      orho1(nR)*r(nR)*r(nR)* ( &
                        &        dTheta1S(lm)*this%AdvtLM(lmPS) -  &
                        &        dTheta1A(lm)*this%AdvtLM(lmPA) +  &
                        &             dPhi(lm)*this%AdvpLM(lmP)  )
                     else if ( l == m ) then
                        dVxVhLM(lm)=      orho1(nR)*r(nR)*r(nR)* ( &
                        &      - dTheta1A(lm)*this%AdvtLM(lmPA) +  &
                        &        dPhi(lm)*this%AdvpLM(lmP)  )
                     end if

                     AdvPol_loc=dLh(lm)*or4(nR)*orho1(nR)*this%AdvrLM(lmP)

                  else

                     AdvPol_loc =zero
                     dVxVhLM(lm)=zero

                  endif

               else ! We don't use the double curl

                  if ( l_corr .and. nBc /= 2 ) then
                     if ( l < l_max .and. l > m ) then
                        CorPol_loc =two*CorFac*or1(nR) * (  &
                        &        dPhi(lm)*dw_Rloc(lm,nR) +  & ! phi-deriv of dw/dr
                        &    dTheta2A(lm)*z_Rloc(lmA,nR) -  & ! sin(theta) dtheta z
                        &    dTheta2S(lm)*z_Rloc(lmS,nR) )
                     else if ( l == l_max ) then
                        CorPol_loc=0.0_cp
                     else if ( l == m ) then
                        CorPol_loc = two*CorFac*or1(nR) * (  &
                        &         dPhi(lm)*dw_Rloc(lm,nR)  + &
                        &     dTheta2A(lm)*z_Rloc(lmA,nR) )
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
                        &       rho0(nR)*s_Rloc(lm,nR)-ViscHeatFac*      &
                        &     (ThExpNb*alpha0(nR)*temp0(nR)+ogrun(nR))*  &
                        &     p_Rloc(lm,nR) )
                     else
                        Buo_temp(lm) =0.0_cp
                     end if
                     if (l_chemical_conv) then
                        Buo_xi(lm) =ChemFac*alpha0(nR)*rgrav(nR)*(       &
                        &     rho0(nR)*xi_Rloc(lm,nR)-ViscHeatFac*       &
                        &     (ThExpNb*alpha0(nR)*temp0(nR)+ogrun(nR))*  &
                        &     p_Rloc(lm,nR) )
                     else
                        Buo_xi(lm) =0.0_cp
                     end if
                  else
                     if (l_heat) then
                        Buo_temp(lm) =BuoFac*rho0(nR)*rgrav(nR)*s_Rloc(lm,nR)
                     else
                        Buo_temp(lm) =0.0_cp
                     end if
                     if (l_chemical_conv) then
                        Buo_xi(lm) =ChemFac*rho0(nR)*rgrav(nR)*xi_Rloc(lm,nR)
                     else
                        Buo_xi(lm)=0.0_cp
                     end if
                  end if

                  if ( l_double_curl ) then
                     ! In that case we have to recompute the Coriolis force
                     ! since we also want the pressure gradient
                     if ( l_corr .and. nBc /= 2 ) then
                        if ( l < l_max .and. l > m ) then
                           CorPol_loc =two*CorFac*or1(nR) * (  &
                           &        dPhi(lm)*dw_Rloc(lm,nR) +  & ! phi-deriv of dw/dr
                           &    dTheta2A(lm)*z_Rloc(lmA,nR) -  & ! sin(theta) dtheta z
                           &    dTheta2S(lm)*z_Rloc(lmS,nR) )
                        else if ( l == l_max ) then
                           CorPol_loc= 0.0_cp
                        else if ( l == m ) then
                           CorPol_loc = two*CorFac*or1(nR) * (  &
                           &         dPhi(lm)*dw_Rloc(lm,nR)  + &
                           &     dTheta2A(lm)*z_Rloc(lmA,nR) )
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
                     dpdr(lm)=Buo_temp(lm)+Buo_xi(lm)-this%dtVrLM(lmP)+         &
                     &        dLh(lm)*or2(nR)*hdif_V(lm)*visc(nR)* (            &
                     &                                        ddw_Rloc(lm,nR)+  &
                     &         (two*dLvisc(nR)-third*beta(nR))*dw_Rloc(lm,nR)-  &
                     &         ( dLh(lm)*or2(nR)+four*third*( dbeta(nR)+        &
                     &         dLvisc(nR)*beta(nR)+(three*dLvisc(nR)+beta(nR))* &
                     &         or1(nR)) )*                      w_Rloc(lm,nR))+ &
                     &        beta(nR)*p_Rloc(lm,nR)+AdvPol_loc+CorPol_loc
                  else
                     dpdr(lm)=dp_Rloc(lm,nR)
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
                     CorTor_loc=          two*CorFac*or2(nR) * (  &
                     &                 dPhi(lm)*z_Rloc(lm,nR)   + &
                     &            dTheta3A(lm)*dw_Rloc(lmA,nR)  + &
                     &    or1(nR)*dTheta4A(lm)* w_Rloc(lmA,nR)  + &
                     &            dTheta3S(lm)*dw_Rloc(lmS,nR)  - &
                     &    or1(nR)*dTheta4S(lm)* w_Rloc(lmS,nR)  )
                  else if ( l == l_max ) then
                     CorTor_loc=0.0_cp
                  else if ( l == m ) then
                     CorTor_loc=  two*CorFac*or2(nR) * (  &
                     &         dPhi(lm)*z_Rloc(lm,nR)   + &
                     &    dTheta3A(lm)*dw_Rloc(lmA,nR)  + &
                     &    or1(nR)*dTheta4A(lm)* w_Rloc(lmA,nR)  )
                  end if
               else
                  CorTor_loc=zero
               end if

               if ( l_conv_nl ) then
                  if ( l > m ) then
                     AdvTor_loc=   -dPhi(lm)*this%AdvtLM(lmP)  + &
                     &          dTheta1S(lm)*this%AdvpLM(lmPS) - &
                     &          dTheta1A(lm)*this%AdvpLM(lmPA)
                  else if ( l == m ) then
                     AdvTor_loc=   -dPhi(lm)*this%AdvtLM(lmP)  - &
                     &          dTheta1A(lm)*this%AdvpLM(lmPA)
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
                        LFTor(lm) =   -dPhi(lm)*this%LFtLM(lmP)  + &
                        &          dTheta1S(lm)*this%LFpLM(lmPS) - &
                        &          dTheta1A(lm)*this%LFpLM(lmPA)
                     else if ( l == m ) then
                        LFTor(lm) =   -dPhi(lm)*this%LFtLM(lmP)  - &
                        &          dTheta1A(lm)*this%LFpLM(lmPA)
                     end if
                     AdvTor(lm)=AdvTor_loc-LFTor(lm)
                  else
                     AdvTor(lm)=AdvTor_loc
                  end if
               end if

            end do
            !$omp end parallel do
            !PERFOFF

            if ( lRmsCalc ) then

               if ( l_conv_nl ) then
                  if ( l_adv_curl ) then
                     do lm=1,lm_max
                        lmP =lm2lmP(lm)
                        AdvPol(lm)=AdvPol(lm)-this%dpkindrLM(lmP)
                     end do
                  end if
                  call hIntRms(AdvPol,nR,1,lm_max,0,Adv2hInt(:,nR),st_map, .false.)
                  call hIntRms(this%Advt2LM,nR,1,lmP_max,1,Adv2hInt(:,nR),st_map, &
                       &       .true.)
                  call hIntRms(this%Advp2LM,nR,1,lmP_max,1,Adv2hInt(:,nR),st_map, &
                       &       .true.)
                  do lm=1,lm_max
                     lmP =lm2lmP(lm)
                     !-- Use Geo as work array
                     Geo(lm)=AdvPol(lm)-this%dtVrLM(lmP)
                  end do
                  call hIntRms(Geo,nR,1,lm_max,0,Iner2hInt(:,nR),st_map, .false.)
                  call hIntRms(this%Advt2LM-this%dtVtLM,nR,1,lmP_max,1, &
                       &       Iner2hInt(:,nR),st_map,.true.)
                  call hIntRms(this%Advp2LM-this%dtVpLM,nR,1,lmP_max,1, &
                       &       Iner2hInt(:,nR),st_map,.true.)
               end if

               if ( l_anelastic_liquid ) then
                  call hIntRms(dpdr,nR,1,lm_max,0,Pre2hInt(:,nR),st_map,.false.)
               else
                  ! rho* grad(p/rho) = grad(p) - beta*p
                  !-- Geo is used to store the pressure Gradient
                  if ( l_adv_curl ) then
                     do lm=1,lm_max
                        lmP =lm2lmP(lm)
                        !-- Use Geo as work array
                        Geo(lm)=dpdr(lm)-beta(nR)*p_Rloc(lm,nR)- &
                        &       this%dpkindrLM(lmP)
                     end do
                  else
                     do lm=1,lm_max
                        !-- Use Geo as work array
                        Geo(lm)=dpdr(lm)-beta(nR)*p_Rloc(lm,nR)
                     end do
                  end if
                  call hIntRms(Geo,nR,1,lm_max,0,Pre2hInt(:,nR),st_map,.false.)
               end if
               call hIntRms(this%PFt2LM,nR,1,lmP_max,1,Pre2hInt(:,nR),st_map,.true.)
               call hIntRms(this%PFp2LM,nR,1,lmP_max,1,Pre2hInt(:,nR),st_map,.true.)

               if ( l_heat ) then
                  call hIntRms(Buo_temp,nR,1,lm_max,0,Buo_temp2hInt(:,nR),&
                       &       st_map,.false.)
               end if
               if ( l_chemical_conv ) then
                  call hIntRms(Buo_xi,nR,1,lm_max,0,Buo_xi2hInt(:,nR),&
                       &       st_map,.false.)
               end if
               if ( l_corr ) then
                  call hIntRms(CorPol,nR,1,lm_max,0,Cor2hInt(:,nR),st_map,.false.)
                  call hIntRms(this%CFt2LM,nR,1,lmP_max,1,Cor2hInt(:,nR), &
                       &       st_map,.true.)
                  calL hIntRms(this%CFp2LM,nR,1,lmP_max,1,Cor2hInt(:,nR), &
                       &       st_map,.true.)
               end if
               if ( l_mag_LF .and. nR>n_r_LCR ) then
                  call hIntRms(LFPol,nR,1,lm_max,0,LF2hInt(:,nR),st_map,.false.)
                  call hIntRms(this%LFt2LM,nR,1,lmP_max,1,LF2hInt(:,nR),st_map,.true.)
                  call hIntRms(this%LFp2LM,nR,1,lmP_max,1,LF2hInt(:,nR),st_map,.true.)
               end if

               do lm=1,lm_max
                  lmP =lm2lmP(lm)
                  Geo(lm)=CorPol(lm)-dpdr(lm)+beta(nR)*p_Rloc(lm,nR)
                  PLF(lm)=LFPol(lm)-dpdr(lm)+beta(nR)*p_Rloc(lm,nR)
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
               call hIntRms(Geo,nR,1,lm_max,0,Geo2hInt(:,nR),st_map,.false.)
               call hIntRms(CLF,nR,1,lm_max,0,CLF2hInt(:,nR),st_map,.false.)
               call hIntRms(PLF,nR,1,lm_max,0,PLF2hInt(:,nR),st_map,.false.)
               call hIntRms(Mag,nR,1,lm_max,0,Mag2hInt(:,nR),st_map,.false.)
               call hIntRms(Arc,nR,1,lm_max,0,Arc2hInt(:,nR),st_map,.false.)
               call hIntRms(ArcMag,nR,1,lm_max,0,ArcMag2hInt(:,nR),st_map,.false.)
               call hIntRms(CIA,nR,1,lm_max,0,CIA2hInt(:,nR),st_map,.false.)

               do lm=1,lm_max
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
               call hIntRms(Geo,nR,1,lm_max,0,Geo2hInt(:,nR),st_map,.true.)
               call hIntRms(CLF,nR,1,lm_max,0,CLF2hInt(:,nR),st_map,.true.)
               call hIntRms(PLF,nR,1,lm_max,0,PLF2hInt(:,nR),st_map,.true.)
               call hIntRms(Mag,nR,1,lm_max,0,Mag2hInt(:,nR),st_map,.true.)
               call hIntRms(Arc,nR,1,lm_max,0,Arc2hInt(:,nR),st_map,.true.)
               call hIntRms(ArcMag,nR,1,lm_max,0,ArcMag2hInt(:,nR),st_map,.true.)
               call hIntRms(CIA,nR,1,lm_max,0,CIA2hInt(:,nR),st_map,.true.)

               do lm=1,lm_max
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
               call hIntRms(Geo,nR,1,lm_max,0,Geo2hInt(:,nR),st_map,.true.)
               call hIntRms(CLF,nR,1,lm_max,0,CLF2hInt(:,nR),st_map,.true.)
               call hIntRms(PLF,nR,1,lm_max,0,PLF2hInt(:,nR),st_map,.true.)
               call hIntRms(Mag,nR,1,lm_max,0,Mag2hInt(:,nR),st_map,.true.)
               call hIntRms(Arc,nR,1,lm_max,0,Arc2hInt(:,nR),st_map,.true.)
               call hIntRms(ArcMag,nR,1,lm_max,0,ArcMag2hInt(:,nR),st_map,.true.)
               call hIntRms(CIA,nR,1,lm_max,0,CIA2hInt(:,nR),st_map,.true.)

            end if

            ! In case double curl is calculated dpdt is useless
            if ( (.not. l_double_curl) .or. lPressNext ) then
            !if ( .true. ) then
               !PERFON('td_cv2')
               !$omp parallel do default(shared) private(lm,l,m,lmS,lmA,lmP) &
               !$omp private(lmPS,AdvPol_loc,CorPol_loc)
               !LIKWID_ON('td_cv2')
               do lm=2,lm_max
                  l   =lm2l(lm)
                  m   =lm2m(lm)
                  lmS =lm2lmS(lm)
                  lmA =lm2lmA(lm)
                  lmP =lm2lmP(lm)
                  lmPS=lmP2lmPS(lmP)
                  lmPA=lmP2lmPA(lmP)

                  !------ Recycle CorPol and AdvPol:
                  if ( l_corr ) then
                     !PERFON('td_cv2c')
                     if ( l < l_max .and. l > m ) then
                        CorPol_loc=           two*CorFac*or2(nR) *  &
                        &           ( -dPhi(lm)  * ( dw_Rloc(lm,nR) &
                        &            +or1(nR)*dLh(lm)*w_Rloc(lm,nR) &
                        &                                         ) &
                        &              +dTheta3A(lm)*z_Rloc(lmA,nR) &
                        &              +dTheta3S(lm)*z_Rloc(lmS,nR) &
                        &           )

                     else if ( l == l_max ) then
                        CorPol_loc=0.0_cp

                     else if ( l == m ) then
                        CorPol_loc=                    two*CorFac*or2(nR) *  &
                        &                    ( -dPhi(lm)  * ( dw_Rloc(lm,nR) &
                        &                     +or1(nR)*dLh(lm)*w_Rloc(lm,nR) &
                        &                                                   )&
                        &                      +dTheta3A(lm)*z_Rloc(lmA,nR)  &
                        &                    )

                     end if
                     !PERFOFF
                  else
                     CorPol_loc=zero
                  end if
                  if ( l_conv_nl ) then
                     !PERFON('td_cv2nl')
                     if ( l > m ) then
                        AdvPol_loc= dTheta1S(lm)*this%AdvtLM(lmPS) - &
                        &           dTheta1A(lm)*this%AdvtLM(lmPA) + &
                        &               dPhi(lm)*this%AdvpLM(lmP)
                     else if ( l == m ) then
                        AdvPol_loc=-dTheta1A(lm)*this%AdvtLM(lmPA) + &
                        &               dPhi(lm)*this%AdvpLM(lmP)
                     end if
                     !PERFOFF
                  else
                     AdvPol_loc=zero
                  end if
                  dpdt(lm)=AdvPol_loc+CorPol_loc

               end do ! lm loop
               !LIKWID_OFF('td_cv2')
               !$omp end parallel do
               !PERFOFF
            end if

         else
            do lm=2,lm_max
               dwdt(lm) =0.0_cp
               dzdt(lm) =0.0_cp
               dpdt(lm) =0.0_cp
            end do
         end if ! l_conv ?

      end if

      if ( nBc == 0 ) then

         if ( l_heat ) then
            dsdt_loc  =epsc*epscProf(nR)!+opr/epsS*divKtemp0(nR)
            dVSrLM(1)=this%VSrLM(1)
            if ( l_anel ) then
               if ( l_anelastic_liquid ) then
                  if ( l_mag_nl ) then
                     dsdt_loc=dsdt_loc+ViscHeatFac*temp0(nR)*this%ViscHeatLM(1)+ &
                     &        OhmLossFac*temp0(nR)*this%OhmLossLM(1)
                  else
                     dsdt_loc=dsdt_loc+ViscHeatFac*temp0(nR)*this%ViscHeatLM(1)
                  end if
               else
                  if ( l_mag_nl ) then
                     dsdt_loc=dsdt_loc+ViscHeatFac*this%ViscHeatLM(1)+ &
                     &        OhmLossFac*this%OhmLossLM(1)
                  else
                     dsdt_loc=dsdt_loc+ViscHeatFac*this%ViscHeatLM(1)
                  end if
               end if
            end if
            dsdt(1)=dsdt_loc

#ifdef WITH_SHTNS
            !$omp parallel do default(shared) private(lm,lmP,dsdt_loc,l,m)
            do lm=2,lm_max
               l   =lm2l(lm)
               m   =lm2m(lm)
               lmP =lm2lmP(lm)
               dVSrLM(lm)=this%VSrLM(lmP)
               dsdt_loc = dLh(lm)*this%VStLM(lmP)
               if ( l_anel ) then
                  if ( l_anelastic_liquid ) then
                     dsdt_loc = dsdt_loc+ViscHeatFac*hdif_V(lo_map%lm2(l,m))* &
                     &          temp0(nR)*this%ViscHeatLM(lmP)
                     if ( l_mag_nl ) then
                        dsdt_loc = dsdt_loc+OhmLossFac*hdif_B(lo_map%lm2(l,m)) * &
                        &          temp0(nR)*this%OhmLossLM(lmP)
                     end if
                  else
                     dsdt_loc = dsdt_loc+ViscHeatFac*hdif_V(lo_map%lm2(l,m)) * &
                     &          this%ViscHeatLM(lmP)
                     if ( l_mag_nl ) then
                        dsdt_loc = dsdt_loc+OhmLossFac*hdif_B(lo_map%lm2(l,m)) * &
                        &          this%OhmLossLM(lmP)
                     end if
                  end if
               end if
               dsdt(lm) = dsdt_loc
            end do
            !$omp end parallel do
#else
            !PERFON('td_heat')
            !$omp parallel do default(shared) private(lm,l,m,lmP,lmPS) &
            !$omp private(lmPA,dsdt_loc)
            !LIKWID_ON('td_heat')
            do lm=2,lm_max
               l   =lm2l(lm)
               m   =lm2m(lm)
               lmP =lm2lmP(lm)
               lmPS=lmP2lmPS(lmP)
               lmPA=lmP2lmPA(lmP)
               !------ This is horizontal heat advection:
               !PERFON('td_h1')

               if ( l > m ) then
                  dsdt_loc= -dTheta1S(lm)*this%VStLM(lmPS) &
                  &         +dTheta1A(lm)*this%VStLM(lmPA) &
                  &         -dPhi(lm)*this%VSpLM(lmP)
               else if ( l == m ) then
                  dsdt_loc=  dTheta1A(lm)*this%VStLM(lmPA) &
                  &          -dPhi(lm)*this%VSpLM(lmP)
               end if
               !PERFOFF
               !PERFON('td_h2')
               if ( l_anel ) then
                  if ( l_anelastic_liquid ) then
                     dsdt_loc = dsdt_loc+ViscHeatFac*hdif_V(lo_map%lm2(l,m))* &
                     &          temp0(nR)*this%ViscHeatLM(lmP)
                     if ( l_mag_nl ) then
                        dsdt_loc = dsdt_loc+OhmLossFac*hdif_B(lo_map%lm2(l,m))*&
                        &          temp0(nR)*this%OhmLossLM(lmP)
                     end if
                  else
                     dsdt_loc = dsdt_loc+ViscHeatFac*hdif_V(lo_map%lm2(l,m))* &
                     &          this%ViscHeatLM(lmP)
                     if ( l_mag_nl ) then
                        dsdt_loc = dsdt_loc+OhmLossFac*hdif_B(lo_map%lm2(l,m))* &
                        &          this%OhmLossLM(lmP)
                     end if
                  end if
               end if
               !PERFOFF
               !-----   simplified form for linear onset !
               !        not ds not saved in the current program form!
               !                 dsdt(lm)=
               !                    -dLh(lm)*w(lm,nR)*or2(nR)*dsR(1)
               dVSrLM(lm)=this%VSrLM(lmP)
               dsdt(lm) = dsdt_loc
            end do
            !LIKWID_OFF('td_heat')
            !$omp end parallel do
            !PERFOFF
#endif
         else
            do lm=2,lm_max
               dsdt(lm)  =0.0_cp
               dVSrLM(lm)=0.0_cp
            end do
         end if

         if ( l_chemical_conv ) then
            dVXirLM(1)=this%VXirLM(1)
            dxidt(1)  =epscXi

#ifdef WITH_SHTNS
            !$omp parallel do default(shared) private(lm,lmP)
            do lm=2,lm_max
               lmP=lm2lmP(lm)
               dVXirLM(lm)=this%VXirLM(lmP)
               dxidt(lm)  =dLh(lm)*this%VXitLM(lmP)
            end do
            !$omp end parallel do
#else
            !$omp parallel do default(shared) private(lm,l,m,lmP,lmPS) &
            !$omp private(lmPA,dxidt_loc)
            do lm=2,lm_max
               l   =lm2l(lm)
               m   =lm2m(lm)
               lmP =lm2lmP(lm)
               lmPS=lmP2lmPS(lmP)
               lmPA=lmP2lmPA(lmP)

               if ( l > m ) then
                  dxidt_loc= -dTheta1S(lm)*this%VXitLM(lmPS) &
                  &          +dTheta1A(lm)*this%VXitLM(lmPA) &
                  &          -dPhi(lm)*this%VXipLM(lmP)
               else if ( l == m ) then
                  dxidt_loc=  dTheta1A(lm)*this%VXitLM(lmPA) &
                  &          -dPhi(lm)*this%VXipLM(lmP)
               end if
               dVXirLM(lm)=this%VXirLM(lmP)
               dxidt(lm) = dxidt_loc
            end do
            !$omp end parallel do
#endif
         end if

         if ( l_mag_nl .or. l_mag_kin  ) then
            !PERFON('td_magnl')

#ifdef WITH_SHTNS
            !$omp parallel do default(shared) private(lm,lmP)
            do lm=1,lm_max
               lmP =lm2lmP(lm)
               dbdt(lm)   = dLh(lm)*this%VxBpLM(lmP)
               dVxBhLM(lm)=-dLh(lm)*this%VxBtLM(lmP)*r(nR)*r(nR)
               djdt(lm)   = dLh(lm)*or4(nR)*this%VxBrLM(lmP)
            end do
            !$omp end parallel do
#else
            !$omp parallel do default(shared) private(lm,l,m,lmP,lmPS,lmPA)
            do lm=1,lm_max
               if (lm == 1) then
                  lmP=1
                  lmPA=lmP2lmPA(lmP)
                  dVxBhLM(lm)= -r(nR)*r(nR)* dTheta1A(lm)*this%VxBtLM(lmPA)
                  dbdt(lm)   = -dTheta1A(lm)*this%VxBpLM(lmPA)
                  djdt(lm)   = zero
                  cycle
               end if
               l   =lm2l(lm)
               m   =lm2m(lm)
               lmP =lm2lmP(lm)
               lmPS=lmP2lmPS(lmP)
               lmPA=lmP2lmPA(lmP)

               !------- This is the radial part of the dynamo terms \curl(VxB)
               !PERFON('td_mnl1')
               if ( l > m ) then
                  dbdt(lm)=  dTheta1S(lm)*this%VxBpLM(lmPS) &
                  &         -dTheta1A(lm)*this%VxBpLM(lmPA) &
                  &         -dPhi(lm)    *this%VxBtLM(lmP)
               else if ( l == m ) then
                  dbdt(lm)= -dTheta1A(lm)*this%VxBpLM(lmPA) &
                  &         -dPhi(lm)    *this%VxBtLM(lmP)
               end if
               !PERFOFF

               !------- Radial component of
               !           \curl\curl(UxB) = \grad\div(UxB) - \laplace(VxB)

               !------- This is the radial part of \laplace (UxB)
               djdt(lm)=dLh(lm)*or4(nR)*this%VxBrLM(lmP)

               !------- This is r^2 * horizontal divergence of (UxB)
               !        Radial derivative performed in get_dr_td
               !PERFON('td_mnl2')
               if ( l > m ) then
                  dVxBhLM(lm)=            r(nR)*r(nR)* ( &
                  &    dTheta1S(lm)*this%VxBtLM(lmPS) -  &
                  &    dTheta1A(lm)*this%VxBtLM(lmPA) +  &
                  &    dPhi(lm)*this%VxBpLM(lmP)  )
               else if ( l == m ) then
                  dVxBhLM(lm)=              r(nR)*r(nR)* ( &
                  &    - dTheta1A(lm)*this%VxBtLM(lmPA) +  &
                  &    dPhi(lm)*this%VxBpLM(lmP)  )
               end if
               !PERFOFF
            end do
            !$omp end parallel do
#endif
            !PERFOFF
         else
            if ( l_mag ) then
               do lm=1,lm_max
                  dbdt(lm)   =zero
                  djdt(lm)   =zero
                  dVxBhLM(lm)=zero
               end do
            end if
         end if

      else   ! boundary !
         !PERFON('td_bnd')
         if ( l_mag_nl .or. l_mag_kin ) then

#ifdef WITH_SHTNS
            dVxBhLM(1)=zero
            dVSrLM(1) =zero
            !$omp parallel do default(shared) private(lm,lmP)
            do lm=2,lm_max
               lmP =lm2lmP(lm)
               dVxBhLM(lm)=-dLh(lm)*this%VxBtLM(lmP)*r(nR)*r(nR)
               dVSrLM(lm) =zero
            end do
            !$omp end parallel do
#else
            !----- Stress free boundary, only nl mag. term for poloidal field needed.
            !      Because the radial derivative will be taken, this will contribute to
            !      the other radial grid points.
            dVxBhLM(1)=zero
            dVSrLM(1) =zero
            !$omp parallel do default(shared) private(lm,lmP,l,m,lmPS,lmPA)
            do lm=2,lm_max
               l   =lm2l(lm)
               m   =lm2m(lm)
               lmP =lm2lmP(lm)
               lmPS=lmP2lmPS(lmP)   ! l-1
               lmPA=lmP2lmPA(lmP)   ! l+1
               if ( l > m ) then
                  dVxBhLM(lm)=r(nR)*r(nR)* (               &
                  &      dTheta1S(lm)*this%VxBtLM(lmPS) -  &
                  &      dTheta1A(lm)*this%VxBtLM(lmPA) +  &
                  &          dPhi(lm)*this%VxBpLM(lmP)  )
               else if ( l == m ) then ! (l-1) not allowed !
                  dVxBhLM(lm)=r(nR)*r(nR)* (               &
                  &    - dTheta1A(lm)*this%VxBtLM(lmPA) +  &
                  &    dPhi(lm)*this%VxBpLM(lmP)  )
               end if
               dVSrLM(lm)=zero
            end do
            !$omp end parallel do
#endif

         else
            do lm=1,lm_max
               if ( l_mag ) dVxBhLM(lm)=zero
               dVSrLM(lm) =zero
            end do
         end if
         if ( l_double_curl ) then
            do lm=1,lm_max
               dVxVhLM(lm)=zero
            end do
         end if
         if ( l_chemical_conv ) then
            do lm=1,lm_max
               dVXirLM(lm)=zero
            end do
         end if
      end if  ! boundary ? lvelo ?

   end subroutine get_td
!-----------------------------------------------------------------------------
end module nonlinear_lm_mod
