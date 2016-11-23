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
       &             l_chemical_conv, l_TP_form, l_single_matrix
   use radial_functions, only: r, or2, or1, beta, rho0, rgrav, epscProf, &
       &                       or4, temp0, alpha0
   use physical_parameters, only: CorFac, ra, epsc, ViscHeatFac, &
       &                          OhmLossFac, n_r_LCR, epscXi,   &
       &                          BuoFac, ThExpNb, ogrun
   use blocking, only: lm2l, lm2m, lm2lmP, lmP2lmPS, lmP2lmPA, lm2lmA, &
       &               lm2lmS, st_map
   use horizontal_data, only: dLh, dTheta1S, dTheta1A, dPhi, dTheta2A, &
       &                      dTheta3A, dTheta4A, dPhi0, dTheta2S,     &
       &                      dTheta3S, dTheta4S, hdif_V, hdif_B
   use RMS, only: Adv2hInt, Pre2hInt, Buo2hInt, Cor2hInt, LF2hInt,  &
       &          Geo2hInt, Mag2hInt, Arc2hInt, CLF2hInt, PLF2hInt, &
       &          CIA2hInt
   use leg_helper_mod, only: leg_helper_t
   use constants, only: zero, two
   use fields, only: w_Rloc,dw_Rloc,z_Rloc
   use RMS_helpers, only: hIntRms
    

   implicit none
 
   type :: nonlinear_lm_t
      !----- Nonlinear terms in lm-space: 
      complex(cp), allocatable :: AdvrLM(:), AdvtLM(:), AdvpLM(:)
      complex(cp), allocatable :: LFrLM(:),  LFtLM(:),  LFpLM(:)
      complex(cp), allocatable :: VxBrLM(:), VxBtLM(:), VxBpLM(:)
      complex(cp), allocatable :: VSrLM(:),  VStLM(:),  VSpLM(:)
      complex(cp), allocatable :: VXirLM(:),  VXitLM(:),  VXipLM(:)
      complex(cp), allocatable :: VPrLM(:)
      complex(cp), allocatable :: ViscHeatLM(:), OhmLossLM(:)
      !----- RMS calculations
      complex(cp), allocatable :: Advt2LM(:), Advp2LM(:)
      complex(cp), allocatable :: LFt2LM(:), LFp2LM(:)
      complex(cp), allocatable :: CFt2LM(:), CFp2LM(:)
      complex(cp), allocatable :: p1LM(:), p2LM(:)
 
   contains
 
      procedure :: initialize
      procedure :: finalize
      procedure :: output
      procedure :: set_zero
      procedure :: get_td
 
   end type nonlinear_lm_t

contains

   subroutine initialize(this,lmP_max)

      class(nonlinear_lm_t) :: this
      integer, intent(in) :: lmP_max

      allocate( this%AdvrLM(lmP_max) )   
      allocate( this%AdvtLM(lmP_max) )   
      allocate( this%AdvpLM(lmP_max) )   
      allocate( this%LFrLM(lmP_max) )    
      allocate( this%LFtLM(lmP_max) )    
      allocate( this%LFpLM(lmP_max) )    
      allocate( this%VxBrLM(lmP_max) )   
      allocate( this%VxBtLM(lmP_max) )   
      allocate( this%VxBpLM(lmP_max) )   
      allocate( this%VSrLM(lmP_max) )    
      allocate( this%VStLM(lmP_max) )    
      allocate( this%VSpLM(lmP_max) )    
      allocate( this%ViscHeatLM(lmP_max) )
      allocate( this%OhmLossLM(lmP_max) )
      bytes_allocated = bytes_allocated + 14*lmP_max*SIZEOF_DEF_COMPLEX

      if ( l_TP_form ) then
         allocate( this%VPrLM(lmP_max) )    
         bytes_allocated = bytes_allocated + lmP_max*SIZEOF_DEF_COMPLEX
      end if
      
      if ( l_chemical_conv ) then
         allocate( this%VXirLM(lmP_max) )    
         allocate( this%VXitLM(lmP_max) )    
         allocate( this%VXipLM(lmP_max) )    
         bytes_allocated = bytes_allocated + 3*lmP_max*SIZEOF_DEF_COMPLEX
      end if

      !-- RMS calculations
      if ( l_RMS ) then
         allocate( this%Advt2LM(lmP_max) )
         allocate( this%Advp2LM(lmP_max) )
         allocate( this%LFt2LM(lmP_max) )
         allocate( this%LFp2LM(lmP_max) )
         allocate( this%CFt2LM(lmP_max) )
         allocate( this%CFp2LM(lmP_max) )
         allocate( this%p1LM(lmP_max) )
         allocate( this%p2LM(lmP_max) )
         bytes_allocated = bytes_allocated + 8*lmP_max*SIZEOF_DEF_COMPLEX
      end if

   end subroutine initialize
!----------------------------------------------------------------------------
   subroutine finalize(this)

      class(nonlinear_lm_t) :: this

      deallocate( this%AdvrLM )   
      deallocate( this%AdvtLM )   
      deallocate( this%AdvpLM )   
      deallocate( this%LFrLM )    
      deallocate( this%LFtLM )    
      deallocate( this%LFpLM )    
      deallocate( this%VxBrLM )   
      deallocate( this%VxBtLM )   
      deallocate( this%VxBpLM )   
      deallocate( this%VSrLM )    
      deallocate( this%VStLM )    
      deallocate( this%VSpLM )    
      deallocate( this%ViscHeatLM )
      deallocate( this%OhmLossLM )

      if ( l_TP_form ) then
         deallocate( this%VPrLM )    
      end if

      if ( l_chemical_conv ) then
         deallocate( this%VXirLM )    
         deallocate( this%VXitLM )    
         deallocate( this%VXipLM )    
      end if

      !-- RMS calculations
      if ( l_RMS ) then
         deallocate( this%Advt2LM )
         deallocate( this%Advp2LM )
         deallocate( this%LFt2LM )
         deallocate( this%LFp2LM )
         deallocate( this%CFt2LM )
         deallocate( this%CFp2LM )
         deallocate( this%p1LM )
         deallocate( this%p2LM )
      end if

   end subroutine finalize
!----------------------------------------------------------------------------
   subroutine set_zero(this)

      class(nonlinear_lm_t) :: this
      
      this%AdvrLM    =zero
      this%AdvtLM    =zero
      this%AdvpLM    =zero
      this%LFrLM     =zero
      this%LFtLM     =zero
      this%LFpLM     =zero
      this%VxBrLM    =zero
      this%VxBtLM    =zero
      this%VxBpLM    =zero
      this%VSrLM     =zero
      this%VStLM     =zero
      this%VSpLM     =zero
      this%ViscHeatLM=zero
      this%OhmLossLM =zero

      if ( l_TP_form ) then
         this%VPrLM =zero
      end if

      if ( l_chemical_conv ) then
         this%VXirLM =zero
         this%VXitLM =zero
         this%VXipLM =zero
      end if

      if ( l_RMS ) then
         this%Advt2LM=zero
         this%Advp2LM=zero
         this%LFp2LM =zero
         this%LFt2LM =zero
         this%CFt2LM =zero
         this%CFp2LM =zero
         this%p1LM   =zero
         this%p2LM   =zero
      end if

   end subroutine set_zero
!----------------------------------------------------------------------------
   subroutine output(this)

      class(nonlinear_lm_t) :: this
      
      write(*,"(A,6ES20.12)") "AdvrLM,AdvtLM,AdvpLM = ",&
           & sum(this%AdvrLM), sum(this%AdvtLM),        &
           & sum(this%AdvpLM)

   end subroutine output
!----------------------------------------------------------------------------
   subroutine get_td(this,nR,nBc,lRmsCalc,dVSrLM,dVPrLM,dVXirLM,dVxBhLM, &
              &      dwdt,dzdt,dpdt,dsdt,dxidt,dbdt,djdt,leg_helper)
      !
      !  Purpose of this to calculate time derivatives
      !  dwdt,dzdt,dpdt,dsdt,dxidt,dbdt,djdt
      !  and auxiliary arrays dVSrLM, dVXirLM and dVxBhLM
      !  from non-linear terms in spectral form,
      !  contained in flmw1-3,flms1-3, flmb1-3 (input)
      !
    
      !-- Input of variables:
      class(nonlinear_lm_t) :: this

      integer,            intent(in) :: nR
      integer,            intent(in) :: nBc ! signifies boundary conditions
      logical,            intent(in) :: lRmsCalc
      type(leg_helper_t), intent(in) :: leg_helper
    
      !-- Output of variables:
      complex(cp), intent(out) :: dwdt(:),dzdt(:)
      complex(cp), intent(out) :: dpdt(:),dsdt(:)
      complex(cp), intent(out) :: dxidt(:)
      complex(cp), intent(out) :: dbdt(:),djdt(:)
      complex(cp), intent(out) :: dVxBhLM(:)
      complex(cp), intent(out) :: dVSrLM(:)
      complex(cp), intent(out) :: dVXirLM(:)
      complex(cp), intent(out) :: dVPrLM(:)
    
      !-- Local variables:
      integer :: l,m,lm,lmS,lmA,lmP,lmPS,lmPA
      complex(cp) :: CorPol(lm_max)
      complex(cp) :: AdvPol(lm_max),AdvTor(lm_max)
      complex(cp) :: LFPol(lm_max),LFTor(lm_max)
      complex(cp) :: Geo(lm_max),CLF(lm_max),PLF(lm_max)
      complex(cp) :: Arc(lm_max),Mag(lm_max),CIA(lm_max)
      complex(cp) :: dpt(lm_max),dpp(lm_max), Buo(lm_max)
      complex(cp) :: AdvPol_loc,CorPol_loc,AdvTor_loc,CorTor_loc
      complex(cp) :: dsdt_loc, dxidt_loc
    
      integer, parameter :: DOUBLE_COMPLEX_PER_CACHELINE=4
    
    
      !write(*,"(I3,A,4ES20.12)") nR,": get_td start: ",SUM(this%AdvrLM), &
      !                                                 SUM(leg_helper%dLHz)
    
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
               CorTor_loc= two*CorFac*or2(nR) * ( &
                    dTheta3A(lm)*dw_Rloc(lmA,nR) + &
                    or1(nR)*dTheta4A(lm)* w_Rloc(lmA,nR) )
            else
               CorPol_loc=zero
               CorTor_loc=zero
            end if
            Buo(lm) =BuoFac*rgrav(nR)*rho0(nR)*leg_helper%sR(lm)
            if ( l_single_matrix ) then
               dwdt(lm)=AdvPol_loc!+CorPol_loc
            else
               dwdt(lm)=AdvPol_loc+CorPol_loc
            end if
            dzdt(lm)=AdvTor_loc+CorTor_loc
            if ( lRmsCalc ) then
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
            end if
    
            !PERFON('td_cv1')
            !$OMP PARALLEL do default(none) &
            !$OMP private(lm,l,m,lmS,lmA,lmP,lmPS,lmPA) &
            !$OMP private(AdvPol_loc,CorPol_loc,AdvTor_loc,CorTor_loc) &
            !$OMP shared(lm2l,lm2m,lm2lmS,lm2lmA,lm2lmP,lmP2lmPS,lmP2lmPA) &
            !$OMP shared(lm_max,l_corr,l_max,l_conv_nl,lRmsCalc,l_mag_LF) &
            !$OMP shared(CorPol,AdvPol,LFPol,AdvTor,LFTor,z_Rloc,nR) &
            !$OMP shared(w_Rloc,this,dw_Rloc,nBc,leg_helper,Buo) &
            !$OMP shared(CorFac,or1,or2,dPhi0,dPhi,dTheta2A,dTheta2S,n_r_LCR) &
            !$OMP shared(dTheta3A,dTheta4A,dTheta3S,dTheta4S,dTheta1S,dTheta1A) &
            !$OMP shared(dwdt,dzdt,rho0,rgrav,BuoFac,l_TP_form,temp0) &
            !$OMP shared(l_anelastic_liquid,alpha0,ThExpNb,ViscHeatFac,ogrun)
            do lm=2,lm_max
               l   =lm2l(lm)
               m   =lm2m(lm)
               lmS =lm2lmS(lm)
               lmA =lm2lmA(lm)
               lmP =lm2lmP(lm)
               lmPS=lmP2lmPS(lmP)
               lmPA=lmP2lmPA(lmP)
    
#if 0
               write(*,"(A,I4,A)") "========== lm = ",lm," =========="
               call print_cache_info_integer("l"//C_NULL_CHAR,l)
               call print_cache_info_dcmplx("dPhi0"//C_NULL_CHAR,dPhi0(lm))
               call print_cache_info_dcmplx("dPhi"//C_NULL_CHAR,dPhi(lm))
               call print_cache_info_dreal("dTheta1A"//C_NULL_CHAR,dTheta1A(lm))
               call print_cache_info_dreal("dTheta1S"//C_NULL_CHAR,dTheta1S(lm))
               call print_cache_info_dreal("dTheta2A"//C_NULL_CHAR,dTheta2A(lm))
               call print_cache_info_dreal("dTheta2S"//C_NULL_CHAR,dTheta2S(lm))
               call print_cache_info_dreal("dTheta3A"//C_NULL_CHAR,dTheta3A(lm))
               call print_cache_info_dreal("dTheta3S"//C_NULL_CHAR,dTheta3S(lm))
               call print_cache_info_dreal("dTheta4A"//C_NULL_CHAR,dTheta4A(lm))
               call print_cache_info_dreal("dTheta4S"//C_NULL_CHAR,dTheta4S(lm))
               call print_cache_info_dcmplx("z_Rloc(lmA)"//C_NULL_CHAR,z_Rloc(lmA,nR))
               call print_cache_info_dcmplx("z_Rloc(lmS)"//C_NULL_CHAR,z_Rloc(lmS,nR))
               call print_cache_info_dcmplx("w_Rloc(lmA)"//C_NULL_CHAR,w_Rloc(lmA,nR))
               call print_cache_info_dcmplx("w_Rloc(lmS)"//C_NULL_CHAR,w_Rloc(lmS,nR))
               call print_cache_info_dcmplx("dw_Rloc(lmA)"//C_NULL_CHAR,dw_Rloc(lmA,nR))
               call print_cache_info_dcmplx("dw_Rloc(lmS)"//C_NULL_CHAR,dw_Rloc(lmS,nR))
               call print_cache_info_dcmplx("dw_Rloc"//C_NULL_CHAR,dw_Rloc(lm,nR))
               call print_cache_info_dcmplx("this%AdvrLM(lmP)"//C_NULL_CHAR, &
                                            this%AdvrLM(lmP))
               call print_cache_info_dcmplx("this%AdvtLM(lmP)"//C_NULL_CHAR, &
                                            this%AdvtLM(lmP))
               call print_cache_info_dcmplx("this%AdvpLM(lmPA)"//C_NULL_CHAR, &
                                            this%AdvpLM(lmPA))
               call print_cache_info_dcmplx("this%AdvpLM(lmPS)"//C_NULL_CHAR, &
                                            this%AdvpLM(lmPS))
               write(*,"(A)") ""
#endif
               if ( l_corr .and. nBc /= 2 ) then
                  if ( l < l_max .and. l > m ) then
                     CorPol_loc =two*CorFac*or1(nR) * (  &
                     &       dPhi0(lm)*dw_Rloc(lm,nR) +  & ! phi-deriv of dw/dr
                     &    dTheta2A(lm)*z_Rloc(lmA,nR) -  & ! sin(theta) dtheta z
                     &    dTheta2S(lm)*z_Rloc(lmS,nR) )
                  else if ( l == l_max ) then
                     CorPol_loc= two*CorFac*or1(nR) * ( &
                     &            dPhi0(lm)*dw_Rloc(lm,nR)  )
                  else if ( l == m ) then
                     CorPol_loc = two*CorFac*or1(nR) * ( &
                     &        dPhi0(lm)*dw_Rloc(lm,nR)  + &
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
               if ( l_TP_form .or. l_anelastic_liquid ) then
                  Buo(lm) =BuoFac*alpha0(nR)*rgrav(nR)*(          &
                  &        rho0(nR)*leg_helper%sR(lm)-ViscHeatFac*&
                  &        (ThExpNb*alpha0(nR)*temp0(nR)+ogrun)*  &
                  &        leg_helper%preR(lm) )
               else
                  Buo(lm) =BuoFac*rho0(nR)*rgrav(nR)*leg_helper%sR(lm)
               end if
               dwdt(lm)=AdvPol_loc+CorPol_loc

               if ( lRmsCalc ) then
                  if ( l_mag_LF .and. nR>n_r_LCR ) then
                     LFPol(lm) =or2(nR)*this%LFrLM(lmP)
                     AdvPol(lm)=AdvPol_loc-LFPol(lm)
                  else
                     AdvPol(lm)=AdvPol_loc
                  end if
                  CorPol(lm)=CorPol_loc
               end if

               if ( l_corr ) then
                  if ( l < l_max .and. l > m ) then
                     CorTor_loc=          two*CorFac*or2(nR) * (  &
                     &                dPhi0(lm)*z_Rloc(lm,nR)   + &
                     &            dTheta3A(lm)*dw_Rloc(lmA,nR)  + &
                     &    or1(nR)*dTheta4A(lm)* w_Rloc(lmA,nR)  + &
                     &            dTheta3S(lm)*dw_Rloc(lmS,nR)  - &
                     &    or1(nR)*dTheta4S(lm)* w_Rloc(lmS,nR)  )
                  else if ( l == l_max ) then
                     CorTor_loc=two*CorFac*or2(nR) * ( &
                     &            dPhi0(lm)*z_Rloc(lm,nR)   )
                  else if ( l == m ) then
                     CorTor_loc=  two*CorFac*or2(nR) * (  &
                     &        dPhi0(lm)*z_Rloc(lm,nR)   + &
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
               ! until here
    
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
            !$OMP END PARALLEL DO
            !PERFOFF
    
            if ( lRmsCalc ) then
    
               if ( l_conv_nl ) then
                  call hIntRms(AdvPol,nR,1,lm_max,0,Adv2hInt(:,nR),st_map)
                  call hIntRms(this%Advt2LM,nR,1,lmP_max,1,Adv2hInt(:,nR),st_map)
                  call hIntRms(this%Advp2LM,nR,1,lmP_max,1,Adv2hInt(:,nR),st_map)
               end if

               lm  =1
               dpt(lm)=zero
               dpp(lm)=zero

               !PERFON('td_cv3')
               !$OMP PARALLEL do default(none) &
               !$OMP private(lm,l,m,lmP,lmPS,lmPA) &
               !$OMP shared(or1,dTheta1A,dTheta1S,this,dpt,dpp,dPhi,nR) &
               !$OMP shared(lm2l,lm2m,lm2lmP,lmP2lmPS,lmP2lmPA,l_max,lm_max)
               do lm=2,lm_max
                  l   =lm2l(lm)
                  m   =lm2m(lm)
                  lmP =lm2lmP(lm)
                  lmPS=lmP2lmPS(lmP)
                  lmPA=lmP2lmPA(lmP)
                  if ( l<l_max .and. l>m ) then
                     dpt(lm)=               or1(nR)* (      &
                     &      -dTheta1A(lm)*this%p1LM(lmPA) + &
                     &       dTheta1S(lm)*this%p1LM(lmPS) - &
                     &                    this%p2LM(lmP)  )                                  
                  else if ( l == l_max ) then
                     dpt(lm)=zero
                  else if ( l == m ) then
                     dpt(lm)=                or1(nR)* (       &
                     &        -dTheta1A(lm)*this%p1LM(lmPA) - &
                     &                      this%p2LM(lmP)  )
                  end if
                  dpp(lm)=or1(nR)*dPhi(lm)*this%p1LM(lmP)
               end do
               !$OMP END PARALLEL DO
               !PERFOFF

               if ( l_TP_form .or. l_anelastic_liquid ) then
                  call hIntRms(leg_helper%dpR,nR,1,lm_max,0, &
                       &       Pre2hInt(:,nR),st_map)
               else
                  call hIntRms(leg_helper%dpR-beta(nR)*leg_helper%preR,&
                       &       nR,1,lm_max,0,Pre2hInt(:,nR),st_map)
               end if
               call hIntRms(dpt,nR,1,lm_max,0,Pre2hInt(:,nR),st_map)
               call hIntRms(dpp,nR,1,lm_max,0,Pre2hInt(:,nR),st_map)

               ! rho* grad(p/rho) = grad(p) - beta*p
               if ( ra /= 0.0_cp ) &
                  call hIntRms(Buo,nR,1,lm_max,0,Buo2hInt(:,nR),st_map)
               if ( l_corr ) then
                  call hIntRms(CorPol,nR,1,lm_max,0,Cor2hInt(:,nR),st_map)
                  call hIntRms(this%CFt2LM,nR,1,lmP_max,1,Cor2hInt(:,nR),st_map)
                  calL hIntRms(this%CFp2LM,nR,1,lmP_max,1,Cor2hInt(:,nR),st_map)
               end if
               if ( l_mag_LF .and. nR>n_r_LCR ) then
                  call hIntRms(LFPol,nR,1,lm_max,0,LF2hInt(:,nR),st_map)
                  call hIntRms(this%LFt2LM,nR,1,lmP_max,1,LF2hInt(:,nR),st_map)
                  call hIntRms(this%LFp2LM,nR,1,lmP_max,1,LF2hInt(:,nR),st_map)
               end if

               do lm=1,lm_max
                  Geo(lm)=CorPol(lm)-leg_helper%dpR(lm)+beta(nR)*leg_helper%preR(lm)
                  CLF(lm)=CorPol(lm)+LFPol(lm)
                  PLF(lm)=LFPol(lm)-leg_helper%dpR(lm)+beta(nR)*leg_helper%preR(lm)
                  Mag(lm)=Geo(lm)+LFPol(lm)
                  Arc(lm)=Mag(lm)+Buo(lm)
                  CIA(lm)=Arc(lm)+AdvPol(lm)
                  !CIA(lm)=CorPol(lm)+Buo(lm)+AdvPol(lm)
               end do
               call hIntRms(Geo,nR,1,lm_max,0,Geo2hInt(:,nR),st_map)
               call hIntRms(CLF,nR,1,lm_max,0,CLF2hInt(:,nR),st_map)
               call hIntRms(PLF,nR,1,lm_max,0,PLF2hInt(:,nR),st_map)
               call hIntRms(Mag,nR,1,lm_max,0,Mag2hInt(:,nR),st_map)
               call hIntRms(Arc,nR,1,lm_max,0,Arc2hInt(:,nR),st_map)
               call hIntRms(CIA,nR,1,lm_max,0,CIA2hInt(:,nR),st_map)

               do lm=1,lm_max
                  lmP =lm2lmP(lm)
                  Geo(lm)=-this%CFt2LM(lmP)-dpt(lm)
                  CLF(lm)=-this%CFt2LM(lmP)+this%LFt2LM(lmP)
                  PLF(lm)=this%LFt2LM(lmP)-dpt(lm)
                  Mag(lm)=Geo(lm)+this%LFt2LM(lmP)
                  Arc(lm)=Mag(lm)
                  CIA(lm)=Arc(lm)+this%Advt2LM(lmP)
                  !CIA(lm)=-this%CFt2LM(lmP)+this%Advt2LM(lmP)
               end do
               call hIntRms(Geo,nR,1,lm_max,0,Geo2hInt(:,nR),st_map)
               call hIntRms(CLF,nR,1,lm_max,0,CLF2hInt(:,nR),st_map)
               call hIntRms(PLF,nR,1,lm_max,0,PLF2hInt(:,nR),st_map)
               call hIntRms(Mag,nR,1,lm_max,0,Mag2hInt(:,nR),st_map)
               call hIntRms(Arc,nR,1,lm_max,0,Arc2hInt(:,nR),st_map)
               call hIntRms(CIA,nR,1,lm_max,0,CIA2hInt(:,nR),st_map)
    
               do lm=1,lm_max
                  lmP =lm2lmP(lm)
                  Geo(lm)=-this%CFp2LM(lmP)-dpp(lm)
                  CLF(lm)=-this%CFp2LM(lmP)+this%LFp2LM(lmP)
                  PLF(lm)=this%LFp2LM(lmP)-dpp(lm)
                  Mag(lm)=Geo(lm)+this%LFp2LM(lmP)
                  Arc(lm)=Mag(lm)
                  CIA(lm)=Arc(lm)+this%Advp2LM(lmP)
                  !CIA(lm)=-this%CFp2LM(lmP)+this%Advp2LM(lmP)
               end do
               call hIntRms(Geo,nR,1,lm_max,0,Geo2hInt(:,nR),st_map)
               call hIntRms(CLF,nR,1,lm_max,0,CLF2hInt(:,nR),st_map)
               call hIntRms(PLF,nR,1,lm_max,0,PLF2hInt(:,nR),st_map)
               call hIntRms(Mag,nR,1,lm_max,0,Mag2hInt(:,nR),st_map)
               call hIntRms(Arc,nR,1,lm_max,0,Arc2hInt(:,nR),st_map)
               call hIntRms(CIA,nR,1,lm_max,0,CIA2hInt(:,nR),st_map)

            end if

            !PERFON('td_cv2')
            !$OMP PARALLEL default(none) &
            !$OMP private(lm,l,m,lmS,lmA,lmP,lmPS,lmPA) &
            !$OMP private(AdvPol_loc,CorPol_loc) &
            !$OMP shared(lm2l,lm2m,lm2lmS,lm2lmA,lm2lmP,lmP2lmpS,lmP2lmPA) &
            !$OMP shared(lm_max,l_max,nR,l_corr,l_conv_nl) &
            !$OMP shared(CorFac,or1,or2,dPhi0,dTheta3A,dTheta3S,dTheta1S,dTheta1A) &
            !$OMP shared(z_Rloc,dPhi,leg_helper,dw_Rloc) &
            !$OMP shared(CorPol,AdvPol,dpdt,this)
            !LIKWID_ON('td_cv2')
            !$OMP DO
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
                     CorPol_loc=                    two*CorFac*or2(nR) * &
                          &               ( -dPhi0(lm) * ( dw_Rloc(lm,nR) &
                          &                  +or1(nR)*leg_helper%dLhw(lm) &
                          &                                             ) &
                          &                  +dTheta3A(lm)*z_Rloc(lmA,nR) &
                          &                  +dTheta3S(lm)*z_Rloc(lmS,nR) &
                          &               )
    
                  else if ( l == l_max ) then
                     CorPol_loc=  two*CorFac*or2(nR) * ( -dPhi0(lm) *  &
                                 ( dw_Rloc(lm,nR) + or1(nR)*leg_helper%dLhw(lm) ) )
    
                  else if ( l == m ) then
                     CorPol_loc=                    two*CorFac*or2(nR) * &
                          &               ( -dPhi0(lm) * ( dw_Rloc(lm,nR) &
                          &                  +or1(nR)*leg_helper%dLhw(lm) &
                          &                                              )&
                          &                 +dTheta3A(lm)*z_Rloc(lmA,nR)  &
                          &               )
    
                  end if
                  !PERFOFF
               else
                  CorPol_loc=zero
               end if
               if ( l_conv_nl ) then
                  !PERFON('td_cv2nl')
                  if ( l > m ) then
                     AdvPol_loc= dTheta1S(lm)*this%AdvtLM(lmPS) - &
                                 dTheta1A(lm)*this%AdvtLM(lmPA) + &
                                     dPhi(lm)*this%AdvpLM(lmP)
                  else if ( l == m ) then
                     AdvPol_loc=-dTheta1A(lm)*this%AdvtLM(lmPA) + &
                                     dPhi(lm)*this%AdvpLM(lmP)
                  end if
                  !PERFOFF
               else
                  AdvPol_loc=zero
               end if
               dpdt(lm)=AdvPol_loc+CorPol_loc
    
            end do ! lm loop
            !$OMP end do 
            !LIKWID_OFF('td_cv2')
            !$OMP END PARALLEL 
            !PERFOFF
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
            if ( l_TP_form ) dVPrLM(1)=this%VPrLM(1)
            if ( l_anel ) then
               if ( l_anelastic_liquid .or. l_TP_form ) then
                  if ( l_mag_nl ) then
                     dsdt_loc=dsdt_loc+                                        &
                     &    ViscHeatFac*hdif_V(1)*temp0(nR)*this%ViscHeatLM(1)+  &
                     &     OhmLossFac*hdif_B(1)*temp0(nR)*this%OhmLossLM(1)
                  else
                     dsdt_loc=dsdt_loc+ &
                     &    ViscHeatFac*hdif_V(1)*temp0(nR)*this%ViscHeatLM(1)
                  end if
               else
                  if ( l_mag_nl ) then
                     dsdt_loc=dsdt_loc+ViscHeatFac*hdif_V(1)*this%ViscHeatLM(1)+ &
                                        OhmLossFac*hdif_B(1)*this%OhmLossLM(1)
                  else
                     dsdt_loc=dsdt_loc+ViscHeatFac*hdif_V(1)*this%ViscHeatLM(1)
                  end if
               end if
            end if
            dsdt(1)=dsdt_loc
    
            !PERFON('td_heat')
            !$OMP PARALLEL DEFAULT(none) &
            !$OMP private(lm,l,m,lmP,lmPS,lmPA,dsdt_loc) &
            !$OMP shared(lm2l,lm2m,lm2lmP,lmP2lmPS,lmP2lmPA) &
            !$OMP shared(lm_max,dsdt,dVSrLM,dTheta1S,dTheta1A,dPhi) &
            !$OMP shared(l_anel,l_anelastic_liquid,l_mag_nl,nR) &
            !$OMP shared(l_TP_form, dVPrLM) &
            !$OMP shared(ViscHeatFac,hdif_V,OhmLossFac,hdif_B,temp0,this)
            !LIKWID_ON('td_heat')
            !$OMP DO
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
                       &    +dTheta1A(lm)*this%VStLM(lmPA) &
                       &    -dPhi(lm)*this%VSpLM(lmP)
               else if ( l == m ) then
                  dsdt_loc=  dTheta1A(lm)*this%VStLM(lmPA) &
                       &     -dPhi(lm)*this%VSpLM(lmP)
               end if
               !PERFOFF
               !PERFON('td_h2')
               if ( l_anel ) then
                  if ( l_anelastic_liquid .or. l_TP_form ) then
                     dsdt_loc = dsdt_loc+ &
                          &     ViscHeatFac*hdif_V(lm)*temp0(nR)*this%ViscHeatLM(lmP)
                     if ( l_mag_nl ) then
                        dsdt_loc = dsdt_loc+ &
                             &     OhmLossFac*hdif_B(lm)*temp0(nR)*this%OhmLossLM(lmP)
                     end if
                  else
                     dsdt_loc = dsdt_loc+ &
                          &     ViscHeatFac*hdif_V(lm)*this%ViscHeatLM(lmP)
                     if ( l_mag_nl ) then
                        dsdt_loc = dsdt_loc+ &
                             &     OhmLossFac*hdif_B(lm)*this%OhmLossLM(lmP)
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
               if ( l_TP_form ) dVPrLM(lm)=this%VPrLM(lmP)
            end do
            !$OMP end do
            !LIKWID_OFF('td_heat')
            !$OMP END PARALLEL
            !PERFOFF
         else
            do lm=2,lm_max
               dsdt(lm)  =0.0_cp
               dVSrLM(lm)=0.0_cp
            end do
         end if

         if ( l_chemical_conv ) then
            dVXirLM(1)=this%VXirLM(1)
            dxidt(1)  =epscXi
    
            !PERFON('td_xi_heat')
            !$OMP PARALLEL DEFAULT(none) &
            !$OMP private(lm,l,m,lmP,lmPS,lmPA,dxidt_loc) &
            !$OMP shared(lm2l,lm2m,lm2lmP,lmP2lmPS,lmP2lmPA) &
            !$OMP shared(lm_max,dxidt,dVXirLM,dTheta1S,dTheta1A,dPhi) &
            !$OMP shared(nR,this) 
            !LIKWID_ON('td_xi_heat')
            !$OMP DO
            do lm=2,lm_max
               l   =lm2l(lm)
               m   =lm2m(lm)
               lmP =lm2lmP(lm)
               lmPS=lmP2lmPS(lmP)
               lmPA=lmP2lmPA(lmP)
               !------ This is horizontal heat advection:
               !PERFON('td_h1')
    
               if ( l > m ) then
                  dxidt_loc= -dTheta1S(lm)*this%VXitLM(lmPS) &
                       &     +dTheta1A(lm)*this%VXitLM(lmPA) &
                       &     -dPhi(lm)*this%VXipLM(lmP)
               else if ( l == m ) then
                  dxidt_loc=  dTheta1A(lm)*this%VXitLM(lmPA) &
                       &     -dPhi(lm)*this%VXipLM(lmP)
               end if
               !PERFOFF
               dVXirLM(lm)=this%VXirLM(lmP)
               dxidt(lm) = dxidt_loc
            end do
            !$OMP end do
            !LIKWID_OFF('td_xi_heat')
            !$OMP END PARALLEL
            !PERFOFF
         end if
    
         if ( l_mag_nl .or. l_mag_kin  ) then
            !PERFON('td_magnl')
    
            !$OMP PARALLEL do default(none) &
            !$OMP private(lm,l,m,lmP,lmPS,lmPA) &
            !$OMP shared(lm_max,lm2l,lm2m,lm2lmP,lmP2lmPS,lmP2lmPA) &
            !$OMP shared(dbdt,djdt,dTheta1S,dTheta1A,dPhi) &
            !$OMP shared(dLh,or4,dVxBhLM,r,nR,this)
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
                       &    -dTheta1A(lm)*this%VxBpLM(lmPA) &
                       &    -dPhi(lm)    *this%VxBtLM(lmP)
               else if ( l == m ) then
                  dbdt(lm)= -dTheta1A(lm)*this%VxBpLM(lmPA) &
                       &    -dPhi(lm)    *this%VxBtLM(lmP)
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
                       dTheta1S(lm)*this%VxBtLM(lmPS) -  &
                       dTheta1A(lm)*this%VxBtLM(lmPA) +  &
                       dPhi(lm)*this%VxBpLM(lmP)  )
               else if ( l == m ) then
                  dVxBhLM(lm)=              r(nR)*r(nR)* ( &
                       - dTheta1A(lm)*this%VxBtLM(lmPA) +  &
                       dPhi(lm)*this%VxBpLM(lmP)  )
               end if
               !PERFOFF
            end do
            !$OMP END PARALLEL DO
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
    
            !----- Stress free boundary, only nl mag. term for poloidal field needed.
            !      Because the radial derivative will be taken, this will contribute to
            !      the other radial grid points.
            dVxBhLM(1)=zero
            dVSrLM(1) =zero
            do lm=2,lm_max
               l   =lm2l(lm)
               m   =lm2m(lm)
               lmP =lm2lmP(lm)
               lmPS=lmP2lmPS(lmP)   ! l-1
               lmPA=lmP2lmPA(lmP)   ! l+1
               if ( l > m ) then
                  dVxBhLM(lm)=r(nR)*r(nR)* (               &
                       & dTheta1S(lm)*this%VxBtLM(lmPS) -  &
                       & dTheta1A(lm)*this%VxBtLM(lmPA) +  &
                       &     dPhi(lm)*this%VxBpLM(lmP)  )
               else if ( l == m ) then ! (l-1) not allowed !
                  dVxBhLM(lm)=r(nR)*r(nR)* (               &
                       - dTheta1A(lm)*this%VxBtLM(lmPA) +  &
                       dPhi(lm)*this%VxBpLM(lmP)  )
               end if
               dVSrLM(lm)=zero
            end do
    
         else
            do lm=1,lm_max
               if ( l_mag ) dVxBhLM(lm)=zero
               dVSrLM(lm) =zero
            end do
         end if
         if ( l_heat ) then
            do lm=1,lm_max
               dVSrLM(lm)=zero
            end do
         end if
         if ( l_chemical_conv ) then
            do lm=1,lm_max
               dVXirLM(lm)=zero
            end do
         end if
         if ( l_TP_form ) then
            do lm=1,lm_max
               dVPrLM(lm)=zero
            end do
         end if
    
      end if  ! boundary ? lvelo ?

   end subroutine get_td
!-----------------------------------------------------------------------------
end module nonlinear_lm_mod
