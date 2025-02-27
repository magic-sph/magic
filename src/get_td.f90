module nonlinear_lm_mod
   !
   ! This module is used to finish the assembly of the nonlinear terms in
   ! :math:`(\ell,m)` space. Derivatives along :math:`\theta` and :math:`\phi`
   ! are handled using recurrence relations.
   !

   use, intrinsic :: iso_c_binding
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: lm_max, lm_maxMag, m_min
   use logic, only : l_anel, l_conv_nl, l_corr, l_heat, l_anelastic_liquid, &
       &             l_mag_nl, l_mag_kin, l_mag_LF, l_conv, l_mag,          &
       &             l_chemical_conv, l_single_matrix, l_double_curl, l_tidal
   use radial_functions, only: r, or2, or1, beta, epscProf, or4, temp0, orho1, l_R
   use physical_parameters, only: CorFac, epsc,  n_r_LCR, epscXi, tidalFac, w_orbit_th !ARS
   use blocking, only: lm2l, lm2lmA, lm2lmS !lm2m ARS
   use horizontal_data, only: dLh, dPhi, dTheta2A, dTheta3A, dTheta4A, dTheta2S, &
       &                      dTheta3S, dTheta4S
   use constants, only: zero, two
   use fields, only: w_Rloc, dw_Rloc, ddw_Rloc, z_Rloc, dz_Rloc, z0v_Rloc, dz0v_Rloc,&
       & wtidal_Rloc, dwtidal_Rloc, ddwtidal_Rloc

   implicit none

   private

   type, public :: nonlinear_lm_t
      !----- Nonlinear terms in lm-space:
      complex(cp), allocatable :: AdvrLM(:), AdvtLM(:), AdvpLM(:)
      complex(cp), allocatable :: VxBrLM(:), VxBtLM(:), VxBpLM(:)
      complex(cp), allocatable :: VStLM(:),  VSpLM(:)
      complex(cp), allocatable :: VXitLM(:),  VXipLM(:)
      complex(cp), allocatable :: heatTermsLM(:)
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: set_zero
      procedure :: get_dsdt
      procedure :: get_dbdt
      procedure :: get_dxidt
      procedure :: get_dwdt_double_curl
      procedure :: get_dwdt
      procedure :: get_dzdt
      procedure :: get_dpdt
   end type nonlinear_lm_t

   integer :: lm_min

contains

   subroutine initialize(this,lmP_max)
      !
      ! Memory allocation of ``get_td`` arrays
      !

      class(nonlinear_lm_t) :: this
      integer, intent(in) :: lmP_max

      allocate( this%AdvrLM(lmP_max), this%AdvtLM(lmP_max), this%AdvpLM(lmP_max))
      allocate( this%VxBrLM(lmP_max), this%VxBtLM(lmP_max), this%VxBpLM(lmP_max))
      bytes_allocated = bytes_allocated + 6*lmP_max*SIZEOF_DEF_COMPLEX

      if ( l_anel) then
         allocate( this%heatTermsLM(lmP_max) )
         bytes_allocated = bytes_allocated+lmP_max*SIZEOF_DEF_COMPLEX
      end if

      if ( l_heat ) then
         allocate(this%VStLM(lmP_max),this%VSpLM(lmP_max))
         bytes_allocated = bytes_allocated + 2*lmP_max*SIZEOF_DEF_COMPLEX
      end if

      if ( l_chemical_conv ) then
         allocate(this%VXitLM(lmP_max),this%VXipLM(lmP_max))
         bytes_allocated = bytes_allocated + 2*lmP_max*SIZEOF_DEF_COMPLEX
      end if

      if ( m_min == 0 ) then
         lm_min = 2
      else
         lm_min = 1
      end if

   end subroutine initialize
!----------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Memory deallocation
      !

      class(nonlinear_lm_t) :: this

      deallocate( this%AdvrLM, this%AdvtLM, this%AdvpLM )
      deallocate( this%VxBrLM, this%VxBtLM, this%VxBpLM )
      if ( l_anel ) deallocate( this%heatTermsLM )
      if ( l_chemical_conv ) deallocate( this%VXitLM, this%VXipLM )
      if ( l_heat ) deallocate( this%VStLM, this%VSpLM )

   end subroutine finalize
!----------------------------------------------------------------------------
   subroutine set_zero(this)
      !
      ! Set all the arrays to zero
      !

      class(nonlinear_lm_t) :: this

      !-- Local variable
      integer :: lm

      !$omp parallel do private(lm)
      do lm=1,lm_max
         this%AdvrLM(lm)=zero
         this%AdvtLM(lm)=zero
         this%AdvpLM(lm)=zero
         this%VxBrLM(lm)=zero
         this%VxBtLM(lm)=zero
         this%VxBpLM(lm)=zero
         if ( l_heat ) then
            this%VStLM(lm)=zero
            this%VSpLM(lm)=zero
         end if
         if ( l_anel ) this%heatTermsLM(lm)=zero
         if ( l_chemical_conv ) then
            this%VXitLM(lm)=zero
            this%VXipLM(lm)=zero
         end if
      end do
      !$omp end parallel do

   end subroutine set_zero
!----------------------------------------------------------------------------
   subroutine get_dwdt(this,nR,nBc,dwdt, eiwt)
      !
      ! This subroutine finishes the assembly of the explicit terms that
      ! enter the equation for the poloidal equation at the radial level
      ! nR.
      !

      !-- Input of variables:
      class(nonlinear_lm_t) :: this

      integer, intent(in) :: nR  ! Radial level
      integer, intent(in) :: nBc ! signifies boundary conditions
      complex(cp), intent(in):: eiwt !for Tidal force

      !-- Output of variables:
      complex(cp), intent(out) :: dwdt(:)

      !-- Local variables:
      integer :: l,lm,lmS,lmA
      complex(cp) :: CorPol_loc

      if (nBc == 0 ) then
         lm =1   ! This is l=0,m=0
         lmA=lm2lmA(lm)
         if ( l_conv_nl ) then
            dwdt(lm)=or2(nR)*this%AdvrLM(lm)
         else
            dwdt(lm)=zero
         end if
         if ( l_corr .and. (.not. l_single_matrix ) ) then
            dwdt(lm)=dwdt(lm)+two*CorFac*or1(nR)*dTheta2A(lm)*(z_Rloc(lmA,nR)+z0v_Rloc(lmA,nR))
         end if

         !$omp parallel default(shared) private(lm,l,lmS,lmA,CorPol_loc)
         if ( l_conv_nl ) then
            !$omp do
            do lm=lm_min,lm_max
               dwdt(lm)=or2(nR)*this%AdvrLM(lm)
            end do
            !$omp end do
         else
            !$omp do
            do lm=lm_min,lm_max
               dwdt(lm)=zero
            end do
            !$omp end do
         end if

         if ( l_corr .and. nBc <2 ) then !ARS
            !$omp do
            do lm=lm_min,lm_max
               l  =lm2l(lm)
               lmS=lm2lmS(lm)
               lmA=lm2lmA(lm)

               if ( l < l_R(nR) ) then
                  CorPol_loc =two*CorFac*or1(nR) * (  &
                       &        dPhi(lm)*(dw_Rloc(lm,nR))+ &!+eiwt*dwtidal_Rloc(lm,nR))+&! phi-deriv of dw/dr
                       &    dTheta2A(lm)*(z_Rloc(lmA,nR)+z0v_Rloc(lmA,nR)) -  & ! sin(theta) dtheta z
                       &    dTheta2S(lm)*(z_Rloc(lmS,nR)+z0v_Rloc(lmS,nR)) )
               else if ( l == l_R(nR) ) then
                  CorPol_loc =two*CorFac*or1(nR) * (  &
                       &        dPhi(lm)*(dw_Rloc(lm,nR))- &!+eiwt*dwtidal_Rloc(lm,nR))-& ! phi-deriv of dw/dr
                       &    dTheta2S(lm)*(z_Rloc(lmS,nR)+z0v_Rloc(lmS,nR)) )
               else
                  CorPol_loc=zero
               end if
               if (l_tidal ) then
                  if ( l < l_R(nR) ) then
                     CorPol_loc =CorPol_loc+two*eiwt*CorFac*or1(nR) * (  &
                          &        dPhi(lm)*dwtidal_Rloc(lm,nR)  ) ! phi-deriv of dw/dr
                  else if ( l == l_R(nR) ) then
                     CorPol_loc =CorPol_loc+two*eiwt*CorFac*or1(nR) * (  &
                          &        dPhi(lm)*dwtidal_Rloc(lm,nR) )   ! phi-deriv of dw/dr
                  end if
               end if
               dwdt(lm)=dwdt(lm)+CorPol_loc
            end do
            !$omp end do
         end if
         !$omp end parallel

      end if

   end subroutine get_dwdt
!----------------------------------------------------------------------------
   subroutine get_dwdt_double_curl(this,nR,nBc,dwdt,dVxVhLM,eiwt)
      !
      ! This subroutine finishes the assembly of the explicit terms that
      ! enter the equation for the poloidal equation in case the double
      ! curl formulation is employed.
      !

      !-- Input of variables:
      class(nonlinear_lm_t) :: this

      integer, intent(in) :: nR  ! Radial level
      integer, intent(in) :: nBc ! signifies boundary conditions
      complex(cp), intent(in):: eiwt !for Tidal force

      !-- Output of variables:
      complex(cp), intent(out) :: dwdt(:)
      complex(cp), intent(out) :: dVxVhLM(:)

      !-- Local variables:
      integer :: l,lm,lmS,lmA
      complex(cp) :: CorPol_loc

      if (nBc == 0 ) then

         lm =1   ! This is l=0,m=0
         lmA=lm2lmA(lm)
         if ( l_conv_nl ) then
            dwdt(lm)=or2(nR)*this%AdvrLM(lm)
         else
            dwdt(lm)=zero
         end if
         if ( l_corr .and. (.not. l_single_matrix ) ) then
            dwdt(lm)=dwdt(lm)+two*CorFac*or1(nR)*dTheta2A(lm)*(z_Rloc(lmA,nR)+z0v_Rloc(lmA,nR))
         end if

         !$omp parallel default(shared) private(lm,l,lmS,lmA,CorPol_loc)
         if ( l_conv_nl ) then
            !$omp do
            do lm=lm_min,lm_max
               dwdt(lm)   =dLh(lm)*or4(nR)*orho1(nR)*this%AdvrLM(lm)
               dVxVhLM(lm)=-orho1(nR)*r(nR)*r(nR)*dLh(lm)*this%AdvtLM(lm)
            end do
            !$omp end do
         else
            !$omp do
            do lm=lm_min,lm_max
               dwdt(lm)   =zero
               dVxVhLM(lm)=zero
            end do
            !$omp end do
         endif

         if ( l_corr ) then
            !$omp do
            do lm=lm_min,lm_max
               l  =lm2l(lm)
               lmS=lm2lmS(lm)
               lmA=lm2lmA(lm)

               if ( l < l_R(nR) ) then
                  CorPol_loc =two*CorFac*or2(nR)*orho1(nR)*(               &
                       &                    dPhi(lm)*(                          &
                       &         -ddw_Rloc(lm,nR)+beta(nR)*dw_Rloc(lm,nR)     + &
                       &             ( beta(nR)*or1(nR)+or2(nR))*               &
                       &                            dLh(lm)*w_Rloc(lm,nR) )   + &
                       &             dTheta3A(lm)*( (dz_Rloc(lmA,nR)+dz0v_Rloc(lmA,nR))-&
                       &                            beta(nR)*(z_Rloc(lmA,nR)+z0v_Rloc(lmA,nR)) ) + &
                       &             dTheta3S(lm)*( (dz_Rloc(lmS,nR)+dz0v_Rloc(lmS,nR))-            &
                       &                            beta(nR)*(z_Rloc(lmS,nR)+z0v_Rloc(lmS,nR)) ) + &
                       &          or1(nR)* (                                    &
                       &             dTheta4A(lm)* (z_Rloc(lmA,nR)+z0v_Rloc(lmA,nR))    &
                       &            -dTheta4S(lm)* (z_Rloc(lmS,nR)+z0v_Rloc(lmS,nR)) ) )
               else if ( l == l_R(nR) ) then
                  CorPol_loc =two*CorFac*or2(nR)*orho1(nR)*(               &
                       &                    dPhi(lm)*(                          &
                       &         -ddw_Rloc(lm,nR)+beta(nR)*dw_Rloc(lm,nR)     + &
                       &             ( beta(nR)*or1(nR)+or2(nR))*               &
                       &                            dLh(lm)*w_Rloc(lm,nR) )   + &
                       &             dTheta3S(lm)*( dz_Rloc(lmS,nR)+dz0v_Rloc(lmS,nR)-            &
                       &                            beta(nR)*(z_Rloc(lmS,nR)+z0v_Rloc(lmS,nR)) ) - &
                       &          or1(nR)* dTheta4S(lm)* (z_Rloc(lmS,nR)+z0v_Rloc(lmS,nR)) )
               else
                  CorPol_loc=zero
               end if
               if ( l_tidal) then
                  if ( l < l_R(nR) ) then
                     CorPol_loc = CorPol_loc+eiwt*two*CorFac*or2(nR)*orho1(nR)*(  &
                          &                    dPhi(lm)*(                          &
                          &         -ddwtidal_Rloc(lm,nR)+beta(nR)*dwtidal_Rloc(lm,nR) + &
                          &             ( beta(nR)*or1(nR)+or2(nR))*               &
                          &                            dLh(lm)*wtidal_Rloc(lm,nR)) )
                  else if ( l == l_R(nR) ) then
                     CorPol_loc =CorPol_loc+two*eiwt*CorFac*or2(nR)*orho1(nR)*( &
                          &                    dPhi(lm)*(                          &
                          &         -ddwtidal_Rloc(lm,nR)+beta(nR)*dwtidal_Rloc(lm,nR) + &
                          &             ( beta(nR)*or1(nR)+or2(nR))*               &
                          &                            dLh(lm)*wtidal_Rloc(lm,nR)) )
                  end if
               end if
               dwdt(lm)=dwdt(lm)+CorPol_loc
            end do
            !$omp end do
         end if
         !$omp end parallel

      else   ! boundary !

         !$omp parallel do
         do lm=1,lm_max
            dVxVhLM(lm)=zero
         end do
         !$omp end parallel do

      end if  ! boundary ? lvelo ?

   end subroutine get_dwdt_double_curl
!----------------------------------------------------------------------------
   subroutine get_dpdt(this,nR,nBc,dpdt,eiwt)
      !
      ! This subroutine finishes the assembly of the explicit terms that
      ! enter the equation for pressure dpdt(:) at the radial level nR.
      !

      !-- Input of variables:
      class(nonlinear_lm_t) :: this

      integer, intent(in) :: nR  ! Radial level
      integer, intent(in) :: nBc ! signifies boundary conditions
      complex(cp), intent(in):: eiwt !for Tidal force

      !-- Output of variables:
      complex(cp), intent(out) :: dpdt(:)

      !-- Local variables:
      integer :: l,lm,lmS,lmA
      complex(cp) :: CorPol_loc

      if (nBc == 0 ) then

         !$omp parallel default(shared) private(lm,l,lmS,lmA,CorPol_loc)
         if ( l_conv_nl ) then
            !$omp do
            do lm=lm_min,lm_max
               dpdt(lm)=-dLh(lm)*this%AdvtLM(lm)
            end do
            !$omp end do
         else
            !$omp do
            do lm=lm_min,lm_max
               dpdt(lm)=zero
            end do
            !$omp end do
         end if

         if ( l_corr ) then
            !$omp do
            do lm=lm_min,lm_max
               l  =lm2l(lm)
               lmS=lm2lmS(lm)
               lmA=lm2lmA(lm)

               if (l_tidal) then
                  if ( l < l_R(nR) ) then
                     CorPol_loc=           two*CorFac*or2(nR) *  &
                          &        ( -dPhi(lm)*((dw_Rloc(lm,nR)+eiwt*dwtidal_Rloc(lm,nR))&
                          &     +or1(nR)*dLh(lm)*(w_Rloc(lm,nR)+eiwt*wtidal_Rloc(lm,nR)) &
                          &                                         ) &
                          &              +dTheta3A(lm)*(z_Rloc(lmA,nR)+z0v_Rloc(lmA,nR)) &
                          &              +dTheta3S(lm)*(z_Rloc(lmS,nR)+z0v_Rloc(lmS,nR)) &
                          &           )
                  else if ( l == l_R(nR) ) then
                     CorPol_loc=           two*CorFac*or2(nR) *  &
                          &   ( -dPhi(lm)  * ((dw_Rloc(lm,nR)+eiwt*dwtidal_Rloc(lm,nR)) &
                          &     +or1(nR)*dLh(lm)*(w_Rloc(lm,nR)+eiwt*wtidal_Rloc(lm,nR))&
                          &                                         ) &
                          &              +dTheta3S(lm)*(z_Rloc(lmS,nR)+z0v_Rloc(lmS,nR)) &
                          &           )
                  else
                     CorPol_loc=zero
                  end if
               else
                  if ( l < l_R(nR) ) then
                     CorPol_loc=           two*CorFac*or2(nR) *  &
                          &           ( -dPhi(lm)  * ((dw_Rloc(lm,nR)) &
                          &            +or1(nR)*dLh(lm)*(w_Rloc(lm,nR)) &
                          &                                         ) &
                          &              +dTheta3A(lm)*(z_Rloc(lmA,nR)+z0v_Rloc(lmA,nR)) &
                          &              +dTheta3S(lm)*(z_Rloc(lmS,nR)+z0v_Rloc(lmS,nR)) &
                          &           )
                  else if ( l == l_R(nR) ) then
                     CorPol_loc=           two*CorFac*or2(nR) *  &
                          &           ( -dPhi(lm)  * ((dw_Rloc(lm,nR)) &
                          &            +or1(nR)*dLh(lm)*(w_Rloc(lm,nR)) &
                          &                                         ) &
                          &              +dTheta3S(lm)*(z_Rloc(lmS,nR)+z0v_Rloc(lmS,nR)) &
                          &           )
                  else
                     CorPol_loc=zero
                  end if
               end if
               dpdt(lm)=dpdt(lm)+CorPol_loc
            end do ! lm loop
            !$omp end do
         end if
         !$omp end parallel
      end if

   end subroutine get_dpdt
!-----------------------------------------------------------------------------
   subroutine get_dzdt(this,nR,nBc,dzdt,eiwt)
      !
      ! This subroutine finishes the assembly of the explicit terms that
      ! enter the toroidal equation dzdt(:) at the radial level nR.
      !

      !-- Input of variables:
      class(nonlinear_lm_t) :: this

      integer, intent(in) :: nR  ! Radial level
      integer, intent(in) :: nBc ! signifies boundary conditions
      complex(cp), intent(in):: eiwt !for Tidal force

      !-- Output of variables:
      complex(cp), intent(out) :: dzdt(:)

      !-- Local variables:
      integer :: l,lm,lmS,lmA
      complex(cp) :: CorTor_loc

      if (nBc == 0 ) then
         lm =1   ! This is l=0,m=0
         lmA=lm2lmA(lm)
         dzdt(lm)=zero!-dTheta1A(lm)*this%AdvpLM(lmA)
         if ( l_corr ) then
            CorTor_loc= two*CorFac*or2(nR) * (                 &
            &                dTheta3A(lm)*dw_Rloc(lmA,nR) +    &
            &        or1(nR)*dTheta4A(lm)* w_Rloc(lmA,nR) )
         else
            CorTor_loc=zero
         end if
         if ( l_tidal .and. l_corr ) then
            CorTor_loc=CorTor_loc+two*eiwt*CorFac*or2(nR) * (       &
                 &                dTheta3A(lm)*dwtidal_Rloc(lmA,nR) +    &
                 &        or1(nR)*dTheta4A(lm)* wtidal_Rloc(lmA,nR) )
         end if
         dzdt(lm)=dzdt(lm)+CorTor_loc

         !$omp parallel default(shared) private(l,lmS,lmA,CorTor_loc)
         if ( l_conv_nl ) then
            !$omp do
            do lm=lm_min,lm_max
               dzdt(lm)=dLh(lm)*this%AdvpLM(lm)
            end do
            !$omp end do
         else
            !$omp do
            do lm=lm_min,lm_max
               dzdt(lm)=zero
            end do
            !$omp end do
         end if

         if ( l_corr ) then
            !$omp do
            do lm=lm_min,lm_max
               l  =lm2l(lm)
               lmS=lm2lmS(lm)
               lmA=lm2lmA(lm)

               if ( l < l_R(nR) ) then
                  CorTor_loc=          two*CorFac*or2(nR) * (  &
                       &                 dPhi(lm)*(z_Rloc(lm,nR)+z0v_Rloc(lm,nR))   + &
                       &            dTheta3A(lm)*(dw_Rloc(lmA,nR))+ &
                       &    or1(nR)*dTheta4A(lm)*(w_Rloc(lmA,nR))  + &
                       &            dTheta3S(lm)*(dw_Rloc(lmS,nR))  - &
                       &    or1(nR)*dTheta4S(lm)*(w_Rloc(lmS,nR))  )
               else if ( l == l_R(nR) ) then
                  CorTor_loc=          two*CorFac*or2(nR) * (  &
                       &                 dPhi(lm)*(z_Rloc(lm,nR)+z0v_Rloc(lm,nR))   + &
                       &            dTheta3S(lm)*(dw_Rloc(lmS,nR)) - &
                       &    or1(nR)*dTheta4S(lm)*(w_Rloc(lmS,nR)) )
               else
                  CorTor_loc=zero
               end if
               if (l_tidal) then
                  if ( l < l_R(nR) ) then
                     CorTor_loc= CorTor_loc +  eiwt*two*CorFac*or2(nR) * (  &
                          &            dTheta3A(lm)*(dwtidal_Rloc(lmA,nR))+ &
                          &    or1(nR)*dTheta4A(lm)*(wtidal_Rloc(lmA,nR))  + &
                          &            dTheta3S(lm)*(dwtidal_Rloc(lmS,nR))  - &
                          &    or1(nR)*dTheta4S(lm)*(wtidal_Rloc(lmS,nR))  )
                  else if ( l == l_R(nR) ) then
                     CorTor_loc=  CorTor_loc +  eiwt*two*CorFac*or2(nR) * (  &
                          &            dTheta3S(lm)*(dwtidal_Rloc(lmS,nR)) - &
                          &    or1(nR)*dTheta4S(lm)*(wtidal_Rloc(lmS,nR)) )
                  end if
               end if
               dzdt(lm)=dzdt(lm)+CorTor_loc
            end do
            !$omp end do
         end if
         !$omp end parallel
      end if

   end subroutine get_dzdt
!-----------------------------------------------------------------------------
   subroutine get_dsdt(this, nR, nBc, dsdt, dVSrLM)
      !
      ! This subroutine finishes the assembly of dsdt(:) and dVSrLM(:)
      ! at the radial level nR.
      !

      class(nonlinear_lm_t) :: this

      !-- Input variables
      integer, intent(in) :: nR  ! Radial level
      integer, intent(in) :: nBc ! Boundary point or not

      !-- Output variables
      complex(cp), intent(out) :: dsdt(:) ! divH(uh*s)
      complex(cp), intent(out) :: dVSrLM(:) ! ur*s

      !-- Local variables
      integer :: lm

      if ( nBc == 0 ) then

         dsdt(1)=epsc*epscProf(nR)!+opr/epsS*divKtemp0(nR)
         if ( l_anel ) then
            if ( l_anelastic_liquid ) then
               dsdt(1)=dsdt(1)+temp0(nR)*this%heatTermsLM(1)
            else
               dsdt(1)=dsdt(1)+this%heatTermsLM(1)
            end if
         end if

         if ( l_anel ) then
            if ( l_anelastic_liquid ) then
               !$omp parallel do
               do lm=lm_min,lm_max
                  dsdt(lm)=dLh(lm)*this%VStLM(lm)+temp0(nR)*this%heatTermsLM(lm)
               end do
               !$omp end parallel do
            else
               !$omp parallel do
               do lm=lm_min,lm_max
                  dsdt(lm)=dLh(lm)*this%VStLM(lm)+this%heatTermsLM(lm)
               end do
               !$omp end parallel do
            end if
         else
            !$omp parallel do
            do lm=lm_min,lm_max
               dsdt(lm)=dLh(lm)*this%VStLM(lm)
            end do
            !$omp end parallel do
         end if

      else   ! boundary !

         !$omp parallel do
         do lm=1,lm_max
            dVSrLM(lm)=zero
         end do
         !$omp end parallel do

      end if  ! boundary ? lvelo ?

   end subroutine get_dsdt
!-----------------------------------------------------------------------------
   subroutine get_dxidt(this, nBc, dxidt, dVXirLM)
      !
      ! This subroutine finishes the assembly of dxidt(:) and dVXirLM(:)
      ! at the radial level nR.
      !

      class(nonlinear_lm_t) :: this

      !-- Input variables
      integer, intent(in) :: nBc ! Boundary point or not

      !-- Output variables
      complex(cp), intent(out) :: dxidt(:) ! divH(uh*xi)
      complex(cp), intent(out) :: dVXirLM(:) ! ur*xi

      !-- Local variables
      integer :: lm

      if ( nBc == 0 ) then
         dxidt(1)  =epscXi
         !$omp parallel do default(shared)
         do lm=lm_min,lm_max
            dxidt(lm)=dLh(lm)*this%VXitLM(lm)
         end do
         !$omp end parallel do
      else   ! boundary !
         !$omp parallel do
         do lm=1,lm_max
            dVXirLM(lm)=zero
         end do
         !$omp end parallel do
      end if  ! boundary ? lvelo ?

   end subroutine get_dxidt
!-----------------------------------------------------------------------------
   subroutine get_dbdt(this, nR, nBc, dbdt, djdt, dVxBhLM)
      !
      ! This subroutine finishes the assembly of dbdt(:), djdt(j)
      ! at the radial level nR.
      !

      class(nonlinear_lm_t) :: this

      !-- Input variables
      integer, intent(in) :: nR  ! Radial level
      integer, intent(in) :: nBc ! Boundary point or not

      !-- Output variables
      complex(cp), intent(out) :: dbdt(:)
      complex(cp), intent(out) :: djdt(:)
      complex(cp), intent(out) :: dVxBhLM(:)

      !-- Local variables
      integer :: lm

      if ( nBc == 0 ) then

         if ( l_mag_nl .or. l_mag_kin  ) then
            !$omp parallel do default(shared)
            do lm=1,lm_max
               dbdt(lm)   = dLh(lm)*this%VxBpLM(lm)
               dVxBhLM(lm)=-dLh(lm)*this%VxBtLM(lm)*r(nR)*r(nR)
               djdt(lm)   = dLh(lm)*or4(nR)*this%VxBrLM(lm)
            end do
            !$omp end parallel do
         else
            if ( l_mag ) then
               !$omp parallel do
               do lm=1,lm_max
                  dbdt(lm)   =zero
                  djdt(lm)   =zero
                  dVxBhLM(lm)=zero
               end do
               !$omp end parallel do
            end if
         end if

      else   ! boundary !

         if ( l_mag_nl .or. l_mag_kin ) then
            dVxBhLM(1)=zero
            !$omp parallel do default(shared)
            do lm=lm_min,lm_max
               dVxBhLM(lm)=-dLh(lm)*this%VxBtLM(lm)*r(nR)*r(nR)
            end do
            !$omp end parallel do
         else if ( l_mag ) then
            !$omp parallel do
            do lm=1,lm_max
               dVxBhLM(lm)=zero
            end do
            !$omp end parallel do
         end if

      end if  ! boundary ? lvelo ?

   end subroutine get_dbdt
!-----------------------------------------------------------------------------
end module nonlinear_lm_mod
