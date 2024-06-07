module nonlinear_lm_mod
   !
   ! This module is used to finish the assembly of the nonlinear terms in
   ! :math:`(\ell,m)` space. Derivatives along :math:`\theta` and :math:`\phi`
   ! are handled using recurrence relations.
   !

   use, intrinsic :: iso_c_binding
   use precision_mod
#ifdef WITH_OMP_GPU
   use mem_alloc, only: bytes_allocated, gpu_bytes_allocated
#else
   use mem_alloc, only: bytes_allocated
#endif
   use truncation, only: lm_max, lm_maxMag
   use logic, only : l_anel, l_conv_nl, l_corr, l_heat, l_anelastic_liquid, &
       &             l_mag_nl, l_mag_kin, l_mag_LF, l_conv, l_mag,          &
       &             l_chemical_conv, l_single_matrix, l_double_curl
   use radial_functions, only: r, or2, or1, beta, epscProf, or4, temp0, orho1, l_R
   use physical_parameters, only: CorFac, epsc,  n_r_LCR, epscXi
   use blocking, only: lm2l, lm2lmA, lm2lmS
   use horizontal_data, only: dLh, dPhi, dTheta2A, dTheta3A, dTheta4A, dTheta2S, &
       &                      dTheta3S, dTheta4S
   use constants, only: zero, two
   use fields, only: w_Rloc, dw_Rloc, ddw_Rloc, z_Rloc, dz_Rloc

   implicit none

   private

   type, public :: nonlinear_lm_t
      !----- Nonlinear terms in lm-space:
      complex(cp), allocatable :: AdvrLM(:), AdvtLM(:), AdvpLM(:)
      complex(cp), allocatable :: VxBrLM(:), VxBtLM(:), VxBpLM(:)
      complex(cp), allocatable :: VStLM(:), VSpLM(:)
      complex(cp), allocatable :: VXitLM(:), VXipLM(:)
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

contains

   subroutine initialize(this, lmP_max)
      !
      ! Memory allocation of ``get_td`` arrays
      !

      class(nonlinear_lm_t) :: this
      integer, intent(in) :: lmP_max

      allocate( this%AdvrLM(lmP_max), this%AdvtLM(lmP_max), this%AdvpLM(lmP_max))
      bytes_allocated = bytes_allocated + 3*lmP_max*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
      gpu_bytes_allocated = gpu_bytes_allocated + 3*lmP_max*SIZEOF_DEF_COMPLEX
#endif

      if ( l_mag ) then
         allocate( this%VxBrLM(lmP_max), this%VxBtLM(lmP_max), this%VxBpLM(lmP_max))
         bytes_allocated = bytes_allocated + 3*lmP_max*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         gpu_bytes_allocated = gpu_bytes_allocated + 3*lmP_max*SIZEOF_DEF_COMPLEX
#endif
      else
         allocate( this%VxBrLM(1), this%VxBtLM(1), this%VxBpLM(1))
      end if

      if ( l_anel) then
         allocate( this%heatTermsLM(lmP_max) )
         bytes_allocated = bytes_allocated+lmP_max*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         gpu_bytes_allocated = gpu_bytes_allocated+lmP_max*SIZEOF_DEF_COMPLEX
#endif
      else
         allocate( this%heatTermsLM(1) )
      end if

      if ( l_heat ) then
         allocate(this%VStLM(lmP_max),this%VSpLM(lmP_max))
         bytes_allocated = bytes_allocated + 2*lmP_max*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         gpu_bytes_allocated = gpu_bytes_allocated + 2*lmP_max*SIZEOF_DEF_COMPLEX
#endif
      else
         allocate(this%VStLM(1),this%VSpLM(1))
      end if

      if ( l_chemical_conv ) then
         allocate(this%VXitLM(lmP_max),this%VXipLM(lmP_max))
         bytes_allocated = bytes_allocated + 2*lmP_max*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         gpu_bytes_allocated = gpu_bytes_allocated + 2*lmP_max*SIZEOF_DEF_COMPLEX
#endif
      else
         allocate(this%VXitLM(1),this%VXipLM(1))
      end if

#ifdef WITH_OMP_GPU
      !$omp target enter data map(alloc: this)
#endif

   end subroutine initialize
!----------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Memory deallocation
      !

      class(nonlinear_lm_t) :: this

#ifdef WITH_OMP_GPU
      !$omp target exit data map(release: this)
#endif

      deallocate( this%AdvrLM, this%AdvtLM, this%AdvpLM )
      deallocate( this%VxBrLM, this%VxBtLM, this%VxBpLM )
      deallocate( this%heatTermsLM, this%VXitLM, this%VXipLM )
      deallocate( this%VStLM, this%VSpLM )

   end subroutine finalize
!----------------------------------------------------------------------------
   subroutine set_zero(this)
      !
      ! Set all the arrays to zero
      !

      class(nonlinear_lm_t) :: this

      !-- Local variable
      integer :: lm

#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do
#else
      !$omp parallel do private(lm)
#endif
      do lm=1,lm_max
         this%AdvrLM(lm)=zero
         this%AdvtLM(lm)=zero
         this%AdvpLM(lm)=zero
         if ( l_mag ) then
            this%VxBrLM(lm)=zero
            this%VxBtLM(lm)=zero
            this%VxBpLM(lm)=zero
         end if
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
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#else
      !$omp end parallel do
#endif

   end subroutine set_zero
!----------------------------------------------------------------------------
   subroutine get_dwdt(this,nR,nBc,dwdt)
      !
      ! This subroutine finishes the assembly of the explicit terms that
      ! enter the equation for the poloidal equation at the radial level
      ! nR.
      !

      !-- Input of variables:
      class(nonlinear_lm_t) :: this

      integer, intent(in) :: nR  ! Radial level
      integer, intent(in) :: nBc ! signifies boundary conditions

      !-- Output of variables:
      complex(cp), intent(out) :: dwdt(:)

      !-- Local variables:
      integer :: l,lm,lmS,lmA
      complex(cp) :: CorPol_loc

      if (nBc == 0 ) then

#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do private(l,lmS,lmA) &
         !$omp private(CorPol_loc)
#else
         !$omp parallel do default(shared) private(l,lmS,lmA) &
         !$omp private(CorPol_loc)
#endif
         do lm=1,lm_max
            l  =lm2l(lm)
            lmS=lm2lmS(lm)
            lmA=lm2lmA(lm)

            if ( l == 0 ) then ! This is l=0,m=0
               if ( l_conv_nl ) then
                  dwdt(lm)=or2(nR)*this%AdvrLM(lm)
               else
                  dwdt(lm)=zero
               end if
               if ( l_corr .and. (.not. l_single_matrix) ) then
                  dwdt(lm)=dwdt(lm)+two*CorFac*or1(nR)*dTheta2A(lm)*z_Rloc(lmA,nR)
               end if
            else
               if ( l_conv_nl ) then
                  dwdt(lm)=or2(nR)*this%AdvrLM(lm)
               else
                  dwdt(lm)=zero
               endif

               if ( l_corr .and. nBc /= 2 ) then
                  if ( l < l_R(nR) ) then
                     CorPol_loc =two*CorFac*or1(nR) * (  &
                     &        dPhi(lm)*dw_Rloc(lm,nR) +  & ! phi-deriv of dw/dr
                     &    dTheta2A(lm)*z_Rloc(lmA,nR) -  & ! sin(theta) dtheta z
                     &    dTheta2S(lm)*z_Rloc(lmS,nR) )
                  else if ( l == l_R(nR) ) then
                     CorPol_loc =two*CorFac*or1(nR) * (  &
                     &        dPhi(lm)*dw_Rloc(lm,nR) -  & ! phi-deriv of dw/dr
                     &    dTheta2S(lm)*z_Rloc(lmS,nR) )
                  else
                     CorPol_loc=zero
                  end if
                  dwdt(lm)=dwdt(lm)+CorPol_loc
               end if
            end if
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#else
         !$omp end parallel do
#endif
      end if

   end subroutine get_dwdt
!----------------------------------------------------------------------------
   subroutine get_dwdt_double_curl(this,nR,nBc,dwdt,dVxVhLM)
      !
      ! This subroutine finishes the assembly of the explicit terms that
      ! enter the equation for the poloidal equation in case the double
      ! curl formulation is employed.
      !

      !-- Input of variables:
      class(nonlinear_lm_t) :: this

      integer, intent(in) :: nR  ! Radial level
      integer, intent(in) :: nBc ! signifies boundary conditions

      !-- Output of variables:
      complex(cp), intent(out) :: dwdt(:)
      complex(cp), intent(out) :: dVxVhLM(:)

      !-- Local variables:
      integer :: l,lm,lmS,lmA
      complex(cp) :: CorPol_loc

      if (nBc == 0 ) then

#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do private(l,lmS,lmA) &
         !$omp private(CorPol_loc)
#else
         !$omp parallel do default(shared) private(l,lmS,lmA) &
         !$omp private(CorPol_loc)
#endif
         do lm=1,lm_max
            l  =lm2l(lm)
            lmS=lm2lmS(lm)
            lmA=lm2lmA(lm)

            if ( l == 0 ) then
               if ( l_conv_nl ) then
                  dwdt(lm)=or2(nR)*this%AdvrLM(lm)
               else
                  dwdt(lm)=zero
               end if
               if ( l_corr .and. (.not. l_single_matrix ) ) then
                  dwdt(lm)=dwdt(lm)+two*CorFac*or1(nR)*dTheta2A(lm)*z_Rloc(lmA,nR)
               end if
            else
               if ( l_conv_nl ) then
                  dwdt(lm)   =dLh(lm)*or4(nR)*orho1(nR)*this%AdvrLM(lm)
                  dVxVhLM(lm)=-orho1(nR)*r(nR)*r(nR)*dLh(lm)*this%AdvtLM(lm)
               else
                  dwdt(lm)   =zero
                  dVxVhLM(lm)=zero
               end if

               if ( l_corr ) then
                  if ( l < l_R(nR) ) then
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
                  else if ( l == l_R(nR) ) then
                     CorPol_loc =two*CorFac*or2(nR)*orho1(nR)*(               &
                     &                    dPhi(lm)*(                          &
                     &         -ddw_Rloc(lm,nR)+beta(nR)*dw_Rloc(lm,nR)     + &
                     &             ( beta(nR)*or1(nR)+or2(nR))*               &
                     &                            dLh(lm)*w_Rloc(lm,nR) )   + &
                     &             dTheta3S(lm)*( dz_Rloc(lmS,nR)-            &
                     &                            beta(nR)*z_Rloc(lmS,nR) ) - &
                     &          or1(nR)* dTheta4S(lm)* z_Rloc(lmS,nR) )
                  else
                     CorPol_loc=zero
                  end if
                  dwdt(lm)=dwdt(lm)+CorPol_loc
               end if
            end if
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#else
         !$omp end parallel do
#endif

      else   ! boundary !

#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do
#else
         !$omp parallel do
#endif
         do lm=1,lm_max
            dVxVhLM(lm)=zero
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#else
         !$omp end parallel do
#endif

      end if  ! boundary ? lvelo ?

   end subroutine get_dwdt_double_curl
!----------------------------------------------------------------------------
   subroutine get_dpdt(this,nR,nBc,dpdt)
      !
      ! This subroutine finishes the assembly of the explicit terms that
      ! enter the equation for pressure dpdt(:) at the radial level nR.
      !

      !-- Input of variables:
      class(nonlinear_lm_t) :: this

      integer, intent(in) :: nR  ! Radial level
      integer, intent(in) :: nBc ! signifies boundary conditions

      !-- Output of variables:
      complex(cp), intent(out) :: dpdt(:)

      !-- Local variables:
      integer :: l,lm,lmS,lmA
      complex(cp) :: CorPol_loc

      if (nBc == 0 ) then

#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do private(l,lmS,lmA) &
         !$omp private(CorPol_loc)
#else
         !$omp parallel do default(shared) private(l,lmS,lmA) &
         !$omp private(CorPol_loc)
#endif
         do lm=1,lm_max
            l  =lm2l(lm)
            lmS=lm2lmS(lm)
            lmA=lm2lmA(lm)

            if ( l > 0 ) then
               if ( l_conv_nl ) then
                  dpdt(lm)=-dLh(lm)*this%AdvtLM(lm)
               else
                  dpdt(lm)=zero
               end if

               if ( l_corr ) then
                  if ( l < l_R(nR) ) then
                     CorPol_loc=           two*CorFac*or2(nR) *  &
                     &           ( -dPhi(lm)  * ( dw_Rloc(lm,nR) &
                     &            +or1(nR)*dLh(lm)*w_Rloc(lm,nR) &
                     &                                         ) &
                     &              +dTheta3A(lm)*z_Rloc(lmA,nR) &
                     &              +dTheta3S(lm)*z_Rloc(lmS,nR) &
                     &           )
                  else if ( l == l_R(nR) ) then
                     CorPol_loc=           two*CorFac*or2(nR) *  &
                     &           ( -dPhi(lm)  * ( dw_Rloc(lm,nR) &
                     &            +or1(nR)*dLh(lm)*w_Rloc(lm,nR) &
                     &                                         ) &
                     &              +dTheta3S(lm)*z_Rloc(lmS,nR) &
                     &           )
                  else
                     CorPol_loc=zero
                  end if
                  dpdt(lm)=dpdt(lm)+CorPol_loc
               end if
            end if
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#else
         !$omp end parallel do
#endif
      end if

   end subroutine get_dpdt
!-----------------------------------------------------------------------------
   subroutine get_dzdt(this,nR,nBc,dzdt)
      !
      ! This subroutine finishes the assembly of the explicit terms that
      ! enter the toroidal equation dzdt(:) at the radial level nR.
      !

      !-- Input of variables:
      class(nonlinear_lm_t) :: this

      integer, intent(in) :: nR  ! Radial level
      integer, intent(in) :: nBc ! signifies boundary conditions

      !-- Output of variables:
      complex(cp), intent(out) :: dzdt(:)

      !-- Local variables:
      integer :: l,lm,lmS,lmA
      complex(cp) :: CorTor_loc

      if (nBc == 0 ) then
#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do private(l,lmS,lmA) &
         !$omp private(CorTor_loc)
#else
         !$omp parallel do default(shared) private(l,lmS,lmA) &
         !$omp private(CorTor_loc)
#endif
         do lm=1,lm_max
            l  =lm2l(lm)
            lmS=lm2lmS(lm)
            lmA=lm2lmA(lm)
            
            if ( l == 0 ) then
               dzdt(lm)=zero!-dTheta1A(lm)*this%AdvpLM(lmA)
               if ( l_corr ) then
                  dzdt(lm)=dzdt(lm)+ two*CorFac*or2(nR) * (          &
                  &                dTheta3A(lm)*dw_Rloc(lmA,nR) +    &
                  &        or1(nR)*dTheta4A(lm)* w_Rloc(lmA,nR) )
               end if
            else
               if ( l_conv_nl ) then
                  dzdt(lm)=dLh(lm)*this%AdvpLM(lm)
               else
                  dzdt(lm)=zero
               end if

               if ( l_corr ) then
                  if ( l < l_R(nR) ) then
                     CorTor_loc=          two*CorFac*or2(nR) * (  &
                     &                 dPhi(lm)*z_Rloc(lm,nR)   + &
                     &            dTheta3A(lm)*dw_Rloc(lmA,nR)  + &
                     &    or1(nR)*dTheta4A(lm)* w_Rloc(lmA,nR)  + &
                     &            dTheta3S(lm)*dw_Rloc(lmS,nR)  - &
                     &    or1(nR)*dTheta4S(lm)* w_Rloc(lmS,nR)  )
                  else if ( l == l_R(nR) ) then
                     CorTor_loc=          two*CorFac*or2(nR) * (  &
                     &                 dPhi(lm)*z_Rloc(lm,nR)   + &
                     &            dTheta3S(lm)*dw_Rloc(lmS,nR)  - &
                     &    or1(nR)*dTheta4S(lm)* w_Rloc(lmS,nR)  )
                  else
                     CorTor_loc=zero
                  end if
                  dzdt(lm)=dzdt(lm)+CorTor_loc
               end if
            end if
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#else
         !$omp end parallel do
#endif
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
      integer :: lm, l

      if ( nBc == 0 ) then

         if ( l_anel ) then
            if ( l_anelastic_liquid ) then
#ifdef WITH_OMP_GPU
               !$omp target teams distribute parallel do private(l)
#else
               !$omp parallel do private(l)
#endif
               do lm=1,lm_max
                  l=lm2l(lm)
                  if ( l == 0 ) then
                     dsdt(lm)=epsc*epscProf(nR)+temp0(nR)*this%heatTermsLM(lm)
                  else
                     dsdt(lm)=dLh(lm)*this%VStLM(lm)+temp0(nR)*this%heatTermsLM(lm)
                  end if
               end do
#ifdef WITH_OMP_GPU
               !$omp end target teams distribute parallel do
#else
               !$omp end parallel do
#endif
            else
#ifdef WITH_OMP_GPU
               !$omp target teams distribute parallel do private(l)
#else
               !$omp parallel do private(l)
#endif
               do lm=1,lm_max
                  l=lm2l(lm)
                  if ( l == 0 ) then
                     dsdt(lm)=epsc*epscProf(nR)+this%heatTermsLM(lm)
                  else
                     dsdt(lm)=dLh(lm)*this%VStLM(lm)+this%heatTermsLM(lm)
                  end if
               end do
#ifdef WITH_OMP_GPU
               !$omp end target teams distribute parallel do
#else
               !$omp end parallel do
#endif
            end if

         else ! Boussinesq

#ifdef WITH_OMP_GPU
            !$omp target teams distribute parallel do private(l)
#else
            !$omp parallel do private(l)
#endif
            do lm=1,lm_max
               l=lm2l(lm)
               if ( l == 0 ) then
                  dsdt(lm)=epsc*epscProf(nR)
               else
                  dsdt(lm)=dLh(lm)*this%VStLM(lm)
               end if
            end do
#ifdef WITH_OMP_GPU
            !$omp end target teams distribute parallel do
#else
            !$omp end parallel do
#endif
         end if

      else   ! boundary !

#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do
#else
         !$omp parallel do
#endif
         do lm=1,lm_max
            dVSrLM(lm)=zero
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#else
         !$omp end parallel do
#endif

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
      integer :: lm, l

      if ( nBc == 0 ) then
#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do private(l)
#else
         !$omp parallel do private(l)
#endif
         do lm=1,lm_max
            l=lm2l(lm)
            if ( l == 0 ) then
               dxidt(lm)=epscXi
            else
               dxidt(lm)=dLh(lm)*this%VXitLM(lm)
            end if
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#else
         !$omp end parallel do
#endif
      else   ! boundary !
#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do
#else
         !$omp parallel do
#endif
         do lm=1,lm_max
            dVXirLM(lm)=zero
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#else
         !$omp end parallel do
#endif
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
      integer :: lm, l

      if ( nBc == 0 ) then

         if ( l_mag_nl .or. l_mag_kin  ) then
#ifdef WITH_OMP_GPU
            !$omp target teams distribute parallel do
#else
            !$omp parallel do default(shared)
#endif
            do lm=1,lm_max
               dbdt(lm)   = dLh(lm)*this%VxBpLM(lm)
               dVxBhLM(lm)=-dLh(lm)*this%VxBtLM(lm)*r(nR)*r(nR)
               djdt(lm)   = dLh(lm)*or4(nR)*this%VxBrLM(lm)
            end do
#ifdef WITH_OMP_GPU
            !$omp end target teams distribute parallel do
#else
            !$omp end parallel do
#endif
         else
            if ( l_mag ) then
#ifdef WITH_OMP_GPU
               !$omp target teams distribute parallel do
#else
               !$omp parallel do
#endif
               do lm=1,lm_max
                  dbdt(lm)   =zero
                  djdt(lm)   =zero
                  dVxBhLM(lm)=zero
               end do
#ifdef WITH_OMP_GPU
               !$omp end target teams distribute parallel do
#else
               !$omp end parallel do
#endif
            end if
         end if

      else   ! boundary !

         if ( l_mag_nl .or. l_mag_kin ) then
#ifdef WITH_OMP_GPU
            !$omp target teams distribute parallel do private(l)
#else
            !$omp parallel do default(shared) private(l)
#endif
            do lm=1,lm_max
               l=lm2l(lm)
               if ( l == 0 ) then
                  dVxBhLM(lm)=zero
               else
                  dVxBhLM(lm)=-dLh(lm)*this%VxBtLM(lm)*r(nR)*r(nR)
               end if
            end do
#ifdef WITH_OMP_GPU
            !$omp end target teams distribute parallel do
#else
            !$omp end parallel do
#endif
         else
            if ( l_mag ) then
#ifdef WITH_OMP_GPU
               !$omp target teams distribute parallel do
#else
               !$omp parallel do
#endif
               do lm=1,lm_max
                  dVxBhLM(lm)=zero
               end do
#ifdef WITH_OMP_GPU
               !$omp end target teams distribute parallel do
#else
               !$omp end parallel do
#endif
            end if
         end if

      end if  ! boundary ? lvelo ?

   end subroutine get_dbdt
!-----------------------------------------------------------------------------
end module nonlinear_lm_mod

module nonlinear_lm_2d_mod
   !
   ! This module is used to finish the assembly of the nonlinear terms in
   ! :math:`(\ell,m)` space. Derivatives along :math:`\theta` and :math:`\phi`
   ! are handled using recurrence relations.
   !

   use, intrinsic :: iso_c_binding
   use precision_mod
#ifdef WITH_OMP_GPU
   use mem_alloc, only: bytes_allocated, gpu_bytes_allocated
#else
   use mem_alloc, only: bytes_allocated
#endif
   use radial_data, only: nRstart, nRstop, n_r_cmb, n_r_icb, nRstartMag, &
       &                  nRstopMag
   use truncation, only: lm_max, lm_maxMag
   use logic, only : l_anel, l_conv_nl, l_corr, l_heat, l_anelastic_liquid, &
       &             l_mag_nl, l_mag_kin, l_mag_LF, l_conv, l_mag,          &
       &             l_chemical_conv, l_single_matrix, l_double_curl,       &
       &             l_adv_curl, l_parallel_solve, l_temperature_diff
   use radial_functions, only: r, or2, or1, beta, epscProf, or4, temp0, orho1, l_R
   use physical_parameters, only: CorFac, epsc,  n_r_LCR, epscXi
   use blocking, only: lm2l, lm2lmA, lm2lmS
   use horizontal_data, only: dLh, dPhi, dTheta2A, dTheta3A, dTheta4A, dTheta2S, &
       &                      dTheta3S, dTheta4S
   use constants, only: zero, two
   use fields, only: w_Rloc, dw_Rloc, ddw_Rloc, z_Rloc, dz_Rloc

   implicit none

   private

   type, public :: nonlinear_lm_2d_t
      !----- Nonlinear terms in lm-space:
      complex(cp), allocatable :: AdvrLM(:,:), AdvtLM(:,:), AdvpLM(:,:)
      complex(cp), allocatable :: VxBrLM(:,:), VxBtLM(:,:), VxBpLM(:,:)
      complex(cp), allocatable :: VStLM(:,:),  VSpLM(:,:)
      complex(cp), allocatable :: VXitLM(:,:), VXipLM(:,:)
      complex(cp), allocatable :: heatTermsLM(:,:)
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: set_zero
      procedure :: get_td
      procedure :: get_dsdt
      procedure :: get_dbdt
      procedure :: get_dxidt
      procedure :: get_dwdt_double_curl
      procedure :: get_dwdt
      procedure :: get_dzdt
      procedure :: get_dpdt
   end type nonlinear_lm_2d_t

contains

   subroutine initialize(this, lmP_max, nRstart, nRstop)
      !
      ! Memory allocation of ``get_td`` arrays
      !

      class(nonlinear_lm_2d_t) :: this
      integer, intent(in) :: lmP_max
      integer, intent(in) :: nRstart
      integer, intent(in) :: nRstop

      allocate( this%AdvrLM(lmP_max,nRstart:nRstop) )
      allocate( this%AdvtLM(lmP_max,nRstart:nRstop) )
      allocate( this%AdvpLM(lmP_max,nRstart:nRstop) )
      bytes_allocated = bytes_allocated + 3*lmP_max*(nRstop-nRstart+1) * &
      &                 SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
      gpu_bytes_allocated = gpu_bytes_allocated + 3*lmP_max*(nRstop-nRstart+1)* &
      &                     SIZEOF_DEF_COMPLEX
#endif

      if ( l_mag ) then
         allocate( this%VxBrLM(lmP_max,nRstart:nRstop) )
         allocate( this%VxBtLM(lmP_max,nRstart:nRstop) )
         allocate( this%VxBpLM(lmP_max,nRstart:nRstop) )
         bytes_allocated = bytes_allocated + 3*lmP_max*(nRstop-nRstart+1) * &
         &                 SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         gpu_bytes_allocated = gpu_bytes_allocated + 3*lmP_max*(nRstop-nRstart+1)* &
         &                     SIZEOF_DEF_COMPLEX
#endif
      else
         allocate( this%VxBrLM(1,1), this%VxBtLM(1,1), this%VxBpLM(1,1))
      end if

      if ( l_anel) then
         allocate( this%heatTermsLM(lmP_max,nRstart:nRstop) )
         bytes_allocated = bytes_allocated + lmP_max*(nRstop-nRstart+1) * &
         &                 SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         gpu_bytes_allocated = gpu_bytes_allocated + lmP_max*(nRstop-nRstart+1)* &
         &                     SIZEOF_DEF_COMPLEX
#endif
      else
         allocate( this%heatTermsLM(1,1) )
      end if

      if ( l_heat ) then
         allocate( this%VStLM(lmP_max,nRstart:nRstop) )
         allocate( this%VSpLM(lmP_max,nRstart:nRstop) )
         bytes_allocated = bytes_allocated + 2*lmP_max*(nRstop-nRstart+1) * &
         &                 SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         gpu_bytes_allocated = gpu_bytes_allocated + 2*lmP_max*(nRstop-nRstart+1)* &
         &                     SIZEOF_DEF_COMPLEX
#endif
      else
         allocate(this%VStLM(1,1),this%VSpLM(1,1))
      end if

      if ( l_chemical_conv ) then
         allocate( this%VXitLM(lmP_max,nRstart:nRstop) )
         allocate( this%VXipLM(lmP_max,nRstart:nRstop) )
         bytes_allocated = bytes_allocated + 2*lmP_max*(nRstop-nRstart+1) * &
         &                 SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         gpu_bytes_allocated = gpu_bytes_allocated + 2*lmP_max*(nRstop-nRstart+1)* &
         &                     SIZEOF_DEF_COMPLEX
#endif
      else
         allocate(this%VXitLM(1,1),this%VXipLM(1,1))
      end if

#ifdef WITH_OMP_GPU
      !$omp target enter data map(alloc: this)
#endif

   end subroutine initialize
!----------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Memory deallocation
      !

      class(nonlinear_lm_2d_t) :: this

#ifdef WITH_OMP_GPU
      !$omp target exit data map(release: this)
#endif

      deallocate( this%AdvrLM, this%AdvtLM, this%AdvpLM )
      deallocate( this%VxBrLM, this%VxBtLM, this%VxBpLM )
      deallocate( this%heatTermsLM, this%VXitLM, this%VXipLM )
      deallocate( this%VStLM, this%VSpLM )

   end subroutine finalize
!----------------------------------------------------------------------------
   subroutine set_zero(this)
      !
      ! Set all the arrays to zero
      !

      class(nonlinear_lm_2d_t) :: this

      !-- Local variable
      integer :: lm, nR

#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do collapse(2)
#else
      !$omp parallel do private(lm)
#endif
      do nR=nRstart,nRstop
         do lm=1,lm_max
            this%AdvrLM(lm,nR)=zero
            this%AdvtLM(lm,nR)=zero
            this%AdvpLM(lm,nR)=zero
            if ( l_mag ) then
               this%VxBrLM(lm,nR)=zero
               this%VxBtLM(lm,nR)=zero
               this%VxBpLM(lm,nR)=zero
            end if
            if ( l_heat ) then
               this%VStLM(lm,nR)=zero
               this%VSpLM(lm,nR)=zero
            end if
            if ( l_anel ) this%heatTermsLM(lm,nR)=zero
            if ( l_chemical_conv ) then
               this%VXitLM(lm,nR)=zero
               this%VXipLM(lm,nR)=zero
            end if
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#else
      !$omp end parallel do
#endif

   end subroutine set_zero
!----------------------------------------------------------------------------
   subroutine get_td(this,lPressNext,dVxVhLM,dVxBhLM,dwdt,dzdt,dpdt,dsdt, &
              &      dxidt,dbdt,djdt)
      !
      !  Purpose of this to calculate time derivatives
      !  ``dwdt``,``dzdt``,``dpdt``,``dsdt``,``dxidt``,``dbdt``,``djdt``
      !  and auxiliary arrays ``dVxBhLM``, ``dVxVhLM``
      !  from non-linear terms in spectral form
      !

      !-- Input of variables:
      class(nonlinear_lm_2d_t) :: this

      logical, intent(in) :: lPressNext

      !-- Output of variables:
      complex(cp), intent(out) :: dwdt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dzdt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dpdt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dsdt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dxidt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dbdt(lm_maxMag,nRstartMag:nRstopMag)
      complex(cp), intent(out) :: djdt(lm_maxMag,nRstartMag:nRstopMag)
      complex(cp), intent(out) :: dVxBhLM(lm_maxMag,nRstartMag:nRstopMag)
      complex(cp), intent(out) :: dVxVhLM(lm_max,nRstart:nRstop)

      !-- Local variables:
      integer :: l,lm,lmS,lmA,nR
      complex(cp) :: AdvPol_loc,CorPol_loc,AdvTor_loc,CorTor_loc,dsdt_loc

      if ( l_conv ) then

#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do collapse(2) &
         !$omp private(l,lmS,lmA) &
         !$omp private(AdvPol_loc,CorPol_loc,AdvTor_loc,CorTor_loc)
#else
         !$omp parallel do default(shared) private(lm,l,lmS,lmA) &
         !$omp private(AdvPol_loc,CorPol_loc,AdvTor_loc,CorTor_loc)
#endif
         do nR=nRstart,nRstop
            do lm=1,lm_max
               l   =lm2l(lm)
               lmS =lm2lmS(lm)
               lmA =lm2lmA(lm)

               if ( l == 0 ) then ! This is l=0,m=0
                  if ( l_conv_nl ) then
                     AdvPol_loc=or2(nR)*this%AdvrLM(lm,nR)
                     AdvTor_loc=zero!-dTheta1A(lm)*this%AdvpLM(lmA,nR)
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
                     dwdt(lm,nR)=AdvPol_loc!+CorPol_loc
                  else
                     dwdt(lm,nR)=AdvPol_loc+CorPol_loc
                  end if

                  dzdt(lm,nR)=AdvTor_loc+CorTor_loc

               else ! l /= 0

                  if ( l_double_curl ) then ! Pressure is not needed

                     if ( l_corr ) then
                        if ( l < l_R(nR) ) then
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
                        else if ( l == l_R(nR) ) then
                           CorPol_loc =two*CorFac*or2(nR)*orho1(nR)*(               &
                           &                    dPhi(lm)*(                          &
                           &         -ddw_Rloc(lm,nR)+beta(nR)*dw_Rloc(lm,nR)     + &
                           &             ( beta(nR)*or1(nR)+or2(nR))*               &
                           &                            dLh(lm)*w_Rloc(lm,nR) )   + &
                           &             dTheta3S(lm)*( dz_Rloc(lmS,nR)-            &
                           &                            beta(nR)*z_Rloc(lmS,nR) ) - &
                           &          or1(nR)* dTheta4S(lm)* z_Rloc(lmS,nR) )
                        else
                           CorPol_loc=zero
                        end if
                     else
                        CorPol_loc=zero
                     end if

                     if ( l_conv_nl ) then
                        AdvPol_loc    =dLh(lm)*or4(nR)*orho1(nR)*this%AdvrLM(lm,nR)
                        dVxVhLM(lm,nR)=-orho1(nR)*r(nR)*r(nR)*dLh(lm)*this%AdvtLM(lm,nR)
                     else
                        AdvPol_loc    =zero
                        dVxVhLM(lm,nR)=zero
                     endif

                  else ! We don't use the double curl

                     if ( l_corr ) then
                        if ( l < l_R(nR) ) then
                           CorPol_loc =two*CorFac*or1(nR) * (  &
                           &        dPhi(lm)*dw_Rloc(lm,nR) +  & ! phi-deriv of dw/dr
                           &    dTheta2A(lm)*z_Rloc(lmA,nR) -  & ! sin(theta) dtheta z
                           &    dTheta2S(lm)*z_Rloc(lmS,nR) )
                        else if ( l == l_R(nR) ) then
                           CorPol_loc =two*CorFac*or1(nR) * (  &
                           &        dPhi(lm)*dw_Rloc(lm,nR) -  & ! phi-deriv of dw/dr
                           &    dTheta2S(lm)*z_Rloc(lmS,nR) )
                        else
                           CorPol_loc=zero
                        end if
                     else
                        CorPol_loc=zero
                     end if

                     if ( l_conv_nl ) then
                        AdvPol_loc=or2(nR)*this%AdvrLM(lm,nR)
                     else
                        AdvPol_loc=zero
                     endif

                  end if ! Double curl or not for the poloidal equation

                  dwdt(lm,nR)=AdvPol_loc+CorPol_loc

                  if ( l_corr ) then
                     if ( l < l_R(nR) ) then
                        CorTor_loc=          two*CorFac*or2(nR) * (  &
                        &                 dPhi(lm)*z_Rloc(lm,nR)   + &
                        &            dTheta3A(lm)*dw_Rloc(lmA,nR)  + &
                        &    or1(nR)*dTheta4A(lm)* w_Rloc(lmA,nR)  + &
                        &            dTheta3S(lm)*dw_Rloc(lmS,nR)  - &
                        &    or1(nR)*dTheta4S(lm)* w_Rloc(lmS,nR)  )
                     else if ( l == l_R(nR) ) then
                        CorTor_loc=          two*CorFac*or2(nR) * (  &
                        &                 dPhi(lm)*z_Rloc(lm,nR)   + &
                        &            dTheta3S(lm)*dw_Rloc(lmS,nR)  - &
                        &    or1(nR)*dTheta4S(lm)* w_Rloc(lmS,nR)  )
                     else
                        CorTor_loc=zero
                     end if
                  else
                     CorTor_loc=zero
                  end if
                  if ( l_conv_nl ) then
                     AdvTor_loc=dLh(lm)*this%AdvpLM(lm,nR)
                  else
                     AdvTor_loc=zero
                  end if
                  dzdt(lm,nR)=CorTor_loc+AdvTor_loc

                  ! In case double curl is calculated dpdt is useless
                  if ( (.not. l_double_curl) .or. lPressNext ) then
                     !------ Recycle CorPol and AdvPol:
                     if ( l_corr ) then
                        if ( l < l_R(nR) ) then
                           CorPol_loc=           two*CorFac*or2(nR) *  &
                           &           ( -dPhi(lm)  * ( dw_Rloc(lm,nR) &
                           &            +or1(nR)*dLh(lm)*w_Rloc(lm,nR) &
                           &                                         ) &
                           &              +dTheta3A(lm)*z_Rloc(lmA,nR) &
                           &              +dTheta3S(lm)*z_Rloc(lmS,nR) &
                           &           )
                        else if ( l == l_R(nR) ) then
                           CorPol_loc=           two*CorFac*or2(nR) *  &
                           &           ( -dPhi(lm)  * ( dw_Rloc(lm,nR) &
                           &            +or1(nR)*dLh(lm)*w_Rloc(lm,nR) &
                           &                                         ) &
                           &              +dTheta3S(lm)*z_Rloc(lmS,nR) &
                           &           )
                        else
                           CorPol_loc=zero
                        end if
                     else
                        CorPol_loc=zero
                     end if
                     if ( l_conv_nl ) then
                        AdvPol_loc=-dLh(lm)*this%AdvtLM(lm,nR)
                     else
                        AdvPol_loc=zero
                     end if
                     dpdt(lm,nR)=AdvPol_loc+CorPol_loc

                  end if ! lPressNext
               end if ! l > 0
            end do
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#else
         !$omp end parallel do
#endif

      else
#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do collapse(2)
#else
         !$omp parallel do default(shared) private(lm)
#endif
         do nR=nRstart,nRstop
            do lm=1,lm_max
               dwdt(lm,nR)=zero
               dzdt(lm,nR)=zero
               dpdt(lm,nR)=zero
            end do
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#else
         !$omp end parallel do
#endif
      end if

      if ( l_heat ) then

#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do collapse(2) &
         !$omp private(dsdt_loc,l)
#else
         !$omp parallel do default(shared) private(lm,dsdt_loc,l)
#endif
         do nR=nRstart,nRstop
            do lm=1,lm_max
               l=lm2l(lm)
               if ( l == 0 ) then
                  dsdt_loc =epsc*epscProf(nR)!+opr/epsS*divKtemp0(nR)
                  if ( l_anel ) then
                     if ( l_anelastic_liquid ) then
                        dsdt_loc=dsdt_loc+temp0(nR)*this%heatTermsLM(lm,nR)
                     else
                        dsdt_loc=dsdt_loc+this%heatTermsLM(lm,nR)
                     end if
                  end if
                  dsdt(lm,nR)=dsdt_loc
               else
                  dsdt_loc     =dLh(lm)*this%VStLM(lm,nR)
                  if ( l_anel ) then
                     if ( l_anelastic_liquid ) then
                        dsdt_loc=dsdt_loc+temp0(nR)*this%heatTermsLM(lm,nR)
                     else
                        dsdt_loc=dsdt_loc+this%heatTermsLM(lm,nR)
                     end if
                  end if
                  dsdt(lm,nR)=dsdt_loc
               end if
            end do
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#else
         !$omp end parallel do
#endif
      end if

      if ( l_chemical_conv ) then

#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do collapse(2) &
         !$omp private(l)
#else
         !$omp parallel do default(shared) private(l,lm)
#endif
         do nR=nRstart,nRstop
            do lm=1,lm_max
               l   =lm2l(lm)
               if ( l == 0 ) then
                  dxidt(lm,nR)=epscXi
               else
                  dxidt(lm,nR)=dLh(lm)*this%VXitLM(lm,nR)
               end if
            end do
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#else
         !$omp end parallel do
#endif
      end if

      if ( l_mag_nl .or. l_mag_kin  ) then
#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do collapse(2)
#else
         !$omp parallel do default(shared) private(lm)
#endif
         do nR=nRstart,nRstop
            do lm=1,lm_max
               dbdt(lm,nR)   = dLh(lm)*this%VxBpLM(lm,nR)
               dVxBhLM(lm,nR)=-dLh(lm)*this%VxBtLM(lm,nR)*r(nR)*r(nR)
               djdt(lm,nR)   = dLh(lm)*or4(nR)*this%VxBrLM(lm,nR)
            end do
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#else
         !$omp end parallel do
#endif
      else
         if ( l_mag ) then
#ifdef WITH_OMP_GPU
            !$omp target teams distribute parallel do collapse(2)
#else
            !$omp parallel do
#endif
            do nR=nRstart,nRstop
               do lm=1,lm_max
                  dbdt(lm,nR)   =zero
                  djdt(lm,nR)   =zero
                  dVxBhLM(lm,nR)=zero
               end do
            end do
#ifdef WITH_OMP_GPU
            !$omp end target teams distribute parallel do
#else
            !$omp end parallel do
#endif
         end if
      end if

      ! boundary values are set to zero except for FD and single matrix with Tdiff
      if ( (.not. l_parallel_solve) .and.  &
      &    (.not. (l_single_matrix .and. l_temperature_diff)) ) then

         do nR=nRstart, nRstop
            if (nR==n_r_cmb .or. nR==n_r_icb ) then

               if ( l_mag_nl .or. l_mag_kin ) then
#ifdef WITH_OMP_GPU
                  !$omp target teams distribute parallel do private(lm,l)
#else
                  !$omp parallel do default(shared) private(lm,l)
#endif
                  do lm=1,lm_max
                     l   =lm2l(lm)
                     if ( l == 0 ) then
                        dVxBhLM(lm,nR)=zero
                     else
                        dVxBhLM(lm,nR)=-dLh(lm)*this%VxBtLM(lm,nR)*r(nR)*r(nR)
                     end if
                  end do
#ifdef WITH_OMP_GPU
                  !$omp end target teams distribute parallel do
#else
                  !$omp end parallel do
#endif
               else
                  if ( l_mag ) then
#ifdef WITH_OMP_GPU
                     !$omp target teams distribute parallel do
#else
                     !$omp parallel do
#endif
                     do lm=1,lm_max
                        dVxBhLM(lm,nR)=zero
                     end do
#ifdef WITH_OMP_GPU
                     !$omp end target teams distribute parallel do
#else
                     !$omp end parallel do
#endif
                  end if
               end if
               if ( l_double_curl ) then
#ifdef WITH_OMP_GPU
                  !$omp target teams distribute parallel do
#else
                  !$omp parallel do
#endif
                  do lm=1,lm_max
                     dVxVhLM(lm,nR)=zero
                  end do
#ifdef WITH_OMP_GPU
                  !$omp end target teams distribute parallel do
#else
                  !$omp end parallel do
#endif
               end if
            end if  ! boundary ? lvelo ?
         end do
      end if ! Do the boundary points need special care?

   end subroutine get_td
!----------------------------------------------------------------------------
   subroutine get_dwdt(this,dwdt)
      !
      ! This subroutine finishes the assembly of the explicit terms that
      ! enter the equation for the poloidal equation for nRstart:nRstop.
      !

      !-- Input of variables:
      class(nonlinear_lm_2d_t) :: this

      !-- Output of variables:
      complex(cp), intent(out) :: dwdt(lm_max,nRstart:nRstop)

      !-- Local variables:
      integer :: l,lm,lmS,lmA,nR
      complex(cp) :: CorPol_loc

#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do collapse(2) private(l,lmS,lmA) &
      !$omp private(CorPol_loc)
#else
      !$omp parallel do default(shared) private(l,lmS,lmA) &
      !$omp private(CorPol_loc)
#endif
      do nR=nRstart,nRstop
         do lm=1,lm_max
            l  =lm2l(lm)
            lmS=lm2lmS(lm)
            lmA=lm2lmA(lm)

            if ( l == 0 ) then ! This is l=0,m=0
               if ( l_conv_nl ) then
                  dwdt(lm,nR)=or2(nR)*this%AdvrLM(lm,nR)
               else
                  dwdt(lm,nR)=zero
               end if
               if ( l_corr .and. (.not. l_single_matrix) ) then
                  dwdt(lm,nR)=dwdt(lm,nR)+two*CorFac*or1(nR)*dTheta2A(lm)*z_Rloc(lmA,nR)
               end if
            else
               if ( l_conv_nl ) then
                  dwdt(lm,nR)=or2(nR)*this%AdvrLM(lm,nR)
               else
                  dwdt(lm,nR)=zero
               endif

               if ( l_corr ) then
                  if ( l < l_R(nR) ) then
                     CorPol_loc =two*CorFac*or1(nR) * (  &
                     &        dPhi(lm)*dw_Rloc(lm,nR) +  & ! phi-deriv of dw/dr
                     &    dTheta2A(lm)*z_Rloc(lmA,nR) -  & ! sin(theta) dtheta z
                     &    dTheta2S(lm)*z_Rloc(lmS,nR) )
                  else if ( l == l_R(nR) ) then
                     CorPol_loc =two*CorFac*or1(nR) * (  &
                     &        dPhi(lm)*dw_Rloc(lm,nR) -  & ! phi-deriv of dw/dr
                     &    dTheta2S(lm)*z_Rloc(lmS,nR) )
                  else
                     CorPol_loc=zero
                  end if
                  dwdt(lm,nR)=dwdt(lm,nR)+CorPol_loc
               end if
            end if
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#else
      !$omp end parallel do
#endif

   end subroutine get_dwdt
!----------------------------------------------------------------------------
   subroutine get_dwdt_double_curl(this,dwdt,dVxVhLM)
      !
      ! This subroutine finishes the assembly of the explicit terms that
      ! enter the equation for the poloidal equation in case the double
      ! curl formulation is employed.
      !

      !-- Input of variables:
      class(nonlinear_lm_2d_t) :: this

      !-- Output of variables:
      complex(cp), intent(out) :: dwdt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dVxVhLM(lm_max,nRstart:nRstop)

      !-- Local variables:
      integer :: l,lm,lmS,lmA,nR
      complex(cp) :: CorPol_loc

#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do collapse(2) private(l,lmS,lmA) &
      !$omp private(CorPol_loc)
#else
      !$omp parallel do default(shared) private(l,lmS,lmA) &
      !$omp private(CorPol_loc)
#endif
      do nR=nRstart,nRstop
         do lm=1,lm_max
            l  =lm2l(lm)
            lmS=lm2lmS(lm)
            lmA=lm2lmA(lm)

            if ( l == 0 ) then
               if ( l_conv_nl ) then
                  dwdt(lm,nR)=or2(nR)*this%AdvrLM(lm,nR)
               else
                  dwdt(lm,nR)=zero
               end if
               if ( l_corr .and. (.not. l_single_matrix ) ) then
                  dwdt(lm,nR)=dwdt(lm,nR)+two*CorFac*or1(nR)*dTheta2A(lm)*z_Rloc(lmA,nR)
               end if
            else
               if ( l_conv_nl ) then
                  dwdt(lm,nR)   =dLh(lm)*or4(nR)*orho1(nR)*this%AdvrLM(lm,nR)
                  dVxVhLM(lm,nR)=-orho1(nR)*r(nR)*r(nR)*dLh(lm)*this%AdvtLM(lm,nR)
               else
                  dwdt(lm,nR)   =zero
                  dVxVhLM(lm,nR)=zero
               end if

               if ( l_corr ) then
                  if ( l < l_R(nR) ) then
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
                  else if ( l == l_R(nR) ) then
                     CorPol_loc =two*CorFac*or2(nR)*orho1(nR)*(               &
                     &                    dPhi(lm)*(                          &
                     &         -ddw_Rloc(lm,nR)+beta(nR)*dw_Rloc(lm,nR)     + &
                     &             ( beta(nR)*or1(nR)+or2(nR))*               &
                     &                            dLh(lm)*w_Rloc(lm,nR) )   + &
                     &             dTheta3S(lm)*( dz_Rloc(lmS,nR)-            &
                     &                            beta(nR)*z_Rloc(lmS,nR) ) - &
                     &          or1(nR)* dTheta4S(lm)* z_Rloc(lmS,nR) )
                  else
                     CorPol_loc=zero
                  end if
                  dwdt(lm,nR)=dwdt(lm,nR)+CorPol_loc
               end if
            end if
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#else
      !$omp end parallel do
#endif

      ! boundary values are set to zero except for FD and single matrix with Tdiff
      if ( (.not. l_parallel_solve) .and.  &
      &    (.not. (l_single_matrix .and. l_temperature_diff)) ) then

         do nR=nRstart,nRstop
            if (nR==n_r_cmb .or. nR==n_r_icb ) then

#ifdef WITH_OMP_GPU
               !$omp target teams distribute parallel do
#else
               !$omp parallel do
#endif
               do lm=1,lm_max
                  dVxVhLM(lm,nR)=zero
               end do
#ifdef WITH_OMP_GPU
               !$omp end target teams distribute parallel do
#else
               !$omp end parallel do
#endif
            end if  ! boundary ? lvelo ?
         end do
      end if ! Do the boundary points need special care?

   end subroutine get_dwdt_double_curl
!----------------------------------------------------------------------------
   subroutine get_dpdt(this,dpdt)
      !
      ! This subroutine finishes the assembly of the explicit terms that
      ! enter the equation for pressure dpdt for nRstart:nRstop
      !

      !-- Input of variables:
      class(nonlinear_lm_2d_t) :: this

      !-- Output of variables:
      complex(cp), intent(out) :: dpdt(lm_max,nRstart:nRstop)

      !-- Local variables:
      integer :: l,lm,lmS,lmA,nR
      complex(cp) :: CorPol_loc

#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do collapse(2) private(l,lmS,lmA) &
         !$omp private(CorPol_loc)
#else
         !$omp parallel do default(shared) private(l,lmS,lmA) &
         !$omp private(CorPol_loc)
#endif
      do nR=nRstart,nRstop
         do lm=1,lm_max
            l  =lm2l(lm)
            lmS=lm2lmS(lm)
            lmA=lm2lmA(lm)

            if ( l > 0 ) then
               if ( l_conv_nl ) then
                  dpdt(lm,nR)=-dLh(lm)*this%AdvtLM(lm,nR)
               else
                  dpdt(lm,nR)=zero
               end if

               if ( l_corr ) then
                  if ( l < l_R(nR) ) then
                     CorPol_loc=           two*CorFac*or2(nR) *  &
                     &           ( -dPhi(lm)  * ( dw_Rloc(lm,nR) &
                     &            +or1(nR)*dLh(lm)*w_Rloc(lm,nR) &
                     &                                         ) &
                     &              +dTheta3A(lm)*z_Rloc(lmA,nR) &
                     &              +dTheta3S(lm)*z_Rloc(lmS,nR) &
                     &           )
                  else if ( l == l_R(nR) ) then
                     CorPol_loc=           two*CorFac*or2(nR) *  &
                     &           ( -dPhi(lm)  * ( dw_Rloc(lm,nR) &
                     &            +or1(nR)*dLh(lm)*w_Rloc(lm,nR) &
                     &                                         ) &
                     &              +dTheta3S(lm)*z_Rloc(lmS,nR) &
                     &           )
                  else
                     CorPol_loc=zero
                  end if
                  dpdt(lm,nR)=dpdt(lm,nR)+CorPol_loc
               end if
            end if
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#else
      !$omp end parallel do
#endif

   end subroutine get_dpdt
!-----------------------------------------------------------------------------
   subroutine get_dzdt(this,dzdt)
      !
      ! This subroutine finishes the assembly of the explicit terms that
      ! enter the toroidal equation dzdt for the range nRstart:nRstop.
      !

      !-- Input of variables:
      class(nonlinear_lm_2d_t) :: this

      !-- Output of variables:
      complex(cp), intent(out) :: dzdt(lm_max,nRstart:nRstop)

      !-- Local variables:
      integer :: l,lm,lmS,lmA,nR
      complex(cp) :: CorTor_loc

#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do collapse(2) private(l,lmS,lmA) &
      !$omp private(CorTor_loc)
#else
      !$omp parallel do default(shared) private(l,lmS,lmA) &
      !$omp private(CorTor_loc)
#endif
      do nR=nRstart,nRstop
         do lm=1,lm_max
            l  =lm2l(lm)
            lmS=lm2lmS(lm)
            lmA=lm2lmA(lm)
            
            if ( l == 0 ) then
               dzdt(lm,nR)=zero!-dTheta1A(lm)*this%AdvpLM(lmA)
               if ( l_corr ) then
                  dzdt(lm,nR)=dzdt(lm,nR)+ two*CorFac*or2(nR) * (    &
                  &                dTheta3A(lm)*dw_Rloc(lmA,nR) +    &
                  &        or1(nR)*dTheta4A(lm)* w_Rloc(lmA,nR) )
               end if
            else
               if ( l_conv_nl ) then
                  dzdt(lm,nR)=dLh(lm)*this%AdvpLM(lm,nR)
               else
                  dzdt(lm,nR)=zero
               end if

               if ( l_corr ) then
                  if ( l < l_R(nR) ) then
                     CorTor_loc=          two*CorFac*or2(nR) * (  &
                     &                 dPhi(lm)*z_Rloc(lm,nR)   + &
                     &            dTheta3A(lm)*dw_Rloc(lmA,nR)  + &
                     &    or1(nR)*dTheta4A(lm)* w_Rloc(lmA,nR)  + &
                     &            dTheta3S(lm)*dw_Rloc(lmS,nR)  - &
                     &    or1(nR)*dTheta4S(lm)* w_Rloc(lmS,nR)  )
                  else if ( l == l_R(nR) ) then
                     CorTor_loc=          two*CorFac*or2(nR) * (  &
                     &                 dPhi(lm)*z_Rloc(lm,nR)   + &
                     &            dTheta3S(lm)*dw_Rloc(lmS,nR)  - &
                     &    or1(nR)*dTheta4S(lm)* w_Rloc(lmS,nR)  )
                  else
                     CorTor_loc=zero
                  end if
                  dzdt(lm,nR)=dzdt(lm,nR)+CorTor_loc
               end if
            end if
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#else
      !$omp end parallel do
#endif

   end subroutine get_dzdt
!-----------------------------------------------------------------------------
   subroutine get_dsdt(this, dsdt, dVSrLM)
      !
      ! This subroutine finishes the assembly of dsdt and dVSrLM
      ! for the whole range of nRstart:nRstop
      !

      class(nonlinear_lm_2d_t) :: this

      !-- Output variables
      complex(cp), intent(out) :: dsdt(lm_max,nRstart:nRstop) ! divH(uh*s)
      complex(cp), intent(out) :: dVSrLM(lm_max,nRstart:nRstop) ! ur*s

      !-- Local variables
      integer :: lm, l, nR

      if ( l_anel ) then
         if ( l_anelastic_liquid ) then
#ifdef WITH_OMP_GPU
            !$omp target teams distribute parallel do collapse(2) private(l)
#else
            !$omp parallel do private(l)
#endif
            do nR=nRstart,nRstop
               do lm=1,lm_max
                  l=lm2l(lm)
                  if ( l == 0 ) then
                     dsdt(lm,nR)=epsc*epscProf(nR)+temp0(nR)*this%heatTermsLM(lm,nR)
                  else
                     dsdt(lm,nR)=dLh(lm)*this%VStLM(lm,nR) + &
                     &           temp0(nR)*this%heatTermsLM(lm,nR)
                  end if
               end do
            end do
#ifdef WITH_OMP_GPU
            !$omp end target teams distribute parallel do
#else
            !$omp end parallel do
#endif
         else
#ifdef WITH_OMP_GPU
            !$omp target teams distribute parallel do collapse(2) private(l)
#else
            !$omp parallel do private(l)
#endif
            do nR=nRstart,nRstop
               do lm=1,lm_max
                  l=lm2l(lm)
                  if ( l == 0 ) then
                     dsdt(lm,nR)=epsc*epscProf(nR)+this%heatTermsLM(lm,nR)
                  else
                     dsdt(lm,nR)=dLh(lm)*this%VStLM(lm,nR)+this%heatTermsLM(lm,nR)
                  end if
               end do
            end do
#ifdef WITH_OMP_GPU
            !$omp end target teams distribute parallel do
#else
            !$omp end parallel do
#endif
         end if

      else ! Boussinesq

#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do collapse(2) private(l)
#else
         !$omp parallel do private(l)
#endif
         do nR=nRstart,nRstop
            do lm=1,lm_max
               l=lm2l(lm)
               if ( l == 0 ) then
                  dsdt(lm,nR)=epsc*epscProf(nR)
               else
                  dsdt(lm,nR)=dLh(lm)*this%VStLM(lm,nR)
               end if
            end do
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#else
         !$omp end parallel do
#endif
      end if

   end subroutine get_dsdt
!-----------------------------------------------------------------------------
   subroutine get_dxidt(this, dxidt, dVXirLM)
      !
      ! This subroutine finishes the assembly of dxidt and dVXirLM
      ! for the whole range nRstart:nRstop
      !

      class(nonlinear_lm_2d_t) :: this

      !-- Output variables
      complex(cp), intent(out) :: dxidt(lm_max,nRstart:nRstop) ! divH(uh*xi)
      complex(cp), intent(out) :: dVXirLM(lm_max,nRstart:nRstop) ! ur*xi

      !-- Local variables
      integer :: lm, l, nR

#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do collapse(2)private(l)
#else
      !$omp parallel do private(l)
#endif
      do nR=nRstart,nRstop
         do lm=1,lm_max
            l=lm2l(lm)
            if ( l == 0 ) then
               dxidt(lm,nR)=epscXi
            else
               dxidt(lm,nR)=dLh(lm)*this%VXitLM(lm,nR)
            end if
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#else
      !$omp end parallel do
#endif

   end subroutine get_dxidt
!-----------------------------------------------------------------------------
   subroutine get_dbdt(this, dbdt, djdt, dVxBhLM)
      !
      ! This subroutine finishes the assembly of dbdt, djdt and dVxBhLM
      ! at the range of radial levels nRstart:nRstop
      !

      class(nonlinear_lm_2d_t) :: this

      !-- Output variables
      complex(cp), intent(out) :: dbdt(lm_maxMag,nRstartMag:nRstopMag)
      complex(cp), intent(out) :: djdt(lm_maxMag,nRstartMag:nRstopMag)
      complex(cp), intent(out) :: dVxBhLM(lm_maxMag,nRstartMag:nRstopMag)

      !-- Local variables
      integer :: lm,l,nR

      if ( l_mag_nl .or. l_mag_kin  ) then
#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do collapse(2)
#else
         !$omp parallel do default(shared)
#endif
         do nR=nRstart,nRstop
            do lm=1,lm_max
               dbdt(lm,nR)   = dLh(lm)*this%VxBpLM(lm,nR)
               dVxBhLM(lm,nR)=-dLh(lm)*this%VxBtLM(lm,nR)*r(nR)*r(nR)
               djdt(lm,nR)   = dLh(lm)*or4(nR)*this%VxBrLM(lm,nR)
            end do
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#else
         !$omp end parallel do
#endif
      else
         if ( l_mag ) then
#ifdef WITH_OMP_GPU
            !$omp target teams distribute parallel do collapse(2)
#else
            !$omp parallel do
#endif
            do nR=nRstart,nRstop
               do lm=1,lm_max
                  dbdt(lm,nR)   =zero
                  djdt(lm,nR)   =zero
                  dVxBhLM(lm,nR)=zero
               end do
            end do
#ifdef WITH_OMP_GPU
            !$omp end target teams distribute parallel do
#else
            !$omp end parallel do
#endif
         end if
      end if

      ! boundary values are set to zero except for FD and single matrix with Tdiff
      if ( (.not. l_parallel_solve) .and.  &
      &    (.not. (l_single_matrix .and. l_temperature_diff)) ) then

         do nR=nRstart,nRstop
            if (nR==n_r_cmb .or. nR==n_r_icb ) then

               if ( l_mag_nl .or. l_mag_kin ) then
#ifdef WITH_OMP_GPU
                  !$omp target teams distribute parallel do private(lm,l)
#else
                  !$omp parallel do default(shared) private(lm,l)
#endif
                  do lm=1,lm_max
                     l   =lm2l(lm)
                     if ( l == 0 ) then
                        dVxBhLM(lm,nR)=zero
                     else
                        dVxBhLM(lm,nR)=-dLh(lm)*this%VxBtLM(lm,nR)*r(nR)*r(nR)
                     end if
                  end do
#ifdef WITH_OMP_GPU
                  !$omp end target teams distribute parallel do
#else
                  !$omp end parallel do
#endif
               else
                  if ( l_mag ) then
#ifdef WITH_OMP_GPU
                     !$omp target teams distribute parallel do
#else
                     !$omp parallel do
#endif
                     do lm=1,lm_max
                        dVxBhLM(lm,nR)=zero
                     end do
#ifdef WITH_OMP_GPU
                     !$omp end target teams distribute parallel do
#else
                     !$omp end parallel do
#endif
                  end if
               end if
            end if  ! boundary ? lvelo ?
         end do
      end if ! Do the boundary points need special care?

   end subroutine get_dbdt
!-----------------------------------------------------------------------------
end module nonlinear_lm_2d_mod
