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
   use truncation, only: lm_max, l_max, lm_maxMag, lmP_max
   use logic, only : l_anel, l_conv_nl, l_corr, l_heat, l_anelastic_liquid, &
       &             l_mag_nl, l_mag_kin, l_mag_LF, l_conv, l_mag,          &
       &             l_chemical_conv, l_single_matrix, l_double_curl,       &
       &             l_adv_curl, l_phase_field
   use radial_functions, only: r, or2, or1, beta, epscProf, or4, temp0, orho1
   use physical_parameters, only: CorFac, epsc,  n_r_LCR, epscXi
   use blocking, only: lm2l, lm2m, lm2lmP, lm2lmA, lm2lmS
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
      complex(cp), allocatable :: VSrLM(:),  VStLM(:),  VSpLM(:)
      complex(cp), allocatable :: VXirLM(:),  VXitLM(:), VXipLM(:)
      complex(cp), allocatable :: heatTermsLM(:), dphidtLM(:)
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: set_zero
      procedure :: get_td
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
         allocate(this%VSrLM(lmP_max),this%VStLM(lmP_max),this%VSpLM(lmP_max))
         bytes_allocated = bytes_allocated + 3*lmP_max*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         gpu_bytes_allocated = gpu_bytes_allocated + 3*lmP_max*SIZEOF_DEF_COMPLEX
#endif
      else
         allocate( this%VSrLM(1), this%VStLM(1),this%VSpLM(1) )
      end if

      if ( l_chemical_conv ) then
         allocate(this%VXirLM(lmP_max),this%VXitLM(lmP_max),this%VXipLM(lmP_max))
         bytes_allocated = bytes_allocated + 3*lmP_max*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         gpu_bytes_allocated = gpu_bytes_allocated + 3*lmP_max*SIZEOF_DEF_COMPLEX
#endif
      else
         allocate(this%VXirLM(1),this%VXitLM(1),this%VXipLM(1))
      end if

      if ( l_phase_field ) then
         allocate(this%dphidtLM(lmP_max))
         bytes_allocated = bytes_allocated + lmP_max*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         gpu_bytes_allocated = gpu_bytes_allocated + lmP_max*SIZEOF_DEF_COMPLEX
#endif
      else
         allocate(this%dphidtLM(1))
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
      deallocate( this%heatTermsLM, this%VXirLM, this%VXitLM, this%VXipLM )
      deallocate( this%VSrLM, this%VStLM, this%VSpLM, this%dphidtLM )

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
      do lm=1,lmP_max
         this%AdvrLM(lm)=zero
         this%AdvtLM(lm)=zero
         this%AdvpLM(lm)=zero
         if ( l_mag ) then
            this%VxBrLM(lm)=zero
            this%VxBtLM(lm)=zero
            this%VxBpLM(lm)=zero
         end if
         if ( l_heat ) then
            this%VSrLM(lm)=zero
            this%VStLM(lm)=zero
            this%VSpLM(lm)=zero
         end if
         if ( l_anel ) this%heatTermsLM(lm)=zero
         if ( l_phase_field ) this%dphidtLM(lm)=zero
         if ( l_chemical_conv ) then
            this%VXirLM(lm)=zero
            this%VXitLM(lm)=zero
            this%VXipLM(lm)=zero
         end if
      end do
      !$omp end parallel do

#ifdef WITH_OMP_GPU
      !$omp target update to(this)
#endif

   end subroutine set_zero
!----------------------------------------------------------------------------
   subroutine get_td(this,nR,nBc,lPressNext,dVSrLM,dVXirLM,dVxVhLM,dVxBhLM, &
              &      dwdt,dzdt,dpdt,dsdt,dxidt,dphidt,dbdt,djdt)
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
      logical, intent(in) :: lPressNext

      !-- Output of variables:
      complex(cp), intent(out) :: dwdt(:),dzdt(:)
      complex(cp), intent(out) :: dpdt(:),dsdt(:)
      complex(cp), intent(out) :: dxidt(:),dphidt(:)
      complex(cp), intent(out) :: dbdt(:),djdt(:)
      complex(cp), intent(out) :: dVxBhLM(:)
      complex(cp), intent(out) :: dVxVhLM(:)
      complex(cp), intent(out) :: dVSrLM(:)
      complex(cp), intent(out) :: dVXirLM(:)

      !-- Local variables:
      integer :: l,m,lm,lmS,lmA,lmP
      complex(cp) :: AdvPol_loc,CorPol_loc,AdvTor_loc,CorTor_loc,dsdt_loc
      !integer, parameter :: DOUBLE_COMPLEX_PER_CACHELINE=4

      !write(*,"(I3,A,4ES20.12)") nR,": get_td start: ",SUM(this%AdvrLM)

      !lm_chunksize=(((lm_max)/nThreads)/DOUBLE_COMPLEX_PER_CACHELINE) * &
      !             & DOUBLE_COMPLEX_PER_CACHELINE
      !lm_chunksize=4
      !write(*,"(A,I4)") "Using a chunksize of ",lm_chunksize

      if (nBc == 0 ) then

         if ( l_conv ) then  ! Convection

#ifdef WITH_OMP_GPU
            !$omp target teams distribute parallel do private(l,m,lmS,lmA,lmP) &
            !$omp& private(AdvPol_loc,CorPol_loc,AdvTor_loc,CorTor_loc)
#else
            !$omp parallel do default(shared) private(lm,l,m,lmS,lmA,lmP) &
            !$omp private(AdvPol_loc,CorPol_loc,AdvTor_loc,CorTor_loc)
#endif
            do lm=1,lm_max
               l   =lm2l(lm)
               m   =lm2m(lm)
               lmS =lm2lmS(lm)
               lmA =lm2lmA(lm)
               lmP =lm2lmP(lm)

               if ( l == 0 ) then ! This is l=0,m=0
                  lmA=lm2lmA(1)
                  lmP=1
                  !lmPA=lmP2lmPA(lmP)
                  if ( l_conv_nl ) then
                     AdvPol_loc=or2(nR)*this%AdvrLM(lmP)
                     AdvTor_loc=zero!-dTheta1A(lm)*this%AdvpLM(lmPA)
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
               else

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
                           CorPol_loc=zero
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
                        AdvPol_loc =dLh(lm)*or4(nR)*orho1(nR)*this%AdvrLM(lmP)
                        dVxVhLM(lm)=-orho1(nR)*r(nR)*r(nR)*dLh(lm)*this%AdvtLM(lmP)
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
                           CorPol_loc=zero
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

                  if ( l_corr ) then
                     if ( l < l_max .and. l > m ) then
                        CorTor_loc=          two*CorFac*or2(nR) * (  &
                        &                 dPhi(lm)*z_Rloc(lm,nR)   + &
                        &            dTheta3A(lm)*dw_Rloc(lmA,nR)  + &
                        &    or1(nR)*dTheta4A(lm)* w_Rloc(lmA,nR)  + &
                        &            dTheta3S(lm)*dw_Rloc(lmS,nR)  - &
                        &    or1(nR)*dTheta4S(lm)* w_Rloc(lmS,nR)  )
                     else if ( l == l_max ) then
                        CorTor_loc=zero
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
                     AdvTor_loc=dLh(lm)*this%AdvpLM(lmP)
                  else
                     AdvTor_loc=zero
                  end if
                  dzdt(lm)=CorTor_loc+AdvTor_loc
               end if

            end do
#ifdef WITH_OMP_GPU
            !$omp end target teams distribute parallel do
#else
            !$omp end parallel do
#endif

            ! In case double curl is calculated dpdt is useless
            if ( (.not. l_double_curl) .or. lPressNext ) then
            !if ( .true. ) then
#ifdef WITH_OMP_GPU
               !$omp target teams distribute parallel do private(lm,l,m,lmS,lmA,lmP) &
               !$omp& private(AdvPol_loc,CorPol_loc)
#else
               !$omp parallel do default(shared) private(lm,l,m,lmS,lmA,lmP) &
               !$omp private(AdvPol_loc,CorPol_loc)
#endif
               do lm=1,lm_max
                  l   =lm2l(lm)
                  if ( l == 0 ) cycle
                  m   =lm2m(lm)
                  lmS =lm2lmS(lm)
                  lmA =lm2lmA(lm)
                  lmP =lm2lmP(lm)

                  !------ Recycle CorPol and AdvPol:
                  if ( l_corr ) then
                     if ( l < l_max .and. l > m ) then
                        CorPol_loc=           two*CorFac*or2(nR) *  &
                        &           ( -dPhi(lm)  * ( dw_Rloc(lm,nR) &
                        &            +or1(nR)*dLh(lm)*w_Rloc(lm,nR) &
                        &                                         ) &
                        &              +dTheta3A(lm)*z_Rloc(lmA,nR) &
                        &              +dTheta3S(lm)*z_Rloc(lmS,nR) &
                        &           )

                     else if ( l == l_max ) then
                        CorPol_loc=zero

                     else if ( l == m ) then
                        CorPol_loc=                    two*CorFac*or2(nR) *  &
                        &                    ( -dPhi(lm)  * ( dw_Rloc(lm,nR) &
                        &                     +or1(nR)*dLh(lm)*w_Rloc(lm,nR) &
                        &                                                   )&
                        &                      +dTheta3A(lm)*z_Rloc(lmA,nR)  &
                        &                    )

                     end if
                  else
                     CorPol_loc=zero
                  end if
                  if ( l_conv_nl ) then
                     AdvPol_loc=-dLh(lm)*this%AdvtLM(lmP)
                  else
                     AdvPol_loc=zero
                  end if
                  dpdt(lm)=AdvPol_loc+CorPol_loc

               end do ! lm loop
#ifdef WITH_OMP_GPU
               !$omp end target teams distribute parallel do
#else
               !$omp end parallel do
#endif
            end if

         else
#ifdef WITH_OMP_GPU
            !$omp target teams distribute parallel do
#endif
            do lm=1,lm_max
               dwdt(lm)=zero
               dzdt(lm)=zero
               dpdt(lm)=zero
            end do
#ifdef WITH_OMP_GPU
            !$omp end target teams distribute parallel do
#endif
         end if ! l_conv ?

      end if

      if ( nBc == 0 ) then

         if ( l_heat ) then

#ifdef WITH_OMP_GPU
            !$omp target teams distribute parallel do private(lmP,dsdt_loc,l)
#else
            !$omp parallel do default(shared) private(lm,lmP,dsdt_loc,l)
#endif
            do lm=1,lm_max
               l   =lm2l(lm)
               lmP =lm2lmP(lm)
               if ( l == 0 ) then
                  dsdt_loc =epsc*epscProf(nR)!+opr/epsS*divKtemp0(nR)
                  dVSrLM(1)=this%VSrLM(1)
                  if ( l_anel ) then
                     if ( l_anelastic_liquid ) then
                        dsdt_loc=dsdt_loc+temp0(nR)*this%heatTermsLM(1)
                     else
                        dsdt_loc=dsdt_loc+this%heatTermsLM(1)
                     end if
                  end if
                  dsdt(1)=dsdt_loc
               else
                  dVSrLM(lm)=this%VSrLM(lmP)
                  dsdt_loc  =dLh(lm)*this%VStLM(lmP)
                  if ( l_anel ) then
                     if ( l_anelastic_liquid ) then
                        dsdt_loc = dsdt_loc+temp0(nR)*this%heatTermsLM(lmP)
                     else
                        dsdt_loc = dsdt_loc+this%heatTermsLM(lmP)
                     end if
                  end if
                  dsdt(lm) = dsdt_loc
               end if
            end do
#ifdef WITH_OMP_GPU
            !$omp end target teams distribute parallel do
#else
            !$omp end parallel do
#endif
         end if

         if ( l_chemical_conv ) then

#ifdef WITH_OMP_GPU
            !$omp target teams distribute parallel do private(lmP,l)
#else
            !$omp parallel do default(shared) private(lm,lmP,l)
#endif
            do lm=1,lm_max
               l   =lm2l(lm)
               lmP =lm2lmP(lm)
               if ( l == 0 ) then
                  dVXirLM(1)=this%VXirLM(1)
                  dxidt(1)  =epscXi
               else
                  dVXirLM(lm)=this%VXirLM(lmP)
                  dxidt(lm)  =dLh(lm)*this%VXitLM(lmP)
               end if
            end do
#ifdef WITH_OMP_GPU
            !$omp end target teams distribute parallel do
#else
            !$omp end parallel do
#endif
         end if

         if ( l_phase_field ) then
#ifdef WITH_OMP_GPU
            !$omp target teams distribute parallel do private(lmP)
#else
            !$omp parallel do default(shared) private(lm,lmP)
#endif
            do lm=1,lm_max
               lmP=lm2lmP(lm)
               dphidt(lm)=this%dphidtLM(lmP)
            end do
#ifdef WITH_OMP_GPU
            !$omp end target teams distribute parallel do
#else
            !$omp end parallel do
#endif
         end if

         if ( l_mag_nl .or. l_mag_kin  ) then
#ifdef WITH_OMP_GPU
            !$omp target teams distribute parallel do private(lmP)
#else
            !$omp parallel do default(shared) private(lm,lmP)
#endif
            do lm=1,lm_max
               lmP =lm2lmP(lm)
               dbdt(lm)   = dLh(lm)*this%VxBpLM(lmP)
               dVxBhLM(lm)=-dLh(lm)*this%VxBtLM(lmP)*r(nR)*r(nR)
               djdt(lm)   = dLh(lm)*or4(nR)*this%VxBrLM(lmP)
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
            !$omp target teams distribute parallel do private(lmP,l)
#else
            !$omp parallel do default(shared) private(lm,lmP,l)
#endif
            do lm=1,lm_max
               l   =lm2l(lm)
               lmP =lm2lmP(lm)
               if ( l == 0 ) then
                  dVxBhLM(1)=zero
                  dVSrLM(1) =zero
               else
                  dVxBhLM(lm)=-dLh(lm)*this%VxBtLM(lmP)*r(nR)*r(nR)
                  dVSrLM(lm) =zero
               end if
            end do
#ifdef WITH_OMP_GPU
            !$omp end target teams distribute parallel do
#else
            !$omp end parallel do
#endif
         else
#ifdef WITH_OMP_GPU
            !$omp target teams distribute parallel do
#else
            !$omp parallel do
#endif
            do lm=1,lm_max
               if ( l_mag ) dVxBhLM(lm)=zero
               dVSrLM(lm) =zero
            end do
#ifdef WITH_OMP_GPU
            !$omp end target teams distribute parallel do
#else
            !$omp end parallel do
#endif
         end if
         if ( l_double_curl ) then
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
#endif
         end if
         if ( l_chemical_conv ) then
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
         end if
      end if  ! boundary ? lvelo ?

   end subroutine get_td
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
   use truncation, only: lm_max, l_max, lm_maxMag, lmP_max
   use logic, only : l_anel, l_conv_nl, l_corr, l_heat, l_anelastic_liquid, &
       &             l_mag_nl, l_mag_kin, l_mag_LF, l_conv, l_mag,          &
       &             l_chemical_conv, l_single_matrix, l_double_curl,       &
       &             l_adv_curl, l_phase_field, l_parallel_solve,           &
       &             l_temperature_diff
   use radial_functions, only: r, or2, or1, beta, epscProf, or4, temp0, orho1
   use physical_parameters, only: CorFac, epsc,  n_r_LCR, epscXi
   use blocking, only: lm2l, lm2m, lm2lmP, lm2lmA, lm2lmS
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
      complex(cp), allocatable :: VSrLM(:,:),  VStLM(:,:),  VSpLM(:,:)
      complex(cp), allocatable :: VXirLM(:,:),  VXitLM(:,:), VXipLM(:,:)
      complex(cp), allocatable :: heatTermsLM(:,:), dphidtLM(:,:)
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: set_zero
      procedure :: get_td
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
         allocate( this%VSrLM(lmP_max,nRstart:nRstop) )
         allocate( this%VStLM(lmP_max,nRstart:nRstop) )
         allocate( this%VSpLM(lmP_max,nRstart:nRstop) )
         bytes_allocated = bytes_allocated + 3*lmP_max*(nRstop-nRstart+1) * &
         &                 SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         gpu_bytes_allocated = gpu_bytes_allocated + 3*lmP_max*(nRstop-nRstart+1)* &
         &                     SIZEOF_DEF_COMPLEX
#endif
      else
         allocate( this%VSrLM(1,1), this%VStLM(1,1),this%VSpLM(1,1) )
      end if

      if ( l_chemical_conv ) then
         allocate( this%VXirLM(lmP_max,nRstart:nRstop) )
         allocate( this%VXitLM(lmP_max,nRstart:nRstop) )
         allocate( this%VXipLM(lmP_max,nRstart:nRstop) )
         bytes_allocated = bytes_allocated + 3*lmP_max*(nRstop-nRstart+1) * &
         &                 SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         gpu_bytes_allocated = gpu_bytes_allocated + 3*lmP_max*(nRstop-nRstart+1)* &
         &                     SIZEOF_DEF_COMPLEX
#endif
      else
         allocate(this%VXirLM(1,1),this%VXitLM(1,1),this%VXipLM(1,1))
      end if

      if ( l_phase_field ) then
         allocate(this%dphidtLM(lmP_max,nRstart:nRstop))
         bytes_allocated = bytes_allocated + lmP_max*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         gpu_bytes_allocated = gpu_bytes_allocated + lmP_max*SIZEOF_DEF_COMPLEX
#endif
      else
         allocate(this%dphidtLM(1,1))
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
      deallocate( this%heatTermsLM, this%VXirLM, this%VXitLM, this%VXipLM )
      deallocate( this%VSrLM, this%VStLM, this%VSpLM, this%dphidtLM )

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
         do lm=1,lmP_max
            this%AdvrLM(lm,nR)=zero
            this%AdvtLM(lm,nR)=zero
            this%AdvpLM(lm,nR)=zero
            if ( l_mag ) then
               this%VxBrLM(lm,nR)=zero
               this%VxBtLM(lm,nR)=zero
               this%VxBpLM(lm,nR)=zero
            end if
            if ( l_heat ) then
               this%VSrLM(lm,nR)=zero
               this%VStLM(lm,nR)=zero
               this%VSpLM(lm,nR)=zero
            end if
            if ( l_anel ) this%heatTermsLM(lm,nR)=zero
            if ( l_phase_field ) this%dphidtLM(lm,nR)=zero
            if ( l_chemical_conv ) then
               this%VXirLM(lm,nR)=zero
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
   subroutine get_td(this,lPressNext,dVSrLM,dVXirLM,dVxVhLM,dVxBhLM, &
              &      dwdt,dzdt,dpdt,dsdt,dxidt,dphidt,dbdt,djdt)
      !
      !  Purpose of this to calculate time derivatives
      !  ``dwdt``,``dzdt``,``dpdt``,``dsdt``,``dxidt``,``dbdt``,``djdt``
      !  and auxiliary arrays ``dVSrLM``, ``dVXirLM`` and ``dVxBhLM``, ``dVxVhLM``
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
      complex(cp), intent(out) :: dphidt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dbdt(lm_maxMag,nRstartMag:nRstopMag)
      complex(cp), intent(out) :: djdt(lm_maxMag,nRstartMag:nRstopMag)
      complex(cp), intent(out) :: dVxBhLM(lm_maxMag,nRstartMag:nRstopMag)
      complex(cp), intent(out) :: dVxVhLM(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dVSrLM(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dVXirLM(lm_max,nRstart:nRstop)

      !-- Local variables:
      integer :: l,m,lm,lmS,lmA,lmP,nR
      complex(cp) :: AdvPol_loc,CorPol_loc,AdvTor_loc,CorTor_loc,dsdt_loc

      if ( l_conv ) then

#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do collapse(2) &
         !$omp private(l,m,lmS,lmA,lmP) &
         !$omp private(AdvPol_loc,CorPol_loc,AdvTor_loc,CorTor_loc)
#else
         !$omp parallel do default(shared) private(lm,l,m,lmS,lmA,lmP) &
         !$omp private(AdvPol_loc,CorPol_loc,AdvTor_loc,CorTor_loc)
#endif
         do nR=nRstart,nRstop
            do lm=1,lm_max
               l   =lm2l(lm)
               m   =lm2m(lm)
               lmS =lm2lmS(lm)
               lmA =lm2lmA(lm)
               lmP =lm2lmP(lm)


               if ( l == 0 ) then ! This is l=0,m=0
                  lmA=lm2lmA(1)
                  lmP=1
                  !lmPA=lmP2lmPA(lmP)
                  if ( l_conv_nl ) then
                     AdvPol_loc=or2(nR)*this%AdvrLM(lmP,nR)
                     AdvTor_loc=zero!-dTheta1A(lm)*this%AdvpLM(lmPA,nR)
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
                           CorPol_loc=zero
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
                        AdvPol_loc    =dLh(lm)*or4(nR)*orho1(nR)*this%AdvrLM(lmP,nR)
                        dVxVhLM(lm,nR)=-orho1(nR)*r(nR)*r(nR)*dLh(lm)*this%AdvtLM(lmP,nR)
                     else
                        AdvPol_loc    =zero
                        dVxVhLM(lm,nR)=zero
                     endif

                  else ! We don't use the double curl

                     if ( l_corr ) then
                        if ( l < l_max .and. l > m ) then
                           CorPol_loc =two*CorFac*or1(nR) * (  &
                           &        dPhi(lm)*dw_Rloc(lm,nR) +  & ! phi-deriv of dw/dr
                           &    dTheta2A(lm)*z_Rloc(lmA,nR) -  & ! sin(theta) dtheta z
                           &    dTheta2S(lm)*z_Rloc(lmS,nR) )
                        else if ( l == l_max ) then
                           CorPol_loc=zero
                        else if ( l == m ) then
                           CorPol_loc = two*CorFac*or1(nR) * (  &
                           &         dPhi(lm)*dw_Rloc(lm,nR)  + &
                           &     dTheta2A(lm)*z_Rloc(lmA,nR) )
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

                  dwdt(lm,nR)=AdvPol_loc+CorPol_loc

                  if ( l_corr ) then
                     if ( l < l_max .and. l > m ) then
                        CorTor_loc=          two*CorFac*or2(nR) * (  &
                        &                 dPhi(lm)*z_Rloc(lm,nR)   + &
                        &            dTheta3A(lm)*dw_Rloc(lmA,nR)  + &
                        &    or1(nR)*dTheta4A(lm)* w_Rloc(lmA,nR)  + &
                        &            dTheta3S(lm)*dw_Rloc(lmS,nR)  - &
                        &    or1(nR)*dTheta4S(lm)* w_Rloc(lmS,nR)  )
                     else if ( l == l_max ) then
                        CorTor_loc=zero
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
                     AdvTor_loc=dLh(lm)*this%AdvpLM(lmP,nR)
                  else
                     AdvTor_loc=zero
                  end if
                  dzdt(lm,nR)=CorTor_loc+AdvTor_loc

                  ! In case double curl is calculated dpdt is useless
                  if ( (.not. l_double_curl) .or. lPressNext ) then
                     !------ Recycle CorPol and AdvPol:
                     if ( l_corr ) then
                        if ( l < l_max .and. l > m ) then
                           CorPol_loc=           two*CorFac*or2(nR) *  &
                           &           ( -dPhi(lm)  * ( dw_Rloc(lm,nR) &
                           &            +or1(nR)*dLh(lm)*w_Rloc(lm,nR) &
                           &                                         ) &
                           &              +dTheta3A(lm)*z_Rloc(lmA,nR) &
                           &              +dTheta3S(lm)*z_Rloc(lmS,nR) &
                           &           )

                        else if ( l == l_max ) then
                           CorPol_loc=zero

                        else if ( l == m ) then
                           CorPol_loc=                    two*CorFac*or2(nR) *  &
                           &                    ( -dPhi(lm)  * ( dw_Rloc(lm,nR) &
                           &                     +or1(nR)*dLh(lm)*w_Rloc(lm,nR) &
                           &                                                   )&
                           &                      +dTheta3A(lm)*z_Rloc(lmA,nR)  &
                           &                    )

                        end if
                     else
                        CorPol_loc=zero
                     end if
                     if ( l_conv_nl ) then
                        AdvPol_loc=-dLh(lm)*this%AdvtLM(lmP,nR)
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
         !$omp private(lmP,dsdt_loc,l)
#else
         !$omp parallel do default(shared) private(lm,lmP,dsdt_loc,l)
#endif
         do nR=nRstart,nRstop
            do lm=1,lm_max
               l   =lm2l(lm)
               lmP =lm2lmP(lm)
               if ( l == 0 ) then
                  dsdt_loc =epsc*epscProf(nR)!+opr/epsS*divKtemp0(nR)
                  dVSrLM(1,nR)=this%VSrLM(1,nR)
                  if ( l_anel ) then
                     if ( l_anelastic_liquid ) then
                        dsdt_loc=dsdt_loc+temp0(nR)*this%heatTermsLM(1,nR)
                     else
                        dsdt_loc=dsdt_loc+this%heatTermsLM(1,nR)
                     end if
                  end if
                  dsdt(1,nR)=dsdt_loc
               else
                  dVSrLM(lm,nR)=this%VSrLM(lmP,nR)
                  dsdt_loc     =dLh(lm)*this%VStLM(lmP,nR)
                  if ( l_anel ) then
                     if ( l_anelastic_liquid ) then
                        dsdt_loc=dsdt_loc+temp0(nR)*this%heatTermsLM(lmP,nR)
                     else
                        dsdt_loc=dsdt_loc+this%heatTermsLM(lmP,nR)
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
         !$omp private(lmP,l)
#else
         !$omp parallel do default(shared) private(lm,lmP,l)
#endif
         do nR=nRstart,nRstop
            do lm=1,lm_max
               l   =lm2l(lm)
               lmP =lm2lmP(lm)
               if ( l == 0 ) then
                  dVXirLM(1,nR)=this%VXirLM(1,nR)
                  dxidt(1,nR)  =epscXi
               else
                  dVXirLM(lm,nR)=this%VXirLM(lmP,nR)
                  dxidt(lm,nR)  =dLh(lm)*this%VXitLM(lmP,nR)
               end if
            end do
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#else
         !$omp end parallel do
#endif
      end if

      if ( l_phase_field ) then
#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do collapse(2) private(lmP)
#else
         !$omp parallel do default(shared) private(lm,lmP)
#endif
         do nR=nRstart,nRstop
            do lm=1,lm_max
               lmP=lm2lmP(lm)
               dphidt(lm,nR)=this%dphidtLM(lmP,nR)
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
         !$omp target teams distribute parallel do collapse(2) private(lmP)
#else
         !$omp parallel do default(shared) private(lm,lmP)
#endif
         do nR=nRstart,nRstop
            do lm=1,lm_max
               lmP =lm2lmP(lm)
               dbdt(lm,nR)   = dLh(lm)*this%VxBpLM(lmP,nR)
               dVxBhLM(lm,nR)=-dLh(lm)*this%VxBtLM(lmP,nR)*r(nR)*r(nR)
               djdt(lm,nR)   = dLh(lm)*or4(nR)*this%VxBrLM(lmP,nR)
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
                  !$omp target teams distribute parallel do private(lmP,l)
#else
                  !$omp parallel do default(shared) private(lm,lmP,l)
#endif
                  do lm=1,lm_max
                     l   =lm2l(lm)
                     lmP =lm2lmP(lm)
                     if ( l == 0 ) then
                        dVxBhLM(1,nR)=zero
                        dVSrLM(1,nR) =zero
                     else
                        dVxBhLM(lm,nR)=-dLh(lm)*this%VxBtLM(lmP,nR)*r(nR)*r(nR)
                        dVSrLM(lm,nR) =zero
                     end if
                  end do
#ifdef WITH_OMP_GPU
                  !$omp end target teams distribute parallel do
#else
                  !$omp end parallel do
#endif
               else
#ifdef WITH_OMP_GPU
                  !$omp target teams distribute parallel do
#else
                  !$omp parallel do
#endif
                  do lm=1,lm_max
                     if ( l_mag ) dVxBhLM(lm,nR)=zero
                     dVSrLM(lm,nR) =zero
                  end do
#ifdef WITH_OMP_GPU
                  !$omp end target teams distribute parallel do
#else
                  !$omp end parallel do
#endif
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
#endif
               end if
               if ( l_chemical_conv ) then
#ifdef WITH_OMP_GPU
                  !$omp target teams distribute parallel do
#else
                  !$omp parallel do
#endif
                  do lm=1,lm_max
                     dVXirLM(lm,nR)=zero
                  end do
#ifdef WITH_OMP_GPU
                  !$omp end target teams distribute parallel do
#else
                  !$omp end parallel do
#endif
               end if
               if ( l_heat ) then
#ifdef WITH_OMP_GPU
                  !$omp target teams distribute parallel do
#else
                  !$omp parallel do
#endif
                  do lm=1,lm_max
                     dVSrLM(lm,nR)=zero
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
!-----------------------------------------------------------------------------
end module nonlinear_lm_2d_mod
