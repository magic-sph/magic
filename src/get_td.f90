module nonlinear_lm_mod
   !
   ! This module is used to finish the assembly of the nonlinear terms in
   ! :math:`(\ell,m)` space. Derivatives along :math:`\theta` and :math:`\phi`
   ! are handled using recurrence relations.
   !

   use, intrinsic :: iso_c_binding
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: lm_max, l_max, lm_maxMag, lmP_max
   use logic, only : l_anel, l_conv_nl, l_corr, l_heat, l_anelastic_liquid, &
       &             l_mag_nl, l_mag_kin, l_mag_LF, l_conv, l_mag,          &
       &             l_chemical_conv, l_single_matrix, l_double_curl,       &
       &             l_adv_curl, l_phase_field
   use radial_functions, only: r, or2, or1, beta, epscProf, or4, temp0, orho1
   use physical_parameters, only: CorFac, epsc,  n_r_LCR, epscXi
   use blocking, only: lm2l, lm2m, lm2lmP, lmP2lmPS, lmP2lmPA, lm2lmA, &
       &               lm2lmS
   use horizontal_data, only: dLh, dPhi, dTheta2A, dTheta3A, dTheta4A, dTheta2S, &
       &                      dTheta3S, dTheta4S
   use constants, only: zero, two
   use fields, only: w_Rloc, dw_Rloc, ddw_Rloc, z_Rloc, dz_Rloc


   implicit none

   type :: nonlinear_lm_t
      !----- Nonlinear terms in lm-space:
      complex(cp), allocatable :: AdvrLM(:), AdvtLM(:), AdvpLM(:)
      complex(cp), allocatable :: VxBrLM(:), VxBtLM(:), VxBpLM(:)
      complex(cp), allocatable :: VSrLM(:),  VStLM(:),  VSpLM(:)
      complex(cp), allocatable :: VXirLM(:),  VXitLM(:),  VXipLM(:)
      complex(cp), allocatable :: heatTermsLM(:), dphidtLM(:)
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: set_zero
      procedure :: get_td
   end type nonlinear_lm_t

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
         allocate(this%VSrLM(lmP_max),this%VStLM(lmP_max),this%VSpLM(lmP_max))
         bytes_allocated = bytes_allocated + 3*lmP_max*SIZEOF_DEF_COMPLEX
      end if

      if ( l_chemical_conv ) then
         allocate(this%VXirLM(lmP_max),this%VXitLM(lmP_max),this%VXipLM(lmP_max))
         bytes_allocated = bytes_allocated + 3*lmP_max*SIZEOF_DEF_COMPLEX
      end if

      if ( l_phase_field ) then
         allocate(this%dphidtLM(lmP_max))
         bytes_allocated = bytes_allocated + lmP_max*SIZEOF_DEF_COMPLEX
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
      if ( l_chemical_conv ) deallocate( this%VXirLM, this%VXitLM, this%VXipLM )
      if ( l_heat ) deallocate( this%VSrLM, this%VStLM, this%VSpLM )
      if ( l_phase_field ) deallocate( this%dphidtLM )

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
         this%VxBrLM(lm)=zero
         this%VxBtLM(lm)=zero
         this%VxBpLM(lm)=zero
         if ( l_heat ) then
            this%VSrLM(lm)=zero
            this%VStLM(lm)=zero
            this%VSpLM(lm)=zero
         end if
         if ( l_anel ) this%heatTermsLM(lm)=zero
         if ( l_chemical_conv ) then
            this%VXirLM(lm)=zero
            this%VXitLM(lm)=zero
            this%VXipLM(lm)=zero
         end if
      end do
      !$omp end parallel do

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
      integer :: l,m,lm,lmS,lmA,lmP,lmPS,lmPA
      complex(cp) :: AdvPol_loc,CorPol_loc,AdvTor_loc,CorTor_loc,dsdt_loc
      !integer, parameter :: DOUBLE_COMPLEX_PER_CACHELINE=4

      !write(*,"(I3,A,4ES20.12)") nR,": get_td start: ",SUM(this%AdvrLM)

      !lm_chunksize=(((lm_max)/nThreads)/DOUBLE_COMPLEX_PER_CACHELINE) * &
      !             & DOUBLE_COMPLEX_PER_CACHELINE
      !lm_chunksize=4
      !write(*,"(A,I4)") "Using a chunksize of ",lm_chunksize

      if (nBc == 0 ) then

         if ( l_conv ) then  ! Convection

            lm =1   ! This is l=0,m=0
            lmA=lm2lmA(lm)
            lmP=1
            lmPA=lmP2lmPA(lmP)
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

            end do
            !$omp end parallel do

            ! In case double curl is calculated dpdt is useless
            if ( (.not. l_double_curl) .or. lPressNext ) then
            !if ( .true. ) then
               !$omp parallel do default(shared) private(lm,l,m,lmS,lmA,lmP) &
               !$omp private(lmPS,AdvPol_loc,CorPol_loc)
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
               !$omp end parallel do
            end if

         else
            do lm=2,lm_max
               dwdt(lm)=zero
               dzdt(lm)=zero
               dpdt(lm)=zero
            end do
         end if ! l_conv ?

      end if

      if ( nBc == 0 ) then

         if ( l_heat ) then
            dsdt_loc  =epsc*epscProf(nR)!+opr/epsS*divKtemp0(nR)
            dVSrLM(1)=this%VSrLM(1)
            if ( l_anel ) then
               if ( l_anelastic_liquid ) then
                  dsdt_loc=dsdt_loc+temp0(nR)*this%heatTermsLM(1)
               else
                  dsdt_loc=dsdt_loc+this%heatTermsLM(1)
               end if
            end if
            dsdt(1)=dsdt_loc

            !$omp parallel do default(shared) private(lm,lmP,dsdt_loc,l)
            do lm=2,lm_max
               l   =lm2l(lm)
               lmP =lm2lmP(lm)
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
            end do
            !$omp end parallel do
         end if

         if ( l_chemical_conv ) then
            dVXirLM(1)=this%VXirLM(1)
            dxidt(1)  =epscXi

            !$omp parallel do default(shared) private(lm,lmP)
            do lm=2,lm_max
               lmP=lm2lmP(lm)
               dVXirLM(lm)=this%VXirLM(lmP)
               dxidt(lm)  =dLh(lm)*this%VXitLM(lmP)
            end do
            !$omp end parallel do
         end if

         if ( l_phase_field ) then
            !$omp parallel do default(shared) private(lm,lmP)
            do lm=1,lm_max
               lmP=lm2lmP(lm)
               dphidt(lm)=this%dphidtLM(lmP)
            end do
            !$omp end parallel do
         end if

         if ( l_mag_nl .or. l_mag_kin  ) then
            !$omp parallel do default(shared) private(lm,lmP)
            do lm=1,lm_max
               lmP =lm2lmP(lm)
               dbdt(lm)   = dLh(lm)*this%VxBpLM(lmP)
               dVxBhLM(lm)=-dLh(lm)*this%VxBtLM(lmP)*r(nR)*r(nR)
               djdt(lm)   = dLh(lm)*or4(nR)*this%VxBrLM(lmP)
            end do
            !$omp end parallel do
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
         if ( l_mag_nl .or. l_mag_kin ) then
            dVxBhLM(1)=zero
            dVSrLM(1) =zero
            !$omp parallel do default(shared) private(lm,lmP)
            do lm=2,lm_max
               lmP =lm2lmP(lm)
               dVxBhLM(lm)=-dLh(lm)*this%VxBtLM(lmP)*r(nR)*r(nR)
               dVSrLM(lm) =zero
            end do
            !$omp end parallel do
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
