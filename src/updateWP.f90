#include "perflib_preproc.cpp"
module updateWP_mod
   !
   ! This module handles the time advance of the poloidal potential w and the pressure p.
   ! It contains the computation of the implicit terms and the linear solves.
   !

   use omp_lib
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, n_r_cmb, n_r_icb, get_openmp_blocks, &
       &                 n_mlo_loc, n_lo_loc, nRstart, nRstop, n_lm_loc
   use radial_functions, only: or1, or2, rho0, rgrav, visc, dLvisc, r, &
       &                       alpha0, temp0, beta, dbeta, ogrun,      &
       &                       rscheme_oc, ddLvisc, ddbeta, orho1
   use physical_parameters, only: kbotv, ktopv, ra, BuoFac, ChemFac,   &
       &                          ViscHeatFac, ThExpNb, ktopp
   use num_param, only: dct_counter, solve_counter
   use horizontal_data, only: hdif_V
   use logic, only: l_update_v, l_chemical_conv, l_RMS, l_double_curl, &
       &            l_fluxProfs, l_finite_diff, l_full_sphere, l_heat
   use RMS, only: DifPol2hInt, DifPolLMr
   use RMS_helpers, only:  hInt2Pol
   use radial_der, only: get_dddr, get_ddr, get_dr, get_dr_Rloc
   use integration, only: rInt_R
   use fields, only: work_LMdist
   use constants, only: zero, one, two, three, four, third, half
   use useful, only: abortRun
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray
   use dense_matrices
   use real_matrices
   use band_matrices
   use LMmapping, only: map_mlo

   implicit none

   private

   !-- Input of recycled work arrays:
   complex(cp), allocatable :: ddddw(:,:)
   complex(cp), allocatable :: dwold(:,:)
   real(cp), allocatable :: work(:)
   complex(cp), allocatable :: Dif(:),Pre(:),Buo(:)
   real(cp), allocatable :: rhs1(:,:,:)
   real(cp), allocatable :: rhs0(:,:,:)
   real(cp), allocatable :: wpMat_fac(:,:,:)
   class(type_realmat), allocatable :: wpMat(:), p0Mat, ellMat(:)
   logical, public, allocatable :: lWPmat(:)
   logical, allocatable :: l_ellMat(:)

   integer :: maxThreads, size_rhs1

   public :: initialize_updateWP, finalize_updateWP, updateWP, assemble_pol, &
   &         finish_exp_pol, get_pol_rhs_imp, finish_exp_pol_Rdist

contains

   subroutine initialize_updateWP(tscheme)
      !
      ! Purpose of this subroutine is to allocate the matrices needed
      ! to time advance the poloidal/pressure equations. Depending on the
      ! radial scheme, it can be either full or band matrices.
      !

      !-- Input variable
      class(type_tscheme), intent(in) :: tscheme ! time scheme

      !-- Local variables:
      integer :: ll, n_bands

      if ( l_finite_diff ) then
         allocate( type_bandmat :: wpMat(n_lo_loc) )

         if ( rscheme_oc%order <= 2 .and. rscheme_oc%order_boundary <= 2 ) then
            n_bands =rscheme_oc%order+3
         else
            n_bands = max(rscheme_oc%order+3,2*rscheme_oc%order_boundary+3)
         end if
         !print*, 'WP', n_bands
         do ll=1,n_lo_loc
            call wpMat(ll)%initialize(n_bands,n_r_max,l_pivot=.true.)
         end do
         allocate( wpMat_fac(n_r_max,2,n_lo_loc) )
         bytes_allocated=bytes_allocated+2*n_r_max*n_lo_loc*    &
         &               SIZEOF_DEF_REAL

         allocate( type_bandmat :: p0Mat )
         n_bands = rscheme_oc%order+1
         call p0Mat%initialize(n_bands,n_r_max,l_pivot=.true.)
      else
         allocate( type_densemat :: wpMat(n_lo_loc) )
         if ( l_double_curl ) then
            do ll=1,n_lo_loc
               call wpMat(ll)%initialize(n_r_max,n_r_max,l_pivot=.true.)
            end do
            allocate( wpMat_fac(n_r_max,2,n_lo_loc) )
            bytes_allocated=bytes_allocated+2*n_r_max*n_lo_loc*    &
            &               SIZEOF_DEF_REAL
         else
            do ll=1,n_lo_loc
               call wpMat(ll)%initialize(2*n_r_max,2*n_r_max,l_pivot=.true.)
            end do
            allocate( wpMat_fac(2*n_r_max,2,n_lo_loc) )
            bytes_allocated=bytes_allocated+4*n_r_max*n_lo_loc*    &
            &               SIZEOF_DEF_REAL
         end if

         allocate( type_densemat :: p0Mat )
         call p0Mat%initialize(n_r_max,n_r_max,l_pivot=.true.)
      end if

      allocate( lWPmat(n_lo_loc) )
      bytes_allocated=bytes_allocated+n_lo_loc*SIZEOF_LOGICAL

      if ( l_double_curl ) then
         allocate( ddddw(n_mlo_loc,n_r_max) )
         bytes_allocated = bytes_allocated+(n_mlo_loc)*n_r_max*SIZEOF_DEF_COMPLEX
         if ( l_RMS .or. l_FluxProfs ) then
            allocate( dwold(n_mlo_loc,n_r_max) )
            bytes_allocated = bytes_allocated+(n_mlo_loc)*n_r_max*SIZEOF_DEF_COMPLEX
            dwold(:,:)=zero
         end if
      end if

      allocate( work(n_r_max) )
      bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_REAL

      allocate( Dif(n_mlo_loc), Pre(n_mlo_loc), Buo(n_mlo_loc) )
      bytes_allocated = bytes_allocated+3*(n_mlo_loc)*SIZEOF_DEF_COMPLEX

      if ( l_double_curl ) then
         size_rhs1 = n_r_max
         allocate( rhs1(n_r_max,2*maxval(map_mlo%n_mi(:)),1) )
         bytes_allocated=bytes_allocated+n_r_max*maxThreads* &
         &               map_mlo%n_mi_max*SIZEOF_DEF_COMPLEX
      else
         size_rhs1 = 2*n_r_max
         allocate( rhs1(2*n_r_max,2*map_mlo%n_mi_max,1) )
         bytes_allocated=bytes_allocated+2*n_r_max*maxThreads* &
         &               map_mlo%n_mi_max*SIZEOF_DEF_COMPLEX
      end if

      if ( tscheme%l_assembly .and. l_double_curl ) then
         allocate( type_bandmat :: ellMat(n_lo_loc) )
         if ( rscheme_oc%order <= 2 .and. rscheme_oc%order_boundary <= 2 .and. &
         &    ktopv /=1 .and. kbotv /=1 ) then
            !n_bands =rscheme_oc%order+1 # should be that but yield matrix singularity?
            n_bands = max(rscheme_oc%order+1,2*rscheme_oc%order_boundary+1)
         else
            n_bands = max(rscheme_oc%order+1,2*rscheme_oc%order_boundary+1)
         end if
         do ll=1,n_lo_loc
            call ellMat(ll)%initialize(n_bands,n_r_max,l_pivot=.true.)
         end do
         allocate( l_ellMat(n_lo_loc) )
         l_ellMat(:) = .false.
         allocate( rhs0(n_r_max,2*map_mlo%n_mi_max,1) )
         rhs0(:,:,:)=zero
         bytes_allocated = bytes_allocated+n_lo_loc*SIZEOF_LOGICAL+&
         &                 n_r_max*maxThreads*2*map_mlo%n_mi_max*SIZEOF_DEF_REAL
      end if

   end subroutine initialize_updateWP
!-----------------------------------------------------------------------------
   subroutine finalize_updateWP(tscheme)
      !
      ! Deallocation of the matrices used to time-advance the poloidal/pressure
      ! equations.
      !

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme ! time scheme

      !-- Local variables:
      integer :: ll

      if ( tscheme%l_assembly .and. l_double_curl ) then
         do ll=1,n_lo_loc
            call ellMat(ll)%finalize()
         end do
         deallocate( l_ellMat, rhs0 )
      end if

      do ll=1,n_lo_loc
         call wpMat(ll)%finalize()
      end do
      call p0Mat%finalize()

      deallocate( wpMat_fac,lWPmat, rhs1, work )
      deallocate( Dif, Pre, Buo )
      if ( l_double_curl ) then
         deallocate( ddddw )
         if ( l_RMS .or. l_FluxProfs ) deallocate( dwold )
      end if

   end subroutine finalize_updateWP
!-----------------------------------------------------------------------------
   subroutine updateWP(s, xi, w, dw, ddw, dwdt, p, dp, dpdt, tscheme, &
              &        lRmsNext, lPressNext)
      !
      !  updates the poloidal velocity potential w, the pressure p, and
      !  their radial derivatives.
      !

      !-- Input/output of scalar fields:
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lRmsNext
      logical,             intent(in) :: lPressNext

      type(type_tarray), intent(inout) :: dpdt
      type(type_tarray), intent(inout) :: dwdt
      complex(cp),       intent(inout) :: s(n_mlo_loc,n_r_max)
      complex(cp),       intent(inout) :: xi(n_mlo_loc,n_r_max)
      complex(cp),       intent(inout) :: w(n_mlo_loc,n_r_max)
      complex(cp),       intent(inout) :: dw(n_mlo_loc,n_r_max)
      complex(cp),       intent(inout) :: ddw(n_mlo_loc,n_r_max)
      complex(cp),       intent(inout) :: p(n_mlo_loc,n_r_max)

      complex(cp),       intent(out) :: dp(n_mlo_loc,n_r_max)

      !-- Local variables:
      integer :: l, m               ! degree and order corresponding
      integer :: lj, mi, i, nRHS    ! l, m and ml counter
      integer :: nR                 ! counts radial grid points
      integer :: n_r_out            ! counts cheb modes
      real(cp) :: rhs(n_r_max)  ! real RHS for l=m=0


      if ( .not. l_update_v ) return
      
      !-- Now assemble the right hand side and store it in work_LMdist
      call tscheme%set_imex_rhs(work_LMdist, dwdt, 1, n_mlo_loc, n_r_max)
      if ( .not. l_double_curl ) then
         call tscheme%set_imex_rhs(ddw, dpdt, 1, n_mlo_loc, n_r_max)
      end if

      call solve_counter%start_count()
      ! Loop over local l
      do lj=1, n_lo_loc
         l = map_mlo%lj2l(lj)

         if ( .not. lWPmat(lj) ) then
            if ( l == 0 ) then
               call get_p0Mat(p0Mat)
            else
               if ( l_double_curl ) then
                  call get_wMat(tscheme,l,hdif_V(l),wpMat(lj),wpMat_fac(:,:,lj))
               else
                  call get_wpMat(tscheme,l,hdif_V(l),wpMat(lj),wpMat_fac(:,:,lj))
               end if
            end if         
            lWPmat(lj)=.true.
         end if

         ! Build RHS
         ! Loop over local m corresponding to current l
         nRHS = map_mlo%n_mi(lj)
         do mi=1,nRHS
            m = map_mlo%milj2m(mi,lj)
            i = map_mlo%milj2i(mi,lj)
            
            if ( l==0 ) then
            
               !-- The integral of rho' r^2 dr vanishes
               if ( ThExpNb*ViscHeatFac /= 0 .and. ktopp==1 ) then
                  if ( rscheme_oc%version == 'cheb' ) then
                     do nR=1,n_r_max
                        work(nR)=ThExpNb*alpha0(nR)*temp0(nR)*rho0(nR)*r(nR)*&
                        &        r(nR)*real(s(map_mlo%m0l0,nR))
                     end do
                     rhs(1)=rInt_R(work,r,rscheme_oc)
                  else
                     rhs(1)=0.0_cp
                  end if
               else
                  rhs(1)=0.0_cp
               end if

               if ( l_chemical_conv ) then
                  do nR=2,n_r_max
                     rhs(nR)=rho0(nR)*BuoFac*rgrav(nR)*real(s(i,nR))+   &
                     &       rho0(nR)*ChemFac*rgrav(nR)*real(xi(i,nR))+ &
                     &       real(dwdt%expl(i,nR,tscheme%istage))
                  end do
               else
                  do nR=2,n_r_max
                     rhs(nR)=rho0(nR)*BuoFac*rgrav(nR)*real(s(i,nR))+  &
                     &       real(dwdt%expl(i,nR,tscheme%istage))
                  end do
               end if

               call p0Mat%solve(rhs)
               
            else ! l /= 0
               rhs1(1,2*mi-1,1)      =0.0_cp
               rhs1(1,2*mi,1)        =0.0_cp
               rhs1(n_r_max,2*mi-1,1)=0.0_cp
               rhs1(n_r_max,2*mi,1)  =0.0_cp
               if ( l_double_curl ) then
                  rhs1(2,2*mi-1,1)        =0.0_cp
                  rhs1(2,2*mi,1)          =0.0_cp
                  rhs1(n_r_max-1,2*mi-1,1)=0.0_cp
                  rhs1(n_r_max-1,2*mi,1)  =0.0_cp
                  do nR=3,n_r_max-2
                     rhs1(nR,2*mi-1,1)= real(work_LMdist(i,nR))
                     rhs1(nR,2*mi,1)  =aimag(work_LMdist(i,nR))
                  end do

                  if ( l_heat ) then
                     do nR=3,n_r_max-2
                        rhs1(nR,2*mi-1,1)=rhs1(nR,2*mi-1,1)+             &
                        &      tscheme%wimp_lin(1)*real(l*(l+1),cp) *    &
                        &      or2(nR)*BuoFac*rgrav(nR)*real(s(i,nR))
                        rhs1(nR,2*mi,1)  =rhs1(nR,2*mi,1)+               &
                        &      tscheme%wimp_lin(1)*real(l*(l+1),cp) *    &
                        &      or2(nR)*BuoFac*rgrav(nR)*aimag(s(i,nR))
                     end do
                  end if

                  if ( l_chemical_conv ) then
                     do nR=3,n_r_max-2
                        rhs1(nR,2*mi-1,1)=rhs1(nR,2*mi-1,1)+          &
                        &      tscheme%wimp_lin(1)*real(l*(l+1),cp) * &
                        &      or2(nR)*ChemFac*rgrav(nR)*real(xi(i,nR))
                        rhs1(nR,2*mi,1)  =rhs1(nR,2*mi,1)+            &
                        &      tscheme%wimp_lin(1)*real(l*(l+1),cp) * &
                        &      or2(nR)*ChemFac*rgrav(nR)*aimag(xi(i,nR))
                     end do
                  end if
               else
                  rhs1(n_r_max+1,2*mi-1,1)=0.0_cp
                  rhs1(n_r_max+1,2*mi,1)  =0.0_cp
                  rhs1(2*n_r_max,2*mi-1,1)=0.0_cp
                  rhs1(2*n_r_max,2*mi,1)  =0.0_cp
                  do nR=2,n_r_max-1
                     rhs1(nR,2*mi-1,1)        = real(work_LMdist(i,nR))
                     rhs1(nR,2*mi,1)          =aimag(work_LMdist(i,nR))
                     rhs1(nR+n_r_max,2*mi-1,1)= real(ddw(i,nR)) ! ddw is a work array
                     rhs1(nR+n_r_max,2*mi,1)  =aimag(ddw(i,nR))
                  end do

                  if ( l_heat ) then
                     do nR=2,n_r_max-1
                        rhs1(nR,2*mi-1,1)=rhs1(nR,2*mi-1,1)+              &
                        &         tscheme%wimp_lin(1)*rho0(nR)*BuoFac*    &
                        &                      rgrav(nR)*real(s(i,nR))
                        rhs1(nR,2*mi,1)  =rhs1(nR,2*mi,1)+                &
                        &         tscheme%wimp_lin(1)*rho0(nR)*BuoFac*    &
                        &                      rgrav(nR)*aimag(s(i,nR))
                     end do
                  end if

                  if ( l_chemical_conv ) then
                     do nR=2,n_r_max-1
                        rhs1(nR,2*mi-1,1)=rhs1(nR,2*mi-1,1)+              &
                        &         tscheme%wimp_lin(1)*rho0(nR)*ChemFac*   &
                        &                      rgrav(nR)*real(xi(i,nR))
                        rhs1(nR,2*mi,1)  =rhs1(nR,2*mi,1)+                &
                        &         tscheme%wimp_lin(1)*rho0(nR)*ChemFac*   &
                        &                      rgrav(nR)*aimag(xi(i,nR))
                     end do
                  end if
                  
               end if
            end if

         end do ! end of loop over local m corresponding to current l
            
         ! Rescale and Solve for all RHS at once
         if ( l>0 ) then

            ! use the mat_fac(:,1) to scale the rhs
            do mi=1,nRHS
               do nR=1,size_rhs1
                  rhs1(nR,2*mi-1,1)=rhs1(nR,2*mi-1,1)*wpMat_fac(nR,1,lj)
                  rhs1(nR,2*mi,1)  =rhs1(nR,2*mi,1)*wpMat_fac(nR,1,lj)
               end do
            end do
            
            call wpMat(lj)%solve(rhs1(:,1:2*nRHS,1),2*nRHS)
            
            ! rescale the solution with mat_fac(:,2)
            do mi=1,nRHS
               do nR=1,size_rhs1
                  rhs1(nR,2*mi-1,1)=rhs1(nR,2*mi-1,1)*wpMat_fac(nR,2,lj)
                  rhs1(nR,2*mi,1)  =rhs1(nR,2*mi,1)*wpMat_fac(nR,2,lj)
               end do
            end do
         end if
         
         
         do mi=1,nRHS
            i = map_mlo%milj2i(mi,lj)
            m = map_mlo%milj2m(mi,lj)
            
            ! Store old dw
            if ( l_double_curl .and. lPressNext .and. tscheme%istage == 1) then
               do nR=1,n_r_max
                  dwold(i,nR)=dw(i,nR)
               end do
            end if
            
            if ( l==0 ) then
               do n_r_out=1,rscheme_oc%n_max
                  p(i,n_r_out)=rhs(n_r_out)
               end do

            else
               if ( l_double_curl ) then
                  if ( m > 0 ) then
                     do n_r_out=1,rscheme_oc%n_max
                        w(i,n_r_out)  =cmplx(rhs1(n_r_out,2*mi-1,1), &
                        &                    rhs1(n_r_out,2*mi,1),cp)
                     end do
                  else
                     do n_r_out=1,rscheme_oc%n_max
                        w(i,n_r_out)  = cmplx(rhs1(n_r_out,2*mi-1,1),0.0_cp,cp)
                     end do
                  end if
               else
                  if ( m > 0 ) then
                     do n_r_out=1,rscheme_oc%n_max
                        w(i,n_r_out)=cmplx(rhs1(n_r_out,2*mi-1,1), &
                        &                  rhs1(n_r_out,2*mi,1),cp)
                        p(i,n_r_out)=cmplx(rhs1(n_r_max+n_r_out,2*mi-1,1), &
                        &                  rhs1(n_r_max+n_r_out,2*mi,1),cp)
                     end do
                  else
                     do n_r_out=1,rscheme_oc%n_max
                        w(i,n_r_out)= cmplx(rhs1(n_r_out,2*mi-1,1),0.0_cp,cp)
                        p(i,n_r_out)= cmplx(rhs1(n_r_max+n_r_out,2*mi-1,1), &
                        &                   0.0_cp,cp)
                     end do
                  end if
               end if
               
            end if
            
         end do

      end do   ! end of loop over local l

      call solve_counter%stop_count()


      !-- set cheb modes > rscheme_oc%n_max to zero (dealiazing)
      do i=1,n_mlo_loc
         do n_r_out=rscheme_oc%n_max+1,n_r_max
            w(i,n_r_out)=zero
            p(i,n_r_out)=zero
         end do
      end do

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dwdt, 1, n_mlo_loc, n_r_max)
      if ( .not. l_double_curl ) call tscheme%rotate_imex(dpdt, 1, n_mlo_loc, n_r_max)

      if ( tscheme%istage == tscheme%nstages ) then
         call get_pol_rhs_imp(s, xi, w, dw, ddw, p, dp, dwdt, dpdt,       &
              &               tscheme, 1, tscheme%l_imp_calc_rhs(1),      &
              &               lPressNext, lRmsNext,                       &
              &               dpdt%expl(:,:,1), l_in_cheb_space=.true.)
         ! dpdt%expl(:,:,1) needed for RMS calc: first stage explicit term
      else
         call get_pol_rhs_imp(s, xi, w, dw, ddw, p, dp, dwdt, dpdt,       &
              &               tscheme, tscheme%istage+1,                  &
              &               tscheme%l_imp_calc_rhs(tscheme%istage+1),   &
              &               lPressNext, lRmsNext,                       &
              &               dpdt%expl(:,:,1), l_in_cheb_space=.true.)
      end if


   end subroutine updateWP
!------------------------------------------------------------------------------
   subroutine get_pol(w, work)
      !
      !  Get the poloidal potential from the solve of an elliptic equation.
      !  Careful: the output is in Chebyshev space!
      !

      !-- Input field
      complex(cp), intent(in) :: work(n_mlo_loc,n_r_max)

      !-- Output field
      complex(cp), intent(out) :: w(n_mlo_loc,n_r_max)

      !-- Local variables:
      integer :: l, m          ! degree and order
      integer :: lj, mi, i     ! l, m and ml counter
      integer :: nR            ! counts radial grid points
      integer :: n_r_out       ! counts cheb modes

      if ( .not. l_update_v ) return

      call solve_counter%start_count()

      ! Loop over local l
      do lj=1, n_lo_loc
         l = map_mlo%lj2l(lj)
         if ( l == 0 ) cycle

         if ( .not. l_ellMat(lj) ) then
            call get_elliptic_mat(l, ellMat(lj))
            l_ellMat(lj) = .true.
         end if

         ! Build RHS
         ! Loop over local m corresponding to current l
         do mi=1,map_mlo%n_mi(lj)
            m = map_mlo%milj2m(mi,lj)
            i = map_mlo%milj2i(mi,lj)

            rhs0(1,2*mi-1,1)        =0.0_cp
            rhs0(1,2*mi,1)          =0.0_cp
            rhs0(2,2*mi-1,1)        =0.0_cp
            rhs0(2,2*mi,1)          =0.0_cp
            rhs0(n_r_max-1,2*mi-1,1)=0.0_cp
            rhs0(n_r_max-1,2*mi,1)  =0.0_cp
            rhs0(n_r_max,2*mi-1,1)  =0.0_cp
            rhs0(n_r_max,2*mi,1)    =0.0_cp
            do nR=3,n_r_max-2
               rhs0(nR,2*mi-1,1)= real(work(i,nR))
               rhs0(nR,2*mi,1)  =aimag(work(i,nR))
            end do
         end do 

         call ellMat(lj)%solve(rhs0(:,1:2*map_mlo%n_mi(lj),1),2*map_mlo%n_mi(lj))

         ! Loop over m corresponding to current l (again)
         do mi=1,map_mlo%n_mi(lj)
            m = map_mlo%milj2m(mi,lj)
            i = map_mlo%milj2i(mi,lj)

            if ( m > 0 ) then
               do n_r_out=1,rscheme_oc%n_max
                  w(i,n_r_out)= cmplx(rhs0(n_r_out,2*mi-1,1), &
                  &                   rhs0(n_r_out,2*mi,1),kind=cp)
               end do
            else
               do n_r_out=1,rscheme_oc%n_max
                  w(i,n_r_out)= cmplx(rhs0(n_r_out,2*mi-1,1),0.0_cp,kind=cp)
               end do
            end if
         end do  ! loop over m (again)
      end do   ! loop over l

      call solve_counter%stop_count()

      !-- set cheb modes > n_cheb_max to zero (dealiazing)
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         do i=1,n_mlo_loc
            w(i,n_r_out)=zero
         end do
      end do

   end subroutine get_pol
!------------------------------------------------------------------------------
   subroutine finish_exp_pol(dVxVhLM, dw_exp_last)

      !-- Input variables
      complex(cp), intent(inout) :: dVxVhLM(n_mlo_loc,n_r_max)

      !-- Output variables
      complex(cp), intent(inout) :: dw_exp_last(n_mlo_loc,n_r_max)

      !-- Local variables
      integer :: n_r, start_lm, stop_lm

      !$omp parallel default(shared) private(start_lm,stop_lm)
      start_lm=1; stop_lm=n_mlo_loc
      call get_openmp_blocks(start_lm,stop_lm)
      call get_dr( dVxVhLM, work_LMdist, n_mlo_loc, start_lm,    &
           &       stop_lm, n_r_max, rscheme_oc, nocopy=.true. )
      !$omp barrier

      !$omp do
      do n_r=1,n_r_max
         dw_exp_last(:,n_r)= dw_exp_last(:,n_r)+or2(n_r)*work_LMdist(:,n_r)
      end do
      !$omp end do
      !$omp end parallel

   end subroutine finish_exp_pol
!------------------------------------------------------------------------------
   subroutine finish_exp_pol_Rdist(dVxVhLM, dw_exp_last)

      !-- Input variables
      complex(cp), intent(inout) :: dVxVhLM(n_lm_loc,nRstart:nRstop)

      !-- Output variables
      complex(cp), intent(inout) :: dw_exp_last(n_lm_loc,nRstart:nRstop)

      !-- Local variables
      complex(cp) :: work_Rloc(n_lm_loc,nRstart:nRstop)
      integer :: n_r

      call get_dr_Rloc(dVxVhLM, work_Rloc, n_lm_loc, nRstart, nRstop, n_r_max, &
           &           rscheme_oc)

      !$omp parallel default(shared)
      !$omp do
      do n_r=nRstart,nRstop
         dw_exp_last(:,n_r)= dw_exp_last(:,n_r)+or2(n_r)*work_Rloc(:,n_r)
      end do
      !$omp end do
      !$omp end parallel

   end subroutine finish_exp_pol_Rdist
!------------------------------------------------------------------------------
   subroutine get_pol_rhs_imp(s, xi, w, dw, ddw, p, dp, dwdt, dpdt, tscheme,     &
              &               istage, l_calc_lin, lPressNext, lRmsNext, dp_expl, &
              &               l_in_cheb_space)
      !
      ! This subroutine computes the derivatives of w and p and assemble the
      ! implicit stage if needed.
      !

      !-- Input variables
      integer,             intent(in) :: istage
      class(type_tscheme), intent(in) :: tscheme
      complex(cp),         intent(in) :: s(n_mlo_loc,n_r_max)
      complex(cp),         intent(in) :: xi(n_mlo_loc,n_r_max)
      logical,             intent(in) :: l_calc_lin
      logical,             intent(in) :: lPressNext
      logical,             intent(in) :: lRmsNext
      logical, optional,   intent(in) :: l_in_cheb_space
      complex(cp),         intent(in) :: dp_expl(n_mlo_loc,n_r_max)

      !-- Output variables
      type(type_tarray), intent(inout) :: dwdt
      type(type_tarray), intent(inout) :: dpdt
      complex(cp),       intent(inout) :: w(n_mlo_loc,n_r_max)
      complex(cp),       intent(inout) :: p(n_mlo_loc,n_r_max)
      complex(cp),       intent(out) :: dp(n_mlo_loc,n_r_max)
      complex(cp),       intent(out) :: dw(n_mlo_loc,n_r_max)
      complex(cp),       intent(out) :: ddw(n_mlo_loc,n_r_max)

      !-- Local variables
      logical :: l_in_cheb
      integer :: n_r_top, n_r_bot, l
      integer :: n_r, lm, start_lm, stop_lm
      real(cp) :: dL

      if ( present(l_in_cheb_space) ) then
         l_in_cheb = l_in_cheb_space
      else
         l_in_cheb = .false.
      end if

      !$omp parallel default(shared)  private(start_lm, stop_lm)
      start_lm=1; stop_lm=n_mlo_loc
      call get_openmp_blocks(start_lm,stop_lm)

      !$omp single
      call dct_counter%start_count()
      !$omp end single
      if ( l_double_curl ) then
         call get_ddr( w, dw, ddw, n_mlo_loc, start_lm, stop_lm, n_r_max, &
              &        rscheme_oc, l_dct_in=.not. l_in_cheb )
         call get_ddr( ddw, work_LMdist, ddddw, n_mlo_loc, start_lm,  &
              &       stop_lm, n_r_max, rscheme_oc )
      else
         call get_dddr( w, dw, ddw, work_LMdist, n_mlo_loc, start_lm, &
              &         stop_lm, n_r_max, rscheme_oc, l_dct_in=.not. l_in_cheb)
         call get_dr( p, dp, n_mlo_loc, start_lm, stop_lm, &
              &       n_r_max, rscheme_oc, l_dct_in=.not. l_in_cheb)
         if ( l_in_cheb ) call rscheme_oc%costf1(p,n_mlo_loc,start_lm,stop_lm)
      end if
      if ( l_in_cheb ) call rscheme_oc%costf1(w,n_mlo_loc,start_lm,stop_lm)
      !$omp barrier
      !$omp single
      call dct_counter%stop_count()
      !$omp end single

      if ( istage == 1 ) then
         if ( l_double_curl ) then
            !$omp do private(n_r,lm,l,dL)
            do n_r=2,n_r_max-1
               do lm=1,n_mlo_loc
                  l = map_mlo%i2l(lm)
                  if ( l == 0 ) cycle
                  dL = real(l*(l+1),cp)
                  dwdt%old(lm,n_r,istage)=dL*or2(n_r)* ( -orho1(n_r)*(  &
                  &                   ddw(lm,n_r)-beta(n_r)*dw(lm,n_r)- &
                  &                            dL*or2(n_r)* w(lm,n_r) ) )
               end do
            end do
            !$omp end do
         else
            !$omp do private(n_r,lm,l,dL)
            do n_r=2,n_r_max-1
               do lm=1,n_mlo_loc
                  l = map_mlo%i2l(lm)
                  if ( l==0 ) cycle ! skips mode (0,0) if it is local
                  dL = real(l*(l+1),cp)
                  dwdt%old(lm,n_r,istage)= dL*or2(n_r)*w(lm,n_r)
                  dpdt%old(lm,n_r,istage)=-dL*or2(n_r)*dw(lm,n_r)
               end do
            end do
            !$omp end do
         end if
      end if

      if ( l_calc_lin .or. (tscheme%istage==tscheme%nstages .and. lRmsNext)) then

         if ( lRmsNext .and. tscheme%istage == tscheme%nstages ) then
            n_r_top=n_r_cmb
            n_r_bot=n_r_icb
         else
            n_r_top=n_r_cmb+1
            n_r_bot=n_r_icb-1
         end if

         !
         !-- Calculate explicit time step part:
         if ( l_double_curl ) then

            if ( lPressNext ) then
               n_r_top=n_r_cmb
               n_r_bot=n_r_icb
            end if

            !$omp do private(n_r,lm,l,Dif,Buo,dL)
            do n_r=n_r_top,n_r_bot
               do lm=1,n_mlo_loc
                  l=map_mlo%i2l(lm)
                  if ( l==0 ) cycle ! skips mode (0,0) if it is local
                  dL=real(l*(l+1),cp)

                  Dif(lm)=     -hdif_V(l)*dL*or2(n_r)*visc(n_r)*orho1(n_r)*      (  &
                  &                                                  ddddw(lm,n_r)  &
                  &           +two*( dLvisc(n_r)-beta(n_r) ) * work_LMdist(lm,n_r)  &
                  &        +( ddLvisc(n_r)-two*dbeta(n_r)+dLvisc(n_r)*dLvisc(n_r)+  &
                  &           beta(n_r)*beta(n_r)-three*dLvisc(n_r)*beta(n_r)-two*  &
                  &           or1(n_r)*(dLvisc(n_r)+beta(n_r))-two*or2(n_r)*dL ) *  &
                  &                                                    ddw(lm,n_r)  &
                  &        +( -ddbeta(n_r)-dbeta(n_r)*(two*dLvisc(n_r)-beta(n_r)+   &
                  &           two*or1(n_r))-ddLvisc(n_r)*(beta(n_r)+two*or1(n_r))+  &
                  &           beta(n_r)*beta(n_r)*(dLvisc(n_r)+two*or1(n_r))-       &
                  &           beta(n_r)*(dLvisc(n_r)*dLvisc(n_r)-two*or2(n_r))-     &
                  &           two*dLvisc(n_r)*or1(n_r)*(dLvisc(n_r)-or1(n_r))+      &
                  &           two*(two*or1(n_r)+beta(n_r)-dLvisc(n_r))*or2(n_r)*dL) &
                  &                                    *                dw(lm,n_r)  &
                  &        + dL*or2(n_r)* ( two*dbeta(n_r)+ddLvisc(n_r)+            &
                  &          dLvisc(n_r)*dLvisc(n_r)-two*third*beta(n_r)*beta(n_r)+ &
                  &          dLvisc(n_r)*beta(n_r)+two*or1(n_r)*(two*dLvisc(n_r)-   &
                  &          beta(n_r)-three*or1(n_r))+dL*or2(n_r) ) *   w(lm,n_r) )

                  Buo(lm) = zero
                  if ( l_heat ) Buo(lm) = BuoFac*dL*or2(n_r)*rgrav(n_r)*s(lm,n_r)
                  if ( l_chemical_conv ) Buo(lm) = Buo(lm)+ChemFac*dL* &
                  &                                     or2(n_r)*rgrav(n_r)*xi(lm,n_r)

                  dwdt%impl(lm,n_r,istage)=Dif(lm)+Buo(lm)

                  if ( lPressNext .and. tscheme%istage==tscheme%nstages) then
                     ! In the double curl formulation, we can estimate the pressure
                     ! if required.
                     p(lm,n_r)=-r(n_r)*r(n_r)/dL*                 dp_expl(lm,n_r)  &
                     &            -one/tscheme%dt(1)*(dw(lm,n_r)-dwold(lm,n_r))+   &
                     &              hdif_V(l)*visc(n_r)* ( work_LMdist(lm,n_r)     &
                     &                       - (beta(n_r)-dLvisc(n_r))*ddw(lm,n_r) &
                     &            - ( dL*or2(n_r)+dLvisc(n_r)*beta(n_r)+dbeta(n_r) &
                     &                  + two*(dLvisc(n_r)+beta(n_r))*or1(n_r)     &
                     &                                              ) * dw(lm,n_r) &
                     &             + dL*or2(n_r)*(two*or1(n_r)+two*third*beta(n_r) &
                     &                     +dLvisc(n_r) )   *            w(lm,n_r) )
                  end if

                  if ( lRmsNext .and. tscheme%istage==tscheme%nstages ) then
                     !-- In case RMS force balance is required, one needs to also
                     !-- compute the classical diffusivity that is used in the non
                     !-- double-curl version
                     Dif(lm) =  hdif_V(l)*dL*or2(n_r)*visc(n_r) *  ( ddw(lm,n_r)   &
                     &        +(two*dLvisc(n_r)-third*beta(n_r))*     dw(lm,n_r)   &
                     &        -( dL*or2(n_r)+four*third*( dbeta(n_r)+dLvisc(n_r)*  &
                     &           beta(n_r)+(three*dLvisc(n_r)+beta(n_r))*or1(n_r)))&
                     &                                         *       w(lm,n_r) )
                  end if
               end do
               if ( lRmsNext .and. tscheme%istage==tscheme%nstages ) then
                  call hInt2Pol(Dif,1,n_mlo_loc,n_r,DifPolLMr(:,n_r), &
                       &        DifPol2hInt(:,n_r))
               end if
            end do
            !$omp end do

         else

            !$omp do private(n_r,lm,l,Dif,Buo,Pre,dL)
            do n_r=n_r_top,n_r_bot
               do lm=1,n_mlo_loc
                  l=map_mlo%i2l(lm)
                  if ( l==0 ) cycle ! skips mode (0,0) if it is local
                  dL=real(l*(l+1),cp)

                  Dif(lm) = hdif_V(l)*dL*or2(n_r)*visc(n_r)*(       ddw(lm,n_r)  & 
                  &       +(two*dLvisc(n_r)-third*beta(n_r))*        dw(lm,n_r)  &
                  &       -( dL*or2(n_r)+four*third*( dbeta(n_r)+dLvisc(n_r)*    &
                  &         beta(n_r)+(three*dLvisc(n_r)+beta(n_r))*or1(n_r)) )* &
                  &                                                   w(lm,n_r)  )
                  Pre(lm) = -dp(lm,n_r)+beta(n_r)*p(lm,n_r)
                  Buo(lm) = zero
                  if ( l_heat )  Buo(lm) = BuoFac*rho0(n_r)*rgrav(n_r)*s(lm,n_r)
                  if ( l_chemical_conv ) Buo(lm) = Buo(lm)+ChemFac*rho0(n_r)* &
                  &                                rgrav(n_r)*xi(lm,n_r)
                  dwdt%impl(lm,n_r,istage)=Pre(lm)+Dif(lm)+Buo(lm)
                  dpdt%impl(lm,n_r,istage)=               dL*or2(n_r)*p(lm,n_r) &
                  &             + hdif_V(l)*visc(n_r)*dL*or2(n_r)               &
                  &                                    * ( -work_LMdist(lm,n_r) &
                  &                       + (beta(n_r)-dLvisc(n_r))*ddw(lm,n_r) &
                  &            + ( dL*or2(n_r)+dLvisc(n_r)*beta(n_r)+dbeta(n_r) &
                  &                  + two*(dLvisc(n_r)+beta(n_r))*or1(n_r)     &
                  &                                           ) *    dw(lm,n_r) &
                  &            - dL*or2(n_r)* ( two*or1(n_r)+two*third*beta(n_r)&
                  &                     +dLvisc(n_r) )   *           w(lm,n_r)  )
               end do
               if ( lRmsNext .and. tscheme%istage==tscheme%nstages ) then
                  call hInt2Pol(Dif,1,n_mlo_loc,n_r,DifPolLMr(:,n_r), &
                       &        DifPol2hInt(:,n_r))
               end if
            end do
            !$omp end do

         end if

      end if

      ! In case pressure is needed in the double curl formulation
      ! we also have to compute the radial derivative of p
      if ( lPressNext .and. l_double_curl ) then
         call get_dr( p, dp,n_mlo_loc, start_lm, stop_lm, n_r_max, rscheme_oc)
         !$omp barrier
      end if

      !$omp end parallel

   end subroutine get_pol_rhs_imp
!------------------------------------------------------------------------------
   subroutine assemble_pol(s, xi, w, dw, ddw, p, dp, dwdt, dpdt, dp_expl, &
              &            tscheme, lPressNext, lRmsNext)
      !
      ! This subroutine is used to assemble w and dw/dr when IMEX RK time schemes
      ! which necessitate an assembly stage are employed. Robin-type boundary
      ! conditions are enforced using Canuto (1986) approach.
      !

      !-- Input variables
      complex(cp),         intent(in) :: s(n_mlo_loc,n_r_max)
      complex(cp),         intent(in) :: xi(n_mlo_loc,n_r_max)
      complex(cp),         intent(in) :: dp_expl(n_mlo_loc,n_r_max)
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lPressNext
      logical,             intent(in) :: lRmsNext

      !-- Output variable
      type(type_tarray), intent(inout) :: dwdt
      type(type_tarray), intent(inout) :: dpdt
      complex(cp),       intent(inout) :: w(n_mlo_loc,n_r_max)
      complex(cp),       intent(out) :: dw(n_mlo_loc,n_r_max)
      complex(cp),       intent(out) :: ddw(n_mlo_loc,n_r_max)
      complex(cp),       intent(inout) :: p(n_mlo_loc,n_r_max)
      complex(cp),       intent(inout) :: dp(n_mlo_loc,n_r_max)

      !-- Local variables
      real(cp) :: fac_top, fac_bot
      integer :: n_r_top, n_r_bot, l, m
      integer :: n_r, lm, start_lm, stop_lm
      real(cp) :: dL

      call tscheme%assemble_imex(work_LMdist, dwdt, 1, n_mlo_loc, n_r_max)
      if ( l_double_curl) then
         call get_pol(w, work_LMdist)
      else
         call tscheme%assemble_imex(ddw, dpdt, 1, n_mlo_loc, n_r_max) ! Use ddw as a work array
      end if

      !$omp parallel default(shared)  private(start_lm, stop_lm)
      start_lm=1; stop_lm=n_mlo_loc
      call get_openmp_blocks(start_lm,stop_lm)

      if ( l_double_curl) then

         !$omp single
            call dct_counter%start_count()
         !$omp end single
         call get_ddr( w, dw, ddw, n_mlo_loc, start_lm, stop_lm, n_r_max, &
              &        rscheme_oc, l_dct_in=.false. )
         call get_ddr( ddw, work_LMdist, ddddw, n_mlo_loc, start_lm, stop_lm, &
              &        n_r_max, rscheme_oc )
         call rscheme_oc%costf1(w, n_mlo_loc, start_lm, stop_lm)
         !$omp barrier
         !$omp single
         call dct_counter%stop_count()
         !$omp end single

         !$omp do private(n_r,lm,l,dL)
         do n_r=2,n_r_max-1
            do lm=1,n_mlo_loc
               l = map_mlo%i2l(lm)
               if ( l == 0 ) cycle
               dL = real(l*(l+1),cp)
               dwdt%old(lm,n_r,1)=-dL*or2(n_r)*orho1(n_r)*(         &
               &                  ddw(lm,n_r)-beta(n_r)*dw(lm,n_r)- &
               &                            dL*or2(n_r)* w(lm,n_r) )
            end do
         end do
         !$omp end do

         if ( tscheme%l_imp_calc_rhs(1) .or. lRmsNext ) then
            if ( lRmsNext ) then
               n_r_top=n_r_cmb
               n_r_bot=n_r_icb
            else
               n_r_top=n_r_cmb+1
               n_r_bot=n_r_icb-1
            end if

            !$omp do private(n_r,lm,l,Dif,Buo,dL)
            do n_r=n_r_top,n_r_bot
               do lm=1,n_mlo_loc
                  l=map_mlo%i2l(lm)
                  if ( l == 0 ) cycle
                  dL=real(l*(l+1),cp)

                  Dif(lm)=-hdif_V(l)*dL*or2(n_r)*visc(n_r)*orho1(n_r)*      (      &
                  &                                                 ddddw(lm,n_r)  &
                  &          +two*( dLvisc(n_r)-beta(n_r) ) * work_LMdist(lm,n_r)  &
                  &       +( ddLvisc(n_r)-two*dbeta(n_r)+dLvisc(n_r)*dLvisc(n_r)+  &
                  &          beta(n_r)*beta(n_r)-three*dLvisc(n_r)*beta(n_r)-two*  &
                  &          or1(n_r)*(dLvisc(n_r)+beta(n_r))-two*or2(n_r)*dL ) *  &
                  &                                                   ddw(lm,n_r)  &
                  &       +( -ddbeta(n_r)-dbeta(n_r)*(two*dLvisc(n_r)-beta(n_r)+   &
                  &          two*or1(n_r))-ddLvisc(n_r)*(beta(n_r)+two*or1(n_r))+  &
                  &          beta(n_r)*beta(n_r)*(dLvisc(n_r)+two*or1(n_r))-       &
                  &          beta(n_r)*(dLvisc(n_r)*dLvisc(n_r)-two*or2(n_r))-     &
                  &          two*dLvisc(n_r)*or1(n_r)*(dLvisc(n_r)-or1(n_r))+      &
                  &          two*(two*or1(n_r)+beta(n_r)-dLvisc(n_r))*or2(n_r)*dL) &
                  &                                   *                dw(lm,n_r)  &
                  &       + dL*or2(n_r)* ( two*dbeta(n_r)+ddLvisc(n_r)+            &
                  &         dLvisc(n_r)*dLvisc(n_r)-two*third*beta(n_r)*beta(n_r)+ &
                  &         dLvisc(n_r)*beta(n_r)+two*or1(n_r)*(two*dLvisc(n_r)-   &
                  &         beta(n_r)-three*or1(n_r))+dL*or2(n_r) ) *   w(lm,n_r) )

                  Buo(lm) = zero
                  if ( l_heat ) Buo(lm) = BuoFac*dL*or2(n_r)*rgrav(n_r)*s(lm,n_r)
                  if ( l_chemical_conv ) Buo(lm) = Buo(lm)+ChemFac*dL*or2(n_r)*&
                  &                                rgrav(n_r)*xi(lm,n_r)

                  dwdt%impl(lm,n_r,1)=Dif(lm)+Buo(lm)

                  if ( lPressNext ) then
                     ! In the double curl formulation, we can estimate the pressure
                     ! if required.
                     p(lm,n_r)=-r(n_r)*r(n_r)/dL*                 dp_expl(lm,n_r)  &
                     &            -one/tscheme%dt(1)*(dw(lm,n_r)-dwold(lm,n_r))+   &
                     &              hdif_V(l)*visc(n_r)* ( work_LMdist(lm,n_r)     &
                     &                       - (beta(n_r)-dLvisc(n_r))*ddw(lm,n_r) &
                     &            - ( dL*or2(n_r)+dLvisc(n_r)*beta(n_r)+dbeta(n_r) &
                     &                  + two*(dLvisc(n_r)+beta(n_r))*or1(n_r)     &
                     &                                              ) * dw(lm,n_r) &
                     &             + dL*or2(n_r)*(two*or1(n_r)+two*third*beta(n_r) &
                     &                     +dLvisc(n_r) )   *            w(lm,n_r) )
                  end if

                  if ( lRmsNext ) then
                     !-- In case RMS force balance is required, one needs to also
                     !-- compute the classical diffusivity that is used in the non
                     !-- double-curl version
                     Dif(lm) =  hdif_V(l)*dL*or2(n_r)*visc(n_r) *  ( ddw(lm,n_r)   &
                     &        +(two*dLvisc(n_r)-third*beta(n_r))*     dw(lm,n_r)   &
                     &        -( dL*or2(n_r)+four*third*( dbeta(n_r)+dLvisc(n_r)*  &
                     &           beta(n_r)+(three*dLvisc(n_r)+beta(n_r))*or1(n_r)))&
                     &                                         *       w(lm,n_r) )
                  end if
               end do
               if ( lRmsNext ) then
                  call hInt2Pol(Dif,1,n_mlo_loc,n_r,DifPolLMr(:,n_r), &
                       &        DifPol2hInt(:,n_r))
               end if
            end do
            !$omp end do

         end if

         ! In case pressure is needed in the double curl formulation
         ! we also have to compute the radial derivative of p
         if ( lPressNext .and. l_double_curl ) then
            call get_dr( p, dp, n_mlo_loc, start_lm, stop_lm, n_r_max, rscheme_oc)
            !$omp barrier
         end if

      else

         !-- Now get the poloidal from the assembly
         !$omp do private(n_r,lm,l,dL,m)
         do n_r=2,n_r_max-1
            do lm=1,n_mlo_loc
               l = map_mlo%i2l(lm)
               if ( l == 0 ) cycle
               m = map_mlo%i2m(lm)
               dL = real(l*(l+1),cp)
               if ( m == 0 ) then
                  w(lm,n_r) = r(n_r)*r(n_r)/dL*cmplx(real(work_LMdist(lm,n_r)),0.0_cp,cp)
                  dw(lm,n_r)=-r(n_r)*r(n_r)/dL*cmplx(real(ddw(lm,n_r)),0.0_cp,cp)
               else
                  w(lm,n_r) = r(n_r)*r(n_r)/dL*work_LMdist(lm,n_r)
                  dw(lm,n_r)=-r(n_r)*r(n_r)/dL*ddw(lm,n_r)
               end if
            end do
         end do
         !$omp end do

         !-- Non-penetration: u_r=0 -> w_lm=0 on both boundaries
         !$omp do private(lm)
         do lm=1,n_mlo_loc
            w(lm,1)      =zero
            w(lm,n_r_max)=zero
         end do
         !$omp end do

         !-- Other boundary condition: stress-free or rigid
         if ( l_full_sphere ) then
            if ( ktopv == 1 ) then ! Stress-free
               fac_top=-two*or1(1)-beta(1)
               !$omp do private(lm,l)
               do lm=1,n_mlo_loc
                  l = map_mlo%i2l(lm)
                  if ( l == 0 ) cycle
                  if ( l == 1 ) then
                     call rscheme_oc%robin_bc(one, fac_top, zero, 0.0_cp, one, zero, dw(lm,:))
                  else
                     call rscheme_oc%robin_bc(one, fac_top, zero, one, 0.0_cp, zero, dw(lm,:))
                  end if
               end do
               !$omp end do
            else
               !$omp do private(lm,l)
               do lm=1,n_mlo_loc
                  l = map_mlo%i2l(lm)
                  if ( l == 0 ) cycle
                  if ( l == 1 ) then
                     dw(lm,1)      =zero
                     dw(lm,n_r_max)=zero
                  else
                     call rscheme_oc%robin_bc(0.0_cp, one, zero, one, 0.0_cp, zero, dw(lm,:))
                  end if
               end do
               !$omp end do
            end if
         else ! Spherical shell
            if ( ktopv /= 1 .and. kbotv /= 1 ) then ! Rigid at both boundaries
               !$omp do private(lm)
               do lm=1,n_mlo_loc
                  dw(lm,1)      =zero
                  dw(lm,n_r_max)=zero
               end do
               !$omp end do
            else if ( ktopv /= 1 .and. kbotv == 1 ) then ! Rigid top/Stress-free bottom
               fac_bot=-two*or1(n_r_max)-beta(n_r_max)
               !$omp do private(lm,l)
               do lm=1,n_mlo_loc
                  l = map_mlo%i2l(lm)
                  if ( l == 0 ) cycle
                  call rscheme_oc%robin_bc(0.0_cp, one, zero, one, fac_bot, zero, dw(lm,:))
               end do
               !$omp end do
            else if ( ktopv == 1 .and. kbotv /= 1 ) then ! Rigid bottom/Stress-free top
               fac_top=-two*or1(1)-beta(1)
               !$omp do private(lm,l)
               do lm=1,n_mlo_loc
                  l = map_mlo%i2l(lm)
                  if ( l == 0 ) cycle
                  call rscheme_oc%robin_bc(one, fac_top, zero, 0.0_cp, one, zero, dw(lm,:))
               end do
               !$omp end do
            else if ( ktopv == 1 .and. kbotv == 1 ) then ! Stress-free at both boundaries
               fac_bot=-two*or1(n_r_max)-beta(n_r_max)
               fac_top=-two*or1(1)-beta(1)
               !$omp do private(lm,l)
               do lm=1,n_mlo_loc
                  l = map_mlo%i2l(lm)
                  if ( l == 0 ) cycle
                  call rscheme_oc%robin_bc(one, fac_top, zero, one, fac_bot, zero, dw(lm,:))
               end do
               !$omp end do
            end if
         end if

         !$omp single
         call dct_counter%start_count()
         !$omp end single
         call get_ddr( dw, ddw, work_LMdist, n_mlo_loc, start_lm, stop_lm, &
              &        n_r_max, rscheme_oc)
         !$omp barrier
         !$omp single
         call dct_counter%stop_count()
         !$omp end single

         !$omp do private(n_r,lm,l,dL)
         do n_r=2,n_r_max-1
            do lm=1,n_mlo_loc
               l = map_mlo%i2l(lm)
               if ( l == 0 )  cycle
               dL = real(l*(l+1),cp)
               dwdt%old(lm,n_r,1)= dL*or2(n_r)*w(lm,n_r)
               dpdt%old(lm,n_r,1)=-dL*or2(n_r)*dw(lm,n_r)
            end do
         end do
         !$omp end do

         if ( tscheme%l_imp_calc_rhs(1) .or. lRmsNext ) then
            if ( lRmsNext ) then
               n_r_top=n_r_cmb
               n_r_bot=n_r_icb
            else
               n_r_top=n_r_cmb+1
               n_r_bot=n_r_icb-1
            end if

            !$omp do private(n_r,lm,l,dL)
            do n_r=n_r_top,n_r_bot
               do lm=1,n_mlo_loc
                  l=map_mlo%i2l(lm)
                  if ( l == 0 ) cycle
                  dL=real(l*(l+1),cp)

                  Dif(lm) =  hdif_V(l)*dL*or2(n_r)*visc(n_r)*(       ddw(lm,n_r)  &
                  &        +(two*dLvisc(n_r)-third*beta(n_r))*        dw(lm,n_r)  &
                  &        -( dL*or2(n_r)+four*third*( dbeta(n_r)+dLvisc(n_r)*    &
                  &          beta(n_r)+(three*dLvisc(n_r)+beta(n_r))*or1(n_r)) )* &
                  &                                                   w(lm,n_r)  )
                  Buo(lm) = zero
                  if ( l_heat )  Buo(lm) = BuoFac*rho0(n_r)*rgrav(n_r)*s(lm,n_r)
                  if ( l_chemical_conv ) Buo(lm) = Buo(lm)+ChemFac*rho0(n_r)* &
                  &                                rgrav(n_r)*xi(lm,n_r)
                  dwdt%impl(lm,n_r,1)=Dif(lm)+Buo(lm)
                  dpdt%impl(lm,n_r,1)= hdif_V(l)*visc(n_r)*dL*or2(n_r)*         &
                  &                                      ( -work_LMdist(lm,n_r) &
                  &                       + (beta(n_r)-dLvisc(n_r))*ddw(lm,n_r) &
                  &            + ( dL*or2(n_r)+dLvisc(n_r)*beta(n_r)+dbeta(n_r) &
                  &                  + two*(dLvisc(n_r)+beta(n_r))*or1(n_r)     &
                  &                                           ) *    dw(lm,n_r) &
                  &            - dL*or2(n_r)* ( two*or1(n_r)+two*third*beta(n_r)&
                  &                     +dLvisc(n_r) )   *            w(lm,n_r) )
               end do

               if ( lRmsNext ) then
                  call hInt2Pol(Dif,1,n_mlo_loc,n_r,DifPolLMr(:,n_r), &
                       &        DifPol2hInt(:,n_r))
               end if
            end do
            !$omp end do
         end if

      end if
      !$omp end parallel

   end subroutine assemble_pol
!------------------------------------------------------------------------------
   subroutine get_wpMat(tscheme,l,hdif,wpMat,wpMat_fac)
      !
      !  Purpose of this subroutine is to contruct the time step matrix
      !  ``wpMat`` for the Navier-Stokes equation.
      !

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme ! time scheme
      real(cp),            intent(in) :: hdif    ! hyperdiffusion
      integer,             intent(in) :: l       ! degree :math:`\ell`

      !-- Output variables:
      class(type_realmat), intent(inout) :: wpMat
      real(cp), intent(out) :: wpMat_fac(2*n_r_max,2)

      !-- local variables:
      integer :: nR,nR_out,nR_p,nR_out_p
      integer :: info
      real(cp) :: dLh

      dLh =real(l*(l+1),kind=cp)

      !-- Now mode l>0

      !----- Boundary conditions, see above:
      do nR_out=1,rscheme_oc%n_max
         nR_out_p=nR_out+n_r_max

         wpMat%dat(1,nR_out)        =rscheme_oc%rnorm*rscheme_oc%rMat(1,nR_out)
         wpMat%dat(1,nR_out_p)      =0.0_cp
         wpMat%dat(n_r_max,nR_out)  =rscheme_oc%rnorm* &
         &                       rscheme_oc%rMat(n_r_max,nR_out)
         wpMat%dat(n_r_max,nR_out_p)=0.0_cp

         if ( ktopv == 1 ) then  ! free slip !
            wpMat%dat(n_r_max+1,nR_out)= rscheme_oc%rnorm * (        &
            &                        rscheme_oc%d2rMat(1,nR_out) -   &
            &    (two*or1(1)+beta(1))*rscheme_oc%drMat(1,nR_out) )
         else                    ! no slip, note exception for l=1,m=0
            wpMat%dat(n_r_max+1,nR_out)=rscheme_oc%rnorm*    &
            &                       rscheme_oc%drMat(1,nR_out)
         end if
         wpMat%dat(n_r_max+1,nR_out_p)=0.0_cp

         if ( l_full_sphere ) then
            if ( l == 1 ) then
               wpMat%dat(2*n_r_max,nR_out)=rscheme_oc%rnorm * &
               &                       rscheme_oc%drMat(n_r_max,nR_out)
            else
               wpMat%dat(2*n_r_max,nR_out)=rscheme_oc%rnorm * &
               &                       rscheme_oc%d2rMat(n_r_max,nR_out)
            end if
         else
            if ( kbotv == 1 ) then  ! free slip !
               wpMat%dat(2*n_r_max,nR_out)=rscheme_oc%rnorm * (       &
               &                  rscheme_oc%d2rMat(n_r_max,nR_out) - &
               &      ( two*or1(n_r_max)+beta(n_r_max))*              &
               &                  rscheme_oc%drMat(n_r_max,nR_out))
            else                 ! no slip, note exception for l=1,m=0
               wpMat%dat(2*n_r_max,nR_out)=rscheme_oc%rnorm * &
               &                       rscheme_oc%drMat(n_r_max,nR_out)
            end if
         end if
         wpMat%dat(2*n_r_max,nR_out_p)=0.0_cp

      end do   !  loop over nR_out

      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme_oc%n_max+1,n_r_max
            nR_out_p=nR_out+n_r_max
            wpMat%dat(1,nR_out)          =0.0_cp
            wpMat%dat(n_r_max,nR_out)    =0.0_cp
            wpMat%dat(n_r_max+1,nR_out)  =0.0_cp
            wpMat%dat(2*n_r_max,nR_out)  =0.0_cp
            wpMat%dat(1,nR_out_p)        =0.0_cp
            wpMat%dat(n_r_max,nR_out_p)  =0.0_cp
            wpMat%dat(n_r_max+1,nR_out_p)=0.0_cp
            wpMat%dat(2*n_r_max,nR_out_p)=0.0_cp
         end do
      end if

      !----- Other points:
      do nR_out=1,n_r_max
         nR_out_p=nR_out+n_r_max
         do nR=2,n_r_max-1
            ! in the BM2 case: visc=1.0,beta=0.0,dLvisc=0.0,dbeta=0.0
            nR_p=nR+n_r_max
            wpMat%dat(nR,nR_out)= rscheme_oc%rnorm *  (                     &
            &               dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out)          &
            &         - tscheme%wimp_lin(1)*hdif*visc(nR)*dLh*or2(nR) * (   &
            &                              rscheme_oc%d2rMat(nR,nR_out)     &
            &        +(two*dLvisc(nR)-third*beta(nR))*                      &
            &                               rscheme_oc%drMat(nR,nR_out)     &
            &        -( dLh*or2(nR)+four*third*( dLvisc(nR)*beta(nR)        &
            &          +(three*dLvisc(nR)+beta(nR))*or1(nR)+dbeta(nR) )     &
            &          )                    *rscheme_oc%rMat(nR,nR_out)  )  )

            wpMat%dat(nR,nR_out_p)= rscheme_oc%rnorm*tscheme%wimp_lin(1)*(  &
            &                            rscheme_oc%drMat(nR,nR_out)        &
            &                  -beta(nR)* rscheme_oc%rMat(nR,nR_out))
            ! the following part gives sometimes very large
            ! matrix entries
            wpMat%dat(nR_p,nR_out)=rscheme_oc%rnorm * (                       &
            &                  -dLh*or2(nR)*rscheme_oc%drMat(nR,nR_out)       &
            &         -tscheme%wimp_lin(1)*hdif*visc(nR)*dLh*or2(nR)      *(  &
            &                                 - rscheme_oc%d3rMat(nR,nR_out)  &
            &          +( beta(nR)-dLvisc(nR) )*rscheme_oc%d2rMat(nR,nR_out)  &
            &          +( dLh*or2(nR)+dbeta(nR)+dLvisc(nR)*beta(nR)           &
            &          +two*(dLvisc(nR)+beta(nR))*or1(nR) )*                  &
            &                                    rscheme_oc%drMat(nR,nR_out)  &
            &          -dLh*or2(nR)*( two*or1(nR)+dLvisc(nR)                  &
            &          +two*third*beta(nR)   )*   rscheme_oc%rMat(nR,nR_out) ) )

            wpMat%dat(nR_p,nR_out_p)= -rscheme_oc%rnorm*tscheme%wimp_lin(1)*  &
            &                          dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out)
         end do
      end do

      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         nR_p=nR+n_r_max
         wpMat%dat(nR,1)          =rscheme_oc%boundary_fac*wpMat%dat(nR,1)
         wpMat%dat(nR,n_r_max)    =rscheme_oc%boundary_fac*wpMat%dat(nR,n_r_max)
         wpMat%dat(nR,n_r_max+1)  =rscheme_oc%boundary_fac*wpMat%dat(nR,n_r_max+1)
         wpMat%dat(nR,2*n_r_max)  =rscheme_oc%boundary_fac*wpMat%dat(nR,2*n_r_max)
         wpMat%dat(nR_p,1)        =rscheme_oc%boundary_fac*wpMat%dat(nR_p,1)
         wpMat%dat(nR_p,n_r_max)  =rscheme_oc%boundary_fac*wpMat%dat(nR_p,n_r_max)
         wpMat%dat(nR_p,n_r_max+1)=rscheme_oc%boundary_fac*wpMat%dat(nR_p,n_r_max+1)
         wpMat%dat(nR_p,2*n_r_max)=rscheme_oc%boundary_fac*wpMat%dat(nR_p,2*n_r_max)
      end do

      ! compute the linesum of each line
      do nR=1,2*n_r_max
         wpMat_fac(nR,1)=one/maxval(abs(wpMat%dat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,2*n_r_max
         wpMat%dat(nR,:) = wpMat%dat(nR,:)*wpMat_fac(nR,1)
      end do

      ! also compute the rowsum of each column
      do nR=1,2*n_r_max
         wpMat_fac(nR,2)=one/maxval(abs(wpMat%dat(:,nR)))
      end do
      ! now divide each row by the rowsum
      do nR=1,2*n_r_max
         wpMat%dat(:,nR) = wpMat%dat(:,nR)*wpMat_fac(nR,2)
      end do

#ifdef MATRIX_CHECK
      block

      integer ::ipiv(2*n_r_max),iwork(2*n_r_max),i,j
      real(cp) :: work(8*n_r_max),anorm,linesum,rcond
      real(cp) :: temp_wpMat(2*n_r_max,2*n_r_max)
      integer, save :: counter=0
      integer :: filehandle
      character(len=100) :: filename
      logical :: first_run=.true.

      write(filename,"(A,I3.3,A,I3.3,A)") "wpMat_",l,"_",counter,".dat"
      open(newunit=filehandle,file=trim(filename))
      counter= counter+1

      do i=1,2*n_r_max
         do j=1,2*n_r_max
            write(filehandle,"(2ES20.12,1X)",advance="no") wpMat%dat(i,j)
         end do
         write(filehandle,"(A)") ""
      end do
      close(filehandle)
      temp_wpMat=wpMat%dat
      anorm = 0.0_cp
      do i=1,2*n_r_max
         linesum = 0.0_cp
         do j=1,2*n_r_max
            linesum = linesum + abs(temp_wpMat(i,j))
         end do
         if (linesum  >  anorm) anorm=linesum
      end do
      write(*,"(A,ES20.12)") "anorm = ",anorm
      ! LU factorization
      call dgetrf(2*n_r_max,2*n_r_max,temp_wpMat,2*n_r_max,ipiv,info)
      ! estimate the condition number
      call dgecon('I',2*n_r_max,temp_wpMat,2*n_r_max,anorm,rcond,work,iwork,info)
      write(*,"(A,I3,A,ES11.3)") "inverse condition number of wpMat for l=",l," is ",rcond

      end block
#endif

      call wpMat%prepare(info)
      if ( info /= 0 ) call abortRun('Singular matrix wpMat!')

   end subroutine get_wpMat
!------------------------------------------------------------------------------
   subroutine get_elliptic_mat(l, ellMat)
      !
      !  Purpose of this subroutine is to contruct the matrix needed
      !  for the derivation of w for the time advance of the poloidal equation
      !  if the double curl form is used.
      !

      !-- Input variables:
      integer, intent(in) :: l       ! degree :math:`\ell`

      !-- Output variables:
      class(type_realmat), intent(inout) :: ellMat

      !-- local variables:
      integer :: nR, nR_out
      integer :: info
      real(cp) :: dLh
      real(cp) :: dat(n_r_max,n_r_max)

      dLh =real(l*(l+1),kind=cp)

      !----- Boundary conditions:
      dat(1,:)      =rscheme_oc%rnorm*rscheme_oc%rMat(1,:)
      dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
      if ( ktopv == 1 ) then
         dat(2,:)=rscheme_oc%rnorm*(rscheme_oc%d2rMat(1,:)-(beta(1)+two*or1(1))* &
         &                          rscheme_oc%drMat(1,:))
      else
         dat(2,:)=rscheme_oc%rnorm*rscheme_oc%drMat(1,:)
      end if
      if ( kbotv == 1 ) then
         dat(n_r_max-1,:)=rscheme_oc%rnorm*(rscheme_oc%d2rMat(n_r_max,:)-&
         &                (beta(n_r_max)+two*or1(n_r_max))*rscheme_oc%drMat(n_r_max,:))
      else
         dat(n_r_max-1,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
      end if

      !----- Bulk points:
      do nR_out=1,n_r_max
         do nR=3,n_r_max-2
            dat(nR,nR_out)= -rscheme_oc%rnorm*dLh*orho1(nR)*or2(nR)* (   &
            &                              rscheme_oc%d2rMat(nR,nR_out)  &
            &                     -beta(nR)*rscheme_oc%drMat(nR,nR_out)  &
            &                   -dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out) )
         end do
      end do

      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         dat(nR,1)      =rscheme_oc%boundary_fac*dat(nR,1)
         dat(nR,n_r_max)=rscheme_oc%boundary_fac*dat(nR,n_r_max)
      end do

      !-- Array copy
      call ellMat%set_data(dat)

      !-- Lu factorisation
      call ellMat%prepare(info)
      if ( info /= 0 ) call abortRun('Singular matrix ellMat!')

   end subroutine get_elliptic_mat
!-----------------------------------------------------------------------------
   subroutine get_wMat(tscheme,l,hdif,wMat,wMat_fac)
      !
      !  Purpose of this subroutine is to contruct the time step matrix
      !  wpmat  for the NS equation. This matrix corresponds here to the
      !  radial component of the double-curl of the Navier-Stokes equation.
      !

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme  ! time scheme
      real(cp),            intent(in) :: hdif     ! hyperdiffusion
      integer,             intent(in) :: l        ! degree :math:`\ell`

      !-- Output variables:
      class(type_realmat), intent(inout) :: wMat
      real(cp), intent(out) :: wMat_fac(n_r_max,2)

      !-- local variables:
      integer :: nR, nR_out
      integer :: info
      real(cp) :: dLh
      real(cp) :: dat(n_r_max,n_r_max)

      dLh =real(l*(l+1),kind=cp)

      !----- Boundary conditions:
      !-- Non-penetration condition at both boundaries
      dat(1,:)      =rscheme_oc%rnorm*rscheme_oc%rMat(1,:)
      dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)

      !-- Second line and n_r_max-1 lines are used for the second BCs
      if ( ktopv == 1 ) then  ! free slip
         dat(2,:)=rscheme_oc%rnorm *(rscheme_oc%d2rMat(1,:)-   &
         &               (two*or1(1)+beta(1))*rscheme_oc%drMat(1,:) )
      else                    ! no slip
         dat(2,:)=rscheme_oc%rnorm*rscheme_oc%drMat(1,:)
      end if

      if ( l_full_sphere ) then
         if ( l == 1 ) then
            dat(n_r_max-1,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         else
            dat(n_r_max-1,:)=rscheme_oc%rnorm*rscheme_oc%d2rMat(n_r_max,:)
         end if
      else
         if ( kbotv == 1 ) then  ! free slip
            dat(n_r_max-1,:)=rscheme_oc%rnorm *(rscheme_oc%d2rMat(n_r_max,:)-  &
            &                      (two*or1(n_r_max)+beta(n_r_max))*           &
            &                                    rscheme_oc%drMat(n_r_max,:) )
         else                 ! no slip
            dat(n_r_max-1,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         end if
      end if

      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme_oc%n_max+1,n_r_max
            dat(1,nR_out)        =0.0_cp
            dat(2,nR_out)        =0.0_cp
            dat(n_r_max-1,nR_out)=0.0_cp
            dat(n_r_max,nR_out)  =0.0_cp
         end do
      end if

      !----- Bulk points:
      do nR_out=1,n_r_max
         do nR=3,n_r_max-2
            dat(nR,nR_out)=rscheme_oc%rnorm* (                              &
            &          -dLh*or2(nR)*orho1(nR)*(                             &
            &                              rscheme_oc%d2rMat(nR,nR_out)     &
            &                     -beta(nR)*rscheme_oc%drMat(nR,nR_out) -   &
            &                    dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out) )   &
            &  + tscheme%wimp_lin(1)*orho1(nR)*hdif*visc(nR)*dLh*or2(nR) * (&
            &                               rscheme_oc%d4rMat(nR,nR_out)    &
            &   +two*(dLvisc(nR)-beta(nR))* rscheme_oc%d3rMat(nR,nR_out)    &
            &    +( ddLvisc(nR)-two*dbeta(nR)+dLvisc(nR)*dLvisc(nR)+        &
            &       beta(nR)*beta(nR)-three*dLvisc(nR)*beta(nR)-            &
            &       two*or1(nR)*(dLvisc(nR)+beta(nR))-two*dLh*or2(nR) ) *   &
            &                               rscheme_oc%d2rMat(nR,nR_out)    &
            &    +( -ddbeta(nR)-dbeta(nR)*(two*dLvisc(nR)-beta(nR)+         &
            &       two*or1(nR))-ddLvisc(nR)*(beta(nR)+two*or1(nR))+        &
            &       beta(nR)*beta(nR)*(dLvisc(nR)+two*or1(nR))-beta(nR)*    &
            &       (dLvisc(nR)*dLvisc(nR)-two*or2(nR))-two*dLvisc(nR)*     &
            &       or1(nR)*(dLvisc(nR)-or1(nR))+two*(two*or1(nR)+          &
            &       beta(nR)-dLvisc(nR))*dLh*or2(nR) ) *                    &
            &                                rscheme_oc%drMat(nR,nR_out)    &
            &    + dLh*or2(nR)*( two*dbeta(nR)+ddLvisc(nR)+dLvisc(nR)*      &
            &      dLvisc(nR)-two*third*beta(nR)*beta(nR)+dLvisc(nR)*       &
            &      beta(nR)+two*or1(nR)*(two*dLvisc(nR)-beta(nR)-three*     &
            &      or1(nR) ) + dLh*or2(nR) ) *                              &
            &                                rscheme_oc%rMat(nR,nR_out) ) )
         end do
      end do

      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         dat(nR,1)      =rscheme_oc%boundary_fac*dat(nR,1)
         dat(nR,n_r_max)=rscheme_oc%boundary_fac*dat(nR,n_r_max)
      end do

      ! compute the linesum of each line
      do nR=1,n_r_max
         wMat_fac(nR,1)=one/maxval(abs(dat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,n_r_max
         dat(nR,:) = dat(nR,:)*wMat_fac(nR,1)
      end do

      ! also compute the rowsum of each column
      do nR=1,n_r_max
         wMat_fac(nR,2)=one/maxval(abs(dat(:,nR)))
      end do
      ! now divide each row by the rowsum
      do nR=1,n_r_max
         dat(:,nR) = dat(:,nR)*wMat_fac(nR,2)
      end do

      !-- Array copy
      call wMat%set_data(dat)

      call wMat%prepare(info)

      if ( info /= 0 ) call abortRun('Singular matrix wMat!')

   end subroutine get_wMat
!-----------------------------------------------------------------------------
   subroutine get_p0Mat(pMat)
      !
      ! This subroutine solves the linear problem of the spherically-symmetric
      ! pressure
      !

      !-- Output variables:
      class(type_realmat), intent(inout) :: pMat ! matrix

      !-- Local variables:
      real(cp) :: dat(n_r_max,n_r_max), delr
      integer :: info, nCheb, nR_out, nR, nCheb_in

      !-- Bulk points
      do nR_out=1,n_r_max
         do nR=2,n_r_max
            dat(nR,nR_out)= rscheme_oc%rnorm * (              &
            &                    rscheme_oc%drMat(nR,nR_out)- &
            &            beta(nR)*rscheme_oc%rMat(nR,nR_out) )
         end do
      end do

      !-- Boundary condition for spherically-symmetric pressure
      !-- The integral of rho' r^2 dr vanishes
      if ( ThExpNb*ViscHeatFac /= 0 .and. ktopp==1 ) then

         work(:) = ThExpNb*ViscHeatFac*ogrun(:)*alpha0(:)*r(:)*r(:)
         call rscheme_oc%costf1(work)
         work(:)      =work(:)*rscheme_oc%rnorm
         work(1)      =rscheme_oc%boundary_fac*work(1)
         work(n_r_max)=rscheme_oc%boundary_fac*work(n_r_max)

         if ( rscheme_oc%version == 'cheb' ) then

            do nCheb=1,rscheme_oc%n_max
               dat(1,nCheb)=0.0_cp
               do nCheb_in=1,rscheme_oc%n_max
                  if ( mod(nCheb+nCheb_in-2,2)==0 ) then
                     dat(1,nCheb)=dat(1,nCheb)+ &
                     &             ( one/(one-real(nCheb_in-nCheb,cp)**2)    + &
                     &               one/(one-real(nCheb_in+nCheb-2,cp)**2) )* &
                     &               work(nCheb_in)*half*rscheme_oc%rnorm
                  end if
               end do
            end do

         else

            !!-- In the finite differences case, we restrict the integral boundary
            !!-- condition to a trapezoidal rule of integration
            !do nR_out=2,rscheme_oc%n_max-1
            !   dat(1,nR_out)=half*work(nR_out)*( r(nR_out+1)-r(nR_out-1) )
            !end do
            !dat(1,1)=half*work(1)*(r(2)-r(1))
            dat(1,:)=rscheme_oc%rMat(1,:)
            delr = r(n_r_max)-r(n_r_max-1)

            !-- First-order on the last point
            dat(n_r_max,n_r_max)    =one/delr-beta(n_r_max)
            dat(n_r_max,n_r_max-1)  =-one/delr
            dat(n_r_max,1:n_r_max-2)=0.0_cp
         end if

      else

         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%rMat(1,:)
         if ( rscheme_oc%version == 'fd' ) then
            delr = r(n_r_max)-r(n_r_max-1)
            !-- First-order on the last point
            dat(n_r_max,n_r_max)  =one/delr-beta(n_r_max)
            dat(n_r_max,n_r_max-1)=-one/delr
         end if

      end if

      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros
         do nR_out=rscheme_oc%n_max+1,n_r_max
            dat(1,nR_out)=0.0_cp
         end do
      end if

      !----- Factors for highest and lowest cheb mode:
      do nR=1,n_r_max
         dat(nR,1)      =rscheme_oc%boundary_fac*dat(nR,1)
         dat(nR,n_r_max)=rscheme_oc%boundary_fac*dat(nR,n_r_max)
      end do

      !-- Array copy
      call pMat%set_data(dat)

      !---- LU decomposition:
      call pMat%prepare(info)
      if ( info /= 0 ) call abortRun('! Singular matrix p0Mat!')

   end subroutine get_p0Mat
!-----------------------------------------------------------------------------
end module updateWP_mod
