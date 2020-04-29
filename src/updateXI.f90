#include "perflib_preproc.cpp"
module updateXi_mod

   use omp_lib
   use precision_mod
   use truncation, only: n_r_max, lm_max, l_max, n_r_icb, n_r_cmb, &
       &                 get_openmp_blocks, n_mlo_loc, n_lo_loc
   use LMmapping, only: map_mlo
   use radial_functions, only: orho1, or1, or2, beta, rscheme_oc, r
   use physical_parameters, only: osc, kbotxi, ktopxi
   use num_param, only: dct_counter, solve_counter
   use init_fields, only: topxi, botxi
   use blocking, only: lo_map, lo_sub_map, llm, ulm
   use horizontal_data, only: hdif_Xi
   use logic, only: l_update_xi, l_finite_diff, l_full_sphere
   use parallel_mod, only: coord_r, chunksize, n_ranks_r
   use radial_der, only: get_ddr, get_dr
   use constants, only: zero, one, two
   use fields, only: work_LMloc, work_LMdist
   use mem_alloc, only: bytes_allocated
   use useful, only: abortRun
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray
   use dense_matrices
   use real_matrices
   use band_matrices


   implicit none

   private

   real(cp), allocatable :: rhs1(:,:,:)
   integer :: maxThreads
   class(type_realmat), pointer :: xiMat(:), xi0Mat
#ifdef WITH_PRECOND_S
   real(cp), allocatable :: xiMat_fac(:,:)
#endif
#ifdef WITH_PRECOND_S0
   real(cp), allocatable :: xi0Mat_fac(:)
#endif
   logical, public, allocatable :: lXimat(:)

   public :: initialize_updateXi, finalize_updateXi, updateXi, assemble_comp, &
   &         finish_exp_comp, get_comp_rhs_imp, get_comp_rhs_imp_dist

contains

   subroutine initialize_updateXi

      integer :: ll, n_bands

      if ( l_finite_diff ) then
         allocate( type_bandmat :: xiMat(n_lo_loc) )
         allocate( type_bandmat :: xi0Mat )

         if ( ktopxi == 1 .and. kbotxi == 1 .and. rscheme_oc%order == 2 &
          &   .and. rscheme_oc%order_boundary <= 2 ) then ! Fixed composition at both boundaries
            n_bands = rscheme_oc%order+1
         else
            n_bands = max(2*rscheme_oc%order_boundary+1,rscheme_oc%order+1)
         end if

         call xi0Mat%initialize(n_bands,n_r_max,l_pivot=.true.)
         do ll=1,n_lo_loc
            call xiMat(ll)%initialize(n_bands,n_r_max,l_pivot=.true.)
         end do
      else
         allocate( type_densemat :: xiMat(n_lo_loc) )
         allocate( type_densemat :: xi0Mat )

         call xi0Mat%initialize(n_r_max,n_r_max,l_pivot=.true.)
         do ll=1,n_lo_loc
            call xiMat(ll)%initialize(n_r_max,n_r_max,l_pivot=.true.)
         end do
      end if

#ifdef WITH_PRECOND_S
      allocate(xiMat_fac(n_r_max,n_lo_loc))
      bytes_allocated = bytes_allocated+n_r_max*n_lo_loc*SIZEOF_DEF_REAL
#endif
#ifdef WITH_PRECOND_S0
      allocate(xi0Mat_fac(n_r_max))
      bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_REAL
#endif
      allocate( lXimat(n_lo_loc) )
      bytes_allocated = bytes_allocated+n_lo_loc*SIZEOF_LOGICAL

#ifdef WITHOMP
      maxThreads=omp_get_max_threads()
#else
      maxThreads=1
#endif
      allocate( rhs1(n_r_max,2*maxval(map_mlo%n_mi(:)), 1) )
      bytes_allocated = bytes_allocated + n_r_max*maxval(map_mlo%n_mi(:))*&
      &                 maxThreads*SIZEOF_DEF_COMPLEX

   end subroutine initialize_updateXi
!------------------------------------------------------------------------------
   subroutine finalize_updateXI

      integer :: ll

      do ll=1,n_lo_loc
         call xiMat(ll)%finalize()
      end do
      call xi0Mat%finalize()

      deallocate( lXimat )
#ifdef WITH_PRECOND_S
      deallocate(xiMat_fac)
#endif
#ifdef WITH_PRECOND_S0
      deallocate(xi0Mat_fac)
#endif
      deallocate( rhs1 )

   end subroutine finalize_updateXI
!------------------------------------------------------------------------------
   subroutine updateXi(xi, dxi, dxidt, tscheme)
      !
      !  updates the entropy field s and its radial derivatives
      !  adds explicit part to time derivatives of s
      !

      !-- Input of variables:
      class(type_tscheme), intent(in) :: tscheme

      !-- Input/output of scalar fields:
      complex(cp),       intent(inout) :: xi(n_mlo_loc,n_r_max)
      type(type_tarray), intent(inout) :: dxidt
      !-- Output: dxi
      complex(cp),       intent(out) :: dxi(n_mlo_loc,n_r_max)

      !-- Local variables:
      integer :: l, m             ! degree and order
      integer :: lj, mi, i        ! l,m and ml counter
      integer :: nR               ! counts radial grid points
      integer :: n_r_out          ! counts cheb modes
      real(cp) ::  rhs(n_r_max) ! real RHS for l=m=0

      if ( .not. l_update_xi ) return

      !-- Now assemble the right hand side and store it in work_LMloc
      call tscheme%set_imex_rhs(work_LMdist, dxidt, 1, n_mlo_loc, n_r_max)

      call solve_counter%start_count()
      do lj=1,n_lo_loc
         l = map_mlo%lj2l(lj)

         if ( .not. lXimat(lj) ) then
            if ( l == 0 ) then
#ifdef WITH_PRECOND_S0
               call get_xi0Mat(tscheme,xi0Mat,xi0Mat_fac)
#else
               call get_xi0Mat(tscheme,xi0Mat)
#endif
            else
#ifdef WITH_PRECOND_S
               call get_xiMat(tscheme,l,hdif_Xi(l), &
                    &         xiMat(lj),xiMat_fac(:,lj))
#else
               call get_xiMat(tscheme,l,hdif_Xi(l),xiMat(lj))
#endif
            end if
            lXimat(lj)=.true.
         end if

         ! Loop over m corresponding to current l
         do mi=1,map_mlo%n_mi(lj)
            m = map_mlo%milj2m(mi,lj)
            i = map_mlo%milj2i(mi,lj)

            if ( l==0 ) then

               rhs(1)      =real(topxi(0,0))
               do nR=2,n_r_max-1
                  rhs(nR)=real(work_LMdist(i,nR))
               end do
               rhs(n_r_max)=real(botxi(0,0))
#ifdef WITH_PRECOND_S0
               rhs(:) = xi0Mat_fac(:)*rhs(:)
#endif
               call xi0Mat%solve(rhs)

            else ! l  /=  0

               ! Builds real RHS
               rhs1(1,2*mi-1,1)    =real(topxi(l,m))
               do nR=2,n_r_max-1
                  rhs1(nR,2*mi-1,1)=real(work_LMdist(i,nR))
               end do
               rhs1(n_r_max,2*mi-1,1)=real(botxi(l,m))

               ! Builds imag RHS
               rhs1(1,2*mi,1)    =aimag(topxi(l,m))
               do nR=2,n_r_max-1
                  rhs1(nR,2*mi,1)=aimag(work_LMdist(i,nR))
               end do
               rhs1(n_r_max,2*mi,1)=aimag(botxi(l,m))

#ifdef WITH_PRECOND_S
               rhs1(:,2*mi-1,1)=xiMat_fac(:,lj)*rhs1(:,2*mi-1,1)
               rhs1(:,2*mi,1)  =xiMat_fac(:,lj)*rhs1(:,2*mi,1)
#endif
            end if
         end do ! loop over m

         ! Solve for all RHS at once
         if ( l>0 ) then
            call xiMat(lj)%solve(rhs1(:,1:2*map_mlo%n_mi(lj),1), &
                 &                 2*map_mlo%n_mi(lj))
         end if

         ! Loop over m corresponding to current l (again)
         do mi=1,map_mlo%n_mi(lj)
            m = map_mlo%milj2m(mi,lj)
            i = map_mlo%milj2i(mi,lj)

            if ( l == 0 ) then
               do n_r_out=1,rscheme_oc%n_max
                  xi(i,n_r_out)=rhs(n_r_out)
               end do
            else
               if ( m > 0 ) then
                  do n_r_out=1,rscheme_oc%n_max
                     xi(i,n_r_out)= cmplx(rhs1(n_r_out,2*mi-1,1), &
                     &                    rhs1(n_r_out,2*mi,1),kind=cp)
                  end do
               else
                  do n_r_out=1,rscheme_oc%n_max
                     xi(i,n_r_out)= cmplx(rhs1(n_r_out,2*mi-1,1),0.0_cp,kind=cp)
                  end do
               end if
            end if
         end do  ! loop over m (again)
      end do     ! loop over l

      call solve_counter%stop_count(l_increment=.false.)

      !-- set cheb modes > rscheme_oc%n_max to zero (dealiazing)
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         do i=1,n_mlo_loc
            xi(i,n_r_out)=zero
         end do
      end do

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dxidt, 1, n_mlo_loc, n_r_max)

      !-- Calculation of the implicit part
      if ( tscheme%istage == tscheme%nstages ) then
         call get_comp_rhs_imp_dist(xi, dxi, dxidt, 1, tscheme%l_imp_calc_rhs(1),  & 
              &                l_in_cheb_space=.true.)
      else
         call get_comp_rhs_imp_dist(xi, dxi, dxidt, tscheme%istage+1,          &
              &                tscheme%l_imp_calc_rhs(tscheme%istage+1),  &
              &                l_in_cheb_space=.true.)
      end if

   end subroutine updateXi
!------------------------------------------------------------------------------
   subroutine finish_exp_comp(dVXirLM, dxi_exp_last)

      !-- Input variables
      complex(cp), intent(inout) :: dVXirLM(n_mlo_loc,n_r_max)

      !-- Output variables
      complex(cp), intent(inout) :: dxi_exp_last(n_mlo_loc,n_r_max)

      !-- Local variables
      integer :: n_r, start_lm, stop_lm

      !$omp parallel default(shared) private(start_lm, stop_lm)
      start_lm=1; stop_lm=n_mlo_loc
      call get_openmp_blocks(start_lm,stop_lm)
      call get_dr( dVXirLM, work_LMdist, n_mlo_loc, start_lm,  &
           &       stop_lm, n_r_max, rscheme_oc, nocopy=.true. )
      !$omp barrier

      !$omp do
      do n_r=1,n_r_max
         dxi_exp_last(:,n_r)=orho1(n_r)*( dxi_exp_last(:,n_r)-   &
         &                        or2(n_r)*work_LMdist(:,n_r) )
      end do
      !$omp end do
      !$omp end parallel

   end subroutine finish_exp_comp
!------------------------------------------------------------------------------
   subroutine get_comp_rhs_imp(xi, dxi, dxidt, istage, l_calc_lin, l_in_cheb_space)

      !-- Input variables
      integer,             intent(in) :: istage
      logical,             intent(in) :: l_calc_lin
      logical, optional,   intent(in) :: l_in_cheb_space

      !-- Output variable
      complex(cp),       intent(inout) :: xi(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: dxi(llm:ulm,n_r_max)
      type(type_tarray), intent(inout) :: dxidt

      !-- Local variables
      logical :: l_in_cheb
      integer :: n_r, lm, start_lm, stop_lm, l1
      real(cp) :: dL
      integer, pointer :: lm2l(:),lm2m(:)

      if ( present(l_in_cheb_space) ) then
         l_in_cheb = l_in_cheb_space
      else
         l_in_cheb = .false.
      end if

      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      !$omp parallel default(shared)  private(start_lm, stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)

      !$omp single
      call dct_counter%start_count()
      !$omp end single
      call get_ddr(xi, dxi, work_LMloc, ulm-llm+1,start_lm-llm+1,  &
           &       stop_lm-llm+1,n_r_max, rscheme_oc, l_dct_in=.not. l_in_cheb)
      if ( l_in_cheb ) call rscheme_oc%costf1(xi,ulm-llm+1,start_lm-llm+1, &
                            &                 stop_lm-llm+1)
      !$omp barrier
      !$omp single
      call dct_counter%stop_count(l_increment=.false.)
      !$omp end single

      if ( istage == 1 ) then
         !$omp do
         do n_r=1,n_r_max
            dxidt%old(:,n_r,istage) = xi(:,n_r)
         end do
         !$omp end do
      end if

      if ( l_calc_lin ) then

         !$omp do private(n_r,lm,l1,dL)
         do n_r=1,n_r_max
            do lm=llm,ulm
               l1 = lm2l(lm)
               dL = real(l1*(l1+1),cp)
               dxidt%impl(lm,n_r,istage)=                    osc*hdif_Xi(l1) *   &
               &     ( work_LMloc(lm,n_r)+(beta(n_r)+two*or1(n_r)) * dxi(lm,n_r) &
               &                                       - dL*or2(n_r)* xi(lm,n_r) )
            end do
         end do
         !$omp end do

      end if

      !$omp end parallel

   end subroutine get_comp_rhs_imp
!------------------------------------------------------------------------------
   subroutine get_comp_rhs_imp_dist(xi, dxi, dxidt, istage, l_calc_lin, l_in_cheb_space)

      !-- Input variables
      integer,             intent(in) :: istage
      logical,             intent(in) :: l_calc_lin
      logical, optional,   intent(in) :: l_in_cheb_space

      !-- Output variable
      complex(cp),       intent(inout) :: xi(n_mlo_loc,n_r_max)
      complex(cp),       intent(out) :: dxi(n_mlo_loc,n_r_max)
      type(type_tarray), intent(inout) :: dxidt

      !-- Local variables
      logical :: l_in_cheb
      integer :: n_r, lm, start_lm, stop_lm, l1
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
      call get_ddr(xi, dxi, work_LMdist, n_mlo_loc, start_lm,  &
           &       stop_lm, n_r_max, rscheme_oc, l_dct_in=.not. l_in_cheb)
      if ( l_in_cheb ) call rscheme_oc%costf1(xi, n_mlo_loc, start_lm, stop_lm)
      !$omp barrier
      !$omp single
      call dct_counter%stop_count(l_increment=.false.)
      !$omp end single

      if ( istage == 1 ) then
         !$omp do
         do n_r=1,n_r_max
            dxidt%old(:,n_r,istage) = xi(:,n_r)
         end do
         !$omp end do
      end if

      if ( l_calc_lin ) then

         !$omp do private(n_r,lm,l1,dL)
         do n_r=1,n_r_max
            do lm=1,n_mlo_loc
               l1 = map_mlo%i2l(lm)
               dL = real(l1*(l1+1),cp)
               dxidt%impl(lm,n_r,istage)=                    osc*hdif_Xi(l1) *   &
               &    ( work_LMdist(lm,n_r)+(beta(n_r)+two*or1(n_r)) * dxi(lm,n_r) &
               &                                       - dL*or2(n_r)* xi(lm,n_r) )
            end do
         end do
         !$omp end do

      end if

      !$omp end parallel

   end subroutine get_comp_rhs_imp_dist
!------------------------------------------------------------------------------
   subroutine assemble_comp(xi, dxi, dxidt, tscheme)
      !
      ! This subroutine is used to assemble the chemical composition when an
      ! IMEX-RK with an assembly stage is employed. Non-Dirichlet boundary
      ! conditions are handled using Canuto (1986) approach.
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variables
      complex(cp),       intent(inout) :: xi(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: dxi(llm:ulm,n_r_max)
      type(type_tarray), intent(inout) :: dxidt

      !-- Local variables
      integer :: lm, l1, m1, n_r
      integer, pointer :: lm2l(:), lm2m(:)

      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      call tscheme%assemble_imex(work_LMloc, dxidt, llm, ulm, n_r_max)

      !$omp parallel default(shared)
      !$omp do private(n_r,lm,m1)
      do n_r=2,n_r_max
         do lm=llm,ulm
            m1 = lm2m(lm)
            if ( m1 == 0 ) then
               xi(lm,n_r)=cmplx(real(work_LMloc(lm,n_r)),0.0_cp,cp)
            else
               xi(lm,n_r)=work_LMloc(lm,n_r)
            end if
         end do
      end do
      !$omp end do

      !-- Boundary conditions
      if ( l_full_sphere) then
         if ( ktopxi == 1 ) then ! Fixed entropy at the outer boundary
            !$omp do private(lm,l1,m1)
            do lm=llm,ulm
               l1 = lm2l(lm)
               m1 = lm2m(lm)
               if ( l1 == 1 ) then
                  call rscheme_oc%robin_bc(0.0_cp, one, topxi(l1,m1), 0.0_cp, one, &
                       &                   botxi(l1,m1), xi(lm,:))
               else
                  call rscheme_oc%robin_bc(0.0_cp, one, topxi(l1,m1), one, 0.0_cp, &
                       &                   botxi(l1,m1), xi(lm,:))
               end if
            end do
            !$omp end do
         else ! Fixed flux at the outer boundary
            !$omp do private(lm,l1,m1)
            do lm=llm,ulm
               l1 = lm2l(lm)
               m1 = lm2m(lm)
               if ( l1 == 1 ) then
                  call rscheme_oc%robin_bc(one, 0.0_cp, topxi(l1,m1), 0.0_cp, one, &
                       &                   botxi(l1,m1), xi(lm,:))
               else
                  call rscheme_oc%robin_bc(one, 0.0_cp, topxi(l1,m1), one, 0.0_cp, &
                       &                   botxi(l1,m1), xi(lm,:))
               end if
            end do
            !$omp end do
         end if

      else ! Spherical shell

         !-- Boundary conditions
         if ( ktopxi==1 .and. kbotxi==1 ) then ! Dirichlet on both sides
            !$omp do private(lm,l1,m1)
            do lm=llm,ulm
               l1 = lm2l(lm)
               m1 = lm2m(lm)
               call rscheme_oc%robin_bc(0.0_cp, one, topxi(l1,m1), 0.0_cp, one, &
                    &                   botxi(l1,m1), xi(lm,:))
            end do
            !$omp end do
         else if ( ktopxi==1 .and. kbotxi /= 1 ) then ! Dirichlet: top and Neumann: bot
            !$omp do private(lm,l1,m1)
            do lm=llm,ulm
               l1 = lm2l(lm)
               m1 = lm2m(lm)
               call rscheme_oc%robin_bc(0.0_cp, one, topxi(l1,m1), one, 0.0_cp, &
                    &                   botxi(l1,m1), xi(lm,:))
            end do
            !$omp end do
         else if ( kbotxi==1 .and. ktopxi /= 1 ) then ! Dirichlet: bot and Neumann: top
            !$omp do private(lm,l1,m1)
            do lm=llm,ulm
               l1 = lm2l(lm)
               m1 = lm2m(lm)
               call rscheme_oc%robin_bc(one, 0.0_cp, topxi(l1,m1), 0.0_cp, one, &
                    &                   botxi(l1,m1), xi(lm,:))
            end do
            !$omp end do
         else if ( kbotxi /=1 .and. kbotxi /= 1 ) then ! Neumann on both sides
            !$omp do private(lm,l1,m1)
            do lm=llm,ulm
               l1 = lm2l(lm)
               m1 = lm2m(lm)
               call rscheme_oc%robin_bc(0.0_cp, one, topxi(l1,m1), 0.0_cp, one, &
                    &                   botxi(l1,m1), xi(lm,:))
            end do
            !$omp end do
         end if
      end if
      !$omp end parallel

      call get_comp_rhs_imp(xi, dxi, dxidt, 1, tscheme%l_imp_calc_rhs(1), .false.)

   end subroutine assemble_comp
!------------------------------------------------------------------------------
#ifdef WITH_PRECOND_S0
   subroutine get_xi0Mat(tscheme,xiMat,xiMat_fac)
#else
   subroutine get_xi0Mat(tscheme,xiMat)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matrix
      !  xiMat0
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme        ! time step

      !-- Output variables
      class(type_realmat), intent(inout) :: xiMat
#ifdef WITH_PRECOND_S0
      real(cp), intent(out) :: xiMat_fac(n_r_max)
#endif

      !-- Local variables:
      real(cp) :: dat(n_r_max,n_r_max)
      integer :: info, nR_out, nR

      !----- Boundary condition:
      if ( ktopxi == 1 ) then
         !--------- Constant chemical composition at CMB:
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%rMat(1,:)
      else
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%drMat(1,:)
      end if

      if ( l_full_sphere ) then
         dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
      else
         if ( kbotxi == 1 ) then
            !--------- Constant composition at ICB:
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
         else
            !--------- Constant flux at ICB:
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         end if
      end if

      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme_oc%n_max+1,n_r_max
            dat(1,nR_out)      =0.0_cp
            dat(n_r_max,nR_out)=0.0_cp
         end do
      end if

      !-- Fill bulk points
      do nR_out=1,n_r_max
         do nR=2,n_r_max-1
            dat(nR,nR_out)= rscheme_oc%rnorm * (                          &
            &                                rscheme_oc%rMat(nR,nR_out) - &
            & tscheme%wimp_lin(1)*osc*(    rscheme_oc%d2rMat(nR,nR_out) + &
            &    (beta(nR)+two*or1(nR))*    rscheme_oc%drMat(nR,nR_out) ) )
         end do
      end do

      !----- Factors for highest and lowest cheb mode:
      do nR=1,n_r_max
         dat(nR,1)      =rscheme_oc%boundary_fac*dat(nR,1)
         dat(nR,n_r_max)=rscheme_oc%boundary_fac*dat(nR,n_r_max)
      end do

#ifdef WITH_PRECOND_S0
      ! compute the linesum of each line
      do nR=1,n_r_max
         xiMat_fac(nR)=one/maxval(abs(dat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,n_r_max
         dat(nR,:) = dat(nR,:)*xiMat_fac(nR)
      end do
#endif

      !-- Array copy
      call xiMat%set_data(dat)

      !---- LU decomposition:
      call xiMat%prepare(info)
      if ( info /= 0 ) call abortRun('! Singular matrix xiMat0!')

   end subroutine get_xi0Mat
!-----------------------------------------------------------------------------
#ifdef WITH_PRECOND_S
   subroutine get_xiMat(tscheme,l,hdif,xiMat,xiMat_fac)
#else
   subroutine get_xiMat(tscheme,l,hdif,xiMat)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matrices
      !  xiMat(i,j) for the equation for the chemical composition.
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme        ! time step
      real(cp),            intent(in) :: hdif
      integer,             intent(in) :: l

      !-- Output variables
      class(type_realmat), intent(inout) :: xiMat
#ifdef WITH_PRECOND_S
      real(cp),intent(out) :: xiMat_fac(n_r_max)
#endif

      !-- Local variables:
      integer :: info, nR_out, nR
      real(cp) :: dLh
      real(cp) :: dat(n_r_max,n_r_max)

      dLh=real(l*(l+1),kind=cp)

      !----- Boundary coditions:
      if ( ktopxi == 1 ) then
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%rMat(1,:)
      else
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%drMat(1,:)
      end if

      if ( l_full_sphere ) then
         !dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         if ( l == 1 ) then
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
         else
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         end if
      else
         if ( kbotxi == 1 ) then
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
         else
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         end if
      end if

      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme_oc%n_max+1,n_r_max
            dat(1,nR_out)      =0.0_cp
            dat(n_r_max,nR_out)=0.0_cp
         end do
      end if

      !----- Bulk points
      do nR_out=1,n_r_max
         do nR=2,n_r_max-1
            dat(nR,nR_out)= rscheme_oc%rnorm * (                           &
            &                                 rscheme_oc%rMat(nR,nR_out) - &
            & tscheme%wimp_lin(1)*osc*hdif*(rscheme_oc%d2rMat(nR,nR_out) + &
            &   ( beta(nR)+two*or1(nR) )*    rscheme_oc%drMat(nR,nR_out) - &
            &        dLh*or2(nR)*             rscheme_oc%rMat(nR,nR_out) ) )
         end do
      end do

      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         dat(nR,1)      =rscheme_oc%boundary_fac*dat(nR,1)
         dat(nR,n_r_max)=rscheme_oc%boundary_fac*dat(nR,n_r_max)
      end do

#ifdef WITH_PRECOND_S
      ! compute the linesum of each line
      do nR=1,n_r_max
         xiMat_fac(nR)=one/maxval(abs(dat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,n_r_max
         dat(nR,:) = dat(nR,:)*xiMat_fac(nR)
      end do
#endif

      !-- Array copy
      call xiMat%set_data(dat)

      !----- LU decomposition:
      call xiMat%prepare(info)
      if ( info /= 0 ) call abortRun('Singular matrix xiMat!')

   end subroutine get_xiMat
!-----------------------------------------------------------------------------
end module updateXi_mod
