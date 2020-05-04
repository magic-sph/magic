#include "perflib_preproc.cpp"
module updateZ_mod

   use init_fields
   use omp_lib
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, l_max, n_r_cmb, n_r_icb, &
       &                 get_openmp_blocks, n_lo_loc, n_mlo_loc
   use LMmapping, only: map_mlo
   use radial_functions, only: visc, or1, or2, rscheme_oc, dLvisc, beta, &
       &                       rho0, r_icb, r_cmb, r, beta, dbeta
   use physical_parameters, only: kbotv, ktopv, prec_angle, po, oek
   use num_param, only: AMstart, dct_counter, solve_counter
   use torsional_oscillations, only: ddzASL
   use horizontal_data, only: hdif_V
   use logic, only: l_rot_ma, l_rot_ic, l_SRMA, l_SRIC, l_z10mat, l_precession, &
       &            l_correct_AMe, l_correct_AMz, l_update_v, l_TO,             &
       &            l_finite_diff, l_full_sphere
   use RMS, only: DifTor2hInt
   use constants, only: c_lorentz_ma, c_lorentz_ic, c_dt_z10_ma, c_dt_z10_ic, &
       &                c_moi_ma, c_moi_ic, c_z10_omega_ma, c_z10_omega_ic,   &
       &                c_moi_oc, y10_norm, y11_norm, zero, one, two, four,   &
       &                pi, third
   use parallel_mod
   use outRot, only: get_angular_moment
   use RMS_helpers, only: hInt2Tor
   use radial_der, only: get_ddr
   use fields, only: work_LMdist
   use useful, only: abortRun
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray, type_tscalar
   use special
   use dense_matrices
   use real_matrices
   use band_matrices

   implicit none

   private

   !-- Input of recycled work arrays:
   real(cp), allocatable :: rhs1(:,:,:) ! RHS for other modes
   complex(cp), allocatable :: Dif(:)
   class(type_realmat), pointer :: zMat(:), z10Mat
#ifdef WITH_PRECOND_Z
   real(cp), allocatable :: zMat_fac(:,:)
#endif
#ifdef WITH_PRECOND_Z10
   real(cp), allocatable :: z10Mat_fac(:)
#endif
   logical, public :: lZ10mat
   logical, public, allocatable :: lZmat(:)

   integer :: maxThreads

   public :: updateZ, initialize_updateZ, finalize_updateZ, get_tor_rhs_imp, &
   &         assemble_tor, finish_exp_tor, get_rot_rates

contains

   subroutine initialize_updateZ

      integer :: ll, n_bands

      if ( l_finite_diff ) then
         allocate( type_bandmat :: zMat(n_lo_loc) )
         allocate( type_bandmat :: z10Mat )

         if ( ktopv /= 1 .and. kbotv /= 1 .and. rscheme_oc%order <= 2  .and. &
         &    rscheme_oc%order_boundary <= 2 ) then ! Rigid at both boundaries
            n_bands = rscheme_oc%order+1
         else
            n_bands = max(2*rscheme_oc%order_boundary+1,rscheme_oc%order+1)
         end if

         do ll=1,n_lo_loc
            call zMat(ll)%initialize(n_bands,n_r_max,l_pivot=.true.)
         end do

         !-- Special care when Inner Core or Mantle is free to rotate
         if ( ktopv /= 1 .and. kbotv /= 1 .and. rscheme_oc%order <= 2  .and. &
         &    rscheme_oc%order_boundary <= 2 .and. (.not. l_rot_ic) .and.    &
         &    (.not. l_rot_ma) ) then ! Rigid at both boundaries
            n_bands = rscheme_oc%order+1
         else
            n_bands = max(2*rscheme_oc%order_boundary+1,rscheme_oc%order+1)
         end if

         call z10Mat%initialize(n_bands,n_r_max,l_pivot=.true.)

      else
         allocate( type_densemat :: zMat(n_lo_loc) )
         allocate( type_densemat :: z10Mat )

         call z10Mat%initialize(n_r_max,n_r_max,l_pivot=.true.)
         do ll=1,n_lo_loc
            call zMat(ll)%initialize(n_r_max,n_r_max,l_pivot=.true.)
         end do
      end if

#ifdef WITH_PRECOND_Z10
      allocate(z10Mat_fac(n_r_max))
      bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_REAL
#endif
#ifdef WITH_PRECOND_Z
      allocate(zMat_fac(n_r_max,n_lo_loc))
      bytes_allocated = bytes_allocated+n_r_max*n_lo_loc*SIZEOF_DEF_REAL
#endif
      allocate( lZmat(n_lo_loc) )
      bytes_allocated = bytes_allocated+n_lo_loc*SIZEOF_LOGICAL

      allocate( Dif(n_mlo_loc) )
      bytes_allocated=bytes_allocated+n_mlo_loc*SIZEOF_DEF_COMPLEX

#ifdef WITHOMP
      maxThreads=omp_get_max_threads()
#else
      maxThreads=1
#endif

      allocate(rhs1(n_r_max,2*maxval(map_mlo%n_mi(:)),1))
      bytes_allocated=bytes_allocated+n_r_max*maxval(map_mlo%n_mi(:))*&
      &               SIZEOF_DEF_COMPLEX

      AMstart=0.0_cp

   end subroutine initialize_updateZ
!-------------------------------------------------------------------------------
   subroutine finalize_updateZ

      integer :: ll

      do ll=1,n_lo_loc
         call zMat(ll)%finalize()
      end do
      call z10Mat%finalize()

#ifdef WITH_PRECOND_Z10
      deallocate( z10Mat_fac )
#endif
#ifdef WITH_PRECOND_Z
      deallocate( zMat_fac )
#endif
      deallocate( rhs1, Dif )

   end subroutine finalize_updateZ
!-------------------------------------------------------------------------------
   subroutine updateZ(time,timeNext,z,dz,dzdt,omega_ma,omega_ic,domega_ma_dt,  &
              &       domega_ic_dt,lorentz_torque_ma_dt, lorentz_torque_ic_dt, &
              &       tscheme,lRmsNext)
      !
      !  updates the toroidal potential z and its radial derivatives
      !  adds explicit part to time derivatives of z
      !

      !-- Input/output of scalar fields:
      class(type_tscheme), intent(in) :: tscheme
      complex(cp), intent(inout) :: z(n_mlo_loc,n_r_max)   ! Toroidal velocity potential z
      type(type_tarray),  intent(inout) :: dzdt
      type(type_tscalar), intent(inout) :: domega_ic_dt
      type(type_tscalar), intent(inout) :: domega_ma_dt
      type(type_tscalar), intent(inout) :: lorentz_torque_ic_dt
      type(type_tscalar), intent(inout) :: lorentz_torque_ma_dt

      !-- Input of other variables:
      real(cp),           intent(in) :: time       ! Current stage time
      real(cp),           intent(in) :: timeNext   ! Next time
      logical,            intent(in) :: lRmsNext   ! Logical for storing update if (l_RMS.and.l_logNext)

      !-- Output variables
      complex(cp), intent(out) :: dz(n_mlo_loc,n_r_max) ! Radial derivative of z
      real(cp),    intent(out) :: omega_ma              ! Calculated OC rotation
      real(cp),    intent(out) :: omega_ic              ! Calculated IC rotation

      !-- local variables:
      integer :: l, m               ! degree and order
      integer :: lj, mi, i          ! l, m and ml counter
      integer :: nR                 ! counts radial grid points
      integer :: n_r_out            ! counts cheb modes
      real(cp) :: rhs(n_r_max)   ! RHS of matrix multiplication
      real(cp) :: prec_fac
      real(cp) :: dom_ma, dom_ic, lo_ma, lo_ic

      if ( l_precession ) then
         prec_fac=sqrt(8.0_cp*pi*third)*po*oek*oek*sin(prec_angle)
      else
         prec_fac = 0.0_cp
      end if

      if ( .not. l_update_v ) return

      if ( l_rot_ma  .and. (.not. l_SRMA) ) then
         if ( ktopv == 1 ) then
            call tscheme%set_imex_rhs_scalar(lo_ma, lorentz_torque_ma_dt)
         else
            call tscheme%set_imex_rhs_scalar(dom_ma, domega_ma_dt)
         end if
      end if

      if ( l_rot_ic .and. (.not. l_SRIC) ) then
         if ( kbotv == 1 ) then
            call tscheme%set_imex_rhs_scalar(lo_ic, lorentz_torque_ic_dt)
         else
            call tscheme%set_imex_rhs_scalar(dom_ic, domega_ic_dt)
         end if
      end if

      !-- Now assemble the right hand side and store it in work_LMdist
      call tscheme%set_imex_rhs(work_LMdist, dzdt, 1, n_mlo_loc, n_r_max)

      call solve_counter%start_count()

      ! Loop over local l
      do lj=1, n_lo_loc
         l = map_mlo%lj2l(lj)

         if ( l == 0 ) cycle

         if ( .not. lZmat(lj) ) then
#ifdef WITH_PRECOND_Z
            call get_zMat(tscheme,l,hdif_V(l),zMat(lj),zMat_fac(:,lj))
#else
            call get_zMat(tscheme,l,hdif_V(l),zMat(lj))
#endif
            lZmat(lj)=.true.
         end if

         ! Loop over m corresponding to current l
         do mi=1,map_mlo%n_mi(lj)
            m = map_mlo%milj2m(mi,lj)
            i = map_mlo%milj2i(mi,lj)

            if ( l_z10mat .and. l==1 .and. m == 0 ) then
               !----- Special treatment of z10 component if ic or mantle
               !      are allowed to rotate about z-axis (l_z10mat=.true.) and
               !      we use no slip boundary condition (ktopv=2,kbotv=2):
               !      Lorentz torque is the explicit part of this time integration
               !      at the boundaries!
               !      Note: no angular momentum correction necessary for this case !
               if ( .not. lZ10mat ) then
#ifdef WITH_PRECOND_Z10
                  call get_z10Mat(tscheme,l,hdif_V(l),z10Mat,z10Mat_fac)
#else
                  call get_z10Mat(tscheme,l,hdif_V(l),z10Mat)
#endif
                  lZ10mat=.true.
               end if

               if ( l_SRMA ) then
                  tOmega_ma1=time+tShift_ma1
                  tOmega_ma2=time+tShift_ma2
                  omega_ma= omega_ma1*cos(omegaOsz_ma1*tOmega_ma1) + &
                  &         omega_ma2*cos(omegaOsz_ma2*tOmega_ma2)
                  rhs(1)=omega_ma
               else if ( ktopv == 2 .and. l_rot_ma ) then  ! time integration
                  rhs(1)=dom_ma
               else
                  rhs(1)=0.0_cp
               end if

               if ( l_SRIC ) then
                  tOmega_ic1=time+tShift_ic1
                  tOmega_ic2=time+tShift_ic2
                  omega_ic= omega_ic1*cos(omegaOsz_ic1*tOmega_ic1) + &
                  &         omega_ic2*cos(omegaOsz_ic2*tOmega_ic2)
                  rhs(n_r_max)=omega_ic
               else if ( kbotv == 2 .and. l_rot_ic ) then  ! time integration
                  rhs(n_r_max)=dom_ic
               else
                  rhs(n_r_max)=0.0_cp
               end if

               !----- This is the normal RHS for the other radial grid points:
               do nR=2,n_r_max-1
                  rhs(nR)=real(work_LMdist(i,nR))
               end do

#ifdef WITH_PRECOND_Z10
               rhs(:) = z10Mat_fac(:)*rhs(:)
#endif

               call z10Mat%solve(rhs)

            else

               rhs1(1,2*mi-1,1)      =0.0_cp
               rhs1(1,2*mi,1)        =0.0_cp
               rhs1(n_r_max,2*mi-1,1)=0.0_cp
               rhs1(n_r_max,2*mi,1)  =0.0_cp

               if (amp_RiIc /= 0.0_cp) then
                  if (l == (m_RiIc + RiSymmIc) .and. m == m_RiIc) then
                     rhs1(n_r_max,2*mi-1,1)=amp_RiIc*cos(omega_RiIc*time)
                     rhs1(n_r_max,2*mi,1)  =amp_RiIc*sin(omega_RiIc*time)
                  end if
               end if

               if (amp_RiMa /= 0.0_cp) then
                  if (l == (m_RiMa + RiSymmMa) .and. m == m_RiMa) then
                     rhs1(1,2*mi-1,1)=amp_RiMa*cos(omega_RiMa*time)
                     rhs1(1,2*mi,1)  =amp_RiMa*sin(omega_RiMa*time)
                  end if
               end if

               do nR=2,n_r_max-1
                  rhs1(nR,2*mi-1,1)=real(work_LMdist(i,nR))
                  rhs1(nR,2*mi,1)  =aimag(work_LMdist(i,nR))
                  if ( l_precession .and. l == 1 .and. m == 1 ) then
                     rhs1(nR,2*mi-1,1)=rhs1(nR,2*mi-1,1)+ tscheme%wimp_lin(1)*    &
                     &                 prec_fac*sin(oek*time)
                     rhs1(nR,2*mi,1)=rhs1(nR,2*mi,1)-tscheme%wimp_lin(1)*prec_fac*&
                     &                    cos(oek*time)
                  end if
               end do

#ifdef WITH_PRECOND_Z
               rhs1(:,2*mi-1,1)=zMat_fac(:,lj)*rhs1(:,2*mi-1,1)
               rhs1(:,2*mi,1)  =zMat_fac(:,lj)*rhs1(:,2*mi,1)
#endif
                  !PERFOFF
            end if
         end do

         call zMat(lj)%solve(rhs1(:,1:2*map_mlo%n_mi(lj),1),2*map_mlo%n_mi(lj))

         ! Loop over m corresponding to current l (again)
         do mi=1,map_mlo%n_mi(lj)
            m = map_mlo%milj2m(mi,lj)
            i = map_mlo%milj2i(mi,lj)

            if ( l_z10mat .and. l==1 .and. m==0 ) then
               do n_r_out=1,rscheme_oc%n_max
                  z(i,n_r_out)=cmplx(rhs(n_r_out),0.0_cp,kind=cp)
               end do
            else
               if ( m > 0 ) then
                  do n_r_out=1,rscheme_oc%n_max
                     z(i,n_r_out)=cmplx(rhs1(n_r_out,2*mi-1,1), &
                     &                  rhs1(n_r_out,2*mi,1),kind=cp)
                  end do
               else
                  do n_r_out=1,rscheme_oc%n_max
                     z(i,n_r_out)=cmplx(rhs1(n_r_out,2*mi-1,1),0.0_cp,cp)
                  end do
               end if
            end if
         end do
      end do       ! end of loop over lm blocks
      call solve_counter%stop_count(l_increment=.false.)

      !-- set cheb modes > rscheme_oc%n_max to zero (dealiazing)
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         do i=1,n_mlo_loc
            z(i,n_r_out)=zero
         end do
      end do

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dzdt, 1, n_mlo_loc, n_r_max)
      call tscheme%rotate_imex_scalar(domega_ma_dt)
      call tscheme%rotate_imex_scalar(domega_ic_dt)
      call tscheme%rotate_imex_scalar(lorentz_torque_ma_dt)
      call tscheme%rotate_imex_scalar(lorentz_torque_ic_dt)

      !-- Calculation of the implicit part
      if (  tscheme%istage == tscheme%nstages ) then
         call get_tor_rhs_imp(timeNext, z, dz, dzdt, domega_ma_dt, domega_ic_dt, &
              &               omega_ic, omega_ma, omega_ic1, omega_ma1, tscheme, &
              &               1, tscheme%l_imp_calc_rhs(1), lRmsNext,            &
              &               l_in_cheb_space=.true.)
         call update_rot_rates(z, lo_ma, lo_ic, lorentz_torque_ma_dt,    &
              &                lorentz_torque_ic_dt, omega_ma,           &
              &                omega_ma1, omega_ic, omega_ic1, 1)
      else
         call get_tor_rhs_imp(time, z, dz, dzdt, domega_ma_dt, domega_ic_dt,     &
              &               omega_ic, omega_ma, omega_ic1, omega_ma1, tscheme, &
              &               tscheme%istage+1, tscheme%l_imp_calc_rhs(          &
              &               tscheme%istage+1), lRmsNext, l_in_cheb_space=.true.)
         call update_rot_rates(z, lo_ma, lo_ic, lorentz_torque_ma_dt,    &
              &                lorentz_torque_ic_dt, omega_ma,           &
              &                omega_ma1, omega_ic, omega_ic1, tscheme%istage+1)
      end if

   end subroutine updateZ
!------------------------------------------------------------------------------
   subroutine get_tor_rhs_imp(time, z, dz, dzdt, domega_ma_dt, domega_ic_dt, &
              &               omega_ic, omega_ma, omega_ic1, omega_ma1,      &
              &               tscheme, istage, l_calc_lin, lRmsNext,         &
              &               l_in_cheb_space)

      !-- Input variables
      real(cp),            intent(in) :: time
      integer,             intent(in) :: istage
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lRmsNext
      logical,             intent(in) :: l_calc_lin
      logical, optional,   intent(in) :: l_in_cheb_space

      !-- Output variable
      type(type_tarray),  intent(inout) :: dzdt
      type(type_tscalar), intent(inout) :: domega_ic_dt
      type(type_tscalar), intent(inout) :: domega_ma_dt
      real(cp),    intent(inout) :: omega_ic
      real(cp),    intent(inout) :: omega_ma
      real(cp),    intent(inout) :: omega_ic1
      real(cp),    intent(inout) :: omega_ma1
      complex(cp), intent(inout) :: z(n_mlo_loc,n_r_max)
      complex(cp), intent(out) :: dz(n_mlo_loc,n_r_max)

      !-- Local variables
      real(cp) :: angular_moment(3)   ! total angular momentum
      real(cp) :: angular_moment_oc(3)! x,y,z component of outer core angular mom.
      real(cp) :: angular_moment_ic(3)! x,y,z component of inner core angular mom.
      real(cp) :: angular_moment_ma(3)! x,y,z component of mantle angular mom.
      complex(cp) :: z10(n_r_max), z11(n_r_max)
      complex(cp) :: corr_l1m0, corr_l1m1
      real(cp) :: r_E_2, nomi, dL, prec_fac
      logical :: l_in_cheb
      integer :: n_r, lm, start_lm, stop_lm, n_r_bot, n_r_top, i
      integer :: l1, m1, l1m0, l1m1
      real(cp) :: ddzASL_loc(l_max+1,n_r_max)

      if ( l_precession ) then
         prec_fac=sqrt(8.0_cp*pi*third)*po*oek*oek*sin(prec_angle)
      else
         prec_fac = 0.0_cp
      end if

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
      call get_ddr( z, dz, work_LMdist, n_mlo_loc, start_lm, &
           &       stop_lm, n_r_max, rscheme_oc, l_dct_in=.not. l_in_cheb)
      if ( l_in_cheb ) call rscheme_oc%costf1(z, n_mlo_loc, start_lm, stop_lm)
      !$omp barrier
      !$omp single
      call dct_counter%stop_count(l_increment=.false.)
      !$omp end single

      l1m0=map_mlo%ml2i(0,1)
      l1m1=map_mlo%ml2i(1,1)

      !--- We correct so that the angular moment about axis in the equatorial plane
      !    vanish and the angular moment about the (planetary) rotation axis
      !    is kept constant.
      !$omp single
      if ( l_correct_AMz .and.  l1m0 > 0 ) then

         z10(:)=z(l1m0,:)
         call get_angular_moment(z10,z11,omega_ic,omega_ma,          &
              &                  angular_moment_oc,angular_moment_ic,&
              &                  angular_moment_ma)
         do i=1,3
            angular_moment(i)=angular_moment_oc(i) + angular_moment_ic(i) + &
            &                 angular_moment_ma(i)
         end do
         if ( ( ktopv == 2 .and. l_rot_ma ) .and. ( kbotv == 2 .and. l_rot_ic ) ) then
            nomi=c_moi_ma*c_z10_omega_ma*r_cmb*r_cmb + &
            &    c_moi_ic*c_z10_omega_ic*r_icb*r_icb + &
            &    c_moi_oc*y10_norm
         else if ( ktopv == 2 .and. l_rot_ma ) then
            nomi=c_moi_ma*c_z10_omega_ma*r_cmb*r_cmb+c_moi_oc*y10_norm
         else if ( kbotv == 2 .and. l_rot_ic ) then
            nomi=c_moi_ic*c_z10_omega_ic*r_icb*r_icb+c_moi_oc*y10_norm
         else
            nomi=c_moi_oc*y10_norm
         end if
         corr_l1m0=cmplx(angular_moment(3)-AMstart,-1.0_cp,kind=cp)/nomi

         !-------- Correct z(2,n_r) and z(l_max+2,n_r) plus the respective
         !         derivatives:
         do n_r=1,n_r_max
            r_E_2=r(n_r)*r(n_r)
            z(l1m0,n_r)  =z(l1m0,n_r)  - rho0(n_r)*r_E_2*corr_l1m0
            dz(l1m0,n_r) =dz(l1m0,n_r) - rho0(n_r)*(          &
            &            two*r(n_r)+r_E_2*beta(n_r))*corr_l1m0
            work_LMdist(l1m0,n_r)=work_LMdist(l1m0,n_r)-rho0(n_r)*( &
            &                 two+four*beta(n_r)*r(n_r) +           &
            &                  dbeta(n_r)*r_E_2 +                   &
            &                beta(n_r)*beta(n_r)*r_E_2 )*corr_l1m0
         end do

         if ( ktopv == 2 .and. l_rot_ma ) &
         &    omega_ma=c_z10_omega_ma*real(z(l1m0,n_r_cmb))
         if ( kbotv == 2 .and. l_rot_ic ) &
         &    omega_ic=c_z10_omega_ic*real(z(l1m0,n_r_icb))
         omega_ic1=omega_ic
         omega_ma1=omega_ma

      end if ! l=1,m=0 contained in lm-block ?

      if ( l_correct_AMe .and.  l1m1 > 0 ) then
         z11(:)=z(l1m1,:)
         call get_angular_moment(z10,z11,omega_ic,omega_ma,          &
              &                  angular_moment_oc,angular_moment_ic,&
              &                  angular_moment_ma)
         do i=1,3
            angular_moment(i)=angular_moment_oc(i) + angular_moment_ic(i) + &
            &                 angular_moment_ma(i)
         end do
         corr_l1m1=cmplx(angular_moment(1),-angular_moment(2),kind=cp) / &
         &         (two*y11_norm*c_moi_oc)

         !-------- Correct z(2,n_r) and z(l_max+2,n_r) plus the respective
         !         derivatives:
         do n_r=1,n_r_max
            r_E_2=r(n_r)*r(n_r)
            z(l1m1,n_r)  =z(l1m1,n_r)  -  rho0(n_r)*r_E_2*corr_l1m1
            dz(l1m1,n_r) =dz(l1m1,n_r) -  rho0(n_r)*(            &
            &            two*r(n_r)+r_E_2*beta(n_r))*corr_l1m1
            work_LMdist(l1m1,n_r)=work_LMdist(l1m1,n_r)-rho0(n_r)*(  &
            &              two+four*beta(n_r)*r(n_r) +               &
            &                          dbeta(n_r)*r_E_2 +            &
            &                beta(n_r)*beta(n_r)*r_E_2 )*corr_l1m1
         end do
      end if ! l=1,m=1 contained in lm-block ?
      !$omp end single


      if ( istage == 1 ) then
         !$omp do private(n_r,lm,l1,dL)
         do n_r=1,n_r_max
            do lm=1,n_mlo_loc
               l1 = map_mlo%i2l(lm)
               dL = real(l1*(l1+1),cp)
               dzdt%old(lm,n_r,istage)=dL*or2(n_r)*z(lm,n_r)
            end do
         end do
         !$omp end do
      end if

      if ( l_calc_lin .or. (tscheme%istage==tscheme%nstages .and. lRmsNext)) then

         if ( lRmsNext ) then
            n_r_top=n_r_cmb
            n_r_bot=n_r_icb
         else
            n_r_top=n_r_cmb+1
            n_r_bot=n_r_icb-1
         end if

         !$omp do private(n_r,lm,Dif,l1,dL)
         do n_r=n_r_top,n_r_bot
            do lm=1,n_mlo_loc
               l1 = map_mlo%i2l(lm)
               m1 = map_mlo%i2m(lm)
               if ( l1 == 0 ) cycle
               dL = real(l1*(l1+1),cp)
               Dif(lm)=hdif_V(l1)*dL*or2(n_r)*visc(n_r)* ( work_LMdist(lm,n_r) + &
               &         (dLvisc(n_r)-beta(n_r))    *              dz(lm,n_r) -  &
               &         ( dLvisc(n_r)*beta(n_r)+two*dLvisc(n_r)*or1(n_r)        &
               &          + dL*or2(n_r)+dbeta(n_r)+two*beta(n_r)*or1(n_r) )*     &
               &                                                    z(lm,n_r) )

               dzdt%impl(lm,n_r,istage)=Dif(lm)
               if ( l_precession .and. l1==1 .and. m1==1 ) then
                  dzdt%impl(lm,n_r,istage)=dzdt%impl(lm,n_r,istage)+prec_fac*cmplx( &
                  &                        sin(oek*time),-cos(oek*time),kind=cp)
               end if
            end do
            if ( lRmsNext .and. tscheme%istage==tscheme%nstages ) then
               call hInt2Tor(Dif,1,n_mlo_loc,n_r,DifTor2hInt(:,n_r))
            end if
         end do
         !$omp end do

      end if

      !$omp end parallel

      if ( (l1m0>0) .and. l_z10mat ) then
         !----- NOTE opposite sign of viscous torque on ICB and CMB:
         if ( .not. l_SRMA .and. ktopv == 2 .and. l_rot_ma ) then
            domega_ma_dt%impl(istage)=visc(n_r_cmb)*( (two*or1(n_r_cmb)+   &
            &                         beta(n_r_cmb))*real(z(l1m0,n_r_cmb))-&
            &                                       real(dz(l1m0,n_r_cmb)) )
            if ( istage == 1 ) domega_ma_dt%old(istage)=c_dt_z10_ma* &
            &                                           real(z(l1m0,n_r_cmb))
         end if
         if ( .not. l_SRIC .and. kbotv == 2 .and. l_rot_ic ) then
            domega_ic_dt%impl(istage)=-visc(n_r_icb)* ( (two*or1(n_r_icb)+   &
            &                          beta(n_r_icb))*real(z(l1m0,n_r_icb))- &
            &                                        real(dz(l1m0,n_r_icb)) )
            if ( istage == 1 ) domega_ic_dt%old(istage)=c_dt_z10_ic* &
            &                                           real(z(l1m0,n_r_icb))
         end if
      end if

      !--- Note: from ddz=work_LMdist only the axisymmetric contributions are needed
      !    beyond this point for the TO calculation.
      !    Parallization note: Very likely, all axisymmetric modes m=0 are
      !    located on the first processor #0.
      if ( l_TO ) then
         !$omp parallel do default(shared) private(n_r,lm,l1,m1)
         do n_r=1,n_r_max
            ddzASL_loc(:,n_r)=0.0_cp
            do lm=1,n_mlo_loc
               l1=map_mlo%i2l(lm)
               m1=map_mlo%i2m(lm)
               if ( m1 == 0 ) ddzASL_loc(l1+1,n_r)=real(work_LMdist(lm,n_r))
            end do
         end do
         !$omp end parallel do

         do n_r=1,n_r_max
#ifdef WITH_MPI
            call MPI_Allreduce(ddzASL_loc(:,n_r), ddzASL(:,n_r), l_max+1, &
                 &             MPI_DEF_REAL, MPI_SUM, comm_r, ierr)
#else
            ddzASL(:,n_r)=ddzASL_loc(:,n_r)
#endif
         end do
      end if


   end subroutine get_tor_rhs_imp
!------------------------------------------------------------------------------
   subroutine assemble_tor(time, z, dz, dzdt, domega_ic_dt, domega_ma_dt,        &
              &            lorentz_torque_ic_dt, lorentz_torque_ma_dt, omega_ic, &
              &            omega_ma, omega_ic1, omega_ma1, lRmsNext, tscheme)
      !
      ! This subroutine is used to assemble the toroidal flow potential when an IMEX
      ! RK time scheme with an assembly stage is employed.
      !

      !-- Input variable
      class(type_tscheme), intent(in) :: tscheme
      real(cp),            intent(in) :: time
      logical,             intent(in) :: lRmsNext

      !-- Output variables
      complex(cp),        intent(inout) :: z(n_mlo_loc,n_r_max)
      complex(cp),        intent(out) :: dz(n_mlo_loc,n_r_max)
      type(type_tarray),  intent(inout) :: dzdt
      type(type_tscalar), intent(inout) :: domega_ic_dt
      type(type_tscalar), intent(inout) :: domega_ma_dt
      type(type_tscalar), intent(inout) :: lorentz_torque_ic_dt
      type(type_tscalar), intent(inout) :: lorentz_torque_ma_dt
      real(cp),           intent(inout) :: omega_ic
      real(cp),           intent(inout) :: omega_ma
      real(cp),           intent(inout) :: omega_ic1
      real(cp),           intent(inout) :: omega_ma1

      !-- Local variables
      complex(cp) :: top_val(n_mlo_loc), bot_val(n_mlo_loc)
      integer :: n_r, lm, l, m, l1m0
      real(cp) :: fac_top, fac_bot, dom_ma, dom_ic, lo_ma, lo_ic
      real(cp) :: dL

      if ( amp_RiIc /= 0.0_cp .or. amp_RiMa /= 0.0_cp ) then
         call abortRun('Not implemented yet in assembly stage of z')
      end if

      if ( l_rot_ma  .and. (.not. l_SRMA) ) then
         if ( ktopv == 1 ) then
            call tscheme%assemble_imex_scalar(lo_ma, lorentz_torque_ma_dt)
         else
            call tscheme%assemble_imex_scalar(dom_ma, domega_ma_dt)
         end if
      end if

      if ( l_rot_ic .and. (.not. l_SRIC) ) then
         if ( kbotv == 1 ) then
            call tscheme%assemble_imex_scalar(lo_ic, lorentz_torque_ic_dt)
         else
            call tscheme%assemble_imex_scalar(dom_ic, domega_ic_dt)
         end if
      end if

      !-- Store the assembled quantity in work_LMdist
      call tscheme%assemble_imex(work_LMdist, dzdt, 1, n_mlo_loc, n_r_max)

      !-- Now get the toroidal potential from the assembly
      !$omp parallel default(shared)
      !$omp do private(n_r,lm,l,m,dL)
      do n_r=2,n_r_max-1
         do lm=1,n_mlo_loc
            l = map_mlo%i2l(lm)
            m = map_mlo%i2m(lm)
            if ( l == 0 ) cycle
            dL = real(l*(l+1),cp)
            if ( m == 0 ) then
               z(lm,n_r) = r(n_r)*r(n_r)/dL*cmplx(real(work_LMdist(lm,n_r)),0.0_cp,cp)
            else
               z(lm,n_r) = r(n_r)*r(n_r)/dL*work_LMdist(lm,n_r)
            end if
         end do
      end do
      !$omp end do

      l1m0 = map_mlo%ml2i(0,1)
      !$omp do private(lm)
      do lm=1,n_mlo_loc
         if ( lm==l1m0 ) then
            if ( l_SRMA ) then
               tOmega_ma1=time+tShift_ma1
               tOmega_ma2=time+tShift_ma2
               omega_ma= omega_ma1*cos(omegaOsz_ma1*tOmega_ma1) + &
               &         omega_ma2*cos(omegaOsz_ma2*tOmega_ma2)
               top_val(lm)=cmplx(omega_ma/c_z10_omega_ma,0.0_cp,kind=cp)
            else if ( ktopv == 2 .and. l_rot_ma ) then  ! time integration
               top_val(lm)=cmplx(dom_ma/c_dt_z10_ma,0.0_cp,kind=cp)
            else
               top_val(lm)=zero
            end if

            if ( l_SRIC ) then
               tOmega_ic1=time+tShift_ic1
               tOmega_ic2=time+tShift_ic2
               omega_ic= omega_ic1*cos(omegaOsz_ic1*tOmega_ic1) + &
               &         omega_ic2*cos(omegaOsz_ic2*tOmega_ic2)
               bot_val(lm)=cmplx(omega_ic/c_z10_omega_ic,0.0_cp,kind=cp)
            else if ( kbotv == 2 .and. l_rot_ic ) then  ! time integration
               bot_val(lm)=cmplx(dom_ic/c_dt_z10_ic,0.0_cp,kind=cp)
            else
               bot_val(lm)=zero
            end if
         else
            top_val(lm)=zero
            bot_val(lm)=zero
         end if
      end do
      !$omp end do

      !-- Boundary conditions
      if ( l_full_sphere ) then
         if ( ktopv /= 1 ) then ! Rigid outer
            !$omp do private(lm,l)
            do lm=1,n_mlo_loc
               l = map_mlo%i2l(lm)
               if ( l == 0 ) cycle
               z(lm,1)      =top_val(lm)
               z(lm,n_r_max)=bot_val(lm)
            end do
            !$omp end do
         else
            fac_top=-two*or1(1)-beta(1)
            !$omp do private(lm,l)
            do lm=1,n_mlo_loc
               l = map_mlo%i2l(lm)
               if ( l == 0 ) cycle
               call rscheme_oc%robin_bc(one, fac_top, zero, 0.0_cp, one, bot_val(lm), z(lm,:))
            end do
            !$omp end do
         end if
      else
         if ( ktopv /= 1 .and. kbotv /= 1 ) then ! Rigid BCs
            !$omp do private(lm,l)
            do lm=1,n_mlo_loc
               l = map_mlo%i2l(lm)
               if ( l == 0 ) cycle
               z(lm,1)      =top_val(lm)
               z(lm,n_r_max)=bot_val(lm)
            end do
            !$omp end do
         else if ( ktopv == 1 .and. kbotv /= 1 ) then ! Stress-free top and rigid bottom
            fac_top=-two*or1(1)-beta(1)
            !$omp do private(lm,l)
            do lm=1,n_mlo_loc
               l = map_mlo%i2l(lm)
               if ( l == 0 ) cycle
               call rscheme_oc%robin_bc(one, fac_top, zero, 0.0_cp, one, bot_val(lm), z(lm,:))
            end do
            !$omp end do
         else if ( kbotv == 1 .and. ktopv /= 1 ) then ! Stress-free bot and rigid top
            fac_bot=-two*or1(n_r_max)-beta(n_r_max)
            !$omp do private(lm,l)
            do lm=1,n_mlo_loc
               l = map_mlo%i2l(lm)
               if ( l == 0 ) cycle
               call rscheme_oc%robin_bc(0.0_cp, one, top_val(lm), one, fac_bot, zero, z(lm,:))
            end do
            !$omp end do
         else if ( ktopv == 1 .and. kbotv == 1 ) then ! Stress-free at both boundaries
            fac_bot=-two*or1(n_r_max)-beta(n_r_max)
            fac_top=-two*or1(1)-beta(1)
            !$omp do private(lm,l)
            do lm=1,n_mlo_loc
               l = map_mlo%i2l(lm)
               if ( l == 0 ) cycle
               call rscheme_oc%robin_bc(one, fac_top, zero, one, fac_bot, zero, z(lm,:))
            end do
            !$omp end do
         end if
      end if
      !$omp end parallel

      call get_tor_rhs_imp(time, z, dz, dzdt, domega_ma_dt, domega_ic_dt, &
           &               omega_ic, omega_ma, omega_ic1, omega_ma1,      &
           &               tscheme, 1, tscheme%l_imp_calc_rhs(1),         &
           &               lRmsNext, .false.)
      call update_rot_rates(z, lo_ma, lo_ic, lorentz_torque_ma_dt,    &
           &                lorentz_torque_ic_dt, omega_ma,           &
           &                omega_ma1, omega_ic, omega_ic1, 1)

   end subroutine assemble_tor
!------------------------------------------------------------------------------
   subroutine update_rot_rates(z, lo_ma, lo_ic, lorentz_torque_ma_dt,       &
              &                lorentz_torque_ic_dt, omega_ma, omega_ma1,   &
              &                omega_ic, omega_ic1, istage)

      !-- Input variables
      complex(cp), intent(in) :: z(n_mlo_loc,n_r_max)
      real(cp),    intent(in) :: lo_ma, lo_ic   ! RHS when stress-free BCs are used
      integer,     intent(in) :: istage

      !-- Output variables
      type(type_tscalar), intent(inout) :: lorentz_torque_ma_dt
      type(type_tscalar), intent(inout) :: lorentz_torque_ic_dt
      real(cp),           intent(out) :: omega_ma, omega_ma1
      real(cp),           intent(out) :: omega_ic, omega_ic1

      !-- Local variables
      integer :: l1m0

      l1m0=map_mlo%ml2i(0,1)
      !--- Update of inner core and mantle rotation:
      if ( ( l1m0 > 0 ) .and. l_z10mat ) then
         if ( l_rot_ma .and. .not. l_SRMA ) then
            if ( ktopv == 1 ) then  ! free slip, explicit time stepping of omega !
               call get_rot_rates(omega_ma, lorentz_torque_ma_dt%old(istage))
               omega_ma=lo_ma
            else if ( ktopv == 2 ) then ! no slip, omega given by z10
               omega_ma=c_z10_omega_ma*real(z(l1m0,n_r_cmb))
            end if
            omega_ma1=omega_ma
         end if
         if ( l_rot_ic .and. .not. l_SRIC ) then
            if ( kbotv == 1 ) then  ! free slip, explicit time stepping of omega !
               call get_rot_rates(omega_ic, lorentz_torque_ic_dt%old(istage))
               omega_ic=lo_ic
            else if ( kbotv == 2 ) then ! no slip, omega given by z10
               omega_ic=c_z10_omega_ic*real(z(l1m0,n_r_icb))
            end if
            omega_ic1=omega_ic
         end if
      end if  ! l=1,m=0 contained in block ?

   end subroutine update_rot_rates
!------------------------------------------------------------------------------
   subroutine get_rot_rates(omega, domega_old)

      !-- Input variables
      real(cp), intent(in) :: omega ! Rotation rate

      !-- Output variable
      real(cp), intent(out) :: domega_old ! Old value of the rotation rate

      domega_old = omega

   end subroutine get_rot_rates
!------------------------------------------------------------------------------
   subroutine finish_exp_tor(lorentz_torque_ma, lorentz_torque_ic, domega_ma_dt_exp, &
              &              domega_ic_dt_exp, lorentz_ma_exp, lorentz_ic_exp)

      !-- Input variables
      real(cp), intent(in) :: lorentz_torque_ma  ! Lorentz torque (for OC rotation)
      real(cp), intent(in) :: lorentz_torque_ic  ! Lorentz torque (for IC rotation)

      !-- Output variables
      real(cp), intent(out) :: domega_ic_dt_exp
      real(cp), intent(out) :: domega_ma_dt_exp
      real(cp), intent(out) :: lorentz_ic_exp
      real(cp), intent(out) :: lorentz_ma_exp

      domega_ic_dt_exp=0.0_cp
      domega_ma_dt_exp=0.0_cp
      lorentz_ic_exp  =0.0_cp
      lorentz_ma_exp  =0.0_cp

      if ( l_rot_ma  .and. (.not. l_SRMA) ) then
         if ( ktopv == 1 ) then !-- Stress-free
            lorentz_ma_exp=lorentz_torque_ma/c_moi_ma
         else
            domega_ma_dt_exp=c_lorentz_ma*lorentz_torque_ma
         end if
      end if

      if ( l_rot_ic .and. (.not. l_SRIC) ) then
         if ( kbotv == 1 ) then !-- Stress-free
            lorentz_ic_exp=lorentz_torque_ic/c_moi_ic
         else
            domega_ic_dt_exp=c_lorentz_ic*lorentz_torque_ic
         end if
      end if

   end subroutine finish_exp_tor
!------------------------------------------------------------------------------
#ifdef WITH_PRECOND_Z10
   subroutine get_z10Mat(tscheme,l,hdif,zMat,zMat_fac)
#else
   subroutine get_z10Mat(tscheme,l,hdif,zMat)
#endif
      !
      !  Purpose of this subroutine is to construct and LU-decompose the
      !  inversion matrix z10mat for the implicit time step of the
      !  toroidal velocity potential z of degree l=1 and order m=0.
      !  This differs from the the normal zmat only if either the ICB or
      !  CMB have no-slip boundary condition and inner core or mantle are
      !  chosen to rotate freely (either kbotv=1 and/or ktopv=1).
      !

      class(type_tscheme), intent(in) :: tscheme ! Time step internal
      real(cp),            intent(in) :: hdif    ! Value of hyperdiffusivity in zMat terms
      integer,             intent(in) :: l       ! Variable to loop over l's

      !-- Output: z10Mat and pivot_z10
      class(type_realmat), intent(inout) :: zMat
#ifdef WITH_PRECOND_Z10
      real(cp), intent(out) :: zMat_fac(n_r_max)     ! Inverse of max(zMat) for inversion
#endif

      !-- local variables:
      integer :: nR,nR_out,info
      real(cp) :: dLh
      real(cp) :: dat(n_r_max,n_r_max)

      dLh=real(l*(l+1),kind=cp)

      !-- Boundary conditions:
      !----- CMB condition:
      !        Note opposite sign of viscous torques (-dz+(2/r+beta) z )*visc
      !        for CMB and ICB!

      if ( ktopv == 1 ) then  ! free slip
         dat(1,:)=                   rscheme_oc%rnorm *     &
         &    ( (two*or1(1)+beta(1))*rscheme_oc%rMat(1,:) - &
         &                          rscheme_oc%drMat(1,:) )
      else if ( ktopv == 2 ) then ! no slip
         if ( l_SRMA ) then
            dat(1,:)= rscheme_oc%rnorm * c_z10_omega_ma* &
            &               rscheme_oc%rMat(1,:)
         else if ( l_rot_ma ) then
            dat(1,:)= rscheme_oc%rnorm *               (        &
            &                c_dt_z10_ma*rscheme_oc%rMat(1,:) - &
            &       tscheme%wimp_lin(1)*visc(1)*(               &
            &       (two*or1(1)+beta(1))*rscheme_oc%rMat(1,:) - &
            &                           rscheme_oc%drMat(1,:) ) )
         else
            dat(1,:)= rscheme_oc%rnorm*rscheme_oc%rMat(1,:)
         end if
      end if

      !----- ICB condition:
      if ( l_full_sphere ) then
         dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
      else
         if ( kbotv == 1 ) then  ! free slip
            dat(n_r_max,:)=rscheme_oc%rnorm *                  &
            &            ( (two*or1(n_r_max)+beta(n_r_max))*   &
            &                   rscheme_oc%rMat(n_r_max,:) -   &
            &                  rscheme_oc%drMat(n_r_max,:) )
         else if ( kbotv == 2 ) then ! no slip
            if ( l_SRIC ) then
               dat(n_r_max,:)= rscheme_oc%rnorm * &
               &             c_z10_omega_ic*rscheme_oc%rMat(n_r_max,:)
            else if ( l_rot_ic ) then     !  time integration of z10
               dat(n_r_max,:)= rscheme_oc%rnorm *             (          &
               &                c_dt_z10_ic*rscheme_oc%rMat(n_r_max,:) + &
               &       tscheme%wimp_lin(1)*visc(n_r_max)*(               &
               &           (two*or1(n_r_max)+beta(n_r_max))*             &
               &                            rscheme_oc%rMat(n_r_max,:) - &
               &                           rscheme_oc%drMat(n_r_max,:) ) )
            else
               dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
            end if
         end if
      end if

      !-- Fill up with zeros:
      do nR_out=rscheme_oc%n_max+1,n_r_max
         dat(1,nR_out)      =0.0_cp
         dat(n_r_max,nR_out)=0.0_cp
      end do

      !----- Other points: (same as zMat)
      do nR_out=1,n_r_max
         do nR=2,n_r_max-1
            dat(nR,nR_out)=rscheme_oc%rnorm * (                         &
            &             dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out) -      &
            &         tscheme%wimp_lin(1)*hdif*dLh*visc(nR)*or2(nR) * ( &
            &                            rscheme_oc%d2rMat(nR,nR_out) + &
            &    (dLvisc(nR)- beta(nR))*  rscheme_oc%drMat(nR,nR_out) - &
            &    ( dLvisc(nR)*beta(nR)+two*dLvisc(nR)*or1(nR)  +        &
            &      dLh*or2(nR)+dbeta(nR)+two*beta(nR)*or1(nR) )*        &
                                           rscheme_oc%rMat(nR,nR_out) ) )
         end do
      end do

      !-- Normalisation
      do nR=1,n_r_max
         dat(nR,1)      =rscheme_oc%boundary_fac*dat(nR,1)
         dat(nR,n_r_max)=rscheme_oc%boundary_fac*dat(nR,n_r_max)
      end do

#ifdef WITH_PRECOND_Z10
      ! compute the linesum of each line
      do nR=1,n_r_max
         zMat_fac(nR)=one/maxval(abs(dat(nR,:)))
         dat(nR,:) = dat(nR,:)*zMat_fac(nR)
      end do
#endif

      !-- Array copy
      call zMat%set_data(dat)

      !-- LU-decomposition of z10mat:
      call zMat%prepare(info)

      if ( info /= 0 ) call abortRun('Error from get_z10Mat: singular matrix!')

   end subroutine get_z10Mat
!-------------------------------------------------------------------------------
#ifdef WITH_PRECOND_Z
   subroutine get_zMat(tscheme,l,hdif,zMat,zMat_fac)
#else
   subroutine get_zMat(tscheme,l,hdif,zMat)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matrices
      !  zmat(i,j) for the NS equation.
      !

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme  ! time step
      integer,             intent(in) :: l        ! Variable to loop over degrees
      real(cp),            intent(in) :: hdif     ! Hyperdiffusivity

      !-- Output variables:
      class(type_realmat), intent(inout) :: zMat
#ifdef WITH_PRECOND_Z
      real(cp), intent(out) :: zMat_fac(n_r_max)     !  Inverse of max(zMat) for the inversion
#endif

      !-- local variables:
      integer :: nR,nR_out
      integer :: info
      real(cp) :: dLh
      real(cp) :: dat(n_r_max,n_r_max)
      character(len=80) :: message
      character(len=14) :: str, str_1

      dLh=real(l*(l+1),kind=cp)

      !----- Boundary conditions, see above:
      if ( ktopv == 1 ) then  ! free slip !
         dat(1,:)=rscheme_oc%rnorm *             (      & 
         &                     rscheme_oc%drMat(1,:) -  &
         & (two*or1(1)+beta(1))*rscheme_oc%rMat(1,:) )
      else                    ! no slip, note exception for l=1,m=0
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%rMat(1,:)
      end if

      if ( l_full_sphere ) then
         dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
      else
         if ( kbotv == 1 ) then  ! free slip !
            dat(n_r_max,:)=rscheme_oc%rnorm *            (   &
            &                  rscheme_oc%drMat(n_r_max,:) - &
            &    (two*or1(n_r_max)+beta(n_r_max))*           &
            &                   rscheme_oc%rMat(n_r_max,:) )
         else                    ! no slip, note exception for l=1,m=0
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
         end if
      end if

      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme_oc%n_max+1,n_r_max
            dat(1,nR_out)      =0.0_cp
            dat(n_r_max,nR_out)=0.0_cp
         end do
      end if

      !----- Bulk points:
      do nR_out=1,n_r_max
         do nR=2,n_r_max-1
            dat(nR,nR_out)=rscheme_oc%rnorm * (                          &
            &               dLh*or2(nR)* rscheme_oc%rMat(nR,nR_out)      &
            &   -tscheme%wimp_lin(1)*hdif*dLh*visc(nR)*or2(nR) * (       &
            &                               rscheme_oc%d2rMat(nR,nR_out) &
            &   + (dLvisc(nR)- beta(nR)) *   rscheme_oc%drMat(nR,nR_out) &
            &      - ( dLvisc(nR)*beta(nR)+two*dLvisc(nR)*or1(nR)        &
            &          +dLh*or2(nR)+dbeta(nR)+two*beta(nR)*or1(nR)       &
            &                             ) * rscheme_oc%rMat(nR,nR_out) ) )
         end do
      end do

      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         dat(nR,1)      =rscheme_oc%boundary_fac*dat(nR,1)
         dat(nR,n_r_max)=rscheme_oc%boundary_fac*dat(nR,n_r_max)
      end do

#ifdef WITH_PRECOND_Z
      ! compute the linesum of each line
      do nR=1,n_r_max
         zMat_fac(nR)=one/maxval(abs(dat(nR,:)))
         dat(nR,:)   =dat(nR,:)*zMat_fac(nR)
      end do
#endif

#ifdef MATRIX_CHECK
      block

      integer :: i,j
      real(cp) :: rcond
      integer ::ipiv(n_r_max),iwork(n_r_max)
      real(cp) :: work(4*n_r_max),anorm,linesum
      real(cp) :: temp_Mat(n_r_max,n_r_max)
      integer, save :: counter=0
      integer :: filehandle
      character(len=100) :: filename

      ! copy the zMat to a temporary variable for modification
      write(filename,"(A,I3.3,A,I3.3,A)") "zMat_",l,"_",counter,".dat"
      open(newunit=filehandle,file=trim(filename))
      counter= counter+1

      do i=1,n_r_max
         do j=1,n_r_max
            write(filehandle,"(2ES20.12,1X)",advance="no") dat(i,j)
         end do
         write(filehandle,"(A)") ""
      end do
      close(filehandle)
      temp_Mat=dat
      anorm = 0.0_cp
      do i=1,n_r_max
         linesum = 0.0_cp
         do j=1,n_r_max
            linesum = linesum + abs(temp_Mat(i,j))
         end do
         if (linesum  >  anorm) anorm=linesum
      end do
      !write(*,"(A,ES20.12)") "anorm = ",anorm
      ! LU factorization
      call dgetrf(n_r_max,n_r_max,temp_Mat,n_r_max,ipiv,info)
      ! estimate the condition number
      call dgecon('I',n_r_max,temp_Mat,n_r_max,anorm,rcond,work,iwork,info)
      write(*,"(A,I3,A,ES11.3)") "inverse condition number of zMat for l=",l," is ",rcond

      end block
#endif

      !-- Array copy
      call zMat%set_data(dat)

      !-- LU decomposition:
      call zMat%prepare(info)

      if ( info /= 0 ) then
         write(str, *) l
         write(str_1, *) info
         message='Singular matrix zmat for l='//trim(adjustl(str))//&
         &       ', info = '//trim(adjustl(str_1))
         call abortRun(message)
      end if

   end subroutine get_zMat
!-------------------------------------------------------------------------------
end module updateZ_mod
