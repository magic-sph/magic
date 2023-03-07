module nonlinear_bcs

   use iso_fortran_env, only: output_unit
   use precision_mod
   use truncation, only: lm_max, n_phi_max, l_max, n_theta_max, nlat_padded
   use radial_data, only: n_r_cmb, n_r_icb, nRstart, nRstop
   use radial_functions, only: r_cmb, r_icb, rho0, orho1, or2
   use physical_parameters, only: sigma_ratio, conductance_ma, prmag, oek
   use horizontal_data, only: cosTheta, sinTheta_E2
   use constants, only: two
   use useful, only: abortRun
   use sht, only: spat_to_sphertor, sht_l_single

   implicit none

   private

   public :: get_br_v_bcs, get_br_v_bcs_batch, get_b_nl_bcs, v_rigid_boundary, v_rigid_boundary_batch

contains

   subroutine get_br_v_bcs(nR,br,vt,vp,omega,br_vt_lm,br_vp_lm)
      !
      !  Purpose of this subroutine is to calculate the nonlinear term
      !  of the magnetic boundary condition for a conducting mantle or
      !  inner core in space (r,lm).
      !  Calculation is performed for the theta block:
      !
      !  .. code-block:: fortran
      !
      !                  n_theta_min<=n_theta<=n_theta_min+n_theta_block-1
      !
      !  On input ``br``, ``vt`` and ``vp`` are given on all phi points and
      !  thetas in the specific block.
      !  On output the contribution of these grid points to all
      !  degree and orders is stored in ``br_vt_lm`` and ``br_vp_lm``.
      !  Output is ``[r/sin(theta)*Br*U]=[(0,br_vt_lm,br_vp_lm)]``
      !

      !-- input:
      real(cp), intent(in) :: br(:,:)      ! :math:`r^2 B_r`
      real(cp), intent(in) :: vt(:,:)      ! :math:`r \sin\theta u_\theta`
      real(cp), intent(in) :: vp(:,:)      ! :math:`r \sin\theta u_\phi`
      real(cp), intent(in) :: omega        ! rotation rate of mantle or IC
      integer,  intent(in) :: nR

      !-- Output variables:
      ! br*vt/(sin(theta)**2*r**2)
      complex(cp), intent(inout) :: br_vt_lm(lm_max)
      ! br*(vp/(sin(theta)**2*r**2)-omega_ma)
      complex(cp), intent(inout) :: br_vp_lm(lm_max)

      !-- Local variables:
      integer :: n_theta, n_phi
      real(cp) :: br_vt(nlat_padded,n_phi_max), br_vp(nlat_padded,n_phi_max)
      real(cp) :: fac          ! 1/( r**2 sin(theta)**2 )

      fac=or2(nR)*orho1(nR)
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do collapse(2)
#else
      !$omp parallel do default(shared) &
      !$omp& private(n_theta,n_phi)
#endif
      do n_phi=1,n_phi_max
         do n_theta=1,n_theta_max
            br_vt(n_theta,n_phi)= fac*br(n_theta,n_phi)*vt(n_theta,n_phi)
            br_vp(n_theta,n_phi)= br(n_theta,n_phi)*( fac*vp(n_theta,n_phi) - &
            &                                  omega*sinTheta_E2(n_theta) )
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#else
      !$omp end parallel do
#endif

      call spat_to_sphertor(sht_l_single, br_vt, br_vp, br_vt_lm, br_vp_lm, l_max)

   end subroutine get_br_v_bcs

   subroutine get_br_v_bcs_batch(nR,br,vt,vp,omega,br_vt_lm,br_vp_lm)
      !
      !  Purpose of this subroutine is to calculate the nonlinear term
      !  of the magnetic boundary condition for a conducting mantle or
      !  inner core in space (r,lm).
      !  Calculation is performed for the theta block:
      !
      !  .. code-block:: fortran
      !
      !                  n_theta_min<=n_theta<=n_theta_min+n_theta_block-1
      !
      !  On input ``br``, ``vt`` and ``vp`` are given on all phi points and
      !  thetas in the specific block.
      !  On output the contribution of these grid points to all
      !  degree and orders is stored in ``br_vt_lm`` and ``br_vp_lm``.
      !  Output is ``[r/sin(theta)*Br*U]=[(0,br_vt_lm,br_vp_lm)]``
      !

      !-- input:
      real(cp), intent(in) :: br(nlat_padded,nRstart:nRstop,n_phi_max)      ! :math:`r^2 B_r`
      real(cp), intent(in) :: vt(nlat_padded,nRstart:nRstop,n_phi_max)      ! :math:`r \sin\theta u_\theta`
      real(cp), intent(in) :: vp(nlat_padded,nRstart:nRstop,n_phi_max)      ! :math:`r \sin\theta u_\phi`
      real(cp), intent(in) :: omega          ! rotation rate of mantle or IC
      integer,  intent(in) :: nR

      !-- Output variables:
      ! br*vt/(sin(theta)**2*r**2)
      complex(cp), intent(inout) :: br_vt_lm(lm_max)
      ! br*(vp/(sin(theta)**2*r**2)-omega_ma)
      complex(cp), intent(inout) :: br_vp_lm(lm_max)

      !-- Local variables:
      integer :: n_theta, n_phi
      real(cp) :: br_vt(nlat_padded,n_phi_max), br_vp(nlat_padded,n_phi_max)
      real(cp) :: fac          ! 1/( r**2 sin(theta)**2 )

      fac=or2(nR)*orho1(nR)
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do
#else
      !$omp parallel do default(shared) &
      !$omp& private(n_theta,n_phi)
#endif
      do n_phi=1,n_phi_max
         do n_theta=1,n_theta_max
            br_vt(n_theta,n_phi)= fac*br(n_theta,nR,n_phi)*vt(n_theta,nR,n_phi)
            br_vp(n_theta,n_phi)= br(n_theta,nR,n_phi)* ( fac*vp(n_theta,nR,n_phi) - &
            &                                  omega*sinTheta_E2(n_theta) )
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#else
      !$omp end parallel do
#endif

      call spat_to_sphertor(sht_l_single, br_vt, br_vp, br_vt_lm, br_vp_lm, l_max)

   end subroutine get_br_v_bcs_batch
!----------------------------------------------------------------------------
   subroutine get_b_nl_bcs(bc,br_vt_lm,br_vp_lm,lm_min_b,lm_max_b,b_nl_bc,aj_nl_bc)
      !
      !  Purpose of this subroutine is to calculate the nonlinear term
      !  of the magnetic boundary condition for a conducting mantle in
      !  physical space (theta,phi), assuming that the conductance
      !  of the mantle is much smaller than that of the core.
      !  Calculation is performed for the theta block:
      !
      !  .. code-block:: fortran
      !
      !                  n_theta_min<=n_theta<=n_theta_min+n_theta_block-1
      !

      !-- Input variables:
      character(len=3), intent(in) :: bc                 ! Distinguishes 'CMB' and 'ICB'
      integer,          intent(in) :: lm_min_b,lm_max_b  ! limits of lm-block
      complex(cp),      intent(in) :: br_vt_lm(lm_max)  ! :math:`B_r u_\theta/(r^2\sin^2\theta)`
      complex(cp),      intent(in) :: br_vp_lm(lm_max)  ! :math:`B_r u_\phi/(r^2\sin^2\theta)`

      !-- Output variables:
      complex(cp), intent(out) :: b_nl_bc(lm_min_b:lm_max_b)  ! nonlinear bc for b
      complex(cp), intent(out) :: aj_nl_bc(lm_min_b:lm_max_b) ! nonlinear bc for aj

      !-- Local variables:
      integer :: lm        ! position of degree and order
      real(cp) :: fac

      if ( bc == 'CMB' ) then

         fac=conductance_ma*prmag
         !$omp parallel do default(shared)
         do lm=lm_min_b,lm_max_b
            b_nl_bc(lm) =-fac * br_vt_lm(lm)
            aj_nl_bc(lm)=-fac * br_vp_lm(lm)
         end do
         !$omp end parallel do

      else if ( bc == 'ICB' ) then

         fac=sigma_ratio*prmag
         !$omp parallel do default(shared) private(lm)
         do lm=lm_min_b,lm_max_b
            aj_nl_bc(lm)=-fac * br_vp_lm(lm)
         end do
         !$omp end parallel do

      else
         call abortRun('Wrong input of bc into get_b_nl_bcs')
      end if

   end subroutine get_b_nl_bcs
!-------------------------------------------------------------------------
   subroutine v_rigid_boundary(nR,omega,lDeriv,vrr,vtr,vpr,cvrr,dvrdtr, &
              &                dvrdpr,dvtdpr,dvpdpr)
      !
      !  Purpose of this subroutine is to set the velocities and their
      !  derivatives at a fixed boundary.
      !  While ``vt`` is zero, since we only allow for rotation about the
      !  :math:`z`-axis, ``vp= r \sin(theta) v_phi = r**2 sin(theta)**2 omega``
      !  and ``cvr= r**2 * radial component of (\curl v) = r**2  2 cos(theta) omega``
      !

      !-- Input of variables:
      integer,  intent(in) :: nR            ! no of radial grid point
      logical,  intent(in) :: lDeriv        ! derivatives required ?

      !-- Input of boundary rotation rate
      real(cp), intent(in) :: omega         ! boundary rotation rate

      !-- output:
      real(cp), intent(out) :: vrr(:,:), vpr(:,:), vtr(:,:)
      real(cp), intent(out) :: cvrr(:,:), dvrdtr(:,:), dvrdpr(:,:)
      real(cp), intent(out) :: dvtdpr(:,:), dvpdpr(:,:)

      !-- Local variables:
      real(cp) :: r2
      integer :: nTheta, nPhi

      if ( nR == n_r_cmb ) then
         r2=r_cmb*r_cmb
      else if ( nR == n_r_icb ) then
         r2=r_icb*r_icb
      else
         write(output_unit,*)
         write(output_unit,*) '! v_rigid boundary called for grid'
         write(output_unit,*) '! point which is not a boundary !  '
         return
      end if

#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do
#else
      !$omp parallel do default(shared) private(nPhi,nTheta)
#endif
      do nPhi=1,n_phi_max
         do nTheta=1,n_theta_max
            vrr(nTheta,nPhi)=0.0_cp
            vtr(nTheta,nPhi)=0.0_cp
            vpr(nTheta,nPhi)=r2*rho0(nR)*sinTheta_E2(nTheta)*omega
            if ( lDeriv ) then
               cvrr(nTheta,nPhi)  =r2*rho0(nR)*two*cosTheta(nTheta)*omega
               dvrdtr(nTheta,nPhi)=0.0_cp
               dvrdpr(nTheta,nPhi)=0.0_cp
               dvtdpr(nTheta,nPhi)=0.0_cp
               dvpdpr(nTheta,nPhi)=0.0_cp
            end if
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#else
      !$omp end parallel do
#endif

   end subroutine v_rigid_boundary

   subroutine v_rigid_boundary_batch(nR,omega,lDeriv,vrr,vtr,vpr,cvrr,dvrdtr, &
              &                dvrdpr,dvtdpr,dvpdpr)
      !
      !  Purpose of this subroutine is to set the velocities and their
      !  derivatives at a fixed boundary.
      !  While ``vt`` is zero, since we only allow for rotation about the
      !  :math:`z`-axis, ``vp= r \sin(theta) v_phi = r**2 sin(theta)**2 omega``
      !  and ``cvr= r**2 * radial component of (\curl v) = r**2  2 cos(theta) omega``
      !

      !-- Input of variables:
      integer,  intent(in) :: nR            ! no of radial grid point
      logical,  intent(in) :: lDeriv        ! derivatives required ?

      !-- Input of boundary rotation rate
      real(cp), intent(in) :: omega         ! boundary rotation rate

      !-- output:
      real(cp), intent(out) :: vrr(:,:,:), vpr(:,:,:), vtr(:,:,:)
      real(cp), intent(out) :: cvrr(:,:,:), dvrdtr(:,:,:), dvrdpr(:,:,:)
      real(cp), intent(out) :: dvtdpr(:,:,:), dvpdpr(:,:,:)

      !-- Local variables:
      real(cp) :: r2
      integer :: nTheta, nPhi, nelem

      if ( nR == n_r_cmb ) then
         r2=r_cmb*r_cmb
      else if ( nR == n_r_icb ) then
         r2=r_icb*r_icb
      else
         write(output_unit,*)
         write(output_unit,*) '! v_rigid boundary called for grid'
         write(output_unit,*) '! point which is not a boundary !  '
         return
      end if

#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do
#else
      !$omp parallel do default(shared) private(nPhi,nTheta)
#endif
      do nPhi=1,n_phi_max
         do nTheta=1,n_theta_max
            vrr(nTheta,nR,nPhi)=0.0_cp
            vtr(nTheta,nR,nPhi)=0.0_cp
            vpr(nTheta,nR,nPhi)=r2*rho0(nR)*sinTheta_E2(nTheta)*omega
            if ( lDeriv ) then
               cvrr(nTheta,nR,nPhi)  =r2*rho0(nR)*two*cosTheta(nTheta)*omega
               dvrdtr(nTheta,nR,nPhi)=0.0_cp
               dvrdpr(nTheta,nR,nPhi)=0.0_cp
               dvtdpr(nTheta,nR,nPhi)=0.0_cp
               dvpdpr(nTheta,nR,nPhi)=0.0_cp
            end if
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#else
      !$omp end parallel do
#endif

   end subroutine v_rigid_boundary_batch
!-------------------------------------------------------------------------
end module nonlinear_bcs
