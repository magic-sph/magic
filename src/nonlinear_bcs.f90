module nonlinear_bcs

   use iso_fortran_env, only: output_unit
   use precision_mod
   use blocking, only: lm2
   use truncation, only: n_phi_max, l_max, nlat_padded, lm_max
   use radial_data, only: n_r_cmb, n_r_icb
   use radial_functions, only: r_cmb, r_icb, rho0
   use physical_parameters, only: sigma_ratio, conductance_ma, prmag, oek
   use horizontal_data, only: cosTheta, sinTheta_E2, phi, sinTheta
   use constants, only: two, y10_norm, y11_norm, zero
   use useful, only: abortRun
   use sht, only: spat_to_sphertor

   implicit none

   private

   public :: get_br_v_bcs, get_b_nl_bcs, v_rigid_boundary, v_center_sphere

contains

   subroutine get_br_v_bcs(br,vt,vp,omega,O_r_E_2,O_rho,br_vt_lm,br_vp_lm)
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
      real(cp), intent(in) :: O_r_E_2      ! :math:`1/r^2`
      real(cp), intent(in) :: O_rho        ! :math:`1/\tilde{\rho}` (anelastic)

      !-- Output variables:
      complex(cp), intent(inout) :: br_vt_lm(:) ! br*vt/(sin(theta)**2*r**2)
      complex(cp), intent(inout) :: br_vp_lm(:) ! br*(vp/(sin(theta)**2*r**2)-omega_ma)

      !-- Local variables:
      integer :: n_phi
      real(cp) :: br_vt(nlat_padded,n_phi_max), br_vp(nlat_padded,n_phi_max)
      real(cp) :: fac          ! 1/( r**2 sin(theta)**2 )

      fac=O_r_E_2*O_rho
      !$omp parallel do default(shared)
      do n_phi=1,n_phi_max
         br_vt(:,n_phi)= fac*br(:,n_phi)*vt(:,n_phi)
         br_vp(:,n_phi)= br(:,n_phi)* ( fac*vp(:,n_phi)-omega*sinTheta_E2(:) )
      end do
      !$omp end parallel do

      call spat_to_sphertor(br_vt, br_vp, br_vt_lm, br_vp_lm, l_max)

   end subroutine get_br_v_bcs
!----------------------------------------------------------------------------
   subroutine get_b_nl_bcs(bc,br_vt_lm,br_vp_lm,b_nl_bc,aj_nl_bc)
      !
      !  Purpose of this subroutine is to calculate the nonlinear term
      !  of the magnetic boundary condition for a conducting mantle in
      !  physical space (theta,phi), assuming that the conductance
      !  of the mantle is much smaller than that of the core.
      !

      !-- Input variables:
      character(len=3), intent(in) :: bc                 ! Distinguishes 'CMB' and 'ICB'
      complex(cp),      intent(in) :: br_vt_lm(:)  ! :math:`B_r u_\theta/(r^2\sin^2\theta)`
      complex(cp),      intent(in) :: br_vp_lm(:)  ! :math:`B_r u_\phi/(r^2\sin^2\theta)`

      !-- Output variables:
      complex(cp), intent(out) :: b_nl_bc(:)  ! nonlinear bc for b
      complex(cp), intent(out) :: aj_nl_bc(:) ! nonlinear bc for aj

      !-- Local variables:
      integer :: lm        ! position of degree and order
      real(cp) :: fac

      if ( bc == 'CMB' ) then

         b_nl_bc(1) =zero                                             
         aj_nl_bc(1)=zero
         fac=conductance_ma*prmag
         !$omp parallel do default(shared)
         do lm=2,lm_max
            b_nl_bc(lm) =-fac * br_vt_lm(lm)
            aj_nl_bc(lm)=-fac * br_vp_lm(lm)
         end do
         !$omp end parallel do

      else if ( bc == 'ICB' ) then

         aj_nl_bc(1)=zero
         fac=sigma_ratio*prmag
         !$omp parallel do default(shared) private(lm)
         do lm=2,lm_max
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
      integer :: nPhi

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

      !$omp parallel do default(shared)
      do nPhi=1,n_phi_max
         vrr(:,nPhi)=0.0_cp
         vtr(:,nPhi)=0.0_cp
         vpr(:,nPhi)=r2*rho0(nR)*sinTheta_E2(:)*omega
         if ( lDeriv ) then
            cvrr(:,nPhi)  =r2*rho0(nR)*two*cosTheta(:)*omega
            dvrdtr(:,nPhi)=0.0_cp
            dvrdpr(:,nPhi)=0.0_cp
            dvtdpr(:,nPhi)=0.0_cp
            dvpdpr(:,nPhi)=0.0_cp
         end if
      end do
      !$omp end parallel do

   end subroutine v_rigid_boundary
!-------------------------------------------------------------------------
   subroutine v_center_sphere(ddw, vrr, vtr, vpr)
      !
      ! This routine is only called for full sphere computations to construct
      ! a vector field at the center of the the sphere. At the center, we have
      ! wlm \propto r^{l+1} and so vr = d2wlm/dr2 for l=1, 0 otherwise
      ! vtheta, vphi = sht(1/l*ddwlm, 0) for l=1, 0 otherwise
      !

      !-- Input variable
      complex(cp), intent(in) :: ddw(:)

      !-- Output variables:
      real(cp), intent(out) :: vrr(:,:), vtr(:,:), vpr(:,:)

      !-- Local variables:
      integer :: nPhi, lm10, lm11

      lm10=lm2(1,0)
      lm11=lm2(1,1)

      if ( lm11 > 0 ) then ! minc = 1
         !$omp parallel do default(shared)
         do nPhi=1,n_phi_max
            vrr(:,nPhi)=y10_norm*real(ddw(lm10))*cosTheta(:)  +       &
            &                two*y11_norm*sinTheta(:)*(               &
            &                        real(ddw(lm11))*cos(phi(nPhi))-  &
            &                       aimag(ddw(lm11))*sin(phi(nPhi)) )
            vtr(:,nPhi)=sinTheta(:)*(                                 &
            &                -y10_norm*real(ddw(lm10))*sinTheta(:) +  &
            &                two*y11_norm*cosTheta(:)*(               &
            &                        real(ddw(lm11))*cos(phi(nPhi))-  &
            &                       aimag(ddw(lm11))*sin(phi(nPhi)) ) )
            vpr(:,nPhi)=-two*y11_norm*sinTheta(:)*(                   &
            &                        real(ddw(lm11))*sin(phi(nPhi))+  &
            &                       aimag(ddw(lm11))*cos(phi(nPhi)) )
         end do
         !$omp end parallel do
      else ! minc /= 1
         !$omp parallel do default(shared)
         do nPhi=1,n_phi_max
            vrr(:,nPhi)=y10_norm*real(ddw(lm10))*cosTheta(:)
            vtr(:,nPhi)=-sinTheta(:)*y10_norm*real(ddw(lm10))*sinTheta(:)
            vpr(:,nPhi)=0.0_cp
         end do
         !$omp end parallel do
      end if

   end subroutine v_center_sphere
!-------------------------------------------------------------------------
end module nonlinear_bcs
