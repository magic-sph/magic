module nonlinear_bcs

   use iso_fortran_env, only: output_unit
   use precision_mod
   use LMmapping, only: map_dist_st
   use truncation, only: nrp, n_phi_max, l_axi, l_max, n_r_cmb, n_r_icb, &
       &                 nThetaStart, nThetaStop, n_lmP_loc, n_lm_loc
   use radial_functions, only: r_cmb, r_icb, rho0
   use physical_parameters, only: sigma_ratio, conductance_ma, prmag, oek
   use horizontal_data, only: dTheta1S_loc, dTheta1A_loc, dPhi_loc, O_sin_theta, &
       &                      dLh_loc, sn2, cosTheta
   use constants, only: two
   use shtns, only: spat_to_SH_dist
   use useful, only: abortRun

   implicit none

   private

   public :: get_br_v_bcs, get_b_nl_bcs, v_rigid_boundary

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
      !     n_theta_min<=n_theta<=n_theta_min+n_theta_block-1
      !
      !  On input br, vt and vp are given on all phi points and
      !  thetas in the specific block.
      !  On output the contribution of these grid points to all
      !  degree and orders is stored in br_vt_lm and br_vp_lm.
      !  Output is [r/sin(theta)*Br*U]=[(0,br_vt_lm,br_vp_lm)]
      !

      !-- input:
      real(cp), intent(in) :: br(nrp,nThetaStart:nThetaStop)      ! r**2 * B_r
      real(cp), intent(in) :: vt(nrp,nThetaStart:nThetaStop)      ! r*sin(theta) U_theta
      real(cp), intent(in) :: vp(nrp,nThetaStart:nThetaStop)      ! r*sin(theta) U_phi
      real(cp), intent(in) :: omega          ! rotation rate of mantle or IC
      real(cp), intent(in) :: O_r_E_2        ! 1/r**2
      real(cp), intent(in) :: O_rho          ! 1/rho0 (anelastic)

      !-- Output variables:
      complex(cp), intent(inout) :: br_vt_lm(n_lmP_loc) ! br*vt/(sin(theta)**2*r**2)
      complex(cp), intent(inout) :: br_vp_lm(n_lmP_loc) ! br*(vp/(sin(theta)**2*r**2)-omega_ma)

      !-- Local variables:
      integer :: n_theta      ! number of theta position
      integer :: n_phi        ! number of longitude
      real(cp) :: br_vt(nrp,nThetaStart:nThetaStop)
      real(cp) :: br_vp(nrp,nThetaStart:nThetaStop)
      real(cp) :: fac          ! 1/( r**2 sin(theta)**2 )

      !$omp parallel do default(shared) private(n_phi,fac,n_theta)
      do n_theta=nThetaStart,nThetaStop
         fac=O_sin_theta(n_theta)*O_sin_theta(n_theta)*O_r_E_2*O_rho

         do n_phi=1,n_phi_max
            br_vt(n_phi,n_theta)= fac*br(n_phi,n_theta)*vt(n_phi,n_theta)

            br_vp(n_phi,n_theta)= br(n_phi,n_theta) * &
            &                        ( fac*vp(n_phi,n_theta) - omega )
         end do
      end do
      !$omp end parallel do

      !-- Fourier transform phi 2 m (real 2 complex!)
      call spat_to_SH_dist(br_vt, br_vt_lm, l_max)
      call spat_to_SH_dist(br_vp, br_vp_lm, l_max)

   end subroutine get_br_v_bcs
!----------------------------------------------------------------------------
   subroutine get_b_nl_bcs(bc,br_vt_lm,br_vp_lm,b_nl_bc,aj_nl_bc)
      !
      !  Purpose of this subroutine is to calculate the nonlinear term
      !  of the magnetic boundary condition for a conducting mantle in
      !  physical space (phi,theta), assuming that the conductance
      !  of the mantle is much smaller than that of the core.
      !  Calculation is performed for the theta block:
      !
      !  .. code-block:: fortran
      !
      !      n_theta_min<=n_theta<=n_theta_min+n_theta_block-1
      !

      !-- Input variables:
      character(len=3), intent(in) :: bc                 ! Distinguishes 'CMB' and 'ICB'
      complex(cp),      intent(in) :: br_vt_lm(n_lmP_loc)  ! [br*vt/(r**2*sin(theta)**2)]
      complex(cp),      intent(in) :: br_vp_lm(n_lmP_loc)  ! [br*vp/(r**2*sin(theta)**2)

      !-- Output variables:
      complex(cp), intent(out) :: b_nl_bc(n_lm_loc)  ! nonlinear bc for b
      complex(cp), intent(out) :: aj_nl_bc(n_lm_loc) ! nonlinear bc for aj

      !-- Local variables:
      integer :: l,m       ! degree and order
      integer :: lm        ! position of degree and order
      integer :: lmP       ! same as lm but for l running to l_max+1
      integer :: lmPS,lmPA ! lmP for l-1 and l+1
      real(cp) :: fac

      if ( bc == 'CMB' ) then

         fac=conductance_ma*prmag

         !$omp parallel do default(shared) private(lm,l,m,lmP,lmPS,lmPA)
         do lm=1,n_lm_loc
            l   =map_dist_st%lm2l(lm)
            if ( l == 0 ) cycle
            m   =map_dist_st%lm2m(lm)
            lmP =map_dist_st%lm2lmP(lm)
            lmPS=map_dist_st%lmP2lmPS(lmP)
            lmPA=map_dist_st%lmP2lmPA(lmP)
            if ( l > m ) then
               b_nl_bc(lm)= fac/dLh_loc(lm) * (  dTheta1S_loc(lm)*br_vt_lm(lmPS)  &
               &                               - dTheta1A_loc(lm)*br_vt_lm(lmPA)  &
               &                               +     dPhi_loc(lm)*br_vp_lm(lmP)   )
               aj_nl_bc(lm)=-fac/dLh_loc(lm) * ( dTheta1S_loc(lm)*br_vp_lm(lmPS)  &
               &                               - dTheta1A_loc(lm)*br_vp_lm(lmPA)  &
               &                               -     dPhi_loc(lm)*br_vt_lm(lmP)    )
            else if ( l == m ) then
               b_nl_bc(lm)= fac/dLh_loc(lm) * ( - dTheta1A_loc(lm)*br_vt_lm(lmPA)  &
               &                                    + dPhi_loc(lm)*br_vp_lm(lmP)   )
               aj_nl_bc(lm)=-fac/dLh_loc(lm) * ( - dTheta1A_loc(lm)*br_vp_lm(lmPA) &
               &                                     - dPhi_loc(lm)*br_vt_lm(lmP)   )
            end if
         end do
         !$omp end parallel do

      else if ( bc == 'ICB' ) then

         fac=sigma_ratio*prmag
         !$omp parallel do default(shared) private(lm,l,m,lmP,lmPS,lmPA)
         do lm=1,n_lm_loc
            l   =map_dist_st%lm2l(lm)
            if ( l == 0 ) cycle
            m   =map_dist_st%lm2m(lm)
            lmP =map_dist_st%lm2lmP(lm)
            lmPS=map_dist_st%lmP2lmPS(lmP)
            lmPA=map_dist_st%lmP2lmPA(lmP)
            if ( l > m ) then
               aj_nl_bc(lm)=-fac/dLh_loc(lm) * ( dTheta1S_loc(lm)*br_vp_lm(lmPS)   &
               &                               - dTheta1A_loc(lm)*br_vp_lm(lmPA)   &
               &                               -     dPhi_loc(lm)*br_vt_lm(lmP)    )
            else if ( l == m ) then
               aj_nl_bc(lm)=-fac/dLh_loc(lm) * (- dTheta1A_loc(lm)*br_vp_lm(lmPA) &
               &                                    - dPhi_loc(lm)*br_vt_lm(lmP)    )
            end if
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
      !  While vt is zero, since we only allow for rotation about the
      !  z-axis, vp= r sin(theta) v_phi = r**2 sin(theta)**2 omega
      !  cvr= r**2 * radial component of (\curl v) =
      !  r**2  2 cos(theta) omega
      !

      !-- Input of variables:
      integer,  intent(in) :: nR            ! no of radial grid point
      logical,  intent(in) :: lDeriv        ! derivatives required ?

      !-- Input of boundary rotation rate
      real(cp), intent(in) :: omega

      !-- output:
      real(cp), intent(out) :: vrr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(out) :: vpr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(out) :: vtr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(out) :: cvrr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(out) :: dvrdtr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(out) :: dvrdpr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(out) :: dvtdpr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(out) :: dvpdpr(nrp,nThetaStart:nThetaStop)

      !-- Local variables:
      real(cp) :: r2
      integer :: nTheta,nThetaNHS
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

      !$omp parallel do default(shared) private(nTheta,nThetaNHS,nPhi)
      do nTheta=nThetaStart,nThetaStop
         nThetaNHS =(nTheta+1)/2 ! northern hemisphere,sn2 has size n_theta_max/2
         do nPhi=1,n_phi_max
            vrr(nPhi,nTheta)=0.0_cp
            vtr(nPhi,nTheta)=0.0_cp
            vpr(nPhi,nTheta)=r2*rho0(nR)*sn2(nThetaNHS)*omega
            if ( lDeriv ) then
               cvrr(nPhi,nTheta)  =r2*rho0(nR)*two*cosTheta(nTheta)*omega
               dvrdtr(nPhi,nTheta)=0.0_cp
               dvrdpr(nPhi,nTheta)=0.0_cp
               dvtdpr(nPhi,nTheta)=0.0_cp
               dvpdpr(nPhi,nTheta)=0.0_cp
            end if
         end do
      end do
      !$omp end parallel do

   end subroutine v_rigid_boundary
!-------------------------------------------------------------------------
end module nonlinear_bcs
