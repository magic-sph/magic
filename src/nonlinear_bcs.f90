!$Id$
module nonlinear_bcs

   use precision_mod
   use truncation, only: nrp, lmP_max, n_phi_max
   use radial_data, only: n_r_cmb, n_r_icb
   use radial_functions, only: r_cmb, r_icb, rho0
   use blocking, only: lm2l, lm2m, lm2lmP, lmP2lmPS, lmP2lmPA, nfs, &
                       sizeThetaB
   use physical_parameters, only: sigma_ratio, conductance_ma, prmag
   use horizontal_data, only: dTheta1S, dTheta1A, dPhi, O_sin_theta, &
                              dLh, sn2, cosTheta
#if (FFTLIB==JW)
   use fft_JW, only: fft_thetab
#elif (FFTLIB==MKL)
   use fft_MKL, only: fft_thetab
#endif
   use legendre_grid_to_spec, only: legTF2
   use const, only: two

   implicit none

   private

   public :: get_br_v_bcs, get_b_nl_bcs, v_rigid_boundary

contains

   subroutine get_br_v_bcs(br,vt,vp,omega,O_r_E_2,O_rho, &
        &                  n_theta_min,n_theta_block, &
        &                  br_vt_lm,br_vp_lm)
      !  +-------------------------------------------------------------------+
      !  |                                                                   |
      !  |  Purpose of this subroutine is to calculate the nonlinear term    |
      !  |  of the magnetic boundary condition for a conducting mantle or    |
      !  |  inner core in space (r,lm).                                      |
      !  |  Calculation is performed for the theta block:                    |
      !  |          n_theta_min<=n_theta<=n_theta_min+n_theta_block-1        |
      !  |  On input br, vt and vp are given on all phi points and           |
      !  |  thetas in the specific block.                                    |
      !  |  On output the contribution of these grid points to all           |
      !  |  degree and orders is stored in br_vt_lm and br_vp_lm.            |
      !  |  Output is [r/sin(theta)*Br*U]=[(0,br_vt_lm,br_vp_lm)]            |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+
    
      !-- input:
      real(cp), intent(in) :: br(nrp,*)      ! r**2 * B_r
      real(cp), intent(in) :: vt(nrp,*)      ! r*sin(theta) U_theta
      real(cp), intent(in) :: vp(nrp,*)      ! r*sin(theta) U_phi
      real(cp), intent(in) :: omega          ! rotation rate of mantle or IC
      real(cp), intent(in) :: O_r_E_2        ! 1/r**2
      real(cp), intent(in) :: O_rho          ! 1/rho0 (anelastic)
      integer,  intent(in) :: n_theta_min    ! start of theta block
      integer,  intent(in) :: n_theta_block  ! size of theta_block
    
      !-- Output variables:
      ! br*vt/(sin(theta)**2*r**2)
      complex(cp), intent(inout) :: br_vt_lm(lmP_max)
      ! br*(vp/(sin(theta)**2*r**2)-omega_ma)
      complex(cp), intent(inout) :: br_vp_lm(lmP_max)
    
      !-- Local variables:
      integer :: n_theta     ! number of theta position
      integer :: n_theta_rel ! number of theta position in block
      integer :: n_phi       ! number of longitude
      real(cp) :: br_vt(nrp,n_theta_block)
      real(cp) :: br_vp(nrp,n_theta_block)
      real(cp) :: fac          ! 1/( r**2 sin(theta)**2 )
    
      n_theta=n_theta_min-1 ! n_theta needed for O_sin_theta_E_2
    
      do n_theta_rel=1,n_theta_block
         n_theta=n_theta+1           ! absolute number of theta
    
         fac=O_sin_theta(n_theta)*O_sin_theta(n_theta)*O_r_E_2*O_rho
    
         do n_phi=1,n_phi_max
            br_vt(n_phi,n_theta_rel)= fac*br(n_phi,n_theta_rel)*vt(n_phi,n_theta_rel)
    
            br_vp(n_phi,n_theta_rel)= br(n_phi,n_theta_rel) * &
                                     ( fac*vp(n_phi,n_theta_rel) - omega )
         end do
      end do
    
      !-- Fourier transform phi 2 m (real 2 complex!)
      call fft_thetab(br_vt,-1)
      call fft_thetab(br_vp,-1)
    
      !-- Legendre transform contribution of thetas in block:
      call legTF2(n_theta_min,br_vt_lm,br_vp_lm,br_vt,br_vp)
    
   end subroutine get_br_v_bcs
!----------------------------------------------------------------------------
   subroutine get_b_nl_bcs(bc,br_vt_lm,br_vp_lm, &
                           lm_min_b,lm_max_b,b_nl_bc,aj_nl_bc)
      !  +-------------------------------------------------------------------+
      !  |                                                                   |
      !  |  Purpose of this subroutine is to calculate the nonlinear term    |
      !  |  of the magnetic boundary condition for a conducting mantle in    |
      !  |  physical space (phi,theta), assuming that the conductance        |
      !  |  of the mantle is much smaller than that of the core.             |
      !  |  Calculation is performed for the theta block:                    |
      !  |          n_theta_min<=n_theta<=n_theta_min+n_theta_block-1        |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+
         
      !-- Input variables:
      character(len=3), intent(in) :: bc                 ! Distinguishes 'CMB' and 'ICB'
      integer,          intent(in) :: lm_min_b,lm_max_b  ! limits of lm-block
      complex(cp),      intent(in) :: br_vt_lm(lmP_max)  ! [br*vt/(r**2*sin(theta)**2)]
      complex(cp),      intent(in) :: br_vp_lm(lmP_max)  ! [br*vp/(r**2*sin(theta)**2)

      !-- Output variables:
      complex(cp), intent(out) :: b_nl_bc(lm_min_b:lm_max_b)  ! nonlinear bc for b
      complex(cp), intent(out) :: aj_nl_bc(lm_min_b:lm_max_b) ! nonlinear bc for aj

      !-- Local variables:
      integer :: l,m       ! degree and order
      integer :: lm        ! position of degree and order
      integer :: lmP       ! same as lm but for l running to l_max+1
      integer :: lmPS,lmPA ! lmP for l-1 and l+1
      real(cp) :: fac

      write(*,"(2A)") "In get_b_nl_bcs with bc=",bc

      if ( bc == 'CMB' ) then

         fac=conductance_ma*prmag

         do lm=lm_min_b,lm_max_b
            l   =lm2l(lm)
            m   =lm2m(lm)
            lmP =lm2lmP(lm)
            lmPS=lmP2lmPS(lmP)
            lmPA=lmP2lmPA(lmP)
            if ( l > m ) then
               b_nl_bc(lm)= fac/dLh(lm) * (  dTheta1S(lm)*br_vt_lm(lmPS)  &
                                           - dTheta1A(lm)*br_vt_lm(lmPA)  &
                                           +     dPhi(lm)*br_vp_lm(lmP)   )
               aj_nl_bc(lm)=-fac/dLh(lm) * ( dTheta1S(lm)*br_vp_lm(lmPS)  &
                                           - dTheta1A(lm)*br_vp_lm(lmPA)  &
                                           -     dPhi(lm)*br_vt_lm(lmP)    )
            else if ( l == m ) then
               b_nl_bc(lm)= fac/dLh(lm) * ( - dTheta1A(lm)*br_vt_lm(lmPA)  &
                                                + dPhi(lm)*br_vp_lm(lmP)   )
               aj_nl_bc(lm)=-fac/dLh(lm) * ( - dTheta1A(lm)*br_vp_lm(lmPA) &
                                                 - dPhi(lm)*br_vt_lm(lmP)   )
            end if
         end do

      else if ( bc == 'ICB' ) then

         fac=sigma_ratio*prmag

         do lm=lm_min_b,lm_max_b
            l   =lm2l(lm)
            m   =lm2m(lm)
            lmP =lm2lmP(lm)
            lmPS=lmP2lmPS(lmP)
            lmPA=lmP2lmPA(lmP)
            if ( l > m ) then
               aj_nl_bc(lm)=-fac/dLh(lm) * ( dTheta1S(lm)*br_vp_lm(lmPS)   &
                                           - dTheta1A(lm)*br_vp_lm(lmPA)   &
                                           -     dPhi(lm)*br_vt_lm(lmP)    )
            else if ( l == m ) then
               aj_nl_bc(lm)=-fac/dLh(lm) * (- dTheta1A(lm)*br_vp_lm(lmPA) &
                                                - dPhi(lm)*br_vt_lm(lmP)    )
            end if
         end do

      else
         write(*,*) 'Wrong input of bc into s_get_b_nl_bcs.f'
         stop
      end if
               
   end subroutine get_b_nl_bcs
!-------------------------------------------------------------------------
   subroutine v_rigid_boundary(nR,omega,lDeriv,vrr,vtr,vpr,      &
            &                  cvrr,dvrdtr,dvrdpr,dvtdpr,dvpdpr, &
            &                  nThetaStart)
      !  +-------------------------------------------------------------------+
      !  |                                                                   |
      !  |  Purpose of this subroutine is to set the velocities and their    |
      !  |  derivatives at a fixed boundary.                                 |
      !  |  While vt is zero, since we only allow for rotation about the     |
      !  |  z-axix, vp= r sin(theta) v_phi = r**2 sin(theta)**2 omega        |
      !  |  cvr= r**2 * radial component of (\curl v) =                      |
      !  |       r**2  2 cos(theta) omega                                    |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+

      !-- Input of variables:
      integer,  intent(in) :: nR            ! no of radial grid point
      logical,  intent(in) :: lDeriv        ! derivatives required ?
      integer,  intent(in) :: nThetaStart   ! no of theta to start with
              
      !-- Input of boundary rotation rate
      real(cp), intent(in) :: omega

      !-- output:
      real(cp), intent(out) :: vrr(nrp,nfs)
      real(cp), intent(out) :: vpr(nrp,nfs)
      real(cp), intent(out) :: vtr(nrp,nfs)
      real(cp), intent(out) :: cvrr(nrp,nfs)
      real(cp), intent(out) :: dvrdtr(nrp,nfs)
      real(cp), intent(out) :: dvrdpr(nrp,nfs)
      real(cp), intent(out) :: dvtdpr(nrp,nfs)
      real(cp), intent(out) :: dvpdpr(nrp,nfs)

      !-- Local variables:
      real(cp) :: r2
      integer :: nThetaCalc,nTheta,nThetaNHS
      integer :: nPhi

      if ( nR == n_r_cmb ) then
         r2=r_cmb*r_cmb
      else if ( nR == n_r_icb ) then
         r2=r_icb*r_icb
      else
         write(*,*)
         write(*,*) '! v_rigid boundary called for a grid'
         write(*,*) '! points which is not a boundary !  '
         return
      end if

      nThetaCalc=nThetaStart-1
      do nTheta=1,sizeThetaB
         nThetaCalc=nThetaCalc+1
         nThetaNHS =(nThetaCalc+1)/2 ! northern hemisphere=odd n_theta
         do nPhi=1,n_phi_max
            vrr(nPhi,nTheta)=0.0_cp
            vtr(nPhi,nTheta)=0.0_cp
            vpr(nPhi,nTheta)=r2*rho0(nR)*sn2(nThetaNHS)*omega
            if ( lDeriv ) then
               cvrr(nPhi,nTheta)  =r2*rho0(nR)*two*cosTheta(nThetaCalc)*omega
               dvrdtr(nPhi,nTheta)=0.0_cp
               dvrdpr(nPhi,nTheta)=0.0_cp
               dvtdpr(nPhi,nTheta)=0.0_cp
               dvpdpr(nPhi,nTheta)=0.0_cp
            end if
         end do
      end do

   end subroutine v_rigid_boundary
!-------------------------------------------------------------------------
end module nonlinear_bcs
