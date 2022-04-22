module outRot
   !
   ! This module handles the writing of several diagnostic files related
   ! to the rotation: angular momentum (AM.TAG), drift (drift.TAG), inner
   ! core and mantle rotations.
   !

   use parallel_mod
   use precision_mod
   use communications, only: allgather_from_rloc, send_lm_pair_to_master
   use truncation, only: n_r_max, n_r_maxMag, minc, n_phi_max, n_theta_max
   use grid_blocking, only: radlatlon2spat
   use radial_data, only: n_r_cmb, n_r_icb, nRstart, nRstop
   use radial_functions, only: r_icb, r_cmb, r, rscheme_oc, beta, visc
   use physical_parameters, only: kbotv, ktopv, LFfac
   use num_param, only: lScale, tScale, vScale
   use blocking, only: lo_map, lm_balance, llm, ulm, llmMag, ulmMag
   use logic, only: l_AM, l_save_out, l_iner, l_SRIC, l_rot_ic, &
       &            l_SRMA, l_rot_ma, l_mag_LF, l_mag, l_drift, &
       &            l_finite_diff, l_full_sphere
   use output_data, only: tag
   use constants, only: c_moi_oc, c_moi_ma, c_moi_ic, pi, y11_norm, &
       &            y10_norm, zero, two, third, four, half
   use integration, only: rInt_R
   use horizontal_data, only: cosTheta, gauss
   use special, only: BIC, lGrenoble
   use useful, only: abortRun

   implicit none

   private

   integer :: n_SRMA_file, n_SRIC_file
   integer :: n_angular_file, n_rot_file
   integer :: n_inerP_file, n_inerT_file
   integer :: n_driftVD_file, n_driftVQ_file
   integer :: n_driftBD_file, n_driftBQ_file
   character(len=72) :: SRMA_file, SRIC_file, rot_file, angular_file
   character(len=72) :: inerP_file, inerT_file
   character(len=72) :: driftVD_file, driftVQ_file
   character(len=72) :: driftBD_file, driftBQ_file

   public :: write_rot, get_viscous_torque, get_angular_moment, get_angular_moment_Rloc, &
   &         get_lorentz_torque, initialize_outRot, finalize_outRot

contains

   subroutine initialize_outRot

      SRIC_file   ='SRIC.'//tag
      SRMA_file   ='SRMA.'//tag
      rot_file    ='rot.'//tag
      angular_file='AM.'//tag
      driftVD_file='driftVD.'//tag
      driftVQ_file='driftVQ.'//tag
      driftBD_file='driftBD.'//tag
      driftBQ_file='driftBQ.'//tag
      inerP_file  ='inerP.'//tag
      inerT_file  ='inerT.'//tag

      if ( rank == 0 .and. (.not. l_save_out) ) then

         if ( l_SRIC ) then
            open(newunit=n_SRIC_file, file=SRIC_file, status='new')
         end if

         if ( l_SRMA ) then
            open(newunit=n_SRMA_file, file=SRMA_file, status='new')
         end if

         if ( .not. l_SRIC .and. .not. l_SRMA ) then
            if ( l_rot_ic .or. l_rot_ma ) then
               open(newunit=n_rot_file, file=rot_file, status='new')
            end if
         end if

         if ( l_AM ) then
            open(newunit=n_angular_file, file=angular_file, status='new')
         end if

         if ( l_drift ) then
            open(newunit=n_driftVD_file, file=driftVD_file, status='new')
            open(newunit=n_driftVQ_file, file=driftVQ_file, status='new')
            if ( l_mag ) then
               open(newunit=n_driftBD_file, file=driftBD_file, status='new')
               open(newunit=n_driftBQ_file, file=driftBQ_file, status='new')
            end if
         end if

         if ( l_iner ) then
            open(newunit=n_inerP_file, file=inerP_file, status='new')
            open(newunit=n_inerT_file, file=inerT_file, status='new')
         end if

      end if

   end subroutine initialize_outRot
!-----------------------------------------------------------------------
   subroutine finalize_outRot

      if ( rank == 0 .and. (.not. l_save_out) ) then
         if ( l_SRIC ) close(n_SRIC_file)
         if ( l_SRMA ) close(n_SRMA_file)
         if ( l_AM ) close(n_angular_file)
         if ( l_iner ) then
            close(n_inerT_file)
            close(n_inerP_file)
         end if
         if ( .not. l_SRIC .and. .not. l_SRMA ) then
            if ( l_rot_ic .or. l_rot_ma ) then
               close(n_rot_file)
            end if
         end if
         if ( l_drift ) then
            close(n_driftVD_file)
            close(n_driftVQ_file)
            if ( l_mag ) then
               close(n_driftBD_file)
               close(n_driftBQ_file)
            end if
         end if
      end if

   end subroutine finalize_outRot
!-----------------------------------------------------------------------
   subroutine write_rot(time,dt,eKinIC,ekinMA,w,z,dz,b,omega_ic,omega_ma, &
              &         lorentz_torque_ic,lorentz_torque_ma)

      !-- Input of variables:
      real(cp),    intent(in) :: omega_ic,omega_ma
      real(cp),    intent(in) :: lorentz_torque_ma,lorentz_torque_ic
      real(cp),    intent(in) :: time,dt
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(in) :: z(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dz(llm:ulm,n_r_max)
      complex(cp), intent(in) :: b(llmMag:ulmMag,n_r_maxMag)

      !-- Output into rot_file
      real(cp), intent(out) :: eKinIC,eKinMA

      !-- Local variables:
      real(cp), parameter :: tolerance=10.0_cp*epsilon(0.0_cp)
      real(cp) :: eKinOC
      integer :: n_r1,n_r2,n_r3,l1m0
      real(cp) :: viscous_torque_ic,viscous_torque_ma
      real(cp) :: AMz,eKinAMz
      real(cp) :: angular_moment_oc(3)
      real(cp) :: angular_moment_ic(3)
      real(cp) :: angular_moment_ma(3)
      complex(cp) :: z10(n_r_max),z11(n_r_max)

      real(cp) :: powerLor,powerVis
      real(cp), save :: AMzLast=0.0_cp,eKinAMzLast=0.0_cp

      integer, pointer :: lm2(:,:)
      integer :: i,l,m,ilm,n_lm_vals
      complex(cp) :: zvals_on_rank0(8,3),bvals_on_rank0(8,3)
      complex(cp) :: vals_on_rank0_1d(21)

      lm2(0:,0:) => lo_map%lm2
      l1m0=lm2(1,0)

      if ( llm <= l1m0 .and. ulm >= l1m0 ) then
         !-- Calculating viscous torques:
         if ( l_rot_ic .and. kbotv == 2 ) then
            call get_viscous_torque(viscous_torque_ic,real(z(l1m0,n_r_max)),    &
                 &                  real(dz(l1m0,n_r_max)),r_icb,beta(n_r_max), &
                 &                  visc(n_r_max))
         else
            viscous_torque_ic=0.0_cp
         end if
         if ( l_rot_ma .and. ktopv == 2 ) then
            call get_viscous_torque(viscous_torque_ma,real(z(l1m0,1)), &
                 &                  real(dz(l1m0,1)),r_cmb,beta(1),visc(1))
         else
            viscous_torque_ma=0.0_cp
         end if
      end if

      call send_lm_pair_to_master(viscous_torque_ic,1,0)
      call send_lm_pair_to_master(viscous_torque_ma,1,0)

      if ( rank == 0 ) then
         if ( l_SRIC ) then
            powerLor=lorentz_torque_ic*omega_IC
            powerVis=viscous_torque_ic*omega_IC
            if ( l_save_out ) then
               open(newunit=n_SRIC_file, file=SRIC_file, status='unknown', &
               &    position='append')
            end if
            write(n_SRIC_file,'(1p,2x,ES20.12,4ES17.6)')    &
            &     time*tScale,omega_ic/tScale,              &
            &     (powerLor+powerVis)*vScale*vScale/tScale, &
            &     powerVis*vScale*vScale/tScale,            &
            &     powerLor*vScale*vScale/tScale
            if ( l_save_out ) close(n_SRIC_file)
         end if
         if ( l_SRMA ) then
            powerLor=lorentz_torque_ma*omega_ma
            powerVis=viscous_torque_ma*omega_ma
            if ( l_save_out ) then
               open(newunit=n_SRMA_file, file=SRMA_file, status='unknown', &
               &    position='append')
            end if
            write(n_SRMA_file,'(1p,2x,ES20.12,4ES17.6)')    &
            &     time*tScale, omega_ma/tScale,             &
            &     (powerLor+powerVis)*vScale*vScale/tScale, &
            &     powerVis*vScale*vScale/tScale,            &
            &     powerLor*vScale*vScale/tScale
            if ( l_save_out ) close(n_SRMA_file)
         end if
      end if

      if ( l_drift ) then
         n_r1=int(third*(n_r_max-1))
         n_r2=int(two*third*(n_r_max-1))
         n_r3=n_r_max-1

         do i=1,4
            call send_lm_pair_to_master(z(:,n_r1),i*minc,i*minc,zvals_on_rank0(i,1))
            call send_lm_pair_to_master(z(:,n_r2),i*minc,i*minc,zvals_on_rank0(i,2))
            call send_lm_pair_to_master(z(:,n_r3),i*minc,i*minc,zvals_on_rank0(i,3))
            call send_lm_pair_to_master(z(:,n_r1),i*minc+1,i*minc,zvals_on_rank0(4+i,1))
            call send_lm_pair_to_master(z(:,n_r2),i*minc+1,i*minc,zvals_on_rank0(4+i,2))
            call send_lm_pair_to_master(z(:,n_r3),i*minc+1,i*minc,zvals_on_rank0(4+i,3))
         end do

         if ( rank == 0 ) then
            if ( l_save_out ) then
               open(newunit=n_driftVD_file, file=driftVD_file, status='unknown', &
               &    position='append')
               open(newunit=n_driftVQ_file, file=driftVQ_file, status='unknown', &
               &    position='append')
            end if
            write(n_driftVD_file,'(1P,2X,ES20.12,24ES12.4)')      &
            &     time*tScale, (zvals_on_rank0(ilm,1),ilm=1,4),   &
            &     (zvals_on_rank0(ilm,2),ilm=1,4),                &
            &     (zvals_on_rank0(ilm,3),ilm=1,4)
            write(n_driftVQ_file,'(1P,2X,ES20.12,24ES12.4)')      &
            &     time*tScale, (zvals_on_rank0(ilm,1),ilm=5,8),   &
            &     (zvals_on_rank0(ilm,2),ilm=5,8),                &
            &     (zvals_on_rank0(ilm,3),ilm=5,8)
            if ( l_save_out ) then
               close(n_driftVD_file)
               close(n_driftVQ_file)
            end if
         end if

         if ( l_mag .or. l_mag_LF ) then
            n_r1=n_r_cmb
            n_r2=n_r_icb
            do i=1,4
               call send_lm_pair_to_master(b(:,n_r1),i*minc,i*minc, &
                    &                      bvals_on_rank0(i,1))
               call send_lm_pair_to_master(b(:,n_r2),i*minc,i*minc, &
                    &                      bvals_on_rank0(i,2))
               call send_lm_pair_to_master(b(:,n_r1),i*minc+1,i*minc,&
                    &                      bvals_on_rank0(4+i,1))
               call send_lm_pair_to_master(b(:,n_r2),i*minc+1,i*minc,&
                    &                      bvals_on_rank0(4+i,2))
            end do

            if ( rank == 0 ) then
               if ( l_save_out ) then
                  open(newunit=n_driftBD_file, file=driftBD_file, status='unknown',&
                  &    position='append')
                  open(newunit=n_driftBQ_file, file=driftBQ_file, status='unknown',&
                  &    position='append')
               end if
               write(n_driftBD_file,'(1P,2X,ES20.12,16ES12.4)')     &
               &     time*tScale, (bvals_on_rank0(ilm,1),ilm=5,8),  &
               &     (bvals_on_rank0(ilm,2),ilm=5,8)
               write(n_driftBQ_file,'(1P,2X,ES20.12,16ES12.4)')     &
               &     time*tScale, (bvals_on_rank0(ilm,1),ilm=1,4),  &
               &     (bvals_on_rank0(ilm,2),ilm=1,4)
               if ( l_save_out ) then
                  close(n_driftBD_file)
                  close(n_driftBQ_file)
               end if
            end if
         end if ! l_mag
      end if

      if ( .not. l_SRIC .and. ( l_rot_ic .or. l_rot_ma ) ) then
         if ( rank == 0 ) then
            if ( l_save_out ) then
               open(newunit=n_rot_file, file=rot_file, status='unknown', &
               &    position='append')
            end if
            write(n_rot_file,'(1P,2X,ES20.12,6ES16.8)')  &
            &     time*tScale, omega_ic/tScale,          &
            &     lScale**2*vScale*lorentz_torque_ic,    &
            &     lScale**2*vScale*viscous_torque_ic,    &
            &     omega_ma/tScale,                       &
            &     lScale**2*vScale*lorentz_torque_ma,    &
            &     -lScale**2*vScale*viscous_torque_ma
            if ( l_save_out ) close(n_rot_file)
         end if
      end if

      if ( l_AM ) then
         call send_lm_pair_to_master(z,1,0,z10)
         call send_lm_pair_to_master(z,1,1,z11)

         if ( rank == 0 ) then
            call get_angular_moment(z10,z11,omega_ic,omega_ma,angular_moment_oc, &
                 &                  angular_moment_ic,angular_moment_ma)
            if ( l_save_out ) then
               open(newunit=n_angular_file, file=angular_file, status='unknown', &
               &    position='append')
            end if
            AMz=angular_moment_oc(3)+angular_moment_ic(3)+angular_moment_ma(3)
            if ( abs(AMz) < tolerance ) AMz=0.0_cp
            if ( l_full_sphere ) then
               eKinAMz=half*(angular_moment_oc(3)**2/c_moi_oc + &
               &             angular_moment_ma(3)**2/c_moi_ma )
            else
               eKinAMz=half*(angular_moment_oc(3)**2/c_moi_oc + &
               &             angular_moment_ic(3)**2/c_moi_ic + &
               &             angular_moment_ma(3)**2/c_moi_ma )
            end if
            if ( abs(eKinAMz) < tolerance ) eKinAMz=0.0_cp
            if ( l_full_sphere ) then
               eKinIC = 0.0_cp
            else
               eKinIC=half*angular_moment_ic(3)**2/c_moi_ic
            end if
            eKinOC=half*angular_moment_oc(3)**2/c_moi_oc
            eKinMA=half*angular_moment_ma(3)**2/c_moi_ma
            if ( AMzLast /= 0.0_cp ) then
               !write(*,"(A,4ES22.15)") "col9 = ",eKinAMz,eKinAMzLast, &
               !     &                  dt,(eKinAMz-eKinAMzLast)
               write(n_angular_file,'(1p,2x,ES20.12,5ES14.6,3ES20.12)', advance='no') &
               &     time*tScale, angular_moment_oc,                                  &
               &     angular_moment_ic(3), angular_moment_ma(3),                      &
               &     AMz,(AMz-AMzLast)/AMzLast/dt,eKinAMz
               if (eKinAMzLast /= 0.0_cp) then
                  write(n_angular_file,'(1ES20.12)', advance='no') &
                  &     (eKinAMz-eKinAMzLast)/eKinAMzLast/dt
               else
                  write(n_angular_file,'(1ES20.12)', advance='no') 0.0
               end if
               write(n_angular_file,'(3ES20.12)') eKinIC,eKinOC,eKinMA
            end if
            if ( l_save_out ) close(n_angular_file)
            AMzLast=AMz
            eKinAMzLast=eKinAMz
         end if
      end if

      if ( l_iner ) then
         ! l_iner can only be .true. for minc=1
         n_r1=int(half*(n_r_max-1))
         ilm=0
         do l=1,6
            do m=1,l
               ilm = ilm + 1
               call send_lm_pair_to_master(w(:,n_r1),l,m,vals_on_rank0_1d(ilm))
            end do
         end do
         n_lm_vals=ilm

         if ( rank == 0 ) then
            if ( l_save_out ) then
               open(newunit=n_inerP_file, file=inerP_file, status='unknown', &
               &    position='append')
            end if
            write(n_inerP_file,'(1P,2X,ES20.12,21ES12.4)') &
            &    time*tScale, ( real(vals_on_rank0_1d(ilm)),ilm=1,n_lm_vals )
            if ( l_save_out ) close(n_inerP_file)
         end if

         n_r1=int(half*(n_r_max-1))
         ilm=0
         do l=1,6
            do m=1,l
               ilm = ilm+1
               call send_lm_pair_to_master(z(:,n_r1),l,m,vals_on_rank0_1d(ilm))
            end do
         end do

         if ( rank == 0 ) then
            if ( l_save_out ) then
               open(newunit=n_inerT_file, file=inerT_file, status='unknown', &
               &    position='append')
            end if
            write(n_inerT_file,'(1P,2X,ES20.12,21ES12.4)') &
            &    time*tScale, ( real(vals_on_rank0_1d(ilm)),ilm=1,n_lm_vals )
            if ( l_save_out ) close(n_inerT_file)
         end if

      end if

   end subroutine write_rot
!-----------------------------------------------------------------------
   subroutine get_viscous_torque(viscous_torque,z10,dz10,r,dLrho,nu)
      !
      !  Purpose of this subroutine is to calculate the viscous torque
      !  on mantle or inner core respectively.
      !
      !  .. math::
      !     \Gamma_\nu=4\sqrt{\pi/3}\nu r\left[ \frac{\partial z_{10}}{\partial r}
      !     -(\frac{2}{r}+\beta)z_{10} \right]
      !  ..
      !

      !-- Input:
      real(cp), intent(in) :: z10,dz10    ! z10 coefficient and its radial deriv.
      real(cp), intent(in) :: r           ! radius (ICB or CMB)
      real(cp), intent(in) :: dLrho       ! dln(rho)/dr
      real(cp), intent(in) :: nu          ! viscosity

      !-- Output:
      real(cp), intent(out) :: viscous_torque

      viscous_torque=-four*sqrt(third*pi)*nu*r*( (two+dLrho*r)*real(z10)- &
      &               r*real(dz10) )

   end subroutine get_viscous_torque
!-----------------------------------------------------------------------
   subroutine get_lorentz_torque(lorentz_torque,br,bp,nR)
      !
      !  Purpose of this subroutine is to calculate the Lorentz torque
      !  on mantle or inner core respectively.
      !
      !  .. note:: ``lorentz_torque`` must be set to zero before loop over
      !            theta blocks is started.
      !
      !  .. warning:: subroutine returns ``-lorentz_torque`` if used at CMB
      !               to calculate torque on mantle because if the inward
      !               surface normal vector.
      !
      !  The Prandtl number is always the Prandtl number of the outer
      !  core. This comes in via scaling of the magnetic field.
      !  Theta alternates between northern and southern hemisphere in
      !  ``br`` and ``bp`` but not in gauss. This has to be cared for, and we
      !  use: ``gauss(latitude)=gauss(-latitude)`` here.
      !

      !-- Input variables:
      real(cp), intent(in) :: br(*)    ! array containing :math:`r^2 B_r`
      real(cp), intent(in) :: bp(*)    ! array containing :math:`r\sin\theta B_\phi`
      integer,  intent(in) :: nR         ! radial level

      real(cp), intent(inout) :: lorentz_torque ! Lorentz torque


      !-- Local variables:
      integer :: nTheta,nPhi,nThetaNHS,nelem
      real(cp) :: fac,b0r

      ! to avoid rounding errors for different theta blocking, we do not
      ! calculate sub sums with lorentz_torque_local, but keep on adding
      ! the contributions to the total lorentz_torque given as argument.

      lorentz_torque=0.0_cp

      fac=LFfac*two*pi/real(n_phi_max,cp) ! 2 pi/n_phi_max

      !$omp parallel do default(shared) &
      !$omp& private(nTheta, nPhi, nThetaNHS, b0r, nelem) &
      !$omp& reduction(+: lorentz_torque)
      do nPhi=1,n_phi_max
         do nTheta=1,n_theta_max
            nelem = radlatlon2spat(nTheta,nPhi,nR)

            nThetaNHS=(nTheta+1)/2 ! northern hemisphere=odd n_theta
            if ( lGrenoble ) then
               if ( r(nR) == r_icb ) then
                  b0r=two*BIC*r_icb**2*cosTheta(nTheta)
               else if ( r(nR) == r_cmb ) then
                  b0r=two*BIC*r_icb**2*cosTheta(nTheta)*(r_icb/r_cmb)
               end if
            else
               b0r=0.0_cp
            end if

            lorentz_torque=lorentz_torque + fac*gauss(nThetaNHS)* &
            &              (br(nelem)-b0r)*bp(nelem)
         end do
      end do
      !$omp end parallel do

   end subroutine get_lorentz_torque
!-----------------------------------------------------------------------
   subroutine get_angular_moment(z10,z11,omega_ic,omega_ma,angular_moment_oc, &
              &                  angular_moment_ic,angular_moment_ma)
      !
      !    Calculates angular momentum of outer core, inner core and
      !    mantle. For outer core we need ``z(l=1|m=0,1|r)``, for
      !    inner core and mantle the respective rotation rates are needed.
      !

      !-- Input of scalar fields:
      complex(cp), intent(in) :: z10(n_r_max),z11(n_r_max)
      real(cp),    intent(in) :: omega_ic,omega_ma

      !-- output:
      real(cp), intent(out) :: angular_moment_oc(:)
      real(cp), intent(out) :: angular_moment_ic(:)
      real(cp), intent(out) :: angular_moment_ma(:)

      !-- local variables:
      integer :: n_r,n
      integer :: l1m1
      real(cp) :: f(n_r_max,3)
      real(cp) :: r_E_2             ! r**2
      real(cp) :: fac

      !----- Construct radial function:
      l1m1=lo_map%lm2(1,1)
      do n_r=1,n_r_max
         r_E_2=r(n_r)*r(n_r)
         if ( l1m1 > 0 ) then
            f(n_r,1)=r_E_2* real(z11(n_r))
            f(n_r,2)=r_E_2*aimag(z11(n_r))
         else
            f(n_r,1)=0.0_cp
            f(n_r,2)=0.0_cp
         end if
         f(n_r,3)=r_E_2*real(z10(n_r))
      end do

      !----- Perform radial integral:
      do n=1,3
         angular_moment_oc(n)=rInt_R(f(:,n),r,rscheme_oc)
      end do

      !----- Apply normalisation factors of chebs and other factors
      !      plus the sign correction for y-component:
      fac=8.0_cp*third*pi
      angular_moment_oc(1)= two*fac*y11_norm * angular_moment_oc(1)
      angular_moment_oc(2)=-two*fac*y11_norm * angular_moment_oc(2)
      angular_moment_oc(3)=     fac*y10_norm * angular_moment_oc(3)

      !----- Now inner core and mantle:
      angular_moment_ic(1)=0.0_cp
      angular_moment_ic(2)=0.0_cp
      angular_moment_ic(3)=c_moi_ic*omega_ic
      angular_moment_ma(1)=0.0_cp
      angular_moment_ma(2)=0.0_cp
      angular_moment_ma(3)=c_moi_ma*omega_ma

   end subroutine get_angular_moment
!-----------------------------------------------------------------------
   subroutine get_angular_moment_Rloc(z10,z11,omega_ic,omega_ma,angular_moment_oc, &
              &                       angular_moment_ic,angular_moment_ma)
      !
      !    Calculates angular momentum of outer core, inner core and
      !    mantle. For outer core we need ``z(l=1|m=0,1|r)``, for
      !    inner core and mantle the respective rotation rates are needed.
      !    This is the version that takes r-distributed arrays as input arrays.
      !

      !-- Input of scalar fields:
      complex(cp), intent(in) :: z10(nRstart:nRstop),z11(nRstart:nRstop)
      real(cp),    intent(in) :: omega_ic,omega_ma

      !-- output:
      real(cp), intent(out) :: angular_moment_oc(:)
      real(cp), intent(out) :: angular_moment_ic(:)
      real(cp), intent(out) :: angular_moment_ma(:)

      !-- local variables:
      integer :: n_r,n
      real(cp) :: f_Rloc(nRstart:nRstop,3), f(n_r_max)
      real(cp) :: r_E_2,fac

      !----- Construct radial function:
      do n_r=nRstart,nRstop
         r_E_2=r(n_r)*r(n_r)
         f_Rloc(n_r,1)=r_E_2* real(z11(n_r))
         f_Rloc(n_r,2)=r_E_2*aimag(z11(n_r))
         f_Rloc(n_r,3)=r_E_2*real(z10(n_r))
      end do

      !----- Perform radial integration: (right now: MPI_Allgather, could be
      ! made local if needed)
      do n=1,3
         call allgather_from_rloc(f_Rloc(:,n), f)
         angular_moment_oc(n)=rInt_R(f,r,rscheme_oc)
      end do

      !----- Apply normalisation factors of chebs and other factors
      !      plus the sign correction for y-component:
      fac=8.0_cp*third*pi
      angular_moment_oc(1)= two*fac*y11_norm * angular_moment_oc(1)
      angular_moment_oc(2)=-two*fac*y11_norm * angular_moment_oc(2)
      angular_moment_oc(3)=     fac*y10_norm * angular_moment_oc(3)

      !----- Now inner core and mantle:
      angular_moment_ic(1)=0.0_cp
      angular_moment_ic(2)=0.0_cp
      angular_moment_ic(3)=c_moi_ic*omega_ic
      angular_moment_ma(1)=0.0_cp
      angular_moment_ma(2)=0.0_cp
      angular_moment_ma(3)=c_moi_ma*omega_ma

   end subroutine get_angular_moment_Rloc
!-----------------------------------------------------------------------
end module outRot
