module outRot

   use parallel_mod
   use precision_mod
   use truncation, only: n_r_max, n_r_maxMag, minc, nrp, n_phi_max
   use radial_data, only: n_r_CMB, n_r_ICB
   use radial_functions, only: r_icb, r_cmb, r, drx, chebt_oc
   use physical_parameters, only: kbotv, ktopv
   use num_param, only: lScale, tScale, vScale
   use blocking, only: lo_map,st_map,lmStartB,lmStopB, lm2
   use logic, only: l_AM, l_save_out, l_iner, l_SRIC, l_rot_ic, &
                    l_SRMA, l_rot_ma, l_mag_LF, l_mag, l_drift
   use output_data, only: tag, angular_file, n_angular_file,    &
                          n_SRIC_file, n_rot_file, n_SRMA_file, &
                          SRMA_file, SRIC_file, rot_file
   use constants, only: c_moi_oc, c_moi_ma, c_moi_ic, pi, y11_norm, &
                    y10_norm, zero, two, third, four, half
   use LMLoop_data, only: llm,ulm,llmMag,ulmMag
   use integration, only: rInt, rInt_R
   use horizontal_data, only: cosTheta, gauss
   use special, only: BIC, lGrenoble

   implicit none

   private

   interface get_viscous_torque
      module procedure get_viscous_torque_real
      module procedure get_viscous_torque_complex
   end interface get_viscous_torque

   public :: write_rot, get_viscous_torque, get_angular_moment, &
             get_lorentz_torque

contains

   subroutine write_rot(time,dt,eKinIC,ekinMA,w,z,dz,b, &
                      & omega_ic,omega_ma,lorentz_torque_ic,lorentz_torque_ma)
    
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
      real(cp), parameter :: tolerance=1e-16
      real(cp) :: eKinOC
      integer :: n_r1,n_r2,n_r3,nR
      integer :: l1m0,l1m1
      real(cp) :: viscous_torque_ic,viscous_torque_ma
      real(cp) :: AMz,eKinAMz
      real(cp) :: angular_moment_oc(3)
      real(cp) :: angular_moment_ic(3)
      real(cp) :: angular_moment_ma(3)
      complex(cp) :: z10(n_r_max),z11(n_r_max)
      character(len=80) :: filename
    
      real(cp) :: powerLor,powerVis
      real(cp), save :: AMzLast=0.0_cp,eKinAMzLast=0.0_cp
    
      integer, pointer :: lm2(:,:)
      integer :: i,l,m,ilm,lm_vals(21),n_lm_vals
      complex(cp) :: zvals_on_rank0(8,3),bvals_on_rank0(8,3)
      complex(cp) :: vals_on_rank0_1d(21)
    
      integer :: sr_tag
#ifdef WITH_MPI
      integer :: status(MPI_STATUS_SIZE),ierr
#endif
      logical :: rank_has_l1m0,rank_has_l1m1
      logical :: DEBUG_OUTPUT=.false.

      ! some arbitrary tag for the send and recv
      sr_tag=12345
    
      lm2(0:,0:) => lo_map%lm2
      l1m0=lm2(1,0)
    
      if ( DEBUG_OUTPUT ) write(*,"(I3,A,3I6)") rank,":lmStartB,lmStopB,l1m0=",& 
                        lmStartB(rank+1),lmStopB(rank+1),l1m0
    
      if ( lmStartB(rank+1) <= l1m0 .and. lmStopB(rank+1) >= l1m0 ) then
         !if (rank /= 0) then
         !   PRINT*,"in s_write_rot, l1m0 not on rank 0"
         !   stop
         !end if
         !-- Calculating viscous torques:
         if ( l_rot_ic .and. kbotv == 2 ) then
            call get_viscous_torque(viscous_torque_ic, &
                 &                  z(l1m0,n_r_max),dz(l1m0,n_r_max),r_icb)
         else
            viscous_torque_ic=0.0_cp
         end if
         if ( l_rot_ma .and. ktopv == 2 ) then
            call get_viscous_torque(viscous_torque_ma,z(l1m0,1),dz(l1m0,1),r_cmb)
         else
            viscous_torque_ma=0.0_cp
         end if
         rank_has_l1m0=.true.
#ifdef WITH_MPI
         if ( rank /= 0 ) then
            ! send viscous_torque_ic and viscous_torque_ma to rank 0 for 
            ! output
            call MPI_Send(viscous_torque_ic,1,MPI_DEF_REAL,0, &
                 &        sr_tag,MPI_COMM_WORLD,ierr)
            call MPI_Send(viscous_torque_ma,1,MPI_DEF_REAL,0, &
                 &        sr_tag+1,MPI_COMM_WORLD,ierr)
         end if
#endif
      else
         rank_has_l1m0=.false.
      end if
    
      if ( rank == 0 ) then
#ifdef WITH_MPI
         if ( .not. rank_has_l1m0 ) then
            call MPI_Recv(viscous_torque_ic,1,MPI_DEF_REAL,MPI_ANY_SOURCE,&
                 &        sr_tag,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(viscous_torque_ma,1,MPI_DEF_REAL,MPI_ANY_SOURCE,&
                 &        sr_tag+1,MPI_COMM_WORLD,status,ierr)
         end if
#endif
         if ( l_SRIC ) then
            powerLor=lorentz_torque_ic*omega_IC
            powerVis=viscous_torque_ic*omega_IC
            open(n_SRIC_file, file=SRIC_file, status="unknown", position='append')
            write(n_SRIC_file,'(1p,2x,ES20.12,4ES17.6)')   &
                 time*tScale,omega_ic/tScale,              &
                 (powerLor+powerVis)*vScale*vScale/tScale, &
                 powerVis*vScale*vScale/tScale,            &
                 powerLor*vScale*vScale/tScale
            close(n_SRIC_file)
         end if
         if ( l_SRMA ) then
            powerLor=lorentz_torque_ma*omega_ma
            powerVis=viscous_torque_ma*omega_ma
            open(n_SRMA_file, file=SRMA_file, status="unknown", position='append')
            write(n_SRMA_file,'(1p,2x,ES20.12,4ES17.6)')   &
                 time*tScale, omega_ma/tScale,             &
                 (powerLor+powerVis)*vScale*vScale/tScale, &
                 powerVis*vScale*vScale/tScale,            &
                 powerLor*vScale*vScale/tScale
            close(n_SRMA_file)
         end if
      end if
    
      if ( l_drift ) then
         do i=1,4
            lm_vals(i)=lm2(i*minc,i*minc)
            lm_vals(4+i)=lm2(i*minc+1,i*minc)
         end do
         n_r1=int(third*(n_r_max-1))
         n_r2=int(two*third*(n_r_max-1))
         n_r3=n_r_max-1
         call sendvals_to_rank0(z,n_r1,lm_vals(1:8),zvals_on_rank0(:,1))
         call sendvals_to_rank0(z,n_r2,lm_vals(1:8),zvals_on_rank0(:,2))
         call sendvals_to_rank0(z,n_r3,lm_vals(1:8),zvals_on_rank0(:,3))
    
         if ( rank == 0 ) then
            filename='driftVD.'//tag
            open(n_SRIC_file, file=filename, status='unknown', position='append')
            write(n_SRIC_file,'(1P,2X,ES20.12,24ES12.4)') &
                 time, (zvals_on_rank0(ilm,1),ilm=1,4),   &
                 (zvals_on_rank0(ilm,2),ilm=1,4),         &
                 (zvals_on_rank0(ilm,3),ilm=1,4)
            close(n_SRIC_file)
            filename='driftVQ.'//tag
            open(n_SRIC_file, file=filename, status='unknown', position='append')
            write(n_SRIC_file,'(1P,2X,ES20.12,24ES12.4)') &
                 time, (zvals_on_rank0(ilm,1),ilm=5,8),   &
                 (zvals_on_rank0(ilm,2),ilm=5,8),         &
                 (zvals_on_rank0(ilm,3),ilm=5,8)
            close(n_SRIC_file)
         end if
         
         if ( l_mag .or. l_mag_LF ) then
            n_r1=n_r_CMB
            n_r2=n_r_ICB
            call sendvals_to_rank0(b,n_r1,lm_vals(1:8),bvals_on_rank0(:,1))
            call sendvals_to_rank0(b,n_r2,lm_vals(1:8),bvals_on_rank0(:,2))
    
            if ( rank == 0 ) then
               filename='driftBD.'//tag
               open(n_SRIC_file, file=filename, status='unknown', position='append')
               write(n_SRIC_file,'(1P,2X,ES20.12,16ES12.4)') &
                    time, (bvals_on_rank0(ilm,1),ilm=5,8),   &
                    (bvals_on_rank0(ilm,2),ilm=5,8)
               close(n_SRIC_file)
               filename='driftBQ.'//tag
               open(n_SRIC_file, file=filename, status='unknown', position='append')
               write(n_SRIC_file,'(1P,2X,ES20.12,16ES12.4)') &
                    time, (bvals_on_rank0(ilm,1),ilm=1,4),   &
                    (bvals_on_rank0(ilm,2),ilm=1,4)
               close(n_SRIC_file)
            end if
         end if ! l_mag
      end if
    
      if ( .not. l_SRIC .and. ( l_rot_ic .or. l_rot_ma ) ) then
         if ( rank == 0 ) then
            if ( l_save_out ) then
               open(n_rot_file, file=rot_file, status='unknown', position='append')
            end if
            write(n_rot_file,'(1P,2X,ES20.12,6ES14.6)') &
                 time*tScale, omega_ic/tScale,          &
                 lScale**2*vScale*lorentz_torque_ic,    &
                 lScale**2*vScale*viscous_torque_ic,    &
                 omega_ma/tScale,                       &
                 lScale**2*vScale*lorentz_torque_ma,    &
                 -lScale**2*vScale*viscous_torque_ma
            if ( l_save_out ) close(n_rot_file)
         end if
      end if
      
      if ( l_AM ) then
         rank_has_l1m0=.false.
         rank_has_l1m1=.false.
         l1m0=lo_map%lm2(1,0)
         l1m1=lo_map%lm2(1,1)
         if ( (lmStartB(rank+1) <= l1m0) .and. (l1m0 <= lmStopB(rank+1)) ) then
            do nR=1,n_r_max
               z10(nR)=z(l1m0,nR)
            end do
            rank_has_l1m0=.true.
#ifdef WITH_MPI
            if (rank /= 0) then
               call MPI_Send(z10,n_r_max,MPI_DEF_COMPLEX,0,sr_tag, & 
                             MPI_COMM_WORLD,ierr)
            end if
#endif
         end if
    
         if ( l1m1 > 0 ) then
            if ( (lmStartB(rank+1) <= l1m1) .and. (l1m1 <= lmStopB(rank+1)) ) then
               do nR=1,n_r_max
                  z11(nR)=z(l1m1,nR)
               end do
               rank_has_l1m1=.true.
#ifdef WITH_MPI
               if ( rank /= 0 ) then
                  call MPI_Send(z11,n_r_max,MPI_DEF_COMPLEX,0, &
                              & sr_tag+1,MPI_COMM_WORLD,ierr)
               end if
#endif
            end if
         else
            do nR=1,n_r_max
               z11(nR)=zero
            end do
         end if
         ! now we have z10 and z11 in the worst case on two different
         ! ranks, which are also different from rank 0
         if ( rank == 0 ) then
#ifdef WITH_MPI
            if ( .not. rank_has_l1m0 ) then
               call MPI_Recv(z10,n_r_max,MPI_DEF_COMPLEX,&
                    &        MPI_ANY_SOURCE,sr_tag,MPI_COMM_WORLD,status,ierr)
            end if
            if ( l1m1 > 0 ) then
               if ( .not. rank_has_l1m1 ) then
                  call MPI_Recv(z11,n_r_max,MPI_DEF_COMPLEX,&
                       &        MPI_ANY_SOURCE,sr_tag+1,MPI_COMM_WORLD,status,ierr)
               end if
            end if
#endif
    
            call get_angular_moment(z10,z11,omega_ic,omega_ma,angular_moment_oc, &
                                    angular_moment_ic,angular_moment_ma)
            if ( l_save_out ) then
               open(n_angular_file, file=angular_file, status='unknown', &
                    position='append')
            end if
            AMz=angular_moment_oc(3)+angular_moment_ic(3)+angular_moment_ma(3)
            if ( abs(AMz) < tolerance ) AMz=0.0_cp
            eKinAMz=half*(angular_moment_oc(3)**2/c_moi_oc + &
                          angular_moment_ic(3)**2/c_moi_ic + &
                          angular_moment_ma(3)**2/c_moi_ma )
            if ( abs(eKinAMz) < tolerance ) eKinAMz=0.0_cp
            eKinIC=half*angular_moment_ic(3)**2/c_moi_ic
            eKinOC=half*angular_moment_oc(3)**2/c_moi_oc
            eKinMA=half*angular_moment_ma(3)**2/c_moi_ma
            if ( AMzLast /= 0.0_cp ) then
               !write(*,"(A,4ES22.15)") "col9 = ",eKinAMz,eKinAMzLast, &
               !     &                  dt,(eKinAMz-eKinAMzLast)
               write(n_angular_file,'(1p,2x,ES20.12,5ES14.6,3ES20.12)', advance='no') &
                    & time*tScale, angular_moment_oc,                                 &
                    & angular_moment_ic(3), angular_moment_ma(3),                     &
                    & AMz,(AMz-AMzLast)/AMzLast/dt,eKinAMz
               if (eKinAMzLast /= 0.0_cp) then
                  write(n_angular_file,'(1ES20.12)', advance='no') &
                    & (eKinAMz-eKinAMzLast)/eKinAMzLast/dt
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
         ! l_iner can only be .TRUE. for minc=1
         n_lm_vals=0
         do l=1,6
            do m=1,l
               n_lm_vals = n_lm_vals + 1
               lm_vals(n_lm_vals)=lm2(l,m)
            end do
         end do
         n_r1=int(half*(n_r_max-1))
         call sendvals_to_rank0(w,n_r1,lm_vals(1:n_lm_vals),vals_on_rank0_1d)
    
         if ( rank == 0 ) then
            filename='inerP.'//tag
            open(n_SRIC_file, file=filename, status='unknown', position='append')
            write(n_SRIC_file,'(1P,2X,ES20.12,21ES12.4)') &
                 time, ( real(vals_on_rank0_1d(ilm)),ilm=1,n_lm_vals )
            close(n_SRIC_file)
         end if
    
         n_r1=int(half*(n_r_max-1))
         call sendvals_to_rank0(z,n_r1,lm_vals(1:n_lm_vals),vals_on_rank0_1d)
    
         if ( rank == 0 ) then
            filename='inerT.'//tag
            open(n_SRIC_file, file=filename, status='unknown', position='append')
            write(n_SRIC_file,'(1P,2X,ES20.12,21ES12.4)') &
                 time, ( real(vals_on_rank0_1d(ilm)),ilm=1,n_lm_vals ) 
            close(n_SRIC_file)
         end if
    
      end if

   end subroutine write_rot
!-----------------------------------------------------------------------
   subroutine get_viscous_torque_real(viscous_torque,z10,dz10,r)
      !
      !  Purpose of this subroutine is to calculate the viscous torque    
      !  on mantle or inner core respectively.                            
      !  NOTE: sign is wrong for torque on mantle!                        
      !

      !-- Input:
      real(cp), intent(in) :: z10,dz10    ! z10 coefficient and its radial deriv.
      real(cp), intent(in) :: r               ! radius (ICB or CMB)

      !-- Output:
      real(cp), intent(out) :: viscous_torque

      viscous_torque=-four*sqrt(third*pi)*r *( two*real(z10) - r*real(dz10) )

   end subroutine get_viscous_torque_real
!-----------------------------------------------------------------------
   subroutine get_viscous_torque_complex(viscous_torque,z10,dz10,r)
      !
      !  Purpose of this subroutine is to calculate the viscous torque    
      !  on mantle or inner core respectively.                            
      !  NOTE: sign is wrong for torque on mantle!                        
      !

      !-- Input:
      complex(cp), intent(in) :: z10,dz10    ! z10 coefficient and its radial deriv.
      real(cp),    intent(in) :: r               ! radius (ICB or CMB)

      !-- Output:
      real(cp), intent(out) :: viscous_torque

      viscous_torque=-four*sqrt(third*pi)*r *( two*real(z10) - r*real(dz10) )

   end subroutine get_viscous_torque_complex
!-----------------------------------------------------------------------
   subroutine get_lorentz_torque(lorentz_torque,nThetaStart, &
                                 sizeThetaB,br,bp,nR)
      !
      !  Purpose of this subroutine is to calculate the lorentz torque    
      !  on mantle or inner core respectively.                            
      !  Blocking in theta can be used to increased performance.          
      !  If no blocking required set n_theta_block=n_theta_max,           
      !  where n_theta_max is the absolut number of thetas used.          
      !
      !  .. note:: Lorentz_torque must be set to zero before loop over        
      !            theta blocks is started.                                   
      !
      !  .. warning:: subroutine returns -lorentz_torque if used at CMB       
      !               to calculate torque on mantle because if the inward        
      !               surface normal vector.
      !
      !  The Prandtl number is always the Prandtl number of the outer     
      !  core. This comes in via scaling of the magnetic field.           
      !  Theta alternates between northern and southern hemisphere in     
      !  br and bp but not in gauss. This has to be cared for, and we     
      !  use: gauss(latitude)=gauss(-latitude) here.                      
      !

      !-- Input variables:
      integer,  intent(in) :: nThetaStart    ! first number of theta in block
      integer,  intent(in) :: sizeThetaB     ! size of theta bloching
      real(cp), intent(in) :: br(nrp,*)      ! array containing
      real(cp), intent(in) :: bp(nrp,*)      ! array containing
      integer,  intent(in) :: nR

      real(cp), intent(inout) :: lorentz_torque ! lorentz_torque for theta(1:n_theta)


      !-- local variables:
      integer :: nTheta,nPhi,nThetaNHS
      integer :: nThetaB
      real(cp) :: fac,b0r

      ! to avoid rounding errors for different theta blocking, we do not
      ! calculate sub sums with lorentz_torque_local, but keep on adding
      ! the contributions to the total lorentz_torque given as argument.

      if ( nThetaStart == 1 ) then
         lorentz_torque=0.0_cp
      end if

      !lorentz_torque_local=0.0_cp
      fac=two*pi/real(n_phi_max,cp) ! 2 pi/n_phi_max

      nTheta=nThetaStart-1
#ifdef WITH_SHTNS
      !$OMP PARALLEL DO default(none) &
      !$OMP& private(nThetaB, nTheta, nPhi, nThetaNHS, b0r) &
      !$OMP& shared(n_phi_max, sizeThetaB, r_icb, r, nR) &
      !$OMP& shared(lGrenoble, nThetaStart, BIC, cosTheta, r_cmb) &
      !$OMP& shared(fac, gauss) &
      !$OMP& reduction(+: lorentz_torque)
#endif
      do nThetaB=1,sizeThetaB
         nTheta=nThetaStart+nThetaB-1
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

         do nPhi=1,n_phi_max
            !lorentz_torque_local=lorentz_torque_local + &
            !                          gauss(nThetaNHS) * &
            !       (br(nPhi,nThetaB)-b0r)*bp(nPhi,nThetaB)
            lorentz_torque=lorentz_torque + fac * gauss(nThetaNHS) * &
                   (br(nPhi,nThetaB)-b0r)*bp(nPhi,nThetaB)
         end do
         !lorentz_torque_local = lorentz_torque_local + gauss(nThetaNHS)*phisum
      end do
#ifdef WITH_SHTNS
      !$OMP END PARALLEL DO
#endif

      !-- normalisation of phi-integration and division by Pm:
      !lorentz_torque=lorentz_torque+fac*lorentz_torque_local
              
   end subroutine get_lorentz_torque
!-----------------------------------------------------------------------

   subroutine get_angular_moment(z10,z11,omega_ic,omega_ma,angular_moment_oc, &
                                 angular_moment_ic,angular_moment_ma)
      !
      !    Calculates angular momentum of outer core, inner core and      
      !    mantle. For outer core we need z(l=1|m=0,1|r), for             
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
      integer :: l1m0,l1m1
      real(cp) :: f(n_r_max,3)
      real(cp) :: r_E_2             ! r**2
      real(cp) :: fac
    
      !----- Construct radial function:
      l1m0=lm2(1,0)
      l1m1=lm2(1,1)
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
         angular_moment_oc(n)=rInt_R(f(1,n),n_r_max,n_r_max,drx,chebt_oc)
      end do
    
      !----- Apply normalisation factors of chebs and other factors
      !      plus the sign correction for y-component:
      fac=8.0_cp*third*pi
      angular_moment_oc(1)= two*fac*y11_norm * angular_moment_oc(1)
      angular_moment_oc(2)=-two*fac*y11_norm * angular_moment_oc(2)
      angular_moment_oc(3)=      fac*y10_norm * angular_moment_oc(3)
    
      !----- Now inner core and mantle:
      angular_moment_ic(1)=0.0_cp
      angular_moment_ic(2)=0.0_cp
      angular_moment_ic(3)=c_moi_ic*omega_ic
      angular_moment_ma(1)=0.0_cp
      angular_moment_ma(2)=0.0_cp
      angular_moment_ma(3)=c_moi_ma*omega_ma

   end subroutine get_angular_moment
!-----------------------------------------------------------------------
   subroutine sendvals_to_rank0(field,n_r,lm_vals,vals_on_rank0)

      !-- Input variables:
      complex(cp), intent(in) :: field(llm:ulm,n_r_max)
      integer,     intent(in) :: n_r
      integer,     intent(in) :: lm_vals(:)

      !-- Output variables:
      complex(cp), intent(out) :: vals_on_rank0(:)

      !-- Local variables:
      integer :: ilm,lm,tag,n_lm_vals
#ifdef WITH_MPI
      integer :: ierr,status(MPI_STATUS_SIZE)
#endif
    
      n_lm_vals=size(lm_vals)
      if ( size(vals_on_rank0) < n_lm_vals ) then
         write(*,"(2(A,I4))") "write_rot: length of vals_on_rank0=",size(vals_on_rank0),&
              &" must be >= size(lm_vals)=",n_lm_vals
#ifdef WITH_MPI
         call mpi_abort(MPI_COMM_WORLD,43,ierr)
#endif
      end if

      do ilm=1,n_lm_vals
         lm=lm_vals(ilm)
         if ( lmStartB(1) <= lm .and. lm <= lmStopB(1) ) then
            ! the value is already on rank 0
            if (rank == 0) vals_on_rank0(ilm)=field(lm,n_r)
         else
            tag=876+ilm
            ! on which process is the lm value?
#ifdef WITH_MPI
            if (lmStartB(rank+1) <= lm .and. lm <= lmStopB(rank+1)) then
               call MPI_Send(field(lm,n_r),1,MPI_DEF_COMPLEX,&
                    & 0,tag,MPI_COMM_WORLD,ierr)
            end if
            if (rank == 0) then
               call MPI_Recv(vals_on_rank0(ilm),1,MPI_DEF_COMPLEX,&
                    & MPI_ANY_SOURCE,tag,MPI_COMM_WORLD,status,ierr)
            end if
#endif
         end if
      end do
   end subroutine sendvals_to_rank0
!-----------------------------------------------------------------------
end module outRot
