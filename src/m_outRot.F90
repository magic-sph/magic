!$Id$
!***********************************************************************
module outRot

   use mpi
   use truncation, only: n_r_max, n_r_maxMag, minc
   use radial_data, only: n_r_CMB, n_r_ICB
   use radial_functions, only: r_icb, r_cmb, r, drx, i_costf_init, &
                               d_costf_init
   use physical_parameters, only: kbotv, ktopv
   use num_param, only: lScale, tScale, vScale
   use blocking, only: lo_map,st_map,lmStartB,lmStopB, lm2
   use logic, only: l_AM, l_save_out, l_iner, l_SRIC, l_rot_ic, &
                    l_SRMA, l_rot_ma, l_mag_LF, l_mag, l_drift
   use output_data, only: tag, angular_file, n_angular_file,    &
                          n_SRIC_file, n_rot_file, n_SRMA_file, &
                          SRMA_file, SRIC_file, rot_file
   use const, only: c_moi_oc, c_moi_ma, c_moi_ic, pi, y11_norm, &
                    y10_norm
   use parallel_mod, only: rank
   use LMLoop_data, only: llm,ulm,llmMag,ulmMag
   use integration, only: rInt, rInt_R

   implicit none

   private

   interface get_viscous_torque
      module procedure get_viscous_torque_real
      module procedure get_viscous_torque_complex
   end interface get_viscous_torque

   public :: write_rot, get_viscous_torque, get_angular_moment

contains

   subroutine write_rot(time,dt,eKinIC,ekinMA,w,z,dz,b, &
                      & omega_ic,omega_ma,lorentz_torque_ic,lorentz_torque_ma)
      !***********************************************************************
    
      !-- Input of variables:
      real(kind=8),    intent(in) :: omega_ic,omega_ma
      real(kind=8),    intent(in) :: lorentz_torque_ma,lorentz_torque_ic
      real(kind=8),    intent(in) :: time,dt
      complex(kind=8), intent(in) :: w(llm:ulm,n_r_max)
      complex(kind=8), intent(in) :: z(llm:ulm,n_r_max)
      complex(kind=8), intent(in) :: dz(llm:ulm,n_r_max)
      complex(kind=8), intent(in) :: b(llmMag:ulmMag,n_r_maxMag)
    
      !-- Output into rot_file
      real(kind=8), intent(out) :: eKinIC,eKinMA
    
      !-- Local variables:
      real(kind=8), parameter :: tolerance=1e-16
      real(kind=8) :: eKinOC
      integer :: n_r1,n_r2,n_r3,nR
      integer :: l1m0,l1m1
      real(kind=8) :: viscous_torque_ic,viscous_torque_ma
      real(kind=8) :: AMz,eKinAMz
      real(kind=8) :: angular_moment_oc(3)
      real(kind=8) :: angular_moment_ic(3)
      real(kind=8) :: angular_moment_ma(3)
      complex(kind=8) :: z10(n_r_max),z11(n_r_max)
      character(len=80) :: filename
    
      real(kind=8) :: powerLor,powerVis
    
      real(kind=8) :: AMzLast,eKinAMzLast
      SAVE AMzLast,eKinAMzLast
    
      integer, pointer :: lm2(:,:)
      integer :: i,l,m,ilm,lm_vals(21),n_lm_vals
      complex(kind=8) :: zvals_on_rank0(8,3),bvals_on_rank0(8,3)
      complex(kind=8) :: vals_on_rank0_1d(21)
    
      integer :: sr_tag,status(MPI_STATUS_SIZE),ierr
      logical :: rank_has_l1m0,rank_has_l1m1
      logical :: DEBUG_OUTPUT=.false.

      ! some arbitrary tag for the send and recv
      sr_tag=12345
    
      lm2(0:,0:) => lo_map%lm2
      l1m0=lm2(1,0)
    
      if (DEBUG_OUTPUT) write(*,"(I3,A,3I6)") rank,":lmStartB,lmStopB,l1m0=",lmStartB(rank+1),lmStopB(rank+1),l1m0
    
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
            viscous_torque_ic=0.d0
         end if
         if ( l_rot_ma .and. ktopv == 2 ) then
            call get_viscous_torque(viscous_torque_ma,z(l1m0,1),dz(l1m0,1),r_cmb)
         else
            viscous_torque_ma=0.d0
         end if
         rank_has_l1m0=.true.
         if ( rank /= 0 ) then
            ! send viscous_torque_ic and viscous_torque_ma to rank 0 for 
            ! output
            call MPI_Send(viscous_torque_ic,1,MPI_DOUBLE_PRECISION,0,&
                 &sr_tag,MPI_COMM_WORLD,ierr)
            call MPI_Send(viscous_torque_ma,1,MPI_DOUBLE_PRECISION,0,&
                 &sr_tag+1,MPI_COMM_WORLD,ierr)
         end if
      else
         rank_has_l1m0=.false.
      end if
    
      if ( rank == 0 ) then
         if ( .not. rank_has_l1m0 ) then
            call MPI_Recv(viscous_torque_ic,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,&
                 &sr_tag,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(viscous_torque_ma,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,&
                 &sr_tag+1,MPI_COMM_WORLD,status,ierr)
         end if
         if ( l_SRIC ) then
            powerLor=lorentz_torque_ic*omega_IC
            powerVis=viscous_torque_ic*omega_IC
            open(n_SRIC_file, file=SRIC_file, status="unknown", position='APPEND')
            write(n_SRIC_file,'(1p,2x,d20.12,4d17.6)')     &
                 time*tScale,omega_ic/tScale,              &
                 (powerLor+powerVis)*vScale*vScale/tScale, &
                 powerVis*vScale*vScale/tScale,            &
                 powerLor*vScale*vScale/tScale
            close(n_SRIC_file)
         end if
         if ( l_SRMA ) then
            powerLor=lorentz_torque_ma*omega_ma
            powerVis=viscous_torque_ma*omega_ma
            open(n_SRMA_file, file=SRMA_file, status="unknown", position='APPEND')
            write(n_SRMA_file,'(1p,2x,d20.12,4d17.6)')     &
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
         n_r1=int(1.D0/3.D0*(n_r_max-1))
         n_r2=int(2.D0/3.D0*(n_r_max-1))
         n_r3=n_r_max-1
         call sendvals_to_rank0(z,n_r1,lm_vals(1:8),zvals_on_rank0(:,1))
         call sendvals_to_rank0(z,n_r2,lm_vals(1:8),zvals_on_rank0(:,2))
         call sendvals_to_rank0(z,n_r3,lm_vals(1:8),zvals_on_rank0(:,3))
    
         if ( rank == 0 ) then
            filename='driftVD.'//tag
            open(n_SRIC_file, file=filename, status='UNKNOWN', position='APPEND')
            write(n_SRIC_file,'(1P,2X,D20.12,24D12.4)') &
                 time, (zvals_on_rank0(ilm,1),ilm=1,4), &
                 (zvals_on_rank0(ilm,2),ilm=1,4),       &
                 (zvals_on_rank0(ilm,3),ilm=1,4)
            close(n_SRIC_file)
            filename='driftVQ.'//tag
            open(n_SRIC_file, file=filename, status='UNKNOWN', position='APPEND')
            write(n_SRIC_file,'(1P,2X,D20.12,24D12.4)') &
                 time, (zvals_on_rank0(ilm,1),ilm=5,8), &
                 (zvals_on_rank0(ilm,2),ilm=5,8),       &
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
               open(n_SRIC_file, file=filename, status='UNKNOWN', position='APPEND')
               write(n_SRIC_file,'(1P,2X,D20.12,16D12.4)') &
                    time, (bvals_on_rank0(ilm,1),ilm=5,8), &
                    (bvals_on_rank0(ilm,2),ilm=5,8)
               close(n_SRIC_file)
               filename='driftBQ.'//tag
               open(n_SRIC_file, file=filename, status='UNKNOWN', position='APPEND')
               write(n_SRIC_file,'(1P,2X,D20.12,16D12.4)') &
                    time, (bvals_on_rank0(ilm,1),ilm=1,4), &
                    (bvals_on_rank0(ilm,2),ilm=1,4)
               close(n_SRIC_file)
            end if
         end if ! l_mag
      end if
    
      if ( .not. l_SRIC .and. ( l_rot_ic .or. l_rot_ma ) ) then
         if ( rank == 0 ) then
            if ( l_save_out ) then
               open(n_rot_file, file=rot_file, status='UNKNOWN', position='APPEND')
            end if
            write(n_rot_file,'(1P,2X,D20.12,6D14.6)') &
                 time*tScale, omega_ic/tScale,        &
                 lScale**2*vScale*lorentz_torque_ic,  &
                 lScale**2*vScale*viscous_torque_ic,  &
                 omega_ma/tScale,                     &
                 lScale**2*vScale*lorentz_torque_ma,  &
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
            if (rank /= 0) then
               call MPI_Send(z10,n_r_max,MPI_DOUBLE_complex,0,sr_tag,MPI_COMM_WORLD,ierr)
            end if
         end if
    
         if ( l1m1 > 0 ) then
            if ( (lmStartB(rank+1) <= l1m1) .and. (l1m1 <= lmStopB(rank+1)) ) then
               do nR=1,n_r_max
                  z11(nR)=z(l1m1,nR)
               end do
               rank_has_l1m1=.true.
               if ( rank /= 0 ) then
                  call MPI_Send(z11,n_r_max,MPI_DOUBLE_complex,0, &
                              & sr_tag+1,MPI_COMM_WORLD,ierr)
               end if
            end if
         else
            do nR=1,n_r_max
               z11(nR)=cmplx(0d0,0d0,kind=kind(0d0))
            end do
         end if
         ! now we have z10 and z11 in the worst case on two different
         ! ranks, which are also different from rank 0
         if ( rank == 0 ) then
            if ( .not. rank_has_l1m0 ) then
               call MPI_Recv(z10,n_r_max,MPI_DOUBLE_complex,&
                    & MPI_ANY_SOURCE,sr_tag,MPI_COMM_WORLD,status,ierr)
            end if
            if ( l1m1 > 0 ) then
               if ( .not. rank_has_l1m1 ) then
                  call MPI_Recv(z11,n_r_max,MPI_DOUBLE_complex,&
                       & MPI_ANY_SOURCE,sr_tag+1,MPI_COMM_WORLD,status,ierr)
               end if
            end if
    
            call get_angular_moment(z10,z11,omega_ic,omega_ma, &
                 angular_moment_oc, angular_moment_ic,angular_moment_ma)
            if ( l_save_out ) then
               open(n_angular_file,file=angular_file,status='UNKNOWN',position='APPEND')
            end if
            AMz=angular_moment_oc(3)+angular_moment_ic(3)+angular_moment_ma(3)
            if ( abs(AMz) < tolerance ) AMz=0.0D0
            eKinAMz=0.5d0*(angular_moment_oc(3)**2/c_moi_oc + &
                           angular_moment_ic(3)**2/c_moi_ic + &
                           angular_moment_ma(3)**2/c_moi_ma )
            if ( abs(eKinAMz) < tolerance ) eKinAMz=0.0D0
            eKinIC=0.5d0*angular_moment_ic(3)**2/c_moi_ic
            eKinOC=0.5d0*angular_moment_oc(3)**2/c_moi_oc
            eKinMA=0.5d0*angular_moment_ma(3)**2/c_moi_ma
            if ( AMzLast /= 0.0D0 ) then
               !write(*,"(A,4ES22.15)") "col9 = ",eKinAMz,eKinAMzLast,dt,(eKinAMz-eKinAMzLast)
               write(n_angular_file,'(1p,2x,d20.12,5d14.6,3d20.12)', advance='no') &
                    time*tScale, angular_moment_oc,               &
                    & angular_moment_ic(3), angular_moment_ma(3), &
                    & AMz,(AMz-AMzLast)/AMzLast/dt,eKinAMz
               if (eKinAMzLast /= 0.0d0) then
                  write(n_angular_file,'(1d20.12)', advance='no') &
                    & (eKinAMz-eKinAMzLast)/eKinAMzLast/dt
               else
                  write(n_angular_file,'(1d20.12)', advance='no') 0.0
               end if
               write(n_angular_file,'(3d20.12)') eKinIC,eKinOC,eKinMA
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
         n_r1=int(1.D0/2.D0*(n_r_max-1))
         call sendvals_to_rank0(w,n_r1,lm_vals(1:n_lm_vals),vals_on_rank0_1d)
    
         if ( rank == 0 ) then
            filename='inerP.'//tag
            open(n_SRIC_file, file=filename, status='UNKNOWN', position='APPEND')
            write(n_SRIC_file,'(1P,2X,D20.12,21D12.4)') &
                 time, ( real(vals_on_rank0_1d(ilm)),ilm=1,n_lm_vals )
            close(n_SRIC_file)
         end if
    
         n_r1=int(1.D0/2.D0*(n_r_max-1))
         call sendvals_to_rank0(z,n_r1,lm_vals(1:n_lm_vals),vals_on_rank0_1d)
    
         if ( rank == 0 ) then
            filename='inerT.'//tag
            open(n_SRIC_file, file=filename, status='UNKNOWN', position='APPEND')
            write(n_SRIC_file,'(1P,2X,D20.12,21D12.4)') &
                 time, ( real(vals_on_rank0_1d(ilm)),ilm=1,n_lm_vals ) 
            close(n_SRIC_file)
         end if
    
      end if

   end subroutine write_rot
!-----------------------------------------------------------------------
   subroutine get_viscous_torque_real(viscous_torque,z10,dz10,r)
      !  +-------------------------------------------------------------------+
      !  |                                                                   |
      !  |  Purpose of this subroutine is to calculate the viscous torque    |
      !  |  on mantle or inner core respectively.                            |
      !  |  NOTE: sign is wrong for torque on mantle!                        |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+

      !-- Input:
      real(kind=8), intent(in) :: z10,dz10    ! z10 coefficient and its radial deriv.
      real(kind=8), intent(in) :: r               ! radius (ICB or CMB)

      !-- Output:
      real(kind=8), intent(out) :: viscous_torque

      !-- Local:
      real(kind=8) :: pi

      viscous_torque=-4.D0*dsqrt(pi/3.D0)*r *( 2.D0*real(z10) - r*real(dz10) )

   end subroutine get_viscous_torque_real
!-----------------------------------------------------------------------
   subroutine get_viscous_torque_complex(viscous_torque,z10,dz10,r)
      !  +-------------------------------------------------------------------+
      !  |                                                                   |
      !  |  Purpose of this subroutine is to calculate the viscous torque    |
      !  |  on mantle or inner core respectively.                            |
      !  |  NOTE: sign is wrong for torque on mantle!                        |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+

      !-- Input:
      complex(kind=8), intent(in) :: z10,dz10    ! z10 coefficient and its radial deriv.
      real(kind=8),    intent(in) :: r               ! radius (ICB or CMB)

      !-- Output:
      real(kind=8), intent(out) :: viscous_torque

      !-- Local:
      real(kind=8) :: pi

      viscous_torque=-4.D0*dsqrt(pi/3.D0)*r *( 2.D0*real(z10) - r*real(dz10) )

   end subroutine get_viscous_torque_complex
!-----------------------------------------------------------------------
   subroutine get_angular_moment(z10,z11,omega_ic,omega_ma,angular_moment_oc, &
                                 angular_moment_ic,angular_moment_ma)
      !  +-------------------------------------------------------------------+
      !  |                                                                   |
      !  |    Calculates angular momentum of outer core, inner core and      |
      !  |    mantle. For outer core we need z(l=1|m=0,1|r), for             |
      !  |    inner core and mantle the respective rotation rates are needed.|
      !  |                                                                   |
      !  +-------------------------------------------------------------------+
    
      !-- Input of scalar fields:
      complex(kind=8), intent(in) :: z10(n_r_max),z11(n_r_max)
      real(kind=8),    intent(in) :: omega_ic,omega_ma
    
      !-- output:
      real(kind=8), intent(out) :: angular_moment_oc(:)
      real(kind=8), intent(out) :: angular_moment_ic(:)
      real(kind=8), intent(out) :: angular_moment_ma(:)
    
      !-- local variables:
      integer :: n_r,n
      integer :: l1m0,l1m1
      real(kind=8) :: f(n_r_max,3)
      real(kind=8) :: r_E_2             ! r**2
      real(kind=8) :: fac
    
      !----- Construct radial function:
      l1m0=lm2(1,0)
      l1m1=lm2(1,1)
      do n_r=1,n_r_max
         r_E_2=r(n_r)*r(n_r)
         if ( l1m1 > 0 ) then
            f(n_r,1)=r_E_2* real(z11(n_r))
            f(n_r,2)=r_E_2*aimag(z11(n_r))
         else
            f(n_r,1)=0.D0
            f(n_r,2)=0.D0
         end if
         f(n_r,3)=r_E_2*real(z10(n_r))
      end do
    
      !----- Perform radial integral:
      do n=1,3
         angular_moment_oc(n)=rInt_R(f(1,n),n_r_max,n_r_max,drx, &
                                     i_costf_init,d_costf_init)
      end do
    
      !----- Apply normalisation factors of chebs and other factors
      !      plus the sign correction for y-component:
      fac=8.d0/3.d0*pi
      angular_moment_oc(1)= 2.d0*fac*y11_norm * angular_moment_oc(1)
      angular_moment_oc(2)=-2.d0*fac*y11_norm * angular_moment_oc(2)
      angular_moment_oc(3)=      fac*y10_norm * angular_moment_oc(3)
    
      !----- Now inner core and mantle:
      angular_moment_ic(1)=0.d0
      angular_moment_ic(2)=0.d0
      angular_moment_ic(3)=c_moi_ic*omega_ic
      angular_moment_ma(1)=0.d0
      angular_moment_ma(2)=0.d0
      angular_moment_ma(3)=c_moi_ma*omega_ma

   end subroutine get_angular_moment
!-----------------------------------------------------------------------
   subroutine sendvals_to_rank0(field,n_r,lm_vals,vals_on_rank0)

      !-- Input variables:
      complex(kind=8), intent(in) :: field(llm:ulm,n_r_max)
      integer,         intent(in) :: n_r
      integer,         intent(in) :: lm_vals(:)

      !-- Output variables:
      complex(kind=8), intent(out) :: vals_on_rank0(:)

      !-- Local variables:
      integer :: ilm,lm,ierr,status(MPI_STATUS_SIZE),tag,n_lm_vals
    
      n_lm_vals=size(lm_vals)
      if ( size(vals_on_rank0) < n_lm_vals ) then
         write(*,"(2(A,I4))") "write_rot: length of vals_on_rank0=",size(vals_on_rank0),&
              &" must be >= size(lm_vals)=",n_lm_vals
         call mpi_abort(MPI_COMM_WORLD,43,ierr)
      end if

      do ilm=1,n_lm_vals
         lm=lm_vals(ilm)
         if ( lmStartB(1) <= lm .and. lm <= lmStopB(1) ) then
            ! the value is already on rank 0
            if (rank == 0) vals_on_rank0(ilm)=field(lm,n_r)
         else
            tag=876+ilm
            ! on which process is the lm value?
            if (lmStartB(rank+1) <= lm .and. lm <= lmStopB(rank+1)) then
               call MPI_Send(field(lm,n_r),1,MPI_DOUBLE_complex,&
                    & 0,tag,MPI_COMM_WORLD,ierr)
            end if
            if (rank == 0) then
               call MPI_Recv(vals_on_rank0(ilm),1,MPI_DOUBLE_complex,&
                    & MPI_ANY_SOURCE,tag,MPI_COMM_WORLD,status,ierr)
            end if
         end if
      end do
   end subroutine sendvals_to_rank0
!-----------------------------------------------------------------------
end module outRot
