module updateZ_mod
   !
   ! This module handles the time advance of the toroidal potential z
   ! It contains the computation of the implicit terms and the linear
   ! solves.
   !

   use init_fields
   use omp_lib
   use precision_mod
   use communications, only: allgather_from_Rloc
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, lm_max, l_max
   use radial_data, only: n_r_cmb, n_r_icb, nRstart, nRstop
   use radial_functions, only: visc, or1, or2, rscheme_oc, dLvisc, beta, &
       &                       rho0, r_icb, r_cmb, r, beta, dbeta
   use physical_parameters, only: kbotv, ktopv, prec_angle, po, oek
   use num_param, only: AMstart, dct_counter, solve_counter
   use torsional_oscillations, only: ddzASL
   use blocking, only: lo_sub_map, lo_map, st_sub_map, llm, ulm, st_map
   use horizontal_data, only: hdif_V
   use logic, only: l_rot_ma, l_rot_ic, l_SRMA, l_SRIC, l_z10mat, l_precession, &
       &            l_correct_AMe, l_correct_AMz, l_update_v, l_TO,             &
       &            l_finite_diff, l_full_sphere, l_parallel_solve
   use RMS, only: DifTor2hInt
   use constants, only: c_lorentz_ma, c_lorentz_ic, c_dt_z10_ma, c_dt_z10_ic, &
       &                c_moi_ma, c_moi_ic, c_z10_omega_ma, c_z10_omega_ic,   &
       &                c_moi_oc, y10_norm, y11_norm, zero, one, two, four,   &
       &                pi, third
   use parallel_mod
   use outRot, only: get_angular_moment, get_angular_moment_Rloc
   use RMS_helpers, only: hInt2Tor
   use radial_der, only: get_ddr, get_ddr_ghost, bulk_to_ghost, exch_ghosts
   use fields, only: work_LMloc
   use useful, only: abortRun
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray, type_tscalar
   use special
   use dense_matrices
   use real_matrices
   use band_matrices
   use parallel_solvers

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
   type(type_tri_par), public :: zMat_FD, z10Mat_FD
   complex(cp), public, allocatable :: z_ghost(:,:), z10_ghost(:)

   integer :: maxThreads

   public :: updateZ, initialize_updateZ, finalize_updateZ, get_tor_rhs_imp,    &
   &         assemble_tor, finish_exp_tor, get_rot_rates, updateZ_FD,           &
   &         get_tor_rhs_imp_ghost, prepareZ_FD, fill_ghosts_Z, assemble_tor_Rloc

contains

   subroutine initialize_updateZ
      !
      ! This subroutine handles the memory allocation of the arrays involved
      ! in the time advance of the toroidal potential.
      !

      integer, pointer :: nLMBs2(:)
      integer :: ll, n_bands

      if ( .not. l_parallel_solve ) then
         nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

         if ( l_finite_diff ) then
            allocate( type_bandmat :: zMat(nLMBs2(1+rank)) )
            allocate( type_bandmat :: z10Mat )

            if ( ktopv /= 1 .and. kbotv /= 1 .and. rscheme_oc%order <= 2  .and. &
            &    rscheme_oc%order_boundary <= 2 ) then ! Rigid at both boundaries
               n_bands = rscheme_oc%order+1
            else
               n_bands = max(2*rscheme_oc%order_boundary+1,rscheme_oc%order+1)
            end if

            do ll=1,nLMBs2(1+rank)
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
            allocate( type_densemat :: zMat(nLMBs2(1+rank)) )
            allocate( type_densemat :: z10Mat )

            call z10Mat%initialize(n_r_max,n_r_max,l_pivot=.true.)
            do ll=1,nLMBs2(1+rank)
               call zMat(ll)%initialize(n_r_max,n_r_max,l_pivot=.true.)
            end do
         end if

#ifdef WITH_PRECOND_Z10
         allocate(z10Mat_fac(n_r_max))
         bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_REAL
#endif
#ifdef WITH_PRECOND_Z
         allocate(zMat_fac(n_r_max,nLMBs2(1+rank)))
         bytes_allocated = bytes_allocated+n_r_max*nLMBs2(1+rank)*SIZEOF_DEF_REAL
#endif

         allocate( Dif(llm:ulm) )
         bytes_allocated=bytes_allocated+(ulm-llm+1)*SIZEOF_DEF_COMPLEX

#ifdef WITHOMP
         maxThreads=omp_get_max_threads()
#else
         maxThreads=1
#endif

         allocate(rhs1(n_r_max,2*lo_sub_map%sizeLMB2max,0:maxThreads-1))
         bytes_allocated=bytes_allocated+n_r_max*maxThreads* &
         &               lo_sub_map%sizeLMB2max*SIZEOF_DEF_COMPLEX
      else
         allocate( Dif(lm_max) )
         bytes_allocated=bytes_allocated+lm_max*SIZEOF_DEF_COMPLEX
         allocate(z_ghost(lm_max,nRstart-1:nRstop+1),z10_ghost(nRstart-1:nRstop+1))
         bytes_allocated=bytes_allocated+(lm_max+1)*(nRstop-nRstart+3)*SIZEOF_DEF_COMPLEX
         z_ghost(:,:)=zero
         z10_ghost(:)=zero

         !-- Create matrix
         call zMat_FD%initialize(1,n_r_max,0,l_max)

         !-- Create matrix
         call z10Mat_FD%initialize(1,n_r_max,1,1)
      end if
      allocate( lZmat(0:l_max) )
      bytes_allocated = bytes_allocated+(l_max+1)*SIZEOF_LOGICAL

      AMstart=0.0_cp

   end subroutine initialize_updateZ
!-------------------------------------------------------------------------------
   subroutine finalize_updateZ
      !
      ! Memory deallocation of arrays associated with time advance of Z
      !

      integer, pointer :: nLMBs2(:)
      integer :: ll

      if ( .not. l_parallel_solve ) then
         nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

         do ll=1,nLMBs2(1+rank)
            call zMat(ll)%finalize()
         end do
         call z10Mat%finalize()

#ifdef WITH_PRECOND_Z10
         deallocate( z10Mat_fac )
#endif
#ifdef WITH_PRECOND_Z
         deallocate( zMat_fac )
#endif
         deallocate( rhs1 )
      else
         deallocate( z_ghost, z10_ghost )
         call zMat_FD%finalize()
         call z10Mat_FD%finalize()
      end if
      deallocate(Dif,lZmat)

   end subroutine finalize_updateZ
!-------------------------------------------------------------------------------
   subroutine updateZ(time,timeNext,z,dz,dzdt,omega_ma,omega_ic,domega_ma_dt,  &
              &       domega_ic_dt,lorentz_torque_ma_dt, lorentz_torque_ic_dt, &
              &       tscheme,lRmsNext)
      !
      !  Updates the toroidal potential z and its radial derivative dz
      !

      !-- Input/output of scalar fields:
      class(type_tscheme), intent(in) :: tscheme
      complex(cp), intent(inout) :: z(llm:ulm,n_r_max)  ! Toroidal velocity potential z
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
      complex(cp), intent(out) :: dz(llm:ulm,n_r_max)   ! Radial derivative of z
      real(cp),    intent(out) :: omega_ma              ! Calculated OC rotation
      real(cp),    intent(out) :: omega_ic              ! Calculated IC rotation

      !-- local variables:
      integer :: l1,m1              ! degree and order
      integer :: lm1,lm,lmB         ! position of (l,m) in array
      integer :: lmStart_00         ! excluding l=0,m=0
      integer :: nLMB2
      integer :: nR                 ! counts radial grid points
      integer :: n_r_out            ! counts cheb modes
      complex(cp) :: rhs(n_r_max)   ! RHS of matrix multiplication
      real(cp) :: prec_fac
      real(cp) :: dom_ma, dom_ic, lo_ma, lo_ic
      integer :: l1m0          ! position of (l=1,m=0) and (l=1,m=1) in lm.
      logical :: l10
      integer :: nLMB

      integer, pointer :: nLMBs2(:),lm2l(:)
      integer, pointer :: sizeLMB2(:,:),lm2(:,:)
      integer, pointer :: lm22lm(:,:,:),lm22l(:,:,:),lm22m(:,:,:)

      integer :: nChunks,iChunk,lmB0,size_of_last_chunk,threadid

      if ( l_precession ) then
         prec_fac=sqrt(8.0_cp*pi*third)*po*oek*oek*sin(prec_angle)
      else
         prec_fac = 0.0_cp
      end if

      if ( .not. l_update_v ) return

      nLMBs2(1:n_procs) => lo_sub_map%nLMBs2
      sizeLMB2(1:,1:) => lo_sub_map%sizeLMB2
      lm22lm(1:,1:,1:) => lo_sub_map%lm22lm
      lm22l(1:,1:,1:) => lo_sub_map%lm22l
      lm22m(1:,1:,1:) => lo_sub_map%lm22m
      lm2(0:,0:) => lo_map%lm2
      lm2l(1:lm_max) => lo_map%lm2l

      nLMB = 1+rank
      lmStart_00 =max(2,llm)
      l1m0       =lm2(1,0)

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

      !-- Now assemble the right hand side and store it in work_LMloc
      call tscheme%set_imex_rhs(work_LMloc, dzdt)

      !$omp parallel default(shared)

      !$omp single
      call solve_counter%start_count()
      !$omp end single
      l10=.false.
      !$omp single
      do nLMB2=1,nLMBs2(nLMB)
         !$omp task default(shared) &
         !$omp firstprivate(nLMB2) &
         !$omp private(lmB,lm,lm1,l1,m1,n_r_out,nR) &
         !$omp private(nChunks,size_of_last_chunk,iChunk) &
         !$omp private(tOmega_ma1,tOmega_ma2,threadid) &
         !$omp private(tOmega_ic1,tOmega_ic2)
         nChunks = (sizeLMB2(nLMB2,nLMB)+chunksize-1)/chunksize
         size_of_last_chunk=chunksize+(sizeLMB2(nLMB2,nLMB)-nChunks*chunksize)

         ! This task treats one l given by l1
         l1=lm22l(1,nLMB2,nLMB)

         if ( l1 /= 0 ) then
            if ( .not. lZmat(l1) ) then
#ifdef WITH_PRECOND_Z
               call get_zMat(tscheme,l1,hdif_V(l1),zMat(nLMB2),zMat_fac(:,nLMB2))
#else
               call get_zMat(tscheme,l1,hdif_V(l1),zMat(nLMB2))
#endif
               lZmat(l1)=.true.
            end if
         end if

         do iChunk=1,nChunks
            !$omp task default(shared) &
            !$omp firstprivate(iChunk) &
            !$omp private(lmB0,lmB,lm,lm1,m1,nR,n_r_out,threadid) &
            !$omp private(tOmega_ma1,tOmega_ma2,tOmega_ic1,tOmega_ic2)
#ifdef WITHOMP
            threadid = omp_get_thread_num()
#else
            threadid = 0
#endif

            lmB0=(iChunk-1)*chunksize
            lmB=lmB0

            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               lm1=lm22lm(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)
               if ( lm1 == l1m0 ) l10= .true.

               if ( l_z10mat .and. lm1 == l1m0 ) then
                  !----- Special treatment of z10 component if ic or mantle
                  !      are allowed to rotate about z-axis (l_z10mat=.true.) and
                  !      we use no slip boundary condition (ktopv=2,kbotv=2):
                  !      Lorentz torque is the explicit part of this time integration
                  !      at the boundaries!
                  !      Note: no angular momentum correction necessary for this case !
                  if ( .not. lZ10mat ) then
#ifdef WITH_PRECOND_Z10
                     call get_z10Mat(tscheme,l1,hdif_V(l1),z10Mat,z10Mat_fac)
#else
                     call get_z10Mat(tscheme,l1,hdif_V(l1),z10Mat)
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
                     rhs(nR)=work_LMloc(lm1,nR)
                  end do

#ifdef WITH_PRECOND_Z10
                  rhs(:) = z10Mat_fac(:)*rhs(:)
#endif

                  call z10Mat%solve(rhs)

               else if ( l1 /= 0 ) then
                  lmB=lmB+1

                  rhs1(1,2*lmB-1,threadid)      =0.0_cp
                  rhs1(1,2*lmB,threadid)        =0.0_cp
                  rhs1(n_r_max,2*lmB-1,threadid)=0.0_cp
                  rhs1(n_r_max,2*lmB,threadid)  =0.0_cp

                  if (amp_RiIc /= 0.0_cp) then
                     if (l1 == (m_RiIc + RiSymmIc) .and. m1 == m_RiIc) then
                        rhs1(n_r_max,2*lmB-1,threadid)=amp_RiIc* &
                        &                              cos(omega_RiIc*time)
                        rhs1(n_r_max,2*lmB,threadid)  =amp_RiIc* &
                        &                              sin(omega_RiIc*time)
                     end if
                  end if

                  if (amp_RiMa /= 0.0_cp) then
                     if (l1 == (m_RiMa + RiSymmMa) .and. m1 == m_RiMa) then
                        rhs1(1,2*lmB-1,threadid)=amp_RiMa*cos(omega_RiMa*time)
                        rhs1(1,2*lmB,threadid)  =amp_RiMa*sin(omega_RiMa*time)
                     end if
                  end if

                  do nR=2,n_r_max-1
                     rhs1(nR,2*lmB-1,threadid)=real(work_LMloc(lm1,nR))
                     rhs1(nR,2*lmB,threadid)  =aimag(work_LMloc(lm1,nR))
                     if ( l_precession .and. l1 == 1 .and. m1 == 1 ) then
                        rhs1(nR,2*lmB-1,threadid)=rhs1(nR,2*lmB-1,threadid)+ &
                        &                         tscheme%wimp_lin(1)*       &
                        &                         prec_fac*sin(oek*time)
                        rhs1(nR,2*lmB,threadid)=rhs1(nR,2*lmB,threadid)-     &
                        &                       tscheme%wimp_lin(1)*prec_fac*&
                        &                       cos(oek*time)
                     end if
                  end do

#ifdef WITH_PRECOND_Z
                  rhs1(:,2*lmB-1,threadid)=zMat_fac(:,nLMB2)*rhs1(:,2*lmB-1,threadid)
                  rhs1(:,2*lmB,threadid)  =zMat_fac(:,nLMB2)*rhs1(:,2*lmB,threadid)
#endif
               end if
            end do

            if ( lmB > lmB0 ) then
               call zMat(nLMB2)%solve(rhs1(:,2*(lmB0+1)-1:2*lmB,threadid),2*(lmB-lmB0))
            end if

            lmB=lmB0
            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               lm1=lm22lm(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)

               if ( l_z10mat .and. lm1 == l1m0 ) then
                  do n_r_out=1,rscheme_oc%n_max
                     z(lm1,n_r_out)=real(rhs(n_r_out))
                  end do
               else if ( l1 /= 0 ) then
                  lmB=lmB+1
                  if ( m1 > 0 ) then
                     do n_r_out=1,rscheme_oc%n_max
                        z(lm1,n_r_out)=cmplx(rhs1(n_r_out,2*lmB-1,threadid), &
                        &                    rhs1(n_r_out,2*lmB,threadid),kind=cp)
                     end do
                  else
                     do n_r_out=1,rscheme_oc%n_max
                        z(lm1,n_r_out)=cmplx(rhs1(n_r_out,2*lmB-1,threadid), &
                        &                    0.0_cp,kind=cp)
                     end do
                  end if
               end if
            end do
            !$omp end task
         end do
         !$omp taskwait
         !$omp end task
      end do       ! end of loop over lm blocks
      !$omp end single
      !$omp taskwait
      !$omp single
      call solve_counter%stop_count(l_increment=.false.)
      !$omp end single

      !-- set cheb modes > rscheme_oc%n_max to zero (dealiazing)
      !$omp do private(n_r_out,lm1) collapse(2)
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         do lm1=llm,ulm
            z(lm1,n_r_out)=zero
         end do
      end do
      !$omp end do
      !$omp end parallel

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dzdt)
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
   subroutine prepareZ_FD(time, tscheme, dzdt, omega_ma, omega_ic, domega_ma_dt, &
              &           domega_ic_dt)

      !-- Input of variables:
      real(cp),            intent(in) :: time
      class(type_tscheme), intent(in) :: tscheme

      !-- Input/output of scalar fields:
      real(cp),           intent(inout) :: omega_ma ! Calculated OC rotation
      real(cp),           intent(inout) :: omega_ic ! Calculated IC rotation
      type(type_tscalar), intent(inout) :: domega_ic_dt
      type(type_tscalar), intent(inout) :: domega_ma_dt
      type(type_tarray),  intent(inout) :: dzdt

      !-- Local variables
      real(cp) :: dom_ma, dom_ic, prec_fac
      integer :: nR, lm_start, lm_stop, lm, l, m, l1m1, l1m0

      if ( .not. l_update_v ) return

      if ( l_precession ) then
         prec_fac=sqrt(8.0_cp*pi*third)*po*oek*oek*sin(prec_angle)
      else
         prec_fac=0.0_cp
      end if

      !-- LU factorisation of the matrix if needed
      if ( .not. lZmat(1) ) then
         call get_zMat_Rdist(tscheme,hdif_V,zMat_FD)
         lZmat(:)=.true.
         if ( l_z10mat ) then
            call get_z10Mat_Rdist(tscheme,hdif_V,z10Mat_FD)
            lZ10mat=.true.
         end if
      end if

      if ( l_rot_ma  .and. (.not. l_SRMA) ) then
         if ( ktopv /= 1 ) then ! Not stress-free
            call tscheme%set_imex_rhs_scalar(dom_ma, domega_ma_dt)
         end if
      end if

      if ( l_rot_ic .and. (.not. l_SRIC) ) then
         if ( kbotv /= 1 ) then ! Not stress-free
            call tscheme%set_imex_rhs_scalar(dom_ic, domega_ic_dt)
         end if
      end if


      !$omp parallel default(shared) private(lm_start,lm_stop,nR, l, m, lm)
      lm_start=1; lm_stop=lm_max
      call get_openmp_blocks(lm_start,lm_stop)
      !$omp barrier

      !-- Assemble the r.h.s.
      call tscheme%set_imex_rhs_ghost(z_ghost, dzdt, lm_start, lm_stop, 1)

      !-- If needed fill z10_ghost
      if ( l_z10mat ) then
         l1m0=st_map%lm2(1,0)
         if ( l1m0 >= lm_start .and. l1m0 <= lm_stop ) then
            do nR=nRstart,nRstop
               z10_ghost(nR)=real(z_ghost(l1m0,nR))
            end do
         end if
      end if

      !-- If precession add one source term
      if ( l_precession ) then
         l1m1=st_map%lm2(1,1)
         if ( l1m1 >= lm_start .and. l1m1 <= lm_stop ) then
            do nR=nRstart,nRstop
               z_ghost(l1m1,nR)=z_ghost(l1m1,nR)+tscheme%wimp_lin(1)*prec_fac* &
               &                cmplx(sin(oek*time),cos(oek*time),cp)
            end do
         end if
      end if

      !-- Boundary conditions
      if ( nRstart==n_r_cmb ) then
         nR = n_r_cmb
         do lm=lm_start,lm_stop
            l = st_map%lm2l(lm)
            if ( l == 0 ) cycle
            m = st_map%lm2m(lm)
            if ( l==1 .and. m==0 .and. l_z10mat ) then
               if ( l_SRMA ) then
                  tOmega_ma1=time+tShift_ma1
                  tOmega_ma2=time+tShift_ma2
                  omega_ma= omega_ma1*cos(omegaOsz_ma1*tOmega_ma1) + &
                  &         omega_ma2*cos(omegaOsz_ma2*tOmega_ma2)
                  z10_ghost(nR)=omega_ma
               else if ( ktopv == 2 .and. l_rot_ma ) then  ! time integration
                  z10_ghost(nR)=dom_ma!/c_dt_z10_ma
               else
                  if ( ktopv == 2 ) z10_ghost(nR)=zero
               end if
               z10_ghost(nR-1)=zero ! Set ghost zone to zero
            end if
            if ( ktopv==2 ) z_ghost(lm,nR)=zero

            if (amp_RiMa /= 0.0_cp .and. l==(m_RiMa+RiSymmMa) .and. m==m_RiMa) then
               z_ghost(lm,nR)=cmplx(amp_RiMa*cos(omega_RiMa*time), &
               &                    amp_RiMa*sin(omega_RiMa*time),cp)
            end if

            z_ghost(lm,nR-1)=zero ! Set ghost zone to zero
         end do
      end if

      if ( nRstop == n_r_icb ) then
         nR=n_r_icb
         do lm=lm_start,lm_stop
            l = st_map%lm2l(lm)
            if ( l == 0 ) cycle
            m = st_map%lm2m(lm)
            if ( l==1 .and. m==0 .and. l_z10mat ) then
               z10_ghost(nR+1)=zero ! Set ghost zone to zero
               if ( l_full_sphere ) then
                  z10_ghost(nR)=zero
               else
                  if ( l_SRIC ) then
                     tOmega_ic1=time+tShift_ic1
                     tOmega_ic2=time+tShift_ic2
                     omega_ic= omega_ic1*cos(omegaOsz_ic1*tOmega_ic1) + &
                     &         omega_ic2*cos(omegaOsz_ic2*tOmega_ic2)
                     z10_ghost(nR)=omega_ic
                  else if ( kbotv == 2 .and. l_rot_ic ) then  ! time integration
                     z10_ghost(nR)=dom_ic!/c_dt_z10_ic
                  else
                     if ( kbotv == 2 ) z10_ghost(nR)=zero
                  end if
               end if
            end if
            if ( l_full_sphere ) then
               z_ghost(lm,nR)=zero
            else
               if ( kbotv==2 ) z_ghost(lm,nR)=zero
            end if

            if (amp_RiIc /= 0.0_cp .and. l==(m_RiIc+RiSymmIc) .and. m==m_RiIc) then
               z_ghost(lm,nR)=cmplx(amp_RiIc*cos(omega_RiIc*time),  &
               &                    amp_RiIc*sin(omega_RiIc*time),cp)
            end if

            z_ghost(lm,nR+1)=zero ! Set ghost zone to zero
         end do
      end if
      !$omp end parallel

   end subroutine prepareZ_FD
!------------------------------------------------------------------------------
   subroutine fill_ghosts_Z(zg)
      !
      ! This subroutine is used to fill the ghosts zones that are located at
      ! nR=n_r_cmb-1 and nR=n_r_icb+1. This is used to properly set the Neuman
      ! boundary conditions. In case Dirichlet BCs are used, a simple first order
      ! extrapolation is employed. This is anyway only used for outputs (like Nusselt
      ! numbers).
      !
      complex(cp), intent(inout) :: zg(lm_max,nRstart-1:nRstop+1)

      !-- Local variables
      integer :: lm, lm_start, lm_stop
      real(cp) :: dr

      if ( .not. l_update_v ) return

      !$omp parallel default(shared) private(lm_start, lm_stop, lm)
      lm_start=1; lm_stop=lm_max
      call get_openmp_blocks(lm_start,lm_stop)

      !-- Handle upper boundary
      if ( nRstart == n_r_cmb ) then
         dr = r(2)-r(1)
         do lm=lm_start,lm_stop
            if ( ktopv == 2 ) then
               zg(lm,nRstart-1)=two*zg(lm,nRstart)-zg(lm,nRstart+1)
            else
               zg(lm,nRstart-1)=zg(lm,nRstart+1)-two*dr*(two*or1(1)+beta(1))* &
               &                zg(lm,nRstart)
            end if
         end do
      end if

      !-- Handle Lower boundary
      if ( nRstop == n_r_icb ) then
         dr = r(n_r_max)-r(n_r_max-1)
         do lm=lm_start,lm_stop
            if ( l_full_sphere ) then
               zg(lm,nRstop+1)=two*zg(lm,nRstop)-zg(lm,nRstop-1)
            else ! Not a full sphere
               if (kbotv == 2) then ! Rigid boundary
                  zg(lm,nRstop+1)=two*zg(lm,nRstop)-zg(lm,nRstop-1)
               else
                  zg(lm,nRstop+1)=zg(lm,nRstop-1)+two*dr* &
                  &               (two*or1(n_r_max)+beta(n_r_max))*zg(lm,nRstop)
               end if
            end if
         end do
      end if
      !$omp end parallel

   end subroutine fill_ghosts_Z
!------------------------------------------------------------------------------
   subroutine updateZ_FD(time, timeNext, z, dz, dzdt, omega_ma, omega_ic,         &
              &          domega_ma_dt, domega_ic_dt, lorentz_torque_ma_dt,        &
              &          lorentz_torque_ic_dt, tscheme,lRmsNext)

      !-- Input of variables:
      real(cp),            intent(in) :: time       ! Current stage time
      real(cp),            intent(in) :: timeNext   ! Next time
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lRmsNext   ! Logical for storing update if (l_RMS.and.l_logNext)

      !-- Input/output of scalar fields:
      type(type_tarray), intent(inout) :: dzdt
      complex(cp),       intent(inout) :: z(lm_max,nRstart:nRstop) ! Toroidal potential

      !-- Output: ds
      complex(cp),        intent(out) :: dz(lm_max,nRstart:nRstop) ! Radial derivative of z
      real(cp),           intent(inout) :: omega_ma        ! Calculated OC rotation
      real(cp),           intent(inout) :: omega_ic        ! Calculated IC rotation
      type(type_tscalar), intent(inout) :: domega_ic_dt
      type(type_tscalar), intent(inout) :: domega_ma_dt
      type(type_tscalar), intent(inout) :: lorentz_torque_ic_dt
      type(type_tscalar), intent(inout) :: lorentz_torque_ma_dt

      !-- Local variables
      real(cp) :: lo_ma, lo_ic
      integer :: lm_start, lm_stop, lm, nR, l

      if ( .not. l_update_v ) return

      if ( l_rot_ma  .and. (.not. l_SRMA) .and. ktopv == 1) then
         call tscheme%set_imex_rhs_scalar(lo_ma, lorentz_torque_ma_dt)
      end if

      if ( l_rot_ic .and. (.not. l_SRIC) .and. kbotv == 1 ) then
         call tscheme%set_imex_rhs_scalar(lo_ic, lorentz_torque_ic_dt)
      end if

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dzdt)
      call tscheme%rotate_imex_scalar(domega_ma_dt)
      call tscheme%rotate_imex_scalar(domega_ic_dt)
      call tscheme%rotate_imex_scalar(lorentz_torque_ma_dt)
      call tscheme%rotate_imex_scalar(lorentz_torque_ic_dt)

      !-- Calculation of the implicit part
      if (  tscheme%istage == tscheme%nstages ) then
         call get_tor_rhs_imp_ghost(timeNext, z_ghost, dz, dzdt, domega_ma_dt,       &
              &                     domega_ic_dt, omega_ic, omega_ma, omega_ic1,     &
              &                     omega_ma1, tscheme, 1, tscheme%l_imp_calc_rhs(1),&
              &                     lRmsNext)
         call update_rot_rates_Rloc(z_ghost, lo_ma, lo_ic, lorentz_torque_ma_dt,   &
              &                     lorentz_torque_ic_dt, omega_ma,                &
              &                     omega_ma1, omega_ic, omega_ic1, 1)
      else
         call get_tor_rhs_imp_ghost(time, z_ghost, dz, dzdt, domega_ma_dt,        &
              &                     domega_ic_dt, omega_ic, omega_ma, omega_ic1,  &
              &                     omega_ma1, tscheme, tscheme%istage+1,         &
              &                     tscheme%l_imp_calc_rhs(tscheme%istage+1), lRmsNext)
         call update_rot_rates_Rloc(z_ghost, lo_ma, lo_ic, lorentz_torque_ma_dt,    &
              &                     lorentz_torque_ic_dt, omega_ma,                 &
              &                     omega_ma1, omega_ic, omega_ic1, tscheme%istage+1)
      end if

      !$omp parallel default(shared) private(lm_start,lm_stop,nR,l,lm)
      lm_start=1; lm_stop=lm_max
      call get_openmp_blocks(lm_start,lm_stop)
      !$omp barrier

      !-- Array copy from z_ghost to z
      do nR=nRstart,nRstop
         do lm=lm_start,lm_stop
            l=st_map%lm2l(lm)
            if ( l==0 ) cycle
            z(lm,nR)=z_ghost(lm,nR)
         end do
      end do
      !$omp end parallel

   end subroutine updateZ_FD
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
      complex(cp), intent(inout) :: z(llm:ulm,n_r_max)
      complex(cp), intent(out) :: dz(llm:ulm,n_r_max)

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
      integer :: lmStart_00, l1, m1, l1m0, l1m1
      integer, pointer :: lm2l(:),lm2m(:), lm2(:,:)
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

      lm2(0:, 0:) => lo_map%lm2
      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m
      lmStart_00 =max(2,llm)

      !$omp parallel default(shared)  private(start_lm, stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)

      !$omp single
      call dct_counter%start_count()
      !$omp end single
      call get_ddr( z, dz, work_LMloc, ulm-llm+1, start_lm-llm+1, &
           &       stop_lm-llm+1, n_r_max, rscheme_oc, l_dct_in=.not. l_in_cheb)
      if ( l_in_cheb ) call rscheme_oc%costf1(z,ulm-llm+1,start_lm-llm+1, &
                            &                 stop_lm-llm+1)
      !$omp barrier
      !$omp single
      call dct_counter%stop_count(l_increment=.false.)
      !$omp end single

      l1m0=lm2(1,0)
      l1m1=lm2(1,1)


      !--- We correct so that the angular moment about axis in the equatorial plane
      !    vanish and the angular moment about the (planetary) rotation axis
      !    is kept constant.
      !$omp single
      if ( l_correct_AMz .and.  l1m0 > 0 .and. lmStart_00 <= l1m0 .and. &
      &    ulm >= l1m0 ) then

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
            work_LMloc(l1m0,n_r)=work_LMloc(l1m0,n_r)-rho0(n_r)*( &
            &               two+four*beta(n_r)*r(n_r) +           &
            &                dbeta(n_r)*r_E_2 +                   &
            &              beta(n_r)*beta(n_r)*r_E_2 )*corr_l1m0
         end do

         if ( ktopv == 2 .and. l_rot_ma ) &
         &    omega_ma=c_z10_omega_ma*real(z(l1m0,n_r_cmb))
         if ( kbotv == 2 .and. l_rot_ic ) &
         &    omega_ic=c_z10_omega_ic*real(z(l1m0,n_r_icb))
         omega_ic1=omega_ic
         omega_ma1=omega_ma

      end if ! l=1,m=0 contained in lm-block ?

      if ( l_correct_AMe .and.  l1m1 > 0 .and. &
      &    lmStart_00 <= l1m1 .and. ulm >= l1m1 ) then

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
            work_LMloc(l1m1,n_r)=work_LMloc(l1m1,n_r)-rho0(n_r)*(    &
            &              two+four*beta(n_r)*r(n_r) +               &
            &                          dbeta(n_r)*r_E_2 +            &
            &                beta(n_r)*beta(n_r)*r_E_2 )*corr_l1m1
         end do
      end if ! l=1,m=1 contained in lm-block ?
      !$omp end single


      if ( istage == 1 ) then
         !$omp do private(n_r,lm,l1,dL)
         do n_r=1,n_r_max
            do lm=llm,ulm
               l1 = lm2l(lm)
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
            do lm=lmStart_00,ulm
               l1 = lm2l(lm)
               m1 = lm2m(lm)
               dL = real(l1*(l1+1),cp)
               Dif(lm)=hdif_V(l1)*dL*or2(n_r)*visc(n_r)* ( work_LMloc(lm,n_r) +  &
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
               call hInt2Tor(Dif,llm,ulm,n_r,lmStart_00,ulm, &
                    &        DifTor2hInt(:,n_r),lo_map)
            end if
         end do
         !$omp end do

      end if

      !$omp end parallel

      if ( ( llm <= l1m0 .and. ulm >= l1m0 ) .and. l_z10mat ) then
         !----- NOTE opposite sign of viscous torque on ICB and CMB:
         if ( .not. l_SRMA .and. ktopv == 2 .and. l_rot_ma ) then
            domega_ma_dt%impl(istage)=visc(1)*( (two*or1(1)+beta(1))* &
            &                         real(z(l1m0,1))-real(dz(l1m0,1)) )
            if ( istage == 1 ) domega_ma_dt%old(istage)=c_dt_z10_ma*real(z(l1m0,1))
         end if
         if ( .not. l_SRIC .and. kbotv == 2 .and. l_rot_ic ) then
            domega_ic_dt%impl(istage)=-visc(n_r_max)* ( (two*or1(n_r_max)+   &
            &                          beta(n_r_max))*real(z(l1m0,n_r_max))- &
            &                          real(dz(l1m0,n_r_max)) )
            if ( istage == 1 ) domega_ic_dt%old(istage)=c_dt_z10_ic* &
            &                                           real(z(l1m0,n_r_max))
         end if
      end if

      !--- Note: from ddz=work_LMloc only the axisymmetric contributions are needed
      !    beyond this point for the TO calculation.
      !    Parallization note: Very likely, all axisymmetric modes m=0 are
      !    located on the first processor #0.
      if ( l_TO ) then
         !$omp parallel do default(shared) private(n_r,lm,l1,m1)
         do n_r=1,n_r_max
            ddzASL_loc(:,n_r)=0.0_cp
            do lm=lmStart_00,ulm
               l1=lm2l(lm)
               m1=lm2m(lm)
               if ( m1 == 0 ) ddzASL_loc(l1+1,n_r)=real(work_LMloc(lm,n_r))
            end do
         end do
         !$omp end parallel do

         do n_r=1,n_r_max
#ifdef WITH_MPI
            call MPI_Allreduce(ddzASL_loc(:,n_r), ddzASL(:,n_r), l_max+1, &
                 &             MPI_DEF_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
            ddzASL(:,n_r)=ddzASL_loc(:,n_r)
#endif
         end do
      end if

   end subroutine get_tor_rhs_imp
!------------------------------------------------------------------------------
   subroutine get_tor_rhs_imp_ghost(time, zg, dz, dzdt, domega_ma_dt, domega_ic_dt,  &
              &                     omega_ic, omega_ma, omega_ic1, omega_ma1,        &
              &                     tscheme, istage, l_calc_lin, lRmsNext)

      !-- Input variables
      real(cp),            intent(in) :: time
      integer,             intent(in) :: istage
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lRmsNext
      logical,             intent(in) :: l_calc_lin

      !-- Output variable
      type(type_tarray),  intent(inout) :: dzdt
      type(type_tscalar), intent(inout) :: domega_ic_dt
      type(type_tscalar), intent(inout) :: domega_ma_dt
      real(cp),    intent(inout) :: omega_ic
      real(cp),    intent(inout) :: omega_ma
      real(cp),    intent(inout) :: omega_ic1
      real(cp),    intent(inout) :: omega_ma1
      complex(cp), intent(inout) :: zg(lm_max,nRstart-1:nRstop+1)
      complex(cp), intent(out) :: dz(lm_max,nRstart:nRstop)

      !-- Local variables
      real(cp) :: angular_moment(3)   ! total angular momentum
      real(cp) :: angular_moment_oc(3)! x,y,z component of outer core angular mom.
      real(cp) :: angular_moment_ic(3)! x,y,z component of inner core angular mom.
      real(cp) :: angular_moment_ma(3)! x,y,z component of mantle angular mom.
      complex(cp) :: z10(nRstart:nRstop), z11(nRstart:nRstop)
      complex(cp) :: corr_l1m0, corr_l1m1
      real(cp) :: r_E_2, nomi, dL, prec_fac
      integer :: n_r, lm, start_lm, stop_lm, i
      integer :: l, m, l1m0, l1m1
      real(cp) :: ddzASL_loc(0:l_max,n_r_max)
      complex(cp) :: work_Rloc(lm_max,nRstart:nRstop)

      if ( l_precession ) then
         prec_fac=sqrt(8.0_cp*pi*third)*po*oek*oek*sin(prec_angle)
      else
         prec_fac = 0.0_cp
      end if

      !$omp parallel default(shared)  private(start_lm, stop_lm, n_r, lm, l, m)
      start_lm=1; stop_lm=lm_max
      call get_openmp_blocks(start_lm,stop_lm)

      !$omp single
      call dct_counter%start_count()
      !$omp end single
      call get_ddr_ghost(zg, dz, work_Rloc, lm_max,start_lm, stop_lm,  nRstart, nRstop, &
           &             rscheme_oc)
      !$omp single
      call dct_counter%stop_count(l_increment=.false.)
      !$omp end single
      !$omp barrier
      !$omp end parallel

      l1m0=st_map%lm2(1,0)
      l1m1=st_map%lm2(1,1)

      !--- We correct so that the angular moment about axis in the equatorial plane
      !    vanish and the angular moment about the (planetary) rotation axis
      !    is kept constant.
      if ( l_correct_AMz ) then

         z10(nRstart:nRstop)=zg(l1m0,nRstart:nRstop)
         call get_angular_moment_Rloc(z10,z11,omega_ic,omega_ma,          &
              &                       angular_moment_oc,angular_moment_ic,&
              &                       angular_moment_ma)
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
         do n_r=nRstart,nRstop
            r_E_2=r(n_r)*r(n_r)
            zg(l1m0,n_r)=zg(l1m0,n_r)  - rho0(n_r)*r_E_2*corr_l1m0
            dz(l1m0,n_r)=dz(l1m0,n_r) - rho0(n_r)*(          &
            &            two*r(n_r)+r_E_2*beta(n_r))*corr_l1m0
            work_Rloc(l1m0,n_r)=work_Rloc(l1m0,n_r)-rho0(n_r)*( &
            &               two+four*beta(n_r)*r(n_r) +         &
            &                dbeta(n_r)*r_E_2 +                 &
            &              beta(n_r)*beta(n_r)*r_E_2 )*corr_l1m0
         end do

         if ( ktopv == 2 .and. l_rot_ma .and. nRstart==n_r_cmb ) &
         &    omega_ma=c_z10_omega_ma*real(zg(l1m0,n_r_cmb))
         if ( kbotv == 2 .and. l_rot_ic .and. nRstop==n_r_icb ) &
         &    omega_ic=c_z10_omega_ic*real(zg(l1m0,n_r_icb))
         omega_ic1=omega_ic
         omega_ma1=omega_ma

      end if

      if ( l_correct_AMe ) then

         z11(nRstart:nRstop)=zg(l1m1,nRstart:nRstop)
         call get_angular_moment_Rloc(z10,z11,omega_ic,omega_ma,          &
              &                       angular_moment_oc,angular_moment_ic,&
              &                       angular_moment_ma)
         do i=1,3
            angular_moment(i)=angular_moment_oc(i) + angular_moment_ic(i) + &
            &                 angular_moment_ma(i)
         end do
         corr_l1m1=cmplx(angular_moment(1),-angular_moment(2),kind=cp) / &
         &         (two*y11_norm*c_moi_oc)

         !-------- Correct z(2,n_r) and z(l_max+2,n_r) plus the respective
         !         derivatives:
         do n_r=nRstart,nRstop
            r_E_2=r(n_r)*r(n_r)
            zg(l1m1,n_r)=zg(l1m1,n_r)  -  rho0(n_r)*r_E_2*corr_l1m1
            dz(l1m1,n_r)=dz(l1m1,n_r) -  rho0(n_r)*(            &
            &            two*r(n_r)+r_E_2*beta(n_r))*corr_l1m1
            work_Rloc(l1m1,n_r)=work_Rloc(l1m1,n_r)-rho0(n_r)*(      &
            &              two+four*beta(n_r)*r(n_r) +               &
            &                          dbeta(n_r)*r_E_2 +            &
            &                beta(n_r)*beta(n_r)*r_E_2 )*corr_l1m1
         end do
      end if ! l=1,m=1 contained in lm-block ?

      !$omp parallel default(shared) private(start_lm, stop_lm, n_r, lm, l, m)
      start_lm=1; stop_lm=lm_max
      call get_openmp_blocks(start_lm,stop_lm)

      if ( istage == 1 ) then
         do n_r=nRstart,nRstop
            do lm=start_lm,stop_lm
               l = st_map%lm2l(lm)
               dL = real(l*(l+1),cp)
               dzdt%old(lm,n_r,istage)=dL*or2(n_r)*zg(lm,n_r)
            end do
         end do
      end if

      if ( l_calc_lin .or. (tscheme%istage==tscheme%nstages .and. lRmsNext)) then

         do n_r=nRstart,nRstop
            do lm=start_lm,stop_lm
               l = st_map%lm2l(lm)
               if ( l == 0 ) cycle
               m = st_map%lm2m(lm)
               dL = real(l*(l+1),cp)
               Dif(lm)=hdif_V(l)*dL*or2(n_r)*visc(n_r)* ( work_Rloc(lm,n_r) +    &
               &         (dLvisc(n_r)-beta(n_r))    *              dz(lm,n_r) -  &
               &         ( dLvisc(n_r)*beta(n_r)+two*dLvisc(n_r)*or1(n_r)        &
               &          + dL*or2(n_r)+dbeta(n_r)+two*beta(n_r)*or1(n_r) )*     &
               &                                                   zg(lm,n_r) )

               dzdt%impl(lm,n_r,istage)=Dif(lm)
               if ( l_precession .and. l==1 .and. m==1 ) then
                  dzdt%impl(lm,n_r,istage)=dzdt%impl(lm,n_r,istage)+prec_fac*cmplx( &
                  &                        sin(oek*time),-cos(oek*time),kind=cp)
               end if
            end do
            if ( lRmsNext .and. tscheme%istage==tscheme%nstages ) then
               call hInt2Tor(Dif,1,lm_max,n_r,start_lm,stop_lm,DifTor2hInt(:,n_r),st_map)
            end if
         end do

      end if
      !$omp end parallel


      if ( l_z10mat ) then
         !----- NOTE opposite sign of viscous torque on ICB and CMB:
         if ( .not. l_SRMA .and. ktopv==2 .and. l_rot_ma .and. nRstart==n_r_cmb ) then
            domega_ma_dt%impl(istage)=visc(1)*( (two*or1(1)+beta(1))* &
            &                         real(zg(l1m0,1))-real(dz(l1m0,1)) )
            if ( istage == 1 ) domega_ma_dt%old(istage)=c_dt_z10_ma*real(zg(l1m0,1))
         end if
         if ( .not. l_SRIC .and. kbotv==2 .and. l_rot_ic .and. nRstop==n_r_icb ) then
            domega_ic_dt%impl(istage)=-visc(n_r_max)* ( (two*or1(n_r_max)+   &
            &                          beta(n_r_max))*real(zg(l1m0,n_r_max))- &
            &                          real(dz(l1m0,n_r_max)) )
            if ( istage == 1 ) domega_ic_dt%old(istage)=c_dt_z10_ic* &
            &                                           real(zg(l1m0,n_r_max))
         end if
      end if

      !--- Note: from ddz=work_LMloc only the axisymmetric contributions are needed
      !    beyond this point for the TO calculation.
      !    Parallization note: Very likely, all axisymmetric modes m=0 are
      !    located on the first processor #0.
      if ( l_TO ) then
         !$omp parallel do default(shared) private(n_r,lm,l,m)
         do n_r=nRstart,nRstop
            ddzASL_loc(:,n_r)=0.0_cp
            do lm=1,lm_max
               l=st_map%lm2l(lm)
               m=st_map%lm2m(lm)
               if ( m == 0 ) ddzASL_loc(l,n_r)=real(work_Rloc(lm,n_r))
            end do
         end do
         !$omp end parallel do

         !-- This is quite ugly here:
         do l=0,l_max
            call allgather_from_rloc(ddzASL_loc(l,:),ddzASL(l+1,:))
         end do
      end if

   end subroutine get_tor_rhs_imp_ghost
!------------------------------------------------------------------------------
   subroutine assemble_tor(time, z, dz, dzdt, domega_ic_dt, domega_ma_dt,        &
              &            lorentz_torque_ic_dt, lorentz_torque_ma_dt, omega_ic, &
              &            omega_ma, omega_ic1, omega_ma1, lRmsNext, tscheme)
      !
      ! This subroutine is used to assemble the toroidal flow potential when an IMEX
      ! RK time scheme with an assembly stage is employed (LM-distributed version).
      !

      !-- Input variable
      class(type_tscheme), intent(in) :: tscheme
      real(cp),            intent(in) :: time
      logical,             intent(in) :: lRmsNext

      !-- Output variables
      complex(cp),        intent(inout) :: z(llm:ulm,n_r_max)
      complex(cp),        intent(out) :: dz(llm:ulm,n_r_max)
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
      complex(cp) :: top_val(llm:ulm), bot_val(llm:ulm)
      integer :: n_r, lm, l1, m1, lmStart_00, l1m0
      real(cp) :: fac_top, fac_bot, dom_ma, dom_ic, lo_ma, lo_ic
      integer, pointer :: lm2l(:),lm2(:,:), lm2m(:)
      real(cp) :: dL

      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m
      lm2(0:,0:) => lo_map%lm2
      lmStart_00 =max(2,llm)
      l1m0       =lm2(1,0)

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

      !-- Store the assembled quantity in work_LMloc
      call tscheme%assemble_imex(work_LMloc, dzdt)

      !-- Now get the toroidal potential from the assembly
      !$omp parallel default(shared)
      !$omp do private(n_r,lm,l1,m1,dL)
      do n_r=2,n_r_max-1
         do lm=lmStart_00,ulm
            l1 = lm2l(lm)
            m1 = lm2m(lm)
            dL = real(l1*(l1+1),cp)
            if ( m1 == 0 ) then
               z(lm,n_r) = r(n_r)*r(n_r)/dL*cmplx(real(work_LMloc(lm,n_r)),0.0_cp,cp)
            else
               z(lm,n_r) = r(n_r)*r(n_r)/dL*work_LMloc(lm,n_r)
            end if
         end do
      end do
      !$omp end do

      !$omp do private(lm)
      do lm=llm,ulm
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
            !$omp do private(lm)
            do lm=lmStart_00,ulm
               z(lm,1)      =top_val(lm)
               z(lm,n_r_max)=bot_val(lm)
            end do
            !$omp end do
         else
            fac_top=-two*or1(1)-beta(1)
            !$omp do private(lm)
            do lm=lmStart_00,ulm
               call rscheme_oc%robin_bc(one, fac_top, zero, 0.0_cp, one, bot_val(lm), z(lm,:))
            end do
            !$omp end do
         end if
      else
         if ( ktopv /= 1 .and. kbotv /= 1 ) then ! Rigid BCs
            !$omp do private(lm)
            do lm=lmStart_00,ulm
               z(lm,1)      =top_val(lm)
               z(lm,n_r_max)=bot_val(lm)
            end do
            !$omp end do
         else if ( ktopv == 1 .and. kbotv /= 1 ) then ! Stress-free top and rigid bottom
            fac_top=-two*or1(1)-beta(1)
            !$omp do private(lm)
            do lm=lmStart_00,ulm
               call rscheme_oc%robin_bc(one, fac_top, zero, 0.0_cp, one, bot_val(lm), z(lm,:))
            end do
            !$omp end do
         else if ( kbotv == 1 .and. ktopv /= 1 ) then ! Stress-free bot and rigid top
            fac_bot=-two*or1(n_r_max)-beta(n_r_max)
            !$omp do private(lm)
            do lm=lmStart_00,ulm
               call rscheme_oc%robin_bc(0.0_cp, one, top_val(lm), one, fac_bot, zero, z(lm,:))
            end do
            !$omp end do
         else if ( ktopv == 1 .and. kbotv == 1 ) then ! Stress-free at both boundaries
            fac_bot=-two*or1(n_r_max)-beta(n_r_max)
            fac_top=-two*or1(1)-beta(1)
            !$omp do private(lm)
            do lm=lmStart_00,ulm
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
   subroutine assemble_tor_Rloc(time, z, dz, dzdt, domega_ic_dt, domega_ma_dt,        &
              &                 lorentz_torque_ic_dt, lorentz_torque_ma_dt, omega_ic, &
              &                 omega_ma, omega_ic1, omega_ma1, lRmsNext, tscheme)
      !
      ! This subroutine is used when a IMEX Runge-Kutta scheme with an assembly
      ! stage is employed (R-distributed version)
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      real(cp),            intent(in) :: time
      logical,             intent(in) :: lRmsNext

      !-- Output variables
      complex(cp),        intent(inout) :: z(lm_max,nRstart:nRstop)
      complex(cp),        intent(out)   :: dz(lm_max,nRstart:nRstop)
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
      complex(cp) :: work_Rloc(lm_max,nRstart:nRstop)
      integer :: start_lm, stop_lm, lm, n_r, l, m
      real(cp) :: dLh, lo_ic, lo_ma, dom_ma, dom_ic, r2

      call tscheme%assemble_imex(work_Rloc, dzdt)
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

      !$omp parallel default(shared) private(start_lm, stop_lm, l, m, dLh)
      start_lm=1; stop_lm=lm_max
      call get_openmp_blocks(start_lm,stop_lm)
      !$omp barrier

      do n_r=nRstart,nRstop
         r2=r(n_r)*r(n_r)
         do lm=start_lm,stop_lm
            l = st_map%lm2l(lm)
            if ( l == 0 ) cycle
            m = st_map%lm2m(lm)
            dLh=real(l*(l+1),cp)
            if ( m == 0 ) then
               z(lm,n_r)=cmplx(real(work_Rloc(lm,n_r)),0.0_cp,cp)*r2/dLh
            else
               z(lm,n_r)=work_Rloc(lm,n_r)*r2/dLh
            end if
         end do
      end do

      !-- Boundary points
      if ( nRstart == n_r_cmb ) then
         n_r=n_r_cmb
         do lm=start_lm,stop_lm
            l=st_map%lm2l(lm)
            if ( l == 0 ) cycle
            m=st_map%lm2m(lm)
            if ( l == 1 .and. m == 0 ) then
               if ( l_SRMA ) then
                  tOmega_ma1=time+tShift_ma1
                  tOmega_ma2=time+tShift_ma2
                  omega_ma= omega_ma1*cos(omegaOsz_ma1*tOmega_ma1) + &
                  &         omega_ma2*cos(omegaOsz_ma2*tOmega_ma2)
                  z(lm,n_r)=cmplx(omega_ma/c_z10_omega_ma,0.0_cp,kind=cp)
               else if ( ktopv == 2 .and. l_rot_ma ) then
                  z(lm,n_r)=cmplx(dom_ma/c_dt_z10_ma,0.0_cp,kind=cp)
               else
                  if ( ktopv == 2 ) z(lm,n_r)=zero
               end if
            else
               if ( ktopv == 2 ) z(lm,n_r)=zero
            end if
         end do
      end if

      if ( nRstop == n_r_icb ) then
         n_r=n_r_icb
         do lm=start_lm,stop_lm
            l=st_map%lm2l(lm)
            if ( l == 0 ) cycle
            m=st_map%lm2m(lm)
            if ( l_full_sphere ) then
               z(lm,n_r)=zero
            else
               if ( l == 1 .and. m == 0 ) then
                  if ( l_SRIC ) then
                     tOmega_ic1=time+tShift_ic1
                     tOmega_ic2=time+tShift_ic2
                     omega_ic= omega_ic1*cos(omegaOsz_ic1*tOmega_ic1) + &
                     &         omega_ic2*cos(omegaOsz_ic2*tOmega_ic2)
                     z(lm,n_r)=cmplx(omega_ic/c_z10_omega_ic,0.0_cp,kind=cp)
                  else if ( kbotv == 2 .and. l_rot_ic ) then  ! time integration
                     z(lm,n_r)=cmplx(dom_ic/c_dt_z10_ic,0.0_cp,kind=cp)
                  else
                     if ( kbotv == 2 ) z(lm,n_r)=zero
                  end if
               else
                  if ( kbotv == 2 ) z(lm,n_r)=zero
               end if
            end if
         end do
      end if

      call bulk_to_ghost(z, z_ghost, 1, nRstart, nRstop, lm_max, start_lm, stop_lm)
      !$omp end parallel

      call exch_ghosts(z_ghost, lm_max, nRstart, nRstop, 1)
      call fill_ghosts_Z(z_ghost)

      !-- Finally call the construction of the implicit terms for the first stage
      !-- of next iteration
      call get_tor_rhs_imp_ghost(time, z_ghost, dz, dzdt, domega_ma_dt,           &
           &                     domega_ic_dt, omega_ic, omega_ma, omega_ic1,     &
           &                     omega_ma1, tscheme, 1, tscheme%l_imp_calc_rhs(1),&
           &                     lRmsNext)

      call update_rot_rates_Rloc(z_ghost, lo_ma, lo_ic, lorentz_torque_ma_dt,   &
           &                     lorentz_torque_ic_dt, omega_ma,                &
           &                     omega_ma1, omega_ic, omega_ic1, 1)

   end subroutine assemble_tor_Rloc
!------------------------------------------------------------------------------
   subroutine update_rot_rates(z, lo_ma, lo_ic, lorentz_torque_ma_dt,       &
              &                lorentz_torque_ic_dt, omega_ma, omega_ma1,   &
              &                omega_ic, omega_ic1, istage)

      !-- Input variables
      complex(cp), intent(in) :: z(llm:ulm,n_r_max)
      real(cp),    intent(in) :: lo_ma, lo_ic   ! RHS when stress-free BCs are used
      integer,     intent(in) :: istage

      !-- Output variables
      type(type_tscalar), intent(inout) :: lorentz_torque_ma_dt
      type(type_tscalar), intent(inout) :: lorentz_torque_ic_dt
      real(cp),           intent(out) :: omega_ma, omega_ma1
      real(cp),           intent(out) :: omega_ic, omega_ic1

      !-- Local variables
      integer :: l1m0, lmStart_00

      lmStart_00 =max(2,llm)
      l1m0=lo_map%lm2(1,0)

      !--- Update of inner core and mantle rotation:
      !$omp single
      !if ( l1m0 > 0 .and. lmStart_00 <= l1m0 .and. ulm >= l1m0 ) then
      if ( ( llm <= l1m0 .and. ulm >= l1m0 ) .and. l_z10mat ) then
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
      !$omp end single

   end subroutine update_rot_rates
!------------------------------------------------------------------------------
   subroutine update_rot_rates_Rloc(z, lo_ma, lo_ic, lorentz_torque_ma_dt,       &
              &                     lorentz_torque_ic_dt, omega_ma, omega_ma1,   &
              &                     omega_ic, omega_ic1, istage)

      !-- Input variables
      complex(cp), intent(in) :: z(lm_max,nRstart-1:nRstop+1)
      real(cp),    intent(in) :: lo_ma, lo_ic   ! RHS when stress-free BCs are used
      integer,     intent(in) :: istage

      !-- Output variables
      type(type_tscalar), intent(inout) :: lorentz_torque_ma_dt
      type(type_tscalar), intent(inout) :: lorentz_torque_ic_dt
      real(cp),           intent(out) :: omega_ma, omega_ma1
      real(cp),           intent(out) :: omega_ic, omega_ic1

      !-- Local variables
      integer :: l1m0

      l1m0=st_map%lm2(1,0)

      !--- Update of inner core and mantle rotation:
      !$omp single
      if ( l_z10mat ) then
         if ( l_rot_ma .and. .not. l_SRMA .and. (nRstart==n_r_cmb) ) then
            if ( ktopv == 1 ) then  ! free slip, explicit time stepping of omega !
               call get_rot_rates(omega_ma, lorentz_torque_ma_dt%old(istage))
               omega_ma=lo_ma
            else if ( ktopv == 2 ) then ! no slip, omega given by z10
               omega_ma=c_z10_omega_ma*real(z(l1m0,n_r_cmb))
            end if
            omega_ma1=omega_ma
         end if
         if ( l_rot_ic .and. .not. l_SRIC .and. (nRstop==n_r_icb) ) then
            if ( kbotv == 1 ) then  ! free slip, explicit time stepping of omega !
               call get_rot_rates(omega_ic, lorentz_torque_ic_dt%old(istage))
               omega_ic=lo_ic
            else if ( kbotv == 2 ) then ! no slip, omega given by z10
               omega_ic=c_z10_omega_ic*real(z(l1m0,n_r_icb))
            end if
            omega_ic1=omega_ic
         end if
      end if  ! l=1,m=0 contained in block ?
      !$omp end single

   end subroutine update_rot_rates_Rloc
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
      !  inversion matrix ``z10mat`` for the implicit time step of the
      !  toroidal velocity potential z of degree :math:`\ell=1` and order :math:`m=0`.
      !  This differs from the the normal zmat only if either the ICB or
      !  CMB have no-slip boundary condition and inner core or mantle are
      !  chosen to rotate freely (either ``kbotv=1`` and/or ``ktopv=1``).
      !

      class(type_tscheme), intent(in) :: tscheme ! Time step internal
      real(cp),            intent(in) :: hdif    ! Value of hyperdiffusivity in zMat terms
      integer,             intent(in) :: l       ! Variable to loop over degrees

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
      !  ``zmat(i,j)`` for the Navier-Stokes equation.
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
!------------------------------------------------------------------------------
   subroutine get_z10Mat_Rdist(tscheme,hdif,zMat)
      !
      ! This subroutine is employed to construct the matrix for the z(l=1,m=0) mode
      ! when the parallel F.D. solver is used.
      !

      class(type_tscheme), intent(in) :: tscheme        ! Time step internal
      real(cp),            intent(in) :: hdif(0:l_max)    ! Value of hyperdiffusivity in zMat terms

      !-- Output variables:
      type(type_tri_par), intent(inout) :: zMat ! Tridiag matrix

      !-- Local variables:
      integer :: nR, l
      real(cp) :: dLh, dr

      l=1 ! This is a matrix for l=1,m=0 only
      dLh=real(l*(l+1),kind=cp)

      !-- Bulk points: we fill all the points: this is then easier to handle
      !-- Neumann boundary conditions
      do nR=1,n_r_max
         zMat%diag(l,nR)=dLh*or2(nR)-tscheme%wimp_lin(1)*hdif(l)*dLh* &
         &     visc(nR)*or2(nR) * (              rscheme_oc%ddr(nR,1) &
         &   + (dLvisc(nR)- beta(nR)) *           rscheme_oc%dr(nR,1) &
         &      - ( dLvisc(nR)*beta(nR)+two*dLvisc(nR)*or1(nR)        &
         &          +dLh*or2(nR)+dbeta(nR)+two*beta(nR)*or1(nR)       &
         &                             ) )
         zMat%low(l,nR)=-tscheme%wimp_lin(1)*hdif(l)*dLh*visc(nR)*or2(nR) * ( &
         &      rscheme_oc%ddr(nR,0)+ (dLvisc(nR)- beta(nR)) *rscheme_oc%dr(nR,0) )
         zMat%up(l,nR) =-tscheme%wimp_lin(1)*hdif(l)*dLh*visc(nR)*or2(nR) * ( &
         &      rscheme_oc%ddr(nR,2)+ (dLvisc(nR)- beta(nR)) *rscheme_oc%dr(nR,2) )
      end do

      !-- Boundary conditions:
      !----- CMB condition:
      !        Note opposite sign of viscous torques (-dz+(2/r+beta) z )*visc
      !        for CMB and ICB!
      if ( ktopv == 1 ) then  ! free slip
         zMat%up(l,1)  =zMat%up(l,1)+zMat%low(l,1)
         zMat%diag(l,1)=zMat%diag(l,1)-two*(r(2)-r(1))*(two*or1(1)+beta(1))* &
         &              zMat%low(l,1)
      else if ( ktopv == 2 ) then ! no slip
         if ( l_SRMA ) then
            zMat%diag(l,1)=c_z10_omega_ma
            zMat%low(l,1) =0.0_cp
            zMat%up(l,1)  =0.0_cp
         else if ( l_rot_ma ) then
            !-- I don't know what to here except approximating via a first order
            !-- approximation. There is no easy way to use ghost zones here.
            !-- Using a three point stencil would be possible but this would
            !-- yield a pentadiagonal matrix here.
            dr = one/(r(2)-r(1))
            zMat%diag(l,1)=c_dt_z10_ma-tscheme%wimp_lin(1)*visc(1)*( &
            &              (two*or1(1)+beta(1))-dr )
            zMat%up(l,1)  =tscheme%wimp_lin(1)*visc(1)*dr
            zMat%low(l,1) =0.0_cp
            !TBD now this is not in place
            call abortRun('in update Z: not implemented yet')
         else
            zMat%diag(l,1)=one
            zMat%low(l,1) =0.0_cp
            zMat%up(l,1)  =0.0_cp
         end if
      end if

      !----- ICB condition:
      if ( l_full_sphere ) then
         zMat%diag(l,n_r_max)=one
         zMat%low(l,n_r_max) =0.0_cp
         zMat%up(l,n_r_max)  =0.0_cp
      else
         if ( kbotv == 1 ) then  ! free slip
            zMat%low(l,n_r_max) =zMat%low(l,n_r_max)+zMat%up(l,n_r_max)
            zMat%diag(l,n_r_max)=zMat%diag(l,n_r_max)+two*(r(n_r_max)-r(n_r_max-1))*&
            &                    (two*or1(n_r_max)+beta(n_r_max))*zMat%up(l,n_r_max)
         else if ( kbotv == 2 ) then ! no slip
            if ( l_SRIC ) then
               zMat%diag(l,n_r_max)=c_z10_omega_ic
               zMat%low(l,n_r_max) =0.0_cp
               zMat%up(l,n_r_max)  =0.0_cp
            else if ( l_rot_ic ) then     !  time integration of z10
               !-- Right now local first order: don't know how to handle it
               !-- otherwise
               dr = one/(r(n_r_max)-r(n_r_max-1))
               zMat%diag(l,n_r_max)=c_dt_z10_ic+tscheme%wimp_lin(1)*visc(n_r_max)*( &
               &                    two*or1(n_r_max)+beta(n_r_max)-dr)
               zMat%low(l,n_r_max) =tscheme%wimp_lin(1)*visc(n_r_max)*dr
               zMat%up(l,n_r_max)  =0.0_cp
            else
               zMat%diag(l,n_r_max)=one
               zMat%low(l,n_r_max) =0.0_cp
               zMat%up(l,n_r_max)  =0.0_cp
            end if
         end if
      end if

      !-- LU-decomposition of z10mat:
      call zMat%prepare_mat()

   end subroutine get_z10Mat_Rdist
!-------------------------------------------------------------------------------
   subroutine get_zMat_Rdist(tscheme,hdif,zMat)
      !
      !  This subroutine is used to construct the z matrices when the parallel
      !  F.D. solver is employed.
      !

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme       ! time step
      real(cp),            intent(in) :: hdif(0:l_max) ! Hyperdiffusivity

      !-- Output variables:
      type(type_tri_par), intent(inout) :: zMat

      !-- local variables:
      integer :: nR, l
      real(cp) :: dLh


      !-- Bulk points: we fill all the points: this is then easier to handle
      !-- Neumann boundary conditions
      do nR=1,n_r_max
         do l=1,l_max
            dLh=real(l*(l+1),kind=cp)
            zMat%diag(l,nR)=dLh*or2(nR)-tscheme%wimp_lin(1)*hdif(l)*dLh* &
            &     visc(nR)*or2(nR) * (              rscheme_oc%ddr(nR,1) &
            &   + (dLvisc(nR)- beta(nR)) *           rscheme_oc%dr(nR,1) &
            &      - ( dLvisc(nR)*beta(nR)+two*dLvisc(nR)*or1(nR)        &
            &          +dLh*or2(nR)+dbeta(nR)+two*beta(nR)*or1(nR)       &
            &                             ) )
            zMat%low(l,nR)=-tscheme%wimp_lin(1)*hdif(l)*dLh*visc(nR)*or2(nR) * ( &
            &      rscheme_oc%ddr(nR,0)+ (dLvisc(nR)- beta(nR)) *rscheme_oc%dr(nR,0) )
            zMat%up(l,nR) =-tscheme%wimp_lin(1)*hdif(l)*dLh*visc(nR)*or2(nR) * ( &
            &      rscheme_oc%ddr(nR,2)+ (dLvisc(nR)- beta(nR)) *rscheme_oc%dr(nR,2) )
         end do
      end do

      !----- Boundary conditions, see above:
      do l=1,l_max
         if ( ktopv == 1 ) then  ! free slip !
            zMat%up(l,1)  =zMat%up(l,1)+zMat%low(l,1)
            zMat%diag(l,1)=zMat%diag(l,1)-two*(r(2)-r(1))*(two*or1(1)+beta(1))* &
            &              zMat%low(l,1)
         else                    ! no slip, note exception for l=1,m=0
            zMat%diag(l,1)=one
            zMat%up(l,1)  =0.0_cp
            zMat%low(l,1) =0.0_cp
         end if

         if ( l_full_sphere ) then
            zMat%diag(l,n_r_max)=one
            zMat%low(l,n_r_max) =0.0_cp
            zMat%up(l,n_r_max)  =0.0_cp
         else
            if ( kbotv == 1 ) then  ! free slip !
               zMat%low(l,n_r_max)=zMat%low(l,n_r_max)+zMat%up(l,n_r_max)
               zMat%diag(l,n_r_max)=zMat%diag(l,n_r_max)+two*(r(n_r_max)-r(n_r_max-1))*&
               &                    (two*or1(n_r_max)+beta(n_r_max))*zMat%up(l,n_r_max)
            else
               zMat%diag(l,n_r_max)=one
               zMat%low(l,n_r_max) =0.0_cp
               zMat%up(l,n_r_max)  =0.0_cp
            end if
         end if
      end do

      !-- LU decomposition:
      call zMat%prepare_mat()

   end subroutine get_zMat_Rdist
!-------------------------------------------------------------------------------
end module updateZ_mod
