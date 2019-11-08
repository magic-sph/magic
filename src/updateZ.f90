#include "perflib_preproc.cpp"
module updateZ_mod

   use init_fields
   use omp_lib
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, lm_max, l_max
   use radial_data, only: n_r_cmb, n_r_icb
   use radial_functions, only: visc, or1, or2, rscheme_oc, dLvisc, beta, &
       &                       rho0, r_icb, r_cmb, r, beta, dbeta
   use physical_parameters, only: kbotv, ktopv, LFfac, prec_angle, po, oek
   use num_param, only: AMstart, dct_counter, solve_counter
   use torsional_oscillations, only: ddzASL
   use blocking, only: lo_sub_map, lo_map, st_map, st_sub_map, llm, ulm
   use horizontal_data, only: dLh, hdif_V
   use logic, only: l_rot_ma, l_rot_ic, l_SRMA, l_SRIC, l_z10mat, l_precession, &
       &            l_correct_AMe, l_correct_AMz, l_update_v, l_TO,             &
       &            l_finite_diff, l_full_sphere
   use RMS, only: DifTor2hInt
   use constants, only: c_lorentz_ma, c_lorentz_ic, c_dt_z10_ma, c_dt_z10_ic, &
       &                c_moi_ma, c_moi_ic, c_z10_omega_ma, c_z10_omega_ic,   &
       &                c_moi_oc, y10_norm, y11_norm, zero, one, two, four,   &
       &                pi, third
   use parallel_mod
   use communications, only:get_global_sum
   use outRot, only: get_angular_moment
   use fieldsLast, only: domega_ic_dt, domega_ma_dt, lorentz_torque_ic_dt, &
       &                 lorentz_torque_ma_dt
   use RMS_helpers, only: hInt2Tor
   use radial_der, only: get_ddr
   use fields, only: work_LMloc
   use useful, only: abortRun
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray
   use special
   use dense_matrices
   use real_matrices
   use band_matrices

   implicit none

   private

   !-- Input of recycled work arrays:
   complex(cp), allocatable :: rhs1(:,:,:) ! RHS for other modes
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

   public :: updateZ, initialize_updateZ, finalize_updateZ, get_tor_rhs_imp

contains

   subroutine initialize_updateZ

      integer, pointer :: nLMBs2(:)
      integer :: ll, n_bands

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

         call z10Mat%initialize(n_bands,n_r_max,l_pivot=.true.)
         do ll=1,nLMBs2(1+rank)
            call zMat(ll)%initialize(n_bands,n_r_max,l_pivot=.true.)
         end do
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
      allocate( lZmat(0:l_max) )
      bytes_allocated = bytes_allocated+(l_max+1)*SIZEOF_LOGICAL

      allocate( Dif(llm:ulm) )
      bytes_allocated=bytes_allocated+(ulm-llm+1)*SIZEOF_DEF_COMPLEX

#ifdef WITHOMP
      maxThreads=omp_get_max_threads()
#else
      maxThreads=1
#endif

      allocate(rhs1(n_r_max,lo_sub_map%sizeLMB2max,0:maxThreads-1))
      bytes_allocated=bytes_allocated+n_r_max*maxThreads* &
      &               lo_sub_map%sizeLMB2max*SIZEOF_DEF_COMPLEX

      AMstart=0.0_cp

   end subroutine initialize_updateZ
!-------------------------------------------------------------------------------
   subroutine finalize_updateZ

      integer, pointer :: nLMBs2(:)
      integer :: ll

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
      deallocate( rhs1, Dif )

   end subroutine finalize_updateZ
!-------------------------------------------------------------------------------
   subroutine updateZ(z,dz,dzdt,time,omega_ma,omega_ic,lorentz_torque_ma,     &
              &       lorentz_torque_ic,tscheme,lRmsNext)
      !
      !  updates the toroidal potential z and its radial derivatives
      !  adds explicit part to time derivatives of z
      !

      !-- Input/output of scalar fields:
      class(type_tscheme), intent(in) :: tscheme
      complex(cp), intent(inout) :: z(llm:ulm,n_r_max)        ! Toroidal velocity potential z
      type(type_tarray), intent(inout) :: dzdt
      !type(type_tscalar), intent(inout) :: domega_ic_dt
      !type(type_tscalar), intent(inout) :: domega_ma_dt
      real(cp),    intent(in) :: lorentz_torque_ma            ! Lorentz torque (for OC rotation)
      real(cp),    intent(in) :: lorentz_torque_ic            ! Lorentz torque (for IC rotation)

      !-- Input of other variables:
      real(cp),    intent(in) :: time       ! Current time
      logical,     intent(in) :: lRmsNext   ! Logical for storing update if (l_RMS.and.l_logNext)

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
      complex(cp) :: z10(n_r_max),z11(n_r_max) ! toroidal flow scalar components
      real(cp) :: angular_moment(3)   ! total angular momentum
      real(cp) :: angular_moment_oc(3)! x,y,z component of outer core angular mom.
      real(cp) :: angular_moment_ic(3)! x,y,z component of inner core angular mom.
      real(cp) :: angular_moment_ma(3)! x,y,z component of mantle angular mom.
      complex(cp) :: corr_l1m0      ! correction factor for z(l=1,m=0)
      complex(cp) :: corr_l1m1      ! correction factor for z(l=1,m=1)
      real(cp) :: r_E_2             ! =r**2
      real(cp) :: nomi              ! nominator for Z10 AM correction
      real(cp) :: prec_fac
      real(cp) :: dom_ma, dom_ic, lo_ma, lo_ic
      integer :: l1m0,l1m1          ! position of (l=1,m=0) and (l=1,m=1) in lm.
      integer :: i                  ! counter
      logical :: l10
      integer :: nLMB
      real(cp) :: ddzASL_loc(l_max+1,n_r_max)

      integer, pointer :: nLMBs2(:),lm2l(:),lm2m(:)
      integer, pointer :: sizeLMB2(:,:),lm2(:,:)
      integer, pointer :: lm22lm(:,:,:),lm22l(:,:,:),lm22m(:,:,:)

      integer :: start_lm, stop_lm
      integer :: nChunks,iChunk,lmB0,size_of_last_chunk,threadid

      if ( l_precession ) then
         prec_fac=sqrt(8.0_cp*pi*third)*po*oek*oek*sin(prec_angle)*tscheme%dt(1)
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
      lm2m(1:lm_max) => lo_map%lm2m


      nLMB = 1+rank
      lmStart_00 =max(2,llm)
      l1m0       =lm2(1,0)

      !-- Calculation of the implicit part
      call get_tor_rhs_imp(z, dz, dzdt%old(:,:,tscheme%istage),    &
           &               dzdt%impl(:,:,tscheme%istage),          &
           &               domega_ma_dt%old(tscheme%istage),       &
           &               domega_ic_dt%old(tscheme%istage),       &
           &               domega_ma_dt%impl(tscheme%istage),      &
           &               domega_ic_dt%impl(tscheme%istage),      &
           &               tscheme,                                &
           &               tscheme%l_imp_calc_rhs(tscheme%istage), &
           &               lRmsNext)

      if ( ktopv == 2 .and. l_rot_ma  .and. (.not. l_SRMA)) then
         domega_ma_dt%expl(tscheme%istage)=LFfac*c_lorentz_ma*lorentz_torque_ma
         call tscheme%set_imex_rhs_scalar(dom_ma, domega_ma_dt)
      end if

      if ( kbotv == 2 .and. l_rot_ic .and. (.not. l_SRIC)) then
         domega_ic_dt%expl(tscheme%istage)=LFfac*c_lorentz_ic*lorentz_torque_ic
         call tscheme%set_imex_rhs_scalar(dom_ic, domega_ic_dt)
      end if

      !-- Now assemble the right hand side and store it in work_LMloc
      call tscheme%set_imex_rhs(work_LMloc, dzdt, llm, ulm, n_r_max)


      !$OMP parallel default(shared) private(start_lm, stop_lm)

      call solve_counter%start_count()
      l10=.false.
      !$OMP SINGLE
      do nLMB2=1,nLMBs2(nLMB)
         !$OMP TASK default(shared) &
         !$OMP firstprivate(nLMB2) &
         !$OMP private(lmB,lm,lm1,l1,m1,n_r_out,nR) &
         !$OMP private(nChunks,size_of_last_chunk,iChunk) &
         !$OMP private(tOmega_ma1,tOmega_ma2,threadid) &
         !$OMP private(tOmega_ic1,tOmega_ic2,rhs_sum)
         nChunks = (sizeLMB2(nLMB2,nLMB)+chunksize-1)/chunksize
         size_of_last_chunk=chunksize+(sizeLMB2(nLMB2,nLMB)-nChunks*chunksize)

         ! This task treats one l given by l1
         l1=lm22l(1,nLMB2,nLMB)
         !write(*,"(3(A,I3),A)") "Launching task for nLMB2=", &
         !     &   nLMB2," (l=",l1,") and scheduling ",nChunks," subtasks."

         if ( l1 /= 0 ) then
            if ( .not. lZmat(l1) ) then
#ifdef WITH_PRECOND_Z
               call get_zMat(tscheme,l1,hdif_V(st_map%lm2(l1,0)), &
                    &        zMat(nLMB2),zMat_fac(:,nLMB2))
#else
               call get_zMat(tscheme,l1,hdif_V(st_map%lm2(l1,0)), &
                    &        zMat(nLMB2))
#endif
               lZmat(l1)=.true.
            !write(*,"(A,I3,A,2ES20.12)") "zMat(",l1,") = ",SUM(zMat(:,:,l1))
            end if
         end if

         do iChunk=1,nChunks
            !$OMP TASK default(shared) &
            !$OMP firstprivate(iChunk) &
            !$OMP private(lmB0,lmB,lm,lm1,m1,nR,n_r_out) &
            !$OMP private(tOmega_ma1,tOmega_ma2) &
            !$OMP private(tOmega_ic1,tOmega_ic2,rhs_sum) &
            !$OMP private(threadid)
#ifdef WITHOMP
            threadid = omp_get_thread_num()
#else
            threadid = 0
#endif

            lmB0=(iChunk-1)*chunksize
            lmB=lmB0

            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               !do lm=1,sizeLMB2(nLMB2,nLMB)
               lm1=lm22lm(lm,nLMB2,nLMB)
               !l1 =lm22l(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)
               if ( lm1 == l1m0 ) l10= .true.

               if ( l_z10mat .and. lm1 == l1m0 ) then
                  !write(*,"(A,3I3)") "l_z10mat and lm1=",lm1,l1,m1
                  !PERFON('upZ_z10')
                  !----- Special treatment of z10 component if ic or mantle
                  !      are allowed to rotate about z-axis (l_z10mat=.true.) and
                  !      we use no slip boundary condition (ktopv=1,kbotv=1):
                  !      Lorentz torque is the explicit part of this time integration
                  !      at the boundaries!
                  !      Note: no angular momentum correction necessary for this case !
                  if ( .not. lZ10mat ) then
#ifdef WITH_PRECOND_Z10
                     call get_z10Mat(tscheme,l1,                              &
                          &          hdif_V(st_map%lm2(lm2l(lm1),lm2m(lm1))), &
                          &          z10Mat,z10Mat_fac)
#else
                     call get_z10Mat(tscheme,l1,hdif_V(                &
                          &          st_map%lm2(lm2l(lm1),lm2m(lm1))), &
                          &          z10Mat)
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
                  !PERFON('upZ_ln0')
                  lmB=lmB+1

                  rhs1(1,lmB,threadid)      =0.0_cp
                  rhs1(n_r_max,lmB,threadid)=0.0_cp

                  if (amp_RiIc /= 0.0_cp) then
                     
                     if (l1 == (m_RiIc + RiSymmIc) .and. m1 == m_RiIc) then
                        rhs1(n_r_max,lmB,threadid) = cmplx(      &
                        &         amp_RiIc*cos(omega_RiIc*time), &
                        &         amp_RiIc*sin(omega_RiIc*time), &
                        &                                  kind=cp)
                     end if
                  end if

                  if (amp_RiMa /= 0.0_cp) then
                     if (l1 == (m_RiMa + RiSymmMa) .and. m1 == m_RiMa) then
                        rhs1(1,lmB,threadid) = cmplx(            &
                        &        amp_RiMa*cos(omega_RiMa*time),  &
                        &        amp_RiMa*sin(omega_RiMa*time),  &
                        &                            kind=cp)
                     end if
                  end if

                  do nR=2,n_r_max-1
                     rhs1(nR,lmB,threadid)=work_LMloc(lm1,nR)
                     if ( l_precession .and. l1 == 1 .and. m1 == 1 ) then
                        rhs1(nR,lmB,threadid)=rhs1(nR,lmB,threadid)+prec_fac* &
                        &    cmplx(sin(oek*time),-cos(oek*time),kind=cp)
                     end if
                  end do

#ifdef WITH_PRECOND_Z
                  do nR=1,n_r_max
                     rhs1(nR,lmB,threadid)=zMat_fac(nR,nLMB2)*rhs1(nR,lmB,threadid)
                  end do
#endif
                  !PERFOFF
               end if
            end do

            !PERFON('upZ_sol')
            if ( lmB > lmB0 ) then
               call zMat(nLMB2)%solve(rhs1(:,lmB0+1:lmB,threadid),lmB-lmB0)
            end if
            !PERFOFF

            lmB=lmB0
            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               !do lm=1,sizeLMB2(nLMB2,nLMB)
               lm1=lm22lm(lm,nLMB2,nLMB)
               !l1 =lm22l(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)

               if ( l_z10mat .and. lm1 == l1m0 ) then
                  do n_r_out=1,rscheme_oc%n_max
                     z(lm1,n_r_out)=real(rhs(n_r_out))
                  end do
               else if ( l1 /= 0 ) then
                  lmB=lmB+1
                  if ( m1 > 0 ) then
                     do n_r_out=1,rscheme_oc%n_max
                        z(lm1,n_r_out)=rhs1(n_r_out,lmB,threadid)
                     end do
                  else
                     do n_r_out=1,rscheme_oc%n_max
                        z(lm1,n_r_out)=cmplx( &
                        &  real(rhs1(n_r_out,lmB,threadid)),0.0_cp,kind=cp)
                     end do
                  end if
               end if
            end do
            !PERFOFF
            !$OMP END TASK
         end do
         !$OMP END TASK
      end do       ! end of loop over lm blocks
      !$OMP END SINGLE
      !$omp single
      call solve_counter%stop_count(l_increment=.false.)
      !$omp end single

      !-- set cheb modes > rscheme_oc%n_max to zero (dealiazing)
      !$omp single
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         do lm1=llm,ulm
            z(lm1,n_r_out)=zero
         end do
      end do
      !$omp end single

      !PERFON('upZ_drv')
      !$omp single
      call dct_counter%start_count()
      !$omp end single
      start_lm=lmStart_00; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)
      call get_ddr(z, dz, work_LMloc, ulm-llm+1, start_lm-llm+1,     &
           &       stop_lm-llm+1,n_r_max, rscheme_oc, l_dct_in=.false.)
      call rscheme_oc%costf1(z,ulm-llm+1,start_lm-llm+1,stop_lm-llm+1)
      !$omp barrier
      !PERFOFF
      !$omp single
      call dct_counter%stop_count(l_increment=.false.)
      !$omp end single

      !$omp end parallel

      !PERFON('upZ_icma')
      !--- Update of inner core and mantle rotation:
      if ( l10 ) then
         if ( l_rot_ma .and. .not. l_SRMA ) then
            if ( ktopv == 1 ) then  ! free slip, explicit time stepping of omega !
               call get_rot_rates(omega_ma, lorentz_torque_ma, c_moi_ma,    &
                    &             lorentz_torque_ma_dt%old(tscheme%istage), &
                    &             lorentz_torque_ma_dt%expl(tscheme%istage) )
               call tscheme%set_imex_rhs_scalar(lo_ma, lorentz_torque_ma_dt)
               omega_ma=lo_ma
            else if ( ktopv == 2 ) then ! no slip, omega given by z10
               omega_ma=c_z10_omega_ma*real(z(l1m0,n_r_cmb))
            end if
            omega_ma1=omega_ma
         end if
         if ( l_rot_ic .and. .not. l_SRIC ) then
            if ( kbotv == 1 ) then  ! free slip, explicit time stepping of omega !
               call get_rot_rates(omega_ic, lorentz_torque_ic, c_moi_ic,    &
                    &             lorentz_torque_ic_dt%old(tscheme%istage), &
                    &             lorentz_torque_ic_dt%expl(tscheme%istage) )
               call tscheme%set_imex_rhs_scalar(lo_ic, lorentz_torque_ic_dt)
               omega_ic=lo_ic
            else if ( kbotv == 2 ) then ! no slip, omega given by z10
               omega_ic=c_z10_omega_ic*real(z(l1m0,n_r_icb))
            end if
            omega_ic1=omega_ic
         end if
      end if  ! l=1,m=0 contained in block ?
      !PERFOFF
      !PERFON('upZ_ang')
      !--- We correct so that the angular moment about axis in the equatorial plane
      !    vanish and the angular moment about the (planetary) rotation axis
      !    is kept constant.
      l1m1=lm2(1,1)
      if ( l_correct_AMz .and.  l1m0 > 0 .and. &
         & lmStart_00 <= l1m0 .and. ulm >= l1m0 ) then

         do nR=1,n_r_max
            z10(nR)=z(l1m0,nR)
         end do
         call get_angular_moment(z10,z11,omega_ic,omega_ma,          &
              &                  angular_moment_oc,angular_moment_ic,&
              &                  angular_moment_ma)
         do i=1,3
            angular_moment(i)=angular_moment_oc(i) + angular_moment_ic(i) + &
            &                 angular_moment_ma(i)
         end do
         if ( ( ktopv == 2 .and. l_rot_ma ) .and. &
              ( kbotv == 2 .and. l_rot_ic ) ) then
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
         corr_l1m0=cmplx(angular_moment(3)-AMstart,0.0_cp,kind=cp)/nomi

         !-------- Correct z(2,nR) and z(l_max+2,nR) plus the respective
         !         derivatives:
         !$OMP PARALLEL do default(shared) PRIVATE(r_E_2,nR)
         do nR=1,n_r_max
            r_E_2=r(nR)*r(nR)
            z(l1m0,nR)  =z(l1m0,nR)  - rho0(nR)*r_E_2*corr_l1m0
            dz(l1m0,nR) =dz(l1m0,nR) - rho0(nR)*( &
            &            two*r(nR)+r_E_2*beta(nR))*corr_l1m0
            work_LMloc(l1m0,nR)=work_LMloc(l1m0,nR)-rho0(nR)*( &
            &              two+four*beta(nR)*r(nR) +           &
            &              dbeta(nR)*r_E_2 +                   &
            &              beta(nR)*beta(nR)*r_E_2 )*corr_l1m0
         end do
         !$OMP END PARALLEL DO
         if ( ktopv == 2 .and. l_rot_ma ) &
         &    omega_ma=c_z10_omega_ma*real(z(l1m0,n_r_cmb))
         if ( kbotv == 2 .and. l_rot_ic ) &
         &    omega_ic=c_z10_omega_ic*real(z(l1m0,n_r_icb))
         omega_ic1=omega_ic
         omega_ma1=omega_ma

      end if ! l=1,m=0 contained in lm-block ?

      if ( l_correct_AMe .and.  l1m1 > 0 .and. &
           lmStart_00 <= l1m1 .and. ulm >= l1m1 ) then

         do nR=1,n_r_max
            z11(nR)=z(l1m1,nR)
         end do
         call get_angular_moment(z10,z11,omega_ic,omega_ma,          &
              &                  angular_moment_oc,angular_moment_ic,&
              &                  angular_moment_ma)
         do i=1,3
            angular_moment(i)=angular_moment_oc(i) + angular_moment_ic(i) + &
            &                 angular_moment_ma(i)
         end do
         corr_l1m1=cmplx(angular_moment(1),-angular_moment(2),kind=cp) / &
         &         (two*y11_norm*c_moi_oc)

         !-------- Correct z(2,nR) and z(l_max+2,nR) plus the respective
         !         derivatives:
         !$OMP PARALLEL do default(shared) private(nR,r_E_2)
         do nR=1,n_r_max
            r_E_2=r(nR)*r(nR)
            z(l1m1,nR)  =z(l1m1,nR)  -  rho0(nR)*r_E_2*corr_l1m1
            dz(l1m1,nR) =dz(l1m1,nR) -  rho0(nR)*( &
            &            two*r(nR)+r_E_2*beta(nR))*corr_l1m1
            work_LMloc(l1m1,nR)=work_LMloc(l1m1,nR)-rho0(nR)*(    &
            &             two+four*beta(nR)*r(nR) +               &
            &                        dbeta(nR)*r_E_2 +            &
            &                beta(nR)*beta(nR)*r_E_2 )*corr_l1m1
         end do
         !$OMP END PARALLEL DO
      end if ! l=1,m=1 contained in lm-block ?

      !-- Calculate explicit time step part:
      !$OMP PARALLEL default(shared) private(nR,lm1,Dif)

      !--- Note: from ddz=work_LMloc only the axisymmetric contributions are needed
      !    beyond this point for the TO calculation.
      !    Parallization note: Very likely, all axisymmetric modes m=0 are
      !    located on the first processor #0.
      if ( l_TO ) then
         !$OMP do private(nR,lm1,l1,m1)
         do nR=1,n_r_max
            ddzASL_loc(:,nR)=0.0_cp
            do lm1=lmStart_00,ulm
               l1=lm2l(lm1)
               m1=lm2m(lm1)
               if ( m1 == 0 ) ddzASL_loc(l1+1,nR)=real(work_LMloc(lm1,nR))
            end do
         end do
         !$OMP end do
      end if
      !$OMP END PARALLEL

      if ( l_TO ) then
         do nR=1,n_r_max
#ifdef WITH_MPI
            call MPI_Allreduce(ddzASL_loc(:,nR), ddzASL(:,nR), l_max+1, &
                 &             MPI_DEF_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
            ddzASL(:,nR)=ddzASL_loc(:,nR)
#endif
         end do
      end if

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dzdt, llm, ulm, n_r_max)
      call tscheme%rotate_imex_scalar(domega_ma_dt)
      call tscheme%rotate_imex_scalar(domega_ic_dt)
      call tscheme%rotate_imex_scalar(lorentz_torque_ma_dt)
      call tscheme%rotate_imex_scalar(lorentz_torque_ic_dt)

   end subroutine updateZ
!-----------------------------------------------------------------------------
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
      !        Note opposite sign of viscous torques (-dz+2 z /r )
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
            &       tscheme%wimp_lin(1)*(                       &
            &                 two*or1(1)*rscheme_oc%rMat(1,:) - &
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
               &       tscheme%wimp_lin(1)*(                             &
               &           two*or1(n_r_max)*rscheme_oc%rMat(n_r_max,:) - &
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
!-----------------------------------------------------------------------
   subroutine get_tor_rhs_imp(z, dz, z_last, dz_imp_last, domega_ma_old,     &
              &               domega_ic_old, domega_ma_last, domega_ic_last, &
              &               tscheme, l_calc_lin_rhs, lRmsNext)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      complex(cp),         intent(in) :: z(llm:ulm,n_r_max)
      logical,             intent(in) :: l_calc_lin_rhs
      logical,             intent(in) :: lRmsNext

      !-- Output variable
      complex(cp), intent(out) :: dz(llm:ulm,n_r_max)
      complex(cp), intent(out) :: z_last(llm:ulm,n_r_max)
      complex(cp), intent(out) :: dz_imp_last(llm:ulm,n_r_max)
      real(cp),    intent(out) :: domega_ma_old
      real(cp),    intent(out) :: domega_ic_old
      real(cp),    intent(out) :: domega_ma_last
      real(cp),    intent(out) :: domega_ic_last

      !-- Local variables
      integer :: n_r, lm, start_lm, stop_lm, n_r_bot, n_r_top
      integer :: lmStart_00, l1, m1, l1m0
      integer, pointer :: lm2l(:),lm2m(:), lm2(:,:)

      lm2(0:, 0:) => lo_map%lm2
      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m
      lmStart_00 =max(2,llm)

      !$omp parallel default(shared)  private(start_lm, stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)

      !$omp do private(n_r,lm,l1,m1) collapse(2)
      do n_r=1,n_r_max
         do lm=llm,ulm
            l1 = lm2l(lm)
            m1 = lm2m(lm)
            z_last(lm,n_r)=dLh(st_map%lm2(l1,m1))*or2(n_r)*z(lm,n_r)
         end do
      end do
      !$omp end do

      if ( l_calc_lin_rhs .or. (tscheme%istage==tscheme%nstages .and. lRmsNext)) then
         call get_ddr( z, dz, work_LMloc, ulm-llm+1,                  &
              &        start_lm-llm+1, stop_lm-llm+1, n_r_max, rscheme_oc)
         !$omp barrier

         if ( lRmsNext ) then
            n_r_top=n_r_cmb
            n_r_bot=n_r_icb
         else
            n_r_top=n_r_cmb+1
            n_r_bot=n_r_icb-1
         end if

         !$omp do default(shared) private(n_r,lm,Dif) collapse(2)
         do n_r=n_r_top,n_r_bot
            do lm=lmStart_00,ulm
               Dif(lm)=hdif_V(st_map%lm2(lm2l(lm),lm2m(lm)))*                     &
               &        dLh(st_map%lm2(lm2l(lm),lm2m(lm)))*or2(n_r)*visc(n_r)*    &
               &      ( work_LMloc(lm,n_r)   +(dLvisc(n_r)-beta(n_r)) *dz(lm,n_r) &
               &        -( dLvisc(n_r)*beta(n_r)+two*dLvisc(n_r)*or1(n_r)         &
               &           + dLh(st_map%lm2(lm2l(lm),lm2m(lm)))*or2(n_r)          &
               &           + dbeta(n_r)+ two*beta(n_r)*or1(n_r) ) * z(lm,n_r) )

               dz_imp_last(lm,n_r)=Dif(lm)
            end do
            if ( lRmsNext .and. tscheme%istage==tscheme%nstages ) then
               call hInt2Tor(Dif,llm,ulm,n_r,lmStart_00,ulm, &
                    &        DifTor2hInt(:,n_r),lo_map)
            end if
         end do
         !$omp end do

      end if
      !$omp end parallel

      l1m0 = lm2(1,0)
      if ( ( llm <= l1m0 .and. ulm >= l1m0 ) .and. l_z10mat ) then
         !----- NOTE opposite sign of viscous torque on ICB and CMB:
         if ( .not. l_SRMA .and. ktopv == 2 .and. l_rot_ma ) then
            domega_ma_last=two*or1(1)*real(z(l1m0,1))-real(dz(l1m0,1))
            domega_ma_old =c_dt_z10_ma*real(z(l1m0,1))
         end if
         if ( .not. l_SRIC .and. kbotv == 2 .and. l_rot_ic ) then
            domega_ic_last=-two*or1(n_r_max)*real(z(l1m0,n_r_max))+ &
            &                  real(dz(l1m0,n_r_max))
            domega_ic_old =c_dt_z10_ic*real(z(l1m0,n_r_max))
         end if
      end if

   end subroutine get_tor_rhs_imp
!-----------------------------------------------------------------------
   subroutine get_rot_rates(omega, lorentz_torque, c_moi, &
              &             domega_old, domega_last)

      real(cp), intent(in) :: omega
      real(cp), intent(in) :: lorentz_torque
      real(cp), intent(in) :: c_moi

      real(cp), intent(out) :: domega_old
      real(cp), intent(out) :: domega_last

      domega_old = omega
      domega_last = LFfac/c_moi * lorentz_torque

   end subroutine get_rot_rates
!-----------------------------------------------------------------------
#ifdef WITH_PRECOND_Z
   subroutine get_zMat(tscheme,l,hdif,zMat,zMat_fac)
#else
   subroutine get_zMat(tscheme,l,hdif,zMat)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matricies
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

#ifdef MATRIX_CHECK
      integer :: i,j
      real(cp) :: rcond
      integer ::ipiv(n_r_max),iwork(n_r_max)
      real(cp) :: work(4*n_r_max),anorm,linesum
      real(cp) :: temp_Mat(n_r_max,n_r_max)
      integer, save :: counter=0
      integer :: filehandle
      character(len=100) :: filename
#endif

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
!-----------------------------------------------------------------------------
end module updateZ_mod
