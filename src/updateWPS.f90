module updateWPS_mod
   !
   ! This module handles the time advance of the poloidal potential w,
   ! the pressure p and the entropy s in one single matrix per degree.
   ! It contains the computation of the implicit terms and the linear
   ! solves.
   !

   use omp_lib
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: lm_max, n_r_max, l_max, m_min, l_axi
   use radial_data, only: n_r_cmb,n_r_icb, nRstart, nRstop
   use radial_functions, only: r_cmb, or1, or2, rho0, rgrav, r, visc, dLvisc,    &
       &                       rscheme_oc, beta, dbeta, dLkappa, dLtemp0, &
       &                       ddLtemp0, alpha0, dLalpha0, ddLalpha0,     &
       &                       ogrun, kappa, orho1, dentropy0, temp0, l_R
   use physical_parameters, only: kbotv, ktopv, ktops, kbots, ra, opr, &
       &                          ViscHeatFac, ThExpNb, BuoFac,        &
       &                          CorFac, ktopp
   use num_param, only: dct_counter, solve_counter
   use init_fields, only: tops, bots
   use special, only: ellip_fac_cmb, ellip_fac_icb, tide_fac20, tide_fac22p, &
       &              tide_fac22n, omega_tide, amp_tide, l_vr_cmb, l_vr_icb, &
       &              l2_h2_ratio
   use init_fields, only: omegaOsz_ma1, tShift_ma1, omegaOsz_ic1, tShift_ic1
   use blocking, only: lo_sub_map, lo_map, llm, ulm, st_map
   use horizontal_data, only: hdif_V, hdif_S
   use logic, only: l_temperature_diff, l_RMS, l_full_sphere
   use RMS, only: DifPol2hInt, DifPolLMr
   use algebra, only: prepare_mat, solve_mat
   use communications, only: get_global_sum
   use parallel_mod, only: chunksize, rank, n_procs, get_openmp_blocks
#ifdef WITH_MPI
   use mpimod
#endif
   use radial_der, only: get_dddr, get_ddr, get_dr, get_dr_Rloc
   use constants, only: zero, one, two, three, four, third, half, pi, osq4pi
   use fields, only: work_LMloc, s_LMloc, ds_LMloc, p_LMloc, dp_LMloc
   use useful, only: abortRun, cc2real
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray

   implicit none

   private

   !-- Input of recycled work arrays:
   complex(cp), allocatable :: workB(:,:), workC(:,:)
   real(cp), allocatable :: rhs1(:,:,:)
   real(cp), allocatable :: ps0Mat(:,:), ps0Mat_fac(:,:)
   integer, allocatable :: ps0Pivot(:)
   real(cp), allocatable :: wpsMat(:,:,:)
   integer, allocatable :: wpsPivot(:,:)
   real(cp), allocatable :: wpsMat_fac(:,:,:)
   real(cp) :: Cor00_fac
   logical, public, allocatable :: lWPSmat(:)

   !-- Boundary values due to the prescribed motion of the boundaries (standard
   !-- lm order, see update_wps_bc)
   complex(cp), allocatable :: wbc_cmb(:), wbc_icb(:)   ! W at the boundaries
   complex(cp), allocatable :: wbc2_cmb(:), wbc2_icb(:) ! r.h.s. of the second b.c. on W
   complex(cp), allocatable :: sbc_cmb(:) ! Tidal perturbation of the outer entropy/temperature b.c.
   logical :: l_sbc_tide ! Tidal perturbation of a fixed entropy/temperature at r_cmb (see update_wps_sbc)

   integer :: maxThreads

   public :: initialize_updateWPS, finalize_updateWPS, updateWPS, finish_exp_smat,&
   &         get_single_rhs_imp, assemble_single, finish_exp_smat_Rdist, &
   &         update_wps_bc

contains

   subroutine initialize_updateWPS

      integer, pointer :: nLMBs2(:)

      nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

      allocate( ps0Mat(2*n_r_max,2*n_r_max) )
      allocate( ps0Mat_fac(2*n_r_max,2) )
      allocate( ps0Pivot(2*n_r_max) )
      bytes_allocated = bytes_allocated+(4*n_r_max+2)*n_r_max*SIZEOF_DEF_REAL &
      &                 +2*n_r_max*SIZEOF_INTEGER
      allocate( wpsMat(3*n_r_max,3*n_r_max,nLMBs2(1+rank)) )
      allocate(wpsMat_fac(3*n_r_max,2,nLMBs2(1+rank)))
      allocate ( wpsPivot(3*n_r_max,nLMBs2(1+rank)) )
      bytes_allocated = bytes_allocated+(9*n_r_max*nLMBs2(1+rank)+6*n_r_max* &
      &                 nLMBs2(1+rank))*SIZEOF_DEF_REAL+3*n_r_max*           &
      &                 nLMBs2(1+rank)*SIZEOF_INTEGER
      allocate( lWPSmat(0:l_max) )
      bytes_allocated = bytes_allocated+(l_max+1)*SIZEOF_LOGICAL

      allocate( workB(llm:ulm,n_r_max) )
      allocate( workC(llm:ulm,n_r_max) )
      bytes_allocated = bytes_allocated+2*(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX

#ifdef WITHOMP
      maxThreads=omp_get_max_threads()
#else
      maxThreads=1
#endif

      allocate( rhs1(3*n_r_max,2*lo_sub_map%sizeLMB2max,0:maxThreads-1) )
      bytes_allocated=bytes_allocated+2*n_r_max*maxThreads* &
                      lo_sub_map%sizeLMB2max*SIZEOF_DEF_COMPLEX

      Cor00_fac=four/sqrt(three)


      allocate( wbc_cmb(lm_max), wbc_icb(lm_max), wbc2_cmb(lm_max), wbc2_icb(lm_max) )
      allocate( sbc_cmb(lm_max) )
      bytes_allocated=bytes_allocated+5*lm_max*SIZEOF_DEF_COMPLEX
      wbc_cmb(:) =zero
      wbc_icb(:) =zero
      wbc2_cmb(:)=zero
      wbc2_icb(:)=zero
      sbc_cmb(:) =zero

      !-- The tidal correction only applies to a fixed value (ktops=1,3): for a
      !-- fixed flux it is O(delta/H) and omitted (see update_wps_sbc)
      l_sbc_tide=amp_tide /= 0.0_cp .and. mod(ktops,2) == 1

   end subroutine initialize_updateWPS
!-----------------------------------------------------------------------------
   subroutine finalize_updateWPS

      deallocate( ps0Mat, ps0Mat_fac, ps0Pivot )
      deallocate( wpsMat, wpsMat_fac, wpsPivot, lWPSmat )
      deallocate( workB, workC, rhs1)

      deallocate( wbc_cmb, wbc_icb, wbc2_cmb, wbc2_icb, sbc_cmb )

   end subroutine finalize_updateWPS
!-----------------------------------------------------------------------------
   subroutine updateWPS(w,dw,ddw,z10,dwdt,p,dp,dpdt,s,ds,dsdt,tscheme,lRmsNext)
      !
      !  updates the poloidal velocity potential w, the pressure p, the entropy
      !  and their radial derivatives.
      !

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme
      real(cp),            intent(in) :: z10(n_r_max)
      logical,             intent(in) :: lRmsNext

      !-- Output variables
      type(type_tarray), intent(inout) :: dwdt
      type(type_tarray), intent(inout) :: dpdt
      type(type_tarray), intent(inout) :: dsdt
      complex(cp),       intent(inout) :: w(llm:ulm,n_r_max)
      complex(cp),       intent(inout) :: dw(llm:ulm,n_r_max)
      complex(cp),       intent(out)   :: ddw(llm:ulm,n_r_max)
      complex(cp),       intent(inout) :: p(llm:ulm,n_r_max)
      complex(cp),       intent(out)   :: dp(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: s(llm:ulm,n_r_max)
      complex(cp),       intent(inout) :: ds(llm:ulm,n_r_max)

      !-- Local variables:
      integer :: l1,m1          ! degree and order
      integer :: lm1,lm         ! position of (l,m) in array
      integer :: lms            ! position of (l,m) in the standard order
      integer :: nLMB2, nLMB
      integer :: nR             ! counts radial grid points
      integer :: n_r_out        ! counts cheb modes
      real(cp) :: rhs(2*n_r_max)  ! real RHS for l=m=0

      integer, pointer :: nLMBs2(:),lm2l(:),lm2m(:)
      integer, pointer :: sizeLMB2(:,:),lm2(:,:)
      integer, pointer :: lm22lm(:,:,:),lm22l(:,:,:),lm22m(:,:,:)

      integer :: nChunks,iChunk,lmB0,size_of_last_chunk,threadid

      nLMBs2(1:n_procs) => lo_sub_map%nLMBs2
      sizeLMB2(1:,1:) => lo_sub_map%sizeLMB2
      lm22lm(1:,1:,1:) => lo_sub_map%lm22lm
      lm22l(1:,1:,1:) => lo_sub_map%lm22l
      lm22m(1:,1:,1:) => lo_sub_map%lm22m
      lm2(0:,0:) => lo_map%lm2
      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      nLMB=1+rank

      !-- Now assemble the right hand side and store it in work_LMloc, dp and ds
      call tscheme%set_imex_rhs(work_LMloc, dwdt)
      call tscheme%set_imex_rhs(dp, dpdt)
      call tscheme%set_imex_rhs(ds, dsdt)

      !$omp parallel default(shared)

      !$omp single
      call solve_counter%start_count()
      !$omp end single
      !$omp single
      ! each of the nLMBs2(nLMB) subblocks have one l value
      do nLMB2=1,nLMBs2(nLMB)

         !$omp task default(shared) &
         !$omp firstprivate(nLMB2) &
         !$omp private(lm,lm1,l1,m1,iChunk,nChunks,size_of_last_chunk,threadid)

         ! determine the number of chunks of m
         ! total number for l1 is sizeLMB2(nLMB2,nLMB)
         ! chunksize is given
         nChunks = (sizeLMB2(nLMB2,nLMB)+chunksize-1)/chunksize
         size_of_last_chunk = chunksize + (sizeLMB2(nLMB2,nLMB)-nChunks*chunksize)

         l1=lm22l(1,nLMB2,nLMB)
         if ( l1 == 0 ) then ! Spherically-symmetric, l=0 mode
            lm1 = lm2(0,0)

            if ( .not. lWPSmat(l1) ) then
               call get_ps0Mat(tscheme,ps0Mat,ps0Pivot,ps0Mat_fac)
               lWPSmat(l1)=.true.
            end if

            do nR=1,n_r_max
               rhs(nR)        =real(ds(lm1,nR))
               rhs(nR+n_r_max)=real(dwdt%expl(lm1,nR,tscheme%istage))+&
               &               Cor00_fac*CorFac*or1(nR)*z10(nR)
            end do
            rhs(1)        =real(tops(0,0))
            rhs(n_r_max)  =real(bots(0,0))
            rhs(n_r_max+1)=0.0_cp

            do nR=1,2*n_r_max
               rhs(nR)=rhs(nR)*ps0Mat_fac(nR,1)
            end do

            call solve_mat(ps0Mat,2*n_r_max,2*n_r_max,ps0Pivot,rhs)

            do nR=1,2*n_r_max
               rhs(nR)=rhs(nR)*ps0Mat_fac(nR,2)
            end do

            do n_r_out=1,rscheme_oc%n_max
               s(lm1,n_r_out)=rhs(n_r_out)
               p(lm1,n_r_out)=rhs(n_r_out+n_r_max)
               w(lm1,n_r_out)=zero
            end do

         else ! l1 /= 0

            if ( .not. lWPSmat(l1) ) then
               call get_wpsMat(tscheme, l1, hdif_V(l1), hdif_S(l1), wpsMat(:,:,nLMB2), &
                    &          wpsPivot(:,nLMB2), wpsMat_fac(:,:,nLMB2))
               lWPSmat(l1)=.true.
            end if

            do iChunk=1,nChunks
               !$omp task default(shared) &
               !$omp firstprivate(iChunk) &
               !$omp private(lmB0,lm,lm1,m1,lms,nR,n_r_out) &
               !$omp private(threadid)

#ifdef WITHOMP
               threadid = omp_get_thread_num()
#else
               threadid = 0
#endif
               lmB0=(iChunk-1)*chunksize
               do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
                  lm1=lm22lm(lm,nLMB2,nLMB)
                  m1 =lm22m(lm,nLMB2,nLMB)

                  lms=st_map%lm2(l1,m1)
                  rhs1(1,2*lm-1,threadid)          = real(wbc_cmb(lms))
                  rhs1(1,2*lm,threadid)            =aimag(wbc_cmb(lms))
                  rhs1(n_r_max,2*lm-1,threadid)    = real(wbc_icb(lms))
                  rhs1(n_r_max,2*lm,threadid)      =aimag(wbc_icb(lms))
                  rhs1(n_r_max+1,2*lm-1,threadid)  = real(wbc2_cmb(lms))
                  rhs1(n_r_max+1,2*lm,threadid)    =aimag(wbc2_cmb(lms))
                  rhs1(2*n_r_max,2*lm-1,threadid)  = real(wbc2_icb(lms))
                  rhs1(2*n_r_max,2*lm,threadid)    =aimag(wbc2_icb(lms))
                  rhs1(2*n_r_max+1,2*lm-1,threadid)= real(tops(l1,m1)+sbc_cmb(lms))
                  rhs1(2*n_r_max+1,2*lm,threadid)  =aimag(tops(l1,m1)+sbc_cmb(lms))
                  rhs1(3*n_r_max,2*lm-1,threadid)  = real(bots(l1,m1))
                  rhs1(3*n_r_max,2*lm,threadid)    =aimag(bots(l1,m1))

                  do nR=2,n_r_max-1
                     !-- dp and ds used as work arrays here
                     rhs1(nR,2*lm-1,threadid)          = real(work_LMloc(lm1,nR))
                     rhs1(nR,2*lm,threadid)            =aimag(work_LMloc(lm1,nR))
                     rhs1(nR+n_r_max,2*lm-1,threadid)  = real(dp(lm1,nR))
                     rhs1(nR+n_r_max,2*lm,threadid)    =aimag(dp(lm1,nR))
                     rhs1(nR+2*n_r_max,2*lm-1,threadid)= real(ds(lm1,nR))
                     rhs1(nR+2*n_r_max,2*lm,threadid)  =aimag(ds(lm1,nR))
                  end do

                  rhs1(:,2*lm-1,threadid)=rhs1(:,2*lm-1,threadid)*wpsMat_fac(:,1,nLMB2)
                  rhs1(:,2*lm,threadid)  =rhs1(:,2*lm,threadid)*wpsMat_fac(:,1,nLMB2)
               end do

               call solve_mat(wpsMat(:,:,nLMB2),3*n_r_max,3*n_r_max,                &
                    &         wpsPivot(:,nLMB2),rhs1(:,2*(lmB0+1)-1:2*(lm-1),threadid),&
                    &         2*(lm-1-lmB0))

               do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
                  lm1=lm22lm(lm,nLMB2,nLMB)
                  m1 =lm22m(lm,nLMB2,nLMB)

                  rhs1(:,2*lm-1,threadid)=rhs1(:,2*lm-1,threadid)*wpsMat_fac(:,2,nLMB2)
                  rhs1(:,2*lm,threadid)  =rhs1(:,2*lm,threadid)*wpsMat_fac(:,2,nLMB2)

                  if ( m1 > 0 ) then
                     do n_r_out=1,rscheme_oc%n_max
                        w(lm1,n_r_out)=cmplx(rhs1(n_r_out,2*lm-1,threadid), &
                        &                    rhs1(n_r_out,2*lm,threadid),cp)
                        p(lm1,n_r_out)=cmplx(rhs1(n_r_max+n_r_out,2*lm-1,threadid),&
                        &                    rhs1(n_r_max+n_r_out,2*lm,threadid),cp)
                        s(lm1,n_r_out)=cmplx(rhs1(2*n_r_max+n_r_out,2*lm-1,  &
                        &                    threadid),rhs1(2*n_r_max+n_r_out,&
                        &                    2*lm,threadid),cp)
                     end do
                  else
                     do n_r_out=1,rscheme_oc%n_max
                        w(lm1,n_r_out)= cmplx(rhs1(n_r_out,2*lm-1,threadid), &
                        &                     0.0_cp,kind=cp)
                        p(lm1,n_r_out)= cmplx(rhs1(n_r_max+n_r_out,2*lm-1,threadid), &
                        &                     0.0_cp,kind=cp)
                        s(lm1,n_r_out)= cmplx(rhs1(2*n_r_max+n_r_out,2*lm-1,threadid), &
                        &                     0.0_cp,kind=cp)
                     end do
                  end if
               end do
               !$omp end task
            end do
         end if
         !$omp taskwait
         !$omp end task
      end do   ! end of loop over l1 subblocks
      !$omp end single
      !$omp taskwait
      !$omp single
      call solve_counter%stop_count()
      !$omp end single

      !-- set cheb modes > rscheme_oc%n_max to zero (dealiazing)
      !$omp do private(n_r_out,lm1) collapse(2)
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         do lm1=llm,ulm
            w(lm1,n_r_out)=zero
            p(lm1,n_r_out)=zero
            s(lm1,n_r_out)=zero
         end do
      end do
      !$omp end do

      !$omp end parallel

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dwdt)
      call tscheme%rotate_imex(dpdt)
      call tscheme%rotate_imex(dsdt)

      if ( tscheme%istage == tscheme%nstages ) then
         call get_single_rhs_imp(s, ds, w, dw, ddw, p, dp, dsdt, dwdt, dpdt, &
              &                  tscheme, 1, tscheme%l_imp_calc_rhs(1),      &
              &                  lRmsNext, l_in_cheb_space=.true.)
      else
         call get_single_rhs_imp(s, ds, w, dw, ddw, p, dp, dsdt, dwdt, dpdt, &
              &                  tscheme, tscheme%istage+1,                  &
              &                  tscheme%l_imp_calc_rhs(tscheme%istage+1),   &
              &                  lRmsNext, l_in_cheb_space=.true.)
      end if

   end subroutine updateWPS
!------------------------------------------------------------------------------
   subroutine finish_exp_smat(dVSrLM, ds_exp_last)

      complex(cp), intent(inout) :: dVSrLM(llm:ulm,n_r_max)

      !-- Output variables
      complex(cp), intent(inout) :: ds_exp_last(llm:ulm,n_r_max)

      !-- Local variables
      integer :: n_r, start_lm, stop_lm, lm, l

      !$omp parallel default(shared) private(start_lm, stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)
      call get_dr( dVSrLM, work_LMloc, ulm-llm+1, start_lm-llm+1,  &
           &       stop_lm-llm+1, n_r_max, rscheme_oc, nocopy=.true. )
      !$omp barrier

      !$omp do private(lm,l)
      do n_r=1,n_r_max
         do lm=llm,ulm
            l=lo_map%lm2l(lm)
            if ( l > l_R(n_r) ) cycle
            ds_exp_last(lm,n_r)=orho1(n_r)*(ds_exp_last(lm,n_r)-   &
            &                      or2(n_r)*work_LMloc(lm,n_r))
         end do
      end do
      !$omp end do

      !$omp end parallel

   end subroutine finish_exp_smat
!------------------------------------------------------------------------------
   subroutine finish_exp_smat_Rdist(dVSrLM, ds_exp_last)

      complex(cp), intent(inout) :: dVSrLM(lm_max,nRstart:nRstop)

      !-- Output variables
      complex(cp), intent(inout) :: ds_exp_last(lm_max,nRstart:nRstop)

      !-- Local variables
      complex(cp) :: work_Rloc(lm_max,nRstart:nRstop)
      integer :: n_r, lm, l, start_lm, stop_lm

      call get_dr_Rloc(dVSrLM, work_Rloc, lm_max, nRstart, nRstop, n_r_max, &
           &           rscheme_oc)

      !$omp parallel default(shared) private(n_r, lm, start_lm, stop_lm, l)
      start_lm=1; stop_lm=lm_max
      call get_openmp_blocks(start_lm, stop_lm)
      !$omp barrier

      do n_r=nRstart,nRstop
         do lm=start_lm,stop_lm
            l=st_map%lm2l(lm)
            if ( l > l_R(n_r) ) cycle
            ds_exp_last(lm,n_r)=orho1(n_r)*(ds_exp_last(lm,n_r) - &
            &                        or2(n_r)*work_Rloc(lm,n_r))
         end do
      end do
      !$omp end parallel

   end subroutine finish_exp_smat_Rdist
!------------------------------------------------------------------------------
   subroutine assemble_single(s, ds, w, dw, ddw, dsdt, dwdt, dpdt, tscheme, lRmsNext)
      !
      ! This routine is used to assemble the solution in case IMEX RK with
      ! an assembly stage are used
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lRmsNext

      !-- Output variables
      type(type_tarray), intent(inout) :: dsdt
      type(type_tarray), intent(inout) :: dwdt
      type(type_tarray), intent(inout) :: dpdt
      complex(cp),       intent(inout) :: s(llm:ulm,n_r_max)
      complex(cp),       intent(inout) :: w(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: ds(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: dw(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: ddw(llm:ulm,n_r_max)

      !-- Local variables:
      complex(cp) :: Dif, Buo
      integer :: n_r, lm, l1, m1, lms, start_lm,stop_lm, lmStart_00
      integer :: n_r_top, n_r_bot
      real(cp) :: dL, a_top, b_top, a_bot, b_bot
      real(cp) :: a2_cmb, a1_cmb, a0_cmb, a2_icb, a1_icb, a0_icb ! Second b.c. on W
      integer, pointer :: lm2l(:), lm2m(:)

      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m
      if ( m_min == 0 ) then
         lmStart_00=max(2,llm)
      else
         lmStart_00=llm
      end if

      if ( l_temperature_diff ) then
         call abortRun('Temperature diff + assembly stage not supported!')
      end if

      !-- The pressure is not assembled here
      if ( ktops > 2 .or. kbots > 2 ) then
         call abortRun('Temperature boundary conditions + assembly stage not supported!')
      end if

      !-- First assemble and store in temporary arrays
      call tscheme%assemble_imex(ddw, dwdt)
      call tscheme%assemble_imex(work_LMloc, dpdt)
      call tscheme%assemble_imex(ds, dsdt)

      !$omp parallel default(shared)  private(start_lm, stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)

      !-- Now get the fields from the assembly
      !$omp do private(n_r,lm,l1,m1,dL)
      do n_r=2,n_r_max-1
         do lm=llm,ulm
            l1 = lm2l(lm)
            m1 = lm2m(lm)
            dL = real(l1*(l1+1),cp)
            if ( m1 == 0 ) then
               if ( l1 > 0 ) then
                  w(lm,n_r) = r(n_r)*r(n_r)/dL*cmplx(real(ddw(lm,n_r)),0.0_cp,cp)
                  dw(lm,n_r)=-r(n_r)*r(n_r)/dL*cmplx(real(work_LMloc(lm,n_r)),0.0_cp,cp)
               end if
               s(lm,n_r) = cmplx(real(ds(lm,n_r)),0.0_cp,cp)
            else
               if ( l1 > 0 ) then
                  w(lm,n_r) = r(n_r)*r(n_r)/dL*ddw(lm,n_r)
                  dw(lm,n_r)=-r(n_r)*r(n_r)/dL*work_LMloc(lm,n_r)
               end if
               s(lm,n_r) = ds(lm,n_r)
            end if
         end do
      end do
      !$omp end do

      !-- Get the entropy boundary points using Canuto (1986) approach
      !$omp do private(lm,l1,m1,a_top,b_top,a_bot,b_bot)
      do lm=llm,ulm
         l1 = lm2l(lm)
         m1 = lm2m(lm)
         call get_wps_sbc_coeffs(n_r_cmb, l1, a_top, b_top)
         call get_wps_sbc_coeffs(n_r_icb, l1, a_bot, b_bot)
         call rscheme_oc%robin_bc(a_top, b_top, tops(l1,m1)+sbc_cmb(st_map%lm2(l1,m1)), &
              &                   a_bot, b_bot, bots(l1,m1), s(lm,:))
      end do
      !$omp end do

      !-- Boundary conditions on W: W and a2*ddW+a1*dW+a0*W at both boundaries
      !$omp do private(lm,l1,lms,a2_cmb,a1_cmb,a0_cmb,a2_icb,a1_icb,a0_icb)
      do lm=lmStart_00,ulm
         l1 =lm2l(lm)
         lms=st_map%lm2(l1,lm2m(lm))
         w(lm,1)      =wbc_cmb(lms)
         w(lm,n_r_max)=wbc_icb(lms)
         call get_bc2_coeffs(n_r_cmb, l1, a2_cmb, a1_cmb, a0_cmb)
         call get_bc2_coeffs(n_r_icb, l1, a2_icb, a1_icb, a0_icb)
         call rscheme_oc%robin_bc(a2_cmb, a1_cmb, wbc2_cmb(lms)-a0_cmb*w(lm,1),     &
              &                   a2_icb, a1_icb, wbc2_icb(lms)-a0_icb*w(lm,n_r_max), &
              &                   dw(lm,:))
      end do
      !$omp end do

      !$omp single
      call dct_counter%start_count()
      !$omp end single
      call get_ddr(s, ds, workB, ulm-llm+1, start_lm-llm+1, stop_lm-llm+1, &
           &       n_r_max, rscheme_oc)
      call get_ddr( dw, ddw, work_LMloc, ulm-llm+1, start_lm-llm+1, &
           &         stop_lm-llm+1, n_r_max, rscheme_oc)
      !$omp barrier
      !$omp single
      call dct_counter%stop_count()
      !$omp end single

      !$omp do private(n_r,lm,l1,dL)
      do n_r=2,n_r_max-1
         do lm=llm,ulm
            l1 = lm2l(lm)
            dL = real(l1*(l1+1),cp)
            dsdt%old(lm,n_r,1)= s(lm,n_r)
            dwdt%old(lm,n_r,1)= dL*or2(n_r)*w(lm,n_r)
            dpdt%old(lm,n_r,1)=-dL*or2(n_r)*dw(lm,n_r)
         end do
      end do
      !$omp end do

      if ( tscheme%l_imp_calc_rhs(1) .or. lRmsNext ) then

         if ( lRmsNext ) then
            n_r_top=n_r_cmb
            n_r_bot=n_r_icb
         else
            n_r_top=n_r_cmb+1
            n_r_bot=n_r_icb-1
         end if

         !$omp do private(n_r,lm,l1,m1,Dif,Buo,dL)
         do n_r=n_r_top,n_r_bot
            do lm=llm,ulm
               l1=lm2l(lm)
               m1=lm2m(lm)
               dL = real(l1*(l1+1),cp)

               Dif = hdif_V(l1)*dL*or2(n_r)*visc(n_r) *  (    ddw(lm,n_r) &
               &    +(two*dLvisc(n_r)-third*beta(n_r))*        dw(lm,n_r) &
               &    -( dL*or2(n_r)+four*third*( dbeta(n_r)+dLvisc(n_r)*   &
               &       beta(n_r)+(three*dLvisc(n_r)+beta(n_r))*or1(n_r)))*&
               &                                                w(lm,n_r) )
               Buo = BuoFac*rho0(n_r)*rgrav(n_r)*s(lm,n_r)
               dwdt%impl(lm,n_r,1)=Buo+Dif
               dpdt%impl(lm,n_r,1)=       hdif_V(l1)* visc(n_r)*dL*or2(n_r) &
               &                                   * ( -work_LMloc(lm,n_r)  &
               &                     + (beta(n_r)-dLvisc(n_r))*ddw(lm,n_r)  &
               &        + ( dL*or2(n_r)+dLvisc(n_r)*beta(n_r)+ dbeta(n_r)   &
               &                   + two*(dLvisc(n_r)+beta(n_r))*or1(n_r)   &
               &                                            ) * dw(lm,n_r)  &
               &        - dL*or2(n_r)*( two*or1(n_r)+two*third*beta(n_r)    &
               &                      +dLvisc(n_r) )   *         w(lm,n_r) )
               dsdt%impl(lm,n_r,1)=               opr*hdif_S(l1)*kappa(n_r)*&
               &        ( workB(lm,n_r) + (beta(n_r)+dLtemp0(n_r)+          &
               &            two*or1(n_r) + dLkappa(n_r) )  * ds(lm,n_r)     &
               &                 - dL*or2(n_r) * s(lm,n_r) ) -dL*or2(n_r)   &
               &              *orho1(n_r)*dentropy0(n_r)*        w(lm,n_r)

               if ( lRmsNext .and. l1 > 0 ) then
                  DifPol2hInt(l1,n_r)=DifPol2hInt(l1,n_r)+r(n_r)**2*cc2real(Dif,m1)
                  DifPolLMr(lm,n_r)  =r(n_r)**2/dL * Dif
               end if

            end do
         end do
         !$omp end do
      end if

      !$omp end parallel

   end subroutine assemble_single
!------------------------------------------------------------------------------
   subroutine get_single_rhs_imp(s, ds, w, dw, ddw, p, dp, dsdt, dwdt, dpdt, &
              &                  tscheme, istage, l_calc_lin, lRmsNext,      &
              &                  l_in_cheb_space)

      !-- Input variables
      integer,             intent(in) :: istage
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: l_calc_lin
      logical,             intent(in) :: lRmsNext
      logical, optional,   intent(in) :: l_in_cheb_space

      !-- Output variables
      type(type_tarray), intent(inout) :: dsdt
      type(type_tarray), intent(inout) :: dwdt
      type(type_tarray), intent(inout) :: dpdt
      complex(cp),       intent(inout) :: s(llm:ulm,n_r_max)
      complex(cp),       intent(inout) :: w(llm:ulm,n_r_max)
      complex(cp),       intent(inout) :: p(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: ds(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: dp(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: dw(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: ddw(llm:ulm,n_r_max)

      !-- Local variables
      complex(cp) :: Pre, Buo, Dif
      logical :: l_in_cheb
      real(cp) :: dL
      integer :: n_r_top, n_r_bot, l1, m1
      integer :: n_r, lm, start_lm, stop_lm
      integer, pointer :: lm2l(:),lm2m(:)

      if ( present(l_in_cheb_space) ) then
         l_in_cheb = l_in_cheb_space
      else
         l_in_cheb = .false.
      end if

      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      !$omp parallel default(shared)  private(start_lm, stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)

      !$omp single
      call dct_counter%start_count()
      !$omp end single
      call get_ddr( s, ds, workB, ulm-llm+1, start_lm-llm+1, stop_lm-llm+1, &
           &        n_r_max, rscheme_oc, l_dct_in=.not. l_in_cheb)
      if ( l_in_cheb ) call rscheme_oc%costf1(s,ulm-llm+1,start_lm-llm+1, &
                            &                 stop_lm-llm+1)
      call get_dddr( w, dw, ddw, work_LMloc, ulm-llm+1, start_lm-llm+1, &
           &         stop_lm-llm+1, n_r_max, rscheme_oc,                &
           &         l_dct_in=.not. l_in_cheb)
      if ( l_in_cheb ) call rscheme_oc%costf1(w,ulm-llm+1,start_lm-llm+1, &
                            &                 stop_lm-llm+1)
      call get_ddr( p, dp, workC, ulm-llm+1, start_lm-llm+1, stop_lm-llm+1, &
           &       n_r_max, rscheme_oc, l_dct_in=.not. l_in_cheb)
      if ( l_in_cheb ) call rscheme_oc%costf1(p,ulm-llm+1,start_lm-llm+1, &
                            &                 stop_lm-llm+1)
      !$omp barrier
      !$omp single
      call dct_counter%stop_count()
      !$omp end single

      if ( istage == 1 ) then
         !$omp do private(n_r,lm,l1,dL)
         do n_r=2,n_r_max-1
            do lm=llm,ulm
               l1 = lm2l(lm)
               dL = real(l1*(l1+1),cp)
               dsdt%old(lm,n_r,istage)= s(lm,n_r)
               dwdt%old(lm,n_r,istage)= dL*or2(n_r)*w(lm,n_r)
               dpdt%old(lm,n_r,istage)=-dL*or2(n_r)*dw(lm,n_r)
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

         !-- Calculate explicit time step part:
         if ( l_temperature_diff ) then
            !$omp do private(n_r,lm,l1,m1,Dif,Pre,Buo,dL)
            do n_r=n_r_top,n_r_bot
               do lm=llm,ulm
                  l1=lm2l(lm)
                  m1=lm2m(lm)
                  dL = real(l1*(l1+1),cp)

                  Dif = hdif_V(l1)*dL*or2(n_r)*visc(n_r) *  (    ddw(lm,n_r) &
                  &    +(two*dLvisc(n_r)-third*beta(n_r))*        dw(lm,n_r) &
                  &    -( dL*or2(n_r)+four*third* (dbeta(n_r)+dLvisc(n_r)*   &
                  &      beta(n_r)+(three*dLvisc(n_r)+beta(n_r))*or1(n_r)))* &
                  &                                                w(lm,n_r) )
                  Pre = -dp(lm,n_r)+beta(n_r)*p(lm,n_r)
                  Buo = BuoFac*rho0(n_r)*rgrav(n_r)*s(lm,n_r)
                  dwdt%impl(lm,n_r,istage)=Pre+Buo+Dif
                  dpdt%impl(lm,n_r,istage)=               dL*or2(n_r)*p(lm,n_r)  &
                  &           + hdif_V(l1)*visc(n_r)*dL*or2(n_r)                 &
                  &                                     * ( -work_LMloc(lm,n_r)  &
                  &                   + (beta(n_r)-dLvisc(n_r))    *ddw(lm,n_r)  &
                  &           + ( dL*or2(n_r)+dLvisc(n_r)*beta(n_r)+dbeta(n_r)   &
                  &               +two*(dLvisc(n_r)+beta(n_r))*or1(n_r) ) *      &
                  &                                                  dw(lm,n_r)  &
                  &           - dL*or2(n_r)*( two*or1(n_r)+two*third*beta(n_r)   &
                  &                      +dLvisc(n_r) )   *           w(lm,n_r) )
                  dsdt%impl(lm,n_r,istage)=opr*hdif_S(l1)* kappa(n_r)*(          &
                  &                                               workB(lm,n_r)  &
                  &          + ( beta(n_r)+two*dLtemp0(n_r)+two*or1(n_r)+        &
                  &              dLkappa(n_r) )                    * ds(lm,n_r)  &
                  &          + ( ddLtemp0(n_r)+ dLtemp0(n_r)*( two*or1(n_r)+     &
                  &              dLkappa(n_r)+dLtemp0(n_r)+beta(n_r))-dL*        &
                  &              or2(n_r) ) *                         s(lm,n_r)  &
                  &          +  alpha0(n_r)*orho1(n_r)*ViscHeatFac*ThExpNb*(     &
                  &                                               workC(lm,n_r)  &
                  &          +  ( dLkappa(n_r)+two*(dLtemp0(n_r)+dLalpha0(n_r))  &
                  &               +two*or1(n_r)-beta(n_r) ) *        dp(lm,n_r)  &
                  &          +  ( (dLkappa(n_r)+dLtemp0(n_r)+dLalpha0(n_r)+      &
                  &               two*or1(n_r))*(dLalpha0(n_r)+dLtemp0(n_r)-     &
                  &               beta(n_r))+ddLtemp0(n_r)+ddLalpha0(n_r)-       &
                  &               dbeta(n_r)-dL*or2(n_r) )*           p(lm,n_r)))&
                  &          - dL*or2(n_r)*orho1(n_r)*dentropy0(n_r)* w(lm,n_r)

                  if ( lRmsNext .and. l1 > 0 ) then
                     DifPol2hInt(l1,n_r)=DifPol2hInt(l1,n_r)+r(n_r)**2*cc2real(Dif,m1)
                     DifPolLMr(lm,n_r)  =r(n_r)**2/dL * Dif
                  end if
               end do
            end do
            !$omp end do

         else ! entropy diffusion

            !$omp do private(n_r,lm,l1,m1,Dif,Pre,Buo,dL)
            do n_r=n_r_top,n_r_bot
               do lm=llm,ulm
                  l1=lm2l(lm)
                  m1=lm2m(lm)
                  dL = real(l1*(l1+1),cp)

                  Dif = hdif_V(l1)*dL*or2(n_r)*visc(n_r) *  (    ddw(lm,n_r) &
                  &    +(two*dLvisc(n_r)-third*beta(n_r))*        dw(lm,n_r) &
                  &    -( dL*or2(n_r)+four*third*( dbeta(n_r)+dLvisc(n_r)*   &
                  &       beta(n_r)+(three*dLvisc(n_r)+beta(n_r))*or1(n_r)))*&
                  &                                                w(lm,n_r) )
                  Pre = -dp(lm,n_r)+beta(n_r)*p(lm,n_r)
                  Buo = BuoFac*rho0(n_r)*rgrav(n_r)*s(lm,n_r)
                  dwdt%impl(lm,n_r,istage)=Pre+Buo+Dif
                  dpdt%impl(lm,n_r,istage)=               dL*or2(n_r)*p(lm,n_r)&
                  &                        + hdif_V(l1)* visc(n_r)*dL*or2(n_r) &
                  &                                   * ( -work_LMloc(lm,n_r)  &
                  &                     + (beta(n_r)-dLvisc(n_r))*ddw(lm,n_r)  &
                  &        + ( dL*or2(n_r)+dLvisc(n_r)*beta(n_r)+ dbeta(n_r)   &
                  &                   + two*(dLvisc(n_r)+beta(n_r))*or1(n_r)   &
                  &                                            ) * dw(lm,n_r)  &
                  &        - dL*or2(n_r)*( two*or1(n_r)+two*third*beta(n_r)    &
                  &                      +dLvisc(n_r) )   *         w(lm,n_r) )
                  dsdt%impl(lm,n_r,istage)=          opr*hdif_S(l1)*kappa(n_r)*&
                  &        ( workB(lm,n_r) + (beta(n_r)+dLtemp0(n_r)+          &
                  &            two*or1(n_r) + dLkappa(n_r) )  * ds(lm,n_r)     &
                  &                 - dL*or2(n_r) * s(lm,n_r) ) -dL*or2(n_r)   &
                  &              *orho1(n_r)*dentropy0(n_r)*        w(lm,n_r)

                  if ( lRmsNext .and. l1 > 0 ) then
                     DifPol2hInt(l1,n_r)=DifPol2hInt(l1,n_r)+r(n_r)**2*cc2real(Dif,m1)
                     DifPolLMr(lm,n_r)  =r(n_r)**2/dL * Dif
                  end if

               end do
            end do
            !$omp end do
         end if

      end if
      !$omp end parallel

   end subroutine get_single_rhs_imp
!------------------------------------------------------------------------------
   subroutine get_wpsMat(tscheme,l,hdif_vel,hdif_s,wpsMat,wpsPivot,wpsMat_fac)
      !
      !  Purpose of this subroutine is to contruct the time step matrix
      !  wpmat  for the NS equation.
      !

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme
      real(cp),            intent(in) :: hdif_vel
      real(cp),            intent(in) :: hdif_s
      integer,             intent(in) :: l

      !-- Output variables:
      real(cp), intent(out) :: wpsMat(3*n_r_max,3*n_r_max)
      real(cp), intent(out) :: wpsMat_fac(3*n_r_max,2)
      integer,  intent(out) :: wpsPivot(3*n_r_max)

      !-- local variables:
      integer :: nR,nR_out,nR_p,nR_s,nR_out_p,nR_out_s,info
      real(cp) :: dLh
      real(cp) :: a2_cmb, a1_cmb, a0_cmb, a2_icb, a1_icb, a0_icb ! Second b.c. on W
      real(cp) :: as_cmb, bs_cmb, ap_cmb, bp_cmb, as_icb, bs_icb, ap_icb, bp_icb ! b.c. on s

      dLh =real(l*(l+1),kind=cp)
      call get_bc2_coeffs(n_r_cmb, l, a2_cmb, a1_cmb, a0_cmb)
      call get_bc2_coeffs(n_r_icb, l, a2_icb, a1_icb, a0_icb)
      call get_entropy_bc_coeffs(n_r_cmb, l, as_cmb, bs_cmb, ap_cmb, bp_cmb)
      call get_entropy_bc_coeffs(n_r_icb, l, as_icb, bs_icb, ap_icb, bp_icb)

      !-- Now mode l>0

      !----- Boundary conditions, see above:
      do nR_out=1,rscheme_oc%n_max
         nR_out_p=nR_out+n_r_max
         nR_out_s=nR_out+2*n_r_max

         wpsMat(1,nR_out)        =rscheme_oc%rnorm*rscheme_oc%rMat(1,nR_out)
         wpsMat(1,nR_out_p)      =0.0_cp
         wpsMat(1,nR_out_s)      =0.0_cp
         wpsMat(n_r_max,nR_out)  =rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,nR_out)
         wpsMat(n_r_max,nR_out_p)=0.0_cp
         wpsMat(n_r_max,nR_out_s)=0.0_cp

         wpsMat(n_r_max+1,nR_out)=rscheme_oc%rnorm*(                &
         &                  a2_cmb*rscheme_oc%d2rMat(1,nR_out)+     &
         &                  a1_cmb*rscheme_oc%drMat(1,nR_out)+      &
         &                  a0_cmb*rscheme_oc%rMat(1,nR_out))
         wpsMat(n_r_max+1,nR_out_p)=0.0_cp
         wpsMat(n_r_max+1,nR_out_s)=0.0_cp

         wpsMat(2*n_r_max,nR_out)=rscheme_oc%rnorm*(                &
         &                  a2_icb*rscheme_oc%d2rMat(n_r_max,nR_out)+ &
         &                  a1_icb*rscheme_oc%drMat(n_r_max,nR_out)+  &
         &                  a0_icb*rscheme_oc%rMat(n_r_max,nR_out))
         wpsMat(2*n_r_max,nR_out_p)=0.0_cp
         wpsMat(2*n_r_max,nR_out_s)=0.0_cp

         wpsMat(2*n_r_max+1,nR_out_s)=rscheme_oc%rnorm*(                    &
         &                            as_cmb*rscheme_oc%drMat(1,nR_out)+ &
         &                            bs_cmb*rscheme_oc%rMat(1,nR_out))
         wpsMat(2*n_r_max+1,nR_out_p)=rscheme_oc%rnorm*(                    &
         &                            ap_cmb*rscheme_oc%drMat(1,nR_out)+ &
         &                            bp_cmb*rscheme_oc%rMat(1,nR_out))
         wpsMat(2*n_r_max+1,nR_out)  =0.0_cp

         wpsMat(3*n_r_max,nR_out_s)=rscheme_oc%rnorm*(                          &
         &                          as_icb*rscheme_oc%drMat(n_r_max,nR_out)+ &
         &                          bs_icb*rscheme_oc%rMat(n_r_max,nR_out))
         wpsMat(3*n_r_max,nR_out_p)=rscheme_oc%rnorm*(                          &
         &                          ap_icb*rscheme_oc%drMat(n_r_max,nR_out)+ &
         &                          bp_icb*rscheme_oc%rMat(n_r_max,nR_out))
         wpsMat(3*n_r_max,nR_out)  =0.0_cp


      end do   !  loop over nR_out

      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme_oc%n_max+1,n_r_max
            nR_out_p=nR_out+n_r_max
            nR_out_s=nR_out+2*n_r_max
            wpsMat(1,nR_out)            =0.0_cp
            wpsMat(n_r_max,nR_out)      =0.0_cp
            wpsMat(n_r_max+1,nR_out)    =0.0_cp
            wpsMat(2*n_r_max,nR_out)    =0.0_cp
            wpsMat(2*n_r_max+1,nR_out)  =0.0_cp
            wpsMat(3*n_r_max,nR_out)    =0.0_cp
            wpsMat(1,nR_out_p)          =0.0_cp
            wpsMat(n_r_max,nR_out_p)    =0.0_cp
            wpsMat(n_r_max+1,nR_out_p)  =0.0_cp
            wpsMat(2*n_r_max,nR_out_p)  =0.0_cp
            wpsMat(2*n_r_max+1,nR_out_p)=0.0_cp
            wpsMat(3*n_r_max,nR_out_p)  =0.0_cp
            wpsMat(1,nR_out_s)          =0.0_cp
            wpsMat(n_r_max,nR_out_s)    =0.0_cp
            wpsMat(n_r_max+1,nR_out_s)  =0.0_cp
            wpsMat(2*n_r_max,nR_out_s)  =0.0_cp
            wpsMat(2*n_r_max+1,nR_out_s)=0.0_cp
            wpsMat(3*n_r_max,nR_out_s)  =0.0_cp
         end do
      end if

      if ( l_temperature_diff ) then ! temperature diffusion

         do nR_out=1,n_r_max
            nR_out_p=nR_out+n_r_max
            nR_out_s=nR_out+2*n_r_max
            do nR=2,n_r_max-1
               nR_p=nR+n_r_max
               nR_s=nR+2*n_r_max

               ! W equation
               wpsMat(nR,nR_out)= rscheme_oc%rnorm *  (                        &
               &                 dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out)        &
               &  - tscheme%wimp_lin(1)*hdif_vel*visc(nR)*dLh*or2(nR) * (      &
               &                                rscheme_oc%d2rMat(nR,nR_out)   &
               &+(two*dLvisc(nR)-third*beta(nR))*rscheme_oc%drMat(nR,nR_out)   &
               &         -( dLh*or2(nR)+four*third*( dLvisc(nR)*beta(nR)       &
               &          +(three*dLvisc(nR)+beta(nR))*or1(nR)+dbeta(nR) )     &
               &          )                      *rscheme_oc%rMat(nR,nR_out) )  )

               ! Buoyancy
               wpsMat(nR,nR_out_s)=-rscheme_oc%rnorm*tscheme%wimp_lin(1)*  &
               &                    BuoFac*rgrav(nR)*rho0(nR)*             &
               &                                 rscheme_oc%rMat(nR,nR_out)

               ! Pressure gradient
               wpsMat(nR,nR_out_p)= rscheme_oc%rnorm*tscheme%wimp_lin(1)*(  &
               &                              rscheme_oc%drMat(nR,nR_out)   &
               &                    -beta(nR)* rscheme_oc%rMat(nR,nR_out) )

               ! P equation
               wpsMat(nR_p,nR_out)= rscheme_oc%rnorm * (                       &
               &             -dLh*or2(nR)*    rscheme_oc%drMat(nR,nR_out)      &
               &   -tscheme%wimp_lin(1)*hdif_vel*visc(nR)*dLh*or2(nR)      *(  &
               &                                 -rscheme_oc%d3rMat(nR,nR_out) &
               &         +( beta(nR)-dLvisc(nR) )*rscheme_oc%d2rMat(nR,nR_out) &
               &             +( dLh*or2(nR)+dbeta(nR)+dLvisc(nR)*beta(nR)      &
               &             +two*(dLvisc(nR)+beta(nR))*or1(nR) )*             &
               &                                   rscheme_oc%drMat(nR,nR_out) &
               &        -dLh*or2(nR)*( two*or1(nR)+dLvisc(nR)                  &
               &           +two*third*beta(nR)   )* rscheme_oc%rMat(nR,nR_out) ) )

               wpsMat(nR_p,nR_out_p)=-rscheme_oc%rnorm*tscheme%wimp_lin(1)*   &
               &                      dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out)

               wpsMat(nR_p,nR_out_s)=0.0_cp

               ! S equation
               wpsMat(nR_s,nR_out_s)= rscheme_oc%rnorm * (                        &
               &                                    rscheme_oc%rMat(nR,nR_out) -  &
               &           tscheme%wimp_lin(1)*opr*hdif_s*kappa(nR)*(             &
               &                                  rscheme_oc%d2rMat(nR,nR_out) +  &
               &      ( beta(nR)+two*dLtemp0(nR)+                                 &
               &        two*or1(nR)+dLkappa(nR) )* rscheme_oc%drMat(nR,nR_out) +  &
               &      ( ddLtemp0(nR)+dLtemp0(nR)*(                                &
               &  two*or1(nR)+dLkappa(nR)+dLtemp0(nR)+beta(nR) )   -              &
               &           dLh*or2(nR) )*           rscheme_oc%rMat(nR,nR_out) ) )

               wpsMat(nR_s,nR_out_p)= -tscheme%wimp_lin(1)*rscheme_oc%rnorm*      &
               &           hdif_s*kappa(nR)*opr*alpha0(nR)*orho1(nR)*             &
               &           ViscHeatFac*ThExpNb*(  rscheme_oc%d2rMat(nR,nR_out) +  &
               &      ( dLkappa(nR)+two*(dLalpha0(nR)+dLtemp0(nR)) -              &
               &        beta(nR) +two*or1(nR) ) *  rscheme_oc%drMat(nR,nR_out) +  &
               & ( (dLkappa(nR)+dLalpha0(nR)+dLtemp0(nR)+two*or1(nR)) *           &
               &        (dLalpha0(nR)+dLtemp0(nR)-beta(nR)) +                     &
               &        ddLalpha0(nR)+ddLtemp0(nR)-dbeta(nR)-                     &
               &        dLh*or2(nR) ) *             rscheme_oc%rMat(nR,nR_out) )


               !Advection of the background entropy u_r * dso/dr
               wpsMat(nR_s,nR_out)=rscheme_oc%rnorm*tscheme%wimp_lin(1)*dLh*      &
               &                   or2(nR)*dentropy0(nR)*orho1(nR)*               &
               &                                       rscheme_oc%rMat(nR,nR_out)

            end do
         end do

      else ! entropy diffusion

         do nR_out=1,n_r_max
            nR_out_p=nR_out+n_r_max
            nR_out_s=nR_out+2*n_r_max
            do nR=2,n_r_max-1
               nR_p=nR+n_r_max
               nR_s=nR+2*n_r_max

               ! W equation
               wpsMat(nR,nR_out)= rscheme_oc%rnorm *  (                          &
               &                         dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out)  &
               &   - tscheme%wimp_lin(1)*hdif_vel*visc(nR)*dLh*or2(nR) * (       &
               &                                   rscheme_oc%d2rMat(nR,nR_out)  &
               &   +(two*dLvisc(nR)-third*beta(nR))*rscheme_oc%drMat(nR,nR_out)  &
               &         -( dLh*or2(nR)+four*third*( dLvisc(nR)*beta(nR)         &
               &          +(three*dLvisc(nR)+beta(nR))*or1(nR)+dbeta(nR) )       &
               &          )                         *rscheme_oc%rMat(nR,nR_out) ) )

               ! Buoyancy
               wpsMat(nR,nR_out_s)=-rscheme_oc%rnorm*tscheme%wimp_lin(1)*BuoFac* &
               &                    rgrav(nR)*rho0(nR)*rscheme_oc%rMat(nR,nR_out)

               ! Pressure gradient
               wpsMat(nR,nR_out_p)= rscheme_oc%rnorm*tscheme%wimp_lin(1)*(  &
               &                                rscheme_oc%drMat(nR,nR_out) &
               &                      -beta(nR)* rscheme_oc%rMat(nR,nR_out) )

               ! P equation
               wpsMat(nR_p,nR_out)= rscheme_oc%rnorm * (                         &
               &                  -     dLh*or2(nR)*rscheme_oc%drMat(nR,nR_out)  &
               &  -tscheme%wimp_lin(1)*hdif_vel*visc(nR)*dLh*or2(nR)  *(         &
                                                  -rscheme_oc%d3rMat(nR,nR_out)  &
               &   +( beta(nR)-dLvisc(nR) )*       rscheme_oc%d2rMat(nR,nR_out)  &
               &          +( dLh*or2(nR)+dbeta(nR)+dLvisc(nR)*beta(nR)           &
               &          +two*(dLvisc(nR)+beta(nR))*or1(nR) )*                  &
               &                                    rscheme_oc%drMat(nR,nR_out)  &
               &          -dLh*or2(nR)*( two*or1(nR)+dLvisc(nR)                  &
               &            +two*third*beta(nR)   )* rscheme_oc%rMat(nR,nR_out)  ) )

               wpsMat(nR_p,nR_out_p)= -rscheme_oc%rnorm*tscheme%wimp_lin(1)*dLh* &
               &                       or2(nR)*rscheme_oc%rMat(nR,nR_out)

               wpsMat(nR_p,nR_out_s)=0.0_cp

               ! S equation
               wpsMat(nR_s,nR_out_s)= rscheme_oc%rnorm * (                       &
               &                                    rscheme_oc%rMat(nR,nR_out) - &
               &              tscheme%wimp_lin(1)*opr*hdif_s*kappa(nR)*(         &
               &                                  rscheme_oc%d2rMat(nR,nR_out)+  &
               &      ( beta(nR)+dLtemp0(nR)+                                    &
               &        two*or1(nR)+dLkappa(nR) )*rscheme_oc%drMat(nR,nR_out) -  &
               &           dLh*or2(nR)*            rscheme_oc%rMat(nR,nR_out) ) )

               wpsMat(nR_s,nR_out_p)=0.0_cp ! temperature diffusion

               !Advection of the background entropy u_r * dso/dr
               wpsMat(nR_s,nR_out)=rscheme_oc%rnorm*tscheme%wimp_lin(1)*dLh* &
               &                   or2(nR)*dentropy0(nR)*orho1(nR)*          &
               &                                   rscheme_oc%rMat(nR,nR_out)

            end do
         end do
      end if

      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         nR_p=nR+n_r_max
         nR_s=nR+2*n_r_max
         wpsMat(nR,1)            =rscheme_oc%boundary_fac*wpsMat(nR,1)
         wpsMat(nR,n_r_max)      =rscheme_oc%boundary_fac*wpsMat(nR,n_r_max)
         wpsMat(nR,n_r_max+1)    =rscheme_oc%boundary_fac*wpsMat(nR,n_r_max+1)
         wpsMat(nR,2*n_r_max)    =rscheme_oc%boundary_fac*wpsMat(nR,2*n_r_max)
         wpsMat(nR,2*n_r_max+1)  =rscheme_oc%boundary_fac*wpsMat(nR,2*n_r_max+1)
         wpsMat(nR,3*n_r_max)    =rscheme_oc%boundary_fac*wpsMat(nR,3*n_r_max)
         wpsMat(nR_p,1)          =rscheme_oc%boundary_fac*wpsMat(nR_p,1)
         wpsMat(nR_p,n_r_max)    =rscheme_oc%boundary_fac*wpsMat(nR_p,n_r_max)
         wpsMat(nR_p,n_r_max+1)  =rscheme_oc%boundary_fac*wpsMat(nR_p,n_r_max+1)
         wpsMat(nR_p,2*n_r_max)  =rscheme_oc%boundary_fac*wpsMat(nR_p,2*n_r_max)
         wpsMat(nR_p,2*n_r_max+1)=rscheme_oc%boundary_fac*wpsMat(nR_p,2*n_r_max+1)
         wpsMat(nR_p,3*n_r_max)  =rscheme_oc%boundary_fac*wpsMat(nR_p,3*n_r_max)
         wpsMat(nR_s,1)          =rscheme_oc%boundary_fac*wpsMat(nR_s,1)
         wpsMat(nR_s,n_r_max)    =rscheme_oc%boundary_fac*wpsMat(nR_s,n_r_max)
         wpsMat(nR_s,n_r_max+1)  =rscheme_oc%boundary_fac*wpsMat(nR_s,n_r_max+1)
         wpsMat(nR_s,2*n_r_max)  =rscheme_oc%boundary_fac*wpsMat(nR_s,2*n_r_max)
         wpsMat(nR_s,2*n_r_max+1)=rscheme_oc%boundary_fac*wpsMat(nR_s,2*n_r_max+1)
         wpsMat(nR_s,3*n_r_max)  =rscheme_oc%boundary_fac*wpsMat(nR_s,3*n_r_max)
      end do

      ! compute the linesum of each line
      do nR=1,3*n_r_max
         wpsMat_fac(nR,1)=one/maxval(abs(wpsMat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,3*n_r_max
         wpsMat(nR,:) = wpsMat(nR,:)*wpsMat_fac(nR,1)
      end do

      ! also compute the rowsum of each column
      do nR=1,3*n_r_max
         wpsMat_fac(nR,2)=one/maxval(abs(wpsMat(:,nR)))
      end do
      ! now divide each row by the rowsum
      do nR=1,3*n_r_max
         wpsMat(:,nR) = wpsMat(:,nR)*wpsMat_fac(nR,2)
      end do

      call prepare_mat(wpsMat,3*n_r_max,3*n_r_max,wpsPivot,info)
      if ( info /= 0 ) then
         call abortRun('Singular matrix wpsMat!')
      end if

   end subroutine get_wpsMat
!-----------------------------------------------------------------------------
   subroutine get_ps0Mat(tscheme,psMat,psPivot,psMat_fac)

      !-- Input variable:
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variables:
      real(cp), intent(out) :: psMat(2*n_r_max,2*n_r_max)
      integer,  intent(out) :: psPivot(2*n_r_max)
      real(cp), intent(out) :: psMat_fac(2*n_r_max,2)

      !-- Local variables:
      integer :: info,nCheb,nR_out,nR_out_p,nR,nR_p,n_cheb_in
      real(cp) :: work(n_r_max),work2(n_r_max)
      real(cp) :: as_cmb, bs_cmb, ap_cmb, bp_cmb, as_icb, bs_icb, ap_icb, bp_icb ! b.c. on s

      if ( l_temperature_diff ) then ! temperature diffusion

         do nR_out=1,n_r_max
            nR_out_p=nR_out+n_r_max
            do nR=1,n_r_max
               nR_p=nR+n_r_max

               psMat(nR,nR_out)= rscheme_oc%rnorm * (                              &
               &                                   rscheme_oc%rMat(nR,nR_out) -    &
               &          tscheme%wimp_lin(1)*opr*kappa(nR)*(                      &
               &                                 rscheme_oc%d2rMat(nR,nR_out) +    &
               &      ( beta(nR)+two*dLtemp0(nR)+                                  &
               &       two*or1(nR)+dLkappa(nR) )* rscheme_oc%drMat(nR,nR_out) +    &
               &      ( ddLtemp0(nR)+dLtemp0(nR)*(                                 &
               &  two*or1(nR)+dLkappa(nR)+dLtemp0(nR)+beta(nR) ) ) *               &
               &                                   rscheme_oc%rMat(nR,nR_out) ) )

               psMat(nR,nR_out_p)= -tscheme%wimp_lin(1)*rscheme_oc%rnorm*kappa(nR)*&
               &        opr*alpha0(nR)*orho1(nR)*ViscHeatFac*ThExpNb*(             &
               &                              rscheme_oc%d2rMat(nR,nR_out) +       &
               &      ( dLkappa(nR)+two*(dLalpha0(nR)+dLtemp0(nR)) -               &
               &    beta(nR) +two*or1(nR) ) *  rscheme_oc%drMat(nR,nR_out) +       &
               & ( (dLkappa(nR)+dLalpha0(nR)+dLtemp0(nR)+two*or1(nR)) *            &
               &        (dLalpha0(nR)+dLtemp0(nR)-beta(nR)) +                      &
               &        ddLalpha0(nR)+ddLtemp0(nR)-dbeta(nR) ) *                   &
               &                                rscheme_oc%rMat(nR,nR_out) )

               psMat(nR_p,nR_out)  = -rscheme_oc%rnorm*rho0(nR)*          &
               &                     BuoFac*rgrav(nR)*rscheme_oc%rMat(nR,nR_out)
               psMat(nR_p,nR_out_p)= rscheme_oc%rnorm *(                  &
               &                             rscheme_oc%drMat(nR,nR_out)- &
               &                     beta(nR)*rscheme_oc%rMat(nR,nR_out) )
            end do
         end do

      else ! entropy diffusion

         do nR_out=1,n_r_max
           nR_out_p=nR_out+n_r_max
            do nR=1,n_r_max
               nR_p=nR+n_r_max

               psMat(nR,nR_out)    = rscheme_oc%rnorm * (                        &
               &                                    rscheme_oc%rMat(nR,nR_out) - &
               &         tscheme%wimp_lin(1)*opr*kappa(nR)*(                     &
               &                                  rscheme_oc%d2rMat(nR,nR_out) + &
               &    (beta(nR)+dLtemp0(nR)+two*or1(nR)+dLkappa(nR))*              &
               &                                   rscheme_oc%drMat(nR,nR_out) ) )
               psMat(nR,nR_out_p)  =0.0_cp ! entropy diffusion

               psMat(nR_p,nR_out)  = -rscheme_oc%rnorm*BuoFac*rho0(nR)* &
               &                     rgrav(nR)*rscheme_oc%rMat(nR,nR_out)
               psMat(nR_p,nR_out_p)= rscheme_oc%rnorm*(rscheme_oc%drMat(nR,nR_out)- &
               &                               beta(nR)*rscheme_oc%rMat(nR,nR_out) )
            end do
         end do

      end if


      !----- Boundary condition:
      call get_entropy_bc_coeffs(n_r_cmb, 0, as_cmb, bs_cmb, ap_cmb, bp_cmb)
      call get_entropy_bc_coeffs(n_r_icb, 0, as_icb, bs_icb, ap_icb, bp_icb)
      do nR_out=1,rscheme_oc%n_max
         nR_out_p=nR_out+n_r_max

         psMat(1,nR_out)        =rscheme_oc%rnorm*(as_cmb*rscheme_oc%drMat(1,nR_out)+ &
         &                                         bs_cmb*rscheme_oc%rMat(1,nR_out))
         psMat(1,nR_out_p)      =rscheme_oc%rnorm*(ap_cmb*rscheme_oc%drMat(1,nR_out)+ &
         &                                         bp_cmb*rscheme_oc%rMat(1,nR_out))
         psMat(n_r_max,nR_out)  =rscheme_oc%rnorm*(                         &
         &                       as_icb*rscheme_oc%drMat(n_r_max,nR_out)+   &
         &                       bs_icb*rscheme_oc%rMat(n_r_max,nR_out))
         psMat(n_r_max,nR_out_p)=rscheme_oc%rnorm*(                         &
         &                       ap_icb*rscheme_oc%drMat(n_r_max,nR_out)+   &
         &                       bp_icb*rscheme_oc%rMat(n_r_max,nR_out))

      end do

      ! In case density perturbations feed back on pressure (non-Boussinesq)
      ! Impose that the integral of (rho' r^2) vanishes
      if ( ViscHeatFac*ThExpNb /= 0.0_cp .and. ktopp == 1 ) then

         work(:)=ThExpNb*ViscHeatFac*ogrun(:)*alpha0(:)*r(:)*r(:)
         call rscheme_oc%costf1(work)
         work         =work*rscheme_oc%rnorm
         work(1)      =rscheme_oc%boundary_fac*work(1)
         work(n_r_max)=rscheme_oc%boundary_fac*work(n_r_max)

         work2(:)=-ThExpNb*alpha0(:)*temp0(:)*rho0(:)*r(:)*r(:)
         call rscheme_oc%costf1(work2)
         work2         =work2*rscheme_oc%rnorm
         work2(1)      =rscheme_oc%boundary_fac*work2(1)
         work2(n_r_max)=rscheme_oc%boundary_fac*work2(n_r_max)

         if ( rscheme_oc%version == 'cheb' ) then

            do nCheb=1,rscheme_oc%n_max
               nR_out_p=nCheb+n_r_max
               psMat(n_r_max+1,nR_out_p)=0.0_cp
               psMat(n_r_max+1,nCheb)   =0.0_cp
               do n_cheb_in=1,rscheme_oc%n_max
                  if (mod(nCheb+n_cheb_in-2,2)==0) then
                     psMat(n_r_max+1,nR_out_p)=psMat(n_r_max+1,nR_out_p)+           &
                     &                     (one/(one-real(n_cheb_in-nCheb,cp)**2)+  &
                     &                     one/(one-real(n_cheb_in+nCheb-2,cp)**2))*&
                     &                       work(n_cheb_in)*half*rscheme_oc%rnorm
                     psMat(n_r_max+1,nCheb)  =psMat(n_r_max+1,nCheb)+               &
                     &                     (one/(one-real(n_cheb_in-nCheb,cp)**2)+  &
                     &                     one/(one-real(n_cheb_in+nCheb-2,cp)**2))*&
                     &                     work2(n_cheb_in)*half*rscheme_oc%rnorm
                  end if
               end do
            end do

         else

            !-- In the finite differences case, we restrict the integral boundary
            !-- condition to a trapezoidal rule of integration
            do nR_out=2,rscheme_oc%n_max-1
               nR_out_p=nR_out+n_r_max
               psMat(n_r_max+1,nR_out)  =half*work2(nR_out)*( r(nR_out+1)-r(nR_out-1) )
               psMat(n_r_max+1,nR_out_p)=half* work(nR_out)*( r(nR_out+1)-r(nR_out-1) )
            end do
            psMat(n_r_max+1,1)        =half*work2(1)*( r(2)-r(1) )
            psMat(n_r_max+1,n_r_max+1)=half* work(1)*( r(2)-r(1) )
            psMat(n_r_max+1,n_r_max)  =half*work2(n_r_max)*( r(n_r_max)-r(n_r_max-1) )
            psMat(n_r_max+1,2*n_r_max)=half* work(n_r_max)*( r(n_r_max)-r(n_r_max-1) )

         end if

      else

         do nR_out=1,rscheme_oc%n_max
            nR_out_p=nR_out+n_r_max
            psMat(n_r_max+1,nR_out)  =0.0_cp
            psMat(n_r_max+1,nR_out_p)=rscheme_oc%rnorm*rscheme_oc%rMat(1,nR_out)
         end do

      end if

      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme_oc%n_max+1,n_r_max
            nR_out_p=nR_out+n_r_max
            psMat(1,nR_out)          =0.0_cp
            psMat(n_r_max,nR_out)    =0.0_cp
            psMat(n_r_max+1,nR_out)  =0.0_cp
            psMat(2*n_r_max,nR_out)  =0.0_cp
            psMat(1,nR_out_p)        =0.0_cp
            psMat(n_r_max,nR_out_p)  =0.0_cp
            psMat(n_r_max+1,nR_out_p)=0.0_cp
         end do
      end if

      !----- Factors for highest and lowest cheb mode:
      do nR=1,n_r_max
         nR_p=nR+n_r_max
         psMat(nR,1)          =rscheme_oc%boundary_fac*psMat(nR,1)
         psMat(nR,n_r_max)    =rscheme_oc%boundary_fac*psMat(nR,n_r_max)
         psMat(nR,n_r_max+1)  =rscheme_oc%boundary_fac*psMat(nR,n_r_max+1)
         psMat(nR,2*n_r_max)  =rscheme_oc%boundary_fac*psMat(nR,2*n_r_max)
         psMat(nR_p,1)        =rscheme_oc%boundary_fac*psMat(nR_p,1)
         psMat(nR_p,n_r_max)  =rscheme_oc%boundary_fac*psMat(nR_p,n_r_max)
         psMat(nR_p,n_r_max+1)=rscheme_oc%boundary_fac*psMat(nR_p,n_r_max+1)
         psMat(nR_p,2*n_r_max)=rscheme_oc%boundary_fac*psMat(nR_p,2*n_r_max)
      end do

      ! compute the linesum of each line
      do nR=1,2*n_r_max
         psMat_fac(nR,1)=one/maxval(abs(psMat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,2*n_r_max
         psMat(nR,:) = psMat(nR,:)*psMat_fac(nR,1)
      end do

      ! also compute the rowsum of each column
      do nR=1,2*n_r_max
         psMat_fac(nR,2)=one/maxval(abs(psMat(:,nR)))
      end do
      ! now divide each row by the rowsum
      do nR=1,2*n_r_max
         psMat(:,nR) = psMat(:,nR)*psMat_fac(nR,2)
      end do


      !---- LU decomposition:
      call prepare_mat(psMat,2*n_r_max,2*n_r_max,psPivot,info)
      if ( info /= 0 ) then
         call abortRun('! Singular matrix ps0Mat!')
      end if

   end subroutine get_ps0Mat
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
   subroutine get_entropy_bc_coeffs(nR, l, a_s, b_s, a_p, b_p)
      !
      ! Boundary condition on the entropy when it is solved together with
      ! W and P (single matrix): a_s*ds/dr+b_s*s+a_p*dp/dr+b_p*p=rhs.
      ! ktops/kbots: 1 fixed entropy, 2 fixed entropy flux (see
      ! get_wps_sbc_coeffs), 3 fixed temperature, 4 fixed temperature flux.
      ! Temperature conditions involve the pressure.
      !

      !-- Input variables:
      integer, intent(in) :: nR         ! n_r_cmb or n_r_icb
      integer, intent(in) :: l

      !-- Output variables:
      real(cp), intent(out) :: a_s, b_s, a_p, b_p

      !-- Local variables:
      integer :: k
      real(cp) :: fac_p

      a_p=0.0_cp
      b_p=0.0_cp
      k=ktops
      if ( nR == n_r_icb ) k=kbots
      if ( k <= 2 .or. (nR == n_r_icb .and. l_full_sphere) ) then
         call get_wps_sbc_coeffs(nR, l, a_s, b_s)
         return
      end if

      fac_p=orho1(nR)*alpha0(nR)*temp0(nR)*ViscHeatFac*ThExpNb
      if ( k == 3 ) then ! Fixed temperature
         a_s=0.0_cp
         b_s=temp0(nR)
         b_p=fac_p
      else ! Fixed temperature flux
         a_s=temp0(nR)
         b_s=temp0(nR)*dLtemp0(nR)
         a_p=fac_p
         b_p=fac_p*(dLalpha0(nR)+dLtemp0(nR)-beta(nR))
      end if

   end subroutine get_entropy_bc_coeffs
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
   subroutine update_wps_bc(time)
      !
      ! Boundary values of the single-matrix solver at a given time
      ! (see update_wps_wbc and update_wps_sbc). It is called once per
      ! stage, before the boundary values are used.
      !

      !-- Input variable:
      real(cp), intent(in) :: time

      call update_wps_wbc(time)
      call update_wps_sbc(time)

   end subroutine update_wps_bc
!-----------------------------------------------------------------------------
   subroutine update_wps_wbc(time)
      !
      ! W and the r.h.s. of the second b.c. on W due to the prescribed motion
      ! of the boundaries (ellipticity, tides, see get_bc2_coeffs).
      !

      !-- Input variable:
      real(cp), intent(in) :: time

      !-- Local variables:
      complex(cp) :: wtide, eiwt
      integer :: lm20, lm22

      if ( .not. (l_vr_cmb .or. l_vr_icb) ) return

      !-- Only the (l=2,m=0) and (l=2,m=2) modes are forced, the other entries
      !-- stay zero. Both modes are in the truncation, except for an
      !-- axisymmetric run which keeps the m=0 tide only (see readNamelists).
      lm20=st_map%lm2(2,0)
      lm22=st_map%lm2(2,2)

      if ( l_vr_cmb ) then
         eiwt=cmplx(cos(omega_tide*time), sin(omega_tide*time), cp)
         wtide=tide_fac20*real(eiwt)
         wbc_cmb(lm20)=wtide
         if ( ktopv /= 1 ) wbc2_cmb(lm20)=6.0_cp*l2_h2_ratio/r_cmb*wtide
         if ( .not. l_axi ) then
            !-- Relative sign from dPhi_2^-2 < 0: with the time origin chosen so that
            !-- the (2,0) flow goes as cos(wt), the sub-main body point stays at phi=0
            wtide=tide_fac22p*eiwt-tide_fac22n*conjg(eiwt)
            wbc_cmb(lm22)=wtide+ellip_fac_cmb*                          &
            &             cmplx(cos(omegaOsz_ma1*(time+tShift_ma1)),    &
            &                   sin(omegaOsz_ma1*(time+tShift_ma1)), cp)
            if ( ktopv /= 1 ) wbc2_cmb(lm22)=6.0_cp*l2_h2_ratio/r_cmb*wtide
         end if
      end if

      if ( l_vr_icb ) then
         wbc_icb(lm22)=ellip_fac_icb*cmplx(cos(omegaOsz_ic1*(time+tShift_ic1)), &
         &                                 sin(omegaOsz_ic1*(time+tShift_ic1)), cp)
      end if

   end subroutine update_wps_wbc
!-----------------------------------------------------------------------------
   subroutine update_wps_sbc(time)
      !
      ! Perturbation of the outer entropy/temperature b.c. due to the tidal
      ! displacement delta of the boundary: a fixed value (ktops=1,3) holds on
      ! r_cmb+delta, which to first order adds -delta*df0/dr to its r.h.s. at
      ! r_cmb (Taylor expansion, see the notes on tidal forcing). f0 is the
      ! current spherically-symmetric entropy (ktops=1) or temperature
      ! (ktops=3). For a fixed flux (ktops=2,4), the correction -delta*d2f0/dr2
      ! is only O(delta/H) at a no-slip wall (see update_sbc) and is omitted.
      !

      !-- Input variable:
      real(cp), intent(in) :: time

      !-- Local variables:
      complex(cp) :: delta20, delta22
      real(cp) :: fac, f0_dr, fac_p, s00, ds00, p00, dp00, dT_s, dT_p
      integer :: lm20, lm22, lm00, nR
#ifdef WITH_MPI
      integer :: ierr
#endif

      if ( .not. l_sbc_tide ) return

      !-- d(f00)/dr at r_cmb, from the derivatives computed with the implicit
      !-- terms (only known on the rank that holds l=m=0)
      f0_dr=0.0_cp
      lm00=lo_map%lm2(0,0)
      if ( llm <= lm00 .and. ulm >= lm00 ) then
         nR=n_r_cmb
         s00 =real(s_LMloc(lm00,nR))
         ds00=real(ds_LMloc(lm00,nR))
         if ( ktops == 1 ) then ! Fixed entropy: f0=s
            f0_dr=ds00
         else ! Fixed temperature: f0=T=temp0*s+fac_p*p (see get_entropy_bc_coeffs)
            p00 =real(p_LMloc(lm00,nR))
            dp00=real(dp_LMloc(lm00,nR))
            fac_p=orho1(nR)*alpha0(nR)*temp0(nR)*ViscHeatFac*ThExpNb
            !-- d(temp0)/dr=temp0*dLtemp0, d(fac_p)/dr=fac_p*(dLalpha0+dLtemp0-beta)
            dT_s=temp0(nR)*(ds00+dLtemp0(nR)*s00)                    ! d(temp0*s)/dr
            dT_p=fac_p*(dp00+(dLalpha0(nR)+dLtemp0(nR)-beta(nR))*p00) ! d(fac_p*p)/dr
            f0_dr=dT_s+dT_p ! dT_p=0 for Boussinesq (fac_p=0)
         end if
      end if
#ifdef WITH_MPI
      call MPI_Allreduce(MPI_IN_PLACE, f0_dr, 1, MPI_DEF_REAL, MPI_SUM, &
           &             MPI_COMM_WORLD, ierr)
#endif

      !-- Radial displacement of the outer boundary by the tide: time
      !-- integral of u_r=l(l+1)/r^2*W
      lm20=st_map%lm2(2,0)
      lm22=st_map%lm2(2,2)
      fac=6.0_cp/r_cmb**2/omega_tide
      delta20=fac*tide_fac20*cmplx(sin(omega_tide*time), 0.0_cp, cp)
      if ( .not. l_axi ) delta22=fac*cmplx( (tide_fac22p-tide_fac22n)*sin(omega_tide*time), &
      &                                    -(tide_fac22p+tide_fac22n)*cos(omega_tide*time), cp)
      fac=-osq4pi*f0_dr
      sbc_cmb(lm20)=fac*delta20 ! Other entries stay zero
      if ( .not. l_axi ) sbc_cmb(lm22)=fac*delta22 ! m=0 tide only if axisymmetric

   end subroutine update_wps_sbc
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
   subroutine get_wps_sbc_coeffs(nR, l, a, b)
      !
      ! Boundary condition on the entropy for ktops/kbots=1 (fixed entropy) or
      ! 2 (fixed flux): a*ds/dr+b*s=rhs. At the centre of a full sphere,
      ! regularity requires ds/dr=0 for l=0 and s=0 otherwise.
      !

      !-- Input variables:
      integer, intent(in) :: nR ! n_r_cmb or n_r_icb
      integer, intent(in) :: l

      !-- Output variables:
      real(cp), intent(out) :: a, b

      !-- Local variable:
      integer :: k

      if ( nR == n_r_icb .and. l_full_sphere ) then
         if ( l == 0 ) then
            k=2
         else
            k=1
         end if
      else if ( nR == n_r_cmb ) then
         k=ktops
      else
         k=kbots
      end if

      if ( k == 1 ) then ! Fixed value
         a=0.0_cp
         b=one
      else ! Fixed flux
         a=one
         b=0.0_cp
      end if

   end subroutine get_wps_sbc_coeffs
!-----------------------------------------------------------------------------
   subroutine get_bc2_coeffs(nR, l, a2, a1, a0)
      !
      ! Second b.c. on W at the boundary nR: a2*d2W/dr2+a1*dW/dr+a0*W=wbc2
      ! (rigid, stress-free with a moving boundary, or regularity at r=0).
      !

      !-- Input variables:
      integer, intent(in) :: nR ! n_r_cmb or n_r_icb
      integer, intent(in) :: l

      !-- Output variables:
      real(cp), intent(out) :: a2, a1, a0

      !-- Local variable:
      integer :: kv

      a2=0.0_cp
      a1=one
      a0=0.0_cp
      if ( nR == n_r_icb .and. l_full_sphere ) then
         if ( l /= 1 ) then ! d2W/dr2=0, dW/dr=0 for l=1
            a2=one
            a1=0.0_cp
         end if
         return
      end if

      kv=ktopv
      if ( nR == n_r_icb ) kv=kbotv
      if ( kv == 1 ) then ! Stress-free
         a2=one
         a1=-(two*or1(nR)+beta(nR))
         a0=real(l*(l+1),cp)*or2(nR)
      end if

   end subroutine get_bc2_coeffs
!-----------------------------------------------------------------------------
end module updateWPS_mod
