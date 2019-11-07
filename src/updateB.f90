#include "perflib_preproc.cpp"
module updateB_mod

#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif

   use omp_lib
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, n_r_tot, n_r_ic_max,             &
       &                 n_cheb_ic_max, n_r_ic_maxMag, n_r_maxMag, &
       &                 n_r_totMag, lm_max, l_maxMag
   use radial_functions, only: chebt_ic,or2,r_cmb,chebt_ic_even, d2cheb_ic,    &
       &                       cheb_norm_ic,dr_fac_ic,lambda,dLlambda,o_r_ic,r,&
       &                       or1, cheb_ic, dcheb_ic,rscheme_oc
   use radial_data, only: n_r_cmb, n_r_icb
   use physical_parameters, only: n_r_LCR,opm,O_sr,kbotb, imagcon, tmagcon, &
       &                         sigma_ratio, conductance_ma, ktopb, kbotb
   use init_fields, only: bpeaktop, bpeakbot
   use num_param, only: solve_counter, dct_counter
   use blocking, only: st_map, lo_map, st_sub_map, lo_sub_map, llmMag, ulmMag
   use horizontal_data, only: dLh, dPhi, hdif_B, D_l, D_lP1
   use logic, only: l_cond_ic, l_LCR, l_rot_ic, l_mag_nl, l_b_nl_icb, &
       &            l_b_nl_cmb, l_update_b, l_RMS, l_finite_diff,     &
       &            l_full_sphere
   use RMS, only: dtBPolLMr, dtBPol2hInt, dtBTor2hInt
   use constants, only: pi, zero, one, two, three, half
   use special
   use parallel_mod, only:  rank, chunksize, n_procs, get_openmp_blocks
   use RMS_helpers, only: hInt2PolLM, hInt2TorLM
   use fields, only: work_LMloc
   use radial_der_even, only: get_ddr_even
   use radial_der, only: get_dr, get_ddr
   use useful, only: abortRun
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray
   use dense_matrices
   use band_matrices
   use real_matrices

   implicit none

   private

   !-- Local work arrays:
   complex(cp), allocatable :: workB(:,:), work1_LMloc(:,:)
   complex(cp), allocatable :: work_ic_LMloc(:,:), work1_ic_LMloc(:,:)
   complex(cp), allocatable :: rhs1(:,:,:),rhs2(:,:,:)
   complex(cp), allocatable :: dtT(:), dtP(:)
   class(type_realmat), pointer :: bMat(:), jMat(:)
#ifdef WITH_PRECOND_BJ
   real(cp), allocatable :: bMat_fac(:,:)
   real(cp), allocatable :: jMat_fac(:,:)
#endif
   logical, public, allocatable :: lBmat(:)

   public :: initialize_updateB, finalize_updateB, updateB, finish_exp_mag, &
   &         get_mag_rhs_imp, get_mag_ic_rhs_imp, finish_exp_mag_ic

contains

   subroutine initialize_updateB

      integer, pointer :: nLMBs2(:)
      integer :: maxThreads, ll, n_bandsJ, n_bandsB

      nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

      if ( l_finite_diff .and. (.not. l_cond_ic) ) then
         allocate( type_bandmat :: jMat(nLMBs2(1+rank)) )
         allocate( type_bandmat :: bMat(nLMBs2(1+rank)) )

         if ( kbotb == 2 .or. ktopb == 2 .or. conductance_ma /= 0 .or. &
         &    rscheme_oc%order  > 2 .or. rscheme_oc%order_boundary > 2 ) then
            !-- Perfect conductor or conducting mantle
            n_bandsJ = max(2*rscheme_oc%order_boundary+1,rscheme_oc%order+1)
         else
            n_bandsJ = rscheme_oc%order+1
         end if

         if ( conductance_ma /= 0 ) then
            !-- Second derivative on the boundary
            n_bandsB = max(2*rscheme_oc%order_boundary+3,rscheme_oc%order+1)
         else
            n_bandsB = max(2*rscheme_oc%order_boundary+1,rscheme_oc%order+1)
         end if

         !print*, 'J, B', n_bandsJ, n_bandsB
         do ll=1,nLMBs2(1+rank)
            call bMat(ll)%initialize(n_bandsB,n_r_tot,l_pivot=.true.)
            call jMat(ll)%initialize(n_bandsJ,n_r_tot,l_pivot=.true.)
         end do
      else
         allocate( type_densemat :: jMat(nLMBs2(1+rank)) )
         allocate( type_densemat :: bMat(nLMBs2(1+rank)) )

         do ll=1,nLMBs2(1+rank)
            call bMat(ll)%initialize(n_r_tot,n_r_tot,l_pivot=.true.)
            call jMat(ll)%initialize(n_r_tot,n_r_tot,l_pivot=.true.)
         end do
      end if

#ifdef WITH_PRECOND_BJ
      allocate(bMat_fac(n_r_tot,nLMBs2(1+rank)))
      allocate(jMat_fac(n_r_tot,nLMBs2(1+rank)))
      bytes_allocated = bytes_allocated+2*n_r_tot*nLMBs2(1+rank)*SIZEOF_DEF_REAL
#endif
      allocate( lBmat(0:l_maxMag) )
      bytes_allocated = bytes_allocated+(l_maxMag+1)*SIZEOF_LOGICAL


      if ( l_RMS ) then
         allocate( workB(llmMag:ulmMag,n_r_max) )
         bytes_allocated = bytes_allocated+(ulmMag-llmMag+1)*n_r_max* &
         &                 SIZEOF_DEF_COMPLEX
      end if

      if ( l_cond_ic ) then
         allocate( work_ic_LMloc(llmMag:ulmMag,n_r_ic_max) )
         allocate( work1_ic_LMloc(llmMag:ulmMag,n_r_ic_max) )
         bytes_allocated = bytes_allocated+2*(ulmMag-llmMag+1)*n_r_ic_max* &
         &                 SIZEOF_DEF_COMPLEX
      end if

      allocate( work1_LMloc(llmMag:ulmMag,n_r_max) )
      bytes_allocated = bytes_allocated+(ulmMag-llmMag+1)*n_r_max* &
      &                 SIZEOF_DEF_COMPLEX

      allocate( dtT(llmMag:ulmMag) )
      allocate( dtP(llmMag:ulmMag) )
      bytes_allocated = bytes_allocated+2*(ulmMag-llmMag+1)*SIZEOF_DEF_COMPLEX

#ifdef WITHOMP
      maxThreads=omp_get_max_threads()
#else
      maxThreads=1
#endif

      allocate(rhs1(2*n_r_max,lo_sub_map%sizeLMB2max,0:maxThreads-1))
      allocate(rhs2(2*n_r_max,lo_sub_map%sizeLMB2max,0:maxThreads-1))
      bytes_allocated=bytes_allocated+4*n_r_max*maxThreads* &
      &               lo_sub_map%sizeLMB2max*SIZEOF_DEF_COMPLEX

   end subroutine initialize_updateB
!-----------------------------------------------------------------------------
   subroutine finalize_updateB

      integer, pointer :: nLMBs2(:)
      integer :: ll

      nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

      do ll=1,nLMBs2(1+rank)
         call jMat(ll)%finalize()
         call bMat(ll)%finalize()
      end do

      if ( l_cond_ic ) deallocate ( work_ic_LMloc, work1_ic_LMloc )
      deallocate( lBmat, work1_LMLoc )

#ifdef WITH_PRECOND_BJ
      deallocate(bMat_fac,jMat_fac)
#endif
      deallocate( dtT, dtP, rhs1, rhs2 )
      if ( l_RMS ) deallocate( workB )

   end subroutine finalize_updateB
!-----------------------------------------------------------------------------
   subroutine updateB(b,db,ddb,aj,dj,ddj,dbdt,djdt,b_ic,db_ic,ddb_ic,aj_ic,  &
              &       dj_ic,ddj_ic,dbdt_ic,djdt_ic,b_nl_cmb,aj_nl_cmb,       &
              &       aj_nl_icb,time,tscheme,lRmsNext)
      !
      !
      !  Calculated update of magnetic field potential and the time
      !  stepping arrays dbdtLast, ...
      !
      !
      !  updates the magnetic field potentials b, aj and
      !  their derivatives,
      !  adds explicit part to time derivatives of b and j

      !  input:  w1 - weight for dbdt-contribution from current time step
      !               (w2=1-w1: weight for contrib. from previous step)
      !          coex - factor depending on weighting alpha of
      !                 implicit contribution
      !          mc_min,mc_max
      !               - range of mca-indices in which field is updated
      !                 (harmonic order m=(mca-1)*minc)
      !          b_nl_cmb = RHS of nonlinear BC for poloidal magnetic field
      !                      potential in the case of stress free CMB
      !          aj_nl_cmb = RHS of nonlinear BC for toroidal magnetic field
      !                      potential in the case of stress free CMB
      !          aj_nl_icb = RHS of nonlinear BC for radial derivative of
      !                      toroidal magnetic field
      !                      potential in the case of stress free ICB
      !

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme
      complex(cp),         intent(in) :: b_nl_cmb(:)  ! nonlinear BC for b at CMB
      complex(cp),         intent(in) :: aj_nl_cmb(:) ! nonlinear BC for aj at CMB
      complex(cp),         intent(in) :: aj_nl_icb(:) ! nonlinear BC for dr aj at ICB
      real(cp),            intent(in) :: time
      logical,             intent(in) :: lRmsNext

      !-- Input/output of scalar potentials and time stepping arrays:
      type(type_tarray), intent(inout) :: dbdt
      type(type_tarray), intent(inout) :: djdt
      type(type_tarray), intent(inout) :: dbdt_ic
      type(type_tarray), intent(inout) :: djdt_ic
      complex(cp),       intent(inout) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(cp),       intent(inout) :: aj(llmMag:ulmMag,n_r_maxMag)
      complex(cp),       intent(inout) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp),       intent(inout) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)

      !-- Output variables:
      complex(cp), intent(out) :: db(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(out) :: ddb(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(out) :: dj(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(out) :: ddj(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(out) :: db_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(out) :: ddb_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(out) :: dj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(out) :: ddj_ic(llmMag:ulmMag,n_r_ic_maxMag)

      !-- Local variables:
      real(cp) :: yl0_norm,prefac!External magnetic field of general l

      integer :: l1,m1               ! degree and order
      integer :: lm1,lm,lmB          ! position of (l,m) in array
      integer :: lmStart_00          ! excluding l=0,m=0
      integer :: nLMB2, nLMB
      integer :: n_r_out             ! No of cheb polynome (degree+1)
      integer :: nR                  ! No of radial grid point

      integer, pointer :: nLMBs2(:),lm2l(:),lm2m(:)
      integer, pointer :: sizeLMB2(:,:),lm2(:,:)
      integer, pointer :: lm22lm(:,:,:),lm22l(:,:,:),lm22m(:,:,:)

      !-- for feedback
      real(cp) :: ff,cimp,aimp,b10max

      real(cp), save :: direction

      integer :: start_lm, stop_lm
      integer :: nChunks,iChunk,lmB0,size_of_last_chunk,threadid

      if ( .not. l_update_b ) RETURN

      nLMBs2(1:n_procs) => lo_sub_map%nLMBs2
      sizeLMB2(1:,1:) => lo_sub_map%sizeLMB2
      lm22lm(1:,1:,1:) => lo_sub_map%lm22lm
      lm22l(1:,1:,1:) => lo_sub_map%lm22l
      lm22m(1:,1:,1:) => lo_sub_map%lm22m
      lm2(0:,0:) => lo_map%lm2
      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      nLMB=1+rank
      lmStart_00=max(2,llmMag)

      call get_mag_rhs_imp(b, db, ddb, aj, dj, ddj,                 &
           &               dbdt%old(:,:,tscheme%istage),            &
           &               djdt%old(:,:,tscheme%istage),            &
           &               dbdt%impl(:,:,tscheme%istage),           &
           &               djdt%impl(:,:,tscheme%istage),           &
           &               tscheme%l_imp_calc_rhs(tscheme%istage),  &
           &               lRmsNext)

      !-- Now assemble the right hand side and store it in work_LMloc
      call tscheme%set_imex_rhs(work_LMloc, dbdt, llmMag, ulmMag, n_r_maxMag)
      call tscheme%set_imex_rhs(work1_LMloc, djdt, llmMag, ulmMag, n_r_maxMag)

      if ( l_cond_ic ) then
         call get_mag_ic_rhs_imp(b_ic, db_ic, ddb_ic, aj_ic, dj_ic, ddj_ic,  &
              &                  dbdt_ic%old(:,:,tscheme%istage),            &
              &                  djdt_ic%old(:,:,tscheme%istage),            &
              &                  dbdt_ic%impl(:,:,tscheme%istage),           &
              &                  djdt_ic%impl(:,:,tscheme%istage),           &
              &                  tscheme%l_imp_calc_rhs(tscheme%istage))

         !-- Now assemble the right hand side and store it in work_LMloc
         call tscheme%set_imex_rhs(work_ic_LMloc, dbdt_ic, llmMag, ulmMag, &
              &                    n_r_ic_max)
         call tscheme%set_imex_rhs(work1_ic_LMloc, djdt_ic, llmMag, ulmMag, &
              &                    n_r_ic_max)
      end if


      !$omp parallel default(shared) private(start_lm,stop_lm)
      start_lm=lmStart_00; stop_lm=ulmMag
      call get_openmp_blocks(start_lm,stop_lm)

      !$omp single
      call solve_counter%start_count()
      !$omp end single
      ! This is a loop over all l values which should be treated on
      ! the actual MPI rank
      !$OMP SINGLE
      do nLMB2=1,nLMBs2(nLMB)
         !$OMP TASK default(shared) &
         !$OMP firstprivate(nLMB2) &
         !$OMP private(lmB,lm,lm1,l1,m1,nR,iChunk,nChunks,size_of_last_chunk) &
         !$OMP private(fac,bpeaktop,ff,cimp,aimp,threadid)

         ! determine the number of chunks of m
         ! total number for l1 is sizeLMB2(nLMB2,nLMB)
         ! chunksize is given
         nChunks = (sizeLMB2(nLMB2,nLMB)+chunksize-1)/chunksize
         size_of_last_chunk = chunksize + (sizeLMB2(nLMB2,nLMB)-nChunks*chunksize)

         l1=lm22l(1,nLMB2,nLMB)
         if ( l1 > 0 ) then
            if ( .not. lBmat(l1) ) then
#ifdef WITH_PRECOND_BJ
               call get_bMat(tscheme,l1,hdif_B(st_map%lm2(l1,0)),   &
                    &        bMat(nLMB2),bMat_fac(:,nLMB2),         &
                    &        jMat(nLMB2),jMat_fac(:,nLMB2))
#else
               call get_bMat(tscheme,l1,hdif_B(st_map%lm2(l1,0)),   &
                    &        bMat(nLMB2),jMat(nLMB2))
#endif
               lBmat(l1)=.true.
            end if
         end if

         do iChunk=1,nChunks
            !$OMP TASK if (nChunks>1) default(shared) &
            !$OMP firstprivate(iChunk) &
            !$OMP private(lmB0,lmB,lm,lm1,m1,nR) &
            !$OMP private(fac,bpeaktop,ff) &
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
               if ( l1 > 0 ) then
                  lmB=lmB+1
                  !-------- Magnetic boundary conditions, outer core:
                  !         Note: the CMB condition is not correct if we assume free slip
                  !         and a conducting mantle (conductance_ma>0).
                  if ( l_b_nl_cmb ) then ! finitely conducting mantle
                     rhs1(1,lmB,threadid) =  b_nl_cmb(st_map%lm2(l1,m1))
                     rhs2(1,lmB,threadid) = aj_nl_cmb(st_map%lm2(l1,m1))
                  else
                     rhs1(1,lmB,threadid)=0.0_cp
                     rhs2(1,lmB,threadid) = 0.0_cp
                  end if

                  rhs1(n_r_max,lmB,threadid)=0.0_cp
                  if ( kbotb == 2 ) rhs1(n_r_max-1,lmB,threadid)=0.0_cp

                  rhs2(n_r_max,lmB,threadid)=0.0_cp
                  if ( m1 == 0 ) then   ! Magnetoconvection boundary conditions
                     if ( imagcon /= 0 .and. tmagcon <= time ) then
                        if ( l_LCR ) then
                           call abortRun('LCR not compatible with imposed field!')
                        end if
                        if ( l1 == 2 .and. imagcon > 0 .and. imagcon  /=  12 ) then
                           rhs2(1,lmB,threadid)      =cmplx(bpeaktop,0.0_cp,kind=cp)
                           rhs2(n_r_max,lmB,threadid)=cmplx(bpeakbot,0.0_cp,kind=cp)
                        else if( l1 == 1 .and. imagcon == 12 ) then
                           rhs2(1,lmB,threadid)      =cmplx(bpeaktop,0.0_cp,kind=cp)
                           rhs2(n_r_max,lmB,threadid)=cmplx(bpeakbot,0.0_cp,kind=cp)
                        else if( l1 == 1 .and. imagcon == -1) then
                           rhs1(n_r_max,lmB,threadid)=cmplx(bpeakbot,0.0_cp,kind=cp)
                        else if( l1 == 1 .and. imagcon == -2) then
                           rhs1(1,lmB,threadid)      =cmplx(bpeaktop,0.0_cp,kind=cp)
                        else if( l1 == 3 .and. imagcon == -10 ) then
                           rhs2(1,lmB,threadid)      =cmplx(bpeaktop,0.0_cp,kind=cp)
                           rhs2(n_r_max,lmB,threadid)=cmplx(bpeakbot,0.0_cp,kind=cp)
                        end if
                     end if

                    if (l_curr .and. (mod(l1,2) /= 0) ) then    !Current carrying loop around equator of sphere, only odd harmonics

                        if ( l_LCR ) then
                           call abortRun('LCR not compatible with imposed field!')
                        end if

                        !General normalization for spherical harmonics of degree l and order 0
                        yl0_norm=half*sqrt((2*l1+1)/pi)

                        !Prefactor for CMB matching condition
                        prefac = real(2*l1+1,kind=cp)/real(l1*(l1+1),kind=cp)


                        bpeaktop=prefac*fac_loop(l1)*amp_curr*r_cmb/yl0_norm

                        rhs1(1,lmB,threadid)=cmplx(bpeaktop,0.0_cp,kind=cp)

                     end if

                     if ( n_imp > 1 .and. l1 == l_imp ) then
                         if ( l_LCR ) then
                            call abortRun('LCR not compatible with imposed field!')
                         end if
                         ! General normalization for degree l and order 0
                         yl0_norm = half*sqrt((2*l1+1)/pi)
                         ! Prefactor for CMB matching condition
                         prefac = real(2*l1+1,kind=cp)/real(l1*(l1+1),kind=cp)

                        if ( n_imp == 2 ) then
                           !  Chose external field coefficient so that amp_imp is
                           !  the amplitude of the external field:
                           bpeaktop=prefac*r_cmb/yl0_norm*amp_imp
                           rhs1(1,lmB,threadid)=cmplx(bpeaktop,0.0_cp,kind=cp)
                        else if ( n_imp == 3 ) then
                           !  Chose external field coefficient so that amp_imp is
                           !  the amplitude of the external field:
                           bpeaktop=prefac*r_cmb/yl0_norm*amp_imp
                           if ( real(b(2,1)) > 1.0e-9_cp ) &
                                direction=real(b(2,1))/abs(real(b(2,1)))
                           rhs1(1,lmB,threadid)=cmplx(bpeaktop,0.0_cp,kind=cp)*direction
                        else if ( n_imp == 4 ) then
                           !  I have forgotten what this was supposed to do:
                           bpeaktop=three/r_cmb*amp_imp*real(b(2,1))**2
                           rhs1(1,lmB,threadid)=cmplx(bpeaktop,0.0_cp,kind=cp)/b(2,1)

                        else

                           ! Special Heyner feedback functions:
                           ! We don't provide an actual external field strength but rather its
                           ! dependence on the internal dipole via a feedback function ff:
                           !              b10_external=ff(b10_internal)
                           ! where b10_external is the external axial dipole field coeff. and
                           ! b10_internal=b(2,1) is the internal axial dipole field coeff.
                           ! Then rhs1 is always given by:
                           !              rhs1 = (2*l+1)/r_cmb * ff
                           ! because of the special formulation of the CMB matching condition!
                           ! Note that
                           !  B_r_internal(r_cmb) = l*(l+1)/r_cmb**2 * b10_internal*y10_norm*cos(theta)
                           !  B_r_external(r)     = l*(l+1)/r_cmb**2 * b10_external*y10_norm*cos(theta)
                           ! This determines the units of ff=b10_external.
                           ! B itself is given in units sqrt(rho*omega/sigma) so that B**2 is
                           ! the Elsasser number. ff is thus given in units L**2*sqrt(rho*omega/sigma)
                           ! with L=(r_cmb-r_icb). Note that the external dipole field does not depend
                           ! on the radius.
                           ! Here amp_imp provides the relative amplitude so that the maximum of
                           ! the external field is   max(b10_external/b10_internal)=amp_imp

                           ff=0.0_cp
                           if ( n_imp == 7 ) then
                              ! Using a feedback function of the form
                              !    ff= aimp * b10_internal**expo_imp / (cimp+b10_internal**expo_imp)
                              ! Here expo_imp is an input parameter and aimp and cimp are determined
                              ! such that the maximum of b10_external/b10_internal is amp_imp
                              ! and is located at b10_internal=bmax_imp.
                              ! amp_imp and bmax_imp are two more input parameters.
                              ! Note that bmax_imp on input hat the dimensionless unit of the magnetic
                              ! field B and we convert is to dimensionless units for dipole
                              ! coefficients b10 first by multiplying with r_cmb**2
                              ! Other factors like  l*(l+1)*y10_norm*mean(cos(theta)) may also be
                              ! considered, but like r_cmb**2 they are all of order one here.
                              ! Since amp_imp is dimensionless aimp and thus ff has the dimension of
                              ! b10 as required.
                              b10max=bmax_imp*r_cmb**2
                              cimp=b10max**expo_imp / (expo_imp-1)
                              aimp=amp_imp*cimp*expo_imp / &
                              &    (cimp*(expo_imp-1))**((expo_imp-1)/expo_imp)

                              ff=  aimp*real(b(2,1))**expo_imp/ &
                              &    (cimp+real(b(2,1))**expo_imp)

                           end if
                           rhs1(1,lmB,threadid)=(2*l1+1)/r_cmb*ff

                        end if
                     end if
                  end if

                  do nR=2,n_r_max-1
                     if ( nR<=n_r_LCR ) then
                        rhs1(nR,lmB,threadid)=0.0_cp
                        rhs2(nR,lmB,threadid)=0.0_cp
                     else
                        rhs1(nR,lmB,threadid)=work_LMloc(lm1,nR)
                        rhs2(nR,lmB,threadid)=work1_LMloc(lm1,nR)
                     end if
                  end do

                  !-- Magnetic boundary conditions, inner core for radial derivatives
                  !         of poloidal and toroidal magnetic potentials:
                  if ( l_cond_ic ) then    ! inner core
                     rhs1(n_r_max+1,lmB,threadid)=0.0_cp
                     if ( l_b_nl_icb ) then
                        rhs2(n_r_max+1,lmB,threadid)=aj_nl_icb(st_map%lm2(l1,m1))
                     else
                        rhs2(n_r_max+1,lmB,threadid)=0.0_cp
                     end if

                     do nR=2,n_r_ic_max
                        rhs1(n_r_max+nR,lmB,threadid)=work_ic_LMloc(lm1,nR)
                        rhs2(n_r_max+nR,lmB,threadid)=work1_ic_LMloc(lm1,nR)
                     end do
                  end if

               end if ! l>0
            end do    ! loop over lm in block

            if ( lmB > lmB0 ) then
               !write(*,"(2(A,I5))") "updateB: Calling cgeslML for l1=",l1," WITH lmB=",lmB
#ifdef WITH_PRECOND_BJ
               do lm=lmB0+1,lmB
                  do nR=1,n_r_tot
                     rhs1(nR,lm,threadid)=rhs1(nR,lm,threadid)*bMat_fac(nR,nLMB2)
                     rhs2(nR,lm,threadid)=rhs2(nR,lm,threadid)*jMat_fac(nR,nLMB2)
                  end do
               end do
#endif

               !LIKWID_ON('upB_sol')
               call bMat(nLMB2)%solve(rhs1(:,lmB0+1:lmB,threadid),lmB-lmB0)
               call jMat(nLMB2)%solve(rhs2(:,lmB0+1:lmB,threadid),lmB-lmB0)
               !LIKWID_OFF('upB_sol')
            end if

            if ( lRmsNext ) then ! Store old b,aj
               do nR=1,n_r_max
                  do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
                     !do lm=1,sizeLMB2(nLMB2,nLMB)
                     lm1=lm22lm(lm,nLMB2,nLMB)
                     work_LMloc(lm1,nR)=b(lm1,nR)
                     workB(lm1,nR)=aj(lm1,nR)
                  end do
               end do
            end if

            !----- Update magnetic field in cheb space:
            !PERFON('upB_set')
            lmB=lmB0
            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               !do lm=1,sizeLMB2(nLMB2,nLMB)
               lm1=lm22lm(lm,nLMB2,nLMB)
               !l1 =lm22l(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)

               if ( l1 > 0 ) then
                  lmB=lmB+1

                  if ( m1 > 0 ) then
                     do n_r_out=1,rscheme_oc%n_max  ! outer core
                        b(lm1,n_r_out) =rhs1(n_r_out,lmB,threadid)
                        aj(lm1,n_r_out)=rhs2(n_r_out,lmB,threadid)
                     end do
                     if ( l_cond_ic ) then   ! inner core
                        do n_r_out=1,n_cheb_ic_max
                           b_ic(lm1,n_r_out) =rhs1(n_r_max+n_r_out,lmB,threadid)
                           aj_ic(lm1,n_r_out)=rhs2(n_r_max+n_r_out,lmB,threadid)
                        end do
                     end if
                  else
                     do n_r_out=1,rscheme_oc%n_max   ! outer core
                        b(lm1,n_r_out) = &
                        &    cmplx(real(rhs1(n_r_out,lmB,threadid)),0.0_cp,kind=cp)
                        aj(lm1,n_r_out)= &
                        &    cmplx(real(rhs2(n_r_out,lmB,threadid)),0.0_cp,kind=cp)
                     end do
                     if ( l_cond_ic ) then    ! inner core
                        do n_r_out=1,n_cheb_ic_max
                           b_ic(lm1,n_r_out)= cmplx( &
                           & real(rhs1(n_r_max+n_r_out,lmB,threadid)),0.0_cp,kind=cp)
                           aj_ic(lm1,n_r_out)= cmplx( &
                           & real(rhs2(n_r_max+n_r_out,lmB,threadid)),0.0_cp,kind=cp)
                        end do
                     end if
                  end if

               end if
            end do
            !$OMP END TASK
         end do
         !$OMP END TASK
      end do      ! end of do loop over lm1
      !$OMP END SINGLE
      !$omp single
      call solve_counter%stop_count(l_increment=.false.)
      !$omp end single

      !-- Set cheb modes > rscheme_oc%n_max to zero (dealiazing)
      !   for inner core modes > 2*n_cheb_ic_max = 0
      !$omp single
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         do lm1=lmStart_00,ulmMag
            b(lm1,n_r_out) =zero
            aj(lm1,n_r_out)=zero
         end do
      end do
      !$omp end single

      if ( l_cond_ic ) then
         !$omp do private(n_r_out, lm1)
         do n_r_out=n_cheb_ic_max+1,n_r_ic_max
            do lm1=lmStart_00,ulmMag
               b_ic(lm1,n_r_out) =zero
               aj_ic(lm1,n_r_out)=zero
            end do
         end do
         !$omp end do
      end if

      !$omp single
      call dct_counter%start_count()
      !$omp end single
      !-- Radial derivatives: dbdtLast and djdtLast used as work arrays
      call get_ddr(b,db,ddb,ulmMag-llmMag+1,start_lm-llmMag+1, &
           &       stop_lm-llmMag+1,n_r_max,rscheme_oc,l_dct_in=.false.)
      call rscheme_oc%costf1(b,ulmMag-llmMag+1,start_lm-llmMag+1, &
           &                 stop_lm-llmMag+1)

      call get_ddr(aj,dj,ddj,ulmMag-llmMag+1,start_lm-llmMag+1, &
           &       stop_lm-llmMag+1,n_r_max,rscheme_oc,l_dct_in=.false.)
      call rscheme_oc%costf1(aj, ulmMag-llmMag+1, start_lm-llmMag+1, &
           &                 stop_lm-llmMag+1)

      !-- Same for inner core:
      if ( l_cond_ic ) then
         call chebt_ic%costf1( b_ic, ulmMag-llmMag+1, start_lm-llmMag+1, &
              &                stop_lm-llmMag+1, work_LMloc)
         call get_ddr_even( b_ic,db_ic,ddb_ic, ulmMag-llmMag+1, &
              &             start_lm-llmMag+1,stop_lm-llmMag+1, &
              &             n_r_ic_max,n_cheb_ic_max, dr_fac_ic,&
              &             work_LMloc,work1_LMloc, chebt_ic, chebt_ic_even )
         call chebt_ic%costf1( aj_ic, ulmMag-llmMag+1, start_lm-llmMag+1, &
              &               stop_lm-llmMag+1, work_LMloc)
         call get_ddr_even( aj_ic,dj_ic,ddj_ic, ulmMag-llmMag+1,  &
              &             start_lm-llmMag+1, stop_lm-llmMag+1,  &
              &             n_r_ic_max,n_cheb_ic_max, dr_fac_ic,  &
              &             work_LMloc,work1_LMloc, chebt_ic, chebt_ic_even )
      end if
      !$omp barrier
      !PERFOFF
      !$omp single
      call dct_counter%stop_count(l_increment=.false.)
      !$omp end single
      !-- We are now back in radial space !

      !PERFON('upB_last')
      if ( l_LCR ) then
         !$omp do private(nR,lm1,l1,m1)
         do nR=n_r_cmb,n_r_icb-1
            if ( nR<=n_r_LCR ) then
               do lm1=lmStart_00,ulmMag
                  l1=lm2l(lm1)
                  m1=lm2m(lm1)

                  b(lm1,nR)=(r(n_r_LCR)/r(nR))**D_l(st_map%lm2(l1,m1))*b(lm1,n_r_LCR)
                  db(lm1,nR)=-real(D_l(st_map%lm2(l1,m1)),kind=cp)*     &
                  &          (r(n_r_LCR))**D_l(st_map%lm2(l1,m1))/      &
                  &          (r(nR))**(D_l(st_map%lm2(l1,m1))+1)*b(lm1,n_r_LCR)
                  ddb(lm1,nR)=real(D_l(st_map%lm2(l1,m1)),kind=cp)*         &
                  &           (real(D_l(st_map%lm2(l1,m1)),kind=cp)+1)      &
                  &           *(r(n_r_LCR))**(D_l(st_map%lm2(l1,m1)))/ &
                  &           (r(nR))**(D_l(st_map%lm2(l1,m1))+2)*b(lm1,n_r_LCR)
                  aj(lm1,nR)=zero
                  dj(lm1,nR)=zero
                  ddj(lm1,nR)=zero
               end do
            end if
         end do
         !$omp end do
      end if

      !$omp end parallel


      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dbdt, llmMag, ulmMag, n_r_maxMag)
      call tscheme%rotate_imex(djdt, llmMag, ulmMag, n_r_maxMag)
      if ( l_cond_ic ) then
         call tscheme%rotate_imex(dbdt_ic, llmMag, ulmMag, n_r_ic_max)
         call tscheme%rotate_imex(djdt_ic, llmMag, ulmMag, n_r_ic_max)
      end if

   end subroutine updateB
!-----------------------------------------------------------------------------	
   subroutine finish_exp_mag_ic(b_ic, aj_ic, omega_ic, db_exp_last, dj_exp_last)

      !-- Input variables
      real(cp),    intent(in) :: omega_ic
      complex(cp), intent(in) :: aj_ic(llmMag:ulmMag,n_r_ic_max)
      complex(cp), intent(in) :: b_ic(llmMag:ulmMag,n_r_ic_max)

      !-- Output variables
      complex(cp), intent(inout) :: db_exp_last(llmMag:ulmMag,n_r_ic_max)
      complex(cp), intent(inout) :: dj_exp_last(llmMag:ulmMag,n_r_ic_max)

      !-- Local variables
      complex(cp) :: fac
      integer :: n_r, lm, l1, m1
      integer, pointer :: lm2l(:),lm2m(:)

      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      if ( omega_ic /= 0.0_cp .and. l_rot_ic .and. l_mag_nl ) then
         !$omp parallel do default(shared) private(lm,n_r,fac)
         do n_r=2,n_r_ic_max-1
            do lm=llmMag, ulmMag
               l1=lm2l(lm)
               m1=lm2m(lm)
               fac = -omega_ic*or2(n_r_max)*dPhi(st_map%lm2(l1,m1))* &
               &    dLh(st_map%lm2(l1,m1))
               db_exp_last(lm,n_r)=fac*b_ic(lm,n_r)
               dj_exp_last(lm,n_r)=fac*aj_ic(lm,n_r)
            end do
         end do
         !$omp end parallel do
      end if

   end subroutine finish_exp_mag_ic
!-----------------------------------------------------------------------------	
   subroutine finish_exp_mag(dVxBhLM, dj_exp_last)


      !-- Input variables
      complex(cp), intent(inout) :: dVxBhLM(llmMag:ulmMag,n_r_maxMag)

      !-- Output variables
      complex(cp), intent(inout) :: dj_exp_last(llmMag:ulmMag,n_r_maxMag)

      !-- Local variables
      integer :: n_r, lm, start_lm, stop_lm, lmStart_00

      lmStart_00 =max(2,llmMag)

      !$omp parallel default(shared) private(start_lm,stop_lm)
      start_lm=lmStart_00; stop_lm=ulmMag
      call get_openmp_blocks(start_lm, stop_lm)

      call get_dr( dVxBhLM,work_LMloc,ulmMag-llmMag+1,start_lm-llmMag+1, &
           &       stop_lm-llmMag+1,n_r_max,rscheme_oc,nocopy=.true. )
      !$omp barrier

      !$omp do private(n_r,lm)
      do n_r=1,n_r_max
         do lm=lmStart_00,ulmMag
            dj_exp_last(lm,n_r)=dj_exp_last(lm,n_r)+or2(n_r)*work_LMloc(lm,n_r)
         end do
      end do
      !$omp end do
      !$omp end parallel

   end subroutine finish_exp_mag
!-----------------------------------------------------------------------------	
   subroutine get_mag_ic_rhs_imp(b_ic, db_ic, ddb_ic, aj_ic, dj_ic, ddj_ic,  &
              &                  b_ic_last, aj_ic_last, db_ic_imp_last,      &
              &                  dj_ic_imp_last, l_calc_lin_rhs)


      !-- Input variables
      complex(cp), intent(in) :: b_ic(llmMag:ulmMag,n_r_ic_max)
      complex(cp), intent(in) :: aj_ic(llmMag:ulmMag,n_r_ic_max)
      logical,     intent(in) :: l_calc_lin_rhs

      !-- Output variable
      complex(cp), intent(out) :: db_ic(llmMag:ulmMag,n_r_ic_max)
      complex(cp), intent(out) :: ddb_ic(llmMag:ulmMag,n_r_ic_max)
      complex(cp), intent(out) :: dj_ic(llmMag:ulmMag,n_r_ic_max)
      complex(cp), intent(out) :: ddj_ic(llmMag:ulmMag,n_r_ic_max)
      complex(cp), intent(out) :: b_ic_last(llmMag:ulmMag,n_r_ic_max)
      complex(cp), intent(out) :: aj_ic_last(llmMag:ulmMag,n_r_ic_max)
      complex(cp), intent(out) :: db_ic_imp_last(llmMag:ulmMag,n_r_ic_max)
      complex(cp), intent(out) :: dj_ic_imp_last(llmMag:ulmMag,n_r_ic_max)

      !-- Local variables 
      integer :: l1, m1, lmStart_00
      integer :: n_r, lm, start_lm, stop_lm
      integer, pointer :: lm2l(:),lm2m(:)

      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m
      lmStart_00 =max(2,llmMag)

      !$omp parallel default(shared)  private(start_lm, stop_lm)
      start_lm=llmMag; stop_lm=ulmMag
      call get_openmp_blocks(start_lm,stop_lm)

      call get_ddr_even( b_ic,db_ic,ddb_ic, ulmMag-llmMag+1, &
           &             start_lm-llmMag+1,stop_lm-llmMag+1, &
           &             n_r_ic_max,n_cheb_ic_max, dr_fac_ic,&
           &             b_ic_last,aj_ic_last, chebt_ic, chebt_ic_even )
      call get_ddr_even( aj_ic,dj_ic,ddj_ic, ulmMag-llmMag+1,  &
           &             start_lm-llmMag+1, stop_lm-llmMag+1,  &
           &             n_r_ic_max,n_cheb_ic_max, dr_fac_ic,  &
           &             b_ic_last,aj_ic_last, chebt_ic, chebt_ic_even )

      !$omp do private(n_r,lm,l1,m1)
      do n_r=1,n_r_ic_max
         do lm=llmMag,ulmMag
            l1 = lm2l(lm)
            m1 = lm2m(lm)
            b_ic_last(lm,n_r) =dLh(st_map%lm2(l1,m1))*or2(n_r)*b_ic(lm,n_r)
            aj_ic_last(lm,n_r)=dLh(st_map%lm2(l1,m1))*or2(n_r)*aj_ic(lm,n_r)
         end do
      end do
      !$omp end do

      if ( l_calc_lin_rhs ) then
         !$omp do private(n_r,lm,l1,m1)
         do n_r=2,n_r_ic_max-1
            do lm=lmStart_00,ulmMag
               l1=lm2l(lm)
               m1=lm2m(lm)
               db_ic_imp_last(lm,n_r)=opm*O_sr*dLh(st_map%lm2(l1,m1))*   &
               &                or2(n_r_max) *  (   ddb_ic(lm,n_r) +    &
               &    two*D_lP1(st_map%lm2(l1,m1))*O_r_ic(n_r)*db_ic(lm,n_r) )
               dj_ic_imp_last(lm,n_r)=opm*O_sr*dLh(st_map%lm2(l1,m1))*   &
               &                or2(n_r_max) *  (   ddj_ic(lm,n_r) +    &
               &    two*D_lP1(st_map%lm2(l1,m1))*O_r_ic(n_r)*dj_ic(lm,n_r) )
            end do
         end do
         !$omp end do
         n_r=n_r_ic_max
         !$omp do private(lm,l1,m1)
         do lm=lmStart_00,ulmMag
            l1=lm2l(lm)
            m1=lm2m(lm)
            db_ic_imp_last(lm,n_r)=opm*O_sr*dLh(st_map%lm2(l1,m1))*   &
            &                                         or2(n_r_max) * &
            &    (one+two*D_lP1(st_map%lm2(l1,m1)))*ddb_ic(lm,n_r)
            dj_ic_imp_last(lm,n_r)=opm*O_sr*dLh(st_map%lm2(l1,m1))*   &
            &                                         or2(n_r_max) * &
            &    (one+two*D_lP1(st_map%lm2(l1,m1)))*ddj_ic(lm,n_r)
         end do
         !$omp end do
      end if

      !$omp end parallel

   end subroutine get_mag_ic_rhs_imp
!-----------------------------------------------------------------------------
   subroutine get_mag_rhs_imp(b, db, ddb, aj, dj, ddj,  b_last, aj_last, &
              &               db_imp_last, dj_imp_last, l_calc_lin_rhs,  &
              &               lRmsNext)


      !-- Input variables
      complex(cp), intent(in) :: b(llmMag:ulmMag,n_r_max)
      complex(cp), intent(in) :: aj(llmMag:ulmMag,n_r_max)
      logical,     intent(in) :: l_calc_lin_rhs
      logical,     intent(in) :: lRmsNext

      !-- Output variable
      complex(cp), intent(out) :: db(llmMag:ulmMag,n_r_max)
      complex(cp), intent(out) :: dj(llmMag:ulmMag,n_r_max)
      complex(cp), intent(out) :: ddj(llmMag:ulmMag,n_r_max)
      complex(cp), intent(out) :: ddb(llmMag:ulmMag,n_r_max)
      complex(cp), intent(out) :: b_last(llmMag:ulmMag,n_r_max)
      complex(cp), intent(out) :: aj_last(llmMag:ulmMag,n_r_max)
      complex(cp), intent(out) :: db_imp_last(llmMag:ulmMag,n_r_max)
      complex(cp), intent(out) :: dj_imp_last(llmMag:ulmMag,n_r_max)

      !-- Local variables 
      integer :: n_r_top, n_r_bot, l1, m1, lmStart_00
      integer :: n_r, lm, start_lm, stop_lm
      integer, pointer :: lm2l(:),lm2m(:)

      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m
      lmStart_00 =max(2,llmMag)

      !$omp parallel default(shared)  private(start_lm, stop_lm)
      start_lm=lmStart_00; stop_lm=ulmMag
      call get_openmp_blocks(start_lm,stop_lm)

      call get_ddr(b,db,ddb,ulmMag-llmMag+1,start_lm-llmMag+1, &
           &       stop_lm-llmMag+1,n_r_max,rscheme_oc)
      call get_ddr(aj,dj,ddj,ulmMag-llmMag+1,start_lm-llmMag+1, &
           &       stop_lm-llmMag+1,n_r_max,rscheme_oc)

      !$omp do private(n_r,lm,l1,m1)
      do n_r=1,n_r_max
         do lm=lmStart_00,ulmMag
            l1 = lm2l(lm)
            m1 = lm2m(lm)
            b_last(lm,n_r) =dLh(st_map%lm2(l1,m1))*or2(n_r)*b(lm,n_r)
            aj_last(lm,n_r)=dLh(st_map%lm2(l1,m1))*or2(n_r)*aj(lm,n_r)
         end do
      end do
      !$omp end do

      if ( l_calc_lin_rhs ) then
         if ( lRmsNext ) then
            n_r_top=n_r_cmb
            n_r_bot=n_r_icb
         else
            n_r_top=n_r_cmb+1
            n_r_bot=n_r_icb-1
         end if

         !$omp do private(n_r,lm,l1,m1,dtP,dtT)
         do n_r=n_r_top,n_r_bot
            do lm=lmStart_00,ulmMag
               l1=lm2l(lm)
               m1=lm2m(lm)
               db_imp_last(lm,n_r)= opm*lambda(n_r)*hdif_B(st_map%lm2(l1,m1))* &
               &                             dLh(st_map%lm2(l1,m1))*or2(n_r) * &
               &    ( ddb(lm,n_r) - dLh(st_map%lm2(l1,m1))*or2(n_r)*b(lm,n_r) )
               dj_imp_last(lm,n_r)= opm*lambda(n_r)*hdif_B(st_map%lm2(l1,m1))* &
               &                             dLh(st_map%lm2(l1,m1))*or2(n_r) * &
               &               ( ddj(lm,n_r) + dLlambda(n_r)*dj(lm,n_r) -      &
               &                  dLh(st_map%lm2(l1,m1))*or2(n_r)*aj(lm,n_r) )
               if ( lRmsNext ) then
                  !dtP(lm)=O_dt*dLh(st_map%lm2(l1,m1))*or2(n_r) &
                  !&             * (  b(lm,n_r)-work_LMloc(lm,n_r) )
                  !dtT(lm)=O_dt*dLh(st_map%lm2(l1,m1))*or2(n_r) &
                  !&             * ( aj(lm,n_r)-workB(lm,n_r) )
               end if
            end do
            if ( lRmsNext ) then
               call hInt2PolLM(dtP,llmMag,ulmMag,n_r,lmStart_00,ulmMag, &
                    &          dtBPolLMr(llmMag:ulmMag,n_r),            &
                    &          dtBPol2hInt(llmMag:ulmMag,n_r),lo_map)
               call hInt2TorLM(dtT,llmMag,ulmMag,n_r,lmStart_00,ulmMag, &
                    &          dtBTor2hInt(llmMag:ulmMag,n_r),lo_map)
            end if
         end do
         !$omp end do

      end if
      !$omp end parallel

   end subroutine get_mag_rhs_imp
!-----------------------------------------------------------------------------
#ifdef WITH_PRECOND_BJ
   subroutine get_bMat(tscheme,l,hdif,bMat,bMat_fac,jMat,jMat_fac)
#else
   subroutine get_bMat(tscheme,l,hdif,bMat,jMat)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matrices
      !  bmat(i,j) and ajmat for the dynamo equations.
      !

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme        ! time step
      integer,             intent(in) :: l
      real(cp),            intent(in) :: hdif

      !-- Output variables:
      class(type_realmat), intent(inout) :: bMat
      class(type_realmat), intent(inout) :: jMat
#ifdef WITH_PRECOND_BJ
      real(cp), intent(out) :: bMat_fac(n_r_totMag),jMat_fac(n_r_totMag)
#endif

      !-- local variables:
      integer :: nR,nCheb,nR_out,nRall
      integer :: info
      real(cp) :: l_P_1
      real(cp) :: dLh
      real(cp) :: rRatio
      real(cp) :: datJmat(n_r_tot,n_r_tot)
      real(cp) :: datBmat(n_r_tot,n_r_tot)

#undef MATRIX_CHECK
#ifdef MATRIX_CHECK
      integer :: i,j
      real(cp) :: rcond
      integer ::ipiv(n_r_tot),iwork(n_r_tot)
      real(cp) :: work(4*n_r_tot),anorm,linesum
      real(cp) :: temp_Mat(n_r_tot,n_r_tot)
      integer, save :: counter=0
      integer :: filehandle
      character(len=100) :: filename
#endif

      nRall=n_r_max
      if ( l_cond_ic ) nRall=nRall+n_r_ic_max
      dLh=real(l*(l+1),kind=cp)

      !-- matricies depend on degree l but not on order m,
      !   we thus have to construct bmat and ajmat for each l:

      l_P_1=real(l+1,kind=cp)

      do nR_out=1,n_r_max
         do nR=2,n_r_max-1
            datBmat(nR,nR_out)=                       rscheme_oc%rnorm * (  &
            &                 dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out) -      &
            &    tscheme%wimp_lin(1)*opm*lambda(nR)*hdif*dLh*or2(nR) * (    &
            &                                rscheme_oc%d2rMat(nR,nR_out) - &
            &                      dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out) ) )

            datJmat(nR,nR_out)=                       rscheme_oc%rnorm * (  &
            &                 dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out) -      &
            &   tscheme%wimp_lin(1)*opm*lambda(nR)*hdif*dLh*or2(nR) * (     &
            &                                rscheme_oc%d2rMat(nR,nR_out) + &
            &                    dLlambda(nR)*rscheme_oc%drMat(nR,nR_out) - &
            &                      dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out) ) )
         end do
      end do

      if  ( l_LCR ) then
         do nR=2,n_r_max-1
            if ( nR<=n_r_LCR ) then
               do nR_out=1,n_r_max
                  datBmat(nR,nR_out)= rscheme_oc%rnorm*(                    &
                  &                           rscheme_oc%drMat(nR,nR_out) + &
                  &    real(l,kind=cp)*or1(nR)*rscheme_oc%rMat(nR,nR_out) )

                  datJmat(nR,nR_out)= rscheme_oc%rnorm*rscheme_oc%rMat(nR,nR_out)
               end do
            end if
         end do
      end if

      !----- boundary conditions for outer core field:
      if ( ktopb == 1 .or. ktopb == 3 ) then
         !-------- at CMB (nR=1):
         !         the internal poloidal field should fit a potential
         !         field (matrix bmat) and the toroidal field has to
         !         vanish (matrix ajmat).

         datBmat(1,1:n_r_max)=    rscheme_oc%rnorm * (   &  
         &                      rscheme_oc%drMat(1,:) +  &
         &     real(l,cp)*or1(1)*rscheme_oc%rMat(1,:) +  &
         &                     conductance_ma* (         &
         &                     rscheme_oc%d2rMat(1,:) -  &
         &            dLh*or2(1)*rscheme_oc%rMat(1,:) ) )

         datJmat(1,1:n_r_max)=    rscheme_oc%rnorm * (   &
         &                       rscheme_oc%rMat(1,:) +  &
         &       conductance_ma*rscheme_oc%drMat(1,:) )
      else if ( ktopb == 2 ) then
         call abortRun('! Boundary condition ktopb=2 not defined!')
      else if ( ktopb == 4 ) then
         !----- pseudo vacuum condition, field has only
         !      a radial component, horizontal components
         !      vanish when aj and db are zero:
         datJmat(1,1:n_r_max)=rscheme_oc%rnorm* rscheme_oc%rMat(1,:)
         datBmat(1,1:n_r_max)=rscheme_oc%rnorm*rscheme_oc%drMat(1,:)
      end if

      !-------- at IC (nR=n_r_max):
      if ( l_full_sphere ) then
         datBmat(n_r_max,1:n_r_max)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
         datJmat(n_r_max,1:n_r_max)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
      else
         if ( kbotb == 1 ) then
            !----------- insulating IC, field has to fit a potential field:
            datJmat(n_r_max,1:n_r_max)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)

            datBmat(n_r_max,1:n_r_max)=rscheme_oc%rnorm * (      &
            &                      rscheme_oc%drMat(n_r_max,:) - &
            &    l_P_1*or1(n_r_max)*rscheme_oc%rMat(n_r_max,:) )
         else if ( kbotb == 2 ) then
            !----------- perfect conducting IC
            datBmat(n_r_max-1,1:n_r_max)=  &
            &        rscheme_oc%rnorm*rscheme_oc%d2rMat(n_r_max,:)
            datJmat(n_r_max,1:n_r_max)=rscheme_oc%rnorm* rscheme_oc%drMat(n_r_max,:)
         else if ( kbotb == 3 ) then
            !---------- finite conducting IC, four boundary conditions:
            !           continuity of b,j, (d b)/(d r) and (d j)/(d r)/sigma.
            !           note: n_r=n_r_max and n_r=n_r_max+1 stand for IC radius
            !           here we set the outer core part of the equations.
            !           the conductivity ratio sigma_ratio is used as
            !           an additional dimensionless parameter.
            datBmat(n_r_max,1:n_r_max)  =rscheme_oc%rnorm* rscheme_oc%rMat(n_r_max,:)
            datBmat(n_r_max+1,1:n_r_max)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
            datJmat(n_r_max,1:n_r_max)  =rscheme_oc%rnorm* rscheme_oc%rMat(n_r_max,:)
            datJmat(n_r_max+1,1:n_r_max)=rscheme_oc%rnorm*sigma_ratio* &
            &                      rscheme_oc%drMat(n_r_max,:)
         else if ( kbotb == 4 ) then
            !----- Pseudovacuum conduction at lower boundary:
            datJmat(n_r_max,1:n_r_max)= rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
            datBmat(n_r_max,1:n_r_max)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         end if
      end if

      !-------- Imposed fields: (overwrites above IC boundary cond.)
      if ( l == 1 .and. ( imagcon == -1 .or. imagcon == -2 ) ) then
         datBmat(n_r_max,1:n_r_max)=rscheme_oc%rnorm * rscheme_oc%rMat(n_r_max,:)
      else if ( l == 3 .and. imagcon == -10 ) then
         if ( l_LCR ) then
            call abortRun('Imposed field not compatible with weak conducting region!')
         end if
         datJmat(1,1:n_r_max)      =rscheme_oc%rnorm * rscheme_oc%rMat(1,:)
         datJmat(n_r_max,1:n_r_max)=rscheme_oc%rnorm * rscheme_oc%rMat(n_r_max,:)
      else if ( n_imp == 1 ) then
         !-- This is Uli Christensen's idea where the external field is
         !   not fixed but compensates the internal field so that the
         !   radial field component vanishes at r/r_cmb=rrMP
         if ( l_LCR ) then
            call abortRun('Imposed field not compatible with weak conducting region!')
         end if
         rRatio=rrMP**real(2*l+1,kind=cp)

         datBmat(1,1:n_r_max)=rscheme_oc%rnorm * (        &
         &                      rscheme_oc%drMat(1,:) +   &
         &     real(l,cp)*or1(1)*rscheme_oc%rMat(1,:) -   &
         &    real(2*l+1,cp)*or1(1)/(1-rRatio) +          &
         &                     conductance_ma* (          &
         &                     rscheme_oc%d2rMat(1,:) -   &
         &            dLh*or2(1)*rscheme_oc%rMat(1,:) ) )
      end if

      !----- fill up with zeros:
      do nR_out=rscheme_oc%n_max+1,n_r_max
         datBmat(1,nR_out)=0.0_cp
         datJmat(1,nR_out)=0.0_cp
         if ( l_LCR ) then
            do nR=2,n_r_LCR
               datBmat(nR,nR_out)=0.0_cp
               datJmat(nR,nR_out)=0.0_cp
            end do
         end if
         if ( kbotb == 1 ) then
            datBmat(n_r_max,nR_out)  =0.0_cp
            datJmat(n_r_max,nR_out)  =0.0_cp
         else if ( kbotb == 2 ) then
            datBmat(n_r_max-1,nR_out)=0.0_cp
            datJmat(n_r_max,nR_out)  =0.0_cp
         else if ( kbotb == 3 ) then
            datBmat(n_r_max,nR_out)  =0.0_cp
            datBmat(n_r_max+1,nR_out)=0.0_cp
            datJmat(n_r_max,nR_out)  =0.0_cp
            datJmat(n_r_max+1,nR_out)=0.0_cp
         else if ( kbotb == 4 ) then
            datBmat(n_r_max,nR_out)  =0.0_cp
            datJmat(n_r_max,nR_out)  =0.0_cp
         end if
      end do

      !----- normalization for highest and lowest Cheb mode:
      do nR=1,n_r_max
         datBmat(nR,1)      =rscheme_oc%boundary_fac*datBmat(nR,1)
         datBmat(nR,n_r_max)=rscheme_oc%boundary_fac*datBmat(nR,n_r_max)
         datJmat(nR,1)      =rscheme_oc%boundary_fac*datJmat(nR,1)
         datJmat(nR,n_r_max)=rscheme_oc%boundary_fac*datJmat(nR,n_r_max)
      end do
      if ( kbotb == 3 ) then
         datBmat(n_r_max+1,1)      =rscheme_oc%boundary_fac*datBmat(n_r_max+1,1)
         datBmat(n_r_max+1,n_r_max)=rscheme_oc%boundary_fac* &
         &                          datBmat(n_r_max+1,n_r_max)
         datJmat(n_r_max+1,1)      =rscheme_oc%boundary_fac*datJmat(n_r_max+1,1)
         datJmat(n_r_max+1,n_r_max)=rscheme_oc%boundary_fac* &
         &                          datJmat(n_r_max+1,n_r_max)
      end if

      !----- Conducting inner core:
      if ( l_cond_ic ) then
         !----- inner core implicit time step matricies for the grid
         !      points n_r=n_r_max+1,...,n_r_max+n_r_ic
         do nCheb=1,n_r_ic_max ! counts even IC cheb modes
            do nR=2,n_r_ic_max-1 ! counts IC radial grid points
               ! n_r=1 represents ICB
               !----------- poloidal field matrix for an inner core field
               !            of the radial form: (r/r_icb)**(l+1)*cheb_ic(r)
               !            where cheb_ic are even chebs only.
               !            NOTE: no hyperdiffusion in inner core !

               datBmat(n_r_max+nR,n_r_max+nCheb) =    &
               &    cheb_norm_ic*dLh*or2(n_r_max) * ( &
               &             cheb_ic(nCheb,nR) -      &
               &     tscheme%wimp_lin(1)*opm*O_sr * ( &
               &                d2cheb_ic(nCheb,nR) + &
               &    two*l_P_1*O_r_ic(nR)*dcheb_ic(nCheb,nR) )   )

               datJmat(n_r_max+nR,n_r_max+nCheb)=datBmat(n_r_max+nR,n_r_max+nCheb)
            end do

            !----- Special treatment for r=0, asymptotic of 1/r dr
            nR=n_r_ic_max
            datBmat(n_r_max+nR,n_r_max+nCheb) =    &
            &    cheb_norm_ic*dLh*or2(n_r_max) * ( &
            &                  cheb_ic(nCheb,nR) - &
            &       tscheme%wimp_lin(1)*opm*O_sr * &
            &    (one+two*l_P_1)*d2cheb_ic(nCheb,nR) )

            datJmat(n_r_max+nR,n_r_max+nCheb)=datBmat(n_r_max+nR,n_r_max+nCheb)
         end do

         !-------- boundary condition at r_icb:
         do nCheb=1,n_cheb_ic_max
            datBmat(n_r_max,n_r_max+nCheb)=-cheb_norm_ic*cheb_ic(nCheb,1)
            datBmat(n_r_max+1,n_r_max+nCheb)=              &
            &    -cheb_norm_ic * ( dcheb_ic(nCheb,1) +     &
            &    l_P_1*or1(n_r_max)*cheb_ic(nCheb,1) )
            datJmat(n_r_max,n_r_max+nCheb)  =datBmat(n_r_max,n_r_max+nCheb)
            datJmat(n_r_max+1,n_r_max+nCheb)=datBmat(n_r_max+1,n_r_max+nCheb)
         end do ! cheb modes

         !-------- fill with zeros:
         do nCheb=n_cheb_ic_max+1,n_r_ic_max
            datBmat(n_r_max,n_r_max+nCheb)  =0.0_cp
            datBmat(n_r_max+1,n_r_max+nCheb)=0.0_cp
            datJmat(n_r_max,n_r_max+nCheb)  =0.0_cp
            datJmat(n_r_max+1,n_r_max+nCheb)=0.0_cp
         end do

         !-------- normalization for lowest Cheb mode:
         do nR=n_r_max,n_r_tot
            datBmat(nR,n_r_max+1)=half*datBmat(nR,n_r_max+1)
            datJmat(nR,n_r_max+1)=half*datJmat(nR,n_r_max+1)
            datBmat(nR,n_r_tot)  =half*datBmat(nR,n_r_tot)
            datJmat(nR,n_r_tot)  =half*datJmat(nR,n_r_tot)
         end do

         !-------- fill matricies up with zeros:
         do nCheb=n_r_max+1,n_r_tot
            do nR=1,n_r_max-1
               datBmat(nR,nCheb)=0.0_cp
               datJmat(nR,nCheb)=0.0_cp
            end do
         end do
         do nCheb=1,n_r_max
            do nR=n_r_max+2,n_r_tot
               datBmat(nR,nCheb)=0.0_cp
               datJmat(nR,nCheb)=0.0_cp
            end do
         end do

      end if  ! conducting inner core ?

#ifdef WITH_PRECOND_BJ
      ! compute the linesum of each line
      do nR=1,nRall
         bMat_fac(nR)=one/maxval(abs(datBmat(nR,1:nRall)))
         datBmat(nR,:) = datBmat(nR,:)*bMat_fac(nR)
      end do
      do nR=1,nRall
         jMat_fac(nR)=one/maxval(abs(datJmat(nR,1:nRall)))
         datJmat(nR,:) = datJmat(nR,:)*jMat_fac(nR)
      end do
#endif

#ifdef MATRIX_CHECK
      ! copy the bMat to a temporary variable for modification
      write(filename,"(A,I3.3,A,I3.3,A)") "bMat_",l,"_",counter,".dat"
      open(newunit=filehandle,file=trim(filename))
      counter= counter+1

      do i=1,n_r_tot
         do j=1,n_r_tot
            write(filehandle,"(2ES20.12,1X)",advance="no") datBmat(i,j)
         end do
         write(filehandle,"(A)") ""
      end do
      close(filehandle)
      temp_Mat=datBmat
      anorm = 0.0_cp
      do i=1,n_r_tot
         linesum = 0.0_cp
         do j=1,n_r_tot
            linesum = linesum + abs(temp_Mat(i,j))
         end do
         if (linesum  >  anorm) anorm=linesum
      end do
      !write(*,"(A,ES20.12)") "anorm = ",anorm
      ! LU factorization
      call dgetrf(n_r_tot,n_r_tot,temp_Mat,n_r_tot,ipiv,info)
      ! estimate the condition number
      call dgecon('I',n_r_tot,temp_Mat,n_r_tot,anorm,rcond,work,iwork,info)
      write(*,"(A,I3,A,ES11.3)") "inverse condition number of bMat for l=",l," is ",rcond

      ! The same computation for jMat.
      ! copy the jMat to a temporary variable for modification
      write(filename,"(A,I3.3,A,I3.3,A)") "jMat_",l,"_",counter,".dat"
      open(newunit=filehandle,file=trim(filename))
      counter= counter+1

      do i=1,n_r_tot
         do j=1,n_r_tot
            write(filehandle,"(2ES20.12,1X)",advance="no") datJmat(i,j)
         end do
         write(filehandle,"(A)") ""
      end do
      close(filehandle)
      temp_Mat=datJmat
      anorm = 0.0_cp
      do i=1,n_r_tot
         linesum = 0.0_cp
         do j=1,n_r_tot
            linesum = linesum + abs(temp_Mat(i,j))
         end do
         if (linesum  >  anorm) anorm=linesum
      end do
      !write(*,"(A,ES20.12)") "anorm = ",anorm
      ! LU factorization
      call dgetrf(n_r_tot,n_r_tot,temp_Mat,n_r_tot,ipiv,info)
      ! estimate the condition number
      call dgecon('I',n_r_tot,temp_Mat,n_r_tot,anorm,rcond,work,iwork,info)
      write(*,"(A,I3,A,ES11.3)") "inverse condition number of jMat for l=",l," is ",rcond
#endif

      !-- Array copy
      call bMat%set_data(datBmat)
      call jMat%set_data(datJmat)

      !----- LU decomposition:
      call bMat%prepare(info)

      if ( info /= 0 ) call abortRun('Singular matrix bMat in get_bmat')

      !----- LU decomposition:
      call jMat%prepare(info)
      if ( info /= 0 ) call abortRun('! Singular matrix jMat in get_bmat!')

   end subroutine get_bMat
!-----------------------------------------------------------------------------
end module updateB_mod
