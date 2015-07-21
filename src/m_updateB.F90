!$Id$
#include "perflib_preproc.cpp"
module updateB_mod
   use omp_lib
   use truncation, only: n_r_max, n_r_tot, n_r_ic_max, n_cheb_max, &
                         n_cheb_ic_max, n_r_ic_maxMag, n_r_maxMag, &
                         n_r_totMag, lm_max
   use radial_functions, only: i_costf_init,d_costf_init,drx,ddrx,or2,r_cmb,&
                            & i_costf1_ic_init,d_costf1_ic_init,            &
                            & i_costf2_ic_init,d_costf2_ic_init,            &
                            & dr_fac_ic,lambda,dLlambda,o_r_ic,r,cheb_norm, &
                            & cheb, dcheb, d2cheb, or1, cheb_ic, dcheb_ic,  &
                            & d2cheb_ic, cheb_norm_ic
   use radial_data, only: n_r_cmb,n_r_icb
   use physical_parameters, only: n_r_LCR,opm,O_sr,kbotb, imagcon, tmagcon, &
                                 sigma_ratio, conductance_ma, ktopb, kbotb
   use init_fields, only: bpeaktop, bpeakbot
   use num_param, only: alpha
   use blocking, only: nLMBs,st_map,lo_map,st_sub_map,lo_sub_map,lmStartB,lmStopB
   use horizontal_data, only: dLh, dPhi, hdif_B, D_l, D_lP1
   use logic, only: l_cond_ic, l_LCR, l_rot_ic, l_mag_nl, l_b_nl_icb, &
                   l_b_nl_cmb, l_update_b
   use matrices, only: bPivot, jPivot, bMat, jMat, &
#ifdef WITH_PRECOND_BJ
                       bMat_fac, jMat_fac,         &
#endif
                       lBmat
   use RMS, only: dtBPolLMr, dtBPol2hInt, dtBPolAs2hInt, dtBTorAs2hInt, &
                 dtBTor2hInt
   use const, only: pi
   use Bext, only: n_imp, l_imp, bmax_imp, expo_imp, amp_imp, rrMP
#ifdef WITH_MKL_LU
   use lapack95, only: getrs, getrf
#else
   use algebra, only: cgeslML, sgefa
#endif
   use LMLoop_data, only: llmMag,ulmMag,llm_realMag,ulm_realMag
   use parallel_mod, only:  rank,chunksize
#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif

   implicit none

   private

   !-- Local work arrays:
   complex(kind=8), allocatable :: workA(:,:),workB(:,:)
   complex(kind=8), allocatable :: rhs1(:,:,:),rhs2(:,:,:)
   integer :: maxThreads

   public :: initialize_updateB,updateB

contains

   subroutine initialize_updateB

      allocate( workA(llmMag:ulmMag,n_r_max) )
      allocate( workB(llmMag:ulmMag,n_r_max) )
#ifdef WITHOMP
      maxThreads=omp_get_max_threads()
#else
      maxThreads=1
#endif

      allocate(rhs1(2*n_r_max,lo_sub_map%sizeLMB2max,0:maxThreads-1))
      allocate(rhs2(2*n_r_max,lo_sub_map%sizeLMB2max,0:maxThreads-1))

   end subroutine initialize_updateB
!-----------------------------------------------------------------------------
   subroutine updateB(b,db,ddb,aj,dj,ddj,dVxBhLM, &
       &             dbdt,dbdtLast,djdt,djdtLast, &
       &             b_ic,db_ic,ddb_ic,aj_ic,dj_ic,ddj_ic, &
       &             dbdt_icLast,djdt_icLast, &
       &             b_nl_cmb,aj_nl_cmb,aj_nl_icb,omega_ic, &
       &             w1,coex,dt,time,nLMB,lRmsNext)
      !  +-------------------------------------------------------------------+
      !  |                                                                   |
      !  |  Calculated update of magnetic field potential and the time       |
      !  |  stepping arrays dbdtLast, ...                                    |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+
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

      !-----------------------------------------------------------------------

      !-- Input variables:
      complex(kind=8), intent(in) :: b_nl_cmb(:)  ! nonlinear BC for b at CMB
      complex(kind=8), intent(in) :: aj_nl_cmb(:) ! nonlinear BC for aj at CMB
      complex(kind=8), intent(in) :: aj_nl_icb(:) ! nonlinear BC for dr aj at ICB
      complex(kind=8), intent(in) :: dVxBhLM(llmMag:ulmMag,n_r_maxMag)
      complex(kind=8), intent(in) :: dbdt(llmMag:ulmMag,n_r_maxMag)
      real(kind=8),    intent(in) :: omega_ic
      real(kind=8),    intent(in) :: w1    ! weight for time step !
      real(kind=8),    intent(in) :: coex  ! factor depending on alpha
      real(kind=8),    intent(in) :: dt
      real(kind=8),    intent(in) :: time
      integer,         intent(in) :: nLMB
      logical,         intent(in) :: lRmsNext

      !-- Input/output of scalar potentials and time stepping arrays:
      complex(kind=8), intent(inout) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(kind=8), intent(inout) :: aj(llmMag:ulmMag,n_r_maxMag)
      complex(kind=8), intent(inout) :: dbdtLast(llmMag:ulmMag,n_r_maxMag)
      complex(kind=8), intent(inout) :: djdt(llmMag:ulmMag,n_r_maxMag)
      complex(kind=8), intent(inout) :: djdtLast(llmMag:ulmMag,n_r_maxMag)
      complex(kind=8), intent(inout) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(kind=8), intent(inout) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(kind=8), intent(inout) :: dbdt_icLast(llmMag:ulmMag,n_r_ic_maxMag)
      complex(kind=8), intent(inout) :: djdt_icLast(llmMag:ulmMag,n_r_ic_maxMag)

      !-- Output variables:
      complex(kind=8), intent(out) :: db(llmMag:ulmMag,n_r_maxMag)
      complex(kind=8), intent(out) :: ddb(llmMag:ulmMag,n_r_maxMag)
      complex(kind=8), intent(out) :: dj(llmMag:ulmMag,n_r_maxMag)
      complex(kind=8), intent(out) :: ddj(llmMag:ulmMag,n_r_maxMag)
      complex(kind=8),intent(out) :: db_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(kind=8),intent(out) :: ddb_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(kind=8),intent(out) :: dj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(kind=8),intent(out) :: ddj_ic(llmMag:ulmMag,n_r_ic_maxMag)

      !-- Local variables:
      real(kind=8) :: w2             ! weight of second time step
      real(kind=8) :: O_dt
      real(kind=8) :: yl0_norm,prefac!External magnetic field of general l

      integer :: l1,m1               ! degree and order
      integer :: lm1,lm,lmB          ! position of (l,m) in array
      integer :: lmStart,lmStop      ! max and min number of orders m
      integer :: lmStart_real        ! range of lm for real array
      integer :: lmStop_real
      integer :: lmStart_00          ! excluding l=0,m=0
      integer :: nLMB2
      integer :: n_cheb              ! No of cheb polynome (degree+1)
      integer :: nR                 ! No of radial grid point
      integer :: n_r_real            ! total number of used grid points

      complex(kind=8) :: fac
      complex(kind=8) :: dbdt_ic,djdt_ic  ! they are calculated here !

      integer, pointer :: nLMBs2(:),lm2l(:),lm2m(:)
      integer, pointer :: sizeLMB2(:,:),lm2(:,:)
      integer, pointer :: lm22lm(:,:,:),lm22l(:,:,:),lm22m(:,:,:)

      !-- for feedback
      real(kind=8) :: ff,cimp,aimp,b10max

      real(kind=8), save :: direction

      integer :: iThread,start_lm,stop_lm,all_lms,per_thread,nThreads,maxThreads
      integer :: nChunks,iChunk,lmB0,size_of_last_chunk,threadid
      
      if ( .not. l_update_b ) RETURN

#ifdef WITHOMP
      maxThreads=omp_get_max_threads()
#else
      maxThreads = 1
#endif

      nLMBs2(1:nLMBs) => lo_sub_map%nLMBs2
      sizeLMB2(1:,1:) => lo_sub_map%sizeLMB2
      lm22lm(1:,1:,1:) => lo_sub_map%lm22lm
      lm22l(1:,1:,1:) => lo_sub_map%lm22l
      lm22m(1:,1:,1:) => lo_sub_map%lm22m
      lm2(0:,0:) => lo_map%lm2
      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      n_r_real=n_r_max
      if ( l_cond_ic ) n_r_real=n_r_max+n_r_ic_max

      lmStart     =lmStartB(nLMB)
      lmStop      =lmStopB(nLMB)
      lmStart_00  =max(2,lmStart)
      lmStart_real=2*lmStart_00-1
      lmStop_real =2*lmStop

      ! output the input 
      !write(*,"(4(A,I4),I4,A,I4)") "nLMB=",nLMB,", from ",lmStart," to ",lmStop,&
      !     &", reals: ",lmStart_real,lmStop_real,", nLMBs2 = ",nLMBs2(nLMB)

      w2  =1.D0-w1
      O_dt=1.D0/dt

      !--- Start with finishing djdt:
      !    dVxBhLM is still in the R-distributed space,
      !    the ouput workA is in the LM-distributed space.
      !if (2*lmStart-1 - llm_realMag+1.NE.1) then
      !   write(*,"(I4,2(A,I6))") rank,": lmStart = ",lmStart,", llm_realMag = ",llm_realMag
      !   STOP
      !end if
      !if (lmStop_real .NE. ulm_realMag) then
      !   write(*,"(I4,A,2I6)") rank,": ",ulm_realMag,lmStop_real
      !   stop
      !end if
      !call get_drNS( dVxBhLM,workA, &
      !     &         ulm_realMag-llm_realMag+1,(2*lmStart-1)-llm_realMag+1,lmStop_real-llm_realMag+1, &
      !     &         n_r_max,n_cheb_max,workB, &
      !     &         i_costf_init,d_costf_init,drx)
      ! simplified interface
      !PRINT*,rank,": computing for ",ulm_realMag-llm_realMag+1," rows, i_costf_init = ",i_costf_init

      !PERFON('upB_fin')
      all_lms=lmStop_real-lmStart_real+1
#ifdef WITHOMP
      if (all_lms < maxThreads) then
         call omp_set_num_threads(all_lms)
      end if
#endif
      !$OMP PARALLEL default(none) &
      !$OMP private(iThread,start_lm,stop_lm) &
      !$OMP shared(all_lms,per_thread,lmStart_real,lmStop_real,lmStart_00,lmStop) &
      !$OMP shared(dVxBhLM,workA,workB,djdt,or2) &
      !$OMP shared(i_costf_init,d_costf_init,drx) &
      !$OMP shared(n_r_max,n_cheb_max,nThreads,llm_realMag,ulm_realMag)
      !$OMP SINGLE
#ifdef WITHOMP
      nThreads=omp_get_num_threads()
#else 
      nThreads=1
#endif
      !-- Get radial derivatives of s: workA,dsdtLast used as work arrays
      per_thread=all_lms/nThreads
      !$OMP END SINGLE
      !$OMP BARRIER
      !$OMP DO
      do iThread=0,nThreads-1
         start_lm=lmStart_real+iThread*per_thread
         stop_lm = start_lm+per_thread-1
         if (iThread == nThreads-1) stop_lm=lmStop_real

         call get_drNS( dVxBhLM,workA, &
              &         ulm_realMag-llm_realMag+1,start_lm-llm_realMag+1,stop_lm-llm_realMag+1, &
              &         n_r_max,n_cheb_max,workB, &
              &         i_costf_init,d_costf_init,drx)

      end do
      !$OMP end do

      !$OMP do private(nR)
      do nR=1,n_r_max
         do lm=lmStart_00,lmStop
            djdt(lm,nR)=djdt(lm,nR)+or2(nR)*workA(lm,nR)
         end do
      end do
      !$OMP end do
      !$OMP END PARALLEL
#ifdef WITHOMP
      call omp_set_num_threads(maxThreads)
#endif
      !PERFOFF

      ! This is a loop over all l values which should be treated on 
      ! the actual MPI rank
      !$OMP PARALLEL default(shared)
      !$OMP SINGLE
      do nLMB2=1,nLMBs2(nLMB)
         !$OMP TASK default(shared) &
         !$OMP firstprivate(nLMB2) &
         !$OMP private(lmB,lm,lm1,l1,m1,nR,iChunk,nChunks,size_of_last_chunk) &
         !$OMP private(dbdt_ic,djdt_ic,fac,bpeaktop,ff,cimp,aimp,threadid)

         ! determine the number of chunks of m
         ! total number for l1 is sizeLMB2(nLMB2,nLMB)
         ! chunksize is given
         nChunks = (sizeLMB2(nLMB2,nLMB)+chunksize-1)/chunksize
         size_of_last_chunk = chunksize + (sizeLMB2(nLMB2,nLMB)-nChunks*chunksize)

         l1=lm22l(1,nLMB2,nLMB)
         if ( l1 > 0 ) then
            if ( .not. lBmat(l1) ) then
#ifdef WITH_PRECOND_BJ
               call get_bMat(dt,l1,hdif_B(st_map%lm2(l1,0)),           &
                             bMat(1,1,l1),bPivot(1,l1), bMat_fac(1,l1),&
                             jMat(1,1,l1),jPivot(1,l1), jMat_fac(1,l1))
#else
               call get_bMat(dt,l1,hdif_B(st_map%lm2(l1,0)), &
                             bMat(1,1,l1),bPivot(1,l1),      &
                             jMat(1,1,l1),jPivot(1,l1) )
#endif
               lBmat(l1)=.TRUE.
            end if
         end if

         do iChunk=1,nChunks
            !$OMP TASK if (nChunks>1) default(shared) &
            !$OMP firstprivate(iChunk) &
            !$OMP private(lmB0,lmB,lm,lm1,m1,nR,n_cheb) &
            !$OMP private(dbdt_ic,djdt_ic,fac,bpeaktop,ff) &
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
                     rhs1(1,lmB,threadid) = 0.D0
                     rhs2(1,lmB,threadid) = 0.D0
                  end if

                  rhs1(n_r_max,lmB,threadid)=0.D0
                  if ( kbotb == 2 ) rhs1(n_r_max-1,lmB,threadid)=0.D0

                  rhs2(n_r_max,lmB,threadid)=0.D0
                  if ( m1 == 0 ) then   ! Magnetoconvection boundary conditions
                     if ( imagcon /= 0 .and. tmagcon <= time ) then
                        if ( l_LCR ) then
                           write(*,*) 'LCR not compatible with imposed field!'
                           stop
                        end if
                        if ( l1 == 2 .and. imagcon > 0 .and. imagcon .NE. 12 ) then
                           rhs2(1,lmB,threadid)      =cmplx(bpeaktop,0.D0,kind=kind(0d0))
                           rhs2(n_r_max,lmB,threadid)=cmplx(bpeakbot,0.D0,kind=kind(0d0))
                        else if( l1 == 1 .and. imagcon == 12 ) then
                           rhs2(1,lmB,threadid)      =cmplx(bpeaktop,0.D0,kind=kind(0d0))
                           rhs2(n_r_max,lmB,threadid)=cmplx(bpeakbot,0.D0,kind=kind(0d0))
                        else if( l1 == 1 .and. imagcon == -1) then
                           rhs1(n_r_max,lmB,threadid)=cmplx(bpeakbot,0.D0,kind=kind(0d0))
                        else if( l1 == 1 .and. imagcon == -2) then
                           rhs1(1,lmB,threadid)      =cmplx(bpeaktop,0.D0,kind=kind(0d0))
                        else if( l1 == 3 .and. imagcon == -10 ) then
                           rhs2(1,lmB,threadid)      =cmplx(bpeaktop,0.D0,kind=kind(0d0))
                           rhs2(n_r_max,lmB,threadid)=cmplx(bpeakbot,0.D0,kind=kind(0d0))
                        end if
                     end if
                     if ( n_imp > 1 .and. l1 == l_imp ) then
                         if ( l_LCR ) then
                            write(*,*) 'LCR not compatible with imposed field!'
                            stop
                         end if
                         ! General normalization for degree l and order 0
                         yl0_norm = 0.5D0*dsqrt((2*l1+1)/pi)   
                         ! Prefactor for CMB matching condition
                         prefac = dble(2*l1+1)/dble(l1*(l1+1)) 

                        if ( n_imp == 2 ) then
                           !  Chose external field coefficient so that amp_imp is
                           !  the amplitude of the external field:
                           bpeaktop=prefac*r_cmb/yl0_norm*amp_imp
                           rhs1(1,lmB,threadid)=cmplx(bpeaktop,0.D0,kind=kind(0d0))
                        else if ( n_imp == 3 ) then
                           !  Chose external field coefficient so that amp_imp is
                           !  the amplitude of the external field:
                           bpeaktop=prefac*r_cmb/yl0_norm*amp_imp
                           if ( real(b(2,1)) > 1.D-9 ) &
                                direction=real(b(2,1))/DABS(real(b(2,1)))
                           rhs1(1,lmB,threadid)=cmplx(bpeaktop,0.D0, &
                                                kind=kind(0d0))*direction
                        else if ( n_imp == 4 ) then
                           !  I have forgotten what this was supposed to do:
                           bpeaktop=3.D0/r_cmb*amp_imp*real(b(2,1))**2
                           rhs1(1,lmB,threadid)=cmplx(bpeaktop,0.D0, &
                                                kind=kind(0d0))/b(2,1)

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

                           ff=0.D0
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
                                   (cimp*(expo_imp-1))**((expo_imp-1)/expo_imp)

                              ff=  aimp*real(b(2,1))**expo_imp/ &
                                   (cimp+real(b(2,1))**expo_imp)

                           end if
                           rhs1(1,lmB,threadid)=(2*l1+1)/r_cmb*ff

                        end if
                     end if
                  end if
                  
                  do nR=2,n_r_max-1
                     if ( nR<=n_r_LCR ) then
                        rhs1(nR,lmB,threadid)=0.D0
                        rhs2(nR,lmB,threadid)=0.D0
                     else
                        rhs1(nR,lmB,threadid)= ( w1*dbdt(lm1,nR)+w2*dbdtLast(lm1,nR) ) &
                          &        + O_dt*dLh(st_map%lm2(l1,m1))*or2(nR)*b(lm1,nR)
                        rhs2(nR,lmB,threadid)= ( w1*djdt(lm1,nR)+w2*djdtLast(lm1,nR) ) &
                          &        + O_dt*dLh(st_map%lm2(l1,m1))*or2(nR)*aj(lm1,nR)
                     end if
                  end do

                  !-- Magnetic boundary conditions, inner core for radial derivatives
                  !         of poloidal and toroidal magnetic potentials:
                  if ( l_cond_ic ) then    ! inner core
                     rhs1(n_r_max+1,lmB,threadid)=0.d0
                     if ( l_b_nl_icb ) then
                        rhs2(n_r_max+1,lmB,threadid)=aj_nl_icb(st_map%lm2(l1,m1))
                     else
                        rhs2(n_r_max+1,lmB,threadid)=0.d0
                     end if

                     do nR=2,n_r_ic_max
                        if ( omega_ic == 0.D0 .or. .not. l_rot_ic .or. &
                             .not. l_mag_nl ) then
                           dbdt_ic=cmplx(0.D0,0.D0,kind=kind(0d0))
                           djdt_ic=cmplx(0.D0,0.D0,kind=kind(0d0))
                        else
                           fac=-omega_ic*or2(n_r_max)*dPhi(st_map%lm2(l1,m1))* &
                                dLh(st_map%lm2(l1,m1))
                           dbdt_ic=fac*b_ic(lm1,nR)
                           djdt_ic=fac*aj_ic(lm1,nR)
                        end if
                        rhs1(n_r_max+nR,lmB,threadid)=                 &
                             & ( w1*dbdt_ic + w2*dbdt_icLast(lm1,nR) ) &
                             & +O_dt*dLh(st_map%lm2(l1,m1))*or2(n_r_max) * b_ic(lm1,nR)
                        rhs2(n_r_max+nR,lmB,threadid)=                 &
                             & ( w1*djdt_ic + w2*djdt_icLast(lm1,nR) ) &
                             & +O_dt*dLh(st_map%lm2(l1,m1))*or2(n_r_max) * aj_ic(lm1,nR)

                        !--------- Store the IC non-linear terms for the usage below:
                        dbdt_icLast(lm1,nR)=dbdt_ic
                        djdt_icLast(lm1,nR)=djdt_ic
                     end do
                  end if

               end if ! l>0
            end do    ! loop over lm in block

            if ( lmB > lmB0 ) then
               !write(*,"(2(A,I5))") "updateB: Calling cgeslML for l1=",l1," WITH lmB=",lmB
#ifdef WITH_PRECOND_BJ
               do lm=lmB0+1,lmB
                  do nR=1,n_r_tot
                     rhs1(nR,lm,threadid)=rhs1(nR,lm,threadid)*bMat_fac(nR,l1)
                     rhs2(nR,lm,threadid)=rhs2(nR,lm,threadid)*jMat_fac(nR,l1)
                  end do
               end do
#endif

               !LIKWID_ON('upB_sol')
#ifdef WITH_MKL_LU
               call getrs(cmplx(bMat(1:n_r_real,1:n_r_real,l1),0.D0,kind=kind(0.D0)), &
                    bPivot(1:n_r_real,l1),rhs1(1:n_r_real,lmB0+1:lmB,threadid))
               call getrs(cmplx(jMat(1:n_r_real,1:n_r_real,l1),0.D0,kind=kind(0.D0)), &
                    jPivot(1:n_r_real,l1),rhs2(1:n_r_real,lmB0+1:lmB,threadid))
#else
               call cgeslML(bMat(:,:,l1),n_r_tot,n_r_real, &
                    bPivot(:,l1),rhs1(:,lmB0+1:lmB,threadid),2*n_r_max,lmB-lmB0)
               call cgeslML(jMat(:,:,l1),n_r_tot,n_r_real, &
                    jPivot(:,l1),rhs2(:,lmB0+1:lmB,threadid),2*n_r_max,lmB-lmB0)
#endif
               !LIKWID_OFF('upB_sol')
            end if

            if ( lRmsNext ) then ! Store old b,aj
               do nR=1,n_r_max
                  do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
                     !do lm=1,sizeLMB2(nLMB2,nLMB)
                     lm1=lm22lm(lm,nLMB2,nLMB)
                     workA(lm1,nR)=b(lm1,nR)
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
                     do n_cheb=1,n_cheb_max  ! outer core
                        b(lm1,n_cheb) =rhs1(n_cheb,lmB,threadid)
                        aj(lm1,n_cheb)=rhs2(n_cheb,lmB,threadid)
                     end do
                     if ( l_cond_ic ) then   ! inner core
                        do n_cheb=1,n_cheb_ic_max
                           b_ic(lm1,n_cheb) = rhs1(n_r_max+n_cheb,lmB,threadid)
                           aj_ic(lm1,n_cheb)= rhs2(n_r_max+n_cheb,lmB,threadid)
                        end do
                     end if
                  else
                     do n_cheb=1,n_cheb_max   ! outer core
                        b(lm1,n_cheb) = &
                             cmplx(real(rhs1(n_cheb,lmB,threadid)),0.D0,kind=kind(0d0))
                        aj(lm1,n_cheb)= &
                             cmplx(real(rhs2(n_cheb,lmB,threadid)),0.D0,kind=kind(0d0))
                     end do
                     if ( l_cond_ic ) then    ! inner core
                        do n_cheb=1,n_cheb_ic_max
                           b_ic(lm1,n_cheb)= cmplx( &
                             real(rhs1(n_r_max+n_cheb,lmB,threadid)),0.D0,kind=kind(0d0))
                           aj_ic(lm1,n_cheb)= cmplx( &
                             real(rhs2(n_r_max+n_cheb,lmB,threadid)),0.D0,kind=kind(0d0))
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
      !$OMP END PARALLEL

      !-- Set cheb modes > n_cheb_max to zero (dealiazing)
      !   for inner core modes > 2*n_cheb_ic_max = 0
      do n_cheb=n_cheb_max+1,n_r_max
         do lm1=lmStart_00,lmStop
            b(lm1,n_cheb) =cmplx(0.D0,0.D0,kind=kind(0d0))
            aj(lm1,n_cheb)=cmplx(0.D0,0.D0,kind=kind(0d0))
         end do
      end do
      if ( l_cond_ic ) then
         do n_cheb=n_cheb_ic_max+1,n_r_ic_max
            do lm1=lmStart_00,lmStop
               b_ic(lm1,n_cheb) =cmplx(0.D0,0.D0,kind=kind(0d0))
               aj_ic(lm1,n_cheb)=cmplx(0.D0,0.D0,kind=kind(0d0))
            end do
         end do
      end if

      !PERFON('upB_drv')
      all_lms=lmStop_real-lmStart_real+1
#ifdef WITHOMP
      if (all_lms < maxThreads) then
         call omp_set_num_threads(all_lms)
      end if
#endif
      !$OMP PARALLEL default(none) &
      !$OMP private(iThread,start_lm,stop_lm) &
      !$OMP shared(all_lms,per_thread,lmStart_real,lmStop_real) &
      !$OMP shared(b,db,ddb,aj,dj,ddj,dbdtLast,djdtLast) &
      !$OMP shared(i_costf_init,d_costf_init,drx,ddrx) &
      !$OMP shared(n_r_max,n_cheb_max,nThreads,llm_realMag,ulm_realMag) &
      !$OMP shared(l_cond_ic,b_ic,db_ic,ddb_ic,aj_ic,dj_ic,ddj_ic) &
      !$OMP shared(i_costf1_ic_init,d_costf1_ic_init,i_costf2_ic_init,d_costf2_ic_init) &
      !$OMP shared(n_r_ic_max,n_cheb_ic_max,dr_fac_ic)
      !$OMP SINGLE
#ifdef WITHOMP
      nThreads=omp_get_num_threads()
#else
      nThreads=1
#endif
      !-- Get radial derivatives of s: workA,dsdtLast used as work arrays
      per_thread=all_lms/nThreads
      !write(*,"(2(A,I5))") "nThreads = ",nThreads,", per_thread = ",per_thread
      !$OMP END SINGLE
      !$OMP BARRIER
      !$OMP DO
      do iThread=0,nThreads-1
         start_lm=lmStart_real+iThread*per_thread
         stop_lm = start_lm+per_thread-1
         if (iThread == nThreads-1) stop_lm=lmStop_real

         !-- Radial derivatives: dbdtLast and djdtLast used as work arrays
         !PERFON('upB_cb')
         call costf1(b,ulm_realMag-llm_realMag+1,&
              &      start_lm-llm_realMag+1,stop_lm-llm_realMag+1, &
              &      dbdtLast,i_costf_init,d_costf_init)
         !PERFOFF
         !PERFON('upB_db')
         call get_ddr(b,db,ddb,ulm_realMag-llm_realMag+1,&
              &       start_lm-llm_realMag+1,stop_lm-llm_realMag+1, &
              &       n_r_max,n_cheb_max,dbdtLast,djdtLast, &
              &       i_costf_init,d_costf_init,drx,ddrx)
         !PERFOFF
         call costf1(aj,ulm_realMag-llm_realMag+1,&
              &      start_lm-llm_realMag+1,stop_lm-llm_realMag+1, &
              &      dbdtLast,i_costf_init,d_costf_init)
         call get_ddr(aj,dj,ddj,ulm_realMag-llm_realMag+1,&
              &       start_lm-llm_realMag+1,stop_lm-llm_realMag+1, &
              &       n_r_max,n_cheb_max,dbdtLast,djdtLast, &
              &       i_costf_init,d_costf_init,drx,ddrx)
         
         !-- Same for inner core:
         if ( l_cond_ic ) then
            call costf1(b_ic,ulm_realMag-llm_realMag+1,&
                 &      start_lm-llm_realMag+1,stop_lm-llm_realMag+1, &
                 &      dbdtLast,i_costf1_ic_init,d_costf1_ic_init)
            call get_ddr_even( b_ic,db_ic,ddb_ic,                     &
                 & ulm_realMag-llm_realMag+1,start_lm-llm_realMag+1,  &
                 & stop_lm-llm_realMag+1, n_r_ic_max,n_cheb_ic_max,   &
                 & dr_fac_ic,dbdtLast,djdtLast, i_costf1_ic_init,     &
                 & d_costf1_ic_init, i_costf2_ic_init,d_costf2_ic_init)
            call costf1(aj_ic,ulm_realMag-llm_realMag+1,&
                 & start_lm-llm_realMag+1,stop_lm-llm_realMag+1, &
                 & dbdtLast,i_costf1_ic_init,d_costf1_ic_init)
            call get_ddr_even( aj_ic,dj_ic,ddj_ic,                   &
                 & ulm_realMag-llm_realMag+1,start_lm-llm_realMag+1, &
                 & stop_lm-llm_realMag+1, n_r_ic_max,n_cheb_ic_max,  &
                 & dr_fac_ic,dbdtLast,djdtLast, i_costf1_ic_init,    &
                 & d_costf1_ic_init, i_costf2_ic_init,d_costf2_ic_init)
         end if
      end do
      !$OMP end do
      !$OMP END PARALLEL
#ifdef WITHOMP
      call omp_set_num_threads(maxThreads)
#endif
      !PERFOFF
      !-- We are now back in radial space !

      !PERFON('upB_last')
      if ( l_LCR ) then
         do nR=n_r_cmb,n_r_icb-1
            if ( nR<=n_r_LCR ) then
               do lm1=lmStart_00,lmStop
                  l1=lm2l(lm1)
                  m1=lm2m(lm1)

                  b(lm1,nR)=(r(n_r_LCR)/r(nR))**D_l(st_map%lm2(l1,m1))*b(lm1,n_r_LCR)
                  db(lm1,nR)=-dble(D_l(st_map%lm2(l1,m1)))*      &
                           (r(n_r_LCR))**D_l(st_map%lm2(l1,m1))/ &
                           (r(nR))**(D_l(st_map%lm2(l1,m1))+1)*b(lm1,n_r_LCR)
                  ddb(lm1,nR)=dble(D_l(st_map%lm2(l1,m1)))*         &
                             (dble(D_l(st_map%lm2(l1,m1)))+1)       &
                           *(r(n_r_LCR))**(D_l(st_map%lm2(l1,m1)))/ &
                           (r(nR))**(D_l(st_map%lm2(l1,m1))+2)*b(lm1,n_r_LCR)
                  aj(lm1,nR)=cmplx(0.D0,0.D0,kind=kind(0d0))
                  dj(lm1,nR)=cmplx(0.D0,0.D0,kind=kind(0d0))
                  ddj(lm1,nR)=cmplx(0.D0,0.D0,kind=kind(0d0))
               end do
            end if
         end do
      end if

      do nR=n_r_cmb,n_r_icb-1
         do lm1=lmStart_00,lmStop
            l1=lm2l(lm1)
            m1=lm2m(lm1)
            dbdtLast(lm1,nR)= dbdt(lm1,nR) -                    &
                 coex*opm*lambda(nR)*hdif_B(st_map%lm2(l1,m1))* &
                               dLh(st_map%lm2(l1,m1))*or2(nR) * &
                 ( ddb(lm1,nR) - dLh(st_map%lm2(l1,m1))*or2(nR)*b(lm1,nR) )
            djdtLast(lm1,nR)= djdt(lm1,nR) -                    &
                 coex*opm*lambda(nR)*hdif_B(st_map%lm2(l1,m1))* &
                               dLh(st_map%lm2(l1,m1))*or2(nR) * &
                 ( ddj(lm1,nR) + dLlambda(nR)*dj(lm1,nR) -      &
                   dLh(st_map%lm2(l1,m1))*or2(nR)*aj(lm1,nR) )
            if ( lRmsNext ) then
               workA(lm1,nR)=O_dt*dLh(st_map%lm2(l1,m1))*or2(nR) &
                             * (  b(lm1,nR)-workA(lm1,nR) )
               workB(lm1,nR)=O_dt*dLh(st_map%lm2(l1,m1))*or2(nR) &
                             * ( aj(lm1,nR)-workB(lm1,nR) )
            end if
         end do
         if ( lRmsNext ) then
            !write(*,"(A,2I3,2ES20.12)") "workA = ",nLMB,nR,SUM( workA(lmStart_00:lmStop,nR) )
            !call hInt2Pol(workA(1,nR),nR,lmStart_00,lmStop,dtBPolLMr, &
            !     dtBPol2hInt(nR,nTh),dtBPolAs2hInt(nR,nTh),lo_map)
            !write(*,"(A,2I3,ES20.13)") "upB: before dtBPol2hInt = ",nLMB,nR,dtBPol2hInt(nR,1)
            !call hInt2Pol(workA(llmMag,nR),ulmMag-llmMag+1,nR,lmStart_00-llmMag+1,lmStop-llmMag+1,&
            !     &dtBPolLMr,dtBPol2hInt(nR,nTh),dtBPolAs2hInt(nR,nTh),lo_map)
            call hInt2Pol(workA(llmMag,nR),llmMag,ulmMag,nR,lmStart_00,lmStop,&
                 &dtBPolLMr,dtBPol2hInt(nR,1),dtBPolAs2hInt(nR,1),lo_map)
            !write(*,"(A,2I3,ES20.13)") "upB: after  dtBPol2hInt = ",nLMB,nR,dtBPol2hInt(nR,1)

            !write(*,"(A,2I3,3ES20.13)") "upB: before dtBTor2hInt = ",nLMB,nR,dtBTor2hInt(nR,1),SUM( workB(lmStart:lmStop,nR) )
            !call hInt2Tor(workB(llmMag,nR),ulmMag-llmMag+1,nR,lmStart_00-llmMag+1,lmStop-llmMag+1, &
            !     &        dtBTor2hInt(nR,nTh),dtBTorAs2hInt(nR,nTh))
            call hInt2Tor(workB(llmMag,nR),llmMag,ulmMag,nR,lmStart_00,lmStop, &
                 &        dtBTor2hInt(nR,1),dtBTorAs2hInt(nR,1),lo_map)
            !write(*,"(A,2I3,ES20.13)") "upB: after  dtBTor2hInt = ",nLMB,nR,dtBTor2hInt(nR,1)
         end if
      end do
      !PERFOFF

      !----- equations for inner core are different:
      !      D_lP1(lm1)=l+1, O_sr=sigma/sigma_ic
      !      NOTE: no hyperdiffusion in inner core !
      if ( l_cond_ic ) then
         !PERFON('upB_ic')
         do nR=2,n_r_ic_max-1
            do lm1=lmStart_00,lmStop
               l1=lm2l(lm1)
               m1=lm2m(lm1)
               dbdt_icLast(lm1,nR)=dbdt_icLast(lm1,nR) -                &
                    coex*opm*O_sr*dLh(st_map%lm2(l1,m1))*or2(n_r_max) * &
                    (                        ddb_ic(lm1,nR) +           &
                    2.d0*D_lP1(st_map%lm2(l1,m1))*O_r_ic(nR)*db_ic(lm1,nR) )
               djdt_icLast(lm1,nR)=djdt_icLast(lm1,nR) -                &
                    coex*opm*O_sr*dLh(st_map%lm2(l1,m1))*or2(n_r_max) * &
                    (                        ddj_ic(lm1,nR) +           &
                    2.d0*D_lP1(st_map%lm2(l1,m1))*O_r_ic(nR)*dj_ic(lm1,nR) )
            end do
         end do
         nR=n_r_ic_max
         do lm1=lmStart_00,lmStop
            l1=lm2l(lm1)
            m1=lm2m(lm1)
            dbdt_icLast(lm1,nR)=dbdt_icLast(lm1,nR) -                &
                 coex*opm*O_sr*dLh(st_map%lm2(l1,m1))*or2(n_r_max) * &
                 (1.d0+2.d0*D_lP1(st_map%lm2(l1,m1)))*ddb_ic(lm1,nR)
            djdt_icLast(lm1,nR)=djdt_icLast(lm1,nR) -                &
                 coex*opm*O_sr*dLh(st_map%lm2(l1,m1))*or2(n_r_max) * &
                 (1.d0+2.d0*D_lP1(st_map%lm2(l1,m1)))*ddj_ic(lm1,nR)
         end do
         !PERFOFF
      end if

   end subroutine updateB
!-----------------------------------------------------------------------------
#ifdef WITH_PRECOND_BJ
   subroutine get_bMat(dt,l,hdif,bMat,bPivot,bMat_fac,jMat,jPivot,jMat_fac)
#else
   subroutine get_bMat(dt,l,hdif,bMat,bPivot,jMat,jPivot)
#endif
      !  +-------------+----------------+------------------------------------+
      !  |                                                                   |
      !  |  Purpose of this subroutine is to contruct the time step matrices |
      !  |  bmat(i,j) and ajmat for the dynamo equations.                    |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+

      !-- Input variables:
      real(kind=8), intent(in) :: dt
      integer,      intent(in) :: l
      real(kind=8), intent(in) :: hdif
    
      !-- Output variables:
      real(kind=8), intent(out) :: bMat(n_r_totMag,n_r_totMag)
      integer,      intent(out) :: bPivot(n_r_totMag)
      real(kind=8), intent(out) :: jMat(n_r_totMag,n_r_totMag)
      integer,      intent(out) :: jPivot(n_r_totMag)
#ifdef WITH_PRECOND_BJ
      real(kind=8), intent(out) :: bMat_fac(n_r_totMag),jMat_fac(n_r_totMag)
#endif
 
      !-- local variables:
      integer :: nR,nCheb,nRall
      integer :: info
      real(kind=8) :: l_P_1
      real(kind=8) :: O_dt,dLh
      real(kind=8) :: rRatio
 
#undef MATRIX_CHECK
#ifdef MATRIX_CHECK
      integer :: i,j
      real(kind=8) :: rcond
      integer ::ipiv(n_r_tot),iwork(n_r_tot)
      real(kind=8) :: work(4*n_r_tot),anorm,linesum
      real(kind=8) :: temp_Mat(n_r_tot,n_r_tot)
      integer,save :: counter=0
      integer :: filehandle
      character(len=100) :: filename
#endif

      nRall=n_r_max
      if ( l_cond_ic ) nRall=nRall+n_r_ic_max
      O_dt=1.D0/dt
      dLh=dble(l*(l+1))
    
      !-- matricies depend on degree l but not on order m,
      !   we thus have to construct bmat and ajmat for each l:
    
      !-- do loop limits introduced to get rid of compiler warning !
    
      l_P_1=dble(l+1)
    
      do nR=2,n_r_max-1
         do nCheb=1,n_r_max
            bMat(nR,nCheb)=                       cheb_norm * ( &
                              O_dt*dLh*or2(nR)*cheb(nCheb,nR) - &
                      alpha*opm*lambda(nR)*hdif*dLh*or2(nR) * ( &
                                             d2cheb(nCheb,nR) - &
                                   dLh*or2(nR)*cheb(nCheb,nR) ) )
    
            jMat(nR,nCheb)=                       cheb_norm * ( &
                              O_dt*dLh*or2(nR)*cheb(nCheb,nR) - &
                      alpha*opm*lambda(nR)*hdif*dLh*or2(nR) * ( &
                                             d2cheb(nCheb,nR) + &
                                 dLlambda(nR)*dcheb(nCheb,nR) - &
                                   dLh*or2(nR)*cheb(nCheb,nR) ) )
         end do
      end do
    
      if  ( l_LCR ) then
         do nR=2,n_r_max-1
            if ( nR<=n_r_LCR ) then
               do nCheb=1,n_r_max
                  bMat(nR,nCheb)= cheb_norm*(     dcheb(nCheb,nR) + &
                                   dble(l)*or1(nR)*cheb(nCheb,nR) ) 
    
                  jMat(nR,nCheb)= cheb_norm*cheb(nCheb,nR)
               end do
            end if
         end do
      end if
    
      !----- boundary conditions for outer core field:
      do nCheb=1,n_cheb_max

         if ( ktopb == 1 .or. ktopb == 3 ) then
         !-------- at CMB (nR=1):
         !         the internal poloidal field should fit a potential
         !         field (matrix bmat) and the toroidal field has to
         !         vanish (matrix ajmat).
            bMat(1,nCheb)=            cheb_norm * ( &
                                   dcheb(nCheb,1) + &
                     dble(l)*or1(1)*cheb(nCheb,1) + &
                                  conductance_ma* ( &
                                  d2cheb(nCheb,1) - &
                         dLh*or2(1)*cheb(nCheb,1) ) )
    
            jMat(1,nCheb)=            cheb_norm * ( &
                                    cheb(nCheb,1) + &
                    conductance_ma*dcheb(nCheb,1) )
         else if ( ktopb == 2 ) then
            write(*,*) '! Boundary condition ktopb=2 not defined!'
            stop
         else if ( ktopb == 4 ) then

            !----- pseudo vacuum condition, field has only
            !      a radial component, horizontal components
            !      vanish when aj and db are zero:
            bMat(1,nCheb)=cheb_norm*dcheb(nCheb,1)
            jMat(1,nCheb)=cheb_norm*cheb(nCheb,1)
         end if

         !-------- at IC (nR=n_r_max):
         if ( kbotb == 1 ) then
            !----------- insulating IC, field has to fit a potential field:
            bMat(n_r_max,nCheb)=       cheb_norm * ( &
                              dcheb(nCheb,n_r_max) - &
                 l_P_1*or1(n_r_max)*cheb(nCheb,n_r_max) )
            jMat(n_r_max,nCheb)=       cheb_norm*cheb(nCheb,n_r_max)
         else if ( kbotb == 2 ) then
            !----------- perfect conducting IC
            bMat(n_r_max-1,nCheb)=cheb_norm*d2cheb(nCheb,n_r_max)
            jMat(n_r_max,nCheb)  =cheb_norm* dcheb(nCheb,n_r_max)
         else if ( kbotb == 3 ) then
            !---------- finite conducting IC, four boundary conditions:
            !           continuity of b,j, (d b)/(d r) and (d j)/(d r)/sigma.
            !           note: n_r=n_r_max and n_r=n_r_max+1 stand for IC radius
            !           here we set the outer core part of the equations.
            !           the conductivity ratio sigma_ratio is used as
            !           an additional dimensionless parameter.
            bMat(n_r_max,nCheb)=  cheb_norm*cheb(nCheb,n_r_max)
            bMat(n_r_max+1,nCheb)=cheb_norm*dcheb(nCheb,n_r_max)
            jMat(n_r_max,nCheb)=  cheb_norm*cheb(nCheb,n_r_max)
            jMat(n_r_max+1,nCheb)=cheb_norm*sigma_ratio*dcheb(nCheb,n_r_max)
         else if ( kbotb == 4 ) then
            !----- Pseudovacuum conduction at lower boundary:
            bMat(n_r_max,nCheb)=cheb_norm*dcheb(nCheb,n_r_max)
            jMat(n_r_max,nCheb)=cheb_norm*cheb(nCheb,n_r_max)
            end if

         !-------- Imposed fields: (overwrites above IC boundary cond.)
         if ( l == 1 .and. ( imagcon == -1 .or. imagcon == -2 ) ) then
            bMat(n_r_max,nCheb)=  cheb_norm * cheb(nCheb,n_r_max)
         else if ( l == 3 .and. imagcon == -10 ) then
            if ( l_LCR ) then
               write(*,*) 'Imposed field not compatible with weak conducting region!'
               stop
            end if
            jMat(1,nCheb)      =  cheb_norm * cheb(nCheb,1)
            jMat(n_r_max,nCheb)=  cheb_norm * cheb(nCheb,n_r_max)
         else if ( n_imp == 1 ) then
            !-- This is the Uli Christensen idea where the external field is
            !   not fixed but compensates the internal field so that the
            !   radial field component vanishes at r/r_cmb=rrMP
            if ( l_LCR ) then
               write(*,*) 'Imposed field not compatible with weak conducting region!'
               stop
            end if
            rRatio=rrMP**dble(2*l+1)
            bMat(1,nCheb)=            cheb_norm * ( &
                                   dcheb(nCheb,1) + &
                     dble(l)*or1(1)*cheb(nCheb,1) - &
                    dble(2*l+1)*or1(1)/(1-rRatio) + &
                                  conductance_ma* ( &
                                  d2cheb(nCheb,1) - &
                         dLh*or2(1)*cheb(nCheb,1) ) )
         end if
    
      end do ! loop over cheb modes !
    
      !----- fill up with zeros:
      do nCheb=n_cheb_max+1,n_r_max
         bMat(1,nCheb)=0.D0
         jMat(1,nCheb)=0.D0
         if ( l_LCR ) then
            do nR=2,n_r_LCR
               bMat(nR,nCheb)=0.D0
               jMat(nR,nCheb)=0.D0
            end do
         end if
         if ( kbotb == 1 ) then
            bMat(n_r_max,nCheb)  =0.D0
            jMat(n_r_max,nCheb)  =0.D0
         else if ( kbotb == 2 ) then
            bMat(n_r_max-1,nCheb)=0.D0
            jMat(n_r_max,nCheb)  =0.D0
         else if ( kbotb == 3 ) then
            bMat(n_r_max,nCheb)  =0.D0
            bMat(n_r_max+1,nCheb)=0.D0
            jMat(n_r_max,nCheb)  =0.D0
            jMat(n_r_max+1,nCheb)=0.D0
         else if ( kbotb == 4 ) then
            bMat(n_r_max,nCheb)  =0.D0
            jMat(n_r_max,nCheb)  =0.D0
         end if
      end do
    
      !----- normalization for highest and lowest Cheb mode:
      do nR=1,n_r_max
         bMat(nR,1)      =0.5D0*bMat(nR,1)
         bMat(nR,n_r_max)=0.5D0*bMat(nR,n_r_max)
         jMat(nR,1)      =0.5D0*jMat(nR,1)
         jMat(nR,n_r_max)=0.5D0*jMat(nR,n_r_max)
      end do
      if ( kbotb == 3 ) then
         bMat(n_r_max+1,1)=0.5D0*bMat(n_r_max+1,1)
         bMat(n_r_max+1,n_r_max)=0.5D0*bMat(n_r_max+1,n_r_max)
         jMat(n_r_max+1,1)=0.5D0*jMat(n_r_max+1,1)
         jMat(n_r_max+1,n_r_max)=0.5D0*jMat(n_r_max+1,n_r_max)
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
    
               bMat(n_r_max+nR,n_r_max+nCheb) =       &
                    cheb_norm_ic*dLh*or2(n_r_max) * ( &
                             O_dt*cheb_ic(nCheb,nR) - &
                                   alpha*opm*O_sr * ( &
                                d2cheb_ic(nCheb,nR) + &
                    2.d0*l_P_1*O_r_ic(nR)*dcheb_ic(nCheb,nR) )   )
    
               jMat(n_r_max+nR,n_r_max+nCheb)=bMat(n_r_max+nR,n_r_max+nCheb)
            end do
    
            !----- Special treatment for r=0, asymptotic of 1/r dr
            nR=n_r_ic_max
            bMat(n_r_max+nR,n_r_max+nCheb) =       &
                 cheb_norm_ic*dLh*or2(n_r_max) * ( &
                          O_dt*cheb_ic(nCheb,nR) - &
                                  alpha*opm*O_sr * &
                 (1.d0+2.d0*l_P_1)*d2cheb_ic(nCheb,nR) )
    
            jMat(n_r_max+nR,n_r_max+nCheb)=bMat(n_r_max+nR,n_r_max+nCheb)
         end do
    
         !-------- boundary condition at r_icb:
         do nCheb=1,n_cheb_ic_max
            bMat(n_r_max,n_r_max+nCheb)=-cheb_norm_ic*cheb_ic(nCheb,1)
            bMat(n_r_max+1,n_r_max+nCheb)=             &
                 -cheb_norm_ic * ( dcheb_ic(nCheb,1) + &
                 l_P_1*or1(n_r_max)*cheb_ic(nCheb,1) )
            jMat(n_r_max,n_r_max+nCheb)=bMat(n_r_max,n_r_max+nCheb)
            jMat(n_r_max+1,n_r_max+nCheb)=bMat(n_r_max+1,n_r_max+nCheb)
         end do ! cheb modes
    
         !-------- fill with zeros:
         do nCheb=n_cheb_ic_max+1,n_r_ic_max
            bMat(n_r_max,n_r_max+nCheb)  =0.D0
            bMat(n_r_max+1,n_r_max+nCheb)=0.D0
            jMat(n_r_max,n_r_max+nCheb)  =0.D0
            jMat(n_r_max+1,n_r_max+nCheb)=0.D0
         end do
    
         !-------- normalization for lowest Cheb mode:
         do nR=n_r_max,n_r_tot
            bMat(nR,n_r_max+1)=0.5D0*bMat(nR,n_r_max+1)
            jMat(nR,n_r_max+1)=0.5D0*jMat(nR,n_r_max+1)
            bMat(nR,n_r_tot)  =0.5D0*bMat(nR,n_r_tot)
            jMat(nR,n_r_tot)  =0.5D0*jMat(nR,n_r_tot)
         end do
    
         !-------- fill matricies up with zeros:
         do nCheb=n_r_max+1,n_r_tot
            do nR=1,n_r_max-1
               bMat(nR,nCheb)=0.D0
               jMat(nR,nCheb)=0.D0
            end do
         end do
         do nCheb=1,n_r_max
            do nR=n_r_max+2,n_r_tot
               bMat(nR,nCheb)=0.D0
               jMat(nR,nCheb)=0.D0
            end do
         end do
    
      end if  ! conducting inner core ?
 
#ifdef WITH_PRECOND_BJ
      ! compute the linesum of each line
      do nR=1,n_r_tot
         bMat_fac(nR)=1.0D0/maxval(abs(bMat(nR,:)))
         bMat(nR,:) = bMat(nR,:)*bMat_fac(nR)
      end do
      do nR=1,n_r_tot
         jMat_fac(nR)=1.0D0/maxval(abs(jMat(nR,:)))
         jMat(nR,:) = jMat(nR,:)*jMat_fac(nR)
      end do
#endif

#ifdef MATRIX_CHECK
      ! copy the bMat to a temporary variable for modification
      write(filename,"(A,I3.3,A,I3.3,A)") "bMat_",l,"_",counter,".dat"
      open(NEWUNIT=filehandle,file=trim(filename))
      counter= counter+1
      
      do i=1,n_r_tot
         do j=1,n_r_tot
            write(filehandle,"(2ES20.12,1X)",advance="no") bMat(i,j)
         end do
         write(filehandle,"(A)") ""
      end do
      close(filehandle)
      temp_Mat=bMat
      anorm = 0.0D0
      do i=1,n_r_tot
         linesum = 0.0D0
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
      open(NEWUNIT=filehandle,file=trim(filename))
      counter= counter+1
      
      do i=1,n_r_tot
         do j=1,n_r_tot
            write(filehandle,"(2ES20.12,1X)",advance="no") jMat(i,j)
         end do
         write(filehandle,"(A)") ""
      end do
      close(filehandle)
      temp_Mat=jMat
      anorm = 0.0D0
      do i=1,n_r_tot
         linesum = 0.0D0
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

      !----- LU decomposition:
#ifdef WITH_MKL_LU
      call getrf(bMat(1:nRall,1:nRall),bPivot(1:nRall),info)
#else
      call sgefa(bMat,n_r_tot,nRall,bPivot,info)
#endif

      if ( info /= 0 ) then
         write(*,*) 'Singular matrix bmat in get_bmat.'
         stop '32'
      end if

      !----- LU decomposition:
#ifdef WITH_MKL_LU
      call getrf(jMat(1:nRall,1:nRall),jPivot(1:nRall),info)
#else
      call sgefa(jMat,n_r_tot,nRall,jPivot,info)
#endif
      if ( info /= 0 ) then
         write(*,*) '! Singular matrix ajmat in get_bmat!'
         stop '33'
      end if

   end subroutine get_bMat
!-----------------------------------------------------------------------------
end module updateB_mod
