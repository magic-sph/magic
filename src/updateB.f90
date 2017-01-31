#include "perflib_preproc.cpp"
module updateB_mod

#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif

   use omp_lib
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, n_r_tot, n_r_ic_max, &
       &                 n_cheb_ic_max, n_r_ic_maxMag, n_r_maxMag, &
       &                 n_r_totMag, lm_max, l_maxMag
   use radial_functions, only: chebt_ic,or2,r_cmb,chebt_ic_even, d2cheb_ic, &
       &                       cheb_norm_ic,dr_fac_ic,lambda,dLlambda,o_r_ic,r,&
       &                       or1, cheb_ic, dcheb_ic,rscheme_oc
   use radial_data, only: n_r_cmb,n_r_icb
   use physical_parameters, only: n_r_LCR,opm,O_sr,kbotb, imagcon, tmagcon, &
                                 sigma_ratio, conductance_ma, ktopb, kbotb
   use init_fields, only: bpeaktop, bpeakbot
   use num_param, only: alpha
   use blocking, only: nLMBs,st_map,lo_map,st_sub_map,lo_sub_map,lmStartB,lmStopB
   use horizontal_data, only: dLh, dPhi, hdif_B, D_l, D_lP1
   use logic, only: l_cond_ic, l_LCR, l_rot_ic, l_mag_nl, l_b_nl_icb, &
       &            l_b_nl_cmb, l_update_b, l_RMS
   use RMS, only: dtBPolLMr, dtBPol2hInt, dtBTor2hInt
   use constants, only: pi, zero, one, two, three, half
   use special
   use algebra, only: cgeslML, sgefa
   use LMLoop_data, only: llmMag,ulmMag
   use parallel_mod, only:  rank,chunksize
   use RMS_helpers, only: hInt2PolLM, hInt2TorLM
   use fields, only: work_LMloc
   use radial_der_even, only: get_ddr_even
   use radial_der, only: get_dr, get_ddr
   
   implicit none

   private

   !-- Local work arrays:
   complex(cp), allocatable :: workB(:,:)
   complex(cp), allocatable :: rhs1(:,:,:),rhs2(:,:,:)
   complex(cp), allocatable :: dtT(:), dtP(:)
   real(cp), allocatable :: bMat(:,:,:)
   real(cp), allocatable :: jMat(:,:,:)
   integer, allocatable :: bPivot(:,:)
   integer, allocatable :: jPivot(:,:)
#ifdef WITH_PRECOND_BJ
   real(cp), allocatable :: bMat_fac(:,:)
   real(cp), allocatable :: jMat_fac(:,:)
#endif
   logical, public, allocatable :: lBmat(:)


   integer :: maxThreads

   public :: initialize_updateB, finalize_updateB, updateB

contains

   subroutine initialize_updateB

      allocate( bMat(n_r_tot,n_r_tot,l_maxMag) )
      allocate( jMat(n_r_tot,n_r_totMag,l_maxMag) )
      bytes_allocated = bytes_allocated+(2*n_r_tot*n_r_tot*l_maxMag) &
      &                 *SIZEOF_DEF_REAL
      allocate( bPivot(n_r_tot,l_maxMag) )
      allocate( jPivot(n_r_tot,l_maxMag) )
      bytes_allocated = bytes_allocated+2*n_r_tot*l_maxMag*SIZEOF_INTEGER

#ifdef WITH_PRECOND_BJ
      allocate(bMat_fac(n_r_tot,l_maxMag))
      allocate(jMat_fac(n_r_tot,l_maxMag))
      bytes_allocated = bytes_allocated+2*n_r_tot*l_maxMag*SIZEOF_DEF_REAL
#endif
      allocate( lBmat(0:l_maxMag) )
      bytes_allocated = bytes_allocated+(l_maxMag+1)*SIZEOF_LOGICAL


      if ( l_RMS ) then
         allocate( workB(llmMag:ulmMag,n_r_max) )
         bytes_allocated = bytes_allocated+(ulmMag-llmMag+1)*n_r_max* & 
         &                 SIZEOF_DEF_COMPLEX
      end if

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
                      lo_sub_map%sizeLMB2max*SIZEOF_DEF_COMPLEX


   end subroutine initialize_updateB
!-----------------------------------------------------------------------------
   subroutine finalize_updateB

      deallocate( bMat, jMat, bPivot, jPivot, lBmat )

#ifdef WITH_PRECOND_BJ
      deallocate(bMat_fac,jMat_fac)
#endif
      deallocate( dtT, dtP, rhs1, rhs2 )
      if ( l_RMS ) deallocate( workB )

   end subroutine finalize_updateB
!-----------------------------------------------------------------------------
   subroutine updateB(b,db,ddb,aj,dj,ddj,dVxBhLM, &
       &             dbdt,dbdtLast,djdt,djdtLast, &
       &             b_ic,db_ic,ddb_ic,aj_ic,dj_ic,ddj_ic, &
       &             dbdt_icLast,djdt_icLast, &
       &             b_nl_cmb,aj_nl_cmb,aj_nl_icb,omega_ic, &
       &             w1,coex,dt,time,nLMB,lRmsNext)
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
      complex(cp), intent(in) :: b_nl_cmb(:)  ! nonlinear BC for b at CMB
      complex(cp), intent(in) :: aj_nl_cmb(:) ! nonlinear BC for aj at CMB
      complex(cp), intent(in) :: aj_nl_icb(:) ! nonlinear BC for dr aj at ICB
      complex(cp), intent(inout) :: dVxBhLM(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: dbdt(llmMag:ulmMag,n_r_maxMag)
      real(cp),    intent(in) :: omega_ic
      real(cp),    intent(in) :: w1    ! weight for time step !
      real(cp),    intent(in) :: coex  ! factor depending on alpha
      real(cp),    intent(in) :: dt
      real(cp),    intent(in) :: time
      integer,     intent(in) :: nLMB
      logical,     intent(in) :: lRmsNext

      !-- Input/output of scalar potentials and time stepping arrays:
      complex(cp), intent(inout) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(inout) :: aj(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(inout) :: dbdtLast(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(inout) :: djdt(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(inout) :: djdtLast(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(inout) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(inout) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(inout) :: dbdt_icLast(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(inout) :: djdt_icLast(llmMag:ulmMag,n_r_ic_maxMag)

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
      real(cp) :: w2             ! weight of second time step
      real(cp) :: O_dt
      real(cp) :: yl0_norm,prefac!External magnetic field of general l

      integer :: l1,m1               ! degree and order
      integer :: lm1,lm,lmB          ! position of (l,m) in array
      integer :: lmStart,lmStop      ! max and min number of orders m
      integer :: lmStart_00          ! excluding l=0,m=0
      integer :: nLMB2
      integer :: n_r_out             ! No of cheb polynome (degree+1)
      integer :: nR                  ! No of radial grid point
      integer :: n_r_real            ! total number of used grid points
      integer :: n_r_top,n_r_bot

      complex(cp) :: fac
      complex(cp) :: dbdt_ic,djdt_ic  ! they are calculated here !

      integer, pointer :: nLMBs2(:),lm2l(:),lm2m(:)
      integer, pointer :: sizeLMB2(:,:),lm2(:,:)
      integer, pointer :: lm22lm(:,:,:),lm22l(:,:,:),lm22m(:,:,:)

      !-- for feedback
      real(cp) :: ff,cimp,aimp,b10max

      real(cp), save :: direction

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

      ! output the input 
      !write(*,"(4(A,I4),I4,A,I4)") "nLMB=",nLMB,", from ",lmStart," to ",lmStop,&
      !     &", reals: ",lmStart_00,lmStop,", nLMBs2 = ",nLMBs2(nLMB)

      w2  =one-w1
      O_dt=one/dt

      !--- Start with finishing djdt:
      !    dVxBhLM is still in the R-distributed space,
      !    the ouput workA is in the LM-distributed space.
      !if (2*lmStart-1 - llmMag+1 /= 1) then
      !   write(*,"(I4,2(A,I6))") rank,": lmStart = ",lmStart,", llmMag = ",llmMag
      !   STOP
      !end if
      !if (lmStop  /=  ulmMag) then
      !   write(*,"(I4,A,2I6)") rank,": ",ulmMag,lmStop
      !   stop
      !end if
      !call get_dr( dVxBhLM,workA, &
      !     &         ulmMag-llmMag+1,(2*lmStart-1)-llmMag+1,lmStop-llmMag+1, &
      !     &         n_r_max,rscheme_oc )
      ! simplified interface
      !PRINT*,rank,": computing for ",ulmMag-llmMag+1," rows, chebt_oc = ",chebt_oc

      !PERFON('upB_fin')
      all_lms=lmStop-lmStart_00+1
#ifdef WITHOMP
      if (all_lms < maxThreads) then
         call omp_set_num_threads(all_lms)
      end if
#endif
      !$OMP PARALLEL &
      !$OMP private(iThread,start_lm,stop_lm) &
      !$OMP shared(all_lms,per_thread,lmStop,lmStart_00) &
      !$OMP shared(dVxBhLM,work_LMloc,djdt,or2,rscheme_oc) &
      !$OMP shared(n_r_max,nThreads,llmMag,ulmMag)
      !$OMP SINGLE
#ifdef WITHOMP
      nThreads=omp_get_num_threads()
#else 
      nThreads=1
#endif
      !-- Get radial derivatives of s: work_LMloc,dsdtLast used as work arrays
      per_thread=all_lms/nThreads
      !$OMP END SINGLE
      !$OMP BARRIER
      !$OMP DO
      do iThread=0,nThreads-1
         start_lm=lmStart_00+iThread*per_thread
         stop_lm = start_lm+per_thread-1
         if (iThread == nThreads-1) stop_lm=lmStop

         call get_dr( dVxBhLM,work_LMloc,ulmMag-llmMag+1,start_lm-llmMag+1, &
              &       stop_lm-llmMag+1,n_r_max,rscheme_oc,nocopy=.true. )

      end do
      !$OMP end do

      !$OMP do private(nR)
      do nR=1,n_r_max
         do lm=lmStart_00,lmStop
            djdt(lm,nR)=djdt(lm,nR)+or2(nR)*work_LMloc(lm,nR)
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
                    &        bMat(1,1,l1),bPivot(1,l1), bMat_fac(1,l1),&
                    &        jMat(1,1,l1),jPivot(1,l1), jMat_fac(1,l1))
#else
               call get_bMat(dt,l1,hdif_B(st_map%lm2(l1,0)), &
                    &        bMat(1,1,l1),bPivot(1,l1),      &
                    &        jMat(1,1,l1),jPivot(1,l1) )
#endif
               lBmat(l1)=.true.
            end if
         end if

         do iChunk=1,nChunks
            !$OMP TASK if (nChunks>1) default(shared) &
            !$OMP firstprivate(iChunk) &
            !$OMP private(lmB0,lmB,lm,lm1,m1,nR) &
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
                     rhs1(1,lmB,threadid) = 0.0_cp
                     rhs2(1,lmB,threadid) = 0.0_cp
                  end if

                  rhs1(n_r_max,lmB,threadid)=0.0_cp
                  if ( kbotb == 2 ) rhs1(n_r_max-1,lmB,threadid)=0.0_cp

                  rhs2(n_r_max,lmB,threadid)=0.0_cp
                  if ( m1 == 0 ) then   ! Magnetoconvection boundary conditions
                     if ( imagcon /= 0 .and. tmagcon <= time ) then
                        if ( l_LCR ) then
                           write(*,*) 'LCR not compatible with imposed field!'
                           stop
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
                           write(*,*) 'LCR not compatible with imposed field!'
                           stop
                        end if                          

                        !General normalization for spherical harmonics of degree l and order 0
                        yl0_norm=half*sqrt((2*l1+1)/pi)

                        !Prefactor for CMB matching condition
                        prefac = real(2*l1+1,kind=cp)/real(l1*(l1+1),kind=cp)     
                        

                        bpeaktop=prefac*fac_loop(l1)*amp_curr*8.0e-1_cp/yl0_norm

                        rhs1(1,lmB,threadid)=cmplx(bpeaktop,0.0_cp,kind=cp)

                     end if
 


                     if ( n_imp > 1 .and. l1 == l_imp ) then
                         if ( l_LCR ) then
                            write(*,*) 'LCR not compatible with imposed field!'
                            stop
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
                        rhs1(nR,lmB,threadid)= ( w1*dbdt(lm1,nR)+w2*dbdtLast(lm1,nR) ) &
                        &          + O_dt*dLh(st_map%lm2(l1,m1))*or2(nR)*b(lm1,nR)
                        rhs2(nR,lmB,threadid)= ( w1*djdt(lm1,nR)+w2*djdtLast(lm1,nR) ) &
                        &          + O_dt*dLh(st_map%lm2(l1,m1))*or2(nR)*aj(lm1,nR)
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
                        if ( omega_ic == 0.0_cp .or. .not. l_rot_ic .or. &
                        &    .not. l_mag_nl ) then
                           dbdt_ic=zero
                           djdt_ic=zero
                        else
                           fac=-omega_ic*or2(n_r_max)*dPhi(st_map%lm2(l1,m1))* &
                           &    dLh(st_map%lm2(l1,m1))
                           dbdt_ic=fac*b_ic(lm1,nR)
                           djdt_ic=fac*aj_ic(lm1,nR)
                        end if
                        rhs1(n_r_max+nR,lmB,threadid)=                 &
                        &      ( w1*dbdt_ic + w2*dbdt_icLast(lm1,nR) ) &
                        &      +O_dt*dLh(st_map%lm2(l1,m1))*or2(n_r_max) * b_ic(lm1,nR)
                        rhs2(n_r_max+nR,lmB,threadid)=                 &
                        &      ( w1*djdt_ic + w2*djdt_icLast(lm1,nR) ) &
                        &      +O_dt*dLh(st_map%lm2(l1,m1))*or2(n_r_max) * aj_ic(lm1,nR)

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
               call cgeslML(bMat(:,:,l1),n_r_tot,n_r_real, &
                    bPivot(:,l1),rhs1(:,lmB0+1:lmB,threadid),lmB-lmB0)
               call cgeslML(jMat(:,:,l1),n_r_tot,n_r_real, &
                    jPivot(:,l1),rhs2(:,lmB0+1:lmB,threadid),lmB-lmB0)
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
                           b_ic(lm1,n_r_out) = rhs1(n_r_max+n_r_out,lmB,threadid)
                           aj_ic(lm1,n_r_out)= rhs2(n_r_max+n_r_out,lmB,threadid)
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
      !$OMP END PARALLEL

      !-- Set cheb modes > rscheme_oc%n_max to zero (dealiazing)
      !   for inner core modes > 2*n_cheb_ic_max = 0
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         do lm1=lmStart_00,lmStop
            b(lm1,n_r_out) =zero
            aj(lm1,n_r_out)=zero
         end do
      end do
      if ( l_cond_ic ) then
         do n_r_out=n_cheb_ic_max+1,n_r_ic_max
            do lm1=lmStart_00,lmStop
               b_ic(lm1,n_r_out) =zero
               aj_ic(lm1,n_r_out)=zero
            end do
         end do
      end if

      !PERFON('upB_drv')
      all_lms=lmStop-lmStart_00+1
#ifdef WITHOMP
      if (all_lms < maxThreads) then
         call omp_set_num_threads(all_lms)
      end if
#endif
      !$OMP PARALLEL &
      !$OMP private(iThread,start_lm,stop_lm) &
      !$OMP shared(all_lms,per_thread,lmStart_00,lmStop) &
      !$OMP shared(b,db,ddb,aj,dj,ddj,dbdtLast,djdtLast) &
      !$OMP shared(rscheme_oc,n_r_max,nThreads,llmMag,ulmMag) &
      !$OMP shared(l_cond_ic,b_ic,db_ic,ddb_ic,aj_ic,dj_ic,ddj_ic) &
      !$OMP shared(chebt_ic,chebt_ic_even) &
      !$OMP shared(n_r_ic_max,n_cheb_ic_max,dr_fac_ic)
      !$OMP SINGLE
#ifdef WITHOMP
      nThreads=omp_get_num_threads()
#else
      nThreads=1
#endif
      !-- Get radial derivatives of s: work_LMloc,dsdtLast used as work arrays
      per_thread=all_lms/nThreads
      !write(*,"(2(A,I5))") "nThreads = ",nThreads,", per_thread = ",per_thread
      !$OMP END SINGLE
      !$OMP BARRIER
      !$OMP DO
      do iThread=0,nThreads-1
         start_lm=lmStart_00+iThread*per_thread
         stop_lm = start_lm+per_thread-1
         if (iThread == nThreads-1) stop_lm=lmStop

         !-- Radial derivatives: dbdtLast and djdtLast used as work arrays
         !PERFON('upB_cb')
         call rscheme_oc%costf1(b,ulmMag-llmMag+1,start_lm-llmMag+1, &
              &                 stop_lm-llmMag+1)
         !PERFOFF
         !PERFON('upB_db')
         call get_ddr(b,db,ddb,ulmMag-llmMag+1,start_lm-llmMag+1, &
              &       stop_lm-llmMag+1,n_r_max,rscheme_oc)
         !PERFOFF
         call rscheme_oc%costf1(aj, ulmMag-llmMag+1, start_lm-llmMag+1, &
              &                 stop_lm-llmMag+1)
         call get_ddr(aj,dj,ddj,ulmMag-llmMag+1,start_lm-llmMag+1, &
              &       stop_lm-llmMag+1,n_r_max,rscheme_oc)
         
         !-- Same for inner core:
         if ( l_cond_ic ) then
            call chebt_ic%costf1( b_ic, ulmMag-llmMag+1, start_lm-llmMag+1, &
                 &                stop_lm-llmMag+1, dbdtLast)
            call get_ddr_even( b_ic,db_ic,ddb_ic, ulmMag-llmMag+1, &
                 &             start_lm-llmMag+1,stop_lm-llmMag+1, &
                 &             n_r_ic_max,n_cheb_ic_max, dr_fac_ic,&
                 &             dbdtLast,djdtLast, chebt_ic, chebt_ic_even )
            call chebt_ic%costf1( aj_ic, ulmMag-llmMag+1, start_lm-llmMag+1, &
                 &               stop_lm-llmMag+1, dbdtLast)
            call get_ddr_even( aj_ic,dj_ic,ddj_ic, ulmMag-llmMag+1,  &
                 &             start_lm-llmMag+1, stop_lm-llmMag+1,  &
                 &             n_r_ic_max,n_cheb_ic_max, dr_fac_ic,  &
                 &             dbdtLast,djdtLast, chebt_ic, chebt_ic_even )
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
      end if

      if ( lRmsNext ) then
         n_r_top=n_r_cmb
         n_r_bot=n_r_icb
      else
         n_r_top=n_r_cmb+1
         n_r_bot=n_r_icb-1
      end if

      do nR=n_r_top,n_r_bot
         do lm1=lmStart_00,lmStop
            l1=lm2l(lm1)
            m1=lm2m(lm1)
            dbdtLast(lm1,nR)= dbdt(lm1,nR) -                    &
            &    coex*opm*lambda(nR)*hdif_B(st_map%lm2(l1,m1))* &
            &                  dLh(st_map%lm2(l1,m1))*or2(nR) * &
            &    ( ddb(lm1,nR) - dLh(st_map%lm2(l1,m1))*or2(nR)*b(lm1,nR) )
            djdtLast(lm1,nR)= djdt(lm1,nR) -                    &
            &    coex*opm*lambda(nR)*hdif_B(st_map%lm2(l1,m1))* &
            &                  dLh(st_map%lm2(l1,m1))*or2(nR) * &
            &    ( ddj(lm1,nR) + dLlambda(nR)*dj(lm1,nR) -      &
            &      dLh(st_map%lm2(l1,m1))*or2(nR)*aj(lm1,nR) )
            if ( lRmsNext ) then
               dtP(lm1)=O_dt*dLh(st_map%lm2(l1,m1))*or2(nR) &
               &             * (  b(lm1,nR)-work_LMloc(lm1,nR) )
               dtT(lm1)=O_dt*dLh(st_map%lm2(l1,m1))*or2(nR) &
               &             * ( aj(lm1,nR)-workB(lm1,nR) )
            end if
         end do
         if ( lRmsNext ) then
            call hInt2PolLM(dtP,llmMag,ulmMag,nR,lmStart_00,lmStop,&
                 &          dtBPolLMr(1,nR),dtBPol2hInt(1,nR,1),lo_map)
            call hInt2TorLM(dtT,llmMag,ulmMag,nR,lmStart_00,lmStop, &
                 &          dtBTor2hInt(1,nR,1),lo_map)
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
               &    coex*opm*O_sr*dLh(st_map%lm2(l1,m1))*or2(n_r_max) * &
               &    (                        ddb_ic(lm1,nR) +           &
               &    two*D_lP1(st_map%lm2(l1,m1))*O_r_ic(nR)*db_ic(lm1,nR) )
               djdt_icLast(lm1,nR)=djdt_icLast(lm1,nR) -                &
               &    coex*opm*O_sr*dLh(st_map%lm2(l1,m1))*or2(n_r_max) * &
               &    (                        ddj_ic(lm1,nR) +           &
               &    two*D_lP1(st_map%lm2(l1,m1))*O_r_ic(nR)*dj_ic(lm1,nR) )
            end do
         end do
         nR=n_r_ic_max
         do lm1=lmStart_00,lmStop
            l1=lm2l(lm1)
            m1=lm2m(lm1)
            dbdt_icLast(lm1,nR)=dbdt_icLast(lm1,nR) -                &
            &    coex*opm*O_sr*dLh(st_map%lm2(l1,m1))*or2(n_r_max) * &
            &    (one+two*D_lP1(st_map%lm2(l1,m1)))*ddb_ic(lm1,nR)
            djdt_icLast(lm1,nR)=djdt_icLast(lm1,nR) -                &
            &    coex*opm*O_sr*dLh(st_map%lm2(l1,m1))*or2(n_r_max) * &
            &    (one+two*D_lP1(st_map%lm2(l1,m1)))*ddj_ic(lm1,nR)
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
      !
      !  Purpose of this subroutine is to contruct the time step matrices 
      !  bmat(i,j) and ajmat for the dynamo equations.                    
      !

      !-- Input variables:
      real(cp), intent(in) :: dt
      integer,  intent(in) :: l
      real(cp), intent(in) :: hdif
    
      !-- Output variables:
      real(cp), intent(out) :: bMat(n_r_totMag,n_r_totMag)
      integer,  intent(out) :: bPivot(n_r_totMag)
      real(cp), intent(out) :: jMat(n_r_totMag,n_r_totMag)
      integer,  intent(out) :: jPivot(n_r_totMag)
#ifdef WITH_PRECOND_BJ
      real(cp), intent(out) :: bMat_fac(n_r_totMag),jMat_fac(n_r_totMag)
#endif
 
      !-- local variables:
      integer :: nR,nCheb,nR_out,nRall
      integer :: info
      real(cp) :: l_P_1
      real(cp) :: O_dt,dLh
      real(cp) :: rRatio
 
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
      O_dt=one/dt
      dLh=real(l*(l+1),kind=cp)
    
      !-- matricies depend on degree l but not on order m,
      !   we thus have to construct bmat and ajmat for each l:
    
      !-- do loop limits introduced to get rid of compiler warning !
    
      l_P_1=real(l+1,kind=cp)
    
      do nR=2,n_r_max-1
         do nR_out=1,n_r_max
            bMat(nR,nR_out)=                       rscheme_oc%rnorm * (     &
            &                 O_dt*dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out) - &
            &         alpha*opm*lambda(nR)*hdif*dLh*or2(nR) * (             &
            &                                rscheme_oc%d2rMat(nR,nR_out) - &
            &                      dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out) ) )
    
            jMat(nR,nR_out)=                       rscheme_oc%rnorm * (     &
            &                 O_dt*dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out) - &
            &         alpha*opm*lambda(nR)*hdif*dLh*or2(nR) * (             &
            &                                rscheme_oc%d2rMat(nR,nR_out) + &
            &                    dLlambda(nR)*rscheme_oc%drMat(nR,nR_out) - &
            &                      dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out) ) )
         end do
      end do
    
      if  ( l_LCR ) then
         do nR=2,n_r_max-1
            if ( nR<=n_r_LCR ) then
               do nR_out=1,n_r_max
                  bMat(nR,nR_out)= rscheme_oc%rnorm*(                       &
                  &                           rscheme_oc%drMat(nR,nR_out) + &
                  &    real(l,kind=cp)*or1(nR)*rscheme_oc%rMat(nR,nR_out) ) 
    
                  jMat(nR,nR_out)= rscheme_oc%rnorm*rscheme_oc%rMat(nR,nR_out)
               end do
            end if
         end do
      end if
    
      !----- boundary conditions for outer core field:
      do nR_out=1,rscheme_oc%n_max

         if ( ktopb == 1 .or. ktopb == 3 ) then
         !-------- at CMB (nR=1):
         !         the internal poloidal field should fit a potential
         !         field (matrix bmat) and the toroidal field has to
         !         vanish (matrix ajmat).
            bMat(1,nR_out)=            rscheme_oc%rnorm * (     &
            &                      rscheme_oc%drMat(1,nR_out) + &
            &     real(l,cp)*or1(1)*rscheme_oc%rMat(1,nR_out) + &
            &                     conductance_ma* (             &
            &                     rscheme_oc%d2rMat(1,nR_out) - &
            &            dLh*or2(1)*rscheme_oc%rMat(1,nR_out) ) )

            jMat(1,nR_out)=            rscheme_oc%rnorm * (     &
            &                       rscheme_oc%rMat(1,nR_out) + &
            &       conductance_ma*rscheme_oc%drMat(1,nR_out) )
         else if ( ktopb == 2 ) then
            write(*,*) '! Boundary condition ktopb=2 not defined!'
            stop
         else if ( ktopb == 4 ) then

            !----- pseudo vacuum condition, field has only
            !      a radial component, horizontal components
            !      vanish when aj and db are zero:
            bMat(1,nR_out)=rscheme_oc%rnorm*rscheme_oc%drMat(1,nR_out)
            jMat(1,nR_out)=rscheme_oc%rnorm* rscheme_oc%rMat(1,nR_out)
         end if

         !-------- at IC (nR=n_r_max):
         if ( kbotb == 1 ) then
            !----------- insulating IC, field has to fit a potential field:
            bMat(n_r_max,nR_out)=rscheme_oc%rnorm * (                 &
            &                      rscheme_oc%drMat(n_r_max,nR_out) - &
            &    l_P_1*or1(n_r_max)*rscheme_oc%rMat(n_r_max,nR_out) )
            jMat(n_r_max,nR_out)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,nR_out)
         else if ( kbotb == 2 ) then
            !----------- perfect conducting IC
            bMat(n_r_max-1,nR_out)=rscheme_oc%rnorm*rscheme_oc%d2rMat(n_r_max,nR_out)
            jMat(n_r_max,nR_out)  =rscheme_oc%rnorm* rscheme_oc%drMat(n_r_max,nR_out)
         else if ( kbotb == 3 ) then
            !---------- finite conducting IC, four boundary conditions:
            !           continuity of b,j, (d b)/(d r) and (d j)/(d r)/sigma.
            !           note: n_r=n_r_max and n_r=n_r_max+1 stand for IC radius
            !           here we set the outer core part of the equations.
            !           the conductivity ratio sigma_ratio is used as
            !           an additional dimensionless parameter.
            bMat(n_r_max,nR_out)  =rscheme_oc%rnorm* rscheme_oc%rMat(n_r_max,nR_out)
            bMat(n_r_max+1,nR_out)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,nR_out)
            jMat(n_r_max,nR_out)  =rscheme_oc%rnorm* rscheme_oc%rMat(n_r_max,nR_out)
            jMat(n_r_max+1,nR_out)=rscheme_oc%rnorm*sigma_ratio* &
            &                      rscheme_oc%drMat(n_r_max,nR_out)
         else if ( kbotb == 4 ) then
            !----- Pseudovacuum conduction at lower boundary:
            bMat(n_r_max,nR_out)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,nR_out)
            jMat(n_r_max,nR_out)= rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,nR_out)
            end if

         !-------- Imposed fields: (overwrites above IC boundary cond.)
         if ( l == 1 .and. ( imagcon == -1 .or. imagcon == -2 ) ) then
            bMat(n_r_max,nR_out)=  rscheme_oc%rnorm * rscheme_oc%rMat(n_r_max,nR_out)
         else if ( l == 3 .and. imagcon == -10 ) then
            if ( l_LCR ) then
               write(*,*) 'Imposed field not compatible with weak conducting region!'
               stop
            end if
            jMat(1,nR_out)      =  rscheme_oc%rnorm * rscheme_oc%rMat(1,nR_out)
            jMat(n_r_max,nR_out)=  rscheme_oc%rnorm * rscheme_oc%rMat(n_r_max,nR_out)
         else if ( n_imp == 1 ) then
            !-- This is the Uli Christensen idea where the external field is
            !   not fixed but compensates the internal field so that the
            !   radial field component vanishes at r/r_cmb=rrMP
            if ( l_LCR ) then
               write(*,*) 'Imposed field not compatible with weak conducting region!'
               stop
            end if
            rRatio=rrMP**real(2*l+1,kind=cp)
            bMat(1,nR_out)=            rscheme_oc%rnorm * (     &
            &                      rscheme_oc%drMat(1,nR_out) + &
            &     real(l,cp)*or1(1)*rscheme_oc%rMat(1,nR_out) - &
            &    real(2*l+1,cp)*or1(1)/(1-rRatio) +             &
            &                     conductance_ma* (             &
            &                     rscheme_oc%d2rMat(1,nR_out) - &
            &            dLh*or2(1)*rscheme_oc%rMat(1,nR_out) ) )
         end if
    
      end do ! loop over cheb modes !
    
      !----- fill up with zeros:
      do nR_out=rscheme_oc%n_max+1,n_r_max
         bMat(1,nR_out)=0.0_cp
         jMat(1,nR_out)=0.0_cp
         if ( l_LCR ) then
            do nR=2,n_r_LCR
               bMat(nR,nR_out)=0.0_cp
               jMat(nR,nR_out)=0.0_cp
            end do
         end if
         if ( kbotb == 1 ) then
            bMat(n_r_max,nR_out)  =0.0_cp
            jMat(n_r_max,nR_out)  =0.0_cp
         else if ( kbotb == 2 ) then
            bMat(n_r_max-1,nR_out)=0.0_cp
            jMat(n_r_max,nR_out)  =0.0_cp
         else if ( kbotb == 3 ) then
            bMat(n_r_max,nR_out)  =0.0_cp
            bMat(n_r_max+1,nR_out)=0.0_cp
            jMat(n_r_max,nR_out)  =0.0_cp
            jMat(n_r_max+1,nR_out)=0.0_cp
         else if ( kbotb == 4 ) then
            bMat(n_r_max,nR_out)  =0.0_cp
            jMat(n_r_max,nR_out)  =0.0_cp
         end if
      end do
    
      !----- normalization for highest and lowest Cheb mode:
      do nR=1,n_r_max
         bMat(nR,1)      =rscheme_oc%boundary_fac*bMat(nR,1)
         bMat(nR,n_r_max)=rscheme_oc%boundary_fac*bMat(nR,n_r_max)
         jMat(nR,1)      =rscheme_oc%boundary_fac*jMat(nR,1)
         jMat(nR,n_r_max)=rscheme_oc%boundary_fac*jMat(nR,n_r_max)
      end do
      if ( kbotb == 3 ) then
         bMat(n_r_max+1,1)=rscheme_oc%boundary_fac*bMat(n_r_max+1,1)
         bMat(n_r_max+1,n_r_max)=rscheme_oc%boundary_fac*bMat(n_r_max+1,n_r_max)
         jMat(n_r_max+1,1)=rscheme_oc%boundary_fac*jMat(n_r_max+1,1)
         jMat(n_r_max+1,n_r_max)=rscheme_oc%boundary_fac*jMat(n_r_max+1,n_r_max)
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
               &    cheb_norm_ic*dLh*or2(n_r_max) * ( &
               &             O_dt*cheb_ic(nCheb,nR) - &
               &                   alpha*opm*O_sr * ( &
               &                d2cheb_ic(nCheb,nR) + &
               &    two*l_P_1*O_r_ic(nR)*dcheb_ic(nCheb,nR) )   )
    
               jMat(n_r_max+nR,n_r_max+nCheb)=bMat(n_r_max+nR,n_r_max+nCheb)
            end do
    
            !----- Special treatment for r=0, asymptotic of 1/r dr
            nR=n_r_ic_max
            bMat(n_r_max+nR,n_r_max+nCheb) =       &
            &    cheb_norm_ic*dLh*or2(n_r_max) * ( &
            &             O_dt*cheb_ic(nCheb,nR) - &
            &                     alpha*opm*O_sr * &
            &    (one+two*l_P_1)*d2cheb_ic(nCheb,nR) )
    
            jMat(n_r_max+nR,n_r_max+nCheb)=bMat(n_r_max+nR,n_r_max+nCheb)
         end do
    
         !-------- boundary condition at r_icb:
         do nCheb=1,n_cheb_ic_max
            bMat(n_r_max,n_r_max+nCheb)=-cheb_norm_ic*cheb_ic(nCheb,1)
            bMat(n_r_max+1,n_r_max+nCheb)=             &
            &    -cheb_norm_ic * ( dcheb_ic(nCheb,1) + &
            &    l_P_1*or1(n_r_max)*cheb_ic(nCheb,1) )
            jMat(n_r_max,n_r_max+nCheb)=bMat(n_r_max,n_r_max+nCheb)
            jMat(n_r_max+1,n_r_max+nCheb)=bMat(n_r_max+1,n_r_max+nCheb)
         end do ! cheb modes
    
         !-------- fill with zeros:
         do nCheb=n_cheb_ic_max+1,n_r_ic_max
            bMat(n_r_max,n_r_max+nCheb)  =0.0_cp
            bMat(n_r_max+1,n_r_max+nCheb)=0.0_cp
            jMat(n_r_max,n_r_max+nCheb)  =0.0_cp
            jMat(n_r_max+1,n_r_max+nCheb)=0.0_cp
         end do
    
         !-------- normalization for lowest Cheb mode:
         do nR=n_r_max,n_r_tot
            bMat(nR,n_r_max+1)=half*bMat(nR,n_r_max+1)
            jMat(nR,n_r_max+1)=half*jMat(nR,n_r_max+1)
            bMat(nR,n_r_tot)  =half*bMat(nR,n_r_tot)
            jMat(nR,n_r_tot)  =half*jMat(nR,n_r_tot)
         end do
    
         !-------- fill matricies up with zeros:
         do nCheb=n_r_max+1,n_r_tot
            do nR=1,n_r_max-1
               bMat(nR,nCheb)=0.0_cp
               jMat(nR,nCheb)=0.0_cp
            end do
         end do
         do nCheb=1,n_r_max
            do nR=n_r_max+2,n_r_tot
               bMat(nR,nCheb)=0.0_cp
               jMat(nR,nCheb)=0.0_cp
            end do
         end do
    
      end if  ! conducting inner core ?
 
#ifdef WITH_PRECOND_BJ
      ! compute the linesum of each line
      do nR=1,nRall
         bMat_fac(nR)=one/maxval(abs(bMat(nR,1:nRall)))
         bMat(nR,:) = bMat(nR,:)*bMat_fac(nR)
      end do
      do nR=1,nRall
         jMat_fac(nR)=one/maxval(abs(jMat(nR,1:nRall)))
         jMat(nR,:) = jMat(nR,:)*jMat_fac(nR)
      end do
#endif

#ifdef MATRIX_CHECK
      ! copy the bMat to a temporary variable for modification
      write(filename,"(A,I3.3,A,I3.3,A)") "bMat_",l,"_",counter,".dat"
      open(newunit=filehandle,file=trim(filename))
      counter= counter+1
      
      do i=1,n_r_tot
         do j=1,n_r_tot
            write(filehandle,"(2ES20.12,1X)",advance="no") bMat(i,j)
         end do
         write(filehandle,"(A)") ""
      end do
      close(filehandle)
      temp_Mat=bMat
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
            write(filehandle,"(2ES20.12,1X)",advance="no") jMat(i,j)
         end do
         write(filehandle,"(A)") ""
      end do
      close(filehandle)
      temp_Mat=jMat
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

      !----- LU decomposition:
      call sgefa(bMat,n_r_tot,nRall,bPivot,info)

      if ( info /= 0 ) then
         write(*,*) 'Singular matrix bmat in get_bmat.'
         stop '32'
      end if

      !----- LU decomposition:
      call sgefa(jMat,n_r_tot,nRall,jPivot,info)
      if ( info /= 0 ) then
         write(*,*) '! Singular matrix ajmat in get_bmat!'
         stop '33'
      end if

   end subroutine get_bMat
!-----------------------------------------------------------------------------
end module updateB_mod
