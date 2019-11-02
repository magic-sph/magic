#include "perflib_preproc.cpp"
module updateXi_mod

   use omp_lib
   use precision_mod
   use truncation, only: n_r_max, lm_max, l_max
   use radial_data, only: n_r_icb, n_r_cmb
   use radial_functions, only: orho1, or1, or2, beta, rscheme_oc, r
   use physical_parameters, only: osc, kbotxi, ktopxi
   use num_param, only: alpha, dct_counter, solve_counter
   use init_fields, only: topxi, botxi
   use blocking, only: st_map, lo_map, lo_sub_map, llm, ulm
   use horizontal_data, only: dLh, hdif_Xi
   use logic, only: l_update_xi, l_finite_diff, l_full_sphere
   use parallel_mod, only: rank, chunksize, n_procs, get_openmp_blocks
   use radial_der, only: get_ddr, get_dr
   use constants, only: zero, one, two
   use fields, only: work_LMloc
   use mem_alloc, only: bytes_allocated
   use useful, only: abortRun
   use dense_matrices
   use real_matrices
   use band_matrices


   implicit none

   private

   !-- Local variables
   complex(cp), allocatable :: rhs1(:,:,:)
   integer :: maxThreads
   class(type_realmat), pointer :: xiMat(:), xi0Mat
#ifdef WITH_PRECOND_S
   real(cp), allocatable :: xiMat_fac(:,:)
#endif
#ifdef WITH_PRECOND_S0
   real(cp), allocatable :: xi0Mat_fac(:)
#endif
   logical, public, allocatable :: lXimat(:)

   public :: initialize_updateXi, finalize_updateXi, updateXi

contains

   subroutine initialize_updateXi

      integer :: ll, n_bands
      integer, pointer :: nLMBs2(:)

      nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

      if ( l_finite_diff ) then
         allocate( type_bandmat :: xiMat(nLMBs2(1+rank)) )
         allocate( type_bandmat :: xi0Mat )

         if ( ktopxi == 1 .and. kbotxi == 1 .and. rscheme_oc%order == 2 &
          &   .and. rscheme_oc%order_boundary <= 2 ) then ! Fixed composition at both boundaries
            n_bands = rscheme_oc%order+1
         else
            n_bands = max(2*rscheme_oc%order_boundary+1,rscheme_oc%order+1)
         end if

         call xi0Mat%initialize(n_bands,n_r_max,l_pivot=.true.)
         do ll=1,nLMBs2(1+rank)
            call xiMat(ll)%initialize(n_bands,n_r_max,l_pivot=.true.)
         end do
      else
         allocate( type_densemat :: xiMat(nLMBs2(1+rank)) )
         allocate( type_densemat :: xi0Mat )

         call xi0Mat%initialize(n_r_max,n_r_max,l_pivot=.true.)
         do ll=1,nLMBs2(1+rank)
            call xiMat(ll)%initialize(n_r_max,n_r_max,l_pivot=.true.)
         end do
      end if

#ifdef WITH_PRECOND_S
      allocate(xiMat_fac(n_r_max,nLMBs2(1+rank)))
      bytes_allocated = bytes_allocated+n_r_max*nLMBs2(1+rank)*SIZEOF_DEF_REAL
#endif
#ifdef WITH_PRECOND_S0
      allocate(xi0Mat_fac(n_r_max))
      bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_REAL
#endif
      allocate( lXimat(0:l_max) )
      bytes_allocated = bytes_allocated+(l_max+1)*SIZEOF_LOGICAL

#ifdef WITHOMP
      maxThreads=omp_get_max_threads()
#else
      maxThreads=1
#endif

      allocate( rhs1(n_r_max,lo_sub_map%sizeLMB2max,0:maxThreads-1) )
      bytes_allocated = bytes_allocated + n_r_max*lo_sub_map%sizeLMB2max*&
      &                 maxThreads*SIZEOF_DEF_COMPLEX

   end subroutine initialize_updateXi
!------------------------------------------------------------------------------
   subroutine finalize_updateXI

      integer, pointer :: nLMBs2(:)
      integer :: ll

      nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

      do ll=1,nLMBs2(1+rank)
         call xiMat(ll)%finalize()
      end do
      call xi0Mat%finalize()

      deallocate( lXimat )
#ifdef WITH_PRECOND_S
      deallocate(xiMat_fac)
#endif
#ifdef WITH_PRECOND_S0
      deallocate(xi0Mat_fac)
#endif
      deallocate( rhs1 )

   end subroutine finalize_updateXI
!------------------------------------------------------------------------------
   subroutine updateXi(xi,dxi,dVXirLM,dxidt,dxidtLast,w1,coex,dt)
      !
      !  updates the entropy field s and its radial derivatives
      !  adds explicit part to time derivatives of s
      !

      !-- Input of variables:
      real(cp),    intent(in) :: w1        ! weight for time step !
      real(cp),    intent(in) :: coex      ! factor depending on alpha
      real(cp),    intent(in) :: dt        ! time step
      complex(cp), intent(inout) :: dVXirLM(llm:ulm,n_r_max)

      !-- Input/output of scalar fields:
      complex(cp), intent(inout) :: xi(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dxidt(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dxidtLast(llm:ulm,n_r_max)
      !-- Output: udpated s,ds,dxidtLast
      complex(cp), intent(out) :: dxi(llm:ulm,n_r_max)

      !-- Local variables:
      real(cp) :: w2            ! weight of second time step
      real(cp) :: O_dt
      integer :: l1,m1              ! degree and order
      integer :: lm1,lmB,lm         ! position of (l,m) in array
      integer :: nLMB2,nLMB
      integer :: nR                 ! counts radial grid points
      integer :: n_r_out             ! counts cheb modes
      real(cp) ::  rhs(n_r_max) ! real RHS for l=m=0

      integer, pointer :: nLMBs2(:),lm2l(:),lm2m(:)
      integer, pointer :: sizeLMB2(:,:),lm2(:,:)
      integer, pointer :: lm22lm(:,:,:),lm22l(:,:,:),lm22m(:,:,:)

      integer :: threadid, start_lm, stop_lm
      integer :: iChunk,nChunks,size_of_last_chunk,lmB0

      if ( .not. l_update_xi ) return

      nLMBs2(1:n_procs) => lo_sub_map%nLMBs2
      sizeLMB2(1:,1:) => lo_sub_map%sizeLMB2
      lm22lm(1:,1:,1:) => lo_sub_map%lm22lm
      lm22l(1:,1:,1:) => lo_sub_map%lm22l
      lm22m(1:,1:,1:) => lo_sub_map%lm22m
      lm2(0:,0:) => lo_map%lm2
      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      nLMB=1+rank
      w2  =one-w1
      O_dt=one/dt


      !PERFON('upS_fin')
      !$omp parallel default(shared) private(start_lm,stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm, stop_lm)
      !--- Finish calculation of dxidt:
      call get_dr( dVXirLM,work_LMloc,ulm-llm+1,start_lm-llm+1,       &
           &       stop_lm-llm+1,n_r_max, rscheme_oc, nocopy=.true. )
      !$omp barrier

      !$omp do private(nR,lm)
      do nR=1,n_r_max
         do lm=llm,ulm
            dxidt(lm,nR)=orho1(nR)*(dxidt(lm,nR)-or2(nR)*work_LMloc(lm,nR))
         end do
      end do
      !$OMP end do
      !PERFOFF

      !$omp single
      call solve_counter%start_count()
      !$omp end single
      ! one subblock is linked to one l value and needs therefore once the matrix
      !$OMP SINGLE
      do nLMB2=1,nLMBs2(nLMB)
         ! this inner loop is in principle over the m values which belong to the
         ! l value
         !$OMP TASK default(shared) &
         !$OMP firstprivate(nLMB2) &
         !$OMP private(lm,lm1,l1,m1,lmB,threadid) &
         !$OMP private(nChunks,size_of_last_chunk,iChunk)
         nChunks = (sizeLMB2(nLMB2,nLMB)+chunksize-1)/chunksize
         size_of_last_chunk = chunksize + (sizeLMB2(nLMB2,nLMB)-nChunks*chunksize)

         ! This task treats one l given by l1
         l1=lm22l(1,nLMB2,nLMB)
         !write(*,"(3(A,I3),A)") "Launching task for nLMB2=",nLMB2," (l=",l1,") and scheduling ",nChunks," subtasks."

         if ( l1 == 0 ) then
            if ( .not. lXimat(l1) ) then
#ifdef WITH_PRECOND_S0
               call get_xi0Mat(dt,xi0Mat,xi0Mat_fac)
#else
               call get_xi0Mat(dt,xi0Mat)
#endif
               lXimat(l1)=.true.
            end if
         else
            if ( .not. lXimat(l1) ) then
#ifdef WITH_PRECOND_S
               call get_xiMat(dt,l1,hdif_Xi(st_map%lm2(l1,0)), &
                    &         xiMat(nLMB2),xiMat_fac(:,nLMB2))
#else
               call get_xiMat(dt,l1,hdif_Xi(st_map%lm2(l1,0)),xiMat(nLMB2))
#endif
                lXimat(l1)=.true.
            end if
         end if

         do iChunk=1,nChunks
            !$OMP TASK default(shared) &
            !$OMP firstprivate(iChunk) &
            !$OMP private(lmB0,lmB,lm,lm1,m1,nR,n_r_out) &
            !$OMP private(threadid)
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

               if ( l1 == 0 ) then
                  rhs(1)      =real(topxi(0,0))
                  rhs(n_r_max)=real(botxi(0,0))
                  do nR=2,n_r_max-1
                     rhs(nR)=real(xi(lm1,nR))*O_dt+ &
                     &     w1*real(dxidt(lm1,nR)) + &
                     &     w2*real(dxidtLast(lm1,nR))
                  end do

#ifdef WITH_PRECOND_S0
                  rhs(:) = xi0Mat_fac(:)*rhs(:)
#endif

                  call xi0Mat%solve(rhs)

               else ! l1  /=  0
                  lmB=lmB+1

                  rhs1(1,lmB,threadid)      =topxi(l1,m1)
                  rhs1(n_r_max,lmB,threadid)=botxi(l1,m1)
                  do nR=2,n_r_max-1
                     rhs1(nR,lmB,threadid)=xi(lm1,nR)*O_dt + &
                     &                    w1*dxidt(lm1,nR) + &
                     &                    w2*dxidtLast(lm1,nR)
                  end do

#ifdef WITH_PRECOND_S
                  do nR=1,n_r_max
                     rhs1(nR,lmB,threadid)=xiMat_fac(nR,nLMB2)*rhs1(nR,lmB,threadid)
                  end do
#endif

               end if
            end do
            !PERFOFF

            !PERFON('upXi_sol')
            if ( lmB  >  lmB0 ) then
               call xiMat(nLMB2)%solve(rhs1(:,lmB0+1:lmB,threadid),lmB-lmB0)
            end if
            !PERFOFF

            lmB=lmB0
            !PERFON('upXi_af')
            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               lm1=lm22lm(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)
               if ( l1 == 0 ) then
                  do n_r_out=1,rscheme_oc%n_max
                     xi(lm1,n_r_out)=rhs(n_r_out)
                  end do
               else
                  lmB=lmB+1
                  if ( m1 > 0 ) then
                     do n_r_out=1,rscheme_oc%n_max
                        xi(lm1,n_r_out)=rhs1(n_r_out,lmB,threadid)
                     end do
                  else
                     do n_r_out=1,rscheme_oc%n_max
                        xi(lm1,n_r_out)= cmplx(real(rhs1(n_r_out,lmB,threadid)), &
                        &                      0.0_cp,kind=cp)
                     end do
                  end if
               end if
            end do
            !PERFOFF
            !$OMP END TASK
         end do
         !$OMP END TASK
      end do     ! loop over lm blocks
      !$OMP END SINGLE
      !$omp single
      call solve_counter%stop_count(l_increment=.false.)
      !$omp end single

      !-- set cheb modes > rscheme_oc%n_max to zero (dealiazing)
      !$omp single
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         do lm1=llm,ulm
            xi(lm1,n_r_out)=zero
         end do
      end do
      !$omp end single

      !PERFON('upXi_drv')
      !$omp single
      call dct_counter%start_count()
      !$omp end single
      call get_ddr(xi, dxi, work_LMloc, ulm-llm+1, start_lm-llm+1, &
           &       stop_lm-llm+1, n_r_max, rscheme_oc, l_dct_in=.false.)
      call rscheme_oc%costf1(xi,ulm-llm+1,start_lm-llm+1,stop_lm-llm+1)
      !$omp barrier


      !-- Calculate explicit time step part:
      !$omp do private(nR,lm1)
      do nR=n_r_cmb+1,n_r_icb-1
         do lm1=llm,ulm
            dxidtLast(lm1,nR)=dxidt(lm1,nR)                                      &
                 & - coex*osc*hdif_Xi(st_map%lm2(lm2l(lm1),lm2m(lm1))) *         &
                 &   ( work_LMloc(lm1,nR)                                        &
                 &     + ( beta(nR)+two*or1(nR) ) * dxi(lm1,nR)                  &
                 &     - dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))*or2(nR)*xi(lm1,nR) &
                 &   )
         end do
      end do
      !$omp end do
      !$omp single
      call dct_counter%stop_count(l_increment=.false.)
      !$omp end single

      !$omp end parallel

   end subroutine updateXi
!------------------------------------------------------------------------------
#ifdef WITH_PRECOND_S0
   subroutine get_xi0Mat(dt,xiMat,xiMat_fac)
#else
   subroutine get_xi0Mat(dt,xiMat)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matrix
      !  xiMat0
      !

      !-- Input variables
      real(cp), intent(in) :: dt

      !-- Output variables
      class(type_realmat), intent(inout) :: xiMat
#ifdef WITH_PRECOND_S0
      real(cp), intent(out) :: xiMat_fac(n_r_max)
#endif

      !-- Local variables:
      real(cp) :: dat(n_r_max,n_r_max)
      integer :: info, nR_out, nR
      real(cp) :: O_dt

      O_dt=one/dt

      !----- Boundary condition:
      if ( ktopxi == 1 ) then
         !--------- Constant chemical composition at CMB:
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%rMat(1,:)
      else
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%drMat(1,:)
      end if

      if ( l_full_sphere ) then
         dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
      else
         if ( kbotxi == 1 ) then
            !--------- Constant composition at ICB:
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
         else
            !--------- Constant flux at ICB:
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         end if
      end if

      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme_oc%n_max+1,n_r_max
            dat(1,nR_out)      =0.0_cp
            dat(n_r_max,nR_out)=0.0_cp
         end do
      end if

      !-- Fill bulk points
      do nR_out=1,n_r_max
         do nR=2,n_r_max-1
            dat(nR,nR_out)= rscheme_oc%rnorm * (                        &
            &                         O_dt*rscheme_oc%rMat(nR,nR_out) - &
            &             alpha*osc*(    rscheme_oc%d2rMat(nR,nR_out) + &
            &  (beta(nR)+two*or1(nR))*    rscheme_oc%drMat(nR,nR_out) ) )
         end do
      end do

      !----- Factors for highest and lowest cheb mode:
      do nR=1,n_r_max
         dat(nR,1)      =rscheme_oc%boundary_fac*dat(nR,1)
         dat(nR,n_r_max)=rscheme_oc%boundary_fac*dat(nR,n_r_max)
      end do

#ifdef WITH_PRECOND_S0
      ! compute the linesum of each line
      do nR=1,n_r_max
         xiMat_fac(nR)=one/maxval(abs(dat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,n_r_max
         dat(nR,:) = dat(nR,:)*xiMat_fac(nR)
      end do
#endif

      !-- Array copy
      call xiMat%set_data(dat)

      !---- LU decomposition:
      call xiMat%prepare(info)
      if ( info /= 0 ) call abortRun('! Singular matrix xiMat0!')

   end subroutine get_xi0Mat
!-----------------------------------------------------------------------------
#ifdef WITH_PRECOND_S
   subroutine get_xiMat(dt,l,hdif,xiMat,xiMat_fac)
#else
   subroutine get_xiMat(dt,l,hdif,xiMat)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matricies
      !  xiMat(i,j) for the equation for the chemical composition.
      !

      !-- Input variables
      real(cp), intent(in) :: dt
      real(cp), intent(in) :: hdif
      integer,  intent(in) :: l

      !-- Output variables
      class(type_realmat), intent(inout) :: xiMat
#ifdef WITH_PRECOND_S
      real(cp),intent(out) :: xiMat_fac(n_r_max)
#endif

      !-- Local variables:
      integer :: info, nR_out, nR
      real(cp) :: O_dt, dLh
      real(cp) :: dat(n_r_max,n_r_max)

      O_dt=one/dt

      dLh=real(l*(l+1),kind=cp)

      !----- Boundary coditions:
      if ( ktopxi == 1 ) then
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%rMat(1,:)
      else
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%drMat(1,:)
      end if

      if ( l_full_sphere ) then
         !dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
      else
         if ( kbotxi == 1 ) then
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
         else
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         end if
      end if

      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme_oc%n_max+1,n_r_max
            dat(1,nR_out)      =0.0_cp
            dat(n_r_max,nR_out)=0.0_cp
         end do
      end if

      !----- Bulk points
      do nR_out=1,n_r_max
         do nR=2,n_r_max-1
            dat(nR,nR_out)= rscheme_oc%rnorm * (                         &
            &                          O_dt*rscheme_oc%rMat(nR,nR_out) - &
            & alpha*osc*hdif*(            rscheme_oc%d2rMat(nR,nR_out) + &
            & ( beta(nR)+two*or1(nR) )*    rscheme_oc%drMat(nR,nR_out) - &
            &      dLh*or2(nR)*             rscheme_oc%rMat(nR,nR_out) ) )
         end do
      end do

      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         dat(nR,1)      =rscheme_oc%boundary_fac*dat(nR,1)
         dat(nR,n_r_max)=rscheme_oc%boundary_fac*dat(nR,n_r_max)
      end do

#ifdef WITH_PRECOND_S
      ! compute the linesum of each line
      do nR=1,n_r_max
         xiMat_fac(nR)=one/maxval(abs(dat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,n_r_max
         dat(nR,:) = dat(nR,:)*xiMat_fac(nR)
      end do
#endif

      !-- Array copy
      call xiMat%set_data(dat)

      !----- LU decomposition:
      call xiMat%prepare(info)
      if ( info /= 0 ) call abortRun('Singular matrix xiMat!')

   end subroutine get_xiMat
!-----------------------------------------------------------------------------
end module updateXi_mod
