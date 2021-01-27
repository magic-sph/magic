#define LM_L_LOOP(lmStart, lmStop, action) do lm=lmStart,lmStop; l=lm2l(lm); action; end do;
module parallel_solvers
   !
   ! This module contains the routines that are used to solve linear banded problems
   ! with R-distributed arrays.
   !

   use precision_mod
   use parallel_mod
   use radial_data, only: n_r_cmb, n_r_icb
   use mem_alloc, only: bytes_allocated
   use constants, only: one
   use blocking, only: lm2l
   use truncation, only: lm_max

   implicit none

   private

   type, public :: type_tri_par   
      integer :: nRMin
      integer :: nRMax
      integer :: lMin 
      integer :: lMax
      real(cp), allocatable :: low(:,:)
      real(cp), allocatable :: diag(:,:)
      real(cp), allocatable :: up(:,:)
   contains
      procedure :: initialize => initialize_3
      procedure :: finalize => finalize_3
      procedure :: prepare_mat => prepare_mat_3
      procedure :: solver_up => solver_up_3
      procedure :: solver_dn => solver_dn_3
      procedure :: solver_single ! Used for one single right hand side
      procedure :: solver_finish => solver_finish_3
   end type

   type, public :: type_penta_par
      integer :: nRMin
      integer :: nRMax
      integer :: lMin
      integer :: lMax
      real(cp), allocatable :: low2(:,:)
      real(cp), allocatable :: low1(:,:)
      real(cp), allocatable :: diag(:,:)
      real(cp), allocatable :: up1(:,:)
      real(cp), allocatable :: up2(:,:)
   contains
      procedure :: initialize => initialize_5
      procedure :: finalize => finalize_5
      procedure :: prepare_mat => prepare_mat_5
      procedure :: solver_up => solver_up_5
      procedure :: solver_dn => solver_dn_5
      procedure :: solver_finish => solver_finish_5
   end type

contains

   subroutine initialize_3(this, nRstart, nRstop, lMin, lMax)
      !
      ! Memory allocation of a parallel tridiagonal matrix
      !

      class(type_tri_par) :: this
      integer, intent(in) :: nRstart
      integer, intent(in) :: nRstop
      integer, intent(in) :: lMin
      integer, intent(in) :: lMax

      this%nRMin=nRstart
      this%nRMax=nRstop
      this%lMin=lMin
      this%lMax=lMax
      allocate( this%low(lMin:lMax,nRstart:nRstop) )
      allocate( this%diag(lMin:lMax,nRstart:nRstop) )
      allocate( this%up(lMin:lMax,nRstart:nRstop) )
      bytes_allocated = bytes_allocated+3*(lMax-lMin+1)*(nRstop-nRstart+1)*SIZEOF_DEF_REAL

      !-- Fill an identity matrix by default
      this%low(:,:) =0.0_cp
      this%diag(:,:)=one
      this%up(:,:)  =0.0_cp

   end subroutine initialize_3
!-------------------------------------------------------------------------------------
   subroutine initialize_5(this, nRstart, nRstop, lMin, lMax)
      !
      ! Memory allocation of a parallel tridiagonal matrix
      !

      class(type_penta_par) :: this
      integer, intent(in) :: nRstart
      integer, intent(in) :: nRstop
      integer, intent(in) :: lMin
      integer, intent(in) :: lMax

      this%nRMin=nRstart
      this%nRMax=nRstop
      this%lMin=lMin
      this%lMax=lMax
      allocate( this%low1(lMin:lMax,nRstart:nRstop), this%low2(lMin:lMax,nRstart:nRstop) )
      allocate( this%diag(lMin:lMax,nRstart:nRstop) )
      allocate( this%up1(lMin:lMax,nRstart:nRstop), this%up2(lMin:lMax,nRstart:nRstop) )
      bytes_allocated = bytes_allocated+5*(lMax-lMin+1)*(nRstop-nRstart+1)*SIZEOF_DEF_REAL

      !-- Fill an identity matrix by default
      this%low2(:,:)=0.0_cp
      this%low1(:,:)=0.0_cp
      this%diag(:,:)=one
      this%up1(:,:) =0.0_cp
      this%up2(:,:) =0.0_cp

   end subroutine initialize_5
!-------------------------------------------------------------------------------------
   subroutine finalize_3(this)
      !
      ! Memory deallocation of a parallel tridiagonal solver
      !
      class(type_tri_par) :: this

      deallocate(this%low, this%diag, this%up)

   end subroutine finalize_3
!-------------------------------------------------------------------------------------
   subroutine finalize_5(this)
      !
      ! Memory deallocation of a parallel pentadiagonal solver
      !
      class(type_penta_par) :: this

      deallocate(this%low1, this%low2, this%diag, this%up1, this%up2)

   end subroutine finalize_5
!-------------------------------------------------------------------------------------
   subroutine prepare_mat_3(this)
      !
      ! LU factorisation of a tridiagonal matrix: the diagonal is overwritten
      !
      class(type_tri_par) :: this

      !-- Local variables
      integer :: l, nR
      real(cp) :: p

      !-- Set 'out-of-bound' values to zero for safety
      do l=this%lMin, this%lMax
         if ( this%nRMin == n_r_cmb ) this%low(l,this%nRMin)=0.0_cp
         if ( this%nRMax == n_r_icb ) this%up(l,this%nRMax) =0.0_cp
      end do

      do nR=this%nRMin,this%nRMax
         do l=this%lMin, this%lMax
            if ( nR == 1 ) p=this%diag(l,nR)
            if ( nR > 1 ) p=this%diag(l,nR)-this%low(l,nR)*this%up(l,nR-1)*this%diag(l,nR-1)
            this%diag(l,nR)=one/p
         end do
      end do

   end subroutine prepare_mat_3
!-------------------------------------------------------------------------------------
   subroutine prepare_mat_5(this)
      !
      ! LU factorisation of a pentadiagonal matrix
      !
      class(type_penta_par) :: this

      !-- Local variables
      integer :: nR, l, start_l, stop_l

      !$omp parallel default(shared) private(start_l,stop_l,l,nR)
      start_l=this%lMin; stop_l=this%lMax
      call get_openmp_blocks(start_l,stop_l)
      !$omp barrier

      !-- Set 'out-of-bound' values to zero for safety
      do l=start_l,stop_l
         if ( this%nRMin == n_r_cmb ) then
            this%low1(l,this%nRMin)  =0.0_cp
            this%low2(l,this%nRMin)  =0.0_cp
            this%low2(l,this%nRMin+1)=0.0_cp
         end if
         if ( this%nRMax == n_r_icb ) then
            this%up1(l,this%nRMax)  =0.0_cp
            this%up2(l,this%nRMax)  =0.0_cp
            this%up2(l,this%nRMax-1)=0.0_cp
         end if
      end do

      !-- Now proper LU factorisation
      nR=2
      do l=start_l,stop_l
         this%up1(l,nR)=this%up1(l,nR)-this%low1(l,nR)*this%up2(l,nR-1)/this%diag(l,nR-1)
         this%diag(l,nR)=this%diag(l,nR)-this%low1(l,nR)*this%up1(l,nR-1)/ &
         &               this%diag(l,nR-1)
      end do

      do nR=3,this%nRMax
         do l=start_l,stop_l
            this%low1(l,nR)=this%low1(l,nR)-this%low2(l,nR)*this%up1(l,nR-2)/ &
            &               this%diag(l,nR-2)
            this%up1(l,nR)=this%up1(l,nR)-this%low1(l,nR)*this%up2(l,nR-1)/ &
            &              this%diag(l,nR-1)
            this%diag(l,nR)=this%diag(l,nR)-this%low1(l,nR)*this%up1(l,nR-1)/ &
            &               this%diag(l,nR-1)-this%low2(l,nR)*this%up2(l,nR-2)/ &
            &               this%diag(l,nR-2)
          end do
      enddo

      do nR=1,this%nRMax
         do l=start_l,stop_l
            this%diag(l,nR)=one/this%diag(l,nR)
            this%up1(l,nR) =this%up1(l,nR)*this%diag(l,nR)
            this%up2(l,nR) =this%up2(l,nR)*this%diag(l,nR)
            this%low1(l,nR)=this%low1(l,nR)*this%diag(l,nR)
            this%low2(l,nR)=this%low2(l,nR)*this%diag(l,nR)
          end do
      enddo
      !$omp end parallel

   end subroutine prepare_mat_5
!-------------------------------------------------------------------------------------
   subroutine solver_single(this, x, nRstart, nRstop)
      !
      ! This routine is used to solve a solve one single linear system that does not
      ! depend on lm. This is for instance used for z10 when l_z10mat is required.
      !
      class(type_tri_par) :: this

      !-- Input variables
      integer, intent(in) :: nRstart, nRstop

      !-- Output variables
      complex(cp), intent(inout) :: x(nRstart-1:nRstop+1)

      !-- Local variables
      integer :: nR0,nR
      integer :: tag

      tag = 53976

      nR0 = nRstart
      if ( nRstart > n_r_cmb ) then ! Not the first block
         nR0 = nR0-1
#ifdef WITH_MPI
         call MPI_Recv(x(nR0), 1, MPI_DEF_COMPLEX, rank-1, tag, MPI_COMM_WORLD, &
              &        MPI_STATUS_IGNORE, ierr)
#endif
      else ! Lower boundary: x -> x - low * x(i-1)
         x(nR0)=x(nR0)-this%low(1,nR0)*x(nR0-1)
      end if

      do nR=nR0+1,nRstop
         x(nR)=x(nR)-this%diag(1,nR-1)*this%low(1,nR)*x(nR-1)
      end do

      if ( nRstop < n_r_icb ) then ! Not the last block
#ifdef WITH_MPI
         !call MPI_ISend(x(nRstop), 1, MPI_DEF_COMPLEX, rank+1, &
         !     &         tag, MPI_COMM_WORLD, array_req(req), ierr)
         !req = req+1
         call MPI_SSend(x(nRstop), 1, MPI_DEF_COMPLEX, rank+1, tag, MPI_COMM_WORLD, ierr)
#endif
      end if

      tag = tag+1
#ifdef WITH_MPI
      if ( nRstop < n_r_icb ) then ! This is not the last chunk
         call MPI_Recv(x(nRstop+1), 1, MPI_DEF_COMPLEX, rank+1, &
              &        tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
#endif

      do nR=nRstop,nRstart,-1
         x(nR)=(x(nR)-this%up(1,nR)*x(nR+1))*this%diag(1,nR)
      end do

#ifdef WITH_MPI
      if ( nRstart > n_r_cmb ) then
         !call MPI_ISend(x(nRstart), 1, MPI_DEF_COMPLEX, rank-1, &
         !     &         tag, MPI_COMM_WORLD, array_req(req), ierr)
         !req=req+1
         call MPI_SSend(x(nRstart), 1, MPI_DEF_COMPLEX, rank-1, tag, MPI_COMM_WORLD, ierr)
      end if
#endif

      tag = tag+1
#ifdef WITH_MPI
      if ( nRstart /= n_r_cmb ) then
         if ( nRstart == n_r_cmb+1 ) then ! send down
            call MPI_Ssend(x(nRstart), 1, MPI_DEF_COMPLEX, rank-1, tag, &
                 &         MPI_COMM_WORLD, ierr)
         end if

         if ( nRstop == n_r_cmb ) then ! send down
            call MPI_Recv(x(nRstop+1), 1, MPI_DEF_COMPLEX, rank+1, tag, &
                 &        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            !req = req+1
         end if
      end if

      tag = tag+1
      if ( nRstart > n_r_cmb ) then ! This is not the first block
         call MPI_Recv(x(nRstart-1), 1, MPI_DEF_COMPLEX, rank-1, &
              &         tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
         !req = req+1
      end if

      if ( nRstop < n_r_icb .and. nRstop >= n_r_cmb ) then ! This is not the last block
         call MPI_Ssend(x(nRstop), 1, MPI_DEF_COMPLEX, rank+1, tag, MPI_COMM_WORLD, ierr)
         !call MPI_Isend(x(nRstop), lmu-lmb+1, MPI_DEF_COMPLEX, rank+1, &
         !     &         tag, MPI_COMM_WORLD, array_req(req), ierr)
         !req = req+1
      end if
#endif

   end subroutine solver_single
!-------------------------------------------------------------------------------------
   subroutine solver_up_3(this, x, lmStart, lmStop, nRstart, nRstop, tag, array_req, &
              &           req, lms_block, nlm_block)
      !
      ! First part of the parallel tridiag solver: forward substitution
      !

      class(type_tri_par) :: this

      !-- Input variables
      integer, intent(in) :: lmStart ! Starting lm (OMP thread dependent)
      integer, intent(in) :: lmStop  ! Stopping lm (OMP thread dependent)
      integer, intent(in) :: nRstart ! Starting nR
      integer, intent(in) :: nRstop  ! Stopping nR
      integer, intent(in) :: lms_block ! Starting block-index of lm
      integer, intent(in) :: nlm_block ! Size of the block
      integer, intent(in) :: tag

      !-- Output variables
      integer,     intent(inout) :: array_req(:)
      complex(cp), intent(inout) :: x(1:lm_max, nRstart-1:nRstop+1)
      integer,     intent(inout) :: req

      !-- Local variables
      integer :: nR, nR0, lm, l, lmb, lmu

      lmb = lms_block
      lmu = lmb+nlm_block-1

      nR0 = nRstart
      if ( nRstart > n_r_cmb ) then ! Not the first block
         nR0 = nR0-1
#ifdef WITH_MPI
         !$omp master
         call MPI_Recv(x(lmb:lmu,nR0), lmu-lmb+1, MPI_DEF_COMPLEX, rank-1, tag, &
              &        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
         !$omp end master
         !$omp barrier
#endif
      else ! Lower boundary: x -> x - low * x(i-1)
         LM_L_LOOP(lmStart, lmStop, x(lm,nR0)=x(lm,nR0)-this%low(l,nR0)*x(lm,nR0-1))
         !$omp barrier
      end if

      do nR=nR0+1,nRstop
         LM_L_LOOP(lmStart, lmStop, x(lm,nR)=x(lm,nR)-this%diag(l,nR-1)*this%low(l,nR)*x(lm,nR-1))
      end do
      !$omp barrier

      if ( nRstop < n_r_icb ) then ! Not the last block
#ifdef WITH_MPI
         !$omp barrier
         !$omp master
         call MPI_ISend(x(lmb:lmu,nRstop), lmu-lmb+1, MPI_DEF_COMPLEX, rank+1, &
              &         tag, MPI_COMM_WORLD, array_req(req), ierr)
         req = req+1
         !call MPI_SSend(x(lmb:lmu,nRstop), lmu-lmb+1, MPI_DEF_COMPLEX, rank+1, &
         !     &         tag, MPI_COMM_WORLD, ierr)
         !$omp end master
#endif
      end if

   end subroutine solver_up_3
!-------------------------------------------------------------------------------------
   subroutine solver_up_5(this, x, lmStart, lmStop, nRstart, nRstop, tag, array_req, &
              &           req, lms_block, nlm_block)
      !
      ! First part of the parallel pentadiag solver: forward substitution
      !
      class(type_penta_par) :: this

      !-- Input variables
      integer, intent(in) :: lmStart ! Starting lm (OMP thread dependent)
      integer, intent(in) :: lmStop  ! Stopping lm (OMP thread dependent)
      integer, intent(in) :: nRstart ! Starting nR
      integer, intent(in) :: nRstop  ! Stopping nR
      integer, intent(in) :: lms_block ! Starting block-index of lm
      integer, intent(in) :: nlm_block ! Size of the block
      integer, intent(in) :: tag

      !-- Output variables
      integer, intent(inout) :: array_req(:)
      complex(cp), intent(inout) :: x(1:lm_max, nRstart-2:nRstop+2)
      integer, intent(inout) :: req

      !-- Local variables
      integer :: nR, nR0, lm, l, recv(2), lb, lu, lmb, lmu

      lmb = lms_block
      lmu = lmb+nlm_block-1
      lb=lmStart
      lu=lmStop

#ifdef WITH_MPI
      if ( nRstart > n_r_cmb ) then ! This is not the first block
         !$omp master
         call MPI_IRecv(x(lmb:lmu,nRstart-2), lmu-lmb+1, MPI_DEF_COMPLEX, rank-1,&
              &        tag, MPI_COMM_WORLD, recv(1), ierr)
         !-- Non blocking receive
         call MPI_IRecv(x(lmb:lmu,nRstart-1), lmu-lmb+1, MPI_DEF_COMPLEX, rank-1,&
              &        tag+1, MPI_COMM_WORLD, recv(2), ierr)
         !-- Non-blocking receive
         call MPI_Waitall(2, recv, MPI_STATUSES_IGNORE, ierr) ! wait to receive
         !$omp end master
         !$omp barrier
      end if
#endif

      do nR=nRstart,nRstop
         LM_L_LOOP(lb,lu,x(lm,nR)=this%diag(l,nR)*x(lm,nR)-this%low1(l,nR)*x(lm,nR-1)-this%low2(l,nR)*x(lm,nR-2))
      end do
      !$omp barrier

#ifdef WITH_MPI
      if ( nRstop < n_r_icb ) then ! This is not the last block
         !$omp barrier
         !$omp master
         call MPI_ISend(x(lmb:lmu,nRstop-1), lmu-lmb+1, MPI_DEF_COMPLEX, rank+1, &
              &         tag, MPI_COMM_WORLD, array_req(req), ierr)
         call MPI_ISend(x(lmb:lmu,nRstop), lmu-lmb+1, MPI_DEF_COMPLEX, rank+1, &
              &         tag+1, MPI_COMM_WORLD, array_req(req+1), ierr)
         req = req+2
         !$omp end master
      end if
#endif

   end subroutine solver_up_5
!-------------------------------------------------------------------------------------
   subroutine solver_dn_3(this, x, lmStart, lmStop, nRstart, nRstop, tag, array_req, &
              &           req, lms_block, nlm_block)
      !
      ! Second part of the parallel tridiag solver: backward substitution
      !

      class(type_tri_par) :: this

      !-- Input variables
      integer, intent(in) :: lmStart ! Starting lm (OMP thread dependent)
      integer, intent(in) :: lmStop  ! Stopping lm (OMP thread dependent)
      integer, intent(in) :: nRstart ! Starting nR
      integer, intent(in) :: nRstop  ! Stopping nR
      integer, intent(in) :: lms_block ! Starting block-index of lm
      integer, intent(in) :: nlm_block ! Size of the block
      integer, intent(in) :: tag

      !-- Output variables
      integer,     intent(inout) :: array_req(:)
      complex(cp), intent(inout) :: x(1:lm_max, nRstart-1:nRstop+1)
      integer,     intent(inout) :: req

      !-- Local variables
      integer :: nR, l, lm, lmb, lmu

      lmb = lms_block
      lmu = lmb+nlm_block-1

#ifdef WITH_MPI
      if ( nRstop < n_r_icb ) then ! This is not the last chunk
         !$omp master
         call MPI_Recv(x(lmb:lmu,nRstop+1), lmu-lmb+1, MPI_DEF_COMPLEX, rank+1, &
              &        tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
         !$omp end master
         !$omp barrier
      end if
#endif

      do nR=nRstop,nRstart,-1
         LM_L_LOOP(lmStart,lmStop,x(lm,nR)=(x(lm,nR)-this%up(l,nR)*x(lm,nR+1))*this%diag(l,nR))
      end do

#ifdef WITH_MPI
      if ( nRstart > n_r_cmb ) then
         !$omp barrier
         !$omp master
         call MPI_ISend(x(lmb:lmu,nRstart), lmu-lmb+1, MPI_DEF_COMPLEX, rank-1, &
              &         tag, MPI_COMM_WORLD, array_req(req), ierr)
         req=req+1
         !call MPI_SSend(x(lmb:lmu,nRstart), lmu-lmb+1, MPI_DEF_COMPLEX, rank-1, &
         !     &         tag, MPI_COMM_WORLD, ierr)
         !$omp end master
      end if
#endif

   end subroutine solver_dn_3
!-------------------------------------------------------------------------------------
   subroutine solver_dn_5(this, x, lmStart, lmStop, nRstart, nRstop, tag, array_req, &
              &           req, lms_block, nlm_block)
      !
      ! Second part of the parallel pentadiagonal solver: backward substitution
      !

      class(type_penta_par) :: this

      !-- Input variables
      integer, intent(in) :: lmStart ! Starting lm (OMP thread dependent)
      integer, intent(in) :: lmStop  ! Stopping lm (OMP thread dependent)
      integer, intent(in) :: nRstart ! Starting nR
      integer, intent(in) :: nRstop  ! Stopping nR
      integer, intent(in) :: lms_block ! Starting block-index of lm
      integer, intent(in) :: nlm_block ! Size of the block
      integer, intent(in) :: tag

      !-- Output variables
      integer,     intent(inout) :: array_req(:)
      complex(cp), intent(inout) :: x(lm_max, nRstart-2:nRstop+2)
      integer,     intent(inout) :: req

      !-- Local variables
      integer :: nR, l, lm, rcv(2), lmb, lmu

      lmb = lms_block
      lmu = lmb+nlm_block-1

#ifdef WITH_MPI
      if ( nRstop < n_r_icb ) then ! This is not the last rank
         !$omp master
         call MPI_IRecv(x(lmb:lmu,nRstop+1), lmu-lmb+1, MPI_DEF_COMPLEX, rank+1, &
              &        tag, MPI_COMM_WORLD, rcv(1), ierr)
         call MPI_IRecv(x(lmb:lmu,nRstop+2), lmu-lmb+1, MPI_DEF_COMPLEX, rank+1, &
              &        tag+1, MPI_COMM_WORLD, rcv(2), ierr)
         call MPI_Waitall(2, rcv, MPI_STATUSES_IGNORE, ierr)
         !$omp end master
         !$omp barrier
      end if
#endif

      do nR=nRstop,nRstart,-1
         LM_L_LOOP(lmStart,lmStop,x(lm,nR)=x(lm,nR)-this%up1(l,nR)*x(lm,nR+1)-this%up2(l,nR)*x(lm,nR+2))
      end do

#ifdef WITH_MPI
      !$omp barrier
      !$omp master
      if ( nRstart > n_r_cmb ) then
         call MPI_Isend(x(lmb:lmu,nRstart), lmu-lmb+1, MPI_DEF_COMPLEX, rank-1, &
              &         tag, MPI_COMM_WORLD, array_req(req), ierr)
         call MPI_Isend(x(lmb:lmu,nRstart+1), lmu-lmb+1, MPI_DEF_COMPLEX, rank-1, &
              &         tag+1, MPI_COMM_WORLD, array_req(req+1), ierr)
         req=req+2
      end if
      !$omp end master
#endif

   end subroutine solver_dn_5
!-------------------------------------------------------------------------------------
   subroutine solver_finish_3(this, x, lms_block, nlm_block, nRstart, nRstop, tag, &
              &               array_req, req)
      !
      ! Last part of the parallel tridiag solver: halo synchronisation
      !

      class(type_tri_par) :: this

      !-- Input variables
      integer, intent(in) :: nRstart  ! Starting index in radius
      integer, intent(in) :: nRstop   ! Stopping index in radius
      integer, intent(in) :: lms_block ! Starting block-index of lm
      integer, intent(in) :: nlm_block ! Size of the block
      integer, intent(in) :: tag

      !-- Output variables
      integer,     intent(inout) :: array_req(:)
      complex(cp), intent(inout) :: x(1:lm_max, nRstart-1:nRstop+1)
      integer,     intent(inout) :: req
      
      !-- Local variables
      integer :: lmb, lmu

      lmb = lms_block
      lmu = lmb+nlm_block-1

#ifdef WITH_MPI
      !$omp master
      if ( nRstart /= n_r_cmb ) then
         if ( nRstart == n_r_cmb+1 ) then ! send down
            call MPI_Isend(x(lmb:lmu,nRstart), lmu-lmb+1, MPI_DEF_COMPLEX, &
                 &         rank-1, tag, MPI_COMM_WORLD, array_req(req), ierr)
            req = req+1
         end if

         if ( nRstop == n_r_cmb ) then ! send down
            call MPI_Irecv(x(lmb:lmu,nRstop+1), lmu-lmb+1, MPI_DEF_COMPLEX, &
                 &         rank+1, tag, MPI_COMM_WORLD, array_req(req), ierr)
            req = req+1
         end if
      end if

      if ( nRstart > n_r_cmb ) then ! This is not the first block
         call MPI_Irecv(x(lmb:lmu,nRstart-1), lmu-lmb+1, MPI_DEF_COMPLEX, rank-1, &
              &         tag, MPI_COMM_WORLD, array_req(req), ierr)
         req = req+1
      end if

      if ( nRstop < n_r_icb .and. nRstop >= n_r_cmb ) then ! This is not the last block
         !call MPI_Ssend(x(lmb:lmu,nRstop), lmu-lmb+1, MPI_DEF_COMPLEX, rank+1, &
         !     &         tag, MPI_COMM_WORLD, ierr)
         call MPI_Isend(x(lmb:lmu,nRstop), lmu-lmb+1, MPI_DEF_COMPLEX, rank+1, &
              &         tag, MPI_COMM_WORLD, array_req(req), ierr)
         req = req+1
      end if
      !$omp end master
#endif

   end subroutine solver_finish_3
!-------------------------------------------------------------------------------------
   subroutine solver_finish_5(this, x, lms_block, nlm_block, nRstart, nRstop, tag, &
              &               array_req, req)
      !
      ! Last part of the parallel pentadiag solver: halo synchronisation
      !

      class(type_penta_par) :: this

      !-- Input variables
      integer, intent(in) :: nRstart  ! Starting index in radius
      integer, intent(in) :: nRstop   ! Stopping index in radius
      integer, intent(in) :: lms_block ! Starting block-index of lm
      integer, intent(in) :: nlm_block ! Size of the block
      integer, intent(in) :: tag

      !-- Output variables
      integer,     intent(inout) :: array_req(:) ! MPI requests
      complex(cp), intent(inout) :: x(lm_max, nRstart-2:nRstop+2) ! Solution
      integer,     intent(inout) :: req

      !-- Local variables
      integer :: istart, istop, lmb, lmu
      logical :: l_send_dn, l_send_up

      lmb = lms_block
      lmu = lmb+nlm_block-1
      l_send_dn = .false.
      l_send_up = .false.
      istart = n_r_cmb
      istop = n_r_icb

#ifdef WITH_MPI
      !$omp master
      if ( nRstart == n_r_cmb+1 ) then ! Send down
         call MPI_Isend(x(lmb:lmu,nRstart), lmu-lmb+1, MPI_DEF_COMPLEX, rank-1, &
              &         tag, MPI_COMM_WORLD, array_req(req), ierr)
         call MPI_Isend(x(lmb:lmu,nRstart+1), lmu-lmb+1, MPI_DEF_COMPLEX, rank-1, &
              &         tag+1, MPI_COMM_WORLD, array_req(req+1), ierr)
         req=req+2
      end if
      if ( nRStop==n_r_cmb ) then ! Receive up
         call MPI_Irecv(x(lmb:lmu,nRstop+1), lmu-lmb+1, MPI_DEF_COMPLEX, rank+1, &
              &         tag, MPI_COMM_WORLD, array_req(req), ierr)
         call MPI_Irecv(x(lmb:lmu,nRstop+2), lmu-lmb+1, MPI_DEF_COMPLEX, rank+1, &
              &         tag+1, MPI_COMM_WORLD, array_req(req+1), ierr)
         req=req+2
      end if

      if ( nRstart > n_r_cmb ) then ! This is not the first block
         l_send_dn=.true.
         istart=nRstart
      end if
      if ( (nRstop < n_r_icb) .and. (nRstop>=istart) ) then
         l_send_up=.true.
         istop=nRstop
      end if

      if ( l_send_up ) then
         !-- Update halo
         call MPI_Isend(x(lmb:lmu,istop), lmu-lmb+1, MPI_DEF_COMPLEX, rank+1, &
              &         tag, MPI_COMM_WORLD, array_req(req), ierr)
         req=req+1
      end if
      if ( l_send_dn ) then
         call MPI_Irecv(x(lmb:lmu,istart-2), lmu-lmb+1, MPI_DEF_COMPLEX, rank-1, &
              &         tag+1, MPI_COMM_WORLD, array_req(req), ierr)
         req=req+1
         if ( (istop > istart) .or. (.not. l_send_up) ) then
            call MPI_Irecv(x(lmb:lmu,istart-1), lmu-lmb+1, MPI_DEF_COMPLEX, &
                 &         rank-1, tag, MPI_COMM_WORLD, array_req(req), ierr)
            req=req+1
         else
            call MPI_Recv(x(lmb:lmu,istart-1), lmu-lmb+1, MPI_DEF_COMPLEX, &
                 &        rank-1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
         end if
      end if
      if ( l_send_up ) then
         call MPI_Isend(x(lmb:lmu,istop-1), lmu-lmb+1, MPI_DEF_COMPLEX, rank+1, &
              &         tag+1, MPI_COMM_WORLD, array_req(req), ierr)
         req=req+1
      end if
      !$omp end master
#endif

   end subroutine solver_finish_5
!-------------------------------------------------------------------------------------
end module parallel_solvers
