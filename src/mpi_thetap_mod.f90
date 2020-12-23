#include "perflib_preproc.cpp"
!----------------------------------------------------------------------------------
   !
module mpi_thetap_mod
   ! 
   ! This module contains the implementation of theta-parallel transpositions
   !

   use constants, only: zero
   use precision_mod
   use parallel_mod
   use mem_alloc
   use truncation
   use blocking, only: lm_balance, lo_map, st_map, llm, ulm
   use mpi_transp, only: type_mpitransp
   use fft, only: fft_many, ifft_many
   use LMmapping

   implicit none

   private

   public :: transform_m2phi, transform_phi2m, transpose_m2th, transpose_th2m 
   public :: transform_new2old, transform_old2new, test_field !TODO: remove!

contains

!----------------------------------------------------------------------------------
   subroutine transpose_m2th(f_mloc, f_thloc, n_fields)

      !-- Input variables:
      integer,     intent(in) :: n_fields
      complex(cp), intent(in) :: f_mloc(n_theta_max,n_m_loc,*)

      !-- Output variable:
      complex(cp), intent(inout) :: f_thloc(nThetaStart:nThetaStop,n_phi_max/2+1,*)
      
      !-- Local variables:
      complex(cp) :: sendbuf(n_m_loc*n_theta_max*n_fields)
      complex(cp) :: recvbuf(n_theta_loc*n_m_max*n_fields)
      integer :: sendcount(0:n_ranks_theta-1),recvcount(0:n_ranks_theta-1)
      integer :: senddispl(0:n_ranks_theta-1),recvdispl(0:n_ranks_theta-1)
      integer :: irank, pos, n_m, n_f, m_idx, n_t, l_t, u_t
      
      PERFON('m2thS')
      do irank=0,n_ranks_theta-1
         recvcount(irank)=dist_m(irank,0) * n_theta_loc * n_fields
         sendcount(irank)=dist_theta(irank,0) * n_m_loc * n_fields
      end do

      senddispl(0)=0
      recvdispl(0)=0
      do irank=1,n_ranks_theta-1
         senddispl(irank)=senddispl(irank-1)+sendcount(irank-1)
         recvdispl(irank)=recvdispl(irank-1)+recvcount(irank-1)
      end do
      PERFOFF

      do irank=0,n_ranks_theta-1
         l_t=dist_theta(irank,1)
         u_t=dist_theta(irank,2)
         n_t=dist_theta(irank,0)
         pos = senddispl(irank)+1
         do n_f=1,n_fields
            do n_m=1, n_m_loc
               sendbuf(pos:pos+n_t-1)=f_mloc(l_t:u_t,n_m,n_f)
               pos = pos+n_t
            end do
         end do
      end do
      
      call mpi_barrier(comm_theta, ierr)
#ifdef WITH_MPI
      PERFON('m2thW')
      call MPI_Alltoallv(sendbuf, sendcount, senddispl, MPI_DEF_COMPLEX, &
           &             recvbuf, recvcount, recvdispl, MPI_DEF_COMPLEX, &
           &             comm_theta, ierr)
      PERFOFF
#endif
      
      PERFON('m2thS')
      do irank=0,n_ranks_theta-1
         pos = recvdispl(irank)+1
         do n_f=1,n_fields
            do n_m=1,dist_m(irank,0)
               m_idx = dist_m(irank,n_m)/minc+1
               f_thloc(:,m_idx,n_f)=recvbuf(pos:pos+n_theta_loc-1)
               pos = pos+n_theta_loc
            end do
         end do
      end do
      PERFOFF

      !f_thloc(n_m_max+1:,:,*)=zero

   end subroutine transpose_m2th
   
!----------------------------------------------------------------------------------
   subroutine transpose_th2m(f_thloc, f_mloc, n_fields)

      !-- Input variables:
      integer,     intent(in) :: n_fields
      complex(cp), intent(in) :: f_thloc(nThetaStart:nThetaStop,n_phi_max/2+1,*)

      !-- Output variable:
      complex(cp), intent(out) :: f_mloc(n_theta_max,n_m_loc,*)
      
      !-- Local variables:
      complex(cp) :: recvbuf(n_m_loc*n_theta_max*n_fields)
      complex(cp) :: sendbuf(n_theta_loc*n_m_max*n_fields)
      integer :: sendcount(0:n_ranks_theta-1),recvcount(0:n_ranks_theta-1)
      integer :: senddispl(0:n_ranks_theta-1),recvdispl(0:n_ranks_theta-1)
      integer :: irank, pos, n_m, n_f, m_idx, n_t, l_t, u_t

      PERFON('th2mS')
      do irank=0,n_ranks_theta-1
         sendcount(irank)=dist_m(irank,0) * n_theta_loc * n_fields
         recvcount(irank)=dist_theta(irank,0) * n_m_loc * n_fields
      end do

      senddispl(0)=0
      recvdispl(0)=0
      do irank=1,n_ranks_theta-1
         senddispl(irank)=senddispl(irank-1)+sendcount(irank-1)
         recvdispl(irank)=recvdispl(irank-1)+recvcount(irank-1)
      end do
      
      do irank=0,n_ranks_theta-1
         pos = senddispl(irank)+1
         do n_f=1,n_fields
            do n_m=1,dist_m(irank,0)
               m_idx = dist_m(irank,n_m)/minc+1
               sendbuf(pos:pos+n_theta_loc-1)=f_thloc(:,m_idx,n_f)
               pos = pos+n_theta_loc
            end do
         end do
      end do
      PERFOFF

      call mpi_barrier(comm_theta, ierr)
#ifdef WITH_MPI
      PERFON('th2mW')
      call MPI_Alltoallv(sendbuf, sendcount, senddispl, MPI_DEF_COMPLEX, &
           &             recvbuf, recvcount, recvdispl, MPI_DEF_COMPLEX, &
           &             comm_theta, ierr)
      PERFOFF
#endif
      
      PERFON('th2mS')
      do irank=0,n_ranks_theta-1
         pos = recvdispl(irank)+1
         l_t=dist_theta(irank,1)
         u_t=dist_theta(irank,2)
         n_t=dist_theta(irank,0)
         do n_f=1,n_fields
            do n_m=1,n_m_loc
               f_mloc(l_t:u_t,n_m,n_f)=recvbuf(pos:pos+n_t-1)
               pos = pos+n_t
            end do
         end do
      end do
      PERFOFF

   end subroutine transpose_th2m
!----------------------------------------------------------------------------------
   subroutine transform_m2phi(fL, f)
      !-- Transforms from (θ,m) space into (φ,θ) space including transpositions 
      !   and FFT. 
      !   
      !   Author: Rafael Lago (MPCDF) April 2020
      !   TODO: some functions might requires this transformation for multiple
      !      fields at once (e.g. vector transform). In that case, a non-blocking 
      !      transpose_theta_m would make more sense. This would make this 
      !      function obsolete.
      !   TODO: there is a lot of room for immprovement here (e.g. in-place, 
      !     use just one intermediate buffer, vector transform, etc)
      !
      
      !-- Input variables
      complex(cp), intent(inout) :: fL(n_theta_max,n_m_loc)
      
      !-- Output variables
      real(cp),    intent(out)   :: f(n_theta_loc,n_phi_max)
      
      !-- Local variables
      complex(cp) :: lF(n_theta_loc,n_m_max)
      complex(cp) :: Ff(n_theta_loc,n_phi_max/2+1)
   
      call transpose_m2th(fL, lF, 1)
      !-- TODO: The FFT must be performed for an array with the dimensions of 
      !   F_loc which may end up paded with zeroes.
      !   Is there any way to tell MKL to perform a "truncated" FFT?
      Ff = zero
      Ff(1:n_theta_loc,1:n_m_max) = lF
      
      call ifft_many(Ff, f)

   end subroutine transform_m2phi
!----------------------------------------------------------------------------------
   subroutine transform_phi2m(f, fL)
      !-- Transforms from (φ,θ) space into (θ,m) space including transpositions 
      !   and FFT. 
      !   
      !   Author: Rafael Lago (MPCDF) April 2020
      !   TODO: some functions might requires this transformation for multiple
      !      fields at once (e.g. vector transform). In that case, a non-blocking 
      !      transpose_theta_m would make more sense. This would make this 
      !      function obsolete.
      !   TODO: there is a lot of room for immprovement here (e.g. in-place, 
      !     use just one intermediate buffer, vector transform, etc)
      !
      
      !-- Input variables
      real(cp),    intent(inout) :: f(n_theta_loc,n_phi_max)
      
      !-- Output variables
      complex(cp), intent(out) :: fL(n_theta_max,n_m_loc)
      
      !-- Local variables
      complex(cp) :: lF(n_theta_loc,n_m_max)
      complex(cp) :: Ff(n_theta_loc,n_phi_max/2+1)
   
      call fft_many(f, Ff)
      lF(1:n_theta_loc,1:n_m_max) = Ff(1:n_theta_loc,1:n_m_max)
      call transpose_th2m(lF, fL, 1)

   end subroutine transform_phi2m
   
   
!----------------------------------------------------------------------------------
   
!  /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!  /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!
!  From here on, are only temporary functions. They are not optimized
!  (actually they are pretty terrible). They should be deleted once the transition 
!  is over
!  
!  LM Loop transposes and Gathers and Etc
!
!-------------------------------------------------------------------------------
   subroutine transform_new2old(Fmlo_new, Fmlo_old, n_r)
      integer,     intent(in) :: n_r ! Needed since it can be either n_r_ic_max, or n_r_max
      complex(cp), intent(in) :: Fmlo_new(n_mlo_loc, n_r)
      complex(cp), intent(inout) :: Fmlo_old(llm:ulm,   n_r)
      
      complex(cp) :: recvbuff(n_r)
      integer :: irank, ierr, lm, l, m, lo
      
      do lm=1,lm_max
         m = map_glbl_st%lm2m(lm)
         l = map_glbl_st%lm2l(lm)
         irank = map_mlo%ml2rnk(m,l)
         recvbuff = 0.0
         if (irank==rank) recvbuff = Fmlo_new(map_mlo%ml2i(m,l),:)
         call mpi_bcast(recvbuff, n_r, MPI_DEF_COMPLEX, irank, mpi_comm_world, ierr)
         lo = lo_map%lm2(l,m)
         if (lo>=llm .and. lo<=ulm) Fmlo_old(lo,:) = recvbuff
      end do
   end subroutine transform_new2old
!----------------------------------------------------------------------------------
   subroutine transform_old2new(Fmlo_old, Fmlo_new, n_r)

      integer,     intent(in) :: n_r ! Needed since it can be either n_r_ic_max, or n_r_max
      complex(cp), intent(in) :: Fmlo_old(llm:ulm,   n_r)
      complex(cp), intent(inout) :: Fmlo_new(n_mlo_loc, n_r)
      
      complex(cp) :: recvbuff(n_r)
      integer :: old2coord(l_max, l_max)
      integer :: irank, ierr, lm, l, m, size_irank, i
      
      do irank=0,n_ranks_r-1
         ! lolololololol
         if (irank==coord_r) then
            call mpi_bcast(ulm-llm+1, 1, MPI_INTEGER, irank, comm_r, ierr)
            do lm=llm,ulm
               m = lo_map%lm2m(lm)
               l = lo_map%lm2l(lm)
               
               call mpi_bcast(m, 1, MPI_INTEGER, irank, comm_r, ierr)
               call mpi_bcast(l, 1, MPI_INTEGER, irank, comm_r, ierr)
               
               call mpi_bcast(Fmlo_old(lm,:), n_r, MPI_DEF_COMPLEX, irank, comm_r, ierr)
               if (map_mlo%ml2rnk(m,l)==rank) Fmlo_new(map_mlo%ml2i(m,l),:) = Fmlo_old(lm,:)
            end do
         else
            call mpi_bcast(size_irank, 1, MPI_INTEGER, irank, comm_r, ierr)
            do i=1, size_irank
               call mpi_bcast(m, 1, MPI_INTEGER, irank, comm_r, ierr)
               call mpi_bcast(l, 1, MPI_INTEGER, irank, comm_r, ierr)
               call mpi_bcast(recvbuff, n_r, MPI_DEF_COMPLEX, irank, comm_r, ierr)
               if (map_mlo%ml2rnk(m,l)==rank) Fmlo_new(map_mlo%ml2i(m,l),:) = recvbuff
            end do
         end if
         
      end do
   end subroutine transform_old2new
!--------------------------------------------------------------------------------
!@> Delete me after conversion!
   subroutine test_field(newfield, oldfield, name, n_r)
      integer,     intent(in) :: n_r ! Needed since it can be either n_r_ic_max, or n_r_max
      character(len=*), intent(in) :: name
      complex(cp), intent(in) :: newfield(n_mlo_loc,n_r)
      complex(cp), intent(in) :: oldfield(llm:ulm,n_r)
      complex(cp) :: test_old(llm:ulm,n_r)
      complex(cp) :: test_new(n_mlo_loc,n_r)
      real(cp)    :: test_norm, error_threshold
      
      error_threshold = 1e-16 !EPSILON(1.0_cp)
      
      call transform_new2old(newfield, test_old, n_r)
      test_norm = ABS(SUM(oldfield - test_old))
      IF (test_norm>error_threshold) print *, "||",name,"n2o|| : ", test_norm

      call transform_old2new(oldfield, test_new, n_r)
      test_norm = ABS(SUM(newfield - test_new))
      IF (test_norm>error_threshold) print *, "||",name,"o2n|| : ", test_norm
      
   end subroutine
   
end module mpi_thetap_mod
