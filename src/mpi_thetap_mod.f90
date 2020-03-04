! !----------------------------------------------------------------------------------
!    !
! module mpi_thetap_mod
!    ! 
!    ! This module contains the implementation of theta-parallel transpositions
!    !
! 
!    use precision_mod
!    use parallel_mod
!    use mem_alloc
!    use truncation, only: lm_max, n_r_max, l_max, minc, l_axi, &
!        &                 nRstart, nRstop, radial_balance,     &
!        &                 nR_per_rank
!    use blocking, only: lm_balance, lo_map, st_map, llm, ulm
!    use mpi_transp, only: type_mpitransp
! 
!    implicit none
! 
!    private
! 
! !    type, public, extends(type_mpitransp) :: type_thetap
! !       integer, allocatable :: rcounts(:)
! !       integer, allocatable :: scounts(:)
! !       integer, allocatable :: rdisp(:)
! !       integer, allocatable :: sdisp(:)
! !       integer :: max_send, max_recv
! !    contains
! !       procedure :: create_comm => create_comm_thetap
! !       procedure :: destroy_comm => destroy_comm_thetap
! !       procedure :: transp_lm2r => transp_lm2r_thetap
! !       procedure :: transp_r2lm => transp_r2lm_thetap
! ! !       procedure :: transp_lm2r => transp_lm2r_thetap_start
! ! !       procedure :: transp_r2lm => transp_r2lm_thetap_start
! ! !       procedure :: transp_lm2r_wait => transp_lm2r_thetap_wait
! ! !       procedure :: transp_r2lm_wait => transp_r2lm_thetap_wait
! !    end type type_thetap
! 
! 
! contains
!    
!    !-- Transposition from (m_loc,θ_glb) to (θ_loc,m_glb).
!    !   
!    !   
!    !   Author: Rafael Lago (MPCDF) August 2017
!    !
!    !-- TODO this with mpi_type to stride the data and check if performance 
!    !--      improves
!    !
!    subroutine transpose_m_theta(f_m_theta, f_theta_m)
!       complex(cp), intent(inout) :: f_m_theta(n_m_max, n_theta_loc)
!       complex(cp), intent(inout) :: f_theta_m(n_theta_max, n_m_loc)
!       
!       complex(cp) :: sendbuf(n_m_max * n_theta_loc)
!       complex(cp) :: recvbuf(n_m_loc, n_theta_max)
!       
!       integer :: sendcount(0:n_ranks_m-1)
!       integer :: recvcount(0:n_ranks_m-1)
!       integer :: senddispl(0:n_ranks_m-1)
!       integer :: recvdispl(0:n_ranks_m-1)
!       integer :: irank, j, itheta, m, pos
!       
!       PERFON('mt_tm')
!       pos = 1
!       do irank=0,n_ranks_m-1
!          !-- Copy each m which belongs to the irank-th rank into the send buffer
!          !   column-wise. That will simplify a lot things later
!          !
!          !-- TODO check performance of this; implementing this with mpi_type
!          !   striding the data will probably be faster
!          senddispl(irank) = pos-1
!          do itheta=1,n_theta_loc
!             do j=1,dist_m(irank,0)
!                m = dist_m(irank,j)/minc
!                sendbuf(pos) = f_m_theta(m+1,itheta)
!                pos = pos + 1
!             end do
!          end do
!          
!          sendcount(irank) = pos - senddispl(irank) - 1
!          recvdispl(irank) = irank*n_m_loc*dist_theta(irank,0)
!          recvcount(irank) =   n_m_loc*dist_theta(irank,0)
!       end do
!       
!       call MPI_ALLTOALLV(sendbuf, sendcount, senddispl, MPI_DOUBLE_COMPLEX, &
!                          recvbuf, recvcount, recvdispl, MPI_DOUBLE_COMPLEX, &
!                          comm_m, irank)
!       f_theta_m = transpose(recvbuf)
!       PERFOFF
!       
!    end subroutine transpose_m_theta
!    
!    !-- Transposition from (θ_loc,m_glb) to (m_loc,θ_glb)
!    !   
!    !   Author: Rafael Lago (MPCDF) August 2017
!    !
!    !-- TODO this with mpi_type to stride the data
!    !
!    subroutine transpose_theta_m(f_theta_m, f_m_theta)
!       complex(cp), intent(inout) :: f_theta_m(n_theta_max, n_m_loc)
!       complex(cp), intent(inout) :: f_m_theta(n_m_max, n_theta_loc)
!       
!       complex(cp) :: sendbuf(n_m_loc * n_theta_max)
!       complex(cp) :: recvbuf(n_theta_loc,  n_m_max)
!       
!       integer :: sendcount(0:n_ranks_theta-1)
!       integer :: recvcount(0:n_ranks_theta-1)
!       integer :: senddispl(0:n_ranks_theta-1)
!       integer :: recvdispl(0:n_ranks_theta-1)
!       integer :: irank, j, pos, n_t, l_t, u_t
!       integer :: m_arr(n_ranks_theta*n_m_array) 
!       
!       PERFON('tm_mt')
!       recvcount = 0
!       pos = 1
!       do irank=0,n_ranks_theta-1
!          !-- Copy each theta chunk so that the send buffer is contiguous
!          !-- TODO check performance of this; implementing this with mpi_type
!          !   striding the data will probably be faster
!          senddispl(irank) = pos-1
!          n_t = dist_theta(irank,0)
!          l_t = dist_theta(irank,1)
!          u_t = dist_theta(irank,2)
!          do j=1, n_m_loc
!             sendbuf(pos:pos + n_t - 1) = f_theta_m(l_t:u_t,j)
!             pos = pos + n_t
!          end do
!          
!          sendcount(irank) = pos - senddispl(irank) - 1
!          recvdispl(irank) = sum(recvcount)
!          recvcount(irank) = dist_m(irank,0) * n_t
!       end do
!       
!       call MPI_ALLTOALLV(sendbuf, sendcount, senddispl, MPI_DOUBLE_COMPLEX, &
!                          recvbuf, recvcount, recvdispl, MPI_DOUBLE_COMPLEX, &
!                          comm_theta, irank)
!       
!       !-- Now we reorder the receiver buffer. If the m distribution looks like:
!       !   rank 0: 0, 4,  8, 12, 16
!       !   rank 1: 1, 5,  9, 13
!       !   rank 2: 2, 6, 10, 14
!       !   rank 3: 3, 7, 11, 15
!       !   then the columns of recvbuf are ordered as 0,4,8,12,16,1,5,9,13(...)
!       !   and so forth. m_arr will contain this ordering (+1):
!       m_arr = reshape(transpose(dist_m(:,1:)), &
!                       (/n_ranks_m*n_m_array/))/minc + 1
!       j = 1
!       do pos = 1, n_ranks_theta*n_m_array
!          if (m_arr(pos) < 1) cycle
!          f_m_theta(m_arr(pos),:) = recvbuf(:,j)
!          j = j + 1
!       end do
!       PERFOFF
!    end subroutine transpose_theta_m
!    
! end module mpi_thetap_mod
