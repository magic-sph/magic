!$Id$
!***************************************************************************
SUBROUTINE mapDataICR(dataR,n_r_ic_old)
  !***************************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/17/02  by JW. -----------

  !--------------------------------------------------------------------
  !
  !  Copy (interpolate) data (read from disc file) from old grid structure 
  !  to new grid. Linear interploation is used in r if the radial grid
  !  structure differs
  !
  !  called in mapdata
  !
  !--------------------------------------------------------------------

  USE truncation
  USE radial_functions
  USE blocking

  IMPLICIT NONE

  !-- input: 
  INTEGER :: n_r_ic_old

  COMPLEX(kind=8) :: dataR(*)  ! old data 


  !-- local:
  INTEGER :: nR

  !----- For cheb transform:
  INTEGER :: n_r_ic_maxL
  !----- Note: n_local_fac allows for work memory that is
  !            n_local_fac times larger than the normal
  !            memory for one field (10*n_r_max)
  INTEGER,PARAMETER :: n_local_fac=10
  !PARAMETER (n_r_ic_maxL=n_local_fac*n_r_ic_max)
  INTEGER,ALLOCATABLE :: i_costf_ic_init_old(:)
  REAL(kind=8),ALLOCATABLE :: d_costf_ic_init_old(:)
  REAL(kind=8),ALLOCATABLE ::  work(:)
  REAL(kind=8) ::  cheb_norm_ic_old,scale


  !-- end of declaration
  !-------------------------------------------------------------------------

  ! allocation of the local arrays
  n_r_ic_maxL=n_local_fac*n_r_ic_max
  allocate( i_costf_ic_init_old(2*n_r_ic_maxL+2) )
  allocate( d_costf_ic_init_old(2*n_r_ic_maxL+5) )
  allocate( work(2*n_r_ic_maxL) )
  ! end of allocation of local arrays

  IF ( n_r_ic_maxL.LT.n_r_ic_old ) THEN
     WRITE(*,*) ' ! DIMENSION N_R_MAX_LOCAL TOO SMALL IN '
     WRITE(*,*) ' ! SUBROUTINE mapDataR !'
     WRITE(*,*) ' ! SHOULD BE AT LEAST n_r_max_old=',n_r_ic_old
     STOP 
  END IF
  IF ( n_r_ic_maxL.LT.n_r_ic_max ) THEN
     WRITE(*,*) ' ! DIMENSION N_R_MAX_LOCAL TOO SMALL IN '
     WRITE(*,*) ' ! SUBROUTINE mapDataR !'
     WRITE(*,*) ' ! SHOULD BE AT LEAST n_r_max_old=',n_r_ic_max
     STOP
  END IF

  !---- Initialize transform to cheb space:
  CALL init_costf1(n_r_ic_old,i_costf_ic_init_old,2*n_r_ic_maxL+2,&
       &                              d_costf_ic_init_old,2*n_r_ic_maxL+5)

  !----- Transform old data to cheb space:
  CALL costf1(dataR,2,1,2,work,                                   &
       &        i_costf_ic_init_old,d_costf_ic_init_old)

  !----- Fill up with zeros:
  IF ( n_r_ic_max.GT.n_r_ic_old ) THEN
     DO nR=n_r_ic_old,n_r_ic_max
        dataR(nR)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
     END DO
  END IF

  !----- Now transform to new radial grid points:
  CALL costf1(dataR,2,1,2,work,i_costf1_ic_init,d_costf1_ic_init)

  !----- Rescale:
  cheb_norm_ic_old=dsqrt(2.d0/dble(n_r_ic_old-1))
  scale=cheb_norm_ic_old/cheb_norm_ic
  DO nR=1,n_r_ic_max
     dataR(nR)=scale*dataR(nR)
  END DO

  ! deallocation of local arrays
  DEALLOCATE( i_costf_ic_init_old )
  DEALLOCATE( d_costf_ic_init_old )
  DEALLOCATE( work )
  ! end of deallocation


END SUBROUTINE mapDataICR

!---------------------------------------------------------------------
