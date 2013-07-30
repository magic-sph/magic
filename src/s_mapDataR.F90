!$Id$
!***************************************************************************
SUBROUTINE mapDataR(dataR,n_r_max_old,lBc)
  !***************************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/17/02  by JW. -----------

  !---------------------------------------------------------------------------
  !
  !  Copy (interpolate) data (read from disc file) from old grid structure 
  !  to new grid. Linear interploation is used in r if the radial grid
  !  structure differs
  !
  !  called in mapdata
  !
  !---------------------------------------------------------------------------

  USE truncation
  USE radial_functions
  USE blocking
  USE horizontal_data

  IMPLICIT NONE

  INTEGER :: n_r_max_old
  LOGICAL :: lBc
  COMPLEX(kind=8) :: dataR(*)  ! old data 

  !-- local:
  INTEGER :: nR

  !----- For cheb transform:
  INTEGER :: n_r_maxL
  INTEGER,PARAMETER :: n_local_fac=10
  !integer :: n_r_maxL=n_local_fac*n_r_max
  !----- Note: n_local_fac allows for a radial grid point number 
  !            that is five times larger than n_r_max.
  INTEGER,ALLOCATABLE :: i_costf_init_old(:)
  REAL(kind=8),ALLOCATABLE :: d_costf_init_old(:)
  REAL(kind=8),ALLOCATABLE ::  work(:)
  REAL(kind=8) :: cheb_norm_old,scale

  !-- end of declaration
  !-------------------------------------------------------------------------
  n_r_maxL=n_local_fac*n_r_max
  allocate( i_costf_init_old(2*n_r_maxL+2) )
  allocate( d_costf_init_old(2*n_r_maxL+5) )
  allocate( work(2*n_r_maxL) )

  IF ( n_r_maxL.LT.n_r_max_old ) THEN
     WRITE(*,*) ' ! DIMENSION N_R_MAXL TOO SMALL IN '
     WRITE(*,*) ' ! SUBROUTINE mapDataR !'
     WRITE(*,*) ' ! SHOULD BE AT LEAST n_r_max_old=',n_r_max_old
     STOP
  END IF
  IF ( n_r_maxL.LT.n_r_max ) THEN
     WRITE(*,*) ' ! DIMENSION N_R_MAXL TOO SMALL IN '
     WRITE(*,*) ' ! SUBROUTINE mapDataR !'
     WRITE(*,*) ' ! SHOULD BE AT LEAST n_r_max_old=',n_r_max
     STOP
  END IF

  !----- Initialize transform to cheb space:
  CALL init_costf1(n_r_max_old,i_costf_init_old,2*n_r_maxL+2,     &
       &                               d_costf_init_old,2*n_r_maxL+5)

  !-- Guess the boundary values, since they have not been stored:
  IF ( lBc ) THEN
     dataR(1)=2.D0*dataR(2)-dataR(3)
     dataR(n_r_max_old)=2.D0*dataR(n_r_max_old-1) -               &
          &                             dataR(n_r_max_old-2)
  END IF

  !----- Transform old data to cheb space:
  !      Note: i_costf_init_old,d_costf_init_old used here!
  CALL costf1(dataR,2,1,2,work,i_costf_init_old,d_costf_init_old)

  !----- Fill up cheb polynomial with zeros:
  IF ( n_r_max.GT.n_r_max_old ) THEN
     DO nR=n_r_max_old+1,n_r_max
        dataR(nR)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
     END DO
  END IF

  !----- Now transform to new radial grid points:
  !      Note: i_costf_init,d_costf_init used here!
  CALL costf1(dataR,2,1,2,work,i_costf_init,d_costf_init)

  !----- Rescale :
  cheb_norm_old=DSQRT(2.D0/DBLE(n_r_max_old-1))
  scale=cheb_norm_old/cheb_norm
  DO nR=1,n_r_max
     dataR(nR)=scale*dataR(nR)
  END DO


  deallocate( i_costf_init_old )
  deallocate( d_costf_init_old )
  deallocate( work )

END SUBROUTINE mapDataR

!---------------------------------------------------------------------
