!$Id$
!***********************************************************************
SUBROUTINE get_zMat(dt,l,hdif,zMat,zPivot)
  !***********************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/17/02  by JW. --------------!

  !  +-------------+----------------+------------------------------------+
  !  |                                                                   |
  !  |  Purpose of this subroutine is to contruct the time step matricies|
  !  |  zmat(i,j) for the NS equation.                                   |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+

  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE blocking
  USE logic
  USE algebra, ONLY: sgefa

  IMPLICIT NONE

  REAL(kind=8),intent(IN) :: dt
  INTEGER,intent(IN) :: l
  REAL(kind=8),intent(IN) :: hdif

  !-- Output:
  REAL(kind=8),intent(OUT) :: zMat(n_r_max,n_r_max)
  INTEGER,intent(OUT) :: zPivot(n_r_max)

  !-- local variables:
  INTEGER :: nR,nCheb
  INTEGER :: info
  REAL(kind=8) :: O_dt,dLh

  !-- End of declaration
  !-----------------------------------------------------------------------

  O_dt=1.D0/dt
  dLh=DBLE(l*(l+1))

  !----- Boundary conditions, see above:
  DO nCheb=1,n_cheb_max
     IF ( ktopv == 1 ) THEN  ! free slip !
        zMat(1,nCheb)=cheb_norm *             ( &
             dcheb(nCheb,1) - &
             (2.d0*or1(1)+beta(1))*cheb(nCheb,1) )
     ELSE                    ! no slip, note exception for l=1,m=0
        zMat(1,nCheb)=cheb_norm*cheb(nCheb,1)
     END IF

     IF ( kbotv == 1 ) THEN  ! free slip !
        zMat(n_r_max,nCheb)= cheb_norm *            ( &
             dcheb(nCheb,n_r_max) - &
             (2.d0*or1(n_r_max)+beta(n_r_max))*cheb(nCheb,n_r_max) )
     ELSE                    ! no slip, note exception for l=1,m=0
        zMat(n_r_max,nCheb)= cheb_norm * cheb(nCheb,n_r_max)
     END IF

  END DO  !  loop over nCheb

  IF ( n_cheb_max < n_r_max ) THEN ! fill with zeros !
     DO nCheb=n_cheb_max+1,n_r_max
        zMat(1,nCheb)      =0.D0
        zMat(n_r_max,nCheb)=0.D0
     END DO
  END IF

  !----- Other points:
  DO nCheb=1,n_r_max
     DO nR=2,n_r_max-1
        zMat(nR,nCheb)=                   cheb_norm * ( &
             O_dt*dLh*or2(nR)*cheb(nCheb,nR) - &
             alpha*hdif*dLh*visc(nR)*or2(nR) * ( &
             d2cheb(nCheb,nR) + &
             (dLvisc(nR)- beta(nR))*      dcheb(nCheb,nR) - &
             (dLvisc(nR)*beta(nR)+2.d0*dLvisc(nR)*or1(nR)  + &
             dLh*or2(nR)+dbeta(nR)+2.D0*beta(nR)*or1(nR))* &
             cheb(nCheb,nR) ) )
     END DO
  END DO

  !----- Factor for highest and lowest cheb:
  DO nR=1,n_r_max
     zMat(nR,1)      =0.5D0*zMat(nR,1)
     zMat(nR,n_r_max)=0.5D0*zMat(nR,n_r_max)
  END DO

  !----- LU decomposition:
  CALL sgefa(zMat,n_r_max,n_r_max,zPivot,info)
  IF ( info /= 0 ) THEN
     WRITE(*,*) 'Singular matrix zmat!'
     STOP '34'
  END IF


  RETURN
end SUBROUTINE get_zMat

!-----------------------------------------------------------------------------
