!$Id$
!***********************************************************************
#ifdef WITH_PRECOND_Z
SUBROUTINE get_zMat(dt,l,hdif,zMat,zPivot,zMat_fac)
#else
SUBROUTINE get_zMat(dt,l,hdif,zMat,zPivot)
#endif
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
#ifdef WITH_MKL_LU
  USE lapack95, ONLY: getrf
#else
  USE algebra, ONLY: sgefa
#endif

  IMPLICIT NONE

  REAL(kind=8),intent(IN) :: dt
  INTEGER,intent(IN) :: l
  REAL(kind=8),intent(IN) :: hdif

  !-- Output:
  REAL(kind=8),intent(OUT) :: zMat(n_r_max,n_r_max)
  INTEGER,intent(OUT) :: zPivot(n_r_max)
#ifdef WITH_PRECOND_Z
  REAL(kind=8),intent(out) :: zMat_fac(n_r_max)
#endif

  !-- local variables:
  INTEGER :: nR,nCheb
  INTEGER :: info
  REAL(kind=8) :: O_dt,dLh

#ifdef MATRIX_CHECK
  INTEGER :: i,j
  real(kind=8) :: rcond
  INTEGER ::ipiv(n_r_max),iwork(n_r_max)
  REAL(kind=8) :: work(4*n_r_max),anorm,linesum
  REAL(kind=8) :: temp_Mat(n_r_max,n_r_max)
  INTEGER,save :: counter=0
  integer :: filehandle
  character(len=100) :: filename
#endif
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
        zMat(nR,nCheb)= cheb_norm * ( &
             &   O_dt*dLh*or2(nR)*cheb(nCheb,nR) &
             &   -alpha*hdif*dLh*visc(nR)*or2(nR) * ( &
             &      d2cheb(nCheb,nR) &
             &      + (dLvisc(nR)- beta(nR)) * dcheb(nCheb,nR) &
             &      - ( dLvisc(nR)*beta(nR)&
             &          +2.d0*dLvisc(nR)*or1(nR) &
             &          +dLh*or2(nR)&
             &          +dbeta(nR)&
             &          +2.D0*beta(nR)*or1(nR) &
             &        ) * cheb(nCheb,nR) ) )
     END DO
  END DO

  !----- Factor for highest and lowest cheb:
  DO nR=1,n_r_max
     zMat(nR,1)      =0.5D0*zMat(nR,1)
     zMat(nR,n_r_max)=0.5D0*zMat(nR,n_r_max)
  END DO

#ifdef WITH_PRECOND_Z
  ! compute the linesum of each line
  DO nR=1,n_r_max
     zMat_fac(nR)=1.0D0/MAXVAL(ABS(zMat(nR,:)))
     zMat(nR,:) = zMat(nR,:)*zMat_fac(nR)
  END DO
#endif

#ifdef MATRIX_CHECK
     ! copy the zMat to a temporary variable for modification
     WRITE(filename,"(A,I3.3,A,I3.3,A)") "zMat_",l,"_",counter,".dat"
     OPEN(NEWUNIT=filehandle,file=TRIM(filename))
     counter= counter+1

     DO i=1,n_r_max
        DO j=1,n_r_max
           WRITE(filehandle,"(2ES20.12,1X)",advance="no") zMat(i,j)
        END DO
        WRITE(filehandle,"(A)") ""
     END DO
     CLOSE(filehandle)
     temp_Mat=zMat
     anorm = 0.0D0
     DO i=1,n_r_max
        linesum = 0.0D0
        DO j=1,n_r_max
           linesum = linesum + ABS(temp_Mat(i,j))
        END DO
        IF (linesum .GT. anorm) anorm=linesum
     END DO
     !WRITE(*,"(A,ES20.12)") "anorm = ",anorm
     ! LU factorization
     CALL dgetrf(n_r_max,n_r_max,temp_Mat,n_r_max,ipiv,info)
     ! estimate the condition number
     CALL dgecon('I',n_r_max,temp_Mat,n_r_max,anorm,rcond,work,iwork,info)
     WRITE(*,"(A,I3,A,ES11.3)") "inverse condition number of zMat for l=",l," is ",rcond
#endif

  !----- LU decomposition:
#ifdef WITH_MKL_LU
  CALL getrf(zMat,zPivot,info)
#else
  CALL sgefa(zMat,n_r_max,n_r_max,zPivot,info)
#endif

  IF ( info /= 0 ) THEN
     WRITE(*,*) 'Singular matrix zmat for l=',l,", info = ",info
     STOP '34'
  END IF


  RETURN
end SUBROUTINE get_zMat

!-----------------------------------------------------------------------------
