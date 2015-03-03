!$Id$
!***********************************************************************
#ifdef WITH_PRECOND_S
  SUBROUTINE get_Smat(dt,l,hdif,sMat,sPivot,sMat_fac)
#else
  SUBROUTINE get_Smat(dt,l,hdif,sMat,sPivot)
#endif
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  |  Purpose of this subroutine is to contruct the time step matricies|
!  |  sMat(i,j) and s0mat for the entropy equation.                    |
!  |                                                                   |
!  +-------------------------------------------------------------------+
      
    USE truncation
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE algebra, ONLY: sgefa
#ifdef WITH_MKL_LU
    use lapack95
#endif
    IMPLICIT NONE

    REAL(kind=8) :: dt
    INTEGER :: l
    REAL(kind=8) :: hdif

!-- Output
    REAL(kind=8),INTENT(OUT) :: sMat(n_r_max,n_r_max)
    INTEGER,INTENT(OUT) :: sPivot(n_r_max)
#ifdef WITH_PRECOND_S
    REAL(kind=8),INTENT(out) :: sMat_fac(n_r_max)
#endif
!-- Local variables:
    INTEGER :: info,nCheb,nR
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

!-- end of declaration
!-----------------------------------------------------------------------

    O_dt=1.D0/dt

    dLh=DBLE(l*(l+1))

!----- Boundary coditions:
    DO nCheb=1,n_cheb_max
        IF ( ktops == 1 ) then
            sMat(1,nCheb)=cheb_norm
        ELSE
            sMat(1,nCheb)=cheb_norm*dcheb(nCheb,1)
        END IF
        IF ( kbots == 1 ) THEN
            sMat(n_r_max,nCheb)=cheb_norm*cheb(nCheb,n_r_max)
        ELSE
            sMat(n_r_max,nCheb)=cheb_norm*dcheb(nCheb,n_r_max)
        END IF
    END DO
    IF ( n_cheb_max < n_r_max ) THEN ! fill with zeros !
        DO nCheb=n_cheb_max+1,n_r_max
            sMat(1,nCheb)      =0.D0
            sMat(n_r_max,nCheb)=0.D0
        END DO
    END IF

!----- Other points:
    DO nCheb=1,n_r_max
        DO nR=2,n_r_max-1
            sMat(nR,nCheb)= cheb_norm * (                    &
                                       O_dt*cheb(nCheb,nR) - &
              alpha*opr*hdif*kappa(nR)*(  d2cheb(nCheb,nR) + &
              ( beta(nR)+otemp1(nR)*dtemp0(nR)+              &
                2.d0*or1(nR)+dLkappa(nR) )*dcheb(nCheb,nR) - &
                   dLh*or2(nR)*             cheb(nCheb,nR) ) )
        END DO
    END DO

!----- Factor for highest and lowest cheb:
    DO nR=1,n_r_max
        sMat(nR,1)      =0.5D0*sMat(nR,1)
        sMat(nR,n_r_max)=0.5D0*sMat(nR,n_r_max)
    END DO

#ifdef WITH_PRECOND_S
    ! compute the linesum of each line
    DO nR=1,n_r_max
       sMat_fac(nR)=1.0D0/MAXVAL(ABS(sMat(nR,:)))
    END DO
    ! now divide each line by the linesum to regularize the matrix
    DO nr=1,n_r_max
       sMat(nR,:) = sMat(nR,:)*sMat_fac(nR)
    END DO
#endif

#ifdef MATRIX_CHECK
     ! copy the sMat to a temporary variable for modification
     WRITE(filename,"(A,I3.3,A,I3.3,A)") "sMat_",l,"_",counter,".dat"
     OPEN(NEWUNIT=filehandle,file=TRIM(filename))
     counter= counter+1

     DO i=1,n_r_max
        DO j=1,n_r_max
           WRITE(filehandle,"(2ES20.12,1X)",advance="no") sMat(i,j)
        END DO
        WRITE(filehandle,"(A)") ""
     END DO
     CLOSE(filehandle)
     temp_Mat=sMat
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
     WRITE(*,"(A,I3,A,ES11.3)") "inverse condition number of sMat for l=",l," is ",rcond
#endif

!----- LU decomposition:
#ifdef WITH_MKL_LU
     CALL getrf(sMat,sPivot,info)
     !CALL dgetrf(n_r_max,n_r_max,sMat,n_r_max,sPivot,info)
#else
    CALL sgefa(sMat,n_r_max,n_r_max,sPivot,info)
#endif
    IF ( info .ne. 0 ) THEN
        WRITE(*,*) 'Singular matrix sMat!'
        STOP '31'
    END IF
            
    RETURN
    end SUBROUTINE get_Smat

!-----------------------------------------------------------------------------
