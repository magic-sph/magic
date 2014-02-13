!$Id$
!***********************************************************************
SUBROUTINE get_wpMat(dt,l,hdif,wpMat,wpPivot,wpMat_fac)
  !***********************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/17/02  by JW. --------------!

  !  +-------------+----------------+------------------------------------+
  !  |                                                                   |
  !  |  Purpose of this subroutine is to contruct the time step matricies|
  !  |  zmat(i,j) and wpmat  for the NS equation.                        |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+

  USE truncation,ONLY: n_cheb_max,n_r_max
  USE radial_functions,ONLY: cheb_norm,cheb,dcheb,d2cheb,d3cheb,&
       & or1,or2,beta,visc,dLvisc,dbeta
  USE physical_parameters, ONLY: ktopv, kbotv
  USE num_param, only: alpha
  USE algebra, ONLY: sgefa

  IMPLICIT NONE

  REAL(kind=8),INTENT(IN) :: dt
  INTEGER,INTENT(IN) :: l
  REAL(kind=8),INTENT(IN) :: hdif

  !-- Output:
  REAL(kind=8),INTENT(OUT) :: wpMat(2*n_r_max,2*n_r_max)
  INTEGER,INTENT(OUT) :: wpPivot(2*n_r_max)
  REAL(kind=8),INTENT(OUT) :: wpMat_fac(2*n_r_max,2)

  !-- local variables:
  INTEGER :: nR,nCheb,nR_p,nCheb_p
  INTEGER :: info
  REAL(kind=8) :: O_dt,dLh

#ifdef MATRIX_CHECK
  INTEGER ::ipiv(2*n_r_max),iwork(2*n_r_max),i,j
  REAL(kind=8) :: work(8*n_r_max),anorm,linesum,rcond
  REAL(kind=8) :: temp_wpMat(2*n_r_max,2*n_r_max)
  INTEGER,save :: counter=0
  integer :: filehandle
  character(len=100) :: filename
  logical :: first_run=.true.
#endif
  !-- end of declaration
  !-----------------------------------------------------------------------

#if 0
  IF (first_run) THEN
     OPEN(NEWUNIT=filehandle,file="cheb.dat")
     DO nR=1,n_r_max
        DO nCheb=1,n_r_max
           WRITE(filehandle,"(ES20.12)",advance='no') cheb(nCheb,nR)
        END DO
        WRITE(filehandle,"(A)") ""
     END DO
     CLOSE(filehandle)
     OPEN(NEWUNIT=filehandle,file="dcheb.dat")
     DO nR=1,n_r_max
        DO nCheb=1,n_r_max
           WRITE(filehandle,"(ES20.12)",advance='no') dcheb(nCheb,nR)
        END DO
        WRITE(filehandle,"(A)") ""
     END DO
     CLOSE(filehandle)
     OPEN(NEWUNIT=filehandle,file="d2cheb.dat")
     DO nR=1,n_r_max
        DO nCheb=1,n_r_max
           WRITE(filehandle,"(ES20.12)",advance='no') d2cheb(nCheb,nR)
        END DO
        WRITE(filehandle,"(A)") ""
     END DO
     CLOSE(filehandle)
     OPEN(NEWUNIT=filehandle,file="d3cheb.dat")
     DO nR=1,n_r_max
        DO nCheb=1,n_r_max
           WRITE(filehandle,"(ES20.12)",advance='no') d3cheb(nCheb,nR)
        END DO
        WRITE(filehandle,"(A)") ""
     END DO
     CLOSE(filehandle)
     first_run=.false.
  END IF
#endif

  O_dt=1.D0/dt
  dLh =DBLE(l*(l+1))

  !-- Now mode l>0

  !----- Boundary conditions, see above:
  DO nCheb=1,n_cheb_max
     nCheb_p=nCheb+n_r_max

     wpMat(1,nCheb)        =cheb_norm*cheb(nCheb,1)
     wpMat(1,nCheb_p)      =0.D0
     wpMat(n_r_max,nCheb)  =cheb_norm*cheb(nCheb,n_r_max)
     wpMat(n_r_max,nCheb_p)=0.D0

     IF ( ktopv == 1 ) then  ! free slip !
        wpMat(n_r_max+1,nCheb)=   cheb_norm * ( &
             d2cheb(nCheb,1) - &
             (2.d0*or1(1)+beta(1))*dcheb(nCheb,1) )
     ELSE                    ! no slip, note exception for l=1,m=0
        wpMat(n_r_max+1,nCheb)=cheb_norm*dcheb(nCheb,1)
     END IF
     wpMat(n_r_max+1,nCheb_p)=0.D0

     IF ( kbotv == 1 ) then  ! free slip !
        wpMat(2*n_r_max,nCheb)=        cheb_norm * ( &
             d2cheb(nCheb,n_r_max) - &
             (2.d0*or1(n_r_max)+beta(n_r_max))*dcheb(nCheb,n_r_max))
     ELSE                 ! no slip, note exception for l=1,m=0
        wpMat(2*n_r_max,nCheb)=cheb_norm * dcheb(nCheb,n_r_max)
     END IF
     wpMat(2*n_r_max,nCheb_p)=0.D0

  END DO   !  loop over nCheb

  IF ( n_cheb_max < n_r_max ) then ! fill with zeros !
     DO nCheb=n_cheb_max+1,n_r_max
        nCheb_p=nCheb+n_r_max
        wpMat(1,nCheb)          =0.D0
        wpMat(n_r_max,nCheb)    =0.D0
        wpMat(n_r_max+1,nCheb)  =0.D0
        wpMat(2*n_r_max,nCheb)  =0.D0
        wpMat(1,nCheb_p)        =0.D0
        wpMat(n_r_max,nCheb_p)  =0.D0
        wpMat(n_r_max+1,nCheb_p)=0.D0
        wpMat(2*n_r_max,nCheb_p)=0.D0
     END DO
  END IF

  !----- Other points:
  DO nCheb=1,n_r_max
     nCheb_p=nCheb+n_r_max
     DO nR=2,n_r_max-1
        !WRITE(*,"(I3,A,6ES11.3)") nR,", visc,beta,dLvisc,dbeta = ",&
        !     & visc(nR),beta(nR),dLvisc(nR),dbeta(nR),hdif,alpha
        ! in the BM2 case: visc=1.0,beta=0.0,dLvisc=0.0,dbeta=0.0
        nR_p=nR+n_r_max
        wpMat(nR,nCheb)= cheb_norm * &
             & ( O_dt*dLh*or2(nR)*cheb(nCheb,nR) &
             &   - alpha*hdif*visc(nR)*dLh*or2(nR) &
             &     * ( d2cheb(nCheb,nR) &
             &         +(2.*dLvisc(nR)-beta(nR)/3.D0)*dcheb(nCheb,nR) &
             &         -( dLh*or2(nR)&
             &            + 4.D0/3.D0*(dLvisc(nR)*beta(nR) &
             &                         +(3.d0*dLvisc(nR)+beta(nR))*or1(nR)&
             &                         +dbeta(nR)) &
             &          )*cheb(nCheb,nR) &
             &       ) &
             & )

        wpMat(nR,nCheb_p)= cheb_norm * &
             alpha*(dcheb(nCheb,nR) &
             -  beta(nR)* cheb(nCheb,nR))
        ! the following part gives sometimes very large 
        ! matrix entries
        wpMat(nR_p,nCheb)= &
             & cheb_norm * ( &
             &  -O_dt*dLh*or2(nR)*dcheb(nCheb,nR) &
             &  -alpha*hdif*visc(nR)*dLh*or2(nR) &
             &    *(- d3cheb(nCheb,nR) &
             &      +( beta(nR)-dLvisc(nR) )*d2cheb(nCheb,nR) &
             &      +( dLh*or2(nR)&
             &         +dbeta(nR)&
             &         +dLvisc(nR)*beta(nR)&
             &         +2.D0*(dLvisc(nR)+beta(nR))*or1(nR) &
             &       )*dcheb(nCheb,nR) &
             &      -dLh*or2(nR)*( 2.d0*or1(nR)&
             &                     +dLvisc(nR)&
             &                     +2.d0/3.d0*beta(nR) &
             &                   )*cheb(nCheb,nR) &
             &     ) &
             & )

        wpMat(nR_p,nCheb_p)= - cheb_norm * &
             alpha*dLh*or2(nR)*cheb(nCheb,nR)
     END DO
  END DO


  !----- Factor for highest and lowest cheb:
  DO nR=1,n_r_max
     nR_p=nR+n_r_max
     wpMat(nR,1)          =0.5D0*wpMat(nR,1)
     wpMat(nR,n_r_max)    =0.5D0*wpMat(nR,n_r_max)
     wpMat(nR,n_r_max+1)  =0.5D0*wpMat(nR,n_r_max+1)
     wpMat(nR,2*n_r_max)  =0.5D0*wpMat(nR,2*n_r_max)
     wpMat(nR_p,1)        =0.5D0*wpMat(nR_p,1)
     wpMat(nR_p,n_r_max)  =0.5D0*wpMat(nR_p,n_r_max)
     wpMat(nR_p,n_r_max+1)=0.5D0*wpMat(nR_p,n_r_max+1)
     wpMat(nR_p,2*n_r_max)=0.5D0*wpMat(nR_p,2*n_r_max)
  END DO

  ! compute the linesum of each line
  DO nR=1,2*n_r_max
     wpMat_fac(nR,1)=1.0D0/MAXVAL(ABS(wpMat(nR,:)))
  END DO
  ! now divide each line by the linesum to regularize the matrix
  DO nr=1,2*n_r_max
     wpMat(nR,:) = wpMat(nR,:)*wpMat_fac(nR,1)
  END DO

  ! also compute the rowsum of each column
  DO nR=1,2*n_r_max
     wpMat_fac(nR,2)=1.0D0/MAXVAL(ABS(wpMat(:,nR)))
  END DO
  ! now divide each row by the rowsum
  DO nR=1,2*n_r_max
     wpMat(:,nR) = wpMat(:,nR)*wpMat_fac(nR,2)
  END DO

#ifdef MATRIX_CHECK
     ! copy the wpMat to a temporary variable for modification
     WRITE(filename,"(A,I3.3,A,I3.3,A)") "wpMat_",l,"_",counter,".dat"
     OPEN(NEWUNIT=filehandle,file=TRIM(filename))
     counter= counter+1
     
     DO i=1,2*n_r_max
        DO j=1,2*n_r_max
           WRITE(filehandle,"(2ES20.12,1X)",advance="no") wpMat(i,j)
        END DO
        WRITE(filehandle,"(A)") ""
     END DO
     CLOSE(filehandle)
     temp_wpMat=wpMat
     anorm = 0.0D0
     DO i=1,2*n_r_max
        linesum = 0.0D0
        DO j=1,2*n_r_max
           linesum = linesum + ABS(temp_wpMat(i,j))
        END DO
        IF (linesum .GT. anorm) anorm=linesum
     END DO
     WRITE(*,"(A,ES20.12)") "anorm = ",anorm
     ! LU factorization
     CALL dgetrf(2*n_r_max,2*n_r_max,temp_wpMat,2*n_r_max,ipiv,info)
     ! estimate the condition number
     CALL dgecon('I',2*n_r_max,temp_wpMat,2*n_r_max,anorm,rcond,work,iwork,info)
     WRITE(*,"(A,I3,A,ES11.3)") "inverse condition number of wpMat for l=",l," is ",rcond
#endif

  ! this is the original code
  CALL sgefa(wpMat,2*n_r_max,2*n_r_max,wpPivot,info)
  IF ( info /= 0 ) THEN
     WRITE(*,*) 'Singular matrix wpmat!'
     STOP '35'
  END IF


  RETURN
end SUBROUTINE get_wpMat

!-----------------------------------------------------------------------------
