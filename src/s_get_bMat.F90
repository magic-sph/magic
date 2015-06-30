!$Id$
!***********************************************************************
#ifdef WITH_PRECOND_BJ
SUBROUTINE get_bMat(dt,l,hdif,bMat,bPivot,bMat_fac,jMat,jPivot,jMat_fac)
#else
SUBROUTINE get_bMat(dt,l,hdif,bMat,bPivot,jMat,jPivot)
#endif
  !***********************************************************************

  !  +-------------+----------------+------------------------------------+
  !  |                                                                   |
  !  |  Purpose of this subroutine is to contruct the time step matricies|
  !  |  bmat(i,j) and ajmat for the dynamo equations.                    |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+
  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE init_fields
  USE blocking
  USE logic
  USE Bext
#ifdef WITH_MKL_LU
  USE lapack95, ONLY: getrf
#else
  USE algebra, ONLY: sgefa
#endif

  IMPLICIT NONE

  REAL(kind=8),INTENT(IN) :: dt
  INTEGER,INTENT(IN) :: l
  REAL(kind=8),INTENT(IN) :: hdif

  !-- Output:
  REAL(kind=8),INTENT(OUT) :: bMat(n_r_totMag,n_r_totMag)
  INTEGER,INTENT(OUT) :: bPivot(n_r_totMag)
  REAL(kind=8),INTENT(OUT) :: jMat(n_r_totMag,n_r_totMag)
  INTEGER,INTENT(OUT) :: jPivot(n_r_totMag)
#ifdef WITH_PRECOND_BJ
  REAL(kind=8),INTENT(OUT) :: bMat_fac(n_r_totMag),jMat_fac(n_r_totMag)
#endif

  !-- local variables:
  INTEGER :: nR,nCheb,nRall
  INTEGER :: info
  REAL(kind=8) :: l_P_1
  REAL(kind=8) :: O_dt,dLh
  REAL(kind=8) :: rRatio

#undef MATRIX_CHECK
#ifdef MATRIX_CHECK
  INTEGER :: i,j
  real(kind=8) :: rcond
  INTEGER ::ipiv(n_r_tot),iwork(n_r_tot)
  REAL(kind=8) :: work(4*n_r_tot),anorm,linesum
  REAL(kind=8) :: temp_Mat(n_r_tot,n_r_tot)
  INTEGER,save :: counter=0
  integer :: filehandle
  character(len=100) :: filename
#endif

  !-- end of declaration
  !-----------------------------------------------------------------------


  nRall=n_r_max
  IF ( l_cond_ic ) nRall=nRall+n_r_ic_max
  O_dt=1.D0/dt
  dLh=DBLE(l*(l+1))

  !-- matricies depend on degree l but not on order m,
  !   we thus have to construct bmat and ajmat for each l:

  !-- do loop limits introduced to get rid of compiler warning !

  l_P_1=DBLE(l+1)

  DO nR=2,n_r_max-1
     DO nCheb=1,n_r_max
        bMat(nR,nCheb)=                       cheb_norm * ( &
                          O_dt*dLh*or2(nR)*cheb(nCheb,nR) - &
                  alpha*opm*lambda(nR)*hdif*dLh*or2(nR) * ( &
                                         d2cheb(nCheb,nR) - &
                               dLh*or2(nR)*cheb(nCheb,nR) ) )

        jMat(nR,nCheb)=                       cheb_norm * ( &
                          O_dt*dLh*or2(nR)*cheb(nCheb,nR) - &
                  alpha*opm*lambda(nR)*hdif*dLh*or2(nR) * ( &
                                         d2cheb(nCheb,nR) + &
                             dLlambda(nR)*dcheb(nCheb,nR) - &
                               dLh*or2(nR)*cheb(nCheb,nR) ) )
     END DO
  END DO

  IF  ( l_LCR ) THEN

     DO nR=2,n_r_max-1
        IF ( nR<=n_r_LCR ) THEN
           DO nCheb=1,n_r_max
              bMat(nR,nCheb)= cheb_norm*(     dcheb(nCheb,nR) + &
                               DBLE(l)*or1(nR)*cheb(nCheb,nR) ) 

              jMat(nR,nCheb)= cheb_norm*cheb(nCheb,nR)
           END DO
        END IF

     END DO

  END IF

  !----- boundary conditions for outer core field:
  DO nCheb=1,n_cheb_max

!-- JW 10.Apr.2014: pseudo vacuum outer boundary condition included as ktopb=4
!option.
     IF ( ktopb == 1 .OR. ktopb == 3 ) THEN

     !-------- at CMB (nR=1):
     !         the internal poloidal field should fit a potential
     !         field (matrix bmat) and the toroidal field has to
     !         vanish (matrix ajmat).

        bMat(1,nCheb)=            cheb_norm * ( &
                               dcheb(nCheb,1) + &
                 DBLE(l)*or1(1)*cheb(nCheb,1) + &
                              conductance_ma* ( &
                              d2cheb(nCheb,1) - &
                     dLh*or2(1)*cheb(nCheb,1) ) )

        jMat(1,nCheb)=            cheb_norm * ( &
                                cheb(nCheb,1) + &
                conductance_ma*dcheb(nCheb,1) )

     ELSE IF ( ktopb == 2 ) THEN

        WRITE(*,*) '! Boundary condition ktopb=2 not defined!'
        STOP

     ELSE IF ( ktopb == 4 ) THEN

!----- pseudo vacuum condition, field has only
!      a radial component, horizontal components
!      vanish when aj and db are zero:
        bMat(1,nCheb)=cheb_norm*dcheb(nCheb,1)

        jMat(1,nCheb)=cheb_norm*cheb(nCheb,1)

     END IF



     !-------- at IC (nR=n_r_max):
     IF ( kbotb == 1 ) THEN

        !----------- insulating IC, field has to fit a potential field:
        bMat(n_r_max,nCheb)=       cheb_norm * ( &
                          dcheb(nCheb,n_r_max) - &
             l_P_1*or1(n_r_max)*cheb(nCheb,n_r_max) )
        jMat(n_r_max,nCheb)=       cheb_norm*cheb(nCheb,n_r_max)

     ELSE IF ( kbotb == 2 ) THEN

        !----------- perfect conducting IC
        bMat(n_r_max-1,nCheb)=cheb_norm*d2cheb(nCheb,n_r_max)
        jMat(n_r_max,nCheb)  =cheb_norm* dcheb(nCheb,n_r_max)

     ELSE IF ( kbotb == 3 ) THEN

        !---------- finite conducting IC, four boundary conditions:
        !           continuity of b,j, (d b)/(d r) and (d j)/(d r)/sigma.
        !           note: n_r=n_r_max and n_r=n_r_max+1 stand for IC radius
        !           here we set the outer core part of the equations.
        !           the conductivity ratio sigma_ratio is used as
        !           an additional dimensionless parameter.

        bMat(n_r_max,nCheb)=  cheb_norm*cheb(nCheb,n_r_max)
        bMat(n_r_max+1,nCheb)=cheb_norm*dcheb(nCheb,n_r_max)
        jMat(n_r_max,nCheb)=  cheb_norm*cheb(nCheb,n_r_max)
        jMat(n_r_max+1,nCheb)=cheb_norm*sigma_ratio*dcheb(nCheb,n_r_max)

!-- JW 10.Apr.2012: pseudo vacuum condition at lower boundary included.
     ELSE IF ( kbotb == 4 ) THEN

!----- Pseudovacuum conduction at lower boundary:
        bMat(n_r_max,nCheb)=cheb_norm*dcheb(nCheb,n_r_max)
        jMat(n_r_max,nCheb)=cheb_norm*cheb(nCheb,n_r_max)

     END IF

     !-------- Imposed fields: (overwrites above IC boundary cond.)
     IF ( l == 1 .AND. ( imagcon == -1 .OR. imagcon == -2 ) ) THEN
        bMat(n_r_max,nCheb)=  cheb_norm * cheb(nCheb,n_r_max)

     ELSE IF ( l == 3 .AND. imagcon == -10 ) THEN
        IF ( l_LCR ) THEN
           WRITE(*,*) 'Imposed field not compatible with weak conducting region!'
           STOP
        END IF
        jMat(1,nCheb)      =  cheb_norm * cheb(nCheb,1)
        jMat(n_r_max,nCheb)=  cheb_norm * cheb(nCheb,n_r_max)

     ELSE IF ( n_imp == 1 ) THEN
        !-- This is the Uli Christensen idea where the external field is
        !   not fixed but compensates the internal field so that the
        !   radial field component vanishes at r/r_cmb=rrMP
        IF ( l_LCR ) THEN
           WRITE(*,*) 'Imposed field not compatible with weak conducting region!'
           STOP
        END IF
        rRatio=rrMP**DBLE(2*l+1)
        bMat(1,nCheb)=            cheb_norm * ( &
                               dcheb(nCheb,1) + &
                 DBLE(l)*or1(1)*cheb(nCheb,1) - &
                DBLE(2*l+1)*or1(1)/(1-rRatio) + &
                              conductance_ma* ( &
                              d2cheb(nCheb,1) - &
                     dLh*or2(1)*cheb(nCheb,1) ) )
     END IF

  END DO ! loop over cheb modes !

  !----- fill up with zeros:
  DO nCheb=n_cheb_max+1,n_r_max
     bMat(1,nCheb)=0.D0
     jMat(1,nCheb)=0.D0
     IF ( l_LCR ) THEN
        DO nR=2,n_r_LCR
           bMat(nR,nCheb)=0.D0
           jMat(nR,nCheb)=0.D0
        END DO
     END IF
     IF ( kbotb == 1 ) THEN
        bMat(n_r_max,nCheb)  =0.D0
        jMat(n_r_max,nCheb)  =0.D0
     ELSE IF ( kbotb == 2 ) THEN
        bMat(n_r_max-1,nCheb)=0.D0
        jMat(n_r_max,nCheb)  =0.D0
     ELSE IF ( kbotb == 3 ) THEN
        bMat(n_r_max,nCheb)  =0.D0
        bMat(n_r_max+1,nCheb)=0.D0
        jMat(n_r_max,nCheb)  =0.D0
        jMat(n_r_max+1,nCheb)=0.D0
     ELSE IF ( kbotb == 4 ) THEN
        bMat(n_r_max,nCheb)  =0.D0
        jMat(n_r_max,nCheb)  =0.D0
     END IF
  END DO

  !----- normalization for highest and lowest Cheb mode:
  DO nR=1,n_r_max
     bMat(nR,1)      =0.5D0*bMat(nR,1)
     bMat(nR,n_r_max)=0.5D0*bMat(nR,n_r_max)
     jMat(nR,1)      =0.5D0*jMat(nR,1)
     jMat(nR,n_r_max)=0.5D0*jMat(nR,n_r_max)
  END DO
  IF ( kbotb == 3 ) THEN
     bMat(n_r_max+1,1)=0.5D0*bMat(n_r_max+1,1)
     bMat(n_r_max+1,n_r_max) = &
          0.5D0*bMat(n_r_max+1,n_r_max)
     jMat(n_r_max+1,1)=0.5D0*jMat(n_r_max+1,1)
     jMat(n_r_max+1,n_r_max)= &
          0.5D0*jMat(n_r_max+1,n_r_max)
  END IF

  !----- Conducting inner core:
  IF ( l_cond_ic ) THEN

     !----- rnner core implicit time step matricies for the grid
     !      points n_r=n_r_max+1,...,n_r_max+n_r_ic

     DO nCheb=1,n_r_ic_max ! counts even IC cheb modes
        DO nR=2,n_r_ic_max-1 ! counts IC radial grid points
           ! n_r=1 represents ICB
           !----------- poloidal field matrix for an inner core field
           !            of the radial form: (r/r_icb)**(l+1)*cheb_ic(r)
           !            where cheb_ic are even chebs only.
           !            NOTE: no hyperdiffusion in inner core !

           bMat(n_r_max+nR,n_r_max+nCheb) = &
                cheb_norm_ic*dLh*or2(n_r_max) * ( &
                O_dt*cheb_ic(nCheb,nR) - &
                alpha*opm*O_sr * ( &
                d2cheb_ic(nCheb,nR) + &
                2.d0*l_P_1*O_r_ic(nR)*dcheb_ic(nCheb,nR) )   )

           jMat(n_r_max+nR,n_r_max+nCheb) = &
                bMat(n_r_max+nR,n_r_max+nCheb)
        END DO

        !----- Special treatment for r=0, asymptotic of 1/r dr
        nR=n_r_ic_max
        bMat(n_r_max+nR,n_r_max+nCheb) = &
             cheb_norm_ic*dLh*or2(n_r_max) * ( &
             O_dt*cheb_ic(nCheb,nR) - &
             alpha*opm*O_sr * &
             (1.d0+2.d0*l_P_1)*d2cheb_ic(nCheb,nR) )

        jMat(n_r_max+nR,n_r_max+nCheb) = &
             bMat(n_r_max+nR,n_r_max+nCheb)

     END DO

     !-------- boundary condition at r_icb:
     DO nCheb=1,n_cheb_ic_max
        bMat(n_r_max,n_r_max+nCheb)= &
             -cheb_norm_ic*cheb_ic(nCheb,1)
        bMat(n_r_max+1,n_r_max+nCheb)= &
             -cheb_norm_ic * ( &
             dcheb_ic(nCheb,1) + &
             l_P_1*or1(n_r_max)*cheb_ic(nCheb,1) )
        jMat(n_r_max,n_r_max+nCheb)= &
             bMat(n_r_max,n_r_max+nCheb)
        jMat(n_r_max+1,n_r_max+nCheb)= &
             bMat(n_r_max+1,n_r_max+nCheb)
     END DO ! cheb modes

     !-------- fill with zeros:
     DO nCheb=n_cheb_ic_max+1,n_r_ic_max
        bMat(n_r_max,n_r_max+nCheb)  =0.D0
        bMat(n_r_max+1,n_r_max+nCheb)=0.D0
        jMat(n_r_max,n_r_max+nCheb)  =0.D0
        jMat(n_r_max+1,n_r_max+nCheb)=0.D0
     END DO

     !-------- normalization for lowest Cheb mode:
     DO nR=n_r_max,n_r_tot
        bMat(nR,n_r_max+1)=0.5D0*bMat(nR,n_r_max+1)
        jMat(nR,n_r_max+1)=0.5D0*jMat(nR,n_r_max+1)
        bMat(nR,n_r_tot)  =0.5D0*bMat(nR,n_r_tot)
        jMat(nR,n_r_tot)  =0.5D0*jMat(nR,n_r_tot)
     END DO

     !-------- fill matricies up with zeros:
     DO nCheb=n_r_max+1,n_r_tot
        DO nR=1,n_r_max-1
           bMat(nR,nCheb)=0.D0
           jMat(nR,nCheb)=0.D0
        END DO
     END DO
     DO nCheb=1,n_r_max
        DO nR=n_r_max+2,n_r_tot
           bMat(nR,nCheb)=0.D0
           jMat(nR,nCheb)=0.D0
        END DO
     END DO

  END IF  ! conducting inner core ?

#ifdef WITH_PRECOND_BJ
  ! compute the linesum of each line
  DO nR=1,n_r_tot
     bMat_fac(nR)=1.0D0/MAXVAL(ABS(bMat(nR,:)))
     bMat(nR,:) = bMat(nR,:)*bMat_fac(nR)
  END DO
  DO nR=1,n_r_tot
     jMat_fac(nR)=1.0D0/MAXVAL(ABS(jMat(nR,:)))
     jMat(nR,:) = jMat(nR,:)*jMat_fac(nR)
  END DO
#endif

#ifdef MATRIX_CHECK
  ! copy the bMat to a temporary variable for modification
  WRITE(filename,"(A,I3.3,A,I3.3,A)") "bMat_",l,"_",counter,".dat"
  OPEN(NEWUNIT=filehandle,file=TRIM(filename))
  counter= counter+1
  
  DO i=1,n_r_tot
     DO j=1,n_r_tot
        WRITE(filehandle,"(2ES20.12,1X)",advance="no") bMat(i,j)
     END DO
     WRITE(filehandle,"(A)") ""
  END DO
  CLOSE(filehandle)
  temp_Mat=bMat
  anorm = 0.0D0
  DO i=1,n_r_tot
     linesum = 0.0D0
     DO j=1,n_r_tot
        linesum = linesum + ABS(temp_Mat(i,j))
     END DO
     IF (linesum .GT. anorm) anorm=linesum
  END DO
  !WRITE(*,"(A,ES20.12)") "anorm = ",anorm
  ! LU factorization
  CALL dgetrf(n_r_tot,n_r_tot,temp_Mat,n_r_tot,ipiv,info)
  ! estimate the condition number
  CALL dgecon('I',n_r_tot,temp_Mat,n_r_tot,anorm,rcond,work,iwork,info)
  WRITE(*,"(A,I3,A,ES11.3)") "inverse condition number of bMat for l=",l," is ",rcond
  
  ! The same computation for jMat.
  ! copy the jMat to a temporary variable for modification
  WRITE(filename,"(A,I3.3,A,I3.3,A)") "jMat_",l,"_",counter,".dat"
  OPEN(NEWUNIT=filehandle,file=TRIM(filename))
  counter= counter+1
  
  DO i=1,n_r_tot
     DO j=1,n_r_tot
        WRITE(filehandle,"(2ES20.12,1X)",advance="no") jMat(i,j)
     END DO
     WRITE(filehandle,"(A)") ""
  END DO
  CLOSE(filehandle)
  temp_Mat=jMat
  anorm = 0.0D0
  DO i=1,n_r_tot
     linesum = 0.0D0
     DO j=1,n_r_tot
        linesum = linesum + ABS(temp_Mat(i,j))
     END DO
     IF (linesum .GT. anorm) anorm=linesum
  END DO
  !WRITE(*,"(A,ES20.12)") "anorm = ",anorm
  ! LU factorization
  CALL dgetrf(n_r_tot,n_r_tot,temp_Mat,n_r_tot,ipiv,info)
  ! estimate the condition number
  CALL dgecon('I',n_r_tot,temp_Mat,n_r_tot,anorm,rcond,work,iwork,info)
  WRITE(*,"(A,I3,A,ES11.3)") "inverse condition number of jMat for l=",l," is ",rcond
#endif

  !----- LU decomposition:
#ifdef WITH_MKL_LU
  CALL getrf(bMat(1:nRall,1:nRall),bPivot(1:nRall),info)
#else
  CALL sgefa(bMat,n_r_tot,nRall,bPivot,info)
#endif

  IF ( info /= 0 ) THEN
     WRITE(*,*) 'Singular matrix bmat in get_bmat.'
     STOP '32'
  END IF

  !----- LU decomposition:
#ifdef WITH_MKL_LU
  CALL getrf(jMat(1:nRall,1:nRall),jPivot(1:nRall),info)
#else
  CALL sgefa(jMat,n_r_tot,nRall,jPivot,info)
#endif
  IF ( info /= 0 ) THEN
     WRITE(*,*) '! Singular matrix ajmat in get_bmat!'
     STOP '33'
  END IF


  RETURN
end SUBROUTINE get_bMat

!-----------------------------------------------------------------------------
