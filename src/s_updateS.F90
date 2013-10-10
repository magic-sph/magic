!$Id$
!*************************************************************************
SUBROUTINE updateS(s,ds,dVSrLM,dsdt,dsdtLast, &
     &             w1,coex,dt,nLMB)
  !*************************************************************************

  !-------------------------------------------------------------------------

  !  updates the entropy field s and its radial derivatives
  !  adds explicit part to time derivatives of s

  !-------------------------------------------------------------------------

  USE truncation
  USE radial_functions,ONLY: i_costf_init,d_costf_init,orho1,or1,or2,n_r_cmb,n_r_icb,beta,drx,ddrx,&
       & kappa,dlkappa
  USE physical_parameters,ONLY: opr,polfac
  USE init_fields,ONLY: tops,bots
  USE blocking,ONLY: nLMBs,st_map,lo_map,lo_sub_map,lmStartB,lmStopB
  USE horizontal_data,ONLY: dLh,hdif_S
  USE logic,only: l_update_s
  USE matrices,ONLY: lSmat,s0Mat,s0Pivot,sMat,sPivot
  USE algebra, ONLY: cgeslML,sgesl
  USE LMLoop_data, ONLY: llm,ulm,llm_real,ulm_real
  USE parallel_mod,only: rank
  IMPLICIT NONE

  !-- Input of variables:
  REAL(kind=8),intent(IN) :: w1        ! weight for time step !
  REAL(kind=8),intent(IN) :: coex      ! factor depending on alpha
  REAL(kind=8),intent(IN) :: dt        ! time step
  INTEGER,intent(IN) :: nLMB

  !-- Input/output of scalar fields:
  COMPLEX(kind=8),intent(INOUT) :: s(llm:ulm,n_r_max)
  COMPLEX(kind=8),INTENT(OUT) :: ds(llm:ulm,n_r_max)
  COMPLEX(kind=8),intent(IN) :: dVSrLM(llm:ulm,n_r_max)
  COMPLEX(kind=8),intent(INOUT) :: dsdt(llm:ulm,n_r_max)
  COMPLEX(kind=8),intent(INOUT) :: dsdtLast(llm:ulm,n_r_max)
  !-- Output: udpated s,ds,dsdtLast

  !-- Input of recycled work arrays:
  !-- Local work arrays
  COMPLEX(kind=8) :: workA(llm:ulm,n_r_max)
  COMPLEX(kind=8) :: workB(llm:ulm,n_r_max)

  !-- Local variables:
  REAL(kind=8) :: w2            ! weight of second time step
  REAL(kind=8) :: O_dt
  INTEGER :: l1,m1              ! degree and order
  INTEGER :: lm1,lmB,lm         ! position of (l,m) in array
  INTEGER :: lmStart_real       ! range of lm for real array
  INTEGER :: lmStop_real        !
  INTEGER :: lmStart,lmStop
  INTEGER :: nLMB2
  INTEGER :: nR                 ! counts radial grid points
  INTEGER :: n_cheb             ! counts cheb modes
  REAL(kind=8) ::  rhs(n_r_max) ! real RHS for l=m=0
  COMPLEX(kind=8) :: rhs1(n_r_max,lo_sub_map%sizeLMB2max) ! complex RHS for l>0

  INTEGER, DIMENSION(:),POINTER :: nLMBs2,lm2l,lm2m
  INTEGER, DIMENSION(:,:),POINTER :: sizeLMB2,lm2
  INTEGER, DIMENSION(:,:,:),POINTER :: lm22lm,lm22l,lm22m

  !-- end of declaration
  !---------------------------------------------------------------------

  IF ( .NOT. l_update_s ) RETURN

  nLMBs2(1:nLMBs) => lo_sub_map%nLMBs2
  sizeLMB2(1:,1:) => lo_sub_map%sizeLMB2
  lm22lm(1:,1:,1:) => lo_sub_map%lm22lm
  lm22l(1:,1:,1:) => lo_sub_map%lm22l
  lm22m(1:,1:,1:) => lo_sub_map%lm22m
  lm2(0:,0:) => lo_map%lm2
  lm2l(1:lm_max) => lo_map%lm2l
  lm2m(1:lm_max) => lo_map%lm2m


  lmStart     =lmStartB(nLMB)
  lmStop      =lmStopB(nLMB)
  lmStart_real=2*lmStart-1
  lmStop_real =2*lmStop
  w2  =1.-w1
  O_dt=1.D0/dt


  !--- Finish calculation of dsdt:
  CALL get_drNS( dVSrLM,workA, &
       &         ulm_real-llm_real+1,lmStart_real-llm_real+1,lmStop_real-llm_real+1, &
       &         n_r_max,n_cheb_max,workB, &
       &         i_costf_init,d_costf_init,drx)

  DO nR=1,n_r_max
     DO lm=lmStart,lmStop
        dsdt(lm,nR)=orho1(nR)*(dsdt(lm,nR)-or2(nR)*workA(lm,nR))
     END DO
  END DO

  !DO l1=0,l_max
  !   DO m1=0,l1,minc
  !      WRITE(*,"(2I3,2ES20.12)") l1,m1,SUM( dsdt(lm2(l1,m1),:) )
  !   END DO
  !END DO
  !WRITE(*,"(A,2ES20.12)") "dsdt: ",SUM(dsdt)

  DO nLMB2=1,nLMBs2(nLMB)
     lmB=0
     DO lm=1,sizeLMB2(nLMB2,nLMB)
        lm1=lm22lm(lm,nLMB2,nLMB)
        l1 =lm22l(lm,nLMB2,nLMB)
        m1 =lm22m(lm,nLMB2,nLMB)
        IF ( l1 == 0 ) THEN
           IF ( .NOT. lSmat(l1) ) THEN
              CALL get_s0Mat(dt,s0Mat,s0Pivot)
              lSmat(l1)=.TRUE.
           END IF
           rhs(1)=      REAL(tops(0,0))
           rhs(n_r_max)=REAL(bots(0,0))
           DO nR=2,n_r_max-1
              rhs(nR)=REAL(s(lm1,nR))*O_dt+ &
                   w1*REAL(dsdt(lm1,nR)) + &
                   w2*REAL(dsdtLast(lm1,nR))
           END DO
           !WRITE(*,"(A,5I3,2ES22.12)") "l1==0, rhs=",nLMB2,lm,lm1,l1,m1,SUM(rhs),SUM(s0Mat)
           CALL sgesl(s0Mat,n_r_max,n_r_max,s0Pivot,rhs)
        ELSE
           IF ( .NOT. lSmat(l1) ) THEN
              CALL get_sMat(dt,l1,hdif_S(st_map%lm2(l1,m1)), &
                   sMat(1,1,l1),sPivot(1,l1))
              lSmat(l1)=.TRUE.
              !WRITE(*,"(A,I3,ES22.14)") "sMat: ",l1,SUM( sMat(:,:,l1) )
           END IF
           lmB=lmB+1
           rhs1(1,lmB)=      tops(l1,m1)
           rhs1(n_r_max,lmB)=bots(l1,m1)

           DO nR=2,n_r_max-1
              rhs1(nR,lmB)=s(lm1,nR)*O_dt + &
                   w1*dsdt(lm1,nR) + &
                   w2*dsdtLast(lm1,nR)
           END DO
        END IF
        !IF (lmB.GT.0) WRITE(*,"(A,3I4,8ES20.12)") "rhs1 = ",lm1,l1,m1,SUM( rhs1(:,lmB) ),&
        !     & SUM( s(lm1,:) ),SUM( dsdt(lm1,:) ),SUM( dsdtLast(lm1,:) )
     END DO
     IF ( lmB > 0 ) THEN
        CALL cgeslML(sMat(1,1,l1),n_r_max,n_r_max, &
             &       sPivot(1,l1),rhs1,n_r_max,lmB)
     END IF
     lmB=0
     DO lm=1,sizeLMB2(nLMB2,nLMB)
        lm1=lm22lm(lm,nLMB2,nLMB)
        l1 =lm22l(lm,nLMB2,nLMB)
        m1 =lm22m(lm,nLMB2,nLMB)
        IF ( l1 == 0 ) THEN
           DO n_cheb=1,n_cheb_max
              s(lm1,n_cheb)=rhs(n_cheb)
           END DO
        ELSE
           lmB=lmB+1
           IF ( m1 > 0 ) THEN
              DO n_cheb=1,n_cheb_max
                 s(lm1,n_cheb)=rhs1(n_cheb,lmB)
              END DO
           ELSE
              DO n_cheb=1,n_cheb_max
                 s(lm1,n_cheb)= CMPLX(REAL(rhs1(n_cheb,lmB)),0.D0,KIND=KIND(0d0))
              END DO
           END IF
        END IF
     END DO

  END DO     ! loop over lm blocks

  !WRITE(*,"(A,2ES22.12)") "s after = ",SUM(s)
  !-- set cheb modes > n_cheb_max to zero (dealiazing)
  DO n_cheb=n_cheb_max+1,n_r_max
     DO lm1=lmStart,lmStop
        s(lm1,n_cheb)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
     END DO
  END DO

  !-- Get radial derivatives of s: workA,dsdtLast used as work arrays
  CALL costf1(s, ulm_real-llm_real+1, lmStart_real-llm_real+1, lmStop_real-llm_real+1, &
       &      dsdtLast, i_costf_init, d_costf_init)
  CALL get_ddr(s, ds, workA, ulm_real-llm_real+1, lmStart_real-llm_real+1, lmStop_real-llm_real+1, &
       &       n_r_max, n_cheb_max, workB, dsdtLast, &
       &       i_costf_init,d_costf_init,drx,ddrx)

  !-- Calculate explicit time step part:
  DO nR=n_r_cmb+1,n_r_icb-1
     DO lm1=lmStart,lmStop
        dsdtLast(lm1,nR)=dsdt(lm1,nR) &
             & - coex*opr*hdif_S(st_map%lm2(lm2l(lm1),lm2m(lm1))) * kappa(nR) * &
             &   ( workA(lm1,nR) &
             &     + ( PolFac*beta(nR) + 2.D0*or1(nR) + dLkappa(nR) ) * ds(lm1,nR) &
             &     - dLh(st_map%lm2(lm2l(lm1),lm2m(lm1))) * or2(nR)   *  s(lm1,nR) &
             &   )
     END DO
  END DO

  !-- workA=dds not needed further after this point, used as work array later


  RETURN
end SUBROUTINE updateS

!-------------------------------------------------------------------------------
