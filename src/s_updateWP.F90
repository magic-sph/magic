!$Id$
!***********************************************************************
SUBROUTINE updateWP(w,dw,ddw,dwdt,dwdtLast, &
     &              p,dp,dpdt,dpdtLast,s, &
     &              w1,coex,dt,nLMB,lRmsNext,nTh)
  !***********************************************************************

  !-----------------------------------------------------------------------

  !  updates the poloidal velocity potential w, the pressure p,  and
  !  their derivatives
  !  adds explicit part to time derivatives of w and p

  !-------------------------------------------------------------------------

  USE truncation
  USE radial_data,ONLY: n_r_cmb,n_r_icb
  USE radial_functions, ONLY: drx,ddrx,dddrx,or1,or2,rho0,agrav,rgrav,&
       &i_costf_init,d_costf_init,&
       &visc,dlvisc,beta,dbeta
  USE physical_parameters
  USE num_param
  USE blocking,only: nLMBs,lo_sub_map,lo_map,st_map,st_sub_map,&
       &lmStartB,lmStopB
  USE horizontal_data
  USE logic
  USE matrices
  USE RMS
  USE algebra, ONLY: cgeslML
  USE LMLoop_data, ONLY:llm,ulm, llm_real,ulm_real
  USE communications, only: get_global_sum
  IMPLICIT NONE

  !-- Input/output of scalar fields:
  COMPLEX(kind=8),INTENT(IN) :: dwdt(llm:ulm,n_r_max)
  COMPLEX(kind=8),intent(IN) :: dpdt(llm:ulm,n_r_max)
  COMPLEX(kind=8),INTENT(IN) :: s(llm:ulm,n_r_max)

  COMPLEX(kind=8),INTENT(INOUT) :: w(llm:ulm,n_r_max)
  COMPLEX(kind=8),intent(INOUT) :: dw(llm:ulm,n_r_max)
  COMPLEX(kind=8),intent(OUT) :: ddw(llm:ulm,n_r_max)
  COMPLEX(kind=8),intent(INOUT) :: dwdtLast(llm:ulm,n_r_max)
  COMPLEX(kind=8),intent(INOUT) :: p(llm:ulm,n_r_max)
  COMPLEX(kind=8),intent(OUT) :: dp(llm:ulm,n_r_max)
  COMPLEX(kind=8),intent(INOUT) :: dpdtLast(llm:ulm,n_r_max)
  !-- Output: updated w,dw,ddw,p,dp,dwdtLast,dpdtLast

  !-- Input of other variables:
  REAL(kind=8),intent(IN) :: w1        ! weight for time step !
  REAL(kind=8),intent(IN) :: coex      ! factor depending on alpha
  REAL(kind=8),intent(IN) :: dt        ! time step
  INTEGER,intent(IN) :: nLMB     ! block number
  LOGICAL,intent(IN) :: lRmsNext
  INTEGER,intent(IN) :: nTh      ! thread number

  !-- Input of recycled work arrays:
  COMPLEX(kind=8) :: workA(llm:ulm,n_r_max)
  COMPLEX(kind=8) :: workB(llm:ulm,n_r_max)

  !-- Local variables:
  REAL(kind=8) :: w2                  ! weight of second time step
  REAL(kind=8) :: O_dt
  INTEGER :: l1,m1              ! degree and order
  INTEGER :: lm1,lm,lmB         ! position of (l,m) in array
  INTEGER :: lmStart,lmStop ! max and min number of orders m
  INTEGER :: lmStart_real      ! range of lm for real array
  INTEGER :: lmStop_real       !
  INTEGER :: lmStart_00        ! excluding l=0,m=0
  INTEGER :: nLMB2
  INTEGER :: nR                ! counts radial grid points
  INTEGER :: n_cheb             ! counts cheb modes

  !COMPLEX(kind=8) :: Dif(lm_max),Pre(lm_max),Buo(lm_max)
  COMPLEX(kind=8) :: Dif(llm:ulm),Pre(llm:ulm),Buo(llm:ulm)

  COMPLEX(kind=8) :: rhs1(2*n_r_max,lo_sub_map%sizeLMB2max)
  COMPLEX(kind=8) :: rhs1_sum,rhs2_sum

  INTEGER, DIMENSION(:),POINTER :: nLMBs2,lm2l,lm2m
  INTEGER, DIMENSION(:,:),POINTER :: sizeLMB2,lm2
  INTEGER, DIMENSION(:,:,:),POINTER :: lm22lm,lm22l,lm22m

  !-- end of declaration
  !---------------------------------------------------------------------

  IF ( .NOT. l_update_v ) RETURN

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
  lmStart_00  =MAX(2,lmStart)
  lmStart_real=2*lmStart_00-1
  lmStop_real =2*lmStop

  w2  =1.D0-w1
  O_dt=1.D0/dt

  DO nLMB2=1,nLMBs2(nLMB)
     lmB=0
     DO lm=1,sizeLMB2(nLMB2,nLMB)
        lm1=lm22lm(lm,nLMB2,nLMB)
        l1 =lm22l(lm,nLMB2,nLMB)
        m1 =lm22m(lm,nLMB2,nLMB)
        IF ( l1 > 0 ) THEN
           IF ( .NOT. lWPmat(l1) ) THEN
#ifdef WITH_PRECOND_WP
              CALL get_wpMat(dt,l1,hdif_V(st_map%lm2(l1,m1)), &
                   wpMat(1,1,l1),wpPivot(1,l1),wpMat_fac(1,1,l1))
#else
              CALL get_wpMat(dt,l1,hdif_V(st_map%lm2(l1,m1)), &
                   wpMat(1,1,l1),wpPivot(1,l1))
#endif
              lWPmat(l1)=.TRUE.
           END IF
           lmB=lmB+1
           rhs1(1,lmB)        =0.D0
           rhs1(n_r_max,lmB)  =0.D0
           rhs1(n_r_max+1,lmB)=0.D0
           rhs1(2*n_r_max,lmB)=0.D0
           DO nR=2,n_r_max-1
              rhs1(nR,lmB)=                         &
                   & O_dt*dLh(st_map%lm2(l1,m1))*or2(nR)*w(lm1,nR) + &
                   & rho0(nR)*agrav(nR)*s(lm1,nR) + &
                   & w1*dwdt(lm1,nR) + &
                   & w2*dwdtLast(lm1,nR)
              rhs1(nR+n_r_max,lmB)=                 &
                   -O_dt*dLh(st_map%lm2(l1,m1))*or2(nR)*dw(lm1,nR) + &
                   w1*dpdt(lm1,nR) + &
                   w2*dpdtLast(lm1,nR)
              !IF ((l1.EQ.60).OR.(l1.EQ.61)) THEN
              !   WRITE(*,"(3I3,3(A,I3),A,12ES20.12)") lm1,l1,m1,": rhs1(",nR,"/",&
              !        & nR+n_r_max,",",lmB,") = ",&
              !        & w(lm1,nR),s(lm1,nR),dwdt(lm1,nR),dwdtLast(lm1,nR),&
              !        & rhs1(nR,lmB),rhs1(nR+n_r_max,lmB)
              !   WRITE(*,"(6X,A,8ES20.12)") "w,s,dwdt,dwdtLast = ",&
              !        & w(lm1,nR),s(lm1,nR),dwdt(lm1,nR),dwdtLast(lm1,nR)
              !END IF
           END DO
        END IF
     END DO
     IF ( lmB > 0 ) THEN
        !IF (lm1.GE.244) THEN
        !rhs1_sum=SUM(rhs1(1:n_r_max,1:lmB))
        !rhs2_sum=SUM(rhs1(n_r_max+1:2*n_r_max,1:lmB))
           
        !WRITE(*,"(2I4,A,2(I4,ES22.15))") nLMB2,lmB,": rhs1 before ",&
        !     & EXPONENT(REAL(rhs1_sum)),FRACTION(REAL(rhs1_sum)),&
        !     & EXPONENT(AIMAG(rhs1_sum)),FRACTION(AIMAG(rhs1_sum))
        !WRITE(*,"(2I4,A,2(I4,ES22.15))") nLMB2,lmB,": rhs2 before ",&
        !     & EXPONENT(REAL(rhs2_sum)),FRACTION(REAL(rhs2_sum)),&
        !     & EXPONENT(AIMAG(rhs2_sum)),FRACTION(AIMAG(rhs2_sum))
        !END IF
        ! use the mat_fac(:,1) to scale the rhs
#ifdef WITH_PRECOND_WP
        DO lm=1,lmB
           DO nR=1,2*n_r_max
              rhs1(nR,lm)=rhs1(nR,lm)*wpMat_fac(nR,1,l1)
           END DO
        END DO
#endif
        CALL cgeslML(wpMat(1,1,l1),2*n_r_max,2*n_r_max,    &
             &       wpPivot(1,l1),rhs1,2*n_r_max,lmB)

#ifdef WITH_PRECOND_WP
        ! rescale the solution with mat_fac(:,2)
        DO lm=1,lmB
           DO nR=1,2*n_r_max
              rhs1(nR,lm)=rhs1(nR,lm)*wpMat_fac(nR,2,l1)
           END DO
        END DO
#endif
        !rhs1_sum=SUM(rhs1(1:n_r_max,1:lmB))
        !rhs2_sum=SUM(rhs1(n_r_max+1:2*n_r_max,1:lmB))
        !WRITE(*,"(2I4,A,2(I4,ES22.15))") nLMB2,lmB,": rhs1 after ",&
        !     & EXPONENT(REAL(rhs1_sum)),FRACTION(REAL(rhs1_sum)),&
        !     & EXPONENT(AIMAG(rhs1_sum)),FRACTION(AIMAG(rhs1_sum))
        !WRITE(*,"(2I4,A,2(I4,ES22.15))") nLMB2,lmB,": rhs2 after ",&
        !     & EXPONENT(REAL(rhs2_sum)),FRACTION(REAL(rhs2_sum)),&
        !     & EXPONENT(AIMAG(rhs2_sum)),FRACTION(AIMAG(rhs2_sum))
        !END IF
     END IF

     IF ( lRmsNext ) THEN ! Store old w
        DO nR=1,n_r_max
           DO lm=1,sizeLMB2(nLMB2,nLMB)
              lm1=lm22lm(lm,nLMB2,nLMB)
              workB(lm1,nR)=w(lm1,nR)
           END DO
        END DO
     END IF

     lmB=0
     DO lm=1,sizeLMB2(nLMB2,nLMB)
        lm1=lm22lm(lm,nLMB2,nLMB)
        l1 =lm22l(lm,nLMB2,nLMB)
        m1 =lm22m(lm,nLMB2,nLMB)
        IF ( l1 > 0 ) THEN
           lmB=lmB+1
           IF ( m1 > 0 ) then
              DO n_cheb=1,n_cheb_max
                 w(lm1,n_cheb)=rhs1(n_cheb,lmB)
                 p(lm1,n_cheb)=rhs1(n_r_max+n_cheb,lmB)
              END DO
           ELSE
              DO n_cheb=1,n_cheb_max
                 w(lm1,n_cheb)= &
                      CMPLX(REAL(rhs1(n_cheb,lmB)),0.D0,KIND=KIND(0d0))
                 p(lm1,n_cheb)= &
                      CMPLX(REAL(rhs1(n_r_max+n_cheb,lmB)),0.D0,KIND=KIND(0d0))
              END DO
           END IF
        END IF
     END DO

  END DO   ! end of loop over lm1 blocks

  !WRITE(*,"(A,I3,4ES22.12)") "w,p after: ",nLMB,get_global_SUM(w),get_global_SUM(p)

  !-- set cheb modes > n_cheb_max to zero (dealiazing)
  DO n_cheb=n_cheb_max+1,n_r_max
     DO lm1=lmStart_00,lmStop
        w(lm1,n_cheb)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        p(lm1,n_cheb)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
     END DO
  END DO

  !-- Transform to radial space and get radial derivatives
  !   using dwdtLast, dpdtLast as work arrays:
  CALL costf1( w, ulm_real-llm_real+1, lmStart_real-llm_real+1, lmStop_real-llm_real+1, &
       &       dwdtLast,i_costf_init,d_costf_init)
  CALL get_dddr( w, dw, ddw, workA, &
       &         ulm_real-llm_real+1, lmStart_real-llm_real+1, lmStop_real-llm_real+1, &
       &         n_r_max,n_cheb_max,dwdtLast,dpdtLast, &
       &         i_costf_init,d_costf_init,drx,ddrx,dddrx)
  CALL costf1( p, ulm_real-llm_real+1, lmStart_real-llm_real+1, lmStop_real-llm_real+1, &
       &       dwdtLast,i_costf_init,d_costf_init)
  !WRITE(*,"(A,I4,A,4I5)") "size(p,1)=",SIZE(p,1),", llm_real,ulm_real,lmStart_real,lmStop_real=",&
  !     &llm_real,ulm_real,lmStart_real,lmStop_real
  CALL get_dr( p, dp, &
       &       ulm_real-llm_real+1, lmStart_real-llm_real+1, lmStop_real-llm_real+1, &
       &       n_r_max,n_cheb_max,dwdtLast,dpdtLast, &
       &       i_costf_init,d_costf_init,drx)
  !WRITE(*,"(A,4ES22.10)") "after get_dr: p,dp = ",get_global_sum(p),get_global_sum(dp)

  !-- Calculate explicit time step part:
  IF ( ra /= 0.D0 ) THEN
     DO nR=n_r_cmb+1,n_r_icb-1
        DO lm1=lmStart_00,lmStop
           l1=lm2l(lm1)
           m1=lm2m(lm1)

           Dif(lm1) = hdif_V(st_map%lm2(l1,m1))*dLh(st_map%lm2(l1,m1))*or2(nR)*visc(nR) * &
                & ( ddw(lm1,nR) &
                &   +(2.d0*dLvisc(nR)-beta(nR)/3.d0)*dw(lm1,nR) &
                &   -( dLh(st_map%lm2(l1,m1))*or2(nR)&
                &      +4.d0/3.d0*( dbeta(nR)&
                &                   +dLvisc(nR)*beta(nR)&
                &                   +(3.d0*dLvisc(nR)+beta(nR))*or1(nR)&
                &                 )&
                &    ) * w(lm1,nR) &
                & )
           Pre(lm1) = -dp(lm1,nR)+beta(nR)*p(lm1,nR)
           Buo(lm1) = rho0(nR)*rgrav(nR)*s(lm1,nR)
           dwdtLast(lm1,nR)=dwdt(lm1,nR) - coex*(Pre(lm1)+Buo(lm1)+Dif(lm1))
           dpdtLast(lm1,nR)=&
                & dpdt(lm1,nR) &
                & - coex*( dLh(st_map%lm2(l1,m1))*or2(nR)*p(lm1,nR) &
                &          + hdif_V(st_map%lm2(l1,m1))*visc(nR)*dLh(st_map%lm2(l1,m1))*or2(nR) &
                &            * ( -workA(lm1,nR) &
                &                + (beta(nR)-dLvisc(nR))*ddw(lm1,nR)&
                &                + ( dLh(st_map%lm2(l1,m1))*or2(nR)&
                &                    + dLvisc(nR)*beta(nR) &
                &                    + dbeta(nR) &
                &                    +2.d0*(dLvisc(nR)+beta(nR))*or1(nR)&
                &                  ) * dw(lm1,nR) &
                &                - dLh(st_map%lm2(l1,m1))*or2(nR) &
                &                  * (2.d0*or1(nR)+2.d0/3.d0*beta(nR)+dLvisc(nR))&
                &                  * w(lm1,nR)    &
                &              )  &
                &        )
           IF ( lRmsNext ) THEN
              workB(lm1,nR)=O_dt*dLh(st_map%lm2(l1,m1))*or2(nR) * &
                   &        ( w(lm1,nR)-workB(lm1,nR) )
              IF ( l_RMStest ) workB(lm1,nR)=workB(lm1,nR)-Dif(lm1)
           END IF
        END DO
        IF ( lRmsNext ) THEN
           CALL hInt2Pol(Dif,llm,ulm,nR,lmStart_00,lmStop,DifPolLMr, &
                DifPol2hInt(nR,nTh),DifPolAs2hInt(nR,nTh),lo_map)
           !WRITE(*,"(A,I4,3ES22.14)") "upWP, work=",nR,SUM(workB(:,nR)),dtVPol2hInt(nR,nTh)
           CALL hInt2Pol(workB(llm,nR),llm,ulm,nR,lmStart_00,lmStop, &
                dtVPolLMr,dtVPol2hInt(nR,nTh),dtVPolAs2hInt(nR,nTh),lo_map)
           !WRITE(*,"(A,2I4,ES22.14)") "upWP: ",nR,nTh,dtVPol2hInt(nR,nTh)
        END IF
     END DO

  ELSE  ! no s-contribution !

     DO nR=n_r_cmb+1,n_r_icb-1
        DO lm1=lmStart_00,lmStop
           l1=lm2l(lm1)
           m1=lm2m(lm1)
           Dif(lm1)=     hdif_V(st_map%lm2(l1,m1))*dLh(st_map%lm2(l1,m1))*or2(nR)*visc(nR)* &
                (                   ddw(lm1,nR)        + &
                (2.d0*dLvisc(nR)-beta(nR)/3.d0)*dw(lm1,nR) - &
                (dLh(st_map%lm2(l1,m1))*or2(nR)+4.d0/3.d0*(dbeta(nR)           + &
                dLvisc(nR)*beta(nR)+ &
                (3.d0*dLvisc(nR)+beta(nR))*or1(nR)))              * &
                w(lm1,nR) )
           Pre(lm1)=-dp(lm1,nR)+beta(nR)*p(lm1,nR)
           dwdtLast(lm1,nR) = dwdt(lm1,nR) - coex*(Pre(lm1)+Dif(lm1))
           dpdtLast(lm1,nR)=        dpdt(lm1,nR) - coex*( &
                dLh(st_map%lm2(l1,m1))*or2(nR)*p(lm1,nR) + &
                hdif_V(st_map%lm2(l1,m1))*visc(nR)*dLh(st_map%lm2(l1,m1))*or2(nR) * ( &
                -workA(lm1,nR)          + &
                (beta(nR)-dLvisc(nR))*ddw(lm1,nR)          + &
                ( dLh(st_map%lm2(l1,m1))*or2(nR)+dLvisc(nR)*beta(nR)    + &
                dbeta(nR) &
                +2.d0*(dLvisc(nR)+beta(nR))*or1(nR))    * &
                dw(lm1,nR)          - &
                dLh(st_map%lm2(l1,m1))*or2(nR)    * &
                (2.d0*or1(nR)+2.d0/3.d0*beta(nR)+dLvisc(nR))* &
                w(lm1,nR)          )  )
           IF ( lRmsNext ) THEN
              workB(lm1,nR)=O_dt*dLh(st_map%lm2(l1,m1))*or2(nR) * &
                   ( w(lm1,nR)-workB(lm1,nR) )
              IF ( l_RMStest ) &
                   workB(lm1,nR)=workB(lm1,nR)-Dif(lm1)
           END IF
        END DO
        IF ( lRmsNext ) THEN
           CALL hInt2Pol(Dif,llm,ulm,nR,lmStart_00,lmStop,DifPolLMr, &
                DifPol2hInt(nR,nTh),DifPolAs2hInt(nR,nTh),lo_map)
           CALL hInt2Pol(workB(llm,nR),llm,ulm,nR,lmStart_00,lmStop, &
                dtVPolLMr, dtVPol2hInt(nR,nTh),dtVpolAs2hInt(nR,nTh),lo_map)
        END IF
     END DO


  END IF


  !  Note: workA=dddw not needed beyond this point!


  RETURN
end SUBROUTINE updateWP

!------------------------------------------------------------------------------
