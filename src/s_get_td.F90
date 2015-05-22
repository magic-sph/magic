!$Id$
#include "perflib_preproc.cpp"
!***********************************************************************
SUBROUTINE get_td(nR,nBc,lRmsCalc,dVSrLM,dVxBhLM, &
     &            dwdt,dzdt,dpdt,dsdt,dbdt,djdt, &
     &            nl_lm,leg_helper)
  !***********************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/17/02  by JW. --------------!

  !  +-------------------------------------------------------------------+
  !  |                                                                   |
  !  |  Purpose of this to calculate time derivatives                    |
  !  |  dwdt,dzdt,dpdt,dsdt,dbdt,djdt                                    |
  !  |  and auxiliary arrays dVSrLM and dVxBhLM                          |
  !  |  from non-linear terms in spectral form,                          |
  !  |  contained in flmw1-3,flms1-3, flmb1-3 (input)                    |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+
  USE,INTRINSIC :: iso_c_binding
  USE truncation
  USE radial_functions,ONLY: r,or2,or1,beta,rho0,rgrav,epscProf,or4,temp0
  USE physical_parameters,ONLY: CorFac,ra,epsc,ViscHeatFac,OhmLossFac,n_r_LCR
  !USE num_param
  USE blocking,ONLY: lm2l,lm2m,lm2lmP,lmP2lmPS,lmP2lmPA,lm2lmA,lm2lmS,&
       & st_map
  USE horizontal_data,ONLY: dLh,dTheta1S,dTheta1A,dPhi,dTheta2A,&
       & dTheta3A,dTheta4A,dPhi0,dTheta2S,dTheta3S,dTheta4S,    &
       & hdif_V,hdif_B
  USE logic
  USE RMS,ONLY: AdvPolLMr,AdvPolAs2hInt,AdvPol2hInt,              &
       & AdvTor2hInt,AdvTorAs2hInt,PreLMr,Pre2hInt,               &
       & PreAs2hInt,BuoLMr,Buo2hInt,BuoAs2hInt,CorPolLMr,         &
       & CorPol2hInt,CorPolAs2hInt,CorTor2hInt,CorTorAs2hInt,     &
       & LFPolLMr,LFPol2hInt,LFPolAs2hInt,LFTor2hInt,LFTorAs2hInt,&
       & GeoLMr, Geo2hInt,GeoAs2hInt,Mag2hInt,MagAs2hInt,MagLMr,  &
       & ArcLMr,Arc2hInt,ArcAs2hInt
  USE leg_helper_mod,only:leg_helper_t
  USE nonlinear_lm_mod,only: nonlinear_lm_t
  USE fields, ONLY: w_Rloc,dw_Rloc,z_Rloc
  !USE cutils
#ifdef WITH_LIKWID
#   include "likwid_f90.h"
#endif
  IMPLICIT NONE

  !-- Input of variables:
  INTEGER,intent(IN) :: nR
  INTEGER,intent(IN) :: nBc ! signifies boundary conditions
  LOGICAL,intent(IN) :: lRmsCalc

  !----- Nonlinear terms:
  TYPE(nonlinear_lm_t),intent(IN) :: nl_lm

  !----- Scalar fields in radial distribute space:
  TYPE(leg_helper_t),INTENT(IN) :: leg_helper

  !-- Output:
  COMPLEX(kind=8),intent(OUT) :: dwdt(lm_max),dzdt(lm_max)
  COMPLEX(kind=8),intent(OUT) :: dpdt(lm_max),dsdt(lm_max)
  COMPLEX(kind=8),intent(OUT) :: dbdt(lm_maxMag),djdt(lm_maxMag)
  COMPLEX(kind=8),intent(OUT) :: dVxBhLM(lm_maxMag)
  COMPLEX(kind=8),intent(OUT) :: dVSrLM(lm_max)

  !-- local:
  INTEGER :: l,m,lm,lmS,lmA,lmP,lmPS,lmPA
  COMPLEX(kind=8) :: CorPol(lm_max),CorTor(lm_max)
  COMPLEX(kind=8) :: AdvPol(lm_max),AdvTor(lm_max)
  COMPLEX(kind=8) :: LFPol(lm_max),LFTor(lm_max)
  !COMPLEX(kind=8) :: zR(lm_max)!,dwR(lm_max)!,wR(lm_max)
  COMPLEX(kind=8) :: AdvPol_loc,CorPol_loc,AdvTor_loc,CorTor_loc
  complex(kind=8) :: dsdt_loc

  !INTEGER,PARAMETER :: lm_chunksize=4
  INTEGER,parameter :: DOUBLE_COMPLEX_PER_CACHELINE=4

  !-- end of declaration
  !-------------------------------------------------------------------------

  !WRITE(*,"(I3,A,4ES20.12)") nR,": get_td start: ",SUM(nl_lm%AdvrLM),SUM(leg_helper%dLHz)

  !lm_chunksize=(((lm_max)/nThreads)/DOUBLE_COMPLEX_PER_CACHELINE) * DOUBLE_COMPLEX_PER_CACHELINE
  !lm_chunksize=4
  !WRITE(*,"(A,I4)") "Using a chunksize of ",lm_chunksize

  IF ( nBc > 0 ) THEN   ! boundary !
     !PERFON('td_bnd')
     IF ( l_mag_nl .OR. l_mag_kin ) THEN

        !----- Stress free boundary, only nl mag. term for poloidal field needed.
        !      Because the radial derivative will be taken, this will contribute to
        !      the other radial grid points.
        dVxBhLM(1)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        dVSrLM(1) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        DO lm=2,lm_max
           l   =lm2l(lm)
           m   =lm2m(lm)
           lmP =lm2lmP(lm)
           lmPS=lmP2lmPS(lmP)   ! l-1
           lmPA=lmP2lmPA(lmP)   ! l+1
           IF ( l > m ) THEN
              dVxBhLM(lm)=r(nR)*r(nR)* (                &
                   & dTheta1S(lm)*nl_lm%VxBtLM(lmPS) - &
                   & dTheta1A(lm)*nl_lm%VxBtLM(lmPA) + &
                   & dPhi(lm)*nl_lm%VxBpLM(lmP)  )
           ELSE IF ( l == m ) THEN ! (l-1) not allowed !
              dVxBhLM(lm)=r(nR)*r(nR)* (                &
                   - dTheta1A(lm)*nl_lm%VxBtLM(lmPA) + &
                   dPhi(lm)*nl_lm%VxBpLM(lmP)  )
           END IF
           dVSrLM(lm)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        END DO

     ELSE
        DO lm=1,lm_max
           IF ( l_mag ) dVxBhLM(lm)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           dVSrLM(lm) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        END DO
     END IF
     IF ( l_heat ) THEN
        DO lm=1,lm_max
           dVSrLM(lm)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        END DO
     END IF
     !PERFOFF
  ELSE

     IF ( l_conv ) THEN  ! Convection

        lm =1   ! This is l=0,m=0
        lmA=lm2lmA(lm)
        lmP=1
        lmPA=lmP2lmPA(lmP)
        IF ( l_conv_nl ) THEN
           AdvPol_loc=      or2(nR)*nl_lm%AdvrLM(lm)
           AdvTor_loc=-dTheta1A(lm)*nl_lm%AdvpLM(lmPA)
        ELSE
           AdvPol_loc=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           AdvTor_loc=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        END IF
        IF ( l_corr ) THEN
           CorPol_loc=2.D0*CorFac*or1(nR) * dTheta2A(lm)* z_Rloc(lmA,nR)
           CorTor_loc= 2.d0*CorFac*or2(nR) * ( &
                dTheta3A(lm)*dw_Rloc(lmA,nR) + &
                or1(nR)*dTheta4A(lm)* w_Rloc(lmA,nR) )
        ELSE
           CorPol_loc=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           CorTor_loc=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        END IF
        dwdt(lm)=AdvPol_loc+CorPol_loc
        dzdt(lm)=AdvTor_loc+CorTor_loc
        IF ( lRmsCalc .AND. l_mag_LF .AND. nR>n_r_LCR ) THEN
           LFPol(lm) =      or2(nR)*nl_lm%LFrLM(lm)
           LFTor(lm) =-dTheta1A(lm)*nl_lm%LFpLM(lmPA)
           AdvPol(lm)=AdvPol_loc-LFPol(lm)
           AdvTor(lm)=AdvTor(lm)-LFTor(lm)
           CorPol(lm) = CorPol_loc
        END IF

        !PERFON('td_cv1')
        !$OMP PARALLEL DO default(none) &
        !$OMP private(lm,l,m,lmS,lmA,lmP,lmPS,lmPA) &
        !$OMP private(AdvPol_loc,CorPol_loc,AdvTor_loc,CorTor_loc) &
        !$OMP shared(lm2l,lm2m,lm2lmS,lm2lmA,lm2lmP,lmP2lmPS,lmP2lmPA) &
        !$OMP shared(lm_max,l_corr,l_max,l_conv_nl,nl_lm,lRmsCalc,l_mag_LF) &
        !$OMP shared(CorPol,AdvPol,LFPol,CorTor,AdvTor,LFTor,z_Rloc,nR,w_Rloc,dw_Rloc) &
        !$OMP shared(CorFac,or1,or2,dPhi0,dPhi,dTheta2A,dTheta2S,n_r_LCR) &
        !$OMP shared(dTheta3A,dTheta4A,dTheta3S,dTheta4S,dTheta1S,dTheta1A) &
        !$OMP shared(dwdt,dzdt)
        DO lm=1,lm_max
           IF (lm.EQ.1) CYCLE
           l   =lm2l(lm)
           m   =lm2m(lm)
           lmS =lm2lmS(lm)
           lmA =lm2lmA(lm)
           lmP =lm2lmP(lm)
           lmPS=lmP2lmPS(lmP)
           lmPA=lmP2lmPA(lmP)

#if 0
           WRITE(*,"(A,I4,A)") "========== lm = ",lm," =========="
           CALL print_cache_info_integer("l"//C_NULL_CHAR,l)
           CALL print_cache_info_dcmplx("dPhi0"//C_NULL_CHAR,dPhi0(lm))
           CALL print_cache_info_dcmplx("dPhi"//C_NULL_CHAR,dPhi(lm))
           CALL print_cache_info_dreal("dTheta1A"//C_NULL_CHAR,dTheta1A(lm))
           CALL print_cache_info_dreal("dTheta1S"//C_NULL_CHAR,dTheta1S(lm))
           CALL print_cache_info_dreal("dTheta2A"//C_NULL_CHAR,dTheta2A(lm))
           CALL print_cache_info_dreal("dTheta2S"//C_NULL_CHAR,dTheta2S(lm))
           CALL print_cache_info_dreal("dTheta3A"//C_NULL_CHAR,dTheta3A(lm))
           CALL print_cache_info_dreal("dTheta3S"//C_NULL_CHAR,dTheta3S(lm))
           CALL print_cache_info_dreal("dTheta4A"//C_NULL_CHAR,dTheta4A(lm))
           CALL print_cache_info_dreal("dTheta4S"//C_NULL_CHAR,dTheta4S(lm))
           CALL print_cache_info_dcmplx("z_Rloc(lmA)"//C_NULL_CHAR,z_Rloc(lmA,nR))
           CALL print_cache_info_dcmplx("z_Rloc(lmS)"//C_NULL_CHAR,z_Rloc(lmS,nR))
           CALL print_cache_info_dcmplx("w_Rloc(lmA)"//C_NULL_CHAR,w_Rloc(lmA,nR))
           CALL print_cache_info_dcmplx("w_Rloc(lmS)"//C_NULL_CHAR,w_Rloc(lmS,nR))
           CALL print_cache_info_dcmplx("dw_Rloc(lmA)"//C_NULL_CHAR,dw_Rloc(lmA,nR))
           CALL print_cache_info_dcmplx("dw_Rloc(lmS)"//C_NULL_CHAR,dw_Rloc(lmS,nR))
           CALL print_cache_info_dcmplx("dw_Rloc"//C_NULL_CHAR,dw_Rloc(lm,nR))
           CALL print_cache_info_dcmplx("nl_lm%AdvrLM(lmP)"//C_NULL_CHAR,nl_lm%AdvrLM(lmP))
           CALL print_cache_info_dcmplx("nl_lm%AdvtLM(lmP)"//C_NULL_CHAR,nl_lm%AdvtLM(lmP))
           CALL print_cache_info_dcmplx("nl_lm%AdvpLM(lmPA)"//C_NULL_CHAR,nl_lm%AdvpLM(lmPA))
           CALL print_cache_info_dcmplx("nl_lm%AdvpLM(lmPS)"//C_NULL_CHAR,nl_lm%AdvpLM(lmPS))
           WRITE(*,"(A)") ""
#endif
           IF ( l_corr ) THEN
              IF ( l < l_max .AND. l > m ) THEN
                 CorPol_loc =2.D0*CorFac*or1(nR) * ( &
                      dPhi0(lm)*dw_Rloc(lm,nR) +  & ! phi-deriv of dw/dr
                      dTheta2A(lm)*z_Rloc(lmA,nR) -  & ! sin(theta) dtheta z
                      dTheta2S(lm)*z_Rloc(lmS,nR) )

                 CorTor_loc=2.d0*CorFac*or2(nR) * ( &
                      dPhi0(lm)*z_Rloc(lm,nR)   + &
                      dTheta3A(lm)*dw_Rloc(lmA,nR)  + &
                      or1(nR)*dTheta4A(lm)* w_Rloc(lmA,nR)  + &
                      dTheta3S(lm)*dw_Rloc(lmS,nR)  - &
                      or1(nR)*dTheta4S(lm)* w_Rloc(lmS,nR)  )
              ELSE IF ( l == l_max ) THEN
                 CorPol_loc= 2.D0*CorFac*or1(nR) * ( &
                      dPhi0(lm)*dw_Rloc(lm,nR)  )
                 CorTor_loc=2.d0*CorFac*or2(nR) * ( &
                      dPhi0(lm)*z_Rloc(lm,nR)   )
              ELSE IF ( l == m ) THEN
                 CorPol_loc =2.D0*CorFac*or1(nR) * ( &
                      dPhi0(lm)*dw_Rloc(lm,nR)  + &
                      dTheta2A(lm)*z_Rloc(lmA,nR) )
                 CorTor_loc=2.d0*CorFac*or2(nR) * ( &
                      dPhi0(lm)*z_Rloc(lm,nR)   + &
                      dTheta3A(lm)*dw_Rloc(lmA,nR)  + &
                      or1(nR)*dTheta4A(lm)* w_Rloc(lmA,nR)  )
              END IF
           ELSE
              CorPol_loc=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              CorTor_loc=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           END IF

           IF ( l_conv_nl ) THEN
              AdvPol_loc=or2(nR)*nl_lm%AdvrLM(lmP)

              IF ( l > m ) THEN
                 AdvTor_loc=   -dPhi(lm)*nl_lm%AdvtLM(lmP)  + &
                      dTheta1S(lm)*nl_lm%AdvpLM(lmPS) - &
                      dTheta1A(lm)*nl_lm%AdvpLM(lmPA)
              ELSE IF ( l == m ) THEN
                 AdvTor_loc=   -dPhi(lm)*nl_lm%AdvtLM(lmP)  - &
                      dTheta1A(lm)*nl_lm%AdvpLM(lmPA)
              END IF
           ELSE
              AdvPol_loc=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              AdvTor_loc=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           END IF

           ! Here false sharing is still possible
           dwdt(lm)=AdvPol_loc+CorPol_loc
           dzdt(lm)=CorTor_loc+AdvTor_loc
           ! until here

           IF ( lRmsCalc .AND. l_mag_LF .AND. nR>n_r_LCR ) THEN
              !------ When RMS values are required, the Lorentz force is treated
              !       separately:
              !--------- FLPol = ( curl(B)xB )_r
              LFPol(lm) =or2(nR)*nl_lm%LFrLM(lmP)
              AdvPol(lm)=AdvPol_loc-LFPol(lm)
              CorPol(lm)=CorPol_loc

              !------ When RMS values are required, the Lorentz force is treated
              !       separately:
              IF ( l > m ) THEN
                 !------- LFTor= 1/(E*Pm) * curl( curl(B) x B )_r
                 LFTor(lm) =   -dPhi(lm)*nl_lm%LFtLM(lmP)  + &
                      dTheta1S(lm)*nl_lm%LFpLM(lmPS) - &
                      dTheta1A(lm)*nl_lm%LFpLM(lmPA)
              ELSE IF ( l == m ) THEN
                 LFTor(lm) =   -dPhi(lm)*nl_lm%LFtLM(lmP)  - &
                      dTheta1A(lm)*nl_lm%LFpLM(lmPA)
              END IF
              AdvTor(lm)=AdvTor_loc-LFTor(lm)
              CorTor(lm)=CorTor_loc
           END IF

        END DO
        !$OMP END PARALLEL DO
        !PERFOFF

        IF ( lRmsCalc ) THEN

           IF ( l_conv_nl ) THEN
              CALL hInt2Pol(AdvPol,1,lm_max,nR,2,lm_max,AdvPolLMr, &
                   AdvPol2hInt(nR),AdvPolAs2hInt(nR),st_map)
              CALL hInt2Tor(AdvTor,1,lm_max,nR,2,lm_max, &
                   AdvTor2hInt(nR),AdvTorAs2hInt(nR),st_map)
           END IF
           ! rho* grad(p/rho) = grad(p) - beta*p
           CALL hInt2Pol(leg_helper%dpR-beta(nR)*leg_helper%preR,1,lm_max,nR,2,lm_max,PreLMr, &
                Pre2hInt(nR),PreAs2hInt(nR),st_map)
           IF ( ra /= 0.D0 ) &
                CALL hInt2Pol(leg_helper%sR,1,lm_max,nR,2,lm_max,BuoLMr, &
                Buo2hInt(nR),BuoAs2hInt(nR),st_map)
           IF ( l_corr ) THEN
              CALL hInt2Pol(CorPol,1,lm_max,nR,2,lm_max,CorPolLMr, &
                   CorPol2hInt(nR),CorPolAs2hInt(nR),st_map)
              CALL hInt2Tor(CorTor,1,lm_max,nR,2,lm_max, &
                   CorTor2hInt(nR),CorTorAs2hInt(nR),st_map)
           END IF
           IF ( l_mag_LF .AND. nR>n_r_LCR ) THEN
              CALL hInt2Pol(LFPol,1,lm_max,nR,2,lm_max,LFPolLMr, &
                   LFPol2hInt(nR),LFPolAs2hInt(nR),st_map)
              CALL hInt2Tor(LFTor,1,lm_max,nR,2,lm_max, &
                   LFTor2hInt(nR),LFTorAs2hInt(nR),st_map)
           END IF

           !----- Calculate balances: recycle AdvPol, AdvTor for this:
           DO lm=2,lm_max
              IF ( l_RMStest ) THEN
                 AdvPol(lm)=AdvPol(lm)+CorPol(lm)+LFPol(lm)- &
                      leg_helper%dpR(lm)+beta(nR)*leg_helper%preR(lm)+ &
                      rho0(nR)*rgrav(nR)*leg_helper%sR(lm)
                 AdvTor(lm)=AdvTor(lm)+CorTor(lm)+LFTor(lm)
                 CorPol(lm)=or2(nR)*leg_helper%dLhw(lm)
              ELSE
                 AdvPol(lm)=CorPol(lm)-leg_helper%dpR(lm)+beta(nR)*leg_helper%preR(lm)
                 AdvTor(lm)=AdvPol(lm)+LFPol(lm)
                 CorPol(lm)=AdvTor(lm)+rho0(nR)*rgrav(nR)*leg_helper%sR(lm)
              END IF
           END DO
           CALL hInt2Pol(AdvPol,1,lm_max,nR,2,lm_max,GeoLMr, &
                Geo2hInt(nR),GeoAs2hInt(nR),st_map)
           IF ( l_RMStest ) THEN
              CALL hInt2Tor(AdvTor,1,lm_max,nR,2,lm_max, &
                   Mag2hInt(nR),MagAs2hInt(nR),st_map)
           ELSE
              CALL hInt2Pol(AdvTor,1,lm_max,nR,2,lm_max,MagLMr, &
                   Mag2hInt(nR),MagAs2hInt(nR),st_map)
           END IF
           CALL hInt2Pol(CorPol,1,lm_max,nR,2,lm_max,ArcLMr, &
                Arc2hInt(nR),ArcAs2hInt(nR),st_map)

        END IF

        !PERFON('td_cv2')
        !$OMP PARALLEL default(none) &
        !$OMP private(lm,l,m,lmS,lmA,lmP,lmPS,lmPA) &
        !$OMP private(AdvPol_loc,CorPol_loc) &
        !$OMP shared(lm2l,lm2m,lm2lmS,lm2lmA,lm2lmP,lmP2lmpS,lmP2lmPA) &
        !$OMP shared(lm_max,l_max,nR,l_corr,l_conv_nl) &
        !$OMP shared(CorFac,or1,or2,dPhi0,dTheta3A,dTheta3S,dTheta1S,dTheta1A) &
        !$OMP shared(z_Rloc,nl_lm,dPhi,leg_helper,dw_Rloc) &
        !$OMP shared(CorPol,AdvPol,dpdt)
        !LIKWID_ON('td_cv2')
        !$OMP DO
        DO lm=2,lm_max
           l   =lm2l(lm)
           m   =lm2m(lm)
           lmS =lm2lmS(lm)
           lmA =lm2lmA(lm)
           lmP =lm2lmP(lm)
           lmPS=lmP2lmPS(lmP)
           lmPA=lmP2lmPA(lmP)

           !------ Recycle CorPol and AdvPol:
           IF ( l_corr ) THEN
              !PERFON('td_cv2c')
              IF ( l < l_max .AND. l > m ) THEN
                 CorPol_loc= 2.d0*CorFac*or2(nR) * &
                      & ( -dPhi0(lm) * ( dw_Rloc(lm,nR) &
                      &                  +or1(nR)*leg_helper%dLhw(lm) &
                      &                ) &
                      &   +dTheta3A(lm)*z_Rloc(lmA,nR) &
                      &   +dTheta3S(lm)*z_Rloc(lmS,nR) &
                      & )

              ELSE IF ( l == l_max ) THEN
                 CorPol_loc=  2.d0*CorFac*or2(nR) * &
                      & ( -dPhi0(lm) * ( dw_Rloc(lm,nR) + or1(nR)*leg_helper%dLhw(lm) ) )

              ELSE IF ( l == m ) THEN
                 CorPol_loc=  2.d0*CorFac*or2(nR) * &
                      & ( -dPhi0(lm) * ( dw_Rloc(lm,nR) &
                      &                  +or1(nR)*leg_helper%dLhw(lm) &
                      &                )&
                      &   +dTheta3A(lm)*z_Rloc(lmA,nR)  &
                      & )

              END IF
              !PERFOFF
           ELSE
              CorPol_loc=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           END IF
           IF ( l_conv_nl ) THEN
              !PERFON('td_cv2nl')
              IF ( l > m ) THEN
                 AdvPol_loc= dTheta1S(lm)*nl_lm%AdvtLM(lmPS) - &
                      dTheta1A(lm)*nl_lm%AdvtLM(lmPA) + &
                      dPhi(lm)*nl_lm%AdvpLM(lmP)
              ELSE IF ( l == m ) THEN
                 AdvPol_loc=-dTheta1A(lm)*nl_lm%AdvtLM(lmPA) + &
                      dPhi(lm)*nl_lm%AdvpLM(lmP)
              END IF
              !PERFOFF
           ELSE
              AdvPol_loc=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           END IF
           dpdt(lm)=AdvPol_loc+CorPol_loc

        END DO ! lm loop
        !$OMP END DO 
        !LIKWID_OFF('td_cv2')
        !$OMP END PARALLEL 
        !PERFOFF
     ELSE
        DO lm=2,lm_max
           dwdt(lm) =0.D0
           dzdt(lm) =0.D0
           dpdt(lm) =0.D0
        END DO
     END IF ! l_conv ?

     IF ( l_heat ) THEN
        dsdt_loc  =epsc*epscProf(nR)!+opr/epsS*divKtemp0(nR)
        dVSrLM(1)=nl_lm%VSrLM(1)
        IF ( l_anel ) THEN
           IF ( l_anelastic_liquid ) THEN
              IF ( l_mag_nl ) THEN
                 dsdt_loc=dsdt_loc+                                        &
                 &    ViscHeatFac*hdif_V(1)*temp0(nR)*nl_lm%ViscHeatLM(1)+ &
                 &     OhmLossFac*hdif_B(1)*temp0(nR)*nl_lm%OhmLossLM(1)
              ELSE
                 dsdt_loc=dsdt_loc+ &
                 &    ViscHeatFac*hdif_V(1)*temp0(nR)*nl_lm%ViscHeatLM(1)
              END IF
           ELSE
              IF ( l_mag_nl ) THEN
                 dsdt_loc=dsdt_loc+ViscHeatFac*hdif_V(1)*nl_lm%ViscHeatLM(1)+ &
                      OhmLossFac*hdif_B(1)*nl_lm%OhmLossLM(1)
              ELSE
                 dsdt_loc=dsdt_loc+ViscHeatFac*hdif_V(1)*nl_lm%ViscHeatLM(1)
              END IF
           END IF
        END IF
        dsdt(1)=dsdt_loc

        !PERFON('td_heat')
        !$OMP PARALLEL DEFAULT(none) &
        !$OMP private(lm,l,m,lmP,lmPS,lmPA,dsdt_loc) &
        !$OMP shared(lm2l,lm2m,lm2lmP,lmP2lmPS,lmP2lmPA) &
        !$OMP shared(lm_max,dsdt,dVSrLM,dTheta1S,nl_lm,dTheta1A,dPhi) &
        !$OMP shared(l_anel,l_anelastic_liquid,l_mag_nl,nR) &
        !$OMP shared(ViscHeatFac,hdif_V,OhmLossFac,hdif_B,temp0)
        !LIKWID_ON('td_heat')
        !$OMP DO
        DO lm=2,lm_max
           l   =lm2l(lm)
           m   =lm2m(lm)
           lmP =lm2lmP(lm)
           lmPS=lmP2lmPS(lmP)
           lmPA=lmP2lmPA(lmP)
           !------ This is horizontal heat advection:
           !PERFON('td_h1')

           IF ( l > m ) THEN
              dsdt_loc= -dTheta1S(lm)*nl_lm%VStLM(lmPS) &
                   &    +dTheta1A(lm)*nl_lm%VStLM(lmPA) &
                   &    -dPhi(lm)*nl_lm%VSpLM(lmP)
           ELSE IF ( l == m ) THEN
              dsdt_loc=  dTheta1A(lm)*nl_lm%VStLM(lmPA) &
                   &     -dPhi(lm)*nl_lm%VSpLM(lmP)
           END IF
           !PERFOFF
           !PERFON('td_h2')
           IF ( l_anel ) THEN
              IF ( l_anelastic_liquid ) THEN
                 dsdt_loc = dsdt_loc+ &
                      &     ViscHeatFac*hdif_V(lm)*temp0(nR)*nl_lm%ViscHeatLM(lmP)
                 IF ( l_mag_nl ) THEN
                    dsdt_loc = dsdt_loc+ &
                         &     OhmLossFac*hdif_B(lm)*temp0(nR)*nl_lm%OhmLossLM(lmP)
                 END IF
              ELSE
                 dsdt_loc = dsdt_loc+ &
                      &     ViscHeatFac*hdif_V(lm)*nl_lm%ViscHeatLM(lmP)
                 IF ( l_mag_nl ) THEN
                    dsdt_loc = dsdt_loc+ &
                         &     OhmLossFac*hdif_B(lm)*nl_lm%OhmLossLM(lmP)
                 END IF
              END IF
           END IF
           !PERFOFF
           !-----   simplified form for linear onset !
           !        not ds not saved in the current program form!
           !                 dsdt(lm)=
           !                    -dLh(lm)*w(lm,nR)*or2(nR)*dsR(1)
           dVSrLM(lm)=nl_lm%VSrLM(lmP)
           dsdt(lm) = dsdt_loc
        END DO
        !$OMP END DO
        !LIKWID_OFF('td_heat')
        !$OMP END PARALLEL
        !PERFOFF
     ELSE
        DO lm=2,lm_max
           dsdt(lm)  =0.D0
           dVSrLM(lm)=0.D0
        END DO
     END IF

     IF ( l_mag_nl .OR. l_mag_kin  ) THEN
        !PERFON('td_magnl')

        !$OMP PARALLEL DO default(none) &
        !$OMP private(lm,l,m,lmP,lmPS,lmPA) &
        !$OMP shared(lm_max,lm2l,lm2m,lm2lmP,lmP2lmPS,lmP2lmPA) &
        !$OMP shared(nl_lm,dbdt,djdt,dTheta1S,dTheta1A,dPhi) &
        !$OMP shared(dLh,or4,dVxBhLM,r,nR)
        DO lm=1,lm_max
           IF (lm.EQ.1) THEN
              lmP=1
              lmPA=lmP2lmPA(lmP)
              dVxBhLM(lm)= -r(nR)*r(nR)* dTheta1A(lm)*nl_lm%VxBtLM(lmPA)
              dbdt(lm)   = -dTheta1A(lm)*nl_lm%VxBpLM(lmPA)
              djdt(lm)   = CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              CYCLE
           END IF
           l   =lm2l(lm)
           m   =lm2m(lm)
           lmP =lm2lmP(lm)
           lmPS=lmP2lmPS(lmP)
           lmPA=lmP2lmPA(lmP)

           !------- This is the radial part of the dynamo terms \curl(VxB)
           !PERFON('td_mnl1')
           IF ( l > m ) THEN
              dbdt(lm)=  dTheta1S(lm)*nl_lm%VxBpLM(lmPS) &
                   &    -dTheta1A(lm)*nl_lm%VxBpLM(lmPA) &
                   &    -dPhi(lm)    *nl_lm%VxBtLM(lmP)
           ELSE IF ( l == m ) THEN
              dbdt(lm)= -dTheta1A(lm)*nl_lm%VxBpLM(lmPA) &
                   &    -dPhi(lm)    *nl_lm%VxBtLM(lmP)
           END IF
           !PERFOFF

           !------- Radial component of
           !           \curl\curl(UxB) = \grad\div(UxB) - \laplace(VxB)

           !------- This is the radial part of \laplace (UxB)
           djdt(lm)=dLh(lm)*or4(nR)*nl_lm%VxBrLM(lmP)

           !------- This is r^2 * horizontal divergence of (UxB)
           !        Radial derivative performed in get_dr_td
           !PERFON('td_mnl2')
           IF ( l > m ) THEN
              dVxBhLM(lm)=          r(nR)*r(nR)* ( &
                   dTheta1S(lm)*nl_lm%VxBtLM(lmPS) - &
                   dTheta1A(lm)*nl_lm%VxBtLM(lmPA) + &
                   dPhi(lm)*nl_lm%VxBpLM(lmP)  )
           ELSE IF ( l == m ) THEN
              dVxBhLM(lm)=          r(nR)*r(nR)* ( &
                   - dTheta1A(lm)*nl_lm%VxBtLM(lmPA) + &
                   dPhi(lm)*nl_lm%VxBpLM(lmP)  )
           END IF
           !PERFOFF
        END DO
        !$OMP END PARALLEL DO
        !PERFOFF
     ELSE
        IF ( l_mag ) THEN
           DO lm=1,lm_max
              dbdt(lm)   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              djdt(lm)   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              dVxBhLM(lm)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           END DO
        END IF
     END IF

  END IF  ! boundary ? lvelo ?


  RETURN
end SUBROUTINE get_td
!-----------------------------------------------------------------------------
