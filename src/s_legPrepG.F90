!$Id$
!***********************************************************************
!SUBROUTINE legPrepG(nR,nBc,lDeriv,lRmsCalc,l_frame, &
!     &              lTOnext,lTOnext2,lTOcalc, &
!     &              dLhw,dLhdw,dLhz,vhG,vhC,dvhdrG,dvhdrC, &
!     &              dLhb,dLhj,bhG,bhC,cbhG,cbhC, &
!     &              sR,dsR,preR,dpR,zAS,dzAS,ddzAS,bCMB,omegaIC,omegaMA)
SUBROUTINE legPrepG(nR,nBc,lDeriv,lRmsCalc,l_frame, &
     &              lTOnext,lTOnext2,lTOcalc, leg_helper)
  !***********************************************************************

  !    !------------ This is release 3 level 1  --------------!
  !    !------------ Created on 12/08/05  by JW. --------------!

  !  +-------------------------------------------------------------------+
  !  |                                                                   |
  !  |  Purpose of this subroutine is to prepare Legendre transforms     |
  !  |  from (r,l,m) space to (r,theta,m) space by calculating           |
  !  |  auxiliary arrays dpdw,dpddw, ....... dLhj which contain          |
  !  |  combinations of harmonic coeffs multiplied with (l,m)-dependend  |
  !  |  factors as well as the radial dependence.                        |
  !  |    nBc  =0 standard inner radial grid point                       |
  !  |    nBc  =1 radial velocity zero, spatial derivs not needed        |
  !  |    nBc  =2 all velocity comp. zero, spatial derivs not needed     |
  !  |   lDeriv=.TRUE. field derivatives required                        |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+

  USE truncation,ONLY: lm_max,l_max
  USE radial_functions,ONLY: n_r_icb,n_r_cmb,or2
  USE torsional_oscillations,only:ddzASL
  USE Grenoble,ONLY: lGrenoble, b0,db0,ddb0
  USE blocking,ONLY: lm2l,lm2m,lm2
  USE horizontal_data,only: dLh
  USE logic,ONLY: l_conv,l_mag_kin,l_heat,l_mag,l_movie_oc,l_mag_LF,l_fluxProfs
  USE fields,ONLY: s_Rloc,ds_Rloc, z_Rloc,dz_Rloc, p_Rloc,dp_Rloc, &
       &           b_Rloc,db_Rloc,ddb_Rloc, aj_Rloc,dj_Rloc,&
       &           w_Rloc,dw_Rloc,ddw_Rloc, omega_ic,omega_ma
  USE const,only: zero
  USE leg_helper_mod, only: leg_helper_t
  IMPLICIT NONE

  !-- Input of variables
  INTEGER,intent(IN) :: nR              ! radial level
  INTEGER,intent(IN) :: nBc             ! boundary condition
  LOGICAL,intent(IN) :: lDeriv          ! get also field derivatives !
  LOGICAL,intent(IN) :: lRmsCalc        ! Rms force balance ?
  LOGICAL,intent(IN) :: l_frame         ! Calculate movie frame?
  LOGICAL,intent(IN) :: lTOnext         ! for TO output
  LOGICAL,intent(IN) :: lTOnext2
  LOGICAL,intent(IN) :: lTOcalc

  !-- Input of scalar fields in LM-distributed space
  !   These have to be collected and stored in the
  !   R-distributed output variables

  !-- Output variables, R-distributed:
  TYPE(leg_helper_t),INTENT(INOUT) :: leg_helper

  !-- local:
  INTEGER :: lm,l,m
  COMPLEX(kind=8) :: dbd

  !-- end of declaration
  !--------------------------------------------------------------------------


  IF ( nR == n_r_icb ) leg_helper%omegaIC=omega_ic
  IF ( nR == n_r_cmb ) leg_helper%omegaMA=omega_ma
  IF ( l_conv .OR. l_mag_kin ) THEN

     IF ( l_heat ) THEN
        DO lm=1,lm_max
           leg_helper%sR(lm) =s_Rloc(lm,nR)   ! used for plotting and Rms
           leg_helper%dsR(lm)=ds_Rloc(lm,nR)  ! used for plotting and Rms
        END DO
     END IF
     IF ( lTOnext .OR. lTOnext2 .OR. lTOCalc ) THEN
        DO lm=1,lm_max
           l=lm2l(lm)
           m=lm2m(lm)
           IF ( l <= l_max .AND. m == 0 ) THEN
              leg_helper%zAS(l+1)  =REAL(z_Rloc(lm,nR))   ! used in TO
              leg_helper%dzAS(l+1) =REAL(dz_Rloc(lm,nR))  ! used in TO (anelastic)
              leg_helper%ddzAS(l+1)=ddzASL(l+1,nR)    ! used in TO
           END IF
        END DO
     END IF
     IF ( lRmsCalc .OR. l_fluxProfs ) THEN
        leg_helper%preR(1)=zero
        leg_helper%dpR(1)=zero
        DO lm=2,lm_max
           leg_helper%preR(lm)=p_Rloc(lm,nR)    ! used for Rms in get_td (anelastic)
           leg_helper%dpR(lm)=dp_Rloc(lm,nR)  ! used for Rms in get_td
        END DO
     END IF
     IF ( l_mag .AND. l_frame .AND. &
          l_movie_oc .AND. nR == n_r_cmb ) THEN
        leg_helper%bCMB(1)=zero ! used in s_store_movie_frame.f
        DO lm=2,lm_max
           leg_helper%bCMB(lm)=b_Rloc(lm,nR)  ! used for movie output of surface field
        END DO
     END IF

     IF ( nBc /= 2 ) THEN ! nBc=2 is flag for fixed boundary
        leg_helper%dLhw(1)=zero
        leg_helper%vhG(1) =zero
        leg_helper%vhC(1) =zero
        DO lm=2,lm_max
           leg_helper%dLhw(lm)=dLh(lm)*w_Rloc(lm,nR)
           leg_helper%vhG(lm) =dw_Rloc(lm,nR) - &
                CMPLX(-AIMAG(z_Rloc(lm,nR)),REAL(z_Rloc(lm,nR)),KIND=KIND(0d0))
           leg_helper%vhC(lm) =dw_Rloc(lm,nR) + &
                CMPLX(-AIMAG(z_Rloc(lm,nR)),REAL(z_Rloc(lm,nR)),KIND=KIND(0d0))
        END DO
     END IF

     IF ( lDeriv ) THEN
        leg_helper%dLhdw(1) =zero
        leg_helper%dLhz(1)  =zero
        leg_helper%dvhdrG(1)=zero
        leg_helper%dvhdrC(1)=zero
        DO lm=2,lm_max
           leg_helper%dLhz(lm)  =dLh(lm)*z_Rloc(lm,nR)
           leg_helper%dLhdw(lm) =dLh(lm)*dw_Rloc(lm,nR)
           leg_helper%dvhdrG(lm)=ddw_Rloc(lm,nR) - &
                CMPLX(-AIMAG(dz_Rloc(lm,nR)),REAL(dz_Rloc(lm,nR)),KIND=KIND(0d0))
           leg_helper%dvhdrC(lm)=ddw_Rloc(lm,nR) + &
                CMPLX(-AIMAG(dz_Rloc(lm,nR)),REAL(dz_Rloc(lm,nR)),KIND=KIND(0d0))
        END DO
     END IF

  END IF

  IF ( l_mag .OR. l_mag_LF ) THEN

     !PRINT*,"aj: ",SUM(ABS(aj(:,nR))),SUM(ABS(dLh))
     !PRINT*,"dj: ",SUM(ABS(dj(:,nR)))
     leg_helper%dLhb(1)=zero
     leg_helper%bhG(1) =zero
     leg_helper%bhC(1) =zero
     DO lm=2,lm_max
        leg_helper%dLhb(lm)=dLh(lm)*b_Rloc(lm,nR)
        leg_helper%bhG(lm) =db_Rloc(lm,nR) - &
             CMPLX(-AIMAG(aj_Rloc(lm,nR)),REAL(aj_Rloc(lm,nR)),KIND=KIND(0d0))
        leg_helper%bhC(lm) =db_Rloc(lm,nR) + &
             CMPLX(-AIMAG(aj_Rloc(lm,nR)),REAL(aj_Rloc(lm,nR)),KIND=KIND(0d0))
     END DO
     IF ( lGrenoble ) THEN ! Add dipole imposed by inner core
        lm=lm2(1,0)
        leg_helper%dLhb(lm)=leg_helper%dLhb(lm)+dLh(lm)*b0(nR)
        leg_helper%bhG(lm) =leg_helper%bhG(lm)+db0(nR)
        leg_helper%bhC(lm) =leg_helper%bhC(lm)+db0(nR)
     END IF
     IF ( lDeriv ) THEN
        leg_helper%dLhj(1)=zero
        leg_helper%cbhG(1)=zero
        leg_helper%cbhC(1)=zero
        DO lm=2,lm_max
           leg_helper%dLhj(lm)=dLh(lm)*aj_Rloc(lm,nR)
           dbd     =or2(nR)*leg_helper%dLhb(lm)-ddb_Rloc(lm,nR)
           leg_helper%cbhG(lm)=dj_Rloc(lm,nR)-CMPLX(-AIMAG(dbd),REAL(dbd),KIND=KIND(0d0))
           leg_helper%cbhC(lm)=dj_Rloc(lm,nR)+CMPLX(-AIMAG(dbd),REAL(dbd),KIND=KIND(0d0))
        END DO
        IF ( lGrenoble ) THEN ! Add dipole imposed by inner core
           lm=lm2(1,0)
           leg_helper%cbhG(lm)=leg_helper%cbhG(lm)+CMPLX(0.D0,ddb0(nR),KIND=KIND(0d0))
           leg_helper%cbhC(lm)=leg_helper%cbhC(lm)-CMPLX(0.D0,ddb0(nR),KIND=KIND(0d0))
        END IF
     END IF

  END IF   ! magnetic terms required ?


  RETURN
end SUBROUTINE legPrepG

!------------------------------------------------------------------------
