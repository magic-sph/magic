!$Id$
!***********************************************************************
    SUBROUTINE legPrepG(nR,nBc,lDeriv,lRmsCalc,l_frame, &
                              lTOnext,lTOnext2,lTOcalc, &
                 dLhw,dLhdw,dLhz,vhG,vhC,dvhdrG,dvhdrC, &
                           dLhb,dLhj,bhG,bhC,cbhG,cbhC, &
     sR,dsR,preR,dpR,zAS,dzAS,ddzAS,bCMB,omegaIC,omegaMA)
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

    USE truncation
    USE radial_functions
    USE num_param
    USE torsional_oscillations
    USE Grenoble
    USE blocking
    USE horizontal_data
    USE logic
    USE fields,ONLY: s_Rloc,ds_Rloc,z_Rloc,dz_Rloc,p_Rloc,dp_Rloc,b_Rloc,db_Rloc,ddb_Rloc,aj_Rloc,dj_Rloc,&
         & w_Rloc,dw_Rloc,ddw_Rloc,omega_ic,omega_ma
    USE const

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
    COMPLEX(kind=8),intent(OUT) :: dLhw(lm_max),dLhdw(lm_max),dLhz(lm_max)
    COMPLEX(kind=8),intent(OUT) :: vhG(lm_max),vhC(lm_max)
    COMPLEX(kind=8),intent(OUT) :: dvhdrG(lm_max),dvhdrC(lm_max)
    COMPLEX(kind=8),intent(OUT) :: dLhb(lm_max),dLhj(lm_max)
    COMPLEX(kind=8),intent(OUT) :: bhG(lm_max),bhC(lm_max)
    COMPLEX(kind=8),intent(OUT) :: cbhG(lm_max),cbhC(lm_max)
    COMPLEX(kind=8),intent(OUT) :: sR(lm_max),dsR(lm_max)
    COMPLEX(kind=8),intent(OUT) :: preR(lm_max),dpR(lm_max)
    COMPLEX(kind=8),intent(OUT) :: bCMB(lm_max)
    REAL(kind=8),intent(OUT) :: zAS(l_max+1),dzAS(l_max+1),ddzAS(l_max+1)
    REAL(kind=8),intent(OUT) :: omegaIC,omegaMA
     
!-- local:
    INTEGER :: lm,l,m
    COMPLEX(kind=8) :: dbd

!-- end of declaration
!--------------------------------------------------------------------------

            
    IF ( nR == n_r_icb ) omegaIC=omega_ic
    IF ( nR == n_r_cmb ) omegaMA=omega_ma
    IF ( l_conv .OR. l_mag_kin ) THEN

        IF ( l_heat ) THEN
            DO lm=1,lm_max
               sR(lm) =s_Rloc(lm,nR)   ! used for plotting and Rms
               dsR(lm)=ds_Rloc(lm,nR)  ! used for plotting and Rms
            END DO
        END IF
        IF ( lTOnext .OR. lTOnext2 .OR. lTOCalc ) THEN
            DO lm=1,lm_max
                l=lm2l(lm)
                m=lm2m(lm)
                IF ( l <= l_max .AND. m == 0 ) THEN
                    zAS(l+1)  =REAL(z_Rloc(lm,nR))   ! used in TO
                    dzAS(l+1) =REAL(dz_Rloc(lm,nR))  ! used in TO (anelastic)
                    ddzAS(l+1)=ddzASL(l+1,nR)    ! used in TO
                END IF
            END DO
        END IF
        IF ( lRmsCalc ) THEN
            preR(1)=zero
            dpR(1)=zero
            DO lm=2,lm_max
                preR(lm)=p_Rloc(lm,nR)    ! used for Rms in get_td (anelastic)
                dpR(lm)=dp_Rloc(lm,nR)  ! used for Rms in get_td
            END DO
        END IF
        IF ( l_mag .AND. l_frame .AND. &
             l_movie_oc .AND. nR == n_r_cmb ) THEN
            bCMB(1)=zero ! used in s_store_movie_frame.f
            DO lm=2,lm_max
                bCMB(lm)=b_Rloc(lm,nR)  ! used for movie output of surface field
            END DO
        END IF

        IF ( nBc /= 2 ) THEN ! nBc=2 is flag for fixed boundary
            dLhw(1)=zero
            vhG(1) =zero
            vhC(1) =zero
            DO lm=2,lm_max
                dLhw(lm)=dLh(lm)*w_Rloc(lm,nR)
                vhG(lm) =dw_Rloc(lm,nR) - &
                     CMPLX(-AIMAG(z_Rloc(lm,nR)),REAL(z_Rloc(lm,nR)),KIND=KIND(0d0))
                vhC(lm) =dw_Rloc(lm,nR) + &
                     CMPLX(-AIMAG(z_Rloc(lm,nR)),REAL(z_Rloc(lm,nR)),KIND=KIND(0d0))
            END DO
        END IF

        IF ( lDeriv ) THEN
            dLhdw(1) =zero
            dLhz(1)  =zero
            dvhdrG(1)=zero
            dvhdrC(1)=zero
            DO lm=2,lm_max
                dLhz(lm)  =dLh(lm)*z_Rloc(lm,nR)
                dLhdw(lm) =dLh(lm)*dw_Rloc(lm,nR)
                dvhdrG(lm)=ddw_Rloc(lm,nR) - &
                     CMPLX(-AIMAG(dz_Rloc(lm,nR)),REAL(dz_Rloc(lm,nR)),KIND=KIND(0d0))
                dvhdrC(lm)=ddw_Rloc(lm,nR) + &
                     CMPLX(-AIMAG(dz_Rloc(lm,nR)),REAL(dz_Rloc(lm,nR)),KIND=KIND(0d0))
            END DO
        END IF

    END IF
     
    IF ( l_mag .OR. l_mag_LF ) THEN
         
       !PRINT*,"aj: ",SUM(ABS(aj(:,nR))),SUM(ABS(dLh))
       !PRINT*,"dj: ",SUM(ABS(dj(:,nR)))
        dLhb(1)=zero
        bhG(1) =zero
        bhC(1) =zero
        DO lm=2,lm_max
            dLhb(lm)=dLh(lm)*b_Rloc(lm,nR)
            bhG(lm) =db_Rloc(lm,nR) - &
                 CMPLX(-AIMAG(aj_Rloc(lm,nR)),REAL(aj_Rloc(lm,nR)),KIND=KIND(0d0))
            bhC(lm) =db_Rloc(lm,nR) + &
                 CMPLX(-AIMAG(aj_Rloc(lm,nR)),REAL(aj_Rloc(lm,nR)),KIND=KIND(0d0))
        END DO
        IF ( lGrenoble ) THEN ! Add dipole imposed by inner core
            lm=lm2(1,0)
            dLhb(lm)=dLhb(lm)+dLh(lm)*b0(nR)
            bhG(lm) =bhG(lm)+db0(nR)
            bhC(lm) =bhC(lm)+db0(nR)
        END IF
        IF ( lDeriv ) THEN
            dLhj(1)=zero
            cbhG(1)=zero
            cbhC(1)=zero
            DO lm=2,lm_max
                dLhj(lm)=dLh(lm)*aj_Rloc(lm,nR)
                dbd     =or2(nR)*dLhb(lm)-ddb_Rloc(lm,nR)
                cbhG(lm)=dj_Rloc(lm,nR)-CMPLX(-AIMAG(dbd),REAL(dbd),KIND=KIND(0d0))
                cbhC(lm)=dj_Rloc(lm,nR)+CMPLX(-AIMAG(dbd),REAL(dbd),KIND=KIND(0d0))
            END DO
            IF ( lGrenoble ) THEN ! Add dipole imposed by inner core
                lm=lm2(1,0)
                cbhG(lm)=cbhG(lm)+CMPLX(0.D0,ddb0(nR),KIND=KIND(0d0))
                cbhC(lm)=cbhC(lm)-CMPLX(0.D0,ddb0(nR),KIND=KIND(0d0))
            END IF
        END IF

    END IF   ! magnetic terms required ?
     

    RETURN
    end SUBROUTINE legPrepG

!------------------------------------------------------------------------
