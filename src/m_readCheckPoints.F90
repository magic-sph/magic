!$Id$
#include "intrinsic_sizes.h"
MODULE readCheckPoints

  IMPLICIT NONE

  PRIVATE

  INTEGER(kind=8) :: bytes_allocated=0

  PUBLIC :: readStartFields

CONTAINS

!***********************************************************************
  SUBROUTINE readStartFields(w,dwdt,z,dzdt,p,dpdt,s,dsdt, &
     &                       b,dbdt,aj,djdt,b_ic,dbdt_ic, &
     &                   aj_ic,djdt_ic,omega_ic,omega_ma, &
     &               lorentz_torque_ic,lorentz_torque_ma, &
     &                    time,dt_old,dt_new,n_time_step)
!***********************************************************************

  !  +-------------+----------------+------------------------------------+
  !  |                                                                   |
  !  |   read initial condition from restart file                        |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+

  USE truncation, ONLY: n_r_max,lm_max,n_r_maxMag,lm_maxMag,n_r_ic_max, &
                        n_r_ic_maxMag,nalias,n_phi_tot,l_max,m_max,    &
                        minc,lMagMem
  USE logic, ONLY: l_rot_ma,l_rot_ic,l_SRIC,l_SRMA,l_cond_ic,l_heat,l_mag, &
                   l_mag_LF
  USE blocking, ONLY: lm2,lmStartB,lmStopB,nLMBs,lm2l,lm2m
  USE init_fields, ONLY: start_file,n_start_file,inform,tOmega_ic1,tOmega_ic2, &
                         tOmega_ma1,tOmega_ma2,omega_ic1,omegaOsz_ic1,         &
                         omega_ic2,omegaOsz_ic2,omega_ma1,omegaOsz_ma1,        &
                         omega_ma2,omegaOsz_ma2,tShift_ic1,tShift_ic2,         &
                         tShift_ma1,tShift_ma2,tipdipole
  USE radial_functions, ONLY: n_r_icb,n_r_CMB,r
  USE physical_parameters, ONLY: ra,ek,pr,prmag,radratio,sigma_ratio,kbotv,ktopv
  USE const, ONLY: c_z10_omega_ic, c_z10_omega_ma, pi

  !-- output:
  REAL(kind=8),INTENT(OUT) :: time,dt_old,dt_new
  INTEGER,INTENT(OUT) :: n_time_step

  !-- Output fields:
  COMPLEX(kind=8),DIMENSION(lm_max,n_r_max),INTENT(OUT) :: w,z,s,p,dwdt,dzdt,dsdt,dpdt
  COMPLEX(kind=8),DIMENSION(lm_maxMag,n_r_maxMag),INTENT(OUT) :: b,aj,dbdt,djdt
  COMPLEX(kind=8),DIMENSION(lm_maxMag,n_r_ic_maxMag),INTENT(OUT) :: b_ic,aj_ic
  COMPLEX(kind=8),DIMENSION(lm_maxMag,n_r_ic_maxMag),INTENT(OUT) :: dbdt_ic,djdt_ic
  REAL(kind=8),INTENT(OUT) :: omega_ic,omega_ma
  REAL(kind=8),INTENT(OUT) :: lorentz_torque_ic,lorentz_torque_ma

  !-- Local:
  INTEGER :: minc_old,n_phi_tot_old,n_theta_max_old,nalias_old
  INTEGER :: l_max_old,n_r_max_old
  INTEGER :: n_r_ic_max_old
  REAL(kind=8) :: pr_old,ra_old,pm_old
  REAL(kind=8) :: ek_old,radratio_old
  REAL(kind=8) :: sigma_ratio_old
  INTEGER :: nLMB,lm,lmStart,lmStop,nR,l1m0
  LOGICAL :: l_mag_old
  LOGICAL :: startfile_does_exist
  LOGICAL :: lreadS
  INTEGER :: informOld,ioerr
  INTEGER :: n_r_maxL,n_r_ic_maxL,n_data_oldP,lm_max_old,n_dataL
  INTEGER,ALLOCATABLE :: lm2lmo(:)

  REAL(kind=8) :: fr
  REAL(kind=8) :: omega_ic1Old,omegaOsz_ic1Old
  REAL(kind=8) :: omega_ic2Old,omegaOsz_ic2Old
  REAL(kind=8) :: omega_ma1Old,omegaOsz_ma1Old
  REAL(kind=8) :: omega_ma2Old,omegaOsz_ma2Old

  COMPLEX(kind=8),ALLOCATABLE :: wo(:),zo(:),po(:),so(:)

  !-- end of declaration
  !----------------------------------------------------------------------------


  INQUIRE(file=start_file,exist=startfile_does_exist)

  IF ( startfile_does_exist ) THEN
     OPEN(n_start_file,FILE=start_file, &
          STATUS='OLD',FORM='UNFORMATTED')
  ELSE
     WRITE(*,*)
     WRITE(*,*) '! The restart file does not exist !'
     STOP
  END IF

  sigma_ratio_old=0.D0  ! assume non conducting inner core !
  IF ( inform == -1 ) THEN ! This is default !
     READ(n_start_file)                                         &
          time,dt_old,ra_old,pr_old,pm_old,ek_old,radratio_old, &
          informOld,n_r_max_old,n_theta_max_old,n_phi_tot_old,  &
          minc_old,nalias_old,n_r_ic_max_old,sigma_ratio_old
     n_time_step=0
  ELSE IF ( inform == 0 ) THEN
     READ(n_start_file)                                         &
          time,dt_old,ra_old,pr_old,pm_old,ek_old,radratio_old, &
          n_time_step,n_r_max_old,n_theta_max_old,n_phi_tot_old,&
          minc_old,nalias_old
  ELSE IF ( inform == 1 ) THEN
     READ(n_start_file)                                         &
          time,dt_old,ra_old,pr_old,pm_old,ek_old,radratio_old, &
          n_time_step,n_r_max_old,n_theta_max_old,n_phi_tot_old,&
          minc_old
     nalias_old=nalias
  ELSE IF ( inform >= 2 ) THEN
     READ(n_start_file)                                         &
          time,dt_old,ra_old,pr_old,pm_old,ek_old,radratio_old, &
          n_time_step,n_r_max_old,n_theta_max_old,n_phi_tot_old,&
          minc_old,nalias_old,n_r_ic_max_old,sigma_ratio_old
  END IF
  IF ( inform == -1 ) inform=informOld

  !---- Compare parameters:
  IF ( ra /= ra_old ) &
       WRITE(*,'(/,'' ! New Rayleigh number (old/new):'',2d16.6)') ra_old,ra
  IF ( ek /= ek_old ) &
       WRITE(*,'(/,'' ! New Ekman number (old/new):'',2d16.6)') ek_old,ek
  IF ( pr /= pr_old ) &
       WRITE(*,'(/,'' ! New Prandtl number (old/new):'',2d16.6)') pr_old,pr
  IF ( prmag /= pm_old )                                         &
       WRITE(*,'(/,'' ! New mag Pr.number (old/new):'',2d16.6)') &
       pm_old,prmag
  IF ( radratio /= radratio_old )                                   &
       WRITE(*,'(/,'' ! New mag aspect ratio (old/new):'',2d16.6)') &
       radratio_old,radratio
  IF ( sigma_ratio /= sigma_ratio_old )                            &
       WRITE(*,'(/,'' ! New mag cond. ratio (old/new):'',2d16.6)') &
       sigma_ratio_old,sigma_ratio

  l_max_old=nalias_old*n_phi_tot_old/60
  l_mag_old=.false.
  IF ( pm_old /= 0.D0 ) l_mag_old= .TRUE. 

  IF ( n_phi_tot_old /= n_phi_tot) &
       WRITE(*,*) '! New n_phi_tot (old,new):',n_phi_tot_old,n_phi_tot
  IF ( nalias_old /= nalias) &
       WRITE(*,*) '! New nalias (old,new)   :',nalias_old,nalias
  IF ( l_max_old /= l_max ) &
       WRITE(*,*) '! New l_max (old,new)    :',l_max_old,l_max


  IF ( inform==6 .OR. inform==7 .OR. inform==9 .OR. inform==11 ) THEN
     lreadS=.FALSE.
  ELSE
     lreadS=.TRUE.
  END IF

  ALLOCATE( lm2lmo(lm_max) )

  CALL getLm2lmO(n_r_max,n_r_max_old,l_max,l_max_old, &
                 m_max,minc,minc_old,inform,lm_max,   &
                 lm_max_old,n_data_oldP,lm2lmo)

  ! allocation of local arrays.
  ! if this becomes a performance bottleneck, one can make a module
  ! and allocate the array only once in the initialization
  !ALLOCATE( wo(n_dataL),zo(n_dataL),po(n_dataL),so(n_dataL) )
  ALLOCATE( wo(n_data_oldP),zo(n_data_oldP),po(n_data_oldP),so(n_data_oldP) )
  bytes_allocated = bytes_allocated + 4*n_data_oldP*SIZEOF_DOUBLE_COMPLEX
  ! end of allocation

  !PERFON('mD_rd')
  IF ( lreadS ) THEN
     !READ(n_start_file) (wo(i),i=1,n_data_oldP),                  &
     !     &             (zo(i),i=1,n_data_oldP),                  &
     !     &             (po(i),i=1,n_data_oldP),                  &
     !     &             (so(i),i=1,n_data_oldP)
     !WRITE(*,"(A,I10,A)") "Reading four fields, each with ",n_data_oldP," double complex entries."
     READ(n_start_file) wo, zo, po, so
  ELSE
     !READ(n_start_file) (wo(i),i=1,n_data_oldP),                  &
     !     &             (zo(i),i=1,n_data_oldP),                  &
     !     &             (po(i),i=1,n_data_oldP)
     READ(n_start_file) wo, zo, po
  END If
  !PERFOFF

  n_r_maxL = MAX(n_r_max,n_r_max_old)

  CALL mapDataHydro( wo,zo,po,so,n_data_oldP,lm2lmo,  &
                    n_r_max_old,lm_max_old,n_r_maxL,  &
                 .FALSE.,.FALSE.,.FALSE.,.FALSE.,w,z,p,s )

  !PERFON('mD_rd_dt')
  IF ( lreadS ) THEN
     !READ(n_start_file) (so(i),i=1,n_data_old),                   &
     !     &                        (wo(i),i=1,n_data_old),                   &
     !     &                        (zo(i),i=1,n_data_old),                   &
     !     &                        (po(i),i=1,n_data_old)
     READ(n_start_file) so,wo,zo,po
  ELSE
     !READ(n_start_file) (wo(i),i=1,n_data_old),                   &
     !     &                        (zo(i),i=1,n_data_old),                   &
     !     &                        (po(i),i=1,n_data_old)
     READ(n_start_file) wo,zo,po
  END IF
  !PERFOFF

  CALL mapDataHydro( wo,zo,po,so,n_data_oldP,lm2lmo,         &
                     n_r_max_old,lm_max_old,n_r_maxL,.TRUE., &
                     .TRUE.,.TRUE.,.TRUE.,dwdt,dzdt,dpdt,dsdt )

  IF ( l_mag_old ) THEN
     READ(n_start_file) so,wo,zo,po

     CALL mapDataMag( wo,zo,po,so,n_data_oldP,n_r_max,n_r_max_old, &
                         lm_max_old,n_r_maxL,lm2lmo,n_r_maxMag,    &
                         .FALSE.,aj,dbdt,djdt,b )
  ELSE
     WRITE(*,*) '! No magnetic data in input file!'
  END IF


  !-- If mapping towards reduced symmetry, add thermal perturbation in
  !   mode (l,m)=(minc,minc) if parameter tipdipole .ne. 0
  IF ( l_heat .AND.                                               &
       &       minc<minc_old .AND. tipdipole>0.D0 ) THEN
     DO nLMB=1,nLMBs ! Blocking of loop over all (l,m)
        lmStart=lmStartB(nLMB)
        lmStop =lmStopB(nLMB)
        lm=l_max+2
        IF ( lmStart<=lm .AND. lmStop>=lm ) THEN
           DO nR=1,n_r_max+1
              fr=dsin(pi*(r(nR)-r(n_r_max)))
              s(lm,nR)=tipdipole*fr
           END DO
        END IF
     END DO
  END IF

  !-- If starting from data file with longitudinal symmetry, add
  !   weak non-axisymmetric dipole component if tipdipole .ne. 0
  IF ( ( l_mag .OR. l_mag_LF )                                    &
       &       .AND. minc==1 .AND. minc_old/=1 .AND.                  &
       &       tipdipole>0.d0 .AND. l_mag_old ) THEN
     DO nLMB=1,nLMBs ! Blocking of loop over all (l,m)
        lmStart=lmStartB(nLMB)
        lmStop =lmStopB(nLMB)
        lm=l_max+2
        IF ( lmStart<=lm .AND. lmStop>=lm ) THEN
           DO nR=1,n_r_max+1
              b(lm,nR)=tipdipole
           END DO
        END IF
     END DO
  END IF

  ! deallocation of the local arrays
  DEALLOCATE( lm2lmo )
  DEALLOCATE( wo,zo,po,so )
  bytes_allocated = bytes_allocated - 4*n_data_oldP*SIZEOF_DOUBLE_COMPLEX

  !CALL mapData(n_r_max_old,l_max_old,minc_old,l_mag_old, &
  !     w,dwdt,z,dzdt,p,dpdt,s,dsdt,b,dbdt,aj,djdt)

  !-- Inner core fields:
  IF ( l_mag_old ) THEN
     IF ( inform >= 2 .AND. sigma_ratio_old /= 0.D0 ) THEN
        ALLOCATE( lm2lmo(lm_max) )
        CALL getLm2lmO(n_r_ic_max,n_r_ic_max_old,l_max,l_max_old, &
                 m_max,minc,minc_old,inform,lm_max,   &
                 lm_max_old,n_data_oldP,lm2lmo)

        n_r_ic_maxL = MAX(n_r_ic_max,n_r_ic_max_old)
        ALLOCATE( wo(n_data_oldP),zo(n_data_oldP),po(n_data_oldP), &
                  so(n_data_oldP) )

        READ(n_start_file) so,wo,zo,po
        CALL mapDataMag( wo,zo,po,so,n_data_oldP,n_r_ic_max,n_r_ic_max_old, &
                         lm_max_old,n_r_ic_maxL,lm2lmo,n_r_ic_maxMag,       &
                         .TRUE.,aj_ic,dbdt_ic,djdt_ic,b_ic )

        DEALLOCATE( lm2lmo )
        DEALLOCATE( wo,zo,po,so )
     ELSE IF ( l_cond_ic ) THEN
        !----- No inner core fields provided by start_file, we thus assume that
        !      simple the inner core field decays like r**(l+1) from
        !      the ICB to r=0:
        WRITE(*,'(/,'' ! USING POTENTIAL IC fields!'')')

        DO nLMB=1,nLMBs ! Blocking of loop over all (l,m)
           lmStart=lmStartB(nLMB)
           lmStop =lmStopB(nLMB)
           DO lm=lmStart,lmStop
              DO nR=1,n_r_ic_max
                 b_ic(lm,nR)   =b(lm,n_r_CMB)
                 aj_ic(lm,nR)  =aj(lm,n_r_CMB)
                 dbdt_ic(lm,nR)=dbdt(lm,n_r_CMB)
                 djdt_ic(lm,nR)=djdt(lm,n_r_CMB)
              END DO
           END DO
        END DO
     END IF
  END IF

  !-- Lorentz-torques:
  !   NOTE: If lMagMem=.false. the memory required to read
  !         magnetic field is not available. The code therefore
  !         cannot read lorentz torques and rotations that
  !         are stored after the magnetic fields.
  !         In this case I set the lorentz torques to zero and
  !         calculate the rotation from the speed at the
  !         boundaries in the case of no slip conditions.
  omega_ic1Old     =0.D0
  omegaOsz_ic1Old  =0.D0
  tOmega_ic1       =0.D0
  omega_ic2Old     =0.D0
  omegaOsz_ic2Old  =0.D0
  tOmega_ic2       =0.D0
  omega_ma1Old     =0.D0
  omegaOsz_ma1Old  =0.D0
  tOmega_ma1       =0.D0
  omega_ma2Old     =0.D0
  omegaOsz_ma2Old  =0.D0
  tOmega_ma2       =0.D0
  dt_new           =dt_old
  IF ( inform == 3 .AND. l_mag_old .AND. lMagMem == 1 ) THEN
     READ(n_start_file,IOSTAT=ioerr) lorentz_torque_ic, &
          lorentz_torque_ma
     IF( ioerr/=0 ) THEN
        WRITE(*,*) '! Could not read last line in input file!'
        WRITE(*,*) '! Data missing or wrong format!'
        WRITE(*,*) '! Change inform accordingly!'
        STOP
     END IF
  ELSE IF ( inform >= 4 .AND. inform <= 6 .AND. lMagMem == 1 )THEN
     READ(n_start_file,IOSTAT=ioerr) lorentz_torque_ic, &
          lorentz_torque_ma,omega_ic,omega_ma
     IF( ioerr/=0 ) THEN
        WRITE(*,*) '! Could not read last line in input file!'
        WRITE(*,*) '! Data missing or wrong format!'
        WRITE(*,*) '! Change inform accordingly!'
        STOP
     END IF
  ELSE IF ( inform == 7 .OR. inform == 8 ) THEN
     READ(n_start_file,IOSTAT=ioerr) lorentz_torque_ic, &
          lorentz_torque_ma, &
          omega_ic1Old,omegaOsz_ic1Old,tOmega_ic1, &
          omega_ic2Old,omegaOsz_ic2Old,tOmega_ic2, &
          omega_ma1Old,omegaOsz_ma1Old,tOmega_ma1, &
          omega_ma2Old,omegaOsz_ma2Old,tOmega_ma2
     IF( ioerr/=0 ) THEN
        WRITE(*,*) '! Could not read last line in input file!'
        WRITE(*,*) '! Data missing or wrong format!'
        WRITE(*,*) '! Change inform accordingly!'
        STOP
     END IF
  ELSE IF ( inform > 8 ) THEN
     READ(n_start_file,IOSTAT=ioerr) lorentz_torque_ic, &
          lorentz_torque_ma, &
          omega_ic1Old,omegaOsz_ic1Old,tOmega_ic1, &
          omega_ic2Old,omegaOsz_ic2Old,tOmega_ic2, &
          omega_ma1Old,omegaOsz_ma1Old,tOmega_ma1, &
          omega_ma2Old,omegaOsz_ma2Old,tOmega_ma2, &
          dt_new
     IF( ioerr/=0 ) THEN
        WRITE(*,*) '! Could not read last line in input file!'
        WRITE(*,*) '! Data missing or wrong format!'
        WRITE(*,*) '! Change inform accordingly!'
        STOP
     END IF
  ELSE
     !-- These could possibly be calcualted from the B-field
     lorentz_torque_ic=0.D0
     lorentz_torque_ma=0.D0
  END IF
  IF ( inform < 11 ) THEN
     lorentz_torque_ic=pm_old*lorentz_torque_ic
     lorentz_torque_ma=pm_old*lorentz_torque_ma
  END IF

  IF ( l_SRIC ) THEN
     IF ( omega_ic1Old /= omega_ic1 )                     &
          WRITE(*,*) '! New IC rotation rate 1 (old/new):', &
          omega_ic1Old,omega_ic1
     IF ( omegaOsz_ic1Old /= omegaOsz_ic1 )                    &
          WRITE(*,*) '! New IC rotation osz. rate 1 (old/new):', &
          omegaOsz_ic1Old,omegaOsz_ic1
     IF ( omega_ic2Old /= omega_ic2 )                     &
          WRITE(*,*) '! New IC rotation rate 2 (old/new):', &
          omega_ic2Old,omega_ic2
     IF ( omegaOsz_ic2Old /= omegaOsz_ic2 )                    &
          WRITE(*,*) '! New IC rotation osz. rate 2 (old/new):', &
          omegaOsz_ic2Old,omegaOsz_ic2
  END IF
  IF ( l_SRMA ) THEN
     IF ( omega_ma1Old /= omega_ma1 )                     &
          WRITE(*,*) '! New MA rotation rate 1 (old/new):', &
          omega_ma1Old,omega_ma1
     IF ( omegaOsz_ma1Old /= omegaOsz_ma1 )                    &
          WRITE(*,*) '! New MA rotation osz. rate 1 (old/new):', &
          omegaOsz_ma1Old,omegaOsz_ma1
     IF ( omega_ma2Old /= omega_ma2 )                     &
          WRITE(*,*) '! New MA rotation rate 2 (old/new):', &
          omega_ma2Old,omega_ma2
     IF ( omegaOsz_ma2Old /= omegaOsz_ma2 )                    &
          WRITE(*,*) '! New MA rotation osz. rate 2 (old/new):', &
          omegaOsz_ma2Old,omegaOsz_ma2
  END IF


  !----- Set IC and mantle rotation rates:
  !      Following cases are covered:
  !       1) Prescribed inner-core rotation omega_ic_pre
  !       2) Rotation has been read above ( inform.ge.4)
  !       3) Rotation calculated from flow field z(l=1,m=0)
  !       4) No rotation
  !       5) Flow driven by prescribed inner core rotation
  !       l_SRIC=.true. (spherical Couette case)
  l1m0=lm2(1,0)
  IF ( l_rot_ic ) THEN
     IF ( l_SRIC .OR. omega_ic1 /= 0.D0 ) THEN
        IF ( tShift_ic1 == 0.D0 ) tShift_ic1=tOmega_ic1-time
        IF ( tShift_ic2 == 0.D0 ) tShift_ic2=tOmega_ic2-time
        tOmega_ic1=time+tShift_ic1
        tOmega_ic2=time+tShift_ic2
        omega_ic=omega_ic1*DCOS(omegaOsz_ic1*tOmega_ic1) + &
             omega_ic2*DCOS(omegaOsz_ic2*tOmega_ic2)
        WRITE(*,*)
        WRITE(*,*) '! I use prescribed inner core rotation rate:'
        WRITE(*,*) '! omega_ic=',omega_ic
        IF ( kbotv == 2 ) &
             z(l1m0,n_r_icb)=CMPLX(omega_ic/c_z10_omega_ic,0d0,KIND=KIND(0d0))
     ELSE IF ( inform >= 7 ) THEN
        omega_ic=omega_ic1Old
     END IF
  ELSE
     omega_ic=0.D0
  END IF

  !----- Mantle rotation, same as for inner core (see above)
  !      exept the l_SRIC case.
  IF ( l_rot_ma ) THEN
     IF ( l_SRMA .OR. omega_ma1 /= 0.D0 ) THEN
        IF ( tShift_ma1 == 0.D0 ) tShift_ma1=tOmega_ma1-time
        IF ( tShift_ma2 == 0.D0 ) tShift_ma2=tOmega_ma2-time
        tOmega_ma1=time+tShift_ma1
        tOmega_ma2=time+tShift_ma2
        omega_ma=omega_ma1*DCOS(omegaOsz_ma1*tOmega_ma1) + &
             omega_ma2*DCOS(omegaOsz_ma2*tOmega_ma2)
        WRITE(*,*)
        WRITE(*,*) '! I use prescribed mantle rotation rate:'
        WRITE(*,*) '! omega_ma =',omega_ma
        WRITE(*,*) '! omega_ma1=',omega_ma1
        IF ( ktopv == 2 ) &
             z(l1m0,n_r_cmb)=CMPLX(omega_ma/c_z10_omega_ma,0d0,KIND=KIND(0d0))
     ELSE IF ( inform >= 7 ) THEN
        omega_ma=omega_ma1Old
     END IF
  ELSE
     omega_ma=0.D0
  END IF



  CLOSE(n_start_file)


  RETURN
  END SUBROUTINE readStartFields
!***********************************************************************
  SUBROUTINE getLm2lmO(n_r_max,n_r_max_old,l_max,l_max_old, &
                       m_max,minc,minc_old,inform,lm_max,   &
                       lm_max_old,n_data_oldP,lm2lmo)

  USE blocking, ONLY: lm2l,lm2m

  !--- Input variables
  INTEGER,INTENT(IN) :: n_r_max,l_max,m_max,minc
  INTEGER,INTENT(IN) :: n_r_max_old,l_max_old,minc_old
  INTEGER,INTENT(IN) :: inform,lm_max

  !--- Output variables
  INTEGER,DIMENSION(lm_max),INTENT(OUT) :: lm2lmo
  INTEGER,INTENT(OUT) :: n_data_oldP
  INTEGER,INTENT(OUT) :: lm_max_old

  !--- Local variables
  INTEGER :: n_data,n_data_old
  INTEGER :: m_max_old
  INTEGER :: l,m,lm,lmo,lo,mo

  !-- Outer core fields:
  n_data  = lm_max*n_r_max
  !-- This allows to increase the number of grid points by 10!

  IF (  l_max.EQ.l_max_old .AND.      &
       & minc.EQ.minc_old .AND.       &
       & n_r_max.EQ.n_r_max_old ) THEN

     !----- Direct reading of fields, grid not changed:
     WRITE(*,'(/,'' ! Reading fields directly.'')')

     n_data_old=n_data
     IF ( inform>2 ) THEN
        n_data_oldP=n_data
     ELSE
        !----- In the past an 'extra' radial grid point has been
        !      stored which was not really necessary
        n_data_oldP=lm_max*(n_r_max+1)
     END IF

     lm_max_old=lm_max
  ELSE

     !----- Mapping onto new grid !
     WRITE(*,'(/,'' ! Mapping onto new grid.'')')

     IF ( MOD(minc_old,minc).NE.0)                                &
          &     WRITE(6,'('' ! Warning: Incompatible old/new minc= '',2i3)')

     m_max_old =(l_max_old/minc_old)*minc_old
     lm_max_old=m_max_old*(l_max_old+1)/minc_old -                &
          &             m_max_old*(m_max_old-minc_old)/(2*minc_old) +        &
          &             l_max_old-m_max_old+1

     n_data_old=lm_max_old*n_r_max_old
     IF ( inform>2 ) THEN
        n_data_oldP=n_data_old
     ELSE
        n_data_oldP=lm_max_old*(n_r_max_old+1)
     END IF

     !-- Write info to STDOUT:
     WRITE(*,'('' ! Old/New  l_max= '',2I4,''  m_max= '',2I4,     &
          &       ''  minc= '',2I3,''  lm_max= '',2I5/)')                    &
          &           l_max_old,l_max,m_max_old,m_max,                       &
          &           minc_old,minc,lm_max_old,lm_max
     IF ( n_r_max_old.NE.n_r_max )                                &
          &        WRITE(*,'('' ! Old/New n_r_max='',2i4)')                  &
          &              n_r_max_old,n_r_max

  END IF

  DO lm=1,lm_max
     l=lm2l(lm)
     m=lm2m(lm)
     lm2lmo(lm)=-1 ! -1 means that there is no data in the startfile
     lmo=0
     DO mo=0,l_max_old,minc_old
        DO lo=mo,l_max_old
           lmo=lmo+1
           IF ( lo.EQ.l .AND. mo.EQ.m ) THEN
              lm2lmo(lm)=lmo ! data found in startfile
              CYCLE
           END IF
        END DO
     END DO
  END DO

  END SUBROUTINE getLm2lmO
!***********************************************************************
  SUBROUTINE mapDataHydro( wo,zo,po,so,n_data_oldP,lm2lmo, &
                          n_r_max_old,lm_max_old,n_r_maxL, &
                          lbc1,lbc2,lbc3,lbc4,w,z,p,s )

  USE logic, ONLY: l_heat
  USE init_fields, ONLY: scale_v,scale_s
  USE truncation, ONLY: n_r_max,lm_max
  USE blocking, ONLY: lmStartB,lmStopB,nLMBs

  !--- Input variables
  INTEGER,INTENT(IN) :: n_r_max_old,lm_max_old
  INTEGER,INTENT(IN) :: n_r_maxL,n_data_oldP
  LOGICAL,INTENT(IN) :: lbc1,lbc2,lbc3,lbc4
  INTEGER,INTENT(IN),DIMENSION(lm_max) :: lm2lmo
  COMPLEX(kind=8),INTENT(IN),DIMENSION(n_data_oldP) :: wo,zo,po,so

  !--- Output variables
  COMPLEX(kind=8),INTENT(OUT),DIMENSION(lm_max,n_r_max) :: w,z,p,s

  !--- Local variables
  INTEGER :: lm,lmo,n,nR,lmStart,lmStop,nLMB
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:) :: woR,zoR
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:) :: poR,soR

  !PRINT*,omp_get_thread_num(),": Before nLMB loop, nLMBs=",nLMBs
  ALLOCATE( woR(n_r_maxL),zoR(n_r_maxL) )
  ALLOCATE( poR(n_r_maxL),soR(n_r_maxL) )
  bytes_allocated = bytes_allocated + 4*n_r_maxL*SIZEOF_DOUBLE_COMPLEX
  WRITE(*,"(A,I12)") "maximal allocated bytes in mapData are ",bytes_allocated

  !PERFON('mD_map')
  DO nLMB=1,nLMBs ! Blocking of loop over all (l,m)
     lmStart=lmStartB(nLMB)
     lmStop =lmStopB(nLMB)

     !PRINT*,nLMB,lmStart,lmStop
     DO lm=lmStart,lmStop
        lmo=lm2lmo(lm)
        IF ( lmo>0 ) THEN
           IF ( n_r_max.NE.n_r_max_old ) THEN
              DO nR=1,n_r_max_old  ! copy on help arrays
                 n=lmo+(nR-1)*lm_max_old
                 woR(nR)=wo(n)
                 zoR(nR)=zo(n)
                 poR(nR)=po(n)
                 IF ( l_heat ) soR(nR)=so(n)
              END DO
              CALL mapDataR(woR,n_r_max,n_r_max_old,n_r_maxL,lBc1,.FALSE.)
              CALL mapDataR(zoR,n_r_max,n_r_max_old,n_r_maxL,lBc2,.FALSE.)
              CALL mapDataR(poR,n_r_max,n_r_max_old,n_r_maxL,lBc3,.FALSE.)
              IF ( l_heat ) CALL mapDataR(soR,n_r_max,n_r_max_old, & 
                                          n_r_maxL,lBc4,.FALSE.)
              DO nR=1,n_r_max
                 IF ( lm>1 ) THEN
                    w(lm,nR)=scale_v*woR(nR)
                    z(lm,nR)=scale_v*zoR(nR)
                    p(lm,nR)=scale_v*poR(nR)
                 END IF
                 IF ( l_heat ) s(lm,nR)=scale_s*soR(nR)
              END DO
           ELSE
              DO nR=1,n_r_max
                 n=lmo+(nR-1)*lm_max_old
                 IF ( lm>1 ) THEN
                    w(lm,nR)=scale_v*wo(n)
                    z(lm,nR)=scale_v*zo(n)
                    p(lm,nR)=scale_v*po(n)
                 END IF
                 IF ( l_heat ) s(lm,nR)=scale_s*so(n)
              END DO
           END IF
        END IF
     END DO
  END DO
  !PERFOFF
  !PRINT*,omp_get_thread_num(),": After nLMB loop"
  DEALLOCATE(woR,zoR,poR,soR)
  bytes_allocated = bytes_allocated - 4*n_r_maxL*SIZEOF_DOUBLE_COMPLEX

  END SUBROUTINE mapDataHydro
!***********************************************************************
  SUBROUTINE mapDataMag( wo,zo,po,so,n_data_oldP,n_rad_tot,n_r_max_old, &
                         lm_max_old,n_r_maxL,lm2lmo,dim1,l_IC,          & 
                         w,z,p,s )

  USE init_fields, ONLY: scale_b
  USE truncation, ONLY: lm_max,lm_maxMag
  USE blocking, ONLY: lmStartB,lmStopB,nLMBs

  !--- Input variables
  INTEGER,INTENT(IN) :: n_rad_tot,n_r_max_old,lm_max_old
  INTEGER,INTENT(IN) :: n_r_maxL,n_data_oldP,dim1
  INTEGER,INTENT(IN),DIMENSION(lm_max) :: lm2lmo
  LOGICAL,INTENT(IN) :: l_IC
  COMPLEX(kind=8),INTENT(IN),DIMENSION(n_data_oldP) :: wo,zo,po,so

  !--- Output variables
  COMPLEX(kind=8),INTENT(OUT),DIMENSION(lm_maxMag,dim1) :: w,z,p,s

  !--- Local variables
  INTEGER :: lm,lmo,n,nR,lmStart,lmStop,nLMB
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:) :: woR,zoR
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:) :: poR,soR

  !PRINT*,omp_get_thread_num(),": Before nLMB loop, nLMBs=",nLMBs
  ALLOCATE( woR(n_r_maxL),zoR(n_r_maxL) )
  ALLOCATE( poR(n_r_maxL),soR(n_r_maxL) )
  bytes_allocated = bytes_allocated + 4*n_r_maxL*SIZEOF_DOUBLE_COMPLEX
  WRITE(*,"(A,I12)") "maximal allocated bytes in mapData are ",bytes_allocated

  !PERFON('mD_map')
  DO nLMB=1,nLMBs ! Blocking of loop over all (l,m)
     lmStart=lmStartB(nLMB)
     lmStop =lmStopB(nLMB)
     lmStart=MAX(2,lmStart)
     DO lm=lmStart,lmStop
        lmo=lm2lmo(lm)
        IF ( lmo>0 ) THEN
           IF ( n_rad_tot/=n_r_max_old ) THEN
              DO nR=1,n_r_max_old  ! copy on help arrays
                 n=lmo+(nR-1)*lm_max_old
                 woR(nR)=wo(n)
                 zoR(nR)=zo(n)
                 poR(nR)=po(n)
                 soR(nR)=so(n)
              END DO
              CALL mapDataR(woR,dim1,n_r_max_old,n_r_maxL,.FALSE.,l_IC)
              CALL mapDataR(zoR,dim1,n_r_max_old,n_r_maxL,.TRUE.,l_IC)
              CALL mapDataR(poR,dim1,n_r_max_old,n_r_maxL,.TRUE.,l_IC)
              CALL mapDataR(soR,dim1,n_r_max_old,n_r_maxL,.FALSE.,l_IC)
              DO nR=1,n_rad_tot
                 w(lm,nR)=scale_b*woR(nR)
                 z(lm,nR)=scale_b*zoR(nR)
                 p(lm,nR)=scale_b*poR(nR)
                 s(lm,nR)=scale_b*soR(nR)
              END DO
           ELSE
              DO nR=1,n_rad_tot
                 n=lmo+(nR-1)*lm_max_old
                 w(lm,nR)=scale_b*wo(n)
                 z(lm,nR)=scale_b*zo(n)
                 p(lm,nR)=scale_b*po(n)
                 s(lm,nR)=scale_b*so(n)
              END DO
           END IF
        END IF
     END DO
  END DO
  !PERFOFF
  !PRINT*,omp_get_thread_num(),": After nLMB loop"
  DEALLOCATE(woR,zoR,poR,soR)
  bytes_allocated = bytes_allocated - 4*n_r_maxL*SIZEOF_DOUBLE_COMPLEX

  END SUBROUTINE mapDataMag
!***************************************************************************
  SUBROUTINE mapDataR(dataR,n_rad_tot,n_r_max_old,n_r_maxL,lBc,l_IC)
!***************************************************************************
  !---------------------------------------------------------------------------
  !
  !  Copy (interpolate) data (read from disc file) from old grid structure 
  !  to new grid. Linear interploation is used in r if the radial grid
  !  structure differs
  !
  !  called in mapdata
  !
  !---------------------------------------------------------------------------

  USE radial_functions, ONLY: i_costf_init,d_costf_init,cheb_norm, &
                              i_costf1_ic_init,d_costf1_ic_init,   &
                              cheb_norm_ic

  !--- Input variables
  INTEGER,INTENT(IN) :: n_r_max_old
  INTEGER,INTENT(IN) :: n_r_maxL,n_rad_tot
  LOGICAL,INTENT(IN) :: lBc,l_IC

  !--- Output variables
  COMPLEX(kind=8),INTENT(OUT),DIMENSION(:) :: dataR  ! old data 

  !-- Local variables
  INTEGER :: nR, n_r_index_start
  INTEGER,     ALLOCATABLE :: i_costf_init_old(:)
  REAL(kind=8),ALLOCATABLE :: d_costf_init_old(:)
  REAL(kind=8),ALLOCATABLE :: work(:)
  REAL(kind=8) :: cheb_norm_old,scale

  !-- end of declaration

  allocate( i_costf_init_old(2*n_r_maxL+2) )
  allocate( d_costf_init_old(2*n_r_maxL+5) )
  allocate( work(2*n_r_maxL) )

  !----- Initialize transform to cheb space:
  CALL init_costf1(n_r_max_old,i_costf_init_old,2*n_r_maxL+2,     &
       &                               d_costf_init_old,2*n_r_maxL+5)

  !-- Guess the boundary values, since they have not been stored:
  IF ( .NOT. l_IC .AND. lBc ) THEN
     dataR(1)=2.D0*dataR(2)-dataR(3)
     dataR(n_r_max_old)=2.D0*dataR(n_r_max_old-1) -               &
          &                             dataR(n_r_max_old-2)
  END IF

  !----- Transform old data to cheb space:
  !      Note: i_costf_init_old,d_costf_init_old used here!
  CALL costf1(dataR,2,1,2,work,i_costf_init_old,d_costf_init_old)

  !----- Fill up cheb polynomial with zeros:
  IF ( n_rad_tot>n_r_max_old ) THEN
     IF ( l_IC) THEN
        n_r_index_start=n_r_max_old
     ELSE 
        n_r_index_start=n_r_max_old+1
     END IF
     DO nR=n_r_index_start,n_rad_tot
        dataR(nR)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
     END DO
  END IF

  !----- Now transform to new radial grid points:
  !      Note: i_costf_init,d_costf_init used here!
  IF ( l_IC ) THEN
     CALL costf1(dataR,2,1,2,work,i_costf1_ic_init,d_costf1_ic_init)
     !----- Rescale :
     cheb_norm_old=DSQRT(2.D0/DBLE(n_r_max_old-1))
     scale=cheb_norm_old/cheb_norm_ic
  ELSE
     CALL costf1(dataR,2,1,2,work,i_costf_init,d_costf_init)
     !----- Rescale :
     cheb_norm_old=DSQRT(2.D0/DBLE(n_r_max_old-1))
     scale=cheb_norm_old/cheb_norm
  END IF
  DO nR=1,n_rad_tot
     dataR(nR)=scale*dataR(nR)
  END DO

  deallocate( i_costf_init_old )
  deallocate( d_costf_init_old )
  deallocate( work )

  END SUBROUTINE mapDataR
!---------------------------------------------------------------------
END MODULE readCheckPoints
