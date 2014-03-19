!$Id$
!***********************************************************************
SUBROUTINE readStartFields(w,dwdt,z,dzdt,p,dpdt,s,dsdt, &
     &                     b,dbdt,aj,djdt, &
     &                     b_ic,dbdt_ic,aj_ic,djdt_ic, &
     &                     omega_ic,omega_ma, &
     &                     lorentz_torque_ic,lorentz_torque_ma, &
     &                     time,dt_old,dt_new,n_time_step)
  !***********************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/17/02  by JW. -----------

  !  +-------------+----------------+------------------------------------+
  !  |                                                                   |
  !  |   read initial condition from restart file                        |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+

  USE truncation
  USE radial_functions
  USE physical_parameters
  USE init_fields
  USE blocking
  USE logic
  USE output_data
  USE const

  IMPLICIT NONE

  !-- output:
  REAL(kind=8),INTENT(OUT) :: time,dt_old,dt_new
  INTEGER,INTENT(OUT) :: n_time_step

  !-- Output fields:
  COMPLEX(kind=8),DIMENSION(lm_max,n_r_max),INTENT(OUT) :: w,z,s,p,dwdt,dzdt,dsdt,dpdt
  COMPLEX(kind=8),DIMENSION(lm_maxMag,n_r_maxMag),INTENT(OUT) :: b,aj,dbdt,djdt
  COMPLEX(kind=8),DIMENSION(lm_maxMag,n_r_ic_maxMag),INTENT(OUT) :: b_ic,aj_ic,dbdt_ic,djdt_ic
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
  INTEGER :: informOld

  REAL(kind=8) :: omega_ic1Old,omegaOsz_ic1Old
  REAL(kind=8) :: omega_ic2Old,omegaOsz_ic2Old
  REAL(kind=8) :: omega_ma1Old,omegaOsz_ma1Old
  REAL(kind=8) :: omega_ma2Old,omegaOsz_ma2Old

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
     READ(n_start_file)                                        &
          time,dt_old,ra_old,pr_old,pm_old,ek_old,radratio_old, &
          informOld,n_r_max_old,n_theta_max_old,n_phi_tot_old, &
          minc_old,nalias_old,n_r_ic_max_old,sigma_ratio_old
     n_time_step=0
  ELSE IF ( inform == 0 ) THEN
     READ(n_start_file)                                        &
          time,dt_old,ra_old,pr_old,pm_old,ek_old,radratio_old, &
          n_time_step,n_r_max_old,n_theta_max_old,n_phi_tot_old, &
          minc_old,nalias_old
  ELSE IF ( inform == 1 ) THEN
     READ(n_start_file)                                        &
          time,dt_old,ra_old,pr_old,pm_old,ek_old,radratio_old, &
          n_time_step,n_r_max_old,n_theta_max_old,n_phi_tot_old, &
          minc_old
     nalias_old=nalias
  ELSE IF ( inform >= 2 ) THEN
     READ(n_start_file)                                        &
          time,dt_old,ra_old,pr_old,pm_old,ek_old,radratio_old, &
          n_time_step,n_r_max_old,n_theta_max_old,n_phi_tot_old, &
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
  IF ( prmag /= pm_old )                                       &
       WRITE(*,'(/,'' ! New mag Pr.number (old/new):'',2d16.6)') &
       pm_old,prmag
  IF ( radratio /= radratio_old )                                 &
       WRITE(*,'(/,'' ! New mag aspect ratio (old/new):'',2d16.6)') &
       radratio_old,radratio
  IF ( sigma_ratio /= sigma_ratio_old )                          &
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

  !-- Outer core fields:
  CALL mapData(n_r_max_old,l_max_old,minc_old,l_mag_old, &
       w,dwdt,z,dzdt,p,dpdt,s,dsdt,b,dbdt,aj,djdt)

  !-- Inner core fields:
  IF ( l_mag_old ) THEN
     IF ( inform >= 2 .AND. sigma_ratio_old /= 0.D0 ) THEN
        CALL mapDataIC(n_r_ic_max_old,l_max_old,minc_old, &
             b_ic,dbdt_ic,aj_ic,djdt_ic)
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
     READ(n_start_file,ERR=100,END=100) lorentz_torque_ic, &
          lorentz_torque_ma
  ELSE IF ( inform >= 4 .AND. inform <= 6 .AND. lMagMem == 1 )THEN
     READ(n_start_file,ERR=100,END=100) lorentz_torque_ic, &
          lorentz_torque_ma,omega_ic,omega_ma
  ELSE IF ( inform == 7 .OR. inform == 8 ) THEN
     READ(n_start_file,ERR=100,END=100) lorentz_torque_ic, &
          lorentz_torque_ma, &
          omega_ic1Old,omegaOsz_ic1Old,tOmega_ic1, &
          omega_ic2Old,omegaOsz_ic2Old,tOmega_ic2, &
          omega_ma1Old,omegaOsz_ma1Old,tOmega_ma1, &
          omega_ma2Old,omegaOsz_ma2Old,tOmega_ma2
  ELSE IF ( inform > 8 ) THEN
     READ(n_start_file,ERR=100,END=100) lorentz_torque_ic, &
          lorentz_torque_ma, &
          omega_ic1Old,omegaOsz_ic1Old,tOmega_ic1, &
          omega_ic2Old,omegaOsz_ic2Old,tOmega_ic2, &
          omega_ma1Old,omegaOsz_ma1Old,tOmega_ma1, &
          omega_ma2Old,omegaOsz_ma2Old,tOmega_ma2, &
          dt_new
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


  GOTO 300
100 WRITE(*,*) '! Could not read last line in input file!'
  WRITE(*,*) '! Data missing or wrong format!'
  WRITE(*,*) '! Change inform accordingly!'
  STOP
300 CONTINUE

  CLOSE(n_start_file)


  RETURN
end SUBROUTINE readStartFields

!-----------------------------------------------------------------------------
