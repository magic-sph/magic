!$Id$
MODULE power
  USE truncation
  IMPLICIT NONE

  REAL(kind=8),DIMENSION(:),ALLOCATABLE :: buoMeanR
  REAL(kind=8),DIMENSION(:),ALLOCATABLE :: curlU2MeanR
  REAL(kind=8),DIMENSION(:),ALLOCATABLE :: ohmDissR

CONTAINS
  
  SUBROUTINE initialize_output_power

    ALLOCATE( buoMeanR(n_r_max) )
    ALLOCATE( ohmDissR(n_r_max) )
    ALLOCATE( curlU2MeanR(n_r_max) )

    buoMeanR=0.d0
    ohmDissR=0.d0
    curlU2MeanR=0.d0

  end subroutine initialize_output_power

  !********************************************************************
  SUBROUTINE get_power(time,timePassed,timeNorm,l_stop_time, &
       &               omega_IC,omega_MA, &
       &               lorentz_torque_IC,lorentz_torque_MA, &
       &               w,ddw,z,dz,s,b,ddb,aj,dj, &
       &               db_ic,ddb_ic,aj_ic,dj_ic, &
       &               viscDiss,ohmDiss)
    !********************************************************************

    !------------ This is release 2 level 6  --------------!
    !------------ Created on 2/25/02  by JW. --------------!

    !--------------------------------------------------------------------

    !  This subroutine calculates power and dissipation of
    !  the core/mantle system.
    !  Energy input into the outer core is by buoyancy and
    !  possibly viscous accelarations at the boundaries if
    !  the rotation rates of inner core or mantle are prescribed
    !  and kept fixed.
    !  The losses are due to Ohmic and viscous dissipation.
    !  If inner core and mantel are allowed to change their
    !  rotation rates due to viscous forces this power is not
    !  lost from the system and has to be respected.

    !  The output is written into a file  power.TAG.
    !  This file contains 10 columns:
    !    column  1: time
    !    column  2: buoyancy
    !    column  3: power due to inner core rotation
    !    column  4: power due to mantle rotation
    !    column  5: viscous dissipation
    !    column  6: Ohmic dissipation in inner and outer core
    !    column  7: mantle power
    !    column  8: inner core power
    !    column  9: power gain (power-diffusion) of core/mantel system
    !    column 10: time integrated power gain

    !--------------------------------------------------------------------

    USE truncation
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE blocking
    USE horizontal_data
    USE logic
    USE movie_data
    USE output_data
    USE usefull, ONLY: cc2real,cc22real
    USE LMLoop_data,ONLY: llm,ulm,llmMag,ulmMag
    USE integration, ONLY: rInt_R,rIntIC
    IMPLICIT NONE

    !-- INPUT of variables:
    LOGICAL,INTENT(IN) :: l_stop_time
    REAL(kind=8),INTENT(IN) :: time,timePassed,timeNorm
    REAL(kind=8),INTENT(IN) :: omega_IC,omega_MA
    REAL(kind=8),INTENT(IN) :: lorentz_torque_IC,lorentz_torque_MA

    !-- Input of scalar fields:
    COMPLEX(kind=8),INTENT(IN) :: w(llm:ulm,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: ddw(llm:ulm,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: z(llm:ulm,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: dz(llm:ulm,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: s(llm:ulm,n_r_max)

    COMPLEX(kind=8),INTENT(IN) :: b(llmMag:ulmMag,n_r_maxMag)
    COMPLEX(kind=8),INTENT(IN) :: ddb(llmMag:ulmMag,n_r_maxMag)
    COMPLEX(kind=8),INTENT(IN) :: aj(llmMag:ulmMag,n_r_maxMag)
    COMPLEX(kind=8),INTENT(IN) :: dj(llmMag:ulmMag,n_r_maxMag)

    COMPLEX(kind=8),INTENT(IN) :: db_ic(llmMag:ulmMag,n_r_ic_maxMag)
    COMPLEX(kind=8),INTENT(IN) :: ddb_ic(llmMag:ulmMag,n_r_ic_maxMag)
    COMPLEX(kind=8),INTENT(IN) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
    COMPLEX(kind=8),INTENT(IN) :: dj_ic(llmMag:ulmMag,n_r_ic_maxMag)

    !-- OUTPUT:
    REAL(kind=8),INTENT(OUT) :: viscDiss,ohmDiss

    !-- local:
    INTEGER :: n_r,lm,l,m,l1m0

    REAL(kind=8) :: r_ratio
    REAL(kind=8) :: curlB2,curlU2,buoy,curlB2_IC
    REAL(kind=8),DIMENSION(n_r_max) :: curlB2_r,curlB2_r_global
    REAL(kind=8),DIMENSION(n_r_max) :: curlU2_r,curlU2_r_global
    REAL(kind=8),DIMENSION(n_r_max) :: buoy_r,buoy_r_global
    REAL(kind=8),DIMENSION(n_r_ic_max) :: curlB2_rIC,curlB2_rIC_global
    REAL(kind=8) :: viscous_torque_ic,viscous_torque_ma

    COMPLEX(kind=8) :: laplace,Bh

    CHARACTER(len=76) :: fileName
    CHARACTER(len=7) :: marker
    REAL(kind=8) :: z10ICB,z10CMB,drz10ICB,drz10CMB
    REAL(kind=8) :: powerIC,powerMA
    REAL(kind=8) :: powerDiff,powerDiffOld,powerDiffT,eDiffInt
    REAL(kind=8) :: tStart
    SAVE   powerDiff,eDiffInt,marker,tStart

    logical :: rank_has_l1m0
    INTEGER :: sr_tag, status(MPI_STATUS_SIZE)
    !-- end of declaration
    !---------------------------------------------------------------------

    IF ( marker /= 'started' ) THEN
       tStart   =time
       powerDiff=0.D0
       eDiffInt =0.D0
    END IF

    DO n_r=1,n_r_max

       IF ( l_conv ) THEN
          curlU2_r(n_r)=0.d0
          !DO lm=2,lm_max
          DO lm=MAX(2,llm),ulm
             l=lo_map%lm2l(lm)
             m=lo_map%lm2m(lm)
             laplace=dLh(st_map%lm2(l,m))*or2(n_r)*w(lm,n_r)-ddw(lm,n_r)
             curlU2_r(n_r)=     curlU2_r(n_r) + dLh(st_map%lm2(l,m)) * ( &
                  dLh(st_map%lm2(l,m))*or2(n_r)*cc2real(z(lm,n_r),m)  + &
                  cc2real(dz(lm,n_r),m) + &
                  cc2real(laplace,m)    )
          END DO
       END IF

       IF ( l_mag ) THEN
          curlB2_r(n_r)=0.d0
          DO lm=MAX(2,llm),ulm
             l=lo_map%lm2l(lm)
             m=lo_map%lm2m(lm)
             laplace=dLh(st_map%lm2(l,m))*or2(n_r)*b(lm,n_r)-ddb(lm,n_r)
             curlB2_r(n_r)=curlB2_r(n_r) + dLh(st_map%lm2(l,m))*lambda(n_r)*( &
                  dLh(st_map%lm2(l,m))*or2(n_r)*cc2real(aj(lm,n_r),m) + &
                  cc2real(dj(lm,n_r),m) + &
                  cc2real(laplace,m)    )
          END DO
       END IF

       IF ( l_heat ) THEN
          buoy_r(n_r)=0.d0
          DO lm=MAX(2,llm),ulm
             l=lo_map%lm2l(lm)
             m=lo_map%lm2m(lm)
             buoy_r(n_r)=buoy_r(n_r) + dLh(st_map%lm2(l,m)) * &
                  rgrav(n_r)*cc22real(w(lm,n_r),s(lm,n_r),m)
          END DO
          !WRITE(*,"(A,I4,2ES22.14)") "buoy_r = ",n_r,buoy_r(n_r),rgrav(n_r)
       END IF

    END DO    ! radial grid points

    IF (l_conv) CALL MPI_Reduce(curlU2_r,curlU2_r_global,n_r_max,&
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if (l_mag) CALL MPI_Reduce(curlB2_r,curlB2_r_global,n_r_max,&
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if (l_heat) CALL MPI_Reduce(buoy_r,buoy_r_global,n_r_max,&
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    IF (rank.EQ.0) THEN
       !-- Transform to cheb space:
       IF ( l_conv ) THEN
          curlU2MeanR=curlU2MeanR+timePassed*curlU2_r_global*eScale
          curlU2=rInt_R(curlU2_r_global,n_r_max,n_r_max,drx, &
               i_costf_init,d_costf_init)
          curlU2=eScale*curlU2
       ELSE
          curlU2=0.D0
       END IF
       IF ( l_mag )  THEN
          ohmDissR=ohmDissR+timePassed*curlB2_r_global*LFfac*opm*eScale
          curlB2=rInt_R(curlB2_r_global,n_r_max,n_r_max,drx, &
               &        i_costf_init,d_costf_init)
          curlB2=LFfac*opm*eScale*curlB2
       ELSE
          curlB2=0.D0
       END IF
       IF ( l_heat ) THEN
          buoMeanR=buoMeanR+timePassed*buoy_r_global*eScale
          buoy=rInt_R(buoy_r_global,n_r_max,n_r_max,drx, &
               i_costf_init,d_costf_init)
          buoy=eScale*buoy
       ELSE
          buoy=0.D0
       END IF
    END IF

    !-- Inner core:

    IF ( l_cond_ic .AND. l_mag ) THEN

       DO n_r=1,n_r_ic_max
          r_ratio=r_ic(n_r)/r_ic(1)
          curlB2_rIC(n_r)=0.d0
          DO lm=MAX(2,llm),ulm
             l=lo_map%lm2l(lm)
             m=lo_map%lm2m(lm)
             Bh=(l+1.D0)*O_r_ic(n_r)*aj_ic(lm,n_r)+dj_ic(lm,n_r)
             laplace=-ddb_ic(lm,n_r) - &
                  2.D0*(l+1.D0)*O_r_ic(n_r)*db_ic(lm,n_r)
             curlB2_rIC(n_r)=curlB2_rIC(n_r) +                   &
                  dLh(st_map%lm2(l,m))*r_ratio**(2*l+2) *  ( &
                  dLh(st_map%lm2(l,m))*O_r_ic2(n_r)*cc2real(aj_ic(lm,n_r),m) + &
                  cc2real(Bh,m) + &
                  cc2real(laplace,m)       )
          END DO

       END DO    ! radial grid points

       CALL MPI_Reduce(curlB2_rIC,curlB2_rIC_global,n_r_ic_max,&
            & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

       IF (rank.EQ.0) THEN
          curlB2_IC=rIntIC(curlB2_rIC_global,n_r_ic_max,dr_fac_ic, &
               i_costf1_ic_init,d_costf1_ic_init)
          curlB2_IC=LFfac*opm*eScale*curlB2_IC
       END IF
    ELSE
       IF (rank.EQ.0) THEN
          curlB2_IC=0.D0
       END IF
    END IF  ! conducting inner core ?

    !-- Add up and correct:
    !  A correction is neccesaary when using the integral over
    !  (curl U)**2 to calculate the dippisated energy and inner core
    !  or mantle are allowed to rotate.
    !  The only flow component here is l=1,m=0, i.e. lm=2
    !  If the rotation rates of inner core of mantel are kept
    !  fixed, the viscous dissipation may actually be a power source if
    !  the radial derivatives drz10 are negative (positive) at the
    !  ICB (CMB)! The energy transfere is described by the very
    !  correction terms.
    l1m0=lo_map%lm2(1,0)
    rank_has_l1m0=.FALSE. ! set default
    sr_tag=46378 !arbitray send-recv tag
    IF (lmStartB(rank+1).LE.l1m0 .AND. lmStopB(rank+1).GE.l1m0 ) THEN
       !IF (rank.NE.0) THEN
       !   PRINT*,"in get_power, l1m0 is not on rank 0!"
       !   stop
       !END IF
       IF ( l_rot_IC ) THEN
          z10ICB  =REAL(z(l1m0,n_r_ICB))
          drz10ICB=REAL(dz(l1m0,n_r_ICB))
       ELSE
          z10ICB  =0.D0
          drz10ICB=0.D0
       END IF
       IF ( l_rot_MA ) THEN
          z10CMB  =REAL(z(l1m0,n_r_CMB))
          drz10CMB=REAL(dz(l1m0,n_r_CMB))
       ELSE
          z10CMB  =0.D0
          drz10CMB=0.D0
       END IF

       IF (rank.NE.0) THEN
          ! send data to rank 0
          CALL MPI_Send(z10ICB,1,MPI_DOUBLE_COMPLEX,0,sr_tag,MPI_COMM_WORLD,ierr)
          CALL MPI_Send(drz10ICB,1,MPI_DOUBLE_COMPLEX,0,sr_tag+1,MPI_COMM_WORLD,ierr)
          CALL MPI_Send(z10CMB,1,MPI_DOUBLE_COMPLEX,0,sr_tag+2,MPI_COMM_WORLD,ierr)
          CALL MPI_Send(drz10CMB,1,MPI_DOUBLE_COMPLEX,0,sr_tag+3,MPI_COMM_WORLD,ierr)
       END IF
       rank_has_l1m0=.TRUE.
    END IF

    IF (rank.EQ.0) THEN
       IF (.NOT.rank_has_l1m0) THEN
          ! receive data from the source ranks
          CALL MPI_Recv(z10ICB,1,MPI_DOUBLE_COMPLEX,&
               & MPI_ANY_SOURCE,sr_tag,MPI_COMM_WORLD,status,ierr)
          CALL MPI_Recv(drz10ICB,1,MPI_DOUBLE_COMPLEX,&
               & MPI_ANY_SOURCE,sr_tag+1,MPI_COMM_WORLD,status,ierr)
          CALL MPI_Recv(z10CMB,1,MPI_DOUBLE_COMPLEX,&
               & MPI_ANY_SOURCE,sr_tag+2,MPI_COMM_WORLD,status,ierr)
          CALL MPI_Recv(drz10CMB,1,MPI_DOUBLE_COMPLEX,&
               & MPI_ANY_SOURCE,sr_tag+3,MPI_COMM_WORLD,status,ierr)
       END IF

       IF ( l_conv ) THEN
          viscDiss= -curlU2
          IF ( l_rot_IC ) viscDiss=viscDiss - 2.D0*z10ICB*drz10ICB
          IF ( l_rot_MA ) viscDiss=viscDiss + 2.D0*z10CMB*drz10CMB
       ELSE
          viscDiss=0.D0
       END IF

       !--- If the inner core or mantle rotation rate is allowed to change due
       !    to viscous drag, the power transfer has to be taken into account:

       !-- Calculating viscous torques:
       IF ( l_rot_ic .AND. kbotv == 2 ) THEN
          CALL get_viscous_torque(viscous_torque_ic, &
               &                  z10ICB,drz10ICB,r_icb)
       ELSE
          viscous_torque_ic=0.d0
       END IF
       IF ( l_rot_ma .AND. ktopv == 2 ) THEN
          CALL get_viscous_torque(viscous_torque_ma, &
               &                  z10CMB,drz10CMB,r_cmb)
       ELSE
          viscous_torque_ma=0.d0
       END IF

       IF ( l_rot_IC .AND. .NOT. l_SRIC ) THEN
          IF ( kbotv == 2 ) THEN
             CALL get_viscous_torque(viscous_torque_ic, &
                  &                  z10ICB,drz10ICB,r_icb)
          ELSE
             viscous_torque_ic=0.d0
          END IF
          powerIC=omega_IC*(viscous_torque_ic+lorentz_torque_ic)
       ELSE
          powerIC=0.D0
       END IF
       IF ( l_rot_MA ) THEN
          IF ( ktopv == 2 ) THEN
             CALL get_viscous_torque(viscous_torque_ma, &
                  &                  z10CMB,drz10CMB,r_cmb)
          ELSE
             viscous_torque_ma=0.d0
          END IF
          powerMA=omega_MA*(-viscous_torque_ma+lorentz_torque_ma)
       ELSE
          powerMA=0.D0
       END IF

       !--- Because the two systems are coupled only the total ohmic dissipation in usefull:
       ohmDiss=-curlB2-curlB2_IC

       powerDiffOld=powerDiff
       powerDiff   =(buoy+powerIC+powerMA+viscDiss+ohmDiss)

       IF ( marker == 'started' ) THEN
          powerDiffT  =1.5D0*powerDiff-0.5D0*powerDiffOld
          eDiffInt=eDiffInt+timePassed*timePassed*powerDiffT
          IF ( l_save_out ) THEN
             OPEN(n_power_file,FILE=power_file,status='unknown', &
                  POSITION='APPEND')
          END IF
          WRITE(n_power_file,'(1P,D20.12,9D16.8)') &
               time*tScale, &
               buoy, &
               -2.D0*z10ICB*drz10ICB, &
               2.D0*z10CMB*drz10CMB, &
               viscDiss, &
               ohmDiss, &
               powerMA, &
               powerIC, &
               powerDiff, &
               eDiffInt/timeNorm
          IF ( l_save_out ) CLOSE(n_power_file)
       ELSE
          marker='started'
       END IF

       IF ( l_stop_time ) THEN
          buoMeanR=buoMeanR/timeNorm
          ohmDissR=ohmDissR/timeNorm
          curlU2MeanR=curlU2MeanR/timeNorm
          fileName='powerR.'//tag
          OPEN(99,FILE=fileName,STATUS='unknown')
          DO n_r=1,n_r_max
             WRITE(99,'(4D20.10)')        &
                  &   r(n_r),             &! 1) radius
                  &   buoMeanR(n_r),      &! 2) Buo power
                  &   curlU2MeanR(n_r),   &! 3) Viscous heating
                  &   ohmDissR(n_r)        ! 4) Ohmic dissipation
          END DO
          CLOSE(99)
       END IF
    END IF
    RETURN

  END SUBROUTINE get_power
  !----------------------------------------------------------------------------
end module power
