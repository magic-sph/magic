!$Id$
MODULE spectrumC_average_mod
  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE blocking
  USE horizontal_data
  USE logic
  USE output_data
  USE const

  IMPLICIT NONE

  !-- output:
  REAL(kind=8),ALLOCATABLE :: T_ave(:)
  REAL(kind=8),ALLOCATABLE :: T_ICB_ave(:)
  REAL(kind=8),ALLOCATABLE :: dT_ICB_ave(:)
  REAL(kind=8),ALLOCATABLE :: T2_ave(:)
  REAL(kind=8),ALLOCATABLE :: T_ICB2_ave(:)
  REAL(kind=8),ALLOCATABLE :: dT_ICB2_ave(:)

contains
  SUBROUTINE initialize_spectrumC_average_mod

    ALLOCATE( T_ave(l_max+1) )
    ALLOCATE( T_ICB_ave(l_max+1) )
    ALLOCATE( dT_ICB_ave(l_max+1) )
    ALLOCATE( T2_ave(l_max+1) )
    ALLOCATE( T_ICB2_ave(l_max+1) )
    ALLOCATE( dT_ICB2_ave(l_max+1) )

  END SUBROUTINE initialize_spectrumC_average_mod

  !********************************************************************
  SUBROUTINE spectrumC_average(n_time_ave,l_stop_time,            &
       &                           time_passed,time_norm,s,ds)
    !********************************************************************

    !    !------------ This is release 2 level 1  --------------!
    !    !------------ Created on 1/17/02  by JW. --------------!

    !--------------------------------------------------------------------
    !
    !--------------------------------------------------------------------

    USE usefull, ONLY: cc2real
    USE integration, ONLY: rInt

    IMPLICIT NONE

    !-- Direct input:
    INTEGER :: n_time_ave
    LOGICAL :: l_stop_time
    REAL(kind=8) time_passed
    REAL(kind=8) time_norm
    COMPLEX(kind=8) :: s(lm_max,n_r_max)
    COMPLEX(kind=8) :: ds(lm_max,n_r_max)


    !-- Local:
    CHARACTER(len=72) :: outFile
    INTEGER :: n_r,lm,l,m,lc
    REAL(kind=8) :: T_temp
    REAL(kind=8) :: dT_temp
    REAL(kind=8) :: surf_ICB
    REAL(kind=8) :: fac,facICB
    REAL(kind=8) :: dt_norm

    REAL(kind=8) :: T_r_l(n_r_max,l_max+1)
    REAL(kind=8) :: T_l(l_max+1)
    REAL(kind=8) :: T_ICB_l(l_max+1)
    REAL(kind=8) :: dT_ICB_l(l_max+1)

    INTEGER :: nOut

    !-- end of declaration
    !---------------------------------------------------------------------

    T_ICB_l=0.0D0
    dT_ICB_l=0.0D0

    DO n_r=1,n_r_max
       DO l=1,l_max+1
          T_r_l(n_r,l)=0.D0
       END DO
       DO lm=1,lm_max
          l =lm2l(lm)
          m =lm2m(lm)
          lc=l+1

          T_temp=DSQRT(cc2real(s(lm,n_r),m))/or2(n_r)
          dT_temp=DSQRT(cc2real(ds(lm,n_r),m))/or2(n_r)
          T_r_l(n_r,lc) =T_r_l(n_r,lc) +  T_temp
          IF ( n_r.EQ.n_r_icb ) THEN
             T_ICB_l(lc) =T_ICB_l(lc) +T_temp
             dT_ICB_l(lc)=dT_ICB_l(lc) +dT_temp
          END IF
       END DO    ! do loop over lms in block 
    END DO    ! radial grid points 

    !-- Radial Integrals:
    surf_ICB=4.D0*pi*r_icb*r_icb
    fac      =1.D0/vol_oc
    facICB   =1.D0/surf_ICB
    DO l=1,l_max+1
       T_l(l)=fac*rInt(T_r_l(1,l),n_r_max,dr_fac,                   &
            &                     i_costf_init,d_costf_init)
       T_ICB_l(l)=facICB*T_ICB_l(l)
       dT_ICB_l(l)=facICB*dT_ICB_l(l)
    END DO

    !-- Averaging:
    IF ( n_time_ave.EQ.1 ) THEN
       DO l=1,l_max+1
          T_ave(l)     =time_passed*T_l(l)
          T_ICB_ave(l) =time_passed*T_ICB_l(l)
          dT_ICB_ave(l) =time_passed*dT_ICB_l(l)
          T2_ave(l)    =time_passed*T_l(l)*T_l(l)
          T_ICB2_ave(l)=time_passed*T_ICB_l(l)*T_ICB_l(l)
          dT_ICB2_ave(l)=time_passed*dT_ICB_l(l)*dT_ICB_l(l)
       END DO
    ELSE
       DO l=1,l_max+1
          T_ave(l)     =T_ave(l)+time_passed*T_l(l)
          T_ICB_ave(l) =T_ICB_ave(l)+time_passed*T_ICB_l(l)
          dT_ICB_ave(l)=dT_ICB_ave(l)+time_passed*dT_ICB_l(l)
          T2_ave(l)    =T2_ave(l)+time_passed*T_l(l)*T_l(l)
          T_ICB2_ave(l)=T_ICB2_ave(l) +                             &
               &                  time_passed*T_ICB_l(l)*T_ICB_l(l)
          dT_ICB2_ave(l)=dT_ICB2_ave(l) +                           &
               &                  time_passed*dT_ICB_l(l)*dT_ICB_l(l)
       END DO
    END IF

    !-- Output:
    IF ( l_stop_time ) THEN

       !------ Normalize:
       dt_norm=1.d0/time_norm
       DO l=1,l_max+1
          T_ave(l)     =dt_norm*T_ave(l)
          T_ICB_ave(l) =dt_norm*T_ICB_ave(l)
          dT_ICB_ave(l)=dt_norm*dT_ICB_ave(l)
          T2_ave(l)    =DSQRT( dt_norm*T2_ave(l) -                  &
               &                             T_ave(l)*T_ave(l) )
          T_ICB2_ave(l)=DSQRT( dt_norm*T_ICB2_ave(l) -              &
               &                         T_ICB_ave(l)*T_ICB_ave(l) )
          dT_ICB2_ave(l)=DSQRT( dt_norm*dT_ICB2_ave(l) -            &
               &                         dT_ICB_ave(l)*dT_ICB_ave(l) )
       END DO

       !------ Output:
       outFile='specAveC.'//TAG
       nOut   =93
       OPEN(nOut,file=outFile,status='UNKNOWN')
       DO l=1,l_max+1
          WRITE(93,'(2X,1P,I4,6D12.4)') l,                          &
               &              T_ave(l),T2_ave(l),                                 &
               &              T_ICB_ave(l),T_ICB2_ave(l),                         &
               &              dT_ICB_ave(l),dT_ICB2_ave(l) 
       END DO
       CLOSE(nOut)

       CALL safeOpen(nLF,log_file)
       WRITE(nLF,'(/,                                               &
            &'' ! TIME AVERAGED T/C SPECTRA STORED IN FILE: '',a)') outFile
       WRITE(nLF,'(                                                 &
            &'' !              No. of averaged spectra: '',I5)') n_time_ave
       CALL safeClose(nLF)

    END IF

  END SUBROUTINE spectrumC_average

  !----------------------------------------------------------------------------
END MODULE spectrumC_average_mod
