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
  USE LMLoop_data,ONLY: llm,ulm

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
       &                       time_passed,time_norm,s,ds)
    !********************************************************************

    !    !------------ This is release 2 level 1  --------------!
    !    !------------ Created on 1/17/02  by JW. --------------!

    !--------------------------------------------------------------------

    USE usefull, ONLY: cc2real
    USE integration, ONLY: rInt

    IMPLICIT NONE

    !-- Direct input:
    INTEGER,INTENT(IN) :: n_time_ave
    LOGICAL,INTENT(IN) :: l_stop_time
    REAL(kind=8),INTENT(IN) :: time_passed
    REAL(kind=8),INTENT(IN) :: time_norm
    COMPLEX(kind=8),INTENT(IN) :: s(llm:ulm,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: ds(llm:ulm,n_r_max)


    !-- Local:
    CHARACTER(len=72) :: outFile
    INTEGER :: n_r,lm,l,m,lc
    REAL(kind=8) :: T_temp
    REAL(kind=8) :: dT_temp
    REAL(kind=8) :: surf_ICB
    REAL(kind=8) :: fac,facICB
    REAL(kind=8) :: dt_norm

    REAL(kind=8),DIMENSION(n_r_max,l_max+1) :: T_r_l,T_r_l_global
    REAL(kind=8),DIMENSION(l_max+1) :: T_l
    REAL(kind=8),DIMENSION(l_max+1) ::  T_ICB_l,  T_ICB_l_global
    REAL(kind=8),DIMENSION(l_max+1) :: dT_ICB_l, dT_ICB_l_global

    REAL(kind=8) :: y,t,comp(l_max+1)
    INTEGER :: nOut

    !-- end of declaration
    !---------------------------------------------------------------------

    T_ICB_l=0.0D0
    dT_ICB_l=0.0D0

    DO n_r=1,n_r_max
       DO l=1,l_max+1
          T_r_l(n_r,l)=0.D0
          comp(l) = 0.0D0
       END DO
       DO lm=llm,ulm
          l =lo_map%lm2l(lm)
          m =lo_map%lm2m(lm)
          lc=l+1

          T_temp=DSQRT(cc2real(s(lm,n_r),m))/or2(n_r)

          !local_sum = 0.0D0
          !c = 0.0D0          !A running compensation for lost low-order bits.
          !DO i=lb,ub
          !   y = arr_local(i) - c    !So far, so good: c is zero.
          !   t = local_sum + y       !Alas, sum is big, y small, so low-order digits of y are lost.
          !   c = (t - local_sum) - y !(t - sum) recovers the high-order part of y; subtracting y recovers -(low part of y)
          !   local_sum = t           !Algebraically, c should always be zero. Beware eagerly optimising compilers!
          !   !Next time around, the lost low part will be added to y in a fresh attempt.
          !END DO
#if 0
          y = T_temp - comp(lc)
          t = T_r_l(n_r,lc) + y
          comp(lc) = (t-T_r_l(n_r,lc)) - y
          T_r_l(n_r,lc) = t
#else
          T_r_l(n_r,lc) =T_r_l(n_r,lc) +  T_temp
#endif

          IF ( n_r.EQ.n_r_icb ) THEN
             dT_temp=DSQRT(cc2real(ds(lm,n_r),m))/or2(n_r)
             T_ICB_l(lc) =  T_ICB_l(lc) +  T_temp
             dT_ICB_l(lc)= dT_ICB_l(lc) + dT_temp
          END IF
       END DO    ! do loop over lms in block 
    END DO    ! radial grid points 

    ! Reduction over all ranks
    CALL MPI_Reduce(T_r_l,T_r_l_global,n_r_max*(l_max+1),&
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(T_ICB_l,T_ICB_l_global,l_max+1,&
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(dT_ICB_l,dT_ICB_l_global,l_max+1,&
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    IF (rank.EQ.0) THEN
       !-- Radial Integrals:
       surf_ICB=4.D0*pi*r_icb*r_icb
       fac      =1.D0/vol_oc
       facICB   =1.D0/surf_ICB
       DO l=1,l_max+1
          T_l(l)=fac*rInt(T_r_l_global(1,l),n_r_max,dr_fac, &
               &          i_costf_init,d_costf_init)
          T_ICB_l(l)=facICB*T_ICB_l_global(l)
          dT_ICB_l(l)=facICB*dT_ICB_l_global(l)
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
             T_ICB2_ave(l)=T_ICB2_ave(l) + time_passed*T_ICB_l(l)*T_ICB_l(l)
             dT_ICB2_ave(l)=dT_ICB2_ave(l) + time_passed*dT_ICB_l(l)*dT_ICB_l(l)
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
             T2_ave(l)    =DSQRT( dt_norm*T2_ave(l) - T_ave(l)**2 ) 
             T_ICB2_ave(l)=DSQRT( dt_norm*T_ICB2_ave(l) - T_ICB_ave(l)**2 )
             dT_ICB2_ave(l)=DSQRT( dt_norm*dT_ICB2_ave(l) - dT_ICB_ave(l)**2 )
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
    END IF
  END SUBROUTINE spectrumC_average

  !----------------------------------------------------------------------------
END MODULE spectrumC_average_mod
