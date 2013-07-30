!$Id$
MODULE spectrum_average_mod
  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE blocking
  USE horizontal_data
  USE logic
  USE output_data
  USE LMLoop_data,ONLY: llm,ulm
  IMPLICIT NONE

  REAL(kind=8),ALLOCATABLE :: e_p_l_ave(:),e_p_m_ave(:)
  REAL(kind=8),ALLOCATABLE :: e_p2_l_ave(:),e_p2_m_ave(:)
  REAL(kind=8),ALLOCATABLE :: e_t_l_ave(:),e_t_m_ave(:)
  REAL(kind=8),ALLOCATABLE :: e_t2_l_ave(:),e_t2_m_ave(:)
  REAL(kind=8),ALLOCATABLE :: e_cmb_l_ave(:),e_cmb_m_ave(:)
  REAL(kind=8),ALLOCATABLE :: e_cmb2_l_ave(:),e_cmb2_m_ave(:)

  REAL(kind=8),ALLOCATABLE :: ek_p_l_ave(:),ek_p_m_ave(:)
  REAL(kind=8),ALLOCATABLE :: ek_p2_l_ave(:),ek_p2_m_ave(:)
  REAL(kind=8),ALLOCATABLE :: ek_t_l_ave(:),ek_t_m_ave(:)
  REAL(kind=8),ALLOCATABLE :: ek_t2_l_ave(:),ek_t2_m_ave(:)

CONTAINS
  SUBROUTINE initialize_spectrum_average_mod

    ALLOCATE( e_p_l_ave(0:l_max),e_p_m_ave(0:l_max) )
    ALLOCATE( e_p2_l_ave(0:l_max),e_p2_m_ave(0:l_max) )
    ALLOCATE( e_t_l_ave(0:l_max),e_t_m_ave(0:l_max) )
    ALLOCATE( e_t2_l_ave(0:l_max),e_t2_m_ave(0:l_max) )
    ALLOCATE( e_cmb_l_ave(0:l_max),e_cmb_m_ave(0:l_max) )
    ALLOCATE( e_cmb2_l_ave(0:l_max),e_cmb2_m_ave(0:l_max) )

    ALLOCATE( ek_p_l_ave(0:l_max),ek_p_m_ave(0:l_max) )
    ALLOCATE( ek_p2_l_ave(0:l_max),ek_p2_m_ave(0:l_max) )
    ALLOCATE( ek_t_l_ave(0:l_max),ek_t_m_ave(0:l_max) )
    ALLOCATE( ek_t2_l_ave(0:l_max),ek_t2_m_ave(0:l_max) )

  END SUBROUTINE initialize_spectrum_average_mod

  !********************************************************************
  SUBROUTINE spectrum_average(n_time_ave,l_stop_time,             &
       &                      time_passed,time_norm,b,aj,db,BV)
    !********************************************************************

    !    !------------ This is release 2 level 1  --------------!
    !    !------------ Created on 1/17/02  by JW. --------------!

    USE usefull, ONLY: cc2real
    USE integration, ONLY: rInt_R

    IMPLICIT NONE

    !-- Direct input:
    INTEGER,INTENT(IN) :: n_time_ave
    LOGICAL,INTENT(IN) :: l_stop_time
    REAL(kind=8),INTENT(IN) :: time_passed
    REAL(kind=8),INTENT(IN) :: time_norm

    COMPLEX(kind=8),INTENT(IN) :: b(llm:ulm,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: aj(llm:ulm,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: db(llm:ulm,n_r_max)
    CHARACTER(len=1),INTENT(IN) :: BV


    !-- output: 
    REAL(kind=8) :: e_p_l(0:l_max),e_t_l(0:l_max)
    REAL(kind=8) :: e_cmb_l(0:l_max)
    REAL(kind=8) :: e_p_m(0:l_max),e_t_m(0:l_max)
    REAL(kind=8) :: e_cmb_m(0:l_max)

    !-- local:
    CHARACTER(len=85) :: outFile
    INTEGER :: nOut
    INTEGER :: nR,lm,l,m

    REAL(kind=8) :: e_p_temp,e_t_temp
    REAL(kind=8) :: fac
    REAL(kind=8) :: dt_norm
    REAL(kind=8) :: SDp_l,SDp_m,SDt_l,SDt_m,SDcmb_l,SDcmb_m

    REAL(kind=8),DIMENSION(n_r_max,0:l_max) :: e_p_r_l,e_p_r_l_global
    REAL(kind=8),DIMENSION(n_r_max,0:l_max) :: e_t_r_l,e_t_r_l_global
    REAL(kind=8),DIMENSION(n_r_max,0:l_max) :: e_p_r_m,e_p_r_m_global
    REAL(kind=8),DIMENSION(n_r_max,0:l_max) :: e_t_r_m,e_t_r_m_global


    !-- end of declaration
    !---------------------------------------------------------------------

    IF ( BV.EQ.'V' ) THEN ! kinetic spectrum (correction of density)

       DO nR=1,n_r_max
          DO l=0,l_max
             e_p_r_l(nR,l)=0.D0
             e_t_r_l(nR,l)=0.D0
             e_p_r_m(nR,l)=0.D0
             e_t_r_m(nR,l)=0.D0
          END DO
          !DO lm=2,lm_max
          DO lm=MAX(llm,2),ulm
             l =lo_map%lm2l(lm)
             m =lo_map%lm2m(lm)
             e_p_temp= orho1(nR) * dLh(st_map%lm2(l,m)) * (           &
                  &      dLh(st_map%lm2(l,m))*or2(nR)*cc2real(b(lm,nR),m) &
                  &      + cc2real(db(lm,nR),m) )
             e_t_temp=orho1(nR)*dLh(st_map%lm2(l,m))*cc2real(aj(lm,nR),m)
             e_p_r_l(nR,l)=e_p_r_l(nR,l)+e_p_temp
             e_t_r_l(nR,l)=e_t_r_l(nR,l)+e_t_temp
             e_p_r_m(nR,m)=e_p_r_m(nR,m)+e_p_temp
             e_t_r_m(nR,m)=e_t_r_m(nR,m)+e_t_temp
          END DO    ! do loop over lms in block 
       END DO    ! radial grid points

    ELSE ! magnetic spectrum

       DO nR=1,n_r_max
          DO l=0,l_max
             e_p_r_l(nR,l)=0.D0
             e_t_r_l(nR,l)=0.D0
             e_p_r_m(nR,l)=0.D0
             e_t_r_m(nR,l)=0.D0
          END DO
          DO lm=MAX(2,llm),ulm
             l =lo_map%lm2l(lm)
             m =lo_map%lm2m(lm)
             e_p_temp=  dLh(st_map%lm2(l,m)) * (             &
                  &       dLh(st_map%lm2(l,m))*or2(nR)*cc2real(b(lm,nR),m) &
                  &       + cc2real(db(lm,nR),m) )
             e_t_temp=dLh(st_map%lm2(l,m))*cc2real(aj(lm,nR),m)
             e_p_r_l(nR,l)=e_p_r_l(nR,l)+e_p_temp
             e_t_r_l(nR,l)=e_t_r_l(nR,l)+e_t_temp
             e_p_r_m(nR,m)=e_p_r_m(nR,m)+e_p_temp
             e_t_r_m(nR,m)=e_t_r_m(nR,m)+e_t_temp
          END DO    ! do loop over lms in block 
       END DO    ! radial grid points

    END IF

    CALL MPI_Reduce(e_p_r_l,e_p_r_l_global,n_r_max*(l_max+1),&
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(e_t_r_l,e_t_r_l_global,n_r_max*(l_max+1),&
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(e_p_r_m,e_p_r_m_global,n_r_max*(l_max+1),&
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(e_t_r_m,e_t_r_m_global,n_r_max*(l_max+1),&
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)


    IF (rank.EQ.0) THEN
       !-- Radial Integrals:
       fac=0.5D0*eScale
       IF ( BV.EQ.'B' ) fac=fac*LFfac
       DO l=0,l_max
          e_p_l(l)  =fac*rInt_R(e_p_r_l_global(1,l),n_r_max,n_r_max,drx,      &
               &                           i_costf_init,d_costf_init)
          e_t_l(l)  =fac*rInt_R(e_t_r_l_global(1,l),n_r_max,n_r_max,drx,      &
               &                           i_costf_init,d_costf_init)
          e_p_m(l)  =fac*rInt_R(e_p_r_m_global(1,l),n_r_max,n_r_max,drx,      &
               &                           i_costf_init,d_costf_init)
          e_t_m(l)  =fac*rInt_R(e_t_r_m_global(1,l),n_r_max,n_r_max,drx,      &
               &                           i_costf_init,d_costf_init)
          IF ( BV.EQ.'B' ) THEN 
             e_cmb_l(l)=fac*e_p_r_l_global(1,l)
             e_cmb_m(l)=fac*e_p_r_m_global(1,l)
          END IF
       END DO


       !-- Averaging:
       IF ( n_time_ave.EQ.1 ) THEN
          DO l=0,l_max
             IF ( BV.EQ.'B' ) THEN
                e_p_l_ave(l)   =time_passed*e_p_l(l)
                e_t_l_ave(l)   =time_passed*e_t_l(l)
                e_p2_l_ave(l)  =time_passed*e_p_l(l)*e_p_l(l)
                e_t2_l_ave(l)  =time_passed*e_t_l(l)*e_t_l(l)
                e_p_m_ave(l)   =time_passed*e_p_m(l)
                e_t_m_ave(l)   =time_passed*e_t_m(l)
                e_p2_m_ave(l)  =time_passed*e_p_m(l)*e_p_m(l)
                e_t2_m_ave(l)  =time_passed*e_t_m(l)*e_t_m(l)
                e_cmb_l_ave(l) =time_passed*e_cmb_l(l)
                e_cmb2_l_ave(l)=time_passed*e_cmb_l(l)*e_cmb_l(l)
                e_cmb_m_ave(l) =time_passed*e_cmb_m(l)
                e_cmb2_m_ave(l)=time_passed*e_cmb_m(l)*e_cmb_m(l)
             ELSE
                ek_p_l_ave(l)   =time_passed*e_p_l(l)
                ek_t_l_ave(l)   =time_passed*e_t_l(l)
                ek_p2_l_ave(l)  =time_passed*e_p_l(l)*e_p_l(l)
                ek_t2_l_ave(l)  =time_passed*e_t_l(l)*e_t_l(l)
                ek_p_m_ave(l)   =time_passed*e_p_m(l)
                ek_t_m_ave(l)   =time_passed*e_t_m(l)
                ek_p2_m_ave(l)  =time_passed*e_p_m(l)*e_p_m(l)
                ek_t2_m_ave(l)  =time_passed*e_t_m(l)*e_t_m(l)
             END IF
          END DO
       ELSE
          DO l=0,l_max
             IF ( BV.EQ.'B' ) THEN
                e_p_l_ave(l)   =e_p_l_ave(l)   +  time_passed*e_p_l(l)
                e_t_l_ave(l)   =e_t_l_ave(l)   +  time_passed*e_t_l(l)
                e_p2_l_ave(l)  =e_p2_l_ave(l)  +  time_passed*e_p_l(l)*e_p_l(l)
                e_t2_l_ave(l)  =e_t2_l_ave(l)  +  time_passed*e_t_l(l)*e_t_l(l)
                e_p_m_ave(l)   =e_p_m_ave(l)   +  time_passed*e_p_m(l)
                e_t_m_ave(l)   =e_t_m_ave(l)   +  time_passed*e_t_m(l)
                e_p2_m_ave(l)  =e_p2_m_ave(l)  +  time_passed*e_p_m(l)*e_p_m(l)
                e_t2_m_ave(l)  =e_t2_m_ave(l)  +  time_passed*e_t_m(l)*e_t_m(l)
                e_cmb_l_ave(l) =e_cmb_l_ave(l) +  time_passed*e_cmb_l(l)
                e_cmb2_l_ave(l)=e_cmb2_l_ave(l)+  time_passed*e_cmb_l(l)*e_cmb_l(l)
                e_cmb_m_ave(l) =e_cmb_m_ave(l) +  time_passed*e_cmb_m(l)
                e_cmb2_m_ave(l)=e_cmb2_m_ave(l)+  time_passed*e_cmb_m(l)*e_cmb_m(l)
             ELSE
                ek_p_l_ave(l)   =ek_p_l_ave(l) +  time_passed*e_p_l(l)
                ek_t_l_ave(l)   =ek_t_l_ave(l)   + time_passed*e_t_l(l)
                ek_p2_l_ave(l)  =ek_p2_l_ave(l)  + time_passed*e_p_l(l)*e_p_l(l)
                ek_t2_l_ave(l)  =ek_t2_l_ave(l)  + time_passed*e_t_l(l)*e_t_l(l)
                ek_p_m_ave(l)   =ek_p_m_ave(l)   + time_passed*e_p_m(l)
                ek_t_m_ave(l)   =ek_t_m_ave(l)   + time_passed*e_t_m(l)
                ek_p2_m_ave(l)  =ek_p2_m_ave(l)  + time_passed*e_p_m(l)*e_p_m(l)
                ek_t2_m_ave(l)  =ek_t2_m_ave(l)  + time_passed*e_t_m(l)*e_t_m(l)
             END IF
          END DO
       END IF


       !-- Output: every 10th averaging step and at end of run
       IF ( l_stop_time .OR. MOD(n_time_ave,10).EQ.0 ) THEN

          !------ Output:
          dt_norm=1.d0/time_norm
          IF ( BV.EQ.'B' ) THEN
             outFile='mag_spec_ave.'//TAG
          ELSE IF ( BV.EQ.'V' ) THEN
             outFile='kin_spec_ave.'//TAG
          ELSE
             WRITE(*,*) 'WRONG BV INPUT TO spectrum_average!'
             STOP
          END IF
          nOut   =93
          OPEN(nOut,file=outFile,status='UNKNOWN')
          IF ( BV.EQ.'B' ) THEN
             DO l=0,l_max
                SDp_l = get_standard_deviation(dt_norm,e_p_l_ave(l),e_p2_l_ave(l))
                SDp_m = get_standard_deviation(dt_norm,e_p_m_ave(l),e_p2_m_ave(l))
                SDt_l = get_standard_deviation(dt_norm,e_t_l_ave(l),e_t2_l_ave(l))
                SDt_m = get_standard_deviation(dt_norm,e_t_m_ave(l),e_t2_m_ave(l))
                SDcmb_l = get_standard_deviation(dt_norm,e_cmb_l_ave(l),e_cmb2_l_ave(l))
                SDcmb_m = get_standard_deviation(dt_norm,e_cmb_m_ave(l),e_cmb2_m_ave(l))
                !WRITE(*,"(A,I4,5ES22.14)") "SDcmb_m = ",l,dt_norm,SDcmb_m,dt_norm*e_cmb2_m_ave(l),&
                !     & (dt_norm*e_cmb_m_ave(l))**2,&
                !     & dt_norm*e_cmb2_m_ave(l) - (dt_norm*e_cmb_m_ave(l))**2
                WRITE(93,'(2X,1P,I4,16D12.4)') l,                      &
                     &              dt_norm*e_p_l_ave(l),   dt_norm*e_p_m_ave(l),       &
                     &              dt_norm*e_t_l_ave(l),   dt_norm*e_t_m_ave(l),       &
                     &              dt_norm*e_cmb_l_ave(l), dt_norm*e_cmb_m_ave(l),     &
                     &              dt_norm*e_p_l_ave(l)+SDp_l,                         &
                     &              dt_norm*e_p_l_ave(l)-SDp_l,                         &
                     &              dt_norm*e_p_m_ave(l)+SDp_m,                         &
                     &              dt_norm*e_p_m_ave(l)-SDp_m,                         &
                     &              dt_norm*e_t_l_ave(l)+SDt_l,                         &
                     &              dt_norm*e_t_l_ave(l)-SDt_l,                         &
                     &              dt_norm*e_t_m_ave(l)+SDt_m,                         &
                     &              dt_norm*e_t_m_ave(l)-SDt_m,                         &
                     &              dt_norm*e_cmb_m_ave(l)+SDcmb_m,                     &
                     &              dt_norm*e_cmb_m_ave(l)-SDcmb_m 
             END DO
          ELSE
             DO l=0,l_max
                WRITE(93,'(2X,1P,I4,8D12.4)') l,                       &
                     &              dt_norm*ek_p_l_ave(l), dt_norm*ek_p_m_ave(l),       &
                     &              dt_norm*ek_t_l_ave(l), dt_norm*ek_t_m_ave(l),       &
                     &              dt_norm*ek_p2_l_ave(l),dt_norm*ek_p2_m_ave(l),      &
                     &              dt_norm*ek_t2_l_ave(l),dt_norm*ek_t2_m_ave(l)
             END DO
          END IF
          CLOSE(nOut)

          IF ( l_stop_time ) THEN
             CALL safeOpen(nLF,log_file)
             WRITE(nLF,'(/,                                            &
                  &'' ! TIME AVERAGED SPECTRA STORED IN FILE: '',a)') outFile
             WRITE(nLF,'(                                              &
                  &'' !              No. of averaged spectra: '',I5)') n_time_ave
             CALL safeClose(nLF)
          END IF

       END IF
    END IF

  END SUBROUTINE spectrum_average

  REAL(kind=8) FUNCTION get_standard_deviation(dt_norm,mean,sum_of_squares) result(stdev)
    REAL(kind=8),INTENT(in) :: dt_norm,mean,sum_of_squares

    REAL(kind=8) :: mean2,variance
    REAL(kind=8),PARAMETER :: tolerance = 5.0D0*EPSILON(1.0D0)

    mean2    = (dt_norm*mean)**2
    variance = ABS(dt_norm*sum_of_squares - mean2 )
    IF (variance/mean2.LT.tolerance) THEN
       stdev = 0.0D0
    ELSE
       stdev = SQRT(variance)
    END IF
  END FUNCTION get_standard_deviation
  !----------------------------------------------------------------------------
END MODULE spectrum_average_mod
