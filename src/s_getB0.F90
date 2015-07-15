!$Id$
!***********************************************************************
    SUBROUTINE getB0
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  |  Purpose of this subroutine is to solve for the Grenoble imposed  |
!  |  magnetic field by the inner core permanent magnet.               |
!  |                                                                   |
!  +-------------------------------------------------------------------+

    USE truncation
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE Grenoble
    USE blocking
    USE horizontal_data
    USE logic
    USE output_data
    USE const
    USE usefull, ONLY: cc2real
    USE parallel_mod,only:rank
    USE integration, ONLY: rInt_R

    IMPLICIT NONE

!-- output:
!       b0,db0,ddb0
!       Written into c_Grenoble.f

!-- local:
    INTEGER :: l,m,lm,nR
    REAL(kind=8) :: ePr(n_r_max),eP

!-- end of declaration
!--------------------------------------------------------------------------


    DO nR=1,n_r_max
        b0(nR)  =BIC*r_icb**2/y10_norm*(r_icb*or1(nR))
        db0(nR) =    -or1(nR)*b0(nR)
        ddb0(nR)=2.D0*or2(nR)*b0(nR)
    END DO

!----- Get magnetic energy:
    DO nR=1,n_r_max
        ePr(nR)=0.D0
        DO l=1,l_max
            m=0
            lm=lm2(l,0)
            ePr(nR)=ePr(nR) +   dLh(lm) * ( &
                dLh(lm)*or2(nR)*b0(nR)**2 + &
                               db0(nR)**2 )
        END DO    ! do loop over lms in block
    END DO    ! radial grid points
    eP=rInt_R(ePr,n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
    eP=0.5D0*LFfac*eScale*eP

    IF (rank == 0) THEN
       WRITE(*,'(1x,1P," ! Energy of imposed IC field:",D16.4)') eP
       IF ( l_save_out ) THEN
          OPEN(nLF,FILE=log_file,STATUS='UNKNOWN',POSITION='APPEND')
       END IF
       WRITE(nLF,'(1x,1P," ! Energy of imposed IC field:",D16.4)') eP
       IF ( l_save_out ) CLOSE(n_log_file)
    END IF
    RETURN
    end SUBROUTINE getB0

!--------------------------------------------------------------------------------
