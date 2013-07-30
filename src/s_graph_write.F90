!$Id$
!******************************************************************************
#include "perflib_preproc.cpp"
  SUBROUTINE graph_write(n_phis,n_thetas,dummy, &
       &                 PRECISION,FORMAT,n_graph_file)
!******************************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!------------------------------------------------------------------------------

!  called in graphout

!  This subroutine writes the data for one theta-band
!  (stored in 'dummy'). Version May, 5, 2000.

!------------------------------------------------------------------------------

    USE truncation

    IMPLICIT NONE

    INTEGER :: n_thetas            ! number of first colatitude value
    INTEGER :: n_phis              ! number of logitudes to be printed
    REAL(kind=4) :: dummy(n_phi_max,*)   ! data
    INTEGER :: precision           ! determins precision if output
    INTEGER :: format              ! formatted/unformatted output
    INTEGER :: n_graph_file        ! output unit

!-- Local:
    INTEGER :: n_phi,n_theta

!-- End of declaration
!------------------------------------------------------------------------------

    PERFON('gwrite')
    DO n_theta=1,n_thetas
         
        IF ( format == 0 ) THEN ! unformatted output
            WRITE(n_graph_file) &
                 (dummy(n_phi,n_theta),n_phi=1,n_phis)
        ELSE                ! formatted output
            IF ( precision == 0 ) THEN
                WRITE(n_graph_file,900) &
                     (dummy(n_phi,n_theta),n_phi=1,n_phis)
            ELSE IF ( precision == 1 ) THEN
                WRITE(n_graph_file,901) &
                     (dummy(n_phi,n_theta),n_phi=1,n_phis)
            END IF
        END IF

    END DO

    900 FORMAT(512(1X,F7.2))
    901 FORMAT(512(1X,F7.3))

    PERFOFF
    RETURN
    end SUBROUTINE graph_write

!------------------------------------------------------------------------------
