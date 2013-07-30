!$Id: s_graphOut.F90 413 2013-02-14 10:03:55Z gastine $
  
#include "perflib_preproc.cpp"
#include "intrinsic_sizes.h"
#define ONE_LARGE_BLOCK
  !*************************************************************************
  ! This routine has the same name as the non-MPI version. Which one is used
  ! is decided at compile time in the Makefile.
  ! At the moment, this is not packed into a module, due to the interface
  ! which is used with real and complex fields.
  !*************************************************************************
  SUBROUTINE graphOut_mpi(time,n_r,which_form,vr,vt,vp,br,bt,bp,sr, &
       &              n_theta_start,n_theta_block_size,lGraphHeader)
    !*************************************************************************

    !    !------------ This is release 2 level 1  --------------!
    !    !------------ Created on 1/17/02  by JW. --------------!

    !-------------------------------------------------------------------------
    !
    !  called in radialLoop
    !
    !  output of components of velocity, magnetic field vector and
    !  entropy for graphics. Version May, 23, 2000.
    !
    !  n_r: (input) for n_r = 0 a header is written
    !               for n_r > 0 values at radial level n_r are written
    !
    !  which_form : (input) = 0: unformatted
    !                       = 1: formatted writing
    !                       =-1: comment lines are included into file for
    !                            easier reading (cannot be used for graphics
    !                            processing in this form)
    !
    !  vr...sr: (input) arrays with grid-point values
    !
    !  n_theta_start : (input) values are written for theta-points :
    !                  n_theta_start <= n_theta <= n_theta_start-1+n_theta_block
    !
    !-------------------------------------------------------------------------

    USE truncation
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE blocking
    USE horizontal_data
    USE logic
    USE parallel_mod
    USE output_data, ONLY: graph_mpi_fh, runid

    IMPLICIT NONE

    REAL(kind=8),INTENT(IN) :: time

    INTEGER,INTENT(IN) :: n_r                          ! radial grod point no.
    INTEGER,INTENT(IN) :: which_form                   ! determins format
    INTEGER,INTENT(IN) :: n_theta_start                ! start theta no.
    INTEGER,INTENT(IN) :: n_theta_block_size           ! size of theta block
    LOGICAL,INTENT(INOUT) :: lGraphHeader

    REAL(kind=8),intent(IN) :: vr(nrp,*),vt(nrp,*),vp(nrp,*)
    REAL(kind=8),intent(IN) :: br(nrp,*),bt(nrp,*),bp(nrp,*)
    REAL(kind=8),intent(IN) :: sr(nrp,*)

    !-- Local:
    INTEGER :: n_phi         ! counter for longitude
    INTEGER :: n_theta       ! counter for colatitude
    INTEGER :: n_theta_cal   ! position of block colat in all colats
    INTEGER :: n_theta_stop  ! end theta no.

    INTEGER :: precision

    REAL(kind=8) :: fac,fac_r
    REAL(kind=4) :: dummy(n_phi_max,nfs)

    CHARACTER(len=20) :: version

    ! MPI related variables
    INTEGER :: info
    INTEGER :: status(MPI_STATUS_SIZE)
    CHARACTER(len=MPI_MAX_ERROR_STRING) :: error_string
    !INTEGER :: count
    INTEGER :: length_of_error,bytes_written
    INTEGER :: size_of_header, size_of_data_per_rank, size_of_data_per_r
    integer :: size_of_data_per_thetaB
    integer(KIND=MPI_OFFSET_KIND) :: disp
    INTEGER :: etype,filetype
    character(len=MPI_MAX_DATAREP_STRING) :: datarep
    ! end of MPI related variables

    !-- End of Declaration
    !----------------------------------------------------------------------
    !PRINT*,"lGraphHeader = ",lGraphHeader
    !WRITE(*,"(A,I5,A,I5)") "Theta block start at ",n_theta_start,", size is ",n_theta_block_size

    IF ( lGraphHeader ) THEN
       !PRINT*,"Setting the View"
       !WRITE(*,"(A,4I5)") "n_phi_max,n_theta_block_size,nThetaBs,nr_per_rank: ",n_phi_max,n_theta_block_size,nThetaBs,nr_per_rank
       size_of_header = 8+LEN(version)+8+LEN(runid)+8+13*SIZEOF_INTEGER+8+n_theta_max*SIZEOF_REAL

#ifdef ONE_LARGE_BLOCK
       size_of_data_per_thetaB = 8+4*SIZEOF_REAL+4*(8+n_phi_max*SIZEOF_REAL*n_theta_block_size)
       IF (l_mag) size_of_data_per_thetaB = size_of_data_per_thetaB + 3*(8+n_phi_max*SIZEOF_REAL*n_theta_block_size)
#else
       size_of_data_per_thetaB = 8+4*SIZEOF_REAL+4*(8+n_phi_max*SIZEOF_REAL)*n_theta_block_size
       IF (l_mag) size_of_data_per_thetaB = size_of_data_per_thetaB + 3*(8+n_phi_max*SIZEOF_REAL)*n_theta_block_size
#endif
       size_of_data_per_r = size_of_data_per_thetaB * nThetaBs
       size_of_data_per_rank = size_of_data_per_r * nr_per_rank

       !PRINT*,"size_of_header = ",size_of_header,", size_of_data/rank = ",size_of_data_per_rank
       IF (rank.EQ.0) THEN
          ! rank zero writes the Header
          disp = 0
          CALL MPI_FILE_SET_VIEW(graph_mpi_fh,disp,MPI_CHARACTER,MPI_CHARACTER,"external32",MPI_INFO_NULL,ierr)
       ELSE
          disp = size_of_header+rank*size_of_data_per_rank
          CALL MPI_FILE_SET_VIEW(graph_mpi_fh,disp,&
               & MPI_CHARACTER,MPI_CHARACTER,"external32",MPI_INFO_NULL,ierr)
       END IF
       !CALL MPI_ERROR_STRING(ierr,error_string,length_of_error,ierr)
       !PRINT*,"MPI_FILE_SET_VIEW returned: ",TRIM(error_string)

       CALL mpi_file_get_view(graph_mpi_fh,disp,etype,filetype,datarep,ierr)
       !PRINT*,"view = ",disp,etype,filetype,datarep

       bytes_written = 0
       !-- Write header & colatitudes for n_r=0:
       IF (rank.EQ.0) THEN
          IF ( which_form /= 0 ) THEN
#if 0
             !----- Formatted output:
             version='Graphout_Version_6'

             !-------- Write parameters:
             WRITE(n_graph_file,'(a)')  version ! identifying header
             ! for IDL magsym5 routine
             WRITE(n_graph_file,'(A64)') runid
             IF ( which_form < 0 ) write(n_graph_file,'(/ &
                  &  "Time,truncation,outputgrid:")')
             WRITE(n_graph_file,'(e12.5,1x,10i6)') &
                  time,n_r_max,n_theta_max,n_phi_tot, &
                  n_r_ic_max-1,minc,nThetaBs
             IF ( which_form < 0 ) write(n_graph_file, &
                  & '(/"Ra, Ek, Pr, Pm, radratio, sigma_ratio:")')
             WRITE(n_graph_file,'(6e14.6)') &
                  ra,ek,pr,prmag,radratio,sigma_ratio

             !-------- Write colatitudes:
             IF ( which_form < 0 ) THEN
                WRITE(n_graph_file, &
                     &   '(/"Colatitudes:", /"point#, theta")')
                DO n_theta=1,n_theta_max/2   ! NHS
                   WRITE(n_graph_file,'(i4, f9.4)') &
                        n_theta,theta_ord(n_theta)
                END DO
                WRITE(n_graph_file,'("--- equator ---")')
                DO n_theta=n_theta_max/2+1,n_theta_max   ! SHS
                   WRITE(n_graph_file,'(i4, f9.4)') &
                        n_theta,theta_ord(n_theta)
                END DO
             ELSE IF ( which_form >= 0 ) then
                WRITE(n_graph_file,'(129(1x,f8.5))') &
                     (theta_ord(n_theta),n_theta=1,n_theta_max)
             END IF
#endif
          ELSE

             !----- Unformatted output:
             version='Graphout_Version_9'

             !-------- Write parameters:
             CALL MPI_FILE_WRITE(graph_mpi_fh,len(version),1,MPI_INTEGER,status,ierr)
             !CALL mpi_get_count(status,MPI_INTEGER,count,ierr)
             !bytes_written = bytes_written + count*SIZEOF_INTEGER
             CALL MPI_FILE_WRITE(graph_mpi_fh,version,LEN(version),MPI_CHARACTER,status,ierr)
             !CALL mpi_get_count(status,MPI_CHARACTER,count,ierr)
             !bytes_written = bytes_written + count*SIZEOF_CHARACTER
             CALL MPI_FILE_WRITE(graph_mpi_fh,len(version),1,MPI_INTEGER,status,ierr)
             !CALL mpi_get_count(status,MPI_INTEGER,count,ierr)
             !bytes_written = bytes_written + count*SIZEOF_INTEGER


             CALL MPI_FILE_WRITE(graph_mpi_fh,LEN(runid),1,MPI_INTEGER,status,ierr)
             !CALL mpi_get_count(status,MPI_INTEGER,count,ierr)
             !bytes_written = bytes_written + count*SIZEOF_INTEGER
             CALL MPI_FILE_WRITE(graph_mpi_fh,runid,len(runid),MPI_CHARACTER,status,ierr)
             !CALL mpi_get_count(status,MPI_CHARACTER,count,ierr)
             !bytes_written = bytes_written + count*SIZEOF_CHARACTER
             CALL MPI_FILE_WRITE(graph_mpi_fh,len(runid),1,MPI_INTEGER,status,ierr)
             !CALL mpi_get_count(status,MPI_INTEGER,count,ierr)
             !bytes_written = bytes_written + count*SIZEOF_INTEGER

             CALL MPI_FILE_WRITE(graph_mpi_fh,13*4,1,MPI_INTEGER,status,ierr)
             !CALL mpi_get_count(status,MPI_INTEGER,count,ierr)
             !bytes_written = bytes_written + count*SIZEOF_INTEGER
             CALL MPI_FILE_WRITE(graph_mpi_fh,SNGL(time),1,MPI_REAL,status,ierr)
             !CALL mpi_get_count(status,MPI_REAL,count,ierr)
             !bytes_written = bytes_written + count*SIZEOF_REAL
             CALL MPI_FILE_WRITE(graph_mpi_fh,FLOAT(n_r_max),1,MPI_REAL,status,ierr)
             !CALL mpi_get_count(status,MPI_REAL,count,ierr)
             !bytes_written = bytes_written + count*SIZEOF_REAL
             CALL MPI_FILE_WRITE(graph_mpi_fh,FLOAT(n_theta_max),1,MPI_REAL,status,ierr)
             !CALL mpi_get_count(status,MPI_REAL,count,ierr)
             !bytes_written = bytes_written + count*SIZEOF_REAL
             CALL MPI_FILE_WRITE(graph_mpi_fh,FLOAT(n_phi_tot),1,MPI_REAL,status,ierr)
             !CALL mpi_get_count(status,MPI_REAL,count,ierr)
             !bytes_written = bytes_written + count*SIZEOF_REAL
             CALL MPI_FILE_WRITE(graph_mpi_fh,FLOAT(n_r_ic_max-1),1,MPI_REAL,status,ierr)
             !CALL mpi_get_count(status,MPI_REAL,count,ierr)
             !bytes_written = bytes_written + count*SIZEOF_REAL
             CALL MPI_FILE_WRITE(graph_mpi_fh,FLOAT(minc),1,MPI_REAL,status,ierr)
             !CALL mpi_get_count(status,MPI_REAL,count,ierr)
             !bytes_written = bytes_written + count*SIZEOF_REAL
             CALL MPI_FILE_WRITE(graph_mpi_fh,FLOAT(nThetaBs),1,MPI_REAL,status,ierr)
             !CALL mpi_get_count(status,MPI_REAL,count,ierr)
             !bytes_written = bytes_written + count*SIZEOF_REAL
             CALL MPI_FILE_WRITE(graph_mpi_fh,SNGL(ra),1,MPI_REAL,status,ierr)
             !CALL mpi_get_count(status,MPI_REAL,count,ierr)
             !bytes_written = bytes_written + count*SIZEOF_REAL
             CALL MPI_FILE_WRITE(graph_mpi_fh,SNGL(ek),1,MPI_REAL,status,ierr)
             !CALL mpi_get_count(status,MPI_REAL,count,ierr)
             !bytes_written = bytes_written + count*SIZEOF_REAL
             CALL MPI_FILE_WRITE(graph_mpi_fh,SNGL(pr),1,MPI_REAL,status,ierr)
             !CALL mpi_get_count(status,MPI_REAL,count,ierr)
             !bytes_written = bytes_written + count*SIZEOF_REAL
             CALL MPI_FILE_WRITE(graph_mpi_fh,SNGL(prmag),1,MPI_REAL,status,ierr)
             !CALL mpi_get_count(status,MPI_REAL,count,ierr)
             !bytes_written = bytes_written + count*SIZEOF_REAL
             CALL MPI_FILE_WRITE(graph_mpi_fh,SNGL(radratio),1,MPI_REAL,status,ierr)
             !CALL mpi_get_count(status,MPI_REAL,count,ierr)
             !bytes_written = bytes_written + count*SIZEOF_REAL
             CALL MPI_FILE_WRITE(graph_mpi_fh,SNGL(sigma_ratio),1,MPI_REAL,status,ierr)
             !CALL mpi_get_count(status,MPI_REAL,count,ierr)
             !bytes_written = bytes_written + count*SIZEOF_REAL
             CALL MPI_FILE_WRITE(graph_mpi_fh,13*4,1,MPI_INTEGER,status,ierr)
             !CALL mpi_get_count(status,MPI_INTEGER,count,ierr)
             !bytes_written = bytes_written + count*SIZEOF_INTEGER

             !-------- Write colatitudes:
             CALL MPI_FILE_WRITE(graph_mpi_fh,n_theta_max*SIZEOF_REAL,1,MPI_INTEGER,status,ierr)
             !CALL mpi_get_count(status,MPI_INTEGER,count,ierr)
             !bytes_written = bytes_written + count*SIZEOF_INTEGER
             DO n_theta=1,n_theta_max
                CALL MPI_FILE_WRITE(graph_mpi_fh,SNGL(theta_ord(n_theta)),1,MPI_REAL,status,ierr)
                !CALL mpi_get_count(status,MPI_REAL,count,ierr)
                !bytes_written = bytes_written + count*SIZEOF_REAL
             END DO
             CALL MPI_FILE_WRITE(graph_mpi_fh,n_theta_max*SIZEOF_REAL,1,MPI_INTEGER,status,ierr)
             !CALL mpi_get_count(status,MPI_INTEGER,count,ierr)
             !bytes_written = bytes_written + count*SIZEOF_INTEGER

          END IF
       END IF
       lGraphHeader=.FALSE.
       !PRINT*,"For the header, we wrote ",bytes_written," bytes."
    ELSE  ! Call not for writing header

       !PERFON('mw_data')
       bytes_written=0

       !-- Determine radius and thetas in this block:
       n_theta_stop=n_theta_start+n_theta_block_size-1
#if 0
       IF ( which_form < 0 ) WRITE(n_graph_file, &
            &   '(/"--- Radial level, radius, i1, i2 ---")')

       IF ( which_form /= 0 ) THEN
          WRITE(n_graph_file,'(I4, F9.5, 2I4)') &
               n_r-1,r(n_r)/r(1), &
               n_theta_start,n_theta_stop
       ELSE
#endif
          CALL MPI_FILE_WRITE(graph_mpi_fh,4*SIZEOF_REAL,1,MPI_INTEGER,status,ierr)
          !CALL mpi_get_count(status,MPI_INTEGER,count,ierr)
          !bytes_written = bytes_written + count*SIZEOF_INTEGER
          CALL MPI_FILE_WRITE(graph_mpi_fh,FLOAT(n_r-1),1,MPI_REAL,status,ierr)
          !CALL mpi_get_count(status,MPI_REAL,count,ierr)
          !bytes_written = bytes_written + count*SIZEOF_REAL
          CALL MPI_FILE_WRITE(graph_mpi_fh,SNGL(r(n_r)/r(1)),1,MPI_REAL,status,ierr)
          !CALL mpi_get_count(status,MPI_REAL,count,ierr)
          !bytes_written = bytes_written + count*SIZEOF_REAL
          CALL MPI_FILE_WRITE(graph_mpi_fh,FLOAT(n_theta_start),1,MPI_REAL,status,ierr)
          !CALL mpi_get_count(status,MPI_REAL,count,ierr)
          !bytes_written = bytes_written + count*SIZEOF_REAL
          CALL MPI_FILE_WRITE(graph_mpi_fh,FLOAT(n_theta_stop),1,MPI_REAL,status,ierr)
          !CALL mpi_get_count(status,MPI_REAL,count,ierr)
          !bytes_written = bytes_written + count*SIZEOF_REAL
          CALL MPI_FILE_WRITE(graph_mpi_fh,4*SIZEOF_REAL,1,MPI_INTEGER,status,ierr)
          !CALL mpi_get_count(status,MPI_INTEGER,count,ierr)
          !bytes_written = bytes_written + count*SIZEOF_INTEGER
#if 0
       END IF
#endif


       !-- Write entropy:
       precision=1
       !IF ( which_form < 0 ) WRITE(n_graph_file,'(/" S: ")')
       DO n_theta=1,n_theta_block_size,2
          DO n_phi=1,n_phi_max ! do loop over phis
             dummy(n_phi,n_theta)  =SNGL(sr(n_phi,n_theta))   ! NHS
             dummy(n_phi,n_theta+1)=SNGL(sr(n_phi,n_theta+1)) ! SHS
          END DO
       END DO
       CALL graph_write_mpi(n_phi_max,n_theta_block_size,dummy, &
            precision,which_form,graph_mpi_fh)

       precision=0

       !-- Calculate and write radial velocity:
       !IF ( which_form < 0 ) WRITE(n_graph_file,'(/" Vr: ")')
       fac=or2(n_r)*vScale*orho1(n_r)
       DO n_theta=1,n_theta_block_size,2
          DO n_phi=1,n_phi_max
             dummy(n_phi,n_theta)  =SNGL(fac*vr(n_phi,n_theta))
             dummy(n_phi,n_theta+1)=SNGL(fac*vr(n_phi,n_theta+1))
          END DO
       END DO
       CALL graph_write_mpi(n_phi_max,n_theta_block_size,dummy, &
            precision,which_form,graph_mpi_fh)


       !-- Calculate and write latitudinal velocity:
       !IF ( which_form < 0 ) WRITE(n_graph_file,'(/" Vt: ")')
       fac_r=or1(n_r)*vScale*orho1(n_r)
       DO n_theta=1,n_theta_block_size,2
          n_theta_cal=n_theta_start+n_theta-1
          fac=fac_r*O_sin_theta(n_theta_cal)
          DO n_phi=1,n_phi_max
             dummy(n_phi,n_theta)  =SNGL(fac*vt(n_phi,n_theta))
             dummy(n_phi,n_theta+1)=SNGL(fac*vt(n_phi,n_theta+1))
          END DO
       END DO
       CALL graph_write_mpi(n_phi_max,n_theta_block_size,dummy, &
            precision,which_form,graph_mpi_fh)

       !-- Calculate and write longitudinal velocity:
       !IF ( which_form < 0 ) WRITE(n_graph_file,'(/" Vp: ")')
       fac_r=or1(n_r)*vScale*orho1(n_r)
       DO n_theta=1,n_theta_block_size,2
          n_theta_cal=n_theta_start+n_theta-1
          fac=fac_r*O_sin_theta(n_theta_cal)
          DO n_phi=1,n_phi_max
             dummy(n_phi,n_theta)  =SNGL(fac*vp(n_phi,n_theta))
             dummy(n_phi,n_theta+1)=SNGL(fac*vp(n_phi,n_theta+1))
          END DO
       END DO
       CALL graph_write_mpi(n_phi_max,n_theta_block_size,dummy, &
            precision,which_form,graph_mpi_fh)

       IF ( l_mag ) THEN

          !-- Calculate and write radial magnetic field:
          !IF ( which_form < 0 ) WRITE(n_graph_file,'(/'' Br: '')')
          fac=or2(n_r)
          DO n_theta=1,n_theta_block_size,2
             DO n_phi=1,n_phi_max
                dummy(n_phi,n_theta)  =SNGL(fac*br(n_phi,n_theta))
                dummy(n_phi,n_theta+1)=SNGL(fac*br(n_phi,n_theta+1))
             END DO
          END DO
          CALL graph_write_mpi(n_phi_max,n_theta_block_size,dummy, &
               precision,which_form,graph_mpi_fh)

          !-- Calculate and write latitudinal magnetic field:
          !IF ( which_form < 0 ) WRITE(n_graph_file,'(/'' Bt: '')')
          DO n_theta=1,n_theta_block_size,2
             n_theta_cal=n_theta_start+n_theta-1
             fac=or1(n_r)*O_sin_theta(n_theta_cal)
             DO n_phi=1,n_phi_max
                dummy(n_phi,n_theta)  =SNGL(fac*bt(n_phi,n_theta))
                dummy(n_phi,n_theta+1)=SNGL(fac*bt(n_phi,n_theta+1))
             END DO
          END DO
          CALL graph_write_mpi(n_phi_max,n_theta_block_size,dummy, &
               precision,which_form,graph_mpi_fh)

          !-- Calculate and write longitudinal magnetic field:
          !IF ( which_form < 0 ) WRITE(n_graph_file,'(/'' Bp: '')')
          DO n_theta=1,n_theta_block_size,2
             n_theta_cal=n_theta_start+n_theta-1
             fac=or1(n_r)*O_sin_theta(n_theta_cal)
             DO n_phi=1,n_phi_max
                dummy(n_phi,n_theta)  =SNGL(fac*bp(n_phi,n_theta))
                dummy(n_phi,n_theta+1)=SNGL(fac*bp(n_phi,n_theta+1))
             END DO
          END DO
          CALL graph_write_mpi(n_phi_max,n_theta_block_size,dummy, &
               precision,which_form,graph_mpi_fh)

       END IF ! l_mag ?

       !WRITE(*,"(A,I8)") "bytes_written = ",bytes_written

       !PERFOFF
    END IF

  END SUBROUTINE graphOut_mpi

  SUBROUTINE graph_write_mpi(n_phis,n_thetas,dummy,PRECISION,which_form,graph_mpi_fh)
    !******************************************************************************

    !    !------------ This is release 2 level 1  --------------!
    !    !------------ Created on 1/17/02  by JW. --------------!

    !------------------------------------------------------------------------------

    !  called in graphout

    !  This subroutine writes the data for one theta-band
    !  (stored in 'dummy'). Version May, 5, 2000.

    !------------------------------------------------------------------------------

    use parallel_mod
    USE truncation
    IMPLICIT NONE

    INTEGER :: n_thetas            ! number of first colatitude value
    INTEGER :: n_phis              ! number of logitudes to be printed
    REAL(kind=4) :: dummy(n_phi_max,*)   ! data
    INTEGER :: precision           ! determins precision if output
    INTEGER :: which_form              ! formatted/unformatted output
    INTEGER :: graph_mpi_fh        ! mpi handle of the mpi file

    !-- Local:
    INTEGER :: n_phi,n_theta

    ! MPI related variables
    INTEGER :: status(MPI_STATUS_SIZE),count
    !-- End of declaration
    !------------------------------------------------------------------------------
#ifdef ONE_LARGE_BLOCK
    !PERFON('gwrite_M')
    ! this is for Graphout_Version_9
    CALL MPI_FILE_WRITE(graph_mpi_fh,n_phis*n_thetas*SIZEOF_REAL,1,MPI_INTEGER,status,ierr)
    CALL MPI_FILE_WRITE(graph_mpi_fh,dummy(:,1:n_thetas),n_phis*n_thetas,MPI_REAL,status,ierr)
    CALL MPI_FILE_WRITE(graph_mpi_fh,n_phis*n_thetas*SIZEOF_REAL,1,MPI_INTEGER,status,ierr)
    !PERFOFF
    !CALL mpi_get_count(status,MPI_REAL,count,ierr)
    !PRINT*,"count = ",count
#else
    !PERFON('gwrite_M')
    DO n_theta=1,n_thetas

       IF ( which_form == 0 ) THEN ! unformatted output
          CALL MPI_FILE_WRITE(graph_mpi_fh,n_phis*SIZEOF_REAL,1,MPI_INTEGER,status,ierr)
          CALL MPI_FILE_WRITE(graph_mpi_fh,dummy(1,n_theta),n_phis,MPI_REAL,status,ierr)
          CALL MPI_FILE_WRITE(graph_mpi_fh,n_phis*SIZEOF_REAL,1,MPI_INTEGER,status,ierr)

          !CALL MPI_FILE_WRITE(n_graph_file) &
          !      (dummy(n_phi,n_theta),n_phi=1,n_phis)
       ELSE                ! formatted output
#if 0  
          IF ( precision == 0 ) THEN
             WRITE(n_graph_file,900) &
                  (dummy(n_phi,n_theta),n_phi=1,n_phis)
          ELSE IF ( precision == 1 ) THEN
             WRITE(n_graph_file,901) &
                  (dummy(n_phi,n_theta),n_phi=1,n_phis)
          END IF
#endif
       END IF

    END DO

900 FORMAT(512(1X,F7.2))
901 FORMAT(512(1X,F7.3))

    !PERFOFF
#endif
  END SUBROUTINE graph_write_mpi

  !------------------------------------------------------------------------------

