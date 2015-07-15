!$Id: s_store.F90 418 2013-02-14 12:45:13Z gastine $
!********************************************************************
#include "intrinsic_sizes.h"
MODULE store_rst
  use truncation
  use output_data
  use parallel_mod
  USE logic, ONLY: l_heat, l_mag, l_cond_ic
  USE num_param, only: tScale
  USE physical_parameters, ONLY: ra, pr, prmag, ek, radratio, sigma_ratio
  USE init_fields, ONLY: omega_ic1, omegaOsz_ic1, tOmega_ic1, &
       &                 omega_ic2, omegaOsz_ic2, tOmega_ic2, &
       &                 omega_ma1, omegaOsz_ma1, tOmega_ma1, &
       &                 omega_ma2, omegaOsz_ma2, tOmega_ma2
  USE fieldsLast, ONLY: dsdtLast,dwdtLast, dzdtLast, dpdtLast, dbdtLast, djdtLast, &
       &dbdt_icLast, djdt_icLast, lorentz_torque_icLast, lorentz_torque_maLast 
  IMPLICIT NONE

  private
  INTEGER :: header_length,size_of_field,per_process,size_of_field_ic

  PUBLIC :: store_mpi
contains
  SUBROUTINE store_mpi(time,dt,dtNew,w,z,p,s,b,aj,b_ic,aj_ic)
    !********************************************************************

    !--------------------------------------------------------------------

    ! *** store results on disc file (restart file)
    !   In addition to the magnetic field and velocity potentials
    !   we store the time derivative terms
    !     djdt(lm,nR),dbdt(lm,nR), ......

    !--------------------------------------------------------------------

    ! File format is at the moment identical to the old one. Can be optimized lateron.
    ! (RM = Fortran record marker), magnetic field depends on the respective switch.
    ! header
    ! RM,w,z,p(,s),RM
    ! RM,(dsdtLast,)dwdtLast,dzdtLast,dpdtLast,RM
    ! RM,b,aj,dbdtLast,djdtLast,RM
    ! RM,b_ic,aj_ic,dbdt_icLast,djdt_icLast,RM
    ! RM,15 REAL8,RM


    IMPLICIT NONE

    !-- Input of variables:
    REAL(kind=8),INTENT(IN) :: time,dt,dtNew

    !-- Input of scalar fields to be stored:
    COMPLEX(kind=8),INTENT(IN) :: w(lm_max,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: z(lm_max,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: p(lm_max,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: s(lm_max,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: b(lm_maxMag,n_r_maxMag)
    COMPLEX(kind=8),INTENT(IN) :: aj(lm_maxMag,n_r_maxMag)
    COMPLEX(kind=8),INTENT(IN) :: b_ic(lm_maxMag,n_r_ic_maxMag)
    COMPLEX(kind=8),INTENT(IN) :: aj_ic(lm_maxMag,n_r_ic_maxMag)

    integer :: inform
    integer :: header_length
    integer(KIND=MPI_OFFSET_KIND) :: disp
    integer :: status(MPI_STATUS_SIZE)
    INTEGER :: number_of_fields,per_process, data_to_write
    INTEGER :: filetype
    !-- end of declaration
    !---------------------------------------------------------------------

    !-- Write parameters:
    IF ( .NOT. l_heat ) THEN
       inform=11
    ELSE
       inform=12
    END IF

    IF ( l_heat ) THEN
       number_of_fields = 4
    ELSE
       number_of_fields = 3
    END IF

    ! the length of the header, which is only written by the first process
    header_length=8*SIZEOF_DOUBLE_PRECISION+7*SIZEOF_INTEGER
    ! also the recordmarker for the following fields w,z,p(,s) is added to the header length
    !header_length = header_length+4
    IF (rank == 0) THEN
       disp = 0
    ELSE
       disp = 2*SIZEOF_INTEGER + header_length
    END IF
    CALL MPI_FILE_SET_VIEW(rst_mpi_fh,disp,MPI_BYTE,MPI_BYTE,"external32",MPI_INFO_NULL,ierr)

    IF (rank == 0) THEN
       ! write the header to file in process 0
       CALL MPI_FILE_WRITE(rst_mpi_fh,header_length,1,MPI_INTEGER,status,ierr)

       CALL MPI_FILE_WRITE(rst_mpi_fh,time*tScale,1,MPI_DOUBLE_PRECISION,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,dt*tScale,1,MPI_DOUBLE_PRECISION,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,ra,1,MPI_DOUBLE_PRECISION,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,pr,1,MPI_DOUBLE_PRECISION,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,prmag,1,MPI_DOUBLE_PRECISION,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,ek,1,MPI_DOUBLE_PRECISION,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,radratio,1,MPI_DOUBLE_PRECISION,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,inform,1,MPI_INTEGER,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,n_r_max,1,MPI_INTEGER,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,n_theta_max,1,MPI_INTEGER,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,n_phi_tot,1,MPI_INTEGER,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,minc,1,MPI_INTEGER,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,nalias,1,MPI_INTEGER,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,n_r_ic_max,1,MPI_INTEGER,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,sigma_ratio,1,MPI_DOUBLE_PRECISION,status,ierr)
       
       CALL MPI_FILE_WRITE(rst_mpi_fh,header_length,1,MPI_INTEGER,status,ierr)

       WRITE(*,"(A,8ES20.13)") "Output in restart file: ",time*tScale,dt*tScale, &
            & ra,pr,prmag,ek,radratio, sigma_ratio

       CALL MPI_FILE_WRITE(rst_mpi_fh,&
            &number_of_fields*lm_max*n_r_max*SIZEOF_DOUBLE_COMPLEX,&
            &1,MPI_INTEGER,status,ierr)
    END IF

    ! define a datatype for the filetype which skips the data from other ranks
    per_process = lm_max*nr_per_rank
    IF (rank == n_procs-1) THEN
       CALL MPI_Type_vector(1,lm_max*(nr_per_rank+1), lm_max*n_r_max, MPI_DOUBLE_COMPLEX, filetype,ierr)
    ELSE
       CALL MPI_Type_vector(1,per_process,            lm_max*n_r_max, MPI_DOUBLE_COMPLEX, filetype,ierr)
    END IF
    CALL MPI_Type_commit(filetype,ierr)

    ! now set the view after the header
    disp = 8+header_length+4+rank*per_process*SIZEOF_DOUBLE_COMPLEX
    CALL MPI_FILE_SET_VIEW(rst_mpi_fh,disp,MPI_DOUBLE_COMPLEX,filetype,"external32",MPI_INFO_NULL,ierr)

    IF (rank < n_procs-1) THEN
       data_to_write = per_process
    ELSE
       ! on the last process, we have to write one radial point more
       data_to_write = lm_max*(nr_per_rank+1)
    END IF

    WRITE(*,"(A,I4,A,2ES20.13)") "w(1,",1+rank*nr_per_rank,") = ",w(1,1+rank*nr_per_rank)
    CALL MPI_FILE_WRITE(rst_mpi_fh,w(:,1+rank*nr_per_rank),data_to_write,MPI_DOUBLE_COMPLEX,status,ierr)
    CALL MPI_FILE_WRITE(rst_mpi_fh,z(:,1+rank*nr_per_rank),data_to_write,MPI_DOUBLE_COMPLEX,status,ierr)
    CALL MPI_FILE_WRITE(rst_mpi_fh,p(:,1+rank*nr_per_rank),data_to_write,MPI_DOUBLE_COMPLEX,status,ierr)
    IF ( l_heat ) THEN
       CALL MPI_FILE_WRITE(rst_mpi_fh,s(:,1+rank*nr_per_rank),data_to_write,MPI_DOUBLE_COMPLEX,status,ierr)
    END IF

    ! now set the view after the first block to write the RM
    disp = 8+header_length+4+number_of_fields*lm_max*n_r_max*SIZEOF_DOUBLE_COMPLEX
    CALL MPI_FILE_SET_VIEW(rst_mpi_fh,disp,MPI_BYTE,MPI_BYTE,"external32",MPI_INFO_NULL,ierr)

    IF (rank  ==  n_procs-1) THEN
       ! write the recordmarker for the first fields
       CALL MPI_FILE_WRITE(rst_mpi_fh,&
            &number_of_fields*lm_max*n_r_max*SIZEOF_DOUBLE_COMPLEX,&
            &1,MPI_INTEGER,status,ierr)
       
       ! write the starting record marker for the ...Last fields
       CALL MPI_FILE_WRITE(rst_mpi_fh,&
            &number_of_fields*lm_max*n_r_max*SIZEOF_DOUBLE_COMPLEX,&
            &1,MPI_INTEGER,status,ierr)
    END IF

    ! now set the view after the RM of the first field block
    disp = 8+header_length+8+number_of_fields*lm_max*n_r_max*SIZEOF_DOUBLE_COMPLEX&
         & + 4 + rank*per_process
    CALL MPI_FILE_SET_VIEW(rst_mpi_fh,disp,MPI_DOUBLE_COMPLEX,filetype,"external32",MPI_INFO_NULL,ierr)

    if ( l_heat ) then
       CALL MPI_FILE_WRITE(rst_mpi_fh,dsdtLast(:,1+rank*nr_per_rank),data_to_write,MPI_DOUBLE_COMPLEX,status,ierr)
    END IF
    CALL MPI_FILE_WRITE(rst_mpi_fh,dwdtLast(:,1+rank*nr_per_rank),data_to_write,MPI_DOUBLE_COMPLEX,status,ierr)
    CALL MPI_FILE_WRITE(rst_mpi_fh,dzdtLast(:,1+rank*nr_per_rank),data_to_write,MPI_DOUBLE_COMPLEX,status,ierr)
    CALL MPI_FILE_WRITE(rst_mpi_fh,dpdtLast(:,1+rank*nr_per_rank),data_to_write,MPI_DOUBLE_COMPLEX,status,ierr)

    ! now set the view after the second block to write the RM
    disp = 8+header_length+8+number_of_fields*lm_max*n_r_max*SIZEOF_DOUBLE_COMPLEX&
         & + 4 + number_of_fields*lm_max*n_r_max*SIZEOF_DOUBLE_COMPLEX
    CALL MPI_FILE_SET_VIEW(rst_mpi_fh,disp,MPI_BYTE,MPI_BYTE,"external32",MPI_INFO_NULL,ierr)

    IF (rank  ==  n_procs-1) THEN
       ! write the recordmarker for the first fields
       CALL MPI_FILE_WRITE(rst_mpi_fh,&
            &number_of_fields*lm_max*n_r_max*SIZEOF_DOUBLE_COMPLEX,&
            &1,MPI_INTEGER,status,ierr)
       
    END IF

    !-- Write magnetic field:
    IF ( l_mag ) THEN
       IF ( rank == n_procs-1 ) THEN
          ! write the starting record marker for the magnetic fields
          CALL MPI_FILE_WRITE(rst_mpi_fh,&
               &4*lm_max*n_r_max*SIZEOF_DOUBLE_COMPLEX,&
               &1,MPI_INTEGER,status,ierr)
       END IF
       ! now set the view after the RM of the first field block
       disp = 8+header_length+8+number_of_fields*lm_max*n_r_max*SIZEOF_DOUBLE_COMPLEX&
            & + 8 + number_of_fields*lm_max*n_r_max*SIZEOF_DOUBLE_COMPLEX + 4
       CALL MPI_FILE_SET_VIEW(rst_mpi_fh,disp,MPI_DOUBLE_COMPLEX,filetype,"external32",MPI_INFO_NULL,ierr)

       CALL MPI_FILE_WRITE(rst_mpi_fh,b(:,1+rank*nr_per_rank),data_to_write,MPI_DOUBLE_COMPLEX,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,aj(:,1+rank*nr_per_rank),data_to_write,MPI_DOUBLE_COMPLEX,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,dbdtLast(:,1+rank*nr_per_rank),data_to_write,MPI_DOUBLE_COMPLEX,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,djdtLast(:,1+rank*nr_per_rank),data_to_write,MPI_DOUBLE_COMPLEX,status,ierr)

       disp = 8+header_length+16+2*number_of_fields*lm_max*n_r_max*SIZEOF_DOUBLE_COMPLEX&
            & + 4 + 4*lm_max*n_r_max*SIZEOF_DOUBLE_COMPLEX
       CALL MPI_FILE_SET_VIEW(rst_mpi_fh,disp,MPI_BYTE,MPI_BYTE,"external32",MPI_INFO_NULL,ierr)

       IF (rank  ==  n_procs-1) THEN
          CALL MPI_FILE_WRITE(rst_mpi_fh,&
               &4*lm_max*n_r_max*SIZEOF_DOUBLE_COMPLEX,&
               &1,MPI_INTEGER,status,ierr)
       END IF
    ENDIF

    !-- Write IC magnetic field:
    IF ( l_mag .AND. l_cond_ic ) THEN
       size_of_field_ic = lm_maxMag*n_r_ic_max
       IF (rank == n_procs-1) THEN

          ! write the starting record marker for the ...Last fields
          CALL MPI_FILE_WRITE(rst_mpi_fh,&
               &4*size_of_field_ic*SIZEOF_DOUBLE_COMPLEX,&
               &1,MPI_INTEGER,status,ierr)
          !WRITE(n_rst_file) b_ic,aj_ic,dbdt_icLast,djdt_icLast
          ! write the IC fields on the last PE
          CALL MPI_FILE_WRITE(rst_mpi_fh,b_ic, size_of_field_ic,MPI_DOUBLE_COMPLEX,status,ierr)
          CALL MPI_FILE_WRITE(rst_mpi_fh,aj_ic,size_of_field_ic,MPI_DOUBLE_COMPLEX,status,ierr)
          CALL MPI_FILE_WRITE(rst_mpi_fh,dbdt_icLast,size_of_field_ic,MPI_DOUBLE_COMPLEX,status,ierr)
          CALL MPI_FILE_WRITE(rst_mpi_fh,djdt_icLast,size_of_field_ic,MPI_DOUBLE_COMPLEX,status,ierr)

          ! write the starting record marker for the ...Last fields
          CALL MPI_FILE_WRITE(rst_mpi_fh,&
               &4*size_of_field_ic*SIZEOF_DOUBLE_COMPLEX,&
               &1,MPI_INTEGER,status,ierr)
       END IF
    ENDIF

    IF (rank  ==  n_procs-1) THEN
       CALL MPI_FILE_WRITE(rst_mpi_fh,&
            &15*SIZEOF_DOUBLE_PRECISION,&
            &1,MPI_INTEGER,status,ierr)

       CALL MPI_FILE_WRITE(rst_mpi_fh,lorentz_torque_icLast,1,MPI_DOUBLE_PRECISION,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,lorentz_torque_maLast,1,MPI_DOUBLE_PRECISION,status,ierr)
       
       CALL MPI_FILE_WRITE(rst_mpi_fh,omega_ic1,1,MPI_DOUBLE_PRECISION,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,omegaOsz_ic1,1,MPI_DOUBLE_PRECISION,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,tOmega_ic1,1,MPI_DOUBLE_PRECISION,status,ierr)
       
       CALL MPI_FILE_WRITE(rst_mpi_fh,omega_ic2,1,MPI_DOUBLE_PRECISION,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,omegaOsz_ic2,1,MPI_DOUBLE_PRECISION,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,tOmega_ic2,1,MPI_DOUBLE_PRECISION,status,ierr)
       
       CALL MPI_FILE_WRITE(rst_mpi_fh,omega_ma1,1,MPI_DOUBLE_PRECISION,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,omegaOsz_ma1,1,MPI_DOUBLE_PRECISION,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,tOmega_ma1,1,MPI_DOUBLE_PRECISION,status,ierr)
       
       CALL MPI_FILE_WRITE(rst_mpi_fh,omega_ma2,1,MPI_DOUBLE_PRECISION,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,omegaOsz_ma2,1,MPI_DOUBLE_PRECISION,status,ierr)
       CALL MPI_FILE_WRITE(rst_mpi_fh,tOmega_ma2,1,MPI_DOUBLE_PRECISION,status,ierr)
       
       CALL MPI_FILE_WRITE(rst_mpi_fh,dtNew,1,MPI_DOUBLE_PRECISION,status,ierr)
       
       CALL MPI_FILE_WRITE(rst_mpi_fh,&
            &15*SIZEOF_DOUBLE_PRECISION,&
            &1,MPI_INTEGER,status,ierr)
    END IF
    CALL MPI_Type_free(filetype,ierr)

  end SUBROUTINE store_mpi

END MODULE store_rst
