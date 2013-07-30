MODULE radial_data
  USE truncation, only: n_r_max
  USE parallel_mod, ONLY: rank, n_procs,nr_per_rank,nr_on_last_rank
  USE logic,ONLY: l_mag, lVerbose
  IMPLICIT NONE

  PRIVATE

  INTEGER,PUBLIC :: nRstart,nRstop,nRstartMag,nRstopMag
  INTEGER,public :: n_r_cmb,n_r_icb

  public :: initialize_radial_data

CONTAINS
  SUBROUTINE initialize_radial_data

    n_r_cmb=1
    n_r_icb=n_r_max

#ifdef WITH_MPI
    nR_per_rank = (n_r_max-1)/n_procs
    nRstart = n_r_cmb + rank*nR_per_rank
    nRstop  = n_r_cmb + (rank+1)*nR_per_rank - 1

    IF (rank.EQ.n_procs-1) THEN
       ! add the last point to the last process, which now has nR_per_rank+1
       ! radial points
       nRstop = nRstop+1
    END IF
    nR_on_last_rank = nR_per_rank+1
#else
    nR_per_rank = n_r_max
    nR_on_last_rank = n_r_max
    nRstart = n_r_cmb
    nRstop  = n_r_icb
#endif
    IF (l_mag) THEN
       nRstartMag = nRstart
       nRstopMag  = nRstop
    ELSE
       nRstartMag = 1
       nRstopMag  = 1
    END IF

    IF (lVerbose) THEN
       WRITE(*,"(4(A,I4))") "On rank ",rank," nR is in (",nRstart,",",nRstop,"), nr_per_rank is ",nr_per_rank
    ENDIF
  END SUBROUTINE initialize_radial_data
END MODULE radial_data

