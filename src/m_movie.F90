!$Id$
!********************************************************************
!  Module containing work arrays
!********************************************************************

!------------ This is release 1 level 1  --------------!
!------------ Created on 1/17/02  by JW. --------------!

MODULE movie_data
  use truncation

  IMPLICIT NONE

  REAL(kind=8) :: movieDipColat,movieDipLon
  REAL(kind=8) :: movieDipStrength,movieDipStrengthGeo
  REAL(kind=8) :: t_movieS(10000)

  !-- Work arrays for storing movie frame:
  INTEGER :: n_frame_work  
  INTEGER :: n_MD
  REAL(kind=8),ALLOCATABLE :: frames(:)
CONTAINS
  SUBROUTINE initialize_movie_data

    n_MD=lMovieMem*n_movie_work*n_theta_max*n_phi_max*n_r_tot
    n_frame_work=MAX0(n_MD,1)

    ALLOCATE( frames(n_frame_work) )

  END SUBROUTINE initialize_movie_data
END MODULE movie_data
