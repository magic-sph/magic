!$Id$
!****************************************************************
!  Common block containing information on start solution
!****************************************************************

!   !------------ This is release 1 level 1  --------------!
!   !------------ Created on 1/17/02  by JW. --------------!

MODULE init_fields
  use truncation

  IMPLICIT NONE

  !-- Initialisation of fields:
  INTEGER :: init_s1,init_s2
  INTEGER :: init_b1,init_v1
  INTEGER :: imagcon
  REAL(kind=8) :: tmagcon

  !----- Entropy amplitudes for initialisation:
  REAL(kind=8) :: amp_s1,amp_s2,amp_v1,amp_b1

  !----- Entropy at CMB and ICB (input):
  INTEGER,PARAMETER :: n_s_bounds=20
  REAL(kind=8) :: s_bot(4*n_s_bounds)  ! input variables for tops,bots
  REAL(kind=8) :: s_top(4*n_s_bounds)
  COMPLEX(kind=8),allocatable :: tops(:,:)
  COMPLEX(kind=8),allocatable :: bots(:,:)

  !----- Peak values for magnetic field:
  REAL(kind=8) :: bpeakbot,bpeaktop

  !----- Initialised IC and mantle rotation rates:
  INTEGER :: nRotMa,nRotIc
  REAL(kind=8) :: omega_ma1,omegaOsz_ma1,tShift_ma1,tOmega_ma1
  REAL(kind=8) :: omega_ma2,omegaOsz_ma2,tShift_ma2,tOmega_ma2
  REAL(kind=8) :: omega_ic1,omegaOsz_ic1,tShift_ic1,tOmega_ic1
  REAL(kind=8) :: omega_ic2,omegaOsz_ic2,tShift_ic2,tOmega_ic2

  !----- About start-file:
  LOGICAL :: l_start_file     ! taking fields from startfile ?
  LOGICAL :: l_reset_t        ! reset time from startfile ?
  INTEGER :: inform           ! format of start_file
  INTEGER :: n_start_file     ! I/O unit of start_file
  CHARACTER(len=72) :: start_file  ! name of start_file           

  !-- Scales for input field:
  REAL(kind=8) :: scale_s
  REAL(kind=8) :: scale_v
  REAL(kind=8) :: scale_b
  REAL(kind=8) :: tipdipole       ! adding to symetric field

  !-----------------------------------------------------------------------------------
CONTAINS
    SUBROUTINE initialize_init_fields

      allocate( tops(0:l_max,0:m_max) )
      allocate( bots(0:l_max,0:m_max) )

    END SUBROUTINE initialize_init_fields
END MODULE init_fields
