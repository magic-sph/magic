!$Id$
!***************************************************************
!  Common blocks containing 'time derivatives' of the fields.
!***************************************************************

MODULE fieldsLast
  use truncation
  implicit none

  !--- The following variables labeled Last are provided
  !    by the restart file for the first time step or
  !    calculated here or by the update routines for the
  !    following time step.
  !    These fields remain in the LM-distributed space 

  COMPLEX(kind=8),ALLOCATABLE :: dwdtLast(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dzdtLast(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dpdtLast(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dsdtLast(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dbdtLast(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: djdtLast(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dbdt_icLast(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: djdt_icLast(:,:)
  !COMMON/dt_fieldLast/dwdtLast,dzdtLast,dpdtLast,                 &
  !     &                      dsdtLast,dbdtLast,djdtLast,                 &
  !     &                      dbdt_icLast,djdt_icLast

  REAL(kind=8) :: d_omega_ma_dtLast,d_omega_ic_dtLast
  REAL(kind=8) :: lorentz_torque_maLast,lorentz_torque_icLast
  !COMMON/dt_field_rotLast/d_omega_ma_dtLast,d_omega_ic_dtLast,    &
  !     &           lorentz_torque_maLast,lorentz_torque_icLast
CONTAINS
  SUBROUTINE initialize_fieldsLast

  ALLOCATE( dwdtLast(lm_max,n_r_max) )
  ALLOCATE( dzdtLast(lm_max,n_r_max) )
  ALLOCATE( dpdtLast(lm_max,n_r_max) )
  ALLOCATE( dsdtLast(lm_max,n_r_max) )
  ALLOCATE( dbdtLast(lm_maxMag,n_r_maxMag) )
  ALLOCATE( djdtLast(lm_maxMag,n_r_maxMag) )
  ALLOCATE( dbdt_icLast(lm_maxMag,n_r_ic_maxMag) )
  ALLOCATE( djdt_icLast(lm_maxMag,n_r_ic_maxMag) )
    
  END SUBROUTINE initialize_fieldsLast
END MODULE fieldsLast
