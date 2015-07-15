!$Id$
!***************************************************************
!  Common blocks containing 'time derivatives' of the fields.
!***************************************************************

MODULE fieldsLast
  use truncation
  USE LMLoop_data,ONLY: llm,ulm,llmMag,ulmMag
  USE parallel_Mod,only: rank
  implicit none

  !--- The following variables labeled Last are provided
  !    by the restart file for the first time step or
  !    calculated here or by the update routines for the
  !    following time step.
  !    These fields remain in the LM-distributed space 

  COMPLEX(kind=8),ALLOCATABLE :: dwdtLast(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dpdtLast(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dwdtLast_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dpdtLast_LMloc(:,:)

  COMPLEX(kind=8),ALLOCATABLE :: dzdtLast(:,:)
  !COMPLEX(kind=8),ALLOCATABLE :: dzdtLast_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dzdtLast_lo(:,:)

  COMPLEX(kind=8),ALLOCATABLE :: dsdtLast(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dsdtLast_LMloc(:,:)

  COMPLEX(kind=8),ALLOCATABLE :: dbdtLast(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: djdtLast(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dbdtLast_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: djdtLast_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dbdt_icLast(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: djdt_icLast(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: dbdt_icLast_LMloc(:,:)
  COMPLEX(kind=8),ALLOCATABLE :: djdt_icLast_LMloc(:,:)

  REAL(kind=8) :: d_omega_ma_dtLast,d_omega_ic_dtLast
  REAL(kind=8) :: lorentz_torque_maLast,lorentz_torque_icLast
CONTAINS
  SUBROUTINE initialize_fieldsLast

    IF (rank == 0) THEN
       ALLOCATE( dwdtLast(lm_max,n_r_max) )
       ALLOCATE( dpdtLast(lm_max,n_r_max) )
       ALLOCATE( dzdtLast(lm_max,n_r_max) )
       ALLOCATE( dsdtLast(lm_max,n_r_max) )
       ALLOCATE( dbdtLast(lm_maxMag,n_r_maxMag) )
       ALLOCATE( djdtLast(lm_maxMag,n_r_maxMag) )
       ALLOCATE( dbdt_icLast(lm_maxMag,n_r_ic_maxMag) )
       ALLOCATE( djdt_icLast(lm_maxMag,n_r_ic_maxMag) )
    ELSE
       ALLOCATE( dwdtLast(1,n_r_max) )
       ALLOCATE( dpdtLast(1,n_r_max) )
       ALLOCATE( dzdtLast(1,n_r_max) )
       ALLOCATE( dsdtLast(1,n_r_max) )
       ALLOCATE( dbdtLast(1,n_r_max) )
       ALLOCATE( djdtLast(1,n_r_max) )
       ALLOCATE( dbdt_icLast(1,n_r_max) )
       ALLOCATE( djdt_icLast(1,n_r_max) )
    END IF
    ALLOCATE( dwdtLast_LMloc(llm:ulm,n_r_max) )
    ALLOCATE( dpdtLast_LMloc(llm:ulm,n_r_max) )
    !ALLOCATE( dzdtLast_LMloc(llm:ulm,n_r_max) )
    ALLOCATE( dzdtLast_lo(llm:ulm,n_r_max) )
    ALLOCATE( dsdtLast_LMloc(llm:ulm,n_r_max) )

    ALLOCATE( dbdtLast_LMloc(llmMag:ulmMag,n_r_maxMag) )
    ALLOCATE( djdtLast_LMloc(llmMag:ulmMag,n_r_maxMag) )
    ALLOCATE( dbdt_icLast_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
    ALLOCATE( djdt_icLast_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
    
  END SUBROUTINE initialize_fieldsLast
END MODULE fieldsLast
