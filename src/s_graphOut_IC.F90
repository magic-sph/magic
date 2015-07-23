!$Id$
!***********************************************************************
SUBROUTINE graphOut_IC(format,b_ic,db_ic,ddb_ic,aj_ic,dj_ic,b)
  !***********************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/17/02  by JW. --------------!

  !  +-------------------------------------------------------------------+
  !  |                                                                   |
  !  |  Purpose of this subroutine is to write inner core magnetic       |
  !  |  field onto graphic output file. If the inner core is             |
  !  |  insulating (l_cond_ic=false) the potential field is calculated   |
  !  |  from the outer core field at r=r_cmb.                            |
  !  |  This version assumes that the fields are fully local on the rank |
  !  |  which is calling this routine (usually rank 0).                  |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+

  USE truncation
  USE radial_functions
  USE num_param
  USE blocking
  USE horizontal_data
  USE logic
  USE output_data
  use parallel_mod
#if (FFTLIB==JW)
  USE fft_JW
#elif (FFTLIB==MKL)
  USE fft_MKL
#endif

  IMPLICIT NONE

  !-- input:
  INTEGER,INTENT(IN) :: FORMAT    ! controls output format

  COMPLEX(kind=8),INTENT(IN) :: b_ic(lm_maxMag,n_r_ic_maxMag)
  COMPLEX(kind=8),INTENT(IN) :: db_ic(lm_maxMag,n_r_ic_maxMag)
  COMPLEX(kind=8),INTENT(IN) :: ddb_ic(lm_maxMag,n_r_ic_maxMag)
  COMPLEX(kind=8),INTENT(IN) :: aj_ic(lm_maxMag,n_r_ic_maxMag)
  COMPLEX(kind=8),INTENT(IN) :: dj_ic(lm_maxMag,n_r_ic_maxMag)
  COMPLEX(kind=8),INTENT(IN) :: b(lm_maxMag,n_r_maxMag)


  !-- local variables:

  INTEGER :: nR
  INTEGER :: nThetaB,nTheta,nThetaStart,nThetaC
  INTEGER :: nPhi

  COMPLEX(kind=8) :: dLhb(lm_max)
  COMPLEX(kind=8) :: bhG(lm_max)
  COMPLEX(kind=8) :: bhC(lm_max)
  COMPLEX(kind=8) :: dLhj(lm_max)
  COMPLEX(kind=8) :: cbhG(lm_max)
  COMPLEX(kind=8) :: cbhC(lm_max)
  REAL(kind=8) :: BrB(nrp,nfs)
  REAL(kind=8) :: BtB(nrp,nfs)
  REAL(kind=8) :: BpB(nrp,nfs)
  REAL(kind=4) :: Br(n_phi_max,n_theta_max)
  REAL(kind=4) :: Bt(n_phi_max,n_theta_max)
  REAL(kind=4) :: Bp(n_phi_max,n_theta_max)


  INTEGER :: precision  ! controls output precision in graph_write

#ifdef WITH_MPI
  ! MPI specific variables
  INTEGER :: status(MPI_STATUS_SIZE)
  ! end MPI variables
#endif

  !-- end of declaration
  !----------------------------------------------------------------------


  precision=1

  !-- Loop over all radial levels:

  DO nR=2,n_r_ic_max  ! nR=1 is ICB

     IF ( l_cond_ic ) THEN
        CALL legPrep_IC(b_ic(1,nR),db_ic(1,nR),ddb_ic(1,nR), &
             aj_ic(1,nR),dj_ic(1,nR), &
             dLh,lm_max,l_max,minc, &
             r_ic(nR),r_ICB,.FALSE.,.TRUE.,l_cond_ic, &
             dLhb,bhG,bhC,dLhj,cbhG,cbhC)
     ELSE
        !CALL legPrep_IC(b(1,n_r_icb),db_ic(1,n_r_icb), &
        !            ddb_ic(1,n_r_icb),aj_ic(1,n_r_icb), &
        !        dj_ic(1,n_r_icb),dLh,lm_max,l_max,minc, &
        !       r_ic(nR),r_ICB,.FALSE.,.TRUE.,l_cond_ic, &
        !                   dLhb,bhG,bhC,dLhj,cbhG,cbhC)
        CALL legPrep_IC(b(1,n_r_icb),db_ic(1,1), &
             ddb_ic(1,1),aj_ic(1,1), &
             dj_ic(1,1),dLh,lm_max,l_max,minc, &
             r_ic(nR),r_ICB,.FALSE.,.TRUE.,l_cond_ic, &
             dLhb,bhG,bhC,dLhj,cbhG,cbhC)
     END IF

     DO nThetaB=1,nThetaBs
        nThetaStart=(nThetaB-1)*sizeThetaB+1

        !------ Preform Legendre transform:
        CALL legTF(dLhb,bhG,bhC,dLhj,cbhG,cbhC, &
             l_max,minc,nThetaStart,sizeThetaB, &
             Plm,dPlm,.TRUE.,.FALSE., &
             BrB,BtB,BpB,BrB,BrB,BrB)

        CALL fft_thetab(BrB,1)
        CALL fft_thetab(BtB,1)
        CALL fft_thetab(BpB,1)

        !------ Copy theta block and calculate real components:
        DO nTheta=1,sizeThetaB
           nThetaC=nThetaStart-1+nTheta
           DO nPhi=1,n_phi_max
              Br(nPhi,nThetaC)=SNGL(BrB(nPhi,nTheta)*O_r_ic2(nR))
              Bt(nPhi,nThetaC)=SNGL(BtB(nPhi,nTheta)*O_r_ic(nR) * &
                   O_sin_theta(nThetaC))
              Bp(nPhi,nThetaC)=SNGL(BpB(nPhi,nTheta)*O_r_ic(nR) * &
                   O_sin_theta(nThetaC))
           END DO
        END DO

     END DO


     !-- Write radius and theta information (all thetas at one go)
     IF ( format == -1 ) WRITE(n_graph_file, &
          &   '(/"--- Radial level IC, radius, i1, i2 ---")')

     IF ( format /= 0 ) then
        WRITE(n_graph_file,'(I4, F9.5, 2I4)') &
             & n_r_max+nR-2,r_ic(nR)/r_cmb, &
             & 1,n_theta_max
     ELSE IF ( format == 0 ) then
        WRITE(n_graph_file) FLOAT(n_r_max+nR-2),SNGL(r_ic(nR)/r_cmb), &
             &              1.E0,FLOAT(n_theta_max)
     END IF


     !-- Write radial magnetic field:
     IF ( format == -1 ) WRITE(n_graph_file,'(/'' Br IC: '')')
     CALL graph_write(n_phi_max,n_theta_max,Br, &
          precision,format,n_graph_file)

     !-- Write latitudinal magnetic field:
     IF ( format == -1 ) WRITE(n_graph_file,'(/'' Bt IC: '')')
     CALL graph_write(n_phi_max,n_theta_max,Bt, &
          precision,format,n_graph_file)

     !-- Write longitudinal magnetic field:
     IF ( format == -1 ) WRITE(n_graph_file,'(/'' Bp IC: '')')
     CALL graph_write(n_phi_max,n_theta_max,Bp, &
          precision,format,n_graph_file)

  END DO  ! Do loop over radial levels nR


  RETURN
end SUBROUTINE graphOut_IC

!-----------------------------------------------------------------------
