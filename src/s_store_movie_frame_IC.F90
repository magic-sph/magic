!$Id$
!*************************************************************************
SUBROUTINE store_movie_frame_IC(b,b_ic,db_ic,ddb_ic,aj_ic,dj_ic)
  !*************************************************************************

  !  +-------------+----------------+------------------------------------+
  !  |                                                                   |
  !  |  Controls storage of IC magnetic field in movie frame.            |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+

  USE truncation
  USE radial_functions
  USE physical_parameters
  USE blocking
  USE horizontal_data
  USE logic
  USE movie_data
  USE output_data
#if (FFTLIB==JW)
  USE fft_JW
#elif (FFTLIB==MKL)
  USE fft_MKL
#endif
  IMPLICIT NONE

  !-- Input of scalar fields:
  COMPLEX(kind=8),INTENT(IN) :: b(lm_maxMag,n_r_maxMag)
  COMPLEX(kind=8),INTENT(IN) :: b_ic(lm_maxMag,n_r_ic_maxMag)
  COMPLEX(kind=8),INTENT(IN) :: db_ic(lm_maxMag,n_r_ic_maxMag)
  COMPLEX(kind=8),INTENT(IN) :: ddb_ic(lm_maxMag,n_r_ic_maxMag)
  COMPLEX(kind=8),INTENT(IN) :: aj_ic(lm_maxMag,n_r_ic_maxMag)
  COMPLEX(kind=8),INTENT(IN) :: dj_ic(lm_maxMag,n_r_ic_maxMag)

  !-- Output: frames stored in array frames(*)

  !-- Local:
  INTEGER :: n_movie        ! No. of movie
  INTEGER :: n_field        ! No. of field
  INTEGER :: n_type         ! Movie type
  INTEGER :: n_surface      ! Surface (1=r,2=theta,3=phi)
  INTEGER :: n_const        ! Gives surface
  INTEGER :: nR
  INTEGER :: n_field_type   ! Numbers field types
  INTEGER :: n_store_last   ! Position i in frame(i) were field starts
  INTEGER :: nTheta,nThetaB,nThetaR,nThetaC,nThetaStart
  INTEGER :: nPhi,nPhi0,nPhi180
  INTEGER :: n_field_size
  INTEGER :: n_fields
  INTEGER :: n_fields_oc
  INTEGER :: n_fields_ic
  INTEGER :: n_o,n_o_r

  COMPLEX(kind=8) :: dLhb(lm_maxMag)
  COMPLEX(kind=8) :: bhG(lm_maxMag)
  COMPLEX(kind=8) :: bhC(lm_maxMag)
  COMPLEX(kind=8) :: dLhj(lm_maxMag)
  COMPLEX(kind=8) :: cbhG(lm_maxMag)
  COMPLEX(kind=8) :: cbhC(lm_maxMag)
  REAL(kind=8) :: BrB(nrp,nfs)
  REAL(kind=8) :: BtB(nrp,nfs)
  REAL(kind=8) :: BpB(nrp,nfs)
  REAL(kind=8) :: cBrB(nrp,nfs)
  REAL(kind=8) :: cBtB(nrp,nfs)
  REAL(kind=8) :: cBpB(nrp,nfs)
  REAL(kind=8) :: fl(nfs),help

  REAL(kind=8) ::  phi_norm

  !-- end of declaration
  !----------------------------------------------------------------------


  phi_norm=1.D0/n_phi_max ! 2 pi /n_phi_max

  DO n_movie=1,n_movies

     n_fields_oc=n_movie_fields(n_movie)
     n_fields_ic=n_movie_fields_ic(n_movie)
     n_fields   =n_fields_oc+n_fields_ic
     n_surface  =n_movie_surface(n_movie)
     n_const    =n_movie_const(n_movie)
     n_type     =n_movie_type(n_movie)

     IF ( n_surface == 0 ) THEN

        DO nR=2,n_r_ic_max-1

           IF ( l_cond_ic ) THEN
              CALL legPrep_IC(b_ic(1,nR),db_ic(1,nR),ddb_ic(1,nR), &
                   aj_ic(1,nR),dj_ic(1,nR), &
                   dLh,lm_max,l_max,minc, &
                   r_ic(nR),r_ICB,.TRUE.,.TRUE.,l_cond_ic, &
                   dLhb,bhG,bhC,dLhj,cbhG,cbhC)
           ELSE
              CALL legPrep_IC(b(1,n_r_icb),db_ic(1,1),ddb_ic(1,1), &
                   aj_ic(1,1),dj_ic(1,1), &
                   dLh,lm_max,l_max,minc, &
                   r_ic(nR),r_ICB,.TRUE.,.TRUE.,l_cond_ic, &
                   dLhb,bhG,bhC,dLhj,cbhG,cbhC)
           END IF

           !------ Calculate magnetic field on grid points:
           DO nThetaB=1,nThetaBs
              nThetaStart=(nThetaB-1)*sizeThetaB+1

              !------ Preform Legendre transform:
              CALL legTF(dLhb,bhG,bhC,dLhj,cbhG,cbhC, &
                   l_max,minc,nThetaStart,sizeThetaB, &
                   Plm,dPlm,lm_max,ncp,.TRUE.,.TRUE., &
                   BrB,BtB,BpB,cBrB,cBtB,cBpB)
              CALL fft_thetab(BrB,1)
              CALL fft_thetab(BtB,1)
              CALL fft_thetab(BpB,1)
              CALL fft_thetab(cBrB,1)
              CALL fft_thetab(cBtB,1)

              DO n_field=n_fields_oc+1,n_fields
                 n_field_type= &
                      n_movie_field_type(n_field,n_movie)
                 n_store_last= &
                      n_movie_field_start(n_field,n_movie)-1

                 IF ( n_store_last >= 0 ) THEN
                    n_o_r=n_store_last + &
                         (nR-2)*n_theta_max*n_phi_max

                    IF ( n_field_type == 1 ) THEN
                       DO nThetaR=1,sizeThetaB
                          nThetaC=nThetaStart-1+nThetaR
                          nTheta=n_theta_cal2ord(nThetaC)
                          n_o=n_o_r+(nTheta-1)*n_phi_max
                          DO nPhi=1,n_phi_max
                             frames(nPhi+n_o)=BrB(nPhi,nThetaR) * &
                                   O_r_ic2(nR)
                          END DO
                       END DO
                    ELSE IF ( n_field_type == 2 ) THEN
                       DO nThetaR=1,sizeThetaB
                          nThetaC=nThetaStart-1+nThetaR
                          nTheta=n_theta_cal2ord(nThetaC)
                          n_o=n_o_r+(nTheta-1)*n_phi_max
                          help=O_r_ic(nR)*O_sin_theta(nThetaC)
                          DO nPhi=1,n_phi_max
                             frames(nPhi+n_o)=help*BtB(nPhi,nThetaR)
                          END DO
                       END DO
                    ELSE IF ( n_field_type == 3 ) THEN
                       DO nThetaR=1,sizeThetaB
                          nThetaC=nThetaStart-1+nThetaR
                          nTheta=n_theta_cal2ord(nThetaC)
                          n_o=n_o_r+(nTheta-1)*n_phi_max
                          help=O_r_ic(nR)*O_sin_theta(nThetaC)
                          DO nPhi=1,n_phi_max
                             frames(nPhi+n_o)=help*BpB(nPhi,nThetaR)
                          END DO
                       END DO
                    ELSE IF ( n_field_type == 54 ) THEN
                       help=LFfac*O_r_ic(nR)*O_r_ic2(nR)
                       DO nThetaR=1,sizeThetaB
                          nThetaC=nThetaStart-1+nThetaR
                          nTheta=n_theta_cal2ord(nThetaC)
                          n_o=n_o_r+(nTheta-1)*n_phi_max
                          DO nPhi=1,n_phi_max
                             frames(nPhi+n_o)= &
                                  help*O_sin_theta(nThetaC) * &
                                  ( cBrB(nPhi,nThetaR)*BtB(nPhi,nThetaR) - &
                                  cBtB(nPhi,nThetaR)*BrB(nPhi,nThetaR) )
                          END DO
                       END DO
                    END IF

                 END IF

              END DO  ! Loop over fields

           END DO     ! Loop over thetas blocks

        END DO        ! Loop over radius

     ELSE IF ( n_surface == 1 ) THEN

        nR=iabs(n_const)

        IF ( l_cond_ic ) THEN
           CALL legPrep_IC(b_ic(1,nR),db_ic(1,nR),ddb_ic(1,nR), &
                aj_ic(1,nR),dj_ic(1,nR), &
                dLh,lm_max,l_max,minc, &
                r_ic(nR),r_ICB,.TRUE.,.TRUE.,l_cond_ic, &
                dLhb,bhG,bhC,dLhj,cbhG,cbhC)
        ELSE
           CALL legPrep_IC(b(1,n_r_icb),db_ic(1,1),ddb_ic(1,1), &
                aj_ic(1,1),dj_ic(1,1), &
                dLh,lm_max,l_max,minc, &
                r_ic(nR),r_ICB,.TRUE.,.TRUE.,l_cond_ic, &
                dLhb,bhG,bhC,dLhj,cbhG,cbhC)
        END IF

        !------ Calculate magnetic field on grid points:
        DO nThetaB=1,nThetaBs
           nThetaStart=(nThetaB-1)*sizeThetaB+1

           !------ Preform Legendre transform:
           CALL legTF(dLhb,bhG,bhC,dLhj,cbhG,cbhC, &
                l_max,minc,nThetaStart,sizeThetaB, &
                Plm,dPlm,lm_max,ncp,.TRUE.,.TRUE., &
                BrB,BtB,BpB,cBrB,cBtB,cBpB)
           CALL fft_thetab(BrB,1)
           CALL fft_thetab(BtB,1)
           CALL fft_thetab(BpB,1)

           DO n_field=n_fields_oc+1,n_fields
              n_field_type= &
                   n_movie_field_type(n_field,n_movie)
              n_store_last= &
                   n_movie_field_start(n_field,n_movie)-1

              IF ( n_store_last >= 0 ) THEN
                 n_o_r=n_store_last

                 IF ( n_field_type == 1 ) THEN
                    DO nThetaR=1,sizeThetaB
                       nThetaC=nThetaStart-1+nThetaR
                       nTheta=n_theta_cal2ord(nThetaC)
                       n_o=n_o_r+(nTheta-1)*n_phi_max
                       DO nPhi=1,n_phi_max
                          frames(nPhi+n_o)=BrB(nPhi,nThetaR) * &
                               O_r_ic2(nR)
                       END DO
                    END DO
                 ELSE IF ( n_field_type == 2 ) THEN
                    DO nThetaR=1,sizeThetaB
                       nThetaC=nThetaStart-1+nThetaR
                       nTheta=n_theta_cal2ord(nThetaC)
                       n_o=n_o_r+(nTheta-1)*n_phi_max
                       help=O_r_ic(nR)*O_sin_theta(nThetaC)
                       DO nPhi=1,n_phi_max
                          frames(nPhi+n_o)=help*BtB(nPhi,nThetaR)
                       END DO
                    END DO
                 ELSE IF ( n_field_type == 3 ) THEN
                    DO nThetaR=1,sizeThetaB
                       nThetaC=nThetaStart-1+nThetaR
                       nTheta=n_theta_cal2ord(nThetaC)
                       n_o=n_o_r+(nTheta-1)*n_phi_max
                       help=O_r_ic(nR)*O_sin_theta(nThetaC)
                       DO nPhi=1,n_phi_max
                          frames(nPhi+n_o)=help*BpB(nPhi,nThetaR)
                       END DO
                    END DO
                 END IF
              END IF

           END DO  ! Loop over fields

        END DO     ! Loop over theta blocks

     ELSE IF ( n_surface == 2 ) THEN

        DO nR=2,n_r_ic_max-1

           IF ( l_cond_ic ) THEN
              CALL legPrep_IC(b_ic(1,nR),db_ic(1,nR),ddb_ic(1,nR), &
                   aj_ic(1,nR),dj_ic(1,nR), &
                   dLh,lm_max,l_max,minc, &
                   r_ic(nR),r_ICB,.TRUE.,.TRUE.,l_cond_ic, &
                   dLhb,bhG,bhC,dLhj,cbhG,cbhC)
           ELSE
              CALL legPrep_IC(b(1,n_r_icb),db_ic(1,1),ddb_ic(1,1), &
                   aj_ic(1,1),dj_ic(1,1), &
                   dLh,lm_max,l_max,minc, &
                   r_ic(nR),r_ICB,.TRUE.,.TRUE.,l_cond_ic, &
                   dLhb,bhG,bhC,dLhj,cbhG,cbhC)
           END IF

           IF ( MOD(n_const,2) == 1 ) THEN
              nTheta=n_const
              nThetaR=1
           ELSE
              nTheta=n_const-1
              nThetaR=2
           END IF

           !------ Preform Legendre transform for 2 theta points
           CALL legTF(dLhb,bhG,bhC,dLhj,cbhG,cbhC, &
                l_max,minc,nTheta,2, &
                Plm,dPlm,lm_max,ncp,.TRUE.,.TRUE., &
                BrB,BtB,BpB,cBrB,cBtB,cBpB)
           CALL fft_thetab(BrB,1)
           CALL fft_thetab(BtB,1)
           CALL fft_thetab(BpB,1)
           CALL fft_thetab(cBtB,1)

           DO n_field=n_fields_oc+1,n_fields

              n_field_type= &
                   n_movie_field_type(n_field,n_movie)
              n_store_last= &
                   n_movie_field_start(n_field,n_movie)-1

              IF ( n_store_last >= 0 ) THEN
                 n_o=n_store_last+(nR-2)*n_phi_max
                 IF ( n_field_type == 1 ) THEN
                    DO nPhi=1,n_phi_max
                       frames(nPhi+n_o)=BrB(nPhi,nThetaR) * &
                            O_r_ic2(nR)
                    END DO
                 ELSE IF ( n_field_type == 2 ) THEN
                    help=O_r_ic(nR)*O_sin_theta(nTheta)
                    DO nPhi=1,n_phi_max
                       frames(nPhi+n_o)=help*BtB(nPhi,nThetaR)
                    END DO
                 ELSE IF ( n_field_type == 3 ) THEN
                    help=O_r_ic(nR)*O_sin_theta(nTheta)
                    DO nPhi=1,n_phi_max
                       frames(nPhi+n_o)=help*BpB(nPhi,nThetaR)
                    END DO
                 ELSE IF ( n_field_type == 13 ) THEN
                    help=-O_r_ic(nR)*O_sin_theta(nTheta)
                    DO nPhi=1,n_phi_max
                       frames(nPhi+n_o)=help*BtB(nPhi,nThetaR)
                    END DO
                 ELSE IF ( n_field_type == 14 ) THEN
                    help=-O_r_ic(nR)*O_sin_theta(nTheta)
                    DO nPhi=1,n_phi_max
                       frames(nPhi+n_o)=help*cBtB(nPhi,nThetaR)
                    END DO
                 END IF

              END IF

           END DO  ! Loop over fields

        END DO     ! Loop over radius


     ELSE IF ( IABS(n_surface) == 3 ) THEN

        DO nR=2,n_r_ic_max-1

           IF ( l_cond_ic ) THEN
              CALL legPrep_IC(b_ic(1,nR),db_ic(1,nR),ddb_ic(1,nR), &
                   aj_ic(1,nR),dj_ic(1,nR), &
                   dLh,lm_max,l_max,minc, &
                   r_ic(nR),r_ICB,.TRUE.,.TRUE.,l_cond_ic, &
                   dLhb,bhG,bhC,dLhj,cbhG,cbhC)
           ELSE
              CALL legPrep_IC(b(1,n_r_icb),db_ic(1,1),ddb_ic(1,1), &
                   aj_ic(1,1),dj_ic(1,1), &
                   dLh,lm_max,l_max,minc, &
                   r_ic(nR),r_ICB,.TRUE.,.TRUE.,l_cond_ic, &
                   dLhb,bhG,bhC,dLhj,cbhG,cbhC)
           END IF

           !------ Get phi no. for left and righty halfspheres:
           nPhi0=n_const
           IF ( MOD(minc,2) == 1 ) THEN
              nPhi180=n_phi_max/2+nPhi0
           ELSE
              nPhi180=nPhi0
           END IF

           !------ Calculate magnetic field on grid points:
           DO nThetaB=1,nThetaBs
              nThetaStart=(nThetaB-1)*sizeThetaB+1

              IF ( n_type == 30 ) THEN
                 !------ get_fl returns field for field line plot:
                 CALL get_fl(fl,nR,nThetaStart,sizeThetaB,.TRUE.)
              ELSE
                 !------ Preform Legendre transform:
                 CALL legTF(dLhb,bhG,bhC,dLhj,cbhG,cbhC, &
                      l_max,minc,nThetaStart,sizeThetaB, &
                      Plm,dPlm,lm_max,ncp,.TRUE.,.TRUE., &
                      BrB,BtB,BpB,cBrB,cBtB,cBpB)
                 CALL fft_thetab(BrB,1)
                 CALL fft_thetab(BtB,1)
                 CALL fft_thetab(BpB,1)
                 CALL fft_thetab(cBrB,1)
                 CALL fft_thetab(cBtB,1)
              END IF

              DO n_field=n_fields_oc+1,n_fields

                 n_field_type= &
                      n_movie_field_type(n_field,n_movie)
                 n_store_last= &
                      n_movie_field_start(n_field,n_movie)-1
                 n_field_size = &
                      ( n_movie_field_stop(n_field,n_movie) - &
                      n_store_last )/2

                 IF ( n_store_last >= 0 ) THEN
                    n_o_r=n_store_last+(nR-2)*n_theta_max

                    IF ( n_field_type == 1 ) THEN
                       DO nThetaR=1,sizeThetaB
                          nThetaC=nThetaStart-1+nThetaR
                          nTheta=n_theta_cal2ord(nThetaC)
                          n_o=n_o_r+nTheta
                          frames(n_o)=BrB(nPhi0,nThetaR) * &
                               O_r_ic2(nR)
                          frames(n_o+n_field_size)= &
                               O_r_ic2(nR)*BrB(nPhi180,nThetaR)
                       END DO
                    ELSE IF ( n_field_type == 2 ) THEN
                       DO nThetaR=1,sizeThetaB
                          nThetaC=nThetaStart-1+nThetaR
                          nTheta=n_theta_cal2ord(nThetaC)
                          n_o=n_o_r+nTheta
                          frames(n_o)=BtB(nPhi0,nThetaR) * &
                               O_r_ic(nR)*O_sin_theta(nThetaC)
                          frames(n_o+n_field_size)= &
                              BtB(nPhi180,nThetaR) * &
                               O_r_ic(nR)*O_sin_theta(nThetaC)
                      END DO
                    ELSE IF ( n_field_type == 3 ) THEN
                        DO nThetaR=1,sizeThetaB
                           nThetaC=nThetaStart-1+nThetaR
                           nTheta=n_theta_cal2ord(nThetaC)
                           n_o=n_o_r+nTheta
                           frames(n_o)=BpB(nPhi0,nThetaR) * &
                               O_r_ic(nR)*O_sin_theta(nThetaC)
                           frames(n_o+n_field_size)= &
                               BpB(nPhi180,nThetaR) * &
                               O_r_ic(nR)*O_sin_theta(nThetaC)
                       END DO
                    ELSE IF ( n_field_type == 8 ) THEN
                       DO nThetaR=1,sizeThetaB
                          nThetaC=nThetaStart-1+nThetaR
                          nTheta=n_theta_cal2ord(nThetaC)
                          n_o=n_o_r+nTheta
                          !WRITE(*,"(A,I5,A)") "store_movie_IC: frames(",n_o,")"
                          frames(n_o)=fl(nThetaR)
                       END DO
                    ELSE IF ( n_field_type == 9 ) then
                       DO nThetaR=1,sizeThetaB
                          nThetaC=nThetaStart-1+nThetaR
                          nTheta=n_theta_cal2ord(nThetaC)
                          n_o=n_o_r+nTheta
                          help=0.D0
                          DO nPhi=1,n_phi_max
                             help=help+BpB(nPhi,nThetaR)
                          END DO
                          frames(n_o)=phi_norm*help * &
                               O_r_ic(nR)*O_sin_theta(nThetaC)
                       END DO
                    ELSE IF ( n_field_type == 54 ) THEN
                       help=LFfac*O_r_ic(nR)*O_r_ic2(nR)
                       DO nThetaR=1,sizeThetaB
                          nThetaC=nThetaStart-1+nThetaR
                          nTheta=n_theta_cal2ord(nThetaC)
                          n_o=n_o_r+nTheta
                          frames(n_o)=help*O_sin_theta(nThetaC) * &
                               ( cBrB(nPhi0,nThetaR)*BtB(nPhi0,nThetaR) - &
                               cBtB(nPhi0,nThetaR)*BrB(nPhi0,nThetaR) )
                          frames(n_o)=help*O_sin_theta(nThetaC) * &
                               ( cBrB(nPhi180,nThetaR)*BtB(nPhi180,nThetaR) - &
                              cBtB(nPhi180,nThetaR)*BrB(nPhi180,nThetaR) )
                       END DO
                    ELSE
                       DO nThetaR=1,sizeThetaB
                          nThetaC=nThetaStart-1+nThetaR
                          nTheta=n_theta_cal2ord(nThetaC)
                          n_o=n_o_r+nTheta
                          frames(n_o)=0.D0
                          frames(n_o+n_field_size)=0.D0
                       END DO
                    END IF
                 END IF

              END DO    ! Loop over fields

           END DO       ! Loop over theta blocks

        END DO          ! Loop over r

     END IF  ! Which surface ?

  END DO     ! Loop over movies


  RETURN
end SUBROUTINE store_movie_frame_IC


!--- End of subroutine  store_movie_frame_ic
!-----------------------------------------------------------------------
