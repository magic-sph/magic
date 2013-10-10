!$Id$
MODULE getDlm_mod
  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE blocking
  USE horizontal_data
  USE const
  USE LMLoop_data,ONLY: llm,ulm
  USE usefull, ONLY: cc2real,cc22real
  USE integration, ONLY: rInt_R
  
  IMPLICIT NONE

CONTAINS

  SUBROUTINE getDlm(w,dw,z,dl,dlR,dm,dlc,dlRc,switch)

    !    !------------ This is release 2 level 1  --------------!
    !    !------------ Created on 1/17/02  by JW. --------------!

    !--------------------------------------------------------------------

    !  calculates energy  = 1/2 Integral(B^2 dV)
    !  integration in theta,phi by summation over harmonic coeffs.
    !  integration in r by Chebycheff integrals

    !  Output:
    !  enbp: Total poloidal        enbt: Total toroidal
    !  apome: Axisym. poloidal     atome: Axisym. toroidal

    !--------------------------------------------------------------------


    COMPLEX(kind=8),INTENT(IN) :: w(llm:ulm,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: dw(llm:ulm,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: z(llm:ulm,n_r_max)
    CHARACTER(len=1),INTENT(IN) :: switch

    !-- Output:
    REAL(kind=8),INTENT(OUT) :: dlR(n_r_max),dlRc(n_r_max)
    REAL(kind=8),INTENT(OUT) :: dl,dlc,dm

    !-- local:
    INTEGER :: nR,lm,l,m,lFirst
    REAL(kind=8) :: e_p,e_t,e_m,e_l
    REAL(kind=8) :: fac
    REAL(kind=8) :: e_lr(n_r_max,l_max),e_lr_c(n_r_max,l_max)
    REAL(kind=8) :: e_lr_global(n_r_max,l_max),e_lr_c_global(n_r_max,l_max)
    REAL(kind=8) :: e_mr(n_r_max,0:l_max)
    REAL(kind=8) :: e_mr_global(n_r_max,0:l_max)
    REAL(kind=8) :: ER(n_r_max),ELR(n_r_max)
    REAL(kind=8) :: E,EL,EM
    REAL(kind=8) :: ERc(n_r_max),ELRc(n_r_max)
    REAL(kind=8) :: Ec,ELc
    REAL(kind=8) :: O_rho ! 1/rho (anelastic)

    !-- end of declaration
    !---------------------------------------------------------------------


    IF ( switch == 'B' ) THEN
       DO nR=1,n_r_max
          DO l=1,l_max
             e_lr(nR,l)=0.D0
             e_lr_c(nR,l)=0.D0
             e_mr(nR,l)=0.D0
          END DO
          DO lm=MAX(2,llm),ulm
             l =lo_map%lm2l(lm)
             m =lo_map%lm2m(lm)

             e_p= dLh(st_map%lm2(l,m)) *  ( &
                  dLh(st_map%lm2(l,m))*or2(nR)*cc2real(w(lm,nR),m) &
                  & + cc2real(dw(lm,nR),m) )
             e_t=dLh(st_map%lm2(l,m))*cc2real(z(lm,nR),m)

             e_lr(nR,l)=e_lr(nR,l) + e_p+e_t
             e_lr_c(nR,l)=0.D0
             e_mr(nR,m)=e_mr(nR,m) + e_p+e_t
          END DO ! do loop over lms in block
          ! We have now a local sum over the local lm in
          ! e_lr(nR,l), e_mr(nR,m)
       END DO    ! radial grid points

       lFirst=2
    ELSE IF ( switch == 'V' ) THEN
       DO nR=1,n_r_max
          O_rho =orho1(nR)
          DO l=1,l_max
             e_lr(nR,l)=0.D0
             e_lr_c(nR,l)=0.D0
             e_mr(nR,l)=0.D0
          END DO
          DO lm=MAX(2,llm),ulm
             l =lo_map%lm2l(lm)
             m =lo_map%lm2m(lm)

             e_p= O_rho * dLh(st_map%lm2(l,m)) *  ( &
                  dLh(st_map%lm2(l,m))*or2(nR)*cc2real(w(lm,nR),m) &
                  & + cc2real(dw(lm,nR),m) )
             e_t=O_rho*dLh(st_map%lm2(l,m))*cc2real(z(lm,nR),m)
             IF ( m /= 0 ) THEN
                e_lr_c(nR,l)=e_lr_c(nR,l) + e_p+e_t
             end if
             e_lr(nR,l)=e_lr(nR,l) + e_p+e_t
             e_mr(nR,m)=e_mr(nR,m) + e_p+e_t
             !IF (l.EQ.2) THEN
             !   WRITE(*,"(A,3I4,4ES20.12)") "e_lr,e_mr,e_p,e_t = ",nR,l,m,&
             !        &e_lr(nR,l),e_mr(nR,m),&
             !        &e_p,e_t
             !END IF
          END DO ! do loop over lms in block
       END DO    ! radial grid points
       lFirst=1
    ELSE
       WRITE(*,*) 'Wrong switch in s_getDlm.f'
       STOP
    END IF

    ! reduce to rank 0
    CALL MPI_Reduce(e_lr,e_lr_global,n_r_max*l_max,MPI_DOUBLE_PRECISION,&
         &MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(e_mr,e_mr_global,n_r_max*(l_max+1),MPI_DOUBLE_PRECISION,&
         &MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(e_lr_c,e_lr_c_global,n_r_max*l_max,MPI_DOUBLE_PRECISION,&
         &MPI_SUM,0,MPI_COMM_WORLD,ierr)
       
    IF (rank.EQ.0) THEN
       !-- Radial Integrals:
       fac=0.5D0*eScale
       E  =0.D0
       EL =0.D0
       Ec =0.D0
       ELc=0.D0

       DO l=lFirst,l_max
          e_l=0.d0
          e_l=fac*rInt_R(e_lr_global(1,l),n_r_max,n_r_max,drx, &
               &         i_costf_init,d_costf_init)
          !WRITE(*,"(A,I5,ES20.12)") "getDlm: l,e_l = ",l,e_l
          E =E+e_l
          EL=EL+DBLE(l)*e_l
          e_l=0.d0
          e_l=fac*rInt_R(e_lr_c_global(1,l),n_r_max,n_r_max,drx, &
               i_costf_init,d_costf_init)
          Ec =Ec+e_l
          ELc=ELc+DBLE(l)*e_l
       END DO
       IF ( EL /= 0d0 ) THEN
          !WRITE(*,"(A,2ES20.12)") "getDlm: E,EL = ",E,EL
          dl=pi*E/EL
       ELSE
          dl=0d0
       END IF
       IF ( switch == 'V' ) THEN
          IF ( ELc /= 0d0 ) THEN
             dlc=pi*Ec/ELc
          ELSE
             dlc=0d0
          END IF
       ELSE IF ( switch == 'B' ) THEN
          dlc=0.d0
       END if
       DO nR=1,n_r_max
          ER(nR)  =0.D0
          ELR(nR) =0.D0
          ERc(nR) =0.D0
          ELRc(nR)=0.D0
          DO l=lFirst,l_max
             e_l=fac*e_lr_global(nR,l)
             ER(nR) =ER(nR)+e_l
             ELR(nR)=ELR(nR)+DBLE(l)*e_l
             IF ( switch == 'V' ) THEN
                e_l=fac*e_lr_c_global(nR,l)
                ERc(nR) =ERc(nR)+e_l
                ELRc(nR)=ELRc(nR)+DBLE(l)*e_l
             END IF
          END DO
          IF ( switch == 'V' ) THEN
             IF ( ELR(nR) /= 0d0 ) THEN
                dlR(nR)=pi*ER(nR)/ELR(nR)
             ELSE
                dlR(nR)=0.d0
             END IF
             IF ( ELRc(nR) /= 0d0 ) THEN
                dlRc(nR)=pi*ERc(nR)/ELRc(nR)
             ELSE
                dlRc(nR)=0.d0
             END IF
          ELSE IF ( switch == 'B' ) THEN
             dlR(nR)=0.D0
             dlRc(nR)=0.D0
          END IF
       END DO
       E =0.D0
       EM=0.D0
       DO m=minc,m_max,minc
          e_m=fac*rInt_R(e_mr_global(1,m),n_r_max,n_r_max,drx, &
               i_costf_init,d_costf_init)
          E =E +e_m
          EM=EM+DBLE(m)*e_m
       END DO
       IF ( EM /= 0d0 ) THEN
          dm=pi*E/EM
       ELSE
          dm=0d0
       END IF
    END IF

    RETURN
  end SUBROUTINE getDlm

END MODULE getDlm_mod
