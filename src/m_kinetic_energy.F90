!$Id$

MODULE kinetic_energy
  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE blocking
  USE horizontal_data
  USE logic
  USE output_data
  USE const
  use parallel_mod
  USE LMLoop_data, ONLY: llm,ulm
  USE communications, only: get_global_sum

  IMPLICIT NONE

  REAL(kind=8),allocatable :: e_pA(:),e_p_asA(:)
  REAL(kind=8),allocatable :: e_tA(:),e_t_asA(:)

CONTAINS
  SUBROUTINE initialize_kinetic_energy

    ALLOCATE( e_pA(n_r_max),e_p_asA(n_r_max) )
    ALLOCATE( e_tA(n_r_max),e_t_asA(n_r_max) )
    
  END SUBROUTINE initialize_kinetic_energy
  !********************************************************************
  SUBROUTINE get_e_kin(time,l_write,l_stop_time,n_e_sets, &
       &               w,dw,z,e_p,e_t,e_p_as,e_t_as,      &
       &               ekinR,ekinRave)
    !********************************************************************

    !--------------------------------------------------------------------
    !
    !  calculates kinetic energy  = 1/2 Integral (v^2 dV)
    !  integration in theta,phi by summation of spherical harmonics
    !  integration in r by using Chebycheff integrals
    !
    !  Output:
    !  e_p: Total poloidal     e_p_as: Total toroidal
    !  e_t: Axisym. poloidal   e_t_as: Axisym. toroidal
    !
    !--------------------------------------------------------------------

    USE integration, ONLY: rInt_R
    USE usefull, ONLY: cc2real

    IMPLICIT NONE

    !-- Input of constant parameters: field via common blocks
    INTEGER,intent(IN) :: n_e_sets

    !-- Input of variables:
    REAL(kind=8),intent(IN) :: time
    LOGICAL,intent(IN) :: l_write
    LOGICAL,intent(IN) :: l_stop_time

    !-- Input of scalar fields:
    COMPLEX(kind=8),intent(IN) :: w(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(IN) :: dw(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(IN) :: z(llm:ulm,n_r_max)

    !-- output: 
    REAL(kind=8),optional :: ekinR(n_r_max)     ! kinetic energy w radius
    REAL(kind=8),optional :: ekinRave(n_r_max)  ! kinetic energy w radius averaged in time
    REAL(kind=8),intent(OUT) :: e_p     ! poloidal energy
    REAL(kind=8),intent(OUT) :: e_t     ! toroidal energy
    REAL(kind=8),intent(OUT) :: e_p_as  ! axisymmetric poloidal energy
    REAL(kind=8),intent(OUT) :: e_t_as  ! axisymmetric toroidal energy

    !-- local:
    REAL(kind=8) :: e_p_es  ! equatorially symmetric poloidal energy 
    REAL(kind=8) :: e_t_es  ! equatorially symmetric toroidal energy
    REAL(kind=8) :: e_p_eas ! equator. & axially symmetric poloidal energy
    REAL(kind=8) :: e_t_eas ! equator. & axially symmetric toroidal energy

    REAL(kind=8) :: e_p_temp,e_t_temp
    REAL(kind=8) :: e_p_r(n_r_max)
    REAL(kind=8) :: e_t_r(n_r_max)
    REAL(kind=8) :: e_p_as_r(n_r_max)
    REAL(kind=8) :: e_t_as_r(n_r_max)
    REAL(kind=8) :: e_p_es_r(n_r_max)
    REAL(kind=8) :: e_t_es_r(n_r_max)
    REAL(kind=8) :: e_p_eas_r(n_r_max)
    REAL(kind=8) :: e_t_eas_r(n_r_max)

    REAL(kind=8),DIMENSION(n_r_max) :: e_p_r_global,e_t_r_global
    REAL(kind=8),DIMENSION(n_r_max) :: e_p_as_r_global,e_t_as_r_global
    REAL(kind=8),DIMENSION(n_r_max) :: e_p_es_r_global,e_t_es_r_global
    REAL(kind=8),DIMENSION(n_r_max) :: e_p_eas_r_global,e_t_eas_r_global

    INTEGER nR,lm,l,m
    REAL(kind=8) :: fac
    REAL(kind=8) :: O_rho ! 1/rho (anelastic)

    !-- time averaging of e(r):
    CHARACTER(len=80) :: filename
    REAL(kind=8) :: timeLast,timeTot,dt,surf
    SAVE timeLast,timeTot

    !-- end of declaration
    !---------------------------------------------------------------------
    
    !WRITE(*,"(A,6ES22.14)") "ekin: w,dw,z = ",get_global_sum( w(llm:ulm,:) ),&
    !     & get_global_sum( dw(llm:ulm,:) ), get_global_sum( z(llm:ulm,:) )

    DO nR=1,n_r_max
       e_p_r(nR)    =0.D0
       e_t_r(nR)    =0.D0
       e_p_as_r(nR) =0.D0
       e_t_as_r(nR) =0.D0
       e_p_es_r(nR) =0.D0
       e_t_es_r(nR) =0.D0
       e_p_eas_r(nR)=0.D0
       e_t_eas_r(nR)=0.D0
       O_rho        =orho1(nR)
       !DO lm=2,lm_max
       DO lm=MAX(2,llm),ulm
          l=lo_map%lm2l(lm)
          m=lo_map%lm2m(lm)

          e_p_temp= O_rho*dLh(st_map%lm2(l,m)) * ( &
               &      dLh(st_map%lm2(l,m))*or2(nR)*cc2real(w(lm,nR),m) &
               &      + cc2real(dw(lm,nR),m) )
          e_t_temp= O_rho*dLh(st_map%lm2(l,m)) * cc2real(z(lm,nR),m)
          !WRITE(*,"(A,3I4,ES22.14)") "e_p_temp = ",nR,l,m,e_p_temp
          IF ( m.EQ.0 ) THEN  ! axisymmetric part
             e_p_as_r(nR) = e_p_as_r(nR) + e_p_temp
             e_t_as_r(nR) = e_t_as_r(nR) + e_t_temp
             IF ( MOD(l,2).EQ.0 ) THEN
                e_p_eas_r(nR)=e_p_eas_r(nR)+e_p_temp
             ELSE
                e_t_eas_r(nR)=e_t_eas_r(nR)+e_t_temp
             END IF
          ELSE
             e_p_r(nR)=e_p_r(nR) + e_p_temp
             e_t_r(nR)=e_t_r(nR) + e_t_temp
          END IF
          IF ( MOD(l+m,2).EQ.0 ) THEN
             e_p_es_r(nR)=e_p_es_r(nR)+e_p_temp
          ELSE
             e_t_es_r(nR)=e_t_es_r(nR)+e_t_temp
          END IF

          !WRITE(*,"(8X,A,4I4,ES22.14)") "e_p_r: ",lm,l,m,nR,e_p_r(nR)
          
       END DO    ! do loop over lms in block
       e_p_r(nR)=e_p_r(nR)+e_p_as_r(nR)
       e_t_r(nR)=e_t_r(nR)+e_t_as_r(nR)
       !WRITE(*,"(4X,A,I4,2ES22.14)") "e_p_r: ",nR,e_p_r(nR),e_p_as_r(nR)
    END DO    ! radial grid points

    ! reduce over the ranks
    CALL MPI_Reduce(e_p_r,    e_p_r_global,     n_r_max, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(e_t_r,    e_t_r_global,     n_r_max, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(e_p_as_r, e_p_as_r_global,  n_r_max, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(e_t_as_r, e_t_as_r_global,  n_r_max, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(e_p_es_r, e_p_es_r_global,  n_r_max, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(e_t_es_r, e_t_es_r_global,  n_r_max, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(e_p_eas_r,e_p_eas_r_global, n_r_max, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(e_t_eas_r,e_t_eas_r_global, n_r_max, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)


    IF (rank.EQ.0) THEN
       !DO nR=1,n_r_max
       !   WRITE(*,"(4X,A,I4,ES22.14)") "e_p_r_global: ",nR,e_p_r_global(nR)
       !END DO
       !-- Radial Integrals:
       e_p    = rInt_R(e_p_r_global,    n_r_max,n_r_max,drx,  i_costf_init,d_costf_init)
       e_t    = rInt_R(e_t_r_global,    n_r_max,n_r_max,drx,  i_costf_init,d_costf_init)
       e_p_as = rInt_R(e_p_as_r_global, n_r_max,n_r_max,drx,  i_costf_init,d_costf_init)
       e_t_as = rInt_R(e_t_as_r_global, n_r_max,n_r_max,drx,  i_costf_init,d_costf_init)
       e_p_es = rInt_R(e_p_es_r_global, n_r_max,n_r_max,drx,  i_costf_init,d_costf_init)
       e_t_es = rInt_R(e_t_es_r_global, n_r_max,n_r_max,drx,  i_costf_init,d_costf_init)
       e_p_eas= rInt_R(e_p_eas_r_global,n_r_max,n_r_max,drx,  i_costf_init,d_costf_init)
       e_t_eas= rInt_R(e_t_eas_r_global,n_r_max,n_r_max,drx,  i_costf_init,d_costf_init)

       fac    =0.5D0*eScale
       e_p    =fac*e_p
       e_t    =fac*e_t
       e_p_as =fac*e_p_as
       e_t_as =fac*e_t_as
       e_p_es =fac*e_p_es
       e_t_es =fac*e_t_es
       e_p_eas=fac*e_p_eas
       e_t_eas=fac*e_t_eas

       !-- OUTPUT:
       IF (PRESENT(ekinR)) THEN
          DO nR=1,n_r_max
             ekinR(nR)=fac*(e_p_r_global(nR)+e_t_r_global(nR))
          END DO
       END IF
       IF ( l_write ) THEN
          IF ( l_save_out ) THEN
             OPEN(n_e_kin_file,FILE=e_kin_file,STATUS='UNKNOWN',       &
                  &             POSITION='APPEND')
          END IF
          WRITE(n_e_kin_file,'(1P,D20.12,8D16.8)')                     &
               & time*tScale, &  ! 1
               & e_p,e_t,       &! 2,3
               & e_p_as,e_t_as, &! 4,5
               & e_p_es,e_t_es, &! 6,7
               & e_p_eas,e_t_eas ! 8,9
          IF ( l_save_out ) CLOSE(n_e_kin_file)
       END IF

       ! NOTE: n_e_sets=0 prevents averaging
       IF ( n_e_sets.EQ.1 ) THEN
          timeTot=1.D0
          e_pA    = e_p_r_global
          e_p_asA = e_p_r_global
          e_tA    = e_t_r_global
          e_t_asA = e_t_r_global
       ELSE IF ( n_e_sets.EQ.2 ) THEN
          dt=time-timeLast
          timeTot=2.D0*dt
          e_pA    = dt*(e_pA   +e_p_r_global   )
          e_p_asA = dt*(e_p_asA+e_p_as_r_global)
          e_tA    = dt*(e_tA   +e_t_r_global   )
          e_t_asA = dt*(e_t_asA+e_t_as_r_global)
       ELSE
          dt=time-timeLast
          timeTot=timeTot+dt
          e_pA    = e_pA    + dt*e_p_r_global
          e_p_asA = e_p_asA + dt*e_p_as_r_global
          e_tA    = e_tA    + dt*e_t_r_global
          e_t_asA = e_t_asA + dt*e_t_as_r_global
       END IF

       !WRITE(*,"(A,2ES22.14)") "e_pA, e_tA = ",SUM( e_pA ),SUM( e_tA )
       IF ( l_stop_time .AND. (n_e_sets.GT.1) ) THEN
          fac=0.5D0*eScale
          filename='eKinR.'//tag
          OPEN(99,FILE=filename,STATUS='UNKNOWN')
          IF (PRESENT(ekinRave)) THEN
             ekinRave(1)      =fac*(e_pA(1)+e_tA(1))/timetot
             ekinRave(n_r_max)=fac*(e_pA(n_r_max)+e_tA(n_r_max))/timetot
             DO nR=1,n_r_max
                ekinRave(nR)  =fac*e_pA(nR)/timetot+fac*e_tA(nR)/timetot
             END DO
          END IF
          DO nR=1,n_r_max
             surf=4.D0*pi*r(nR)**2
             WRITE(99,'(2x,9D12.4)',advance='no') r(nR),                           &
                  &               fac*e_pA(nR)/timetot,                              &
                  &               fac*e_p_asA(nR)/timetot,                           &
                  &               fac*e_tA(nR)/timetot,                              &
                  &               fac*e_t_asA(nR)/timetot,                           &
                  &               fac*e_pA(nR)/timetot/surf,                         &
                  &               fac*e_p_asA(nR)/timetot/surf,                      &
                  &               fac*e_tA(nR)/timetot/surf,                         &
                  &               fac*e_t_asA(nR)/timetot/surf

             IF (PRESENT(ekinR)) THEN
                WRITE(99,'(D12.4)',advance='no') ekinR(nR)
             ELSE
                WRITE(99,'(A)') ""
             END IF
             IF (PRESENT(ekinRave)) THEN
                WRITE(99,'(D12.4)') ekinRave(nR)
             ELSE
                WRITE(99,'(A)') ""
             END IF
          END DO
          CLOSE(99)
       END IF
       timeLast=time
    END IF

    ! broadcast the output arguments of the function to have them on all ranks
    ! e_p,e_t,e_p_as,e_t_as
    CALL MPI_Bcast(e_p,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Bcast(e_t,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Bcast(e_p_as,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Bcast(e_t_as,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    IF (PRESENT(ekinR)) THEN
       CALL MPI_Bcast(ekinR,n_r_max,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    END IF
    IF (PRESENT(ekinRave)) THEN
       CALL MPI_Bcast(ekinRave,n_r_max,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    END IF

    RETURN
    !-----------------------------------------------------------------------------
  END SUBROUTINE get_e_kin
END MODULE kinetic_energy
