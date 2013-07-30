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

  IMPLICIT NONE

  REAL(kind=8),allocatable :: e_pA(:),e_p_asA(:)
  REAL(kind=8),allocatable :: e_tA(:),e_t_asA(:)

CONTAINS
  SUBROUTINE initialize_kinetic_energy

    ALLOCATE( e_pA(n_r_max),e_p_asA(n_r_max) )
    ALLOCATE( e_tA(n_r_max),e_t_asA(n_r_max) )
    
  END SUBROUTINE initialize_kinetic_energy
  !********************************************************************
  SUBROUTINE get_e_kin(time,l_write,l_stop_time,n_e_sets,         &
       &                            w,dw,z,e_p,e_t,e_p_as,e_t_as, &
       &                                          ekinR)
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
    COMPLEX(kind=8),intent(IN) :: w(lm_max,n_r_max)
    COMPLEX(kind=8),intent(IN) :: dw(lm_max,n_r_max)
    COMPLEX(kind=8),intent(IN) :: z(lm_max,n_r_max)

    !-- output: 
    REAL(kind=8),optional :: ekinR(n_r_max)     ! kinetic energy w radius
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

    INTEGER nR,lm,l,m
    REAL(kind=8) :: fac
    REAL(kind=8) :: O_rho ! 1/rho (anelastic)

    !-- time averaging of e(r):
    CHARACTER(len=80) :: filename
    REAL(kind=8) :: timeLast,timeTot,dt,surf
    SAVE timeLast,timeTot

    !-- end of declaration
    !---------------------------------------------------------------------


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
       DO lm=2,lm_max
          l=lm2l(lm)
          m=lm2m(lm)

          e_p_temp= O_rho*dLh(lm) * ( dLh(lm)*or2(nR)*cc2real(w(lm,nR),m) + cc2real(dw(lm,nR),m) )
          e_t_temp= O_rho*dLh(lm)*cc2real(z(lm,nR),m)

          IF ( m.EQ.0 ) THEN  ! axisymmetric part
             e_p_as_r(nR)=e_p_as_r(nR) + e_p_temp
             e_t_as_r(nR)=e_t_as_r(nR)+ e_t_temp
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

       END DO    ! do loop over lms in block
       e_p_r(nR)=e_p_r(nR)+e_p_as_r(nR)
       e_t_r(nR)=e_t_r(nR)+e_t_as_r(nR)
    END DO    ! radial grid points

    !-- Radial Integrals:
    e_p    = rInt_R(e_p_r,    n_r_max,n_r_max,drx,  i_costf_init,d_costf_init)
    e_t    = rInt_R(e_t_r,    n_r_max,n_r_max,drx,  i_costf_init,d_costf_init)
    e_p_as = rInt_R(e_p_as_r, n_r_max,n_r_max,drx,  i_costf_init,d_costf_init)
    e_t_as = rInt_R(e_t_as_r, n_r_max,n_r_max,drx,  i_costf_init,d_costf_init)
    e_p_es = rInt_R(e_p_es_r, n_r_max,n_r_max,drx,  i_costf_init,d_costf_init)
    e_t_es = rInt_R(e_t_es_r, n_r_max,n_r_max,drx,  i_costf_init,d_costf_init)
    e_p_eas= rInt_R(e_p_eas_r,n_r_max,n_r_max,drx,  i_costf_init,d_costf_init)
    e_t_eas= rInt_R(e_t_eas_r,n_r_max,n_r_max,drx,  i_costf_init,d_costf_init)

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
    DO nR=1,n_r_max
       ekinR(nR)=fac*(e_p_r(nR)+e_t_r(nR))
    end do
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
       DO nR=1,n_r_max
          e_pA(nR)   =e_p_r(nR)
          e_p_asA(nR)=e_p_r(nR)
          e_tA(nR)   =e_t_r(nR)
          e_t_asA(nR)=e_t_r(nR)
       END DO
    ELSE IF ( n_e_sets.EQ.2 ) THEN
       dt=time-timeLast
       timeTot=2.D0*dt
       DO nR=1,n_r_max
          e_pA(nR)   =dt*(e_pA(nR)   +e_p_r(nR)   )
          e_p_asA(nR)=dt*(e_p_asA(nR)+e_p_as_r(nR))
          e_tA(nR)   =dt*(e_tA(nR)   +e_t_r(nR)   )
          e_t_asA(nR)=dt*(e_t_asA(nR)+e_t_as_r(nR))
       END DO
    ELSE
       dt=time-timeLast
       timeTot=timeTot+dt
       DO nR=1,n_r_max
          e_pA(nR)   =e_pA(nR)   +dt*e_p_r(nR)
          e_p_asA(nR)=e_p_asA(nR)+dt*e_p_as_r(nR)
          e_tA(nR)   =e_tA(nR)   +dt*e_t_r(nR)
          e_t_asA(nR)=e_t_asA(nR)+dt*e_t_as_r(nR)
       END DO
    END IF
    IF ( l_stop_time .and. (n_e_sets.gt.1) ) THEN
       fac=0.5D0*eScale
       filename='eKinR.'//tag
       OPEN(99,FILE=filename,STATUS='UNKNOWN')
       DO nR=1,n_r_max
          surf=4.D0*pi*r(nR)**2
          WRITE(99,'(2x,10D12.4)') r(nR),                           &
               &               fac*e_pA(nR)/timetot,                              &
               &               fac*e_p_asA(nR)/timetot,                           &
               &               fac*e_tA(nR)/timetot,                              &
               &               fac*e_t_asA(nR)/timetot,                           &
               &               fac*e_pA(nR)/timetot/surf,                         &
               &               fac*e_p_asA(nR)/timetot/surf,                      &
               &               fac*e_tA(nR)/timetot/surf,                         &
               &               fac*e_t_asA(nR)/timetot/surf,                      &
               &               ekinR(nR)
       END DO
       CLOSE(99)
    END IF
    timeLast=time


    RETURN
    !-----------------------------------------------------------------------------
  END SUBROUTINE get_e_kin
END MODULE kinetic_energy
