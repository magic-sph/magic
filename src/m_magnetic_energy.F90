!$Id$
MODULE magnetic_energy
  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE blocking
  USE horizontal_data
  USE logic
  USE movie_data
  USE output_data
  USE const
  USE Bext
  USE LMLoop_data, ONLY: llmMag, ulmMag
  IMPLICIT NONE

  REAL(kind=8),ALLOCATABLE :: e_dipA(:)
  REAL(kind=8),ALLOCATABLE :: e_pA(:),e_p_asA(:)
  REAL(kind=8),ALLOCATABLE :: e_tA(:),e_t_asA(:)
  
contains
  SUBROUTINE initialize_magnetic_energy

    ALLOCATE( e_dipA(n_r_max) )
    ALLOCATE( e_pA(n_r_max),e_p_asA(n_r_max) )
    ALLOCATE( e_tA(n_r_max),e_t_asA(n_r_max) )
    
  END SUBROUTINE initialize_magnetic_energy

  !********************************************************************
  SUBROUTINE get_e_mag(time,l_write,l_stop_time,n_e_sets,         &
       &               b,db,aj,b_ic,db_ic,aj_ic,         &
       &               e_p,e_t,e_p_as,e_t_as,         &
       &               e_p_ic,e_t_ic,e_p_as_ic,e_t_as_ic,         &
       &               e_p_os,e_p_as_os,e_cmb,Dip,DipCMB,         &
       &               elsAnel)
    !********************************************************************

    !------------ This is release 2 level 6  --------------!
    !------------ Created on 2/25/02  by JW. --------------!

    !--------------------------------------------------------------------
    !
    !  calculates magnetic energy  = 1/2 Integral(B^2 dV)
    !  integration in theta,phi by summation over harmonic coeffs.
    !  integration in r by Chebycheff integrals
    !
    !  Output:
    !  enbp: Total poloidal        enbt: Total toroidal
    !  apome: Axisym. poloidal     atome: Axisym. toroidal
    !
    !--------------------------------------------------------------------
    USE integration, ONLY: rInt_R,rIntIC
    USE usefull, ONLY: cc2real,cc22real

    IMPLICIT NONE

    !-- Input of constant parameters: 
    INTEGER,intent(IN) :: n_e_sets

    !-- Input of variables:
    REAL(kind=8),intent(IN) :: time
    LOGICAL,intent(IN) :: l_write
    LOGICAL,intent(IN) :: l_stop_time

    !-- Input of scalar fields:
    COMPLEX(kind=8),intent(IN) :: b(llmMag:ulmMag,n_r_maxMag)
    COMPLEX(kind=8),intent(IN) :: db(llmMag:ulmMag,n_r_maxMag)
    COMPLEX(kind=8),intent(IN) :: aj(llmMag:ulmMag,n_r_maxMag)
    COMPLEX(kind=8),intent(IN) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
    COMPLEX(kind=8),intent(IN) :: db_ic(llmMag:ulmMag,n_r_ic_maxMag)
    COMPLEX(kind=8),intent(IN) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)

    !-- output:
    REAL(kind=8),intent(OUT) :: e_p,e_t            ! poloidal, toroidal energy
    REAL(kind=8),intent(OUT) :: e_p_as,e_t_as      ! axisymmetric poloidal, toroidal energy
    REAL(kind=8),intent(OUT) :: e_p_ic,e_t_ic   
    REAL(kind=8),intent(OUT) :: e_p_as_ic,e_t_as_ic
    REAL(kind=8),intent(OUT) :: e_p_os,e_p_as_os
    REAL(kind=8),intent(OUT) :: e_cmb
    REAL(kind=8),intent(OUT) :: Dip,DipCMB
    REAL(kind=8),intent(OUT) :: elsAnel

    !-- local:
    INTEGER :: nR,lm,l,m,l1m0,l1m1
    INTEGER :: l_geo

    REAL(kind=8),DIMENSION(n_r_max) :: e_p_r, e_p_r_global
    REAL(kind=8),DIMENSION(n_r_max) :: e_t_r, e_t_r_global
    REAL(kind=8),DIMENSION(n_r_max) :: els_r, els_r_global
    REAL(kind=8),DIMENSION(n_r_max) :: e_p_as_r, e_p_as_r_global
    REAL(kind=8),DIMENSION(n_r_max) :: e_t_as_r, e_t_as_r_global
    REAL(kind=8),DIMENSION(n_r_max) :: e_p_es_r, e_p_es_r_global,e_t_es_r, e_t_es_r_global
    REAL(kind=8),DIMENSION(n_r_max) :: e_p_eas_r, e_p_eas_r_global,e_t_eas_r, e_t_eas_r_global
    REAL(kind=8),DIMENSION(n_r_max) :: e_dipole_r, e_dipole_r_global
    REAL(kind=8),DIMENSION(n_r_max) :: e_dipole_ax_r, e_dipole_ax_r_global

    REAL(kind=8),DIMENSION(n_r_ic_max) :: e_p_ic_r, e_p_ic_r_global
    REAL(kind=8),DIMENSION(n_r_ic_max) :: e_t_ic_r, e_t_ic_r_global   
    REAL(kind=8),DIMENSION(n_r_ic_max) :: e_p_as_ic_r, e_p_as_ic_r_global
    REAL(kind=8),DIMENSION(n_r_ic_max) :: e_t_as_ic_r, e_t_as_ic_r_global

    REAL(kind=8) :: e_geo,e_es_geo,e_as_geo,e_eas_geo
    REAL(kind=8) :: e_geo_global,e_es_geo_global,e_as_geo_global,e_eas_geo_global
    REAL(kind=8) :: e_p_ic_global, e_p_as_ic_global, e_p_os_global, e_p_as_os_global
    REAL(kind=8) :: e_p_e,e_p_as_e
    REAL(kind=8) :: e_p_e_global,e_p_as_e_global

    REAL(kind=8) :: r_ratio,fac
    REAL(kind=8) :: e_p_temp,e_t_temp
    REAL(kind=8) :: e_dipole, e_dipole_ax, e_dipole_ax_cmb, e_dipole_e, e_dipole_e_global
    REAL(kind=8) :: e_p_e_ratio
    REAL(kind=8) :: O_r_icb_E_2,rad
    REAL(kind=8) :: e_p_es,e_t_es,e_es_cmb,e_as_cmb
    REAL(kind=8) :: e_p_eas,e_t_eas,e_eas_cmb

    REAL(kind=8) :: e_dip_cmb,eTot,eDR
    REAL(kind=8) :: theta_dip,phi_dip

    COMPLEX(kind=8) :: r_dr_b,b10,b11

    !-- time averaging of e(r):
    CHARACTER(len=80) :: filename
    REAL(kind=8) :: timeLast,timeTot,dt,surf
    LOGICAL :: rank_has_l1m0,rank_has_l1m1
    INTEGER :: status(MPI_STATUS_SIZE),sr_tag,request1,request2
    SAVE timeLast,timeTot

    !-- end of declaration
    !---------------------------------------------------------------------

    ! some arbitrary send recv tag
    sr_tag=18654

    l_geo=11   ! max degree for geomagnetic field seen on Earth  
    ! surface

    e_p      =0.D0
    e_t      =0.D0
    e_p_as   =0.D0
    e_t_as   =0.D0
    e_p_ic   =0.D0
    e_t_ic   =0.D0
    e_p_as_ic=0.D0
    e_t_as_ic=0.D0
    e_p_os   =0.D0
    e_p_as_os=0.D0
    e_geo    =0.D0
    e_es_geo =0.D0
    e_as_geo =0.D0
    e_eas_geo=0.D0
    Dip      =0.D0
    DipCMB   =0.D0

    IF ( .NOT.( l_mag .OR. l_mag_LF ) ) RETURN

    DO nR=1,n_r_max

       e_p_r(nR)     =0.D0
       e_t_r(nR)     =0.D0
       e_p_as_r(nR)  =0.D0
       e_t_as_r(nR)  =0.D0
       e_p_es_r(nR)  =0.D0
       e_t_es_r(nR)  =0.D0
       e_p_eas_r(nR) =0.D0
       e_t_eas_r(nR) =0.D0
       e_dipole_r(nR)=0.D0
       e_dipole_ax_r(nR)=0.D0

       !DO lm=2,lm_max
       DO lm=MAX(2,llmMag),ulmMag
          l=lo_map%lm2l(lm)
          m=lo_map%lm2m(lm)

          e_p_temp= dLh(st_map%lm2(l,m)) * ( dLh(st_map%lm2(l,m))*or2(nR)*cc2real(b(lm,nR),m) &
               &                             + cc2real(db(lm,nR),m) )
          e_t_temp= dLh(st_map%lm2(l,m)) * cc2real(aj(lm,nR),m)

          IF ( m.EQ.0 ) THEN  ! axisymmetric part 
             e_p_as_r(nR)=e_p_as_r(nR) + e_p_temp
             e_t_as_r(nR)=e_t_as_r(nR) + e_t_temp
             IF ( MOD(l,2).EQ.1 ) THEN
                e_p_eas_r(nR)=e_p_eas_r(nR)+e_p_temp
             ELSE
                e_t_eas_r(nR)=e_t_eas_r(nR)+e_t_temp
             END IF
          ELSE
             e_p_r(nR)=e_p_r(nR) + e_p_temp
             e_t_r(nR)=e_t_r(nR) + e_t_temp
          END IF
          IF ( MOD(l+m,2).EQ.1 ) THEN
             e_p_es_r(nR)=e_p_es_r(nR) + e_p_temp
          ELSE
             e_t_es_r(nR)=e_t_es_r(nR) + e_t_temp
          END IF
          IF ( l.LE.l_geo .AND. nR.EQ.n_r_cmb ) THEN     
             e_geo    =e_geo    +e_p_temp
             IF ( MOD(l+m,2).EQ.1 ) e_es_geo = e_es_geo + e_p_temp
             IF ( m.EQ.0 )          e_as_geo = e_as_geo + e_p_temp
             IF ( MOD(l+m,2).EQ.1 .AND. m.EQ.0 )                    &
                  &                 e_eas_geo=e_eas_geo+e_p_temp
          END IF

          if ( l.eq.1 .AND. m.eq.0 ) e_dipole_ax_r(nR)=e_p_temp
          IF ( l.EQ.1 )  e_dipole_r(nR)=e_dipole_r(nR)+e_p_temp

       end do    ! do loop over lms in block 

       e_p_r(nR)=e_p_r(nR)+e_p_as_r(nR)
       e_t_r(nR)=e_t_r(nR)+e_t_as_r(nR)

       ! In anelastic models it is also interesting to have Lambda
       els_r(nR)=(e_p_r(nR)+e_t_r(nR))*orho1(nR)*sigma(nR)

    end do    ! radial grid points 

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
    CALL MPI_Reduce(e_dipole_ax_r,e_dipole_ax_r_global, n_r_max, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(e_dipole_r,   e_dipole_r_global,    n_r_max, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(els_r,   els_r_global,    n_r_max, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    ! reduce some scalars
    CALL MPI_Reduce(e_geo,     e_geo_global,    1, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(e_es_geo,  e_es_geo_global, 1, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(e_as_geo,  e_as_geo_global, 1, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(e_eas_geo, e_eas_geo_global,1, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    IF (rank.EQ.0) THEN
       !-- Get Values at CMB:
       e_cmb          =e_p_r_global(n_r_cmb)+e_t_r_global(n_r_cmb)
       e_dip_cmb      =e_dipole_r_global(n_r_cmb)
       e_dipole_ax_cmb=e_dipole_ax_r_global(n_r_cmb)
       e_es_cmb       =e_p_es_r_global(n_r_cmb)
       e_as_cmb       =e_p_as_r_global(n_r_cmb)
       e_eas_cmb      =e_p_eas_r_global(n_r_cmb)

       ! NOTE: n_e_sets=0 prevents averaging
       IF ( n_e_sets.EQ.1 ) THEN
          timeTot=1.D0
          DO nR=1,n_r_max
             e_dipA(nR) =e_dipole_r_global(nR)
             e_pA(nR)   =e_p_r_global(nR)
             e_p_asA(nR)=e_p_r_global(nR)
             e_tA(nR)   =e_t_r_global(nR)
             e_t_asA(nR)=e_t_r_global(nR)
          END DO
       ELSE IF ( n_e_sets.EQ.2 ) THEN
          dt=time-timeLast
          timeTot=2.D0*dt
          DO nR=1,n_r_max
             e_dipA(nR) =dt*(e_dipA(nR) +e_dipole_r_global(nR))
             e_pA(nR)   =dt*(e_pA(nR)   +e_p_r_global(nR)     )
             e_p_asA(nR)=dt*(e_p_asA(nR)+e_p_as_r_global(nR)  )
             e_tA(nR)   =dt*(e_tA(nR)   +e_t_r_global(nR)     )
             e_t_asA(nR)=dt*(e_t_asA(nR)+e_t_as_r_global(nR)  )
          END DO
       ELSE
          dt=time-timeLast
          timeTot=timeTot+dt
          DO nR=1,n_r_max
             e_dipA(nR) =e_dipA(nR) +dt*e_dipole_r_global(nR)
             e_pA(nR)   =e_pA(nR)   +dt*e_p_r_global(nR)
             e_p_asA(nR)=e_p_asA(nR)+dt*e_p_as_r_global(nR)
             e_tA(nR)   =e_tA(nR)   +dt*e_t_r_global(nR)
             e_t_asA(nR)=e_t_asA(nR)+dt*e_t_as_r_global(nR)
          END DO
       END IF
       IF ( l_stop_time ) THEN
          fac=0.5D0*LFfac*eScale
          filename='eMagR.'//tag
          OPEN(99,FILE=filename,STATUS='UNKNOWN')
          DO nR=1,n_r_max
             eTot=e_pA(nR)+e_tA(nR)
             IF ( e_dipA(nR) .LT. 1.D-4*eTot ) THEN
                eDR=0.D0
             ELSE
                eDR=e_dipA(nR)/eTot
             END IF
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
                  &               eDR
          END DO
          CLOSE(99)
       END IF
       timeLast=time


       !-- Radial integrals:
       e_p        =rInt_R(e_p_r_global        ,n_r_max,n_r_max,drx, i_costf_init,d_costf_init)
       e_t        =rInt_R(e_t_r_global        ,n_r_max,n_r_max,drx, i_costf_init,d_costf_init)
       e_p_as     =rInt_R(e_p_as_r_global     ,n_r_max,n_r_max,drx, i_costf_init,d_costf_init)
       e_t_as     =rInt_R(e_t_as_r_global     ,n_r_max,n_r_max,drx, i_costf_init,d_costf_init)
       e_p_es     =rInt_R(e_p_es_r_global     ,n_r_max,n_r_max,drx, i_costf_init,d_costf_init)
       e_t_es     =rInt_R(e_t_es_r_global     ,n_r_max,n_r_max,drx, i_costf_init,d_costf_init)
       e_p_eas    =rInt_R(e_p_eas_r_global    ,n_r_max,n_r_max,drx, i_costf_init,d_costf_init)
       e_t_eas    =rInt_R(e_t_eas_r_global    ,n_r_max,n_r_max,drx, i_costf_init,d_costf_init)
       e_dipole   =rInt_R(e_dipole_r_global   ,n_r_max,n_r_max,drx, i_costf_init,d_costf_init)
       e_dipole_ax=rInt_R(e_dipole_ax_r_global,n_r_max,n_r_max,drx, i_costf_init,d_costf_init)
       elsAnel    =rInt_R(els_r_global        ,n_r_max,n_r_max,drx, i_costf_init,d_costf_init)

       fac=0.5*LFfac*eScale
       e_p        =fac*e_p
       e_t        =fac*e_t
       e_p_as     =fac*e_p_as
       e_t_as     =fac*e_t_as
       e_p_es     =fac*e_p_es
       e_t_es     =fac*e_t_es
       e_p_eas    =fac*e_p_eas
       e_t_eas    =fac*e_t_eas
       e_dipole   =fac*e_dipole
       e_dipole_ax=fac*e_dipole_ax
       elsAnel    =eScale*elsAnel

       e_cmb          =fac*e_cmb
       e_dip_cmb      =fac*e_dip_cmb
       e_dipole_ax_cmb=fac*e_dipole_ax_cmb
       e_es_cmb       =fac*e_es_cmb
       e_as_cmb       =fac*e_as_cmb
       e_eas_cmb      =fac*e_eas_cmb
       e_geo          =fac*e_geo_global
       e_es_geo       =fac*e_es_geo_global
       e_as_geo       =fac*e_as_geo_global
       e_eas_geo      =fac*e_eas_geo_global
    END IF

    !-- Inner core:

    IF ( l_cond_ic ) THEN 

       O_r_icb_E_2=1.d0/(r(n_r_max)*r(n_r_max))

       DO nR=1,n_r_ic_max

          r_ratio=r_ic(nR)/r_ic(1)

          e_p_ic_r(nR)   =0.D0
          e_t_ic_r(nR)   =0.D0
          e_p_as_ic_r(nR)=0.D0
          e_t_as_ic_r(nR)=0.D0

          !DO lm=2,lm_max
          DO lm=MAX(2,llmMag),ulmMag
             l=lo_map%lm2l(lm)
             m=lo_map%lm2m(lm)
             r_dr_b=r_ic(nR)*db_ic(lm,nR)

             e_p_temp=     dLh(st_map%lm2(l,m))*O_r_icb_E_2*r_ratio**(2*l) * (   &
                  &             DBLE((l+1)*(2*l+1))*cc2real(b_ic(lm,nR),m)      +    &
                  &                DBLE(2*(l+1))*cc22real(b_ic(lm,nR),r_dr_b,m) +    &
                  &                                 cc2real(r_dr_b,m)            )
             e_t_temp=  dLh(st_map%lm2(l,m))*r_ratio**(2*l+2) *                  &
                  &                                 cc2real(aj_ic(lm,nR),m)

             IF ( m.EQ.0 ) then  ! axisymmetric part
                e_p_as_ic_r(nR)=e_p_as_ic_r(nR) + e_p_temp
                e_t_as_ic_r(nR)=e_t_as_ic_r(nR) + e_t_temp
             ELSE
                e_p_ic_r(nR)   =e_p_ic_r(nR) + e_p_temp
                e_t_ic_r(nR)   =e_t_ic_r(nR) + e_t_temp
             END IF

          END DO    ! do loop over lms in block

          e_p_ic_r(nR)=e_p_ic_r(nR)+e_p_as_ic_r(nR)
          e_t_ic_r(nR)=e_t_ic_r(nR)+e_t_as_ic_r(nR)

       END DO    ! radial grid points

       ! reduce over the ranks
       CALL MPI_Reduce(e_p_ic_r,    e_p_ic_r_global,    n_r_ic_max, &
            & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
       CALL MPI_Reduce(e_t_ic_r,    e_t_ic_r_global,    n_r_ic_max, &
            & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
       CALL MPI_Reduce(e_p_as_ic_r, e_p_as_ic_r_global, n_r_ic_max, &
            & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
       CALL MPI_Reduce(e_t_as_ic_r, e_t_as_ic_r_global, n_r_ic_max, &
            & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

       IF (rank.EQ.0) THEN
          e_p_ic   =rIntIC(e_p_ic_r_global,n_r_ic_max,dr_fac_ic,              &
               &                      i_costf1_ic_init,d_costf1_ic_init)
          e_t_ic   =rIntIC(e_t_ic_r_global,n_r_ic_max,dr_fac_ic,              &
               &                      i_costf1_ic_init,d_costf1_ic_init)
          e_p_as_ic=rIntIC(e_p_as_ic_r_global,n_r_ic_max,dr_fac_ic,           &
               &                      i_costf1_ic_init,d_costf1_ic_init)
          e_t_as_ic=rIntIC(e_t_as_ic_r_global,n_r_ic_max,dr_fac_ic,           &
               &                      i_costf1_ic_init,d_costf1_ic_init)
          fac=LFfac*eScale/2.D0
          e_p_ic   =fac*e_p_ic
          e_t_ic   =fac*e_t_ic
          e_p_as_ic=fac*e_p_as_ic
          e_t_as_ic=fac*e_t_as_ic
       END IF

    ELSE

       !DO lm=2,lm_max
       DO lm=MAX(2,llmMag),ulmMag
          l=lo_map%lm2l(lm)
          m=lo_map%lm2m(lm)
          fac=DBLE(l*(l+1)*(l+1))
          e_p_temp=fac*cc2real(b(lm,n_r_max),m)
          e_p_ic=e_p_ic + e_p_temp
          IF ( m.EQ.0 ) e_p_as_ic=e_p_as_ic+e_p_temp
       END DO    ! do loop over lms in block

       CALL MPI_Reduce(e_p_ic,    e_p_ic_global,   1, &
            & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
       CALL MPI_Reduce(e_p_as_ic, e_p_as_ic_global,1, &
            & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

       IF (rank.EQ.0) THEN
          fac      =0.5D0*LFfac/r_icb*eScale
          e_p_ic   =fac*e_p_ic_global
          e_t_ic   =0.D0
          e_p_as_ic=fac*e_p_as_ic_global
          e_t_as_ic=0.D0
       END IF

    END IF  ! conducting inner core ?


    !-- Outside energy:
    nR=n_r_cmb
    e_p_os   =0.D0
    e_p_as_os=0.D0
    !DO lm=2,lm_max
    DO lm=MAX(2,llmMag),ulmMag
       l=lo_map%lm2l(lm)
       m=lo_map%lm2m(lm)
       fac=DBLE( l*l*(l+1) )
       e_p_temp=fac*cc2real(b(lm,nR),m)
       e_p_os  =e_p_os + e_p_temp
       IF ( m.EQ.0 ) e_p_as_os=e_p_as_os + e_p_temp
    END DO

    CALL MPI_Reduce(e_p_os,    e_p_os_global,   1, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(e_p_as_os, e_p_as_os_global,1, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    IF (rank.EQ.0) THEN
       fac      =0.5D0*LFfac/r_cmb*eScale
       e_p_os   =fac*e_p_os_global
       e_p_as_os=fac*e_p_as_os_global
    END IF

    !-- External potential field energy in Uli case (n_imp=1)
    e_p_e     =0.D0
    e_p_as_e  =0.D0
    e_dipole_e=0.D0
    IF ( n_imp.EQ.1 ) THEN
       !DO lm=2,lm_max
       DO lm=MAX(2,llmMag),ulmMag
          l=lo_map%lm2l(lm)
          m=lo_map%lm2m(lm)
          fac=DBLE(l*(l+1)**2*(2*l+1)) *                            &
               &              1.D0/(rrMP**(2*l+1)-1.D0)
          e_p_temp=fac*cc2real(b(lm,nR),m)
          e_p_e   =e_p_e  + e_p_temp
          IF ( m.EQ.0 ) e_p_as_e =e_p_as_e  + e_p_temp
          IF ( l.EQ.1 ) e_dipole_e=e_dipole_e+e_p_temp
       END DO

       CALL MPI_Reduce(e_p_e,    e_p_e_global,   1, &
            & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
       CALL MPI_Reduce(e_p_as_e, e_p_as_e_global,1, &
            & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
       CALL MPI_Reduce(e_dipole_e, e_dipole_e_global,1, &
            & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
       
       IF (rank.EQ.0) THEN
          fac       =0.5D0*LFfac/r_cmb**2*eScale
          e_p_e     =fac*e_p_e_global 
          e_p_as_e  =fac*e_p_as_e_global
          e_dipole_e=fac*e_dipole_e_global
       END IF
    END IF


    !-- Output of OC and outside energies:
    IF (rank.EQ.0) THEN
       IF ( l_write ) THEN
          IF ( l_save_out ) THEN
             OPEN(n_e_mag_oc_file,FILE=e_mag_oc_file,STATUS='UNKNOWN', &
                  &             POSITION='APPEND')
          END IF
          WRITE(n_e_mag_oc_file,'(1P,D20.12,12D16.8)')                 &
               &                             time*tScale,                         &! 1
               &                             e_p,e_t,                             &! 2,3
               &                             e_p_as,e_t_as,                       &! 4,5
               &                             e_p_os,e_p_as_os,                    &! 6,7
               &                             e_p_es,e_t_es,                       &! 8,9
               &                             e_p_eas,e_t_eas,                     &! 10,11
               &                             e_p_e,e_p_as_e    ! 12,13
          IF ( l_save_out ) CLOSE(n_e_mag_oc_file)
       END IF

       !-- Output of IC energies:
       IF ( l_write ) THEN
          IF ( l_save_out ) THEN
             OPEN(n_e_mag_ic_file,FILE=e_mag_ic_file,STATUS='UNKNOWN', &
                  &             POSITION='APPEND')
          END IF
          WRITE(n_e_mag_ic_file,'(1P,D20.12,4D16.8)')                  &
               &                       time*tScale,                               &
               &                       e_p_ic,e_t_ic,                             &
               &                       e_p_as_ic,e_t_as_ic
          IF ( l_save_out ) CLOSE(n_e_mag_ic_file)
       END IF
    END IF

    l1m0=lo_map%lm2(1,0)
    l1m1=lo_map%lm2(1,1)
    !IF ((l1m0.LT.lmStartB(1)).OR.(l1m0.GT.lmStopB(1))) THEN
    !   IF (rank.EQ.0) PRINT*,"in get_e_mag, dipole part: l1m0 not on rank 0"
    !   stop
    !END IF
    !IF (  ( l1m1.GT.0).AND.&
    !     &( (l1m1.LT.lmStartB(1)).OR.(l1m1.GT.lmStopB(1)) ) ) THEN
    !   IF (rank.EQ.0) PRINT*,"in get_e_mag, dipole part: l1m1 not on rank 0:"
    !   stop
    !END IF
    rank_has_l1m0=.FALSE.
    rank_has_l1m1=.FALSE.
    !WRITE(*,"(I5,A,2I5,A,2I5)") rank,": l1m0,l1m1 = ",l1m0,l1m1,&
    !     & ", lm block: ",lmStartB(rank+1),lmStopB(rank+1)
    IF ( (l1m0.GE.lmStartB(rank+1)) .AND. (l1m0.LE.lmStopB(rank+1)) ) THEN
       b10=b(l1m0,n_r_cmb)
       IF (rank.NE.0) THEN
          CALL MPI_Send(b10,1,MPI_DOUBLE_COMPLEX,0,sr_tag,MPI_COMM_WORLD,ierr)
       END IF
       rank_has_l1m0=.TRUE.
    END IF
    IF (l1m1.GT.0) THEN
       IF ( (l1m1.GE.lmStartB(rank+1)) .AND. (l1m1.LE.lmStopB(rank+1)) ) THEN
          b11=b(l1m1,n_r_cmb)
          IF (rank.NE.0) THEN
             CALL MPI_Send(b11,1,MPI_DOUBLE_COMPLEX,0,sr_tag+1,MPI_COMM_WORLD,ierr)
          END IF
          rank_has_l1m1=.TRUE.
       END IF
    ELSE
       b11=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
       rank_has_l1m1=.TRUE.
    END IF

       
    IF (rank.EQ.0) THEN
       !-- Calculate pole position:
       rad =45.D0/DATAN(1.D0)   ! 180/pi
       IF (.NOT.rank_has_l1m0) THEN
          CALL MPI_IRecv(b10,1,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE,&
               &sr_tag,MPI_COMM_WORLD,request1, ierr)
       END IF
       IF (.NOT.rank_has_l1m1) THEN
          CALL MPI_IRecv(b11,1,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE,&
               &sr_tag+1,MPI_COMM_WORLD,request2,ierr)
       END IF
       IF (.NOT.rank_has_l1m0) THEN
          CALL MPI_Wait(request1,status,ierr)
       END IF
       IF (.NOT.rank_has_l1m1) THEN
          CALL MPI_Wait(request2,status,ierr)
       END IF

       !print*, "------------", b11
       theta_dip= rad*DATAN2(DSQRT(2.D0)*ABS(b11),REAL(b10))
       IF ( theta_dip.LT.0.D0 ) theta_dip=180.D0+theta_dip
       phi_dip  =-rad*DATAN2(AIMAG(b11),REAL(b11))
       Dip      =e_dipole_ax/(e_p+e_t)
       DipCMB   =e_dipole_ax_cmb/e_cmb

       !-- Output of pole position:
       IF ( l_write ) THEN
          IF ( l_save_out ) THEN
             OPEN(n_dipole_file,FILE=dipole_file,STATUS='UNKNOWN',     &
                  &             POSITION='APPEND')
          END IF
          IF ( e_p_e.EQ.0 ) THEN
             e_p_e_ratio=0.D0
          ELSE
             e_p_e_ratio=e_dipole_e/e_p_e
          END IF
          !WRITE(*,"(A,4ES25.17)") "e_cmb = ",e_cmb,e_es_cmb,e_cmb-e_es_cmb,(e_cmb-e_es_cmb)/e_cmb
          ! There are still differences in field 17 of the dipole file. These
          ! differences are due to the summation for e_es_cmb and are only of the order
          ! of machine accuracy.
          WRITE(n_dipole_file,'(1P,D20.12,19D12.4)')                   &
               &       time*tScale,                            &! 1
               &       theta_dip,phi_dip,                      &! 2,3
               &       Dip,                                    &! 4  
               &       DipCMB,                                 &! 5
               &       e_dipole_ax_cmb/e_geo,                  &! 6
               &       e_dipole/(e_p+e_t),                     &! 7
               &       e_dip_cmb/e_cmb,                        &! 8
               &       e_dip_cmb/e_geo,                        &! 9
               &       e_dip_cmb,e_dipole_ax_cmb,              &! 10,11
               &       e_dipole,e_dipole_ax,                   &! 12,13
               &       e_cmb,e_geo,                            &! 14,15
               &       e_p_e_ratio,                            &! 16
               &       (e_cmb-e_es_cmb)/e_cmb,                 &! 17
               &       (e_cmb-e_as_cmb)/e_cmb,                 &! 18
               &       (e_geo-e_es_geo)/e_geo,                 &! 19
               &       (e_geo-e_as_geo)/e_geo        ! 20
          IF ( l_save_out ) CLOSE(n_dipole_file)
       END IF
       ! Store values needed for movie output:
       movieDipColat      =theta_dip
       movieDipLon        =phi_dip
       movieDipStrength   =e_dip_cmb/e_cmb
       movieDipStrengthGeo=e_dip_cmb/e_geo
    END IF


  END SUBROUTINE get_e_mag

  !----------------------------------------------------------------------------
END MODULE magnetic_energy
