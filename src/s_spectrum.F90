!$Id$
!********************************************************************
SUBROUTINE spectrum(time,n_spec,w,dw,z, &
     &              b,db,aj,b_ic,db_ic,aj_ic)
  !********************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/17/02  by JW. --------------!

  !--------------------------------------------------------------------

  !  calculates magnetic energy  = 1/2 Integral(B^2 dV)
  !  integration in theta,phi by summation over harmonic coeffs.
  !  integration in r by Chebycheff integrals

  !  Output:
  !  enbp: Total poloidal        enbt: Total toroidal
  !  apome: Axisym. poloidal     atome: Axisym. toroidal

  !--------------------------------------------------------------------
  use parallel_mod
  USE truncation
  USE radial_data,ONLY: n_r_cmb
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE blocking
  USE horizontal_data
  USE logic
  USE output_data
  USE charmanip, ONLY: length_to_blank
  USE usefull, ONLY: cc2real,cc22real
  USE LMLoop_data,ONLY: llm,ulm,llmMag,ulmMag
  USE integration,ONLY: rInt_R,rIntIC
  IMPLICIT NONE

  !-- Input of variables:
  INTEGER,intent(IN) :: n_spec     ! number of spectrum/call, file
  REAL(kind=8),intent(IN) :: time

  !-- Input of fields:
  COMPLEX(kind=8),intent(IN) :: w(llm:ulm,n_r_max)
  COMPLEX(kind=8),intent(IN) :: dw(llm:ulm,n_r_max)
  COMPLEX(kind=8),intent(IN) :: z(llm:ulm,n_r_max)
  COMPLEX(kind=8),intent(IN) :: b(llmMag:ulmMag,n_r_maxMag)
  COMPLEX(kind=8),intent(IN) :: db(llmMag:ulmMag,n_r_maxMag)
  COMPLEX(kind=8),intent(IN) :: aj(llmMag:ulmMag,n_r_maxMag)
  COMPLEX(kind=8),intent(IN) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
  COMPLEX(kind=8),intent(IN) :: db_ic(llmMag:ulmMag,n_r_ic_maxMag)
  COMPLEX(kind=8),intent(IN) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)

  !-- Output:
  REAL(kind=8) :: b_rms
  REAL(kind=8) :: e_mag_p_l(l_max),e_mag_t_l(l_max)
  REAL(kind=8) :: e_kin_p_l(l_max),e_kin_t_l(l_max)
  REAL(kind=8) :: e_mag_p_ic_l(l_max),e_mag_t_ic_l(l_max)
  REAL(kind=8) :: u2_p_l(l_max),u2_t_l(l_max)

  REAL(kind=8) :: e_mag_p_m(l_max+1),e_mag_t_m(l_max+1)
  REAL(kind=8) :: e_kin_p_m(l_max+1),e_kin_t_m(l_max+1)
  REAL(kind=8) :: e_mag_p_ic_m(l_max+1),e_mag_t_ic_m(l_max+1)
  REAL(kind=8) :: u2_p_m(l_max+1),u2_t_m(l_max+1)

  REAL(kind=8) :: e_mag_cmb_l(l_max)
  REAL(kind=8) :: e_mag_cmb_m(l_max+1)
  REAL(kind=8) :: e_kin_nearSurf_l(l_max)
  REAL(kind=8) :: e_kin_nearSurf_m(l_max+1)

  REAL(kind=8) :: eCMB(l_max),eCMB_global(l_max)

  !-- local:
  CHARACTER(len=14) :: string
  CHARACTER(len=72) :: mag_spec_file,kin_spec_file,u2_spec_file
  INTEGER :: n_r,lm,ml,l,mc,m,n_const

  REAL(kind=8) :: r_ratio,O_r_icb_E_2
  REAL(kind=8) :: e_mag_p_temp,e_mag_t_temp
  REAL(kind=8) :: e_kin_p_temp,e_kin_t_temp
  REAL(kind=8) :: u2_p_temp,u2_t_temp
  REAL(kind=8) :: O_surface,pi
  REAL(kind=8) :: fac_mag,fac_kin
  REAL(kind=8) :: nearSurfR

  REAL(kind=8),DIMENSION(n_r_max,l_max) :: e_mag_p_r_l,e_mag_p_r_l_global
  REAL(kind=8),DIMENSION(n_r_max,l_max) :: e_mag_t_r_l,e_mag_t_r_l_global
  REAL(kind=8),DIMENSION(n_r_max,l_max) :: e_kin_p_r_l,e_kin_p_r_l_global
  REAL(kind=8),DIMENSION(n_r_max,l_max) :: e_kin_t_r_l,e_kin_t_r_l_global
  REAL(kind=8),DIMENSION(n_r_max,l_max) :: u2_p_r_l,u2_p_r_l_global
  REAL(kind=8),DIMENSION(n_r_max,l_max) :: u2_t_r_l,u2_t_r_l_global
  REAL(kind=8),DIMENSION(n_r_max,l_max+1) :: e_mag_p_r_m,e_mag_p_r_m_global
  REAL(kind=8),DIMENSION(n_r_max,l_max+1) :: e_mag_t_r_m,e_mag_t_r_m_global
  REAL(kind=8),DIMENSION(n_r_max,l_max+1) :: e_kin_p_r_m,e_kin_p_r_m_global
  REAL(kind=8),DIMENSION(n_r_max,l_max+1) :: e_kin_t_r_m,e_kin_t_r_m_global
  REAL(kind=8),DIMENSION(n_r_max,l_max+1) :: u2_p_r_m,u2_p_r_m_global
  REAL(kind=8),DIMENSION(n_r_max,l_max+1) :: u2_t_r_m,u2_t_r_m_global

  REAL(kind=8),DIMENSION(n_r_ic_max,l_max) :: e_mag_p_ic_r_l,e_mag_p_ic_r_l_global
  REAL(kind=8),DIMENSION(n_r_ic_max,l_max) :: e_mag_t_ic_r_l,e_mag_t_ic_r_l_global
  REAL(kind=8),DIMENSION(n_r_ic_max,l_max+1) :: e_mag_p_ic_r_m,e_mag_p_ic_r_m_global
  REAL(kind=8),DIMENSION(n_r_ic_max,l_max+1) :: e_mag_t_ic_r_m,e_mag_t_ic_r_m_global

  COMPLEX(kind=8) :: r_dr_b

  !-- end of declaration
  !---------------------------------------------------------------------

  pi=4.d0*DATAN(1.d0)
  eCMB=0.0d0

  DO n_r=1,n_r_max

     DO l=1,l_max
        IF ( l_mag ) THEN
           e_mag_p_r_l(n_r,l)=0.d0
           e_mag_t_r_l(n_r,l)=0.d0
        END IF
        IF ( l_anel ) THEN
           u2_p_r_l(n_r,l)=0.d0
           u2_t_r_l(n_r,l)=0.d0
        END IF
        e_kin_p_r_l(n_r,l)=0.d0
        e_kin_t_r_l(n_r,l)=0.d0
     END DO
     DO mc=1,l_max+1
        IF ( l_mag ) THEN
           e_mag_p_r_m(n_r,mc)=0.d0
           e_mag_t_r_m(n_r,mc)=0.d0
        END IF
        IF ( l_anel ) THEN
           u2_p_r_m(n_r,mc)=0.d0
           u2_t_r_m(n_r,mc)=0.d0
        END IF
        e_kin_p_r_m(n_r,mc)=0.d0
        e_kin_t_r_m(n_r,mc)=0.d0
     END DO

     !DO lm=2,lm_max
     DO lm=MAX(llm,2),ulm

        l  =lo_map%lm2l(lm)
        m  =lo_map%lm2m(lm)
        mc=m+1

        IF ( l_mag ) THEN
           e_mag_p_temp= dLh(st_map%lm2(l,m)) * ( &
                &          dLh(st_map%lm2(l,m))*or2(n_r)*cc2real(b(lm,n_r),m) + &
                &          cc2real(db(lm,n_r),m) )
           e_mag_t_temp=dLh(st_map%lm2(l,m))*cc2real(aj(lm,n_r),m)
        END IF
        IF ( l_anel ) THEN
           u2_p_temp=  orho2(n_r)*dLh(st_map%lm2(l,m)) *  ( &
                &        dLh(st_map%lm2(l,m))*or2(n_r)*cc2real(w(lm,n_r),m) + &
                &        cc2real(dw(lm,n_r),m) )
           u2_t_temp=orho2(n_r)*dLh(st_map%lm2(l,m))*cc2real(z(lm,n_r),m)
        END IF
        e_kin_p_temp= orho1(n_r)*dLh(st_map%lm2(l,m)) *  ( &
             &          dLh(st_map%lm2(l,m))*or2(n_r)*cc2real(w(lm,n_r),m) + &
             &          cc2real(dw(lm,n_r),m) )
        e_kin_t_temp=orho1(n_r)*dLh(st_map%lm2(l,m))*cc2real(z(lm,n_r),m)

        !----- l-spectra:
        IF ( l_mag ) THEN
           e_mag_p_r_l(n_r,l) = e_mag_p_r_l(n_r,l) + e_mag_p_temp
           e_mag_t_r_l(n_r,l) = e_mag_t_r_l(n_r,l) + e_mag_t_temp
           IF ( m == 0 .AND. n_r == n_r_cmb ) eCMB(l)=e_mag_p_temp
        END IF
        IF ( l_anel ) THEN
           u2_p_r_l(n_r,l) = u2_p_r_l(n_r,l) + u2_p_temp
           u2_t_r_l(n_r,l) = u2_t_r_l(n_r,l) + u2_t_temp
        END IF
        e_kin_p_r_l(n_r,l) = e_kin_p_r_l(n_r,l) + e_kin_p_temp
        e_kin_t_r_l(n_r,l) = e_kin_t_r_l(n_r,l) + e_kin_t_temp

        !----- m-spectra:
        IF ( l_mag ) THEN
           e_mag_p_r_m(n_r,mc) = e_mag_p_r_m(n_r,mc) + e_mag_p_temp
           e_mag_t_r_m(n_r,mc) = e_mag_t_r_m(n_r,mc) + e_mag_t_temp
        END IF
        IF ( l_anel ) THEN
           u2_p_r_m(n_r,mc) = u2_p_r_m(n_r,mc) + u2_p_temp
           u2_t_r_m(n_r,mc) = u2_t_r_m(n_r,mc) + u2_t_temp
        END IF
        e_kin_p_r_m(n_r,mc)=e_kin_p_r_m(n_r,mc) + e_kin_p_temp
        e_kin_t_r_m(n_r,mc)=e_kin_t_r_m(n_r,mc) + e_kin_t_temp

     END DO    ! do loop over lms in block

  END DO    ! radial grid points

  ! ----------- We need a reduction here ----------------
  ! first the l-spectra
  IF (l_mag) THEN
     CALL MPI_Reduce(e_mag_p_r_l, e_mag_p_r_l_global, n_r_max*l_max,&
          &MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     CALL MPI_Reduce(e_mag_t_r_l, e_mag_t_r_l_global, n_r_max*l_max,&
          &MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  END IF
  IF (l_anel) THEN
     CALL MPI_Reduce(u2_p_r_l, u2_p_r_l_global, n_r_max*l_max,&
          &MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     CALL MPI_Reduce(u2_t_r_l, u2_t_r_l_global, n_r_max*l_max,&
          &MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  END IF
  CALL MPI_Reduce(e_kin_p_r_l, e_kin_p_r_l_global, n_r_max*l_max,&
       &MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  CALL MPI_Reduce(e_kin_t_r_l, e_kin_t_r_l_global, n_r_max*l_max,&
       &MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

  ! then the m-spectra
  IF (l_mag) THEN
     CALL MPI_Reduce(e_mag_p_r_m, e_mag_p_r_m_global, n_r_max*(l_max+1),&
          &MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     CALL MPI_Reduce(e_mag_t_r_m, e_mag_t_r_m_global, n_r_max*(l_max+1),&
          &MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     CALL MPI_Reduce(eCMB, eCMB_global, l_max,&
          &MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  END IF
  IF (l_anel) THEN
     CALL MPI_Reduce(u2_p_r_m, u2_p_r_m_global, n_r_max*(l_max+1),&
          &MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     CALL MPI_Reduce(u2_t_r_m, u2_t_r_m_global, n_r_max*(l_max+1),&
          &MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  END IF
  CALL MPI_Reduce(e_kin_p_r_m, e_kin_p_r_m_global, n_r_max*(l_max+1),&
       &MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  CALL MPI_Reduce(e_kin_t_r_m, e_kin_t_r_m_global, n_r_max*(l_max+1),&
       &MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

  ! now switch to rank 0 for the postprocess
  

  ! Getting appropriate radius index for e_kin_nearSurf spectra
  nearSurfR = r_icb+0.99
  do n_r=2,n_r_max
     if ( r(n_r-1) > nearSurfR .AND. &
          r(n_r)  <= nearSurfR ) then
        if ( r(n_r-1)-nearSurfR < &
             nearSurfR-r(n_r) ) then
           n_const=n_r-1
        else
           n_const=n_r
        end if
     end if
  end do

  IF (rank.EQ.0) THEN
     !-- Save CMB energy spectra:
     O_surface=1.d0/(4.d0*pi*r(1)*r(1))

     IF ( l_mag ) THEN
        b_rms=0.d0
        DO l=1,l_max
           e_mag_cmb_l(l)=e_mag_p_r_l_global(1,l)
           b_rms=b_rms + e_mag_cmb_l(l)
        END DO
        b_rms=dsqrt(b_rms*O_surface)
        DO mc=1,l_max+1
           e_mag_cmb_m(mc)=e_mag_p_r_m_global(1,mc)
        END DO
     END IF

     !-- Save nearSurf kin energy spectra:
     DO l=1,l_max
        e_kin_nearSurf_l(l)=e_kin_p_r_l_global(n_const,l)
     END DO
     DO mc=1,l_max+1
        e_kin_nearSurf_m(mc)=e_kin_p_r_m_global(n_const,mc)
     END DO


     !-- Radial Integrals:
     fac_mag=0.5*LFfac*eScale
     fac_kin=0.5*eScale
     DO l=1,l_max
        IF ( l_mag ) THEN
           e_mag_p_l(l)=                                            &
                fac_mag*rInt_R(e_mag_p_r_l_global(1,l),n_r_max,n_r_max,drx, &
                i_costf_init,d_costf_init)
           e_mag_t_l(l)=                                            &
                fac_mag*rInt_R(e_mag_t_r_l_global(1,l),n_r_max,n_r_max,drx, &
                i_costf_init,d_costf_init)
           e_mag_cmb_l(l)=fac_mag*e_mag_cmb_l(l)
        END IF
        IF ( l_anel ) THEN
           u2_p_l(l)  =                                          &
                fac_kin*rInt_R(u2_p_r_l_global(1,l),n_r_max,n_r_max,drx, &
                i_costf_init,d_costf_init)
           u2_t_l(l)  =                                          &
                fac_kin*rInt_R(u2_t_r_l_global(1,l),n_r_max,n_r_max,drx, &
                i_costf_init,d_costf_init)

        END IF
        e_kin_p_l(l)  =                                          &
             fac_kin*rInt_R(e_kin_p_r_l_global(1,l),n_r_max,n_r_max,drx, &
             i_costf_init,d_costf_init)
        e_kin_t_l(l)  =                                          &
             fac_kin*rInt_R(e_kin_t_r_l_global(1,l),n_r_max,n_r_max,drx, &
             i_costf_init,d_costf_init)
        e_kin_nearSurf_l(l)=fac_kin*e_kin_nearSurf_l(l)
     END DO
     DO m=1,l_max+1 ! Note: counter m is actual order+1
        IF ( l_mag )  THEN
           e_mag_p_m(m)=                                            &
                fac_mag*rInt_R(e_mag_p_r_m_global(1,m),n_r_max,n_r_max,drx, &
                i_costf_init,d_costf_init)
           e_mag_t_m(m)=                                            &
                fac_mag*rInt_R(e_mag_t_r_m_global(1,m),n_r_max,n_r_max,drx, &
                i_costf_init,d_costf_init)
           e_mag_cmb_m(m)=fac_mag*e_mag_cmb_m(m)
        END IF
        IF ( l_anel ) THEN
           u2_p_m(m)  =                                          &
                fac_kin*rInt_R(u2_p_r_m_global(1,m),n_r_max,n_r_max,drx, &
                i_costf_init,d_costf_init)
           u2_t_m(m)  =                                          &
                fac_kin*rInt_R(u2_t_r_m_global(1,m),n_r_max,n_r_max,drx, &
                i_costf_init,d_costf_init)
        END IF
        e_kin_p_m(m)  =                                          &
             fac_kin*rInt_R(e_kin_p_r_m_global(1,m),n_r_max,n_r_max,drx, &
             i_costf_init,d_costf_init)
        e_kin_t_m(m)  =                                          &
             fac_kin*rInt_R(e_kin_t_r_m_global(1,m),n_r_max,n_r_max,drx, &
             i_costf_init,d_costf_init)
        e_kin_nearSurf_m(m)=fac_kin*e_kin_nearSurf_m(m)
     END DO
  END IF

  !-- inner core:

  IF ( l_cond_ic ) THEN

     O_r_icb_E_2=1.d0/(r_ic(1)*r_ic(1))

     DO n_r=1,n_r_ic_max

        r_ratio=r_ic(n_r)/r_ic(1)

        DO mc=1,l_max+1
           e_mag_p_ic_r_m(n_r,mc)=0.d0
           e_mag_t_ic_r_m(n_r,mc)=0.d0
        END DO
        DO l=1,l_max
           e_mag_p_ic_r_l(n_r,l)=0.d0
           e_mag_t_ic_r_l(n_r,l)=0.d0
        END DO

        !DO lm=2,lm_max
        DO lm=MAX(llm,2),ulm

           l =lo_map%lm2l(lm)
           m =lo_map%lm2m(lm)
           mc=m+1
           r_dr_b=r_ic(n_r)*db_ic(lm,n_r)

           e_mag_p_temp=                                        &
                dLh(st_map%lm2(l,m))*O_r_icb_E_2*r_ratio**(2*l) * ( &
                DBLE((2*l+1)*(l+1))*cc2real(b_ic(lm,n_r),m)        + &
                DBLE(2*(l+1))*cc22real(b_ic(lm,n_r),r_dr_b,m) + &
                cc2real(r_dr_b,m) )
           e_mag_t_temp= dLh(st_map%lm2(l,m))*r_ratio**(2*l+2) * &
                cc2real(aj_ic(lm,n_r),m)

           e_mag_p_ic_r_l(n_r,l)=e_mag_p_ic_r_l(n_r,l) + &
                e_mag_p_temp
           e_mag_t_ic_r_l(n_r,l)=e_mag_t_ic_r_l(n_r,l) + &
                e_mag_t_temp
           e_mag_p_ic_r_m(n_r,mc)=e_mag_p_ic_r_m(n_r,mc) + &
                e_mag_p_temp
           e_mag_t_ic_r_m(n_r,mc)=e_mag_t_ic_r_m(n_r,mc) + &
                e_mag_t_temp

        END DO  ! loop over lm's

     END DO ! loop over radial levels

     CALL MPI_Reduce(e_mag_p_ic_r_l, e_mag_p_ic_r_l_global, n_r_ic_max*l_max,&
          &MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     CALL MPI_Reduce(e_mag_t_ic_r_l, e_mag_t_ic_r_l_global, n_r_ic_max*l_max,&
          &MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     CALL MPI_Reduce(e_mag_p_ic_r_m, e_mag_p_ic_r_m_global, n_r_ic_max*(l_max+1),&
          &MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     CALL MPI_Reduce(e_mag_t_ic_r_m, e_mag_t_ic_r_m_global, n_r_ic_max*(l_max+1),&
          &MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)


     IF (rank.EQ.0) THEN
        !----- Radial Integrals:
        fac_mag=LFfac*eScale/2.D0
        DO l=1,l_max
           e_mag_p_ic_l(l)=fac_mag*rIntIC(e_mag_p_ic_r_l_global(1,l), &
                n_r_ic_max,dr_fac_ic,i_costf1_ic_init,d_costf1_ic_init)
           e_mag_t_ic_l(l)=fac_mag*rIntIC(e_mag_t_ic_r_l_global(1,l), &
                n_r_ic_max,dr_fac_ic,i_costf1_ic_init,d_costf1_ic_init)
        END DO
        DO m=1,l_max+1
           e_mag_p_ic_m(m)=fac_mag*rIntIC(e_mag_p_ic_r_m_global(1,m), &
                n_r_ic_max,dr_fac_ic,i_costf1_ic_init,d_costf1_ic_init)
           e_mag_t_ic_m(m)=fac_mag*rIntIC(e_mag_t_ic_r_m_global(1,m), &
                n_r_ic_max,dr_fac_ic,i_costf1_ic_init,d_costf1_ic_init)
        END DO
     END IF
  ELSE
     DO l=1,l_max
        e_mag_p_ic_l(l)=0.d0
        e_mag_t_ic_l(l)=0.d0
     END DO
     DO mc=1,l_max+1
        e_mag_p_ic_m(mc)=0.d0
        e_mag_t_ic_m(mc)=0.d0
     END DO
  END IF  ! conducting inner core ?

  IF (rank.EQ.0) THEN
     !-- Output into files:
     IF ( l_mag ) THEN
        WRITE(string, *) n_spec
        mag_spec_file='mag_spec_'//TRIM(ADJUSTL(string))//'.'//tag
        OPEN(n_mag_spec_file,FILE=mag_spec_file,STATUS='UNKNOWN')
        IF ( n_spec == 0 ) THEN
           WRITE(n_mag_spec_file,'(1x, &
                &      ''Magnetic energy spectra of time averaged field:'')')
        ELSE
           WRITE(n_mag_spec_file,'(1x, &
                &      ''Magnetic energy spectra at time:'', &
                &      d20.12)') time*tScale
        END IF
        WRITE(n_mag_spec_file,'(1p,i4,11d12.4)')         &
             0,0.d0,e_mag_p_m(1)   ,0.d0,e_mag_t_m(1), &
             0.d0,e_mag_p_ic_m(1),0.d0,e_mag_t_ic_m(1), &
             0.d0,e_mag_cmb_m(1),0.d0
        DO ml=1,l_max
           WRITE(n_mag_spec_file,'(1p,i4,11d12.4)') &
                ml,e_mag_p_l(ml),   e_mag_p_m(ml+1), &
                e_mag_t_l(ml),   e_mag_t_m(ml+1), &
                e_mag_p_ic_l(ml),e_mag_p_ic_m(ml+1), &
                e_mag_t_ic_l(ml),e_mag_t_ic_m(ml+1), &
                e_mag_cmb_l(ml), e_mag_cmb_m(ml+1), &
                eCMB_global(ml)
        END DO
        CLOSE(n_mag_spec_file)
     END IF

     IF ( l_anel ) THEN
        WRITE(string, *) n_spec
        u2_spec_file='u2_spec_'//trim(adjustl(string))//'.'//tag
        OPEN(n_u2_spec_file,FILE=u2_spec_file,STATUS='UNKNOWN')
        IF ( n_spec == 0 ) THEN
           WRITE(n_u2_spec_file,'(1x, &
                &     ''Velocity square spectra of time averaged field:'')')
        ELSE
           WRITE(n_u2_spec_file,'(1x,                &
                &     ''Velocity square spectra at time:'', &
                &     d20.12)') time*tScale
        END IF
        WRITE(n_u2_spec_file,'(1p,i4,4d12.4)') &
             &   0,0.d0,u2_p_m(1),0.d0,u2_t_m(1)
        DO ml=1,l_max
           WRITE(n_u2_spec_file,'(1p,i4,4d12.4)') &
                &       ml,u2_p_l(ml),u2_p_m(ml+1),u2_t_l(ml),u2_t_m(ml+1)
        END DO
        CLOSE(n_u2_spec_file)
     END IF

     WRITE(string, *) n_spec
     kin_spec_file='kin_spec_'//TRIM(ADJUSTL(string))//'.'//tag
     OPEN(n_kin_spec_file,FILE=kin_spec_file,STATUS='UNKNOWN')
     IF ( n_spec == 0 ) THEN
        WRITE(n_kin_spec_file,'(1x, &
             &      ''Kinetic energy spectra of time averaged field:'')')
     ELSE
        WRITE(n_kin_spec_file,'(1x,                &
             &      ''Kinetic energy spectra at time:'', &
             &      d20.12)') time*tScale
     END IF
     WRITE(n_kin_spec_file,'(1p,i4,6d12.4)')        &
          0,0.d0,e_kin_p_m(1),0.d0,e_kin_t_m(1), &
          0.d0, e_kin_nearSurf_m(1)
     DO ml=1,l_max
        WRITE(n_kin_spec_file,'(1p,i4,6d12.4)')    &
             ml,e_kin_p_l(ml),e_kin_p_m(ml+1), &
             e_kin_t_l(ml),e_kin_t_m(ml+1), &
             e_kin_nearSurf_l(ml), e_kin_nearSurf_m(ml+1)
     END DO
     CLOSE(n_kin_spec_file)
  END IF

  RETURN
end SUBROUTINE spectrum

!----------------------------------------------------------------------------
