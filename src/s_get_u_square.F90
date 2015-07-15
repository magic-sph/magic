! Id: s_get_u_square.f 400 2013-02-13 12:59:49Z gastine $
!********************************************************************
    SUBROUTINE get_u_square(time,w,dw,z,RolR,dlR,dlRc)
!********************************************************************

!--------------------------------------------------------------------

!  calculates square velocity  = 1/2 Integral (v^2 dV)
!  integration in theta,phi by summation of spherical harmonics
!  integration in r by using Chebychef integrals

!  Write the different contributions in u_square.TAG file

!--------------------------------------------------------------------

    USE truncation
    use parallel_mod
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE blocking
    USE horizontal_data
    USE logic
    USE output_data
    USE const
    USE usefull, ONLY: cc2real
    USE integration, ONLY: rInt_R
    USE LMLoop_data,ONLY: llm,ulm
    IMPLICIT NONE

!-- Input of scalar fields:
    REAL(kind=8),intent(IN) :: time
    COMPLEX(kind=8),intent(IN) :: w(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(IN) :: dw(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(IN) :: z(llm:ulm,n_r_max)

!-- Output:
    REAL(kind=8),intent(OUT) :: dlR(n_r_max)
    REAL(kind=8),intent(OUT) :: dlRc(n_r_max)
    REAL(kind=8),intent(OUT) :: RolR(n_r_max)

!-- u**2
    REAL(kind=8) :: e_p     ! poloidal u**2
    REAL(kind=8) :: e_t     ! toroidal u**2
    REAL(kind=8) :: e_p_as  ! axisymmetric poloidal u**2
    REAL(kind=8) :: e_t_as  ! axisymmetric toroidal u**2
    REAL(kind=8) :: e_kin   ! total u**2

!-- local:
    REAL(kind=8) :: e_p_temp,e_t_temp
    REAL(kind=8),DIMENSION(n_r_max) :: e_p_r,e_p_r_global
    REAL(kind=8),DIMENSION(n_r_max) :: e_t_r,e_t_r_global
    REAL(kind=8),DIMENSION(n_r_max) :: e_p_as_r,e_p_as_r_global
    REAL(kind=8),DIMENSION(n_r_max) :: e_t_as_r,e_t_as_r_global
    REAL(kind=8),DIMENSION(n_r_max,l_max) :: e_lr, e_lr_global,e_lr_c,e_lr_c_global
    REAL(kind=8) :: ER(n_r_max),ELR(n_r_max),ReR(n_r_max),RoR(n_r_max)
    REAL(kind=8) :: ERc(n_r_max),ELRc(n_r_max)
    REAL(kind=8) :: ekinR(n_r_max)
    REAL(kind=8) :: RmR(n_r_max)
    REAL(kind=8) :: e_l,E,EL,Ec,ELc

    INTEGER :: nR,lm,l,m
    REAL(kind=8) :: fac
    REAL(kind=8) :: O_rho ! 1/rho**2 (anelastic)

!-- property parameters
    REAL(kind=8) :: Re,Rm,Ro,Rol,dl,dlc
    REAL(kind=8) :: ReConv,RoConv,RolC


!-- end of declaration
!---------------------------------------------------------------------
     
    DO nR=1,n_r_max
        e_p_r(nR)    =0.D0
        e_t_r(nR)    =0.D0
        e_p_as_r(nR) =0.D0
        e_t_as_r(nR) =0.D0
        O_rho        =orho2(nR) ! divided by rho**2
        DO l=1,l_max
            e_lr(nR,l)=0.D0
            e_lr_c(nR,l)=0.D0
        END DO

        DO lm=MAX(2,llm),ulm
           l=lo_map%lm2l(lm)
           m=lo_map%lm2m(lm)
           !DO lm=2,lm_max
           !  l=lm2l(lm)
           !  m=lm2m(lm)

            e_p_temp= O_rho*dLh(st_map%lm2(l,m)) * ( &
                 &      dLh(st_map%lm2(l,m))*or2(nR)*cc2real(w(lm,nR),m) &
                 &      + cc2real(dw(lm,nR),m) )
            e_t_temp= O_rho*dLh(st_map%lm2(l,m))*cc2real(z(lm,nR),m)
            IF ( m == 0 ) THEN  ! axisymmetric part
                e_p_as_r(nR)=e_p_as_r(nR)+ e_p_temp
                e_t_as_r(nR)=e_t_as_r(nR)+ e_t_temp
            ELSE
                e_p_r(nR)=e_p_r(nR) + e_p_temp
                e_t_r(nR)=e_t_r(nR) + e_t_temp
                e_lr_c(nR,l)=e_lr_c(nR,l) + e_p_temp + e_t_temp
            END IF

            e_lr(nR,l)=e_lr(nR,l) + e_p_temp + e_t_temp
        END DO    ! do loop over lms in block
        e_p_r(nR)=e_p_r(nR)+e_p_as_r(nR)
        e_t_r(nR)=e_t_r(nR)+e_t_as_r(nR)
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
    CALL MPI_Reduce(e_lr_c, e_lr_c_global,  n_r_max*l_max, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(e_lr, e_lr_global,  n_r_max*l_max, &
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    IF (rank == 0) THEN
       !-- Radial Integrals:
       e_p    =rInt_R(e_p_r_global,    n_r_max,n_r_max,drx, &
            i_costf_init,d_costf_init)
       e_t    =rInt_R(e_t_r_global,    n_r_max,n_r_max,drx, &
            i_costf_init,d_costf_init)
       e_p_as =rInt_R(e_p_as_r_global, n_r_max,n_r_max,drx, &
            i_costf_init,d_costf_init)
       e_t_as =rInt_R(e_t_as_r_global, n_r_max,n_r_max,drx, &
            i_costf_init,d_costf_init)
       fac    =0.5D0*eScale
       e_p    =fac*e_p
       e_t    =fac*e_t
       e_p_as =fac*e_p_as
       e_t_as =fac*e_t_as

       e_kin  =e_t+e_p
       DO nR=1,n_r_max
          ekinR(nR)=fac*(e_p_r_global(nR)+e_t_r_global(nR))
       END DO

       !-- Rossby number
       Re=DSQRT(2.D0*e_kin/vol_oc)
       ReConv=DSQRT(2.D0*(e_kin-e_p_as-e_t_as)/vol_oc)
       IF ( l_non_rot ) THEN
          Ro=0.D0
          RoConv=0.D0
       ELSE
          Ro=Re*ek
          RoConv=ReConv*ek
       END IF

       !-- Length Scale
       E  =0.D0
       EL =0.D0
       Ec =0.D0
       ELc=0.D0
       DO l=1,l_max
          e_l=fac*rInt_R(e_lr_global(1,l),n_r_max,n_r_max,drx, &
               i_costf_init,d_costf_init)
          E =E+e_l
          EL=EL+DBLE(l)*e_l
          e_l=fac*rInt_R(e_lr_c_global(1,l),n_r_max,n_r_max,drx, &
               i_costf_init,d_costf_init)
          Ec =Ec+e_l
          ELc=ELc+DBLE(l)*e_l
       END DO
       IF ( EL /= 0d0 ) THEN
          dl=pi*E/EL
          dlc=pi*Ec/ELc
       ELSE
          dl=0d0
          dlc=0d0
       END IF
       DO nR=1,n_r_max
          ER(nR)  =0.D0
          ELR(nR) =0.D0
          ERc(nR) =0.D0
          ELRc(nR)=0.D0
          DO l=1,l_max
             e_l=fac*e_lr_global(nR,l)
             ER(nR) =ER(nR)+e_l
             ELR(nR)=ELR(nR)+DBLE(l)*e_l
             e_l=fac*e_lr_c_global(nR,l)
             ERc(nR) =ERc(nR)+e_l
             ELRc(nR)=ELRc(nR)+DBLE(l)*e_l
          END DO
          IF ( ELR(nR) /= 0d0 ) THEN
             dlR(nR)=pi*ER(nR)/ELR(nR)
             dlRc(nR)=pi*ERc(nR)/ELRc(nR)
          ELSE
             dlR(nR)=0d0
             dlRc(nR)=0d0
          END IF
       END DO

       !-- Local Rossby number
       IF ( dl/=0d0 ) THEN
          Rol = Ro/dl
          RolC = RoConv/dlc
       ELSE
          Rol = Ro
          RolC = RoConv
       END IF
       DO nR=1,n_r_max
          ReR(nR)=DSQRT(2.D0*ekinR(nR)*or2(nR)/(4*pi))
          RoR(nR)=ReR(nR)*ek
          IF ( dlR(nR) /= 0d0 ) THEN
             RolR(nR)=RoR(nR)/dlR(nR)
          ELSE
             RolR(nR)=RoR(nR)
          END IF
          RmR(nR)=ReR(nR)*prmag*sigma(nR)*r(nR)*r(nR)
       END DO

       !-- Magnetic reynolds number
       IF ( prmag /= 0 .AND. nVarCond > 0 ) THEN
          Rm=0.d0
          Rm=rInt_R(RmR,n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
          Rm=Rm*3/(r_cmb**3-r_icb**3)
       ELSEIF ( prmag /= 0 ) THEN
          Rm=Re*prmag
       ELSE
          Rm=Re
       END IF

       !-- Output
       IF ( l_save_out ) THEN
          OPEN(n_u_square_file,FILE=u_square_file,STATUS='UNKNOWN', &
               POSITION='APPEND')
       END IF
       WRITE(n_u_square_file,'(1P,D20.12,10D16.8)') &
            &  time*tScale,     & ! 1
            &      e_p,e_t,     & ! 2,3
            &e_p_as,e_t_as,     & ! 4,5
            &        Ro,Rm,     & ! 6,7
            &       Rol,dl,     & ! 8,9
            &     RolC,dlc        ! 10,11
       IF ( l_save_out ) CLOSE(n_u_square_file)
    END IF
    RETURN
    end SUBROUTINE get_u_square
!-----------------------------------------------------------------------------
