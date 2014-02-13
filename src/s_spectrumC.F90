!$Id$
!********************************************************************
    subroutine spectrumC(time,n_spec,s,ds)
!********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!--------------------------------------------------------------------

!  calculates spectra of temperature and composition

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
    USE usefull, ONLY: cc2real,cc22real
    USE LMLoop_data,ONLY: llm,ulm
    USE integration, ONLY: rInt,rIntIC
    IMPLICIT NONE

!-- Input of variables:
    INTEGER :: n_spec     ! number of spectrum/call, file
    REAL(kind=8) :: time

!-- Input of fields:
    COMPLEX(kind=8) :: s(llm:ulm,n_r_max)
    COMPLEX(kind=8) :: ds(llm:ulm,n_r_max)

!-- output:
    REAL(kind=8) :: T_l(l_max+1)
    REAL(kind=8) :: T_m(l_max+1)
    REAL(kind=8),DIMENSION(l_max+1) :: T_ICB_l,T_ICB_l_global
    REAL(kind=8),DIMENSION(l_max+1) :: T_ICB_m,T_ICB_m_global
    REAL(kind=8),DIMENSION(l_max+1) :: dT_ICB_l,dT_ICB_l_global
    REAL(kind=8),DIMENSION(l_max+1) :: dT_ICB_m,dT_ICB_m_global

!-- local:
    CHARACTER(len=14) :: string
    CHARACTER(len=72) :: spec_file
    INTEGER :: n_r,lm,ml,l,mc,m,lc
    REAL(kind=8) :: T_temp
    REAL(kind=8) :: dT_temp
    REAL(kind=8) :: surf_ICB
    REAL(kind=8) :: fac,facICB

    REAL(kind=8),DIMENSION(n_r_max,l_max+1) :: T_r_l,T_r_l_global
    REAL(kind=8),DIMENSION(n_r_max,l_max+1) :: T_r_m,T_r_m_global

!-- end of declaration
!---------------------------------------------------------------------

    DO l=1,l_max+1
        T_l(l)=0.D0
        T_ICB_l(l)=0.D0
        dT_ICB_l(l)=0.D0
        T_m(l)=0.D0
        T_ICB_m(l)=0.D0
        dT_ICB_m(l)=0.D0
    END DO

    DO n_r=1,n_r_max

        DO l=1,l_max+1
            T_r_l(n_r,l)=0.D0
            T_ICB_l(l)  =0.D0
            dT_ICB_l(l) =0.D0
            T_r_m(n_r,l)=0.D0
            T_ICB_m(l)  =0.D0
            dT_ICB_m(l) =0.D0
        END DO

        !DO lm=1,lm_max
        DO lm=llm,ulm

            l =lo_map%lm2l(lm)
            m =lo_map%lm2m(lm)
            lc=l+1
            mc=m+1

            T_temp=DSQRT(cc2real(s(lm,n_r),m))/or2(n_r)
            dT_temp=DSQRT(cc2real(ds(lm,n_r),m))/or2(n_r)
        !----- l-spectra:
            T_r_l(n_r,lc) =T_r_l(n_r,lc) +  T_temp
        !----- m-spectra:
            T_r_m(n_r,mc)=T_r_m(n_r,mc) + T_temp

        !----- ICB spectra:
            IF ( n_r == n_r_icb ) THEN
                T_ICB_l(lc) =T_ICB_l(lc) +T_temp
                T_ICB_m(mc)=T_ICB_m(mc)+T_temp
                dT_ICB_l(lc) =dT_ICB_l(lc) +dT_temp
                dT_ICB_m(mc)=dT_ICB_m(mc)+dT_temp
            END IF

        END DO

    END DO

    ! reduction over all ranks
    CALL MPI_Reduce(T_r_l,T_r_l_global,n_r_max*(l_max+1),&
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(T_r_m,T_r_m_global,n_r_max*(l_max+1),&
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(T_ICB_l,T_ICB_l_global,l_max+1,&
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(T_ICB_m,T_ICB_m_global,l_max+1,&
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(dT_ICB_l,dT_ICB_l_global,l_max+1,&
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Reduce(dT_ICB_m,dT_ICB_m_global,l_max+1,&
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    IF (rank.EQ.0) THEN
       !-- Radial Integrals:
       surf_ICB=4.D0*pi*r_icb*r_icb
       fac      =1.D0/vol_oc
       facICB   =1.D0/surf_ICB
       DO l=1,l_max+1
          T_l(l)=fac*rInt(T_r_l_global(1,l),n_r_max,dr_fac, &
               i_costf_init,d_costf_init)
          T_ICB_l(l)=facICB*T_ICB_l_global(l)
          dT_ICB_l(l)=facICB*dT_ICB_l_global(l)
       END DO
       DO m=1,l_max+1 ! Note: counter m is actual order+1
          T_m(m)=fac*rInt(T_r_m_global(1,m),n_r_max,dr_fac, &
               i_costf_init,d_costf_init)
          T_ICB_m(m)=facICB*T_ICB_m_global(m)
          dT_ICB_m(m)=facICB*dT_ICB_m_global(m)
       END DO

       !-- Output into files:
       WRITE(string, *) n_spec
       spec_file='TC_spec_'//TRIM(ADJUSTL(string))//'.'//tag
       OPEN(98,file=spec_file,status='UNKNOWN')
       WRITE(98,'(1x,''TC spectra at time:'', &
            &    D20.12)') time*tScale
       DO ml=1,l_max+1
          WRITE(98,'(1P,I4,6D12.4)') &
               ml-1,T_l(ml),T_m(ml), &
               T_ICB_l(ml),T_ICB_m(ml), &
               dT_ICB_l(ml),dT_ICB_m(ml)
       END DO
       CLOSE(98)

    END IF

    RETURN
    end subroutine spectrumC

!----------------------------------------------------------------------------
