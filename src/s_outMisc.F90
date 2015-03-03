!$Id$
!*************************************************************************
SUBROUTINE outMisc(timeScaled,HelLMr,Hel2LMr,HelnaLMr,Helna2LMr, &
     &             nLogs,w,dw,ddw,z,dz,s,ds,p,Geos,dpFlow,dzFlow)
  !*************************************************************************

  use mpi
  USE truncation
  USE radial_functions, only: botcond
  USE physical_parameters
  USE num_param
  USE blocking
  USE horizontal_data
  USE logic
  USE output_data
  use Egeos_mod
  USE const
  USE parallel_mod,only: nr_per_rank
  USE usefull, ONLY: cc2real
  USE integration, ONLY: rInt,rInt_R
  USE LMLoop_data,ONLY: llm,ulm
  IMPLICIT NONE

  !-- Input of variables:
  REAL(kind=8),INTENT(IN) :: timeScaled
  REAL(kind=8),INTENT(IN) :: HelLMr(l_max+1,nRstart:nRstop)
  REAL(kind=8),INTENT(IN) :: Hel2LMr(l_max+1,nRstart:nRstop)
  REAL(kind=8),INTENT(IN) :: HelnaLMr(l_max+1,nRstart:nRstop)
  REAL(kind=8),INTENT(IN) :: Helna2LMr(l_max+1,nRstart:nRstop)
  INTEGER,INTENT(IN) :: nLogs

  !-- Input of scalar fields:
  COMPLEX(kind=8),INTENT(IN) :: s(llm:ulm,n_r_max)
  COMPLEX(kind=8),INTENT(IN) :: ds(llm:ulm,n_r_max)
  !---- Fields transfered to getEgeos:
  COMPLEX(kind=8),INTENT(IN) :: w(llm:ulm,n_r_max)
  COMPLEX(kind=8),INTENT(IN) :: dw(llm:ulm,n_r_max)
  COMPLEX(kind=8),INTENT(IN) :: ddw(llm:ulm,n_r_max)
  COMPLEX(kind=8),INTENT(IN) :: z(llm:ulm,n_r_max)
  COMPLEX(kind=8),INTENT(IN) :: dz(llm:ulm,n_r_max)
  COMPLEX(kind=8),INTENT(IN) :: p(llm:ulm,n_r_max)

  !-- Output: (and stuff written in misc.TAG files)
  REAL(kind=8),INTENT(OUT) :: Geos
  REAL(kind=8),INTENT(OUT) :: dpFlow,dzFlow

  !-- Local stuff:
  INTEGER :: nTheta,nThetaStart,nThetaBlock,nThetaNHS,n,lm44
  logical :: lm44_is_local
  REAL(kind=8) :: pplot_global(n_r_max),pplot(n_r_max)
  REAL(kind=8),DIMENSION(nRstart:nRstop) :: HelNr,HelSr,HelnaNr,HelnaSr
  REAL(kind=8),DIMENSION(nRstart:nRstop) :: Helna2Nr, Helna2Sr,Hel2Nr, Hel2Sr
  REAL(kind=8),DIMENSION(nRstart:nRstop) :: HelEAr
  REAL(kind=8),DIMENSION(n_r_max) :: HelNr_global,HelSr_global
  REAL(kind=8),DIMENSION(n_r_max) :: HelnaNr_global,HelnaSr_global
  REAL(kind=8),DIMENSION(n_r_max) :: Helna2Nr_global,Helna2Sr_global
  REAL(kind=8),DIMENSION(n_r_max) :: Hel2Nr_global,Hel2Sr_global,HelEAr_global
  COMPLEX(kind=8),dimension(n_r_max) :: p44_local
  REAL(kind=8) :: Hel(nfs),Hel2(nfs),Helna(nfs),Helna2(nfs),r2
  REAL(kind=8) :: HelN,HelS
  REAL(kind=8) :: HelnaN,HelnaS
  REAL(kind=8) :: HelnaRMSN,HelnaRMSS
  REAL(kind=8) :: HelRMSN,HelRMSS,HelEA,HelRMS,HelnaRMS
  REAL(kind=8) :: Egeos,EkNTC,EkSTC,Ekin
  REAL(kind=8) :: CVzOTC,CVorOTC,CHelOTC
  REAL(kind=8) :: topnuss,botnuss
  REAL(kind=8) :: topflux,botflux

  INTEGER :: n_r,m,lm,mytag,status(MPI_STATUS_SIZE)
  REAL(kind=8) :: osq4pi
  INTEGER :: i,sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1),ierr

  CHARACTER(len=76) :: filename2

  !-- end of declaration
  !---------------------------------------------------------------------


  IF ( l_hel )  THEN

     !------ Integration of Helicity, on input the Helicity is
     !       already axisymmetric !
     DO n_r=nRstart,nRstop
        r2=r(n_r)*r(n_r)
        HelNr(n_r) =0.D0
        HelSr(n_r) =0.D0
        HelnaNr(n_r) =0.D0
        HelnaSr(n_r) =0.D0
        HelEAr(n_r)=0.D0
        Hel2Nr(n_r) =0.D0
        Hel2Sr(n_r) =0.D0
        Helna2Nr(n_r) =0.D0
        Helna2Sr(n_r) =0.D0

        DO n=1,nThetaBs ! Loop over theta blocks
           nTheta=(n-1)*sizeThetaB
           nThetaStart=nTheta+1
           CALL lmAS2pt(HelLMr(1,n_r),Hel,       &
                &       nThetaStart,sizeThetaB)
           CALL lmAS2pt(Hel2LMr(1,n_r),Hel2,     &
                &       nThetaStart,sizeThetaB)
           CALL lmAS2pt(HelnaLMr(1,n_r),Helna,   &
                &       nThetaStart,sizeThetaB)
           CALL lmAS2pt(Helna2LMr(1,n_r),Helna2, &
                &       nThetaStart,sizeThetaB)
           DO nThetaBlock=1,sizeThetaB
              nTheta=nTheta+1
              nThetaNHS=(nTheta+1)/2

              !------ Integration over theta:
              IF ( MOD(nTheta,2) == 1 ) THEN ! NHS
                 Hel2Nr(n_r)=Hel2Nr(n_r)+gauss(nThetaNHS)* &
                      r2*Hel2(nThetaBlock)
                 Helna2Nr(n_r)=Helna2Nr(n_r)+gauss(nThetaNHS)* &
                      r2*Helna2(nThetaBlock)
                 HelEAr(n_r)=HelEAr(n_r)+gauss(nThetaNHS)* &
                      r2*Hel(nThetaBlock)
                 HelNr(n_r) =HelNr(n_r)+gauss(nThetaNHS)* &
                      r2*Hel(nThetaBlock)
                 HelnaNr(n_r) =HelnaNr(n_r)+gauss(nThetaNHS)* &
                      r2*Helna(nThetaBlock)
              ELSE
                 Hel2Sr(n_r)=Hel2Sr(n_r)+gauss(nThetaNHS)* &
                      r2*Hel2(nThetaBlock)
                 Helna2Sr(n_r)=Helna2Sr(n_r)+gauss(nThetaNHS)* &
                      r2*Helna2(nThetaBlock)
                 HelEAr(n_r)=HelEAr(n_r)-gauss(nThetaNHS)* &
                      r2*Hel(nThetaBlock)
                 HelSr(n_r) =HelSr(n_r)+gauss(nThetaNHS)* &
                      r2*Hel(nThetaBlock)
                 HelnaSr(n_r)=HelnaSr(n_r)+gauss(nThetaNHS)* &
                      r2*Helna(nThetaBlock)
              END IF
           END DO
        END DO

     END DO

     ! Now we have to gather the results on rank 0 for
     ! the arrays: Hel2Nr,Helna2Nr,HelEAr,HelNr,HelnaNr
     ! Hel2Sr,Helna2Sr,HelSr,HelnaSr

     sendcount  = (nRstop-nRstart+1)
     recvcounts = nr_per_rank
     recvcounts(n_procs-1) = (nr_per_rank+1)
     DO i=0,n_procs-1
        displs(i) = i*nr_per_rank
     END DO
     CALL MPI_GatherV(Hel2Nr,sendcount,MPI_DOUBLE_PRECISION,&
          &           Hel2Nr_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
          &           0,MPI_COMM_WORLD,ierr)
     CALL MPI_GatherV(Helna2Nr,sendcount,MPI_DOUBLE_PRECISION,&
          &           Helna2Nr_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
          &           0,MPI_COMM_WORLD,ierr)
     CALL MPI_GatherV(HelEAr,sendcount,MPI_DOUBLE_PRECISION,&
          &           HelEAr_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
          &           0,MPI_COMM_WORLD,ierr)
     CALL MPI_GatherV(HelNr,sendcount,MPI_DOUBLE_PRECISION,&
          &           HelNr_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
          &           0,MPI_COMM_WORLD,ierr)
     CALL MPI_GatherV(HelnaNr,sendcount,MPI_DOUBLE_PRECISION,&
          &           HelnaNr_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
          &           0,MPI_COMM_WORLD,ierr)
     CALL MPI_GatherV(HelSr,sendcount,MPI_DOUBLE_PRECISION,&
          &           HelSr_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
          &           0,MPI_COMM_WORLD,ierr)
     CALL MPI_GatherV(Helna2Sr,sendcount,MPI_DOUBLE_PRECISION,&
          &           Helna2Sr_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
          &           0,MPI_COMM_WORLD,ierr)
     CALL MPI_GatherV(Hel2Sr,sendcount,MPI_DOUBLE_PRECISION,&
          &           Hel2Sr_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
          &           0,MPI_COMM_WORLD,ierr)
     CALL MPI_GatherV(HelnaSr,sendcount,MPI_DOUBLE_PRECISION,&
          &           HelnaSr_global,recvcounts,displs,MPI_DOUBLE_PRECISION,&
          &           0,MPI_COMM_WORLD,ierr)

     IF (rank.EQ.0) THEN
        !------ Integration over r without the boundaries and normalization:
        HelN  =rInt(HelNr_global,n_r_max,dr_fac,i_costf_init,d_costf_init)
        HelS  =rInt(HelSr_global,n_r_max,dr_fac,i_costf_init,d_costf_init)
        HelnaN=rInt(HelnaNr_global,n_r_max,dr_fac,i_costf_init,d_costf_init)
        HelnaS=rInt(HelnaSr_global,n_r_max,dr_fac,i_costf_init,d_costf_init)
        HelEA =rInt(HelEAr_global,n_r_max,dr_fac,i_costf_init,d_costf_init)
        HelRMSN=rInt(Hel2Nr_global,n_r_max,dr_fac,i_costf_init,d_costf_init)
        HelRMSS=rInt(Hel2Sr_global,n_r_max,dr_fac,i_costf_init,d_costf_init)
        HelnaRMSN=rInt(Helna2Nr_global,n_r_max,dr_fac,i_costf_init,d_costf_init)
        HelnaRMSS=rInt(Helna2Sr_global,n_r_max,dr_fac,i_costf_init,d_costf_init)

        HelN  =2.D0*pi*HelN/(vol_oc/2) ! Note integrated over half spheres only !
        HelS  =2.D0*pi*HelS/(vol_oc/2) ! Factor 2*pi is from phi integration
        HelnaN  =2.D0*pi*HelnaN/(vol_oc/2) ! Note integrated over half spheres only !
        HelnaS  =2.D0*pi*HelnaS/(vol_oc/2) ! Factor 2*pi is from phi integration
        HelEA =2.D0*pi*HelEA/vol_oc
        HelRMSN=DSQRT(2.D0*pi*HelRMSN/(vol_oc/2))
        HelRMSS=DSQRT(2.D0*pi*HelRMSS/(vol_oc/2))
        HelnaRMSN=DSQRT(2.D0*pi*HelnaRMSN/(vol_oc/2))
        HelnaRMSS=DSQRT(2.D0*pi*HelnaRMSS/(vol_oc/2))
        HelRMS=HelRMSN+HelRMSS
        HelnaRMS=HelnaRMSN+HelnaRMSS

        IF ( HelnaRMS /= 0 ) THEN
           HelnaN =HelnaN/HelnaRMSN
           HelnaS =HelnaS/HelnaRMSS
        ELSE
           HelnaN =0.D0
           HelnaS =0.D0
        END IF
        IF ( HelRMS /= 0 ) THEN
           HelN =HelN/HelRMSN
           HelS =HelS/HelRMSS
           HelEA=HelEA/HelRMS
        ELSE
           HelN =0.D0
           HelS =0.D0
           HelEA=0.D0
        END IF
     END IF
  ELSE
     HelN     =0.D0
     HelS     =0.D0
     HelEA    =0.D0
     HelRMSN  =0.D0
     HelRMSS  =0.D0
     HelnaN   =0.D0
     HelnaS   =0.D0
     HelnaRMSN=0.D0
     HelnaRMSS=0.D0
  ENDIF

  IF ( l_par ) THEN
     CALL getEgeos(timeScaled,nLogs,w,dw,ddw,z,dz, &
          &        Egeos,EkNTC,EkSTC,Ekin,         &
          &        dpFlow,dzFlow,CVzOTC,CVorOTC,CHelOTC)
     Geos=Egeos/Ekin ! Output, relative geostrophic kinetic Energy
  ELSE
     Egeos  =0.D0
     EkNTC  =0.D0
     EkSTC  =0.D0
     Ekin   =-1.D0 ! Only used for ratio, musst thus be non-zero
     dpFlow =0.D0
     dzFlow =0.D0
     Geos   =0.D0
     CVzOTC =0.D0
     CVorOTC=0.D0
     CHelOTC=0.D0
  ENDIF

  IF (rank.EQ.0) THEN
     !-- Evaluate nusselt numbers (boundary heat flux density):
     osq4pi =1.D0/DSQRT(4.D0*pi)
     IF (topcond/=0.D0 .AND. l_heat) THEN
        botnuss=-osq4pi/botcond*REAL(ds(1,n_r_icb))/lScale
        topnuss=-osq4pi/topcond*REAL(ds(1,n_r_cmb))/lScale
     ELSE
        botnuss=0.D0
        topnuss=0.D0
     END IF
     botflux=-rho0(n_r_max)*temp0(n_r_max)*REAL(ds(1,n_r_max))/lScale* &
          r_icb**2*DSQRT(4.D0*pi)*kappa(n_r_max)
     topflux=-rho0(1)*temp0(1)*REAL(ds(1,1))/lScale*r_cmb**2* &
          DSQRT(4.D0*pi)*kappa(1)

     IF ( l_save_out ) THEN
        OPEN(n_misc_file,file=misc_file,status='unknown', &
             POSITION='APPEND')
     ENDIF
     WRITE(n_misc_file,'(1P,D20.12,21D16.8)')     &
          & timeScaled,botnuss,topnuss, &
          & REAL(s(1,n_r_icb)),REAL(s(1,n_r_cmb)), &
          & HelN,HelS,HelRMSN,HelRMSS, &
          & Egeos/Ekin,EkNTC/Ekin,EkSTC/Ekin,Ekin, & !10-13
          & CVzOTC,CVorOTC,CHelOTC,                & !14-16   
          & HelnaN,HelnaS,HelnaRMSN,HelnaRMSS, &
          & botflux,topflux
     IF ( l_save_out ) CLOSE(n_misc_file)
     !--- NOTE: Ekin can be compared with energy in e_kin.TAG to
     !    get an idea of the precision of cylindrical integration in getEgeos.
  END IF

  IF ( l_prms ) THEN
     DO n_r=1,n_r_max
        pplot(n_r)=0.D0
        DO lm=llm,ulm
           m=lo_map%lm2m(lm)
           pplot(n_r)=pplot(n_r)+cc2real(p(lm,n_r),m)
        end do
     END DO
     CALL MPI_Reduce(pplot,pplot_global,n_r_max,MPI_DOUBLE_PRECISION,&
          & MPI_SUM,0,MPI_COMM_WORLD,ierr)
     ! Send the p(4,4) value to rank 0
     lm44=lo_map%lm2(4,4)
     lm44_is_local=(llm.LE.lm44).AND.(lm44.LE.ulm)
     mytag=120
     IF (lm44_is_local.and.(rank.ne.0)) THEN
        ! copy one row of p into a vector to send
        ! it to rank 0
        p44_local=p(lm44,:)
        CALL MPI_Send(p44_local,n_r_max,MPI_DOUBLE_COMPLEX,0,mytag,MPI_COMM_WORLD,ierr)
     END IF
     IF (rank.EQ.0) THEN
        IF (.NOT.lm44_is_local) THEN
           CALL MPI_Recv(p44_local,n_r_max,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE,mytag,&
                & MPI_COMM_WORLD,status,ierr)
        ELSE
           p44_local=p(lm44,:)
        END IF
        filename2='p.'//TAG
        OPEN(94,FILE=filename2,STATUS='UNKNOWN')
        DO n_r=1,n_r_max
           pplot_global(n_r)=SQRT(pplot_global(n_r)/lm_max)
           WRITE(94,*) r(n_r),pplot_global(n_r),REAL(p44_local(n_r))
        END DO
        CLOSE(94)
     END IF
  END IF

  RETURN
end SUBROUTINE outMisc

!---------------------------------------------------------------------------
