!$Id$
!*************************************************************************
    SUBROUTINE outMisc(timeScaled,HelLMr,Hel2LMr,HelnaLMr,Helna2LMr, &
                       nLogs,w,dw,ddw,z,dz,s,ds,p,Geos,dpFlow,dzFlow)
!*************************************************************************

    USE truncation
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE blocking
    USE horizontal_data
    USE logic
    USE output_data
    use Egeos_mod
    USE const
    USE usefull, ONLY: cc2real
    USE integration, ONLY: rInt,rInt_R

    IMPLICIT NONE

!-- Input of variables:
    REAL(kind=8) :: timeScaled
    REAL(kind=8) :: HelLMr(l_max+1,n_r_max)
    REAL(kind=8) :: Hel2LMr(l_max+1,n_r_max)
    REAL(kind=8) :: HelnaLMr(l_max+1,n_r_max)
    REAL(kind=8) :: Helna2LMr(l_max+1,n_r_max)
    INTEGER :: nLogs

!-- Input of scalar fields:
    COMPLEX(kind=8) :: s(lm_max,n_r_max)
    COMPLEX(kind=8) :: ds(lm_max,n_r_max)
!---- Fields transfered to getEgeos:
    COMPLEX(kind=8) :: w(lm_max,n_r_max)
    COMPLEX(kind=8) :: dw(lm_max,n_r_max)
    COMPLEX(kind=8) :: ddw(lm_max,n_r_max)
    COMPLEX(kind=8) :: z(lm_max,n_r_max)
    COMPLEX(kind=8) :: dz(lm_max,n_r_max)
    COMPLEX(kind=8) :: p(lm_max,n_r_max)
    REAL(kind=8) :: pplot(n_r_max)


!-- Output: (and stuff written in misc.TAG files)
    REAL(kind=8) :: Geos
    REAL(kind=8) :: dpFlow,dzFlow

!-- Local stuff:
    INTEGER :: nTheta,nThetaStart,nThetaBlock,nThetaNHS,n
    REAL(kind=8) :: Hel(nfs),Hel2(nfs),Helna(nfs),Helna2(nfs),r2
    REAL(kind=8) :: HelNr(n_r_max),HelSr(n_r_max),HelN,HelS
    REAL(kind=8) :: HelnaNr(n_r_max),HelnaSr(n_r_max),HelnaN,HelnaS
    REAL(kind=8) :: Helna2Nr(n_r_max),Helna2Sr(n_r_max),HelnaRMSN,HelnaRMSS
    REAL(kind=8) :: Hel2Nr(n_r_max),Hel2Sr(n_r_max),HelEAr(n_r_max)
    REAL(kind=8) :: HelRMSN,HelRMSS,HelEA,HelRMS,HelnaRMS
    REAL(kind=8) :: Egeos,EkNTC,EkSTC,Ekin
    REAL(kind=8) :: CVzOTC,CVorOTC,CHelOTC
    REAL(kind=8) :: topnuss,botnuss
    REAL(kind=8) :: topflux,botflux

    INTEGER :: n_r,m
    REAL(kind=8) :: osq4pi

    CHARACTER(len=76) :: filename2

!-- end of declaration
!---------------------------------------------------------------------


    IF ( l_hel )  THEN

    !------ Integration of Helicity, on input the Helicity is
    !       already axisymmetric !
        DO n_r=1,n_r_max
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
                             nThetaStart,sizeThetaB)
                CALL lmAS2pt(Hel2LMr(1,n_r),Hel2,     &
                             nThetaStart,sizeThetaB)
                CALL lmAS2pt(HelnaLMr(1,n_r),Helna,   &
                             nThetaStart,sizeThetaB)
                CALL lmAS2pt(Helna2LMr(1,n_r),Helna2, &
                             nThetaStart,sizeThetaB)
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

    !------ Integration over r without the boundries and normalization:
        HelN  =rInt(HelNr,n_r_max,dr_fac,i_costf_init,d_costf_init)
        HelS  =rInt(HelSr,n_r_max,dr_fac,i_costf_init,d_costf_init)
        HelnaN=rInt(HelnaNr,n_r_max,dr_fac,i_costf_init,d_costf_init)
        HelnaS=rInt(HelnaSr,n_r_max,dr_fac,i_costf_init,d_costf_init)
        HelEA =rInt(HelEAr,n_r_max,dr_fac,i_costf_init,d_costf_init)
        HelRMSN=rInt(Hel2Nr,n_r_max,dr_fac,i_costf_init,d_costf_init)
        HelRMSS=rInt(Hel2Sr,n_r_max,dr_fac,i_costf_init,d_costf_init)
        HelnaRMSN=rInt(Helna2Nr,n_r_max,dr_fac,i_costf_init,d_costf_init)
        HelnaRMSS=rInt(Helna2Sr,n_r_max,dr_fac,i_costf_init,d_costf_init)
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
                              Egeos,EkNTC,EkSTC,Ekin, &
                dpFlow,dzFlow,CVzOTC,CVorOTC,CHelOTC)
        Geos=Egeos/Ekin ! Output, realtive geostrophic kinetic Energy
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

!-- Evaluate nusselt numbers (boundary heat flux density):
    osq4pi =1.D0/DSQRT(4.D0*pi)
    IF (topcond.NE.0.D0) THEN
        botnuss=-osq4pi/botcond*REAL(ds(1,n_r_icb))/lScale
        topnuss=-osq4pi/topcond*REAL(ds(1,n_r_cmb))/lScale
    ELSE
        botnuss=0.D0
        topnuss=0.D0
    END IF
    botflux=-rho0(n_r_max)**PolFac*REAL(ds(1,n_r_max))/lScale* &
             r_icb**2*DSQRT(4.D0*pi)*kappa(n_r_max)
    topflux=-rho0(1)**PolFac*REAL(ds(1,1))/lScale*r_cmb**2* &
             DSQRT(4.D0*pi)*kappa(1)

    IF ( l_save_out ) THEN
        OPEN(n_misc_file,file=misc_file,status='unknown', &
             POSITION='APPEND')
    ENDIF
    WRITE(n_misc_file,'(1P,D20.12,21D16.8)')     &
                     timeScaled,botnuss,topnuss, &
        REAL(s(1,n_r_icb)),REAL(s(1,n_r_cmb)), &
                      HelN,HelS,HelRMSN,HelRMSS, &
          Egeos/Ekin,EkNTC/Ekin,EkSTC/Ekin,Ekin, &
                         CVzOTC,CVorOTC,CHelOTC, &
              HelnaN,HelnaS,HelnaRMSN,HelnaRMSS, &
                                botflux,topflux
    IF ( l_save_out ) CLOSE(n_misc_file)
!--- NOTE: Ekin can be compared with energy in e_kin.TAG to
!    get an idea of the precission of cylindrical integration in getEgeos.

    if ( l_prms ) then
        filename2='p.'//TAG
        OPEN(94,FILE=filename2,STATUS='UNKNOWN')
        DO n_r=1,n_r_max
            pplot(n_r)=0.D0
            do n=1,lm_max
                m=lm2m(n)
                pplot(n_r)=pplot(n_r)+cc2real(p(n,n_r),m)
            end do
            pplot(n_r)=dsqrt(pplot(n_r)/lm_max)
            WRITE(94,*) r(n_r),pplot(n_r),REAL(p(lm2(4,4),n_r))
        END DO
        CLOSE(94)
    end if

    RETURN
    end SUBROUTINE outMisc

!---------------------------------------------------------------------------
