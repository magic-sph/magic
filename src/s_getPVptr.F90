!$Id$
!****************************************************************************
    SUBROUTINE getPVptr(w,dw,ddw,z,dz,rMin,rMax,rS, &
                    nZmax,nZmaxA,PlmS,dPlmS,OsinTS, &
                           VrS,VpS,VtS,VorS,dpVorS)
!****************************************************************************
!  This subroutine calculated the three flow conponents VrS,VtS,VpS at
!  (r,theta,all phis) and (r,pi-theta, all phis). Here r=rS, PlmS=Plm(theta),
!  dPlmS=sin(theta)*dTheta Plm(theta), and OsinTS=1/sin(theta).
!  The flow is calculated for all n_phi_max azimuthal points used in the code,
!  and for corresponding latitudes north and south of the equator.
!  For lDeriv=.TRUE. the subroutine also calculates dpEk and dzEk which
!  are phi averages of (d Vr/d phi)**2 + (d Vtheta/ d phi)**2 + (d Vphi/ d phi)**2
!  and (d Vr/d z)**2 + (d Vtheta/ d z)**2 + (d Vphi/ d z)**2, respectively.
!  These two quantities are used ot calculate z and phi scale of the flow in
!  s_getEgeos.f
!  NOTE: on input w=l*(l+1)*w
!---------------------------------------------------------------------------------

        USE truncation
        USE radial_functions
        USE blocking
        USE horizontal_data
        USE const
#if (FFTLIB==JW)
        USE fft_JW
#elif (FFTLIB==MKL)
        USE fft_MKL
#endif
    IMPLICIT NONE

!--- Input:
    COMPLEX(kind=8) :: w(lm_max,n_r_max)
    COMPLEX(kind=8) :: dw(lm_max,n_r_max)
    COMPLEX(kind=8) :: ddw(lm_max,n_r_max)
    COMPLEX(kind=8) :: z(lm_max,n_r_max)
    COMPLEX(kind=8) :: dz(lm_max,n_r_max)
    REAL(kind=8) ::     rMin,rMax  ! radial bounds
    INTEGER ::    nZmax,nZmaxA ! number of (r,theta) point to be transformed
    REAL(kind=8) ::     rS(nZmaxA)
    REAL(kind=8) ::     PlmS(lm_max,nZmaxA/2+1)
    REAL(kind=8) ::     dPlmS(lm_max,nZmaxA/2+1)
    REAL(kind=8) ::     OsinTS(nZmaxA/2+1)

!--- Output: function on azimuthal grid points defined by FT!
    COMPLEX(kind=8) :: VrS(ncp,nZmaxA)
    COMPLEX(kind=8) :: VtS(ncp,nZmaxA)
    COMPLEX(kind=8) :: VpS(ncp,nZmaxA)
    COMPLEX(kind=8) :: VorS(ncp,nZmaxA)
    COMPLEX(kind=8) :: dpVorS(ncp,nZmaxA)

!--- Local:
    INTEGER,PARAMETER :: nChebMaxL=150
    REAL(kind=8) ::    chebS(nChebMaxL)
    INTEGER ::   nS,nN,mc,lm,l,m,nCheb
    REAL(kind=8) ::    x,phiNorm,mapFac,OS,cosT,sinT,Or_e1,Or_e2
    COMPLEX(kind=8) :: Vr,Vt,Vt1,Vt2,Vp1,Vp2,Vor,Vot1,Vot2
    COMPLEX(kind=8) :: VotS(ncp,nZmaxA)
    COMPLEX(kind=8) :: wSr,dwSr,ddwSr,zSr,dzSr

    INTEGER,PARAMETER :: nZmaxAL=1000
    INTEGER,PARAMETER :: nZmaxL=MAX0(nZmaxAL,2)  ! Data along z !
    !REAL(kind=8) :: work(nrp,nZmaxL)
    COMPLEX(kind=8) :: dp

!----------------------------------------------------------------------------

    IF ( n_r_max > nChebMaxL ) THEN
        WRITE(*,*) '! nChebMaxL too small in getDVptr!'
        WRITE(*,*) '! Should be at least:',n_r_max
        STOP
    END IF
    IF ( nZmax > nZmaxL) THEN
        WRITE(*,*) '! nZmaxL too small in getDVptr!'
        WRITE(*,*) '! Should be at least:',2*nZmax
        STOP
    END IF

    mapFac=2.D0/(rMax-rMin)
    phiNorm=2.D0*pi/n_phi_max

    DO nS=1,nZmax
        DO mc=1,ncp
            VrS(mc,nS) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            VtS(mc,nS) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            VpS(mc,nS) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            VorS(mc,nS)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            VotS(mc,nS)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        END DO
    END DO

    DO nN=1,nZmax/2    ! Loop over all (r,theta) points in NHS
        nS=nZmax-nN+1   ! Southern counterpart !

    !------ Calculate Chebs:
    !------ Map r to cheb intervall [-1,1]:
    !       and calculate the cheb polynomia:
    !       Note: the factor cheb_norm is needed
    !       for renormalisation. Its not needed if one used
    !       costf1 for the back transform.
        x=2.D0*(rS(nN)-0.5D0*(rMin+rMax))/(rMax-rMin)
        chebS(1) =1.D0*cheb_norm ! Extra cheb_norm cheap here
        chebS(2) =x*cheb_norm
        DO nCheb=3,n_r_max
            chebS(nCheb)=2.D0*x*chebS(nCheb-1)-chebS(nCheb-2)
        END DO
        chebS(1)      =0.5D0*chebS(1)
        chebS(n_r_max)=0.5D0*chebS(n_r_max)
        Or_e2=1.D0/rS(nN)**2

        DO lm=1,lm_max     ! Sum over lms
            l =lm2l(lm)
            m =lm2m(lm)
            mc=lm2mc(lm)
            wSr  =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            dwSr =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            ddwSr=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            zSr  =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            dzSr =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            DO nCheb=1,n_r_max
                wSr  =  wSr+  w(lm,nCheb)*chebS(nCheb)
                dwSr = dwSr+ dw(lm,nCheb)*chebS(nCheb)
                ddwSr=ddwSr+ddw(lm,nCheb)*chebS(nCheb)
                zSr  =  zSr+  z(lm,nCheb)*chebS(nCheb)
                dzSr = dzSr+ dz(lm,nCheb)*chebS(nCheb)
            END DO
            Vr  =  wSr* PlmS(lm,nN)
            Vt1 = dwSr*dPlmS(lm,nN)
            Vt2 =  zSr* PlmS(lm,nN)*dPhi(lm)
            Vp1 = dwSr* PlmS(lm,nN)*dPhi(lm)
            Vp2 = -zSr*dPlmS(lm,nN)
            Vor =  zSr* PlmS(lm,nN)*dLh(lm)
            Vot1= dzSr*dPlmS(lm,nN)
            Vot2= (wSr*Or_e2-ddwSr) * &
                   PlmS(lm,nN)*dPhi(lm)
            VrS(mc,nN) =VrS(mc,nN) +Vr
            VtS(mc,nN) =VtS(mc,nN) +Vt1+Vt2
            VpS(mc,nN) =VpS(mc,nN) +Vp1+Vp2
            VorS(mc,nN)=VorS(mc,nN)+Vor
            VotS(mc,nN)=VotS(mc,nN)+Vot1+Vot2
            IF ( MOD(l+m,2) == 0 ) THEN
                VrS(mc,nS) =VrS(mc,nS) +Vr
                VtS(mc,nS) =VtS(mc,nS) -Vt1+Vt2
                VpS(mc,nS) =VpS(mc,nS) +Vp1-Vp2
                VorS(mc,nS)=VorS(mc,nS)+Vor
                VotS(mc,nS)=VotS(mc,nS)-Vot1+Vot2
            ELSE
                VrS(mc,nS) =VrS(mc,nS) -Vr
                VtS(mc,nS) =VtS(mc,nS) +Vt1-Vt2
                VpS(mc,nS) =VpS(mc,nS) -Vp1+Vp2
                VorS(mc,nS)=VorS(mc,nS)-Vor
                VotS(mc,nS)=VotS(mc,nS)+Vot1-Vot2
            END IF
        END DO

    END DO

    IF ( MOD(nZmax,2) == 1 ) THEN ! Remaining equatorial point
        nS=(nZmax+1)/2

        x=2.D0*(rS(nS)-0.5D0*(rMin+rMax))/(rMax-rMin)
        chebS(1)=1.D0*cheb_norm ! Extra cheb_norm cheap here
        chebS(2)=x*cheb_norm
        DO nCheb=3,n_r_max
            chebS(nCheb)=2.D0*x*chebS(nCheb-1)-chebS(nCheb-2)
        END DO
        chebS(1)      =0.5D0*chebS(1)
        chebS(n_r_max)=0.5D0*chebS(n_r_max)
        Or_e2=1.D0/rS(nS)**2

        DO lm=1,lm_max     ! Sum over lms
            l =lm2l(lm)
            m =lm2m(lm)
            mc=lm2mc(lm)
            wSr  =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            dwSr =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            ddwSr=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            zSr  =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            dzSr =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            DO nCheb=1,n_r_max
                wSr  =  wSr+  w(lm,nCheb)*chebS(nCheb)
                dwSr = dwSr+ dw(lm,nCheb)*chebS(nCheb)
                ddwSr=ddwSr+ddw(lm,nCheb)*chebS(nCheb)
                zSr  =  zSr+  z(lm,nCheb)*chebS(nCheb)
                dzSr = dzSr+ dz(lm,nCheb)*chebS(nCheb)
            END DO
            Vr  =  wSr* PlmS(lm,nS)
            Vt1 = dwSr*dPlmS(lm,nS)
            Vt2 =  zSr* PlmS(lm,nS)*dPhi(lm)
            Vp1 = dwSr* PlmS(lm,nS)*dPhi(lm)
            Vp2 = -zSr*dPlmS(lm,nS)
            Vor =  zSr* PlmS(lm,nS)*dLh(lm)
            Vot1= dzSr*dPlmS(lm,nS)
            Vot2= (wSr*Or_e2-ddwSr) * &
                   PlmS(lm,nS)*dPhi(lm)
            VrS(mc,nS) =VrS(mc,nS) +Vr
            VtS(mc,nS) =VtS(mc,nS) +Vt1+Vt2
            VpS(mc,nS) =VpS(mc,nS) +Vp1+Vp2
            VorS(mc,nS)=VorS(mc,nS)+Vor
            VotS(mc,nS)=VotS(mc,nS)+Vot1+Vot2
        END DO

    END IF ! Equatorial point ?

!--- Extra factors, contructing z-vorticity:
    DO nS=1,(nZmax+1)/2 ! North HS
        OS   =OsinTS(nS)
        sinT =1.D0/OS
        cosT =SQRT(1.D0-sinT**2)
        Or_e1=1.D0/rS(nS)
        Or_e2=Or_e1*Or_e1
    !           VrS(1,nS) =0.D0
    !           VpS(1,nS) =0.D0
    !           VtS(1,nS) =0.D0
    !           VorS(1,nS)=0.D0
        DO mc=1,n_m_max
            Vr=Or_e2*VrS(mc,nS)
            Vt=Or_e1*OS*VtS(mc,nS)
            VrS(mc,nS) =sinT*Vr+cosT*Vt ! this is now Vs
            VpS(mc,nS) =Or_e1*OS*VpS(mc,nS)
            VtS(mc,nS) =cosT*Vr-sinT*Vt ! this is now Vz
            VorS(mc,nS)=cosT*Or_e2*VorS(mc,nS) - &
                        Or_e1*VotS(mc,nS)
        END DO
        DO mc=n_m_max+1,ncp
            VrS(mc,nS) =0.D0
            VpS(mc,nS) =0.D0
            VtS(mc,nS) =0.D0
            VorS(mc,nS)=0.D0
        END DO
    END DO

    DO nS=(nZmax+1)/2+1,nZmax ! South HS
        OS   =OsinTS(nZmax-nS+1)
        sinT =1.D0/OS
        cosT =-SQRT(1.D0-sinT**2)
        Or_e1=1.D0/rS(nZmax-nS+1)
        Or_e2=Or_e1*Or_e1
    !           VrS(1,nS) =0.D0
    !           VpS(1,nS) =0.D0
    !           VtS(1,nS) =0.D0
    !           VorS(1,nS)=0.D0
        DO mc=1,n_m_max
            Vr=Or_e2*VrS(mc,nS)
            Vt=Or_e1*OS*VtS(mc,nS)
            VrS(mc,nS) =sinT*Vr+cosT*Vt ! this is now Vs
            VpS(mc,nS) =Or_e1*OS*VpS(mc,nS)
            VtS(mc,nS) =cosT*Vr-sinT*Vt ! this is now Vz
            VorS(mc,nS)=cosT*Or_e2*VorS(mc,nS) - &
                        Or_e1*VotS(mc,nS)
        END DO
        DO mc=n_m_max+1,ncp
            VrS(mc,nS) =0.D0
            VpS(mc,nS) =0.D0
            VtS(mc,nS) =0.D0
            VorS(mc,nS)=0.D0
        END DO
    END DO

    DO nS=1,nZmax
    !          dpVorS(1,nS)=0.D0
        DO mc=1,n_m_max
            dp=CMPLX(0.D0,1.D0,KIND=KIND(0d0))*DBLE((mc-1)*minc)  ! - i m
            dpVorS(mc,nS)=dp*VorS(mc,nS) ! first test with Vs!!!
        END DO
        DO mc=n_m_max+1,ncp
            dpVorS(mc,nS)=0.D0
        END DO
    END DO


!----- Transform m 2 phi for flow field:
        call fft_to_real(VrS,nrp,nZmax)
        call fft_to_real(VtS,nrp,nZmax)
        call fft_to_real(VpS,nrp,nZmax)
        call fft_to_real(VorS,nrp,nZmax)
        call fft_to_real(dpVorS,nrp,nZmax)
!        cALL fftJW(VrS,nrp,n_phi_max,1,nZmax,work,nrp,nZmaxL,
!     &                                 i_fft_init,d_fft_init)
!        CALL fftJW(VtS,nrp,n_phi_max,1,nZmax,work,nrp,nZmaxL,
!     &                                 i_fft_init,d_fft_init)
!        CALL fftJW(VpS,nrp,n_phi_max,1,nZmax,work,nrp,nZmaxL,
!     &                                 i_fft_init,d_fft_init)
!        CALL fftJW(VorS,nrp,n_phi_max,1,nZmax,work,nrp,nZmaxL,
!     &                                  i_fft_init,d_fft_init)
!        CALL fftJW(dpVorS,nrp,n_phi_max,1,nZmax,work,nrp,nZmaxL,
!     &                                    i_fft_init,d_fft_init)


    RETURN
    end SUBROUTINE getPVptr

!----------------------------------------------------------------------
