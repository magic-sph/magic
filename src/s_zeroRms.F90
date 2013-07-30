!$Id$
!********************************************************************
    SUBROUTINE zeroRms
!********************************************************************

!--------------------------------------------------------------------

!  Zeros integrals that are set in s_get_td, s_update_z,
!  s_update_wp, s_update_b, s_dtVrms and s_dtBrms

!--------------------------------------------------------------------

    USE truncation
    USE blocking
    USE RMS

    IMPLICIT NONE

    INTEGER :: n,nR,lm

!-------------------------------------------------------------------


    DO n=1,nThreadsMax
        DO nR=1,n_r_max
            dtVPol2hInt(nR,n)  =0.D0
            dtVTor2hInt(nR,n)  =0.D0
            dtVPolAs2hInt(nR,n)=0.D0
            dtVTorAs2hInt(nR,n)=0.D0
            DifPol2hInt(nR,n)  =0.D0
            DifTor2hInt(nR,n)  =0.D0
            DifPolAs2hInt(nR,n)=0.D0
            DifTorAs2hInt(nR,n)=0.D0
        END DO
        DO nR=1,n_r_maxMag
            dtBPol2hInt(nR,n)  =0.D0
            dtBTor2hInt(nR,n)  =0.D0
            dtBPolAs2hInt(nR,n)=0.D0
            dtBTorAs2hInt(nR,n)=0.D0
        END DO
    END DO

    DO nR=1,n_r_max
        AdvPol2hInt(nR)  =0.D0
        AdvTor2hInt(nR)  =0.D0
        AdvPolAs2hInt(nR)=0.D0
        AdvTorAs2hInt(nR)=0.D0
        CorPol2hInt(nR)  =0.D0
        CorTor2hInt(nR)  =0.D0
        CorPolAs2hInt(nR)=0.D0
        CorTorAs2hInt(nR)=0.D0
        LFPol2hInt(nR)   =0.D0
        LFTor2hInt(nR)   =0.D0
        LFPolAs2hInt(nR) =0.D0
        LFTorAs2hInt(nR) =0.D0
        Buo2hInt(nR)     =0.D0
        BuoAs2hInt(nR)   =0.D0
        Pre2hInt(nR)     =0.D0
        PreAs2hInt(nR)   =0.D0
        Geo2hInt(nR)     =0.D0
        GeoAs2hInt(nR)   =0.D0
        Mag2hInt(nR)     =0.D0
        MagAs2hInt(nR)   =0.D0
        Arc2hInt(nR)     =0.D0
        ArcAs2hInt(nR)   =0.D0
    END DO
    DO nR=1,n_r_max
        DO lm=1,lm_max
            dtVPolLMr(lm,nR)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            AdvPolLMr(lm,nR)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            CorPolLMr(lm,nR)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            DifPolLMr(lm,nR)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            BuoLMr(lm,nR)   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            PreLMr(lm,nR)   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            GeoLMr(lm,nR)   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            ArcLMr(lm,nR)   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        END DO
    END DO
    DO nR=1,n_r_maxMag
        DO lm=1,lm_maxMag
            dtBPolLMr(lm,nR)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            LFPolLMr(lm,nR) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        END DO
    END DO


    RETURN
    end SUBROUTINE zeroRms

!----------------------------------------------------------------
