!$Id$
!********************************************************************
    SUBROUTINE dtVrms(time,nRMS_sets)
!********************************************************************

!  +-------------+----------------+------------------------------------+
!  |  For testing RMS balance and determining the necessart values of  |
!  |  rDea and rCut one can use the l_RMStest=true options.            |
!  |  In this case the poloidal and toroidal RMS dtV which also        |
!  |  include the diffusion effects should be identical (as close as   |
!  |  desired) to the RMS sum of forces stored in the Geo value and    |
!  |  the Mag value, respectively.                                     |
!  |  An additional tests is the Arc value which should be identical   |
!  |  to the poloidal kinetic energy.                                  |
!  |                                                                   |
!  |  N.B. The second test with the Arc value cannot work so easily in |
!  |       the anelastic version as the density enters now in the      |
!  |       force balance. The first test should work, though.          |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+

    USE truncation
    USE radial_functions
    USE physical_parameters
    USE blocking
    USE logic
    USE RMS
    USE output_data
    USE const
    USE parallel_mod
    USE integration, ONLY: rInt_R

    IMPLICIT NONE

!-- Input variable:
    REAL(kind=8) :: time
    INTEGER :: nRMS_sets

!-- Output:
    REAL(kind=8) :: dtVPolRms,dtVPolAsRms
    REAL(kind=8) :: dtVTorRms,dtVTorAsRms
    REAL(kind=8) :: CorPolRms,CorPolAsRms
    REAL(kind=8) :: CorTorRms,CorTorAsRms
    REAL(kind=8) :: AdvPolRms,AdvPolAsRms
    REAL(kind=8) :: AdvTorRms,AdvTorAsRms
    REAL(kind=8) :: LFPolRms, LFPolAsRms
    REAL(kind=8) :: LFTorRms, LFTorAsRms
    REAL(kind=8) :: DifPolRms,DifPolAsRms
    REAL(kind=8) :: DifTorRms,DifTorAsRms
    REAL(kind=8) :: BuoRms,   BuoAsRms
    REAL(kind=8) :: PreRms,   PreAsRms
    REAL(kind=8) :: GeoRms,   GeoAsRms
    REAL(kind=8) :: MagRms,   MagAsRms
    REAL(kind=8) :: ArcRms,   ArcAsRms

!-- Local:
    INTEGER :: nR,n,nRC
    REAL(kind=8) :: volC
    REAL(kind=8) :: bal1,bal2,bal3

    COMPLEX(kind=8) :: workA(lm_max,n_r_max),workB(lm_max,n_r_max)

!-- end of declaration
!---------------------------------------------------------------------

    IF ( nRMS_sets == 0 ) THEN
        nRMS_sets=nRMS_sets+1
    !--- Initialize new cut-back grid:
        CALL init_rNB(r,n_r_max,n_cheb_max,rCut,rDea, &
                        rC,n_r_maxC,n_cheb_maxC,nCut, &
                    dr_facC,i_costf_initC,nDi_costf1, &
                            d_costf_initC,nDd_costf1)
    END IF
    nRC=nCut+1
    volC=4.D0/3.D0*pi*(r(1+nCut)**3-r(n_r_max-nCut)**3)

    IF ( l_conv ) THEN

    !------ Coriolis force
        IF ( l_corr ) THEN
            CALL get_drNS_R(     CorPolLMr(1,nRC),workA(1,nRC), &
                lm_max_real,1,lm_max_real,n_r_maxC,n_cheb_maxC, &
                     dr_facC,workB,i_costf_initC,d_costf_initC)
            DO nR=1,n_r_maxC
                CALL hInt2dPol(workA(1,nR+nCut),2,lm_max, &
                               CorPol2hInt(nR+nCut),CorPolAs2hInt(nR+nCut))
            END DO
            CorPolRms=rInt_R(CorPol2hInt(nRC),n_r_maxC, &
                                   n_cheb_maxC,dr_facC, &
                           i_costf_initC,d_costf_initC)
            CorPolAsRms=rInt_R(CorPolAs2hInt(nRC),n_r_maxC, &
                                       n_cheb_maxC,dr_facC, &
                               i_costf_initC,d_costf_initC)
            CorTorRms  =rInt_R(CorTor2hInt(nRC),n_r_maxC, &
                                     n_cheb_maxC,dr_facC, &
                             i_costf_initC,d_costf_initC)
            CorTorAsRms=rInt_R(CorTorAs2hInt(nRC),n_r_maxC, &
                                       n_cheb_maxC,dr_facC, &
                               i_costf_initC,d_costf_initC)
            CorPolRms  =DSQRT(CorPolRms  /volC)
            CorPolAsRms=DSQRT(CorPolAsRms/volC)
            CorTorRms  =DSQRT(CorTorRms  /volC)
            CorTorAsRms=DSQRT(CorTorAsRms/volC)
        END IF

    !------ Advection:
        IF ( l_conv_nl ) THEN
            CALL get_drNS_R(    AdvPolLMr(1,nRC),workA(1,nRC), &
               lm_max_real,1,lm_max_real,n_r_maxC,n_cheb_maxC, &
                    dr_facC,workB,i_costf_initC,d_costf_initC)
            DO nR=1,n_r_maxC
                CALL hInt2dPol(workA(1,nR+nCut),2,lm_max, &
                               AdvPol2hInt(nR+nCut),AdvPolAs2hInt(nR+nCut))
            END DO
            AdvPolRms  =rInt_R(AdvPol2hInt(nRC),n_r_maxC, &
                                     n_cheb_maxC,dr_facC, &
                             i_costf_initC,d_costf_initC)
            AdvPolAsRms=rInt_R(AdvPolAs2hInt(nRC),n_r_maxC, &
                                       n_cheb_maxC,dr_facC, &
                               i_costf_initC,d_costf_initC)
            AdvTorRms  =rInt_R(AdvTor2hInt(nRC),n_r_maxC, &
                                     n_cheb_maxC,dr_facC, &
                             i_costf_initC,d_costf_initC)
            AdvTorAsRms=rInt_R(AdvTorAs2hInt(nRC),n_r_maxC, &
                                       n_cheb_maxC,dr_facC, &
                               i_costf_initC,d_costf_initC)
            AdvPolRms  =DSQRT(AdvPolRms  /volC)
            AdvPolAsRms=DSQRT(AdvPolAsRms/volC)
            AdvTorRms  =DSQRT(AdvTorRms  /volC)
            AdvTorAsRms=DSQRT(AdvTorAsRms/volC)
        END IF

    !------ Lorentz force:
        IF ( l_mag_LF ) THEN
            CALL get_drNS_R(     LFPolLMr(1,nRC),workA(1,nRC), &
               lm_max_real,1,lm_max_real,n_r_maxC,n_cheb_maxC, &
                    dr_facC,workB,i_costf_initC,d_costf_initC)
            DO nR=1,n_r_maxC
                CALL hInt2dPol(     workA(1,nR+nCut),2,lm_max, &
                      LFPol2hInt(nR+nCut),LFPolAs2hInt(nR+nCut))
            END DO
            LFPolRms  =rInt_R(LFPol2hInt(nRC),n_r_maxC, &
                                   n_cheb_maxC,dr_facC, &
                           i_costf_initC,d_costf_initC)
            LFPolAsRms=rInt_R(LFPolAs2hInt(nRC),n_r_maxC, &
                                     n_cheb_maxC,dr_facC, &
                             i_costf_initC,d_costf_initC)
            LFTorRms  =rInt_R(LFTor2hInt(nRC),n_r_maxC, &
                                   n_cheb_maxC,dr_facC, &
                           i_costf_initC,d_costf_initC)
            LFTorAsRms=rInt_R(LFTorAs2hInt(nRC),n_r_maxC, &
                                     n_cheb_maxC,dr_facC, &
                             i_costf_initC,d_costf_initC)
            LFPolRms  =DSQRT(LFPolRms  /volC)
            LFPolAsRms=DSQRT(LFPolAsRms/volC)
            LFTorRms  =DSQRT(LFTorRms  /volC)
            LFTorAsRms=DSQRT(LFTorAsRms/volC)
        ELSE
            LFPolRms  =0.D0
            LFPolAsRms=0.D0
            LFTorRms  =0.D0
            LFTorAsRms=0.D0
        END IF

    !------ Buoyancy:
        IF ( l_heat ) THEN
            CALL get_drNS_R(      BuoLMr(1,nRC),workA(1,nRC), &
              lm_max_real,1,lm_max_real,n_r_maxC,n_cheb_maxC, &
                   dr_facC,workB,i_costf_initC,d_costf_initC)
            DO nR=1,n_r_maxC
                CALL hInt2dPol(    workA(1,nR+nCut),2,lm_max, &
                       Buo2hInt(nR+nCut),BuoAs2hInt(nR+nCut))
                Buo2hInt(nR)  =Buo2hInt(nR)*rgrav(nR)**2
                BuoAs2hInt(nR)=BuoAs2hInt(nR)*rgrav(nR)**2
            END DO
            BuoRms  =rInt_R(Buo2hInt(nRC),n_r_maxC, &
                               n_cheb_maxC,dr_facC, &
                       i_costf_initC,d_costf_initC)
            BuoAsRms=rInt_R(BuoAs2hInt(nRC),n_r_maxC, &
                                 n_cheb_maxC,dr_facC, &
                         i_costf_initC,d_costf_initC)
            BuoRms  =DSQRT(BuoRms  /volC)
            BuoAsRms=DSQRT(BuoAsRms/volC)
        ELSE
            BuoRms=0.D0
            BuoAsRms=0.D0
        END IF

    !------ Pressure gradient:
        CALL get_drNS_R(       PreLMr(1,nRC),workA(1,nRC), &
           lm_max_real,1,lm_max_real,n_r_maxC,n_cheb_maxC, &
                dr_facC,workB,i_costf_initC,d_costf_initC)
        DO nR=1,n_r_maxC
            CALL hInt2dPol(     workA(1,nR+nCut),2,lm_max, &
                    Pre2hInt(nR+nCut),PreAs2hInt(nR+nCut))
        END DO
        PreRms  =rInt_R(Pre2hInt(nRC),n_r_maxC, &
                           n_cheb_maxC,dr_facC, &
                   i_costf_initC,d_costf_initC)
        PreAsRms=rInt_R(PreAs2hInt(nRC),n_r_maxC, &
                             n_cheb_maxC,dr_facC, &
                     i_costf_initC,d_costf_initC)
        PreRms  =DSQRT(PreRms  /volC)
        PreAsRms=DSQRT(PreAsRms/volC)

    !------ Geostrophic balance:
        CALL get_drNS_R(       GeoLMr(1,nRC),workA(1,nRC), &
           lm_max_real,1,lm_max_real,n_r_maxC,n_cheb_maxC, &
                dr_facC,workB,i_costf_initC,d_costf_initC)
        DO nR=1,n_r_maxC
            CALL hInt2dPol(     workA(1,nR+nCut),2,lm_max, &
                    Geo2hInt(nR+nCut),GeoAs2hInt(nR+nCut))
        END DO
        GeoRms  =rInt_R(Geo2hInt(nRC),n_r_maxC, &
                           n_cheb_maxC,dr_facC, &
                   i_costf_initC,d_costf_initC)
        GeoAsRms=rInt_R(GeoAs2hInt(nRC),n_r_maxC, &
                             n_cheb_maxC,dr_facC, &
                     i_costf_initC,d_costf_initC)
        GeoRms  =DSQRT(GeoRms  /volC)
        GeoAsRms=DSQRT(GeoAsRms/volC)

    !------ Magnetostrophic balance:
        IF ( .NOT. l_RMStest ) THEN
            CALL get_drNS_R(       MagLMr(1,nRC),workA(1,nRC), &
               lm_max_real,1,lm_max_real,n_r_maxC,n_cheb_maxC, &
                    dr_facC,workB,i_costf_initC,d_costf_initC)
            DO nR=1,n_r_maxC
                CALL hInt2dPol(     workA(1,nR+nCut),2,lm_max, &
                        Mag2hInt(nR+nCut),MagAs2hInt(nR+nCut))
            END DO
        END IF
        MagRms  =rInt_R(Mag2hInt(nRC),n_r_maxC, &
                           n_cheb_maxC,dr_facC, &
                   i_costf_initC,d_costf_initC)
        MagAsRms=rInt_R(MagAs2hInt(nRC),n_r_maxC, &
                             n_cheb_maxC,dr_facC, &
                     i_costf_initC,d_costf_initC)
        MagRms  =DSQRT(MagRms  /volC)
        MagAsRms=DSQRT(MagAsRms/volC)

    !------ Archemidian balance:
        CALL get_drNS_R(        ArcLMr(1,nRC),workA(1,nRC), &
            lm_max_real,1,lm_max_real,n_r_maxC,n_cheb_maxC, &
                 dr_facC,workB,i_costf_initC,d_costf_initC)
        DO nR=1,n_r_maxC
            CALL hInt2dPol(      workA(1,nR+nCut),2,lm_max, &
                     Arc2hInt(nR+nCut),ArcAs2hInt(nR+nCut))
        END DO
        ArcRms  =rInt_R(Arc2hInt(nRC),n_r_maxC, &
                           n_cheb_maxC,dr_facC, &
                   i_costf_initC,d_costf_initC)
        ArcAsRms=rInt_R(ArcAs2hInt(nRC),n_r_maxC, &
                             n_cheb_maxC,dr_facC, &
                     i_costf_initC,d_costf_initC)
        IF ( l_RMStest ) THEN
            ArcRms  =ArcRms/2.D0
            ArcAsRms=ArcAsRms/2.D0
        ELSE
            ArcRms  =DSQRT(ArcRms  /volC)
            ArcAsRms=DSQRT(ArcAsRms/volC)
        END IF

    !------ Diffusion:
        CALL get_drNS_R(     DifPolLMr(1,nRC),workA(1,nRC), &
            lm_max_real,1,lm_max_real,n_r_maxC,n_cheb_maxC, &
                 dr_facC,workB,i_costf_initC,d_costf_initC)
        DO nR=1,n_r_maxC
            CALL hInt2dPol(      workA(1,nR+nCut),2,lm_max, &
             DifPol2hInt(nR+nCut,1),DifPolAs2hInt(nR+nCut,1))
            DO n=2,nThreadsLMmax
                DifPol2hInt(nR+nCut,1)  =DifPol2hInt(nR+nCut,1) + &
                                         DifPol2hInt(nR+nCut,n)
                DifPolAs2hInt(nR+nCut,1)=DifPolAs2hInt(nR+nCut,1) + &
                                         DifPolAs2hInt(nR+nCut,n)
                DifTor2hInt(nR+nCut,1)  =DifTor2hInt(nR+nCut,1) + &
                                         DifTor2hInt(nR+nCut,n)
                DifTorAs2hInt(nR+nCut,1)=DifTorAs2hInt(nR+nCut,1) + &
                                         DifTorAs2hInt(nR+nCut,n)
            END DO
        END DO
        DifPolRms  =rInt_R(DifPol2hInt(nRC,1),n_r_maxC, &
                                   n_cheb_maxC,dr_facC, &
                           i_costf_initC,d_costf_initC)
        DifPolAsRms=rInt_R(DifPolAs2hInt(nRC,1),n_r_maxC, &
                                     n_cheb_maxC,dr_facC, &
                             i_costf_initC,d_costf_initC)
        DifTorRms  =rInt_R(DifTor2hInt(nRC,1),n_r_maxC, &
                                   n_cheb_maxC,dr_facC, &
                           i_costf_initC,d_costf_initC)
        DifTorAsRms=rInt_R(DifTorAs2hInt(nRC,1),n_r_maxC, &
                                     n_cheb_maxC,dr_facC, &
                             i_costf_initC,d_costf_initC)
        DifPolRms  =DSQRT(DifPolRms  /volC)
        DifPolAsRms=DSQRT(DifPolAsRms/volC)
        DifTorRms  =DSQRT(DifTorRms  /volC)
        DifTorAsRms=DSQRT(DifTorAsRms/volC)

    !------ Flow changes: Intertia - Advection
        CALL get_drNS_R(     dtVPolLMr(1,nRC),workA(1,nRC), &
            lm_max_real,1,lm_max_real,n_r_maxC,n_cheb_maxC, &
                 dr_facC,workB,i_costf_initC,d_costf_initC)
        DO nR=1,n_r_maxC
            CALL hInt2dPol(         workA(1,nR+nCut),2,lm_max, &
              dtVPol2hInt(nR+nCut,1),dtVPolAs2hInt(nR+nCut,1))
            DO n=2,nThreadsLMmax
                dtVPol2hInt(nR+nCut,1)  =dtVPol2hInt(nR+nCut,1) + &
                                         dtVPol2hInt(nR+nCut,n)
                dtVPolAs2hInt(nR+nCut,1)=dtVPolAs2hInt(nR+nCut,1) + &
                                         dtVPolAs2hInt(nR+nCut,n)
                dtVTor2hInt(nR+nCut,1)  =dtVTor2hInt(nR+nCut,1) + &
                                         dtVTor2hInt(nR+nCut,n)
                dtVTorAs2hInt(nR+nCut,1)=dtVTorAs2hInt(nR+nCut,1) + &
                                         dtVTorAs2hInt(nR+nCut,n)
            END DO
        END DO
        dtVPolRms  =rInt_R(dtVPol2hInt(nRC,1),n_r_maxC, &
                                   n_cheb_maxC,dr_facC, &
                           i_costf_initC,d_costf_initC)
        dtVPolAsRms=rInt_R(dtVPolAs2hInt(nRC,1),n_r_maxC, &
                                     n_cheb_maxC,dr_facC, &
                             i_costf_initC,d_costf_initC)
        dtVTorRms  =rInt_R(dtVTor2hInt(nRC,1),n_r_maxC, &
                                   n_cheb_maxC,dr_facC, &
                           i_costf_initC,d_costf_initC)
        dtVTorAsRms=rInt_R(dtVTorAs2hInt(nRC,1),n_r_maxC, &
                                     n_cheb_maxC,dr_facC, &
                             i_costf_initC,d_costf_initC)
        dtVPolRms  =DSQRT(dtVPolRms  /volC)
        dtVPolAsRms=DSQRT(dtVPolAsRms/volC)
        dtVTorRms  =DSQRT(dtVTorRms  /volC)
        dtVTorAsRms=DSQRT(dtVTorAsRms/volC)

    !----- Output:
        IF ( l_save_out) THEN
            OPEN(n_dtvrms_file,FILE=dtvrms_file,FORM='FORMATTED', &
                 STATUS='UNKNOWN',POSITION='APPEND')
        END IF
        IF ( l_RMStest ) THEN
            WRITE(n_dtvrms_file,'(1P,D20.12,15D16.8)') &
                                                 time, &
                                  dtVPolRms,dtVTorRms, &
                                  CorPolRms,CorTorRms, &
                                    LFPolRms,LFTorRms, &
                                  AdvPolRms,AdvTorRms, &
                                  DifPolRms,DifTorRms, &
                                        BuoRms,PreRms, &
                                               GeoRms, &
                                               MagRms, &
                                               ArcRms
        ELSE
            WRITE(n_dtvrms_file,'(1P,D20.12,15D16.8)') &
                                                 time, &
                                  dtVPolRms,dtVTorRms, &
                                  CorPolRms,CorTorRms, &
                                    LFPolRms,LFTorRms, &
                                  AdvPolRms,AdvTorRms, &
                                  DifPolRms,DifTorRms, &
                                        BuoRms,PreRms, &
                            GeoRms/(CorPolRms+PreRms), &
                   MagRms/(CorPolRms+PreRms+LFPolRms), &
            ArcRms/(CorPolRms+PreRms+LFPolRms+BuoRms)
        END IF
        IF ( l_save_out) THEN
            CLOSE(n_dtvrms_file)
        END IF
        IF ( l_save_out) THEN
            OPEN(n_dtvasrms_file,FILE=dtvasrms_file,FORM='FORMATTED', &
                 STATUS='UNKNOWN',POSITION='APPEND')
        END IF
        IF ( l_RMStest ) THEN
            WRITE(n_dtvasrms_file,'(1P,D20.12,15D16.8)') &
                                                   time, &
                                dtVPolAsRms,dtVTorAsRms, &
                                CorPolAsRms,CorTorAsRms, &
                                  LFPolAsRms,LFTorAsRms, &
                                AdvPolAsRms,AdvTorAsRms, &
                                DifPolAsRms,DifTorAsRms, &
                                      BuoAsRms,PreAsRms, &
                                               GeoAsRms, &
                                               MagAsRms, &
                                               ArcAsRms
        ELSE
            IF ( PreAsRms/= 0d0 ) THEN
               bal1=GeoAsRms/(CorPolAsRms+PreAsRms)
               bal2=MagAsRms/(CorPolAsRms+PreAsRms+LFPolAsRms)
               bal3=ArcAsRms/(CorPolAsRms+PreAsRms+LFPolAsRms+BuoAsRms)
            ELSE
               bal1=0d0
               bal2=0d0
               bal3=0d0
            END IF
            WRITE(n_dtvasrms_file,'(1P,D20.12,15D16.8)') &
                                                   time, &
                                dtVPolAsRms,dtVTorAsRms, &
                                CorPolAsRms,CorTorAsRms, &
                                  LFPolAsRms,LFTorAsRms, &
                                AdvPolAsRms,AdvTorAsRms, &
                                DifPolAsRms,DifTorAsRms, &
                                      BuoAsRms,PreAsRms, &
                                         bal1,bal2,bal3
        END IF
        IF ( l_save_out) THEN
            CLOSE(n_dtvasrms_file)
        END IF

    END IF


    RETURN
    end SUBROUTINE dtVrms

!-----------------------------------------------------------------------------
