!$Id$
!*************************************************************************
    SUBROUTINE dtBrms(time)
!*************************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!-------------------------------------------------------------------------

!-------------------------------------------------------------------------

    USE truncation
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE blocking
    USE horizontal_data
    USE logic
    USE RMS
    USE dtB_mod
    USE output_data
    USE const
    USE parallel_mod
    USE integration, ONLY: rInt_R

    IMPLICIT NONE

!-- Input of variables:
    REAL(kind=8) :: time

!-- Local
    INTEGER :: nR,n,l1m0,lm
    CHARACTER(len=80) :: fileName
            
    REAL(kind=8) :: dtBPolRms,dtBPolAsRms
    REAL(kind=8) :: dtBTorRms,dtBTorAsRms
    REAL(kind=8) :: DstrRms
    REAL(kind=8) :: DadvRms
    REAL(kind=8) :: DdifRms
    REAL(kind=8) :: DdynRms
    REAL(kind=8) :: PdynRms,PdynAsRms
    REAL(kind=8) :: TdynRms,TdynAsRms
    REAL(kind=8) :: dummy1,dummy2,dummy3

!----- Note: five full additional fields needed here!
!      This may cause memory problems.
    COMPLEX(kind=8) :: PdynLM(lm_max_dtB,n_r_max_dtB)
    COMPLEX(kind=8) :: drPdynLM(lm_max_dtB,n_r_max_dtB)
    COMPLEX(kind=8) :: TdynLM(lm_max_dtB,n_r_max_dtB)
    COMPLEX(kind=8) :: workA(lm_max_dtB,n_r_max_dtB)
    COMPLEX(kind=8) :: workB(lm_max_dtB,n_r_max_dtB)

!-- For new movie output
    INTEGER :: nField,nFields,nFieldSize
    INTEGER :: nTheta,nThetaN,nThetaS,nThetaStart
    INTEGER :: nPos
    REAL(kind=8) :: dumm(12),rS
    REAL(kind=8) :: fOut(n_theta_max*n_r_max)
    REAL(kind=8) :: outBlock(nfs)
    CHARACTER(len=80) :: version
    LOGICAL :: lRmsMov


!-- end of declaration
!----------------------------------------------------------------------

    l1m0=lm2(1,0)

!--- Stretching
    CALL get_drNS(PstrLM,workA,lm_max_real,1,lm_max_real, &
                                n_r_max,n_cheb_max,workB, &
                           i_costf_init,d_costf_init,drx)
!--- Finalize rms poloidal and toroidal stretching:
    CALL get_PolTorRms(PstrLM,workA,TstrLM, &
                       PstrRms,TstrRms,PstrAsRms,TstrAsRms)
!--- Calculate dipole stretching and copy for total dynamo term:
    DO nR=1,n_r_max
        DO lm=1,lm_max
            PdynLM(lm,nR)  =PstrLM(lm,nR)
            drPdynLM(lm,nR)=workA(lm,nR)
            IF ( lm /= l1m0 ) THEN
                PstrLM(lm,nR)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                workA(lm,nR) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            END IF
        END DO
    END DO
!--- Get dipols stretching
    CALL get_PolTorRms(PstrLM,workA,TstrLM, &
                       DstrRms,dummy1,dummy2,dummy3)

!--- Finalize advection
    CALL get_drNS(PadvLM,workA,lm_max_real,1,lm_max_real, &
                                n_r_max,n_cheb_max,workB, &
                           i_costf_init,d_costf_init,drx)
    CALL get_PolTorRms(PadvLM,workA,TadvLM, &
                       PadvRms,TadvRms,PadvAsRms,TadvAsRms)
    DO nR=1,n_r_max
        DO lm=1,lm_max
            PstrLM(lm,nR)  =PdynLM(lm,nR)
            PdynLM(lm,nR)  =PdynLM(lm,nR)-PadvLM(lm,nR)
            drPdynLM(lm,nR)=drPdynLM(lm,nR)+workA(lm,nR)
            TdynLM(lm,nR)  =TstrLM(lm,nR)-TadvLM(lm,nR)
            IF ( lm /= l1m0 ) THEN
                PadvLM(lm,nR)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                workA(lm,nR) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            END IF
        END DO
    END DO
!--- Get dipole advection and total dynamo terms:
    CALL get_PolTorRms(PadvLM,workA,TadvLM, &
                       DadvRms,dummy1,dummy2,dummy3)
    CALL get_PolTorRms(PdynLM,drPdynLM,TdynLM, &
                       PdynRms,TdynRms,PdynAsRms,TdynAsRms)
    DO nR=1,n_r_max
        DO lm=1,lm_max
            PadvLM(lm,nR)=PdynLM(lm,nR)
            IF ( lm /= l1m0 ) THEN
                PdynLM(lm,nR)  =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                drPdynLM(lm,nR)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            END IF
        END DO
    END DO
!--- Get dipole dynamo terms:
    CALL get_PolTorRms(PdynLM,drPdynLM,TdynLM, &
                       DdynRms,dummy1,dummy2,dummy3)

!--- Diffusion:
    CALL get_drNS(PdifLM,workA,lm_max_real,1,lm_max_real, &
                                n_r_max,n_cheb_max,workB, &
                           i_costf_init,d_costf_init,drx)
    CALL get_PolTorRms(PdifLM,workA,TdifLM, &
                       PdifRms,TdifRms,PdifAsRms,TdifAsRms)
    DO nR=1,n_r_max
        DO lm=1,lm_max
            PdynLM(lm,nR)=PdifLM(lm,nR)
            IF ( lm /= l1m0 ) THEN
                PdifLM(lm,nR)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                workA(lm,nR) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            END IF
        END DO
    END DO
    CALL get_PolTorRms(PdifLM,workA,TdifLM, &
                       DdifRms,dummy1,dummy2,dummy3)


!--- Omega effect:
    CALL get_PolTorRms(PdifLM,workA,TomeLM, &
                       dummy1,TomeRms,dummy2,TomeAsRms)


!--- B changes:
    CALL get_drNS(dtBPolLMr,workA,lm_max_real,1,lm_max_real, &
                                   n_r_max,n_cheb_max,workB, &
                              i_costf_init,d_costf_init,drx)
    DO nR=1,n_r_max
        CALL hInt2dPol(workA(1,nR),2,lm_max, &
                       dtBPol2hInt(nR,1),dtBPolAs2hInt(nR,1))
        DO n=2,nThreadsLMmax
            dtBPol2hInt(nR,1)  =dtBPol2hInt(nR,1)   + &
                                dtBPol2hInt(nR,n)
            dtBPolAs2hInt(nR,1)=dtBPolAs2hInt(nR,1) + &
                                dtBPolAs2hInt(nR,n)
            dtBTor2hInt(nR,1)  =dtBTor2hInt(nR,1)   + &
                                dtBTor2hInt(nR,n)
            dtBTorAs2hInt(nR,1)=dtBTorAs2hInt(nR,1) + &
                                dtBTorAs2hInt(nR,n)
        END DO
    END DO
    dtBPolRms  =rInt_R(dtBPol2hInt(1,1),n_r_max, &
                       n_r_max,drx,i_costf_init,d_costf_init)
    dtBPolAsRms=rInt_R(dtBPolAs2hInt(1,1),n_r_max, &
                       n_r_max,drx,i_costf_init,d_costf_init)
    dtBTorRms  =rInt_R(dtBTor2hInt(1,1),n_r_max, &
                       n_r_max,drx,i_costf_init,d_costf_init)
    dtBTorAsRms=rInt_R(dtBTorAs2hInt(1,1),n_r_max, &
                       n_r_max,drx,i_costf_init,d_costf_init)
    dtBPolRms  =DSQRT(dtBPolRms  /vol_oc)
    dtBPolAsRms=DSQRT(dtBPolAsRms/vol_oc)
    dtBTorRms  =DSQRT(dtBTorRms  /vol_oc)
    dtBTorAsRms=DSQRT(dtBTorAsRms/vol_oc)

!-- OUTPUT:
    IF ( l_save_out) THEN
        OPEN(n_dtbrms_file,FILE=dtbrms_file,FORM='FORMATTED', &
             STATUS='UNKNOWN',POSITION='APPEND')
    END IF
    WRITE(n_dtbrms_file,'(1P,D20.12,12D16.8)') &
                                         time, &
                          dtBPolRms,dtBTorRms, &
                              PstrRms,TstrRms, &
                              PadvRms,TadvRms, &
                              PdifRms,TdifRms, &
                      TomeRms/TstrRms,TomeRms, &
                              PdynRms,TdynRms
    IF ( l_save_out) THEN
        CLOSE(n_dtbrms_file)
    END IF

    IF ( l_save_out) THEN
        OPEN(n_dtdrms_file,FILE=dtdrms_file,FORM='FORMATTED', &
             STATUS='UNKNOWN',POSITION='APPEND')
    END IF
    WRITE(n_dtdrms_file,'(1P,D20.12,3D16.8)') &
          time,DstrRms,DadvRms,DdifRms
    IF ( l_save_out) THEN
        CLOSE(n_dtdrms_file)
    END IF

!-- Output of movie files for axisymmetric toroidal field changes:
!   Tstr,Tome,Tdyn=Tstr+Tadv,
    lRmsMov=.FALSE.
    IF ( lRmsMov ) THEN
                   
        nFieldSize=n_theta_max*n_r_max
        nFields=7
        fileName='dtTas_mov.'//TAG
        OPEN(90,file=fileName,STATUS='UNKNOWN',FORM='UNFORMATTED')

    !------ Write header
        version='JW_Movie_Version_2'
        WRITE(90) version
        dumm(1)=112           ! type of input
        dumm(2)=3             ! marker for constant phi plane
        dumm(3)=0.D0          ! surface constant
        dumm(4)=nFields       ! no of fields
        WRITE(90) (real(dumm(n),4),n=1,4)

    !------ Define marker for output fields stored in movie field
        dumm(1)=101           ! Field marker for AS Br stretching
        dumm(2)=102           ! Field marker for AS Br dynamo term
        dumm(3)=103           ! Field marker for AS Br diffusion
        dumm(4)=104           ! Field marker for AS Bp stretching
        dumm(5)=105           ! Field marker for AS Bp dynamo term
        dumm(6)=106           ! Field marker for AS Bp omega effect
        dumm(7)=107           ! Field marker for AS Bp diffusion
        WRITE(90) (SNGL(dumm(n)),n=1,nFields)

    !------ Now other info about grid and parameters:
        WRITE(90) runid        ! run identifier
        dumm( 1)=n_r_max          ! total number of radial points
        dumm( 2)=n_r_max          ! no of radial point in outer core
        dumm( 3)=n_theta_max      ! no. of theta points
        dumm( 4)=n_phi_max        ! no. of phi points
        dumm( 5)=minc             ! imposed symmetry
        dumm( 6)=ra               ! control parameters
        dumm( 7)=ek               ! (for information only)
        dumm( 8)=pr               !      -"-
        dumm( 9)=prmag            !      -"-
        dumm(10)=radratio         ! ratio of inner / outer core
        dumm(11)=tScale           ! timescale
        WRITE(90) (SNGL(dumm(n)),     n=1,11)
        WRITE(90) (SNGL(r(n)/r_CMB),  n=1,n_r_max)
        WRITE(90) (SNGL(theta_ord(n)),n=1,n_theta_max)
        WRITE(90) (SNGL(phi(n)),      n=1,n_phi_max)

        dumm(1)=1    ! time frame number for movie
        dumm(2)=0.D0 ! time
        dumm(3)=0.D0
        dumm(4)=0.D0
        dumm(5)=0.D0
        dumm(6)=0.D0
        dumm(7)=0.D0
        dumm(8)=0.D0
        WRITE(90) (SNGL(dumm(n)),n=1,8)

    !------ Loop over different output field:
        Do nField=1,nFields

        !------ Loop over r and theta:
            DO nR=1,n_r_max ! Loop over radial points
                rS=r(nR)
                DO n=1,nThetaBs ! Loop over theta blocks
                    nThetaStart=(n-1)*sizeThetaB+1

                !------ Convert from lm to theta block and store in outBlock:
                    IF ( nField == 1 ) THEN
                        CALL get_RAS(PstrLM(1,nR),outBlock, &
                                     rS,nThetaStart,sizeThetaB)
                    ELSE IF ( nField == 2 ) THEN
                    ! Note that PadvLM stores PdynLM=PstrLM+PadvLM at this point!
                        CALL get_RAS(PadvLM(1,nR),outBlock, &
                                     rS,nThetaStart,sizeThetaB)
                    ELSE IF ( nField == 3 ) THEN
                    ! Note that PdynLM stores PdifLM at this point!
                        CALL get_RAS(PdynLM(1,nR),outBlock, &
                                     rS,nThetaStart,sizeThetaB)
                    ELSE IF ( nField == 4 ) THEN
                        CALL get_PASLM(TstrLM(1,nR),outBlock, &
                                       rS,nThetaStart,sizeThetaB)
                    ELSE IF ( nField == 5 ) THEN
                        CALL get_PASLM(TdynLM(1,nR),outBlock, &
                                       rS,nThetaStart,sizeThetaB)
                    ELSE IF ( nField == 6 ) THEN
                        CALL get_PASLM(TomeLM(1,nR),outBlock, &
                                       rS,nThetaStart,sizeThetaB)
                    ELSE IF ( nField == 7 ) THEN
                        CALL get_PASLM(TdifLM(1,nR),outBlock, &
                                       rS,nThetaStart,sizeThetaB)
                    END IF

                !------ Storage of field in fout for theta block
                    DO nTheta=1,sizeThetaB,2
                    !--------- Convert to correct order in theta grid points and store of fOut:
                        nThetaN=(nThetaStart+nTheta)/2
                        nPos=(nR-1)*n_theta_max+nThetaN
                        fOut(nPos)=outBlock(nTheta)
                        nThetaS=n_theta_max-nThetaN+1
                        nPos=(nR-1)*n_theta_max+nThetaS
                        fOut(nPos)=outBlock(nTheta+1)
                    END DO ! Loop over thetas in block

                END DO ! Loop over theta blocks

            END DO ! Loop over R

        !------ Output of field:
            WRITE(90) (SNGL(fOut(nPos)),nPos=1,nFieldSize)

        END DO ! Loop over different fields

    END IF ! output of mov fields ?


    RETURN
    end SUBROUTINE dtBrms

!----------------------------------------------------------------------------
