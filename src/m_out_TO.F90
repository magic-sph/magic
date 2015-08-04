!$Id$
module outTO_mod

   use truncation, only: n_r_max, n_r_maxStr, n_theta_maxStr, l_max, &
                         n_theta_max, n_phi_max, minc, lStressMem
   use radial_functions, only: r_ICB, i_costf_init, d_costf_init, r, &
                               r_CMB, orho1, drx
   use physical_parameters, only: ra, ek, pr, prmag, radratio, LFfac
   use torsional_oscillations, only: V2AS, Bs2AS, BspAS, BszAS, BpzAS, &
                                     BspdAS, BpsdAS, BzpdAS, BpzdAS,   &
                                     dzdVpLMr, dzddVpLMr, dzRstrLMr,   &
                                     dzAstrLMr, dzStrLMr, dzLFLMr,     &
                                     dzCorLMr
   use num_param, only: tScale
   use blocking, only: nThetaBs, sizeThetaB, nfs, lo_map
   use horizontal_data, only: phi, sinTheta, theta_ord, gauss
   use logic, only: lVerbose, l_save_out
   use output_data, only: sDens, zDens, TAG, log_file, runid, n_log_file
   use const, only: pi, vol_oc
   use LMLoop_data, only: llm, ulm
   use charmanip, only: dble2str
   use integration, only: rInt_R
   use plms_theta, only: plm_theta
   use TO_helpers, only: getPAStr, get_PAS, getAStr
   use useful, only: logWrite
 
   implicit none 

   private
   
   integer :: lmMaxS
   integer :: nZmaxA,nZmaxL
   integer :: nSmaxA,nSmaxL
 
   !-- Plms: Plm,sin
   real(kind=8), allocatable :: PlmS(:,:,:)
   real(kind=8), allocatable :: dPlmS(:,:,:)
   real(kind=8), allocatable :: OsinTS(:,:)
   real(kind=4), allocatable :: VpM(:,:), LFM(:,:), dVpM(:,:), AstrM(:,:)
   real(kind=4), allocatable :: RstrM(:,:), CorM(:,:), StrM(:,:), CLM(:,:)
   real(kind=8), allocatable :: d_costf_initZ(:,:), zZ(:,:), rZ(:,:)
   integer, allocatable :: i_costf_initZ(:,:), nZmaxS(:)

   public :: initialize_outTO_mod, outTO

contains

   subroutine initialize_outTO_mod

      nZmaxL=lStressMem*722
      nZmaxA=max0(2,nZmaxL)
      nSmaxL=lStressMem*625
      nSmaxA=max0(3,nSmaxL)
      lmMaxS = l_max+1

      allocate( PlmS(lmMaxS,nZmaxA/2+1,nSmaxA) )
      allocate( dPlmS(lmMaxS,nZmaxA/2+1,nSmaxA) )
      allocate( OsinTS(nZmaxA/2+1,nSmaxA) )
      allocate( vpM(nZmaxA/2,nSmaxA) )
      allocate( LFM(nZmaxA/2,nSmaxA) )
      allocate( dVpM(nZmaxA/2,nSmaxA) )
      allocate( AstrM(nZmaxA/2,nSmaxA) )
      allocate( RstrM(nZmaxA/2,nSmaxA) )
      allocate( CorM(nZmaxA/2,nSmaxA) )
      allocate( StrM(nZmaxA/2,nSmaxA) )
      allocate( CLM(nZmaxA/2,nSmaxA) )
      allocate( i_costf_initZ(2*nZmaxA+2,nSmaxA) )
      allocate( d_costf_initZ(2*nZmaxA+5,nSmaxA) )
      allocate( zZ(nZmaxA,nSmaxA) )
      allocate( rZ(nZmaxA/2+1,nSmaxA) )
      allocate( nZmaxS(nSmaxA) )

   end subroutine initialize_outTO_mod
!----------------------------------------------------------------------------
   subroutine outTO(time,n_time_step,                   &
        &           eKin,eKinTAS,nOutFile,nOutFile2,    &
        &           TOfileNhs,TOfileShs,movFile,tayFile,&
        &           nTOsets,nTOmovSets,nTOrmsSets,      &
        &           lTOmov,lTOrms,lTOZwrite,            &
        &           z,omega_ic,omega_ma)
      !-----------------------------------------------------------------------
      !   Output of axisymmetric zonal flow, its relative strength,
      !   its time variation, and all forces acting on it.
      !   The slowest part in the TO process is the repetitious calculation
      !   of plms by subroutine plm_theta. They are needed in getAStr and 
      !   getPAStr when I transform on the cylindrical grid. 
      !   The necessary plms could simply be calculated one and then 
      !   be stored for later use! See s_outTOnew.f.
      !-----------------------------------------------------------------------

      !-- Input of variables:
      real(kind=8),     intent(in) :: time
      integer,          intent(in) :: n_time_step
      real(kind=8),     intent(in) :: eKin, eKinTAS
      integer,          intent(in) :: nOutFile, nOutFile2
      character(len=*), intent(in) :: TOfileNhs, TOfileShs, movFile, tayFile
      logical,          intent(in) :: lTOmov
      complex(kind=8),  intent(in) :: z(llm:ulm,n_r_max)
      real(kind=8),     intent(in) :: omega_ic, omega_ma
      integer,          intent(inout) :: nTOsets, nTOmovSets, nTOrmsSets
      logical,          intent(inout) :: lTOrms, lTOZwrite

      !-- Output field:
      real(kind=8) :: fOut(n_theta_maxStr*n_r_maxStr)

      !-- Local variables:
      logical :: lTC,lStopRun

      !-- (l,r) Representation of the different contributions
      real(kind=8) :: dzVpLMr(lmMaxS,n_r_maxStr)
      real(kind=8) :: V2LMr(lmMaxS,n_r_maxStr)
      real(kind=8) :: Bs2LMr(lmMaxS,n_r_maxStr)
      real(kind=8) :: BszLMr(lmMaxS,n_r_maxStr)
      real(kind=8) :: BspLMr(lmMaxS,n_r_maxStr)
      real(kind=8) :: BpzLMr(lmMaxS,n_r_maxStr)
      real(kind=8) :: BspdLMr(lmMaxS,n_r_maxStr)
      real(kind=8) :: BpsdLMr(lmMaxS,n_r_maxStr)
      real(kind=8) :: BzpdLMr(lmMaxS,n_r_maxStr)
      real(kind=8) :: BpzdLMr(lmMaxS,n_r_maxStr)

      !---- Work array:
      complex(kind=8) :: workA(lmMaxS,n_r_maxStr)

      integer :: lm,l ! counter for degree

      integer :: nSmax,nS,nSI
      real(kind=8) :: zNorm  ! Norm z intervall
      integer :: nNorm  ! No. of grid points for norm intervall
      real(kind=8) :: zMin,zMax!,help ! integration boundarie, help variable
      logical :: lAS    ! .true. if axisymmetric (m=0) functions
      real(kind=8) :: sZ(nSmaxA),dsZ ! cylindrical radius s and s-step
      real(kind=8) :: h(nSmaxA),Oh(nSmaxA)
      real(kind=8) :: Os2(nSmaxA)

      integer :: nR     ! counter for radial grid point
      integer :: n      ! counter for theta blocks
      integer :: nOut,nFields ! counter for output fields
      integer :: nTheta ! counter for all thetas
      integer :: nThetaBlock ! counter for thetas in block
      integer :: nThetaOrd ! counter for ordered thetas
      integer :: nThetaNHS
      integer :: nThetaStart

      integer :: nZ,nZmax,nZmaxNS,nZmaxH!,nZP
      real(kind=8) :: VpS(nZmaxA)      
      real(kind=8) :: dVpS(nZmaxA)      
      real(kind=8) :: ddVpS(nZmaxA)      
      real(kind=8) :: V2S(nZmaxA)
      real(kind=8) :: LFS(nZmaxA)
      real(kind=8) :: CorS(nZmaxA)
      real(kind=8) :: RstrS(nZmaxA)
      real(kind=8) :: AstrS(nZmaxA)
      real(kind=8) :: StrS(nZmaxA)
      real(kind=8) :: Bs2S(nZmaxA)
      real(kind=8) :: BspS(nZmaxA)
      real(kind=8) :: BspdS(nZmaxA)
      real(kind=8) :: BpsdS(nZmaxA)
      real(kind=8) :: TayS(nZmaxA)
      real(kind=8) :: TayRS(nZmaxA)
      real(kind=8) :: TayVS(nZmaxA)

      real(kind=8) :: VpIntN(nSmaxA)  ,VpIntS(nSmaxA)    ! integration results
      real(kind=8) :: dVpIntN(nSmaxA) ,dVpIntS(nSmaxA)   ! integration results
      real(kind=8) :: ddVpIntN(nSmaxA),ddVpIntS(nSmaxA)  ! integration results
      real(kind=8) :: VpRIntN(nSmaxA) ,VpRIntS(nSmaxA)   ! for different s and 
      real(kind=8) :: V2IntS(nSmaxA)  ,V2IntN(nSmaxA)
      real(kind=8) :: LFIntN(nSmaxA)  ,LFIntS(nSmaxA)   
      real(kind=8) :: RstrIntN(nSmaxA),RstrIntS(nSmaxA) 
      real(kind=8) :: AstrIntN(nSmaxA),AstrIntS(nSmaxA) 
      real(kind=8) :: StrIntN(nSmaxA) ,StrIntS(nSmaxA) 
      real(kind=8) :: Bs2IntN(nSmaxA) ,Bs2IntS(nSmaxA)
      real(kind=8) :: BspIntN(nSmaxA) ,BspIntS(nSmaxA)
      real(kind=8) :: BspdIntN(nSmaxA),BspdIntS(nSmaxA)
      real(kind=8) :: BpsdIntN(nSmaxA),BpsdIntS(nSmaxA)
      real(kind=8) :: TayIntN(nSmaxA) ,TayIntS(nSmaxA) 
      real(kind=8) :: TayRIntN(nSmaxA),TayRIntS(nSmaxA) 
      real(kind=8) :: TayVIntN(nSmaxA),TayVIntS(nSmaxA) 
      real(kind=8) :: SVpIntN(nSmaxA) ,SVpIntS(nSmaxA)   ! help arrays and values for 
      real(kind=8) :: SBs2IntN(nSmaxA),SBs2IntS(nSmaxA)  ! differentiation in s   
      real(kind=8) :: SBspIntN(nSmaxA),SBspIntS(nSmaxA)
      real(kind=8) :: dSVpIntN, dSVpIntS
      real(kind=8) :: d2SVpIntN,d2SVpIntS
      real(kind=8) :: dSBspIntN,dSBspIntS
      real(kind=8) :: dSBs2IntN,dSBs2IntS
      real(kind=8) :: TauN(nSmaxA),TauS(nSmaxA)          ! Taylor integral
      real(kind=8) :: TauBN(nSmaxA),TauBS(nSmaxA)       
      real(kind=8) :: dTauBN(nSmaxA),dTauBS(nSmaxA)    
      real(kind=8) :: dTTauN(nSmaxA),dTTauS(nSmaxA)      ! time change of Tau...
      real(kind=8) :: dTTauBN(nSmaxA),dTTauBS(nSmaxA)   

      !-- For integration along z:
      real(kind=8) :: zALL(2*nZmaxA)
      real(kind=8) :: thetaZ
      real(kind=8) :: fac
      real(kind=8) :: vSF,fSF,f1,f2

      !-- For boundaries:
      real(kind=8) :: rB(2),BspB(2),BspdB(2),BpsdB(2)
      real(kind=8) :: Bs2B(2),BszB(2),BpzB(2),BzpdB(2),BpzdB(2)

      !-- For sphere integration:
      real(kind=8) :: CorR(n_r_max)
      real(kind=8) :: RstrR(n_r_max)
      real(kind=8) :: AstrR(n_r_max)
      real(kind=8) :: StrR(n_r_max)
      real(kind=8) :: VpR(n_r_max),dVpR(n_r_max)
      real(kind=8) :: LFR(n_r_max),LFABSR(n_r_max)
      real(kind=8) :: TayR(n_r_max),TayRMSR(n_r_max)
      real(kind=8) :: TayRMS,TaySRMS
      real(kind=8) :: TayRRMS,TayVRMS
      real(kind=8) :: VpRMS,VRMS,VgRMS
      real(kind=8) :: rS,sS
      real(kind=8) :: outBlock(nfs)
      real(kind=8) :: timeLast!,tNorm
      real(kind=4) :: timeAve,dt
      SAVE timeLast,timeAve

      character(len=255) :: message
      character(len=64) :: version,fileName
      integer :: nFieldSize,nPos

      !integer :: nLines
      real(kind=8),EXTERNAL :: chebInt

      real(kind=8) :: dumm(12)

      !-- Huge arrays for time average ....
      logical :: l_TOZave

      !-- For TOZ output files:
      integer,SAVE :: nTOZfile!,length
      character(len=10) :: string

      if ( lVerbose ) write(*,*) '! Starting outTO!'

      nTOsets=nTOsets+1

      l_TOZave=.true.

      !--- Rescaling for rotation time scale and planetary radius
      !    length scale, for the velocity I use the Rossby number
      !       vSF=ek/r_CMB                   
      !       fSF=ek*ek/(4.D0*pi**2*r_CMB)  
      vSF=1.D0
      fSF=1.D0

      nFieldSize=n_theta_maxStr*n_r_maxStr

      !-- Start with calculating advection due to axisymmetric flows:

      zNorm=1.D0               ! This is r_CMB-r_ICB
      nNorm=int(zDens*n_r_max) ! Covered with nNorm  points !
      nSmax=n_r_max+int(r_ICB*dble(n_r_max))
      nSmax=int(sDens*nSmax)
      if ( nSmax > nSmaxA ) then
         write(*,*) 'Increase nSmaxA in ouTO!'
         write(*,*) 'Should be at least nSmax=',nSmax
         stop
      end if
      lAS=.true.

      !--- Transform to lm-space for all radial grid points:

      do nR=1,n_r_max
         do n=1,nThetaBs
            nThetaStart=(n-1)*sizeThetaB+1
            call legtfAS(V2LMr(1,nR),V2AS(nThetaStart,nR),            &
                 &               l_max+1,nThetaStart,sizeThetaB)
            call legtfAS2(Bs2LMr(1,nR),BszLMr(1,nR),                            &
                 &               Bs2AS(nThetaStart,nR),BszAS(nThetaStart,nR),   &
                 &               l_max+1,nThetaStart,sizeThetaB)
            call legtfAS2(BspLMr(1,nR),BpzLMr(1,nR),                            &
                 &               BspAS(nThetaStart,nR),BpzAS(nThetaStart,nR),   &
                 &               l_max+1,nThetaStart,sizeThetaB)
            call legtfAS2(BspdLMr(1,nR),BpsdLMr(1,nR),                          &
                 &               BspdAS(nThetaStart,nR),BpsdAS(nThetaStart,nR), &
                 &               l_max+1,nThetaStart,sizeThetaB)
            call legtfAS2(BzpdLMr(1,nR),BpzdLMr(1,nR),                          &
                 &               BzpdAS(nThetaStart,nR),BpzdAS(nThetaStart,nR), &
                 &               l_max+1,nThetaStart,sizeThetaB)
         end do
      end do ! Loop over radial grid points

      do nR=1,n_r_max
         do l=1,l_max
            lm=lo_map%lm2(l,0)
            dzVpLMr(l+1,nR)=real(z(lm,nR))
         end do
      end do

      !---- Transform the contributions to cheb space for z-integral:
      call costf1(dzVpLMr,lmMaxS,1,lmMaxS,workA,i_costf_init,d_costf_init)
      call costf1(V2LMr,lmMaxS,1,lmMaxS,workA,i_costf_init,d_costf_init)
      call costf1(dzdVpLMr,lmMaxS,1,lmMaxS,workA,i_costf_init,d_costf_init)
      call costf1(dzddVpLMr,lmMaxS,1,lmMaxS,workA,i_costf_init,d_costf_init)
      call costf1(Bs2LMr,lmMaxS,1,lmMaxS,workA,i_costf_init,d_costf_init)
      call costf1(BszLMr,lmMaxS,1,lmMaxS,workA,i_costf_init,d_costf_init)
      call costf1(BspLMr,lmMaxS,1,lmMaxS,workA,i_costf_init,d_costf_init)
      call costf1(BpzLMr,lmMaxS,1,lmMaxS,workA,i_costf_init,d_costf_init)
      call costf1(dzRstrLMr,lmMaxS,1,lmMaxS,workA,i_costf_init,d_costf_init)
      call costf1(dzAstrLMr,lmMaxS,1,lmMaxS,workA,i_costf_init,d_costf_init)
      call costf1(dzStrLMr,lmMaxS,1,lmMaxS,workA,i_costf_init,d_costf_init)
      call costf1(dzLFLMr,lmMaxS,1,lmMaxS,workA,i_costf_init,d_costf_init)
      call costf1(dzCorLMr,lmMaxS,1,lmMaxS,workA,i_costf_init,d_costf_init)
      call costf1(BspdLMr,lmMaxS,1,lmMaxS,workA,i_costf_init,d_costf_init)
      call costf1(BpsdLMr,lmMaxS,1,lmMaxS,workA,i_costf_init,d_costf_init)
      call costf1(BpzdLMr,lmMaxS,1,lmMaxS,workA,i_costf_init,d_costf_init)
      call costf1(BzpdLMr,lmMaxS,1,lmMaxS,workA,i_costf_init,d_costf_init)

      dsZ   =r_CMB/dble(nSmax)  ! Step in s controlled by nSmax
      nSI   =0
      do nS=1,nSmax
         sZ(nS)=(nS-0.5D0)*dsZ
         if ( sZ(nS) < r_ICB .and. nS > nSI ) nSI=nS
      end do

      if ( nTOsets == 1 ) nTOZfile=0
      if ( lTOZwrite ) then
         nTOZfile=nTOZfile+1
         call dble2str(dble(nTOZfile),string)
         fileName='TOZ_'//trim(adjustl(string))//'.'//TAG
         open(95, file=fileName, form='unformatted', status='unknown')
         write(95) sngl(time), float(nSmax), sngl(omega_ic), sngl(omega_ma)
         write(95) (sngl(sZ(nS)),nS=1,nSmax)
      end if
      if ( nTOsets > 1 .and. l_TOZave ) then
         fileName='TOZM.'//TAG
         open(96,file=fileName, form='unformatted', status='unknown')
         write(96) float(nSmax), sngl(omega_ic), sngl(omega_ma)
         write(96) (sngl(sZ(nS)),nS=1,nSmax)
      end if

      lStopRun=.false.
      do nS=1,nSmax

         !------ Get integral boundaries for this s:
         zMax=dsqrt(r_CMB*r_CMB-sZ(nS)*sZ(nS))
         if ( sZ(nS) < r_ICB ) then
            lTC=.true.
            zMin=dsqrt(r_ICB*r_ICB-sZ(nS)*sZ(nS))
         else
            lTC=.false.
            zMin=-zMax
         end if
         h(nS) =zMax-zMin
         Oh(nS)=1.D0/h(nS)

         if ( nTOsets == 1 ) then
            !------ Initialize integration for NHS:
            !       Each processor calculates Cheb transform data
            !       for HIS nS and the Plms along the Cylinder
            !       chebIntInit returns zZ,nZmaxS,i_costf_initZ and
            !       d_costfInitZ:
            !       Note that this returns z in the MAGIC way, 
            !       starting with zMax, ending with zMin
            !       z(1,nS)=zMin, z(nZmax,nS)=zMax
            call chebIntInit(zMin,zMax,zNorm,nNorm,                   &
                 &                    nZmaxA,zZ(1,nS),nZmaxS(nS),     &
                 &       i_costf_initZ(1,nS),d_costf_initZ(1,nS))

            !--- Points in nothers halfsphere
            if ( lTC ) then
               nZmax=nZmaxS(nS)  ! nZmax point in each polar region
               if ( 2*nZmax > nZmaxA ) then 
                  write(*,*) '! nZmaxA too small in outTO!'
                  write(*,*) '! Should be at least:',2*nZmax
                  lStopRun=.true.
               end if
            else
               nZmax=(nZmaxS(nS)-1)/2+1 ! Odd point number !
               ! all together nZmaxS(nS) from
               ! south to north including equator
               if ( nZmaxS(nS) > nZmaxA ) then 
                  write(*,*) '! nZmaxA too small in outTO!'
                  write(*,*) '! Should be at least:',nZmaxS(nS)
                  lStopRun=.true.
               end if
            end if
            if ( lStopRun ) GOTO 99
            do nZ=1,nZmax
               rZ(nZ,nS)    =dsqrt(zZ(nZ,nS)**2+sZ(nS)**2)
               thetaZ       =DATAN2(sZ(nS),zZ(nZ,nS))
               OsinTS(nZ,nS)=1.D0/DSIN(thetaZ)
               call plm_theta(thetaZ,l_max,0,minc,                    &
                    &            PlmS(1,nZ,nS),dPlmS(1,nZ,nS),lmMaxS,2)
            end do
         end if

         !--------- Get the flow components for all northern and
         !          southern thetas and all phis:

         nZmax=nZmaxS(nS)
         if ( lTC ) then
            nZmaxNS=2*nZmax
            nZmaxH =nZmax
            do nZ=1,nZmax
               zALL(nZ)=zZ(nZ,nS)
               zALL(nZmaxNS-nZ+1)=-zZ(nZ,nS)
            end do
         else
            nZmaxNS=nZmax
            nZmaxH =(nZmax-1)/2+1
            do nZ=1,nZmax
               zALL(nZ)=zZ(nZ,nS)
            end do
         end if

         call getPAStr(VpS,dzVpLMr,nZmaxNS,nZmaxA,lmMaxS,             &
              &                      l_max,r_ICB,r_CMB,n_r_max,       &
              &                 rZ(1,nS),dPlmS(1,1,nS),OsinTS(1,nS))
         call getPAStr(dVpS,dzdVpLMr,nZmaxNS,nZmaxA,lmMaxS,           &
              &                        l_max,r_ICB,r_CMB,n_r_max,     &
              &                   rZ(1,nS),dPlmS(1,1,nS),OsinTS(1,nS))
         call getPAStr(ddVpS,dzddVpLMr,nZmaxNS,nZmaxA,lmMaxS,         &
              &                          l_max,r_ICB,r_CMB,n_r_max,   &
              &                     rZ(1,nS),dPlmS(1,1,nS),OsinTS(1,nS))
         call getPAStr(RstrS,dzRstrLMr,nZmaxNS,nZmaxA,lmMaxS,         &
              &                          l_max,r_ICB,r_CMB,n_r_max,   &
              &                     rZ(1,nS),dPlmS(1,1,nS),OsinTS(1,nS))
         call getPAStr(AstrS,dzAstrLMr,nZmaxNS,nZmaxA,lmMaxS,         &
              &                          l_max,r_ICB,r_CMB,n_r_max,   &
              &                     rZ(1,nS),dPlmS(1,1,nS),OsinTS(1,nS))
         call getPAStr(StrS,dzStrLMr,nZmaxNS,nZmaxA,lmMaxS,           &
              &                        l_max,r_ICB,r_CMB,n_r_max,     &
              &                   rZ(1,nS),dPlmS(1,1,nS),OsinTS(1,nS))
         call getPAStr(LFS,dzLFLMr,nZmaxNS,nZmaxA,lmMaxS,             &
              &                      l_max,r_ICB,r_CMB,n_r_max,       &
              &                 rZ(1,nS),dPlmS(1,1,nS),OsinTS(1,nS))
         call getPAStr(CorS,dzCorLMr,nZmaxNS,nZmaxA,lmMaxS,           &
              &                        l_max,r_ICB,r_CMB,n_r_max,     &
              &                   rZ(1,nS),dPlmS(1,1,nS),OsinTS(1,nS))
         do nZ=1,nZmaxNS
            TayS(nZ) =dabs(LFS(nZ))
            TayRS(nZ)=dabs(RstrS(nZ))
            TayVS(nZ)=dabs(StrS(nZ))
         end do

         call getAStr(V2S,V2LMr,nZmaxNS,nZmaxA,lmMaxS,                &
              &                   l_max,r_ICB,r_CMB,n_r_max,          &
              &                            rZ(1,nS),PlmS(1,1,nS))
         call getAStr(Bs2S,Bs2LMr,nZmaxNS,nZmaxA,lmMaxS,              &
              &                     l_max,r_ICB,r_CMB,n_r_max,        &
              &                              rZ(1,nS),PlmS(1,1,nS))
         call getAStr(BspS,BspLMr,nZmaxNS,nZmaxA,lmMaxS,              &
              &                     l_max,r_ICB,r_CMB,n_r_max,        &
              &                              rZ(1,nS),PlmS(1,1,nS))
         call getAStr(BspdS,BspdLMr,nZmaxNS,nZmaxA,lmMaxS,            &
              &                       l_max,r_ICB,r_CMB,n_r_max,      &
              &                                rZ(1,nS),PlmS(1,1,nS))
         call getAStr(BpsdS,BpsdLMr,nZmaxNS,nZmaxA,lmMaxS,            &
              &                       l_max,r_ICB,r_CMB,n_r_max,      &
              &                                rZ(1,nS),PlmS(1,1,nS))

         if ( l_TOZave ) then
            if ( nTOsets == 1 ) then
               timeAve=1.E0
               do nZ=1,nZmaxNS
                  VpM(nZ,nS)  =sngl(VpS(nZ))
                  dVpM(nZ,nS) =sngl(dVpS(nZ))
                  LFM(nZ,nS)  =sngl(LFfac*LFS(nZ))
                  RstrM(nZ,nS)=sngl(RstrS(nZ))
                  AstrM(nZ,nS)=sngl(AstrS(nZ))
                  StrM(nZ,nS) =sngl(StrS(nZ))
                  CorM(nZ,nS) =sngl(CorS(nZ))
                  CLM(nZ,nS)  =sngl(CorS(nZ)+LFfac*LFS(nZ))
               end do
            else if ( nTOsets == 2 ) then
               dt=sngl(time-timeLast)
               timeAve=dt
               do nZ=1,nZmaxNS
                  VpM(nZ,nS)  =dt*(VpM(nZ,nS)  +sngl(VpS(nZ)))
                  dVpM(nZ,nS) =dt*(dVpM(nZ,nS) +sngl(dVpS(nZ)))
                  LFM(nZ,nS)  =dt*(LFM(nZ,nS)  +sngl(LFfac*LFS(nZ)))
                  RstrM(nZ,nS)=dt*(RstrM(nZ,nS)+sngl(RstrS(nZ)))
                  AstrM(nZ,nS)=dt*(AstrM(nZ,nS)+sngl(AstrS(nZ)))
                  StrM(nZ,nS) =dt*(StrM(nZ,nS) +sngl(StrS(nZ)))
                  CorM(nZ,nS) =dt*(CorM(nZ,nS) +sngl(CorS(nZ)))
                  CLM(nZ,nS)  =dt*(CLM(nZ,nS)  +sngl(CorS(nZ)+LFfac*LFS(nZ)))
               end do
            else
               dt=sngl(time-timeLast)
               timeAve=timeAve+dt
               do nZ=1,nZmaxNS
                  VpM(nZ,nS)  =VpM(nZ,nS)  +dt*sngl(VpS(nZ))
                  dVpM(nZ,nS) =dVpM(nZ,nS) +dt*sngl(dVpS(nZ))
                  LFM(nZ,nS)  =LFM(nZ,nS)  +dt*sngl(LFfac*LFS(nZ))
                  RstrM(nZ,nS)=RstrM(nZ,nS)+dt*sngl(RstrS(nZ))
                  AstrM(nZ,nS)=AstrM(nZ,nS)+dt*sngl(AstrS(nZ))
                  StrM(nZ,nS) =StrM(nZ,nS) +dt*sngl(StrS(nZ))
                  CorM(nZ,nS) =CorM(nZ,nS) +dt*sngl(CorS(nZ))
                  CLM(nZ,nS)  =CLM(nZ,nS)  +dt*sngl(CorS(nZ)+LFfac*LFS(nZ))
               end do
            end if
         end if

         if ( l_TOZave .and. nTOsets > 1 ) then
            write(96) float(nZmaxNS)
            write(96) (sngl(zALL(nZ))      ,nZ=1,nZmaxNS), &
                 &    (VpM(nZ,nS)/timeAve  ,nZ=1,nZmaxNS), &
                 &    (dVpM(nZ,nS)/timeAve ,nZ=1,nZmaxNS), &
                 &    (RstrM(nZ,nS)/timeAve,nZ=1,nZmaxNS), &
                 &    (AstrM(nZ,nS)/timeAve,nZ=1,nZmaxNS), &
                 &    (LFM(nZ,nS)/timeAve  ,nZ=1,nZmaxNS), &
                 &    (StrM(nZ,nS)/timeAve ,nZ=1,nZmaxNS), &
                 &    (CorM(nZ,nS)/timeAve ,nZ=1,nZmaxNS), &
                 &    (CLM(nZ,nS)/timeAve  ,nZ=1,nZmaxNS)
         end if
         if ( lTOZwrite ) then
            write(95) float(nZmaxNS)
            write(95) (sngl(zALL(nZ)) ,nZ=1,nZmaxNS),      &
                 &    (sngl(VpS(nZ))  ,nZ=1,nZmaxNS),      &
                 &    (sngl(dVpS(nZ)) ,nZ=1,nZmaxNS),      &
                 &    (sngl(RstrS(nZ)),nZ=1,nZmaxNS),      &
                 &    (sngl(AstrS(nZ)),nZ=1,nZmaxNS),      &
                 &    (sngl(LFfac*LFS(nZ)),nZ=1,nZmaxNS),  &
                 &    (sngl(StrS(nZ)) ,nZ=1,nZmaxNS),      &
                 &    (sngl(CorS(nZ)) ,nZ=1,nZmaxNS)
         end if

         !--- Z-integrals:
         VpIntN(nS)  =chebInt(VpS,zMin,zMax,nZmax,nZmaxA,             &
              &             i_costf_initZ(1,nS),d_costf_initZ(1,nS))
         dVpIntN(nS) =chebInt(dVpS,zMin,zMax,nZmax,nZmaxA,            &
              &              i_costf_initZ(1,nS),d_costf_initZ(1,nS))
         ddVpIntN(nS)=chebInt(ddVpS,zMin,zMax,nZmax,nZmaxA,           &
              &               i_costf_initZ(1,nS),d_costf_initZ(1,nS))
         LFIntN(nS)  =chebInt(LFS,zMin,zMax,nZmax,nZmaxA,             &
              &             i_costf_initZ(1,nS),d_costf_initZ(1,nS))
         TayIntN(nS) =chebInt(TayS,zMin,zMax,nZmax,nZmaxA,            &
              &              i_costf_initZ(1,nS),d_costf_initZ(1,nS))
         TayRIntN(nS)=chebInt(TayRS,zMin,zMax,nZmax,nZmaxA,           &
              &              i_costf_initZ(1,nS),d_costf_initZ(1,nS))
         TayVIntN(nS)=chebInt(TayVS,zMin,zMax,nZmax,nZmaxA,           &
              &              i_costf_initZ(1,nS),d_costf_initZ(1,nS))
         RstrIntN(nS)=chebInt(RstrS,zMin,zMax,nZmax,nZmaxA,           &
              &               i_costf_initZ(1,nS),d_costf_initZ(1,nS))
         AstrIntN(nS)=chebInt(AstrS,zMin,zMax,nZmax,nZmaxA,           &
              &               i_costf_initZ(1,nS),d_costf_initZ(1,nS))
         StrIntN(nS) =chebInt(StrS,zMin,zMax,nZmax,nZmaxA,            &
              &              i_costf_initZ(1,nS),d_costf_initZ(1,nS))
         V2IntN(nS)  =chebInt(V2S,zMin,zMax,nZmax,nZmaxA,             &
              &             i_costf_initZ(1,nS),d_costf_initZ(1,nS))
         Bs2IntN(nS) =chebInt(Bs2S,zMin,zMax,nZmax,nZmaxA,            &
              &              i_costf_initZ(1,nS),d_costf_initZ(1,nS))
         BspIntN(nS) =chebInt(BspS,zMin,zMax,nZmax,nZmaxA,            &
              &              i_costf_initZ(1,nS),d_costf_initZ(1,nS))
         BspdIntN(nS)=chebInt(BspdS,zMin,zMax,nZmax,nZmaxA,           &
              &               i_costf_initZ(1,nS),d_costf_initZ(1,nS))
         BpsdIntN(nS)=chebInt(BpsdS,zMin,zMax,nZmax,nZmaxA,           &
              &               i_costf_initZ(1,nS),d_costf_initZ(1,nS))

         if ( V2IntN(nS) < 0.D0 ) then
            VpRIntN(nS)=1.D0
         else
            VpRIntN(nS)=dabs(VpIntN(nS))/dsqrt(V2IntN(nS))
         end if
         VpRIntN(nS) =Dmin1(1.D0,VpRIntN(nS))
         TayIntN(nS) =LFIntN(nS)/TayIntN(nS)
         TayRIntN(nS)=RstrIntN(nS)/TayRIntN(nS)
         TayVIntN(nS)=StrIntN(nS)/TayVIntN(nS)

         !--- Z-integration inside northern TC:
         if ( lTC ) then
            VpIntS(nS)  =chebInt(VpS(nZmax+1),zMin,zMax,nZmax,        &
                 &           nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
            dVpIntS(nS) =chebInt(dVpS(nZmax+1),zMin,zMax,nZmax,       &
                 &            nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
            ddVpIntS(nS)=chebInt(ddVpS(nZmax+1),zMin,zMax,nZmax,      &
                 &             nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
            LFIntS(nS)  =chebInt(LFS(nZmax+1),zMin,zMax,nZmax,        &
                 &           nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
            TayIntS(nS) =chebInt(TayS(nZmax+1),zMin,zMax,nZmax,       &
                 &            nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
            TayRIntS(nS)=chebInt(TayRS(nZmax+1),zMin,zMax,nZmax,      &
                 &            nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
            TayVIntS(nS)=chebInt(TayVS(nZmax+1),zMin,zMax,nZmax,      &
                 &            nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
            RstrIntS(nS)=chebInt(RstrS(nZmax+1),zMin,zMax,nZmax,      &
                 &             nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
            AstrIntS(nS)=chebInt(AstrS(nZmax+1),zMin,zMax,nZmax,      &
                 &             nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
            StrIntS(nS) =chebInt(StrS(nZmax+1),zMin,zMax,nZmax,       &
                 &            nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
            V2IntS(nS)  =chebInt(V2S(nZmax+1),zMin,zMax,nZmax,        &
                 &           nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
            Bs2IntS(nS) =chebInt(Bs2S(nZmax+1),zMin,zMax,nZmax,       &
                 &            nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
            BspIntS(nS) =chebInt(BspS(nZmax+1),zMin,zMax,nZmax,       &
                 &            nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
            BspdIntS(nS)=chebInt(BspdS(nZmax+1),zMin,zMax,nZmax,      &
                 &             nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
            BpsdIntS(nS)=chebInt(BpsdS(nZmax+1),zMin,zMax,nZmax,      &
                 &             nZmaxA,i_costf_initZ(1,nS),d_costf_initZ(1,nS))
            if ( V2IntS(nS) < 0.D0 ) then
               VpRIntS(nS)=1.D0
            else
               VpRIntS(nS)=dabs(VpIntS(nS))/dsqrt(V2IntS(nS))
            end if
            VpRIntS(nS) =dmin1(1.D0,VpRIntS(nS))
            TayIntS(nS) =LFIntS(nS)/TayIntS(nS)
            TayRIntS(nS)=RstrIntS(nS)/TayRIntS(nS)
            TayVIntS(nS)=StrIntS(nS)/TayVIntS(nS)
         else
            VpIntS(nS)  =VpIntN(nS)
            dVpIntS(nS) =dVpIntN(nS)
            ddVpIntS(nS)=ddVpIntN(nS)
            LFIntS(nS)  =LFIntN(nS)
            RstrIntS(nS)=RstrIntN(nS)
            AstrIntS(nS)=AstrIntN(nS)
            StrIntS(nS) =StrIntN(nS)
            V2IntS(nS)  =V2IntN(nS)
            Bs2IntS(nS) =Bs2IntN(nS)
            BspIntS(nS) =BspIntN(nS)
            BspdIntS(nS)=BspdIntN(nS)
            BpsdIntS(nS)=BpsdIntN(nS)
            VpRIntS(nS) =VpRIntN(nS)
            TayIntS(nS) =TayIntN(nS)                 
            TayRIntS(nS)=TayRIntN(nS)                 
            TayVIntS(nS)=TayVIntN(nS)                 
         end if

         !------ Boundary Values:
         if ( lTC ) then ! inside TC
            rB(1)=r_CMB
            zMax = dsqrt(r_CMB**2-sZ(nS)**2)
            zMin =-Zmax
            BspB(1) =BspS(1)
            BspB(2) =BspS(nZmaxNS)
            BspdB(1)=BspdS(1)
            BspdB(2)=BspdS(nZmaxNS)
            BpsdB(1)=BpsdS(1)
            BpsdB(2)=BpsdS(nZmaxNS)
            Bs2B(1) =Bs2S(1)
            Bs2B(2) =Bs2S(nZmaxNS)
            call getAStr(BszB,BszLMr,2,2,lmMaxS,l_max,           &
                 &                   r_ICB,r_CMB,n_r_max,rB,PlmS(1,1,nS))
            call getAStr(BpzB,BpzLMr,2,2,lmMaxS,l_max,           &
                 &                   r_ICB,r_CMB,n_r_max,rB,PlmS(1,1,nS))
            call getAStr(BzpdB,BzpdLMr,2,2,lmMaxS,l_max,         &
                 &                     r_ICB,r_CMB,n_r_max,rB,PlmS(1,1,nS))
            call getAStr(BpzdB,BpzdLMr,2,2,lmMaxS,l_max,         &
                 &                     r_ICB,r_CMB,n_r_max,rB,PlmS(1,1,nS))

            TauBS(nS)  =-(BpzB(2)+sZ(nS)/zMin*BspB(2))
            dTauBS(nS) =-(BpzdB(2)+BzpdB(2) + sZ(nS)/zMin*(BspdB(2)+BpsdB(2)))
            dTTauBS(nS)=-(BszB(2)+sZ(nS)/zMin*Bs2B(2))
            TauBN(nS)  =  BpzB(1)+sZ(nS)/zMax*BspB(1)
            dTauBN(nS) =  BpzdB(1)+BzpdB(1) + sZ(nS)/zMax*(BspdB(1)+BpsdB(1))
            dTTauBN(nS)=  BszB(1)+sZ(nS)/zMax*Bs2B(1)

            rB(1)=r_ICB 
            zMax = dsqrt(r_ICB**2-sZ(nS)**2)
            zMin =-zMax
            BspB(1) =BspS(nZmax)
            BspB(2) =BspS(nZmax+1)
            BspdB(1)=BspdS(nZmax)
            BspdB(2)=BspdS(nZmax+1)
            BpsdB(1)=BpsdS(nZmax)
            BpsdB(2)=BpsdS(nZmax+1)
            Bs2B(1) =Bs2S(nZmax)
            Bs2B(2) =Bs2S(nZmax+1)
            call getAStr(BszB,BszLMr,2,2,lmMaxS,l_max,           &
                 &               r_ICB,r_CMB,n_r_max,rB,PlmS(1,nZmax,nS))
            call getAStr(BpzB,BpzLMr,2,2,lmMaxS,l_max,           &
                 &               r_ICB,r_CMB,n_r_max,rB,PlmS(1,nZmax,nS))
            call getAStr(BzpdB,BzpdLMr,2,2,lmMaxS,l_max,         &
                 &                 r_ICB,r_CMB,n_r_max,rB,PlmS(1,nZmax,nS))
            call getAStr(BpzdB,BpzdLMr,2,2,lmMaxS,l_max,         &
                 &                 r_ICB,r_CMB,n_r_max,rB,PlmS(1,nZmax,nS))

            TauBS(nS)  =TauBS(nS)  +(BpzB(2)+sZ(nS)/zMin*BspB(2))
            dTauBS(nS) =dTauBS(nS) +(BpzdB(2)+BzpdB(2) +              &
                 &                    sZ(nS)/zMin*(BspdB(2)+BpsdB(2)))
            dTTauBS(nS)=dTTauBS(nS)+(BszB(2)+sZ(nS)/zMin*Bs2B(2))
            TauBN(nS)  =TauBN(nS)  -(BpzB(1) +sZ(nS)/zMax*BspB(1))
            dTauBN(nS) =dTauBN(nS) -(BpzdB(1)+BzpdB(1) +              &
                 &                    sZ(nS)/zMax*(BspdB(1)+BpsdB(1)))    
            dTTauBN(nS)=dTTauBN(nS)-(BszB(1)+sZ(nS)/zMax*Bs2B(1))

         else ! outside TC

            rB(1)=r_CMB
            zMax = dsqrt(r_CMB**2-sZ(nS)**2)
            zMin =-Zmax
            BspB(1)=BspS(1)
            BspB(2)=BspS(nZmaxNS)
            BspdB(1)=BspdS(1)
            BspdB(2)=BspdS(nZmaxNS)
            BpsdB(1)=BpsdS(1)
            BpsdB(2)=BpsdS(nZmaxNS)
            Bs2B(1) =Bs2S(1)
            Bs2B(2) =Bs2S(nZmaxNS)
            call getAStr(BpzB,BpzLMr,2,2,lmMaxS,l_max,           &
                 &                   r_ICB,r_CMB,n_r_max,rB,PlmS(1,1,nS))
            call getAStr(BzpdB,BzpdLMr,2,2,lmMaxS,l_max,         &
                 &                     r_ICB,r_CMB,n_r_max,rB,PlmS(1,1,nS))
            call getAStr(BpzdB,BpzdLMr,2,2,lmMaxS,l_max,         &
                 &                     r_ICB,r_CMB,n_r_max,rB,PlmS(1,1,nS))

            TauBS(nS)  = BpzB(1)+sZ(nS)/zMax*BspB(1) - BpzB(2)-sZ(nS)/zMin*BspB(2)
            dTauBS(nS) = BpzdB(1)+BzpdB(1) + sZ(nS)/zMax*(BspdB(1)+BpsdB(1)) -   &
                 &                     BpzdB(2)-BzpdB(2) -                       &
                 &                     sZ(nS)/zMin*(BspdB(2)+BpsdB(2))
            dTTauBS(nS)= BszB(1)+sZ(nS)/zMax*Bs2B(1) - BszB(2)-sZ(nS)/zMin*Bs2B(2)
            TauBN(nS)  =TauBS(nS)
            dTauBN(nS) =dTauBS(nS)
            dTTauBN(nS)=dTTauBS(nS)

         end if

99     continue

      end do  ! Loop over s 
      ! Integration finished

      close (95)
      if ( l_TOZave .and. nTOsets > 1 ) close (96)

      if ( lStopRun ) stop

      !--- Integrate Geostrophic azumithal flow energy and Taylor measure:
      VgRMS  =0.D0
      TayRMS =0.D0
      TayRRMS=0.D0
      TayVRMS=0.D0
      !--- Integrals changed to not represent the different volumes
      !    represented by each cylinder on Nov. 2 2007:
      do nS=1,nSmax
         ! Old form used again for VgRMS for make is comparable with
         ! VpRMS below.
         VgRMS=VgRMS + 2.D0*pi*h(nS)*sZ(nS)*dsZ * VpIntN(nS)*VpIntN(nS)
         TayRMS =TayRMS +dsZ*dabs(TayIntN(nS))
         TayRRMS=TayRRMS+dsZ*dabs(TayRIntN(nS))
         TayVRMS=TayVRMS+dsZ*dabs(TayVIntN(nS))
         if ( nS <= nSI ) then
            VgRMS=VgRMS + 2.D0*pi*h(nS)*sZ(nS)*dsZ * VpIntS(nS)*VpIntS(nS)
            TayRMS =TayRMS +dsZ*dabs(TayIntS(nS))
            TayRRMS=TayRRMS+dsZ*dabs(TayRIntS(nS))
            TayVRMS=TayVRMS+dsZ*dabs(TayVIntS(nS))
         end if
      end do

      !--- s-derivatives:
      !------ Create arrays to be differentiated:
      do nS=1,nSmax
         Os2(nS)=1.D0/(sZ(nS)*sZ(nS))
         SVpIntN(nS) =VpIntN(nS)/sZ(nS)
         SBspIntN(nS)=h(nS)*sZ(nS)*sZ(nS)*BspIntN(nS)
         SBs2IntN(nS)=h(nS)*sZ(nS)**3*Bs2IntN(nS)
         if ( sZ(nS) < r_ICB ) then ! inside TC
            SVpIntS(nS) =VpIntS(nS)/sZ(nS)
            SBspIntS(nS)=h(nS)*sZ(nS)*sZ(nS)*BspIntS(nS)
            SBs2IntS(nS)=h(nS)*sZ(nS)**3*Bs2IntS(nS)
         end if
      end do
      !------ Create Derivatives:
      f1=1.D0/(12.D0*dsZ)
      f2=f1/dsZ
      do nS=3,nSmax-2
         dSVpIntN =f1*(     SVpIntN(nS-2)-8.D0*SVpIntN(nS-1) +        &
              &                   8.D0*SVpIntN(nS+1)-     SVpIntN(nS+2) )
         d2SVpIntN=f2*(    -SVpIntN(nS-2)+16.D0*SVpIntN(nS-1) -       &
              &       30.D0*SVpIntN(nS)+16.D0*SVpIntN(nS+1)-SVpIntN(nS+2))
         dSBspIntN=f1*(     SBspIntN(nS-2)-8.D0*SBspIntN(nS-1) +      &
              &                   8.D0*SBspIntN(nS+1)-     SBspIntN(nS+2) )
         dSBs2IntN=f1*(     SBs2IntN(nS-2)-8.D0*SBs2IntN(nS-1) +      &
              &                   8.D0*SBs2IntN(nS+1)-     SBs2IntN(nS+2) )
         TauN(nS)  =Oh(nS)*(Os2(nS)*dSBspIntN+TauBN(nS))
         dTTauN(nS)=sZ(nS)*d2SVpIntN*Bs2IntN(nS) +                    &
              &               Oh(nS)*(Os2(nS)*dSVpIntN*dSBs2IntN +    &
              &                        sZ(nS)*dSVpIntN*dTTauBN(nS) )
      end do
      TauN(1)        =TauN(3)
      TauN(2)        =TauN(3)
      TauN(nSmax-1)  =TauN(nSmax-2)
      TauN(nSmax)    =TauN(nSmax-2)
      dTTauN(1)      =dTTauN(3)
      dTTauN(2)      =dTTauN(3)
      dTTauN(nSmax-1)=dTTauN(nSmax-2)
      dTTauN(nSmax)  =dTTauN(nSmax-2)

      !------ South hemisphere:
      do nS=1,nSmax
         if ( sZ(nS) < r_ICB .and. nS > 2 ) then
            dSVpIntS =f1*(SVpIntS(nS-2)-8.D0*SVpIntS(nS-1) +          &
                 &                          8.D0*SVpIntS(nS+1)-SVpIntS(nS+2) )
            d2SVpIntS=f2*(-SVpIntS(nS-2)+16.D0*SVpIntS(nS-1) -        &
                 &          30.D0*SVpIntS(nS)+16.D0*SVpIntS(nS+1)-SVpIntS(nS+2))
            dSBspIntS=f1*(SBspIntS(nS-2)-8.D0*SBspIntS(nS-1) +        &
                 &                      8.D0*SBspIntS(nS+1)-SBspIntS(nS+2) )
            dSBs2IntS=f1*(SBs2IntS(nS-2)-8.D0*SBs2IntS(nS-1) +        &
                 &                      8.D0*SBs2IntS(nS+1)-SBs2IntS(nS+2) )
            TauS(nS) =Oh(nS)*(Os2(nS)*dSBspIntS+TauBS(nS))
            dTTauS(nS)=sZ(nS)*d2SVpIntS*Bs2IntS(nS) +                 &
                 &               Oh(nS)*(Os2(nS)*dSVpIntS*dSBs2IntS + &
                 &                        sZ(nS)*dSVpIntS*dTTauBS(nS) )
         else
            TauS(nS)  =TauN(nS)
            dTTauS(nS)=dTTauN(nS)
         end if
      end do

      !--- Output of z-integral:
      open(nOutFile, file=TOfileNhs, status='unknown',    &
           &       form='unformatted', position='append')
      IF( nTOsets == 1 ) then

         write(message,'(" ! TO: No. of s-values:",i4)') int(nSmax/r_cmb)
         call logWrite(message)
         write(message,'(" ! TO: No. of z-values:",i4)') nNorm
         call logWrite(message)

         write(nOutFile) float(nSmax)                     ! 1
         write(nOutFile) (sngl(sZ(nS)),nS=1,nSmax)        ! 2
      end if
      write(nOutFile)  sngl(time),                                 &! 3
           &          (sngl(VpIntN(nS))             ,nS=1,nSmax),  &! 4
           &          (sngl(dVpIntN(nS))            ,nS=1,nSmax),  &! 5
           &          (sngl(ddVpIntN(nS))           ,nS=1,nSmax),  &! 6
           &          (sngl(VpRIntN(nS))            ,nS=1,nSmax),  &! 7
           &          (sngl(RstrIntN(nS))           ,nS=1,nSmax),  &! 8
           &          (sngl(AstrIntN(nS))           ,nS=1,nSmax),  &! 9
           &          (sngl(LFfac*LFIntN(nS))       ,nS=1,nSmax),  &! 10
           &          (sngl(StrIntN(nS))            ,nS=1,nSmax),  &! 11
           &          (sngl(TayIntN(nS))            ,nS=1,nSmax),  &! 12
           &          (sngl(LFfac*TauN(nS))         ,nS=1,nSmax),  &! 13
           &          (sngl(LFfac*TauBN(nS)/h(nS))  ,nS=1,nSmax),  &! 14
           &          (sngl(LFfac*BspdIntN(nS))     ,nS=1,nSmax),  &! 15 For first part of dTau
           &          (sngl(LFfac*BpsdIntN(nS))     ,nS=1,nSmax),  &! 16 For second part of dTau
           &          (sngl(LFfac*dTauBN(nS)/h(nS)) ,nS=1,nSmax),  &! 17 Boundary contribution !
           &          (sngl(LFfac*dTTauN(nS))       ,nS=1,nSmax),  &! 18
           &          (sngl(LFfac*dTTauBN(nS)/h(nS)),nS=1,nSmax),  &! 19
           &          (sngl(LFfac*Bs2IntN(nS))      ,nS=1,nSmax)       ! 20 
      close(nOutFile)

      open(nOutFile, file=TOfileShs, status='unknown',       &
           &       form='unformatted', position='append')
      IF( nTOsets == 1 ) then
         write(nOutFile) float(nSmax)
         write(nOutFile) (sngl(sZ(nS)),nS=1,nSmax)
      end if
      write(nOutFile)  sngl(time),                                 &
           &          (sngl(VpIntS(nS))             ,nS=1,nSmax),  &
           &          (sngl(dVpIntS(nS))            ,nS=1,nSmax),  &
           &          (sngl(ddVpIntS(nS))           ,nS=1,nSmax),  &
           &          (sngl(VpRIntS(nS))            ,nS=1,nSmax),  &
           &          (sngl(RstrIntS(nS))           ,nS=1,nSmax),  &
           &          (sngl(AstrIntS(nS))           ,nS=1,nSmax),  &
           &          (sngl(LFfac*LFIntS(nS))       ,nS=1,nSmax),  &
           &          (sngl(StrIntS(nS))            ,nS=1,nSmax),  &
           &          (sngl(TayIntS(nS))            ,nS=1,nSmax),  &
           &          (sngl(LFfac*TauS(nS))         ,nS=1,nSmax),  &
           &          (sngl(LFfac*TauBS(nS)/h(nS))  ,nS=1,nSmax),  &
           &          (sngl(LFfac*BspdIntS(nS))     ,nS=1,nSmax),  &! For first part of dTau
           &          (sngl(LFfac*BpsdIntS(nS))     ,nS=1,nSmax),  &! For second part of dTau
           &          (sngl(LFfac*dTauBS(nS)/h(nS)) ,nS=1,nSmax),  &! Boundary contribution !
           &          (sngl(LFfac*dTTauS(nS))       ,nS=1,nSmax),  &
           &          (sngl(LFfac*dTTauBS(nS)/h(nS)),nS=1,nSmax),  &
           &          (sngl(LFfac*Bs2IntS(nS))      ,nS=1,nSmax)
      close(nOutFile)

      ! Note: integral changed on Nov. 2 2007 (see above)
      fac=1.D0/sZ(nSmax) ! I integrate (roughly) from s=0 to s=r_cmb.
      VgRMS  =vSF*dsqrt(VgRMS/vol_oc)
      TayRMS =fac*TayRMS
      TayRRMS=fac*TayRRMS
      TayVRMS=fac*TayVRMS

      lTOrms=.true.
      if ( lTOmov .or. lTOrms ) then

         !--- Output of graphic file:
         !--- Transform back to radial space:
         do nR=1,n_r_max
            do l=1,l_max
               lm=lo_map%lm2(l,0)
               !              dzRstrLMr(lm,nR)=dzRstrLMr(lm,nR) +
               !     +                         dzAstrLMr(lm,nR)
               dzVpLMr(l+1,nR)=orho1(nR)*real(z(lm,nR)) ! instead of transform copy again
            end do
         end do

         call costf1(dzdVpLMr,lmMaxS,1,lmMaxS,workA,i_costf_init,d_costf_init)
         call costf1(dzRstrLMr,lmMaxS,1,lmMaxS,workA,i_costf_init,d_costf_init)
         call costf1(dzAstrLMr,lmMaxS,1,lmMaxS,workA,i_costf_init,d_costf_init)
         call costf1(dzStrLMr,lmMaxS,1,lmMaxS,workA,i_costf_init,d_costf_init)
         call costf1(dzLFLMr,lmMaxS,1,lmMaxS,workA,i_costf_init,d_costf_init)
         call costf1(dzCorLMr,lmMaxS,1,lmMaxS,workA,i_costf_init,d_costf_init)

         !--- Open output file
         nFields=7
         if ( lTOmov ) then

            open(nOutFile, file=movFile, status='unknown',              &
                 &             form='unformatted', position='append')

            !--- Write header into output file:
            if ( nTOmovSets == 0 ) then
               version='JW_Movie_Version_2'
               write(nOutFile) version
               dumm(1)=102           ! type of input
               dumm(2)=3             ! marker for constant phi plane
               dumm(3)=0.D0          ! surface constant
               dumm(4)=nFields       ! no of fields 
               write(nOutFile) (real(dumm(n),4),n=1,4)

               dumm(1)=11.0          ! Field marker for AS vPhi 
               dumm(2)=61.0          ! Field marker for Reynolds Force 
               dumm(3)=62.0          ! Field marker for Advective Force  
               dumm(4)=63.0          ! Field marker for Viscous Force 
               dumm(5)=64.0          ! Field marker for Lorentz Force
               dumm(6)=65.0          ! Field marker for Coriolis force
               dumm(7)=66.0          ! Field marker for dtVp
               write(nOutFile) (sngl(dumm(n)),n=1,nFields)

               !------ Now other info about grid and parameters:
               write(nOutFile) runid        ! run identifier 
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
               write(nOutFile) (sngl(dumm(n)),     n=1,11)  
               write(nOutFile) (sngl(r(n)/r_CMB),  n=1,n_r_max)
               write(nOutFile) (sngl(theta_ord(n)),n=1,n_theta_max)
               write(nOutFile) (sngl(phi(n)),      n=1,n_phi_max)

            end if ! Write Header ?

            nTOmovSets=nTOmovSets+1

            dumm(1)=nTOmovSets        ! time frame number for movie
            dumm(2)=time              ! time      
            dumm(3)=0.D0      
            dumm(4)=0.D0     
            dumm(5)=0.D0          
            dumm(6)=0.D0          
            dumm(7)=0.D0         
            dumm(8)=0.D0        
            write(nOutFile) (sngl(dumm(n)),n=1,8)

         end if ! Produce movie output ?

         if ( lTOrms ) then

            do nR=1,n_r_max
               VpR(nR)   =0.D0
               dVpR(nR)  =0.D0
               RstrR(nR) =0.D0
               AstrR(nR) =0.D0
               LFR(nR)   =0.D0
               LFABSR(nR)=0.D0
               StrR(nR)  =0.D0
               CorR(nR)  =0.D0
            end do

            nTOrmsSets=nTOrmsSets+1
            open(nOutFile2, file=tayFile, form='unformatted',           &
                 &             status='unknown', position='append')
            if ( nTOrmsSets == 1 ) then
               write(nOutFile2) float(n_r_max)
               write(nOutFile2) (sngl(r(nR)),nR=1,n_r_max)
            end if

         end if

         do nOut=1,nFields ! Loop over four output fields

            do nR=1,n_r_max ! Loop over radial points
               rS=r(nR)

               do n=1,nThetaBs ! Loop over theta blocks
                  nThetaStart=(n-1)*sizeThetaB+1

                  !------ Convert from lm to theta block: 
                  if ( nOut == 1 ) then
                     call get_PAS(dzVpLMr(1,nR),outBlock,rS,nThetaStart,sizeThetaB)
                  else if ( nOut == 2 ) then
                     call get_PAS(dzRstrLMr(1,nR),outBlock,rS,nThetaStart,sizeThetaB)
                  else if ( nOut == 3 ) then
                     call get_PAS(dzAstrLMr(1,nR),outBlock,rS,nThetaStart,sizeThetaB)
                  else if ( nOut == 4 ) then
                     call get_PAS(dzStrLMr(1,nR),outBlock,rS,nThetaStart,sizeThetaB)
                  else if ( nOut == 5 ) then
                     call get_PAS(dzLFLMr(1,nR),outBlock,rS,nThetaStart,sizeThetaB)
                  else if ( nOut == 6 ) then
                     call get_PAS(dzCorLMr(1,nR),outBlock,rS,nThetaStart,sizeThetaB)
                  else if ( nOut == 7 ) then
                     call get_PAS(dzdVpLMr(1,nR),outBlock,rS,nThetaStart,sizeThetaB)
                  end if

                  !------ Storage of field in fout for theta block,
                  !       integration and Max/Min values
                  nTheta=(n-1)*sizeThetaB
                  do nThetaBlock=1,sizeThetaB
                     nTheta=nTheta+1
                     nThetaNHS=(nTheta+1)/2
                     sS       =rS*sinTheta(nTheta)
                     !--------- Convert to correct order in theta grid points:
                     if ( mod(nTheta,2) == 1 ) then
                        nThetaOrd=(nTheta+1)/2
                     else
                        nThetaOrd=n_theta_max-nTheta/2+1
                     end if
                     nPos=(nR-1)*n_theta_max+nThetaOrd

                     if ( nOut == 1 ) then
                        !--------------- Zonal flow:
                        fOut(nPos)=vSF*outBlock(nThetaBlock)
                        VpR(nR)=VpR(nR) + gauss(nThetaNHS)*fOut(nPos)/sS
                     else if ( nOut == 2 ) then
                        !--------------- Reynolds force:
                        fOut(nPos)=fSF*outBlock(nThetaBlock)
                        RstrR(nR)=RstrR(nR) + gauss(nThetaNHS)*fOut(nPos)/sS
                     else if ( nOut == 3 ) then
                        !--------------- Advective force:
                        fOut(nPos)=fSF*outBlock(nThetaBlock)
                        AstrR(nR)=AstrR(nR) + gauss(nThetaNHS)*fOut(nPos)/sS
                     else if ( nOut == 4 ) then
                        !--------------- Viscous force:
                        fOut(nPos)=fSF*outBlock(nThetaBlock)
                        StrR(nR)=StrR(nR) + gauss(nThetaNHS)*fOut(nPos)/sS           
                     else if ( nOut == 5 ) then
                        !--------------- Lorentz force:
                        fOut(nPos)=fSF*LFfac*outBlock(nThetaBlock)
                        LFR(nR)=LFR(nR) + gauss(nThetaNHS)*fOut(nPos)/sS
                        LFABSR(nR)=LFABSR(nR) + gauss(nThetaNHS)*dabs(fOut(nPos))/sS
                     else if ( nOut == 6 ) then
                        !--------------- Corriolis force:
                        fOut(nPos)=fSF*outBlock(nThetaBlock)
                        CorR(nR)=CorR(nR) + gauss(nThetaNHS)*fOut(nPos)/sS
                     else if ( nOut == 7 ) then
                        !--------------- dtVp:
                        fOut(nPos)=vSF*outBlock(nThetaBlock)
                        dVpR(nR)  =dVpR(nR) + gauss(nThetaNHS)*fOut(nPos)/sS
                     end if

                  end do ! Loop over thetas in block
               end do    ! Loop over theta blocks

            end do ! Loop over R

            !------ Output of stress contributions for one radial grid point:
            !       Write all fields into movie style file
            if ( lTOmov ) write(nOutFile) (sngl(fOut(nPos)),nPos=1,nFieldSize)

         end do ! Loop over output functions

         if ( l_save_out ) then
            open(n_log_file, file=log_file, status='unknown', position='append')
         end if
         if ( lTOmov ) then 
            close(nOutFile)
            write(n_log_file,'(1p,/,A,I8,A,i8,D16.6)')                    &
                 &           ' ! WRITING TO MOVIE FRAME NO ',nTOmovSets,  &
                 &           ' AT STEP/TIME :',n_time_step,time*tScale  
         end if

         !--- Realtive importance of azimuthal and geostrophic flow
         VRMS   =dsqrt(2.D0*eKin/vol_oc)
         VpRMS  =dsqrt(2.D0*eKinTAS/vol_oc)

         if ( VRMS /= 0d0) then
            VpRMS  =VpRMS/VRMS
            VgRMS  =VgRMS/VRMS
         end if

         do nR=1,n_r_max
            fac=1.D0/2.D0 ! Cancel factor 2 from theta integral
            VpR(nR)    =fac*VpR(nR)
            dVpR(nR)   =fac*dVpR(nR)
            RstrR(nR)  =fac*RstrR(nR)
            AstrR(nR)  =fac*AstrR(nR)
            StrR(nR)   =fac*StrR(nR)
            LFR(nR)    =fac*LFR(nR)
            CorR(nR)   =fac*CorR(nR)
            if ( LFABSR(nR) /= 0d0 ) then
               TayR(nR)   =dabs(LFR(nR))/(fac*LFABSR(nR))
            else
               TayR(nR)   =0d0
            end if
            !              TayRMSR(nR)=rS*rS*TayR(nR)
            TayRMSR(nR)=TayR(nR)
         end do

         !--- Now perform the radial integral: ( not tested )
         TaySRMS=rInt_R(TayRMSR,n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
         !--- And finally calculate the mean value, the factor 4*pi comes from
         !    the fact that the surface integral has already been cared for
         !    NOTE: Integral for RMS Taylorisation changed to not respect the 
         !    different volumes represented by each shell on 2 Nov 2007.
         !          fac    =4.D0*pi/vol_oc
         fac    =1.D0/(r_cmb-r_icb) ! =1
         TaySRMS=fac*TaySRMS

         write(nOutFile2)                                             &
              &          sngl(time),sngl(VpRMS**2),sngl(VgRMS**2),    &
              &          sngl(TayRMS),sngl(TaySRMS),sngl(TayRRMS),    &
              &          sngl(TayVRMS),sngl(eKin),                    &! 3
              &          (sngl(VpR(nR))  ,nR=1,n_r_max),              &! 4
              &          (sngl(dVpR(nR)) ,nR=1,n_r_max),              &! 5
              &          (sngl(RstrR(nR)),nR=1,n_r_max),              &! 6
              &          (sngl(AstrR(nR)),nR=1,n_r_max),              &! 7
              &          (sngl(LFR(nR))  ,nR=1,n_r_max),              &! 8
              &          (sngl(StrR(nR)) ,nR=1,n_r_max),              &! 9
              &          (sngl(CorR(nR)) ,nR=1,n_r_max),              &! 10
              &          (sngl(TayR(nR)) ,nR=1,n_r_max)  ! 11
         close(nOutFile2)

      end if

      timeLast=time
      if ( lVerbose ) write(*,*) '! End of outTO!'

      lTOZwrite=.false.

   end subroutine outTO
!----------------------------------------------------------------------------
end module outTO_mod
