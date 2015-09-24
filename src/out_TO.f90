module outTO_mod

   use parallel_mod, only: rank
   use precision_mod
   use truncation, only: n_r_max, n_r_maxStr, n_theta_maxStr, l_max, &
                         n_theta_max, n_phi_max, minc, lStressMem,   &
                         lm_max
   use radial_functions, only: r_ICB, i_costf_init, d_costf_init, r, &
                               r_CMB, orho1, drx
   use physical_parameters, only: ra, ek, pr, prmag, radratio, LFfac
   use torsional_oscillations, only: V2AS, Bs2AS, BspAS, BszAS, BpzAS, &
                                     BspdAS, BpsdAS, BzpdAS, BpzdAS,   &
                                     dzdVpLMr, dzddVpLMr, dzRstrLMr,   &
                                     dzAstrLMr, dzStrLMr, dzLFLMr,     &
                                     dzCorLMr, TO_gather_Rloc_on_rank0, &
                                     dzRstrLMr_Rloc
   use num_param, only: tScale
   use blocking, only: nThetaBs, sizeThetaB, nfs, st_map
   use horizontal_data, only: phi, sinTheta, theta_ord, gauss
   use logic, only: lVerbose, l_save_out
   use output_data, only: sDens, zDens, TAG, log_file, runid, n_log_file
   use const, only: pi, vol_oc, one, two, half, four
   use LMLoop_data, only: llm, ulm
   use charmanip, only: dble2str
   use integration, only: rInt_R
   use plms_theta, only: plm_theta
   use TO_helpers, only: getPAStr, get_PAS, getAStr
   use useful, only: logWrite
   use legendre_grid_to_spec, only: legTFAS, legTFAS2
   use chebInt_mod, only: chebInt, chebIntInit
   use cosine_transform, only: costf1
   use communications, only: gather_all_from_lo_to_rank0,gt_OC

   implicit none 

   private
   
   integer :: lmMaxS
   integer :: nZmaxA,nZmaxL
   integer :: nSmaxA,nSmaxL
 
   !-- Plms: Plm,sin
   real(cp), allocatable :: PlmS(:,:,:)
   real(cp), allocatable :: dPlmS(:,:,:)
   real(cp), allocatable :: OsinTS(:,:)
   real(outp), allocatable :: VpM(:,:), LFM(:,:), dVpM(:,:), AstrM(:,:)
   real(outp), allocatable :: RstrM(:,:), CorM(:,:), StrM(:,:), CLM(:,:)
   real(cp), allocatable :: d_costf_initZ(:,:), zZ(:,:), rZ(:,:)
   integer, allocatable :: i_costf_initZ(:,:), nZmaxS(:)

   public :: initialize_outTO_mod, outTO

contains

   subroutine initialize_outTO_mod

      nZmaxL=lStressMem*722
      nZmaxA=max(2,nZmaxL)
      nSmaxL=lStressMem*625
      nSmaxA=max(3,nSmaxL)
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
      !
      !   Output of axisymmetric zonal flow, its relative strength,
      !   its time variation, and all forces acting on it.
      !   The slowest part in the TO process is the repetitious calculation
      !   of plms by subroutine plm_theta. They are needed in getAStr and 
      !   getPAStr when I transform on the cylindrical grid. 
      !   The necessary plms could simply be calculated one and then 
      !   be stored for later use!
      !

      !-- Input of variables:
      real(cp),         intent(in) :: time
      integer,          intent(in) :: n_time_step
      real(cp),         intent(in) :: eKin, eKinTAS
      integer,          intent(in) :: nOutFile, nOutFile2
      character(len=*), intent(in) :: TOfileNhs, TOfileShs, movFile, tayFile
      logical,          intent(in) :: lTOmov
      complex(cp),      intent(in) :: z(llm:ulm,n_r_max)
      real(cp),         intent(in) :: omega_ic, omega_ma
      integer,          intent(inout) :: nTOsets, nTOmovSets, nTOrmsSets
      logical,          intent(inout) :: lTOrms, lTOZwrite

      !-- Output field:
      real(cp) :: fOut(n_theta_maxStr*n_r_maxStr)

      !-- Local variables:
      logical :: lTC,lStopRun

      !-- (l,r) Representation of the different contributions
      real(cp) :: dzVpLMr(lmMaxS,n_r_maxStr)
      real(cp) :: V2LMr(lmMaxS,n_r_maxStr)
      real(cp) :: Bs2LMr(lmMaxS,n_r_maxStr)
      real(cp) :: BszLMr(lmMaxS,n_r_maxStr)
      real(cp) :: BspLMr(lmMaxS,n_r_maxStr)
      real(cp) :: BpzLMr(lmMaxS,n_r_maxStr)
      real(cp) :: BspdLMr(lmMaxS,n_r_maxStr)
      real(cp) :: BpsdLMr(lmMaxS,n_r_maxStr)
      real(cp) :: BzpdLMr(lmMaxS,n_r_maxStr)
      real(cp) :: BpzdLMr(lmMaxS,n_r_maxStr)
      complex(cp) :: zS(lm_max,n_r_maxStr)

      !---- Work array:
      real(cp) :: workA(lmMaxS,n_r_maxStr)

      integer :: lm,l ! counter for degree

      integer :: nSmax,nS,nSI
      real(cp) :: zNorm  ! Norm z interval
      integer :: nNorm  ! No. of grid points for norm interval
      real(cp) :: zMin,zMax!,help ! integration boundarie, help variable
      logical :: lAS    ! .true. if axisymmetric (m=0) functions
      real(cp) :: sZ(nSmaxA),dsZ ! cylindrical radius s and s-step
      real(cp) :: h(nSmaxA),Oh(nSmaxA)
      real(cp) :: Os2(nSmaxA)

      integer :: nR     ! counter for radial grid point
      integer :: n      ! counter for theta blocks
      integer :: nOut,nFields ! counter for output fields
      integer :: nTheta ! counter for all thetas
      integer :: nThetaBlock ! counter for thetas in block
      integer :: nThetaOrd ! counter for ordered thetas
      integer :: nThetaNHS
      integer :: nThetaStart

      integer :: nZ,nZmax,nZmaxNS,nZmaxH!,nZP
      real(cp) :: VpS(nZmaxA)      
      real(cp) :: dVpS(nZmaxA)      
      real(cp) :: ddVpS(nZmaxA)      
      real(cp) :: V2S(nZmaxA)
      real(cp) :: LFS(nZmaxA)
      real(cp) :: CorS(nZmaxA)
      real(cp) :: RstrS(nZmaxA)
      real(cp) :: AstrS(nZmaxA)
      real(cp) :: StrS(nZmaxA)
      real(cp) :: Bs2S(nZmaxA)
      real(cp) :: BspS(nZmaxA)
      real(cp) :: BspdS(nZmaxA)
      real(cp) :: BpsdS(nZmaxA)
      real(cp) :: TayS(nZmaxA)
      real(cp) :: TayRS(nZmaxA)
      real(cp) :: TayVS(nZmaxA)

      real(cp) :: VpIntN(nSmaxA)  ,VpIntS(nSmaxA)    ! integration results
      real(cp) :: dVpIntN(nSmaxA) ,dVpIntS(nSmaxA)   ! integration results
      real(cp) :: ddVpIntN(nSmaxA),ddVpIntS(nSmaxA)  ! integration results
      real(cp) :: VpRIntN(nSmaxA) ,VpRIntS(nSmaxA)   ! for different s and 
      real(cp) :: V2IntS(nSmaxA)  ,V2IntN(nSmaxA)
      real(cp) :: LFIntN(nSmaxA)  ,LFIntS(nSmaxA)   
      real(cp) :: RstrIntN(nSmaxA),RstrIntS(nSmaxA) 
      real(cp) :: AstrIntN(nSmaxA),AstrIntS(nSmaxA) 
      real(cp) :: StrIntN(nSmaxA) ,StrIntS(nSmaxA) 
      real(cp) :: Bs2IntN(nSmaxA) ,Bs2IntS(nSmaxA)
      real(cp) :: BspIntN(nSmaxA) ,BspIntS(nSmaxA)
      real(cp) :: BspdIntN(nSmaxA),BspdIntS(nSmaxA)
      real(cp) :: BpsdIntN(nSmaxA),BpsdIntS(nSmaxA)
      real(cp) :: TayIntN(nSmaxA) ,TayIntS(nSmaxA) 
      real(cp) :: TayRIntN(nSmaxA),TayRIntS(nSmaxA) 
      real(cp) :: TayVIntN(nSmaxA),TayVIntS(nSmaxA) 
      real(cp) :: SVpIntN(nSmaxA) ,SVpIntS(nSmaxA)   ! help arrays and values for 
      real(cp) :: SBs2IntN(nSmaxA),SBs2IntS(nSmaxA)  ! differentiation in s   
      real(cp) :: SBspIntN(nSmaxA),SBspIntS(nSmaxA)
      real(cp) :: dSVpIntN, dSVpIntS
      real(cp) :: d2SVpIntN,d2SVpIntS
      real(cp) :: dSBspIntN,dSBspIntS
      real(cp) :: dSBs2IntN,dSBs2IntS
      real(cp) :: TauN(nSmaxA),TauS(nSmaxA)          ! Taylor integral
      real(cp) :: TauBN(nSmaxA),TauBS(nSmaxA)       
      real(cp) :: dTauBN(nSmaxA),dTauBS(nSmaxA)    
      real(cp) :: dTTauN(nSmaxA),dTTauS(nSmaxA)      ! time change of Tau...
      real(cp) :: dTTauBN(nSmaxA),dTTauBS(nSmaxA)   

      !-- For integration along z:
      real(cp) :: zALL(2*nZmaxA)
      real(cp) :: thetaZ
      real(cp) :: fac
      real(cp) :: vSF,fSF,f1,f2

      !-- For boundaries:
      real(cp) :: rB(2),BspB(2),BspdB(2),BpsdB(2)
      real(cp) :: Bs2B(2),BszB(2),BpzB(2),BzpdB(2),BpzdB(2)

      !-- For sphere integration:
      real(cp) :: CorR(n_r_max)
      real(cp) :: RstrR(n_r_max)
      real(cp) :: AstrR(n_r_max)
      real(cp) :: StrR(n_r_max)
      real(cp) :: VpR(n_r_max),dVpR(n_r_max)
      real(cp) :: LFR(n_r_max),LFABSR(n_r_max)
      real(cp) :: TayR(n_r_max),TayRMSR(n_r_max)
      real(cp) :: TayRMS,TaySRMS
      real(cp) :: TayRRMS,TayVRMS
      real(cp) :: VpRMS,VRMS,VgRMS
      real(cp) :: rS,sS
      real(cp) :: outBlock(nfs)
      real(cp), save :: timeLast!,tNorm
      real(outp), save :: timeAve
      real(outp) :: dt

      character(len=255) :: message
      character(len=64) :: version,fileName
      integer :: nFieldSize,nPos

      !integer :: nLines
      real(cp) :: dumm(12)

      !-- Huge arrays for time average ....
      logical :: l_TOZave

      !-- For TOZ output files:
      integer, save :: nTOZfile!,length
      character(len=10) :: string

      if ( lVerbose ) write(*,*) '! Starting outTO!'

      call TO_gather_Rloc_on_rank0

      call gather_all_from_lo_to_rank0(gt_OC,z,zS)

      if ( rank == 0 ) then

         nTOsets=nTOsets+1

         l_TOZave=.true.

         !--- Rescaling for rotation time scale and planetary radius
         !    length scale, for the velocity I use the Rossby number
         !       vSF=ek/r_CMB                   
         !       fSF=ek*ek/(four*pi**2*r_CMB)  
         vSF=one
         fSF=one

         nFieldSize=n_theta_maxStr*n_r_maxStr

         !-- Start with calculating advection due to axisymmetric flows:

         zNorm=one               ! This is r_CMB-r_ICB
         nNorm=int(zDens*n_r_max) ! Covered with nNorm  points !
         nSmax=n_r_max+int(r_ICB*real(n_r_max,cp))
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
               call legTFAS(V2LMr(1,nR),V2AS(nThetaStart,nR),            &
                    &               l_max+1,nThetaStart,sizeThetaB)
               call legTFAS2(Bs2LMr(1,nR),BszLMr(1,nR),                            &
                    &               Bs2AS(nThetaStart,nR),BszAS(nThetaStart,nR),   &
                    &               l_max+1,nThetaStart,sizeThetaB)
               call legTFAS2(BspLMr(1,nR),BpzLMr(1,nR),                            &
                    &               BspAS(nThetaStart,nR),BpzAS(nThetaStart,nR),   &
                    &               l_max+1,nThetaStart,sizeThetaB)
               call legTFAS2(BspdLMr(1,nR),BpsdLMr(1,nR),                          &
                    &               BspdAS(nThetaStart,nR),BpsdAS(nThetaStart,nR), &
                    &               l_max+1,nThetaStart,sizeThetaB)
               call legTFAS2(BzpdLMr(1,nR),BpzdLMr(1,nR),                          &
                    &               BzpdAS(nThetaStart,nR),BpzdAS(nThetaStart,nR), &
                    &               l_max+1,nThetaStart,sizeThetaB)
            end do
         end do ! Loop over radial grid points

         do nR=1,n_r_max
            do l=1,l_max
               lm=st_map%lm2(l,0)
               dzVpLMr(l+1,nR)=real(zS(lm,nR))
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

         dsZ   =r_CMB/real(nSmax,cp)  ! Step in s controlled by nSmax
         nSI   =0
         do nS=1,nSmax
            sZ(nS)=(nS-half)*dsZ
            if ( sZ(nS) < r_ICB .and. nS > nSI ) nSI=nS
         end do

         if ( nTOsets == 1 ) nTOZfile=0
         if ( lTOZwrite ) then
            nTOZfile=nTOZfile+1
            call dble2str(real(nTOZfile,cp),string)
            fileName='TOZ_'//trim(adjustl(string))//'.'//TAG
            open(95, file=fileName, form='unformatted', status='unknown')
            write(95) real(time,kind=outp), real(nSmax,kind=outp), &
                    & real(omega_ic,kind=outp), real(omega_ma,kind=outp)
            write(95) (real(sZ(nS),kind=outp),nS=1,nSmax)
         end if
         if ( nTOsets > 1 .and. l_TOZave ) then
            fileName='TOZM.'//TAG
            open(96,file=fileName, form='unformatted', status='unknown')
            write(96) real(nSmax,kind=outp), real(omega_ic,kind=outp), &
                    & real(omega_ma,kind=outp)
            write(96) (real(sZ(nS),kind=outp),nS=1,nSmax)
         end if

         lStopRun=.false.
         do nS=1,nSmax

            !------ Get integral boundaries for this s:
            zMax=sqrt(r_CMB*r_CMB-sZ(nS)*sZ(nS))
            if ( sZ(nS) < r_ICB ) then
               lTC=.true.
               zMin=sqrt(r_ICB*r_ICB-sZ(nS)*sZ(nS))
            else
               lTC=.false.
               zMin=-zMax
            end if
            h(nS) =zMax-zMin
            Oh(nS)=one/h(nS)

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
                  rZ(nZ,nS)    =sqrt(zZ(nZ,nS)**2+sZ(nS)**2)
                  thetaZ       =atan2(sZ(nS),zZ(nZ,nS))
                  OsinTS(nZ,nS)=one/sin(thetaZ)
                  call plm_theta(thetaZ,l_max,0,minc,                    &
                       &         PlmS(1,nZ,nS),dPlmS(1,nZ,nS),lmMaxS,2)
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
               TayS(nZ) =abs(LFS(nZ))
               TayRS(nZ)=abs(RstrS(nZ))
               TayVS(nZ)=abs(StrS(nZ))
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
                  timeAve=1.e0_outp
                  do nZ=1,nZmaxNS
                     VpM(nZ,nS)  =real(VpS(nZ),kind=outp)
                     dVpM(nZ,nS) =real(dVpS(nZ),kind=outp)
                     LFM(nZ,nS)  =real(LFfac*LFS(nZ),kind=outp)
                     RstrM(nZ,nS)=real(RstrS(nZ),kind=outp)
                     AstrM(nZ,nS)=real(AstrS(nZ),kind=outp)
                     StrM(nZ,nS) =real(StrS(nZ),kind=outp)
                     CorM(nZ,nS) =real(CorS(nZ),kind=outp)
                     CLM(nZ,nS)  =real(CorS(nZ)+LFfac*LFS(nZ),kind=outp)
                  end do
               else if ( nTOsets == 2 ) then
                  dt=real(time-timeLast,kind=outp)
                  timeAve=dt
                  do nZ=1,nZmaxNS
                     VpM(nZ,nS)  =dt*(VpM(nZ,nS)  +real(VpS(nZ),kind=outp))
                     dVpM(nZ,nS) =dt*(dVpM(nZ,nS) +real(dVpS(nZ),kind=outp))
                     LFM(nZ,nS)  =dt*(LFM(nZ,nS)  +real(LFfac*LFS(nZ),kind=outp))
                     RstrM(nZ,nS)=dt*(RstrM(nZ,nS)+real(RstrS(nZ),kind=outp))
                     AstrM(nZ,nS)=dt*(AstrM(nZ,nS)+real(AstrS(nZ),kind=outp))
                     StrM(nZ,nS) =dt*(StrM(nZ,nS) +real(StrS(nZ),kind=outp))
                     CorM(nZ,nS) =dt*(CorM(nZ,nS) +real(CorS(nZ),kind=outp))
                     CLM(nZ,nS)  =dt*(CLM(nZ,nS)  +real(CorS(nZ)+LFfac*LFS(nZ),kind=outp))
                  end do
               else
                  dt=real(time-timeLast,kind=outp)
                  timeAve=timeAve+dt
                  do nZ=1,nZmaxNS
                     VpM(nZ,nS)  =VpM(nZ,nS)  +dt*real(VpS(nZ),kind=outp)
                     dVpM(nZ,nS) =dVpM(nZ,nS) +dt*real(dVpS(nZ),kind=outp)
                     LFM(nZ,nS)  =LFM(nZ,nS)  +dt*real(LFfac*LFS(nZ),kind=outp)
                     RstrM(nZ,nS)=RstrM(nZ,nS)+dt*real(RstrS(nZ),kind=outp)
                     AstrM(nZ,nS)=AstrM(nZ,nS)+dt*real(AstrS(nZ),kind=outp)
                     StrM(nZ,nS) =StrM(nZ,nS) +dt*real(StrS(nZ),kind=outp)
                     CorM(nZ,nS) =CorM(nZ,nS) +dt*real(CorS(nZ),kind=outp)
                     CLM(nZ,nS)  =CLM(nZ,nS)  +dt*real(CorS(nZ)+LFfac*LFS(nZ),kind=outp)
                  end do
               end if
            end if

            if ( l_TOZave .and. nTOsets > 1 ) then
               write(96) real(nZmaxNS,kind=outp)
               write(96) (real(zALL(nZ),kind=outp),nZ=1,nZmaxNS), &
                    &    (VpM(nZ,nS)/timeAve  ,nZ=1,nZmaxNS),     &
                    &    (dVpM(nZ,nS)/timeAve ,nZ=1,nZmaxNS),     &
                    &    (RstrM(nZ,nS)/timeAve,nZ=1,nZmaxNS),     &
                    &    (AstrM(nZ,nS)/timeAve,nZ=1,nZmaxNS),     &
                    &    (LFM(nZ,nS)/timeAve  ,nZ=1,nZmaxNS),     &
                    &    (StrM(nZ,nS)/timeAve ,nZ=1,nZmaxNS),     &
                    &    (CorM(nZ,nS)/timeAve ,nZ=1,nZmaxNS),     &
                    &    (CLM(nZ,nS)/timeAve  ,nZ=1,nZmaxNS)
            end if
            if ( lTOZwrite ) then
               write(95) real(nZmaxNS)
               write(95) (real(zALL(nZ),kind=outp) ,nZ=1,nZmaxNS),      &
                    &    (real(VpS(nZ),kind=outp)  ,nZ=1,nZmaxNS),      &
                    &    (real(dVpS(nZ),kind=outp) ,nZ=1,nZmaxNS),      &
                    &    (real(RstrS(nZ),kind=outp),nZ=1,nZmaxNS),      &
                    &    (real(AstrS(nZ),kind=outp),nZ=1,nZmaxNS),      &
                    &    (real(LFfac*LFS(nZ),kind=outp),nZ=1,nZmaxNS),  &
                    &    (real(StrS(nZ),kind=outp) ,nZ=1,nZmaxNS),      &
                    &    (real(CorS(nZ),kind=outp) ,nZ=1,nZmaxNS)
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

            if ( V2IntN(nS) < 0.0_cp ) then
               VpRIntN(nS)=one
            else
               VpRIntN(nS)=abs(VpIntN(nS))/sqrt(V2IntN(nS))
            end if
            VpRIntN(nS) =min(one,VpRIntN(nS))
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
               if ( V2IntS(nS) < 0.0_cp ) then
                  VpRIntS(nS)=one
               else
                  VpRIntS(nS)=abs(VpIntS(nS))/sqrt(V2IntS(nS))
               end if
               VpRIntS(nS) =min(one,VpRIntS(nS))
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
               zMax = sqrt(r_CMB**2-sZ(nS)**2)
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
               zMax = sqrt(r_ICB**2-sZ(nS)**2)
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
               zMax = sqrt(r_CMB**2-sZ(nS)**2)
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
         VgRMS  =0.0_cp
         TayRMS =0.0_cp
         TayRRMS=0.0_cp
         TayVRMS=0.0_cp
         !--- Integrals changed to not represent the different volumes
         !    represented by each cylinder on Nov. 2 2007:
         do nS=1,nSmax
            ! Old form used again for VgRMS for make is comparable with
            ! VpRMS below.
            VgRMS=VgRMS + two*pi*h(nS)*sZ(nS)*dsZ * VpIntN(nS)*VpIntN(nS)
            TayRMS =TayRMS +dsZ*abs(TayIntN(nS))
            TayRRMS=TayRRMS+dsZ*abs(TayRIntN(nS))
            TayVRMS=TayVRMS+dsZ*abs(TayVIntN(nS))
            if ( nS <= nSI ) then
               VgRMS=VgRMS + two*pi*h(nS)*sZ(nS)*dsZ * VpIntS(nS)*VpIntS(nS)
               TayRMS =TayRMS +dsZ*abs(TayIntS(nS))
               TayRRMS=TayRRMS+dsZ*abs(TayRIntS(nS))
               TayVRMS=TayVRMS+dsZ*abs(TayVIntS(nS))
            end if
         end do

         !--- s-derivatives:
         !------ Create arrays to be differentiated:
         do nS=1,nSmax
            Os2(nS)=one/(sZ(nS)*sZ(nS))
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
         f1=one/(12.0_cp*dsZ)
         f2=f1/dsZ
         do nS=3,nSmax-2
            dSVpIntN =f1*(     SVpIntN(nS-2)-8.0_cp*SVpIntN(nS-1) +        &
                 &                   8.0_cp*SVpIntN(nS+1)-     SVpIntN(nS+2) )
            d2SVpIntN=f2*(    -SVpIntN(nS-2)+16.0_cp*SVpIntN(nS-1) -       &
                 &       30.0_cp*SVpIntN(nS)+16.0_cp*SVpIntN(nS+1)-SVpIntN(nS+2))
            dSBspIntN=f1*(     SBspIntN(nS-2)-8.0_cp*SBspIntN(nS-1) +      &
                 &                   8.0_cp*SBspIntN(nS+1)-     SBspIntN(nS+2) )
            dSBs2IntN=f1*(     SBs2IntN(nS-2)-8.0_cp*SBs2IntN(nS-1) +      &
                 &                   8.0_cp*SBs2IntN(nS+1)-     SBs2IntN(nS+2) )
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
               dSVpIntS =f1*(SVpIntS(nS-2)-8.0_cp*SVpIntS(nS-1) +          &
                    &                          8.0_cp*SVpIntS(nS+1)-SVpIntS(nS+2) )
               d2SVpIntS=f2*(-SVpIntS(nS-2)+16.0_cp*SVpIntS(nS-1) -        &
                    &          30.0_cp*SVpIntS(nS)+16.0_cp*SVpIntS(nS+1)-SVpIntS(nS+2))
               dSBspIntS=f1*(SBspIntS(nS-2)-8.0_cp*SBspIntS(nS-1) +        &
                    &                      8.0_cp*SBspIntS(nS+1)-SBspIntS(nS+2) )
               dSBs2IntS=f1*(SBs2IntS(nS-2)-8.0_cp*SBs2IntS(nS-1) +        &
                    &                      8.0_cp*SBs2IntS(nS+1)-SBs2IntS(nS+2) )
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
         if ( nTOsets == 1 ) then

            write(message,'(" ! TO: No. of s-values:",i4)') int(nSmax/r_cmb)
            call logWrite(message)
            write(message,'(" ! TO: No. of z-values:",i4)') nNorm
            call logWrite(message)

            write(nOutFile) real(nSmax,kind=outp)                        ! 1
            write(nOutFile) (real(sZ(nS),kind=outp),nS=1,nSmax)          ! 2
         end if
         write(nOutFile)  real(time,kind=outp),                         &! 3
              &  (real(VpIntN(nS),kind=outp)             ,nS=1,nSmax),  &! 4
              &  (real(dVpIntN(nS),kind=outp)            ,nS=1,nSmax),  &! 5
              &  (real(ddVpIntN(nS),kind=outp)           ,nS=1,nSmax),  &! 6
              &  (real(VpRIntN(nS),kind=outp)            ,nS=1,nSmax),  &! 7
              &  (real(RstrIntN(nS),kind=outp)           ,nS=1,nSmax),  &! 8
              &  (real(AstrIntN(nS),kind=outp)           ,nS=1,nSmax),  &! 9
              &  (real(LFfac*LFIntN(nS),kind=outp)       ,nS=1,nSmax),  &! 10
              &  (real(StrIntN(nS),kind=outp)            ,nS=1,nSmax),  &! 11
              &  (real(TayIntN(nS),kind=outp)            ,nS=1,nSmax),  &! 12
              &  (real(LFfac*TauN(nS),kind=outp)         ,nS=1,nSmax),  &! 13
              &  (real(LFfac*TauBN(nS)/h(nS),kind=outp)  ,nS=1,nSmax),  &! 14
              &  (real(LFfac*BspdIntN(nS),kind=outp)     ,nS=1,nSmax),  &! 15 For first part of dTau
              &  (real(LFfac*BpsdIntN(nS),kind=outp)     ,nS=1,nSmax),  &! 16 For second part of dTau
              &  (real(LFfac*dTauBN(nS)/h(nS),kind=outp) ,nS=1,nSmax),  &! 17 Boundary contribution !
              &  (real(LFfac*dTTauN(nS),kind=outp)       ,nS=1,nSmax),  &! 18
              &  (real(LFfac*dTTauBN(nS)/h(nS),kind=outp),nS=1,nSmax),  &! 19
              &  (real(LFfac*Bs2IntN(nS),kind=outp)      ,nS=1,nSmax)    ! 20 
         close(nOutFile)

         open(nOutFile, file=TOfileShs, status='unknown',       &
              &       form='unformatted', position='append')
         if ( nTOsets == 1 ) then
            write(nOutFile) real(nSmax,kind=outp)
            write(nOutFile) (real(sZ(nS),kind=outp),nS=1,nSmax)
         end if
         write(nOutFile)  real(time,kind=outp),                         &
              &  (real(VpIntS(nS),kind=outp)             ,nS=1,nSmax),  &
              &  (real(dVpIntS(nS),kind=outp)            ,nS=1,nSmax),  &
              &  (real(ddVpIntS(nS),kind=outp)           ,nS=1,nSmax),  &
              &  (real(VpRIntS(nS),kind=outp)            ,nS=1,nSmax),  &
              &  (real(RstrIntS(nS),kind=outp)           ,nS=1,nSmax),  &
              &  (real(AstrIntS(nS),kind=outp)           ,nS=1,nSmax),  &
              &  (real(LFfac*LFIntS(nS),kind=outp)       ,nS=1,nSmax),  &
              &  (real(StrIntS(nS),kind=outp)            ,nS=1,nSmax),  &
              &  (real(TayIntS(nS),kind=outp)            ,nS=1,nSmax),  &
              &  (real(LFfac*TauS(nS),kind=outp)         ,nS=1,nSmax),  &
              &  (real(LFfac*TauBS(nS)/h(nS),kind=outp)  ,nS=1,nSmax),  &
              &  (real(LFfac*BspdIntS(nS),kind=outp)     ,nS=1,nSmax),  &
              &  (real(LFfac*BpsdIntS(nS),kind=outp)     ,nS=1,nSmax),  &
              &  (real(LFfac*dTauBS(nS)/h(nS),kind=outp) ,nS=1,nSmax),  &
              &  (real(LFfac*dTTauS(nS),kind=outp)       ,nS=1,nSmax),  &
              &  (real(LFfac*dTTauBS(nS)/h(nS),kind=outp),nS=1,nSmax),  &
              &  (real(LFfac*Bs2IntS(nS),kind=outp)      ,nS=1,nSmax)
         close(nOutFile)

         ! Note: integral changed on Nov. 2 2007 (see above)
         fac=one/sZ(nSmax) ! I integrate (roughly) from s=0 to s=r_cmb.
         VgRMS  =vSF*sqrt(VgRMS/vol_oc)
         TayRMS =fac*TayRMS
         TayRRMS=fac*TayRRMS
         TayVRMS=fac*TayVRMS

         lTOrms=.true.
         if ( lTOmov .or. lTOrms ) then

            !--- Output of graphic file:
            !--- Transform back to radial space:
            do nR=1,n_r_max
               do l=1,l_max
                  lm=st_map%lm2(l,0)
                  dzVpLMr(l+1,nR)=orho1(nR)*real(zS(lm,nR)) ! instead of transform copy again
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
                  dumm(3)=0.0_cp          ! surface constant
                  dumm(4)=nFields       ! no of fields 
                  write(nOutFile) (real(dumm(n),kind=outp),n=1,4)

                  dumm(1)=11.0          ! Field marker for AS vPhi 
                  dumm(2)=61.0          ! Field marker for Reynolds Force 
                  dumm(3)=62.0          ! Field marker for Advective Force  
                  dumm(4)=63.0          ! Field marker for Viscous Force 
                  dumm(5)=64.0          ! Field marker for Lorentz Force
                  dumm(6)=65.0          ! Field marker for Coriolis force
                  dumm(7)=66.0          ! Field marker for dtVp
                  write(nOutFile) (real(dumm(n),kind=outp),n=1,nFields)

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
                  write(nOutFile) (real(dumm(n),kind=outp),     n=1,11)  
                  write(nOutFile) (real(r(n)/r_CMB,kind=outp),  n=1,n_r_max)
                  write(nOutFile) (real(theta_ord(n),kind=outp),n=1,n_theta_max)
                  write(nOutFile) (real(phi(n),kind=outp),      n=1,n_phi_max)

               end if ! Write Header ?

               nTOmovSets=nTOmovSets+1

               dumm(1)=nTOmovSets        ! time frame number for movie
               dumm(2)=time              ! time      
               dumm(3)=0.0_cp      
               dumm(4)=0.0_cp     
               dumm(5)=0.0_cp          
               dumm(6)=0.0_cp          
               dumm(7)=0.0_cp         
               dumm(8)=0.0_cp        
               write(nOutFile) (real(dumm(n),kind=outp),n=1,8)

            end if ! Produce movie output ?

            if ( lTOrms ) then

               do nR=1,n_r_max
                  VpR(nR)   =0.0_cp
                  dVpR(nR)  =0.0_cp
                  RstrR(nR) =0.0_cp
                  AstrR(nR) =0.0_cp
                  LFR(nR)   =0.0_cp
                  LFABSR(nR)=0.0_cp
                  StrR(nR)  =0.0_cp
                  CorR(nR)  =0.0_cp
               end do

               nTOrmsSets=nTOrmsSets+1
               open(nOutFile2, file=tayFile, form='unformatted',           &
                    &             status='unknown', position='append')
               if ( nTOrmsSets == 1 ) then
                  write(nOutFile2) real(n_r_max,kind=outp)
                  write(nOutFile2) (real(r(nR),kind=outp),nR=1,n_r_max)
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
                           LFABSR(nR)=LFABSR(nR) + gauss(nThetaNHS)*abs(fOut(nPos))/sS
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
               if ( lTOmov ) write(nOutFile) (real(fOut(nPos),kind=outp),nPos=1,nFieldSize)

            end do ! Loop over output functions

            if ( l_save_out ) then
               open(n_log_file, file=log_file, status='unknown', position='append')
            end if
            if ( lTOmov ) then 
               close(nOutFile)
               call logWrite(' ')
               write(message,'(1p,A,I8,A,ES16.6,I8)')              &
                    & "! WRITING TO MOVIE FRAME NO ",nTOmovSets,   &
                    & "       at time/step",time*tScale,n_time_step
               call logWrite(message)
            end if

            !--- Realtive importance of azimuthal and geostrophic flow
            VRMS   =sqrt(two*eKin/vol_oc)
            VpRMS  =sqrt(two*eKinTAS/vol_oc)

            if ( VRMS /= 0.0_cp) then
               VpRMS  =VpRMS/VRMS
               VgRMS  =VgRMS/VRMS
            end if

            do nR=1,n_r_max
               fac=half ! Cancel factor 2 from theta integral
               VpR(nR)    =fac*VpR(nR)
               dVpR(nR)   =fac*dVpR(nR)
               RstrR(nR)  =fac*RstrR(nR)
               AstrR(nR)  =fac*AstrR(nR)
               StrR(nR)   =fac*StrR(nR)
               LFR(nR)    =fac*LFR(nR)
               CorR(nR)   =fac*CorR(nR)
               if ( LFABSR(nR) /= 0.0_cp ) then
                  TayR(nR)   =abs(LFR(nR))/(fac*LFABSR(nR))
               else
                  TayR(nR)   =0.0_cp
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
            !          fac    =four*pi/vol_oc
            fac    =one/(r_cmb-r_icb) ! =1
            TaySRMS=fac*TaySRMS

            write(nOutFile2)                                                       &
                 &          real(time,kind=outp),real(VpRMS**2,kind=outp),         &
                 &          real(VgRMS**2,kind=outp),real(TayRMS,kind=outp),       &
                 &          real(TaySRMS,kind=outp),real(TayRRMS,kind=outp),       &
                 &          real(TayVRMS,kind=outp),real(eKin,kind=outp),          &! 3
                 &          (real(VpR(nR),kind=outp)  ,nR=1,n_r_max),              &! 4
                 &          (real(dVpR(nR),kind=outp) ,nR=1,n_r_max),              &! 5
                 &          (real(RstrR(nR),kind=outp),nR=1,n_r_max),              &! 6
                 &          (real(AstrR(nR),kind=outp),nR=1,n_r_max),              &! 7
                 &          (real(LFR(nR),kind=outp)  ,nR=1,n_r_max),              &! 8
                 &          (real(StrR(nR),kind=outp) ,nR=1,n_r_max),              &! 9
                 &          (real(CorR(nR),kind=outp) ,nR=1,n_r_max),              &! 10
                 &          (real(TayR(nR),kind=outp) ,nR=1,n_r_max)                ! 11
            close(nOutFile2)

         end if

         timeLast=time

         lTOZwrite=.false.

      end if ! Rank 0

      if ( lVerbose ) write(*,*) '! End of outTO!'

   end subroutine outTO
!----------------------------------------------------------------------------
end module outTO_mod
