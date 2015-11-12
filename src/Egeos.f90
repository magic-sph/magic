module Egeos_mod
 
   use precision_mod
   use truncation, only: n_r_max, lm_max, n_m_max, n_phi_max, nrpGeos, &
                         n_r_maxGeos, lm_maxGeos, minc, l_max, m_max
   use parallel_mod, only: rank
   use radial_functions, only: cheb_norm, r_ICB, r_CMB, chebt_oc
   use physical_parameters, only: ra, ek, pr, prmag, g0, g1, g2, &
                                  radratio, polind, strat
   use num_param, only: tScale
   use blocking, only: lm2l, lm2m, lm2mc
   use horizontal_data, only: dLh, phi, dPhi
   use logic, only: lVerbose, l_corrMov, l_anel
   use output_data, only: sDens, zDens, tag, runid, nSmaxA, nZmaxA
   use constants, only: pi, zero, ci, one, two, half
   use LMLoop_data, only: llm,ulm
   use communications, only: gather_all_from_lo_to_rank0,gt_OC
   use plms_theta, only: plm_theta
   use fft, only: fft_to_real
   use cosine_transform_odd
   use chebInt_mod

   implicit none 
 
   private
 
   real(cp), allocatable :: PlmS(:,:,:)  ! This is huge !
   real(cp), allocatable :: dPlmS(:,:,:) ! This is huge !
   real(cp), allocatable :: OsinTS(:,:)

   type(costf_odd_t), allocatable :: chebt_Z(:)
   integer, allocatable :: nZmaxS(:)
   real(cp), allocatable :: zZ(:,:), rZ(:,:)
   real(cp), parameter :: eps = 10.0_cp*epsilon(one)
 
   public :: initialize_Egeos_mod, getEgeos

contains

   subroutine initialize_Egeos_mod

      allocate( OsinTS(nZmaxA/2+1,nSmaxA) )
      allocate( PlmS(lm_maxGeos,nZmaxA/2+1,nSmaxA) )  ! This is huge !
      allocate( dPlmS(lm_maxGeos,nZmaxA/2+1,nSmaxA) ) ! This is huge !
      allocate( chebt_Z(nSmaxA) )
      allocate( nZmaxS(nSmaxA) )
      allocate( zZ(nZmaxA,nSmaxA) )
      allocate( rZ(nZmaxA,nSmaxA) )

   end subroutine initialize_Egeos_mod
!----------------------------------------------------------------------------
   subroutine getEgeos(time,nGeosSets,w,dw,ddw,z,dz,         &
        &              Egeos,EkNTC,EkSTC,Ekin,               &
        &              dpFlow,dzFlow,CVzOTC,CVorOTC,CHelOTC)
      !
      !   Output of axisymmetric zonal flow, its relative strength,
      !   its time variation, and all forces acting on it.
      !   The slowest part in the TO process is the repitions calculation
      !   of Plms by subroutine plm_theta. They are needed in getDVptr 
      !   when I transform on the cylindrical grid. 
      !   The necessary plms could simply be calculated one and then 
      !   be stored for later use! See s_outTOnew.f.
      !

      !-- Input of variables:
      real(cp),    intent(in) :: time
      integer,     intent(in) :: nGeosSets 
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dw(llm:ulm,n_r_max)
      complex(cp), intent(in) :: ddw(llm:ulm,n_r_max)
      complex(cp), intent(in) :: z(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dz(llm:ulm,n_r_max)

      !-- Output variables:
      real(cp), intent(out) :: Egeos,EkNTC,EkSTC,Ekin
      real(cp), intent(out) :: dpFlow,dzFlow ! RMS length scales
      real(cp), intent(out) :: CVzOTC,CVorOTC,CHelOTC

      !-- Local variables:
      logical :: lDeriv
      integer :: nSmax,nS,nS_ICB
      real(cp) :: zNorm          ! Norm z interval
      integer :: nNorm           ! No. of grid points for norm interval
      real(cp) :: zMin,zMax,help ! integration boundarie, help variable
      logical :: lAS             ! .true. if axisymmetric (m=0) functions
      real(cp) :: sZ(nSmaxA),dsZ ! cylindrical radius s and s-step
      integer :: nPhi,nI
      real(cp) :: phiNorm
      logical :: lTC

      !-- Local field copies to avoid changes by back and forth transform:
      complex(cp) :: wS(lm_maxGeos,n_r_maxGeos)
      complex(cp) :: dwS(lm_maxGeos,n_r_maxGeos)
      complex(cp) :: ddwS(lm_maxGeos,n_r_maxGeos)
      complex(cp) :: zS(lm_maxGeos,n_r_maxGeos)
      complex(cp) :: dzS(lm_maxGeos,n_r_maxGeos)
      !---- Additional work array
      complex(cp) :: workA(lm_maxGeos,n_r_maxGeos)

      !-- Representation in (phi,z):
      real(cp) :: VrS(nrpGeos,nZmaxA),VrInt(nZmaxA),VrIntS
      real(cp) :: VtS(nrpGeos,nZmaxA),VtInt(nZmaxA),VtIntS
      real(cp) :: VpS(nrpGeos,nZmaxA),VpInt(nZmaxA),VpIntS
      real(cp) :: VozS(nrpGeos,nZmaxA)
      real(cp) :: sinT,cosT
      integer :: nInt,nInts   ! index for NHS and SHS integral
      integer :: nZ,nZmax,nZS,nZN
      integer :: lm,nR
      real(cp) :: EkInt(nZmaxA),EkIntS
      real(cp) :: EkSTC_s(nSmaxA),EkNTC_s(nSmaxA),EkOTC_s(nSmaxA)
      real(cp) :: Egeos_s(nSmaxA)
      real(cp) :: EkOTC
      real(cp) :: dpEkInt(nZmaxA),dpEkIntS
      real(cp) :: dzEkInt(nZmaxA),dzEkIntS
      real(cp) :: dpEk_s(nSmaxA),dzEk_s(nSmaxA)
      real(cp) :: thetaZ

      logical :: lStopRun

      !-- Correlation (new Oct. 4 2007)
      logical :: lCorrel
      real(cp) :: VzS,VzN,VorS,VorN,surf,delz
      real(cp) :: VzSN,VzSS,VzNN,VorSN,VorSS,VorNN,HelZZ,VZZ,VorZZ
      real(cp) :: CVz_I,CVor_I,CHel_I
      real(cp) :: CVz_s(nSmaxA),CVor_s(nSmaxA),CHel_S(nSmaxA)

      !-- Movie output
      integer :: nOutFile,n
      character(len=66) :: version,movFile
      integer :: nFields,nFieldSize
      real(cp) :: dumm(40)
      real(outp) :: CVz(nrpGeos,nZmaxA)
      real(outp) :: CVor(nrpGeos,nZmaxA)
      real(outp) :: CHel(nrpGeos,nZmaxA)


      if ( lVerbose ) write(*,*) '! Starting outGeos!'

      call gather_all_from_lo_to_rank0(gt_OC,w,wS)
      call gather_all_from_lo_to_rank0(gt_OC,dw,dwS)
      call gather_all_from_lo_to_rank0(gt_OC,ddw,ddwS)
      call gather_all_from_lo_to_rank0(gt_OC,z,zS)
      call gather_all_from_lo_to_rank0(gt_OC,dz,dzS)

      if ( rank == 0 ) then
         lCorrel=.true. ! Calculate Vz and Vorz north/south correlation

         phiNorm=two*pi/n_phi_max
         lStopRun=.false.
         lDeriv=.true.

         !---- Get resolution in s and z for z-integral:
         zNorm=one               ! This is r_CMB-r_ICB
         nNorm=int(zDens*n_r_max) ! Covered with nNorm  points !
         nSmax=n_r_max+int(r_ICB*real(n_r_max,cp)) 
         nSmax=int(sDens*nSmax)
         dsZ  =r_CMB/real(nSmax,cp)  ! Step in s controlled by nSmax
         if ( nSmax > nSmaxA ) then
            write(*,*) 'Increase nSmaxA in getGeos!'
            write(*,*) 'Should be at least nSmax=',nSmax
            stop
         end if
         lAS=.false.

         do nS=1,nSmaxA
            EkSTC_s(nS)=0.0_cp
            EkNTC_s(nS)=0.0_cp
            EkOTC_s(nS)=0.0_cp
            Egeos_s(nS)=0.0_cp
            dpEk_s(nS) =0.0_cp
            dzEk_s(nS) =0.0_cp
            CVz_s(nS)  =0.0_cp
            CVor_s(nS) =0.0_cp
            CHel_s(nS) =0.0_cp
         end do

         !---- Copy for following costf so that no      
         !     back transform is needed that would change
         !     the field (just a little bit) anyway...
         do nR=1,n_r_max
            do lm=1,lm_max
               wS(lm,nR)=wS(lm,nR)*dLh(lm)
            end do
         end do

         call chebt_oc%costf1(wS,lm_max,1,lm_max,workA)
         call chebt_oc%costf1(dwS,lm_max,1,lm_max,workA)
         call chebt_oc%costf1(ddwS,lm_max,1,lm_max,workA)
         call chebt_oc%costf1(zS,lm_max,1,lm_max,workA)
         call chebt_oc%costf1(dzS,lm_max,1,lm_max,workA)

         !---- Contributions are now in fully spectral space!
         !---- Do the z-integral:
         nI=0

         do nS=1,nSmax
            sZ(nS)=(nS-half)*dsZ
            if ( sZ(nS) < r_ICB ) then
               lTC=.true.
            else
               lTC=.false.
            end if
            if ( nS > 1 ) then
               if ( sZ(nS-1) < r_ICB .and. sZ(ns) >= r_ICB ) nS_ICB=nS
            end if

            !------ Get integral boundaries for this s:
            zMin=-sqrt(r_CMB*r_CMB-sZ(nS)*sZ(nS))
            if ( lTC ) then
               zMax=-sqrt(r_ICB*r_ICB-sZ(nS)*sZ(nS))
            else
               zMax=-zMin
            end if

            if ( nGeosSets == 1 ) then
               !------ Initialize integration for NHS:
               !       Each processor calculates Cheb transform data 
               !       for HIS nS and the Plms along the Cylinder
               !       chebIntInit returns zZ,nZmaxS,chebt_Z
               call chebIntInit(zMin,zMax,zNorm,nNorm,                   &
                    &           nZmaxA,zZ(1,nS),nZmaxS(nS),chebt_Z(nS))
               !------ Calculate and store 1/sin(theta) and Plms,dPlms for 
               !       southern HS:

               if ( lTC ) then 
                  nZmax=nZmaxS(nS)  ! nZmax point in each polar region
               else
                  nZmax=(nZmaxS(nS)-1)/2+1 ! Odd point number !    
                  ! all together nZmaxS(nS) from
                  ! south to north including equator
               end if
               if ( 2*nZmax > nZmaxA ) then ! TC case most critical
                  write(*,*) '! nZmaxA too small in getEgeos!'
                  write(*,*) '! Should be at least:',2*nZmax
                  lStopRun=.true.
               end if
               do nZ=1,nZmax
                  rZ(nZ,nS)    =sqrt(zZ(nZ,nS)**2+sZ(nS)**2)
                  thetaZ       =atan2(sZ(nS),zZ(nZ,nS))
                  OsinTS(nZ,nS)=one/sin(thetaZ)
                  call plm_theta(thetaZ,l_max,m_max,minc,                &
                       &            PlmS(1,nZ,nS),dPlmS(1,nZ,nS),lm_max,2)
               end do
            end if

            !--------- Get the flow components for all northern and
            !          southern thetas and all phis:
            if ( lTC ) then
               nZmax=2*nZmaxS(nS) ! north and south points
               ! calculated in one go
            else
               nZmax=nZmaxS(nS)
            end if

            call getDVptr(wS,dwS,ddwS,zS,dzS,r_ICB,r_CMB,rZ(1,nS),            &
                 &      nZmax,nZmaxA,PlmS(1,1,nS),dPlmS(1,1,nS),OsinTS(1,nS), &
                 &                           lDeriv,VrS,VtS,VpS,VozS,dpEkInt)

            nZmax=nZmaxS(nS)

            !------- Perform z-integral(s) for all phis:
            do nPhi=1,n_phi_max

               !------- Two seperate integrals inside TC:
               if ( lTC ) then
                  nInts=2 ! separate north and south integral
               else
                  nInts=1
               end if
               do nInt=1,nInts
                  if ( nInt == 1 ) then
                     do nZ=1,nZmax ! Copy on simpler array
                        VrInt(nZ)=VrS(nPhi,nZ)
                        VtInt(nZ)=VtS(nPhi,nZ)
                        VpInt(nZ)=VpS(nPhi,nZ)
                        EkInt(nZ)=(VrInt(nZ)**2+VtInt(nZ)**2+VpInt(nZ)**2)
                     end do
                  else if ( nInt == 2 ) then
                     help=zMax
                     zMax=-zMin
                     zMin=-help
                     do nZ=1,nZmax
                        VrInt(nZ)=VrS(nPhi,nZ+nZmax)
                        VtInt(nZ)=VtS(nPhi,nZ+nZmax)
                        VpInt(nZ)=VpS(nPhi,nZ+nZmax)
                        EkInt(nZ)=(VrInt(nZ)**2+VtInt(nZ)**2+VpInt(nZ)**2)
                     end do
                  end if

                  !-------- NOTE: chebIntD replaces VrInt with z-derivative
                  !               for lDeriv=.true.
                  VrIntS=chebIntD(VrInt,lDeriv,zMin,zMax,nZmax,nZmaxA,chebt_Z(nS))
                  VtIntS=chebIntD(VtInt,lDeriv,zMin,zMax,nZmax,nZmaxA,chebt_Z(nS))
                  VpIntS=chebIntD(VpInt,lDeriv,zMin,zMax,nZmax,nZmaxA,chebt_Z(nS))
                  EkIntS=chebIntD(EkInt,.false.,zMin,zMax,nZmax,nZmaxA,chebt_Z(nS))

                  !-------- Get volume integral of energies:
                  Egeos_s(nS)=Egeos_s(nS) +                              &
                       &      half*phiNorm*(zMax-zMin)*sZ(nS)*dsZ *     &
                       &      (VrIntS**2+VtIntS**2+VpIntS**2)
                  if ( lTC ) then
                     if ( nInt == 1 ) then
                        EkSTC_s(nS)=EkSTC_s(nS) +                        &
                             &      half*phiNorm*(zMax-zMin)*sZ(nS)*dsZ*EkIntS
                     else if ( nInt == 2 ) then
                        EkNTC_s(nS)=EkNTC_s(nS) +                        &
                             &      half*phiNorm*(zMax-zMin)*sZ(nS)*dsZ*EkIntS
                     end if
                  else
                     EkOTC_s(nS)=EkOTC_s(nS) +                           &
                          &      half*phiNorm*(zMax-zMin)*sZ(nS)*dsZ*EkIntS
                  end if

                  !-------- Note: chebIntD returns the z derivative
                  if ( lDeriv ) then
                     do nZ=1,nZmax
                        dzEkInt(nZ)=VrInt(nZ)*VrInt(nZ) +                &
                             &      VtInt(nZ)*VtInt(nZ) + VpInt(nZ)*VpInt(nZ)
                     end do
                     dzEkIntS=chebInt(dzEkInt,zMin,zMax,nZmax,nZmaxA,chebt_Z(nS))
                     dzEk_s(nS)=dzEk_s(nS) +                             &
                          &     half*phiNorm*(zMax-zMin)*sZ(nS)*dsZ*dzEkIntS
                  end if

               end do ! Loop over north/south integral


               ! --- All the stuff for North/South correlation only outside TC
               !     and only 0.1 away from boundaries:
               !     Here only integrals over one hemisphere are required. I thus
               !     copy the southern points to the northern hermisphere.
               if ( lCorrel .and. .not.lTC ) then
                  VzSN =0.0_cp
                  VzSS =0.0_cp
                  VzNN =0.0_cp
                  VorSN=0.0_cp
                  VorSS=0.0_cp
                  VorNN=0.0_cp
                  HelZZ=0.0_cp
                  VZZ  =0.0_cp
                  VorZZ=0.0_cp
                  do nZN=1,nZmax/2 ! Dont use equatorial point
                     nZS  =nZmax-nZN+1
                     delz =zZ(nZN,nS)-zZ(nZN+1,nS)
                     sinT =one/OsinTS(nZN,nS)
                     cosT =sqrt(one-sinT**2)
                     VzN  = cosT*VrS(nPhi,nZN)-sinT*VtS(nPhi,nZN)
                     VzS  =-cosT*VrS(nPhi,nZS)-sinT*VtS(nPhi,nZS)
                     VorN =VozS(nPhi,nZN)
                     VorS =VozS(nPhi,nZS)
                     VzSN =VzSN+delz*VzS*VzN
                     VzNN =VzNN+delz*VzN*VzN
                     VzSS =VzSS+delz*VzS*VzS
                     VorSN=VorSN+delz*VorS*VorN
                     VorNN=VorNN+delz*VorN*VorN
                     VorSS=VorSS+delz*VorS*VorS
                     HelZZ=HelZZ+delz*(VorN*VzN+VorS*VzS)**2
                     VZZ  =VZZ  +delz*(VzN*VzN+VzS*VzS)
                     VorZZ=VorZZ+delz*(VorN*VorN+VorS*VorS)
                  end do

                  !------------ Calculate the integrals from zMin+0.1 to zMax-0.1:
                  if ( VzSS /= 0.0_cp .and. VzNN /= 0.0_cp ) then
                     CVz_I = VzSN/sqrt(VzSS*VzNN)
                  else
                     CVz_I = 0.0_cp
                  end if
                  if ( VorSS /= 0.0_cp .and. VorNN /= 0.0_cp ) then
                     CVor_I=VorSN/sqrt(VorSS*VorNN)
                  else
                     CVor_I=0.0_cp
                  end if
                  if ( VZZ /= 0.0_cp .and. VorZZ /= 0.0_cp .and. HelZZ > eps ) then
                     CHel_I=sqrt(HelZZ/(VZZ*VorZZ))
                  else
                     CHel_I=0.0_cp
                  end if
                  CVz(nPhi,nSmax-nS+1) =real(CVz_I,kind=outp)
                  CVor(nPhi,nSmax-nS+1)=real(CVor_I,kind=outp)
                  CHel(nPhi,nSmax-nS+1)=real(CHel_I,kind=outp)
                  CVz_s(nS) =CVz_s(nS) +phiNorm*sZ(nS)*dsZ*CVz_I  
                  CVor_s(nS)=CVor_s(nS)+phiNorm*sZ(nS)*dsZ*CVor_I
                  CHel_s(nS)=CHel_s(nS)+phiNorm*sZ(nS)*dsZ*CHel_I

               end if ! lCorrel ?

            end do  ! Loop over phi

            if ( lDeriv ) then
               !--------- dpEkInt treated differently cause phi intergral has 
               !          already been preformed by getDVptr
               if ( lTC ) then
                  nInts=2 ! separate north and south integral
               else
                  nInts=1
               end if
               do nInt=1,nInts
                  if ( nInt == 2 ) then
                     help=zMax
                     zMax=-zMin
                     zMin=-help
                     do nZ=1,nZmax
                        dpEkInt(nZ)=dpEkInt(nZ+nZmax)
                     end do
                  end if
                  dpEkIntS=chebInt(dpEkInt,zMin,zMax,nZmax,nZmaxA,chebt_Z(nS))
                  dpEk_s(nS)=dpEk_s(nS) +                                &
                       &     half*(zMax-zMin)*sZ(nS)*dsZ*dpEkIntS /     &
                       &     (sZ(nS)**2) ! Convert angle to length
               end do
            end if

         end do ! big loop over S

         !--- Collect:
         EkSTC  =0.0_cp
         EkNTC  =0.0_cp
         EkOTC  =0.0_cp
         Egeos  =0.0_cp
         CVzOTC =0.0_cp
         CVorOTC=0.0_cp
         CHelOTC=0.0_cp
         do nS=1,nSmax
            Egeos=Egeos+Egeos_s(nS)
            if ( sZ(ns) < r_ICB ) then
               EkSTC=EkSTC+EkSTC_s(nS)
               EkNTC=EkNTC+EkNTC_s(nS)
            else
               EkOTC  =EkOTC  +EkOTC_s(nS)
               if ( lCorrel ) then
                  CVzOTC =CVzOTC +CVz_s(nS)
                  CVorOTC=CVorOTC+CVor_s(nS)
                  CHelOTC=CHelOTC+CHel_s(nS)
               end if
            end if
         end do
         if ( lCorrel ) then
            surf=0.0_cp
            do nS=1,nSmax
               if ( sZ(nS) >= r_ICB ) surf=surf+sZ(nS)*dsZ
            end do
            surf   =two*pi*surf
            CVzOTC =CVzOTC/surf
            CVorOTC=CVorOTC/surf
            CHelOTC=CHelOTC/surf
         end if
         Ekin=EkSTC+EkNTC+EKOTC ! Can be used for testing 

         dpFlow=0.0_cp
         dzFlow=0.0_cp
         if ( lDeriv ) then
            do nS=1,nSmax
               dpFlow=dpFlow+dpEk_s(nS)
               dzFlow=dzFlow+dzEk_s(nS)
            end do
            if ( dpFlow /= 0.0_cp ) then
               dpFlow = sqrt(Ekin/dpFlow)
            else
               dpFlow = 0.0_cp
            end if
            if ( dzFlow /= 0.0_cp ) then
               dzFlow = sqrt(Ekin/dzFlow)
            else
               dzFlow = 0.0_cp
            end if
         end if

         !        write(99,*) 'E:',EkSTC,EkNTC,EkOTC
         !        write(99,*) 'Ekin:',Ekin
         !        write(99,*) 'Egeos:',Egeos

         !--- Write correlation movie:
         if ( l_corrMov ) then

            !--- Determine s used for correl
            n=0
            do nS=1,nSmax
               if ( sZ(nS) >= r_ICB+0.1_cp .and. sZ(nS) <= r_CMB-0.1_cp ) then
                  n=n+1
               else 
                  do nPhi=1,n_phi_max
                     CVz(nPhi,nSmax-nS+1)=0.0_cp
                     CVor(nPhi,nSmax-nS+1)=0.0_cp
                  end do
               end if
            end do
            !           ndS=nSmax-n
            !          nSmax=n

            nOutFile=93
            movFile ='CVorz_mov.'//tag
            open(nOutFile, file=movFile, status='unknown',   &
            &    form='unformatted', position='append')

            !--- Write header into output file:
            if ( nGeosSets == 1 ) then

               nFields=3
               nFieldSize=(nSmax-nS_ICB+1)*n_phi_max
               version='JW_Movie_Version_2'
               write(nOutFile) version
               dumm(1)=111           ! type of input
               dumm(2)=2             ! marker for constant theta plane
               dumm(3)=90.0_cp         ! surface constant
               dumm(4)=nFields       ! no of fields
               write(nOutFile) (real(dumm(n),kind=outp),n=1,4)

               dumm(1)=92.0          ! Field marker for AS vPhi
               dumm(2)=93.0          ! Field marker for Reynolds Force
               dumm(3)=94.0          ! Field marker for Reynolds Force
               write(nOutFile) (real(dumm(n),kind=outp),n=1,nFields)

               !------ Now other info about grid and parameters:
               write(nOutFile) runid     ! run identifier
               dumm( 1)=nSmax            ! total number of radial points
               dumm( 2)=nSmax            ! no of radial point in outer core
               dumm( 3)=1                ! no. of theta points
               dumm( 4)=n_phi_max        ! no. of phi points
               dumm( 5)=minc             ! imposed symmetry
               dumm( 6)=ra               ! control parameters
               dumm( 7)=ek               ! (for information only)
               dumm( 8)=pr               !      -"-
               dumm( 9)=prmag            !      -"-
               dumm(10)=radratio         ! ratio of inner / outer core
               dumm(11)=tScale           ! timescale
               write(nOutFile) (real(dumm(n),kind=outp),   n=1,11)
               write(nOutFile) (real(sZ(nSmax-n+1)/r_CMB,kind=outp),n=1,nSmax)
               write(nOutFile)  90.0_outp
               write(nOutFile) (real(phi(n),kind=outp), n=1,n_phi_max)

            end if ! Write Header ?

            dumm(1)=nGeosSets          ! time frame number for movie
            dumm(2)=time              ! time
            dumm(3)=0.0_cp
            dumm(4)=0.0_cp
            dumm(5)=0.0_cp
            dumm(6)=0.0_cp
            dumm(7)=0.0_cp
            dumm(8)=0.0_cp
            write(nOutFile) (real(dumm(n),kind=outp),n=1,8)

            write(nOutFile) ((CVz(nPhi,nS) ,nPhi=1,n_phi_max), nS=1,nSmax)
            write(nOutFile) ((CVor(nPhi,nS),nPhi=1,n_phi_max), nS=1,nSmax)
            write(nOutFile) ((CHel(nPhi,nS),nPhi=1,n_phi_max), nS=1,nSmax)

         end if

      end if
      if ( lVerbose ) write(*,*) '! End of getGeos!'

   end subroutine getEgeos
!----------------------------------------------------------------------------
   subroutine getDVptr(wS,dwS,ddwS,zS,dzS,rMin,rMax,rS,       &
        &              nZmax,nZmaxA,PlmS,dPlmS,OsinTS,lDeriv, &
        &              VrS,VtS,VpS,VorS,dpEk)
      !
      !  This subroutine calculates the three flow components VrS,VtS,VpS at
      !  a (r,theta,all phis) and (t,pi=theta, all phis). Here r=rS, PlmS=Plm(theta),
      !  dPlmS=sin(theta)*dTheta Plm(theta), and OsinTS=1/sin(theta).
      !  The flow is calculated for all n_phi_max azimuthal points used in the code,
      !  and for corresponding latitudes north and south of the equator.
      !  For lDeriv=.true. the subroutine also calculates dpEk and dzEk which
      !  are phi averages of (d Vr/d phi)**2 + (d Vtheta/ d phi)**2 + (d Vphi/ d phi)**2
      !  and (d Vr/d z)**2 + (d Vtheta/ d z)**2 + (d Vphi/ d z)**2, respectively.
      !  These two quantities are used to calculate z and phi scale of the flow in
      !  s_getEgeos.f
      !  
      !  .. note:: on input wS=w/r^2, dwS=dw/r, ddwS=ddw/r, zS=z/r
      !

      !--- Input variables:
      complex(cp), intent(in) :: wS(lm_max,n_r_max)
      complex(cp), intent(in) :: dwS(lm_max,n_r_max)
      complex(cp), intent(in) :: ddwS(lm_maxGeos,n_r_maxGeos)
      complex(cp), intent(in) :: zS(lm_maxGeos,n_r_maxGeos)
      complex(cp), intent(in) :: dzS(lm_maxGeos,n_r_maxGeos)
      real(cp),    intent(in) :: rMin,rMax  ! radial bounds
      integer,     intent(in) :: nZmax,nZmaxA ! number of (r,theta) points
      real(cp),    intent(in) :: rS(nZmaxA)
      real(cp),    intent(in) :: PlmS(lm_maxGeos,nZmaxA/2+1)
      real(cp),    intent(in) :: dPlmS(lm_maxGeos,nZmaxA/2+1)
      real(cp),    intent(in) :: OsinTS(nZmaxA/2+1)
      logical,     intent(in) ::  lDeriv
    
      !--- Output: function on azimuthal grid points defined by FT!
      real(cp), intent(out) :: VrS(nrpGeos,nZmaxA)
      real(cp), intent(out) :: VtS(nrpGeos,nZmaxA)
      real(cp), intent(out) :: VpS(nrpGeos,nZmaxA)
      real(cp), intent(out) :: VorS(nrpGeos,nZmaxA)
      real(cp), intent(out) :: dpEk(nZmaxA)
    
      !--- Local variables:
      real(cp) :: chebS(n_r_max)
      integer :: nS,nN,mc,lm,l,m,nCheb,nPhi,n
      real(cp) :: x,phiNorm,mapFac,OS,cosT,sinT,Or_e1,Or_e2
      complex(cp) :: Vr,Vt1,Vt2,Vp1,Vp2,Vor,Vot1,Vot2
      real(cp) :: VotS(nrpGeos,nZmaxA)
      complex(cp) :: wSr,dwSr,ddwSr,zSr,dzSr
    
      real(cp) :: dV(nrpGeos,2*nZmax)
      complex(cp) :: dp
    
    
      mapFac=two/(rMax-rMin)
      phiNorm=two*pi/n_phi_max
    
      do nS=1,nZmax
         do mc=1,nrpGeos
            VrS(mc,nS) =0.0_cp
            VtS(mc,nS) =0.0_cp
            VpS(mc,nS) =0.0_cp
            VorS(mc,nS)=0.0_cp
            VotS(mc,nS)=0.0_cp
         end do
      end do
    
      do nN=1,nZmax/2    ! Loop over all (r,theta) points in NHS
         nS=nZmax-nN+1   ! Southern counterpart !
    
         !------ Calculate Chebs:
         !------ Map r to cheb interval [-1,1]:
         !       and calculate the cheb polynomia:
         !       Note: the factor cheb_norm is needed
         !       for renormalisation. Its not needed if one used
         !       costf1 for the back transform.
         x=two*(rS(nN)-half*(rMin+rMax))/(rMax-rMin)
         chebS(1) =one*cheb_norm ! Extra cheb_norm cheap here
         chebS(2) =x*cheb_norm
         do nCheb=3,n_r_max
            chebS(nCheb)=two*x*chebS(nCheb-1)-chebS(nCheb-2)
         end do
         chebS(1)      =half*chebS(1)
         chebS(n_r_max)=half*chebS(n_r_max)
         Or_e2=one/rS(nN)**2
    
         do lm=1,lm_max     ! Sum over lms
            l =lm2l(lm)
            m =lm2m(lm)
            mc=lm2mc(lm)
            wSr  =zero
            dwSr =zero
            ddwSr=zero
            zSr  =zero
            dzSr =zero
            do nCheb=1,n_r_max
               wSr  =  wSr+  wS(lm,nCheb)*chebS(nCheb)
               dwSr = dwSr+ dwS(lm,nCheb)*chebS(nCheb)
               ddwSr=ddwSr+ddwS(lm,nCheb)*chebS(nCheb)
               zSr  =  zSr+  zS(lm,nCheb)*chebS(nCheb)
               dzSr = dzSr+ dzS(lm,nCheb)*chebS(nCheb)
            end do
            Vr  =  wSr* PlmS(lm,nN)
            Vt1 = dwSr*dPlmS(lm,nN)
            Vt2 =  zSr* PlmS(lm,nN)*dPhi(lm)
            Vp1 = dwSr* PlmS(lm,nN)*dPhi(lm)
            Vp2 = -zSr*dPlmS(lm,nN)
            Vor =  zSr* PlmS(lm,nN)*dLh(lm)
            Vot1= dzSr*dPlmS(lm,nN)
            Vot2= (wSr*Or_e2-ddwSr) * PlmS(lm,nN)*dPhi(lm)
            VrS(2*mc-1,nN) =VrS(2*mc-1,nN) + real(Vr)
            VrS(2*mc  ,nN) =VrS(2*mc  ,nN) +aimag(Vr)
            VtS(2*mc-1,nN) =VtS(2*mc-1,nN) + real(Vt1+Vt2)
            VtS(2*mc  ,nN) =VtS(2*mc  ,nN) +aimag(Vt1+Vt2)
            VpS(2*mc-1,nN) =VpS(2*mc-1,nN) + real(Vp1+Vp2)
            VpS(2*mc  ,nN) =VpS(2*mc  ,nN) +aimag(Vp1+Vp2)
            VorS(2*mc-1,nN)=VorS(2*mc-1,nN)+ real(Vor)
            VorS(2*mc  ,nN)=VorS(2*mc  ,nN)+aimag(Vor)
            VotS(2*mc-1,nN)=VotS(2*mc-1,nN)+ real(Vot1+Vot2)
            VotS(2*mc  ,nN)=VotS(2*mc  ,nN)+aimag(Vot1+Vot2)
            if ( mod(l+m,2) == 0 ) then
               VrS(2*mc-1,nS) =VrS(2*mc-1,nS) + real(Vr)
               VrS(2*mc  ,nS) =VrS(2*mc  ,nS) +aimag(Vr)
               VtS(2*mc-1,nS) =VtS(2*mc-1,nS) + real(Vt2-Vt1)
               VtS(2*mc  ,nS) =VtS(2*mc  ,nS) +aimag(Vt2-Vt1)
               VpS(2*mc-1,nS) =VpS(2*mc-1,nS) + real(Vp1-Vp2)
               VpS(2*mc  ,nS) =VpS(2*mc  ,nS) +aimag(Vp1-Vp2)
               VorS(2*mc-1,nS)=VorS(2*mc-1,nS)+ real(Vor)
               VorS(2*mc  ,nS)=VorS(2*mc  ,nS)+aimag(Vor)
               VotS(2*mc-1,nS)=VotS(2*mc-1,nS)+ real(Vot2-Vot1)
               VotS(2*mc  ,nS)=VotS(2*mc  ,nS)+aimag(Vot2-Vot1)
            else
               VrS(2*mc-1,nS) =VrS(2*mc-1,nS) - real(Vr)
               VrS(2*mc  ,nS) =VrS(2*mc  ,nS) -aimag(Vr)
               VtS(2*mc-1,nS) =VtS(2*mc-1,nS) + real(Vt1-Vt2)
               VtS(2*mc  ,nS) =VtS(2*mc  ,nS) +aimag(Vt1-Vt2)
               VpS(2*mc-1,nS) =VpS(2*mc-1,nS) + real(Vp2-Vp1)
               VpS(2*mc  ,nS) =VpS(2*mc  ,nS) +aimag(Vp2-Vp1)
               VorS(2*mc-1,nS)=VorS(2*mc-1,nS)- real(Vor)
               VorS(2*mc  ,nS)=VorS(2*mc  ,nS)-aimag(Vor)
               VotS(2*mc-1,nS)=VotS(2*mc-1,nS)+ real(Vot1-Vot2)
               VotS(2*mc  ,nS)=VotS(2*mc  ,nS)+aimag(Vot1-Vot2)
            end if
    
         end do
    
      end do
    
      if ( mod(nZmax,2) == 1 ) then ! Remaining equatorial point
         nS=(nZmax+1)/2
    
         x=two*(rS(nS)-half*(rMin+rMax))/(rMax-rMin)
         chebS(1)=one*cheb_norm ! Extra cheb_norm cheap here
         chebS(2)=x*cheb_norm
         do nCheb=3,n_r_max
            chebS(nCheb)=two*x*chebS(nCheb-1)-chebS(nCheb-2)
         end do
         chebS(1)      =half*chebS(1)
         chebS(n_r_max)=half*chebS(n_r_max)
         Or_e2=one/rS(nS)**2
    
         do lm=1,lm_max     ! Sum over lms
            l =lm2l(lm)
            m =lm2m(lm)
            mc=lm2mc(lm)
            wSr  =zero
            dwSr =zero
            ddwSr=zero
            zSr  =zero
            dzSr =zero
            do nCheb=1,n_r_max
               wSr  =  wSr+  wS(lm,nCheb)*chebS(nCheb)
               dwSr = dwSr+ dwS(lm,nCheb)*chebS(nCheb)
               ddwSr=ddwSr+ddwS(lm,nCheb)*chebS(nCheb)
               zSr  =  zSr+  zS(lm,nCheb)*chebS(nCheb)
               dzSr = dzSr+ dzS(lm,nCheb)*chebS(nCheb)
            end do
            Vr  =  wSr* PlmS(lm,nS)
            Vt1 = dwSr*dPlmS(lm,nS)
            Vt2 =  zSr* PlmS(lm,nS)*dPhi(lm)
            Vp1 = dwSr* PlmS(lm,nS)*dPhi(lm)
            Vp2 = -zSr*dPlmS(lm,nS)
            Vor =  zSr* PlmS(lm,nS)*dLh(lm)
            Vot1= dzSr*dPlmS(lm,nS)
            Vot2= (wSr*Or_e2-ddwSr) * PlmS(lm,nS)*dPhi(lm)
    
            VrS(2*mc-1, nN)=VrS(2*mc-1,nN) + real(Vr)
            VrS(2*mc  , nN)=VrS(2*mc  ,nN) +aimag(Vr)
            VtS(2*mc-1, nN)=VtS(2*mc-1,nN) + real(Vt1+Vt2)
            VtS(2*mc  , nN)=VtS(2*mc  ,nN) +aimag(Vt1+Vt2)
            VpS(2*mc-1, nN)=VpS(2*mc-1,nN) + real(Vp1+Vp2)
            VpS(2*mc  , nN)=VpS(2*mc  ,nN) +aimag(Vp1+Vp2)
            VorS(2*mc-1,nN)=VorS(2*mc-1,nN)+ real(Vor)
            VorS(2*mc  ,nN)=VorS(2*mc  ,nN)+aimag(Vor)
            VotS(2*mc-1,nN)=VotS(2*mc-1,nN)+ real(Vot1+Vot2)
            VotS(2*mc  ,nN)=VotS(2*mc  ,nN)+aimag(Vot1+Vot2)
         end do
    
      end if ! Equatorial point ?
    
      !--- Extra factors, contructing z-vorticity:
      do nS=1,(nZmax+1)/2 ! North HS
         OS   =OsinTS(nS)
         sinT =one/OS
         cosT =sqrt(one-sinT**2)
         Or_e1=one/rS(nS)
         Or_e2=Or_e1*Or_e1
         do mc=1,n_m_max
            VrS(2*mc-1,nS) =Or_e2*VrS(2*mc-1,nS)
            VrS(2*mc  ,nS) =Or_e2*VrS(2*mc  ,nS)
            VtS(2*mc-1,nS) =Or_e1*OS*VtS(2*mc-1,nS)
            VtS(2*mc  ,nS) =Or_e1*OS*VtS(2*mc  ,nS)
            VpS(2*mc-1,nS) =Or_e1*OS*VpS(2*mc-1,nS)
            VpS(2*mc  ,nS) =Or_e1*OS*VpS(2*mc  ,nS)
            VorS(2*mc-1,nS)=cosT*Or_e2*VorS(2*mc-1,nS) - Or_e1*VotS(2*mc-1,nS)
            VorS(2*mc  ,nS)=cosT*Or_e2*VorS(2*mc  ,nS) - Or_e1*VotS(2*mc  ,nS)
         end do
      end do
    
      do nS=(nZmax+1)/2+1,nZmax ! South HS
         OS   =OsinTS(nZmax-nS+1)
         sinT =one/OS
         cosT =-sqrt(one-sinT**2)
         Or_e1=one/rS(nZmax-nS+1)
         Or_e2=Or_e1*Or_e1
         do mc=1,n_m_max
            VrS(2*mc-1,nS) =Or_e2*VrS(2*mc-1,nS)
            VrS(2*mc  ,nS) =Or_e2*VrS(2*mc  ,nS)
            VtS(2*mc-1,nS) =Or_e1*OS*VtS(2*mc-1,nS)
            VtS(2*mc  ,nS) =Or_e1*OS*VtS(2*mc  ,nS)
            VpS(2*mc-1,nS) =Or_e1*OS*VpS(2*mc-1,nS)
            VpS(2*mc  ,nS) =Or_e1*OS*VpS(2*mc  ,nS)
            VorS(2*mc-1,nS)=cosT*Or_e2*VorS(2*mc-1,nS) - Or_e1*VotS(2*mc-1,nS)
            VorS(2*mc  ,nS)=cosT*Or_e2*VorS(2*mc  ,nS) - Or_e1*VotS(2*mc  ,nS)
         end do
      end do
    
      if ( lDeriv ) then
         !--- Calculate phi derivative in lm-space:
         do n=1,3
    
            do nS=1,nZmax
               if ( n == 1 ) then
                  dpEk(nS)=0.0_cp
                  do mc=1,nrpGeos
                     dp=ci*real((mc-1)*minc,cp) ! - i m
                     dV(2*mc-1,nS)= real(dp)*VrS(2*mc-1,nS)-aimag(dp)*VrS(2*mc,nS)
                     dV(2*mc  ,nS)=aimag(dp)*VrS(2*mc-1,nS)+ real(dp)*VrS(2*mc,nS)
                  end do
               else if ( n == 2 ) then
                  do mc=1,nrpGeos
                     dp=ci*real((mc-1)*minc,cp) ! - i m
                     dV(2*mc-1,nS)= real(dp)*VtS(2*mc-1,nS)-aimag(dp)*VtS(2*mc,nS)
                     dV(2*mc  ,nS)=aimag(dp)*VtS(2*mc-1,nS)+ real(dp)*VtS(2*mc,nS)
                  end do
               else if ( n == 3 ) then
                  do mc=1,nrpGeos
                     dp=ci*real((mc-1)*minc,cp) ! - i m
                     dV(2*mc-1,nS)= real(dp)*VpS(2*mc-1,nS)-aimag(dp)*VpS(2*mc,nS)
                     dV(2*mc  ,nS)=aimag(dp)*VpS(2*mc-1,nS)+ real(dp)*VpS(2*mc,nS)
                  end do
               end if
            end do
    
            !--- Transform m 2 phi for phi-derivative
            call fft_to_real(dV,nrpGeos,nZmax)
    
            !--- Phi average
            do nS=1,nZmax
               do nPhi=1,n_phi_max
                  if ( mod(nPhi,2) == 1 ) then
                     mc=(nPhi+1)/2
                     dpEk(nS)=dpEk(nS)+dV(2*mc-1,1)**2 ! Real part
                  else
                     mc=nPhi/2
                     dpEk(nS)=dpEk(nS)+dV(2*mc,1)**2 ! Imaginary part
                  end if
               end do
               if ( n == 3 ) dpEk(nS)=phiNorm*dpEk(nS)
            end do
    
         end do ! Loop over components
    
      end if
    
      !----- Transform m 2 phi for flow field:
      call fft_to_real(VrS,nrpGeos,nZmax)
      call fft_to_real(VtS,nrpGeos,nZmax)
      call fft_to_real(VpS,nrpGeos,nZmax)
      call fft_to_real(VorS,nrpGeos,nZmax)
    
   end subroutine getDVptr
!----------------------------------------------------------------------------
end module Egeos_mod
