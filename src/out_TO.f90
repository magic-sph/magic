module outTO_mod

   use iso_fortran_env, only: output_unit
   use parallel_mod
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, n_r_maxStr, n_theta_maxStr, l_max, &
       &                 n_theta_max, n_phi_max, minc, lStressMem,   &
       &                 lm_max
   use radial_functions, only: r_ICB, rscheme_oc, r, r_CMB, orho1, rscheme_oc
   use radial_data, only: nRstart, nRstop, radial_balance
   use physical_parameters, only: ra, ek, pr, prmag, radratio, LFfac
   use torsional_oscillations, only: BpzAS_Rloc, BspdAS_Rloc, BpsdAS_Rloc, &
       &                             BzpdAS_Rloc, BpzdAS_Rloc, dzCorLMr,   &
       &                             dzdVpLMr, dzddVpLMr, dzRstrLMr,       &
       &                             dzAstrLMr, dzStrLMr, dzLFLMr,         &
       &                             dzRstrLMr_Rloc, Bs2AS_Rloc,           &
       &                             V2AS_Rloc, BspAS_Rloc, BszAS_Rloc,    &
       &                             dzAstrLMr_Rloc, dzCorLMr_Rloc,        &
       &                             dzddVpLMr_Rloc, dzdVpLMr_Rloc,        &
       &                             dzLFLMr_Rloc, dzStrLMr_Rloc
   use num_param, only: tScale
   use blocking, only: lo_map, llm, ulm
   use horizontal_data, only: phi, sinTheta, theta_ord, gauss
   use logic, only: lVerbose, l_save_out
   use output_data, only: sDens, zDens, tag, log_file, runid, n_log_file
   use constants, only: pi, vol_oc, one, two, half, four
   use integration, only: rInt_R
   use plms_theta, only: plm_theta
   use TO_helpers, only: getPAStr, get_PAS, getAStr
   use useful, only: logWrite, abortRun
   use sht, only: spat_to_SH_axi
   use chebInt_mod, only: chebInt, chebIntInit
   use cosine_transform_odd

   implicit none

   private

   integer :: nTOZfile
   integer :: nSmax, nZmaxA
   integer :: nSstart, nSstop
   type(load), allocatable :: cyl_balance(:)

   real(cp) :: timeLast!,tNorm
   real(outp) :: timeAve

   !-- s-distributed arrays
   integer, allocatable :: nZmaxS_Sloc(:)
   real(cp), allocatable :: zZ_Sloc(:,:)
   real(outp), allocatable :: VpM_Sloc(:,:), dVpM_Sloc(:,:), LFM_Sloc(:,:)
   real(outp), allocatable :: AstrM_Sloc(:,:), RstrM_Sloc(:,:), CorM_Sloc(:,:)
   real(outp), allocatable :: StrM_Sloc(:,:), CLM_Sloc(:,:)

   !-- global arrays for outputs
   integer, allocatable :: nZmaxS(:)
   real(cp), allocatable :: zZ(:,:)
   real(outp), allocatable :: VpM(:,:), dVpM(:,:), LFM(:,:), AstrM(:,:)
   real(outp), allocatable :: RstrM(:,:), CorM(:,:), StrM(:,:), CLM(:,:)

   !-- Plms: Plm,sin
   real(cp), allocatable :: PlmS(:,:,:)
   real(cp), allocatable :: dPlmS(:,:,:)
   real(cp), allocatable :: OsinTS(:,:)
   real(cp), allocatable :: rZ(:,:)
   type(costf_odd_t), allocatable :: chebt_Z(:)

   !--
      !-- (l,r) Representation of the different contributions
   real(cp), allocatable :: dzVpLMr_loc(:,:)
   real(cp), allocatable :: V2LMr_Rloc(:,:)
   real(cp), allocatable :: Bs2LMr_Rloc(:,:)
   real(cp), allocatable :: BszLMr_Rloc(:,:)
   real(cp), allocatable :: BspLMr_Rloc(:,:)
   real(cp), allocatable :: BpzLMr_Rloc(:,:)
   real(cp), allocatable :: BspdLMr_Rloc(:,:)
   real(cp), allocatable :: BpsdLMr_Rloc(:,:)
   real(cp), allocatable :: BzpdLMr_Rloc(:,:)
   real(cp), allocatable :: BpzdLMr_Rloc(:,:)
   real(cp), allocatable :: dzVpLMr(:,:)
   real(cp), allocatable :: V2LMr(:,:)
   real(cp), allocatable :: Bs2LMr(:,:)
   real(cp), allocatable :: BszLMr(:,:)
   real(cp), allocatable :: BspLMr(:,:)
   real(cp), allocatable :: BpzLMr(:,:)
   real(cp), allocatable :: BspdLMr(:,:)
   real(cp), allocatable :: BpsdLMr(:,:)
   real(cp), allocatable :: BzpdLMr(:,:)
   real(cp), allocatable :: BpzdLMr(:,:)

   !-- Output files
   character(len=64) :: TOfileNhs,TOfileShs,movFile
   character(len=66) :: tayFile

   public :: initialize_outTO_mod, finalize_outTO_mod, outTO

contains

   subroutine initialize_outTO_mod

      !-- R-distributed arrays
      allocate( V2LMr_Rloc(l_max+1,nRstart:nRstop) )
      allocate( Bs2LMr_Rloc(l_max+1,nRstart:nRstop) )
      allocate( BszLMr_Rloc(l_max+1,nRstart:nRstop) )
      allocate( BspLMr_Rloc(l_max+1,nRstart:nRstop) )
      allocate( BpzLMr_Rloc(l_max+1,nRstart:nRstop) )
      allocate( BspdLMr_Rloc(l_max+1,nRstart:nRstop) )
      allocate( BpsdLMr_Rloc(l_max+1,nRstart:nRstop) )
      allocate( BzpdLMr_Rloc(l_max+1,nRstart:nRstop) )
      allocate( BpzdLMr_Rloc(l_max+1,nRstart:nRstop) )
      bytes_allocated=bytes_allocated+9*(nRstop-nRstart+1)*(l_max+1)* &
      &               SIZEOF_DEF_REAL

      !-- Global arrays
      allocate( dzVpLMr_loc(l_max+1,n_r_max) )
      allocate( dzVpLMr(l_max+1,n_r_max) )
      allocate( V2LMr(l_max+1,n_r_max) )
      allocate( Bs2LMr(l_max+1,n_r_max) )
      allocate( BszLMr(l_max+1,n_r_max) )
      allocate( BspLMr(l_max+1,n_r_max) )
      allocate( BpzLMr(l_max+1,n_r_max) )
      allocate( BspdLMr(l_max+1,n_r_max) )
      allocate( BpsdLMr(l_max+1,n_r_max) )
      allocate( BzpdLMr(l_max+1,n_r_max) )
      allocate( BpzdLMr(l_max+1,n_r_max) )
      bytes_allocated=bytes_allocated+11*n_r_max*(l_max+1)*SIZEOF_DEF_REAL

      nSmax=n_r_max+int(r_ICB*real(n_r_max,cp))
      nSmax=int(sDens*nSmax)
      nZmaxA=2*nSmax

      !-- Distribute over the ranks
      allocate ( cyl_balance(0:n_procs-1) )
      call getBlocks(cyl_balance, nSmax, n_procs)
      nSstart = cyl_balance(rank)%nStart
      nSstop = cyl_balance(rank)%nStop

      !-- Allocate s-distributed arrays
      allocate( PlmS(l_max+1,nZmaxA/2+1,nSstart:nSstop) )
      allocate( dPlmS(l_max+1,nZmaxA/2+1,nSstart:nSstop) )
      bytes_allocated = bytes_allocated + &
                        2*(l_max+1)*(nZmaxA/2+1)*(nSstop-nSstart+1)*SIZEOF_DEF_REAL
      allocate( OsinTS(nZmaxA/2+1,nSstart:nSstop) )
      bytes_allocated = bytes_allocated + (nZmaxA/2+1)*(nSstop-nSstart+1)* &
      &                 SIZEOF_DEF_REAL

      allocate( VpM_Sloc(nZmaxA,nSstart:nSstop) )
      allocate( dVpM_Sloc(nZmaxA,nSstart:nSstop) )
      allocate( LFM_Sloc(nZmaxA,nSstart:nSstop) )
      allocate( AstrM_Sloc(nZmaxA,nSstart:nSstop) )
      allocate( RstrM_Sloc(nZmaxA,nSstart:nSstop) )
      allocate( CorM_Sloc(nZmaxA,nSstart:nSstop) )
      allocate( StrM_Sloc(nZmaxA,nSstart:nSstop) )
      allocate( CLM_Sloc(nZmaxA,nSstart:nSstop) )
      bytes_allocated = bytes_allocated + 8*nZmaxA*(nSstop-nSstart+1)* &
      &                 SIZEOF_OUT_REAL

      allocate( chebt_Z(nSstart:nSstop) )
      allocate( zZ_Sloc(nZmaxA,nSstart:nSstop) )
      bytes_allocated = bytes_allocated + nZmaxA*(nSstop-nSstart+1)* &
      &                 SIZEOF_DEF_REAL
      allocate( rZ(nZmaxA/2+1,nSstart:nSstop) )
      bytes_allocated = bytes_allocated + (nZmaxA/2+1)*(nSstop-nSstart+1)* &
      &                 SIZEOF_DEF_REAL
      allocate( nZmaxS_Sloc(nSstart:nSstop) )
      bytes_allocated = bytes_allocated+(nSstop-nSstart+1)*SIZEOF_INTEGER

      !-- Global arrays for outputs
      if ( rank == 0 ) then
         allocate ( nZmaxS(nSmax) )
         allocate( zZ(nZmaxA,nSmax) )
         allocate( VpM(nZmaxA,nSmax), dVpM(nZmaxA,nSmax) )
         allocate( LFM(nZmaxA,nSmax), AstrM(nZmaxA,nSmax) )
         allocate( RstrM(nZmaxA,nSmax), CorM(nZmaxA,nSmax) )
         allocate( StrM(nZmaxA,nSmax), CLM(nZmaxA,nSmax) )

         bytes_allocated = bytes_allocated+8*nZmaxA*nSmax*SIZEOF_OUT_REAL+&
         &                 nZmaxA*nSmax*SIZEOF_DEF_REAL+nSmax*SIZEOF_INTEGER
      else
         allocate ( nZmaxS(1) )
         allocate( zZ(1,1), VpM(1,1), dVpM(1,1), LFM(1,1), AstrM(1,1) )
         allocate( RstrM(1,1), CorM(1,1), StrM(1,1), CLM(1,1) )

      end if

      TOfileNhs='TOnhs.'//tag
      TOfileShs='TOshs.'//tag
      movFile  ='TO_mov.'//tag
      tayFile  ='TaySphere.'//tag

   end subroutine initialize_outTO_mod
!----------------------------------------------------------------------------
   subroutine finalize_outTO_mod
      !
      ! Memory deallocation
      !

      deallocate( PlmS, dPlmS, OsinTS, VpM_Sloc, LFM_Sloc, dVpM_Sloc)
      deallocate( AstrM_Sloc, RstrM_Sloc, CorM_Sloc )
      deallocate( StrM_Sloc, CLM_Sloc, zZ_Sloc, rZ, nZmaxS_Sloc )

      deallocate( nZmaxS, zZ, VpM, dVpM, AstrM, RstrM, CorM, LFM, StrM )

      deallocate( cyl_balance )

   end subroutine finalize_outTO_mod
!----------------------------------------------------------------------------
   subroutine outTO(time,n_time_step,eKin,eKinTAS,nTOsets,nTOmovSets, &
              &     nTOrmsSets,lTOmov,lTOrms,lTOZwrite,z,omega_ic,omega_ma)
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
      logical,          intent(in) :: lTOmov
      complex(cp),      intent(in) :: z(llm:ulm,n_r_max)
      real(cp),         intent(in) :: omega_ic, omega_ma
      integer,          intent(inout) :: nTOsets, nTOmovSets, nTOrmsSets
      logical,          intent(inout) :: lTOrms, lTOZwrite

      !-- Output field:
      real(cp) :: fOut(n_theta_maxStr*n_r_maxStr)

      !-- Local variables:
      logical :: lTC,lStopRun

      !---- Work array:
      real(cp) :: workA(l_max+1,n_r_maxStr)

      integer :: lm,l,m ! counter for degree and order
      integer :: nOutFile, nOutFile2

      integer :: nS,nSI
      real(cp) :: zNorm  ! Norm z interval
      integer :: nNorm  ! No. of grid points for norm interval
      real(cp) :: zMin,zMax!,help ! integration boundarie, help variable
      logical :: lAS    ! .true. if axisymmetric (m=0) functions
      real(cp) :: sZ(nSmax),dsZ ! cylindrical radius s and s-step
      real(cp) :: h(nSmax),Oh(nSmax)
      real(cp) :: Os2(nSmax)

      integer :: nR     ! counter for radial grid point
      integer :: n      ! counter for theta blocks
      integer :: nOut,nFields ! counter for output fields
      integer :: nTheta ! counter for all thetas
      integer :: nThetaOrd ! counter for ordered thetas
      integer :: nThetaNHS

      integer :: nZ,nZmax,nZmaxNS!,nZP

      !-- S-distributed arrays
      real(cp) :: VpS_Sloc(nZmaxA,nSstart:nSstop), dVpS_Sloc(nZmaxA,nSstart:nSstop)
      real(cp) :: LFS_Sloc(nZmaxA,nSstart:nSstop), CorS_Sloc(nZmaxA,nSstart:nSstop)
      real(cp) :: RstrS_Sloc(nZmaxA,nSstart:nSstop), AstrS_Sloc(nZmaxA,nSstart:nSstop)
      real(cp) :: StrS_Sloc(nZmaxA,nSstart:nSstop)
      real(cp) :: VpIntN_Sloc(nSstart:nSstop), VpIntS_Sloc(nSstart:nSstop)
      real(cp) :: SVpIntN_Sloc(nSstart:nSstop), SVpIntS_Sloc(nSstart:nSstop)
      real(cp) :: SBspIntN_Sloc(nSstart:nSstop), SBspIntS_Sloc(nSstart:nSstop)
      real(cp) :: SBs2IntN_Sloc(nSstart:nSstop), SBs2IntS_Sloc(nSstart:nSstop)
      real(cp) :: TauBN_Sloc(nSstart:nSstop),TauBS_Sloc(nSstart:nSstop)
      real(cp) :: dTauBN_Sloc(nSstart:nSstop),dTauBS_Sloc(nSstart:nSstop)
      real(cp) :: dTTauBN_Sloc(nSstart:nSstop),dTTauBS_Sloc(nSstart:nSstop)
      real(cp) :: Bs2IntN_Sloc(nSstart:nSstop) ,Bs2IntS_Sloc(nSstart:nSstop)
      real(cp) :: dVpIntN_Sloc(nSstart:nSstop) ,dVpIntS_Sloc(nSstart:nSstop)
      real(cp) :: ddVpIntN_Sloc(nSstart:nSstop),ddVpIntS_Sloc(nSstart:nSstop)
      real(cp) :: VpRIntN_Sloc(nSstart:nSstop) ,VpRIntS_Sloc(nSstart:nSstop)
      real(cp) :: LFIntN_Sloc(nSstart:nSstop)  ,LFIntS_Sloc(nSstart:nSstop)
      real(cp) :: RstrIntN_Sloc(nSstart:nSstop),RstrIntS_Sloc(nSstart:nSstop)
      real(cp) :: AstrIntN_Sloc(nSstart:nSstop),AstrIntS_Sloc(nSstart:nSstop)
      real(cp) :: StrIntN_Sloc(nSstart:nSstop) ,StrIntS_Sloc(nSstart:nSstop)
      real(cp) :: TayIntN_Sloc(nSstart:nSstop) ,TayIntS_Sloc(nSstart:nSstop)
      real(cp) :: BspdIntN_Sloc(nSstart:nSstop),BspdIntS_Sloc(nSstart:nSstop)
      real(cp) :: BpsdIntN_Sloc(nSstart:nSstop),BpsdIntS_Sloc(nSstart:nSstop)

      !-- S-distributed arrays that don't need to be gathered
      real(cp) :: V2IntS(nSstart:nSstop)  ,V2IntN(nSstart:nSstop)
      real(cp) :: BspIntN(nSstart:nSstop) ,BspIntS(nSstart:nSstop)
      real(cp) :: TayRIntN(nSstart:nSstop),TayRIntS(nSstart:nSstop)
      real(cp) :: TayVIntN(nSstart:nSstop),TayVIntS(nSstart:nSstop)

      !-- Global arrays
      real(cp) :: VpS(nZmaxA,nSmax), dVpS(nZmaxA,nSmax), RstrS(nZmaxA,nSmax)
      real(cp) :: AstrS(nZmaxA,nSmax), LFS(nZmaxA,nSmax), CorS(nZmaxA,nSmax)
      real(cp) :: StrS(nZmaxA,nSmax)
      real(cp) :: ddVpS(nZmaxA)
      real(cp) :: V2S(nZmaxA), Bs2S(nZmaxA), BspS(nZmaxA), BspdS(nZmaxA)
      real(cp) :: BpsdS(nZmaxA), TayS(nZmaxA), TayRS(nZmaxA), TayVS(nZmaxA)
      real(cp) :: VpIntN(nSmax), VpIntS(nSmax)    ! integration results
      real(cp) :: dVpIntN(nSmax) ,dVpIntS(nSmax)   ! integration results
      real(cp) :: ddVpIntN(nSmax),ddVpIntS(nSmax)  ! integration results
      real(cp) :: VpRIntN(nSmax) ,VpRIntS(nSmax)   ! for different s and
      real(cp) :: LFIntN(nSmax)  ,LFIntS(nSmax)
      real(cp) :: RstrIntN(nSmax),RstrIntS(nSmax)
      real(cp) :: AstrIntN(nSmax),AstrIntS(nSmax)
      real(cp) :: StrIntN(nSmax) ,StrIntS(nSmax)
      real(cp) :: TayIntN(nSmax) ,TayIntS(nSmax)
      real(cp) :: BspdIntN(nSmax),BspdIntS(nSmax)
      real(cp) :: BpsdIntN(nSmax),BpsdIntS(nSmax)
      real(cp) :: Bs2IntN(nSmax) ,Bs2IntS(nSmax)
      real(cp) :: SVpIntN(nSmax) ,SVpIntS(nSmax)   ! help arrays and values for
      real(cp) :: SBs2IntN(nSmax),SBs2IntS(nSmax)  ! differentiation in s
      real(cp) :: SBspIntN(nSmax),SBspIntS(nSmax)
      real(cp) :: dSVpIntN, dSVpIntS
      real(cp) :: d2SVpIntN,d2SVpIntS
      real(cp) :: dSBspIntN,dSBspIntS
      real(cp) :: dSBs2IntN,dSBs2IntS
      real(cp) :: TauN(nSmax),TauS(nSmax)          ! Taylor integral
      real(cp) :: TauBN(nSmax),TauBS(nSmax)
      real(cp) :: dTauBN(nSmax),dTauBS(nSmax)
      real(cp) :: dTTauN(nSmax),dTTauS(nSmax)      ! time change of Tau...
      real(cp) :: dTTauBN(nSmax),dTTauBS(nSmax)

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
      real(cp) :: outBlock(n_theta_max)
      real(outp) :: dt

      character(len=255) :: message
      character(len=64) :: version,fileName
      integer :: nFieldSize,nPos
      integer :: n_toz_file, n_tozm_file

      !integer :: nLines
      real(cp) :: dumm(12)

      !-- Huge arrays for time average ....
      logical :: l_TOZave

      !-- For TOZ output files:
      character(len=14) :: string

#ifdef WITH_MPI
      real(cp) :: global_sum
      integer :: i,sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1)
#endif

      if ( lVerbose ) write(output_unit,*) '! Starting outTO!'

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
      lAS=.true.

      !--- Transform to lm-space for all radial grid points:

      do nR=nRstart,nRstop
         call spat_to_SH_axi(V2AS_Rloc(:,nR),V2LMr_Rloc(:,nR))
         call spat_to_SH_axi(Bs2AS_Rloc(:,nR),Bs2LMr_Rloc(:,nR))
         call spat_to_SH_axi(BszAS_Rloc(:,nR),BszLMr_Rloc(:,nR))
         call spat_to_SH_axi(BspAS_Rloc(:,nR),BspLMr_Rloc(:,nR))
         call spat_to_SH_axi(BpzAS_Rloc(:,nR),BpzLMr_Rloc(:,nR))
         call spat_to_SH_axi(BspdAS_Rloc(:,nR),BspdLMr_Rloc(:,nR))
         call spat_to_SH_axi(BpsdAS_Rloc(:,nR),BpsdLMr_Rloc(:,nR))
         call spat_to_SH_axi(BzpdAS_Rloc(:,nR),BzpdLMr_Rloc(:,nR))
         call spat_to_SH_axi(BpzdAS_Rloc(:,nR),BpzdLMr_Rloc(:,nR))
      end do ! Loop over radial grid points


      !-- All gather everything at this stage
      call outTO_allgather_Rloc()

      do nR=1,n_r_max
         dzVpLMr_loc(:,nR)=0.0_cp
         do lm=llm,ulm
            l=lo_map%lm2l(lm)
            m=lo_map%lm2m(lm)
            if ( m == 0 ) dzVpLMr_loc(l+1,nR)=orho1(nR)*real(z(lm,nR))
         end do
#ifdef WITH_MPI
         call MPI_Allreduce(dzVpLMr_loc(:,nR), dzVPLMr(:,nR), l_max+1, &
              &             MPI_DEF_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
         dzVPLMr(:,nR)=dzVpLMr_loc(:,nR)
#endif
      end do


      !---- Transform the contributions to cheb space for z-integral:
      call rscheme_oc%costf1(dzVpLMr,l_max+1,1,l_max+1,workA)
      call rscheme_oc%costf1(V2LMr,l_max+1,1,l_max+1,workA)
      call rscheme_oc%costf1(dzdVpLMr,l_max+1,1,l_max+1,workA)
      call rscheme_oc%costf1(dzddVpLMr,l_max+1,1,l_max+1,workA)
      call rscheme_oc%costf1(Bs2LMr,l_max+1,1,l_max+1,workA)
      call rscheme_oc%costf1(BszLMr,l_max+1,1,l_max+1,workA)
      call rscheme_oc%costf1(BspLMr,l_max+1,1,l_max+1,workA)
      call rscheme_oc%costf1(BpzLMr,l_max+1,1,l_max+1,workA)
      call rscheme_oc%costf1(dzRstrLMr,l_max+1,1,l_max+1,workA)
      call rscheme_oc%costf1(dzAstrLMr,l_max+1,1,l_max+1,workA)
      call rscheme_oc%costf1(dzStrLMr,l_max+1,1,l_max+1,workA)
      call rscheme_oc%costf1(dzLFLMr,l_max+1,1,l_max+1,workA)
      call rscheme_oc%costf1(dzCorLMr,l_max+1,1,l_max+1,workA)
      call rscheme_oc%costf1(BspdLMr,l_max+1,1,l_max+1,workA)
      call rscheme_oc%costf1(BpsdLMr,l_max+1,1,l_max+1,workA)
      call rscheme_oc%costf1(BpzdLMr,l_max+1,1,l_max+1,workA)
      call rscheme_oc%costf1(BzpdLMr,l_max+1,1,l_max+1,workA)


      !-- First loop to fill some help arrays
      dsZ   =r_CMB/real(nSmax,cp)  ! Step in s controlled by nSmax
      nSI   =0
      do nS=1,nSmax
         sZ(nS)=(nS-half)*dsZ
         if ( sZ(nS) < r_ICB .and. nS > nSI ) nSI=nS

         zMax=sqrt(r_CMB*r_CMB-sZ(nS)*sZ(nS))
         if ( sZ(nS) < r_ICB ) then
            zMin=sqrt(r_ICB*r_ICB-sZ(nS)*sZ(nS))
         else
            zMin=-zMax
         end if
         h(nS) =zMax-zMin
         Oh(nS)=one/h(nS)
         Os2(nS)=one/(sZ(nS)*sZ(nS))
      end do

#ifdef WITH_MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

      lStopRun=.false.
      outer: do nS=nSstart,nSstop

         !------ Get integral boundaries for this s:
         zMax=sqrt(r_CMB*r_CMB-sZ(nS)*sZ(nS))
         if ( sZ(nS) < r_ICB ) then
            lTC=.true.
            zMin=sqrt(r_ICB*r_ICB-sZ(nS)*sZ(nS))
         else
            lTC=.false.
            zMin=-zMax
         end if

         if ( nTOsets == 1 ) then
            !------ Initialize integration for NHS:
            !       Each processor calculates Cheb transform data
            !       for HIS nS and the Plms along the Cylinder
            !       chebIntInit returns zZ_Sloc,nZmaxS,chebt_Z
            !       Note that this returns z in the MAGIC way,
            !       starting with zMax, ending with zMin
            !       z(:,nS)=zMin, z(nZmax,nS)=zMax
            call chebIntInit(zMin,zMax,zNorm,nNorm,                   &
                 &            nZmaxA,zZ_Sloc(:,nS),nZmaxS_Sloc(nS),chebt_Z(nS))

            !--- Points in nothers halfsphere
            if ( lTC ) then
               nZmax=nZmaxS_Sloc(nS)  ! nZmax point in each polar region
               if ( 2*nZmax > nZmaxA ) then
                  write(output_unit,*) '! nZmaxA too small in outTO!'
                  write(output_unit,*) '! Should be at least:',2*nZmax
                  lStopRun=.true.
               end if
            else
               nZmax=(nZmaxS_Sloc(nS)-1)/2+1 ! Odd point number !
               ! all together nZmaxS_Sloc(nS) from
               ! south to north including equator
               if ( nZmaxS_Sloc(nS) > nZmaxA ) then
                  write(output_unit,*) '! nZmaxA too small in outTO!'
                  write(output_unit,*) '! Should be at least:',nZmaxS_Sloc(nS)
                  lStopRun=.true.
               end if
            end if
            if ( lStopRun ) cycle outer
            do nZ=1,nZmax
               rZ(nZ,nS)    =sqrt(zZ_Sloc(nZ,nS)**2+sZ(nS)**2)
               thetaZ       =atan2(sZ(nS),zZ_Sloc(nZ,nS))
               OsinTS(nZ,nS)=one/sin(thetaZ)
               call plm_theta(thetaZ,l_max,0,minc,                    &
                    &         PlmS(:,nZ,nS),dPlmS(:,nZ,nS),l_max+1,2)
            end do
         end if

         !--------- Get the flow components for all northern and
         !          southern thetas and all phis:

         nZmax=nZmaxS_Sloc(nS)
         if ( lTC ) then
            nZmaxNS=2*nZmax
         else
            nZmaxNS=nZmax
         end if

         call getPAStr(VpS_Sloc(:,nS),dzVpLMr,nZmaxNS,nZmaxA,l_max,    &
              &        r_ICB,r_CMB,n_r_max,rZ(:,nS),dPlmS(:,:,nS),OsinTS(:,nS))
         call getPAStr(dVpS_Sloc(:,nS),dzdVpLMr,nZmaxNS,nZmaxA,l_max,  &
              &        r_ICB,r_CMB,n_r_max,rZ(:,nS),dPlmS(:,:,nS),OsinTS(:,nS))
         call getPAStr(ddVpS,dzddVpLMr,nZmaxNS,nZmaxA,l_max,r_ICB,     &
              &        r_CMB,n_r_max,rZ(:,nS),dPlmS(:,:,nS),OsinTS(:,nS))
         call getPAStr(RstrS_Sloc(:,nS),dzRstrLMr,nZmaxNS,nZmaxA,l_max,&
              &        r_ICB,r_CMB,n_r_max,rZ(:,nS),dPlmS(:,:,nS),OsinTS(:,nS))
         call getPAStr(AstrS_Sloc(:,nS),dzAstrLMr,nZmaxNS,nZmaxA,l_max,&
              &        r_ICB,r_CMB,n_r_max,rZ(:,nS),dPlmS(:,:,nS),OsinTS(:,nS))
         call getPAStr(StrS_Sloc(:,nS),dzStrLMr,nZmaxNS,nZmaxA,l_max,  &
              &        r_ICB,r_CMB,n_r_max,rZ(:,nS),dPlmS(:,:,nS),OsinTS(:,nS))
         call getPAStr(LFS_Sloc(:,nS),dzLFLMr,nZmaxNS,nZmaxA,l_max,    &
              &        r_ICB,r_CMB,n_r_max,rZ(:,nS),dPlmS(:,:,nS),OsinTS(:,nS))
         call getPAStr(CorS_Sloc(:,nS),dzCorLMr,nZmaxNS,nZmaxA,l_max,  &
              &        r_ICB,r_CMB,n_r_max,rZ(:,nS),dPlmS(:,:,nS),OsinTS(:,nS))
         do nZ=1,nZmaxNS
            TayS(nZ) =abs(LFS_Sloc(nZ,nS))
            TayRS(nZ)=abs(RstrS_Sloc(nZ,nS))
            TayVS(nZ)=abs(StrS_Sloc(nZ,nS))
         end do

         call getAStr(V2S,V2LMr,nZmaxNS,nZmaxA,l_max,r_ICB,r_CMB,     &
              &       n_r_max,rZ(:,nS),PlmS(:,:,nS))
         call getAStr(Bs2S,Bs2LMr,nZmaxNS,nZmaxA,l_max,r_ICB,r_CMB,   &
              &       n_r_max,rZ(:,nS),PlmS(:,:,nS))
         call getAStr(BspS,BspLMr,nZmaxNS,nZmaxA,l_max,r_ICB,r_CMB,   &
              &       n_r_max,rZ(:,nS),PlmS(:,:,nS))
         call getAStr(BspdS,BspdLMr,nZmaxNS,nZmaxA,l_max,r_ICB,r_CMB, &
              &       n_r_max,rZ(:,nS),PlmS(:,:,nS))
         call getAStr(BpsdS,BpsdLMr,nZmaxNS,nZmaxA,l_max,r_ICB,r_CMB, &
              &       n_r_max,rZ(:,nS),PlmS(:,:,nS))

         if ( l_TOZave ) then
            if ( nTOsets == 1 ) then
               timeAve=1.e0_outp
               do nZ=1,nZmaxNS
                  VpM_Sloc(nZ,nS)  =real(VpS_Sloc(nZ,nS),kind=outp)
                  dVpM_Sloc(nZ,nS) =real(dVpS_Sloc(nZ,nS),kind=outp)
                  LFM_Sloc(nZ,nS)  =real(LFfac*LFS_Sloc(nZ,nS),kind=outp)
                  RstrM_Sloc(nZ,nS)=real(RstrS_Sloc(nZ,nS),kind=outp)
                  AstrM_Sloc(nZ,nS)=real(AstrS_Sloc(nZ,nS),kind=outp)
                  StrM_Sloc(nZ,nS) =real(StrS_Sloc(nZ,nS),kind=outp)
                  CorM_Sloc(nZ,nS) =real(CorS_Sloc(nZ,nS),kind=outp)
                  CLM_Sloc(nZ,nS)  =real(CorS_Sloc(nZ,nS)+ &
                  &                 LFfac*LFS_Sloc(nZ,nS),kind=outp)
               end do
            else
               dt=real(time-timeLast,kind=outp)
               if ( nS == nSstart ) timeAve=timeAve+dt
               do nZ=1,nZmaxNS
                  VpM_Sloc(nZ,nS)  =VpM_Sloc(nZ,nS)  +               &
                  &                 dt*real(VpS_Sloc(nZ,nS),kind=outp)
                  dVpM_Sloc(nZ,nS) =dVpM_Sloc(nZ,nS) +               &
                  &                 dt*real(dVpS_Sloc(nZ,nS),kind=outp)
                  LFM_Sloc(nZ,nS)  =LFM_Sloc(nZ,nS)  +               &
                  &                 dt*real(LFfac*LFS_Sloc(nZ,nS),kind=outp)
                  RstrM_Sloc(nZ,nS)=RstrM_Sloc(nZ,nS)+               &
                  &                 dt*real(RstrS_Sloc(nZ,nS),kind=outp)
                  AstrM_Sloc(nZ,nS)=AstrM_Sloc(nZ,nS)+               &
                  &                 dt*real(AstrS_Sloc(nZ,nS),kind=outp)
                  StrM_Sloc(nZ,nS) =StrM_Sloc(nZ,nS) +               &
                  &                 dt*real(StrS_Sloc(nZ,nS),kind=outp)
                  CorM_Sloc(nZ,nS) =CorM_Sloc(nZ,nS) +               &
                  &                 dt*real(CorS_Sloc(nZ,nS),kind=outp)
                  CLM_Sloc(nZ,nS)  =CLM_Sloc(nZ,nS)  +               &
                  &                 dt*real(CorS_Sloc(nZ,nS)+        &
                  &                 LFfac*LFS_Sloc(nZ,nS),kind=outp)
               end do
            end if
         end if

         !--- Z-integrals:
         VpIntN_Sloc(nS)  =chebInt(VpS_Sloc(:,nS),zMin,zMax,nZmax,nZmaxA,chebt_Z(nS))
         dVpIntN_Sloc(nS) =chebInt(dVpS_Sloc(:,nS),zMin,zMax,nZmax,nZmaxA,chebt_Z(nS))
         ddVpIntN_Sloc(nS)=chebInt(ddVpS,zMin,zMax,nZmax,nZmaxA,chebt_Z(nS))
         LFIntN_Sloc(nS)  =chebInt(LFS_Sloc(:,nS),zMin,zMax,nZmax,nZmaxA,chebt_Z(nS))
         TayIntN_Sloc(nS) =chebInt(TayS,zMin,zMax,nZmax,nZmaxA,chebt_Z(nS))
         TayRIntN(nS)     =chebInt(TayRS,zMin,zMax,nZmax,nZmaxA,chebt_Z(nS))
         TayVIntN(nS)     =chebInt(TayVS,zMin,zMax,nZmax,nZmaxA,chebt_Z(nS))
         RstrIntN_Sloc(nS)=chebInt(RstrS_Sloc(:,nS),zMin,zMax,nZmax,nZmaxA,chebt_Z(nS))
         AstrIntN_Sloc(nS)=chebInt(AstrS_Sloc(:,nS),zMin,zMax,nZmax,nZmaxA,chebt_Z(nS))
         StrIntN_Sloc(nS) =chebInt(StrS_Sloc(:,nS),zMin,zMax,nZmax,nZmaxA,chebt_Z(nS))
         V2IntN(nS)       =chebInt(V2S,zMin,zMax,nZmax,nZmaxA,chebt_Z(nS))
         Bs2IntN_Sloc(nS) =chebInt(Bs2S,zMin,zMax,nZmax,nZmaxA,chebt_Z(nS))
         BspIntN(nS)      =chebInt(BspS,zMin,zMax,nZmax,nZmaxA,chebt_Z(nS))
         BspdIntN_Sloc(nS)=chebInt(BspdS,zMin,zMax,nZmax,nZmaxA,chebt_Z(nS))
         BpsdIntN_Sloc(nS)=chebInt(BpsdS,zMin,zMax,nZmax,nZmaxA,chebt_Z(nS))

         if ( V2IntN(nS) < 0.0_cp ) then
            VpRIntN_Sloc(nS)=one
         else
            VpRIntN_Sloc(nS)=abs(VpIntN_Sloc(nS))/sqrt(V2IntN(nS))
         end if
         VpRIntN_Sloc(nS) =min(one,VpRIntN_Sloc(nS))
         TayIntN_Sloc(nS) =LFIntN_Sloc(nS)/TayIntN_Sloc(nS)
         TayRIntN(nS)=RstrIntN_Sloc(nS)/TayRIntN(nS)
         TayVIntN(nS)=StrIntN_Sloc(nS)/TayVIntN(nS)

         !--- Z-integration inside northern TC:
         if ( lTC ) then
            VpIntS_Sloc(nS)  =chebInt(VpS_Sloc(nZmax+1,nS),zMin,zMax,nZmax,&
            &                         nZmaxA,chebt_Z(nS))
            dVpIntS_Sloc(nS) =chebInt(dVpS_Sloc(nZmax+1,nS),zMin,zMax,nZmax,&
            &                         nZmaxA,chebt_Z(nS))
            ddVpIntS_Sloc(nS)=chebInt(ddVpS(nZmax+1),zMin,zMax,nZmax,nZmaxA,&
            &                         chebt_Z(nS))
            LFIntS_Sloc(nS)  =chebInt(LFS_Sloc(nZmax+1,nS),zMin,zMax,nZmax, &
            &                         nZmaxA,chebt_Z(nS))
            TayIntS_Sloc(nS) =chebInt(TayS(nZmax+1),zMin,zMax,nZmax,nZmaxA, &
            &                         chebt_Z(nS))
            TayRIntS(nS)     =chebInt(TayRS(nZmax+1),zMin,zMax,nZmax,nZmaxA,&
            &                         chebt_Z(nS))
            TayVIntS(nS)     =chebInt(TayVS(nZmax+1),zMin,zMax,nZmax,nZmaxA,&
            &                         chebt_Z(nS))
            RstrIntS_Sloc(nS)=chebInt(RstrS_Sloc(nZmax+1,nS),zMin,zMax,nZmax,&
            &                         nZmaxA,chebt_Z(nS))
            AstrIntS_Sloc(nS)=chebInt(AstrS_Sloc(nZmax+1,nS),zMin,zMax,nZmax,&
            &                         nZmaxA,chebt_Z(nS))
            StrIntS_Sloc(nS) =chebInt(StrS_Sloc(nZmax+1,nS),zMin,zMax,nZmax, &
            &                         nZmaxA,chebt_Z(nS))
            V2IntS(nS)       =chebInt(V2S(nZmax+1),zMin,zMax,nZmax,nZmaxA, &
            &                         chebt_Z(nS))
            Bs2IntS_Sloc(nS) =chebInt(Bs2S(nZmax+1),zMin,zMax,nZmax,nZmaxA,&
            &                         chebt_Z(nS))
            BspIntS(nS)      =chebInt(BspS(nZmax+1),zMin,zMax,nZmax,nZmaxA, &
            &                         chebt_Z(nS))
            BspdIntS_Sloc(nS)=chebInt(BspdS(nZmax+1),zMin,zMax,nZmax,nZmaxA,&
            &                         chebt_Z(nS))
            BpsdIntS_Sloc(nS)=chebInt(BpsdS(nZmax+1),zMin,zMax,nZmax,nZmaxA,&
            &                         chebt_Z(nS))
            if ( V2IntS(nS) < 0.0_cp ) then
               VpRIntS_Sloc(nS)=one
            else
               VpRIntS_Sloc(nS)=abs(VpIntS_Sloc(nS))/sqrt(V2IntS(nS))
            end if
            VpRIntS_Sloc(nS) =min(one,VpRIntS_Sloc(nS))
            TayIntS_Sloc(nS) =LFIntS_Sloc(nS)/TayIntS_Sloc(nS)
            TayRIntS(nS)     =RstrIntS_Sloc(nS)/TayRIntS(nS)
            TayVIntS(nS)     =StrIntS_Sloc(nS)/TayVIntS(nS)
         else
            VpIntS_Sloc(nS)  =VpIntN_Sloc(nS)
            dVpIntS_Sloc(nS) =dVpIntN_Sloc(nS)
            ddVpIntS_Sloc(nS)=ddVpIntN_Sloc(nS)
            LFIntS_Sloc(nS)  =LFIntN_Sloc(nS)
            RstrIntS_Sloc(nS)=RstrIntN_Sloc(nS)
            AstrIntS_Sloc(nS)=AstrIntN_Sloc(nS)
            StrIntS_Sloc(nS) =StrIntN_Sloc(nS)
            V2IntS(nS)       =V2IntN(nS)
            Bs2IntS_Sloc(nS) =Bs2IntN_Sloc(nS)
            BspIntS(nS)      =BspIntN(nS)
            BspdIntS_Sloc(nS)=BspdIntN_Sloc(nS)
            BpsdIntS_Sloc(nS)=BpsdIntN_Sloc(nS)
            VpRIntS_Sloc(nS) =VpRIntN_Sloc(nS)
            TayIntS_Sloc(nS) =TayIntN_Sloc(nS)
            TayRIntS(nS)     =TayRIntN(nS)
            TayVIntS(nS)     =TayVIntN(nS)
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
            call getAStr(BszB,BszLMr,2,2,l_max,r_ICB,r_CMB,n_r_max,rB,PlmS(:,:,nS))
            call getAStr(BpzB,BpzLMr,2,2,l_max,r_ICB,r_CMB,n_r_max,rB,PlmS(:,:,nS))
            call getAStr(BzpdB,BzpdLMr,2,2,l_max,r_ICB,r_CMB,n_r_max,rB,PlmS(:,:,nS))
            call getAStr(BpzdB,BpzdLMr,2,2,l_max,r_ICB,r_CMB,n_r_max,rB,PlmS(:,:,nS))

            TauBS_Sloc(nS)  =-(BpzB(2)+sZ(nS)/zMin*BspB(2))
            dTauBS_Sloc(nS) =-(BpzdB(2)+BzpdB(2) + sZ(nS)/zMin*(BspdB(2)+BpsdB(2)))
            dTTauBS_Sloc(nS)=-(BszB(2)+sZ(nS)/zMin*Bs2B(2))
            TauBN_Sloc(nS)  =  BpzB(1)+sZ(nS)/zMax*BspB(1)
            dTauBN_Sloc(nS) =  BpzdB(1)+BzpdB(1) + sZ(nS)/zMax*(BspdB(1)+BpsdB(1))
            dTTauBN_Sloc(nS)=  BszB(1)+sZ(nS)/zMax*Bs2B(1)

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
            call getAStr(BszB,BszLMr,2,2,l_max,r_ICB,r_CMB,n_r_max,rB,PlmS(1,nZmax,nS))
            call getAStr(BpzB,BpzLMr,2,2,l_max,r_ICB,r_CMB,n_r_max,rB,PlmS(1,nZmax,nS))
            call getAStr(BzpdB,BzpdLMr,2,2,l_max,r_ICB,r_CMB,n_r_max, &
                 &       rB,PlmS(1,nZmax,nS))
            call getAStr(BpzdB,BpzdLMr,2,2,l_max,r_ICB,r_CMB,n_r_max, &
                 &       rB,PlmS(1,nZmax,nS))

            TauBS_Sloc(nS)  =TauBS_Sloc(nS)  +(BpzB(2)+sZ(nS)/zMin*BspB(2))
            dTauBS_Sloc(nS) =dTauBS_Sloc(nS) +(BpzdB(2)+BzpdB(2) +              &
            &                         sZ(nS)/zMin*(BspdB(2)+BpsdB(2)))
            dTTauBS_Sloc(nS)=dTTauBS_Sloc(nS)+(BszB(2)+sZ(nS)/zMin*Bs2B(2))
            TauBN_Sloc(nS)  =TauBN_Sloc(nS)  -(BpzB(1) +sZ(nS)/zMax*BspB(1))
            dTauBN_Sloc(nS) =dTauBN_Sloc(nS) -(BpzdB(1)+BzpdB(1) +              &
            &                         sZ(nS)/zMax*(BspdB(1)+BpsdB(1)))
            dTTauBN_Sloc(nS)=dTTauBN_Sloc(nS)-(BszB(1)+sZ(nS)/zMax*Bs2B(1))

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
            call getAStr(BpzB,BpzLMr,2,2,l_max,r_ICB,r_CMB,n_r_max,rB,PlmS(:,:,nS))
            call getAStr(BzpdB,BzpdLMr,2,2,l_max,r_ICB,r_CMB,n_r_max,rB,PlmS(:,:,nS))
            call getAStr(BpzdB,BpzdLMr,2,2,l_max,r_ICB,r_CMB,n_r_max,rB,PlmS(:,:,nS))

            TauBS_Sloc(nS)  = BpzB(1)+sZ(nS)/zMax*BspB(1) - BpzB(2)-sZ(nS)/zMin*BspB(2)
            dTauBS_Sloc(nS) = BpzdB(1)+BzpdB(1) + sZ(nS)/zMax*(BspdB(1)+BpsdB(1)) - &
            &                 BpzdB(2)-BzpdB(2) - sZ(nS)/zMin*(BspdB(2)+BpsdB(2))
            dTTauBS_Sloc(nS)= BszB(1)+sZ(nS)/zMax*Bs2B(1) - BszB(2)-sZ(nS)/zMin*Bs2B(2)
            TauBN_Sloc(nS)  =TauBS_Sloc(nS)
            dTauBN_Sloc(nS) =dTauBS_Sloc(nS)
            dTTauBN_Sloc(nS)=dTTauBS_Sloc(nS)

         end if

      end do outer ! Loop over s
      ! Integration finished

      !--- Integrate Geostrophic azumithal flow energy and Taylor measure:
      VgRMS  =0.0_cp
      TayRMS =0.0_cp
      TayRRMS=0.0_cp
      TayVRMS=0.0_cp
      !--- Integrals changed to not represent the different volumes
      !    represented by each cylinder on Nov. 2 2007:
      do nS=nSstart,nSstop
         ! Old form used again for VgRMS for make is comparable with
         ! VpRMS below.
         VgRMS=VgRMS + two*pi*h(nS)*sZ(nS)*dsZ * VpIntN_Sloc(nS)*VpIntN_Sloc(nS)
         TayRMS =TayRMS +dsZ*abs(TayIntN_Sloc(nS))
         TayRRMS=TayRRMS+dsZ*abs(TayRIntN(nS))
         TayVRMS=TayVRMS+dsZ*abs(TayVIntN(nS))
         if ( nS <= nSI ) then
            VgRMS=VgRMS + two*pi*h(nS)*sZ(nS)*dsZ * VpIntS_Sloc(nS)*VpIntS_Sloc(nS)
            TayRMS =TayRMS +dsZ*abs(TayIntS_Sloc(nS))
            TayRRMS=TayRRMS+dsZ*abs(TayRIntS(nS))
            TayVRMS=TayVRMS+dsZ*abs(TayVIntS(nS))
         end if
      end do

      !--- s-derivatives:
      !------ Create arrays to be differentiated:
      do nS=nSstart,nSstop
         SVpIntN_Sloc(nS) =VpIntN_Sloc(nS)/sZ(nS)
         SBspIntN_Sloc(nS)=h(nS)*sZ(nS)*sZ(nS)*BspIntN(nS)
         SBs2IntN_Sloc(nS)=h(nS)*sZ(nS)**3*Bs2IntN_Sloc(nS)
         if ( sZ(nS) < r_ICB ) then ! inside TC
            SVpIntS_Sloc(nS) =VpIntS_Sloc(nS)/sZ(nS)
            SBspIntS_Sloc(nS)=h(nS)*sZ(nS)*sZ(nS)*BspIntS(nS)
            SBs2IntS_Sloc(nS)=h(nS)*sZ(nS)**3*Bs2IntS_Sloc(nS)
         end if
      end do

      !------
      !-- Here starts the communication stack
      !------

#ifdef WITH_MPI
      call MPI_Reduce(VgRMS,global_sum,1,MPI_DEF_REAL,MPI_SUM,0, &
           &          MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) VgRMS = global_sum
      call MPI_Reduce(TayRMS,global_sum,1,MPI_DEF_REAL,MPI_SUM,0, &
           &          MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) TayRMS = global_sum
      call MPI_Reduce(TayRRMS,global_sum,1,MPI_DEF_REAL,MPI_SUM,0, &
           &          MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) TayRRMS = global_sum
      call MPI_Reduce(TayVRMS,global_sum,1,MPI_DEF_REAL,MPI_SUM,0, &
           &          MPI_COMM_WORLD,ierr)
      if ( rank == 0 ) TayVRMS = global_sum

      sendcount  = (nSstop-nSstart+1)*nZmaxA
      do i=0,n_procs-1
         recvcounts(i) = cyl_balance(i)%n_per_rank*nZmaxA
      end do
      displs(0)=0
      do i=1,n_procs-1
         displs(i)=displs(i-1)+recvcounts(i-1)
      end do
      call MPI_GatherV(zZ_Sloc, sendcount, MPI_DEF_REAL,        &
           &           zZ, recvcounts, displs, MPI_DEF_REAL,    &
           &           0, MPI_COMM_WORLD, ierr)

      if ( lTOZwrite ) then
         call MPI_GatherV(VpS_Sloc, sendcount, MPI_DEF_REAL,       &
              &           VpS, recvcounts, displs, MPI_DEF_REAL,   &
              &           0, MPI_COMM_WORLD, ierr)
         call MPI_GatherV(dVpS_Sloc, sendcount, MPI_DEF_REAL,      &
              &           dVpS, recvcounts, displs, MPI_DEF_REAL,  &
              &           0, MPI_COMM_WORLD, ierr)
         call MPI_GatherV(RstrS_Sloc, sendcount, MPI_DEF_REAL,     &
              &           RstrS, recvcounts, displs, MPI_DEF_REAL, &
              &           0, MPI_COMM_WORLD, ierr)
         call MPI_GatherV(AstrS_Sloc, sendcount, MPI_DEF_REAL,     &
              &           AstrS, recvcounts, displs, MPI_DEF_REAL, &
              &           0, MPI_COMM_WORLD, ierr)
         call MPI_GatherV(StrS_Sloc, sendcount, MPI_DEF_REAL,      &
              &           StrS, recvcounts, displs, MPI_DEF_REAL,  &
              &           0, MPI_COMM_WORLD, ierr)
         call MPI_GatherV(CorS_Sloc, sendcount, MPI_DEF_REAL,      &
              &           CorS, recvcounts, displs, MPI_DEF_REAL,  &
              &           0, MPI_COMM_WORLD, ierr)
         call MPI_GatherV(LFS_Sloc, sendcount, MPI_DEF_REAL,       &
              &           LFS, recvcounts, displs, MPI_DEF_REAL,   &
              &           0, MPI_COMM_WORLD, ierr)
      end if

      if ( l_TOZave .and. nTOsets > 1 ) then
         call MPI_GatherV(VpM_Sloc, sendcount, MPI_OUT_REAL,       &
              &           VpM, recvcounts, displs, MPI_OUT_REAL,   &
              &           0, MPI_COMM_WORLD, ierr)
         call MPI_GatherV(dVpM_Sloc, sendcount, MPI_OUT_REAL,      &
              &           dVpM, recvcounts, displs, MPI_OUT_REAL,  &
              &           0, MPI_COMM_WORLD, ierr)
         call MPI_GatherV(RstrM_Sloc, sendcount, MPI_OUT_REAL,     &
              &           RstrM, recvcounts, displs, MPI_OUT_REAL, &
              &           0, MPI_COMM_WORLD, ierr)
         call MPI_GatherV(AstrM_Sloc, sendcount, MPI_OUT_REAL,     &
              &           AstrM, recvcounts, displs, MPI_OUT_REAL, &
              &           0, MPI_COMM_WORLD, ierr)
         call MPI_GatherV(StrM_Sloc, sendcount, MPI_OUT_REAL,      &
              &           StrM, recvcounts, displs, MPI_OUT_REAL,  &
              &           0, MPI_COMM_WORLD, ierr)
         call MPI_GatherV(CorM_Sloc, sendcount, MPI_OUT_REAL,      &
              &           CorM, recvcounts, displs, MPI_OUT_REAL,  &
              &           0, MPI_COMM_WORLD, ierr)
         call MPI_GatherV(LFM_Sloc, sendcount, MPI_OUT_REAL,       &
              &           LFM, recvcounts, displs, MPI_OUT_REAL,   &
              &           0, MPI_COMM_WORLD, ierr)
         call MPI_GatherV(CLM_Sloc, sendcount, MPI_OUT_REAL,       &
              &           CLM, recvcounts, displs, MPI_OUT_REAL,   &
              &           0, MPI_COMM_WORLD, ierr)
      end if

      sendcount  = (nSstop-nSstart+1)
      do i=0,n_procs-1
         recvcounts(i) = cyl_balance(i)%n_per_rank
      end do
      displs(0)=0
      do i=1,n_procs-1
         displs(i)=displs(i-1)+recvcounts(i-1)
      end do

      call MPI_GatherV(nZmaxS_Sloc, sendcount, MPI_INTEGER,         &
           &           nZmaxS, recvcounts, displs, MPI_INTEGER,     &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(VpIntN_Sloc, sendcount, MPI_DEF_REAL,        &
           &           VpIntN, recvcounts, displs, MPI_DEF_REAL,    &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(VpIntS_Sloc, sendcount, MPI_DEF_REAL,        &
           &           VpIntS, recvcounts, displs, MPI_DEF_REAL,    &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(SVpIntN_Sloc, sendcount, MPI_DEF_REAL,       &
           &           SVpIntN, recvcounts, displs, MPI_DEF_REAL,   &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(SBspIntN_Sloc, sendcount, MPI_DEF_REAL,      &
           &           SBspIntN, recvcounts, displs, MPI_DEF_REAL,  &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(SBs2IntN_Sloc, sendcount, MPI_DEF_REAL,      &
           &           SBs2IntN, recvcounts, displs, MPI_DEF_REAL,  &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(SVpIntS_Sloc, sendcount, MPI_DEF_REAL,       &
           &           SVpIntS, recvcounts, displs, MPI_DEF_REAL,   &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(SBspIntS_Sloc, sendcount, MPI_DEF_REAL,      &
           &           SBspIntS, recvcounts, displs, MPI_DEF_REAL,  &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(SBs2IntS_Sloc, sendcount, MPI_DEF_REAL,      &
           &           SBs2IntS, recvcounts, displs, MPI_DEF_REAL,  &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(TauBN_Sloc, sendcount, MPI_DEF_REAL,         &
           &           TauBN, recvcounts, displs, MPI_DEF_REAL,     &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(TauBS_Sloc, sendcount, MPI_DEF_REAL,         &
           &           TauBS, recvcounts, displs, MPI_DEF_REAL,     &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(dTauBN_Sloc, sendcount, MPI_DEF_REAL,        &
           &           dTauBN, recvcounts, displs, MPI_DEF_REAL,    &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(dTauBS_Sloc, sendcount, MPI_DEF_REAL,        &
           &           dTauBS, recvcounts, displs, MPI_DEF_REAL,    &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(dTTauBN_Sloc, sendcount, MPI_DEF_REAL,       &
           &           dTTauBN, recvcounts, displs, MPI_DEF_REAL,   &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(dTTauBS_Sloc, sendcount, MPI_DEF_REAL,       &
           &           dTTauBS, recvcounts, displs, MPI_DEF_REAL,   &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(Bs2IntS_Sloc, sendcount, MPI_DEF_REAL,       &
           &           Bs2IntS, recvcounts, displs, MPI_DEF_REAL,   &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(Bs2IntN_Sloc, sendcount, MPI_DEF_REAL,       &
           &           Bs2IntN, recvcounts, displs, MPI_DEF_REAL,   &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(dVpIntN_Sloc, sendcount, MPI_DEF_REAL,       &
           &           dVpIntN, recvcounts, displs, MPI_DEF_REAL,   &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(dVpIntS_Sloc, sendcount, MPI_DEF_REAL,       &
           &           dVpIntS, recvcounts, displs, MPI_DEF_REAL,   &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(ddVpIntN_Sloc, sendcount, MPI_DEF_REAL,      &
           &           ddVpIntN, recvcounts, displs, MPI_DEF_REAL,  &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(ddVpIntS_Sloc, sendcount, MPI_DEF_REAL,      &
           &           ddVpIntS, recvcounts, displs, MPI_DEF_REAL,  &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(VpRIntN_Sloc, sendcount, MPI_DEF_REAL,       &
           &           VpRIntN, recvcounts, displs, MPI_DEF_REAL,   &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(VpRIntS_Sloc, sendcount, MPI_DEF_REAL,       &
           &           VpRIntS, recvcounts, displs, MPI_DEF_REAL,   &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(LFIntN_Sloc, sendcount, MPI_DEF_REAL,        &
           &           LFIntN, recvcounts, displs, MPI_DEF_REAL,    &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(LFIntS_Sloc, sendcount, MPI_DEF_REAL,        &
           &           LFIntS, recvcounts, displs, MPI_DEF_REAL,    &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(RstrIntN_Sloc, sendcount, MPI_DEF_REAL,      &
           &           RstrIntN, recvcounts, displs, MPI_DEF_REAL,  &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(RstrIntS_Sloc, sendcount, MPI_DEF_REAL,      &
           &           RstrIntS, recvcounts, displs, MPI_DEF_REAL,  &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(AstrIntN_Sloc, sendcount, MPI_DEF_REAL,      &
           &           AstrIntN, recvcounts, displs, MPI_DEF_REAL,  &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(AstrIntS_Sloc, sendcount, MPI_DEF_REAL,      &
           &           AstrIntS, recvcounts, displs, MPI_DEF_REAL,  &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(StrIntN_Sloc, sendcount, MPI_DEF_REAL,       &
           &           StrIntN, recvcounts, displs, MPI_DEF_REAL,   &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(StrIntS_Sloc, sendcount, MPI_DEF_REAL,       &
           &           StrIntS, recvcounts, displs, MPI_DEF_REAL,   &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(TayIntN_Sloc, sendcount, MPI_DEF_REAL,       &
           &           TayIntN, recvcounts, displs, MPI_DEF_REAL,   &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(TayIntS_Sloc, sendcount, MPI_DEF_REAL,       &
           &           TayIntS, recvcounts, displs, MPI_DEF_REAL,   &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(BspdIntN_Sloc, sendcount, MPI_DEF_REAL,      &
           &           BspdIntN, recvcounts, displs, MPI_DEF_REAL,  &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(BspdIntS_Sloc, sendcount, MPI_DEF_REAL,      &
           &           BspdIntS, recvcounts, displs, MPI_DEF_REAL,  &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(BpsdIntN_Sloc, sendcount, MPI_DEF_REAL,      &
           &           BpsdIntN, recvcounts, displs, MPI_DEF_REAL,  &
           &           0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(BpsdIntS_Sloc, sendcount, MPI_DEF_REAL,      &
           &           BpsdIntS, recvcounts, displs, MPI_DEF_REAL,  &
           &           0, MPI_COMM_WORLD, ierr)
#else
      nZmaxS(:)  =nZmaxS_Sloc(:)
      zZ(:,:)    =zZ_Sloc(:,:)
      VpIntN(:)  =VpIntN_Sloc(:)
      VpIntS(:)  =VpIntS_Sloc(:)
      SVpIntN(:) =SVpIntN_Sloc(:)
      SBspIntN(:)=SBspIntN_Sloc(:)
      SBs2IntN(:)=SBs2IntN_Sloc(:)
      SVpIntS(:) =SVpIntS_Sloc(:)
      SBspIntS(:)=SBspIntS_Sloc(:)
      SBs2IntS(:)=SBs2IntS_Sloc(:)
      TauBN(:)   =TauBN_Sloc(:)
      TauBS(:)   =TauBS_Sloc(:)
      dTauBN(:)  =dTauBN_Sloc(:)
      dTauBS(:)  =dTauBS_Sloc(:)
      dTTauBN(:) =dTTauBN_Sloc(:)
      dTTauBS(:) =dTTauBS_Sloc(:)
      Bs2IntS(:) =Bs2IntS_Sloc(:)
      Bs2IntN(:) =Bs2IntN_Sloc(:)
      dVpIntN(:) =dVpIntN_Sloc(:)
      dVpIntS(:) =dVpIntS_Sloc(:)
      ddVpIntN(:)=ddVpIntN_Sloc(:)
      ddVpIntS(:)=ddVpIntS_Sloc(:)
      VpRIntN(:) =VpRIntN_Sloc(:)
      VpRIntS(:) =VpRIntS_Sloc(:)
      LFIntN(:)  =LFIntN_Sloc(:)
      LFIntS(:)  =LFIntS_Sloc(:)
      RstrIntN(:)=RstrIntN_Sloc(:)
      RstrIntS(:)=RstrIntS_Sloc(:)
      AstrIntN(:)=AstrIntN_Sloc(:)
      AstrIntS(:)=AstrIntS_Sloc(:)
      StrIntN(:) =StrIntN_Sloc(:)
      StrIntS(:) =StrIntS_Sloc(:)
      TayIntN(:) =TayIntN_Sloc(:)
      TayIntS(:) =TayIntS_Sloc(:)
      BspdIntN(:)=BspdIntN_Sloc(:)
      BspdIntS(:)=BspdIntS_Sloc(:)
      BpsdIntN(:)=BpsdIntN_Sloc(:)
      BpsdIntS(:)=BpsdIntS_Sloc(:)

      if ( lTOZwrite ) then
         VpS(:,:)  =VpS_Sloc(:,:)
         dVpS(:,:) =dVpS_Sloc(:,:)
         RstrS(:,:)=RstrS_Sloc(:,:)
         AstrS(:,:)=AstrS_Sloc(:,:)
         StrS(:,:) =StrS_Sloc(:,:)
         CorS(:,:) =CorS_Sloc(:,:)
         LFS(:,:)  =LFS_Sloc(:,:)
      end if

      if ( l_TOZave .and. nTOsets > 1 ) then
         VpM(:,:)  =VpM_Sloc(:,:)
         dVpM(:,:) =dVpM_Sloc(:,:)
         RstrM(:,:)=RstrM_Sloc(:,:)
         AstrM(:,:)=AstrM_Sloc(:,:)
         StrM(:,:) =StrM_Sloc(:,:)
         CorM(:,:) =CorM_Sloc(:,:)
         LFM(:,:)  =LFM_Sloc(:,:)
         CLM(:,:)  =CLM_Sloc(:,:)
      end if
#endif


      if ( rank == 0 ) then
         !-- Outputs
         if ( nTOsets == 1 ) nTOZfile=0
         if ( lTOZwrite ) then
            nTOZfile=nTOZfile+1
            write(string, *) nTOZfile
            fileName='TOZ_'//trim(adjustl(string))//'.'//tag
            open(newunit=n_toz_file, file=fileName, form='unformatted', &
            &    status='unknown')
            write(n_toz_file) real(time,kind=outp), real(nSmax,kind=outp), &
            &                 real(omega_ic,kind=outp),                    &
            &                 real(omega_ma,kind=outp)
            write(n_toz_file) (real(sZ(nS),kind=outp),nS=1,nSmax)
         end if
         if ( nTOsets > 1 .and. l_TOZave ) then
            fileName='TOZ_ave.'//tag
            open(newunit=n_tozm_file,file=fileName, form='unformatted', &
            &    status='unknown')
            write(n_tozm_file) real(nSmax,kind=outp),real(omega_ic,kind=outp), &
            &                  real(omega_ma,kind=outp)
            write(n_tozm_file) (real(sZ(nS),kind=outp),nS=1,nSmax)
         end if

         do nS=1,nSmax

            if ( sZ(nS) < r_ICB ) then
               lTC=.true.
            else
               lTC=.false.
            end if

            nZmax=nZmaxS(nS)
            if ( lTC ) then
               nZmaxNS=2*nZmax
               do nZ=1,nZmax
                  zALL(nZ)=zZ(nZ,nS)
                  zALL(nZmaxNS-nZ+1)=-zZ(nZ,nS)
               end do
            else
               nZmaxNS=nZmax
               do nZ=1,nZmax
                  zALL(nZ)=zZ(nZ,nS)
               end do
            end if

            if ( l_TOZave .and. nTOsets > 1 ) then
               write(n_tozm_file) real(nZmaxNS,kind=outp)
               write(n_tozm_file) (real(zALL(nZ),kind=outp),nZ=1,nZmaxNS), &
               &                  (VpM(nZ,nS)/timeAve  ,nZ=1,nZmaxNS),     &
               &                  (dVpM(nZ,nS)/timeAve ,nZ=1,nZmaxNS),     &
               &                  (RstrM(nZ,nS)/timeAve,nZ=1,nZmaxNS),     &
               &                  (AstrM(nZ,nS)/timeAve,nZ=1,nZmaxNS),     &
               &                  (LFM(nZ,nS)/timeAve  ,nZ=1,nZmaxNS),     &
               &                  (StrM(nZ,nS)/timeAve ,nZ=1,nZmaxNS),     &
               &                  (CorM(nZ,nS)/timeAve ,nZ=1,nZmaxNS),     &
               &                  (CLM(nZ,nS)/timeAve  ,nZ=1,nZmaxNS)
            end if
            if ( lTOZwrite ) then
               write(n_toz_file) real(nZmaxNS)
               write(n_toz_file) (real(zALL(nZ),kind=outp) ,nZ=1,nZmaxNS),       &
               &                 (real(VpS(nZ,nS),kind=outp)  ,nZ=1,nZmaxNS),    &
               &                 (real(dVpS(nZ,nS),kind=outp) ,nZ=1,nZmaxNS),    &
               &                 (real(RstrS(nZ,nS),kind=outp),nZ=1,nZmaxNS),    &
               &                 (real(AstrS(nZ,nS),kind=outp),nZ=1,nZmaxNS),    &
               &                 (real(LFfac*LFS(nZ,nS),kind=outp),nZ=1,nZmaxNS),&
               &                 (real(StrS(nZ,nS),kind=outp) ,nZ=1,nZmaxNS),    &
               &                 (real(CorS(nZ,nS),kind=outp) ,nZ=1,nZmaxNS)
            end if
         end do

         if ( lTOZwrite ) close(n_toz_file)
         if ( l_TOZave .and. nTOsets > 1 ) close (n_tozm_file)

         if ( lStopRun ) call abortRun('Stop run in out_TO')

         !------ Create Derivatives:
         f1=one/(12.0_cp*dsZ)
         f2=f1/dsZ
         do nS=3,nSmax-2
            dSVpIntN =f1*(     SVpIntN(nS-2)-8.0_cp*SVpIntN(nS-1) +        &
            &                        8.0_cp*SVpIntN(nS+1)-     SVpIntN(nS+2) )
            d2SVpIntN=f2*(    -SVpIntN(nS-2)+16.0_cp*SVpIntN(nS-1) -       &
            &            30.0_cp*SVpIntN(nS)+16.0_cp*SVpIntN(nS+1)-SVpIntN(nS+2))
            dSBspIntN=f1*(     SBspIntN(nS-2)-8.0_cp*SBspIntN(nS-1) +      &
            &                        8.0_cp*SBspIntN(nS+1)-     SBspIntN(nS+2) )
            dSBs2IntN=f1*(     SBs2IntN(nS-2)-8.0_cp*SBs2IntN(nS-1) +      &
            &                        8.0_cp*SBs2IntN(nS+1)-     SBs2IntN(nS+2) )
            TauN(nS)  =Oh(nS)*(Os2(nS)*dSBspIntN+TauBN(nS))
            dTTauN(nS)=sZ(nS)*d2SVpIntN*Bs2IntN(nS) +                    &
            &                    Oh(nS)*(Os2(nS)*dSVpIntN*dSBs2IntN +    &
            &                             sZ(nS)*dSVpIntN*dTTauBN(nS) )
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
               &                               8.0_cp*SVpIntS(nS+1)-SVpIntS(nS+2) )
               d2SVpIntS=f2*(-SVpIntS(nS-2)+16.0_cp*SVpIntS(nS-1) -        &
               &               30.0_cp*SVpIntS(nS)+16.0_cp*SVpIntS(nS+1)-SVpIntS(nS+2))
               dSBspIntS=f1*(SBspIntS(nS-2)-8.0_cp*SBspIntS(nS-1) +        &
               &                           8.0_cp*SBspIntS(nS+1)-SBspIntS(nS+2) )
               dSBs2IntS=f1*(SBs2IntS(nS-2)-8.0_cp*SBs2IntS(nS-1) +        &
               &                           8.0_cp*SBs2IntS(nS+1)-SBs2IntS(nS+2) )
               TauS(nS) =Oh(nS)*(Os2(nS)*dSBspIntS+TauBS(nS))
               dTTauS(nS)=sZ(nS)*d2SVpIntS*Bs2IntS(nS) +                 &
               &                    Oh(nS)*(Os2(nS)*dSVpIntS*dSBs2IntS + &
               &                             sZ(nS)*dSVpIntS*dTTauBS(nS) )
            else
               TauS(nS)  =TauN(nS)
               dTTauS(nS)=dTTauN(nS)
            end if
         end do

         !--- Output of z-integral:
         open(newunit=nOutFile, file=TOfileNhs, status='unknown',    &
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

         open(newunit=nOutFile, file=TOfileShs, status='unknown',       &
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

            !--- Transform back to radial space:
            call rscheme_oc%costf1(dzVpLMr,l_max+1,1,l_max+1,workA)
            call rscheme_oc%costf1(dzdVpLMr,l_max+1,1,l_max+1,workA)
            call rscheme_oc%costf1(dzRstrLMr,l_max+1,1,l_max+1,workA)
            call rscheme_oc%costf1(dzAstrLMr,l_max+1,1,l_max+1,workA)
            call rscheme_oc%costf1(dzStrLMr,l_max+1,1,l_max+1,workA)
            call rscheme_oc%costf1(dzLFLMr,l_max+1,1,l_max+1,workA)
            call rscheme_oc%costf1(dzCorLMr,l_max+1,1,l_max+1,workA)

            !--- Open output file
            nFields=7
            if ( lTOmov ) then

               open(newunit=nOutFile, file=movFile, status='unknown',  &
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
               open(newunit=nOutFile2, file=tayFile, form='unformatted',  &
               &    status='unknown', position='append')
               if ( nTOrmsSets == 1 ) then
                  write(nOutFile2) real(n_r_max,kind=outp)
                  write(nOutFile2) (real(r(nR),kind=outp),nR=1,n_r_max)
               end if

            end if

            do nOut=1,nFields ! Loop over four output fields

               do nR=1,n_r_max ! Loop over radial points
                  rS=r(nR)

                  !------ Convert from lm to theta block:
                  if ( nOut == 1 ) then
                     call get_PAS(dzVpLMr(:,nR),outBlock,rS)
                  else if ( nOut == 2 ) then
                     call get_PAS(dzRstrLMr(:,nR),outBlock,rS)
                  else if ( nOut == 3 ) then
                     call get_PAS(dzAstrLMr(:,nR),outBlock,rS)
                  else if ( nOut == 4 ) then
                     call get_PAS(dzStrLMr(:,nR),outBlock,rS)
                  else if ( nOut == 5 ) then
                     call get_PAS(dzLFLMr(:,nR),outBlock,rS)
                  else if ( nOut == 6 ) then
                     call get_PAS(dzCorLMr(:,nR),outBlock,rS)
                  else if ( nOut == 7 ) then
                     call get_PAS(dzdVpLMr(:,nR),outBlock,rS)
                  end if

                  !------ Storage of field in fout for theta block,
                  !       integration and Max/Min values
                  do nTheta=1,n_theta_max
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
                        fOut(nPos)=vSF*outBlock(nTheta)
                        VpR(nR)=VpR(nR) + gauss(nThetaNHS)*fOut(nPos)/sS
                     else if ( nOut == 2 ) then
                        !--------------- Reynolds force:
                        fOut(nPos)=fSF*outBlock(nTheta)
                        RstrR(nR)=RstrR(nR) + gauss(nThetaNHS)*fOut(nPos)/sS
                     else if ( nOut == 3 ) then
                        !--------------- Advective force:
                        fOut(nPos)=fSF*outBlock(nTheta)
                        AstrR(nR)=AstrR(nR) + gauss(nThetaNHS)*fOut(nPos)/sS
                     else if ( nOut == 4 ) then
                        !--------------- Viscous force:
                        fOut(nPos)=fSF*outBlock(nTheta)
                        StrR(nR)=StrR(nR) + gauss(nThetaNHS)*fOut(nPos)/sS
                     else if ( nOut == 5 ) then
                        !--------------- Lorentz force:
                        fOut(nPos)=fSF*LFfac*outBlock(nTheta)
                        LFR(nR)=LFR(nR) + gauss(nThetaNHS)*fOut(nPos)/sS
                        LFABSR(nR)=LFABSR(nR) + gauss(nThetaNHS)*abs(fOut(nPos))/sS
                     else if ( nOut == 6 ) then
                        !--------------- Corriolis force:
                        fOut(nPos)=fSF*outBlock(nTheta)
                        CorR(nR)=CorR(nR) + gauss(nThetaNHS)*fOut(nPos)/sS
                     else if ( nOut == 7 ) then
                        !--------------- dtVp:
                        fOut(nPos)=vSF*outBlock(nTheta)
                        dVpR(nR)  =dVpR(nR) + gauss(nThetaNHS)*fOut(nPos)/sS
                     end if

                  end do ! Loop over thetas in block

               end do ! Loop over R

               !------ Output of stress contributions for one radial grid point:
               !       Write all fields into movie style file
               if ( lTOmov ) write(nOutFile) (real(fOut(nPos),kind=outp),nPos=1,nFieldSize)

            end do ! Loop over output functions

            if ( l_save_out ) then
               open(newunit=n_log_file, file=log_file, status='unknown', &
               &    position='append')
            end if
            if ( lTOmov ) then
               close(nOutFile)
               call logWrite(' ')
               write(message,'(1p,A,I8,A,ES16.6,I8)')              &
                    & "! WRITING TO MOVIE FRAME NO ",nTOmovSets,   &
                    & "       at time/step",time*tScale,n_time_step
               call logWrite(message)
            end if
            if ( l_save_out ) close(n_log_file)

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
                  TayR(nR)=abs(LFR(nR))/(fac*LFABSR(nR))
               else
                  TayR(nR)=0.0_cp
               end if
               !              TayRMSR(nR)=rS*rS*TayR(nR)
               TayRMSR(nR)=TayR(nR)
            end do

            !--- Now perform the radial integral: ( not tested )
            TaySRMS=rInt_R(TayRMSR,r,rscheme_oc)
            !--- And finally calculate the mean value, the factor 4*pi comes from
            !    the fact that the surface integral has already been cared for
            !    NOTE: Integral for RMS Taylorisation changed to not respect the
            !    different volumes represented by each shell on 2 Nov 2007.
            !          fac    =four*pi/vol_oc
            fac    =one/(r_cmb-r_icb) ! =1
            TaySRMS=fac*TaySRMS

            write(nOutFile2)                                                  &
                 &          real(time,kind=outp),real(VpRMS**2,kind=outp),    &
                 &          real(VgRMS**2,kind=outp),real(TayRMS,kind=outp),  &
                 &          real(TaySRMS,kind=outp),real(TayRRMS,kind=outp),  &
                 &          real(TayVRMS,kind=outp),real(eKin,kind=outp),     &! 3
                 &          (real(VpR(nR),kind=outp)  ,nR=1,n_r_max),         &! 4
                 &          (real(dVpR(nR),kind=outp) ,nR=1,n_r_max),         &! 5
                 &          (real(RstrR(nR),kind=outp),nR=1,n_r_max),         &! 6
                 &          (real(AstrR(nR),kind=outp),nR=1,n_r_max),         &! 7
                 &          (real(LFR(nR),kind=outp)  ,nR=1,n_r_max),         &! 8
                 &          (real(StrR(nR),kind=outp) ,nR=1,n_r_max),         &! 9
                 &          (real(CorR(nR),kind=outp) ,nR=1,n_r_max),         &! 10
                 &          (real(TayR(nR),kind=outp) ,nR=1,n_r_max)           ! 11
            close(nOutFile2)

         end if



      end if ! Rank 0

      timeLast=time

      lTOZwrite=.false.

      if ( lVerbose ) write(output_unit,*) '! End of outTO!'

   end subroutine outTO
!----------------------------------------------------------------------------
   subroutine outTO_allgather_Rloc
      !
      ! MPI communicators for TO outputs
      !

#ifdef WITH_MPI
      !-- Local variables:
      integer :: sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1)
      integer :: p,ierr

      sendcount  = nR_per_rank*(l_max+1)
      do p=0,n_procs-1
         recvcounts(p)=radial_balance(p)%n_per_rank*(l_max+1)
      end do
      displs(0)=0
      do p=1,n_procs-1
         displs(p)=displs(p-1)+recvcounts(p-1)
      end do
      call MPI_AllGatherV(dzStrLMr_Rloc,sendcount,MPI_DEF_REAL,    &
           &              dzStrLMr,recvcounts,displs,MPI_DEF_REAL, &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllGatherV(dzRstrLMr_Rloc,sendcount,MPI_DEF_REAL,   &
           &              dzRstrLMr,recvcounts,displs,MPI_DEF_REAL,&
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllGatherV(dzAstrLMr_Rloc,sendcount,MPI_DEF_REAL,   &
           &              dzAstrLMr,recvcounts,displs,MPI_DEF_REAL,&
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllGatherV(dzCorLMr_Rloc,sendcount,MPI_DEF_REAL,    &
           &              dzCorLMr,recvcounts,displs,MPI_DEF_REAL, &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllGatherV(dzLFLMr_Rloc,sendcount,MPI_DEF_REAL,     &
           &              dzLFLMr,recvcounts,displs,MPI_DEF_REAL,  &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllGatherV(dzdVpLMr_Rloc,sendcount,MPI_DEF_REAL,    &
           &              dzdVpLMr,recvcounts,displs,MPI_DEF_REAL, &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllGatherV(dzddVpLMr_Rloc,sendcount,MPI_DEF_REAL,   &
           &              dzddVpLMr,recvcounts,displs,MPI_DEF_REAL,&
           &              MPI_COMM_WORLD,ierr)

      call MPI_AllGatherV(V2LMr_Rloc,sendcount,MPI_DEF_REAL,       &
           &              V2LMr,recvcounts,displs,MPI_DEF_REAL,    &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllGatherV(Bs2LMr_Rloc,sendcount,MPI_DEF_REAL,      &
           &              Bs2LMr,recvcounts,displs,MPI_DEF_REAL,   &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllGatherV(BszLMr_Rloc,sendcount,MPI_DEF_REAL,      &
           &              BszLMr,recvcounts,displs,MPI_DEF_REAL,   &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllGatherV(BspLMr_Rloc,sendcount,MPI_DEF_REAL,      &
           &              BspLMr,recvcounts,displs,MPI_DEF_REAL,   &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllGatherV(BpzLMr_Rloc,sendcount,MPI_DEF_REAL,      &
           &              BpzLMr,recvcounts,displs,MPI_DEF_REAL,   &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllGatherV(BspdLMr_Rloc,sendcount,MPI_DEF_REAL,     &
           &              BspdLMr,recvcounts,displs,MPI_DEF_REAL,  &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllGatherV(BpsdLMr_Rloc,sendcount,MPI_DEF_REAL,     &
           &              BpsdLMr,recvcounts,displs,MPI_DEF_REAL,  &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllGatherV(BzpdLMr_Rloc,sendcount,MPI_DEF_REAL,     &
           &              BzpdLMr,recvcounts,displs,MPI_DEF_REAL,  &
           &              MPI_COMM_WORLD,ierr)
      call MPI_AllGatherV(BpzdLMr_Rloc,sendcount,MPI_DEF_REAL,     &
           &              BpzdLMr,recvcounts,displs,MPI_DEF_REAL,  &
           &              MPI_COMM_WORLD,ierr)
#else
     dzStrLMr(:,:) =dzStrLMr_Rloc(:,:)
     dzRstrLMr(:,:)=dzRstrLMr_Rloc(:,:)
     dzAstrLMr(:,:)=dzAstrLMr_Rloc(:,:)
     dzCorLMr(:,:) =dzCorLMr_Rloc(:,:)
     dzLFLMr(:,:)  =dzLFLMr_Rloc(:,:)
     dzdVpLMr(:,:) =dzdVpLMr_Rloc(:,:)
     dzddVpLMr(:,:)=dzddVpLMr_Rloc(:,:)
     V2LMr(:,:)    =V2LMr_Rloc(:,:)
     Bs2LMr(:,:)   =Bs2LMr_Rloc(:,:)
     BszLMr(:,:)   =BszLMr_Rloc(:,:)
     BspLMr(:,:)   =BspLMr_Rloc(:,:)
     BpzLMr(:,:)   =BpzLMr_Rloc(:,:)
     BspdLMr(:,:)  =BspdLMr_Rloc(:,:)
     BpsdLMr(:,:)  =BpsdLMr_Rloc(:,:)
     BzpdLMr(:,:)  =BzpdLMr_Rloc(:,:)
     BpzdLMr(:,:)  =BpzdLMr_Rloc(:,:)
#endif

   end subroutine outTO_allgather_Rloc
!----------------------------------------------------------------------------
end module outTO_mod
