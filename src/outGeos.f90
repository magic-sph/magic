module geos_mod
   !
   ! This module is used to calculate the outputs when either l_par = .true.
   !

   use precision_mod
   use parallel_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, lm_max, n_m_max, n_phi_max, nrp,     &
       &                 minc, l_max, m_max, l_axi, n_phi_tot
   use radial_functions, only: r_ICB, r_CMB, rscheme_oc, orho1
   use physical_parameters, only: ra, ek, pr, prmag, radratio
   use num_param, only: tScale
   use blocking, only: lm2l, lm2m, lm2mc, lo_map, st_map, llm, ulm
   use horizontal_data, only: dLh, phi, dPhi
   use logic, only: lVerbose, l_corrMov, l_anel, l_save_out, l_SRIC
   use output_data, only: sDens, zDens, tag, runid
   use constants, only: pi, zero, ci, one, two, three, four,  half
   use communications, only: gather_all_from_lo_to_rank0,gt_OC
   use plms_theta, only: plm_theta
   use fft, only: fft_to_real
   use TO_helpers, only: getPAStr
   use cosine_transform_odd, only: costf_odd_t
   use chebInt_mod

   implicit none

   private

   real(cp) :: timeOld
   real(cp), allocatable :: PlmS(:,:,:), PlmZ(:,:,:)
   real(cp), allocatable :: dPlmS(:,:,:), dPlmZ(:,:,:)
   real(cp), allocatable :: sinTS(:,:)

   complex(cp), allocatable :: wS_global(:,:), dwS_global(:,:), ddwS_global(:,:)
   complex(cp), allocatable :: zS_global(:,:), dzS_global(:,:)

   type(costf_odd_t), allocatable :: chebt_Z(:)
   integer, allocatable :: nZmaxS(:), nZC(:), nZC_Sloc(:), nZ2(:,:)
   real(cp), allocatable :: zZ(:,:), rZS(:,:)
   real(cp), allocatable :: VorOld(:,:,:)
   real(cp), parameter :: eps = 10.0_cp*epsilon(one)

   integer :: n_geos_file, nrp_geos
   integer :: nSmax, nZmaxA, nZmaxAH, nZnorm, nSstart, nSstop
   real(cp) :: zNorm
   type(load), allocatable :: cyl_balance(:)
   character(len=72) :: geos_file

   public :: initialize_geos_mod, getEgeos, finalize_geos_mod

contains

   subroutine initialize_geos_mod(l_geos)
      !
      ! Memory allocation in Egeos
      !
      logical, intent(in) :: l_geos

#ifdef WITH_SHTNS
      nrp_geos=nrp+2 ! One has to include the 2 extra points again
#else
      nrp_geos=nrp
#endif

      ! Number of grid point in s:
      nSmax=n_r_max+int(r_ICB*real(n_r_max,cp))
      nSmax=int(sDens*nSmax) ! scale up with sDens
      ! Reference number of grid point in z covering 
      ! one hemisphere at s=r_ICB.
      ! Maximum number used in z
      zNorm=two*sqrt(1-radratio**2)*r_CMB
      nZnorm=2*int(real(nSmax,cp)*zNorm/r_CMB)
      nZmaxA=nZnorm+51 ! Add safety margin
      nZmaxAH=nZnorm+25

      !-- Distribute over the ranks
      allocate(cyl_balance(0:n_procs-1))
      call getBlocks(cyl_balance, nSmax, n_procs)
      nSstart = cyl_balance(rank)%nStart
      nSstop = cyl_balance(rank)%nStop

      !-- The following global arrays are required in getDVptr
      allocate( wS_global(lm_max,n_r_max), dwS_global(lm_max,n_r_max) )
      allocate( ddwS_global(lm_max,n_r_max), zS_global(lm_max,n_r_max) )
      allocate( dzS_global(lm_max,n_r_max) )
      bytes_allocated = bytes_allocated+ 5*lm_max*n_r_max*SIZEOF_DEF_COMPLEX

      if ( l_geos ) then
         allocate( sinTS(nZmaxAH,nSstart:nSstop) )
         allocate( PlmS(lm_max,nZmaxAH,nSstart:nSstop) )
         allocate( dPlmS(lm_max,nZmaxAH,nSstart:nSstop) )
         allocate( zZ(nZmaxA,nSstart:nSstop) )
         allocate( rZS(nZmaxA,nSstart:nSstop) )
         bytes_allocated = bytes_allocated+ (nZmaxAH*(nSstop-nSstart+1)* &
         &                (1+2*lm_max))*SIZEOF_DEF_REAL + &
         &                2*nZmaxA*(nSstop-nSstart+1)*SIZEOF_DEF_REAL

         allocate( nZmaxS(nSstart:nSstop) )
         allocate( chebt_Z(nSstart:nSstop) )
         bytes_allocated = bytes_allocated+(nSstop-nSstart+1)*SIZEOF_INTEGER

         geos_file='geos.'//tag
         if ( rank == 0 .and. (.not. l_save_out) ) then
            open(newunit=n_geos_file, file=geos_file, status='new')
         end if
      end if

   end subroutine initialize_geos_mod

!----------------------------------------------------------------------------

   subroutine finalize_geos_mod(l_geos)
      !
      ! Memory deallocation
      !
      logical, intent(in) :: l_geos

      deallocate( wS_global, dwS_global, ddwS_global, zS_global, dzS_global)

      if ( l_geos ) then
         deallocate( sinTS, PlmS, dPlmS, zZ, rZS )
         deallocate( chebt_Z, nZmaxS )
         if ( rank == 0 .and. (.not. l_save_out) ) close(n_geos_file)
      end if

      deallocate( cyl_balance )

   end subroutine finalize_geos_mod
 
!----------------------------------------------------------------------------

   subroutine getEgeos(time,nGeosSets,w,dw,ddw,z,dz,Geos,GeosA,GeosZ,GeosM,GeosNA, &
   &                   GeosNAP,dpFlow,dzFlow,volume,Ekin)
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
      real(cp), intent(out) :: Geos    ! ratio of geostrophic to total energy
      real(cp), intent(out) :: GeosA   ! ratio of geostr. axisymm. energy
      real(cp), intent(out) :: GeosZ   ! ratio of geostr. zonal energy
      real(cp), intent(out) :: GeosM   ! ratio of geostr. meridional energy
      real(cp), intent(out) :: GeosNA  ! ration of geostr. non-axisymmetric energy
      real(cp), intent(out) :: GeosNAP ! ration of geostr. non-axisymmetric perp energy
      real(cp), intent(out) :: dpFlow  ! RMS lengths scale
      real(cp), intent(out) :: dzFlow ! RMS lengths scale
      real(cp), intent(out) :: Ekin

      real(cp), intent(out) :: volume
      real(cp) :: volume_s,volumeOTC_s,volumeOTC

      !-- Local variables:
      real(cp) :: Egeos,EkNTC,EkSTC
      real(cp) :: EgeosA,EgeosZ,EgeosM,EgeosNA,EgeosNAP
      real(cp) :: CVzOTC,CVorOTC,CHelOTC
      integer :: nS,nS_ICB,kindCalc
      real(cp) :: zMinN,zMaxN,zMin,zMax ! integration boundarie, help variable
      real(cp) :: sZ(nSmax),dsZ,s1,s2 ! cylindrical radius s and s-step
      integer :: nPhi,nI
      real(cp) :: phiNorm
      logical :: lTC


      !-- Representation in (phi,z):
      real(cp) :: VrS(nrp_geos,nZmaxA),Vr(nZmaxA),VrIntS
      real(cp) :: VtS(nrp_geos,nZmaxA),Vt(nZmaxA),VtIntS
      real(cp) :: VpS(nrp_geos,nZmaxA),Vp(nZmaxA),VpIntS
      real(cp) :: VozS(nrp_geos,nZmaxA),V2(nZmaxA)
      real(cp) :: VrAS(nZmaxA),VrAIntS
      real(cp) :: VtAS(nZmaxA),VtAIntS
      real(cp) :: VpAS(nZmaxA),VpAIntS
      real(cp) :: VrNA(nZmaxA),VrNAIntS
      real(cp) :: VtNA(nZmaxA),VtNAIntS
      real(cp) :: VpNA(nZmaxA),VpNAIntS
      real(cp) :: VsNA(nZmaxA),VsNAIntS
      real(cp) :: VzNA(nZmaxA),VzNAIntS
      real(cp) :: VA2(nZmaxA),VZ2(nZmaxA),VM2(nZmaxA)
      real(cp) :: VNA2(nZmaxA),VNAP2(nZmaxA)
      real(cp) :: rZ(nZmaxA),sinT(nZmaxA),cosT(nZmaxA)
      integer :: nInt,nInts   ! index for NHS and SHS integral
      integer :: nZ,nZmaxH,nZmaxI,nZmaxV,nZS,nZN
      real(cp) :: EA_s,EAIntS,EA,EZ_s,EZIntS,EZ,EM_s,EMIntS,EM
      real(cp) :: ENA_s,ENAIntS,ENA,EkSTC_s,EkNTC_s
      real(cp) :: ENAP_s,ENAPIntS,ENAP
      real(cp) :: Egeos_s,EkOTC_s,EgeosA_s,EgeosZ_s,EgeosM_s
      real(cp) :: EgeosNA_s,EgeosNAP_s,EkOTC,EkIntS
      real(cp) :: dpEkInt(nZmaxA),dpEkIntS,dzEkInt(nZmaxA),dzEkIntS
      real(cp) :: dpEk_s,dzEk_s,thetaZ,wZ,wZP,threehalf
      real(cp) :: dpEk,dzEk

      !-- Correlation (new Oct. 4 2007)
      logical :: lCorrel
      real(cp) :: VzS,VzN,VorS,VorN,surf,delz
      real(cp) :: VzSN,VzSS,VzNN,VorSN,VorSS,VorNN,HelZZ,VZZ,VorZZ
      real(cp) :: CVz_I,CVor_I,CHel_I
      real(cp) :: CVz_s,CVor_s,CHel_s

      !-- Movie output
      integer :: nOutFile,n
      character(len=66) :: movFile
      character(len=64) :: version
      integer :: nFields,nFieldSize
      real(cp) :: dumm(40)
      real(outp) :: CVz(nrp_geos,nSmax)
      real(outp) :: CVor(nrp_geos,nSmax)
      real(outp) :: CHel(nrp_geos,nSmax)
      real(outp) :: CVz_Sloc(nrp_geos,nSstart:nSstop)
      real(outp) :: CVor_Sloc(nrp_geos,nSstart:nSstop)
      real(outp) :: CHel_Sloc(nrp_geos,nSstart:nSstop)

#ifdef WITH_MPI
      integer :: i,sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1)
#endif

      if ( lVerbose ) write(*,*) '! Starting outGeos!'

      call costf_arrays(w,dw,ddw,z,dz)

      lCorrel=.true. ! Calculate Vz and Vorz north/south correlation
      phiNorm=two*pi/n_phi_max
      ! if lDeriv=true chebIntD returns the z-derivative ...
      kindCalc=1

      !---- Get resolution in s, z resolution is defined above
      dsZ  =r_CMB/real(nSmax,cp)  ! Step in s controlled by nSmax

      !-- Construct new cylindrical grid so that first s=dsz/2 and
      !   last s is r_CMB-dsZ/2 
      !   The z grid points are calculated in chebIntInit
      do nS=1,nSmax
         sZ(nS)=(real(nS,cp)-half)*dsZ
      end do

      EkSTC    =0.0_cp
      EkNTC    =0.0_cp
      EkOTC    =0.0_cp
      Egeos    =0.0_cp
      EgeosA   =0.0_cp
      EgeosZ   =0.0_cp
      EgeosM   =0.0_cp
      EgeosNA  =0.0_cp
      EgeosNAP =0.0_cp
      EA       =0.0_cp
      EZ       =0.0_cp
      EM       =0.0_cp
      ENA      =0.0_cp
      ENAP     =0.0_cp
      CVzOTC   =0.0_cp
      CVorOTC  =0.0_cp
      CHelOTC  =0.0_cp
      volume   =0.0_cp
      volumeOTC=0.0_cp
      dpEk     =0.0_cp
      dzEk     =0.0_cp


      !---- Contributions are now in fully spectral space!
      !---- Do the z-integral:
      nI=0

      !--- Start the big loop over s grid points. This loop is 
      !    parallized in MPI.
      do nS=nSstart,nSstop 

         if ( sZ(nS) < r_ICB ) then
            lTC=.true.
         else
            lTC=.false.
         end if
         if ( nS > nSstart ) then
            if ( sZ(nS-1) < r_ICB .and. sZ(nS) >= r_ICB ) nS_ICB=nS
         end if


      ! zero contributions for this s:
         volume_s   =0.0_cp
         volumeOTC_s=0.0_cp
         EkSTC_s    =0.0_cp
         EkNTC_s    =0.0_cp
         EkOTC_s    =0.0_cp
         EA_s       =0.0_cp
         EZ_s       =0.0_cp
         EM_s       =0.0_cp
         ENA_s      =0.0_cp
         ENAP_s     =0.0_cp
         Egeos_s    =0.0_cp
         EgeosA_s   =0.0_cp
         EgeosZ_s   =0.0_cp
         EgeosM_s   =0.0_cp
         EgeosNA_s  =0.0_cp
         EgeosNAP_s =0.0_cp
         dpEk_s     =0.0_cp
         dzEk_s     =0.0_cp
         CVz_s      =0.0_cp
         CVor_s     =0.0_cp
         CHel_s     =0.0_cp


         !--- Get lower integral boundary for this s:
         zMinN=-sqrt(r_CMB*r_CMB-sZ(nS)*sZ(nS))
         if ( lTC ) then
            zMaxN=-sqrt(r_ICB*r_ICB-sZ(nS)*sZ(nS))
         else
            zMaxN=-zMinN
         end if

         !--- Get the z-grid for the integral, only done once
         if ( nGeosSets==1 ) then
            !--- Prepare z-integration
            !---  Inside TC this covers only one hemiphere. 
            !     Output of chebIntInit is the z-grid zZ and the 
            !     number of points nZmaxS. 
            call chebIntInit(zMinN,zMaxN,zNorm,nZnorm,                   &
                 &           nZmaxA,zZ(:,nS),nZmaxS(nS),chebt_Z(nS))
         endif

         ! Number of z points:
         !    nZmaxH = points in one hemisphere (plus equator outside TC)
         !    nZmaxV = all points for one s, used in getDPV ...
         !    nZmaxI = points used for one intergal
         if ( lTC ) then
            nInts=2
            nZmaxH=nZmaxS(nS)
            nZmaxV=2*nZmaxH
            nZmaxI=nZmaxH
         else
            nInts=1
            nZmaxH=(nZmaxS(nS)-1)/2+1 ! Hemisphere plus equator
            nZmaxV=nZmaxS(nS)
            nZmaxI=nZmaxV
         end if

         !--- Calculate the Plm, sin(theta) for the southern hermisphere
         !    Only done once. 
         !    Each processor only calculates the plms for his s!
         if ( nGeosSets==1 ) then
            do nZ=1,nZmaxH ! This is only for one HS
               rZS(nZ,nS)   =sqrt(zZ(nZ,nS)**2+sZ(nS)**2)
               thetaZ       =atan2(sZ(nS),zZ(nZ,nS))
               sinTS(nZ,nS)=sin(thetaZ)
               call plm_theta(thetaZ,l_max,m_max,minc,                &
               &              PlmS(:,nZ,nS),dPlmS(:,nZ,nS),lm_max,2)
            end do
         end if 

         !--- Get sin,cos for the whole sphere (nZmaxV points
         sinT(1:nZmaxH)=sinTS(1:nZmaxH,nS)
         cosT(1:nZmaxH)=sqrt(one-sinT(1:nZmaxH))
         rZ(1:nZmaxH)=rZS(1:nZmaxH,nS)
         if ( lTC ) then 
            ! Northern HS
            sinT(nZmaxH+1:nZmaxV)= sinT(nZmaxH:1:-1)
            cosT(nZmaxH+1:nZmaxV)=-cosT(nZmaxH:1:-1)
            rZ(nZmaxH+1:nZmaxV)  =rZ(nZmaxH:1:-1)
         else
            ! Nothern HS, mind the equator !
            nZ=nZmaxH-1
            sinT(nZmaxH+1:nZmaxV)=sinT(nZ:1:-1)
            cosT(nZmaxH+1:nZmaxV)=-cosT(nZ:1:-1)
            rZ(nZmaxH+1:nZmaxV)  =rZ(nZ:1:-1)
         end if

         ! Weight for individual z-integral contribution,
         ! NOTE: another two pi is multiplied later
         s1=sZ(nS)-half*dsZ
         s2=sZ(nS)+half*dsZ
         threehalf=three*half
         if ( lTC ) then
            wZ=one/three*( (r_ICB**2-s2**2)**threehalf - &
            &              (r_ICB**2-s1**2)**threehalf - &
            &              (r_CMB**2-s2**2)**threehalf + &
            &              (r_CMB**2-s1**2)**threehalf )
         else
            wZ=two/three*( (r_CMB**2-s1**2)**threehalf - &
            &              (r_CMB**2-s2**2)**threehalf ) 
         end if
         wZP=wZ*phiNorm

         !--------- Get the flow components for all northern and
         !          southern thetas and all phis:
         call getDVptr(wS_global,dwS_global,ddwS_global,zS_global,dzS_global, &
              &        r_ICB,r_CMB,rZ,nZmaxV,nZmaxA,PlmS(:,:,nS),       &
              &        dPlmS(:,:,nS),sinT,cosT,kindCalc,VrS,VtS,VpS,VozS,  &
              &        VrAS,VtAS,VpAS,dpEkInt)

         !------- Perform z-integrals for axisymmetric stuff:
         do nInt=1,nInts
            VZ2(:)=0.0_cp
            if ( nInt == 1 ) then
               zMax=zMaxN
               zMin=zMinN
               Vr(1:nZmaxI)=VrAS(1:nZmaxI)
               Vt(1:nZmaxI)=VtAS(1:nZmaxI)
               Vp(1:nZmaxI)=VpAS(1:nZmaxI)
            else if ( nInt == 2 ) then ! second integral inside TC
               zMax=-zMinN
               zMin=-zMaxN
               Vr(1:nZmaxI)=VrAS(nZmaxI+1:nZmaxV)
               Vt(1:nZmaxI)=VtAS(nZmaxI+1:nZmaxV)
               Vp(1:nZmaxI)=VpAS(nZmaxI+1:nZmaxV)
            end if
            VA2(1:nZmaxI)=Vr(1:nZmaxI)**2+Vt(1:nZmaxI)**2+Vp(1:nZmaxI)**2
            VZ2(1:nZmaxI)=Vp(1:nZmaxI)**2
            VM2(1:nZmaxI)=Vr(1:nZmaxI)**2+Vt(1:nZmaxI)**2
            !--- Perform z-integral:
            VrAIntS=chebIntD(Vr,.false.,zMin,zMax,nZmaxI,nZmaxA,chebt_Z(nS))
            VtAIntS=chebIntD(Vt,.false.,zMin,zMax,nZmaxI,nZmaxA,chebt_Z(nS))
            VpAIntS=chebIntD(Vp,.false.,zMin,zMax,nZmaxI,nZmaxA,chebt_Z(nS))
            EAIntS=chebIntD(VA2,.false.,zMin,zMax,nZmaxI,nZmaxA,chebt_Z(nS))
            EZIntS=chebIntD(VZ2,.false.,zMin,zMax,nZmaxI,nZmaxA,chebt_Z(nS))
            EMIntS=chebIntD(VM2,.false.,zMin,zMax,nZmaxI,nZmaxA,chebt_Z(nS))
            ! Note: the factor is 2*pi*wZ(V^2/2) = pi*wZ*v^2
            EgeosA_s=EgeosA_s+pi*wZ*(VrAIntS**2+VtAIntS**2+VpAIntS**2)
            EgeosZ_s=EgeosZ_s+pi*wZ*VpAIntS**2
            EgeosM_s=EgeosM_s+pi*wZ*(VrAIntS**2+VtAIntS**2)
            EA_s=EA_s+pi*wZ*EAIntS
            EZ_s=EZ_s+pi*wZ*EZIntS
            EM_s=EM_s+pi*wZ*EMIntS
            volume_s=volume_s+two*pi*wZ
         end do

         !------- Perform z-integral(s) for all phis:

         do nPhi=1,n_phi_max

            do nInt=1,nInts
               if ( nInt == 1 ) then
                  zMax=zMaxN
                  zMin=zMinN
                  Vr(1:nZmaxI)=VrS(nPhi,1:nZmaxI)
                  Vt(1:nZmaxI)=VtS(nPhi,1:nZmaxI)
                  Vp(1:nZmaxI)=VpS(nPhi,1:nZmaxI)
                  VrNA(1:nZmaxI)=Vr(1:nZmaxI)-VrAS(1:nZmaxI)
                  VtNA(1:nZmaxI)=Vt(1:nZmaxI)-VtAS(1:nZmaxI)
                  VpNA(1:nZmaxI)=Vp(1:nZmaxI)-VpAS(1:nZmaxI)
               else if ( nInt == 2 ) then
                  zMax=-zMinN
                  zMin=-zMaxN
                  Vr(1:nZmaxI)  =VrS(nPhi,nZmaxI+1:2*nZmaxI)
                  Vt(1:nZmaxI)  =VtS(nPhi,nZmaxI+1:2*nZmaxI)
                  Vp(1:nZmaxI)  =VpS(nPhi,nZmaxI+1:2*nZmaxI)
                  VrNA(1:nZmaxI)=Vr(1:nZmaxI)-VrAS(nZmaxI+1:2*nZmaxI)
                  VtNA(1:nZmaxI)=Vt(1:nZmaxI)-VtAS(nZmaxI+1:2*nZmaxI)
                  VpNA(1:nZmaxI)=Vp(1:nZmaxI)-VpAS(nZmaxI+1:2*nZmaxI)
               end if
               VsNA(1:nZmaxI)=sinT(1:nZmaxI)*VrNA(1:nZmaxI) + &
               &              cosT(1:nZmaxI)*VtNA(1:nZmaxI)
               VzNA(1:nZmaxI)=cosT(1:nZmaxI)*VrNA(1:nZmaxI) - & 
               &              sinT(1:nZmaxI)*VtNA(1:nZmaxI)
               V2(:)=Vr(:)**2+Vt(:)**2+Vp(:)**2
               VNA2(:)=VrNA(:)**2+VtNA(:)**2+VpNA(:)**2 
               VNA2(:)=VrNA(:)**2+VtNA(:)**2+VpNA(:)**2 
               VNAP2(:)=VsNA(:)**2+VpNA(:)**2 

               !--- Perform z-integral:
               !-------- NOTE: chebIntD replaces VrInt with z-derivative
               !               for lDeriv=.true. (second parameter)
               VrIntS=chebIntD(Vr,.true.,zMin,zMax,nZmaxI,nZmaxA,chebt_Z(nS))
               VtIntS=chebIntD(Vt,.true.,zMin,zMax,nZmaxI,nZmaxA,chebt_Z(nS))
               VpIntS=chebIntD(Vp,.true.,zMin,zMax,nZmaxI,nZmaxA,chebt_Z(nS))
               EkIntS=chebIntD(V2,.false.,zMin,zMax,nZmaxI,nZmaxA,chebt_Z(nS))
               VrNAIntS=chebIntD(VrNA,.false.,zMin,zMax,nZmaxI,nZmaxA,chebt_Z(nS))
               VtNAIntS=chebIntD(VtNA,.false.,zMin,zMax,nZmaxI,nZmaxA,chebt_Z(nS))
               VpNAIntS=chebIntD(VpNA,.false.,zMin,zMax,nZmaxI,nZmaxA,chebt_Z(nS))
               VsNAIntS=chebIntD(VsNA,.false.,zMin,zMax,nZmaxI,nZmaxA,chebt_Z(nS))
               ENAIntS=chebIntD(VNA2,.false.,zMin,zMax,nZmaxI,nZmaxA,chebt_Z(nS))
               ENAPIntS=chebIntD(VNAP2,.false.,zMin,zMax,nZmaxI,nZmaxA,chebt_Z(nS))

               !-------- Get volume integral of geostrophic energies:
               Egeos_s=Egeos_s+(VrIntS**2+VtIntS**2+VpIntS**2)
               EgeosNA_s =EgeosNA_s +(VrNAIntS**2+VtNAIntS**2+VpNAIntS**2) 
               EgeosNAP_s=EgeosNAP_s+(VsNAIntS**2+VpNAIntS**2) 
               ENA_s=ENA_s+ENAIntS
               ENAP_s=ENAP_s+ENAPIntS
               if ( lTC ) then
                  if ( nInt == 1 ) then
                     EkSTC_s=EkSTC_s+EkIntS
                  else if ( nInt == 2 ) then
                     EkNTC_s=EkNTC_s+EkIntS
                  end if
               else
                  EkOTC_s=EkOTC_s+EkIntS
               end if

               !-------- Note: chebIntD above returns the z derivative 
               dzEkInt(:)=Vr(:)**2+Vt(:)**2+Vp(:)**2
               dzEkIntS=chebInt(dzEkInt,zMin,zMax,nZmaxI,nZmaxA,chebt_Z(nS))
               dzEk_s=dzEk_s+dzEkIntS

            end do ! Loop over north/south integral


            ! --- Here I calculate the Peason correlation coeff between
            !     z-integrated stuff in northern and southern HS.
            !     Considered are Vz,z-vorticity Vor, and axial helicity Vz*Vor.
            !     The correlation is averaged over the shell.
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
               do nZN=1,nZmaxH 
                  nZS  =nZmaxI-nZN+1
                  delz =zZ(nZN,nS)-zZ(nZN+1,nS)
                  VzN  = cosT(nZN)*VrS(nPhi,nZN)-sinT(nZN)*VtS(nPhi,nZN)
                  VzS  =-cosT(nZN)*VrS(nPhi,nZS)-sinT(nZN)*VtS(nPhi,nZS)
                  VorN =VozS(nPhi,nZN)
                  VorS =VozS(nPhi,nZS)
                  VzSN =VzSN+delz*VzS*VzN
                  VzNN =VzNN+delz*VzN*VzN
                  VzSS =VzSS+delz*VzS*VzS
                  VorSN=VorSN+delz*VorS*VorN
                  VorNN=VorNN+delz*VorN*VorN
                  VorSS=VorSS+delz*VorS*VorS
                  HelZZ=HelZZ+delz*(VorN*VzN-VorS*VzS)
                  VZZ  =VZZ  +delz*(VzN*VzN+VzS*VzS)
                  VorZZ=VorZZ+delz*(VorN*VorN+VorS*VorS)
               end do

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
               if ( ( VZZ /= 0.0_cp .and. VorZZ /= 0.0_cp ) .and. HelZZ > eps ) then
                  CHel_I=HelZZ/sqrt(VZZ*VorZZ)
               else
                  CHel_I=0.0_cp
               end if
               CVz_Sloc(nPhi,nS) =real(CVz_I,kind=outp)
               CVor_Sloc(nPhi,nS)=real(CVor_I,kind=outp)
               CHel_Sloc(nPhi,nS)=real(CHel_I,kind=outp)

               CVz_s =CVz_s +CVz_I
               CVor_s=CVor_s+CVor_I
               CHel_s=CHel_s+CHel_I

            end if ! lCorrel ?

         end do  ! Loop over phi

         ! Note: the extral factor half is needed to convert v^2 to energy
         ! Factor wZP is the weight for the specific s.
         Egeos_s   =wZP*half*Egeos_s
         EgeosNA_s =wZP*half*EgeosNA_s
         EgeosNAP_s=wZP*half*EgeosNAP_s
         ENA_s     =wZP*half*ENA_s
         ENAP_s    =wZP*half*ENAP_s
         if ( lTC ) then
            EkSTC_s=wZP*half*EkSTC_s
            EkNTC_s=wZP*half*EkNTC_s
         else
            EkOTC_s=wZP*half*EkOTC_s
         end if
         dzEk_s=wZP*half*dzEk_s

         if ( lCorrel .and. .not.lTC ) then
            CVz_s =wZP*CVz_s
            CVor_s=wZP*CVor_s
            CHel_s=wZP*CHel_s
            volumeOTC_s=volumeOTC_s+two*pi*wZ
         end if 

         !--------- dpEkInt treated differently cause phi integral has
         !          already been performed by getDVptr
         do nInt=1,nInts
            if ( nInt == 1 ) then
               zMax=-zMaxN
               zMin=-zMinN
            elseif ( nInt == 2 ) then
               zMax=-zMinN
               zMin=-zMaxN
               do nZ=1,nZmaxI
                  dpEkInt(nZ)=dpEkInt(nZ+nZmaxI)
               end do
            end if
            ! Perform the z-integral with cheb integration:
            dpEkIntS=chebInt(dpEkInt,zMin,zMax,nZmaxI,nZmaxA,chebt_Z(nS))
            ! Multiply with weight wZ for s and also apply 1/2*1/(2*pi*s)**2
            ! to convert into energy and azimuthal instead of phi derivative
            ! NOTE: the factor 1/2 is contained in wZ
            dpEk_s=dpEk_s+wZ*dpEkIntS/((two*pi*sZ(nS))**2) 
         end do


   ! Collect, integrate in s:
         Egeos   =Egeos   +Egeos_s
         EgeosA  =EgeosA  +EgeosA_s
         EgeosZ  =EgeosZ  +EgeosZ_s
         EgeosM  =EgeosM  +EgeosM_s
         EgeosNA =EgeosNA +EgeosNA_s
         EgeosNAP=EgeosNAP+EgeosNAP_s
         EA      =EA      +EA_s
         EZ      =EZ      +EZ_s
         EM      =EM      +EM_s
         ENA     =ENA     +ENA_s
         ENAP    =ENAP    +ENA_s
         dpEk    =dpEk    +dpEk_s
         dzEk    =dzEk    +dzEk_s
         volume  =volume  +volume_s
         if ( lTC ) then
            EkSTC=EkSTC+EkSTC_s
            EkNTC=EkNTC+EkNTC_s
         else
            EkOTC=EkOTC+EkOTC_s
            if ( lCorrel ) then
               CVzOTC =CVzOTC +CVz_s
               CVorOTC=CVorOTC+CVor_s
               CHelOTC=CHelOTC+CHel_s
               volumeOTC=volumeOTC+volumeOTC_s
            end if
         end if

      end do ! big loop over s

#ifdef WITH_MPI
      call MPI_Allreduce(MPI_IN_PLACE, volume, 1, MPI_DEF_REAL, MPI_SUM,  &
           &             MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, EkSTC, 1, MPI_DEF_REAL, MPI_SUM,  &
           &             MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, EkNTC, 1, MPI_DEF_REAL, MPI_SUM,  &
           &             MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, EkOTC, 1, MPI_DEF_REAL, MPI_SUM,  &
           &             MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, Egeos, 1, MPI_DEF_REAL, MPI_SUM,  &
           &             MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, EgeosA, 1, MPI_DEF_REAL, MPI_SUM,  &
           &             MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, EgeosZ, 1, MPI_DEF_REAL, MPI_SUM,  &
           &             MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, EgeosM, 1, MPI_DEF_REAL, MPI_SUM,  &
           &             MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, EgeosNA,1, MPI_DEF_REAL, MPI_SUM,  &
           &             MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, EgeosNAP,1, MPI_DEF_REAL, MPI_SUM,  &
           &             MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, EA, 1, MPI_DEF_REAL, MPI_SUM,  &
           &             MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, EZ, 1, MPI_DEF_REAL, MPI_SUM,  &
           &             MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, EM, 1, MPI_DEF_REAL, MPI_SUM,  &
           &             MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, ENA, 1, MPI_DEF_REAL, MPI_SUM,  &
           &             MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, ENAP, 1, MPI_DEF_REAL, MPI_SUM,  &
           &             MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, dpEk, 1, MPI_DEF_REAL, MPI_SUM,  &
           &             MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, dzEk, 1, MPI_DEF_REAL, MPI_SUM,  &
           &             MPI_COMM_WORLD, ierr)
      if ( lCorrel ) then
         call MPI_Allreduce(MPI_IN_PLACE, volumeOTC, 1, MPI_DEF_REAL, MPI_SUM,  &
              &             MPI_COMM_WORLD, ierr)
         call MPI_Allreduce(MPI_IN_PLACE, CVzOTC, 1, MPI_DEF_REAL, MPI_SUM, &
              &             MPI_COMM_WORLD, ierr)
         call MPI_Allreduce(MPI_IN_PLACE, CVorOTC, 1, MPI_DEF_REAL, MPI_SUM,&
              &             MPI_COMM_WORLD, ierr)
         call MPI_Allreduce(MPI_IN_PLACE, CHelOTC, 1, MPI_DEF_REAL, MPI_SUM,&
              &             MPI_COMM_WORLD, ierr)
      end if
#endif

      if ( lCorrel ) then
         CVzOTC =CVzOTC/volumeOTC
         CVorOTC=CVorOTC/volumeOTC
         CHelOTC=CHelOTC/volumeOTC
      end if
 
      ! Calculate total energy, can be used for testing:
      Ekin=EkSTC+EkNTC+EkOTC 

      ! Calculate flow length scales based on kinetic energy:
      if ( dpEk/= 0.0_cp ) then
         dpFlow = sqrt(Ekin/dpEk)
      else
         dpFlow = 0.0_cp
      end if
      if ( dzEk /= 0.0_cp ) then
         dzFlow = sqrt(Ekin/dzEk)
      else
         dzFlow = 0.0_cp
      end if

      Geos   = 0.0_cp
      GeosA  = 0.0_cp
      GeosZ  = 0.0_cp
      GeosM  = 0.0_cp
      GeosNA = 0.0_cp
      GeosNAP= 0.0_cp
      if ( Ekin>0.0_cp ) Geos   =Egeos/Ekin     ! rel. geostrophic kinetic energy
      if ( EA>0.0_cp )   GeosA  =EgeosA/EA      ! rel. geostr. axisymmtric kinetic energy
      if ( EZ>0.0_cp )   GeosZ  = EgeosZ/EZ     ! rel. geostr. zonal kinetic energy
      if ( EM>0.0_cp )   GeosM  = EgeosM/EM     ! rel. geostr. meridional kinetic energy
      if ( ENA>0.0_cp )  GeosNA = EgeosNA/ENA   ! rel. geostr. non-axisymmetric kinetic energy
      if ( ENAP>0.0_cp ) GeosNAP= EgeosNAP/ENAP ! rel. geostr. non-axisymmetric kinetic energy
                                                ! in flow perpendicular to rotation axis

      if ( l_corrMov ) then

         !--- Determine s used for correl
         n=0
         do nS=nSstart,nSstop
            if ( sZ(nS) >= r_ICB+0.1_cp .and. sZ(nS) <= r_CMB-0.1_cp ) then
               n=n+1
            else
               do nPhi=1,n_phi_max
                  CVz_Sloc(nPhi,nS) =0.0_cp
                  CVor_Sloc(nPhi,nS)=0.0_cp
                  CHel_Sloc(nPhi,nS)=0.0_cp
               end do
            end if
         end do

#ifdef WITH_MPI
         sendcount  = (nSstop-nSstart+1)*nrp_geos
         do i=0,n_procs-1
            recvcounts(i)=cyl_balance(i)%n_per_rank
         end do
         displs(0)=0
         do i=1,n_procs-1
            displs(i) = displs(i-1)+recvcounts(i-1)
         end do

         call MPI_GatherV(CVZ_Sloc,sendcount,MPI_OUT_REAL,     &
              &           CVZ,recvcounts,displs,MPI_OUT_REAL,  &
              &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(CVor_Sloc,sendcount,MPI_OUT_REAL,    &
              &           CVor,recvcounts,displs,MPI_OUT_REAL, &
              &           0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(CHel_Sloc,sendcount,MPI_OUT_REAL,    &
              &           CHel,recvcounts,displs,MPI_OUT_REAL, &
              &           0,MPI_COMM_WORLD,ierr)
#else
         CVz(:,:) =CVz_Sloc(:,:)
         CVor(:,:)=CVor_Sloc(:,:)
         CHel(:,:)=CHel_Sloc(:,:)
#endif

      end if

      if ( rank == 0 ) then

         if ( l_save_out ) then
            open(newunit=n_geos_file, file=geos_file, status='unknown', &
            &    position='append')
         end if

         ! JW28April2020: GeosNA replaced with GeosNAP (non-axi flow perp. to
         ! rotation axis (Vs,Vp)., GeosNA could be added to output...
         write(n_geos_file,'(1P,ES20.12,11ES16.8)') time, Geos, EkNTC/Ekin, &
         &     EkSTC/Ekin, Ekin, CVzOTC, CVorOTC, CHelOTC, GeosA, GeosZ,    &
         &     GeosM, GeosNAP

         if ( l_save_out ) close(n_geos_file)
         !--- NOTE: Ekin can be compared with energy in e_kin.TAG to
         !    get an idea of the precision of cylindrical integration in getEgeos.

         !        write(99,*) 'E:',EkSTC,EkNTC,EkOTC
         !        write(99,*) 'Ekin:',Ekin
         !        write(99,*) 'Egeos:',Egeos

         !--- Write correlation movie:
         if ( l_corrMov ) then

            movFile ='CVorz_mov.'//tag
            open(newunit=nOutFile, file=movFile, status='unknown',   &
            &    form='unformatted', position='append')

            !--- Write header into output file:
            if ( nGeosSets == 1 ) then

               nFields=3
               nFieldSize=(nSmax-nS_ICB+1)*n_phi_max
               version='JW_Movie_Version_2'
               write(nOutFile) version
               dumm(1)=111           ! type of input
               dumm(2)=2             ! marker for constant theta plane
               dumm(3)=90.0_cp       ! surface constant
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

            write(nOutFile) ((CVz(nPhi,nS) ,nPhi=1,n_phi_max), nS=nSmax,1,-1)
            write(nOutFile) ((CVor(nPhi,nS),nPhi=1,n_phi_max), nS=nSmax,1,-1)
            write(nOutFile) ((CHel(nPhi,nS),nPhi=1,n_phi_max), nS=nSmax,1,-1)

            close(nOutFile)

         end if ! l_corrMov

      end if ! rank == 0

      if ( lVerbose ) write(*,*) '! End of getGeos!'

   end subroutine getEgeos
!---------------------------------------------------------------------------------
   subroutine getDVptr(w,dw,ddw,z,dz,rMin,rMax,rS,nZmax,nZmaxA,PlmS,dPlmS, &
              &        sinT,cosT,kindCalc,VrS,VtS,VpS,VorS,VrAS,VtAS,VpAS,dpEk)
      !
      !  This subroutine calculates the three flow components VrS,VtS,VpS at
      !  r=rS, all phis, and a list of nZmax theta values defined by
      !  PlmS=Plm(theta), dPlmS=sin(theta)*dTheta Plm(theta), and sinT=sin(theta).
      !  Also calculated are the axisymmetric flow contributions VpAS,VtAS,VpAS
      !  and the z-vorticity VorS (z-component of curl V).
      !  NOTE: For kindCalc=2 the cylindrical flow contributions are returned:
      !  VrS=VsS,VtS=VzS,VrAS=VsAS,VtAS=VzAS.
      !  Also calculated is dpEk, the  phi average of 
      !    (d Vr/d phi)**2 + (d Vtheta/ d phi)**2 + (d Vphi/ d phi)**2
      !  used to calculate the phi length scale.
      !  .. note:: on input wS=w/r^2, dwS=dw/r, ddwS=ddw/r, zS=z/r
      !

      !--- Input variables:
      complex(cp), intent(in) :: w(lm_max,n_r_max)
      complex(cp), intent(in) :: dw(lm_max,n_r_max)
      complex(cp), intent(in) :: ddw(lm_max,n_r_max)
      complex(cp), intent(in) :: z(lm_max,n_r_max)
      complex(cp), intent(in) :: dz(lm_max,n_r_max)
      real(cp),    intent(in) :: rMin,rMax  ! radial bounds
      integer,     intent(in) :: nZmax,nZmaxA ! number of (r,theta) points
      real(cp),    intent(in) :: rS(nZmaxA)
      real(cp),    intent(in) :: PlmS(lm_max,nZmaxAH)
      real(cp),    intent(in) :: dPlmS(lm_max,nZmaxAH)
      real(cp),    intent(in) :: sinT(nZmaxA)
      real(cp),    intent(in) :: cosT(nZmaxA)
      integer,     intent(in) :: kindCalc

      !--- Output: function on azimuthal grid points defined by FT!
      real(cp), intent(out) :: VrS(nrp_geos,nZmaxA)
      real(cp), intent(out) :: VtS(nrp_geos,nZmaxA)
      real(cp), intent(out) :: VpS(nrp_geos,nZmaxA)
      real(cp), intent(out) :: VrAS(nZmaxA)
      real(cp), intent(out) :: VtAS(nZmaxA)
      real(cp), intent(out) :: VpAS(nZmaxA)
      real(cp), intent(out) :: VorS(nrp_geos,nZmaxA)
      real(cp), optional, intent(out) :: dpEk(nZmaxA)

      !--- Local variables:
      real(cp) :: chebS(n_r_max)
      integer :: nSmaxH,nS,nN,mc,lm,l,m,nCheb,nPhi,n
      integer :: nZmaxH,nEquator
      real(cp) :: x,phiNorm,mapFac,OS,Or_e1,Or_e2
      complex(cp) :: Vr,Vt1,Vt2,Vp1,Vp2,Vor,Vot1,Vot2
      real(cp) :: VotS(nrp_geos,nZmaxA)
      complex(cp) :: wSr,dwSr,ddwSr,zSr,dzSr
      real(cp) :: phi_norm

      real(cp) :: dV(nrp_geos,nZmaxA)
      complex(cp) :: dp


      mapFac=two/(rMax-rMin)
      phiNorm=two*pi/n_phi_max

      VrS(:,:) =zero   
      VtS(:,:) =zero   
      VpS(:,:) =zero   
      VorS(:,:)=zero   
      VotS(:,:)=zero   
      if ( mod(nZmax,2)==0 ) then
         nZmaxH=nZmax/2
      else
         nZmaxH=(nZmax-1)/2+1 ! Number of points including equator
         nEquator=nZmaxH 
      end if

      do nN=1,nZmaxH    ! Loop over all (r,theta) points in NHS
         nS=nZmax-nN+1   ! Southern counterpart !

         !------ Calculate Chebs:
         !------ Map r to cheb interval [-1,1]:
         !       and calculate the cheb polynomial:
         !       Note: the factor cheb_norm is needed
         !       for renormalisation. Its not needed if one used
         !       costf1 for the back transform.
         x=two*(rS(nN)-half*(rMin+rMax))/(rMax-rMin)
         chebS(1) =one*rscheme_oc%rnorm ! Extra cheb_norm cheap here
         chebS(2) =x*rscheme_oc%rnorm
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
               wSr  =  wSr+  w(lm,nCheb)*chebS(nCheb)
               dwSr = dwSr+ dw(lm,nCheb)*chebS(nCheb)
               ddwSr=ddwSr+ddw(lm,nCheb)*chebS(nCheb)
               zSr  =  zSr+  z(lm,nCheb)*chebS(nCheb)
               dzSr = dzSr+ dz(lm,nCheb)*chebS(nCheb)
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

            if ( nN /= nEquator ) then  ! The other hemisphere 
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
            end if

         end do

      end do

      !--- Extra factors, contructing z-vorticity:
      do nS=1,nZmax 
         OS   =one/sinT(nS)
         Or_e1=one/rS(nS)
         Or_e2=Or_e1*Or_e1
         do mc=1,n_m_max
            if ( kindCalc == 1 ) then
               VrS(2*mc-1,nS) =Or_e2*VrS(2*mc-1,nS)
               VrS(2*mc  ,nS) =Or_e2*VrS(2*mc  ,nS)
               VtS(2*mc-1,nS) =Or_e1*OS*VtS(2*mc-1,nS)
               VtS(2*mc  ,nS) =Or_e1*OS*VtS(2*mc  ,nS)
            else if ( kindCalc == 2 ) then
               !-- This is now Vs
               VrS(2*mc-1,nS)=sinT(nS)*Or_e2*VrS(2*mc-1,nS)+cosT(nS)*Or_e1*OS*VtS(2*mc-1,nS)
               VrS(2*mc  ,nS)=sinT(nS)*Or_e2*VrS(2*mc  ,nS)+cosT(nS)*Or_e1*OS*VtS(2*mc  ,nS)
               !-- This is now Vz
               VtS(2*mc-1,nS)=cosT(nS)*Or_e2*VrS(2*mc-1,nS)-sinT(nS)*Or_e1*OS*VtS(2*mc-1,nS)
               VtS(2*mc  ,nS)=cosT(nS)*Or_e2*VrS(2*mc  ,nS)-sinT(nS)*Or_e1*OS*VtS(2*mc  ,nS)
            end if
            VpS(2*mc-1,nS) =Or_e1*OS*VpS(2*mc-1,nS)
            VpS(2*mc  ,nS) =Or_e1*OS*VpS(2*mc  ,nS)
            VorS(2*mc-1,nS)=cosT(nS)*Or_e2*VorS(2*mc-1,nS) - Or_e1*VotS(2*mc-1,nS)
            VorS(2*mc  ,nS)=cosT(nS)*Or_e2*VorS(2*mc  ,nS) - Or_e1*VotS(2*mc  ,nS)
         end do
      end do

      if ( present(dpEk) ) then

         !--- Calculate phi derivative in lm-space:
         dpEk(:)=0.0_cp
         do n=1,3
            do nS=1,nZmax
               if ( n == 1 ) then
                  do mc=1,n_m_max
                     dp=ci*real((mc-1)*minc,cp) ! - i m
                     dV(2*mc-1,nS)= real(dp)*VrS(2*mc-1,nS)-aimag(dp)*VrS(2*mc,nS)
                     dV(2*mc  ,nS)=aimag(dp)*VrS(2*mc-1,nS)+ real(dp)*VrS(2*mc,nS)
                  end do
               else if ( n == 2 ) then
                  do mc=1,n_m_max
                     dp=ci*real((mc-1)*minc,cp) ! - i m
                     dV(2*mc-1,nS)= real(dp)*VtS(2*mc-1,nS)-aimag(dp)*VtS(2*mc,nS)
                     dV(2*mc  ,nS)=aimag(dp)*VtS(2*mc-1,nS)+ real(dp)*VtS(2*mc,nS)
                  end do
               else if ( n == 3 ) then
                  do mc=1,n_m_max
                     dp=ci*real((mc-1)*minc,cp) ! - i m
                     dV(2*mc-1,nS)= real(dp)*VpS(2*mc-1,nS)-aimag(dp)*VpS(2*mc,nS)
                     dV(2*mc  ,nS)=aimag(dp)*VpS(2*mc-1,nS)+ real(dp)*VpS(2*mc,nS)
                  end do
               end if
               do mc=2*n_m_max+1,nrp_geos
                  dV(mc,nS)=0.0_cp
               end do
            end do

            !--- Transform m 2 phi for phi-derivative
            if ( .not. l_axi ) call fft_to_real(dV,nrp_geos,nZmax)

            !--- Phi average
            do nS=1,nZmax
               do nPhi=1,n_phi_max
                  if ( mod(nPhi,2) == 1 ) then
                     mc=(nPhi+1)/2
                     dpEk(nS)=dpEk(nS)+dV(2*mc-1,nS)**2 ! Real part
                  else
                     mc=nPhi/2
                     dpEk(nS)=dpEk(nS)+dV(2*mc,nS)**2 ! Imaginary part
                  end if
               end do
            end do

         end do ! Loop over components
         dpEk(:)=phiNorm*dpEk(:) ! Now this is the phi integrated (d V/d phi)**2 
                                 ! for one s and all z.

      end if

      !----- Transform m 2 phi for flow field:
      if ( .not. l_axi ) then
         call fft_to_real(VrS,nrp_geos,nZmax)
         call fft_to_real(VtS,nrp_geos,nZmax)
         call fft_to_real(VpS,nrp_geos,nZmax)
         call fft_to_real(VorS,nrp_geos,nZmax)
         VrAS(:)=zero  
         VtAS(:)=zero  
         VpAS(:)=zero  
         do nS=1,nZmax
            do nPhi=1,n_phi_max
                VrAS(nS)=VrAS(nS)+VrS(nPhi,nS) 
                VtAS(nS)=VtAS(nS)+VtS(nPhi,nS) 
                VpAS(nS)=VpAS(nS)+VpS(nPhi,nS) 
            end do
         end do
         phi_norm=one/n_phi_max
         VrAS(:)=phi_norm*VrAS(:)
         VtAS(:)=phi_norm*VtAS(:)
         VpAS(:)=phi_norm*VpAS(:)
      end if

   end subroutine getDVptr

!----------------------------------------------------------------------------

   subroutine costf_arrays(w,dw,ddw,z,dz)
      !
      ! This subroutine performs local copy of LM-distributed arrays, do
      ! a cosine transform and then allgather them
      !

      !-- Input variables
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dw(llm:ulm,n_r_max)
      complex(cp), intent(in) :: ddw(llm:ulm,n_r_max)
      complex(cp), intent(in) :: z(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dz(llm:ulm,n_r_max)

      !-- Local variables
      complex(cp) :: wS(llm:ulm,n_r_max)
      complex(cp) :: dwS(llm:ulm,n_r_max)
      complex(cp) :: ddwS(llm:ulm,n_r_max)
      complex(cp) :: zS(llm:ulm,n_r_max)
      complex(cp) :: dzS(llm:ulm,n_r_max)

      integer :: nR, lm, l, m

      do nR=1,n_r_max
         do lm=llm,ulm
            l = lo_map%lm2l(lm)
            m = lo_map%lm2m(lm)
            wS(lm,nR)  =orho1(nR)*w(lm,nR)*dLh(st_map%lm2(l,m))
            dwS(lm,nR) =orho1(nR)*dw(lm,nR)
            ddwS(lm,nR)=orho1(nR)*ddw(lm,nR)
            zS(lm,nR)  =orho1(nR)*z(lm,nR)
            dzS(lm,nR) =orho1(nR)*dz(lm,nR)
         end do
      end do

      call rscheme_oc%costf1(wS,ulm-llm+1,1,ulm-llm+1)
      call rscheme_oc%costf1(dwS,ulm-llm+1,1,ulm-llm+1)
      call rscheme_oc%costf1(ddwS,ulm-llm+1,1,ulm-llm+1)
      call rscheme_oc%costf1(zS,ulm-llm+1,1,ulm-llm+1)
      call rscheme_oc%costf1(dzS,ulm-llm+1,1,ulm-llm+1)

      !-- Unfortunately they need to be broadcasted...
      call gather_all_from_lo_to_rank0(gt_OC,wS,wS_global)
      call gather_all_from_lo_to_rank0(gt_OC,dwS,dwS_global)
      call gather_all_from_lo_to_rank0(gt_OC,ddwS,ddwS_global)
      call gather_all_from_lo_to_rank0(gt_OC,zS,zS_global)
      call gather_all_from_lo_to_rank0(gt_OC,dzS,dzS_global)
#ifdef WITH_MPI
      call MPI_Bcast(wS_global,n_r_max*lm_max, MPI_DEF_COMPLEX, 0,  &
           &         MPI_COMM_WORLD, ierr)
      call MPI_Bcast(dwS_global,n_r_max*lm_max, MPI_DEF_COMPLEX, 0, &
           &         MPI_COMM_WORLD, ierr)
      call MPI_Bcast(ddwS_global,n_r_max*lm_max, MPI_DEF_COMPLEX, 0,&
           &         MPI_COMM_WORLD, ierr)
      call MPI_Bcast(zS_global,n_r_max*lm_max, MPI_DEF_COMPLEX, 0,  &
           &         MPI_COMM_WORLD, ierr)
      call MPI_Bcast(dzS_global,n_r_max*lm_max, MPI_DEF_COMPLEX, 0, &
           &         MPI_COMM_WORLD, ierr)
#endif
   end subroutine costf_arrays
!------------------------------------------------------------------------------
end module geos_mod
