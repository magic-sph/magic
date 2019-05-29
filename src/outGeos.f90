module geos_mod
   !
   ! This module is used to calculate the outputs when either l_par
   ! or l_PV are set to .true.
   !

   use precision_mod
   use parallel_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, lm_max, n_m_max, n_phi_max, nrp,     &
       &                 minc, l_max, m_max, l_axi
   use radial_functions, only: r_ICB, r_CMB, rscheme_oc, orho1
   use physical_parameters, only: ra, ek, pr, prmag, radratio
   use num_param, only: tScale
   use blocking, only: lm2l, lm2m, lm2mc, lo_map, st_map
   use horizontal_data, only: dLh, phi, dPhi
   use logic, only: lVerbose, l_corrMov, l_anel, l_save_out, l_SRIC
   use output_data, only: sDens, zDens, tag, runid
   use constants, only: pi, zero, ci, one, two, half
   use LMLoop_data, only: llm,ulm
   use communications, only: gather_all_from_lo_to_rank0,gt_OC
   use plms_theta, only: plm_theta
   use fft, only: fft_to_real
   use TO_helpers, only: getPAStr
   use cosine_transform_odd, only: costf_odd_t
   use chebInt_mod

   implicit none

   private

   real(cp) :: timeOld
   real(cp), allocatable :: PlmS(:,:,:), PlmZ(:,:,:), PlmSPV(:,:,:)
   real(cp), allocatable :: dPlmS(:,:,:), dPlmZ(:,:,:), dPlmSPV(:,:,:)
   real(cp), allocatable :: OsinTS(:,:), OsinTSPV(:,:)

   complex(cp), allocatable :: wS_global(:,:), dwS_global(:,:), ddwS_global(:,:)
   complex(cp), allocatable :: zS_global(:,:), dzS_global(:,:)

   type(costf_odd_t), allocatable :: chebt_Z(:)
   integer, allocatable :: nZmaxS(:), nZC(:), nZC_Sloc(:), nZ2(:,:)
   real(cp), allocatable :: zZ(:,:), rZ(:,:), rZPV(:,:)
   real(cp), allocatable :: VorOld(:,:,:)
   real(cp), parameter :: eps = 10.0_cp*epsilon(one)

   integer :: n_geos_file, nrp_geos
   integer :: nSmax, nZmaxA, nSstart, nSstop
   type(load), allocatable :: cyl_balance(:)
   character(len=72) :: geos_file

   public :: initialize_geos_mod, getEgeos, finalize_geos_mod, outPV

contains

   subroutine initialize_geos_mod(l_geos, l_PV)
      !
      ! Memory allocation in Egeos
      !
      logical, intent(in) :: l_geos
      logical, intent(in) :: l_PV

#ifdef WITH_SHTNS
      nrp_geos=nrp+2 ! One has to include the 2 extra points again
#else
      nrp_geos=nrp
#endif

      nSmax=n_r_max+int(r_ICB*real(n_r_max,cp))
      nSmax=int(sDens*nSmax)

      nZmaxA=2*nSmax

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

      if ( l_PV ) then
         allocate( OsinTSPV(nZmaxA/2+1,nSstart:nSstop) )
         allocate( PlmSPV(lm_max,nZmaxA/2+1,nSstart:nSstop) )
         allocate( dPlmSPV(lm_max,nZmaxA/2+1,nSstart:nSstop) )
         allocate( rZPV(nZmaxA,nSstart:nSstop) )
         bytes_allocated = bytes_allocated+ ((nZmaxA/2+1)*(nSstop-nSstart+1)* &
         &                (1+2*lm_max))*SIZEOF_DEF_REAL + &
         &                nZmaxA*(nSstop-nSstart+1)*SIZEOF_DEF_REAL
         allocate( PlmZ(l_max+1,nZmaxA/2+1,nSstart:nSstop) )
         allocate( dPlmZ(l_max+1,nZmaxA/2+1,nSstart:nSstop) )
         bytes_allocated = bytes_allocated + 2*(l_max+1)*(nZmaxA/2+1)* &
         &                 (nSstop-nSstart+1)*SIZEOF_DEF_REAL
         allocate( VorOld(nrp_geos,nZmaxA,nSstart:nSstop) )
         bytes_allocated = bytes_allocated + nrp_geos*nZmaxA*(nSstop-nSstart+1)* &
         &                 SIZEOF_DEF_REAL

         allocate( nZC_Sloc(nSstart:nSstop),nZ2(nZmaxA,nSstart:nSstop) )
         bytes_allocated = bytes_allocated + (nSstop-nSstart+1)*(1+nZmaxA)* &
         &                 SIZEOF_INTEGER
         if ( rank == 0 ) then
            allocate (nZC(nSmax))
            bytes_allocated = bytes_allocated+nSmax*SIZEOF_INTEGER
         else
            allocate (nZC(1))
         end if
      end if

      if ( l_geos ) then
         allocate( OsinTS(nZmaxA/2+1,nSstart:nSstop) )
         allocate( PlmS(lm_max,nZmaxA/2+1,nSstart:nSstop) )
         allocate( dPlmS(lm_max,nZmaxA/2+1,nSstart:nSstop) )
         allocate( zZ(nZmaxA,nSstart:nSstop) )
         allocate( rZ(nZmaxA,nSstart:nSstop) )
         bytes_allocated = bytes_allocated+ ((nZmaxA/2+1)*(nSstop-nSstart+1)* &
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
   subroutine finalize_geos_mod(l_geos, l_PV)
      !
      ! Memory deallocation
      !
      logical, intent(in) :: l_geos
      logical, intent(in) :: l_PV

      deallocate( wS_global, dwS_global, ddwS_global, zS_global, dzS_global)

      if ( l_geos ) then
         deallocate( OsinTS, PlmS, dPlmS, zZ, rZ )
         deallocate( chebt_Z, nZmaxS )
         if ( rank == 0 .and. (.not. l_save_out) ) close(n_geos_file)
      end if

      if ( l_PV ) then
         deallocate( PlmZ, dPlmZ, VorOld, nZC_Sloc, nZC, nZ2 )
         deallocate( OsinTSPV, PlmSPV, dPlmSPV, rZPV )
      end if

      deallocate( cyl_balance )

   end subroutine finalize_geos_mod
!----------------------------------------------------------------------------
   subroutine getEgeos(time,nGeosSets,w,dw,ddw,z,dz,Geos,dpFlow,dzFlow)
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
      real(cp), intent(out) :: Geos ! Degree of geostrophy
      real(cp), intent(out) :: dpFlow ! RMS lengths scale
      real(cp), intent(out) :: dzFlow ! RMS lengths scale

      !-- Local variables:
      real(cp) :: Egeos,EkNTC,EkSTC,Ekin
      real(cp) :: CVzOTC,CVorOTC,CHelOTC
      logical :: lDeriv
      integer :: nS,nS_ICB,kindCalc
      real(cp) :: zNorm          ! Norm z interval
      integer :: nNorm           ! No. of grid points for norm interval
      real(cp) :: zMin,zMax,help ! integration boundarie, help variable
      real(cp) :: sZ(nSmax),dsZ ! cylindrical radius s and s-step
      integer :: nPhi,nI
      real(cp) :: phiNorm
      logical :: lTC


      !-- Representation in (phi,z):
      real(cp) :: VrS(nrp_geos,nZmaxA),VrInt(nZmaxA),VrIntS
      real(cp) :: VtS(nrp_geos,nZmaxA),VtInt(nZmaxA),VtIntS
      real(cp) :: VpS(nrp_geos,nZmaxA),VpInt(nZmaxA),VpIntS
      real(cp) :: VozS(nrp_geos,nZmaxA)
      real(cp) :: sinT,cosT
      integer :: nInt,nInts   ! index for NHS and SHS integral
      integer :: nZ,nZmax,nZS,nZN
      real(cp) :: EkInt(nZmaxA),EkIntS
      real(cp) :: EkSTC_s(nSstart:nSstop),EkNTC_s(nSstart:nSstop)
      real(cp) :: Egeos_s(nSstart:nSstop),EkOTC_s(nSstart:nSstop)
      real(cp) :: EkOTC
      real(cp) :: dpEkInt(nZmaxA),dpEkIntS
      real(cp) :: dzEkInt(nZmaxA),dzEkIntS
      real(cp) :: dpEk_s(nSstart:nSstop),dzEk_s(nSstart:nSstop)
      real(cp) :: thetaZ

      !-- Correlation (new Oct. 4 2007)
      logical :: lCorrel
      real(cp) :: VzS,VzN,VorS,VorN,surf,delz
      real(cp) :: VzSN,VzSS,VzNN,VorSN,VorSS,VorNN,HelZZ,VZZ,VorZZ
      real(cp) :: CVz_I,CVor_I,CHel_I
      real(cp) :: CVz_s(nSstart:nSstop),CVor_s(nSstart:nSstop),CHel_S(nSstart:nSstop)

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
      lDeriv=.true.
      kindCalc=1

      !---- Get resolution in s and z for z-integral:
      zNorm=one               ! This is r_CMB-r_ICB
      nNorm=int(zDens*n_r_max) ! Covered with nNorm  points !
      dsZ  =r_CMB/real(nSmax,cp)  ! Step in s controlled by nSmax

      !-- Global array
      do nS=1,nSmax
         sZ(nS)=(nS-half)*dsZ
      end do

      do nS=nSstart,nSstop
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

      !---- Contributions are now in fully spectral space!
      !---- Do the z-integral:
      nI=0

      do nS=nSstart,nSstop
         if ( sZ(nS) < r_ICB ) then
            lTC=.true.
         else
            lTC=.false.
         end if
         if ( nS > nSstart ) then
            if ( sZ(nS-1) < r_ICB .and. sZ(nS) >= r_ICB ) nS_ICB=nS
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
                 &           nZmaxA,zZ(:,nS),nZmaxS(nS),chebt_Z(nS))
            !------ Calculate and store 1/sin(theta) and Plms,dPlms for
            !       southern HS:

            if ( lTC ) then
               nZmax=nZmaxS(nS)  ! nZmax point in each polar region
            else
               nZmax=(nZmaxS(nS)-1)/2+1 ! Odd point number !
               ! all together nZmaxS(nS) from
               ! south to north including equator
            end if
            do nZ=1,nZmax
               rZ(nZ,nS)    =sqrt(zZ(nZ,nS)**2+sZ(nS)**2)
               thetaZ       =atan2(sZ(nS),zZ(nZ,nS))
               OsinTS(nZ,nS)=one/sin(thetaZ)
               call plm_theta(thetaZ,l_max,m_max,minc,                &
                    &            PlmS(:,nZ,nS),dPlmS(:,nZ,nS),lm_max,2)
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

         call getDVptr(wS_global,dwS_global,ddwS_global,zS_global,dzS_global, &
              &        r_ICB,r_CMB,rZ(:,nS),nZmax,nZmaxA,PlmS(:,:,nS),        &
              &        dPlmS(:,:,nS),OsinTS(:,nS),kindCalc,VrS,VtS,VpS,VozS,  &
              &        dpEkInt)

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
               Egeos_s(nS)=Egeos_s(nS) +                             &
               &           half*phiNorm*(zMax-zMin)*sZ(nS)*dsZ *     &
               &           (VrIntS**2+VtIntS**2+VpIntS**2)
               if ( lTC ) then
                  if ( nInt == 1 ) then
                     EkSTC_s(nS)=EkSTC_s(nS) +                        &
                     &           half*phiNorm*(zMax-zMin)*sZ(nS)*dsZ*EkIntS
                  else if ( nInt == 2 ) then
                     EkNTC_s(nS)=EkNTC_s(nS) +                        &
                     &           half*phiNorm*(zMax-zMin)*sZ(nS)*dsZ*EkIntS
                  end if
               else
                  EkOTC_s(nS)=EkOTC_s(nS) +                           &
                  &           half*phiNorm*(zMax-zMin)*sZ(nS)*dsZ*EkIntS
               end if

               !-------- Note: chebIntD returns the z derivative
               if ( lDeriv ) then
                  do nZ=1,nZmax
                     dzEkInt(nZ)=VrInt(nZ)*VrInt(nZ) +                &
                     &           VtInt(nZ)*VtInt(nZ) + VpInt(nZ)*VpInt(nZ)
                  end do
                  dzEkIntS=chebInt(dzEkInt,zMin,zMax,nZmax,nZmaxA,chebt_Z(nS))
                  dzEk_s(nS)=dzEk_s(nS) +                             &
                  &          half*phiNorm*(zMax-zMin)*sZ(nS)*dsZ*dzEkIntS
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
               CVz_Sloc(nPhi,nS) =real(CVz_I,kind=outp)
               CVor_Sloc(nPhi,nS)=real(CVor_I,kind=outp)
               CHel_Sloc(nPhi,nS)=real(CHel_I,kind=outp)
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
               dpEk_s(nS)=dpEk_s(nS) +                               &
               &          half*(zMax-zMin)*sZ(nS)*dsZ*dpEkIntS /     &
               &          (sZ(nS)**2) ! Convert angle to length
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
      do nS=nSstart,nSstop
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

#ifdef WITH_MPI
      call MPI_Allreduce(MPI_IN_PLACE, EkSTC, 1, MPI_DEF_REAL, MPI_SUM,  &
           &             MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, EkNTC, 1, MPI_DEF_REAL, MPI_SUM,  &
           &             MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, EkOTC, 1, MPI_DEF_REAL, MPI_SUM,  &
           &             MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, Egeos, 1, MPI_DEF_REAL, MPI_SUM,  &
           &             MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, CVzOTC, 1, MPI_DEF_REAL, MPI_SUM, &
           &             MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, CVorOTC, 1, MPI_DEF_REAL, MPI_SUM,&
           &             MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, CHelOTC, 1, MPI_DEF_REAL, MPI_SUM,&
           &             MPI_COMM_WORLD, ierr)
#endif

      if ( lCorrel ) then
         surf=0.0_cp
         do nS=nSstart,nSstop
            if ( sZ(nS) >= r_ICB ) surf=surf+sZ(nS)*dsZ
         end do
#ifdef WITH_MPI
         call MPI_Allreduce(MPI_IN_PLACE, surf, 1, MPI_DEF_REAL, MPI_SUM,  &
              &             MPI_COMM_WORLD, ierr)
#endif

         surf   =two*pi*surf
         CVzOTC =CVzOTC/surf
         CVorOTC=CVorOTC/surf
         CHelOTC=CHelOTC/surf
      end if
      Ekin=EkSTC+EkNTC+EKOTC ! Can be used for testing

      dpFlow=0.0_cp
      dzFlow=0.0_cp
      if ( lDeriv ) then
         do nS=nSstart,nSstop
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
#ifdef WITH_MPI
      call MPI_Allreduce(MPI_IN_PLACE, dpFlow, 1, MPI_DEF_REAL, MPI_SUM,  &
           &             MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, dzFlow, 1, MPI_DEF_REAL, MPI_SUM,  &
           &             MPI_COMM_WORLD, ierr)
#endif

      if ( Ekin > 0.0_cp ) then
         Geos = Egeos/Ekin ! relative geostrophic kinetic energy
      else
         Geos = 0.0_cp
         Ekin = -one
      end if

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

         write(n_geos_file,'(1P,ES20.12,7ES16.8)')       &
         &     time, Geos, EkNTC/Ekin, EkSTC/Ekin, Ekin, &
         &     CVzOTC, CVorOTC, CHelOTC

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

            write(nOutFile) ((CVz(nPhi,nS) ,nPhi=1,n_phi_max), nS=nSmax,1,-1)
            write(nOutFile) ((CVor(nPhi,nS),nPhi=1,n_phi_max), nS=nSmax,1,-1)
            write(nOutFile) ((CHel(nPhi,nS),nPhi=1,n_phi_max), nS=nSmax,1,-1)

            close(nOutFile)

         end if ! l_corrMov

      end if ! rank == 0

      if ( lVerbose ) write(*,*) '! End of getGeos!'

   end subroutine getEgeos
!----------------------------------------------------------------------------
   subroutine outPV(time,l_stop_time,nPVsets,w,dw,ddw,z,dz,omega_IC,omega_MA)
      !
      !   Output of z-integrated axisymmetric rotation rate Vp/s
      !   and s derivatives
      !

      !-- Input of variables:
      real(cp),    intent(in) :: time
      real(cp),    intent(in) :: omega_IC,omega_MA
      logical,     intent(in) :: l_stop_time
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dw(llm:ulm,n_r_max)
      complex(cp), intent(in) :: ddw(llm:ulm,n_r_max)
      complex(cp), intent(in) :: z(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dz(llm:ulm,n_r_max)

      integer, intent(inout) :: nPVsets

      !-- (l,r) Representation of the different contributions
      real(cp) :: dzVpLMr(l_max+1,n_r_max), dzVpLMr_loc(l_max+1,n_r_max)

      !--- Work array:
      real(cp) :: workAr(lm_max,n_r_max)

      integer :: lm, l, m, kindCalc

      real(cp) :: fac

      !--- define Grid
      integer :: nS
      real(cp) :: sZ(nSmax),dsZ ! cylindrical radius s and s-step

      integer :: nZ,nZmaxNS
      integer, save :: nZS
      real(cp) :: zZ(nZmaxA),zstep!,zZC
      real(cp) :: VpAS(nZmaxA),omS(nZmaxA,nSmax)

      !-- Local arrays
      real(cp) :: omS_Sloc(nZmaxA,nSstart:nSstop)
      real(outp) :: frame_Sloc(5,n_phi_max*nZmaxA,nSstart:nSstop)

      !-- Plms: Plm,sin
      integer :: nC,nR,nPhi
      real(cp) :: thetaZ,rZS!,sinT,cosT

      !-- For PV output files:
      character(len=80) :: fileName

      !-- Output of all three field components:
      real(cp) :: VsS(nrp_geos,nZmaxA)
      real(cp) :: VpS(nrp_geos,nZmaxA)
      real(cp) :: VzS(nrp_geos,nZmaxA)
      real(cp) :: VorS(nrp_geos,nZmaxA)
      real(outp) :: frame(5,n_phi_max*nZmaxA,nSmax)

      integer :: n_pvz_file, n_vcy_file

#ifdef WITH_MPI
      integer :: i,sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1)
#endif

      if ( lVerbose ) write(*,*) '! Starting outPV!'

      kindCalc = 2

      call costf_arrays(w,dw,ddw,z,dz)

      if ( l_stop_time ) then
         if ( l_SRIC  .and. omega_IC /= 0 ) then
            fac=one/omega_IC
         else
            fac=one
         end if
         do nR=1,n_r_max
            dzVpLMr_loc(:,nR)=0.0_cp
            do lm=llm,ulm
               l=lo_map%lm2l(lm)
               m=lo_map%lm2m(lm)
               if ( m == 0 ) dzVpLMr_loc(l+1,nR)=fac*orho1(nR)*real(z(lm,nR))
            end do
#ifdef WITH_MPI
            call MPI_Allreduce(dzVpLMr_loc(:,nR), dzVPLMr(:,nR), l_max+1, &
                 &             MPI_DEF_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
            dzVPLMr(:,nR)=dzVpLMr_loc(:,nR)
#endif
         end do

         !---- Transform the contributions to cheb space:
         call rscheme_oc%costf1(dzVpLMr,l_max+1,1,l_max+1,workAr)

      end if


      nPVsets=nPVsets+1

      !-- Global arrays
      dsZ=r_CMB/real(nSmax,kind=cp)  ! Step in s controlled by nSmax
      do nS=1,nSmax
         sZ(nS)=(nS-half)*dsZ
      end do
      zstep=2*r_CMB/real(nZmaxA-1,kind=cp)
      do nZ=1,nZmaxA
         zZ(nZ)=r_CMB-(nZ-1)*zstep
      end do


      do nS=nSstart,nSstop

         !------ Get r,theta,Plm,dPlm for northern hemishere:
         if ( nPVsets == 1 ) then ! do this only for the first call !
            nZC_Sloc(nS)=0 ! Points within shell
            do nZ=1,nZmaxA
               rZS=sqrt(zZ(nZ)**2+sZ(nS)**2)
               if ( rZS >= r_ICB .and. rZS <= r_CMB ) then
                  nZC_Sloc(nS)=nZC_Sloc(nS)+1  ! Counts all z within shell
                  nZ2(nZ,nS)=nZC_Sloc(nS) ! No of point within shell
                  if ( zZ(nZ) > 0 ) then ! Onl north hemisphere
                     rZPV(nZC_Sloc(nS),nS)=rZS
                     thetaZ=atan2(sZ(nS),zZ(nZ))
                     OsinTSPV(nZC_Sloc(nS),nS)=one/sin(thetaZ)
                     call plm_theta(thetaZ,l_max,0,minc,PlmZ(:,nZC_Sloc(nS),nS),&
                          &         dPlmZ(:,nZC_Sloc(nS),nS),l_max+1,2)
                     call plm_theta(thetaZ,l_max,m_max,minc,            &
                          &         PlmSPV(:,nZC_Sloc(nS),nS),          &
                          &         dPlmSPV(:,nZC_Sloc(nS),nS),lm_max,2)
                  end if
               else
                  nZ2(nZ,nS)=-1 ! No z found within shell !
               end if
            end do
         end if

         !-- Get azimuthal flow component in the shell
         nZmaxNS=nZC_Sloc(nS) ! all z points within shell
         if ( l_stop_time ) then
            call getPAStr(VpAS,dzVpLMr,nZmaxNS,nZmaxA,l_max,r_ICB,r_CMB,n_r_max, &
                 &        rZPV(:,nS),dPlmZ(:,:,nS),OsinTSPV(:,nS))

            !-- Copy to array with all z-points
            do nZ=1,nZmaxA
               rZS=sqrt(zZ(nZ)**2+sZ(nS)**2)
               nZS=nZ2(nZ,nS)
               if ( nZS > 0 ) then
                  omS_Sloc(nZ,nS)=VpAS(nZS)/sZ(nS)
               else
                  if ( rZS <= r_ICB ) then
                     omS_Sloc(nZ,nS)=one
                  else
                     omS_Sloc(nZ,nS)=fac*omega_MA
                  end if
               end if
            end do
         end if

         !-- Get all three components in the shell
         call getDVptr(wS_global,dwS_global,ddwS_global,zS_global,dzS_global,   &
              &        r_ICB,r_CMB,rZPV(:,nS),nZmaxNS,nZmaxA,PlmSPV(:,:,nS),    &
              &        dPlmSPV(:,:,nS),OsinTSPV(:,nS),kindCalc,VsS,VzS,VpS,VorS)

         if ( l_stop_time ) then
            nC=0
            do nZ=1,nZmaxNS
               do nPhi=1,n_phi_max
                  nC=nC+1
                  frame_Sloc(1,nC,nS)=real(VsS(nPhi,nZ),kind=outp) ! Vs
                  frame_Sloc(2,nC,nS)=real(VpS(nPhi,nZ),kind=outp) ! Vphi
                  frame_Sloc(3,nC,nS)=real(VzS(nPhi,nZ),kind=outp) ! Vz
                  frame_Sloc(4,nC,nS)=real(VorS(nPhi,nZ),kind=outp)
                  frame_Sloc(5,nC,nS)=(real(VorS(nPhi,nZ)-             &
                  &                    VorOld(nPhi,nZ,nS),kind=outp))/ &
                  &                   (real(time-timeOld,kind=outp))
               end do
            end do
         else
            timeOld=time
            do nZ=1,nZmaxNS
               do nPhi=1,n_phi_max
                  VorOld(nPhi,nZ,nS)=VorS(nPhi,nZ)
               end do
            end do
         end if

      end do  ! Loop over s

      if ( l_stop_time ) then

#ifdef WITH_MPI
         sendcount  = (nSstop-nSstart+1)*nZmaxA
         do i=0,n_procs-1
            recvcounts(i) = cyl_balance(i)%n_per_rank*nZmaxA
         end do
         displs(0)=0
         do i=1,n_procs-1
            displs(i) = displs(i-1)+recvcounts(i-1)
         end do
         call MPI_GatherV(omS_Sloc, sendcount, MPI_DEF_REAL,      &
              &           omS, recvcounts, displs, MPI_DEF_REAL,  &
              &           0, MPI_COMM_WORLD, ierr)

         sendcount  = (nSstop-nSstart+1)*nZmaxA*n_phi_max*5
         do i=0,n_procs-1
            recvcounts(i) = cyl_balance(i)%n_per_rank*nZmaxA*n_phi_max*5
         end do
         displs(0)=0
         do i=1,n_procs-1
            displs(i) = displs(i-1)+recvcounts(i-1)
         end do
         call MPI_GatherV(frame_Sloc, sendcount, MPI_OUT_REAL,    &
              &           frame, recvcounts, displs, MPI_OUT_REAL,&
              &           0, MPI_COMM_WORLD, ierr)

         sendcount  = (nSstop-nSstart+1)
         do i=0,n_procs-1
            recvcounts(i) = cyl_balance(i)%n_per_rank
         end do
         displs(0)=0
         do i=1,n_procs-1
            displs(i) = displs(i-1)+recvcounts(i-1)
         end do
         call MPI_GatherV(nZC_Sloc, sendcount, MPI_INTEGER,       &
              &           nZC, recvcounts, displs, MPI_INTEGER,   &
              &           0, MPI_COMM_WORLD, ierr)
#else
         omS(:,:)    =omS_Sloc(:,:)
         frame(:,:,:)=frame_Sloc(:,:,:)
         nZC(:)      =nZC_Sloc(:)
#endif
         !-- Write output only at the final timestep
         if ( rank == 0 ) then

            !--- Open file for output:
            fileName='PVZ.'//tag
            open(newunit=n_pvz_file, file=fileName, form='unformatted', &
            &    status='unknown')
            write(n_pvz_file) real(time,kind=outp), real(nSmax,kind=outp), &
            &     real(nZmaxA,kind=outp), real(omega_IC,kind=outp),        &
            &     real(omega_ma,kind=outp)
            write(n_pvz_file) (real(sZ(nS),kind=outp),nS=1,nSmax)
            write(n_pvz_file) (real(zZ(nZ),kind=outp),nZ=1,nZmaxA)


            !--- Open file for the three flow components:
            fileName='Vcy.'//tag
            open(newunit=n_vcy_file, file=fileName,form='unformatted', &
            &    status='unknown')
            write(n_vcy_file) real(time,kind=outp), real(nSmax,kind=outp),&
            &     real(nZmaxA,kind=outp), real(n_phi_max,kind=outp),      &
            &     real(omega_IC,kind=outp), real(omega_ma,kind=outp),     &
            &     real(radratio,kind=outp), real(minc,kind=outp)
            write(n_vcy_file) (real(sZ(nS),kind=outp),nS=1,nSmax)
            write(n_vcy_file) (real(zZ(nZ),kind=outp),nZ=1,nZmaxA)

            do nS=1,nSmax

               nZmaxNS=nZC(nS) ! all z points within shell
               nC = n_phi_max*nZmaxNS

               write(n_pvz_file) (real(omS(nZ,nS),kind=outp),nZ=1,nZmaxA)

               write(n_vcy_file) real(nZmaxNS,kind=outp)
               write(n_vcy_file) (frame(1,nZ,nS),nZ=1,nC)
               write(n_vcy_file) (frame(2,nZ,nS),nZ=1,nC)
               write(n_vcy_file) (frame(3,nZ,nS),nZ=1,nC)
               write(n_vcy_file) (frame(4,nZ,nS),nZ=1,nC)
               write(n_vcy_file) (frame(5,nZ,nS),nZ=1,nC)
            end do

            close(n_pvz_file)
            close(n_vcy_file)

         end if ! Rank 0

      end if ! l_stop_time

   end subroutine outPV
!---------------------------------------------------------------------------------
   subroutine getDVptr(w,dw,ddw,z,dz,rMin,rMax,rS,nZmax,nZmaxA,PlmS,dPlmS, &
              &        OsinTS,kindCalc,VrS,VtS,VpS,VorS,dpEk)
      !
      !  This subroutine calculates the three flow components VrS,VtS,VpS at
      !  a (r,theta,all phis) and (t,pi=theta, all phis). Here r=rS, PlmS=Plm(theta),
      !  dPlmS=sin(theta)*dTheta Plm(theta), and OsinTS=1/sin(theta).
      !  The flow is calculated for all n_phi_max azimuthal points used in the code,
      !  and for corresponding latitudes north and south of the equator.
      !  When kindCalc=1, the subroutine also calculates dpEk and dzEk which
      !  are phi averages of (d Vr/d phi)**2 + (d Vtheta/ d phi)**2 + (d Vphi/ d phi)**2
      !  and (d Vr/d z)**2 + (d Vtheta/ d z)**2 + (d Vphi/ d z)**2, respectively.
      !
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
      real(cp),    intent(in) :: PlmS(lm_max,nZmaxA/2+1)
      real(cp),    intent(in) :: dPlmS(lm_max,nZmaxA/2+1)
      real(cp),    intent(in) :: OsinTS(nZmaxA/2+1)
      integer,     intent(in) :: kindCalc

      !--- Output: function on azimuthal grid points defined by FT!
      real(cp), intent(out) :: VrS(nrp_geos,nZmaxA)
      real(cp), intent(out) :: VtS(nrp_geos,nZmaxA)
      real(cp), intent(out) :: VpS(nrp_geos,nZmaxA)
      real(cp), intent(out) :: VorS(nrp_geos,nZmaxA)
      real(cp), optional, intent(out) :: dpEk(nZmaxA)

      !--- Local variables:
      real(cp) :: chebS(n_r_max)
      integer :: nS,nN,mc,lm,l,m,nCheb,nPhi,n
      real(cp) :: x,phiNorm,mapFac,OS,cosT,sinT,Or_e1,Or_e2
      complex(cp) :: Vr,Vt1,Vt2,Vp1,Vp2,Vor,Vot1,Vot2
      real(cp) :: VotS(nrp_geos,nZmaxA)
      complex(cp) :: wSr,dwSr,ddwSr,zSr,dzSr

      real(cp) :: dV(nrp_geos,nZmaxA)
      complex(cp) :: dp


      mapFac=two/(rMax-rMin)
      phiNorm=two*pi/n_phi_max

      do nS=1,nZmax
         do mc=1,nrp_geos
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
         chebS(1)=one*rscheme_oc%rnorm ! Extra cheb_norm cheap here
         chebS(2)=x*rscheme_oc%rnorm
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
               wSr  =  wSr+  w(lm,nCheb)*chebS(nCheb)
               dwSr = dwSr+ dw(lm,nCheb)*chebS(nCheb)
               ddwSr=ddwSr+ddw(lm,nCheb)*chebS(nCheb)
               zSr  =  zSr+  z(lm,nCheb)*chebS(nCheb)
               dzSr = dzSr+ dz(lm,nCheb)*chebS(nCheb)
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
            if ( kindCalc == 1 ) then
               VrS(2*mc-1,nS) =Or_e2*VrS(2*mc-1,nS)
               VrS(2*mc  ,nS) =Or_e2*VrS(2*mc  ,nS)
               VtS(2*mc-1,nS) =Or_e1*OS*VtS(2*mc-1,nS)
               VtS(2*mc  ,nS) =Or_e1*OS*VtS(2*mc  ,nS)
            else if ( kindCalc == 2 ) then
               !-- This is now Vs
               VrS(2*mc-1,nS)=sinT*Or_e2*VrS(2*mc-1,nS)+cosT*Or_e1*OS*VtS(2*mc-1,nS)
               VrS(2*mc  ,nS)=sinT*Or_e2*VrS(2*mc  ,nS)+cosT*Or_e1*OS*VtS(2*mc  ,nS)
               !-- This is now Vz
               VtS(2*mc-1,nS)=cosT*Or_e2*VrS(2*mc-1,nS)-sinT*Or_e1*OS*VtS(2*mc-1,nS)
               VtS(2*mc  ,nS)=cosT*Or_e2*VrS(2*mc  ,nS)-sinT*Or_e1*OS*VtS(2*mc  ,nS)
            end if
            VpS(2*mc-1,nS) =Or_e1*OS*VpS(2*mc-1,nS)
            VpS(2*mc  ,nS) =Or_e1*OS*VpS(2*mc  ,nS)
            VorS(2*mc-1,nS)=cosT*Or_e2*VorS(2*mc-1,nS) - Or_e1*VotS(2*mc-1,nS)
            VorS(2*mc  ,nS)=cosT*Or_e2*VorS(2*mc  ,nS) - Or_e1*VotS(2*mc  ,nS)
         end do
         do mc=2*n_m_max+1,nrp_geos
            VrS(mc,nS) =0.0_cp
            VpS(mc,nS) =0.0_cp
            VtS(mc,nS) =0.0_cp
            VorS(mc,nS)=0.0_cp
         end do
      end do

      do nS=(nZmax+1)/2+1,nZmax ! South HS
         OS   =OsinTS(nZmax-nS+1)
         sinT =one/OS
         cosT =-sqrt(one-sinT**2)
         Or_e1=one/rS(nZmax-nS+1)
         Or_e2=Or_e1*Or_e1
         do mc=1,n_m_max
            if ( kindCalc == 1) then
               VrS(2*mc-1,nS) =Or_e2*VrS(2*mc-1,nS)
               VrS(2*mc  ,nS) =Or_e2*VrS(2*mc  ,nS)
               VtS(2*mc-1,nS) =Or_e1*OS*VtS(2*mc-1,nS)
               VtS(2*mc  ,nS) =Or_e1*OS*VtS(2*mc  ,nS)
            else if ( kindCalc == 2 ) then
               !-- This is now Vs
               VrS(2*mc-1,nS)=sinT*Or_e2*VrS(2*mc-1,nS)+cosT*Or_e1*OS*VtS(2*mc-1,nS)
               VrS(2*mc  ,nS)=sinT*Or_e2*VrS(2*mc  ,nS)+cosT*Or_e1*OS*VtS(2*mc  ,nS)
               !-- This is now Vz
               VtS(2*mc-1,nS)=cosT*Or_e2*VrS(2*mc-1,nS)-sinT*Or_e1*OS*VtS(2*mc-1,nS)
               VtS(2*mc  ,nS)=cosT*Or_e2*VrS(2*mc  ,nS)-sinT*Or_e1*OS*VtS(2*mc  ,nS)
            end if
            VpS(2*mc-1,nS) =Or_e1*OS*VpS(2*mc-1,nS)
            VpS(2*mc  ,nS) =Or_e1*OS*VpS(2*mc  ,nS)
            VorS(2*mc-1,nS)=cosT*Or_e2*VorS(2*mc-1,nS) - Or_e1*VotS(2*mc-1,nS)
            VorS(2*mc  ,nS)=cosT*Or_e2*VorS(2*mc  ,nS) - Or_e1*VotS(2*mc  ,nS)
         end do
         do mc=2*n_m_max+1,nrp_geos
            VrS(mc,nS) =0.0_cp
            VpS(mc,nS) =0.0_cp
            VtS(mc,nS) =0.0_cp
            VorS(mc,nS)=0.0_cp
         end do
      end do

      if ( present(dpEk) ) then
         !--- Calculate phi derivative in lm-space:
         do n=1,3
            do nS=1,nZmax
               if ( n == 1 ) then
                  dpEk(nS)=0.0_cp
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
      if ( .not. l_axi ) then
         call fft_to_real(VrS,nrp_geos,nZmax)
         call fft_to_real(VtS,nrp_geos,nZmax)
         call fft_to_real(VpS,nrp_geos,nZmax)
         call fft_to_real(VorS,nrp_geos,nZmax)
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
