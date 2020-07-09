module torsional_oscillations
   !
   !  This module contains information for TO calculation and output
   !

   use iso_fortran_env, only: output_unit
   use precision_mod
   use parallel_mod
   use mem_alloc, only: bytes_allocated
   use LMmapping, only: map_glbl_st, map_dist_st
   use truncation, only: nrp, n_phi_maxStr, n_r_maxStr, l_max, n_theta_maxStr, &
       &                 n_r_cmb, nRstart, nRstop, nThetaStart, nThetaStop,    &
       &                 n_lm_loc, m_tsid
   use radial_functions, only: r, or1, or2, or3, or4, beta, orho1, dbeta
   use physical_parameters, only: CorFac, kbotv, ktopv
   use horizontal_data, only: sinTheta, cosTheta, hdif_V, dTheta1A, dTheta1S
   use constants, only: one, two
   use logic, only: lVerbose, l_mag
   use sht, only: spat_to_SH_axi

   implicit none

   private

   real(cp), public, allocatable :: dzStrLMr(:,:)
   real(cp), public, allocatable :: dzRstrLMr(:,:)
   real(cp), public, allocatable :: dzAstrLMr(:,:)
   real(cp), public, allocatable :: dzCorLMr(:,:)
   real(cp), public, allocatable :: dzLFLMr(:,:)
   real(cp), public, allocatable :: dzdVpLMr(:,:)
   real(cp), public, allocatable :: dzddVpLMr(:,:)
   real(cp), public, allocatable :: ddzASL(:,:)
   real(cp), allocatable :: dzASL(:)
   real(cp), allocatable :: zASL(:)

   real(cp), public, allocatable :: dzStrLMr_Rloc(:,:)
   real(cp), public, allocatable :: dzRstrLMr_Rloc(:,:)
   real(cp), public, allocatable :: dzAStrLMr_Rloc(:,:)
   real(cp), public, allocatable :: dzCorLMr_Rloc(:,:)
   real(cp), public, allocatable :: dzLFLMr_Rloc(:,:)
   real(cp), public, allocatable :: dzdVpLMr_Rloc(:,:)
   real(cp), public, allocatable :: dzddVpLMr_Rloc(:,:)
   real(cp), public, allocatable :: V2AS_Rloc(:,:)
   real(cp), public, allocatable :: Bs2AS_Rloc(:,:)
   real(cp), public, allocatable :: BszAS_Rloc(:,:)
   real(cp), public, allocatable :: BspAS_Rloc(:,:)
   real(cp), public, allocatable :: BpzAS_Rloc(:,:)
   real(cp), public, allocatable :: BspdAS_Rloc(:,:)
   real(cp), public, allocatable :: BpsdAS_Rloc(:,:)
   real(cp), public, allocatable :: BzpdAS_Rloc(:,:)
   real(cp), public, allocatable :: BpzdAS_Rloc(:,:)

   real(cp), allocatable :: BsLast(:,:,:), BpLast(:,:,:), BzLast(:,:,:)

   public :: initialize_TO, prep_TO_axi, getTO, getTOnext, getTOfinish, finalize_TO

contains

   subroutine initialize_TO
      !
      ! Allocate the memory needed
      !
      allocate( dzStrLMr(l_max+1,n_r_maxStr) )
      allocate( dzRstrLMr(l_max+1,n_r_maxStr) )
      allocate( dzAstrLMr(l_max+1,n_r_maxStr) )
      allocate( dzCorLMr(l_max+1,n_r_maxStr) )
      allocate( dzLFLMr(l_max+1,n_r_maxStr) )
      allocate( dzdVpLMr(l_max+1,n_r_maxStr) )
      allocate( dzddVpLMr(l_max+1,n_r_maxStr) )
      bytes_allocated = bytes_allocated+7*n_r_maxStr*(l_max+1)*SIZEOF_DEF_REAL

      allocate( dzStrLMr_Rloc(l_max+1,nRstart:nRstop) )
      allocate( dzRstrLMr_Rloc(l_max+1,nRstart:nRstop) )
      allocate( dzAStrLMr_Rloc(l_max+1,nRstart:nRstop) )
      allocate( dzCorLMr_Rloc(l_max+1,nRstart:nRstop) )
      allocate( dzLFLMr_Rloc(l_max+1,nRstart:nRstop) )
      allocate( dzdVpLMr_Rloc(l_max+1,nRstart:nRstop) )
      allocate( dzddVpLMr_Rloc(l_max+1,nRstart:nRstop) )
      bytes_allocated = bytes_allocated+ &
      &                 7*(nRstop-nRstart+1)*(l_max+1)*SIZEOF_DEF_REAL

      allocate( V2AS_Rloc(nThetaStart:nThetaStop,nRstart:nRstop) )
      allocate( Bs2AS_Rloc(nThetaStart:nThetaStop,nRstart:nRstop) )
      allocate( BszAS_Rloc(nThetaStart:nThetaStop,nRstart:nRstop) )
      allocate( BspAS_Rloc(nThetaStart:nThetaStop,nRstart:nRstop) )
      allocate( BpzAS_Rloc(nThetaStart:nThetaStop,nRstart:nRstop) )
      allocate( BspdAS_Rloc(nThetaStart:nThetaStop,nRstart:nRstop) )
      allocate( BpsdAS_Rloc(nThetaStart:nThetaStop,nRstart:nRstop) )
      allocate( BzpdAS_Rloc(nThetaStart:nThetaStop,nRstart:nRstop) )
      allocate( BpzdAS_Rloc(nThetaStart:nThetaStop,nRstart:nRstop) )
      bytes_allocated = bytes_allocated+9*(nRstop-nRstart+1)* &
      &                 (nThetaStop-nThetaStart+1)*SIZEOF_DEF_REAL

      allocate( ddzASL(l_max+1,n_r_maxStr), dzASL(l_max+1), zASL(l_max+1) )
      bytes_allocated = bytes_allocated+ 3*(l_max+1)*n_r_maxStr*SIZEOF_DEF_REAL
      ddzASL(:,:)=0.0_cp
      dzASL(:)   =0.0_cp
      zASL(:)    =0.0_cp

      allocate( BsLast(n_phi_maxStr,nThetaStart:nThetaStop,nRstart:nRstop) )
      allocate( BpLast(n_phi_maxStr,nThetaStart:nThetaStop,nRstart:nRstop) )
      allocate( BzLast(n_phi_maxStr,nThetaStart:nThetaStop,nRstart:nRstop) )
      bytes_allocated = bytes_allocated+3*n_phi_maxStr*               &
      &                 (nThetaStop-nThetaStart+1)*(nRstop-nRstart+1)*&
      &                 SIZEOF_DEF_REAL

   end subroutine initialize_TO
!-----------------------------------------------------------------------------
   subroutine finalize_TO
      !
      ! Deallocate the memory
      !

      deallocate( ddzASL, BpzdAS_Rloc, BzpdAS_Rloc, BpsdAS_Rloc, BspdAS_Rloc )
      deallocate( BpzAS_Rloc, BspAS_Rloc, BszAS_Rloc, Bs2AS_Rloc, V2AS_Rloc )
      deallocate( dzddVpLMr_Rloc, dzdVpLMr_Rloc, dzLFLMr_Rloc, dzCorLMr_Rloc )
      deallocate( dzAStrLMr_Rloc, dzRstrLMr_Rloc, dzStrLMr_Rloc )
      deallocate( dzStrLMr, dzRstrLMr, dzAstrLMr, dzCorLMr, dzLFLMr, dzdVpLMr )
      deallocate( dzddVpLMr, zASL, dzASL,BsLast, BpLast, BzLast )

   end subroutine finalize_TO
!-----------------------------------------------------------------------------
   subroutine prep_TO_axi(z,dz)

      !-- Input variables
      complex(cp), intent(in) :: z(n_lm_loc)
      complex(cp), intent(in) :: dz(n_lm_loc)

      !-- Local variables
      integer :: l, lm, Rq(2)

      if ( map_dist_st%has_m0 ) then
         do l=0,l_max
            lm = map_dist_st%lm2(l,0)
            zASL(l+1)  =real(z(lm))   ! used in TO
            dzASL(l+1) =real(dz(lm))  ! used in TO (anelastic)
         end do
      end if

      ! m_tsid(0) is the rank which has m=0 
      ! funfact: m_tsid = dist_m written backwards! Terrible nomeclature, pardon me
#ifdef WITH_MPI
      call MPI_IBcast(zASL, l_max+1, MPI_DEF_REAL,  m_tsid(0), &
           &          comm_theta, Rq(1), ierr)
      call MPI_IBcast(dzASL, l_max+1, MPI_DEF_REAL,  m_tsid(0), &
           &          comm_theta, Rq(2), ierr)
      call MPI_WaitAll(2, Rq, MPI_STATUSES_IGNORE, ierr)
#endif

   end subroutine prep_TO_axi
!-----------------------------------------------------------------------------
   subroutine getTO(vr,vt,vp,cvr,dvpdr,br,bt,bp,cbr,cbt, &
              &     dzRstrLM,dzAstrLM,dzCorLM,dzLFLM,dtLast,nR)
      !
      !  This program calculates various axisymmetric linear
      !  and nonlinear variables for a radial grid point nR and
      !  a theta-block.
      !  Input are the fields vr,vt,vp,cvr,dvpdr
      !  Output are linear azimuthally averaged  field VpAS (flow phi component),
      !  VpAS2 (square of flow phi component), V2AS (V*V),
      !  and Coriolis force Cor. These are give in (r,theta)-space.
      !  Also in (r,theta)-space are azimuthally averaged correlations of
      !  non-axisymmetric flow components and the respective squares:
      !  Vsp=Vs*Vp,Vzp,Vsz,VspC,VzpC,VszC. These are used to calulcate
      !  the respective correlations and Reynolds stress.
      !  In addition three output field are given in (lm,r) space:
      !  dzRstrLMr,dzAstrLMr,dzCorLM,dzLFLM.
      !
      !  These are used to calculate the total Reynolds stress,
      !  advection and viscous stress later. Their calculation
      !  retraces the calculations done in the time-stepping part
      !  of the code.
      !

      !-- Input of variables
      real(cp), intent(in) :: dtLast              ! last time step
      integer,  intent(in) :: nR                 ! radial grid point
      real(cp), intent(in) :: vr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: vt(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: vp(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: cvr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: dvpdr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: br(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: bt(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: bp(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: cbr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: cbt(nrp,nThetaStart:nThetaStop)

      !-- Output of arrays needing further treatment in getTOfinish:
      real(cp), intent(out) :: dzRstrLM(l_max+2),dzAstrLM(l_max+2)
      real(cp), intent(out) :: dzCorLM(l_max+2),dzLFLM(l_max+2)

      !-- Local variables:
      integer :: nTheta,nPhi
      real(cp) :: VrMean,VtMean,VpMean,Vr2Mean,Vt2Mean,Vp2Mean
      real(cp) :: LFmean,cvrMean,dvpdrMean,VrdVpdrMean,VtcVrMean
      real(cp) :: Bs2Mean,BszMean,BspMean,BpzMean,BspdMean,BpsdMean
      real(cp) :: BzpdMean,BpzdMean
      real(cp) :: Rmean(nThetaStart:nThetaStop),Amean(nThetaStart:nThetaStop)
      real(cp) :: dzCorMean(nThetaStart:nThetaStop),dzLFmean(nThetaStart:nThetaStop)
      real(cp) :: sinT,Osin,Osin2,cosT,cosOsin2,phiNorm
      real(cp) :: BsL,BzL,BpL,Bs2F1,Bs2F2,Bs2F3,BspF1,BspF2
      real(cp) :: BpzF1,BpzF2,BszF1,BszF2,BszF3
      real(cp) :: BsF1,BsF2,BpF1,BzF1,BzF2

      if ( lVerbose ) write(output_unit,*) '! Starting getTO!'

      phiNorm=one/real(n_phi_maxStr, kind=cp)

      !-- Big loop over thetas in block:
      do nTheta=nThetaStart,nThetaStop
         sinT =sinTheta(nTheta)
         cosT =cosTheta(nTheta)
         Osin =one/sinT
         Osin2=Osin*Osin
         cosOsin2=cosT*Osin2
         Bs2F1=sinT*sinT*or4(nR)
         Bs2F2=cosT*cosT*Osin2*or2(nR)
         Bs2F3=two*cosT*or3(nR)
         BspF1=or3(nR)
         BspF2=cosT*Osin2*or2(nR)
         BpzF1=cosT*Osin*or3(nR)
         BpzF2=or2(nR)*Osin
         BszF1=sinT*cosT*or4(nR)
         BszF2=(two*cosT*cosT-one)*Osin*or3(nR)
         BszF3=cosT*Osin*or2(nR)
         BsF1 =sinT*or2(nR)
         BsF2 =cosT*Osin*or1(nR)
         BpF1 =Osin*or1(nR)
         BzF1 =cosT*or2(nR)
         BzF2 =or1(nR)

         !--- Get zonal means of velocity and derivatives:
         VrMean     =0.0_cp
         VtMean     =0.0_cp
         VpMean     =0.0_cp
         Vr2Mean    =0.0_cp
         Vt2Mean    =0.0_cp
         Vp2Mean    =0.0_cp
         dVpdrMean  =0.0_cp
         cVrMean    =0.0_cp
         VrdVpdrMean=0.0_cp
         VtcVrMean  =0.0_cp
         LFmean     =0.0_cp
         Bs2Mean    =0.0_cp
         BspMean    =0.0_cp
         BpzMean    =0.0_cp
         BszMean    =0.0_cp
         BspdMean   =0.0_cp
         BpsdMean   =0.0_cp
         BzpdMean   =0.0_cp
         BpzdMean   =0.0_cp
         do nPhi=1,n_phi_maxStr
            VrMean =VrMean +vr(nPhi,nTheta)
            VtMean =VtMean +vt(nPhi,nTheta)
            VpMean =VpMean +vp(nPhi,nTheta)
            Vr2Mean=Vr2Mean+vr(nPhi,nTheta)*vr(nPhi,nTheta)
            Vt2Mean=Vt2Mean+vt(nPhi,nTheta)*vt(nPhi,nTheta)
            Vp2Mean=Vp2Mean+vp(nPhi,nTheta)*vp(nPhi,nTheta)
            dVpdrMean=dVpdrMean    + orho1(nR)*( dvpdr(nPhi,nTheta) - &
            &                              beta(nR)*vp(nPhi,nTheta) )
            cVrMean  =cVrMean  +cvr(nPhi,nTheta)
            VrdVpdrMean=VrdVpdrMean + orho1(nR)*vr(nPhi,nTheta)*(  &
            &           dvpdr(nPhi,nTheta)-beta(nR)*vp(nPhi,nTheta) )
            VtcVrMean=VtcVrMean    + orho1(nR)*vt(nPhi,nTheta)*cvr(nPhi,nTheta)
            if ( l_mag ) then
               LFmean=LFmean + cbr(nPhi,nTheta)*bt(nPhi,nTheta) - &
               &               cbt(nPhi,nTheta)*br(nPhi,nTheta)
               Bs2Mean=Bs2Mean + Bs2F1*br(nPhi,nTheta)*br(nPhi,nTheta) + &
               &                 Bs2F2*bt(nPhi,nTheta)*bt(nPhi,nTheta) + &
               &                 Bs2F3*br(nPhi,nTheta)*bt(nPhi,nTheta)
               BspMean=BspMean + BspF1*br(nPhi,nTheta)*bp(nPhi,nTheta) + &
               &                 BspF2*bt(nPhi,nTheta)*bp(nPhi,nTheta)
               BpzMean=BpzMean + BpzF1*br(nPhi,nTheta)*bp(nPhi,nTheta) - &
               &                 BpzF2*bt(nPhi,nTheta)*bp(nPhi,nTheta)
               BszMean=BszMean + BszF1*br(nPhi,nTheta)*br(nPhi,nTheta) + &
               &                 BszF2*br(nPhi,nTheta)*bt(nPhi,nTheta) - &
               &                 BszF3*bt(nPhi,nTheta)*bt(nPhi,nTheta)
               BsL=BsF1*br(nPhi,nTheta) + BsF2*bt(nPhi,nTheta)
               BpL=BpF1*bp(nPhi,nTheta)
               BzL=BzF1*br(nPhi,nTheta) - BzF2*bt(nPhi,nTheta)
               BspdMean=BspdMean+BsL*(BpL-BpLast(nPhi,nTheta,nR))
               BpsdMean=BpsdMean+BpL*(BsL-BsLast(nPhi,nTheta,nR))
               BzpdMean=BzpdMean+BzL*(BpL-BpLast(nPhi,nTheta,nR))
               BpzdMean=BpzdMean+BpL*(BzL-BzLast(nPhi,nTheta,nR))
            end if
         end do

         Vr2Mean=phiNorm*or4(nR)*Vr2Mean
         Vt2Mean=phiNorm*or2(nR)*Osin2*Vt2Mean
         Vp2Mean=phiNorm*or2(nR)*Osin2*Vp2Mean
         if ( nR == n_r_CMB ) then
            VrMean =0.0_cp
            Vr2Mean=0.0_cp
            if ( ktopv == 2 ) then
               VtMean =0.0_cp
               Vt2Mean=0.0_cp
               VpMean =0.0_cp
               Vp2Mean=0.0_cp
            end if
         end if
         if ( nR == n_r_CMB ) then
            VrMean =0.0_cp
            Vr2Mean=0.0_cp
            if ( kbotv == 2 ) then
               VtMean =0.0_cp
               Vt2Mean=0.0_cp
               VpMean =0.0_cp
               Vp2Mean=0.0_cp
            end if
         end if
         V2AS_Rloc(nTheta,nR)=Vr2Mean+Vt2Mean+Vp2Mean
         VpMean =phiNorm*or1(nR)*Osin*VpMean
         !--- This is Coriolis force / r*sin(theta)
         dzCorMean(nTheta)=phiNorm*two*CorFac * &
         &                      (or3(nR)*VrMean+or2(nR)*cosOsin2*VtMean)
         if ( l_mag ) then
            !--- This is Lorentz force/ r*sin(theta)
            dzLFmean(nTheta)=phiNorm*or4(nR)*Osin2*LFmean
            Bs2AS_Rloc(nTheta,nR) =phiNorm*Bs2Mean
            BspAS_Rloc(nTheta,nR) =phiNorm*BspMean
            BpzAS_Rloc(nTheta,nR) =phiNorm*BpzMean
            BszAS_Rloc(nTheta,nR) =phiNorm*BszMean
            BspdAS_Rloc(nTheta,nR)=phiNorm*(BspdMean/dtLast)
            BpsdAS_Rloc(nTheta,nR)=phiNorm*(BpsdMean/dtLast)
            BzpdAS_Rloc(nTheta,nR)=phiNorm*(BzpdMean/dtLast)
            BpzdAS_Rloc(nTheta,nR)=phiNorm*(BpzdMean/dtLast)
         end if

         ! dVpdrMean, VtcVrMean and VrdVpdrMean are already divided by rho
         Rmean(nTheta)=phiNorm * or4(nR)*Osin2* (VrdVpdrMean - &
         &                          phiNorm*VrMean*dVpdrMean + &
         &                                         VtcVrMean - &
         &                  orho1(nR)*phiNorm*VtMean*cVrMean )
         Amean(nTheta)=phiNorm*or4(nR)*Osin2*phiNorm*(VrMean*dVpdrMean + &
         &                                  orho1(nR)*VtMean*cVrMean )

      end do ! Loop over thetas in block !

      !--- Transform and Add to LM-Space:
      !------ Add contribution from thetas in block:
      !       note legtfAS2 returns modes l=0 -- l=l_max+1
      call spat_to_SH_axi(Rmean,dzRstrLM)
      call spat_to_SH_axi(Amean,dzAstrLM)
      call spat_to_SH_axi(dzCorMean,dzCorLM)
      call spat_to_SH_axi(dZLFmean,dzLFLM)

      if ( lVerbose ) write(output_unit,*) '! End of getTO!'

   end subroutine getTO
!-----------------------------------------------------------------------------
   subroutine getTOnext(br,bt,bp,lTONext,lTONext2,dt,dtLast,nR)
      !
      !  Preparing TO calculation by storing flow and magnetic field
      !  contribution to build time derivative.
      !

      !-- Input of variables:
      real(cp), intent(in) :: dt,dtLast
      integer,  intent(in) :: nR
      logical,  intent(in) :: lTONext,lTONext2
      real(cp), intent(in) :: br(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: bt(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: bp(nrp,nThetaStart:nThetaStop)

      !-- Local variables:
      integer :: l, lm, nTheta, nPhi
      real(cp) :: sinT,cosT,BsF1,BsF2,BpF1,BzF1,BzF2

      if ( lVerbose ) write(output_unit,*) '! Starting getTOnext!',dtLast

      if ( lTONext2 ) then

         dzddVpLMr_Rloc(1,nR)=0.0_cp
         do l=1,l_max
            lm=map_glbl_st%lm2(l,0)
            dzddVpLMr_Rloc(l+1,nR)=zASL(l+1)
         end do

      else if ( lTOnext ) then

         !$omp parallel do default(shared)                     &
         !$omp& private(nTheta,nPhi,sinT,cosT,BsF1,BsF2,BpF1,BzF1,BzF2)
         do nTheta=nThetaStart,nThetaStop
            sinT =sinTheta(nTheta)
            cosT =cosTheta(nTheta)
            BsF1 =sinT*or2(nR)
            BsF2 =cosT/sinT*or1(nR)
            BpF1 =or1(nR)/sinT
            BzF1 =cosT*or2(nR)
            BzF2 =or1(nR)

            !--- Get zonal means of velocity and derivatives:
            do nPhi=1,n_phi_maxStr
               BsLast(nPhi,nTheta,nR)=BsF1*br(nPhi,nTheta) + BsF2*bt(nPhi,nTheta)
               BpLast(nPhi,nTheta,nR)=BpF1*bp(nPhi,nTheta)
               BzLast(nPhi,nTheta,nR)=BzF1*br(nPhi,nTheta) - BzF2*bt(nPhi,nTheta)
            end do

         end do ! Loop over thetas in block !
         !$omp end parallel do

         dzdVpLMr_Rloc(1,nR) =0.0_cp
         dzddVpLMr_Rloc(1,nR)=0.0_cp
         do l=1,l_max
            lm=map_glbl_st%lm2(l,0)
            dzdVpLMr_Rloc(l+1,nR) = zASL(l+1)
            dzddVpLMr_Rloc(l+1,nR)= ( dzddVpLMr_Rloc(l+1,nR) - &
            &                  ((dtLast+dt)/dt)*zASL(l+1) )/dtLast
         end do

      end if

      if ( lVerbose ) write(output_unit,*) '! End of getTOnext!'

   end subroutine getTOnext
!-----------------------------------------------------------------------------
   subroutine getTOfinish(nR,dtLast,dzRstrLM,dzAstrLM,dzCorLM,dzLFLM)
      !
      !  This program was previously part of getTO(...)
      !  It has now been separated to get it out of the theta-block loop.
      !

      !-- Input of variables:
      integer,  intent(in) :: nR
      real(cp), intent(in) :: dtLast
      real(cp), intent(in) :: dzRstrLM(l_max+2),dzAstrLM(l_max+2)
      real(cp), intent(in) :: dzCorLM(l_max+2),dzLFLM(l_max+2)

      !-- Local variables:
      real(cp) :: dLh
      integer :: l,lS,lA,lm

      !------ When all thetas are done calculate viscous stress in LM space:
      dzStrLMr_Rloc(1,nR) =0.0_cp
      dzRstrLMr_Rloc(1,nR)=0.0_cp
      dzAstrLMr_Rloc(1,nR)=0.0_cp
      dzCorLMr_Rloc(1,nR) =0.0_cp
      dzLFLMr_Rloc(1,nR)  =0.0_cp
      dzdVpLMr_Rloc(1,nR) =0.0_cp
      dzddVpLMr_Rloc(1,nR)=0.0_cp
      do l=1,l_max
         lS=(l-1)+1
         lA=(l+1)+1
         lm=map_glbl_st%lm2(l,0)
         dLh = real(l*(l+1),cp)
         dzStrLMr_Rloc(l+1,nR)= hdif_V(l) * (  ddzASL(l+1,nR) - &
         &                               beta(nR)* dzASL(l+1) - &
         &  (dLh*or2(nR)+dbeta(nR)+two*beta(nR)*or1(nR))*  &
         &                                          zASL(l+1) )
         !---- -r**2/(l(l+1)) 1/sin(theta) dtheta sin(theta)**2
         !     minus sign to bring stuff on the RHS of NS equation !
         dzRstrLMr_Rloc(l+1,nR)=-r(nR)*r(nR)/dLh * ( &
         &                   dTheta1S(lm)*dzRstrLM(lS) - &
         &                   dTheta1A(lm)*dzRstrLM(lA) )
         dzAstrLMr_Rloc(l+1,nR)=-r(nR)*r(nR)/dLh * ( &
         &                   dTheta1S(lm)*dzAstrLM(lS) - &
         &                   dTheta1A(lm)*dzAstrLM(lA) )
         dzCorLMr_Rloc(l+1,nR) =-r(nR)*r(nR)/dLh * ( &
         &                    dTheta1S(lm)*dzCorLM(lS) - &
         &                    dTheta1A(lm)*dzCorLM(lA) )
         dzLFLMr_Rloc(l+1,nR)  = r(nR)*r(nR)/dLh * ( &
         &                     dTheta1S(lm)*dzLFLM(lS) - &
         &                     dTheta1A(lm)*dzLFLM(lA) )
         dzdVpLMr_Rloc(l+1,nR) =(zASL(l+1)-dzdVpLMr_Rloc(l+1,nR))/dtLast
         dzddVpLMr_Rloc(l+1,nR)=(zASL(l+1)/dtLast+dzddVpLMr_Rloc(l+1,nR))/dtLast
      end do

   end subroutine getTOfinish
!-----------------------------------------------------------------------------
end module torsional_oscillations
