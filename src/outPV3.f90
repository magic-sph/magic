module outPV3

   use precision_mod
   use parallel_mod, only: rank
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_m_max, n_phi_max, n_r_max, nrp, lm_max, &
       &                 l_max, minc, m_max, l_axi
   use radial_functions, only: cheb_norm, r_ICB, rscheme_oc, r_CMB, &
       &                       rscheme_oc, chebt_oc
   use physical_parameters, only: radratio
   use communications, only: gather_all_from_lo_to_rank0,gt_OC
   use blocking, only: lm2, lm2m, lm2l, lm2mc, st_map, lo_map
   use horizontal_data, only: dLh, dPhi
   use logic, only: lVerbose, l_SRIC
   use output_data, only: tag, sDens, nSmaxA, nZmaxA
   use LMLoop_data, only: llm, ulm
   use plms_theta, only: plm_theta
   use constants, only: pi, zero, one, two, half, ci
   use fft, only: fft_to_real
   use TO_helpers, only: getPAStr
   use cosine_transform_odd
 
   implicit none 
 
   private
 
   real(cp), allocatable :: rZ(:,:)
   real(cp), allocatable :: PlmS(:,:,:)
   real(cp), allocatable :: dPlmS(:,:,:)
   real(cp), allocatable :: PlmZ(:,:,:)
   real(cp), allocatable :: dPlmZ(:,:,:)
   real(cp), allocatable :: OsinTS(:,:)
   real(cp), allocatable :: VorOld(:,:,:)
 
   public :: initialize_outPV3, finalize_outPV3, outPV
  
contains

   subroutine initialize_outPV3

      allocate( rZ(nZmaxA/2+1,nSmaxA) )
      allocate( OsinTS(nZmaxA/2+1,nSmaxA) )
      bytes_allocated = bytes_allocated + 2*(nZmaxA/2+1)*nSmaxA*SIZEOF_DEF_REAL
      allocate( PlmS(l_max+1,nZmaxA/2+1,nSmaxA) )
      allocate( dPlmS(l_max+1,nZmaxA/2+1,nSmaxA) )
      bytes_allocated = bytes_allocated + &
                        2*(l_max+1)*(nZmaxA/2+1)*nSmaxA*SIZEOF_DEF_REAL
      allocate( PlmZ(lm_max,nZmaxA/2+1,nSmaxA) )
      allocate( dPlmZ(lm_max,nZmaxA/2+1,nSmaxA) )
      bytes_allocated = bytes_allocated + &
                        2*(lm_max)*(nZmaxA/2+1)*nSmaxA*SIZEOF_DEF_REAL
      allocate( VorOld(nrp,nZmaxA,nSmaxA) )
      bytes_allocated = bytes_allocated + nrp*nZmaxA*nSmaxA*SIZEOF_DEF_REAL

   end subroutine initialize_outPV3
!---------------------------------------------------------------------------------
   subroutine finalize_outPV3

      deallocate( rZ, OsinTS, PlmS, dPlmS, PlmZ, dPlmZ, VorOld )

   end subroutine finalize_outPV3
!---------------------------------------------------------------------------------
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
      real(cp) :: dzVpLMr(l_max+1,n_r_max)

      !--- Work array:
      real(cp) :: workAr(lm_max,n_r_max)

      integer :: lm, l, m

      real(cp) :: fac

      !--- define Grid
      integer :: nSmax,nS,nSI
      real(cp) ::  sZ(nSmaxA),dsZ ! cylindrical radius s and s-step

      integer :: nZ,nZmax,nZmaxNS
      integer, save :: nZC(nSmaxA),nZ2(nZmaxA,nSmaxA)
      integer, save :: nZS
      real(cp) :: zZ(nZmaxA),zstep!,zZC
      real(cp) :: VpAS(nZmaxA),omS(nZmaxA)

      !-- Plms: Plm,sin
      integer :: nR,nPhi,nC
      real(cp) :: thetaZ,rZS!,sinT,cosT

      !-- For PV output files: 
      character(len=80) :: fileName

      !-- Output of all three field components:
      real(cp) :: VsS(nrp,nZmaxA)
      real(cp) :: VpS(nrp,nZmaxA)
      real(cp) :: VzS(nrp,nZmaxA)
      real(cp) :: VorS(nrp,nZmaxA)
      real(cp) :: dpVorS(nrp,nZmaxA)
      real(outp) :: out1(n_phi_max*nZmaxA)
      real(outp) :: out2(n_phi_max*nZmaxA)
      real(outp) :: out3(n_phi_max*nZmaxA)
      real(outp) :: out4(n_phi_max*nZmaxA)
      real(outp) :: out5(n_phi_max*nZmaxA)
      real(cp), save :: timeOld

      complex(cp) :: wP(llm:ulm,n_r_max)
      complex(cp) :: dwP(llm:ulm,n_r_max)
      complex(cp) :: ddwP(llm:ulm,n_r_max)
      complex(cp) :: zP(llm:ulm,n_r_max)
      complex(cp) :: dzP(llm:ulm,n_r_max)

      !-- This may be deleted later:
      complex(cp), allocatable :: wP_global(:,:), dwP_global(:,:), ddwP_global(:,:)
      complex(cp), allocatable :: zP_global(:,:), dzP_global(:,:)

      integer :: n_pvz_file, n_vcy_file


      if ( lVerbose ) write(*,*) '! Starting outPV!'

      do nR=1,n_r_max
         do lm=llm,ulm
            l = lo_map%lm2l(lm)
            m = lo_map%lm2m(lm)
            wP(lm,nR) =w(lm,nR)*dLh(st_map%lm2(l,m))
            dwP(lm,nR)=dw(lm,nR)
            ddwP(lm,nR)=ddw(lm,nR)
            zP(lm,nR)=z(lm,nR)
            dzP(lm,nR)=dz(lm,nR)
         end do
      end do

      call rscheme_oc%costf1(wP,ulm-llm+1,1,ulm-llm+1)
      call rscheme_oc%costf1(dwP,ulm-llm+1,1,ulm-llm+1)
      call rscheme_oc%costf1(ddwP,ulm-llm+1,1,ulm-llm+1)
      call rscheme_oc%costf1(zP,ulm-llm+1,1,ulm-llm+1)
      call rscheme_oc%costf1(dzP,ulm-llm+1,1,ulm-llm+1)

      if ( rank == 0 ) then
         allocate( wP_global(1:lm_max,1:n_r_max) )
         allocate( dwP_global(1:lm_max,1:n_r_max) )
         allocate( ddwP_global(1:lm_max,1:n_r_max) )
         allocate( zP_global(1:lm_max,1:n_r_max) )
         allocate( dzP_global(1:lm_max,1:n_r_max) )
      else
         allocate( wP_global(1,1) )
         allocate( dwP_global(1,1) )
         allocate( ddwP_global(1,1) )
         allocate( zP_global(1,1) )
         allocate( dzP_global(1,1) )
      end if

      call gather_all_from_lo_to_rank0(gt_OC,wP,wP_global)
      call gather_all_from_lo_to_rank0(gt_OC,dwP,dwP_global)
      call gather_all_from_lo_to_rank0(gt_OC,ddwP,ddwP_global)
      call gather_all_from_lo_to_rank0(gt_OC,zP,zP_global)
      call gather_all_from_lo_to_rank0(gt_OC,dzP,dzP_global)

      if ( rank == 0 ) then

         nPVsets=nPVsets+1

         !-- Start with calculating advection due to axisymmetric flows:

         nSmax=n_r_max+int(r_ICB*real(n_r_max,kind=cp))
         nSmax=int(sDens*nSmax)
         if ( nSmax > nSmaxA ) then
            write(*,*) 'Increase nSmaxA in outPV!'
            write(*,*) 'Should be at least nSmax=',nSmax
            write(*,*) 'But is only=',nSmaxA
            stop
         end if
         nZmax=2*nSmax

         if ( l_stop_time ) then
            if ( l_SRIC  .and. omega_IC /= 0 ) then
               fac=one/omega_IC
            else
               fac=one
            end if
            do nR=1,n_r_max
               do l=1,l_max
                  lm=lm2(l,0)
                  dzVpLMr(l+1,nR)=fac*real(z(lm,nR))
               end do
            end do

            !---- Transform the contributions to cheb space:
            call chebt_oc%costf1(dzVpLMr,l_max+1,1,l_max+1,workAr)
         end if

         dsZ=r_CMB/real(nSmax,kind=cp)  ! Step in s controlled by nSmax
         nSI=0                  ! Inner core position
         do nS=1,nSmax
            sZ(nS)=(nS-half)*dsZ
            if ( sZ(nS) < r_ICB .and. nS > nSI ) nSI=nS
         end do
         zstep=2*r_CMB/real(nZmax-1,kind=cp)
         do nZ=1,nZmax
            zZ(nZ)=r_CMB-(nZ-1)*zstep
         end do

         !--- Open file for output:
         if ( l_stop_time ) then
            fileName='PVZ.'//tag
            open(newunit=n_pvz_file, file=fileName, form='unformatted', &
            &    status='unknown')
            write(n_pvz_file) real(time,kind=outp), real(nSmax,kind=outp), &
            &     real(nZmax,kind=outp), real(omega_IC,kind=outp),         &
            &     real(omega_ma,kind=outp)
            write(n_pvz_file) (real(sZ(nS),kind=outp),nS=1,nSmax)
            write(n_pvz_file) (real(zZ(nZ),kind=outp),nZ=1,nZmax)


            !--- Open file for the three flow components:
            fileName='Vcy.'//tag
            open(newunit=n_vcy_file, file=fileName,form='unformatted', &
            &    status='unknown')
            write(n_vcy_file) real(time,kind=outp), real(nSmax,kind=outp),&
            &     real(nZmax,kind=outp), real(n_phi_max,kind=outp),       &
            &     real(omega_IC,kind=outp), real(omega_ma,kind=outp),     &
            &     real(radratio,kind=outp), real(minc,kind=outp)
            write(n_vcy_file) (real(sZ(nS),kind=outp),nS=1,nSmax)
            write(n_vcy_file) (real(zZ(nZ),kind=outp),nZ=1,nZmax)
         end if



         do nS=1,nSmax

            !------ Get r,theta,Plm,dPlm for northern hemishere:
            if ( nPVsets == 1 ) then ! do this only for the first call !
               nZC(nS)=0 ! Points within shell
               do nZ=1,nZmax
                  rZS=sqrt(zZ(nZ)**2+sZ(nS)**2)
                  if ( rZS >= r_ICB .and. rZS <= r_CMB ) then
                     nZC(nS)=nZC(nS)+1  ! Counts all z within shell
                     nZ2(nZ,nS)=nZC(nS) ! No of point within shell
                     if ( zZ(nZ) > 0 ) then ! Onl north hemisphere
                        rZ(nZC(nS),nS)=rZS
                        thetaZ=atan2(sZ(nS),zZ(nZ))
                        OsinTS(nZC(nS),nS)=one/sin(thetaZ)
                        call plm_theta(thetaZ,l_max,0,minc,              &
                             &    PlmS(1,nZC(nS),nS),dPlmS(1,nZC(nS),nS),l_max+1,2)
                        call plm_theta(thetaZ,l_max,m_max,minc,          &
                             &        PlmZ(1,nZC(nS),nS),dPlmZ(1,nZC(nS),nS),lm_max,2)
                     end if
                  else
                     nZ2(nZ,nS)=-1 ! No z found within shell !
                  end if
               end do
            end if

            !-- Get azimuthal flow component in the shell
            nZmaxNS=nZC(nS) ! all z points within shell
            if ( l_stop_time ) then
               call getPAStr(VpAS,dzVpLMr,nZmaxNS,nZmaxA,l_max+1,      &
                    &        l_max,r_ICB,r_CMB,n_r_max,                &
                    &        rZ(1,nS),dPlmS(1,1,nS),OsinTS(1,nS))

               !-- Copy to array with all z-points
               do nZ=1,nZmax
                  rZS=sqrt(zZ(nZ)**2+sZ(nS)**2)
                  nZS=nZ2(nZ,nS)
                  if ( nZS > 0 ) then
                     omS(nZ)=VpAS(nZS)/sZ(nS)
                  else
                     if ( rZS <= r_ICB ) then
                        omS(nZ)=one
                     else
                        omS(nZ)=fac*omega_MA
                     end if
                  end if
               end do
            end if

            !-- Get all three components in the shell
            call getPVptr(wP,dwP,ddwP,zP,dzP,r_ICB,r_CMB,rZ(1,nS),                 &
                 &        nZmaxNS,nZmaxA,PlmZ(1,1,nS),dPlmZ(1,1,nS),OsinTS(1,nS),  &
                 &        VsS,VpS,VzS,VorS,dpVorS)

            if ( l_stop_time ) then
               write(n_pvz_file) (real(omS(nZ),kind=outp),nZ=1,nZmax)
               write(n_vcy_file) real(nZmaxNS,kind=outp)
               nC=0
               do nZ=1,nZmaxNS
                  do nPhi=1,n_phi_max
                     nC=nC+1
                     out1(nC)=real(VsS(nPhi,nZ),kind=outp) ! Vs
                     out2(nC)=real(VpS(nPhi,nZ),kind=outp) ! Vphi
                     out3(nC)=real(VzS(nPhi,nZ),kind=outp) ! Vz
                     out4(nC)=real(VorS(nPhi,nZ),kind=outp)
                     out5(nC)=(real(VorS(nPhi,nZ)-VorOld(nPhi,nZ,nS),kind=outp))/ &
                              (real(time-timeOld,kind=outp))
                  end do
               end do
               write(n_vcy_file) (out1(nZ),nZ=1,nC)
               write(n_vcy_file) (out2(nZ),nZ=1,nC)
               write(n_vcy_file) (out3(nZ),nZ=1,nC)
               write(n_vcy_file) (out4(nZ),nZ=1,nC)
               write(n_vcy_file) (out5(nZ),nZ=1,nC)
            else
               timeOld=time
               do nZ=1,nZmaxNS
                  do nPhi=1,n_phi_max
                     VorOld(nPhi,nZ,nS)=VorS(nPhi,nZ)
                  end do
               end do
            end if

         end do  ! Loop over s 

         if ( l_stop_time ) close(n_pvz_file)
         if ( l_stop_time ) close(n_vcy_file)

      end if ! Rank 0

      deallocate( wP_global, dwP_global, ddwP_global, zP_global, dzP_global )

   end subroutine outPV
!---------------------------------------------------------------------------------
   subroutine getPVptr(w,dw,ddw,z,dz,rMin,rMax,rS, &
                   nZmax,nZmaxA,PlmS,dPlmS,OsinTS, &
                          VrS,VpS,VtS,VorS,dpVorS)
      !
      !  This subroutine calculates the three flow conponents VrS,VtS,VpS at
      !  (r,theta,all phis) and (r,pi-theta, all phis). Here r=rS, PlmS=Plm(theta),
      !  dPlmS=sin(theta)*dTheta Plm(theta), and OsinTS=1/sin(theta).
      !  The flow is calculated for all n_phi_max azimuthal points used in the code,
      !  and for corresponding latitudes north and south of the equator.
      !  For lDeriv=.true. the subroutine also calculates dpEk and dzEk which
      !  are phi averages of (d Vr/d phi)**2 + (d Vtheta/ d phi)**2 + (d Vphi/ d phi)**2
      !  and (d Vr/d z)**2 + (d Vtheta/ d z)**2 + (d Vphi/ d z)**2, respectively.
      !  These two quantities are used ot calculate z and phi scale of the flow in
      !  s_getEgeos.f
      !  NOTE: on input w=l*(l+1)*w
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

      !--- Output: function on azimuthal grid points defined by FT!
      real(cp), intent(out) :: VrS(nrp,nZmaxA)
      real(cp), intent(out) :: VtS(nrp,nZmaxA)
      real(cp), intent(out) :: VpS(nrp,nZmaxA)
      real(cp), intent(out) :: VorS(nrp,nZmaxA)
      real(cp), intent(out) :: dpVorS(nrp,nZmaxA)

      !--- Local:
      real(cp) :: chebS(n_r_max)
      integer :: nS,nN,mc,lm,l,m,nCheb
      real(cp) :: x,phiNorm,mapFac,OS,cosT,sinT,Or_e1,Or_e2
      complex(cp) :: Vr,Vt,Vt1,Vt2,Vp1,Vp2,Vor,Vot1,Vot2
      real(cp) :: VotS(nrp,nZmaxA)
      complex(cp) :: wSr,dwSr,ddwSr,zSr,dzSr
      complex(cp) :: dp

      mapFac=two/(rMax-rMin)
      phiNorm=two*pi/n_phi_max

      do nS=1,nZmax
         do mc=1,nrp
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
            Vot2= (wSr*Or_e2-ddwSr)*PlmS(lm,nN)*dPhi(lm)
            VrS(2*mc-1, nN)=VrS(2*mc-1, nN)+ real(Vr)
            VrS(2*mc  , nN)=VrS(2*mc  , nN)+aimag(Vr)
            VtS(2*mc-1, nN)=VtS(2*mc-1, nN)+ real(Vt1+Vt2)
            VtS(2*mc  , nN)=VtS(2*mc  , nN)+aimag(Vt1+Vt2)
            VpS(2*mc-1, nN)=VpS(2*mc-1, nN)+ real(Vp1+Vp2)
            VpS(2*mc  , nN)=VpS(2*mc  , nN)+aimag(Vp1+Vp2)
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
            VrS(2*mc-1,nS)=sinT*Or_e2*VrS(2*mc-1,nS)+cosT*Or_e1*OS*VtS(2*mc-1,nS)
            VrS(2*mc  ,nS)=sinT*Or_e2*VrS(2*mc  ,nS)+cosT*Or_e1*OS*VtS(2*mc  ,nS)
            VpS(2*mc-1,nS)=Or_e1*OS*VpS(2*mc-1,nS)
            VpS(2*mc  ,nS)=Or_e1*OS*VpS(2*mc  ,nS)
            VtS(2*mc-1,nS)=cosT*Or_e2*VrS(2*mc-1,nS)-sinT*Or_e1*OS*VtS(2*mc-1,nS)
            VtS(2*mc  ,nS)=cosT*Or_e2*VrS(2*mc  ,nS)-sinT*Or_e1*OS*VtS(2*mc  ,nS)
            VorS(2*mc-1,nS)=cosT*Or_e2*VorS(2*mc-1,nS)-Or_e1*VotS(2*mc-1,nS)
            VorS(2*mc  ,nS)=cosT*Or_e2*VorS(2*mc  ,nS)-Or_e1*VotS(2*mc  ,nS)
         end do
         do mc=2*n_m_max+1,nrp
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
            Vr=cmplx(Or_e2*VrS(2*mc-1,nS), Or_e2*VrS(2*mc,nS),kind=cp)
            Vt=cmplx(Or_e1*OS*VtS(2*mc-1,nS), Or_e1*OS*VtS(2*mc,nS),kind=cp)
            VrS(2*mc-1,nS) =sinT* real(Vr)+cosT* real(Vt) ! this is now Vs
            VrS(2*mc  ,nS) =sinT*aimag(Vr)+cosT*aimag(Vt) ! this is now Vs
            VpS(2*mc-1,nS) =Or_e1*OS*VpS(2*mc-1,nS)
            VpS(2*mc  ,nS) =Or_e1*OS*VpS(2*mc  ,nS)
            VtS(2*mc-1,nS) =cosT* real(Vr)-sinT* real(Vt) ! this is now Vz
            VtS(2*mc  ,nS) =cosT*aimag(Vr)-sinT*aimag(Vt) ! this is now Vz
            VorS(2*mc-1,nS)=cosT*Or_e2*VorS(2*mc-1,nS)-Or_e1*VotS(2*mc-1,nS)
            VorS(2*mc  ,nS)=cosT*Or_e2*VorS(2*mc  ,nS)-Or_e1*VotS(2*mc  ,nS)
         end do
         do mc=2*n_m_max+1,nrp
            VrS(mc,nS) =0.0_cp
            VpS(mc,nS) =0.0_cp
            VtS(mc,nS) =0.0_cp
            VorS(mc,nS)=0.0_cp
         end do
      end do

      do nS=1,nZmax
         do mc=1,n_m_max
            dp=ci*real((mc-1)*minc,kind=cp)  ! - i m
            dpVorS(2*mc-1,nS)= real(dp)*VorS(2*mc-1,nS)-aimag(dp)*VorS(2*mc,nS)
            dpVorS(2*mc  ,nS)=aimag(dp)*VorS(2*mc-1,nS)+ real(dp)*VorS(2*mc,nS)
         end do
         do mc=2*n_m_max+1,nrp
            dpVorS(mc,nS)=0.0_cp
         end do
      end do

      !----- Transform m 2 phi for flow field:
      if ( .not. l_axi ) then
         call fft_to_real(VrS,nrp,nZmax)
         call fft_to_real(VtS,nrp,nZmax)
         call fft_to_real(VpS,nrp,nZmax)
         call fft_to_real(VorS,nrp,nZmax)
         call fft_to_real(dpVorS,nrp,nZmax)
      end if

   end subroutine getPVptr
!---------------------------------------------------------------------------------
end module outPV3
