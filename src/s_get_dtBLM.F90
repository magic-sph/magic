!$Id$
!***********************************************************************
    SUBROUTINE  get_dtBLM(nR,vr,vt,vp,br,bt,bp, &
                   n_theta_start,n_theta_block, &
                   BtVrLM,BpVrLM,BrVtLM,BrVpLM, &
                   BtVpLM,BpVtLM,BrVZLM,BtVZLM, &
                 BtVpCotLM,BpVtCotLM,BtVZcotLM, &
                 BtVpSn2LM,BpVtSn2LM,BtVZsn2LM)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!-----------------------------------------------------------------------

!  calculates non-linear products in grid-space for radial
!  level n_r and returns them in arrays wnlr1-3, snlr1-3, bnlr1-3

!  if lvelo >0 velocities are zero only the (vxB)
!  contributions to bnlr2-3 need to be calculated

!  vr...sr: (input) velocity, magnetic field comp. and derivs, entropy
!                   on grid points
!  i1: (input) range of points in theta for which calculation is done

!-----------------------------------------------------------------------

    USE truncation
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE blocking
    USE horizontal_data
    USE logic
#if (FFTLIB==JW)
  USE fft_JW
#elif (FFTLIB==MKL)
  USE fft_MKL
#endif

    IMPLICIT NONE

!-- input:
! include 'truncation.f'
! include 'c_horizontal.f'
! include 'c_num_param.f'
! include 'c_phys_param.f'
! include 'c_radial.f'
! include 'c_logic.f'
! include 'c_blocking.f'

    integer :: n_theta_start,n_theta_block,nR

!-- field components:
    real(kind=8) :: vr(nrp,nfs),vt(nrp,nfs),vp(nrp,nfs)
    real(kind=8) :: br(nrp,nfs),bt(nrp,nfs),bp(nrp,nfs)

!-- output:
    complex(kind=8) :: BtVrLM(*),BpVrLM(*)
    complex(kind=8) :: BrVtLM(*),BrVpLM(*)
    complex(kind=8) :: BtVpLM(*),BpVtLM(*)
    complex(kind=8) :: BrVZLM(*),BtVZLM(*)
    complex(kind=8) :: BpVtCotLM(*),BtVpCotLM(*),BtVZcotLM(*)
    complex(kind=8) :: BtVpSn2LM(*),BpVtSn2LM(*)
    complex(kind=8) :: BtVZsn2LM(*)

!-- local:
    integer :: n_theta,n_phi,n_theta_nhs
    real(kind=8) :: fac,facCot
    real(kind=8) :: BtVr(nrp,nfs),BpVr(nrp,nfs)
    real(kind=8) :: BrVt(nrp,nfs),BrVp(nrp,nfs)
    real(kind=8) :: BtVp(nrp,nfs),BpVt(nrp,nfs)
    real(kind=8) :: BrVZ(nrp,nfs),BtVZ(nrp,nfs)
    real(kind=8) :: BpVtCot(nrp,nfs),BtVpCot(nrp,nfs)
    real(kind=8) :: BpVtSn2(nrp,nfs),BtVpSn2(nrp,nfs)
    real(kind=8) :: BtVZcot(nrp,nfs),BtVZsn2(nrp,nfs)
    real(kind=8) :: vpAS

    INTEGER :: lm,l

!-- end of declaration
!---------------------------------------------------------------------------


    DO n_theta=1,n_theta_block ! loop over ic-points, alternating north/south

        n_theta_nhs=(n_theta_start+n_theta)/2
        fac=osn2(n_theta_nhs)
        facCot=cosn2(n_theta_nhs)*osn1(n_theta_nhs)
        IF ( MOD(n_theta,2) == 0 ) facCot=-facCot  ! SHS

        DO n_phi=1,n_phi_max
            BtVr(n_phi,n_theta)= &
                fac*orho1(nR)*bt(n_phi,n_theta)*vr(n_phi,n_theta)
            BpVr(n_phi,n_theta)= &
                fac*orho1(nR)*bp(n_phi,n_theta)*vr(n_phi,n_theta)
        END DO
        BtVr(n_phi_max+1,n_theta)=0.d0
        BtVr(n_phi_max+2,n_theta)=0.d0
        BpVr(n_phi_max+1,n_theta)=0.d0
        BpVr(n_phi_max+2,n_theta)=0.d0

        do n_phi=1,n_phi_max
            BrVt(n_phi,n_theta)= &
                fac*orho1(nR)*vt(n_phi,n_theta)*br(n_phi,n_theta)
            BrVp(n_phi,n_theta)= &
                fac*orho1(nR)*vp(n_phi,n_theta)*br(n_phi,n_theta)
        end do
        BrVt(n_phi_max+1,n_theta)=0.d0
        BrVt(n_phi_max+2,n_theta)=0.d0
        BrVp(n_phi_max+1,n_theta)=0.d0
        BrVp(n_phi_max+2,n_theta)=0.d0

        vpAS=0.d0
        do n_phi=1,n_phi_max
            BtVp(n_phi,n_theta)= &
                fac*orho1(nR)*bt(n_phi,n_theta)*vp(n_phi,n_theta)
            BpVt(n_phi,n_theta)= &
                fac*orho1(nR)*bp(n_phi,n_theta)*vt(n_phi,n_theta)
            vpAS=vpAS+orho1(nR)*vp(n_phi,n_theta)
        end do
        BtVp(n_phi_max+1,n_theta)=0.d0
        BtVp(n_phi_max+2,n_theta)=0.d0
        BpVt(n_phi_max+1,n_theta)=0.d0
        BpVt(n_phi_max+2,n_theta)=0.d0
        vpAS=vpAS/dble(n_phi_max)

    !---- For toroidal terms that cancel:
        do n_phi=1,n_phi_max
            BpVtCot(n_phi,n_theta)= &
                facCot*orho1(nR)*bp(n_phi,n_theta)*vt(n_phi,n_theta)
            BtVpCot(n_phi,n_theta)= &
                facCot*orho1(nR)*bt(n_phi,n_theta)*vp(n_phi,n_theta)
            BpVtSn2(n_phi,n_theta)= &
                fac*fac*orho1(nR)*bp(n_phi,n_theta)*vt(n_phi,n_theta)
            BtVpSn2(n_phi,n_theta)= &
                fac*fac*orho1(nR)*bt(n_phi,n_theta)*vp(n_phi,n_theta)
        end do
        BpVtCot(n_phi_max+1,n_theta)=0.d0
        BpVtCot(n_phi_max+2,n_theta)=0.d0
        BtVpCot(n_phi_max+1,n_theta)=0.d0
        BtVpCot(n_phi_max+2,n_theta)=0.d0
        BpVtSn2(n_phi_max+1,n_theta)=0.d0
        BpVtSn2(n_phi_max+2,n_theta)=0.d0
        BtVpSn2(n_phi_max+1,n_theta)=0.d0
        BtVpSn2(n_phi_max+2,n_theta)=0.d0

    !---- For omega effect:
        do n_phi=1,n_phi_max
            BrVZ(n_phi,n_theta)= &
                fac*br(n_phi,n_theta)*vpAS
            BtVZ(n_phi,n_theta)= &
                fac*bt(n_phi,n_theta)*vpAS
            BtVZcot(n_phi,n_theta)= &
                facCot*bt(n_phi,n_theta)*vpAS
            BtVZsn2(n_phi,n_theta)= &
                fac*fac*bt(n_phi,n_theta)*vpAS
        end do
        BrVZ(n_phi_max+1,n_theta)=0.d0
        BrVZ(n_phi_max+2,n_theta)=0.d0
        BtVZ(n_phi_max+1,n_theta)=0.d0
        BtVZ(n_phi_max+2,n_theta)=0.d0
        BtVZCot(n_phi_max+1,n_theta)=0.d0
        BtVZCot(n_phi_max+2,n_theta)=0.d0
        BtVZsn2(n_phi_max+1,n_theta)=0.d0
        BtVZsn2(n_phi_max+2,n_theta)=0.d0

    end do

!-- Fourier transform phi2m
    CALL fft_thetab(BtVpSn2,-1)
    call fft_thetab(BpVtSn2,-1)
    call fft_thetab(BtVpCot,-1)
    call fft_thetab(BpVtCot,-1)
    call fft_thetab(BtVr,-1)
    call fft_thetab(BpVr,-1)
    call fft_thetab(BrVt,-1)
    call fft_thetab(BrVp,-1)
    call fft_thetab(BtVp,-1)
    call fft_thetab(BpVt,-1)
    call fft_thetab(BrVZ,-1)
    call fft_thetab(BtVZ,-1)
    call fft_thetab(BtVZcot,-1)
    call fft_thetab(BtVZsn2,-1)

     
!-- Legendre transform: theta2l
    call legTF3(n_theta_start, &
                BtVrLM,BpVrLM,BrVtLM,BtVr,BpVr,BrVt)
    call legTF3(n_theta_start, &
                BrVpLM,BtVpLM,BpVtLM,BrVp,BtVp,BpVt)
    call legTF3(n_theta_start,                 &
                BtVpCotLM,BpVtCotLM,BtVZcotLM, &
                BtVpCot,BpVtCot,BtVZcot)
    call legTF3(n_theta_start,           &
                BrVZLM,BtVZLM,BtVZsn2LM, &
                BrVZ,BtVZ,BtVZsn2)
    call legTF2(n_theta_start, &
                BtVpSn2LM,BpVtSn2LM,BtVpSn2,BpVtSn2)

    DO l=0,l_max
        lm=l2lmAS(l)
        BtVrLM(lm)   =CMPLX(REAL(BtVrLM(lm)),   0.D0,KIND=KIND(0.d0))
        BpVrLM(lm)   =CMPLX(REAL(BpVrLM(lm)),   0.D0,KIND=KIND(0.d0))
        BrVtLM(lm)   =CMPLX(REAL(BrVtLM(lm)),   0.D0,KIND=KIND(0.d0))
        BrVpLM(lm)   =CMPLX(REAL(BrVpLM(lm)),   0.D0,KIND=KIND(0.d0))
        BtVpLM(lm)   =CMPLX(REAL(BtVpLM(lm)),   0.D0,KIND=KIND(0.d0))
        BpVtLM(lm)   =CMPLX(REAL(BpVtLM(lm)),   0.D0,KIND=KIND(0.d0))
        BtVpCotLM(lm)=CMPLX(REAL(BtVpCotLM(lm)),0.D0,KIND=KIND(0.d0))
        BpVtCotLM(lm)=CMPLX(REAL(BpVtCotLM(lm)),0.D0,KIND=KIND(0.d0))
        BtVZCotLM(lm)=CMPLX(REAL(BtVZCotLM(lm)),0.D0,KIND=KIND(0.d0))
        BrVZLM(lm)   =CMPLX(REAL(BrVZLM(lm))   ,0.D0,KIND=KIND(0.d0))
        BtVZLM(lm)   =CMPLX(REAL(BtVZLM(lm))   ,0.D0,KIND=KIND(0.d0))
        BtVZSn2LM(lm)=CMPLX(REAL(BtVZSn2LM(lm)),0.D0,KIND=KIND(0.d0))
        BtVpSn2LM(lm)=CMPLX(REAL(BtVpSn2LM(lm)),0.D0,KIND=KIND(0.d0))
        BpVtSn2LM(lm)=CMPLX(REAL(BpVtSn2LM(lm)),0.D0,KIND=KIND(0.d0))
    END DO


!-- The different parts of the magnetic fieldline stretching
!   and fieldline advection are now in (l,m) space

    RETURN
    END SUBROUTINE get_dtBLM
!------------------------------------------------------------------------------------
