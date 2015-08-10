!$Id$
#include "perflib_preproc.cpp"
module legendre_grid_to_spec

   use truncation, only: n_m_max, nrp, lmP_max, n_theta_max
   use blocking, only: nfs, sizeThetaB
   use horizontal_data, only: lStartP, wPlm, lmOddP, lStopP

   implicit none

   private

   public :: legTF1, legTF2, legTF3, legTFAS, legTFAS2

contains

   subroutine legTF1(nThetaStart,f1LM,f1TM)
      !  +-------------+----------------+------------------------------------+
      !  |                                                                   |
      !  |  Legendre transform (n_r,n_theta,m) to (n_r,l,m)                  |
      !  |  [grid to spectral] for 2 arrays                                  |
      !  |  f1TM (input) to f1LM (output)                                    |
      !  |  One call to this routine does part of the transform              |
      !  |  by summation over theta points in one theta block:               |
      !  |      nThetaStart,..,nThetaStart+n_theta_block-1                   |
      !  |                                                                   |
      !--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+

      !-- Input variables:
      integer,      intent(in) :: nThetaStart ! First no of theta on block
      real(kind=8), intent(in) :: f1TM(nrp,nfs)

      !-- Output variable:
      complex(kind=8), intent(inout) :: f1LM(lmP_max)

      !-- Local variables:
      integer :: nThetaN     ! No. of theta in NHS
      integer :: nThetaS     ! No. of theta in SHS
      integer :: nThetaNHS   ! No. of thetas in one HS only
      integer :: nThetaB1    ! No. of theta in block
      integer :: nThetaB2    ! No. of theta in block
      integer :: nTheta1     ! No. of theta (in one HS)
      integer :: nTheta2     ! No. of theta (in one HS)
      integer :: nThetaMin   ! Where to start in block

      integer :: mc          ! counter of spherical order
      integer :: lmS,lm      ! counter of spherical mode

      complex(kind=8) :: f1ES(n_m_max,nfs/2),f1ES1,f1ES2
      complex(kind=8) :: f1EA(n_m_max,nfs/2),f1EA1,f1EA2

      !PERFON('legTF1')
      !-- Unscrambles equatorially symmetric and antisymmetric contributions:
      nThetaNHS=0
      do nThetaN=1,sizeThetaB,2 ! thetas in NHS
         nThetaS=nThetaN+1      ! thetas in SHS
         nThetaNHS=nThetaNHS+1  ! thetas counted in NHS only
         do mc=1,n_m_max        ! counts spherical harmonic orders
            f1ES(mc,nThetaNHS)=cmplx(f1TM(2*mc-1,nThetaN),f1TM(2*mc,nThetaN),&
                                     kind=kind(0.d0)) +                      &
                               cmplx(f1TM(2*mc-1,nThetaS),f1TM(2*mc,nThetaS),&
                                     kind=kind(0.d0))                        
            f1EA(mc,nThetaNHS)=cmplx(f1TM(2*mc-1,nThetaN),f1TM(2*mc,nThetaN),&
                                     kind=kind(0.d0)) -                      &
                               cmplx(f1TM(2*mc-1,nThetaS),f1TM(2*mc,nThetaS),&
                                     kind=kind(0.d0))                       
         end do
      end do

      !- Start with first two thetas for first theta block:
      if ( nThetaStart == 1 ) then

         do mc=1,n_m_max
            lmS=lStopP(mc)
            f1ES1=f1ES(mc,1)
            f1ES2=f1ES(mc,2)
            f1EA1=f1EA(mc,1)
            f1EA2=f1EA(mc,2)
            do lm=lStartP(mc),lmS-1,2
               f1LM(lm)  =f1ES1*wPlm(lm,1) +f1ES2*wPlm(lm,2)
               f1LM(lm+1)=f1EA1*wPlm(lm+1,1)+f1EA2*wPlm(lm+1,2)
            end do
            if ( lmOddP(mc) ) then
               f1LM(lmS) =f1ES1*wPlm(lmS,1)+f1ES2*wPlm(lmS,2)
            end if
         end do
           
         if ( sizeThetaB <= 4 ) return ! return

      end if

       
      !-- Loop over half of the thetas with step 2 unrolling:
      nThetaMin=1
      if ( nThetaStart == 1 ) nThetaMin=3
      nTheta1=(nThetaStart-1)/2+nThetaMin-2 ! NHS thetas covered before
      do nThetaB1=nThetaMin,sizeThetaB/2,2
         nThetaB2=nThetaB1+1
         nTheta1 =nTheta1+2
         nTheta2 =nTheta1+1

         do mc=1,n_m_max
            lmS=lStopP(mc)
            f1ES1=f1ES(mc,nThetaB1)
            f1ES2=f1ES(mc,nThetaB2)
            f1EA1=f1EA(mc,nThetaB1)
            f1EA2=f1EA(mc,nThetaB2)
            do lm=lStartP(mc),lmS-1,2
               f1LM(lm)  =f1LM(lm)   + f1ES1*wPlm(lm,nTheta1) + f1ES2*wPlm(lm,nTheta2)
               f1LM(lm+1)=f1LM(lm+1) + f1EA1*wPlm(lm+1,nTheta1) + f1EA2*wPlm(lm+1,nTheta2)
            end do
            if ( lmOddP(mc) ) then
               f1LM(lmS)=f1LM(lmS) + f1ES1*wPlm(lmS,nTheta1) + f1ES2*wPlm(lmS,nTheta2)
            end if
        end do

      end do  !  loop over theta in block

      !PERFOFF

   end subroutine legTF1
!------------------------------------------------------------------------
   subroutine legTF2(nThetaStart,f1LM,f2LM,f1TM,f2TM)
      !  +-------------+----------------+------------------------------------+
      !  |                                                                   |
      !  |  Legendre transform (n_r,n_theta,m) to (n_r,l,m)                  |
      !  |  [grid to spectral] for 2 arrays                                  |
      !  |  f1TM,f2TM (input) to f1LM,f2LM (output)                          |
      !  |  One call to this routine does part of the transform              |
      !  |  by summation over theta points in on theta block:                |
      !  |      nThetaStart,..,nThetaStart+n_theta_block-1                   |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+

      !-- Input variables:
      integer,      intent(in) :: nThetaStart ! First no of theta on block
      real(kind=8), intent(in) :: f1TM(nrp,nfs),f2TM(nrp,nfs)

      !-- Output variables:
      complex(kind=8), intent(inout) :: f1LM(lmP_max),f2LM(lmP_max)

      !-- Local variables:
      integer :: nThetaN     ! No. of theta in NHS
      integer :: nThetaS     ! No. of theta in SHS
      integer :: nThetaNHS   ! No. of thetas in one HS only
      integer :: nThetaB1    ! No. of theta in block
      integer :: nThetaB2    ! No. of theta in block
      integer :: nTheta1     ! No. of theta (in one HS)
      integer :: nTheta2     ! No. of theta (in one HS)
      integer :: nThetaMin   ! Where to start in block

      integer :: mc          ! counter of spherical order
      integer :: lmS,lm      ! counter of spherical mode

      complex(kind=8) :: f1ES(n_m_max,nfs/2),f1ES1,f1ES2
      complex(kind=8) :: f1EA(n_m_max,nfs/2),f1EA1,f1EA2
      complex(kind=8) :: f2ES(n_m_max,nfs/2),f2ES1,f2ES2
      complex(kind=8) :: f2EA(n_m_max,nfs/2),f2EA1,f2EA2

      !PERFON('legTF2')
      !-- Unscrambles equatorially symmetric and antisymmetric contributions:
      nThetaNHS=0
      do nThetaN=1,sizeThetaB,2 ! thetas in NHS
         nThetaS=nThetaN+1      ! thetas in SHS
         nThetaNHS=nThetaNHS+1  ! thetas counted in NHS only
         do mc=1,n_m_max        ! counts spherical harmonic orders
            f1ES(mc,nThetaNHS)=cmplx(f1TM(2*mc-1,nThetaN),f1TM(2*mc,nThetaN),&
                                     kind=kind(0.d0)) +                      &
                               cmplx(f1TM(2*mc-1,nThetaS),f1TM(2*mc,nThetaS),&
                                     kind=kind(0.d0))
            f1EA(mc,nThetaNHS)=cmplx(f1TM(2*mc-1,nThetaN),f1TM(2*mc,nThetaN),&
                                     kind=kind(0.d0)) -                      &
                               cmplx(f1TM(2*mc-1,nThetaS),f1TM(2*mc,nThetaS),&
                                     kind=kind(0.d0))
            f2ES(mc,nThetaNHS)=cmplx(f2TM(2*mc-1,nThetaN),f2TM(2*mc,nThetaN),&
                                     kind=kind(0.d0)) +                      &
                               cmplx(f2TM(2*mc-1,nThetaS),f2TM(2*mc,nThetaS),&
                                     kind=kind(0.d0))
            f2EA(mc,nThetaNHS)=cmplx(f2TM(2*mc-1,nThetaN),f2TM(2*mc,nThetaN),&
                                     kind=kind(0.d0)) -                      &
                               cmplx(f2TM(2*mc-1,nThetaS),f2TM(2*mc,nThetaS),&
                                     kind=kind(0.d0))
         end do
      end do

      !- Start with first two thetas for first theta block:
      if ( nThetaStart == 1 ) then

         do mc=1,n_m_max
            lmS=lStopP(mc)
            f1ES1=f1ES(mc,1)
            f1ES2=f1ES(mc,2)
            f2ES1=f2ES(mc,1)
            f2ES2=f2ES(mc,2)
            do lm=lStartP(mc),lmS-1,2
               f1LM(lm)  =f1ES1*wPlm(lm,1) +f1ES2*wPlm(lm,2)
               f2LM(lm)  =f2ES1*wPlm(lm,1) +f2ES2*wPlm(lm,2)
            end do
            if ( lmOddP(mc) ) then
               f1LM(lmS) =f1ES1*wPlm(lmS,1)+f1ES2*wPlm(lmS,2)
               f2LM(lmS) =f2ES1*wPlm(lmS,1)+f2ES2*wPlm(lmS,2)
            end if
         end do
         do mc=1,n_m_max
            lmS=lStopP(mc)
            f1EA1=f1EA(mc,1)
            f1EA2=f1EA(mc,2)
            f2EA1=f2EA(mc,1)
            f2EA2=f2EA(mc,2)
            do lm=lStartP(mc),lmS-1,2
               f1LM(lm+1)=f1EA1*wPlm(lm+1,1)+f1EA2*wPlm(lm+1,2)
               f2LM(lm+1)=f2EA1*wPlm(lm+1,1)+f2EA2*wPlm(lm+1,2)
            end do
         end do
           
         if ( sizeThetaB <= 4 ) return !RETURN

       end if

      !-- Loop over half of the thetas with step 2 unrolling:
      nThetaMin=1
      if ( nThetaStart == 1 ) nThetaMin=3
      nTheta1=(nThetaStart-1)/2+nThetaMin-2 ! NHS thetas covered before
      do nThetaB1=nThetaMin,sizeThetaB/2,2
         nThetaB2=nThetaB1+1
         nTheta1 =nTheta1+2
         nTheta2 =nTheta1+1

         do mc=1,n_m_max
            lmS=lStopP(mc)
            f1ES1=f1ES(mc,nThetaB1)
            f1ES2=f1ES(mc,nThetaB2)
            f1EA1=f1EA(mc,nThetaB1)
            f1EA2=f1EA(mc,nThetaB2)
            do lm=lStartP(mc),lmS-1,2
               f1LM(lm)  =f1LM(lm)   + f1ES1*wPlm(lm,nTheta1) + f1ES2*wPlm(lm,nTheta2)
               f1LM(lm+1)=f1LM(lm+1) + f1EA1*wPlm(lm+1,nTheta1) + f1EA2*wPlm(lm+1,nTheta2)
            end do
            if ( lmOddP(mc) ) then
               f1LM(lmS)=f1LM(lmS) + f1ES1*wPlm(lmS,nTheta1) + f1ES2*wPlm(lmS,nTheta2)
            end if
         end do
         do mc=1,n_m_max
            lmS=lStopP(mc)
            f2ES1=f2ES(mc,nThetaB1)
            f2ES2=f2ES(mc,nThetaB2)
            f2EA1=f2EA(mc,nThetaB1)
            f2EA2=f2EA(mc,nThetaB2)
            do lm=lStartP(mc),lmS-1,2
               f2LM(lm)  =f2LM(lm)   + f2ES1*wPlm(lm,nTheta1) + f2ES2*wPlm(lm,nTheta2)
               f2LM(lm+1)=f2LM(lm+1) + f2EA1*wPlm(lm+1,nTheta1) + f2EA2*wPlm(lm+1,nTheta2)
            end do
            if ( lmOddP(mc) ) then
               f2LM(lmS)=f2LM(lmS) + f2ES1*wPlm(lmS,nTheta1) + f2ES2*wPlm(lmS,nTheta2)
            end if
         end do

      end do  !  loop over theta in block

      !PERFOFF

   end subroutine legTF2
!------------------------------------------------------------------------
   subroutine legTF3(nThetaStart,f1LM,f2LM,f3LM,f1TM,f2TM,f3TM)
      !  +-------------+----------------+------------------------------------+
      !  |                                                                   |
      !  |  Legendre transform (n_r,n_theta,m) to (n_r,l,m)                  |
      !  |  [grid to spectral] for 2 arrays                                  |
      !  |  ancl1/2/3 (input) to flm1/2/3 (output)                               |
      !  |  One call to this routine does part of the transform              |
      !  |  by summation over theta points in on theta block:                |
      !  |      nThetaStart,..,nThetaStart+n_theta_block-1                   |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+

      !-- Input variables:
      integer,      intent(in) :: nThetaStart ! First no of theta on block
      real(kind=8), intent(in) :: f1TM(nrp,nfs),f2TM(nrp,nfs),f3TM(nrp,nfs)

      !-- Output variables:
      complex(kind=8), intent(inout) :: f1LM(lmP_max),f2LM(lmP_max),f3LM(lmP_max)

      !-- Local variables:
      integer :: nThetaN     ! No. of theta in NHS
      integer :: nThetaS     ! No. of theta in SHS
      integer :: nThetaNHS   ! No. of thetas in one HS only
      integer :: nThetaB1    ! No. of theta in block
      integer :: nThetaB2    ! No. of theta in block
      integer :: nTheta1     ! No. of theta (in one HS)
      integer :: nTheta2     ! No. of theta (in one HS)
      integer :: nThetaMin   ! Where to start in block

      integer :: mc          ! counter of spherical order
      integer :: lmS,lm      ! counter of spherical mode

      complex(kind=8) :: f1ES(n_m_max,nfs/2),f1ES1,f1ES2
      complex(kind=8) :: f1EA(n_m_max,nfs/2),f1EA1,f1EA2
      complex(kind=8) :: f2ES(n_m_max,nfs/2),f2ES1,f2ES2
      complex(kind=8) :: f2EA(n_m_max,nfs/2),f2EA1,f2EA2
      complex(kind=8) :: f3ES(n_m_max,nfs/2),f3ES1,f3ES2
      complex(kind=8) :: f3EA(n_m_max,nfs/2),f3EA1,f3EA2

      !PERFON('legTF3')
      !-- Unscrambles equatorially symmetric and antisymmetric contributions:
      nThetaNHS=0
      do nThetaN=1,sizeThetaB,2 ! thetas in NHS
         nThetaS=nThetaN+1      ! thetas in SHS
         nThetaNHS=nThetaNHS+1  ! thetas counted in NHS only
         do mc=1,n_m_max        ! counts spherical harmonic orders
            f1ES(mc,nThetaNHS)=cmplx(f1TM(2*mc-1,nThetaN),f1TM(2*mc,nThetaN),&
                                     kind=kind(0.d0)) +                      &
                               cmplx(f1TM(2*mc-1,nThetaS),f1TM(2*mc,nThetaS),&
                                     kind=kind(0.d0))
            f1EA(mc,nThetaNHS)=cmplx(f1TM(2*mc-1,nThetaN),f1TM(2*mc,nThetaN),&
                                     kind=kind(0.d0)) -                      &
                               cmplx(f1TM(2*mc-1,nThetaS),f1TM(2*mc,nThetaS),&
                                     kind=kind(0.d0))
            f2ES(mc,nThetaNHS)=cmplx(f2TM(2*mc-1,nThetaN),f2TM(2*mc,nThetaN),&
                                     kind=kind(0.d0)) +                      &
                               cmplx(f2TM(2*mc-1,nThetaS),f2TM(2*mc,nThetaS),&
                                     kind=kind(0.d0))
            f2EA(mc,nThetaNHS)=cmplx(f2TM(2*mc-1,nThetaN),f2TM(2*mc,nThetaN),&
                                     kind=kind(0.d0)) -                      &
                               cmplx(f2TM(2*mc-1,nThetaS),f2TM(2*mc,nThetaS),&
                                     kind=kind(0.d0))
            f3ES(mc,nThetaNHS)=cmplx(f3TM(2*mc-1,nThetaN),f3TM(2*mc,nThetaN),&
                                     kind=kind(0.d0)) +                      &
                               cmplx(f3TM(2*mc-1,nThetaS),f3TM(2*mc,nThetaS),&
                                     kind=kind(0.d0))
            f3EA(mc,nThetaNHS)=cmplx(f3TM(2*mc-1,nThetaN),f3TM(2*mc,nThetaN),&
                                     kind=kind(0.d0)) -                      &
                               cmplx(f3TM(2*mc-1,nThetaS),f3TM(2*mc,nThetaS),&
                                     kind=kind(0.d0))
         end do
      end do

      !- Start with first two thetas for first theta block:

      if ( nThetaStart == 1 ) then

         do mc=1,n_m_max
            lmS=lStopP(mc)
            f1ES1=f1ES(mc,1)
            f1ES2=f1ES(mc,2)
            f2ES1=f2ES(mc,1)
            f2ES2=f2ES(mc,2)
            f3ES1=f3ES(mc,1)
            f3ES2=f3ES(mc,2)
            do lm=lStartP(mc),lmS-1,2
               f1LM(lm)  =f1ES1*wPlm(lm,1) +f1ES2*wPlm(lm,2)
               f2LM(lm)  =f2ES1*wPlm(lm,1) +f2ES2*wPlm(lm,2)
               f3LM(lm)  =f3ES1*wPlm(lm,1) +f3ES2*wPlm(lm,2)
            end do
            if ( lmOddP(mc) ) then
               f1LM(lmS) =f1ES1*wPlm(lmS,1)+f1ES2*wPlm(lmS,2)
               f2LM(lmS) =f2ES1*wPlm(lmS,1)+f2ES2*wPlm(lmS,2)
               f3LM(lmS) =f3ES1*wPlm(lmS,1)+f3ES2*wPlm(lmS,2)
            end if
         end do
         do mc=1,n_m_max
            lmS=lStopP(mc)
            f1EA1=f1EA(mc,1)
            f1EA2=f1EA(mc,2)
            f2EA1=f2EA(mc,1)
            f2EA2=f2EA(mc,2)
            f3EA1=f3EA(mc,1)
            f3EA2=f3EA(mc,2)
            do lm=lStartP(mc),lmS-1,2
               f1LM(lm+1)=f1EA1*wPlm(lm+1,1)+f1EA2*wPlm(lm+1,2)
               f2LM(lm+1)=f2EA1*wPlm(lm+1,1)+f2EA2*wPlm(lm+1,2)
               f3LM(lm+1)=f3EA1*wPlm(lm+1,1)+f3EA2*wPlm(lm+1,2)
            end do
         end do

         if ( sizeThetaB <= 4 ) return !RETURN

      end if


      !-- Loop over half of the thetas with step 2 unrolling:
      nThetaMin=1
      if ( nThetaStart == 1 ) nThetaMin=3
      nTheta1=(nThetaStart-1)/2+nThetaMin-2 ! NHS thetas covered before
      do nThetaB1=nThetaMin,sizeThetaB/2,2
         nThetaB2=nThetaB1+1
         nTheta1 =nTheta1+2
         nTheta2 =nTheta1+1

         do mc=1,n_m_max
            lmS=lStopP(mc)
            f1ES1=f1ES(mc,nThetaB1)
            f1ES2=f1ES(mc,nThetaB2)
            f1EA1=f1EA(mc,nThetaB1)
            f1EA2=f1EA(mc,nThetaB2)
            f2ES1=f2ES(mc,nThetaB1)
            f2ES2=f2ES(mc,nThetaB2)
            do lm=lStartP(mc),lmS-1,2
               f1LM(lm)  =f1LM(lm)   + f1ES1*wPlm(lm,nTheta1) + &
                    f1ES2*wPlm(lm,nTheta2)
               f1LM(lm+1)=f1LM(lm+1) + f1EA1*wPlm(lm+1,nTheta1) + &
                    f1EA2*wPlm(lm+1,nTheta2)
               f2LM(lm)  =f2LM(lm)   + f2ES1*wPlm(lm,nTheta1) + &
                    f2ES2*wPlm(lm,nTheta2)
            end do
            if ( lmOddP(mc) ) then
               f1LM(lmS)=f1LM(lmS) + f1ES1*wPlm(lmS,nTheta1) + &
                    f1ES2*wPlm(lmS,nTheta2)
               f2LM(lmS)=f2LM(lmS) + f2ES1*wPlm(lmS,nTheta1) + &
                    f2ES2*wPlm(lmS,nTheta2)
            end if
         end do
         do mc=1,n_m_max
            lmS=lStopP(mc)
            f3ES1=f3ES(mc,nThetaB1)
            f3ES2=f3ES(mc,nThetaB2)
            f3EA1=f3EA(mc,nThetaB1)
            f3EA2=f3EA(mc,nThetaB2)
            f2EA1=f2EA(mc,nThetaB1)
            f2EA2=f2EA(mc,nThetaB2)
            do lm=lStartP(mc),lmS-1,2
               f3LM(lm)  =f3LM(lm)   + f3ES1*wPlm(lm,nTheta1) + &
                    f3ES2*wPlm(lm,nTheta2)
               f3LM(lm+1)=f3LM(lm+1) + f3EA1*wPlm(lm+1,nTheta1) + &
                    f3EA2*wPlm(lm+1,nTheta2)
               f2LM(lm+1)=f2LM(lm+1) + f2EA1*wPlm(lm+1,nTheta1) + &
                    f2EA2*wPlm(lm+1,nTheta2)
            end do
            if ( lmOddP(mc) ) &
                 f3LM(lmS)=f3LM(lmS) + f3ES1*wPlm(lmS,nTheta1) + &
                 f3ES2*wPlm(lmS,nTheta2)
         end do


      end do  !  loop over theta in block

      !PERFOFF
   end subroutine legTF3
!------------------------------------------------------------------------
   subroutine legTFAS(flm1,ft1,lmMax,nThetaStart,sizeThetaB)
      !  +-------------+----------------+------------------------------------+
      !  |                                                                   |
      !  |  Legendre transform (n_r,n_theta,m) to (n_r,l,m)                  |
      !  |  [grid to spectral] for 2 arrays                                  |
      !  |  ancl1/2 (input) to flm1/2 (output)                               |
      !  |  One call to this routine does part of the transform              |
      !  |  by summation over theta points in on theta block:                |
      !  |      n_theta_min,..,n_theta_min+n_theta_block-1                   |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+

      implicit none

      !-- Input variables:
      integer,      intent(in) :: lmMax          ! Number of modes to be processed
      integer,      intent(in) :: nThetaStart    ! First no of theta on block
      integer,      intent(in) :: sizeThetaB     ! Size of theta block
      real(kind=8), intent(in) :: ft1(*)

      !-- Output: transformed arrays anlc1,anlc2
      real(kind=8), intent(out) :: flm1(*)

      !-- Local variables:
      integer :: nThetaN,nThetaS
      integer :: nThetaHS
      integer :: nTheta1,nTheta2
      integer :: nThetaC1,nThetaC2
      integer :: nThetaMin
      integer :: lm1,lm2

      real(kind=8) :: f1p(n_theta_max/2),f1m(n_theta_max/2)

      !-- Prepare arrays of sums and differences:
      nThetaHS=0
      do nThetaN=1,sizeThetaB,2 ! thetas in NHS
         nThetaS =nThetaN+1         ! thetas in SHS
         nThetaHS=nThetaHS+1        ! thetas in one HS
         f1p(nThetaHS)=ft1(nThetaN)+ft1(nThetaS) ! Symm
         f1m(nThetaHS)=ft1(nThetaN)-ft1(nThetaS) ! ASymm
      end do

      !-- Start with first two thetas in first theta block:
      !   This initalizes the flm1 and flm2
      if ( nThetaStart == 1 ) then
         do lm1=1,lmMax-1,2
            lm2=lm1+1
            flm1(lm1)=f1p(1)*wPlm(lm1,1)+f1p(2)*wPlm(lm1,2)
            flm1(lm2)=f1m(1)*wPlm(lm2,1)+f1m(2)*wPlm(lm2,2)
         end do
         if ( lm2 < lmMax ) then
            lm1=lmMax
            flm1(lm1)=f1p(1)*wPlm(lm1,1)+f1p(2)*wPlm(lm1,2)
         end if
          
         if ( sizeThetaB <= 4 ) return
           
      end if

      !-- Following values of n_theta=3,4,5... when n_theta_min=1
      !   or all values of n_theta when n_theta_min > 0
       
      nThetaMin=1
      If ( nThetaStart == 1 ) nThetaMin=3 ! 2 already done

      !-- Calculating position for all thetas:
      !     (nThetaStart-1)/2 are previous blocks,
      !     (nThetaMin-1) is what has been done already
      nThetaC1=(nThetaStart-1)/2+nThetaMin-2

      !-- Loop over half of the thetas with step 2 unrolling:
      do nTheta1=nThetaMin,sizeThetaB/2,2
         nTheta2=nTheta1+1
         nThetaC1=nThetaC1+2
         nThetaC2=nThetaC1+1

         do lm1=1,lmMax-1,2
            lm2=lm1+1
            flm1(lm1)=flm1(lm1) + f1p(nTheta1)*wPlm(lm1,nThetaC1) + &
                                  f1p(nTheta2)*wPlm(lm1,nThetaC2)
            flm1(lm2)=flm1(lm2) + f1m(nTheta1)*wPlm(lm2,nThetaC1) + &
                                  f1m(nTheta2)*wPlm(lm2,nThetaC2)
         end do
         if ( lm2 < lmMax ) then
            lm1=lmMax
            flm1(lm1)=flm1(lm1) + f1p(nTheta1)*wPlm(lm1,nThetaC1) + &
                                  f1p(nTheta2)*wPlm(lm1,nThetaC2)
         end if
           
      end do
        

   end subroutine legTFAS
!------------------------------------------------------------------------
   subroutine legTFAS2(flm1,flm2,ft1,ft2,lmMax,nThetaStart,sizeThetaB)
      !  +-------------+----------------+------------------------------------+
      !  |                                                                   |
      !  |  Legendre transform (n_r,n_theta,m) to (n_r,l,m)                  |
      !  |  [grid to spectral] for 2 arrays                                  |
      !  |  ancl1/2 (input) to flm1/2 (output)                               |
      !  |  One call to this routine does part of the transform              |
      !  |  by summation over theta points in on theta block:                |
      !  |      n_theta_min,..,n_theta_min+n_theta_block-1                   |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+

      !-- Input variables
      integer,      intent(in) :: lmMax          ! Number of modes to be processed
      integer,      intent(in) :: nThetaStart    ! First no of theta on block
      integer,      intent(in) :: sizeThetaB     ! Size of theta block
      real(kind=8), intent(in) :: ft1(*),ft2(*)

      !-- Output: transformed arrays anlc1,anlc2
      real(kind=8), intent(out) :: flm1(*),flm2(*)

      !-- Local variables:
      integer :: nThetaN,nThetaS
      integer :: nThetaHS
      integer :: nTheta1,nTheta2
      integer :: nThetaC1,nThetaC2
      integer :: nThetaMin
      integer :: lm1,lm2

      real(kind=8) :: f1p(n_theta_max/2),f1m(n_theta_max/2)
      real(kind=8) :: f2p(n_theta_max/2),f2m(n_theta_max/2)

      !-- Prepare arrays of sums and differences:
      nThetaHS=0
      do nThetaN=1,sizeThetaB,2 ! thetas in NHS
         nThetaS =nThetaN+1         ! thetas in SHS
         nThetaHS=nThetaHS+1        ! thetas in one HS
         f1p(nThetaHS)=ft1(nThetaN)+ft1(nThetaS) ! Symm
         f1m(nThetaHS)=ft1(nThetaN)-ft1(nThetaS) ! ASymm
         f2p(nThetaHS)=ft2(nThetaN)+ft2(nThetaS)
         f2m(nThetaHS)=ft2(nThetaN)-ft2(nThetaS)
      end do

      !-- Start with first two thetas in first theta block:
      !   This initalizes the flm1 and flm2
      if ( nThetaStart == 1 ) then
         do lm1=1,lmMax-1,2
            lm2=lm1+1
            flm1(lm1)=f1p(1)*wPlm(lm1,1)+f1p(2)*wPlm(lm1,2)
            flm1(lm2)=f1m(1)*wPlm(lm2,1)+f1m(2)*wPlm(lm2,2)
         end do
         do lm1=1,lmMax-1,2
            lm2=lm1+1
            flm2(lm1)=f2p(1)*wPlm(lm1,1)+f2p(2)*wPlm(lm1,2)
            flm2(lm2)=f2m(1)*wPlm(lm2,1)+f2m(2)*wPlm(lm2,2)
         end do
         if ( lm2 < lmMax ) then
            lm1=lmMax
            flm1(lm1)=f1p(1)*wPlm(lm1,1)+f1p(2)*wPlm(lm1,2)
            flm2(lm1)=f2p(1)*wPlm(lm1,1)+f2p(2)*wPlm(lm1,2)
         end if
          
         if ( sizeThetaB <= 4 ) return
           
      end if

      !-- Following values of n_theta=3,4,5... when n_theta_min=1
      !   or all values of n_theta when n_theta_min > 0
      nThetaMin=1
      if ( nThetaStart == 1 ) nThetaMin=3 ! 2 already done

      !-- Calculating position for all thetas:
      !     (nThetaStart-1)/2 are previous blocks,
      !     (nThetaMin-1) is what has been done already
      nThetaC1=(nThetaStart-1)/2+nThetaMin-2

      !-- Loop over half of the thetas with step 2 unrolling:
      do nTheta1=nThetaMin,sizeThetaB/2,2
         nTheta2=nTheta1+1
         nThetaC1=nThetaC1+2
         nThetaC2=nThetaC1+1

         do lm1=1,lmMax-1,2
            lm2=lm1+1
            flm1(lm1)=flm1(lm1) + f1p(nTheta1)*wPlm(lm1,nThetaC1) + &
                                  f1p(nTheta2)*wPlm(lm1,nThetaC2)
            flm1(lm2)=flm1(lm2) + f1m(nTheta1)*wPlm(lm2,nThetaC1) + &
                                  f1m(nTheta2)*wPlm(lm2,nThetaC2)
         end do
         do lm1=1,lmMax-1,2
            lm2=lm1+1
            flm2(lm1)=flm2(lm1) + f2p(nTheta1)*wPlm(lm1,nThetaC1) + &
                                  f2p(nTheta2)*wPlm(lm1,nThetaC2)
            flm2(lm2)=flm2(lm2) + f2m(nTheta1)*wPlm(lm2,nThetaC1) + &
                                  f2m(nTheta2)*wPlm(lm2,nThetaC2)
         end do
         if ( lm2 < lmMax ) then
            lm1=lmMax
            flm1(lm1)=flm1(lm1) + f1p(nTheta1)*wPlm(lm1,nThetaC1) + &
                                  f1p(nTheta2)*wPlm(lm1,nThetaC2)
            flm2(lm1)=flm2(lm1) + f2p(nTheta1)*wPlm(lm1,nThetaC1) + &
                                  f2p(nTheta2)*wPlm(lm1,nThetaC2)
         end if
           
      end do
        
   end subroutine legTFAS2
!------------------------------------------------------------------------
end module legendre_grid_to_spec
