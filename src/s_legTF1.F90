!$Id$
!***********************************************************************
#include "perflib_preproc.cpp"
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

   use truncation, only: n_m_max, nrp, lmP_max
   use blocking, only: nfs, sizeThetaB
   use horizontal_data, only: lStartP, wPlm, lmOddP, lStopP

   implicit none

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
