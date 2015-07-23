!$Id$
#include "intrinsic_sizes.h"
module general_arrays_mod
 
   implicit none
 
   private
 
   type, public, abstract :: general_arrays_t
 
   end type general_arrays_t
 
end module general_arrays_mod
!----------------------------------------------------------------------------
module grid_space_arrays_mod

   use general_arrays_mod
   use truncation, only: nrp
   use blocking, only: nfs

   implicit none

   type, extends(general_arrays_t) :: grid_space_arrays_t
      !----- Nonlinear terms in phi/theta space: 
      real(kind=8), allocatable :: Advr(:,:), Advt(:,:), Advp(:,:)
      real(kind=8), allocatable :: LFr(:,:),  LFt(:,:),  LFp(:,:)
      real(kind=8), allocatable :: VxBr(:,:), VxBt(:,:), VxBp(:,:)
      real(kind=8), allocatable :: VSr(:,:),  VSt(:,:),  VSp(:,:)
      real(kind=8), allocatable :: ViscHeat(:,:), OhmLoss(:,:)

      !----- Fields calculated from these help arrays by legtf:
      real(kind=8), allocatable :: vrc(:,:), vtc(:,:), vpc(:,:)
      real(kind=8), allocatable :: dvrdrc(:,:), dvtdrc(:,:), dvpdrc(:,:)
      real(kind=8), allocatable :: cvrc(:,:), sc(:,:), drSc(:,:)
      real(kind=8), allocatable :: dvrdtc(:,:), dvrdpc(:,:)
      real(kind=8), allocatable :: dvtdpc(:,:), dvpdpc(:,:)
      real(kind=8), allocatable :: brc(:,:), btc(:,:), bpc(:,:)
      real(kind=8), allocatable :: cbrc(:,:), cbtc(:,:), cbpc(:,:)
      real(kind=8), allocatable :: pc(:,:)
      real(kind=8), allocatable :: dsdtc(:,:), dsdpc(:,:)

   contains

      procedure :: initialize
      procedure :: finalize
      procedure :: output
      procedure :: output_nl_input

   end type grid_space_arrays_t

contains

   subroutine initialize(this)

      class(grid_space_arrays_t) :: this
      integer :: size_in_bytes

      allocate( this%Advr(nrp,nfs) )
      allocate( this%Advt(nrp,nfs) )
      allocate( this%Advp(nrp,nfs) )
      allocate( this%LFr(nrp,nfs) )
      allocate( this%LFt(nrp,nfs) )
      allocate( this%LFp(nrp,nfs) )
      allocate( this%VxBr(nrp,nfs) )
      allocate( this%VxBt(nrp,nfs) )
      allocate( this%VxBp(nrp,nfs) )
      allocate( this%VSr(nrp,nfs) )
      allocate( this%VSt(nrp,nfs) )
      allocate( this%VSp(nrp,nfs) )
      allocate( this%ViscHeat(nrp,nfs) )
      allocate( this%OhmLoss(nrp,nfs) )
      size_in_bytes=14*nrp*nfs*SIZEOF_DOUBLE_PRECISION

      !----- Fields calculated from these help arrays by legtf:
      allocate( this%vrc(nrp,nfs),this%vtc(nrp,nfs),this%vpc(nrp,nfs) )
      allocate( this%dvrdrc(nrp,nfs),this%dvtdrc(nrp,nfs) )
      allocate( this%dvpdrc(nrp,nfs),this%cvrc(nrp,nfs) )
      allocate( this%dvrdtc(nrp,nfs),this%dvrdpc(nrp,nfs) )
      allocate( this%dvtdpc(nrp,nfs),this%dvpdpc(nrp,nfs) )
      allocate( this%brc(nrp,nfs),this%btc(nrp,nfs),this%bpc(nrp,nfs) )
      this%btc=1.0d50
      this%bpc=1.0d50
      allocate( this%cbrc(nrp,nfs),this%cbtc(nrp,nfs),this%cbpc(nrp,nfs) )
      allocate( this%sc(nrp,nfs),this%drSc(nrp,nfs) )
      allocate( this%pc(nrp,nfs) )
      allocate( this%dsdtc(nrp,nfs),this%dsdpc(nrp,nfs) )
      size_in_bytes=size_in_bytes + 21*nrp*nfs*SIZEOF_DOUBLE_PRECISION
      !write(*,"(A,I15,A)") "grid_space_arrays: allocated ",size_in_bytes,"B."

   end subroutine initialize
!----------------------------------------------------------------------------
   subroutine finalize(this)

      class(grid_space_arrays_t) :: this
      integer :: size_in_bytes

      deallocate( this%Advr )
      deallocate( this%Advt )
      deallocate( this%Advp )
      deallocate( this%LFr )
      deallocate( this%LFt )
      deallocate( this%LFp )
      deallocate( this%VxBr )
      deallocate( this%VxBt )
      deallocate( this%VxBp )
      deallocate( this%VSr )
      deallocate( this%VSt )
      deallocate( this%VSp )
      deallocate( this%ViscHeat )
      deallocate( this%OhmLoss )
      size_in_bytes=14*nrp*nfs*SIZEOF_DOUBLE_PRECISION

      !----- Fields calculated from these help arrays by legtf:
      deallocate( this%vrc,this%vtc,this%vpc )
      deallocate( this%dvrdrc,this%dvtdrc )
      deallocate( this%dvpdrc,this%cvrc )
      deallocate( this%dvrdtc,this%dvrdpc )
      deallocate( this%dvtdpc,this%dvpdpc )
      deallocate( this%brc,this%btc,this%bpc )
      deallocate( this%cbrc,this%cbtc,this%cbpc )
      deallocate( this%sc,this%drSc )
      deallocate( this%pc )
      deallocate( this%dsdtc, this%dsdpc )
      size_in_bytes=size_in_bytes + 21*nrp*nfs*SIZEOF_DOUBLE_PRECISION
      write(*,"(A,I15,A)") "grid_space_arrays: deallocated ",size_in_bytes,"B."

   end subroutine finalize
!----------------------------------------------------------------------------
   subroutine output(this)

      class(grid_space_arrays_t) :: this
   
      write(*,"(A,3ES20.12)") "Advr,Advt,Advp = ",sum(this%Advr), &
                                   sum(this%Advt),sum(this%Advp)

   end subroutine output
!----------------------------------------------------------------------------
   subroutine output_nl_input(this)

      class(grid_space_arrays_t) :: this
   
      write(*,"(A,6ES20.12)") "vr,vt,vp = ",sum(this%vrc),sum(this%vtc), &
                                            sum(this%vpc)

   end subroutine output_nl_input
!----------------------------------------------------------------------------
end module grid_space_arrays_mod
