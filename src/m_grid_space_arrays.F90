#include "intrinsic_sizes.h"
MODULE general_arrays_mod
  implicit none

  TYPE,PUBLIC,ABSTRACT :: general_arrays_t
  END TYPE general_arrays_t
END MODULE general_arrays_mod

MODULE grid_space_arrays_mod
  use general_arrays_mod
  USE truncation,ONLY: nrp,ncp
  USE blocking, only: nfs
  implicit none

  TYPE,extends(general_arrays_t) :: grid_space_arrays_t
     !----- Nonlinear terms in phi/theta space: 
     REAL(kind=8),DIMENSION(:,:),ALLOCATABLE :: Advr, Advt, Advp
     REAL(kind=8),DIMENSION(:,:),ALLOCATABLE :: LFr,  LFt,  LFp
     REAL(kind=8),DIMENSION(:,:),ALLOCATABLE :: VxBr, VxBt, VxBp
     REAL(kind=8),DIMENSION(:,:),ALLOCATABLE :: VSr,  VSt,  VSp
     REAL(kind=8),DIMENSION(:,:),ALLOCATABLE :: ViscHeat, OhmLoss

     !----- Fields calculated from these help arrays by legtf:
     COMPLEX(kind=8),DIMENSION(:,:),ALLOCATABLE :: vrc, vtc, vpc
     COMPLEX(kind=8),DIMENSION(:,:),ALLOCATABLE :: dvrdrc, dvtdrc
     COMPLEX(kind=8),DIMENSION(:,:),ALLOCATABLE :: dvpdrc, cvrc
     COMPLEX(kind=8),DIMENSION(:,:),ALLOCATABLE :: dvrdtc, dvrdpc
     COMPLEX(kind=8),DIMENSION(:,:),ALLOCATABLE :: dvtdpc, dvpdpc
     COMPLEX(kind=8),DIMENSION(:,:),ALLOCATABLE :: brc, btc, bpc
     COMPLEX(kind=8),DIMENSION(:,:),ALLOCATABLE :: cbrc, cbtc, cbpc
     COMPLEX(kind=8),DIMENSION(:,:),ALLOCATABLE :: sc, drSc
     COMPLEX(kind=8),DIMENSION(:,:),ALLOCATABLE :: pc
     COMPLEX(kind=8),DIMENSION(:,:),ALLOCATABLE :: dsdtc, dsdpc
   CONTAINS
     procedure :: initialize
     procedure :: finalize
     procedure :: output
     procedure :: output_nl_input
  END TYPE grid_space_arrays_t
CONTAINS
  SUBROUTINE initialize(this)
    class(grid_space_arrays_t) :: this

    INTEGER :: size_in_bytes
    ALLOCATE( this%Advr(nrp,nfs) )
    ALLOCATE( this%Advt(nrp,nfs) )
    ALLOCATE( this%Advp(nrp,nfs) )
    ALLOCATE( this%LFr(nrp,nfs) )
    ALLOCATE( this%LFt(nrp,nfs) )
    ALLOCATE( this%LFp(nrp,nfs) )
    ALLOCATE( this%VxBr(nrp,nfs) )
    ALLOCATE( this%VxBt(nrp,nfs) )
    ALLOCATE( this%VxBp(nrp,nfs) )
    ALLOCATE( this%VSr(nrp,nfs) )
    ALLOCATE( this%VSt(nrp,nfs) )
    ALLOCATE( this%VSp(nrp,nfs) )
    ALLOCATE( this%ViscHeat(nrp,nfs) )
    ALLOCATE( this%OhmLoss(nrp,nfs) )
    size_in_bytes=14*nrp*nfs*SIZEOF_DOUBLE_PRECISION

    !----- Fields calculated from these help arrays by legtf:
    ALLOCATE( this%vrc(ncp,nfs),this%vtc(ncp,nfs),this%vpc(ncp,nfs) )
    ALLOCATE( this%dvrdrc(ncp,nfs),this%dvtdrc(ncp,nfs) )
    ALLOCATE( this%dvpdrc(ncp,nfs),this%cvrc(ncp,nfs) )
    ALLOCATE( this%dvrdtc(ncp,nfs),this%dvrdpc(ncp,nfs) )
    ALLOCATE( this%dvtdpc(ncp,nfs),this%dvpdpc(ncp,nfs) )
    ALLOCATE( this%brc(ncp,nfs),this%btc(ncp,nfs),this%bpc(ncp,nfs) )
    this%btc=1.0d50
    this%bpc=1.0d50
    ALLOCATE( this%cbrc(ncp,nfs),this%cbtc(ncp,nfs),this%cbpc(ncp,nfs) )
    ALLOCATE( this%sc(ncp,nfs),this%drSc(ncp,nfs) )
    ALLOCATE( this%pc(ncp,nfs) )
    ALLOCATE( this%dsdtc(ncp,nfs),this%dsdpc(ncp,nfs) )
    size_in_bytes=size_in_bytes + 21*ncp*nfs*SIZEOF_DOUBLE_COMPLEX
    !WRITE(*,"(A,I15,A)") "grid_space_arrays: allocated ",size_in_bytes,"B."
  END SUBROUTINE initialize

  SUBROUTINE finalize(this)
    class(grid_space_arrays_t) :: this

    INTEGER :: size_in_bytes
    DEALLOCATE( this%Advr )
    DEALLOCATE( this%Advt )
    DEALLOCATE( this%Advp )
    DEALLOCATE( this%LFr )
    DEALLOCATE( this%LFt )
    DEALLOCATE( this%LFp )
    DEALLOCATE( this%VxBr )
    DEALLOCATE( this%VxBt )
    DEALLOCATE( this%VxBp )
    DEALLOCATE( this%VSr )
    DEALLOCATE( this%VSt )
    DEALLOCATE( this%VSp )
    DEALLOCATE( this%ViscHeat )
    DEALLOCATE( this%OhmLoss )
    size_in_bytes=14*nrp*nfs*SIZEOF_DOUBLE_PRECISION

    !----- Fields calculated from these help arrays by legtf:
    DEALLOCATE( this%vrc,this%vtc,this%vpc )
    DEALLOCATE( this%dvrdrc,this%dvtdrc )
    DEALLOCATE( this%dvpdrc,this%cvrc )
    DEALLOCATE( this%dvrdtc,this%dvrdpc )
    DEALLOCATE( this%dvtdpc,this%dvpdpc )
    DEALLOCATE( this%brc,this%btc,this%bpc )
    DEALLOCATE( this%cbrc,this%cbtc,this%cbpc )
    DEALLOCATE( this%sc,this%drSc )
    DEALLOCATE( this%pc )
    DEALLOCATE( this%dsdtc, this%dsdpc )
    size_in_bytes=size_in_bytes + 21*ncp*nfs*SIZEOF_DOUBLE_COMPLEX
    WRITE(*,"(A,I15,A)") "grid_space_arrays: deallocated ",size_in_bytes,"B."
  END SUBROUTINE finalize

  SUBROUTINE output(this)
    class(grid_space_arrays_t) :: this
    
    WRITE(*,"(A,3ES20.12)") "Advr,Advt,Advp = ",SUM(this%Advr),SUM(this%Advt),SUM(this%Advp)
  END SUBROUTINE output

  SUBROUTINE output_nl_input(this)
    class(grid_space_arrays_t) :: this
    
    WRITE(*,"(A,6ES20.12)") "vr,vt,vp = ",SUM(this%vrc),SUM(this%vtc),SUM(this%vpc)
  END SUBROUTINE output_nl_input
END MODULE grid_space_arrays_mod
