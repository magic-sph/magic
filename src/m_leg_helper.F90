module leg_helper_mod

   implicit none
   
   type :: leg_helper_t
      !----- Help arrays for Legendre transform calculated in legPrepG:
      !      Parallelizatio note: these are the R-distributed versions
      !      of the field scalars.
      complex(kind=8), allocatable :: dLhw(:), dLhdw(:), dLhz(:), dLhb(:), dLhj(:)
      complex(kind=8), allocatable :: vhG(:), vhC(:), dvhdrG(:), dvhdrC(:)
      complex(kind=8), allocatable :: bhG(:), bhC(:), cbhG(:), cbhC(:)
      !----- R-distributed versions of scalar fields (see c_fields.f):
      complex(kind=8), allocatable :: sR(:), dsR(:), preR(:), dpR(:)
      real(kind=8), allocatable :: zAS(:), dzAS(:), ddzAS(:) ! used in TO
      real(kind=8) :: omegaIC,omegaMA
      complex(kind=8), allocatable :: bCMB(:)
 
   contains
 
      procedure :: initialize
 
   end type leg_helper_t

contains

   subroutine initialize(this,lm_max,lm_maxMag,l_max)

      class(leg_helper_t) :: this
      integer,intent(in) :: lm_max,lm_maxMag,l_max

      allocate( this%dLhw(lm_max) )
      allocate( this%dLhdw(lm_max) )
      allocate( this%dLhz(lm_max) )
      allocate( this%dLhb(lm_max) )
      allocate( this%dLhj(lm_max) )
      allocate( this%vhG(lm_max) )
      allocate( this%vhC(lm_max) )
      allocate( this%dvhdrG(lm_max) )
      allocate( this%dvhdrC(lm_max) )
      allocate( this%bhG(lm_maxMag) )
      allocate( this%bhC(lm_maxMag) )
      allocate( this%cbhG(lm_maxMag) )
      allocate( this%cbhC(lm_maxMag) )
      !----- R-distributed versions of scalar fields (see c_fields.f):
      allocate( this%sR(lm_max),this%dsR(lm_max) )
      allocate( this%preR(lm_max),this%dpR(lm_max) )
      allocate( this%zAS(l_max+1),this%dzAS(l_max+1),&
           & this%ddzAS(l_max+1) ) ! used in TO

      allocate( this%bCMB(lm_maxMag) )

   end subroutine initialize
!------------------------------------------------------------------------------
end module leg_helper_mod
