module LMmapping

   use precision_mod
   use truncation, only: l_axi
   use mem_alloc, only: bytes_allocated

   implicit none
 
   private
 
   type, public :: mappings
      integer :: l_max,m_max,lm_max,lmP_max
      integer, allocatable :: lm2(:,:),lm2l(:),lm2m(:)
      integer, allocatable :: lm2mc(:),l2lmAS(:)      
      integer, allocatable :: lm2lmS(:),lm2lmA(:)     
                                                     
      integer, allocatable :: lmP2(:,:),lmP2l(:),lmP2m(:)
      integer, allocatable :: lmP2lmPS(:),lmP2lmPA(:) 
                                                     
      integer, allocatable :: lm2lmP(:),lmP2lm(:)     
 
   end type mappings
 
   type, public :: subblocks_mappings
      integer :: nLMBs,l_max,m_max,sizeLMB2max
      integer, allocatable :: nLMBs2(:)
      integer, allocatable :: sizeLMB2(:,:)
      integer, allocatable :: lm22lm(:,:,:)
      integer, allocatable :: lm22l(:,:,:)
      integer, allocatable :: lm22m(:,:,:)
   end type subblocks_mappings
 
   public :: allocate_mappings, deallocate_mappings, &
   &         allocate_subblocks_mappings, deallocate_subblocks_mappings

contains

   subroutine allocate_mappings(self,l_max,lm_max,lmP_max)

      type(mappings) :: self
      integer, intent(in) :: l_max, lm_max, lmP_max

      self%l_max = l_max
      if ( .not. l_axi ) then
         self%m_max = l_max
      else
         self%m_max = 0
      end if
      self%lm_max = lm_max
      self%lmP_max = lmP_max

      allocate( self%lm2(0:l_max,0:l_max),self%lm2l(lm_max),self%lm2m(lm_max) )
      allocate( self%lm2mc(lm_max),self%l2lmAS(0:l_max) )
      allocate( self%lm2lmS(lm_max),self%lm2lmA(lm_max) )
      bytes_allocated = bytes_allocated + &
                        ((l_max+1)*(l_max+1)+5*lm_max+l_max+1)*SIZEOF_INTEGER

      allocate( self%lmP2(0:l_max+1,0:l_max+1),self%lmP2l(lmP_max) )
      allocate( self%lmP2m(lmP_max) )
      allocate( self%lmP2lmPS(lmP_max),self%lmP2lmPA(lmP_max) )
      allocate( self%lm2lmP(lm_max),self%lmP2lm(lmP_max) )
      bytes_allocated = bytes_allocated + &
                        ((l_max+2)*(l_max+2)+5*lmP_max+lm_max)*SIZEOF_INTEGER

   end subroutine allocate_mappings
!-------------------------------------------------------------------------------
   subroutine deallocate_mappings(self)

      type(mappings) :: self

      deallocate( self%lm2, self%lm2l, self%lm2m, self%lm2mc, self%l2lmAS )
      deallocate( self%lm2lmS, self%lm2lmA, self%lmP2, self%lmP2l )
      deallocate( self%lmP2m, self%lmP2lmPS, self%lmP2lmPA, self%lm2lmP )
      deallocate( self%lmP2lm )

   end subroutine deallocate_mappings
!-------------------------------------------------------------------------------
   subroutine allocate_subblocks_mappings(self,map,nLMBs,l_max,lmStartB,lmStopB)

      !-- Input variables
      type(subblocks_mappings) :: self
      type(mappings), intent(in) :: map
      integer,        intent(in) :: nLMBs, l_max
      integer,        intent(in) :: lmStartB(nLMBs), lmStopB(nLMBs)

      !-- Local variables
      integer :: nLMB,lm1,l1,max_size_of_subblock,lmStart,lmStop
      integer :: counter(0:l_max)

      self%nLMBs = nLMBs
      self%l_max = l_max

      if ( .not. l_axi ) then
         self%m_max = l_max
      else
         self%m_max = 0
      end if

      ! now determine the maximal size of a subblock (was sizeLMB2max parameter
      ! in former versions).
      max_size_of_subblock=0
      do nLMB=1,nLMBs
         lmStart=lmStartB(nLMB)
         lmStop =lmStopB(nLMB)
         counter=0
         do lm1=lmStart,lmStop
            l1=map%lm2l(lm1)
            counter(l1) = counter(l1) + 1
         end do
         if (maxval(counter) > max_size_of_subblock) then
            max_size_of_subblock=maxval(counter)
         end if
      end do
      self%sizeLMB2max = max_size_of_subblock
      !write(*,"(A,I5)") "Using sizeLMB2max = ",self%sizeLMB2max

      allocate( self%nLMBs2(nLMBs),self%sizeLMB2(l_max+1,nLMBs) )
      allocate( self%lm22lm(self%sizeLMB2max,l_max+1,nLMBs) )
      allocate( self%lm22l(self%sizeLMB2max,l_max+1,nLMBs) )
      allocate( self%lm22m(self%sizeLMB2max,l_max+1,nLMBs) )
      bytes_allocated = bytes_allocated + &
                        (nLMBs+(l_max+1)*nLMBs+ &
                        3*(l_max+1)*nLMBS*self%sizeLMB2max)*SIZEOF_INTEGER

   end subroutine allocate_subblocks_mappings
!-------------------------------------------------------------------------------
   subroutine deallocate_subblocks_mappings(self)

      type(subblocks_mappings) :: self

      deallocate( self%nLMBs2, self%sizeLMB2, self%lm22lm, self%lm22l )
      deallocate( self%lm22m )

   end subroutine deallocate_subblocks_mappings
!-------------------------------------------------------------------------------
end module LMmapping
