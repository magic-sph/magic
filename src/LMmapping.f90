module LMmapping

   use precision_mod
   use truncation, only: l_axi
   use mem_alloc, only: bytes_allocated
   use parallel_mod, only: load

   implicit none
 
   private
 
   type, public :: mappings
      integer :: l_max, m_max, lm_max ,m_min
      integer, allocatable :: lm2(:,:),lm2l(:),lm2m(:)
      integer, allocatable :: lm2lmS(:), lm2lmA(:)
   end type mappings
 
   type, public :: subblocks_mappings
      integer :: nLMBs,l_max,m_max,sizeLMB2max,m_min
      integer, allocatable :: nLMBs2(:)
      integer, allocatable :: sizeLMB2(:,:)
      integer, allocatable :: lm22lm(:,:,:)
      integer, allocatable :: lm22l(:,:,:)
      integer, allocatable :: lm22m(:,:,:)
   end type subblocks_mappings
 
   public :: allocate_mappings, deallocate_mappings, &
   &         allocate_subblocks_mappings, deallocate_subblocks_mappings

contains

   subroutine allocate_mappings(self,l_max,m_min,m_max,lm_max)

      type(mappings) :: self
      integer, intent(in) :: l_max, m_min, lm_max, m_max

      self%l_max = l_max
      if ( .not. l_axi ) then
         self%m_max = m_max
         self%m_min = m_min
      else
         self%m_max = 0
         self%m_min = 0
      end if
      self%lm_max = lm_max

      allocate( self%lm2(0:l_max,0:l_max),self%lm2l(lm_max),self%lm2m(lm_max) )
      allocate( self%lm2lmS(lm_max),self%lm2lmA(lm_max) )
      bytes_allocated = bytes_allocated + &
      &                 ((l_max+1)*(l_max+1)+4*lm_max)*SIZEOF_INTEGER

   end subroutine allocate_mappings
!-------------------------------------------------------------------------------
   subroutine deallocate_mappings(self)

      type(mappings) :: self

      deallocate( self%lm2, self%lm2l, self%lm2m, self%lm2lmS, self%lm2lmA )

   end subroutine deallocate_mappings
!-------------------------------------------------------------------------------
   subroutine allocate_subblocks_mappings(self,map,nLMBs,l_max,m_min,m_max,lm_balance)

      !-- Input variables
      type(subblocks_mappings) :: self
      type(mappings), intent(in) :: map
      integer,        intent(in) :: nLMBs, l_max, m_min, m_max
      type(load),     intent(in) :: lm_balance(0:nLMBs-1)

      !-- Local variables
      integer :: n_proc,lm1,l1,max_size_of_subblock,lmStart,lmStop
      integer :: counter(0:l_max)

      self%nLMBs = nLMBs
      self%l_max = l_max

      if ( .not. l_axi ) then
         self%m_max = m_max
         self%m_min = m_min
      else
         self%m_max = 0
         self%m_min = 0
      end if

      ! now determine the maximal size of a subblock (was sizeLMB2max parameter
      ! in former versions).
      max_size_of_subblock=0
      do n_proc=0,nLMBs-1
         lmStart=lm_balance(n_proc)%nStart
         lmStop =lm_balance(n_proc)%nStop
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
      bytes_allocated = bytes_allocated +       &
      &                 (nLMBs+(l_max+1)*nLMBs+ &
      &                 3*(l_max+1)*nLMBS*self%sizeLMB2max)*SIZEOF_INTEGER

   end subroutine allocate_subblocks_mappings
!-------------------------------------------------------------------------------
   subroutine deallocate_subblocks_mappings(self)

      type(subblocks_mappings) :: self

      deallocate( self%nLMBs2, self%sizeLMB2, self%lm22lm, self%lm22l )
      deallocate( self%lm22m )

   end subroutine deallocate_subblocks_mappings
!-------------------------------------------------------------------------------
end module LMmapping
