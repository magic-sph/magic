MODULE LMmapping
  implicit none

  PRIVATE
  TYPE,PUBLIC :: mappings
     INTEGER :: l_max,lm_max,lmP_max
     INTEGER,ALLOCATABLE :: lm2(:,:),lm2l(:),lm2m(:)
     INTEGER,ALLOCATABLE :: lm2mc(:),l2lmAS(:)      
     INTEGER,ALLOCATABLE :: lm2lmS(:),lm2lmA(:)     
                                                    
     INTEGER,ALLOCATABLE :: lmP2(:,:),lmP2l(:)      
     INTEGER,ALLOCATABLE :: lmP2lmPS(:),lmP2lmPA(:) 
                                                    
     INTEGER,ALLOCATABLE :: lm2lmP(:),lmP2lm(:)     

  END TYPE mappings

  TYPE,PUBLIC :: subblocks_mappings
     INTEGER :: nLMBs,l_max,sizeLMB2max
     INTEGER,ALLOCATABLE :: nLMBs2(:)
     INTEGER,ALLOCATABLE :: sizeLMB2(:,:)
     INTEGER,ALLOCATABLE :: lm22lm(:,:,:)
     INTEGER,ALLOCATABLE :: lm22l(:,:,:)
     INTEGER,ALLOCATABLE :: lm22m(:,:,:)
  END TYPE subblocks_mappings

  !INTEGER :: sizeLMB2max

  PUBLIC :: allocate_mappings,allocate_subblocks_mappings!,sizeLMB2max
contains
  SUBROUTINE allocate_mappings(self,l_max,lm_max,lmP_max)
    TYPE(mappings) :: self
    INTEGER :: l_max,lm_max,lmP_max

    self%l_max = l_max
    self%lm_max = lm_max
    self%lmP_max = lmP_max

    ALLOCATE( self%lm2(0:l_max,0:l_max),self%lm2l(lm_max),self%lm2m(lm_max) )
    ALLOCATE( self%lm2mc(lm_max),self%l2lmAS(0:l_max) )
    ALLOCATE( self%lm2lmS(lm_max),self%lm2lmA(lm_max) )

    ALLOCATE( self%lmP2(0:l_max+1,0:l_max+1),self%lmP2l(lmP_max) )
    ALLOCATE( self%lmP2lmPS(lmP_max),self%lmP2lmPA(lmP_max) )

    ALLOCATE( self%lm2lmP(lm_max),self%lmP2lm(lmP_max) )

  END SUBROUTINE allocate_mappings

  SUBROUTINE allocate_subblocks_mappings(self,map,nLMBs,l_max,lmStartB,lmStopB)
    type(subblocks_mappings) :: self
    type(mappings) :: map
    INTEGER :: nLMBs,l_max
    integer,dimension(nLMBs) :: lmStartB,lmStopB

    integer :: nLMB,lm1,l1,max_size_of_subblock,lmStart,lmStop
    integer :: counter(0:l_max)

    self%nLMBs = nLMBs
    self%l_max = l_max

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

    ALLOCATE( self%nLMBs2(nLMBs),self%sizeLMB2(l_max+1,nLMBs) )
    ALLOCATE( self%lm22lm(self%sizeLMB2max,l_max+1,nLMBs) )
    ALLOCATE( self%lm22l(self%sizeLMB2max,l_max+1,nLMBs) )
    ALLOCATE( self%lm22m(self%sizeLMB2max,l_max+1,nLMBs) )
  END SUBROUTINE allocate_subblocks_mappings

END MODULE LMmapping
