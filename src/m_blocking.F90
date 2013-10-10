!$Id$
!********************************************************************
!  Common block containing blocking information
!********************************************************************

!    !------------ This is release 1 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

MODULE blocking
  USE logic
  USE parallel_mod
  USE truncation
  USE output_data
  USE LMmapping,ONLY: mappings,allocate_mappings,&
       & allocate_subblocks_mappings,subblocks_mappings!,sizeLMB2max
  IMPLICIT NONE

  !------------------------------------------------------------------------
  !  Dividing loops over all spherical harmonics into blocks that
  !  contain approx. nChunk harmonics . The number of blocks, nLMBs,   
  !  is a multiple of nThreadsUse (the number of processors used).
  !  Parameter nChunk controls the number (and size) of blocks.
  !  When nThreadUse is large, the size of the blocks may be 
  !  considerably smaller than the chosen nChunk,
  !  since nLMBs must be a multiple of nThreadsUse!

  INTEGER,PARAMETER :: nChunk=512
  INTEGER :: nThreadsMax
  ! nthreads > 1
  INTEGER,POINTER :: lm2(:,:),lm2l(:),lm2m(:)
  INTEGER,POINTER :: lm2mc(:),l2lmAS(:)
  INTEGER,POINTER :: lm2lmS(:),lm2lmA(:)

  INTEGER,POINTER :: lmP2(:,:),lmP2l(:)
  INTEGER,POINTER :: lmP2lmPS(:),lmP2lmPA(:)

  INTEGER,POINTER :: lm2lmP(:),lmP2lm(:)

  
  TYPE(mappings),target :: st_map
  TYPE(mappings),target :: lo_map
  !TYPE(mappings),TARGET :: sn_map

  !INTEGER :: nLMBsMax
  INTEGER :: nLMBs,sizeLMB

  INTEGER,ALLOCATABLE :: lmStartB(:),lmStopB(:)
  !INTEGER,PARAMETER :: sizeLMB2max=201

  INTEGER,POINTER :: nLMBs2(:),sizeLMB2(:,:)
  INTEGER,POINTER :: lm22lm(:,:,:)
  INTEGER,POINTER :: lm22l(:,:,:)
  INTEGER,POINTER :: lm22m(:,:,:)

  TYPE(subblocks_mappings),TARGET :: st_sub_map, lo_sub_map,sn_sub_map

  INTEGER :: sizeRB


  !------------------------------------------------------------------------
  !  Following divides loops over points in theta-direction (index ic) into
  !  blocks. Enhances performance by trying to decrease memory access
  !  but is not relevant for SMP parallel processing.
  !  It is tried to devide the theta loop into parts whos data
  !  fit into cache. Thus ideal block sizes (nfs) highly depend on the
  !  used computer.
  !  The value nfs used here has been determined experientally
  !  by uli Christensen for an IBM SP2. It represents an upper bound.
  !  The real block size sizeThetaB used in the code is determined in
  !  s_prep.f by decreasing the size starting with nfs,
  !  until n_theta_max is a multiple of sizeThetaB. The maximum number of
  !  blocks is limited to 8 here (see s_prep.f).
  !  To find a new blocking for a different computer I suggest to
  !  vary the number of blocks nThetaBs (see s_prep.f) for a fixed
  !  resolution. If the ideal no. of blocks nThetaBsI has been found
  !  this determins the ideal block size:
  !         sizeThetaBI=((n_theta_max-1)/nThetaBsI+1)*n_phi_max
  !  This can than be rescaled for different resolutions n_phi_max:
  !         nfs=sizeThetaBI/n_phi_max+1
  !  The resulting block number can be kept down by multiplying this
  !  with an integer number nBDown, Uli has used K=8 here.
  !  For the memory access this is equivalent to inceasing the
  !  number of blocks.
  !  I dont know exactly why Uli has the additional 16. It may just
  !  decrease the block size to be on the safe side, which is a good
  !  idea, since you dont want to end up with blocks that are just
  !  a little bit larger than the cache.
  !  Thus it seems a good idea to use
  !        nfs=sizeThetaBI/(n_phi_max+nBSave)+1)*nBDown
  !

  INTEGER :: nfs
  INTEGER,PARAMETER :: sizeThetaBI=284,nBSave=16,nBDown=8

  INTEGER :: nThetaBs,sizeThetaB
contains
  SUBROUTINE initialize_blocking
    integer :: nThreadsAva

    REAL(kind=8) :: load
    INTEGER :: iLoad
    INTEGER :: n
    integer :: LMB_with_l1m0,l1m0,irank

    character(len=255) :: message
    !--- End of declaration

    CALL allocate_mappings(st_map,l_max,lm_max,lmP_max)
    CALL allocate_mappings(lo_map,l_max,lm_max,lmP_max)
    !CALL allocate_mappings(sn_map,l_max,lm_max,lmP_max)

    IF ( (rank.EQ.0).AND.l_save_out ) THEN
       OPEN(nLF,FILE=log_file,STATUS='UNKNOWN',POSITION='APPEND')
    END IF

    !--- Setting and checking thread number
    !    The number of threads can be selected by the input
    !    variable nThreadsRun. If nThreadsRun=0 I check
    !    use the number of processors available.
#ifdef WITHOMP
    nThreadsAva=OMP_GET_NUM_PROCS()
#else
    nThreadsAva=1
#endif
    IF ( nThreadsRun == 0 .OR. nThreadsRun > nThreadsAva ) THEN
       nThreads=nThreadsAva
    ELSE
       nThreads=nThreadsRun
    END IF
#ifdef WITHOMP
    CALL OMP_SET_NUM_THREADS(nThreads)
    WRITE(*,*)
    WRITE(*,*) '! Max thread number available :',nThreadsAva
    WRITE(*,*) '! Max thread number demanded  :',nThreadsRun
    WRITE(*,*) '! Number of threads I will use:',nThreads
    WRITE(*,*)
    WRITE(nLF,*)
    WRITE(nLF,*) '! Max thread number available :',nThreadsAva
    WRITE(nLF,*) '! Max thread number demanded  :',nThreadsRun
    WRITE(nLF,*) '! Number of threads I will use:',nThreads
    WRITE(nLF,*)
#endif


    IF (rank.EQ.0) THEN
       WRITE(message,*) '! Number of ranks I will use:',n_procs
       call logWrite(message)
    END IF

    nLMBs = n_procs
    nLMBs_per_rank = nLMBs/n_procs
    nThreadsMax = 1
    !PRINT*,"nLMBs first = ",nLMBs
    sizeLMB=(lm_max-1)/nLMBs+1
    !nLMBs = nLMBs - (nLMBs*sizeLMB-lm_max)/sizeLMB

    IF ( nLMBs*sizeLMB > lm_max ) THEN
       WRITE(message,*) '! Uneven load balancing in LM blocks!'
       call logWrite(message)
       load=DBLE(lm_max-(nLMBs-1)*sizeLMB)/sizeLMB
       WRITE(message,*) '! Load percentage of last block:',load*100.D0
       call logWrite(message)
       iLoad=INT(load,4)
       IF ( iLoad >= 1 ) THEN
          WRITE(*,*) '! No. of redundant blocks:',iLoad
          nLMBs=nLMBs-iLoad
       END IF
    END IF
    !PRINT*,"nLMBs final = ",nLMBs
    ALLOCATE( lmStartB(nLMBs),lmStopB(nLMBs) )
    !ALLOCATE( nLMBs2(nLMBs),sizeLMB2(l_max+1,nLMBs) )
    !ALLOCATE( lm22lm(sizeLMB2max,l_max+1,nLMBs) )
    !ALLOCATE( lm22l(sizeLMB2max,l_max+1,nLMBs) )
    !ALLOCATE( lm22m(sizeLMB2max,l_max+1,nLMBs) )

    nfs=(sizeThetaBI/(n_phi_tot+nBSave)+1) * nBDown


    IF ( l_RMS .AND. nThreads > nThreadsMax ) THEN
       IF (rank.EQ.0) THEN
          WRITE(*,*) '! Too small value of nThreadsMax !'
          WRITE(*,*) '! for calculating RMS forces!'
          WRITE(*,*) '! See c_RMS.f!'
          WRITE(*,*) '! Increase nThreadsMax in m_blocking.F90!'
       END IF
       STOP
    END IF

    !--- Get radial blocking
    IF ( MOD(n_r_max-1,nThreads) /= 0 ) THEN
       IF (rank.EQ.0) THEN
          WRITE(*,*) 'Number of threads has to be multiple of n_r_max-1!'
          WRITE(*,*) 'nThreads :',nThreads
          WRITE(*,*) 'n_r_max-1:',n_r_max-1
       END IF
       STOP
    END IF
    sizeRB=(n_r_max-1)/nThreads

    !--- Calculate lm and ml blocking:
    DO n=1,nLMBs
       lmStartB(n)=(n-1)*sizeLMB+1
       lmStopB(n) =MIN(n*sizeLMB,lm_max)
       IF ( lmStopB(n) == lm_max ) EXIT
    END DO

    CALL get_standard_lm_blocking(st_map,minc)
    !CALL get_standard_lm_blocking(lo_map,minc)
    IF (n_procs.LE.l_max/2) THEN
       !better load balancing, but only up to l_max/2
       CALL get_snake_lm_blocking(lo_map,minc)
       WRITE(message,*) "Using snake ordering."
    ELSE
       CALL get_lorder_lm_blocking(lo_map,minc)
       WRITE(message,*) "Using lorder ordering."
    END IF
    CALL logWrite(message)

    DO n=1,nLMBs
       WRITE(message,*) n,lmStartB(n),lmStopB(n),lmStopB(n)-lmStartB(n)+1
       CALL logWrite(message)
       IF ( lmStopB(n) == lm_max ) EXIT
    END DO


    !--- Get the block (rank+1) with the l1m0 mode
    l1m0 = lo_map%lm2(1,0)
    DO n=1,nLMBs
       IF ( (l1m0.GE.lmStartB(n)) .AND. (l1m0.LE.lmStopB(n)) ) THEN
          LMB_with_l1m0=n
          exit
       END IF
    END DO

    ! which rank does have the LMB with LMB_with_l1m0?
    do irank=0,n_procs-1
       IF ((LMB_with_l1m0-1.GE.irank*nLMBs_per_rank).AND.&
            &(LMB_with_l1m0-1.LE.(irank+1)*nLMBs_per_rank-1)) THEN
          rank_with_l1m0 = irank
       end if
    end do
    WRITE(message,"(2(A,I4))") "rank no ",rank_with_l1m0," has l1m0 in block ",LMB_with_l1m0
    CALL logWrite(message)
       
    ! set the standard ordering as default
    lm2(0:,0:) => st_map%lm2
    lm2l(1:lm_max) => st_map%lm2l
    lm2m(1:lm_max) => st_map%lm2m
    lm2mc(1:lm_max)=> st_map%lm2mc
    l2lmAS(0:l_max)=> st_map%l2lmAS
    lm2lmS(1:lm_max) => st_map%lm2lmS
    lm2lmA(1:lm_max) => st_map%lm2lmA
    lmP2(0:,0:) => st_map%lmP2
    lmP2l(1:lmP_max) => st_map%lmP2l
    lmP2lmPS(1:lmP_max) => st_map%lmP2lmPS
    lmP2lmPA(1:lmP_max) => st_map%lmP2lmPA
    lm2lmP(1:lm_max) => st_map%lm2lmP
    lmP2lm(1:lmP_max) => st_map%lmP2lm

    CALL allocate_subblocks_mappings(st_sub_map,st_map,nLMBs,l_max,lmStartB,lmStopB)
    CALL allocate_subblocks_mappings(lo_sub_map,lo_map,nLMBs,l_max,lmStartB,lmStopB)
    !CALL allocate_subblocks_mappings(sn_sub_map,sn_map,nLMBs,l_max,lmStartB,lmStopB)

    !--- Getting lm sub-blocks:
    !PRINT*," ---------------- Making the standard subblocks -------------- "
    CALL get_subblocks(st_map, st_sub_map) !nLMBs,lm2l,lm2m, nLMBs2,sizeLMB2,lm22lm,lm22l,lm22m)
    !PRINT*," ---------------- Making the lorder subblocks ---------------- "
    CALL get_subblocks(lo_map, lo_sub_map)
    !PRINT*," ---------------- Making the snake order subblocks ----------- "
    !CALL get_subblocks(sn_map, sn_sub_map)

    ! default mapping
    nLMBs2(1:nLMBs) => st_sub_map%nLMBs2
    sizeLMB2(1:,1:) => st_sub_map%sizeLMB2
    lm22lm(1:,1:,1:) => st_sub_map%lm22lm
    lm22l(1:,1:,1:) => st_sub_map%lm22l
    lm22m(1:,1:,1:) => st_sub_map%lm22m

    !-- Calculate blocking parameters for blocking loops over theta:
    !   This is not relevant for parallelisation so far.
    !   The desired block size is nfs. If n_theta_max is not a
    !   multiple of nfs I use nThetaBs=n_theta_max/sizeThetaB as a first
    !   guess for the number of blocks and then increase this
    !   number adjusting the block size sizeThetaB.
    !        nThetaBs  =n_theta_max/sizeThetaB
    !        DO n=1,100
    !           IF ( nThetaBs*sizeThetaB.EQ.n_theta_max ) GOTO 50  ! done
    !           nThetaBs=nThetaBs+1
    !           sizeThetaB=n_theta_max/nThetaBs
    !        END DO
    !        IF ( nThetaBs*sizeThetaB.EQ.n_theta_max ) GOTO 50
    !        WRITE(*,*)
    !        WRITE(*,*) '! No proper blocking for theta-loops found!'
    !        STOP
    ! 0      CONTINUE ! jump point for succesfull theta blocking

    ! JW 6 Mar 2012: There are problems when nfs.NE.sizeThetaB and
    !                I don't know why!
    !sizeThetaB=nfs
    sizeThetaB=MIN(n_theta_max,nfs)
    nThetaBs  =n_theta_max/sizeThetaB
    IF ( nThetaBs*sizeThetaB /= n_theta_max ) THEN
       WRITE(*,*)
       WRITE(*,*) '! n_theta_max is not multiple of nfs!'
       WRITE(*,*) '! n_theta_max    =',n_theta_max
       WRITE(*,*) '! nfs            =',nfs
       WRITE(*,*) '! n_theta_max/nfs=',n_theta_max/nfs
       WRITE(*,*) '! Please decrease sizeThetaBI or nBDown in m_blocking.F90!'
       STOP
    END IF



    IF (rank.EQ.0) THEN
       WRITE(*,*) '!-- Blocking information:'
       WRITE(*,*)
       WRITE(*,*) '!    Number of LM-blocks:',nLMBs
       WRITE(*,*) '!    Size   of LM-blocks:',sizeLMB
       WRITE(*,*) '!               nChunk  :',nChunk
       WRITE(*,*) '!               nThreads:',nThreads
       WRITE(*,*)
       WRITE(*,*) '! Number of theta blocks:',nThetaBs
       WRITE(*,*) '!   size of theta blocks:',sizeThetaB
       WRITE(*,*) '!       ideal size (nfs):',nfs
       WRITE(nLF,*) '!-- Blocking information:'
       WRITE(nLF,*)
       WRITE(nLF,*) '!    Number of LM-blocks:',nLMBs
       WRITE(nLF,*) '!    Size   of LM-blocks:',sizeLMB
       WRITE(nLF,*) '!               nChunk  :',nChunk
       WRITE(nLF,*) '!               nThreads:',nThreads
       WRITE(nLF,*)
       WRITE(nLF,*) '! Number of theta blocks:',nThetaBs
       WRITE(nLF,*) '!   size of theta blocks:',sizeThetaB
       WRITE(nLF,*) '!       ideal size (nfs):',nfs

       IF ( l_save_out ) CLOSE(nLF)
    END IF



    RETURN
  END SUBROUTINE initialize_blocking

  !------------------------------------------------------------------------

  ! Uses also the module variables lmStartB, lmStopB
  !SUBROUTINE get_subblocks(number_of_blocks,lm2l,lm2m,&
  !     &                   number_of_subblocks,size_of_subblocks,lm22lm,lm22l,lm22m)
  SUBROUTINE get_subblocks(map,sub_map) 
    TYPE(mappings),intent(IN) :: map
    TYPE(subblocks_mappings),intent(INOUT) :: sub_map

    
    ! Local variables
    INTEGER :: number_of_blocks
    logical :: lStop
    integer :: size
    INTEGER :: nB2,n,n2,n3
    INTEGER :: help,help1(lm_max),help2(lm_max),help3(lm_max)
    INTEGER :: lm,l,m
    INTEGER :: check(0:l_max,0:l_max)

    LOGICAL :: DEBUG_OUTPUT=.false.

    number_of_blocks=sub_map%nLMBs
    
    check = 0
    lStop=.FALSE.
    size=0
    nB2=0
    DO n=1,nLMBs
       sub_map%nLMBs2(n)=1
       lm=lmStartB(n)
       !PRINT*,n,": lm = ",lm
       !------ Start first sub-block:
       sub_map%sizeLMB2(1,n) =1
       sub_map%lm22lm(1,1,n) =lm
       sub_map%lm22l(1,1,n)  =map%lm2l(lm)
       sub_map%lm22m(1,1,n)  =map%lm2m(lm)
       DO lm=lmStartB(n)+1,lmStopB(n)
          !WRITE(*,"(4X,A,I4)") "lm = ",lm
          DO n2=1,sub_map%nLMBs2(n)
             !WRITE(*,"(8X,A,I4)") "n2 = ",n2
             IF ( sub_map%lm22l(1,n2,n) == map%lm2l(lm) ) THEN
                !------ Add to old block
                sub_map%sizeLMB2(n2,n)=sub_map%sizeLMB2(n2,n)+1
                GOTO 20
             END IF
          END DO
          !------ Start new l-block:
          n2 = sub_map%nLMBs2(n)+1
          sub_map%nLMBs2(n)     =n2
          sub_map%sizeLMB2(n2,n)=1
20        CONTINUE
          sub_map%lm22lm(sub_map%sizeLMB2(n2,n),n2,n)=lm
          sub_map%lm22l(sub_map%sizeLMB2(n2,n),n2,n) =map%lm2l(lm)
          sub_map%lm22m(sub_map%sizeLMB2(n2,n),n2,n) =map%lm2m(lm)
       END DO

       !------ Resort:
       IF ( sub_map%nLMBs2(n) > 1 ) THEN
          DO n2=1,sub_map%nLMBs2(n)
             DO n3=n2+1,sub_map%nLMBs2(n)
                IF  ( sub_map%lm22m(1,n2,n) > sub_map%lm22m(1,n3,n) ) THEN
                   help=sub_map%sizeLMB2(n2,n)
                   DO lm=1,help
                      help1(lm)=sub_map%lm22l(lm,n2,n)
                      help2(lm)=sub_map%lm22m(lm,n2,n)
                      help3(lm)=sub_map%lm22lm(lm,n2,n)
                   END DO
                   sub_map%sizeLMB2(n2,n)=sub_map%sizeLMB2(n3,n)
                   DO lm=1,sub_map%sizeLMB2(n2,n)
                      sub_map%lm22l(lm,n2,n) =sub_map%lm22l(lm,n3,n)
                      sub_map%lm22m(lm,n2,n) =sub_map%lm22m(lm,n3,n)
                      sub_map%lm22lm(lm,n2,n)=sub_map%lm22lm(lm,n3,n)
                   END DO
                   sub_map%sizeLMB2(n3,n)=help
                   DO lm=1,help
                      sub_map%lm22l(lm,n3,n) =help1(lm)
                      sub_map%lm22m(lm,n3,n) =help2(lm)
                      sub_map%lm22lm(lm,n3,n)=help3(lm)
                   END DO
                END IF
             END DO
          END DO
       END IF

       nB2=nB2+sub_map%nLMBs2(n)
       DO n2=1,sub_map%nLMBs2(n)
          IF ( sub_map%sizeLMB2(n2,n) > sub_map%sizeLMB2max ) THEN
             lStop=.TRUE.
             size=MAX(size,sub_map%sizeLMB2(n2,n))
          END IF
          DO n3=1,sub_map%sizeLMB2(n2,n)
             l=sub_map%lm22l(n3,n2,n)
             m=sub_map%lm22m(n3,n2,n)
             check(l,m)=check(l,m)+1
          END DO
          !             WRITE(99,*) n,n2,sub_map%sizeLMB2(n2,n)
       END DO
       if (DEBUG_OUTPUT) then
          IF (rank.EQ.0) THEN
             WRITE(*,"(4X,2(A,I4))") "Subblocks of Block ",n,"/",nLMBs
             DO n2=1,sub_map%nLMBs2(n)
                WRITE(*,"(8X,3(A,I4))") "subblock no. ",n2,", of ",sub_map%nLMBs2(n)," with size ",sub_map%sizeLMB2(n2,n)
                DO n3=1,sub_map%sizeLMB2(n2,n)
                   WRITE(*,"(10X,2(A,I4),2I4)") "local lm is ",n3,&
                        &" translates into global lm,l,m : ",sub_map%lm22lm(n3,n2,n),sub_map%lm22l(n3,n2,n),sub_map%lm22m(n3,n2,n)
                END DO
             END DO
          END IF
       endif

    END DO

    IF ( lStop ) THEN
       WRITE(*,*) '! Increase sizeLMB2max in m_blocking.F90!'
       WRITE(*,*) '! to at least:',size
       STOP
    END IF

    DO m=0,l_max,minc
       DO l=m,l_max
          IF ( check(l,m) == 0 ) THEN
             WRITE(*,*) 'Warning, forgotten l,m:',l,m,map%lm2(l,m)
             STOP
          ELSE IF ( check(l,m) > 1 ) THEN
             WRITE(*,*) 'Warning, too much l,m:',l,m,check(l,m)
             STOP
          END IF
       END DO
    END DO

  END SUBROUTINE get_subblocks

  SUBROUTINE get_standard_lm_blocking(map,minc)
    type(mappings) :: map
    INTEGER,INTENT(IN) :: minc
    
    ! Local variables
    INTEGER :: m,l,lm,lmP,mc
    INTEGER,dimension(lmP_max) :: lmP2m
    
    DO m=0,map%l_max
       DO l=m,map%l_max
          map%lm2(l,m)  =-1
          map%lmP2(l,m) =-1
          !check(l,m)=0
       END DO
       l=map%l_max+1
       map%lmP2(l,m)=-1
    END DO

    lm =0
    lmP=0
    mc =0
    DO m=0,map%l_max,minc
       mc=mc+1
       !m2mc(m)=mc
       DO l=m,map%l_max
          lm         =lm+1
          map%lm2l(lm)   =l
          map%lm2m(lm)   =m
          map%lm2mc(lm)  =mc
          map%lm2(l,m)   =lm
          IF ( m == 0 ) map%l2lmAS(l)=lm
          lmP        =lmP+1
          map%lmP2l(lmP) = l
          lmP2m(lmP) = m
          map%lmP2(l,m)  =lmP
          !IF ( m == 0 ) l2lmPAS(l)=lmP
          map%lm2lmP(lm) =lmP
          map%lmP2lm(lmP)=lm
       END DO
       l=map%l_max+1    ! Extra l for lmP
       lmP=lmP+1
       map%lmP2l(lmP) =l
       lmP2m(lmP) = m
       map%lmP2(l,m)  =lmP
       !IF ( m == 0 ) l2lmPAS(l)=lmP
       map%lmP2lm(lmP)=-1
    END DO
    IF ( lm /= map%lm_max ) THEN
       WRITE(*,"(2(A,I6))") 'Wrong lm=',lm," != map%lm_max = ",map%lm_max
       STOP
    END IF
    IF ( lmP /= map%lmP_max ) THEN
       WRITE(*,*) 'Wrong lmP!'
       STOP
    END IF
    DO lm=1,map%lm_max
       l=map%lm2l(lm)
       m=map%lm2m(lm)
       IF ( l > 0 .AND. l > m ) THEN
          map%lm2lmS(lm)=map%lm2(l-1,m)
       ELSE
          map%lm2lmS(lm)=-1
       END IF
       IF ( l < map%l_max ) THEN
          map%lm2lmA(lm)=map%lm2(l+1,m)
       ELSE
          map%lm2lmA(lm)=-1
       END IF
    END DO
    DO lmP=1,map%lmP_max
       l=map%lmP2l(lmP)
       m=lmP2m(lmP)
       IF ( l > 0 .AND. l > m ) THEN
          map%lmP2lmPS(lmP)=map%lmP2(l-1,m)
       ELSE
          map%lmP2lmPS(lmP)=-1
       END IF
       IF ( l < map%l_max+1 ) THEN
          map%lmP2lmPA(lmP)=map%lmP2(l+1,m)
       ELSE
          map%lmP2lmPA(lmP)=-1
       END IF
    END DO
  END SUBROUTINE get_standard_lm_blocking

  SUBROUTINE get_lorder_lm_blocking(map,minc)
    type(mappings) :: map
    INTEGER,INTENT(IN) :: minc
    
    ! Local variables
    INTEGER :: m,l,lm,lmP,mc
    INTEGER,dimension(lmP_max) :: lmP2m
    
    DO m=0,map%l_max
       DO l=m,map%l_max
          map%lm2(l,m)  =-1
          map%lmP2(l,m) =-1
          !check(l,m)=0
       END DO
       l=map%l_max+1
       map%lmP2(l,m)=-1
    END DO

    lm =0
    lmP=0
    DO l=0,map%l_max
       mc =0
       ! set l2lmAS for m==0
       map%l2lmAS(l)=lm
       DO m=0,l,minc
          mc=mc+1

          lm         =lm+1
          map%lm2l(lm)   =l
          map%lm2m(lm)   =m
          map%lm2mc(lm)  =mc
          map%lm2(l,m)   =lm

          lmP        =lmP+1
          map%lmP2l(lmP) = l
          lmP2m(lmP) = m
          map%lmP2(l,m)  =lmP
          !IF ( m == 0 ) l2lmPAS(l)=lmP
          map%lm2lmP(lm) =lmP
          map%lmP2lm(lmP)=lm
       END DO
    END DO
    l=map%l_max+1    ! Extra l for lmP
    mc =0
    DO m=0,map%l_max,minc
       mc=mc+1

       lmP=lmP+1
       map%lmP2l(lmP) =l
       lmP2m(lmP) = m
       map%lmP2(l,m)  =lmP
       map%lmP2lm(lmP)=-1
    END DO

    IF ( lm /= map%lm_max ) THEN
       WRITE(*,"(2(A,I6))") 'get_lorder_lm_blocking: Wrong lm = ',lm," != map%lm_max = ",map%lm_max
       STOP
    END IF
    IF ( lmP /= map%lmP_max ) THEN
       WRITE(*,*) 'Wrong lmP!'
       STOP
    END IF
    DO lm=1,map%lm_max
       l=map%lm2l(lm)
       m=map%lm2m(lm)
       IF ( l > 0 .AND. l > m ) THEN
          map%lm2lmS(lm)=map%lm2(l-1,m)
       ELSE
          map%lm2lmS(lm)=-1
       END IF
       IF ( l < map%l_max ) THEN
          map%lm2lmA(lm)=map%lm2(l+1,m)
       ELSE
          map%lm2lmA(lm)=-1
       END IF
    END DO
    DO lmP=1,map%lmP_max
       l=map%lmP2l(lmP)
       m=lmP2m(lmP)
       IF ( l > 0 .AND. l > m ) THEN
          map%lmP2lmPS(lmP)=map%lmP2(l-1,m)
       ELSE
          map%lmP2lmPS(lmP)=-1
       END IF
       IF ( l < map%l_max+1 ) THEN
          map%lmP2lmPA(lmP)=map%lmP2(l+1,m)
       ELSE
          map%lmP2lmPA(lmP)=-1
       END IF
    END DO
  END SUBROUTINE get_lorder_lm_blocking

  SUBROUTINE get_snake_lm_blocking(map,minc)
    type(mappings) :: map
    INTEGER,INTENT(IN) :: minc
    
    ! Local variables
    INTEGER :: l,proc,lm,m,i_l,lmP,mc
    logical :: Ascending
    INTEGER,DIMENSION(map%lmP_max) :: lmP2m
    INTEGER,DIMENSION(0:n_procs-1,map%l_max+1) :: l_list
    INTEGER,DIMENSION(0:n_procs-1) :: l_counter
    INTEGER :: temp_l_counter,l0proc,pc,src_proc,temp
    INTEGER,DIMENSION(map%l_max+1) :: temp_l_list

    LOGICAL,PARAMETER :: DEBUG_OUTPUT=.false.
    ! First we loop over all l values and jump for each
    ! new l value to the next process in a snake like fashion.
    proc=0
    Ascending=.true.
    l_counter=1
    DO l=map%l_max,0,-1
       ! this l block is distributed to the actual proc
       l_list(proc,l_counter(proc))=l
       !WRITE(*,"(A,3I3)") "l,l_list,l_counter=",l,l_list(proc,l_counter(proc)),l_counter(proc)
       l_counter(proc) = l_counter(proc)+1
       if (l.eq.0) l0proc=proc
       ! now determine on which proc to put the next l value
       IF (Ascending) THEN
          IF (proc.LT.n_procs-1) THEN
             proc=proc+1
          ELSE IF (proc.EQ.n_procs-1) THEN
             Ascending=.FALSE.
          END IF
       ELSE
          IF (proc.GT.0) THEN
             proc=proc-1
          ELSE IF (proc.EQ.0) THEN
             Ascending=.True.
          END IF
       END IF
    END DO

    IF (DEBUG_OUTPUT) THEN
       DO proc=0,n_procs-1
          IF (proc.EQ.l0proc) THEN
             WRITE(*,"(A,I4,A)") "==== proc ",proc," has l=0 ===="
          ELSE
             WRITE(*,"(A,I4,A)") "---- proc ",proc," ----"
          END IF
          DO i_l=1,l_counter(proc)-1
             WRITE(*,"(I4,2X)",advance="no") l_list(proc,i_l)
          END DO
          WRITE(*,"(A)") ""
       END DO
    END IF

    ! Now distribution is as equal as possible. We rotate the distribution
    ! now to have the l0proc as first process. 
    IF (l0proc.NE.0) THEN
       temp_l_list=l_list(0,:)
       temp_l_counter=l_counter(0)
       pc = 0
       do while (.true.)
          src_proc=MODULO(l0proc+pc,n_procs)
          IF (src_proc.NE.0) THEN
             l_list(pc,:)=l_list(src_proc,:)
             l_counter(pc)=l_counter(src_proc)
          ELSE
             l_list(pc,:)=temp_l_list
             l_counter(pc)=temp_l_counter
             EXIT
          END IF
          ! now we can overwrite src_proc
          pc=src_proc
       END DO
       l0proc=0
    END IF

    ! Last step in preparation is to put the l=0 on process 0
    ! as the first l in the list
    DO i_l=1,l_counter(0)-1
       IF (l_list(0,i_l).EQ.0) THEN
          !WRITE(*,"(A,I3)") "i_l = ",i_l
          temp=l_list(0,1)
          l_list(0,1)=l_list(0,i_l)
          l_list(0,i_l)=temp
          EXIT
       END IF
    END DO

    IF (DEBUG_OUTPUT) THEN
       WRITE(*,"(A)") "Ordering after the l0proc reordering:"
       DO proc=0,n_procs-1
          IF (proc.EQ.l0proc) THEN
             WRITE(*,"(A,I4,A)") "==== proc ",proc," has l=0 ===="
          ELSE
             WRITE(*,"(A,I4,A)") "---- proc ",proc," ----"
          END IF
          DO i_l=1,l_counter(proc)-1
             WRITE(*,"(I4,2X)",advance="no") l_list(proc,i_l)
          END DO
          WRITE(*,"(A)") ""
       END DO
    END IF

    lm=1
    lmP=1
    DO proc=0,n_procs-1
       lmStartB(proc+1)=lm
       DO i_l=1,l_counter(proc)-1
          l=l_list(proc,i_l)
          mc = 0
          !WRITE(*,"(3I3)") i_l,proc,l
          DO m=0,l,minc
             mc = mc+1
             map%lm2(l,m)=lm
             map%lm2l(lm)=l
             map%lm2m(lm)=m
             map%lm2mc(lm)=mc

             map%lmP2(l,m)=lmP
             map%lmP2l(lmP)=l
             lmP2m(lmP) = m
             map%lm2lmP(lm)=lmP
             map%lmP2lm(lmP)=lm

             lm = lm+1
             lmP = lmP+1
          END DO
       END DO
       lmStopB(proc+1)=lm-1
       WRITE(*,"(I3,2(A,I6))") proc,": lmStartB=",lmStartB(proc+1),", lmStopB=",lmStopB(proc+1)
    END DO
    
    IF ( lm-1 .NE. map%lm_max ) THEN
       WRITE(*,"(2(A,I6))") 'get_snake_lm_blocking: Wrong lm-1 = ',lm-1,&
            & " != map%lm_max = ",map%lm_max
       STOP
    END IF

    l=map%l_max+1    ! Extra l for lmP
    mc =0
    DO m=0,map%l_max,minc
       mc=mc+1

       map%lmP2l(lmP) =l
       lmP2m(lmP) = m
       map%lmP2(l,m)  =lmP
       map%lmP2lm(lmP)=-1
       lmP=lmP+1
    END DO

    IF ( lmP-1 .ne. map%lmP_max ) THEN
       WRITE(*,*) 'Wrong lmP!'
       STOP
    END IF

    DO lm=1,map%lm_max
       l=map%lm2l(lm)
       m=map%lm2m(lm)
       IF ( l > 0 .AND. l > m ) THEN
          map%lm2lmS(lm)=map%lm2(l-1,m)
       ELSE
          map%lm2lmS(lm)=-1
       END IF
       IF ( l < map%l_max ) THEN
          map%lm2lmA(lm)=map%lm2(l+1,m)
       ELSE
          map%lm2lmA(lm)=-1
       END IF
    END DO
    DO lmP=1,map%lmP_max
       l=map%lmP2l(lmP)
       m=lmP2m(lmP)
       IF ( l > 0 .AND. l > m ) THEN
          map%lmP2lmPS(lmP)=map%lmP2(l-1,m)
       ELSE
          map%lmP2lmPS(lmP)=-1
       END IF
       IF ( l < map%l_max+1 ) THEN
          map%lmP2lmPA(lmP)=map%lmP2(l+1,m)
       ELSE
          map%lmP2lmPA(lmP)=-1
       END IF
    END DO
  END SUBROUTINE get_snake_lm_blocking

END MODULE blocking
