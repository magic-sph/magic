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
  INTEGER,ALLOCATABLE :: lm2(:,:),lm2l(:),lm2m(:)
  INTEGER,ALLOCATABLE :: lm2mc(:),l2lmAS(:)
  !INTEGER,ALLOCATABLE :: m2mc(:)
  INTEGER,ALLOCATABLE :: lm2lmS(:),lm2lmA(:)

  INTEGER,ALLOCATABLE :: lmP2(:,:),lmP2l(:)!,lmP2m(:)
  !INTEGER,ALLOCATABLE :: l2lmPAS(:)
  INTEGER,ALLOCATABLE :: lmP2lmPS(:),lmP2lmPA(:)

  INTEGER,ALLOCATABLE :: lm2lmP(:),lmP2lm(:)

  !INTEGER :: nLMBsMax
  INTEGER :: nLMBs,sizeLMB

  INTEGER,ALLOCATABLE :: lmStartB(:),lmStopB(:)
  INTEGER,ALLOCATABLE :: nLMBs2(:),sizeLMB2(:,:)
  INTEGER,PARAMETER :: sizeLMB2max=201
  INTEGER,ALLOCATABLE :: lm22lm(:,:,:)
  INTEGER,ALLOCATABLE :: lm22l(:,:,:)
  INTEGER,ALLOCATABLE :: lm22m(:,:,:)

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

    ! nthreads > 1
    ALLOCATE( lm2(0:l_max,0:l_max),lm2l(lm_max),lm2m(lm_max) )
    ALLOCATE( lm2mc(lm_max),l2lmAS(0:l_max) )
    !ALLOCATE( m2mc(0:l_max) )
    ALLOCATE( lm2lmS(lm_max),lm2lmA(lm_max) )

    ALLOCATE( lmP2(0:l_max+1,0:l_max+1),lmP2l(lmP_max) )!,lmP2m(lmP_max) )
    !ALLOCATE( l2lmPAS(0:l_max+1) )
    ALLOCATE( lmP2lmPS(lmP_max),lmP2lmPA(lmP_max) )

    ALLOCATE( lm2lmP(lm_max),lmP2lm(lmP_max) )

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
#endif
    WRITE(*,*)
    WRITE(*,*) '! Max thread number available :',nThreadsAva
    WRITE(*,*) '! Max thread number demanded  :',nThreadsRun
    WRITE(*,*) '! Number of threads I will use:',nThreads
    WRITE(*,*)

#ifdef WITH_MPI
    nLMBs=2*((lm_max-1)/(nChunk*n_procs)+1) * n_procs
#else
#ifdef WITHOMP
    nThreadsMax = omp_get_max_threads()
#else
    nThreadsMax = 1
#endif
    nLMBs=2*((lm_max-1)/(nChunk*nThreads)+1) * nThreads
#endif

    ALLOCATE( lmStartB(nLMBs),lmStopB(nLMBs) )
    ALLOCATE( nLMBs2(nLMBs),sizeLMB2(l_max+1,nLMBs) )
    ALLOCATE( lm22lm(sizeLMB2max,l_max+1,nLMBs) )
    ALLOCATE( lm22l(sizeLMB2max,l_max+1,nLMBs) )
    ALLOCATE( lm22m(sizeLMB2max,l_max+1,nLMBs) )

    nfs=(sizeThetaBI/(n_phi_tot+nBSave)+1) * nBDown

  END SUBROUTINE initialize_blocking

  !***********************************************************************
  SUBROUTINE getBlocking

    INTEGER :: lm,l,m,n,n2,n3
    INTEGER :: lmP,mc
    INTEGER :: nThreadsAva,nB2

    INTEGER :: check(0:l_max,0:l_max)
    INTEGER :: help,help1(lm_max),help2(lm_max),help3(lm_max)

    integer :: lmP2m(lmP_max)
    REAL(kind=8) :: load
    INTEGER :: iLoad

    LOGICAL :: lStop
    INTEGER :: size, necessary_nLMBs

    !--- End of declaration
    !----------------------------------------------------------------------

    IF ( l_save_out ) THEN
       OPEN(nLF,FILE=log_file,STATUS='UNKNOWN',POSITION='APPEND')
    END IF

#ifdef WITHOMP
    nThreadsAva=OMP_GET_NUM_PROCS()
#else
    nThreadsAva=1
#endif
    WRITE(nLF,*)
    WRITE(nLF,*) '! Max thread number available :',nThreadsAva
    WRITE(nLF,*) '! Max thread number demanded  :',nThreadsRun
    WRITE(nLF,*) '! Number of threads I will use:',nThreads
    WRITE(nLF,*)



    IF ( l_RMS .AND. nThreads > nThreadsMax ) THEN
       WRITE(*,*) '! Too small value of nThreadsMax !'
       WRITE(*,*) '! for calculating RMS forces!'
       WRITE(*,*) '! See c_RMS.f!'
       WRITE(*,*) '! Increase nThreadsMax in m_blocking.F90!'
       STOP
    END IF

    !--- Get radial blocking
    IF ( MOD(n_r_max-1,nThreads) /= 0 ) THEN
       WRITE(*,*) 'Number of threads has to be multiple of n_r_max-1!'
       WRITE(*,*) 'nThreads :',nThreads
       WRITE(*,*) 'n_r_max-1:',n_r_max-1
       STOP
    END IF
    sizeRB=(n_r_max-1)/nThreads

    !--- Calculate lm and ml blocking:

    DO m=0,l_max
       DO l=m,l_max
          lm2(l,m)  =-1
          lmP2(l,m) =-1
          check(l,m)=0
       END DO
       l=l_max+1
       lmP2(l,m)=-1
    END DO

    lm =0
    lmP=0
    mc =0
    DO m=0,l_max,minc
       mc=mc+1
       !m2mc(m)=mc
       DO l=m,l_max
          lm         =lm+1
          lm2l(lm)   =l
          lm2m(lm)   =m
          lm2mc(lm)  =mc
          lm2(l,m)   =lm
          IF ( m == 0 ) l2lmAS(l)=lm
          lmP        =lmP+1
          lmP2l(lmP) =l
          lmP2m(lmP) =m
          lmP2(l,m)  =lmP
          !IF ( m == 0 ) l2lmPAS(l)=lmP
          lm2lmP(lm) =lmP
          lmP2lm(lmP)=lm
       END DO
       l=l_max+1    ! Extra l for lmP
       lmP=lmP+1
       lmP2l(lmP) =l
       lmP2m(lmP) =m
       lmP2(l,m)  =lmP
       !IF ( m == 0 ) l2lmPAS(l)=lmP
       lmP2lm(lmP)=-1
    END DO
    IF ( lm /= lm_max ) THEN
       WRITE(*,*) 'Wrong lm!'
       STOP
    END IF
    IF ( lmP /= lmP_max ) THEN
       WRITE(*,*) 'Wrong lmP!'
       STOP
    END IF
    DO lm=1,lm_max
       l=lm2l(lm)
       m=lm2m(lm)
       IF ( l > 0 .AND. l > m ) THEN
          lm2lmS(lm)=lm2(l-1,m)
       ELSE
          lm2lmS(lm)=-1
       END IF
       IF ( l < l_max ) THEN
          lm2lmA(lm)=lm2(l+1,m)
       ELSE
          lm2lmA(lm)=-1
       END IF
    END DO
    DO lm=1,lmP_max
       l=lmP2l(lm)
       m=lmP2m(lm)
       IF ( l > 0 .AND. l > m ) THEN
          lmP2lmPS(lm)=lmP2(l-1,m)
       ELSE
          lmP2lmPS(lm)=-1
       END IF
       IF ( l < l_max+1 ) THEN
          lmP2lmPA(lm)=lmP2(l+1,m)
       ELSE
          lmP2lmPA(lm)=-1
       END IF
    END DO

    !--- Check whether we need less than nLMBsMax blocks:
    !#ifdef WITH_MPI
    !  nLMBs  =((lm_max-1)/(nChunk*n_procs)+1)*n_procs
    !#else
    !  nLMBs  =((lm_max-1)/(nChunk*nThreads)+1)*nThreads
    !#endif
    !       nLMBs=nThreads
    sizeLMB=(lm_max-1)/nLMBs+1
    !  IF ( nLMBs > nLMBsMax ) THEN
    !     WRITE(*,*) '! Too small nLMBsMax in m_blocking.F90!'
    !     WRITE(*,*) '! nLMBs=',nLMBs
    !     WRITE(*,*) '! nLMBsMax=',nLMBsMax
    !     WRITE(*,*) '! nThreads=',nThreads
    !     WRITE(*,*) '! nThreadsMax=',nThreadsMax
    !     WRITE(*,*) '! Increase nThreadsMax!'
    !     STOP
    !  END IF
    !PRINT*,"lm_max = ",lm_max,", nLMBs = ",nLMBs,", nThreads = ",nThreads
    necessary_nLMBs = (lm_max+(sizeLMB-1))/sizeLMB
    IF ( necessary_nLMBs .LT. nLMBs ) THEN
       WRITE(*,*) '!  Not all lm-Blocks necessary!'
       WRITE(*,*) '!  No of redundant blocks:',nLMBs-necessary_nLMBs
       nLMBs = necessary_nLMBs
    END IF
    IF ( nLMBs*sizeLMB > lm_max ) THEN
       WRITE(*,*) '! Uneven load balancing in LM blocks!'
       load=DBLE(lm_max-(nLMBs-1)*sizeLMB)/sizeLMB
       WRITE(*,*) '! Load percentage of last block:',load*100.D0
       iLoad=INT(load,4)
       IF ( iLoad >= 1 ) THEN
          WRITE(*,*) '! No. of redundant blocks:',iLoad
          nLMBs=nLMBs-iLoad
       END IF
    END IF
    !--- Get lm start and stop for main blocks:
    DO n=1,nLMBs
       lmStartB(n)=(n-1)*sizeLMB+1
       lmStopB(n) =MIN(n*sizeLMB,lm_max)
       WRITE(*,*) n,lmStartB(n),lmStopB(n),lmStopB(n)-lmStartB(n)+1
       IF ( lmStopB(n) == lm_max ) EXIT
    END DO

    !--- Getting lm sub-blocks:
    lStop=.FALSE.
    size=0
    nB2=0
    DO n=1,nLMBs
       nLMBs2(n)=1
       lm=lmStartB(n)
       !------ Start first sub-block:
       sizeLMB2(1,n) =1
       lm22lm(1,1,n) =lm
       lm22l(1,1,n)  =lm2l(lm)
       lm22m(1,1,n)  =lm2m(lm)
       IF ( lmStartB(n) < lmStopB(n) ) THEN
          DO lm=lmStartB(n)+1,lmStopB(n)
             DO n2=1,nLMBs2(n)
                IF ( lm22l(1,n2,n) == lm2l(lm) ) THEN
                   !------ Add to old block
                   sizeLMB2(n2,n)=sizeLMB2(n2,n)+1
                   GOTO 20
                END IF
             END DO
             !------ Start new l-block:
             n2            =nLMBs2(n)+1
             nLMBs2(n)     =n2
             sizeLMB2(n2,n)=1
20           CONTINUE
             lm22lm(sizeLMB2(n2,n),n2,n)=lm
             lm22l(sizeLMB2(n2,n),n2,n) =lm2l(lm)
             lm22m(sizeLMB2(n2,n),n2,n) =lm2m(lm)
          END DO
       END IF

       !------ Resort:
       IF ( nLMBs2(n) > 1 ) THEN
          n2=1
30        DO n3=n2+1,nLMBs2(n)
             IF  ( lm22m(1,n2,n) > lm22m(1,n3,n) ) THEN
                help=sizeLMB2(n2,n)
                DO lm=1,help
                   help1(lm)=lm22l(lm,n2,n)
                   help2(lm)=lm22m(lm,n2,n)
                   help3(lm)=lm22lm(lm,n2,n)
                END DO
                sizeLMB2(n2,n)=sizeLMB2(n3,n)
                DO lm=1,sizeLMB2(n2,n)
                   lm22l(lm,n2,n) =lm22l(lm,n3,n)
                   lm22m(lm,n2,n) =lm22m(lm,n3,n)
                   lm22lm(lm,n2,n)=lm22lm(lm,n3,n)
                END DO
                sizeLMB2(n3,n)=help
                DO lm=1,help
                   lm22l(lm,n3,n) =help1(lm)
                   lm22m(lm,n3,n) =help2(lm)
                   lm22lm(lm,n3,n)=help3(lm)
                END DO
                GOTO 30
             END IF
          END DO
          n2=n2+1
          IF ( n2 == nLMBs2(n) ) GOTO 40
          GOTO 30
       END IF
40     CONTINUE

       nB2=nB2+nLMBs2(n)
       DO n2=1,nLMBs2(n)
          IF ( sizeLMB2(n2,n) > sizeLMB2max ) THEN
             lStop=.TRUE.
             size=MAX(size,sizeLMB2(n2,n))
          END IF
          DO n3=1,sizeLMB2(n2,n)
             l=lm22l(n3,n2,n)
             m=lm22m(n3,n2,n)
             check(l,m)=check(l,m)+1
          END DO
          !             WRITE(99,*) n,n2,sizeLMB2(n2,n)
       END DO

    END DO
    IF ( lStop ) THEN
       WRITE(*,*) '! Increase sizeLMB2max in m_blocking.F90!'
       WRITE(*,*) '! to at least:',size
       STOP
    END IF


    DO m=0,l_max,minc
       DO l=m,l_max
          IF ( check(l,m) == 0 ) THEN
             WRITE(*,*) 'Warning, forgotten l,m:',l,m,lm2(l,m)
             STOP
          ELSE IF ( check(l,m) > 1 ) THEN
             WRITE(*,*) 'Warning, too much l,m:',l,m,check(l,m)
             STOP
          END IF
       END DO
    END DO

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
       WRITE(*,*) '! Please decrease sizeThetaBI or &
            &          nBDown in m_blocking.F90!'
       STOP
    END IF



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


    IF ( l_save_out ) CLOSE(n_log_file)

    RETURN
  end SUBROUTINE getBlocking

  !------------------------------------------------------------------------

END MODULE blocking
