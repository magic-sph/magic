!$Id$
!***********************************************************************
SUBROUTINE mapDataIC(n_r_ic_max_old,l_max_old,minc_old,         &
     &                              b_ic,dbdt_ic,aj_ic,djdt_ic)
  !***********************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/17/02  by JW. -----------

  !  +-------------+----------------+------------------------------------+
  !  |                                                                   |
  !  |   map data from input file with different grid structure          |
  !  |   or different longitudinal symmetry to the present grid          |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+
  !  |  ruler                                                            |
  !  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
  !--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+

  USE truncation
  USE radial_functions
  USE init_fields
  USE blocking

  IMPLICIT NONE

  INTEGER :: n_r_ic_max_old ! number of radial levels of old grid
  INTEGER :: l_max_old      ! max. harmonic degree for old grid
  INTEGER :: minc_old       ! longitudinal symmetry of old grid

  !-- output:
  COMPLEX(kind=8) :: b_ic(lm_maxMag,n_r_ic_maxMag)
  COMPLEX(kind=8) :: aj_ic(lm_maxMag,n_r_ic_maxMag)
  COMPLEX(kind=8) :: dbdt_ic(lm_maxMag,n_r_ic_maxMag)
  COMPLEX(kind=8) :: djdt_ic(lm_maxMag,n_r_ic_maxMag)

  !-- local:
  INTEGER :: n_data,n_dataL,n_r_ic_maxL,n_map_fac
  INTEGER :: m_max_old,lm_max_old,n_data_old,i
  INTEGER :: nR,l,m,lm,lmo,n,lmStart,lmStop,nLMB
  INTEGER :: lo,mo
  INTEGER,ALLOCATABLE :: lm2lmo(:)

  COMPLEX(kind=8),ALLOCATABLE :: b_ic_old(:),aj_ic_old(:)
  COMPLEX(kind=8),ALLOCATABLE :: dbdt_ic_old(:),djdt_ic_old(:)
  COMPLEX(kind=8),ALLOCATABLE :: b_ic_oldR(:),aj_ic_oldR(:)
  COMPLEX(kind=8),ALLOCATABLE :: dbdt_ic_oldR(:),djdt_ic_oldR(:)

  !-- end of declaration
  !-----------------------------------------------------------------------

  n_data=lm_max*n_r_ic_max
  !-- This allows to increase the number of grip points by 10!
  n_r_ic_maxL=10*n_r_ic_max
  ALLOCATE( lm2lmo(lm_max) )

  IF ( l_max.EQ.l_max_old .AND.                                   &
       &       minc.EQ.minc_old .AND.                                     &
       &       n_r_ic_max.EQ.n_r_ic_max_old) THEN

     n_data_old=n_data
     lm_max_old=lm_max
     n_map_fac=1
  ELSE 

     m_max_old =(l_max_old/minc_old)*minc_old
     lm_max_old=m_max_old*(l_max_old+1)/minc_old -                &
          &             m_max_old*(m_max_old-minc_old)/(2*minc_old) +        &
          &             l_max_old-m_max_old+1

     n_data_old=lm_max_old*n_r_ic_max_old

     !-- Write info to STDOUT:
     IF ( n_r_ic_max_old.NE.n_r_ic_max )                          &
          &        WRITE(*,'('' ! Old/New n_r_ic_max='',2i4)')               &
          &              n_r_ic_max_old,n_r_ic_max

     n_map_fac=INT(n_data_old/DBLE(n_data))+1
  END IF

  ! allocation of local arrays
  n_dataL=n_map_fac*lm_max*(n_r_ic_max+1)
  ALLOCATE( b_ic_old(n_dataL),aj_ic_old(n_dataL) )
  ALLOCATE( dbdt_ic_old(n_dataL),djdt_ic_old(n_dataL) )
  ! end of allocation of local arrays

  DO lm=1,lm_max
     l=lm2l(lm)
     m=lm2m(lm)
     lm2lmo(lm)=-1 ! -1 means that there is no data in the startfile
     lmo=0
     DO mo=0,l_max_old,minc_old
        DO lo=mo,l_max_old
           lmo=lmo+1
           IF ( lo.EQ.l .AND. mo.EQ.m ) THEN
              lm2lmo(lm)=lmo ! data found in startfile
              GOTO 10
           END IF
        END DO
     END DO
10   CONTINUE
  END DO

  READ(n_start_file) (b_ic_old(i),i=1,n_data_old),                &
       &                     (aj_ic_old(i),i=1,n_data_old),               &
       &                     (dbdt_ic_old(i),i=1,n_data_old),             &
       &                     (djdt_ic_old(i),i=1,n_data_old)

  ALLOCATE( b_ic_oldR(n_r_ic_maxL),aj_ic_oldR(n_r_ic_maxL) )
  ALLOCATE( dbdt_ic_oldR(n_r_ic_maxL),djdt_ic_oldR(n_r_ic_maxL) )

  DO nLMB=1,nLMBs ! Blocking of loop over all (l,m)
     lmStart=lmStartB(nLMB)
     lmStop =lmStopB(nLMB)
     lmStart=MAX(2,lmStart)
     DO lm=lmStart,lmStop
        lmo=lm2lmo(lm)
        IF ( lmo.GT.0 ) THEN
           IF ( n_r_ic_max.NE.n_r_ic_max_old ) THEN
              DO nR=1,n_r_ic_max_old  ! copy on help arrays
                 n=lmo+(nR-1)*lm_max_old
                 b_ic_oldR(nR)   =b_ic_old(n)
                 aj_ic_oldR(nR)  =aj_ic_old(n)
                 dbdt_ic_oldR(nR)=dbdt_ic_old(n)
                 djdt_ic_oldR(nR)=djdt_ic_old(n)
              END DO
              CALL mapDataICR(b_ic_oldR,n_r_ic_max_old)
              CALL mapDataICR(aj_ic_oldR,n_r_ic_max_old)
              CALL mapDataICR(dbdt_ic_oldR,n_r_ic_max_old)
              CALL mapDataICR(djdt_ic_oldR,n_r_ic_max_old)
              DO nR=1,n_r_ic_max
                 b_ic(lm,nR)   =b_ic_oldR(nR)
                 aj_ic(lm,nR)  =aj_ic_oldR(nR)
                 dbdt_ic(lm,nR)=dbdt_ic_oldR(nR)
                 djdt_ic(lm,nR)=djdt_ic_oldR(nR)
              END DO
           ELSE
              DO nR=1,n_r_ic_max
                 n=lmo+(nR-1)*lm_max_old
                 b_ic(lm,nR)   =b_ic_old(n)
                 aj_ic(lm,nR)  =aj_ic_old(n)
                 dbdt_ic(lm,nR)=dbdt_ic_old(n)
                 djdt_ic(lm,nR)=djdt_ic_old(n)
              END DO
           END IF
        END IF
     END DO
  END DO
  DEALLOCATE( b_ic_oldR,aj_ic_oldR )
  DEALLOCATE( dbdt_ic_oldR,djdt_ic_oldR )
  DEALLOCATE( lm2lmo )
  DEALLOCATE( b_ic_old,aj_ic_old )
  DEALLOCATE( dbdt_ic_old,djdt_ic_old )

END SUBROUTINE mapDataIC

!-------------------------------------------------------------------------------
