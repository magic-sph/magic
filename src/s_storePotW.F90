!$Id$
!********************************************************************
    SUBROUTINE storePotW(time,b,aj,b_ic,aj_ic,workA,workB,workC, &
                         nPotSets,root,omega_ma,omega_ic)
!********************************************************************

!--------------------------------------------------------------------

!--------------------------------------------------------------------

    USE truncation
    USE radial_functions
    USE physical_parameters
    USE horizontal_data
    USE logic
    USE output_data
    USE charmanip, only: length_to_char
    USE LMLoop_data,ONLY: llm,ulm,llm_real,ulm_real
    USE parallel_mod,only: rank
    USE communications, only: gather_from_lo_to_rank0
    IMPLICIT NONE

    REAL(kind=8),INTENT(IN) :: time

!-- Input of fields to be stored:
    COMPLEX(kind=8),INTENT(IN) :: b(llm:ulm,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: aj(llm:ulm,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: b_ic(llm:ulm,n_r_ic_max)
    COMPLEX(kind=8),INTENT(IN) :: aj_ic(llm:ulm,n_r_ic_max)
    INTEGER,INTENT(INOUT) :: nPotSets
    CHARACTER(len=9),INTENT(IN) :: root
    REAL(kind=8),INTENT(IN) :: omega_ma,omega_ic

!-- Input work arrays:
    COMPLEX(kind=8),INTENT(OUT) :: workA(llm:ulm,n_r_max)
    COMPLEX(kind=8),INTENT(OUT) :: workB(llm:ulm,n_r_max)
    COMPLEX(kind=8),INTENT(OUT) :: workC(llm:ulm,n_r_max)

    COMPLEX(kind=8),DIMENSION(:,:),ALLOCATABLE :: workA_global,workB_global
    CHARACTER(len=80) :: string
    INTEGER :: n_r,lm,n_cheb
    INTEGER :: lengthR
    CHARACTER(len=80) :: fileName
    LOGICAL :: lVB
    REAL(kind=8) :: chebNorm
     
!-- end of declaration
!---------------------------------------------------------------------
            
    nPotSets=nPotSets+1
    lengthR=length_to_char(root,'.')
    lVB=.FALSE.
    IF ( root(1:1) /= 'T' ) lVB= .TRUE. 

!--- Copy:
    DO n_r=1,n_r_max
        DO lm =llm,ulm
            workA(lm,n_r)= b(lm,n_r)
            IF ( lVB ) workB(lm,n_r)=aj(lm,n_r)
        END DO
    END DO

!--- Transform to Cheb-space:
    CALL costf1(workA,ulm_real-llm_real+1,1,ulm_real-llm_real+1, &
                workC,i_costf_init,d_costf_init)
    IF ( lVB ) CALL costf1(workB,ulm_real-llm_real+1,1,ulm_real-llm_real+1, &
         &                 workC,i_costf_init,d_costf_init)

!--- Correct amplitude:
    chebNorm=DSQRT(2.D0/(n_r_max-1))
    DO n_cheb=1,n_cheb_max
        DO lm=llm,ulm
            IF ( n_cheb == 1 .OR. n_cheb == n_r_max ) THEN
                workA(lm,n_cheb)=chebNorm*0.5D0*workA(lm,n_cheb)
                IF ( lVB) workB(lm,n_cheb)=chebNorm*0.5D0*workB(lm,n_cheb)
            ELSE
                workA(lm,n_cheb)=chebNorm*workA(lm,n_cheb)
                IF ( lVB) workB(lm,n_cheb)=chebNorm*workB(lm,n_cheb)
            END IF
        END DO
    END DO

    ! now gather the fields on rank 0 and write them to file
    ! it would be nicer to write the fields with MPI IO in parallel
    ! but then presumably the file format will change
    IF (rank.EQ.0) THEN
       ALLOCATE(workA_global(1:lm_max,1:n_cheb_max))
       ALLOCATE(workB_global(1:lm_max,1:n_cheb_max))
#ifdef WITH_DEBUG
    ELSE
       ALLOCATE(workA_global(1,1:n_cheb_max))
       ALLOCATE(workB_global(1,1:n_cheb_max))
#endif
    END IF

    DO n_cheb=1,n_cheb_max
       CALL gather_from_lo_to_rank0(workA(llm,n_cheb),workA_global(1,n_cheb))
       CALL gather_from_lo_to_rank0(workB(llm,n_cheb),workB_global(1,n_cheb))
    END DO


    IF (rank.EQ.0) THEN
       !--- Write:
       IF ( nPotSets == 0 ) THEN ! nPotSets=-1 on call
          fileName=root(1:lengthR)//'.'//tag
       ELSE
          !------ Names including the time:
          !           IF ( l_graph_time ) THEN
          !              CALL dble2string(time,'_',6,string,length)
          !              fileName=root(1:lengthR)//'_t='//string(1:length)//'.'//tag
          !           ELSE
          !------ Numbered names:
          write(string, *) nPotSets
          fileName=root(1:lengthR)//'_'//trim(adjustl(string))//'.'//tag
          !         END IF
       END IF

       OPEN(99,FILE=fileName,FORM='UNFORMATTED',STATUS='UNKNOWN')

       WRITE(99) l_max,n_cheb_max,n_cheb_ic_max,minc,lm_max
       WRITE(99) SNGL(ra),SNGL(ek),SNGL(pr),SNGL(prmag), &
            SNGL(radratio),SNGL(sigma_ratio), &
            SNGL(omega_ma),SNGL(omega_ic)
       WRITE(99) SNGL(time), &
            ((CMPLX(REAL(workA_global(lm,n_cheb)),AIMAG(workA_global(lm,n_cheb)),KIND=KIND(0d0)), &
            lm=1,lm_max),n_cheb=1,n_cheb_max)
       IF ( lVB ) &
            WRITE(99) SNGL(time), &
            ((CMPLX(REAL(workB_global(lm,n_cheb)),AIMAG(workB_global(lm,n_cheb)),KIND=KIND(0d0)), &
            lm=1,lm_max),n_cheb=1,n_cheb_max)
    END IF
    
!-- Now inner core field
    IF ( root(1:1) == 'B' .AND. l_cond_ic ) THEN

        DO n_r=1,n_r_ic_max
            DO lm =llm,ulm
                workA(lm,n_r)= b_ic(lm,n_r)
                workB(lm,n_r)=aj_ic(lm,n_r)
            END DO
        END DO

        CALL costf1(workA,ulm_real-llm_real+1,1,ulm_real-llm_real+1, &
                    workC,i_costf1_ic_init,d_costf1_ic_init)
        CALL costf1(workB,ulm_real-llm_real+1,1,ulm_real-llm_real+1, &
                    workC,i_costf1_ic_init,d_costf1_ic_init)

        chebNorm=DSQRT(2.D0/(n_r_ic_max-1))
        DO n_cheb=1,n_cheb_ic_max
           DO lm=llm,ulm
                IF ( n_cheb == 1 .OR. n_cheb == n_r_ic_max ) THEN
                    workA(lm,n_cheb)=chebNorm*0.5D0*workA(lm,n_cheb)
                    workB(lm,n_cheb)=chebNorm*0.5D0*workB(lm,n_cheb)
                ELSE
                    workA(lm,n_cheb)=chebNorm*workA(lm,n_cheb)
                    workB(lm,n_cheb)=chebNorm*workB(lm,n_cheb)
                END IF
            END DO
        END DO

        DO n_cheb=1,n_cheb_ic_max
           CALL gather_from_lo_to_rank0(workA(llm,n_cheb),workA_global(1,n_cheb))
           CALL gather_from_lo_to_rank0(workB(llm,n_cheb),workB_global(1,n_cheb))
        END DO

        IF (rank.EQ.0) THEN
           WRITE(*,*) 'WRITING IC DATA INTO FILE:',fileName

           WRITE(99) SNGL(time), &
                ((CMPLX(REAL(workA_global(lm,n_cheb)),AIMAG(workA_global(lm,n_cheb)),KIND=KIND(0d0)), &
                lm=1,lm_max),n_cheb=1,n_cheb_ic_max)
           WRITE(99) SNGL(time), &
                ((CMPLX(REAL(workB_global(lm,n_cheb)),AIMAG(workB_global(lm,n_cheb)),KIND=KIND(0d0)), &
                lm=1,lm_max),n_cheb=1,n_cheb_ic_max)
        END IF

    END IF

    CLOSE(99)

    if (rank.eq.0) DEALLOCATE(workA_global,workB_global)

    RETURN
    end SUBROUTINE storePotW

!----------------------------------------------------------------------
