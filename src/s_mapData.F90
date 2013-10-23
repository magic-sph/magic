!$Id$
!***************************************************************************
SUBROUTINE mapData(n_r_max_old,l_max_old,minc_old,l_mag_old,      &
     &                   w,dwdt,z,dzdt,p,dpdt,s,dsdt,b,dbdt,aj,djdt)  
  !***************************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/17/02  by JW. -----------

  !--------------------------------------------------------------------------

  USE truncation
  USE radial_functions
  USE init_fields
  USE blocking
  USE logic
  use omp_lib

  IMPLICIT NONE

  !-- input:
  INTEGER :: n_r_max_old
  INTEGER :: l_max_old
  INTEGER :: minc_old
  LOGICAL :: l_mag_old

  !-- output:
  COMPLEX(kind=8) :: w(lm_max,n_r_max),dwdt(lm_max,n_r_max)
  COMPLEX(kind=8) :: z(lm_max,n_r_max),dzdt(lm_max,n_r_max)
  COMPLEX(kind=8) :: p(lm_max,n_r_max),dpdt(lm_max,n_r_max)
  COMPLEX(kind=8) :: s(lm_max,n_r_max),dsdt(lm_max,n_r_max)
  COMPLEX(kind=8) :: b(lm_maxMag,n_r_maxMag),dbdt(lm_maxMag,n_r_maxMag)
  COMPLEX(kind=8) :: aj(lm_maxMag,n_r_maxMag),djdt(lm_maxMag,n_r_maxMag)

  !-- local:
  INTEGER :: n_data,n_dataL,n_r_maxL,n_map_fac
  INTEGER :: n_data_old,n_data_oldP
  INTEGER :: i,m_max_old,lm_max_old,nR
  INTEGER :: lmStart,lmStop,lm,lmo,l,m,nLMB
  INTEGER :: lo,mo,n
  INTEGER,ALLOCATABLE :: lm2lmo(:)

  REAL(kind=8) :: pi,fr
  LOGICAL :: lreadS

  COMPLEX(kind=8),ALLOCATABLE :: wo(:),zo(:),po(:),so(:)
  COMPLEX(kind=8),ALLOCATABLE :: woR(:),zoR(:)
  COMPLEX(kind=8),ALLOCATABLE :: poR(:),soR(:)

  !-- end of declaration
  !-----------------------------------------------------------------------

  n_data=lm_max*n_r_max
  n_r_maxL=10*n_r_max
  !-- This allows to increase the number of grid points by 10!
  ALLOCATE( lm2lmo(lm_max) )

  IF ( inform.EQ.6 .OR. inform.EQ.7 .OR. inform.EQ.9 .OR.         &
       &       inform.EQ.11 ) THEN
     lreadS=.FALSE.
  ELSE
     lreadS=.TRUE.
  END IF

  IF ( l_max.EQ.l_max_old .AND.                                   &
       &       minc.EQ.minc_old .AND.                                     &
       &       n_r_max.EQ.n_r_max_old ) THEN

     !----- Direct reading of fields, grid not changed:
     WRITE(*,'(/,'' ! Reading fields directly.'')')

     n_data_old=n_data
     IF ( inform.GT.2 ) THEN
        n_data_oldP=n_data
     ELSE
        !----- In the past an 'extra' radial grid point has been
        !      stored which was not really necessary
        n_data_oldP=lm_max*(n_r_max+1)
     END IF

     lm_max_old=lm_max
     n_map_fac=1

  ELSE 

     !----- Mapping onto new grid !
     WRITE(*,'(/,'' ! Mapping onto new grid.'')')

     IF ( MOD(minc_old,minc).NE.0)                                &
          &     WRITE(6,'('' ! Warning: Incompatible old/new minc= '',2i3)')

     m_max_old =(l_max_old/minc_old)*minc_old
     lm_max_old=m_max_old*(l_max_old+1)/minc_old -                &
          &             m_max_old*(m_max_old-minc_old)/(2*minc_old) +        &
          &             l_max_old-m_max_old+1

     n_data_old=lm_max_old*n_r_max_old
     IF ( inform.GT.2 ) THEN
        n_data_oldP=n_data_old
     ELSE
        n_data_oldP=lm_max_old*(n_r_max_old+1)
     END IF

     !-- Write info to STDOUT:
     WRITE(*,'('' ! Old/New  l_max= '',2I4,''  m_max= '',2I4,     &
          &       ''  minc= '',2I3,''  lm_max= '',2I5/)')                    &
          &           l_max_old,l_max,m_max_old,m_max,                       &
          &           minc_old,minc,lm_max_old,lm_max
     IF ( n_r_max_old.NE.n_r_max )                                &
          &        WRITE(*,'('' ! Old/New n_r_max='',2i4)')                  &
          &              n_r_max_old,n_r_max

     n_map_fac=INT(n_data_old/DBLE(n_data))+1
  END IF

  ! allocation of local arrays.
  ! if this becomes a performance bottleneck, one can make a module
  ! and allocate the array only once in the initialization
  n_dataL=n_map_fac*lm_max*(n_r_max+1)
  ALLOCATE( wo(n_dataL),zo(n_dataL),po(n_dataL),so(n_dataL) )
  ! end of allocation

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
              CYCLE
           END IF
        END DO
     END DO
  END DO

  IF ( lreadS ) THEN
     READ(n_start_file) (wo(i),i=1,n_data_oldP),                  &
          &                        (zo(i),i=1,n_data_oldP),                  &
          &                        (po(i),i=1,n_data_oldP),                  &
          &                        (so(i),i=1,n_data_oldP)
  ELSE
     READ(n_start_file) (wo(i),i=1,n_data_oldP),                  &
          &                        (zo(i),i=1,n_data_oldP),                  &
          &                        (po(i),i=1,n_data_oldP)
  END If

  !-- Select only the spherical harmonic modes you need:

  !PRINT*,omp_get_thread_num(),": Before nLMB loop, nLMBs=",nLMBs
  ALLOCATE( woR(n_r_maxL),zoR(n_r_maxL) )
  ALLOCATE( poR(n_r_maxL),soR(n_r_maxL) )
  
  DO nLMB=1,nLMBs ! Blocking of loop over all (l,m)
     lmStart=lmStartB(nLMB)
     lmStop =lmStopB(nLMB)

     !PRINT*,nLMB,lmStart,lmStop

     DO lm=lmStart,lmStop
        lmo=lm2lmo(lm) 
        IF ( lmo.GT.0 ) THEN
           IF ( n_r_max.NE.n_r_max_old ) THEN
              DO nR=1,n_r_max_old  ! copy on help arrays
                 n=lmo+(nR-1)*lm_max_old
                 woR(nR)=wo(n)
                 zoR(nR)=zo(n)
                 poR(nR)=po(n)
                 IF ( l_heat ) soR(nR)=so(n)
              END DO
              CALL mapDataR(woR,n_r_max_old,.FALSE.)
              CALL mapDataR(zoR,n_r_max_old,.FALSE.)
              CALL mapDataR(poR,n_r_max_old,.FALSE.)
              IF ( l_heat ) CALL mapDataR(soR,n_r_max_old,.FALSE.)
              DO nR=1,n_r_max
                 IF ( lm.GT.1 ) THEN
                    w(lm,nR)=woR(nR)
                    z(lm,nR)=zoR(nR)
                    p(lm,nR)=poR(nR)
                 END IF
                 IF ( l_heat ) s(lm,nR)=soR(nR)
              END DO
           ELSE  
              DO nR=1,n_r_max
                 n=lmo+(nR-1)*lm_max_old
                 IF ( lm.GT.1 ) THEN
                    w(lm,nR)=wo(n)
                    z(lm,nR)=zo(n)
                    p(lm,nR)=po(n)
                 END IF
                 IF ( l_heat ) s(lm,nR)=so(n)
              END DO
           END IF
        END IF
     END DO
  END DO
  !PRINT*,omp_get_thread_num(),": After nLMB loop"
  DEALLOCATE(woR,zoR,poR,soR)

  !---- Read time derivatives:
  IF ( lreadS ) THEN
     READ(n_start_file) (so(i),i=1,n_data_old),                   &
          &                        (wo(i),i=1,n_data_old),                   &
          &                        (zo(i),i=1,n_data_old),                   &
          &                        (po(i),i=1,n_data_old)
  ELSE
     READ(n_start_file) (wo(i),i=1,n_data_old),                   &
          &                        (zo(i),i=1,n_data_old),                   &
          &                        (po(i),i=1,n_data_old)
  END IF


  ALLOCATE( woR(n_r_maxL),zoR(n_r_maxL) )
  ALLOCATE( poR(n_r_maxL),soR(n_r_maxL) )

  DO nLMB=1,nLMBs ! Blocking of loop over all (l,m)
     lmStart=lmStartB(nLMB)
     lmStop =lmStopB(nLMB)
     DO lm=lmStart,lmStop
        lmo=lm2lmo(lm) ! Where is this in the input file?
        IF ( lmo.GT.0 ) THEN
           IF ( n_r_max.NE.n_r_max_old ) THEN
              DO nR=1,n_r_max_old
                 n=lmo+(nR-1)*lm_max_old
                 IF ( l_heat ) soR(nR)=so(n)
                 woR(nR)=wo(n)
                 zoR(nR)=zo(n)
                 poR(nR)=po(n)
              END DO
              IF ( l_heat ) CALL mapDataR(soR,n_r_max_old,.TRUE.)
              CALL mapDataR(woR,n_r_max_old,.TRUE.)
              CALL mapDataR(zoR,n_r_max_old,.TRUE.)
              CALL mapDataR(poR,n_r_max_old,.TRUE.)
              DO nR=1,n_r_max
                 IF ( l_heat ) dsdt(lm,nR)=soR(nR)
                 IF ( lm.GT.1 ) THEN
                    dwdt(lm,nR)=woR(nR)
                    dzdt(lm,nR)=zoR(nR)
                    dpdt(lm,nR)=poR(nR)
                 END IF
              END DO
           ELSE
              DO nR=1,n_r_max
                 n=lmo+(nR-1)*lm_max_old
                 IF ( l_heat ) dsdt(lm,nR)=so(n)
                 IF ( lm.GT.1 ) THEN
                    dwdt(lm,nR)=wo(n)
                    dzdt(lm,nR)=zo(n)
                    dpdt(lm,nR)=po(n)
                 END IF
              END DO
           END IF
        END IF
     END DO
  END DO
  DEALLOCATE( woR,zoR,poR,soR )


  !--- Read magnetic field
  IF ( l_mag_old ) THEN

     READ(n_start_file) (so(i),i=1,n_data_oldP),               &
             &                           (wo(i),i=1,n_data_oldP),               &
             &                           (zo(i),i=1,n_data_old),                &
             &                           (po(i),i=1,n_data_old)

     ALLOCATE( woR(n_r_maxL),zoR(n_r_maxL) )
     ALLOCATE( poR(n_r_maxL),soR(n_r_maxL) )

     DO nLMB=1,nLMBs ! Blocking of loop over all (l,m)
        lmStart=lmStartB(nLMB)
        lmStop =lmStopB(nLMB)
        lmStart=MAX(2,lmStart)
        DO lm=lmStart,lmStop
           lmo=lm2lmo(lm) ! Where is this in the input file?
           IF ( lmo.GT.0 ) THEN
              IF ( n_r_max.NE.n_r_max_old ) THEN
                 DO nR=1,n_r_max_old
                    n=lmo+(nR-1)*lm_max_old
                    soR(nR)=so(n)
                    woR(nR)=wo(n)
                    zoR(nR)=zo(n)
                    poR(nR)=po(n)
                 END DO
                 CALL mapDataR(soR,n_r_max_old,.FALSE.)
                 CALL mapDataR(woR,n_r_max_old,.FALSE.)
                 CALL mapDataR(zoR,n_r_max_old,.TRUE.)
                 CALL mapDataR(poR,n_r_max_old,.TRUE.)
                 DO nR=1,n_r_max
                    b(lm,nR)   =soR(nR)
                    aj(lm,nR)  =woR(nR)
                    dbdt(lm,nR)=zoR(nR)
                    djdt(lm,nR)=poR(nR)
                 END DO
              ELSE
                 DO nR=1,n_r_max
                    n=lmo+(nR-1)*lm_max_old
                    b(lm,nR)   =so(n)
                    aj(lm,nR)  =wo(n)
                    dbdt(lm,nR)=zo(n)
                    djdt(lm,nR)=po(n)
                 END DO
              END IF
           END IF
        END DO
     END DO
     DEALLOCATE( woR,zoR,poR,soR )
  ELSE
     WRITE(*,*) '! No magnetic data in input file!'
  END IF


  !-- If mapping towards reduced symmetry, add thermal perturbation in
  !   mode (l,m)=(minc,minc) if parameter tipdipole .ne. 0
  IF ( l_heat .AND.                                               &
       &       minc.LT.minc_old .AND. tipdipole.GT.0.D0 ) THEN
     DO nLMB=1,nLMBs ! Blocking of loop over all (l,m)
        lmStart=lmStartB(nLMB)
        lmStop =lmStopB(nLMB)
        lm=l_max+2
        IF ( lmStart.LE.lm .AND. lmStop.GE.lm ) THEN
           pi=4.d0*DATAN(1.D0)
           DO nR=1,n_r_max+1
              fr=dsin(pi*(r(nR)-r(n_r_max)))
              s(lm,nR)=tipdipole*fr
           END DO
        END IF
     END DO
  END IF

  !-- If starting from data file with longitudinal symmetry, add
  !   weak non-axisymmetric dipole component if tipdipole .ne. 0
  IF ( ( l_mag .OR. l_mag_LF )                                    &
       &       .AND. minc.EQ.1 .AND. minc_old.NE.1 .AND.                  &
       &       tipdipole.GT.0.d0 .AND. l_mag_old ) THEN
     DO nLMB=1,nLMBs ! Blocking of loop over all (l,m)
        lmStart=lmStartB(nLMB)
        lmStop =lmStopB(nLMB)
        lm=l_max+2
        IF ( lmStart.LE.lm .AND. lmStop.GE.lm ) THEN
           DO nR=1,n_r_max+1
              b(lm,nR)=tipdipole
           END DO
        END IF
     END DO
  END IF

  ! deallocation of the local arrays
  DEALLOCATE( lm2lmo )
  DEALLOCATE( wo,zo,po,so )
  ! end of deallocation of local arrays

END SUBROUTINE mapData

!------------------------------------------------------------------------
