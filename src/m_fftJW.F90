!$Id$
!*******************************************************************************
#include "perflib_preproc.cpp"

MODULE fft_JW
  use truncation
  use blocking
  use parallel_mod
  IMPLICIT NONE

  PRIVATE
  !-- For Fourier transform (calculated in init_fft)
  INTEGER,PARAMETER :: ni=100
  INTEGER :: nd
  INTEGER :: i_fft_init(ni)
  REAL(kind=8),ALLOCATABLE ::  d_fft_init(:)
  !$OMP THREADPRIVATE( i_fft_init,d_fft_init )

  INTERFACE fft_thetab
     module procedure fft_thetab_real
     module procedure fft_thetab_cmplx
  END INTERFACE

  PUBLIC :: fft_thetab, init_fft, fft_to_real

CONTAINS
  
  !SUBROUTINE init_fft(n,i_fft_init,ni,d_fft_init,nd)
  SUBROUTINE init_fft(n)
    !***********************************************************************
    
    !  +-------------+----------------+------------------------------------+
    !  |  Purpose of this subroutine is to calculate and store several     |
    !  |  values that will be needed for a fast fft transform.             |
    !  |  The actual transform is performed by the subroutine fftJW.       |
    !  +-------------------------------------------------------------------+

    use usefull, only: factorise
    USE const, only: pi

    IMPLICIT NONE
    !-- Input
    INTEGER :: n     ! Dimension of problem, number of grid points
    !INTEGER :: ni,nd ! Dimensions reserved for output arrays

    !-- Output
    !INTEGER :: i_fft_init(ni) ! Contrains factorization
    !REAL(kind=8) :: d_fft_init(nd)  ! Contrains sin and cos

    INTEGER :: i,j,nFacs,nFactors,help
    INTEGER, PARAMETER :: nFacsA=100
    INTEGER :: fac(nFacsA),factor(nFacsA)
    REAL(kind=8) :: phi,dPhi

    !-- end of declaration
    !-------------------------------------------------------------------
    !$OMP MASTER
    nd = 3*(n/2)
    !$OMP END MASTER
    !$OMP BARRIER

    !-- For Fourier transform (calculated in init_fft)
    ALLOCATE( d_fft_init(nd) )



    !-- Checking number of datapoints:
    IF ( n <= 3 ) then
       WRITE(*,*) '! Message from subroutine init_fft:'
       WRITE(*,*) '! Sorry, I need more than 3 grid points!'
       STOP
    END IF
    IF ( MOD(n,4) /= 0 ) THEN
       WRITE(*,*) '! Note from subroutine init_fft:'
       WRITE(*,*) '! Number of data points has to be'
       WRITE(*,*) '! a mutiple of 4!'
       STOP
    END IF

    nFacs=4  ! factors to be tried
    fac(1)=4
    fac(2)=2
    fac(3)=3
    fac(4)=5

    CALL factorise(n/2,nFacs,fac,nFactors,factor)
    IF ( nFactors > nFacsA ) then
       WRITE(*,*) '! Message from subroutine init_fft:'
       WRITE(*,*) '! Please increas nFacsA!'
       WRITE(*,*) '! Should be >= ',nFactors
       STOP
    END IF

    !-- Sort in ascending order:
    IF ( nFactors > 1 ) THEN
       DO i=2,nFactors
          help=factor(i)
          DO j=i-1,1,-1
             IF ( factor(j) <= help ) GOTO 15
             factor(j+1)=factor(j) ! SHIFT UP
          END DO
          j=0
15        factor(j+1)=help ! INSERT
       END DO
    END IF

    !-- Store:
    IF ( ni <= nFactors+1 ) THEN
       WRITE(*,*) '! Message from subroutine init_fft:'
       WRITE(*,*) '! Increase dimension of array i_fft_init'
       WRITE(*,*) '! in calling routine.'
       WRITE(*,*) '! Should be at least:',nFactors+1
       STOP
    END IF
    i_fft_init(1)=nFactors
    DO i=1,nFactors
       i_fft_init(i+1)=factor(i)
    END DO

    !-- Calculate trigonometric functions:
    j=n/2
    IF ( nd < n+j ) THEN
       WRITE(*,*) '! Message from subroutine init_fft:'
       WRITE(*,*) '! Increase dimension of array d_fft_init'
       WRITE(*,*) '! in calling routine.'
       WRITE(*,*) '! Should be at least:',n+j
       WRITE(*,*) '! But is only       :',nd
       STOP
    END IF
    !pi=4.D0*DATAN(1.D0)
    dPhi=2.D0*pi/DBLE(n)
    DO i=1,n,2
       phi=DBLE(i-1)*dPhi
       d_fft_init(i)  =DCOS(phi)
       d_fft_init(i+1)=DSIN(phi)
    END DO
    dPhi=0.5D0*dPhi
    DO i=1,j,2
       phi=DBLE(i-1)*dPhi
       d_fft_init(n+i)  =DCOS(phi)
       d_fft_init(n+i+1)=DSIN(phi)
    END DO

  end SUBROUTINE init_fft

  SUBROUTINE fft_to_real(f,ld_f,nrep)
    INTEGER,INTENT(IN) :: ld_f,nrep
    COMPLEX(kind=8), INTENT(INOUT) :: f(ld_f/2,nrep)

    REAL(kind=8) :: work(ld_f,nrep)

    PERFON('fft2r')
    CALL fftJW(f,ld_f,n_phi_max,1,nrep,work,ld_f,nrep,i_fft_init,d_fft_init)
    PERFOFF
  END SUBROUTINE fft_to_real

  !-------------------------------------------------------------------
  SUBROUTINE fft_thetab_cmplx(f,dir)
    COMPLEX(kind=8),intent(INOUT) :: f(nrp/2,nfs)

    INTEGER,intent(IN) :: dir            ! back or forth transform

    REAL(kind=8) :: work(nrp,nfs)

    !---------------------------------------------------------------

    PERFON('fft_thc')
    !WRITE(*,"(6(A,I6))") "fft_thetab_cmplx: nrp = ",nrp,", nfs = ", nfs,&
    !     &", size(f) = ",SIZE(f),", n_phi_max = ",n_phi_max,", dir = ",dir,", sizeThetaB = ",sizeThetaB
    !flush(6)
    CALL fftJW(f,nrp,n_phi_max,dir,sizeThetaB,work,nrp,nfs,i_fft_init,d_fft_init)
    PERFOFF
  END SUBROUTINE fft_thetab_cmplx

  SUBROUTINE fft_thetab_real(f,dir)
    REAL(kind=8) :: f(nrp,nfs)

    INTEGER :: dir            ! back or forth transform

    REAL(kind=8) :: work(nrp,nfs)

    PERFON('fft_thr')
    !---------------------------------------------------------------
    CALL fftJW(f,nrp,n_phi_max,dir,sizeThetaB, &
         work,nrp,nfs,i_fft_init,d_fft_init)
    PERFOFF
  END SUBROUTINE fft_thetab_real

END MODULE fft_JW
  

  SUBROUTINE fftJW(a,ld_a,n,isign,nsize,wrk,wd1,wd2, &
       i_fft_init,d_fft_init)
    !*******************************************************************************
    !  This file contains also the subroutines called by fftJW:
    !             fft99a, fft99b, wpass2, wpass3, wpass4 and wpass5
    !  The routines has been adopted by Gary Glatzmaier and has
    !  subsequently been modified by Uli Christensen and Johannes Wicht

    ! purpose      perform a number of simultaneous real/half-complex
    !              periodic fast fourier transforms or corresponding inverse
    !              transforms, using ordinary spatial order of
    !              gridpoint values.  given a set
    !              of real data vectors, the package returns a set of
    !              "half-complex" fourier coefficient vectors, or vice
    !              versa.  the length of the transforms must be an even
    !              number that has no other factors except possibly powers
    !              of 2, 3, and 5.  this version of fft991 is
    !              optimized for use on the cray-1.

    ! on input     a(ld_a,*)
    !               an array of length (ld_a,nsize) containing the input data
    !               or coefficient vectors.  This array is overwritten by
    !               the results.

    !              n
    !               the length of each transform (see definition of
    !               transforms, below).

    !              nsize
    !               the number of transforms to be done simultaneously.

    !              isign
    !               = +1 for a transform from fourier coefficients to
    !                    gridpoint values.
    !               = -1 for a transform from gridpoint values to fourier
    !                    coefficients.

    ! on output    a
    !               if isign = +1, and n_theta_max coefficient vectors are supplied
    !               each containing the sequence

    !               a(0),b(0),a(1),b(1),...,a(n/2),b(n/2)  (n+2 values)

    !               then the result consists of n_theta_max data vectors each
    !               containing the corresponding n+2 gridpoint values

    !               for fft991, x(0), x(1), x(2),...,x(n-1),0,0.
    !                    (n+2) real values with x(n)=x(n+1)=0

    !               when isign = +1, the transform is defined by
    !                 x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
    !                 where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
    !                 and i=sqrt (-1)
    !                    for k=0,...,n/2    i.e., (n/2+1) complex values
    !                    with c(0) = c(n) = a(0) and c(n/2)=a(n/2)=0

    !               if isign = -1, and n_theta_max data vectors are supplied each
    !               containing a sequence of gridpoint values x(j) as
    !               defined above, then the result consists of n_theta_max vectors
    !               each containing the corresponding fourier cofficients
    !               a(k), b(k), 0 .le. k .le n/2.

    !               when isign = -1, the inverse transform is defined by
    !                 c(k)=(1/n)*sum(j=0,...,n-1)(x(j)*exp(-2*i*j*k*pi/n))
    !                 where c(k)=a(k)+i*b(k) and i=sqrt(-1)
    !                 for k=0,...,n/2

    !               a call with isign=+1 followed by a call with isign=-1
    !               (or vice versa) returns the original data.

    !               note the fact that the gridpoint values x(j) are real
    !               implies that b(0)=b(n/2)=0.  for a call with isign=+1,
    !               it is not actually necessary to supply these zeros.
    !               note starting from grid with x(n)=x(n+1)=0
    !               then transforming to spectral (sign=-1)
    !               then c(n/2)=a(n/2) is not necessarily 0
    !               unless there is no aliasing.
    !---------------------------------------------------------------------------------------------

    !-- input:
    INTEGER :: ld_a           ! leading dimension of a
    REAL(kind=8) :: a(ld_a,*)       ! fields to be transformed
    INTEGER :: n              ! dimension of problem
    INTEGER :: isign          ! back/forth transtorm for isign=1/-1
    INTEGER :: nsize          ! number of fields for be transformed (second dim of a)

    INTEGER :: wd1,wd2
    REAL(kind=8) :: wrk(wd1,wd2)    ! work array

    INTEGER :: i_fft_init(*)  ! factorization information from init_fft
    REAL(kind=8) ::  d_fft_init(*)  ! trigonometric functions from init_fft

    !-- output: transformed a(ld_a,*)

    !-- local:
    LOGICAL :: tofro
    INTEGER :: nFactors,nodd,njap1,nrp
    INTEGER :: i,i0,i1,i2,ia,ic,id,j
    INTEGER :: la,k,fac

    !-- end of declaration
    !--------------------------------------------------------------------------------------------

    IF ( ld_a /= wd1 ) THEN
       WRITE(*,*) 'ERROR IN fftJW'
       WRITE(*,*) 'NOTE: first dim of work array has to be'
       WRITE(*,*) '      indentical to first dim of a!'
       STOP
    END IF
    IF ( nsize > wd2 ) THEN
       WRITE(*,*) 'TOO SMALL WORK ARRAY IN fftJW!'
       WRITE(*,*) 'SECOND DIM SHOULD BE >=:',nsize
       STOP
    END IF

    nFactors=i_fft_init(1)
    nodd    =MOD(nFactors,2)
    njap1   =n+1
    nrp     =n+2

    IF ( nrp < ld_a ) THEN
       WRITE(*,*) 'ERROR IN fftJW'
       WRITE(*,*) 'NOTE: first dim of input array a has to be'
       WRITE(*,*) '      at least n+2!'
       STOP
    END IF
    IF ( nrp < wd1 ) THEN
       WRITE(*,*) 'ERROR IN fftJW'
       WRITE(*,*) 'NOTE: first dim of work array a has to be'
       WRITE(*,*) '      at least n+2!'
       STOP
    END IF


    i0=0
    i1=i0+1
    i2=nsize

    IF ( isign == +1 ) THEN   ! Preprocessing
       CALL fft99aJW(a(1,i1),wrk,d_fft_init,nrp,nsize)
       tofro=.FALSE.
    ELSE
       tofro=.TRUE.
       IF ( nodd == 0 ) THEN
          DO i=1,nsize,2
             ia=i+1
             ic=i0+i
             id=ic+1
             DO j=1,n
                wrk(j,i) =a(j,ic)
                wrk(j,ia)=a(j,id)
             END DO
          END DO
          tofro=.FALSE.
       END IF
    END IF

    !--- Complex transform

    la=1
    DO k=1,nFactors
       fac=i_fft_init(k+1)
       IF ( tofro ) then
          IF ( fac == 2 ) THEN
             CALL wpass2JW(a(1,i1),a(2,i1),wrk,wrk(2,1), &
                  d_fft_init,nrp,nsize)
          ELSE IF ( fac == 3 ) THEN
             CALL wpass3JW(a(1,i1),a(2,i1),wrk,wrk(2,1), &
                  d_fft_init,nrp,la,nsize)
          ELSE IF ( fac == 4 ) THEN
             CALL wpass4JW(a(1,i1),a(2,i1),wrk,wrk(2,1), &
                  d_fft_init,nrp,la,nsize)
          ELSE IF ( fac == 5 ) THEN
             CALL wpass5JW(a(1,i1),a(2,i1),wrk,wrk(2,1), &
                  d_fft_init,nrp,la,nsize)
          END IF
       ELSE
          IF ( fac == 2 ) THEN
             CALL wpass2JW(wrk,wrk(2,1),a(1,i1),a(2,i1), &
                  d_fft_init,nrp,nsize)
          ELSE IF ( fac == 3 ) THEN
             CALL wpass3JW(wrk,wrk(2,1),a(1,i1),a(2,i1), &
                  d_fft_init,nrp,la,nsize)
          ELSE IF ( fac == 4 ) THEN
             CALL wpass4JW(wrk,wrk(2,1),a(1,i1),a(2,i1), &
                  d_fft_init,nrp,la,nsize)
          ELSE IF( fac == 5 ) THEN
             CALL wpass5JW(wrk,wrk(2,1),a(1,i1),a(2,i1), &
                  d_fft_init,nrp,la,nsize)
          END IF
       END IF
       la=la*fac
       tofro=.NOT.tofro
    END DO ! k

    IF ( isign == -1 ) THEN
       CALL fft99bJW(wrk,a(1,i1),d_fft_init,nrp,nsize)
    ELSE
       
       IF ( nodd == 0 ) THEN
          DO i=1,nsize,2
             ia=i+1
             ic=i0+i
             id=ic+1
             DO j=1,n
                a(j,ic)=wrk(j,i)
                a(j,id)=wrk(j,ia)
             END DO
          END DO
       END IF
       
       DO ic=i1,i2  ! Fill zeros into 2 extra elements
          a(njap1,ic)=0.D0
          a(nrp,ic)  =0.D0
       END DO
    END IF

  END SUBROUTINE fftJW

  !-----------------------------------------------------------------------


  !  **************        NEW SUBROUTINE        *************************

  subroutine fft99aJW(a,work,trigs,nrp,nsize)

    implicit none

    !-- input/output:
    integer :: nrp,nsize
    REAL(kind=8) :: a(*),work(*),trigs(*)

    !-- local:
    integer :: nja,njah,njap1
    integer :: ic,ia0,ib0,ia,ib,k,kk,kkmax
    REAL(kind=8) :: c,s

    !-- end of declaration
    !------------------------------------------------------------------------

    !     preprocessing step (isign=+1)
    !     (spectral to gridpoint transform)

    !     called in fftJW

    nja=nrp-2
    njah=nja/2
    njap1=nja+1
    kkmax=njah/2

    do ic=1,nsize
       ia0=(ic-1)*nrp
       ib0=ic*nrp
       ia=ia0+1
       ib=ib0-1
       work(ia)=a(ia)+a(ib)
       work(ia+1)=a(ia)-a(ib)
       k=1
       do kk=2,kkmax
          k=k+2
          c=trigs(nja+k)
          s=trigs(njap1+k)
          ia=ia0+k
          ib=ib0-k
          work(ia)=(a(ia)+a(ib))- &
               (s*(a(ia)-a(ib))+c*(a(ia+1)+a(ib+1)))
          work(ib)=(a(ia)+a(ib))+ &
               (s*(a(ia)-a(ib))+c*(a(ia+1)+a(ib+1)))
          work(ia+1)=(c*(a(ia)-a(ib))-s*(a(ia+1)+a(ib+1)))+ &
               (a(ia+1)-a(ib+1))
          work(ib+1)=(c*(a(ia)-a(ib))-s*(a(ia+1)+a(ib+1)))- &
               (a(ia+1)-a(ib+1))
       enddo
       ia=ia0+njah+1
       work(ia)=2.0*a(ia)
       work(ia+1)=-2.0*a(ia+1)
    enddo

    return
  end subroutine fft99aJW

  !  **************        NEW SUBROUTINE        *************************
  subroutine fft99bJW(work,a,trigsf,nrp,nsize)
    !-----------------------------------------------------------------------

    !     postprocessing step (isign=-1)
    !     (gridpoint to spectral transform)

    !     called in fftJW

    !-----------------------------------------------------------------------

    implicit none

    !-- input/output:
    integer :: nrp,nsize
    REAL(kind=8) :: work(*),a(*),trigsf(*)

    !-- local:
    integer :: nja,njah,kkmax,k,kk
    integer :: ia,ib,ic,ia0,ib0
    REAL(kind=8) :: scal1,scal2,s,c

    !-- end of declaration
    !-----------------------------------------------------------------------

    !     postprocessing step (isign=-1)
    !     (gridpoint to spectral transform)

    !     called in fftJW

    nja=nrp-2
    njah=nja/2
    kkmax=njah/2
    scal1=1.d0/dble(nja)
    scal2=0.5d0*scal1

    do ic=1,nsize
       ia0=(ic-1)*nrp
       ib0=ic*nrp
       ia=ia0+1
       ib=ib0-1
       a(ia)=scal1*(work(ia)+work(ia+1))
       a(ib)=scal1*(work(ia)-work(ia+1))
       k=1
       do kk=2,kkmax
          k=k+2
          c=trigsf(nja+k)
          s=trigsf(nja+k+1)
          ia=ia0+k
          ib=ib0-k
          a(ia)=scal2*((work(ia)+work(ib)) &
               +(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
          a(ib)=scal2*((work(ia)+work(ib)) &
               -(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
          a(ia+1)=scal2*((c*(work(ia)-work(ib)) &
               -s*(work(ia+1)+work(ib+1)))+(work(ib+1)-work(ia+1)))
          a(ib+1)=scal2*((c*(work(ia)-work(ib)) &
               -s*(work(ia+1)+work(ib+1)))-(work(ib+1)-work(ia+1)))
       enddo ! k
       ia=ia0+njah+1
       a(ia)=scal1*work(ia)
       a(ia+1)=-scal1*work(ia+1)
    enddo   ! ic

    return
  end subroutine fft99bJW

  !  **************        NEW SUBROUTINE        *************************
  subroutine wpass2JW(a,b,c,d,trigs,nrp,nsize)
    !-----------------------------------------------------------------------

    !     called in fftJW
    !     reduction for factor 2

    !     if(la.ne.1) stop 'call to wpass2 with la .ne. 1'

    !-----------------------------------------------------------------------

    implicit none

    !-- input/ouput:
    integer :: nrp,nsize
    REAL(kind=8) :: a(*),b(*),c(*),d(*),trigs(*)

    !-- local:
    integer :: i,j,ijk,iadd,in,n
    REAL(kind=8) :: c1,s1,an,bn

    !-- end of declaration
    !-----------------------------------------------------------------------

    n=(nrp-2)/2

    do ijk=1,nsize
       iadd=(ijk-1)*nrp - 1
       do in=2,n,2
          i=iadd+in
          j=i+in-2
          c1=trigs(in-1)
          s1=trigs(in)
          an=a(i)-a(n+i)
          bn=b(i)-b(n+i)
          c(j)=a(i)+a(n+i)
          d(j)=b(i)+b(n+i)
          c(2+j)=c1*an-s1*bn
          d(2+j)=s1*an+c1*bn
       enddo
    enddo

    return
  end subroutine wpass2JW

  !  **************        NEW SUBROUTINE        *************************
  subroutine wpass3JW(a,b,c,d,trigs,nrp,la,nsize)
    !-----------------------------------------------------------------------

    !     called in fftJW

    !-----------------------------------------------------------------------

    implicit none

    !-- input/output:
    integer :: nrp,la,nsize
    REAL(kind=8) :: a(*),b(*),c(*),d(*),trigs(*)

    !-- local
    integer :: n,m,iink,jink,jump,ib,jb,ic,jc
    integer :: ims,ijk,iadd,l,kc,kk,i,j,im

    integer, parameter :: mdim=4*180
    integer :: iindex(mdim),jindex(mdim)
    REAL(kind=8) :: c1(mdim),c2(mdim)
    REAL(kind=8) :: s1(mdim),s2(mdim)

    REAL(kind=8) :: sin60

    !-- end of declaration
    !------------------------------------------------------------------------


    sin60=0.866025403784437D0

    n=(nrp-2)/2
    m=n/3
    iink=m*2
    jink=la*2
    jump=2*jink

    !     coding for factor 3

    ib=iink
    jb=jink
    ic=ib+iink
    jc=jb+jink

    ims=1
    if(la < m .AND. la < 16) go to 65
    do ijk=1,nsize
       iadd=(ijk-1)*nrp-1
       do l=1,la
          i=l*2+iadd
          c(  i)=a(  i)+(a(ib+i)+a(ic+i))
          d(  i)=b(  i)+(b(ib+i)+b(ic+i))
          c(jb+i)=(a(  i)-0.5*(a(ib+i)+a(ic+i))) &
               -(sin60*(b(ib+i)-b(ic+i)))
          c(jc+i)=(a(  i)-0.5*(a(ib+i)+a(ic+i))) &
               +(sin60*(b(ib+i)-b(ic+i)))
          d(jb+i)=(b(  i)-0.5*(b(ib+i)+b(ic+i))) &
               +(sin60*(a(ib+i)-a(ic+i)))
          d(jc+i)=(b(  i)-0.5*(b(ib+i)+b(ic+i))) &
               -(sin60*(a(ib+i)-a(ic+i)))
       enddo
    enddo

    if (la == m) return
    ims=la+1

65  DO im=ims,m
       IF ( im > mdim ) THEN
          WRITE(*,*) 'Please increase mdim in wpass3!'
          WRITE(*,*) 'Should be at least:',m
          STOP
       END IF
       kc=(im-1)/la
       kk=kc*la
       c1(im)=trigs(2*kk+1)
       s1(im)=trigs(2*kk+2)
       c2(im)=trigs(4*kk+1)
       s2(im)=trigs(4*kk+2)
       iindex(im)=im*2
       jindex(im)=im*2+jump*kc
    enddo

    do ijk=1,nsize
       iadd=(ijk-1)*nrp - 1
       do im=ims,m
          i= iindex(im)+iadd
          j= jindex(im)+iadd
          c(  j)=a(  i)+(a(ib+i)+a(ic+i))
          d(  j)=b(  i)+(b(ib+i)+b(ic+i))
          c(jb+j)= &
               c1(im)*((a(  i)-0.5*(a(ib+i)+a(ic+i))) &
               -(sin60*(b(ib+i)-b(ic+i)))) &
               -s1(im)*((b(  i)-0.5*(b(ib+i)+b(ic+i))) &
               +(sin60*(a(ib+i)-a(ic+i))))
          d(jb+j)= &
               s1(im)*((a(  i)-0.5*(a(ib+i)+a(ic+i))) &
               -(sin60*(b(ib+i)-b(ic+i)))) &
               +c1(im)*((b(  i)-0.5*(b(ib+i)+b(ic+i))) &
               +(sin60*(a(ib+i)-a(ic+i))))
          c(jc+j)= &
               c2(im)*((a(  i)-0.5*(a(ib+i)+a(ic+i))) &
               +(sin60*(b(ib+i)-b(ic+i)))) &
               -s2(im)*((b(  i)-0.5*(b(ib+i)+b(ic+i))) &
               -(sin60*(a(ib+i)-a(ic+i))))
          d(jc+j)= &
               s2(im)*((a(  i)-0.5*(a(ib+i)+a(ic+i))) &
               +(sin60*(b(ib+i)-b(ic+i)))) &
               +c2(im)*((b(  i)-0.5*(b(ib+i)+b(ic+i))) &
               -(sin60*(a(ib+i)-a(ic+i))))
       enddo
    enddo

    return
  end subroutine wpass3JW

  !  **************        NEW SUBROUTINE        *************************


  subroutine wpass4JW(a,b,c,d,trigs,nrp,la,nsize)
    !-----------------------------------------------------------------------

    !     called in fftJW
    !     reduction for factor 4

    !-----------------------------------------------------------------------

    implicit none

    !-- input/output:
    integer :: nrp,la,nsize
    REAL(kind=8) :: a(*),b(*),c(*),d(*),trigs(*)

    !-- local:
    integer :: n,m,iink,jink,jump,ib,jb,ic,jc,id,jd
    integer :: ims,ijk,iadd,l,kc,kk,i,j,im
    REAL(kind=8) :: abdm,aacm,aac,abd
    REAL(kind=8) :: bbdm,bacm,bbd,bac

    integer, parameter :: mdim=4*135
    integer :: jindex(mdim)
    REAL(kind=8) :: c1(mdim),c2(mdim),c3(mdim)
    REAL(kind=8) :: s1(mdim),s2(mdim),s3(mdim)

    !-- end of declaration
    !------------------------------------------------------------------------

    n=(nrp-2)/2
    m=n/4
    iink=m*2
    jink=la*2
    jump=3*jink

    ib=iink
    jb=jink
    ic=ib+iink
    jc=jb+jink
    id=ic+iink
    jd=jc+jink

    ims=1
    if(la < m .AND. la < 64) go to 105
    do ijk=1,nsize
       iadd=(ijk-1)*nrp-1
       do l=1,la
          i=l*2+iadd
          aac=a(i)+a(ic+i)
          abd=a(ib+i)+a(id+i)
          bac=b(i)+b(ic+i)
          bbd=b(ib+i)+b(id+i)
          aacm=a(i)-a(ic+i)
          abdm=a(ib+i)-a(id+i)
          bacm=b(i)-b(ic+i)
          bbdm=b(ib+i)-b(id+i)
          c(   i)=aac+abd
          c(jc+i)=aac-abd
          d(   i)=bac+bbd
          d(jc+i)=bac-bbd
          c(jb+i)=aacm-bbdm
          c(jd+i)=aacm+bbdm
          d(jb+i)=bacm+abdm
          d(jd+i)=bacm-abdm
       enddo
    enddo

    if (la == m) return
    ims=la+1

105 do im=ims,m
       IF ( im > mdim ) THEN
          WRITE(*,*) 'Please increase mdim in wpass4!'
          WRITE(*,*) 'Should be at least:',m
          STOP
       END IF
       kc=(im-1)/la
       kk=kc*la
       c1(im)=trigs(2*kk+1)
       s1(im)=trigs(2*kk+2)
       c2(im)=trigs(4*kk+1)
       s2(im)=trigs(4*kk+2)
       c3(im)=trigs(6*kk+1)
       s3(im)=trigs(6*kk+2)
       jindex(im)=im*2+jump*kc
    enddo

    do ijk=1,nsize
       iadd=(ijk-1)*nrp - 1
       i=iadd+ims+ims-2
       do im=ims,m
          i=i+2
          j= jindex(im)+iadd
          aac=a(i)+a(ic+i)
          abd=a(ib+i)+a(id+i)
          bac=b(i)+b(ic+i)
          bbd=b(ib+i)+b(id+i)
          aacm=a(i)-a(ic+i)
          abdm=a(ib+i)-a(id+i)
          bacm=b(i)-b(ic+i)
          bbdm=b(ib+i)-b(id+i)
          c(   j)=aac+abd
          d(   j)=bac+bbd
          c(jc+j)=c2(im)*(aac-abd)-s2(im)*(bac-bbd)
          d(jc+j)=s2(im)*(aac-abd)+c2(im)*(bac-bbd)
          c(jb+j)=c1(im)*(aacm-bbdm)-s1(im)*(bacm+abdm)
          d(jb+j)=s1(im)*(aacm-bbdm)+c1(im)*(bacm+abdm)
          c(jd+j)=c3(im)*(aacm+bbdm)-s3(im)*(bacm-abdm)
          d(jd+j)=s3(im)*(aacm+bbdm)+c3(im)*(bacm-abdm)
       enddo
    enddo
    return

  end subroutine wpass4JW

  !  **************        NEW SUBROUTINE        *************************
  subroutine wpass5JW(a,b,c,d,trigs,nrp,la,nsize)
    !-----------------------------------------------------------------------

    !     called in fftJW
    !     reduction for factor 5

    !-----------------------------------------------------------------------

    implicit none

    !-- input/output:
    integer :: nrp,la,nsize
    REAL(kind=8) :: a(*),b(*),c(*),d(*),trigs(*)

    !-- local:
    integer :: n,m,iink,jink,jump
    integer :: ib,jb,ic,jc,id,jd,ie,je
    integer :: ims,ijk,iadd,l,kc,kk,i,j,im

    integer, parameter :: mdim=4*108
    integer :: iindex(mdim),jindex(mdim)
    REAL(kind=8) :: c1(mdim),c2(mdim),c3(mdim),c4(mdim)
    REAL(kind=8) :: s1(mdim),s2(mdim),s3(mdim),s4(mdim)

    REAL(kind=8) :: sin36,cos36,sin72,cos72

    !-- end of declaration
    !------------------------------------------------------------------------

    sin36=0.587785252292473D0
    cos36=0.809016994374947D0
    sin72=0.951056516295154D0
    cos72=0.309016994374947D0

    n=(nrp-2)/2
    m=n/5
    iink=m*2
    jink=la*2
    jump=4*jink

    ib= iink
    jb= jink
    ic=ib+iink
    jc=jb+jink
    id=ic+iink
    jd=jc+jink
    ie=id+iink
    je=jd+jink

    ims=1
    if(la < m .AND. la < 16) go to 145
    do ijk=1,nsize
       iadd=(ijk-1)*nrp-1
       do l=1,la
          i=l*2+iadd
          c( i)=a(i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
          d( i)=b(i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
          c(jb+i)=(a(i)+cos72*(a(ib+i)+a(ie+i)) &
               -cos36*(a(ic+i)+a(id+i))) &
               -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
          c(je+i)=(a(i)+cos72*(a(ib+i)+a(ie+i)) &
               -cos36*(a(ic+i)+a(id+i))) &
               +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
          d(jb+i)=(b(i)+cos72*(b(ib+i)+b(ie+i)) &
               -cos36*(b(ic+i)+b(id+i))) &
               +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
          d(je+i)=(b(i)+cos72*(b(ib+i)+b(ie+i)) &
               -cos36*(b(ic+i)+b(id+i))) &
               -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
          c(jc+i)=(a(i)-cos36*(a(ib+i)+a(ie+i)) &
               +cos72*(a(ic+i)+a(id+i))) &
               -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
          c(jd+i)=(a(i)-cos36*(a(ib+i)+a(ie+i)) &
               +cos72*(a(ic+i)+a(id+i))) &
               +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
          d(jc+i)=(b(i)-cos36*(b(ib+i)+b(ie+i)) &
               +cos72*(b(ic+i)+b(id+i))) &
               +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
          d(jd+i)=(b(i)-cos36*(b(ib+i)+b(ie+i)) &
               +cos72*(b(ic+i)+b(id+i))) &
               -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
       enddo
    enddo

    if (la == m) return
    ims=la+1

145 do im=ims,m
       IF ( im > mdim ) THEN
          WRITE(*,*) 'Please increase mdim in wpass5!'
          WRITE(*,*) 'Should be at least:',m
          STOP
       END IF
       kc=(im-1)/la
       kk=kc*la
       c1(im)=trigs(2*kk+1)
       s1(im)=trigs(2*kk+2)
       c2(im)=trigs(4*kk+1)
       s2(im)=trigs(4*kk+2)
       c3(im)=trigs(6*kk+1)
       s3(im)=trigs(6*kk+2)
       c4(im)=trigs(8*kk+1)
       s4(im)=trigs(8*kk+2)
       iindex(im)=im*2
       jindex(im)=im*2+jump*kc
    enddo

    do ijk=1,nsize
       iadd=(ijk-1)*nrp - 1
       do im=ims,m
          i= iindex(im)+iadd
          j= jindex(im)+iadd
          c( j)=a(i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
          d( j)=b(i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
          c(jb+j)= &
               c1(im)*((a(i)+cos72*(a(ib+i)+a(ie+i)) &
               -cos36*(a(ic+i)+a(id+i))) &
               -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))) &
               -s1(im)*((b(i)+cos72*(b(ib+i)+b(ie+i)) &
               -cos36*(b(ic+i)+b(id+i))) &
               +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
          d(jb+j)= &
               s1(im)*((a(i)+cos72*(a(ib+i)+a(ie+i)) &
               -cos36*(a(ic+i)+a(id+i))) &
               -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))) &
               +c1(im)*((b(i)+cos72*(b(ib+i)+b(ie+i)) &
               -cos36*(b(ic+i)+b(id+i))) &
               +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
          c(je+j)= &
               c4(im)*((a(i)+cos72*(a(ib+i)+a(ie+i)) &
               -cos36*(a(ic+i)+a(id+i))) &
               +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))) &
               -s4(im)*((b(i)+cos72*(b(ib+i)+b(ie+i)) &
               -cos36*(b(ic+i)+b(id+i))) &
               -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
          d(je+j)= &
               s4(im)*((a(i)+cos72*(a(ib+i)+a(ie+i)) &
               -cos36*(a(ic+i)+a(id+i))) &
               +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))) &
               +c4(im)*((b(i)+cos72*(b(ib+i)+b(ie+i)) &
               -cos36*(b(ic+i)+b(id+i))) &
               -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
          c(jc+j)= &
               c2(im)*((a(i)-cos36*(a(ib+i)+a(ie+i)) &
               +cos72*(a(ic+i)+a(id+i))) &
               -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))) &
               -s2(im)*((b(i)-cos36*(b(ib+i)+b(ie+i)) &
               +cos72*(b(ic+i)+b(id+i))) &
               +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
          d(jc+j)= &
               s2(im)*((a(i)-cos36*(a(ib+i)+a(ie+i)) &
               +cos72*(a(ic+i)+a(id+i))) &
               -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))) &
               +c2(im)*((b(i)-cos36*(b(ib+i)+b(ie+i)) &
               +cos72*(b(ic+i)+b(id+i))) &
               +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
          c(jd+j)= &
               c3(im)*((a(i)-cos36*(a(ib+i)+a(ie+i)) &
               +cos72*(a(ic+i)+a(id+i))) &
               +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))) &
               -s3(im)*((b(i)-cos36*(b(ib+i) &
               +b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
               -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
          d(jd+j)= &
               s3(im)*((a(i)-cos36*(a(ib+i)+a(ie+i)) &
               +cos72*(a(ic+i)+a(id+i))) &
               +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))) &
               +c3(im)*((b(i)-cos36*(b(ib+i)+b(ie+i)) &
               +cos72*(b(ic+i)+b(id+i))) &
               -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
       enddo
    enddo

    return
  end subroutine wpass5JW
