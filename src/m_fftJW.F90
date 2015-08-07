!$Id$
#include "perflib_preproc.cpp"

module fft_JW

   use useful, only: factorise
   use const, only: pi, sin36, sin60, sin72, cos36, cos72

   use truncation
   use blocking
   use parallel_mod
 
   implicit none
 
   private
 
   !-- For Fourier transform (calculated in init_fft)
   integer, parameter :: ni=100
   integer :: nd
   integer :: i_fft_init(ni)
   real(kind=8), allocatable ::  d_fft_init(:)
 
   public :: fft_thetab, init_fft, fft_to_real

contains
  
   subroutine init_fft(n)
      !  +-------------+----------------+------------------------------------+
      !  |  Purpose of this subroutine is to calculate and store several     |
      !  |  values that will be needed for a fast fft transform.             |
      !  |  The actual transform is performed by the subroutine fftJW.       |
      !  +-------------------------------------------------------------------+


      !-- Input variable:
      integer, intent(in) :: n     ! Dimension of problem, number of grid points


      !-- Local variables:
      integer :: i,j,nFacs,nFactors,help
      integer, parameter :: nFacsA=100
      integer :: fac(nFacsA),factor(nFacsA)
      real(kind=8) :: phi,dPhi

      nd = 3*(n/2)

      !-- For Fourier transform (calculated in init_fft)
      allocate( d_fft_init(nd) )



      !-- Checking number of datapoints:
      if ( n <= 3 ) then
         write(*,*) '! Message from subroutine init_fft:'
         write(*,*) '! Sorry, I need more than 3 grid points!'
         stop
      end if
      if ( mod(n,4) /= 0 ) then
         write(*,*) '! Note from subroutine init_fft:'
         write(*,*) '! Number of data points has to be'
         write(*,*) '! a mutiple of 4!'
         stop
      end if

      nFacs=4  ! factors to be tried
      fac(1)=4
      fac(2)=2
      fac(3)=3
      fac(4)=5

      call factorise(n/2,nFacs,fac,nFactors,factor)
      if ( nFactors > nFacsA ) then
         write(*,*) '! Message from subroutine init_fft:'
         write(*,*) '! Please increas nFacsA!'
         write(*,*) '! Should be >= ',nFactors
         stop
      end if

      !-- Sort in ascending order:
      if ( nFactors > 1 ) then
         do i=2,nFactors
            help=factor(i)
            do j=i-1,1,-1
               if ( factor(j) <= help ) goto 15
               factor(j+1)=factor(j) ! SHIFT UP
            end do
            j=0
15          factor(j+1)=help ! INSERT
         end do
      end if

      !-- Store:
      if ( ni <= nFactors+1 ) then
         write(*,*) '! Message from subroutine init_fft:'
         write(*,*) '! Increase dimension of array i_fft_init'
         write(*,*) '! in calling routine.'
         write(*,*) '! Should be at least:',nFactors+1
         stop
      end if
      i_fft_init(1)=nFactors
      do i=1,nFactors
         i_fft_init(i+1)=factor(i)
      end do

      !-- Calculate trigonometric functions:
      j=n/2
      if ( nd < n+j ) then
         write(*,*) '! Message from subroutine init_fft:'
         write(*,*) '! Increase dimension of array d_fft_init'
         write(*,*) '! in calling routine.'
         write(*,*) '! Should be at least:',n+j
         write(*,*) '! But is only       :',nd
         stop
      end if
      dPhi=2.D0*pi/dble(n)
      do i=1,n,2
         phi=dble(i-1)*dPhi
         d_fft_init(i)  =dcos(phi)
         d_fft_init(i+1)=dsin(phi)
      end do
      dPhi=0.5D0*dPhi
      do i=1,j,2
         phi=dble(i-1)*dPhi
         d_fft_init(n+i)  =dcos(phi)
         d_fft_init(n+i+1)=dsin(phi)
      end do

   end subroutine init_fft
!------------------------------------------------------------------------------
   subroutine fft_to_real(f,ld_f,nrep)

      integer,      intent(in) :: ld_f, nrep
      real(kind=8), intent(inout) :: f(ld_f, nrep)

      real(kind=8) :: work(ld_f,nrep)

      !PERFON('fft2r')
      call fftJW(f,ld_f,n_phi_max,1,nrep,work,ld_f,nrep,i_fft_init,d_fft_init)
      !PERFOFF

   end subroutine fft_to_real
!------------------------------------------------------------------------------
   subroutine fft_thetab(f,dir)

      real(kind=8), intent(inout) :: f(nrp,nfs)
      integer,      intent(in) :: dir            ! back or forth transform

      real(kind=8) :: work(nrp,nfs)

      !PERFON('fft_thr')
      call fftJW(f,nrp,n_phi_max,dir,sizeThetaB,work,nrp,nfs,i_fft_init,d_fft_init)
      !PERFOFF

   end subroutine fft_thetab
!------------------------------------------------------------------------------
   subroutine fftJW(a,ld_a,n,isign,nsize,wrk,wd1,wd2,i_fft_init,d_fft_init)
      !-------------------------------------------------------------------------
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
      !               a(k), b(k), 0  <=  k .le n/2.
  
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
      !-------------------------------------------------------------------------

      !-- Input variables:
      integer,      intent(in) :: ld_a         ! leading dimension of a
      real(kind=8), intent(inout) :: a(ld_a,*) ! fields to be transformed
      integer,      intent(in) :: n            ! dimension of problem
      integer,      intent(in) :: isign        ! back/forth transtorm for isign=1/-1
      integer,      intent(in) :: nsize        ! number of fields for 
                                               ! be transformed (second dim of a)
  
      integer,      intent(in) :: wd1,wd2
      real(kind=8), intent(inout) :: wrk(wd1,wd2)    ! work array
  
      integer,      intent(in) :: i_fft_init(*) ! factorization information from init_fft
      real(kind=8), intent(in) :: d_fft_init(*) ! trigonometric functions from init_fft
  
  
      !-- Local variables:
      logical :: tofro
      integer :: nFactors,nodd,njap1,nrp
      integer :: i,i0,i1,i2,ia,ic,id,j
      integer :: la,k,fac
  
  
      if ( ld_a /= wd1 ) then
         write(*,*) 'ERROR IN fftJW'
         write(*,*) 'NOTE: first dim of work array has to be'
         write(*,*) '      indentical to first dim of a!'
         stop
      end if
      if ( nsize > wd2 ) then
         write(*,*) 'TOO SMALL WORK ARRAY IN fftJW!'
         write(*,*) 'SECOND DIM SHOULD BE >=:',nsize
         stop
      end if
  
      nFactors=i_fft_init(1)
      nodd    =mod(nFactors,2)
      njap1   =n+1
      nrp     =n+2
  
      if ( nrp < ld_a ) then
         write(*,*) 'ERROR IN fftJW'
         write(*,*) 'NOTE: first dim of input array a has to be'
         write(*,*) '      at least n+2!'
         stop
      end if
      if ( nrp < wd1 ) then
         write(*,*) 'ERROR IN fftJW'
         write(*,*) 'NOTE: first dim of work array a has to be'
         write(*,*) '      at least n+2!'
         stop
      end if
  
  
      i0=0
      i1=i0+1
      i2=nsize
  
      if ( isign == +1 ) then   ! Preprocessing
         call fft99aJW(a(1,i1),wrk,d_fft_init,nrp,nsize)
         tofro=.false.
      else
         tofro=.true.
         if ( nodd == 0 ) then
            do i=1,nsize,2
               ia=i+1
               ic=i0+i
               id=ic+1
               do j=1,n
                  wrk(j,i) =a(j,ic)
                  wrk(j,ia)=a(j,id)
               end do
            end do
            tofro=.false.
         end if
      end if
  
      !--- Complex transform
  
      la=1
      do k=1,nFactors
         fac=i_fft_init(k+1)
         if ( tofro ) then
            if ( fac == 2 ) then
               call wpass2JW(a(1,i1),a(2,i1),wrk,wrk(2,1), &
                             d_fft_init,nrp,nsize)
            else if ( fac == 3 ) then
               call wpass3JW(a(1,i1),a(2,i1),wrk,wrk(2,1), &
                             d_fft_init,nrp,la,nsize)
            else if ( fac == 4 ) then
               call wpass4JW(a(1,i1),a(2,i1),wrk,wrk(2,1), &
                             d_fft_init,nrp,la,nsize)
            else if ( fac == 5 ) then
               call wpass5JW(a(1,i1),a(2,i1),wrk,wrk(2,1), &
                             d_fft_init,nrp,la,nsize)
            end if
         else
            if ( fac == 2 ) then
               call wpass2JW(wrk,wrk(2,1),a(1,i1),a(2,i1), &
                             d_fft_init,nrp,nsize)
            else if ( fac == 3 ) then
               call wpass3JW(wrk,wrk(2,1),a(1,i1),a(2,i1), &
                             d_fft_init,nrp,la,nsize)
            else if ( fac == 4 ) then
               call wpass4JW(wrk,wrk(2,1),a(1,i1),a(2,i1), &
                             d_fft_init,nrp,la,nsize)
            else if( fac == 5 ) then
               call wpass5JW(wrk,wrk(2,1),a(1,i1),a(2,i1), &
                             d_fft_init,nrp,la,nsize)
            end if
         end if
         la=la*fac
         tofro=.not.tofro
      end do ! k
  
      if ( isign == -1 ) then
         call fft99bJW(wrk,a(1,i1),d_fft_init,nrp,nsize)
      else
         
         if ( nodd == 0 ) then
            do i=1,nsize,2
               ia=i+1
               ic=i0+i
               id=ic+1
               do j=1,n
                  a(j,ic)=wrk(j,i)
                  a(j,id)=wrk(j,ia)
               end do
            end do
         end if
         
         do ic=i1,i2  ! Fill zeros into 2 extra elements
            a(njap1,ic)=0.D0
            a(nrp,ic)  =0.D0
         end do
      end if
  
    end subroutine fftJW
!------------------------------------------------------------------------------
   subroutine fft99aJW(a,work,trigs,nrp,nsize)

      !-- Input/output:
      integer,      intent(in) :: nrp,nsize
      real(kind=8), intent(in) :: trigs(*)
      real(kind=8), intent(inout) :: a(*),work(*)
  
      !-- Local variables:
      integer :: nja,njah,njap1
      integer :: ic,ia0,ib0,ia,ib,k,kk,kkmax
      real(kind=8) :: c,s
  
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
            work(ia)=(a(ia)+a(ib))- (s*(a(ia)-a(ib))+c*(a(ia+1)+a(ib+1)))
            work(ib)=(a(ia)+a(ib))+ (s*(a(ia)-a(ib))+c*(a(ia+1)+a(ib+1)))
            work(ia+1)=(c*(a(ia)-a(ib))-s*(a(ia+1)+a(ib+1)))+ (a(ia+1)-a(ib+1))
            work(ib+1)=(c*(a(ia)-a(ib))-s*(a(ia+1)+a(ib+1)))- (a(ia+1)-a(ib+1))
         enddo
         ia=ia0+njah+1
         work(ia)=2.0*a(ia)
         work(ia+1)=-2.0*a(ia+1)
      enddo

   end subroutine fft99aJW
!------------------------------------------------------------------------------
   subroutine fft99bJW(work,a,trigsf,nrp,nsize)
      !-----------------------------------------------------------------------
      !     postprocessing step (isign=-1)
      !     (gridpoint to spectral transform)
      !
      !     called in fftJW
      !
      !-----------------------------------------------------------------------
  
      !-- input/output:
      integer,      intent(in) :: nrp,nsize
      real(kind=8), intent(in) :: trigsf(*)
      real(kind=8), intent(inout) :: work(*),a(*)
  
      !-- Local variables:
      integer :: nja,njah,kkmax,k,kk
      integer :: ia,ib,ic,ia0,ib0
      real(kind=8) :: scal1,scal2,s,c
  
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

   end subroutine fft99bJW
!------------------------------------------------------------------------------
   subroutine wpass2JW(a,b,c,d,trigs,nrp,nsize)
      !-----------------------------------------------------------------------
  
      !     called in fftJW
      !     reduction for factor 2
  
      !     if(la /= 1) stop 'call to wpass2 with la  /=  1'
  
      !-----------------------------------------------------------------------
  
      !-- input/ouput:
      integer,      intent(in) :: nrp,nsize
      real(kind=8), intent(in) :: a(*),b(*),trigs(*)
      real(kind=8), intent(out) :: c(*),d(*)
  
      !-- Local variables:
      integer :: i,j,ijk,iadd,in,n
      real(kind=8) :: c1,s1,an,bn
  
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
  
   end subroutine wpass2JW
!------------------------------------------------------------------------------
   subroutine wpass3JW(a,b,c,d,trigs,nrp,la,nsize)
      !-----------------------------------------------------------------------
      !     called in fftJW
      !-----------------------------------------------------------------------
  
      !-- input/output:
      integer,      intent(in) :: nrp,la,nsize
      real(kind=8), intent(in) :: a(*),b(*),trigs(*)
      real(kind=8), intent(out) :: c(*),d(*)
  
      !-- Local variables:
      integer :: n,m,iink,jink,jump,ib,jb,ic,jc
      integer :: ims,ijk,iadd,l,kc,kk,i,j,im
  
      integer, parameter :: mdim=4*180
      integer :: iindex(mdim),jindex(mdim)
      real(kind=8) :: c1(mdim),c2(mdim)
      real(kind=8) :: s1(mdim),s2(mdim)
  
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
      if(la < m .and. la < 16) go to 65
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
  
65    do im=ims,m
         if ( im > mdim ) then
            write(*,*) 'Please increase mdim in wpass3!'
            write(*,*) 'Should be at least:',m
            stop
         end if
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
  
   end subroutine wpass3JW
!------------------------------------------------------------------------------
   subroutine wpass4JW(a,b,c,d,trigs,nrp,la,nsize)
      !-----------------------------------------------------------------------
      !     called in fftJW
      !     reduction for factor 4
      !-----------------------------------------------------------------------
  
      !-- input/output:
      integer,      intent(in) :: nrp,la,nsize
      real(kind=8), intent(in) :: a(*),b(*),trigs(*)
      real(kind=8), intent(out) :: c(*),d(*)
  
      !-- Local variables:
      integer :: n,m,iink,jink,jump,ib,jb,ic,jc,id,jd
      integer :: ims,ijk,iadd,l,kc,kk,i,j,im
      real(kind=8) :: abdm,aacm,aac,abd
      real(kind=8) :: bbdm,bacm,bbd,bac
  
      integer, parameter :: mdim=4*135
      integer :: jindex(mdim)
      real(kind=8) :: c1(mdim),c2(mdim),c3(mdim)
      real(kind=8) :: s1(mdim),s2(mdim),s3(mdim)
  
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
      if(la < m .and. la < 64) go to 105
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
  
105   do im=ims,m
         if ( im > mdim ) then
            write(*,*) 'Please increase mdim in wpass4!'
            write(*,*) 'Should be at least:',m
            stop
         end if
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

   end subroutine wpass4JW
!------------------------------------------------------------------------------
   subroutine wpass5JW(a,b,c,d,trigs,nrp,la,nsize)
      !-----------------------------------------------------------------------
      !     called in fftJW
      !     reduction for factor 5
      !-----------------------------------------------------------------------
  
      !-- input/output:
      integer,      intent(in) :: nrp,la,nsize
      real(kind=8), intent(in) :: a(*),b(*),trigs(*)
      real(kind=8), intent(out) :: c(*),d(*)
  
      !-- local:
      integer :: n,m,iink,jink,jump
      integer :: ib,jb,ic,jc,id,jd,ie,je
      integer :: ims,ijk,iadd,l,kc,kk,i,j,im
  
      integer, parameter :: mdim=4*108
      integer :: iindex(mdim),jindex(mdim)
      real(kind=8) :: c1(mdim),c2(mdim),c3(mdim),c4(mdim)
      real(kind=8) :: s1(mdim),s2(mdim),s3(mdim),s4(mdim)
  
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
      if(la < m .and. la < 16) go to 145
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

145   do im=ims,m
         if ( im > mdim ) then
            write(*,*) 'Please increase mdim in wpass5!'
            write(*,*) 'Should be at least:',m
            stop
         end if
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

   end subroutine wpass5JW
!------------------------------------------------------------------------------
end module fft_JW
