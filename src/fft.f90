module fft
   !
   ! This module contains the native subroutines used to compute FFT's.
   ! They are based on the FFT99 package from Temperton:
   ! http://www.cesm.ucar.edu/models/cesm1.2/cesm/cesmBbrowser/html_code/cam/fft99.F90.html
   ! I simply got rid of the 'go to' and Fortran legacy statements
   ! Those transforms only work for prime decomposition that involve factors
   ! of 2, 3, 5
   !

   use iso_fortran_env, only: output_unit
   use precision_mod
   use useful, only: factorise, abortRun
   use constants, only: pi, sin36, sin60, sin72, cos36, cos72, one, two, half
   use truncation, only: n_phi_max, nlat_padded
   use parallel_mod
 
   implicit none
 
   private
 
   !-- For Fourier transform (calculated in init_fft)
   integer, parameter :: ni=100
   integer :: nd
   integer :: i_fft_init(ni)
   real(cp), allocatable ::  d_fft_init(:)
 
   public :: init_fft, fft_to_real, finalize_fft, fft_many, ifft_many

contains
  
   subroutine init_fft(n)
      !
      ! Purpose of this subroutine is to calculate and store several
      ! values that will be needed for a fast Fourier transforms.
      !

      !-- Input variable:
      integer, intent(in) :: n     ! Dimension of problem, number of grid points

      !-- Local variables:
      integer :: i,j,nFacs,nFactors,help
      integer, parameter :: nFacsA=100
      integer :: fac(nFacsA),factor(nFacsA)
      logical :: lAfter
      real(cp) :: phi,dPhi

      nd = 3*(n/2)+1

      !-- For Fourier transform (calculated in init_fft)
      allocate( d_fft_init(nd) )
      d_fft_init(:) = 0.0_cp

      !-- Checking number of datapoints:
      if ( n <= 3 ) then
         write(output_unit,*) '! Message from subroutine init_fft:'
         write(output_unit,*) '! Sorry, I need more than 3 grid points!'
         call abortRun('Stop run in fft')
      end if
      if ( mod(n,4) /= 0 ) then
         write(output_unit,*) '! Note from subroutine init_fft:'
         write(output_unit,*) '! Number of data points has to be'
         write(output_unit,*) '! a mutiple of 4!'
         call abortRun('Stop run in fft')
      end if

      nFacs=4  ! factors to be tried
      fac(1)=4
      fac(2)=2
      fac(3)=3
      fac(4)=5

      call factorise(n/2,nFacs,fac,nFactors,factor)
      if ( nFactors > nFacsA ) then
         write(output_unit,*) '! Message from subroutine init_fft:'
         write(output_unit,*) '! Please increas nFacsA!'
         write(output_unit,*) '! Should be >= ',nFactors
         call abortRun('Stop run in fft')
      end if

      !-- Sort in ascending order:
      if ( nFactors > 1 ) then
         do i=2,nFactors
            help=factor(i)
            do j=i-1,1,-1
               if ( factor(j) <= help ) then
                  factor(j+1)=help ! INSERT
                  lAfter=.false.
                  exit
               else
                  factor(j+1)=factor(j) ! SHIFT UP
                  lAfter=.true.
               end if
            end do
            if ( lAfter ) factor(1)=help ! INSERT
         end do
      end if

      !-- Store:
      if ( ni <= nFactors+1 ) then
         write(output_unit,*) '! Message from subroutine init_fft:'
         write(output_unit,*) '! Increase dimension of array i_fft_init'
         write(output_unit,*) '! in calling routine.'
         write(output_unit,*) '! Should be at least:',nFactors+1
         call abortRun('Stop run in fft')
      end if
      i_fft_init(1)=nFactors
      do i=1,nFactors
         i_fft_init(i+1)=factor(i)
      end do

      !-- Calculate trigonometric functions:
      j=n/2
      if ( nd < n+j ) then
         write(output_unit,*) '! Message from subroutine init_fft:'
         write(output_unit,*) '! Increase dimension of array d_fft_init'
         write(output_unit,*) '! in calling routine.'
         write(output_unit,*) '! Should be at least:',n+j
         write(output_unit,*) '! But is only       :',nd
         call abortRun('Stop run in fft')
      end if
      dPhi=two*pi/real(n,cp)
      do i=1,n,2
         phi=real(i-1,cp)*dPhi
         d_fft_init(i)  =cos(phi)
         d_fft_init(i+1)=sin(phi)
      end do
      dPhi=half*dPhi
      do i=1,j,2
         phi=real(i-1,cp)*dPhi
         d_fft_init(n+i)  =cos(phi)
         d_fft_init(n+i+1)=sin(phi)
      end do

   end subroutine init_fft
!------------------------------------------------------------------------------
   subroutine finalize_fft
      ! Memory deallocation of FFT help arrays

      deallocate(d_fft_init)

   end subroutine finalize_fft
!------------------------------------------------------------------------------
   subroutine fft_to_real(f,ld_f,nrep)

      !-- Input variable
      integer,  intent(in) :: ld_f, nrep

      !-- In/Out variable
      real(cp), intent(inout) :: f(ld_f, nrep)

      !-- Local variable
      real(cp) :: tmp(ld_f+2,nrep),work(ld_f+1,nrep)

      tmp(:ld_f,:)=f(:,:)
      !-- FFT along the first axis
      call fft991(tmp,work,d_fft_init,i_fft_init,1,ld_f+2,ld_f,nrep,1) 
      f(:,:)=tmp(:ld_f,:)

   end subroutine fft_to_real
!------------------------------------------------------------------------------
   subroutine fft_many(g,f)
      !
      ! Fourier transform:  f(nlat,nlon) -> fhat(nlat,nlon/2+1)
      !

      !-- Input variable
      real(cp), intent(in) :: g(:,:)

      !-- Output variable
      complex(cp), intent(out) :: f(nlat_padded,n_phi_max/2+1)

      !-- Local variables
      integer :: nt, np, nThStart, nThStop, nThLoc
      real(cp) :: tmp(nlat_padded,n_phi_max+2),work(nlat_padded,n_phi_max+1)

      !$omp parallel default(shared) private(nThStart, nThStop, nThLoc, tmp, work)
      nThStart=1; nThStop=nlat_padded
      call get_openmp_blocks(nThStart,nThStop)
      nThLoc = nThStop-nThStart+1

      !-- Copy in a larger array with the two trailing zeroes
      do np=1,n_phi_max
         do nt=1,nThLoc
            tmp(nt,np) =g(nt+nThStart-1,np)
         end do
      end do
      do np=n_phi_max+1,n_phi_max+2
         do nt=1,nThLoc
            tmp(nt,np)=0.0_cp
         end do
      end do

      !-- Phi first
      !call fft991(tmp,work,d_fft_init,i_fft_init,1,n_phi_max+2,n_phi_max,nlat_padded,-1) 
      !-- Theta first
      call fft991(tmp,work,d_fft_init,i_fft_init,nlat_padded,1,n_phi_max,nThLoc,-1) 

      !-- Real to complex
      do np=1,n_phi_max/2+1
         do nt=nThStart,nThStop
            f(nt,np)=cmplx(tmp(nt-nThStart+1,2*np-1),tmp(nt-nThstart+1,2*np),cp)
         end do
      end do
      !$omp end parallel

   end subroutine fft_many
!------------------------------------------------------------------------------
   subroutine ifft_many(f,g)
      !
      ! Inverse Fourier transform: fhat(nlat, nlon/2+1) -> f(nlat,nlon)
      !

      !-- Input variable
      complex(cp), intent(in) :: f(:,:)

      !-- Output variable
      real(cp), intent(out) :: g(nlat_padded,n_phi_max)

      !-- Local variables
      integer :: nt, np, nThStart, nThStop, nThloc
      real(cp) :: tmp(nlat_padded,n_phi_max+2),work(nlat_padded,n_phi_max+1)

      !$omp parallel default(shared) private(nThStart, nThStop, nThLoc, tmp, work)
      nThStart=1; nThStop=nlat_padded
      call get_openmp_blocks(nThStart,nThStop)
      nThLoc = nThStop-nThStart+1

      !-- Complex to real
      do np=1,n_phi_max/2+1
         do nt=1,nThLoc
            tmp(nt,2*np-1)=real (f(nt+nThStart-1,np))
            tmp(nt,2*np)  =aimag(f(nt+nThStart-1,np))
         end do
      end do

      !-- Phi-first
      !call fft991(tmp,work,d_fft_init,i_fft_init,1,n_phi_max+2,n_phi_max,nlat_padded,1) 
      !-- Theta first
      call fft991(tmp,work,d_fft_init,i_fft_init,nlat_padded,1,n_phi_max,nThloc,1)

      !-- Copy the solution into g, removing the two last trailing
      do np=1,n_phi_max
         do nt=nThStart,nThStop
            g(nt,np)=tmp(nt-nThStart+1,np)
         end do
      end do
      !$omp end parallel

   end subroutine ifft_many
!------------------------------------------------------------------------------
   subroutine fft991(a,work,trigs,ifax,inc,jump,n,lot,isign)
      !
      !     subroutine "fft991" - multiple real/half-complex periodic
      !     fast fourier transform
      !
      !     same as fft99 except that ordering of data corresponds to
      !     that in mrfft2
      !
      !     procedure used to convert to half-length complex transform
      !     is given by cooley, lewis and welch (j. sound vib., vol. 12
      !     (1970), 315-337)
      !
      !     Definition of transforms:
      !
      !     isign=+1: x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
      !     where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
      !
      !     isign=-1: a(k)=(1/n)*sum(j=0,...,n-1)(x(j)*cos(2*j*k*pi/n))
      !     b(k)=-(1/n)*sum(j=0,...,n-1)(x(j)*sin(2*j*k*pi/n))
      !

      !-- Input variables
      integer,  intent(in)    :: inc      ! increment within each data 'vector'
      integer,  intent(in)    :: jump     ! increment between the start of each data vector
      integer,  intent(in)    :: n        ! length of the data vectors
      integer,  intent(in)    :: lot      ! number of data vectors
      integer,  intent(in)    :: isign    ! sign of the FFT
      integer,  intent(inout) :: ifax(:)  ! previously prepared list of factors of n/2
      real(cp), intent(in)    :: trigs(:) ! previously prepared list of trig function values

      !-- Output variables
      real(cp), intent(inout) :: a(*) ! array containing input and output data
      real(cp), intent(inout) :: work((n+1)*lot) ! array of size (n+1)*lot

      !-- Local variables
      integer :: nfax, nx, nh, ink, igo, ibase, jbase
      integer :: i, j, k, L, m, ia, la, ib

      nfax=ifax(1)
      nx=n+1
      nh=n/2
      ink=inc+inc
      if (isign /= 1) then
         !-- If necessary, transfer data to work area
         igo=50
         if (mod(nfax,2) /= 1) then
            ibase=1
            jbase=1
            do L=1,lot
               i=ibase
               j=jbase
               do m=1,n
                  work(j)=a(i)
                  i=i+inc
                  j=j+1
               end do
               ibase=ibase+jump
               jbase=jbase+nx
            end do
       
            igo=60
         end if
      else
         !  preprocessing (isign=+1)
         !  ------------------------
         call fft99a(a,work,trigs,inc,jump,n,lot)
         igo=60
      end if

      !   complex transform
      !   -----------------
      ia=1
      la=1
      do k=1,nfax
         if (igo /= 60) then
            call vpassm(a(ia),a(ia+inc),work(1),work(2),trigs, &
                 &      ink,2,jump,nx,lot,nh,ifax(k+1),la)
            igo=60
         else
            call vpassm(work(1),work(2),a(ia),a(ia+inc),trigs, &
                 &      2,ink,nx,jump,lot,nh,ifax(k+1),la)
            igo=50
         end if
         la=la*ifax(k+1)
      end do

      if (isign == -1) then
         !     postprocessing (isign=-1):
         !     --------------------------
         call fft99b(work,a,trigs,inc,jump,n,lot)
         return
      end if

      !-- if necessary, transfer data from work area
      if (mod(nfax,2) /= 1) then
         ibase=1
         jbase=1
         do L=1,lot
            i=ibase
            j=jbase
            do m=1,n
               a(j)=work(i)
               i=i+1
               j=j+inc
            end do
            ibase=ibase+nx
            jbase=jbase+jump
         end do
      end if

      !-- Fill in zeros at end
      ib=n*inc+1
      do L=1,lot
         a(ib)=0.0_cp
         a(ib+inc)=0.0_cp
         ib=ib+jump
      end do

   end subroutine fft991
!-------------------------------------------------------------------------------
   subroutine fft99a(a,work,trigs,inc,jump,n,lot)

      !-- Input variables:
      integer,  intent(in)    :: inc,jump,n,lot
      real(cp), intent(in)    :: trigs(:)

      !-- Output variable:
      real(cp),    intent(inout) :: a(*),work(*)

      !-- Local variables
      integer :: nh, nx, ink, k, L
      integer :: ia, ib, ja, jb, iabase, ibbase, jabase, jbbase
      real(cp) :: c, s

      nh=n/2
      nx=n+1
      ink=inc+inc

      !--   a(0) and a(n/2)
      ia=1
      ib=n*inc+1
      ja=1
      jb=2
      do L=1,lot
         work(ja)=a(ia)+a(ib)
         work(jb)=a(ia)-a(ib)
         ia=ia+jump
         ib=ib+jump
         ja=ja+nx
         jb=jb+nx
      end do
    
      !-- Remaining wavenumbers
      iabase=2*inc+1
      ibbase=(n-2)*inc+1
      jabase=3
      jbbase=n-1

      do k=3,nh,2
         ia=iabase
         ib=ibbase
         ja=jabase
         jb=jbbase
         c=trigs(n+k)
         s=trigs(n+k+1)
         do L=1,lot
            work(ja)=(a(ia)+a(ib))-(s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
            work(jb)=(a(ia)+a(ib))+(s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
            work(ja+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))+(a(ia+inc)-a(ib+inc))
            work(jb+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))-(a(ia+inc)-a(ib+inc))
            ia=ia+jump
            ib=ib+jump
            ja=ja+nx
            jb=jb+nx
         end do
         iabase=iabase+ink
         ibbase=ibbase-ink
         jabase=jabase+2
         jbbase=jbbase-2
       end do

      !--  Wavenumber n/4 (if it exists)
      if (iabase == ibbase) then
         ia=iabase
         ja=jabase
         do L=1,lot
            work(ja)  =two*a(ia)
            work(ja+1)=-two*a(ia+inc)
            ia=ia+jump
            ja=ja+nx
         end do
      end if

   end subroutine fft99a
!------------------------------------------------------------------------------
   subroutine fft99b(work,a,trigs,inc,jump,n,lot)

      !-- Input variables
      integer,  intent(in)    :: inc,jump,n,lot
      real(cp), intent(in)    :: trigs(:)

      !-- Output variables
      real(cp), intent(inout) :: a(*),work(*)

      !-- Local variables
      integer :: nh, nx, ink, k, L
      integer :: ia, ib, ja, jb, iabase, ibbase, jabase, jbbase
      real(cp) :: sca, c, s

      nh=n/2
      nx=n+1
      ink=inc+inc

      !--   a(0) and a(n/2)
      sca=one/real(n)
      ia=1
      ib=2
      ja=1
      jb=n*inc+1
      do L=1,lot
         a(ja)=sca*(work(ia)+work(ib))
         a(jb)=sca*(work(ia)-work(ib))
         a(ja+inc)=0.0_cp
         a(jb+inc)=0.0_cp
         ia=ia+nx
         ib=ib+nx
         ja=ja+jump
         jb=jb+jump
      end do

      !-- Remaining wavenumbers
      sca=half*sca
      iabase=3
      ibbase=n-1
      jabase=2*inc+1
      jbbase=(n-2)*inc+1

      do k=3,nh,2
         ia=iabase
         ib=ibbase
         ja=jabase
         jb=jbbase
         c=trigs(n+k)
         s=trigs(n+k+1)
         do L=1,lot
            a(ja)=sca*((work(ia)+work(ib)) &
            &  +(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
            a(jb)=sca*((work(ia)+work(ib)) &
            &  -(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
            a(ja+inc)=sca*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1))) &
            &   +(work(ib+1)-work(ia+1)))
            a(jb+inc)=sca*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1))) &
            &   -(work(ib+1)-work(ia+1)))
            ia=ia+nx
            ib=ib+nx
            ja=ja+jump
            jb=jb+jump
         end do
         iabase=iabase+2
         ibbase=ibbase-2
         jabase=jabase+ink
         jbbase=jbbase-ink
      end do

      !-- Wavenumber n/4 (if it exists)
      if (iabase == ibbase) then
         ia=iabase
         ja=jabase
         sca=two*sca
         do L=1,lot
            a(ja)=sca*work(ia)
            a(ja+inc)=-sca*work(ia+1)
            ia=ia+nx
            ja=ja+jump
         end do
      end if

   end subroutine fft99b
!-------------------------------------------------------------------------------
   subroutine vpassm(a,b,c,d,trigs,inc1,inc2,inc3,inc4,lot,n,ifac,la)

      !-- Input variables
      integer,  intent(in)  :: inc1, inc2, inc3, inc4, lot, n, ifac, la
      real(cp), intent(in)  :: a(*),b(*),trigs(*)

      !-- Output variables
      real(cp), intent(out) :: c(*),d(*)

      !-- Local variables
      integer :: i, j, k, L, m, iink, jink, jump, ibase, jbase, igo, ijk, la1
      integer :: ia, ja, ib, jb, kb, ic, jc, kc, id, jd, kd, ie, je, ke
      real(cp) :: c1, s1, c2, s2, c3, s3, c4, s4

      m=n/ifac
      iink=m*inc1
      jink=la*inc2
      jump=(ifac-1)*jink
      ibase=0
      jbase=0
      igo=ifac-1
      if (igo > 4) return

      select case (igo)

      !-- Coding for factor 2
      case (1)
         ia=1
         ja=1
         ib=ia+iink
         jb=ja+jink
         do L=1,la
            i=ibase
            j=jbase
            do ijk=1,lot
               c(ja+j)=a(ia+i)+a(ib+i)
               d(ja+j)=b(ia+i)+b(ib+i)
               c(jb+j)=a(ia+i)-a(ib+i)
               d(jb+j)=b(ia+i)-b(ib+i)
               i=i+inc3
               j=j+inc4
            end do
            ibase=ibase+inc1
            jbase=jbase+inc2
         end do
         if (la == m) return
         la1=la+1
         jbase=jbase+jump
         do k=la1,m,la
            kb=k+k-2
            c1=trigs(kb+1)
            s1=trigs(kb+2)
            do L=1,la
               i=ibase
               j=jbase
               do ijk=1,lot
                  c(ja+j)=a(ia+i)+a(ib+i)
                  d(ja+j)=b(ia+i)+b(ib+i)
                  c(jb+j)=c1*(a(ia+i)-a(ib+i))-s1*(b(ia+i)-b(ib+i))
                  d(jb+j)=s1*(a(ia+i)-a(ib+i))+c1*(b(ia+i)-b(ib+i))
                  i=i+inc3
                  j=j+inc4
               end do
               ibase=ibase+inc1
               jbase=jbase+inc2
            end do
            jbase=jbase+jump
         end do

      !-- Coding for factor 3
      case (2)
         ia=1
         ja=1
         ib=ia+iink
         jb=ja+jink
         ic=ib+iink
         jc=jb+jink
         do L=1,la
            i=ibase
            j=jbase
            do ijk=1,lot
               c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
               d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
               c(jb+j)=(a(ia+i)-half*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))
               c(jc+j)=(a(ia+i)-half*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))
               d(jb+j)=(b(ia+i)-half*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i)))
               d(jc+j)=(b(ia+i)-half*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i)))
               i=i+inc3
               j=j+inc4
            end do
            ibase=ibase+inc1
            jbase=jbase+inc2
         end do
         if (la == m) return
         la1=la+1
         jbase=jbase+jump
         do k=la1,m,la
            kb=k+k-2
            kc=kb+kb
            c1=trigs(kb+1)
            s1=trigs(kb+2)
            c2=trigs(kc+1)
            s2=trigs(kc+2)
            do L=1,la
               i=ibase
               j=jbase
               do ijk=1,lot
                  c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
                  d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
                  c(jb+j)=                                                            &
                  &   c1*((a(ia+i)-half*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))) &
                  &  -s1*((b(ia+i)-half*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
                  d(jb+j)=                                                            &
                  &   s1*((a(ia+i)-half*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))) &
                  &  +c1*((b(ia+i)-half*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
                  c(jc+j)=                                                            &
                  &   c2*((a(ia+i)-half*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))) &
                  &  -s2*((b(ia+i)-half*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
                  d(jc+j)=                                                            &
                  &   s2*((a(ia+i)-half*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))) &
                  &  +c2*((b(ia+i)-half*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
                  i=i+inc3
                  j=j+inc4
               end do
               ibase=ibase+inc1
               jbase=jbase+inc2
            end do
            jbase=jbase+jump
         end do

      !-- Coding for factor 4
      case (3)
         ia=1
         ja=1
         ib=ia+iink
         jb=ja+jink
         ic=ib+iink
         jc=jb+jink
         id=ic+iink
         jd=jc+jink
         do L=1,la
            i=ibase
            j=jbase
            do ijk=1,lot
               c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
               c(jc+j)=(a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))
               d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
               d(jc+j)=(b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i))
               c(jb+j)=(a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))
               c(jd+j)=(a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))
               d(jb+j)=(b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i))
               d(jd+j)=(b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i))
               i=i+inc3
               j=j+inc4
            end do
            ibase=ibase+inc1
            jbase=jbase+inc2
         end do
         if (la == m) return
         la1=la+1
         jbase=jbase+jump
         do k=la1,m,la
            kb=k+k-2
            kc=kb+kb
            kd=kc+kb
            c1=trigs(kb+1)
            s1=trigs(kb+2)
            c2=trigs(kc+1)
            s2=trigs(kc+2)
            c3=trigs(kd+1)
            s3=trigs(kd+2)
            do L=1,la
               i=ibase
               j=jbase
               do ijk=1,lot
                  c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
                  d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
                  c(jc+j)=                                     &
                  &   c2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) &
                  &  -s2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
                  d(jc+j)=                                     &
                  &   s2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) &
                  &  +c2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
                  c(jb+j)=                                     &
                  &   c1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))) &
                  &  -s1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
                  d(jb+j)=                                     &
                  &   s1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))) &
                  &  +c1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
                  c(jd+j)=                                     &
                  &   c3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))) &
                  &  -s3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
                  d(jd+j)=                                     &
                  &   s3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))) &
                  &  +c3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
                  i=i+inc3
                  j=j+inc4
               end do
               ibase=ibase+inc1
               jbase=jbase+inc2
            end do
            jbase=jbase+jump
         end do

      !-- Coding for factor 5
      case (4)
         ia=1
         ja=1
         ib=ia+iink
         jb=ja+jink
         ic=ib+iink
         jc=jb+jink
         id=ic+iink
         jd=jc+jink
         ie=id+iink
         je=jd+jink
         do L=1,la
            i=ibase
            j=jbase
            do ijk=1,lot
               c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
               d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
               c(jb+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
               & -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
               c(je+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
               & +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
               d(jb+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
               & +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
               d(je+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
               & -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
               c(jc+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
               & -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
               c(jd+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
               & +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
               d(jc+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
               & +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
               d(jd+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
               & -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
               i=i+inc3
               j=j+inc4
            end do
            ibase=ibase+inc1
            jbase=jbase+inc2
         end do
         if (la == m) return
         la1=la+1
         jbase=jbase+jump
         do k=la1,m,la
            kb=k+k-2
            kc=kb+kb
            kd=kc+kb
            ke=kd+kb
            c1=trigs(kb+1)
            s1=trigs(kb+2)
            c2=trigs(kc+1)
            s2=trigs(kc+2)
            c3=trigs(kd+1)
            s3=trigs(kd+2)
            c4=trigs(ke+1)
            s4=trigs(ke+2)
            do L=1,la
               i=ibase
               j=jbase
               do ijk=1,lot
                  c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
                  d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
                  c(jb+j)=                                                          &
                  &   c1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
                  &     -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))         &
                  &  -s1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
                  &     +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
                  d(jb+j)=                                                          &
                  &   s1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
                  &     -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))         &
                  &  +c1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
                  &     +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
                  c(je+j)=                                                          &
                  &   c4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
                  &     +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))         &
                  &  -s4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
                  &     -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
                  d(je+j)=                                                          &
                  &   s4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
                  &     +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))         &
                  &  +c4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
                  &     -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
                  c(jc+j)=                                                          &
                  &   c2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
                  &     -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))         &
                  &  -s2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
                  &     +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
                  d(jc+j)=                                                          &
                  &   s2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
                  &     -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))         &
                  &  +c2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
                  &     +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
                  c(jd+j)=                                                          &
                  &   c3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
                  &     +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))         &
                  &  -s3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
                  &     -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
                  d(jd+j)=                                                          &
                  &   s3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
                  &     +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))         &
                  &  +c3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
                  &     -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
                  i=i+inc3
                  j=j+inc4
               end do
               ibase=ibase+inc1
               jbase=jbase+inc2
            end do
            jbase=jbase+jump
         end do

      end select

   end subroutine vpassm
!------------------------------------------------------------------------------
end module fft
