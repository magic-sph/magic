module useful
   !
   !  library with several useful subroutines
   !

   use precision_mod
   use parallel_mod
   use output_data, only: n_log_file, log_file
   use logic, only: l_save_out
   use constants, only: half, one, two

   implicit none

   private

   public :: l_correct_step, random, factorise, cc2real, cc22real, &
   &         logWrite, getMSD2, polynomial_interpolation, abortRun

contains

   logical function l_correct_step(n,t,t_last,n_max,n_step,n_intervals, &
                                   n_ts,times,n_eo)
      !
      ! Suppose we have a (loop) maximum of n_max steps!
      ! If n_intervals times in these steps a certain action should be carried out
      ! this can be invoked by l_correct_step=true if on input n_intervals>0
      ! and n_step=0.
      ! Alternatively the action can be invoked every n_step steps if
      ! on input n_intervals=0 and n_step>0.
      ! In both cases l_correct_step=true for n=n_max.
      !
      ! The argument controls whether in addition n should be
      !       even: n_eo=2
      !    or  odd: n_eo=1
      !

      !-- Input variables:
      integer,  intent(in) :: n            ! current step
      real(cp), intent(in) :: t            ! time at current step
      real(cp), intent(in) :: t_last       ! last time at current step
      integer,  intent(in) :: n_max        ! max number of steps
      integer,  intent(in) :: n_step       ! action interval
      integer,  intent(in) :: n_intervals ! number of actions
      integer,  intent(in) :: n_ts         ! number of times t
      real(cp), intent(in) :: times(*)     ! times where l_correct_step == true
      integer,  intent(in) :: n_eo         ! even/odd controller
  
      !-- Local variables:
      integer :: n_delta      ! corrector for even/odd n
      integer :: n_offset     ! offset with no action
      integer :: n_t          ! counter for times
      integer :: n_steps      ! local step width
  
  
      if ( n_step /= 0 .and. n_intervals /= 0 ) then
         write(*,*) '! Error message from function l_correct_step:'
         write(*,*) '! Either n_step or n_interval have to be zero!'
         call abortRun('Stop run in l_correct_step')
      end if
  
      l_correct_step=.false.
  
      if ( n_intervals /= 0 ) then
         n_steps=n_max/n_intervals
         if ( ( n_eo == 2 .and. mod(n_step,2) /= 0 ) .or. &
              ( n_eo == 1 .and. mod(n_step,2) /= 1 ) ) then
            n_steps=n_steps+1
         end if
  
         n_offset=n_max-n_steps*n_intervals
  
         if ( n > n_offset .and. mod(n-n_offset,n_steps) == 0 ) l_correct_step= .true. 
      else if ( n_step /= 0 ) then
         n_delta=0
         if ( ( n_eo == 1 .and. mod(n,2) == 0 ) .or. &
              ( n_eo == 2 .and. mod(n,2) == 1 ) ) then
            n_delta=1
         end if
  
         if ( n == n_max .or. mod(n-n_delta,n_step) == 0 ) l_correct_step= .true. 
      end if
  
      if ( n_ts >= 1 ) then
         do n_t=1,n_ts
            if ( times(n_t) < t .and. times(n_t) >= t_last ) then
               l_correct_step=.true.
               exit
            end if
         end do
      end if
               
   end function l_correct_step
!----------------------------------------------------------------------------
   real(cp) function random(r)
      !
      !     random number generator
      !
      !     if ( r == 0 ) then
      !        random(r) = next random number (between 0. and 1.)
      !     if ( r < 0 ) then
      !        random(r) = previous random number
      !     if ( r > 0 ) then
      !        random(r) = a new sequence of random numbers is started
      !                  with seed r mod 1
      !                  note: r must be a non-integer to get a different seq
      !
      !     called in sinit
      !

      !-- Input variables:
      real(cp), intent(in) :: r
  
      !-- Local variables:
      integer :: ia1, ia0, ia1ma0, ic
      integer :: iy0, iy1
      integer, save :: ix1, ix0
  
      ia1   =1536
      ia0   =1029
      ia1ma0=507
      ic    =1731
  
      if ( r == 0.0_cp ) then
         iy0 = ia0*ix0
         iy1 = ia1*ix1 + ia1ma0*(ix0-ix1) + iy0
         iy0 = iy0 + ic
  
         ix0 = mod(iy0,2048)
         iy1 = iy1 + (iy0-ix0)/2048
         ix1 = mod(iy1,2048)
      else if ( r > 0.0_cp ) then
         ix1 = int(mod(r,one)*4194304.0_cp + half)
         ix0 = mod(ix1,2048)
         ix1 = (ix1-ix0)/2048
      end if
       
      random = ix1*2048 + ix0
      random = random / 4194304.0_cp
  
   end function random
!----------------------------------------------------------------------------
   subroutine factorise(n,n_facs,fac,n_factors,factor)
      !
      !  Purpose of this subroutine is factorize n into a number of      
      !  given factors fac(i).                                           
      !
      

      !-- Input variables:
      integer, intent(in) :: n         ! number to be factorised    !
      integer, intent(in) :: fac(*)    ! list of fators to be tried !
      integer, intent(in) :: n_facs    ! number of facs to be tried !
  
      !-- Output:
      integer, intent(out) :: n_factors ! number of factors used
      integer, intent(out) :: factor(*) ! list of factors used
  
      !-- Local variables:
      character(len=14) :: str
      integer :: n_rest,n_fac
      integer :: factor_tot,factor_test
  
      if ( n < 1 )  then
         write(*,*) '! Error message from factorise:'
         write(*,*) '! n should be larger than 0.'
         call abortRun('Stop run in factorise')
      else if ( n == 1 ) then
         n_factors=0
         factor(1)=0
         return
      end if
  
      n_factors=0    ! number of factors
      factor_tot=1   ! total factor
      n_rest=n       ! remaining
  
      do n_fac=1,n_facs  ! use new factor !
         factor_test=fac(n_fac)
         do while ( mod(n_rest,factor_test) == 0 )
            n_rest=n_rest/factor_test
            n_factors=n_factors+1
            factor(n_factors)=factor_test
            factor_tot=factor_tot*factor_test
            if ( n_rest == 1 ) return
         end do
      end do
  
      if ( n_rest /= 1 ) then
         write(str,*) n
         call abortRun('Sorry, no factorisation possible of:'//trim(adjustl(str)))
      end if

   end subroutine factorise
!----------------------------------------------------------------------------
   real(cp)  function cc2real(c,m)

      !-- Input variables:
      complex(cp), intent(in) :: c
      integer,     intent(in) :: m
  
      if ( m == 0 ) then
         cc2real=real(c)*real(c)
      else
         cc2real=two*(real(c)*real(c)+aimag(c)*aimag(c))
      end if

   end function cc2real
!----------------------------------------------------------------------------
   real(cp) function cc22real(c1,c2,m)

      !-- Input variables:
      complex(cp), intent(in) :: c1,c2
      integer,         intent(in) :: m
  
      if ( m == 0 ) then
         cc22real=real(c1)*real(c2)
      else
         cc22real=two*(real(c1)*real(c2)+aimag(c1)*aimag(c2))
      end if

   end function cc22real
!----------------------------------------------------------------------------
   subroutine logWrite(message)

      !-- Input variable:
      character(len=*), intent(in) :: message

       if ( rank == 0 ) then
          if ( l_save_out ) then
             open(newunit=n_log_file, file=log_file, status='unknown', &
             &    position='append')
          end if
          write(n_log_file,*) trim(message)
          write(*,*)          trim(message)
          if ( l_save_out ) close(n_log_file)
       end if

   end subroutine logWrite
!----------------------------------------------------------------------------
   subroutine getMSD2(mean,SD,x,n,dt,totalTime)
      ! This subroutine computes the mean and standard deviation according 
      ! to a method introduced by Donald Knuth (1962). I rederived his formulas
      ! for a variable time step. On output SD still needs to be normalized with 
      ! the totalTime and then you have to take the square root!!
      ! The input integer counts the number of calls. For n=1 initialisation
      ! is necessary.

      !-- Input variables:
      real(cp), intent(in) :: x         ! quantity to be averaged
      real(cp), intent(in) :: dt        ! time since last averaging step
      real(cp), intent(in) :: totalTime ! total averaging time up to now
      integer,  intent(in) :: n         ! number of calls( only n=1 needed)

      !-- Output variables:
      real(cp), intent(out) :: mean     ! Time-average
      real(cp), intent(out) :: SD       ! Standard-deviation

      !-- Local variable:
      real(cp) :: delta

      if ( n == 1) then
         mean=x
         sd=0.0_cp
      else
         delta=x-mean
         mean=mean+delta*dt/totalTime
         SD=SD+dt*delta*(x-mean)
      end if

   end subroutine getMSD2
!----------------------------------------------------------------------------
   subroutine polynomial_interpolation(xold, yold, xnew ,ynew)

      !-- Input variables
      real(cp),    intent(in) :: xold(4) 
      complex(cp), intent(in) :: yold(4)
      real(cp),    intent(in) :: xnew

      !-- Output variables
      complex(cp), intent(out) :: ynew 

      !-- Local variables
      real(cp) :: yold_real(4), yold_imag(4)
      real(cp) :: ynew_real, ynew_imag

      yold_real= real(yold)
      yold_imag=aimag(yold)

      call polynomial_interpolation_real(xold, yold_real, xnew, ynew_real)
      call polynomial_interpolation_real(xold, yold_imag, xnew, ynew_imag)

      ynew = cmplx(ynew_real, ynew_imag, kind=cp)

   end subroutine polynomial_interpolation
!----------------------------------------------------------------------------
   subroutine polynomial_interpolation_real(xold,yold,xnew,ynew)

      !-- Input variables:
      real(cp), intent(in) :: xold(:)
      real(cp), intent(in) :: yold(:)
      real(cp), intent(in) :: xnew

      !-- Output variables:
      real(cp), intent(out) :: ynew

      !-- Local variables:
      integer :: n_stencil
      integer :: n_st, n_st_out, n_s
      real(cp) :: diff, diff_tmp, dy
      real(cp) :: ho, hp, den, work_diff
      real(cp), allocatable :: work1(:), work2(:)

      n_stencil=size(xold)
      allocate( work1(n_stencil), work2(n_stencil) )

      n_s=1
      diff=abs(xnew-xold(1))
      do n_st=1,n_stencil
         diff_tmp=abs(xnew-xold(n_st))
         if ( diff_tmp < diff ) then
            n_s =n_st
            diff=diff_tmp
         end if
         work1(n_st)=yold(n_st)
         work2(n_st)=yold(n_st)
      end do
      ynew=yold(n_s)

      n_s=n_s-1
      do n_st_out=1,n_stencil-1
         do n_st=1,n_stencil-n_st_out
            ho       =xold(n_st)-xnew
            hp       =xold(n_st+n_st_out)-xnew
            work_diff=work1(n_st+1)-work2(n_st)
            den=ho-hp
            if ( den == 0.0_cp ) call abortRun('Stop in polynomial interpolation')
            den        =work_diff/den
            work2(n_st)=hp*den
            work1(n_st)=ho*den
         end do
         if ( 2*n_s < n_stencil-n_st_out )then
            dy=work1(n_s+1)
         else
            dy=work2(n_s)
            n_s=n_s-1
         end if
         ynew=ynew+dy
      end do

      deallocate( work1, work2 )

   end subroutine polynomial_interpolation_real
!----------------------------------------------------------------------------
   subroutine abortRun(message)
      !
      ! This routine properly terminates a run
      !

      !-- Input variable
      character(len=*), intent(in) :: message

      !-- Local variables:
      integer :: code

      code = 32

      write(*,*) 
      write(*,*) 
      write(*,*) 
      write(*,*) '! Something went wrong, MagIC will stop now'
      write(*,*) '! See below the error message:'
      write(*,*) 
      write(*,*) message
      write(*,*) 

#ifdef WITH_MPI
      call MPI_Abort(MPI_COMM_WORLD, code, ierr)
#else
      stop
#endif

   end subroutine abortRun
!----------------------------------------------------------------------------
end module useful
