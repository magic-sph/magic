!$Id$
module useful
   !------------------------------------------------------------------------
   !  library with several useful shit
   !------------------------------------------------------------------------

   use precision_mod, only: cp
   use parallel_mod, only: rank
   use output_data, only: n_log_file, log_file
   use logic, only: l_save_out
   use const, only: half, one, two

   implicit none

   private

   public :: l_correct_step, random, factorise, cc2real, cc22real, &
             safeOpen, safeClose, logWrite

contains

   logical function l_correct_step(n,t,t_last,n_max,n_step,n_intervalls, &
                                   n_ts,times,n_eo)
      !------------------------------------------------------------------------
      ! Suppose we have a (loop) maximum of n_max steps!
      ! If n_intervalls times in these steps a certain action should be carried out
      ! this can be invoked by l_correct_step=true if on input n_intervalls>0
      ! and n_step=0.
      ! Alternatively the action can be invoked every n_step steps if
      ! on input n_intervalls=0 and n_step>0.
      ! In both cases l_correct_step=true for n=n_max.

      ! The argument controlls whether in addition n should be
      !       even: n_eo=2
      !    or  odd: n_eo=1
      !------------------------------------------------------------------------

      !-- Input variables:
      integer,  intent(in) :: n            ! current step
      real(cp), intent(in) :: t            ! time at current step
      real(cp), intent(in) :: t_last       ! last time at current step
      integer,  intent(in) :: n_max        ! max number of steps
      integer,  intent(in) :: n_step       ! action intervall
      integer,  intent(in) :: n_intervalls ! number of actions
      integer,  intent(in) :: n_ts         ! number of times t
      real(cp), intent(in) :: times(*)     ! times where l_correct_step == true
      integer,  intent(in) :: n_eo         ! even/odd controller
  
      !-- Local variables:
      integer :: n_delta      ! corrector for even/odd n
      integer :: n_offset     ! offset with no action
      integer :: n_t          ! counter for times
      integer :: n_steps      ! local step width
  
  
      if ( n_step /= 0 .and. n_intervalls /= 0 ) then
         write(*,*) '! ERROR MESSAGE FROM FUNCTION L_CORRECT_STEP:'
         write(*,*) '! EITHER N_STEP OR N_INTERVALL HAVE TO BE ZERO!'
         stop
      end if
  
      l_correct_step=.false.
  
      if ( n_intervalls /= 0 ) then
         n_steps=n_max/n_intervalls
         if ( ( n_eo == 2 .and. mod(n_step,2) /= 0 ) .or. &
              ( n_eo == 1 .and. mod(n_step,2) /= 1 ) ) then
            n_steps=n_steps+1
         end if
  
         n_offset=n_max-n_steps*n_intervalls
  
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
      !----------------------------------------------------------------------
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
      !-----------------------------------------------------------------------

      !-- Input variables:
      real(cp), intent(in) :: r
  
      !-- Local variables:
      integer :: ia1, ia0, ia1ma0, ic, ix1, ix0
      integer :: iy0, iy1
  
      save ix1, ix0
  
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
      !  +-------------+----------------+------------------------------------+
      !  |  Purpose of this subroutine is factorize n into a number of       |
      !  |  given factors fac(i).                                            |
      !  +-------------------------------------------------------------------+
      

      !-- Input variables:
      integer, intent(in) :: n         ! number to be factorised    !
      integer, intent(in) :: fac(*)    ! list of fators to be tried !
      integer, intent(in) :: n_facs    ! number of facs to be tried !
  
      !-- Output:
      integer, intent(out) :: n_factors ! number of factors used
      integer, intent(out) :: factor(*) ! list of factors used
  
      !-- Local variables:
      integer :: n_rest,n_fac
      integer :: factor_tot,factor_test
  
      if ( n < 1 )  then
         write(*,*) '! Error message from factorise:'
         write(*,*) '! n should be larger than 0.'
         stop
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
         write(*,*) 'Sorry, no factorisation possible of:',n
         stop
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
   subroutine safeOpen(nf, file_name)

      integer,          intent(in) :: nf
      character(len=*), intent(in) :: file_name
  
      if ( l_save_out ) then
         open(nf, file=file_name, status='unknown', position='append')
      end if

   end subroutine safeOpen
!----------------------------------------------------------------------------
   subroutine safeClose(nf)

      integer, intent(in) :: nf
  
      if ( l_save_out ) close(nf)

   end subroutine safeClose
!----------------------------------------------------------------------------
   subroutine logWrite(message)

      !-- Input variable:
      character(len=*), intent(in) :: message

       if ( rank == 0 ) then
          if ( l_save_out ) then
             open(n_log_file, file=log_file, status='unknown', position='append')
          end if
          write(n_log_file,*) trim(message)
          write(*,*)          trim(message)
          if ( l_save_out ) close(n_log_file)
       end if

   end subroutine logWrite
!----------------------------------------------------------------------------
end module useful
