module multistep_schemes

   use precision_mod
   use iso_fortran_env, only: output_unit
   use parallel_mod
   use num_param, only: alpha
   use constants, only: one, half, two, ci, zero
   use mem_alloc, only: bytes_allocated
   use useful, only: abortRun
   use output_data, only: log_file
   use logic, only: l_save_out
   use radial_functions, only: or1
   use time_schemes, only: type_tscheme
   use time_array

   implicit none

   private

   type, public, extends(type_tscheme) :: type_multistep
      real(cp), allocatable :: wimp(:)
      real(cp), allocatable :: wexp(:)
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: set_weights
      procedure :: set_dt_array
      procedure :: set_imex_rhs
      procedure :: set_imex_rhs_scalar
      procedure :: rotate_imex
      procedure :: rotate_imex_scalar
      procedure :: bridge_with_cnab2
      procedure :: start_with_ab1
      procedure :: get_time_stage
   end type type_multistep

contains

   subroutine initialize(this, time_scheme, courfac_nml, intfac_nml, alffac_nml)

      class(type_multistep) :: this

      !-- Input/output variables
      real(cp),          intent(in) :: courfac_nml
      real(cp),          intent(in) :: intfac_nml
      real(cp),          intent(in) :: alffac_nml
      character(len=72), intent(inout) :: time_scheme

      !-- Local variables
      real(cp) :: courfac_loc, alffac_loc, intfac_loc

      !-- Number of stages per iteration is always one in this case
      this%nstages = 1
      this%istage = 1
      this%family = 'MULTISTEP'

      allocate( this%l_exp_calc(1) )
      allocate( this%l_imp_calc_rhs(1) )
      this%l_exp_calc(1) = .true.
      this%l_imp_calc_rhs(1) = .true.
      bytes_allocated = bytes_allocated+2*SIZEOF_LOGICAL

      if ( index(time_scheme, 'CNAB2') /= 0 ) then
         this%time_scheme = 'CNAB2'
         this%norder_imp_lin = 2
         this%nold = 1
         this%norder_exp = 2
         courfac_loc = 2.5_cp
         alffac_loc  = 1.0_cp
         intfac_loc  = 0.15_cp
      else if ( index(time_scheme, 'MODCNAB') /= 0 ) then
         this%time_scheme = 'MODCNAB'
         this%nold = 2
         this%norder_imp_lin = 3
         this%norder_exp = 2
         courfac_loc = 2.5_cp
         alffac_loc  = 1.0_cp
         intfac_loc  = 0.15_cp
      else if ( index(time_scheme, 'CNLF') /= 0 ) then
         this%time_scheme = 'CNLF'
         this%nold = 2
         this%norder_imp_lin = 3
         this%norder_exp = 2
         courfac_loc = 2.5_cp
         alffac_loc  = 1.0_cp
         intfac_loc  = 0.15_cp
      else if ( index(time_scheme, 'SBDF2') /= 0 ) then
         this%time_scheme = 'SBDF2'
         this%nold = 2
         this%norder_imp_lin = 2 ! it should be one but we need to restart
         this%norder_exp = 2
         this%l_imp_calc_rhs(1) = .false.
         courfac_loc = 2.5_cp
         alffac_loc  = 1.0_cp
         intfac_loc  = 0.15_cp
      else if ( index(time_scheme, 'SBDF3') /= 0 ) then
         this%time_scheme = 'SBDF3'
         this%nold = 3
         this%norder_imp_lin = 2 ! it should be one but we need to restart
         this%norder_exp = 3
         this%l_imp_calc_rhs(1) = .false.
         courfac_loc = 4.0_cp
         alffac_loc  = 1.6_cp
         intfac_loc  = 0.09_cp
      else if ( index(time_scheme, 'TVB33') /= 0 ) then
         this%time_scheme = 'TVB33'
         this%nold = 3
         this%norder_imp_lin = 4 ! it should be one but we need to restart
         this%norder_exp = 3
         courfac_loc = 4.0_cp
         alffac_loc  = 1.6_cp
         intfac_loc  = 0.09_cp
      else if ( index(time_scheme, 'SBDF4') /= 0 ) then
         this%time_scheme = 'SBDF4'
         this%nold = 4
         this%norder_imp_lin = 2 ! it should be one but we need to restart
         this%norder_exp = 4
         this%l_imp_calc_rhs(1) = .false.
         courfac_loc = 5.5_cp
         alffac_loc  = 2.2_cp
         intfac_loc  = 0.065_cp
      end if

      if ( abs(courfac_nml) >= 1.0e3_cp ) then
         this%courfac=courfac_loc
      else
         this%courfac=courfac_nml
      end if

      if ( abs(alffac_nml) >= 1.0e3_cp ) then
         this%alffac=alffac_loc
      else
         this%alffac=alffac_nml
      end if

      if ( abs(intfac_nml) >= 1.0e3_cp ) then
         this%intfac=intfac_loc
      else
         this%intfac=intfac_nml
      end if

      allocate ( this%dt(this%norder_exp) )
      allocate ( this%wimp(this%nold) )
      allocate ( this%wimp_lin(this%norder_imp_lin) )
      allocate ( this%wexp(this%norder_exp) )

      this%dt(:)       = 0.0_cp
      this%wimp(:)     = 0.0_cp
      this%wimp_lin(:) = 0.0_cp
      this%wexp(:)     = 0.0_cp

      bytes_allocated = bytes_allocated+(2*this%norder_exp+this%nold+&
      &                 this%norder_imp_lin)*SIZEOF_DEF_REAL

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this)

      class(type_multistep) :: this

      deallocate( this%l_exp_calc, this%dt, this%wimp, this%wimp_lin, this%wexp )
      deallocate( this%l_imp_calc_rhs )

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine set_weights(this, lMatNext)

      class(type_multistep) :: this
      logical, intent(inout) :: lMatNext

      !-- Local variables
      real(cp) :: delta, delta_n, delta_n_1, delta_n_2
      real(cp) :: a0, a1, a2, a3, a4, b0, b1, b2, b3, c0, c1, c2, c3
      real(cp) :: gam, theta, c 
      real(cp) :: wimp_old

      wimp_old = this%wimp_lin(1)

      select case ( this%time_scheme )
         case ('CNAB2') 
            this%wimp(1)    =one
            this%wimp_lin(1)=alpha*this%dt(1)
            this%wimp_lin(2)=(1-alpha)*this%dt(1)

            this%wexp(1)=(one+half*this%dt(1)/this%dt(2))*this%dt(1)
            this%wexp(2)=-half*this%dt(1)*this%dt(1)/this%dt(2)
         case ('CNLF')
            delta = this%dt(1)/this%dt(2)
            this%wimp(1)    =(one-delta)*(one+delta)
            this%wimp(2)    =delta*delta
            this%wimp_lin(1)=half*(one+delta)/delta*this%dt(1)
            this%wimp_lin(2)=half*(one+delta)*(delta-one)/delta*this%dt(1)
            this%wimp_lin(3)=half*(one+delta)*this%dt(1)

            this%wexp(1)=(one+delta)*this%dt(1)
            this%wexp(2)=0.0_cp
         case ('MODCNAB') 
            delta = this%dt(1)/this%dt(2)
            this%wimp(1)    =one
            this%wimp(2)    =0.0_cp
            this%wimp_lin(1)=(half+1.0_cp/delta/16.0_cp)*this%dt(1)
            this%wimp_lin(2)=(7.0_cp/16.0_cp-1.0_cp/delta/16.0_cp)*this%dt(1)
            this%wimp_lin(3)=1.0_cp/16.0_cp*this%dt(1)

            this%wexp(1)=(one+half*delta)*this%dt(1)
            this%wexp(2)=-half*delta*this%dt(1)
         case ('SBDF2')
            delta = this%dt(1)/this%dt(2)
            this%wimp_lin(1)=(one+delta)/(one+two*delta)*this%dt(1)
            this%wimp_lin(2)=0.0_cp
            this%wimp(1)=(one+delta)*(one+delta)/(one+two*delta)
            this%wimp(2)=-delta*delta/(one+two*delta)

            this%wexp(1)=(one+delta)*(one+delta)*this%dt(1)/(one+two*delta)
            this%wexp(2)=-delta*(one+delta)*this%dt(1)/(one+two*delta)
         case ('SBDF3')
            delta_n   = this%dt(2)/this%dt(1)
            delta_n_1 = this%dt(3)/this%dt(1)
            a0 = one+one/(one+delta_n)+one/(one+delta_n+delta_n_1)
            a1 = (one+delta_n)*(one+delta_n+delta_n_1)/ &
            &    (delta_n*(delta_n+delta_n_1))
            a2 = -(one+delta_n+delta_n_1)/(delta_n*delta_n_1* &
            &    (one+delta_n))
            a3 = (one+delta_n)/(delta_n_1*(delta_n+delta_n_1)* &
            &     (one+delta_n+delta_n_1))
            b0 = (one+delta_n)*(one+delta_n+delta_n_1)/(delta_n* &
            &    (delta_n+delta_n_1))
            b1 = -(one+delta_n+delta_n_1)/(delta_n*delta_n_1)
            b2 = (one+delta_n)/(delta_n_1*(delta_n+delta_n_1))

            this%wimp_lin(1)=one/a0 * this%dt(1)
            this%wimp_lin(2)=0.0_cp

            this%wimp(1)=a1/a0
            this%wimp(2)=a2/a0
            this%wimp(3)=a3/a0

            this%wexp(1)=b0/a0 * this%dt(1)
            this%wexp(2)=b1/a0 * this%dt(1)
            this%wexp(3)=b2/a0 * this%dt(1)

         case ('TVB33')
            gam = (2619.0_cp-sqrt(5995281.0_cp))/345.0_cp
            theta = (3525.0_cp*gam*gam-3118.0_cp*gam-5084.0_cp)/5468.0_cp
            c = (2178.0_cp*theta-959.0_cp*gam*gam+130.0_cp*gam+726.0_cp)/4096.0_cp

            delta_n = this%dt(1)/this%dt(2)
            delta_n_1 = this%dt(2)/this%dt(3)

            a3 = - delta_n_1**3 * delta_n**2 * (3.0_cp*gam**2*delta_n+two*gam* &
            &      (one-delta_n)-one)/(one+delta_n_1)/(one+delta_n_1*(one+delta_n))
            a2 = - delta_n**2*(gam*delta_n*delta_n_1*(two-3.0_cp*gam)+(one-two*gam)*&
            &      (one+delta_n_1))/(one+delta_n)
            a1 = - (delta_n_1*(one+gam*delta_n)*(1+delta_n*(3.0_cp*gam-two))+ &
            &       delta_n*(two*gam-one)+one)/(one+delta_n_1)-theta
            a0 = theta+(one+two*gam*delta_n+delta_n_1*(one+gam*delta_n)*&
            &    (one+3.0_cp*gam*delta_n))/(one+delta_n)/               &
            &    (one+delta_n_1*(one+delta_n))

            b2 = delta_n_1**2*delta_n*(6.0_cp*gam*(one+gam*delta_n)+theta* &
            &    (3.0_cp+two*delta_n))/(6.0_cp*(one+delta_n_1))
            b1 = -gam*delta_n*(one+delta_n_1*(one+gam*delta_n))-theta/6.0_cp* &
            &    delta_n*(3.0_cp+delta_n_1*(3.0_cp+two*delta_n))
            b0 = (one+gam*delta_n)*(one+delta_n_1*(one+gam*delta_n))/(one+delta_n_1)&
            &    +theta*(one+0.5_cp*delta_n+delta_n_1*delta_n*(3.0_cp+two*delta_n)/ &
            &    (6.0_cp*(one+delta_n_1)))

            c3 = theta*delta_n_1**2*delta_n*(3.0_cp+two*delta_n)/(6.0_cp* &
            &    (one+delta_n_1))-c
            c2 = (c*(one+delta_n_1)*(one+delta_n_1*(one+delta_n))-delta_n_1**2* &
            &    delta_n**2*gam*(one-gam))/(delta_n_1**2*(one+delta_n))-theta*  &
            &    delta_n*(3.0_cp+delta_n_1*(3.0_cp+2.0_cp*delta_n))/6.0_cp
            c1 = (delta_n_1**2*delta_n*(one-gam)*(one+gam*delta_n)-c*(one+delta_n_1*&
            &    (one+delta_n)))/delta_n_1**2/delta_n+theta*(one+0.5_cp*delta_n+    &
            &    delta_n_1*delta_n*(3.0_cp+two*delta_n)/(6.0_cp*(one+delta_n_1)))
            c0 = (delta_n_1**2*delta_n*gam*(one+gam*delta_n)+c*(one+delta_n_1)) / &
            &    (delta_n_1**2*delta_n*(one+delta_n))

            this%wimp_lin(1)=c0/a0 * this%dt(1)
            this%wimp_lin(2)=c1/a0 * this%dt(1)
            this%wimp_lin(3)=c2/a0 * this%dt(1)
            this%wimp_lin(4)=c3/a0 * this%dt(1)
            !this%wimp_lin(5)=0.0_cp

            this%wimp(1)=-a1/a0
            this%wimp(2)=-a2/a0
            this%wimp(3)=-a3/a0

            this%wexp(1)=b0/a0 * this%dt(1)
            this%wexp(2)=b1/a0 * this%dt(1)
            this%wexp(3)=b2/a0 * this%dt(1)

         case ('SBDF4')
            delta_n = this%dt(1)/this%dt(2)
            delta_n_1 = this%dt(2)/this%dt(3)
            delta_n_2 = this%dt(3)/this%dt(4)

            c1 = one+delta_n_2*(one+delta_n_1)
            c2 = one+delta_n_1*(one+delta_n)
            c3 = one+delta_n_2*c2

            a4 = (one+delta_n)/(one+delta_n_2) * c2/c1/c3 *(delta_n_2**4* &
            &    delta_n_1**3*delta_n**2)
            a3 = -delta_n_1**3*delta_n**2*(one+delta_n)/(one+delta_n_1) * c3/c2
            a2 = delta_n*(delta_n/(one+delta_n)+delta_n_1*delta_n * (c3+delta_n_2)/&
            &    (one+delta_n_2))
            a1 = -one-delta_n*(one+delta_n_1*(one+delta_n)*(one+delta_n_2*c2/c1)/ &
            &    (one+delta_n_1))
            a0 = one + delta_n/(one+delta_n)+delta_n_1*delta_n/c2+delta_n_2* &
            &    delta_n_1*delta_n/c3

            b3 = -delta_n_2**3*delta_n_1**2*delta_n*(one+delta_n)/(one+delta_n_2)* &
            &    c2/c1
            b2 = delta_n_1**2*delta_n*(one+delta_n)/(one+delta_n_1)*c3
            b1 = -c2*c3 * delta_n/(one+delta_n_2)
            b0 = delta_n_1*(one+delta_n)/(one+delta_n_1) * ((one+delta_n) * &
            &    (c3+delta_n_2)+(one+delta_n_2)/delta_n_1)/c1

            this%wimp_lin(1)=one/a0 * this%dt(1)
            this%wimp_lin(2)=0.0_cp

            this%wimp(1)=-a1/a0
            this%wimp(2)=-a2/a0
            this%wimp(3)=-a3/a0
            this%wimp(4)=-a4/a0

            this%wexp(1)=b0/a0 * this%dt(1)
            this%wexp(2)=b1/a0 * this%dt(1)
            this%wexp(3)=b2/a0 * this%dt(1)
            this%wexp(4)=b3/a0 * this%dt(1)
      end select

      if ( this%wimp_lin(1) /= wimp_old ) lMatNext = .true.

   end subroutine set_weights
!------------------------------------------------------------------------------
   subroutine set_dt_array(this, dt_new, dt_min, time, n_log_file,  &
              &            n_time_step, l_new_dtNext)
      !
      ! This subroutine adjusts the time step
      !

      class(type_multistep) :: this

      !-- Input variables
      real(cp), intent(in) :: dt_new
      real(cp), intent(in) :: dt_min
      real(cp), intent(in) :: time
      integer,  intent(in) :: n_log_file
      integer,  intent(in) :: n_time_step
      logical,  intent(in) :: l_new_dtNext

      !-- Local variables
      real(cp) :: dt_old

      dt_old = this%dt(1)

      !-- First roll the dt array
      this%dt   =cshift(this%dt,shift=this%norder_exp-1)
      !-- Then overwrite the first element by the new timestep
      this%dt(1)=dt_new

      !----- Stop if time step has become too small:
      if ( dt_new < dt_min ) then
         if ( rank == 0 ) then
            write(output_unit,'(1p,/,A,ES14.4,/,A)')   &
            &    " ! Time step too small, dt=",dt_new, &
            &    " ! I thus stop the run !"
            if ( l_save_out ) then
               open(n_log_file, file=log_file, status='unknown', &
               &    position='append')
            end if
            write(n_log_file,'(1p,/,A,ES14.4,/,A)')    &
            &    " ! Time step too small, dt=",dt_new, &
            &    " ! I thus stop the run !"
            if ( l_save_out ) close(n_log_file)
         end if
         call abortRun('Stop run in steptime!')
      end if

      if ( l_new_dtNext ) then
         !------ Writing info and getting new weights:
         if ( rank == 0 ) then
            write(output_unit,'(1p,/,A,ES18.10,/,A,i9,/,A,ES15.8,/,A,ES15.8)')  &
            &    " ! Changing time step at time=",(time+this%dt(1)),            &
            &    "                 time step no=",n_time_step,                  &
            &    "                      last dt=",dt_old,                       &
            &    "                       new dt=",dt_new
            if ( l_save_out ) then
               open(n_log_file, file=log_file, status='unknown', &
               &    position='append')
            end if
            write(n_log_file,                                         &
            &    '(1p,/,A,ES18.10,/,A,i9,/,A,ES15.8,/,A,ES15.8)')     &
            &    " ! Changing time step at time=",(time+this%dt(1)),  &
            &    "                 time step no=",n_time_step,        &
            &    "                      last dt=",dt_old,             &
            &    "                       new dt=",dt_new
            if ( l_save_out ) close(n_log_file)
         end if
      end if

   end subroutine set_dt_array
!------------------------------------------------------------------------------
   subroutine set_imex_rhs(this, rhs, dfdt, lmStart, lmStop, len_rhs)
      !
      ! This subroutine assembles the right-hand-side of an IMEX scheme
      !

      class(type_multistep) :: this

      !-- Input variables:
      integer,     intent(in) :: lmStart
      integer,     intent(in) :: lmStop
      integer,     intent(in) :: len_rhs
      type(type_tarray), intent(in) :: dfdt

      !-- Output variable
      complex(cp), intent(out) :: rhs(lmStart:lmStop,len_rhs)

      !-- Local variables
      integer :: n_o, n_r, lm, startR, stopR

      !$omp parallel default(shared) private(startR, stopR,lm,n_r)
      startR=1; stopR=len_rhs
      call get_openmp_blocks(startR,stopR)
      
      do n_o=1,this%nold
         if ( n_o == 1 ) then
            do n_r=startR,stopR
               do lm=lmStart,lmStop
                  rhs(lm,n_r)=this%wimp(n_o)*dfdt%old(lm,n_r,n_o)
               end do
            end do
         else
            do n_r=startR,stopR
               do lm=lmStart,lmStop
                  rhs(lm,n_r)=rhs(lm,n_r)+this%wimp(n_o)*dfdt%old(lm,n_r,n_o)
               end do
            end do
         end if
      end do

      do n_o=1,this%norder_imp_lin-1
         do n_r=startR,stopR
            do lm=lmStart,lmStop
               rhs(lm,n_r)=rhs(lm,n_r)+this%wimp_lin(n_o+1)*dfdt%impl(lm,n_r,n_o)
            end do
         end do
      end do

      do n_o=1,this%norder_exp
         do n_r=startR,stopR
            do lm=lmStart,lmStop
               rhs(lm,n_r)=rhs(lm,n_r)+this%wexp(n_o)*dfdt%expl(lm,n_r,n_o)
            end do
         end do
      end do

      !$omp end parallel

   end subroutine set_imex_rhs
!------------------------------------------------------------------------------
   subroutine set_imex_rhs_scalar(this, rhs, dfdt)
      !
      ! This subroutine assembles the right-hand-side of an IMEX scheme
      !

      class(type_multistep) :: this

      !-- Input variables:
      type(type_tscalar), intent(in) :: dfdt

      !-- Output variable
      real(cp), intent(out) :: rhs

      !-- Local variables
      integer :: n_o

      do n_o=1,this%nold
         if ( n_o == 1 ) then
            rhs=this%wimp(n_o)*dfdt%old(n_o)
         else
            rhs=rhs+this%wimp(n_o)*dfdt%old(n_o)
         end if
      end do

      do n_o=1,this%norder_imp_lin-1
         rhs=rhs+this%wimp_lin(n_o+1)*dfdt%impl(n_o)
      end do

      do n_o=1,this%norder_exp
         rhs=rhs+this%wexp(n_o)*dfdt%expl(n_o)
      end do

   end subroutine set_imex_rhs_scalar
!------------------------------------------------------------------------------
   subroutine rotate_imex(this, dfdt, lmStart, lmStop, n_r_max)
      !
      ! This subroutine is used to roll the time arrays from one time step
      !

      class(type_multistep) :: this

      !-- Input variables:
      integer,     intent(in) :: lmStart
      integer,     intent(in) :: lmStop
      integer,     intent(in) :: n_r_max

      !-- Output variables:
      type(type_tarray), intent(inout) :: dfdt

      !-- Local variables:
      integer :: n_o, lm, n_r, startR, stopR

      !$omp parallel default(shared) private(startR, stopR,lm,n_r)
      startR=1; stopR=n_r_max
      call get_openmp_blocks(startR,stopR)

      do n_o=this%norder_exp,2,-1
         do n_r=startR,stopR
            do lm=lmStart,lmStop
               dfdt%expl(lm,n_r,n_o)=dfdt%expl(lm,n_r,n_o-1)
            end do
         end do
      end do

      do n_o=this%nold,2,-1
         do n_r=startR,stopR
            do lm=lmStart,lmStop
               dfdt%old(lm,n_r,n_o)=dfdt%old(lm,n_r,n_o-1)
            end do
         end do
      end do

      do n_o=this%norder_imp_lin-1,2,-1
         do n_r=startR,stopR
            do lm=lmStart,lmStop
               dfdt%impl(lm,n_r,n_o)=dfdt%impl(lm,n_r,n_o-1)
            end do
         end do
      end do

      !$omp end parallel

   end subroutine rotate_imex
!------------------------------------------------------------------------------
   subroutine rotate_imex_scalar(this, dfdt)
      !
      ! This subroutine is used to roll the time arrays from one time step
      !

      class(type_multistep) :: this

      !-- Output variables:
      type(type_tscalar), intent(inout) :: dfdt

      !-- Local variables:
      integer :: n_o

      do n_o=this%norder_exp,2,-1
         dfdt%expl(n_o)=dfdt%expl(n_o-1)
      end do

      do n_o=this%nold,2,-1
         dfdt%old(n_o)=dfdt%old(n_o-1)
      end do

      do n_o=this%norder_imp_lin-1,2,-1
         dfdt%impl(n_o)=dfdt%impl(n_o-1)
      end do


   end subroutine rotate_imex_scalar
!------------------------------------------------------------------------------
   subroutine bridge_with_cnab2(this)

      class(type_multistep) :: this

      !-- Local variables
      logical :: lMatNext
      character(len=8) :: old_scheme
      integer :: old_order

      if (rank == 0 ) write(output_unit,*) '! Crank-Nicolson for this time-step'

      old_order=this%norder_imp_lin
      this%norder_imp_lin=2

      old_scheme         =this%time_scheme
      this%time_scheme='CNAB2'
      call this%set_weights(lMatNext)
      !-- Since CN has only two coefficients, one has to set the remainings to zero
      this%wimp(2:size(this%wimp))=0.0_cp
      this%wimp_lin(3:size(this%wimp_lin))=0.0_cp
      this%wexp(3:size(this%wexp))=0.0_cp
      this%time_scheme   =old_scheme
      this%norder_imp_lin=old_order

   end subroutine bridge_with_cnab2
!------------------------------------------------------------------------------
   subroutine start_with_ab1(this)

      class(type_multistep) :: this

      if (rank == 0 ) write(output_unit,*) &
      &                  '! 1st order Adams-Bashforth for 1st time step'
      this%wexp(1)=this%dt(1) ! Instead of one
      this%wexp(2:this%norder_exp)=0.0_cp

   end subroutine start_with_ab1
!------------------------------------------------------------------------------
   subroutine get_time_stage(this, tlast, tstage)

      class(type_multistep) :: this
      real(cp), intent(in) :: tlast
      real(cp), intent(out) :: tstage

      tstage = tlast+this%dt(1)

   end subroutine get_time_stage
!------------------------------------------------------------------------------
end module multistep_schemes
