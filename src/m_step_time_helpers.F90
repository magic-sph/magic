!$Id$
module step_time_helpers

   use parallel_mod
   use num_param, only: dtMin, tScale
   use logic, only: l_save_out
   use output_data, only: n_time_hits, t_rst, t_graph, t_log,      &
                          t_spec, t_cmb, t_movie, t_TO, t_TOmovie, &
                          n_log_file, log_file
   use radial_data, only: nRstart, nRstop


   implicit none

   private

   public :: check_time_hits, dt_courant

contains

   subroutine check_time_hits(l_new_dt,time,dt,dt_new)
      !  +-------------+----------------+------------------------------------+
      !  |                                                                   |
      !  |  Checks whether a certain dt is required to hit a                 |
      !  |  specific output-time.                                            |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+

      !-- Output: ev. modified dt
      logical, intent(out) :: l_new_dt ! signfies change of dt !
      real(kind=8), intent(inout) :: time,dt,dt_new
       
      !-- Local variables:
      integer :: n_dt_hit
      integer, parameter :: n_dt_hit_max=10
      real(kind=8) ::  dt_hit(n_dt_hit_max) ! dt for different hit times
      integer :: n                          ! counter
      real(kind=8) ::  time_new             ! Next time step

      time_new=time+dt
      l_new_dt=.false.

      n_dt_hit=7

      do n=1,n_dt_hit
         dt_hit(n)=0.D0
      end do

      do n=1,n_time_hits
         if ( t_rst(n) > time .and. t_rst(n) < time_new ) &
              dt_hit(1)=t_rst(n)-time
         if ( t_graph(n) > time .and. t_graph(n) < time_new ) &
              dt_hit(2)=t_graph(n)-time
         if ( t_log(n) > time .and. t_log(n) < time_new ) &
              dt_hit(3)=t_log(n)-time
         if ( t_spec(n) > time .and. t_spec(n) < time_new ) &
              dt_hit(4)=t_spec(n)-time
         if ( t_cmb(n) > time .and. t_cmb(n) < time_new ) &
              dt_hit(5)=t_cmb(n)-time
         if ( t_movie(n) > time .and. t_movie(n) < time_new ) &
              dt_hit(6)=t_movie(n)-time
         if ( t_TO(n) > time .and. t_TO(n) < time_new ) &
              dt_hit(7)=t_TO(n)-time
         if ( t_TOmovie(n) > time .and. t_TOmovie(n) < time_new ) &
              dt_hit(7)=t_TOmovie(n)-time
      end do

      do n=1,n_dt_hit
         if ( dt_hit(n) /= 0.D0 .and. dt_hit(n) < dt_new ) then
            l_new_dt=.true.
            dt_new=dt_hit(n)
         end if
      end do

      if ( l_new_dt ) then
         if ( dt_new < dtMin ) dt_new=dtMin
         time_new=time+dt_new
         write(*, '(/," ! TIME STEP CHANGED TO HIT TIME:",1p,2d16.6)') &
         &     time_new*tScale,time*tScale
         if ( rank == 0 ) then
            if ( l_save_out ) then
               open(n_log_file, file=log_file, status='unknown', position='append')
               write(n_log_file, &
                    &     '(/," ! TIME STEP CHANGED TO HIT TIME:",1p,2d16.6)') &
                    &     time_new*tScale,time*tScale
               close(n_log_file)
            else
               write(n_log_file, &
                    &    '(/," ! TIME STEP CHANGED TO HIT TIME:",1p,2d16.6)') &
                    &    time_new*tScale,time*tScale
            end if
         end if
      end if

   end subroutine check_time_hits
!------------------------------------------------------------------------------
   subroutine dt_courant(dt_r,dt_h,l_new_dt,dt,dt_new,dtMax,dtrkc,dthkc)
      !---------------------------------------------------------------------------

      ! *** Check if Courant criterion based on combined
      ! *** fluid and Alfven velocity is satisfied
      ! *** Returns new value of time step dtnew

      !     dtr,dth: (output) radial/horizontal Courant time step
      !     n_time_step: (input) time step number
      !     l_new_dt: (output) flag indicating that time step is changed (=1) or not (=0)
      !     dt: (input) old time step
      !     dtnew: (output) new time step
      !     dtMin: (input) lower limit for time step (termination if dtnew < dtMin)
      !     dtMax: (input) upper limit for time step
      !     dtrkc: (input) radial Courant time step as function of radial level
      !     dthkc: (input) horizontal Courant time step as function of radial level

      !---------------------------------------------------------------------------

      !-- Input variables:
      real(kind=8), intent(in) :: dt
      real(kind=8), intent(in) :: dtMax
      real(kind=8), intent(in) :: dtrkc(nRstart:nRstop),dthkc(nRstart:nRstop)
    
      !-- Output variables:
      logical,      intent(out) :: l_new_dt
      real(kind=8), intent(out) :: dt_new
      real(kind=8), intent(out) :: dt_r,dt_h
    
      !-- Local:
      integer :: n_r
      real(kind=8) :: dt_rh,dt_2
      real(kind=8) :: dt_fac
    
      character(len=200) :: message
    
    
      dt_fac=2.D0
      dt_r  =1000.D0*dtMax
      dt_h  =dt_r
      do n_r=nRstart,nRstop
         dt_r=min(dtrkc(n_r),dt_r)
         dt_h=min(dthkc(n_r),dt_h)
      end do
      call MPI_Allreduce(MPI_IN_PLACE,dt_r,1,MPI_DOUBLE_PRECISION, &
                         MPI_min,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(MPI_IN_PLACE,dt_h,1,MPI_DOUBLE_PRECISION, &
                         MPI_min,MPI_COMM_WORLD,ierr)
    
      dt_rh=min(dt_r,dt_h)
      dt_2 =min(0.5D0*(1.D0/dt_fac+1.D0)*dt_rh,dtMax)

      if ( dt > dtMax ) then
    
         l_new_dt=.true.
         dt_new=dtMax
         write(message,'(1P," ! COURANT: dt=dtMax =",d12.4,A)') dtMax,&
              &" ! Think about changing dtMax !"
         call logWrite(message)
    
      else if ( dt > dt_rh ) then
    
         l_new_dt=.true.
         dt_new  =dt_2
         write(message,'(1P," ! COURANT: dt=",D11.4," > dt_r=",D12.4, &
              &       " and dt_h=",D12.4)') dt,dt_r,dt_h
         call logWrite(message)
    
      else if ( dt_fac*dt < dt_rh .and. dt < dtMax ) then
    
         l_new_dt=.true.
         dt_new=dt_2
         write(message,'(" ! COURANT: ",F4.1,1P,"*dt=",D11.4, &
              &     " < dt_r=",D12.4," and dt_h=",D12.4)') &
              &     dt_fac,dt_fac*dt,dt_r,dt_h
         call logWrite(message)
    
      end if
    
      if ( dt == dt_new ) l_new_dt= .false. 
       
   end subroutine dt_courant
!------------------------------------------------------------------------------
end module step_time_helpers
