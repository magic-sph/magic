module mean_sd
   !
   ! This module contains a small type that simply handles two arrays (mean and SD)
   ! This type is used for time-averaged outputs (and their standard deviations).
   !

   use mem_alloc
   use precision_mod
   use parallel_mod, only: rank
   use constants, only: osq4pi
   use useful, only: cc2real
   use blocking, only: lo_map
   use communications, only: reduce_radial

   implicit none

   private

   type, public :: mean_sd_type
      real(cp), allocatable :: mean(:)
      real(cp), allocatable :: SD(:)
      integer :: nstart
      integer :: nstop
   contains
      procedure :: initialize => initialize_1D
      procedure :: compute_1D_real
      procedure :: compute_1D_complex
      generic :: compute => compute_1D_real, compute_1D_complex
      procedure :: finalize_SD => finalize_SD_1D
      procedure :: finalize => finalize_1D
   end type mean_sd_type

   type, public :: mean_sd_2D_type
      real(cp), allocatable :: mean(:,:)
      real(cp), allocatable :: SD(:,:)
      logical :: l_SD
      integer :: n_start_col
      integer :: n_stop_col
      integer :: n_start_row
      integer :: n_stop_row
   contains
      procedure :: initialize => initialize_2D
      procedure :: compute_2D_1D_input
      procedure :: compute_2D_2D_input
      procedure :: finalize_SD => finalize_SD_2D
      procedure :: finalize => finalize_2D
      generic   :: compute => compute_2D_1D_input, compute_2D_2D_input
   end type mean_sd_2D_type

contains

   subroutine initialize_1D(this, n_start, n_stop)
      !
      ! Memory allocation
      !
      class(mean_sd_type) :: this

      !-- Input variables:
      integer, intent(in) :: n_start
      integer, intent(in) :: n_stop

      this%nstart = n_start
      this%nstop = n_stop

      allocate( this%mean(n_start:n_stop), this%SD(n_start:n_stop) )
      bytes_allocated=bytes_allocated+2*(n_stop-n_start+1)*SIZEOF_DEF_REAL

      this%mean(:)=0.0_cp
      this%SD(:)  =0.0_cp

   end subroutine initialize_1D
!------------------------------------------------------------------------------
   subroutine compute_1D_real(this, input_data, n_ave, dt, totalTime)

      class(mean_sd_type)  :: this
      real(cp), intent(in) :: input_data(:)
      real(cp), intent(in) :: dt
      real(cp), intent(in) :: totalTime
      integer,  intent(in) :: n_ave

      !-- Local variable:
      real(cp) :: delta(this%nstart:this%nstop)

      if ( n_ave == 1) then
         this%mean(:)=input_data(:)
         this%SD(:)  =0.0_cp
      else
         delta(:)    =input_data(:)-this%mean(:)
         this%mean(:)=this%mean(:)+delta(:)*dt/totalTime
         this%SD(:)  =this%SD(:)+dt*delta(:)*(input_data(:)-this%mean(:))
      end if

   end subroutine compute_1D_real
!------------------------------------------------------------------------------
   subroutine compute_1D_complex(this, scal_lm, llm, ulm, n_r_max, n_ave, dt, totalTime)
      !
      ! This computes the mean and the spatial-variance of a scalar quantity
      !

      !-- Input variables
      class(mean_sd_type)  :: this
      integer,     intent(in) :: llm
      integer,     intent(in) :: ulm
      integer,     intent(in) :: n_r_max
      complex(cp), intent(in) :: scal_lm(llm:ulm,n_r_max)
      real(cp),    intent(in) :: dt
      real(cp),    intent(in) :: totalTime
      integer,     intent(in) :: n_ave

      !-- Local variables
      real(cp) :: scal_square(n_r_max), scal_square_glob(n_r_max)
      real(cp) :: scal(n_r_max), SD(n_r_max)
      integer :: lm, m, n_r, lm00

      !-- Compute the square of the scalar quantity
      scal_square(:)=0.0_cp
      do n_r=1,n_r_max
         do lm=llm,ulm
            m=lo_map%lm2m(lm)
            scal_square(n_r)=scal_square(n_r)+cc2real(scal_lm(lm,n_r), m)
         end do
      end do

      !-- Gather data on rank==0
      call reduce_radial(scal_square, scal_square_glob, 0)

      lm00 = lo_map%lm2(0,0) ! Spherically-symmetric mode
      if ( llm <= lm00 .and. ulm >= lm00 ) then
         scal_square_glob(:)=scal_square_glob(:)*osq4pi*osq4pi
         scal(:)=real(scal_lm(lm00,:))*osq4pi
         SD(:)  =scal_square_glob(:)-scal(:)*scal(:)

         if ( n_ave == 1) then
            this%mean(:)=scal(:)
            this%SD(:)  =SD(:)
         else
            !-- Time-average of mean and STD
            this%mean(:)=this%mean(:) + (scal(:)-this%mean(:))*dt/totalTime
            this%SD(:)  =this%SD(:) + (SD(:)-this%SD(:))*dt/totalTime
         end if
      end if

   end subroutine compute_1D_complex
!------------------------------------------------------------------------------
   subroutine finalize_SD_1D(this,totalTime)
      !
      ! Finish computation of standard-deviation
      !
      class(mean_sd_type)  :: this
      real(cp), intent(in) :: totalTime

      this%SD(:)=sqrt(this%SD(:)/totalTime)

   end subroutine finalize_SD_1D
!------------------------------------------------------------------------------
   subroutine finalize_1D(this)
      !
      ! Memory deallocation
      !
      class(mean_sd_type) :: this

      deallocate( this%mean, this%SD )

   end subroutine finalize_1D
!------------------------------------------------------------------------------
   subroutine initialize_2D(this, n_start_row, n_stop_row, n_start_col, &
              &             n_stop_col, l_SD)
      !
      ! Memory allocation
      !
      class(mean_sd_2D_type) :: this

      !-- Input variables:
      integer, intent(in) :: n_start_row
      integer, intent(in) :: n_stop_row
      integer, intent(in) :: n_start_col
      integer, intent(in) :: n_stop_col
      logical, intent(in) :: l_SD

      this%l_SD=l_SD
      this%n_start_col=n_start_col
      this%n_stop_col =n_stop_col
      this%n_start_row=n_start_row
      this%n_stop_row =n_stop_row

      allocate( this%mean(n_start_row:n_stop_row,n_start_col:n_stop_col) )
      bytes_allocated=bytes_allocated+(n_stop_row-n_start_row+1)* &
      &               (n_stop_col-n_start_col+1)*SIZEOF_DEF_REAL
      this%mean(:,:)=0.0_cp

      if ( l_SD ) then
         allocate( this%SD(n_start_row:n_stop_row,n_start_col:n_stop_col) )
         bytes_allocated=bytes_allocated+(n_stop_row-n_start_row+1)* &
         &               (n_stop_col-n_start_col+1)*SIZEOF_DEF_REAL
         this%SD(:,:)=0.0_cp
      end if

   end subroutine initialize_2D
!------------------------------------------------------------------------------
   subroutine compute_2D_1D_input(this, input_data, n_ave, dt, totalTime, ind)

      class(mean_sd_2D_type) :: this
      real(cp), intent(in)   :: input_data(:)
      real(cp), intent(in)   :: dt
      real(cp), intent(in)   :: totalTime
      integer,  intent(in)   :: n_ave, ind

      !-- Local variable:
      integer :: n_input
      real(cp), allocatable :: delta(:)

      n_input = size(input_data)
      allocate( delta(n_input) )

      if ( n_ave == 1 ) then
         this%mean(ind,:)=input_data(:)
         if ( this%l_SD ) then
            this%SD(ind,:)=0.0_cp
         end if
      else
         delta(:)        =input_data(:)-this%mean(ind,:)
         this%mean(ind,:)=this%mean(ind,:)+delta(:)*dt/totalTime
         if ( this%l_SD ) then
            this%SD(ind,:)  =this%SD(ind,:)+dt*delta(:)*(input_data(:)-this%mean(ind,:))
         end if
      end if

      deallocate( delta )

   end subroutine compute_2D_1D_input
!------------------------------------------------------------------------------
   subroutine compute_2D_2D_input(this, input_data, n_ave, dt, totalTime)

      class(mean_sd_2D_type) :: this
      real(cp), intent(in)   :: input_data(this%n_start_row:this%n_stop_row, &
                                &          this%n_start_col:this%n_stop_col)
      real(cp), intent(in)   :: dt
      real(cp), intent(in)   :: totalTime
      integer,  intent(in)   :: n_ave

      !-- Local variables
      real(cp) :: delta
      integer  :: n,m

      if ( n_ave == 1 ) then
         this%mean(:,:)=input_data(:,:)
         if ( this%l_SD ) then
            this%SD(:,:)  =0.0_cp
         end if
      else
         do m=this%n_start_col,this%n_stop_col
            do n=this%n_start_row,this%n_stop_row
               delta         =input_data(n,m)-this%mean(n,m)
               this%mean(n,m)=this%mean(n,m)+delta*dt/totalTime
               if ( this%l_SD ) then
                  this%SD(n,m)=this%SD(n,m)+dt*delta*(input_data(n,m)- &
                  &            this%mean(n,m))
               end if
            end do
         end do
      end if

   end subroutine compute_2D_2D_input
!------------------------------------------------------------------------------
   subroutine finalize_SD_2D(this,totalTime)
      !
      ! Finish computation of standard-deviation
      !
      class(mean_sd_2D_type) :: this
      real(cp), intent(in)   :: totalTime

      this%SD(:,:)=sqrt(this%SD(:,:)/totalTime)

   end subroutine finalize_SD_2D
!------------------------------------------------------------------------------
   subroutine finalize_2D(this)
      !
      ! Memory deallocation
      !
      class(mean_sd_2D_type) :: this

      deallocate( this%mean )
      if ( this%l_SD ) deallocate( this%SD )

   end subroutine finalize_2D
!------------------------------------------------------------------------------
end module mean_sd
