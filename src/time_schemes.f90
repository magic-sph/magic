module time_schemes
   !
   ! This module defines an abstract class type_tscheme which is employed for
   ! the time advance of the code.
   !

   use iso_fortran_env, only: output_unit
   use time_array
   use logic, only: l_save_out
   use output_data, only: log_file
   use precision_mod

   implicit none

   private

   type, abstract, public :: type_tscheme

      character(len=10) :: family ! Family of the time integrator: MULTISTEP, DIRK
      integer :: nstages ! Total number of stages
      integer :: istage ! index of the current stage
      integer :: nexp ! Number of explicit terms
      integer :: nold ! Number of old state
      integer :: nimp ! Number of implicit terms
      character(len=8) :: time_scheme ! Name of the time scheme
      logical :: l_assembly
      real(cp), allocatable :: dt(:) ! Array that contains the timesteps
      real(cp), pointer :: wimp_lin(:) ! Weighting factor
      logical,  allocatable :: l_exp_calc(:) ! Array of booleans to specify the calculation of an explicit stage
      logical, allocatable :: l_imp_calc_rhs(:) ! Array of booleans to specify the calculation of an implicit stage
      real(cp) :: courfac ! Courant factor
      real(cp) :: alffac ! Courant factor for Alfven waves
      real(cp) :: intfac ! Coriolis factor

   contains

      procedure(initialize_if), deferred :: initialize
      procedure(finalize_if),  deferred :: finalize
      procedure(set_weights_if), deferred :: set_weights
      procedure(set_dt_array_if), deferred :: set_dt_array
      procedure(assemble_imex_if),  deferred :: assemble_imex
      procedure(assemble_imex_scalar_if),  deferred :: assemble_imex_scalar
      procedure(set_imex_rhs_if),  deferred :: set_imex_rhs
      procedure(set_imex_rhs_ghost_if),  deferred :: set_imex_rhs_ghost
      procedure(set_imex_rhs_scalar_if),  deferred :: set_imex_rhs_scalar
      procedure(rotate_imex_if), deferred :: rotate_imex
      procedure(rotate_imex_scalar_if), deferred :: rotate_imex_scalar
      procedure(bridge_with_cnab2_if), deferred :: bridge_with_cnab2
      procedure(start_with_ab1_if), deferred :: start_with_ab1
      procedure(get_time_stage_if), deferred :: get_time_stage
      procedure :: print_info

   end type type_tscheme

   interface

      subroutine initialize_if(this, time_scheme, courfac_nml, intfac_nml, &
                 &             alffac_nml)
         import
         class(type_tscheme) :: this
         real(cp),          intent(in)    :: courfac_nml
         real(cp),          intent(in)    :: intfac_nml
         real(cp),          intent(in)    :: alffac_nml
         character(len=72), intent(inout) :: time_scheme
      end subroutine initialize_if

      subroutine finalize_if(this)
         import
         class(type_tscheme) :: this
      end subroutine finalize_if

      subroutine set_weights_if(this, lMatNext)
         import
         class(type_tscheme) :: this
         logical, intent(inout) :: lMatNext
      end subroutine set_weights_if

      subroutine set_dt_array_if(this, dt_new, dt_min, time, n_log_file,  &
                 &                n_time_step, l_new_dtNext)
         import
         class(type_tscheme) :: this
         real(cp), intent(in) :: dt_new
         real(cp), intent(in) :: dt_min
         real(cp), intent(in) :: time
         integer,  intent(inout) :: n_log_file
         integer,  intent(in) :: n_time_step
         logical,  intent(in) :: l_new_dtNext
      end subroutine set_dt_array_if

      subroutine set_imex_rhs_if(this, rhs, dfdt)
         import
         class(type_tscheme) :: this
         type(type_tarray), intent(in) :: dfdt
         complex(cp), intent(inout) :: rhs(dfdt%llm:dfdt%ulm,dfdt%nRstart:dfdt%nRstop)
      end subroutine set_imex_rhs_if

      subroutine set_imex_rhs_ghost_if(this, rhs, dfdt, start_lm, stop_lm, ng)
         import
         class(type_tscheme) :: this
         type(type_tarray), intent(in) :: dfdt
         integer,           intent(in) :: start_lm ! Starting index
         integer,           intent(in) :: stop_lm  ! Stopping index
         integer,           intent(in) :: ng       ! Number of ghosts zones
         complex(cp), intent(out) :: rhs(dfdt%llm:dfdt%ulm,dfdt%nRstart-ng:dfdt%nRstop+ng)
      end subroutine set_imex_rhs_ghost_if

      subroutine assemble_imex_if(this, rhs, dfdt)
         import
         class(type_tscheme) :: this
         type(type_tarray), intent(in) :: dfdt
         complex(cp), intent(out) :: rhs(dfdt%llm:dfdt%ulm,dfdt%nRstart:dfdt%nRstop)
      end subroutine assemble_imex_if

      subroutine assemble_imex_scalar_if(this, rhs, dfdt)
         import
         class(type_tscheme) :: this
         type(type_tscalar), intent(in) :: dfdt
         real(cp), intent(out) :: rhs
      end subroutine assemble_imex_scalar_if

      subroutine set_imex_rhs_scalar_if(this, rhs, dfdt)
         import
         class(type_tscheme) :: this
         type(type_tscalar), intent(in) :: dfdt
         real(cp), intent(out) :: rhs
      end subroutine set_imex_rhs_scalar_if

      subroutine rotate_imex_if(this, dfdt)
         import
         class(type_tscheme) :: this
         type(type_tarray), intent(inout) :: dfdt
      end subroutine rotate_imex_if

      subroutine rotate_imex_scalar_if(this, dfdt)
         import
         class(type_tscheme) :: this
         type(type_tscalar), intent(inout) :: dfdt
      end subroutine rotate_imex_scalar_if

      subroutine bridge_with_cnab2_if(this)
         import
         class(type_tscheme) :: this
      end subroutine bridge_with_cnab2_if

      subroutine start_with_ab1_if(this)
         import
         class(type_tscheme) :: this
      end subroutine start_with_ab1_if

      subroutine get_time_stage_if(this, tlast, tstage)
         import
         class(type_tscheme) :: this
         real(cp), intent(in) :: tlast
         real(cp), intent(out) :: tstage
      end subroutine get_time_stage_if

   end interface

contains

   subroutine print_info(this, n_log_file)
      !
      ! This subroutine prints some informations about the time stepping scheme
      !

      class(type_tscheme) :: this

      !-- Input variable
      integer, intent(inout) :: n_log_file ! File unit of the log.TAG file

      !-- Local variables
      integer :: n, n_out

      do n=1,2
         if ( n == 1 ) n_out=output_unit
         if ( n == 2 ) n_out=n_log_file
         if ( n == 2 .and. l_save_out ) then
            open(newunit=n_log_file, file=log_file, status='unknown', &
            &    position='append')
            n_out=n_log_file
         end if
         write(n_out,*) ''
         write(n_out, '('' ! Time integrator   :'',1p,A10)') this%time_scheme
         write(n_out, '('' ! CFL (flow) value  :'',es12.4)') this%courfac
         write(n_out, '('' ! CFL (Alfven) value:'',es12.4)') this%alffac
         write(n_out, '('' ! CFL (Ekman) value :'',es12.4)') this%intfac
         write(n_out,*) ''
         if ( n == 2 .and. l_save_out ) close(n_log_file)
      end do

   end subroutine print_info
!------------------------------------------------------------------------------
end module time_schemes
