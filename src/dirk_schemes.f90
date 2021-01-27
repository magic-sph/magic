module dirk_schemes
   !
   ! This module defines the type_dirk which inherits from the
   ! abstract type_tscheme. It actually implements all the routine required
   ! to time-advance an diagonally implicit Runge-Kutta scheme. It makes use
   ! of Butcher tables to construct the right-hand-sides.
   !

   use iso_fortran_env, only: output_unit
   use precision_mod
   use parallel_mod
   use num_param, only: alpha
   use constants, only: one, half, two, zero
   use mem_alloc, only: bytes_allocated
   use useful, only: abortRun
   use logic, only: l_save_out
   use output_data, only: log_file
   use time_schemes, only: type_tscheme
   use time_array

   implicit none

   private

   type, public, extends(type_tscheme) :: type_dirk
      real(cp), allocatable :: butcher_imp(:,:) ! Implicit Butcher table
      real(cp), allocatable :: butcher_exp(:,:) ! Explicit Butcher table
      real(cp), allocatable :: butcher_c(:)     ! Stage time Butcher vector
      real(cp), allocatable :: butcher_ass_imp(:) ! Implicit Assembly stage
      real(cp), allocatable :: butcher_ass_exp(:) ! Explicit Assembly stage
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: set_weights
      procedure :: set_dt_array
      procedure :: set_imex_rhs
      procedure :: set_imex_rhs_ghost
      procedure :: set_imex_rhs_scalar
      procedure :: rotate_imex
      procedure :: rotate_imex_scalar
      procedure :: bridge_with_cnab2
      procedure :: start_with_ab1
      procedure :: get_time_stage
      procedure :: assemble_imex
      procedure :: assemble_imex_scalar
   end type type_dirk

contains

   subroutine initialize(this, time_scheme, courfac_nml, intfac_nml, alffac_nml)

      class(type_dirk) :: this

      !-- Input/output variables
      real(cp),          intent(in) :: courfac_nml
      real(cp),          intent(in) :: intfac_nml
      real(cp),          intent(in) :: alffac_nml
      character(len=72), intent(inout) :: time_scheme

      !-- Local variables
      integer :: sizet
      real(cp) :: courfac_loc, intfac_loc, alffac_loc

      allocate ( this%dt(1) )
      this%dt(:)=0.0_cp
      allocate ( this%wimp_lin(1) )
      this%wimp_lin(1)=0.0_cp

      this%family='DIRK'

      this%nold = 1
      if ( index(time_scheme, 'ARS222') /= 0 ) then
         this%time_scheme = 'ARS222'
         this%nimp = 2
         this%nexp = 2
         this%nstages = 2
         this%istage = 1
         courfac_loc = 1.35_cp
         alffac_loc  = 0.54_cp
         intfac_loc  = 0.28_cp
         this%l_assembly = .false.
      else if ( index(time_scheme, 'ARS232') /= 0 ) then
         this%time_scheme = 'ARS232'
         this%nimp = 3
         this%nexp = 3
         this%nstages = 3
         this%istage = 1
         courfac_loc = 1.35_cp
         alffac_loc  = 1.0_cp
         intfac_loc  = 0.28_cp
         this%l_assembly = .true.
      else if ( index(time_scheme, 'ARS233') /= 0 ) then
         this%time_scheme = 'ARS233'
         this%nimp = 3
         this%nexp = 3
         this%nstages = 3
         this%istage = 1
         courfac_loc = 1.35_cp
         alffac_loc  = 1.0_cp
         intfac_loc  = 0.28_cp
         this%l_assembly = .true.
      else if ( index(time_scheme, 'ARS343') /= 0 ) then
         this%time_scheme = 'ARS343'
         this%nimp = 4
         this%nexp = 4
         this%nstages = 4
         this%istage = 1
         courfac_loc = 0.5_cp
         alffac_loc  = 0.2_cp
         intfac_loc  = 1.0_cp
         this%l_assembly = .true.
      else if ( index(time_scheme, 'MARS343') /= 0 ) then
         this%time_scheme = 'MARS343'
         this%nimp = 4
         this%nexp = 4
         this%nstages = 4
         this%istage = 1
         courfac_loc = 0.5_cp
         alffac_loc  = 0.2_cp
         intfac_loc  = 1.0_cp
         this%l_assembly = .true.
      else if ( index(time_scheme, 'CB3') /= 0 ) then
         this%time_scheme = 'CB3'
         this%nimp = 4
         this%nexp = 4
         this%nstages = 4
         this%istage = 1
         courfac_loc = 0.5_cp
         alffac_loc  = 0.2_cp
         intfac_loc  = 1.0_cp
         this%l_assembly = .true.
      else if ( index(time_scheme, 'KC343') /= 0 ) then
         this%time_scheme = 'KC343'
         this%nimp = 4
         this%nexp = 4
         this%nstages = 4
         this%istage = 1
         courfac_loc = 0.5_cp
         alffac_loc  = 1.0_cp
         intfac_loc  = 1.0_cp
         this%l_assembly = .true.
      else if ( index(time_scheme, 'KC564') /= 0 ) then
         this%time_scheme = 'KC564'
         this%nimp = 6
         this%nexp = 6
         this%nstages = 6
         this%istage = 1
         courfac_loc = 0.5_cp
         alffac_loc  = 0.2_cp
         intfac_loc  = 1.0_cp
         this%l_assembly = .true.
      else if ( index(time_scheme, 'KC674') /= 0 ) then
         this%time_scheme = 'KC674'
         this%nimp = 7
         this%nexp = 7
         this%nstages = 7
         this%istage = 1
         courfac_loc = 0.5_cp
         alffac_loc  = 0.2_cp
         intfac_loc  = 1.0_cp
         this%l_assembly = .true.
      else if ( index(time_scheme, 'KC785') /= 0 ) then
         this%time_scheme = 'KC785'
         this%nimp = 8
         this%nexp = 8
         this%nstages = 8
         this%istage = 1
         courfac_loc = 0.5_cp
         alffac_loc  = 0.2_cp
         intfac_loc  = 1.0_cp
         this%l_assembly = .true.
      else if ( index(time_scheme, 'BHR553') /= 0 ) then
         this%time_scheme = 'BHR553'
         this%nimp = 5
         this%nexp = 5
         this%nstages = 5
         this%istage = 1
         courfac_loc = 0.5_cp
         alffac_loc  = 0.2_cp
         intfac_loc  = 1.0_cp
         this%l_assembly = .true.
      else if ( index(time_scheme, 'LZ232') /= 0 ) then
         this%time_scheme = 'LZ232'
         this%nimp = 2
         this%nexp = 2
         this%nstages = 2
         this%istage = 1
         courfac_loc = 1.25_cp
         alffac_loc  = 0.5_cp
         intfac_loc  = 0.3_cp
         this%l_assembly = .false.
      else if ( index(time_scheme, 'CK232') /= 0 ) then
         this%time_scheme = 'CK232'
         this%nimp = 2
         this%nexp = 2
         this%nstages = 2
         this%istage = 1
         courfac_loc = 1.25_cp
         alffac_loc  = 0.5_cp
         intfac_loc  = 0.3_cp
         this%l_assembly = .false.
      else if ( index(time_scheme, 'ARS443') /= 0 ) then
         this%time_scheme = 'ARS443'
         this%nimp = 4
         this%nexp = 4
         this%nstages = 4
         this%istage = 1
         courfac_loc = 1.0_cp
         alffac_loc  = 0.4_cp
         intfac_loc  = 0.35_cp
         this%l_assembly = .false.
      else if ( index(time_scheme, 'DBM453') /= 0 ) then
         this%time_scheme = 'DBM453'
         this%nimp = 5
         this%nexp = 5
         this%nstages = 5
         this%istage = 1
         courfac_loc = 0.5_cp
         alffac_loc  = 0.2_cp
         intfac_loc  = 1.0_cp
         this%l_assembly = .true.
      else if ( index(time_scheme, 'LZ453') /= 0 ) then
         this%time_scheme = 'LZ453'
         this%nimp = 4
         this%nexp = 4
         this%nstages = 4
         this%istage = 1
         courfac_loc = 1.15_cp
         alffac_loc  = 0.45_cp
         intfac_loc  = 0.3_cp
         this%l_assembly = .false.
      else if ( index(time_scheme, 'BPR353') /= 0 ) then
         this%time_scheme = 'BPR353'
         this%nimp = 4
         this%nexp = 4
         this%nstages = 4
         this%istage = 1
         courfac_loc = 0.8_cp
         alffac_loc  = 0.35_cp
         intfac_loc  = 0.46_cp
         this%l_assembly = .false.
      else if ( index(time_scheme, 'PC2') /= 0 ) then
         this%time_scheme = 'PC2'
         this%nimp = 3
         this%nexp = 3
         this%nstages = 3
         this%istage = 1
         courfac_loc = 0.7_cp
         alffac_loc  = 0.28_cp
         intfac_loc  = 0.5_cp
         this%l_assembly = .false.
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

      if ( .not. this%l_assembly ) then
         sizet = this%nstages+1
      else
         sizet = this%nstages
      end if
      allocate( this%butcher_imp(sizet, sizet), this%butcher_exp(sizet, sizet) )
      this%butcher_imp(:,:)=0.0_cp
      this%butcher_exp(:,:)=0.0_cp
      bytes_allocated=bytes_allocated+2*sizet*sizet*SIZEOF_DEF_REAL

      allocate( this%butcher_ass_imp(sizet), this%butcher_ass_exp(sizet) )
      this%butcher_ass_imp(:)=0.0_cp
      this%butcher_ass_exp(:)=0.0_cp
      bytes_allocated=bytes_allocated+2*sizet*SIZEOF_DEF_REAL

      allocate( this%l_exp_calc(this%nstages) )
      allocate( this%l_imp_calc_rhs(this%nstages) )
      this%l_exp_calc(:) = .true.
      this%l_imp_calc_rhs(:) = .true.
      bytes_allocated=bytes_allocated+2*this%nstages*SIZEOF_LOGICAL

      allocate( this%butcher_c(sizet) )
      bytes_allocated=bytes_allocated+sizet*SIZEOF_LOGICAL

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this)

      class(type_dirk) :: this

      deallocate( this%butcher_ass_exp, this%butcher_ass_imp )
      deallocate( this%dt, this%wimp_lin, this%butcher_exp, this%butcher_imp )
      deallocate( this%l_exp_calc, this%l_imp_calc_rhs, this%butcher_c )

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine set_weights(this, lMatNext)

      class(type_dirk) :: this
      logical, intent(inout) :: lMatNext

      !-- Local variables
      real(cp) :: wimp_old, del, gam, b1, b2, b3, b4, c3, c4
      real(cp) :: aa31, aa32, aa41, aa43, aa51, aa52, aa53, aa54, a31, a32
      real(cp) :: a41, a42, a43
      wimp_old = this%wimp_lin(1)

      select case ( this%time_scheme )
         case ('ARS222')
            gam = half * (two-sqrt(two))
            del = one-half/gam
            this%wimp_lin(1) = gam
            this%butcher_imp(:,:) = reshape([  0.0_cp,  0.0_cp, 0.0_cp,  &
                                    &          0.0_cp,  gam   , 0.0_cp,  &
                                    &          0.0_cp, one-gam, gam    ],&
                                    &          [3,3],order=[2,1])
            this%butcher_ass_imp(:) = [0.0_cp, one-gam, gam]
            this%butcher_exp(:,:) = reshape([  0.0_cp,  0.0_cp, 0.0_cp,  &
                                    &             gam,  0.0_cp, 0.0_cp,  &
                                    &             del, one-del, 0.0_cp], &
                                    &          [3,3],order=[2,1])
            this%butcher_ass_exp(:) = [del, one-del, 0.0_cp]
            this%l_imp_calc_rhs(1)=.false.
            this%butcher_c(:) = [gam, one]
         case ('ARS343')
            gam = 0.435866521508459_cp
            b1 = -1.5_cp*gam*gam+4.0_cp*gam-0.25_cp
            b2 = 1.5_cp*gam*gam-5.0_cp*gam+1.25_cp
            this%wimp_lin(1) = gam
            this%butcher_imp(:,:) = reshape([                           &
            &                  0.0_cp,         0.0_cp, 0.0_cp,  0.0_cp, &
            &                  0.0_cp,            gam, 0.0_cp,  0.0_cp, &
            &                  0.0_cp, half*(one-gam),    gam,  0.0_cp, &
            &                  0.0_cp,             b1,     b2,    gam], &
            &                               [4,4], order=[2,1])
            this%butcher_ass_imp(:) = [0.0_cp, b1, b2, gam]
            this%butcher_ass_exp(:) = this%butcher_ass_imp(:)
            this%butcher_exp(:,:) = reshape([                                           &
            &                0.0_cp,               0.0_cp,               0.0_cp,0.0_cp, &
            &                   gam,               0.0_cp,               0.0_cp,0.0_cp, &
            & 0.3212788860286278_cp,0.3966543747256017_cp,               0.0_cp,0.0_cp, &
            &-0.1058582960718797_cp,0.5529291480359398_cp,0.5529291480359398_cp,0.0_cp],&
            &                               [4,4], order=[2,1])
            this%l_imp_calc_rhs(1)=.false.
            this%butcher_c(:) = [gam, half*(gam+one), one]
         case ('BHR553')
            gam = 0.435866521508482_cp
            aa31 = gam
            aa32 = aa31
            aa41 = -475883375220285986033264.0_cp/594112726933437845704163.0_cp
            aa43 = 1866233449822026827708736.0_cp/594112726933437845704163.0_cp
            aa51 = 62828845818073169585635881686091391737610308247.0_cp/ &
            &      176112910684412105319781630311686343715753056000.0_cp
            aa52 = -302987763081184622639300143137943089.0_cp/ &
            &       1535359944203293318639180129368156500.0_cp
            aa53 = 262315887293043739337088563996093207.0_cp/  &
            &      297427554730376353252081786906492000.0_cp
            aa54 = -987618231894176581438124717087.0_cp/ &
            &       23877337660202969319526901856000.0_cp
            a31 = gam
            a32 = -31733082319927313.0_cp/455705377221960889379854647102.0_cp
            a41 = -3012378541084922027361996761794919360516301377809610.0_cp/ &
            &     45123394056585269977907753045030512597955897345819349.0_cp
            a42 = -62865589297807153294268.0_cp/ &
            &      102559673441610672305587327019095047.0_cp
            a43 = 418769796920855299603146267001414900945214277000.0_cp/ &
            &     212454360385257708555954598099874818603217167139.0_cp
            b1 = 487698502336740678603511.0_cp/1181159636928185920260208.0_cp
            b2 = 0.0_cp
            b3 = -aa52
            b4 = -105235928335100616072938218863.0_cp/ &
            &      2282554452064661756575727198000.0_cp
            c3 = 902905985686.0_cp/1035759735069.0_cp
            c4 = 2684624.0_cp/1147171.0_cp

            this%wimp_lin(1) = gam

            this%butcher_imp(:,:) = reshape([ 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, &
            &                      gam, gam, 0.0_cp, 0.0_cp, 0.0_cp,                  &
            &                      a31, a32, gam, 0.0_cp, 0.0_cp,                     &
            &                      a41, a42, a43, gam, 0.0_cp,                        &
            &                      b1, 0.0_cp, b3, b4, gam], [5,5], order=[2,1])

            this%butcher_exp(:,:) = reshape([0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, &
            &                        two*gam, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,        &
            &                        aa31, aa32, 0.0_cp, 0.0_cp, 0.0_cp,             &
            &                        aa41, 0.0_cp, aa43, 0.0_cp, 0.0_cp,             &
            &                        aa51, aa52, aa53, aa54, 0.0_cp], [5,5], order=[2,1])

            this%butcher_ass_imp(:) = [b1, 0.0_cp, b3, b4, gam]
            this%butcher_ass_exp(:) = this%butcher_ass_imp(:)
            this%butcher_c(:) = [two*gam, c3, c4, one]
         case ('MARS343')
            gam = 0.435866521508458_cp
            b1 =  1.20849664917601_cp
            b2 = -0.64436317068446_cp
            this%wimp_lin(1) = gam
            this%butcher_imp(:,:) = reshape([                           &
                            &  0.0_cp,         0.0_cp, 0.0_cp,  0.0_cp, &
                            &  0.0_cp,            gam, 0.0_cp,  0.0_cp, &
                            &  0.0_cp, half*(one-gam),    gam,  0.0_cp, &
                            &  0.0_cp,             b1,     b2,    gam], &
                            &          [4,4],order=[2,1])
            this%butcher_ass_imp(:) = [0.0_cp, b1, b2, gam]
            this%butcher_ass_exp(:) = this%butcher_ass_imp(:)
            this%butcher_exp(:,:) = reshape([                                    &
                   &          0.0_cp,          0.0_cp,         0.0_cp,  0.0_cp,  &
                   & 0.877173304301691_cp,     0.0_cp,         0.0_cp,  0.0_cp,  &
                   & 0.535396540307354_cp, 0.182536720446875_cp, 0.0_cp,  0.0_cp,&
                   & 0.63041255815287_cp, -0.83193390106308_cp, 1.20152134291021_cp,&
                   & 0.0_cp], [4,4], order=[2,1])
            this%l_imp_calc_rhs(1)=.false.
            this%butcher_c(:) = [gam, half*(gam+one), one]
         case ('CB3')
            gam = half+sqrt(3.0_cp)/6.0_cp
            b1 = -sqrt(3.0_cp)/3.0_cp
            b2 = -0.64436317068446_cp
            this%wimp_lin(1) = gam
            this%butcher_imp(:,:) = reshape([                           &
                            &  0.0_cp,         0.0_cp, 0.0_cp,  0.0_cp, &
                            &  0.0_cp,            gam, 0.0_cp,  0.0_cp, &
                            &  0.0_cp,             b1,    gam,  0.0_cp, &
                            &  0.0_cp,         0.0_cp, 0.0_cp,    gam], &
                            &          [4,4],order=[2,1])
            this%butcher_ass_imp(:) = [0.0_cp, 0.0_cp, half, half]
            this%butcher_ass_exp(:) = this%butcher_ass_imp(:)
            this%butcher_exp(:,:) = reshape([                               &
                   &    0.0_cp,                   0.0_cp, 0.0_cp,  0.0_cp,  &
                   &       gam,                   0.0_cp, 0.0_cp,  0.0_cp,  &
                   &    0.0_cp, half-sqrt(3.0_cp)/6.0_cp, 0.0_cp,  0.0_cp,  &
                   &    0.0_cp,                   0.0_cp,    gam,  0.0_cp], &
                   &    [4,4], order=[2,1])
            this%l_imp_calc_rhs(1)=.false.
            this%butcher_c(:) = [gam, half-sqrt(3.0_cp)/6.0_cp, gam]
         case ('ARS233')
            gam = (3.0_cp+sqrt(3.0_cp))/6.0_cp
            this%wimp_lin(1) = gam
            this%butcher_imp(:,:) = reshape([  0.0_cp,  0.0_cp,    0.0_cp,  &
                                    &          0.0_cp,  gam   ,    0.0_cp,  &
                                    &          0.0_cp, one-two*gam,   gam], &
                                    &          [3,3], order=[2,1])
            this%butcher_ass_imp(:) = [0.0_cp, half, half]
            this%butcher_exp(:,:) = reshape([  0.0_cp,        0.0_cp, 0.0_cp,  &
                                    &             gam,        0.0_cp, 0.0_cp,  &
                                    &         gam-one, two*(one-gam), 0.0_cp], &
                                    &          [3,3], order=[2,1])
            this%butcher_ass_exp(:) = [0.0_cp, half, half]
            this%l_imp_calc_rhs(1)=.false.
            this%butcher_c(:) = [gam, one-gam]
         case ('ARS232')
            gam = half * (two-sqrt(two))
            del = -two*sqrt(two)/3.0_cp
            this%wimp_lin(1) = gam
            this%butcher_imp(:,:) = reshape([  0.0_cp,  0.0_cp, 0.0_cp,  &
                                    &          0.0_cp,  gam   , 0.0_cp,  &
                                    &          0.0_cp, one-gam,    gam], &
                                    &          [3,3], order=[2,1])
            this%butcher_ass_imp(:) = [0.0_cp, one-gam, gam]
            this%butcher_exp(:,:) = reshape([  0.0_cp,  0.0_cp, 0.0_cp,  &
                                    &             gam,  0.0_cp, 0.0_cp,  &
                                    &             del, one-del, 0.0_cp], &
                                    &          [3,3], order=[2,1])
            this%butcher_ass_exp(:) = [0.0_cp, one-gam, gam]
            this%l_imp_calc_rhs(1)=.false.
            this%butcher_c(:) = [gam, one]
         case ('KC343')
            gam = 1767732205903.0_cp/4055673282236.0_cp
            this%wimp_lin(1) = gam

            this%butcher_imp(:,:) = reshape([0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,        &
            &    gam, gam, 0.0_cp, 0.0_cp, 2746238789719.0_cp/10658868560708.0_cp,  &
            &    -640167445237.0_cp/6845629431997.0_cp,  gam, 0.0_cp,               &
            &    1471266399579.0_cp/7840856788654.0_cp, -4482444167858.0_cp/        &
            &    7529755066697.0_cp, 11266239266428.0_cp/11593286722821.0_cp, gam], &
            &                     [4,4],  order=[2,1])

            this%butcher_ass_imp(:) = [1471266399579.0_cp/7840856788654.0_cp,  &
            &                          -4482444167858.0_cp/7529755066697.0_cp, &
            &                          11266239266428.0_cp/11593286722821.0_cp,&
            &                          gam]

            this%butcher_exp(:,:) = reshape([0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,    &
            &    1767732205903.0_cp/2027836641118.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, &
            &    5535828885825.0_cp/10492691773637.0_cp, 788022342437.0_cp/     &
            &    10882634858940.0_cp,0.0_cp, 0.0_cp, 6485989280629.0_cp/        &
            &    16251701735622.0_cp,-4246266847089.0_cp/9704473918619.0_cp,    &
            &    10755448449292.0_cp/10357097424841.0_cp, 0.0_cp], [4,4],       &
            &    order=[2,1])
            this%butcher_ass_exp(:)=this%butcher_ass_imp(:)
            this%butcher_c(:)=[two*gam, 0.6_cp, one]
         case ('KC564')
            this%wimp_lin(1) = 0.25_cp

            this%butcher_imp(:,:) = reshape([0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,   &
            &    0.0_cp, 0.0_cp, 0.25_cp, 0.25_cp, 0.0_cp, 0.0_cp, 0.0_cp,     &
            &    0.0_cp, 8611.0_cp/62500.0_cp, -1743.0_cp/31250.0_cp, 0.25_cp, &
            &    0.0_cp, 0.0_cp, 0.0_cp, 5012029.0_cp/34652500.0_cp,           &
            &    -654441.0_cp/2922500.0_cp, 174375.0_cp/388108.0_cp, 0.25_cp,  &
            &    0.0_cp, 0.0_cp, 15267082809.0_cp/155376265600.0_cp,           &
            &    -71443401.0_cp/120774400.0_cp, 730878875.0_cp/902184768.0_cp, &
            &    2285395.0_cp/8070912.0_cp, 0.25_cp, 0.0_cp, 82889.0_cp/       &
            &    524892.0_cp, 0.0_cp, 15625.0_cp/83664.0_cp, 69875.0_cp/       &
            &    102672.0_cp, -2260.0_cp/8211.0_cp, 0.25_cp], [6,6], order=[2,1])

            this%butcher_ass_imp(:) = [82889.0_cp/524892.0_cp, 0.0_cp, &
            &                          15625.0_cp/83664.0_cp,          &
            &                          69875.0_cp/102672.0_cp,         &
            &                          -2260.0_cp/8211.0_cp, 0.25_cp]

            this%butcher_exp(:,:) = reshape([ 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,    &
            &    0.0_cp, 0.0_cp, 0.5_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, &
            &    13861.0_cp/62500.0_cp,6889.0_cp/62500.0_cp, 0.0_cp, 0.0_cp,     &
            &    0.0_cp, 0.0_cp, -116923316275.0_cp/2393684061468.0_cp,          &
            &    -2731218467317.0_cp/15368042101831.0_cp, 9408046702089.0_cp/    &
            &    11113171139209.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, -451086348788.0_cp/&
            &    2902428689909.0_cp, -2682348792572.0_cp/7519795681897.0_cp,     &
            &    12662868775082.0_cp/11960479115383.0_cp, 3355817975965.0_cp/    &
            &    11060851509271.0_cp, 0.0_cp, 0.0_cp, 647845179188.0_cp/         &
            &    3216320057751.0_cp, 73281519250.0_cp/8382639484533.0_cp,        &
            &    552539513391.0_cp/3454668386233.0_cp, 3354512671639.0_cp/       &
            &    8306763924573.0_cp, 4040.0_cp/17871.0_cp, 0.0_cp], [6,6],       &
            &    order=[2,1])

            this%butcher_ass_exp(:)=this%butcher_ass_imp(:)
            this%butcher_c(:)=[half, 0.332_cp, 0.62_cp, 0.85_cp, one]

         case ('KC674')
            this%wimp_lin(1) = 0.1235_cp

            this%butcher_imp(:,:) = reshape([0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,        &
            &    0.0_cp, 0.0_cp, 0.0_cp, 0.1235_cp, 0.1235_cp, 0.0_cp, 0.0_cp,      &
            &    0.0_cp, 0.0_cp, 0.0_cp, 624185399699.0_cp/4186980696204.0_cp,      &
            &    624185399699.0_cp/4186980696204.0_cp, 0.1235_cp, 0.0_cp, 0.0_cp,   &
            &    0.0_cp, 0.0_cp, 1258591069120.0_cp/10082082980243.0_cp,            &
            &    1258591069120.0_cp/10082082980243.0_cp, -322722984531.0_cp/        &
            &    8455138723562.0_cp, 0.1235_cp, 0.0_cp, 0.0_cp, 0.0_cp,             &
            &    -436103496990.0_cp/5971407786587.0_cp, -436103496990.0_cp/         &
            &    5971407786587.0_cp, -2689175662187.0_cp/11046760208243.0_cp,       &
            &    4431412449334.0_cp/12995360898505.0_cp, 0.1235_cp, 0.0_cp, 0.0_cp, &
            &    -2207373168298.0_cp/14430576638973.0_cp, -2207373168298.0_cp/      &
            &    14430576638973.0_cp, 242511121179.0_cp/3358618340039.0_cp,         &
            &    3145666661981.0_cp/7780404714551.0_cp,  5882073923981.0_cp/        &
            &    14490790706663.0_cp, 0.1235_cp, 0.0_cp, 0.0_cp, 0.0_cp,            &
            &    9164257142617.0_cp/17756377923965.0_cp, -10812980402763.0_cp/      &
            &    74029279521829.0_cp, 1335994250573.0_cp/5691609445217.0_cp,        &
            &    2273837961795.0_cp/8368240463276.0_cp, 0.1235_cp], [7,7], order=[2,1])

            this%butcher_ass_imp(:) = [0.0_cp, 0.0_cp, 9164257142617.0_cp/        &
            &                          17756377923965.0_cp, -10812980402763.0_cp/ &
            &                          74029279521829.0_cp, 1335994250573.0_cp/   &
            &                          5691609445217.0_cp,  2273837961795.0_cp/   &
            &                          8368240463276.0_cp, 0.1235_cp]

            this%butcher_exp(:,:) = reshape([ 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,      &
            &    0.0_cp, 0.0_cp, 0.0_cp, 0.247_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, &
            &    0.0_cp, 0.0_cp, 0.06175_cp, 2694949928731.0_cp/7487940209513.0_cp,&
            &    0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 464650059369.0_cp/        &
            &    8764239774964.0_cp, 878889893998.0_cp/2444806327765.0_cp,         &
            &    -952945855348.0_cp/12294611323341.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,   &
            &    0.0_cp, 476636172619.0_cp/8159180917465.0_cp, -1271469283451.0_cp/&
            &    7793814740893.0_cp, -859560642026.0_cp/4356155882851.0_cp,        &
            &    1723805262919.0_cp/4571918432560.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,    &
            &    6338158500785.0_cp/11769362343261.0_cp, -4970555480458.0_cp/      &
            &    10924838743837.0_cp, 3326578051521.0_cp/2647936831840.0_cp,       &
            &    -880713585975.0_cp/1841400956686.0_cp, -1428733748635.0_cp/       &
            &    8843423958496.0_cp, 0.0_cp, 0.0_cp, 760814592956.0_cp/            &
            &    3276306540349.0_cp, 760814592956.0_cp/3276306540349.0_cp,         &
            &    -47223648122716.0_cp/6934462133451.0_cp, 71187472546993.0_cp/     &
            &    9669769126921.0_cp, -13330509492149.0_cp/9695768672337.0_cp,      &
            &    11565764226357.0_cp/8513123442827.0_cp, 0.0_cp], [7,7], order=[2,1])

            this%butcher_ass_exp(:)=this%butcher_ass_imp(:)
            this%butcher_c(:)=[0.247_cp, 4276536705230.0_cp/10142255878289.0_cp, &
            &                  0.335_cp, 0.075_cp, 0.7_cp, one]

         case ('KC785')
            this%wimp_lin(1) = 2.0_cp/9.0_cp

            this%butcher_imp(:,:) = reshape([0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,          &
            &  0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 2.0_cp/9.0_cp, 2.0_cp/9.0_cp,          &
            &  0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 2366667076620.0_cp/    &
            &  8822750406821.0_cp, 2366667076620.0_cp/8822750406821.0_cp,             &
            &  2.0_cp/9.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,                 &
            &  -257962897183.0_cp/4451812247028.0_cp, -257962897183.0_cp/             &
            &  4451812247028.0_cp, 128530224461.0_cp/14379561246022.0_cp,             &
            &  2.0_cp/9.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, -486229321650.0_cp/     &
            &  11227943450093.0_cp, -486229321650.0_cp/11227943450093.0_cp,           &
            &  -225633144460.0_cp/6633558740617.0_cp, 1741320951451.0_cp/             &
            &  6824444397158.0_cp, 2.0_cp/9.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,             &
            &  621307788657.0_cp/4714163060173.0_cp, 621307788657.0_cp/               &
            &  4714163060173.0_cp, -125196015625.0_cp/3866852212004.0_cp,             &
            &  940440206406.0_cp/7593089888465.0_cp, 961109811699.0_cp/               &
            &  6734810228204.0_cp, 2.0_cp/9.0_cp, 0.0_cp, 0.0_cp, 2036305566805.0_cp/ &
            &  6583108094622.0_cp, 2036305566805.0_cp/6583108094622.0_cp,             &
            &  -3039402635899.0_cp/4450598839912.0_cp, -1829510709469.0_cp/           &
            &  31102090912115.0_cp, -286320471013.0_cp/6931253422520.0_cp,            &
            &  8651533662697.0_cp/9642993110008.0_cp, 2.0_cp/9.0_cp, 0.0_cp, 0.0_cp,  &
            &  0.0_cp, 3517720773327.0_cp/20256071687669.0_cp, 4569610470461.0_cp/    &
            &  17934693873752.0_cp,  2819471173109.0_cp/11655438449929.0_cp,          &
            &  3296210113763.0_cp/10722700128969.0_cp, -1142099968913.0_cp/           &
            &  5710983926999.0_cp, 2.0_cp/9.0_cp], [8,8], order=[2,1])

            this%butcher_ass_imp(:) = [0.0_cp, 0.0_cp, 3517720773327.0_cp/    &
            &    20256071687669.0_cp, 4569610470461.0_cp/17934693873752.0_cp, &
            &    2819471173109.0_cp/11655438449929.0_cp, 3296210113763.0_cp/  &
            &    10722700128969.0_cp, -1142099968913.0_cp/5710983926999.0_cp, &
            &    2.0_cp/9.0_cp]

            this%butcher_exp(:,:) = reshape([ 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,  &
            &  0.0_cp, 0.0_cp, 0.0_cp, 4.0_cp/9.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,  &
            &  0.0_cp, 0.0_cp, 0.0_cp, 1.0_cp/9.0_cp, 1183333538310.0_cp/              &
            &  1827251437969.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,     &
            &  895379019517.0_cp/9750411845327.0_cp, 477606656805.0_cp/                &
            &  13473228687314.0_cp, -112564739183.0_cp/9373365219272.0_cp, 0.0_cp,     &
            &  0.0_cp, 0.0_cp, 0.0_cp,0.0_cp, -4458043123994.0_cp/13015289567637.0_cp, &
            &  -2500665203865.0_cp/9342069639922.0_cp, 983347055801.0_cp/              &
            &  8893519644487.0_cp, 2185051477207.0_cp/2551468980502.0_cp, 0.0_cp,      &
            &  0.0_cp, 0.0_cp, 0.0_cp, -167316361917.0_cp/17121522574472.0_cp,         &
            &  1605541814917.0_cp/7619724128744.0_cp, 991021770328.0_cp/               &
            &  13052792161721.0_cp,  2342280609577.0_cp/11279663441611.0_cp,           &
            &  3012424348531.0_cp/12792462456678.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,         &
            &  6680998715867.0_cp/14310383562358.0_cp, 5029118570809.0_cp/             &
            &  3897454228471.0_cp, 2415062538259.0_cp/6382199904604.0_cp,              &
            &  -3924368632305.0_cp/6964820224454.0_cp, -4331110370267.0_cp/            &
            &  15021686902756.0_cp, -3944303808049.0_cp/11994238218192.0_cp, 0.0_cp,   &
            &  0.0_cp, 2193717860234.0_cp/3570523412979.0_cp, 2193717860234.0_cp/      &
            &  3570523412979.0_cp, 5952760925747.0_cp/18750164281544.0_cp,             &
            &  -4412967128996.0_cp/6196664114337.0_cp, 4151782504231.0_cp/             &
            &  36106512998704.0_cp,  572599549169.0_cp/6265429158920.0_cp,             &
            &  -457874356192.0_cp/11306498036315.0_cp, 0.0_cp], [8,8], order=[2,1])

            this%butcher_ass_exp(:)=this%butcher_ass_imp(:)
            this%butcher_c(:)=[4.0_cp/9.0_cp, 6456083330201.0_cp/8509243623797.0_cp,  &
            &   1632083962415.0_cp/14158861528103.0_cp, 6365430648612.0_cp/           &
            &   17842476412687.0_cp,  0.72_cp, 0.955_cp, one]

         case ('LZ232')
            this%wimp_lin(1) = half
            this%butcher_imp(:,:) = reshape([  0.0_cp,  0.0_cp, 0.0_cp,  &
                                    &        -0.25_cp,    half, 0.0_cp,  &
                                    &            half,  0.0_cp, half ],  &
                                    &          [3,3],order=[2,1])
            this%butcher_ass_imp(:) = [ half,  0.0_cp, half ]
            this%butcher_exp(:,:) = reshape([  0.0_cp, 0.0_cp, 0.0_cp,  &
                                    &         0.25_cp, 0.0_cp, 0.0_cp,  &
                                    &            -one,    two, 0.0_cp], &
                                    &          [3,3],order=[2,1])
            this%butcher_ass_exp(:) = [ -one,    two, 0.0_cp ]
            this%butcher_c(:) = [0.25_cp, one]
         case ('CK232')
            gam = one-half*sqrt(two)
            this%wimp_lin(1) = gam
            this%butcher_imp(:,:) = reshape([                                     &
            &                         0.0_cp,                     0.0_cp, 0.0_cp, &
            &  -1.0_cp/3.0_cp+half*sqrt(two),                        gam, 0.0_cp, &
            &      0.75_cp-0.25_cp*sqrt(two), -0.75_cp+0.75_cp*sqrt(two),    gam],&
            &                               [3,3],order=[2,1])
            this%butcher_ass_imp(:) = [ 0.75_cp-0.25_cp*sqrt(two), &
            &                          -0.75_cp+0.75_cp*sqrt(two), gam ]
            this%butcher_exp(:,:) = reshape([        0.0_cp,  0.0_cp, 0.0_cp,  &
                                    &         2.0_cp/3.0_cp,  0.0_cp, 0.0_cp,  &
                                    &               0.25_cp, 0.75_cp, 0.0_cp], &
                                    &        [3,3],order=[2,1])
            this%butcher_ass_exp(:) = [ 0.25_cp, 0.75_cp, 0.0_cp ]
            this%butcher_c(:) = [2.0_cp/3.0_cp, one]
         case ('ARS443')
            this%wimp_lin(1) = half
            this%butcher_imp(:,:) = reshape(                               &
            &            [ 0.0_cp,       0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,   &
            &              0.0_cp,         half, 0.0_cp, 0.0_cp, 0.0_cp,   &
            &              0.0_cp,1.0_cp/6.0_cp,   half, 0.0_cp, 0.0_cp,   &
            &              0.0_cp,        -half,   half,   half, 0.0_cp,   &
            &              0.0_cp,       1.5_cp,-1.5_cp,   half,   half],  &
            &              [5,5],order=[2,1])
            this%butcher_ass_imp(:) = [ 0.0_cp, 1.5_cp, -1.5_cp, half, half ]
            this%butcher_exp(:,:) = reshape(                                     &
            &       [           0.0_cp,       0.0_cp,  0.0_cp,   0.0_cp, 0.0_cp, &
            &                    half,        0.0_cp,  0.0_cp,   0.0_cp, 0.0_cp, &
            &         11.0_cp/18.0_cp,1.0_cp/18.0_cp,  0.0_cp,   0.0_cp, 0.0_cp, &
            &           5.0_cp/6.0_cp,-5.0_cp/6.0_cp,    half,   0.0_cp, 0.0_cp, &
            &                 0.25_cp,       1.75_cp, 0.75_cp, -1.75_cp, 0.0_cp],&
            &         [5,5],order=[2,1])
            this%butcher_ass_exp(:) = [ 0.25_cp, 1.75_cp, 0.75_cp, -1.75_cp, 0.0_cp ]
            this%l_imp_calc_rhs(1)=.false.
            this%butcher_c(:) = [half, 2.0_cp/3.0_cp, half, one]
         case ('DBM453')
            this%wimp_lin(1) =  0.32591194130117247_cp
            this%butcher_imp(:,:) = reshape(                                     &
            &       [         0.0_cp,        0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,     &
            &  -0.22284985318525410_cp, 0.32591194130117247_cp, 0.0_cp, 0.0_cp,  &
            &          0.0_cp, -0.46801347074080545_cp, 0.86349284225716961_cp,  &
            &   0.32591194130117247_cp, 0.0_cp, 0.0_cp, -0.46509906651927421_cp, &
            &   0.81063103116959553_cp, 0.61036726756832357_cp,                  &
            &   0.32591194130117247_cp, 0.0_cp, 0.87795339639076675_cp,          &
            &  -0.72692641526151547_cp, 0.75204137157372720_cp,                  &
            &  -0.22898029400415088_cp, 0.32591194130117247_cp], [5,5],order=[2,1])
            this%butcher_exp(:,:) = reshape(                                     &
            &  [ 0.0_cp,       0.0_cp,  0.0_cp,   0.0_cp, 0.0_cp,                &
            &    0.10306208811591838_cp,  0.0_cp,  0.0_cp,   0.0_cp, 0.0_cp,     &
            &   -0.94124866143519894_cp, 1.6626399742527356_cp, 0.0_cp, 0.0_cp,  &
            &    0.0_cp, -1.3670975201437765_cp,  1.3815852911016873_cp,         &
            &    1.2673234025619065_cp, 0.0_cp, 0.0_cp, -0.81287582068772448_cp, &
            &    0.81223739060505738_cp, 0.90644429603699305_cp,                 &
            &    0.094194134045674111_cp, 0.0_cp], [5,5],order=[2,1])
            this%butcher_c(:) = [0.10306208811591838_cp, 0.7213913128175367_cp, &
            &                    1.2818111735198172_cp, one]

            this%butcher_ass_exp(:)=[0.87795339639076672_cp, -0.72692641526151549_cp, &
                                     0.7520413715737272_cp, -0.22898029400415090_cp,  &
                                     0.32591194130117247_cp]
            this%butcher_ass_imp(:)=this%butcher_ass_exp(:)
         case ('BPR353')
            this%wimp_lin(1) = half
            this%butcher_imp(:,:) = reshape(                                   &
            &       [         0.0_cp,        0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,   &
            &                   half,          half, 0.0_cp, 0.0_cp, 0.0_cp,   &
            &         5.0_cp/18.0_cp,-1.0_cp/9.0_cp,   half, 0.0_cp, 0.0_cp,   &
            &                   half,        0.0_cp, 0.0_cp,   half, 0.0_cp,   &
            &                0.25_cp,        0.0_cp,0.75_cp,  -half,   half],  &
            &              [5,5],order=[2,1])
            this%butcher_ass_imp(:) = [ 0.25_cp, 0.0_cp, 0.75_cp, -half, half ]
            this%butcher_exp(:,:) = reshape(                                   &
            &       [        0.0_cp,       0.0_cp,  0.0_cp,   0.0_cp, 0.0_cp,  &
            &                   one,       0.0_cp,  0.0_cp,   0.0_cp, 0.0_cp,  &
            &         4.0_cp/9.0_cp,2.0_cp/9.0_cp,  0.0_cp,   0.0_cp, 0.0_cp,  &
            &               0.25_cp,       0.0_cp, 0.75_cp,   0.0_cp, 0.0_cp,  &
            &               0.25_cp,       0.0_cp, 0.75_cp,   0.0_cp, 0.0_cp], &
            &         [5,5],order=[2,1])
            this%butcher_ass_exp(:) = [ 0.25_cp, 0.0_cp, 0.75_cp, 0.0_cp, 0.0_cp ]
            this%l_exp_calc(4)=.false. ! No need to calculte the explicit solve
            this%butcher_c(:) = [one, 2.0_cp/3.0_cp, one, one]
         case ('LZ453')
            this%wimp_lin(1) = 1.2_cp
            this%butcher_imp(:,:) = reshape(                                                        &
            &[            0.0_cp,           0.0_cp,            0.0_cp,            0.0_cp, 0.0_cp,   &
            &   -44.0_cp/45.0_cp,           1.2_cp,            0.0_cp,            0.0_cp, 0.0_cp,   &
            &  -47.0_cp/300.0_cp,         -0.71_cp,            1.2_cp,            0.0_cp, 0.0_cp,   &
            &           3.375_cp,         -3.25_cp, -59.0_cp/120.0_cp,            1.2_cp, 0.0_cp,   &
            &    89.0_cp/50.0_cp,-486.0_cp/55.0_cp,            8.9_cp,-562.0_cp/275.0_cp, 1.2_cp],  &
            &[5,5],order=[2,1])
            this%butcher_ass_imp(:) = [ 89.0_cp/50.0_cp, -486.0_cp/55.0_cp,    &
            &                           8.9_cp, -562.0_cp/275.0_cp, 1.2_cp ]
            this%butcher_exp(:,:) = reshape(                                   &
            & [            0.0_cp,            0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,  &
            &       2.0_cp/9.0_cp,            0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,  &
            &    71.0_cp/420.0_cp,  23.0_cp/140.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,  &
            &  -281.0_cp/336.0_cp, 187.0_cp/112.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,  &
            &             0.1_cp,             0.0_cp, 0.5_cp, 0.4_cp, 0.0_cp], &
            &  [5,5],order=[2,1])
            this%butcher_ass_exp(:) = [ 0.1_cp, 0.0_cp, 0.5_cp, 0.4_cp, 0.0_cp ]
            this%butcher_c(:) = [2.0_cp/9.0_cp, 1.0_cp/3.0_cp, 5.0_cp/6.0_cp,  &
            &                    one]
         case ('PC2')
            this%wimp_lin(1) = half
            this%butcher_imp(:,:) = reshape([ 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, &
            &                                   half,   half, 0.0_cp, 0.0_cp, &
            &                                   half, 0.0_cp,   half, 0.0_cp, &
            &                                   half, 0.0_cp, 0.0_cp,   half],&
            &                               [4,4],order=[2,1])
            this%butcher_ass_imp(:) = [ half, 0.0_cp, 0.0_cp, half]
            this%butcher_exp(:,:) = reshape([ 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, &
            &                                    one, 0.0_cp, 0.0_cp, 0.0_cp, &
            &                                   half,   half, 0.0_cp, 0.0_cp, &
            &                                   half, 0.0_cp,   half, 0.0_cp],&
            &                               [4,4],order=[2,1])
            this%butcher_ass_exp(:) = [ half, 0.0_cp, half, 0.0_cp ]
            this%butcher_c(:) = [one, one, one]
      end select

      this%wimp_lin(1)      = this%dt(1)*this%wimp_lin(1)
      this%butcher_imp(:,:) = this%dt(1)*this%butcher_imp(:,:)
      this%butcher_exp(:,:) = this%dt(1)*this%butcher_exp(:,:)

      this%butcher_ass_imp(:) = this%dt(1)*this%butcher_ass_imp(:)
      this%butcher_ass_exp(:) = this%dt(1)*this%butcher_ass_exp(:)
         
   end subroutine set_weights
!------------------------------------------------------------------------------
   subroutine set_dt_array(this, dt_new, dt_min, time, n_log_file,  &
              &            n_time_step, l_new_dtNext)
      !
      ! This subroutine adjusts the time step
      !

      class(type_dirk) :: this

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
   subroutine set_imex_rhs(this, rhs, dfdt)
      !
      ! This subroutine assembles the right-hand-side of an IMEX scheme
      !

      class(type_dirk) :: this

      !-- Input variables:
      type(type_tarray), intent(in) :: dfdt

      !-- Output variable
      complex(cp), intent(out) :: rhs(dfdt%llm:dfdt%ulm,dfdt%nRstart:dfdt%nRstop)

      !-- Local variables
      integer :: n_stage, n_r, startR, stopR

      !$omp parallel default(shared) private(startR,stopR,n_r)
      startR=dfdt%nRstart; stopR=dfdt%nRstop
      call get_openmp_blocks(startR,stopR)

      do n_r=startR,stopR
         rhs(dfdt%llm:dfdt%ulm,n_r)=dfdt%old(dfdt%llm:dfdt%ulm,n_r,1)
      end do

      do n_stage=1,this%istage
         do n_r=startR,stopR
            rhs(dfdt%llm:dfdt%ulm,n_r)=rhs(dfdt%llm:dfdt%ulm,n_r) +          &
            &                       this%butcher_exp(this%istage+1,n_stage)* &
            &                       dfdt%expl(dfdt%llm:dfdt%ulm,n_r,n_stage)
         end do
      end do

      do n_stage=1,this%istage
         do n_r=startR,stopR
            rhs(dfdt%llm:dfdt%ulm,n_r)=rhs(dfdt%llm:dfdt%ulm,n_r) +          &
            &                       this%butcher_imp(this%istage+1,n_stage)* &
            &                       dfdt%impl(dfdt%llm:dfdt%ulm,n_r,n_stage)
         end do
      end do
      !$omp end parallel

   end subroutine set_imex_rhs
!------------------------------------------------------------------------------
   subroutine set_imex_rhs_ghost(this, rhs, dfdt, start_lm, stop_lm, ng)
      !
      ! This subroutine assembles the right-hand-side of an IMEX scheme in case
      ! an array with ghosts zones is provided
      !

      class(type_dirk) :: this

      !-- Input variables:
      type(type_tarray), intent(in) :: dfdt
      integer,           intent(in) :: start_lm ! Starting lm index
      integer,           intent(in) :: stop_lm  ! Stopping lm index
      integer,           intent(in) :: ng       ! Number of ghost zones

      !-- Output variable
      complex(cp), intent(out) :: rhs(dfdt%llm:dfdt%ulm,dfdt%nRstart-ng:dfdt%nRstop+ng)

      !-- Local variables
      integer :: n_stage, n_r

      do n_r=dfdt%nRstart,dfdt%nRstop
         rhs(start_lm:stop_lm,n_r)=dfdt%old(start_lm:stop_lm,n_r,1)
      end do

      do n_stage=1,this%istage
         do n_r=dfdt%nRstart,dfdt%nRstop
            rhs(start_lm:stop_lm,n_r)=rhs(start_lm:stop_lm,n_r) +            &
            &                       this%butcher_exp(this%istage+1,n_stage)* &
            &                       dfdt%expl(start_lm:stop_lm,n_r,n_stage)
         end do
      end do

      do n_stage=1,this%istage
         do n_r=dfdt%nRstart,dfdt%nRstop
            rhs(start_lm:stop_lm,n_r)=rhs(start_lm:stop_lm,n_r) +            &
            &                       this%butcher_imp(this%istage+1,n_stage)* &
            &                       dfdt%impl(start_lm:stop_lm,n_r,n_stage)
         end do
      end do

   end subroutine set_imex_rhs_ghost
!------------------------------------------------------------------------------
   subroutine assemble_imex(this, rhs, dfdt)
      !
      ! This subroutine performs the assembly stage of an IMEX-RK scheme
      !

      class(type_dirk) :: this

      !-- Input variables:
      type(type_tarray), intent(in) :: dfdt

      !-- Output variable
      complex(cp), intent(out) :: rhs(dfdt%llm:dfdt%ulm,dfdt%nRstart:dfdt%nRstop)

      !-- Local variables
      integer :: n_stage, n_r, startR, stopR

      !$omp parallel default(shared) private(startR,stopR,n_r)
      startR=dfdt%nRstart; stopR=dfdt%nRstop
      call get_openmp_blocks(startR,stopR)

      do n_r=startR,stopR
         rhs(dfdt%llm:dfdt%ulm,n_r)=dfdt%old(dfdt%llm:dfdt%ulm,n_r,1)
      end do

      do n_stage=1,this%nstages
         do n_r=startR,stopR
            rhs(dfdt%llm:dfdt%ulm,n_r)=rhs(dfdt%llm:dfdt%ulm,n_r) +        &
            &                       this%butcher_ass_exp(n_stage)*         &
            &                       dfdt%expl(dfdt%llm:dfdt%ulm,n_r,n_stage)
         end do
      end do

      do n_stage=1,this%nstages
         do n_r=startR,stopR
            rhs(dfdt%llm:dfdt%ulm,n_r)=rhs(dfdt%llm:dfdt%ulm,n_r) +       &
            &                       this%butcher_ass_imp(n_stage)*        &
            &                       dfdt%impl(dfdt%llm:dfdt%ulm,n_r,n_stage)
         end do
      end do
      !$omp end parallel

   end subroutine assemble_imex
!------------------------------------------------------------------------------
   subroutine set_imex_rhs_scalar(this, rhs, dfdt)
      !
      ! This subroutine assembles the right-hand-side of an IMEX scheme
      !

      class(type_dirk) :: this

      !-- Input variables:
      type(type_tscalar), intent(in) :: dfdt

      !-- Output variable
      real(cp), intent(out) :: rhs

      !-- Local variables
      integer :: n_stage

      rhs=dfdt%old(1)

      do n_stage=1,this%istage
         rhs=rhs+this%butcher_exp(this%istage+1,n_stage)*dfdt%expl(n_stage)
      end do

      do n_stage=1,this%istage
         rhs=rhs+this%butcher_imp(this%istage+1,n_stage)*dfdt%impl(n_stage)
      end do

   end subroutine set_imex_rhs_scalar
!------------------------------------------------------------------------------
   subroutine assemble_imex_scalar(this, rhs, dfdt)
      !
      ! This subroutine performs the assembly stage of an IMEX-RK scheme
      ! for a scalar quantity
      !
      class(type_dirk) :: this

      !-- Input variable
      type(type_tscalar), intent(in) :: dfdt

      !-- Output variable
      real(cp), intent(out) :: rhs

      !-- Local variable
      integer :: n_stage

      rhs=dfdt%old(1)
      do n_stage=1,this%nstages
         rhs=rhs + this%butcher_ass_exp(n_stage)*dfdt%expl(n_stage)
      end do

      do n_stage=1,this%nstages
         rhs=rhs + this%butcher_ass_imp(n_stage)*dfdt%impl(n_stage)
      end do

   end subroutine assemble_imex_scalar
!------------------------------------------------------------------------------
   subroutine rotate_imex(this, dfdt)
      !
      ! This subroutine is used to roll the time arrays from one time step
      !

      class(type_dirk) :: this

      !-- Output variables:
      type(type_tarray), intent(inout) :: dfdt

   end subroutine rotate_imex
!------------------------------------------------------------------------------
   subroutine rotate_imex_scalar(this, dfdt)
      !
      ! This subroutine is used to roll the time arrays from one time step
      !

      class(type_dirk) :: this

      !-- Output variables:
      type(type_tscalar), intent(inout) :: dfdt

   end subroutine rotate_imex_scalar
!------------------------------------------------------------------------------
   subroutine bridge_with_cnab2(this)

      class(type_dirk) :: this

   end subroutine bridge_with_cnab2
!------------------------------------------------------------------------------
   subroutine start_with_ab1(this)

      class(type_dirk) :: this

   end subroutine start_with_ab1
!------------------------------------------------------------------------------
   subroutine get_time_stage(this, tlast, tstage)

      class(type_dirk) :: this
      real(cp), intent(in) :: tlast
      real(cp), intent(out) :: tstage

      tstage = tlast+this%dt(1)*this%butcher_c(this%istage)

   end subroutine get_time_stage
!------------------------------------------------------------------------------
end module dirk_schemes
