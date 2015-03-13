MODULE Namelists
  USE const
  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE torsional_oscillations
  USE init_fields
  USE Grenoble
  USE logic
  USE output_data
  USE parallel_mod
  USE Bext
  USE movie_data,ONLY: movie,n_movies
  USE charmanip, ONLY: length_to_blank
  USE blocking, only: cacheblock_size_in_B
  IMPLICIT NONE

  private
  PUBLIC :: readNamelists, writeNamelists

CONTAINS
    !***********************************************************************

    !    !------------ This is release 2 level 5  --------------!
    !    !------------ Created on 1/18/02  by JW. -----------

    !  +-------------+----------------+------------------------------------+
    !  |                                                                   |
    !  |  Purpose of this subroutine is to read the input namelists.       |
    !  |  This program also determins logical parameters that are stored   |
    !  |  in c_logic.f.                                                    |
    !  |                                                                   |
    !  +-------------------------------------------------------------------+
  SUBROUTINE readNamelists
    !-- Local stuff
    INTEGER :: n
    INTEGER :: nCounts
    REAL(kind=8) :: sCMB(4*n_impS_max),rad ! cmb heat boundary condition
    logical :: log_does_exist
    INTEGER :: length
    integer :: argument_count
    character(len=100) :: input_filename
    integer :: nThreadsRun
#ifdef WITH_MPI
  character(len=100) :: new_tag
#endif

    !-- Name lists:
    INTEGER :: runHours,runMinutes,runSeconds
    NAMELIST/grid/n_r_max,n_cheb_max,n_phi_tot, &
         n_r_ic_max,n_cheb_ic_max, &
         minc,nalias

    NAMELIST/control/                                &
         mode,tag,n_time_steps, &
         n_tScale,n_lScale,alpha,enscale, &
         l_update_v,l_update_b,l_update_s, &
         dtstart,dtMax,courfac,alffac,intfac,n_cour_step, &
         difnu,difeta,difkap,ldif,ldifexp,l_correct_AMe, &
         l_correct_AMz,tEND,l_non_rot,l_isothermal, &
         l_newmap,alph1,alph2,l_interior_model, &
         runHours,runMinutes,runSeconds,nThreadsRun,&
         & cacheblock_size_in_B
    

    NAMELIST/phys_param/                             &
         ra,pr,prmag,ek,epsc0,radratio, &
         ktops,kbots,ktopv,kbotv,ktopb,kbotb, &
         s_top,s_bot,impS,sCMB, &
         nVarCond,con_DecRate,con_RadRatio,con_LambdaMatch, &
         con_LambdaOut,con_FuncWidth, &
         strat,polind,g0,g1,g2,r_cut_model, &
         nVarDiff,nVarVisc,difExp,nVarEps

    NAMELIST/B_external/ &
         rrMP,amp_imp,expo_imp,bmax_imp,n_imp

    NAMELIST/start_field/                            &
         l_start_file,start_file,inform, &
         l_reset_t,scale_s,scale_b,scale_v,tipdipole, &
         init_s1,init_s2,init_v1,init_b1,imagcon,tmagcon, &
         amp_s1,amp_s2,amp_v1,amp_b1

    NAMELIST/output_control/                         &
         n_graph_step,n_graphs,t_graph, &
         t_graph_start,t_graph_stop,dt_graph, &
         l_graph_time, &
         n_stores,n_rst_step,n_rsts,t_rst, &
         t_rst_start,t_rst_stop,dt_rst, &
         n_log_step,n_logs,t_log, &
         t_log_start,t_log_stop,dt_log, &
         n_p_step,n_ps,t_p, &
         t_p_start,t_p_stop,dt_p, &
         n_spec_step,n_specs,t_spec, &
         t_spec_start,t_spec_stop,dt_spec, &
         n_cmb_step,n_cmbs,t_cmb, &
         t_cmb_start,t_cmb_stop,dt_cmb, &
         n_r_field_step,n_r_fields,t_r_field, &
         t_r_field_start,t_r_field_stop,dt_r_field, &
         n_Bpot_step,n_Bpots,t_Bpot, &
         t_Bpot_start,t_Bpot_stop,dt_Bpot, &
         n_Vpot_step,n_Vpots,t_Vpot, &
         t_Vpot_start,t_Vpot_stop,dt_Vpot, &
         n_Tpot_step,n_Tpots,t_Tpot, &
         t_Tpot_start,t_Tpot_stop,dt_Tpot, &
         n_pot_step,n_pots,t_pot, &
         t_pot_start,t_pot_stop,dt_pot, &
         ngform,runid, &
         movie,n_movie_step,n_movie_frames,t_movie, &
         t_movie_start,t_movie_stop,dt_movie, &
         n_TO_step,n_TOs,t_TO, &
         t_TO_start,t_TO_stop,dt_TO, &
         n_TOZ_step,n_TOZs,t_TOZ, &
         t_TOZ_start,t_TOZ_stop,dt_TOZ, &
         n_TOmovie_step,n_TOmovie_frames,t_TOmovie, &
         t_TOmovie_start,t_TOmovie_stop,dt_TOmovie, &
         l_movie,l_average,l_save_out,l_true_time, &
         l_cmb_field,l_rMagSpec,l_DTrMagSpec, &
         l_dt_cmb_field,l_max_cmb, &
         l_r_field,l_r_fieldT,n_r_step,l_max_r,n_r_array, &
         l_TO,l_TOmovie,l_hel,lVerbose, &
         l_AM,l_power,l_drift,l_storeBpot,l_storeVpot, &
         l_storeTpot,l_storePot,sDens,zDens,l_RMS, &
         l_RMStest,l_par,l_corrMov,rCut,rDea,l_prms, &
         l_plotmap,l_PV,l_iner,l_viscBcCalc,l_fluxProfs

    NAMELIST/mantle/conductance_ma,nRotMa,rho_ratio_ma, &
         omega_ma1,omegaOsz_ma1,tShift_ma1, &
         omega_ma2,omegaOsz_ma2,tShift_ma2

    NAMELIST/inner_core/sigma_ratio,nRotIc,rho_ratio_ic, &
         omega_ic1,omegaOsz_ic1,tShift_ic1, &
         omega_ic2,omegaOsz_ic2,tShift_ic2, &
         BIC

    !-- end of declaration
    !----------------------------------------------------------------------

    DO n=1,4*n_impS_max
       sCMB(n)=0.D0
    END DO

    runHours  =0
    runMinutes=0
    runSeconds=0

    !-- Set default values of control parameters:
    CALL defaultNamelists


    ! get the filename of the input file as first argument from the command line
    argument_count = command_argument_count()
    IF (argument_count.EQ.0) THEN
       WRITE(*,"(A,/,A)") "The filename of the input file must be provided as first argument.",&
            &"Aborting!"
       STOP
    ELSE
       CALL get_command_argument(1,input_filename)

       OPEN(105,file=trim(input_filename))

       !-- Reading control parameters from namelists in STDIN:
       if (rank.eq.0) write(*,*) '!  Reading grid parameters!'
       READ(105,grid)

       !-- Reading control parameters from namelists in STDIN:
       if (rank.eq.0) write(*,*) '!  Reading control parameters!'
       READ(105,control)

       !-- Reading physical parameters from namelists in STDIN:
       if (rank.eq.0) write(*,*) '!  Reading physical parameters!'
       READ(105,phys_param)

       !-- Reading external field parameters for feedback:
       !           WRITE(*,*) '!  Reading B external parameters!'
       !           READ(105,B_external)

       !-- Reading start field info from namelists in STDIN:
       if (rank.eq.0) write(*,*) '!  Reading start information!'
       READ(105,start_field)

       !-- Reading output parameters from namelists in STDIN:
       if (rank.eq.0) write(*,*) '!  Reading output information!'
       READ(105,output_control)

       !-- Reading mantle parameters from namelists in STDIN:
       if (rank.eq.0) write(*,*) '!  Reading mantle information!'
       READ(105,mantle)

       !-- Reading inner-core parameter from namelists in STDIN:
       if (rank.eq.0) write(*,*) '!  Reading inner core information!'
       READ(105,inner_core)

       close(105)
       !-- Correcting some parameters:
    END IF


#ifdef WITH_MPI
    tag_wo_rank=tag
    !WRITE(new_tag,"(I4.4,A,A)") rank,".",TRIM(tag)
    !tag=new_tag
#endif
    !-- Does log-file already exist?
    log_file='log.'//tag
    IF (rank.EQ.0) THEN
       INQUIRE(file=TRIM(log_file),exist=log_does_exist)
    END IF
    CALL MPI_Bcast(log_does_exist,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    IF (log_does_exist) THEN
       IF (rank.EQ.0) THEN 
          WRITE(*,*)
          WRITE(*,*) '! The log-file exists already !'
          WRITE(*,*) '! I add _BIS to the tag and create new files!'
       END IF
       length=length_to_blank(tag)
       tag=tag(1:length)//'_BIS'
    END IF

    n_stores=MAX(n_stores,n_rsts)

    l_TOmovie=l_TOmovie.AND.l_TO

    !-- Determine what has to be calculated depending on mode:

    lMagMem  =1
    lAveMem  =0
    ldtBMem  =0
    lStressMem=0
    lMovieMem=0
    l_conv   =.TRUE.
    l_conv_nl=.TRUE.
    l_mag    =.TRUE.
    l_mag_nl =.TRUE.
    l_mag_kin=.FALSE.
    l_mag_LF =.TRUE.
    l_heat   =.TRUE.
    l_heat_nl=.TRUE.
    l_SRIC   =.FALSE.
    l_SRMA   =.FALSE.

    IF ( mode == 1 ) THEN

       !----- Only convection:
       l_mag   =.FALSE.
       l_mag_nl=.FALSE.
       l_mag_LF=.FALSE.

    ELSE IF ( mode == 2 ) THEN

       !----- Kinematic dynamo:
       l_conv=.FALSE.
       l_conv_nl=.FALSE.
       l_mag_kin=.TRUE.
       l_mag_LF =.FALSE.
       l_heat   =.FALSE.
       l_heat_nl=.FALSE.

    ELSE IF ( mode == 3 ) THEN

       !----- Magnetic decay modes:
       l_conv   =.FALSE.
       l_conv_nl=.FALSE.
       l_mag_nl =.FALSE.
       l_mag_LF =.FALSE.
       l_heat   =.FALSE.
       l_heat_nl=.FALSE.

    ELSE IF ( mode == 4 ) THEN

       !----- Magneto convection:
       l_mag    =.FALSE.
       l_mag_nl =.FALSE.
       l_mag_LF =.FALSE.
       l_rot_ic =.FALSE.
       l_rot_ma =.FALSE.

    ELSE IF ( mode == 5 ) THEN

       !----- Onset of convection (linear):
       l_conv_nl=.FALSE.
       l_mag    =.FALSE.
       l_mag_nl =.FALSE.
       l_mag_LF =.FALSE.
       l_rot_ic =.FALSE.
       l_rot_ma =.FALSE.

    ELSE IF ( mode == 6 ) THEN

       !----- Selfconsistent dynamo, but no Lorentz Force
       l_mag_LF=.FALSE.

    ELSE IF ( mode == 7 ) THEN

       !----- Super-rotating IC or MA, no convection, no magnetic field
       l_heat   =.FALSE.
       l_mag    =.FALSE.
       l_mag_nl =.FALSE.
       l_mag_LF =.FALSE.

    ELSE IF ( mode == 8 ) THEN

       !----- Super-rotating IC or MA, no convection, plus magnetic field
       l_heat   =.FALSE.

    ELSE IF ( mode == 9 ) THEN

       !----- Super-rotating IC or MA, no convection, magnetic field, no LF
       l_heat   =.FALSE.
       l_mag_LF =.FALSE.


    ELSE IF ( mode == 10 ) THEN

       !----- Super-rotating IC or MA, no convection, no magnetic field, no LF, no advection
       l_heat   =.FALSE.
       l_conv_nl=.FALSE.
       l_mag    =.FALSE.
       l_mag_nl =.FALSE.
       l_mag_LF =.FALSE.

    END IF

    IF ( mode == 7 .OR. mode == 8 .OR. &
         mode == 9 .OR. mode == 10 ) THEN
       !kbotv=2
       !ktopv=2
       l_correct_AMz=.FALSE.
       l_heat       =.FALSE.
    END IF

    IF ( nRotIc > 0 ) THEN
       l_rot_ic=.TRUE.
    ELSE IF ( nRotIc == 0 ) THEN
       l_rot_ic=.FALSE.
    ELSE IF ( nRotIc == -1 ) THEN
       l_rot_ic=.TRUE.
       l_SRIC  =.TRUE.
    END IF

    IF ( nRotMa > 1 ) THEN
       l_rot_ma=.TRUE.
    ELSE IF ( nRotMa == 0 ) THEN
       l_rot_ma=.FALSE.
    ELSE IF ( nRotMa == -1 ) THEN
       l_rot_ma=.TRUE.
       l_SRMA  =.TRUE.
    END IF

    IF ( ek < 0.D0 ) l_non_rot= .TRUE. 
    IF ( l_non_rot ) THEN
       l_corr=.FALSE.
       ek=-1.D0 ! used as a flag, not used for the calculation
    ELSE
       l_corr=.TRUE.
    END IF

    IF ( strat > 0.D0 ) l_anel= .TRUE. 

    IF ( l_interior_model ) l_anel= .TRUE. 

    IF ( prmag == 0.D0 ) THEN
       l_mag   =.FALSE.
       l_mag_nl=.FALSE.
       l_mag_LF=.FALSE.
    END IF
    IF ( .NOT. l_conv ) THEN
       l_conv_nl=.FALSE.
       l_mag_LF =.FALSE.
    END IF

!-- JW 10.Apr.2014: new checking of magnetic boundary condition.
    IF ( kbotb > 4 ) THEN
       WRITE(*,*) '! Only outer boundary conditions kbotb<=4 implemented!'
       STOP
    END IF
    IF ( sigma_ratio == 0.D0 ) THEN
       l_cond_ic=.FALSE.
       IF ( kbotb == 3 ) THEN
          WRITE(*,*) '! For an insulating  IC with sigma_ratio=0   !'
          WRITE(*,*) '! boundary condition kbotb=3 is not appropriate!'
          STOP
       END IF
    ELSE
       l_cond_ic=.TRUE.      ! tell the code to use a conducting inner core
       IF ( kbotb .NE. 3 ) THEN
          WRITE(*,*) '! For a conducting IC with sigma_ratio>0   !'
          WRITE(*,*) '! boundary condition kbotb=3 is appropriate!'
          STOP
       END IF
    END IF

    IF ( ktopb > 4 ) THEN
       WRITE(*,*) '! Only outer boundary conditions ktopb<=4 implemented!'
       STOP
    END IF
    IF ( conductance_ma == 0.D0 ) THEN
       l_cond_ma=.FALSE.
       IF ( ktopb == 3 ) THEN
          WRITE(*,*) '! For an insulating mantle with conductance_ma=0 !'
          WRITE(*,*) '! boundary condition ktopb=3 is not appropriate!'
          STOP
       END IF
    ELSE
       l_cond_ma=.TRUE.      ! tell the code to use a conducting mantle
       IF ( ktopb .NE. 3 ) THEN
          WRITE(*,*) '! For a conducting mantle with conductance_ma>0   !'
          WRITE(*,*) '! boundary condition ktopb=3 is appropriate!'
          STOP
       END IF
    END IF

    IF ( .NOT. l_mag ) THEN
       prmag    =0.D0
       l_mag_nl =.FALSE.
       l_cond_ic=.FALSE.
       lMagMem  =0
    END IF

    IF ( l_corrMov ) l_par= .TRUE. 

    !--- Stuff for the radial non-linear mapping
    !     alph1 can be any positive number, above 0
    !     alph2 has to be a value between -1 and 1 (interval in Chebyshev space)
    if ( (alph1 == 0.D0) .OR. &
         (alph2 < -1.D0) .OR. (alph2 > 1.D0) ) then
       l_newmap=.FALSE.
    elseif ( l_newmap .AND. (alph1 < 0.D0) ) then
       alph1=dabs(alph1)
    end if

    !--- Stuff for spherical magnetosphere boundary: rrMP=r(magnetosphere)/r_core
    IF ( n_imp /= 0 ) imagcon=0
    IF ( n_imp == 1 .AND. rrMP <= 1.D0 ) THEN
       WRITE(*,*) '! For runs with n_imp=1!'
       WRITE(*,*) '! please provide rrMP>1!'
       STOP
    ELSE IF ( n_imp >= 2 .AND. amp_imp == 0.D0 ) THEN
       WRITE(*,*) '! For runs with n_imp=2!'
       WRITE(*,*) '! please provide amp_imp!'
       STOP
    END IF
    IF ( imagcon /= 0 .AND. tmagcon == 0 ) tmagcon=1.D18

    IF ( imagcon == -10 ) THEN
       !----- This is used to test variable conductivity
       !      with an analytical solution, see s_initB.f
       l_conv   =.FALSE.
       l_conv_nl=.FALSE.
       l_mag_kin=.TRUE.
       l_mag_nl =.FALSE.
       l_mag_LF =.FALSE.
       l_heat   =.FALSE.
       l_heat_nl=.FALSE.
       l_cond_ic=.FALSE.
       l_cond_ma=.FALSE.
       kbotv    =2
       ktopv    =2
       kbotb    =1
       ktopb    =1
       nRotMa   =0
       nRotIc   =0
       radratio =0.5
    END IF

    IF ( l_rot_ma ) THEN
       WRITE(*,*)
       WRITE(*,*) '! I ALLOW FOR ROTATING MANTLE.'
       IF ( ktopv == 1 .AND. .NOT. l_cond_ma ) THEN
          WRITE(*,*)
          WRITE(*,*) '! No torques on mantle!'
          WRITE(*,*) '! I dont update mantle rotation omega_ma.'
          l_rot_ma=.FALSE.
       END IF
    END IF

    IF ( .NOT. l_mag ) THEN
       l_cmb_field   =.FALSE.
       l_dt_cmb_field=.FALSE.
       l_rMagSpec    =.FALSE.
       l_DTrMagSpec  =.FALSE.
       l_storeBpot   =.FALSE.
    END IF

    l_b_nl_icb=.FALSE.
    IF ( l_mag_nl .AND. kbotv == 1 .AND. l_cond_ic ) THEN
       WRITE(*,*)
       WRITE(*,*) '! Nonlinear magnetic BC required at ICB!'
       l_b_nl_icb=.TRUE.
    END IF

    l_b_nl_cmb=.FALSE.
    IF ( l_mag_nl .AND. ktopv == 1 .AND. l_cond_ma ) THEN
       WRITE(*,*)
       WRITE(*,*) '! Nonlinear magnetic BC required at CMB!'
       l_b_nl_cmb=.TRUE.
    END IF

    !-- Special matrix for z(l=1,m=0) which is the solid body rotation:
    l_z10mat=.FALSE.
    IF ( ( l_rot_ma .AND. ktopv == 2 ) .OR. &
         ( l_rot_ic .AND. kbotv == 2 )      ) l_z10mat= .TRUE. 

    !-- Check Courant criteria at even time steps:
    IF ( MOD(n_cour_step,2) /= 0 ) n_cour_step=n_cour_step+1

    !-- Check whether memory has been reserved:
    IF ( l_average ) lAveMem=1
    IF ( l_TO ) lStressMem=1
    IF ( l_RMS .OR. l_DTrMagSpec ) ldtBMem=1
    IF ( l_movie .OR. l_TOmovie ) lMovieMem=1

    !-- Output of angular moment?
    l_AM=l_AM .OR. l_correct_AMe .OR. l_correct_AMz
    l_AM=l_AM .OR. l_power

    !-- Heat boundary condition
    IF ( impS /= 0 ) THEN
       rad=4.D0*DATAN(1.D0)/180
       n_impS=0
       DO n=1,4*n_impS_max,4
          IF ( sCMB(n) /= 0.D0 ) THEN
             n_impS=n_impS+1
             peakS(n_impS) =sCMB(n)
             thetaS(n_impS)=rad*sCMB(n+1)
             phiS(n_impS)  =rad*sCMB(n+2)
             widthS(n_impS)=rad*sCMB(n+3)
          END IF
       END DO
    END IF

    !-- Grenoble stuff:
    lGrenoble=.FALSE.
    IF ( BIC /= 0.D0 .AND. l_mag ) THEN
       lGrenoble=.TRUE.
       WRITE(*,*)
       WRITE(*,*) '! Running the Grenoble case !'
    END IF


    ! Setting up truncation is required to set up ldif and l_max_r
    CALL initialize_truncation

    !-- Coeffs at radial levels:
    IF ( l_r_fieldT ) l_r_field=.TRUE.

    IF ( l_r_field ) THEN
       IF ( n_r_step == 0 ) n_r_step=2
       IF ( l_max_r == 0 .OR. l_max_r > l_max ) l_max_r=l_max
       nCounts = COUNT(n_r_array>0)

       IF ( nCounts > 0 ) THEN
           n_coeff_r_max=nCounts
       ELSE
           n_coeff_r_max=5
       END IF
    END IF

    !-- ldif determins at which l hyperdiffusivity starts:
    ldif=MAX(1,ldif)
    ldif=MIN(l_max,ldif)

    !-- Coeffs at CMB:
    l_max_cmb=MIN(l_max_cmb,l_max)



    !-- Maximum run time specified?
    runTimeLimit(1)=runHours
    runTimeLimit(2)=runMinutes
    runTimeLimit(3)=runSeconds
    runTimeLimit(4)=0
    l_runTimeLimit =.FALSE.
    IF ( runHours+runMinutes+runSeconds > 0 ) &
         l_runTimeLimit= .TRUE. 


    RETURN
  end SUBROUTINE readNamelists
  
  !***********************************************************************
  
  !     !------------ This is release 2 level 1  --------------!
  !     !------------ Created on 1/17/02  by JW. --------------!
  
  !  +-------------+----------------+------------------------------------+
  !  |                                                                   |
  !  |  Purpose of this subroutine is to write the namelist to           |
  !  |  file unit n_out. This file has to be open before calling this    |
  !  |  routine.                                                         |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+

  SUBROUTINE writeNamelists(n_out)

    !-- input:
    INTEGER :: n_out

    !-- local:
    INTEGER :: l,m,n,i
    REAL(kind=8) ::  rad,sq4pi
    INTEGER :: length

    !-- end of declaration
    !----------------------------------------------------------------------


    !-- Output of name lists:

    write(n_out,*) "&control"
    write(n_out,'(''  mode        ='',i3,'','')') mode
    length=length_to_blank(tag)
    write(n_out,*) " tag         = """,tag(1:length),""","
    write(n_out,'(''  n_time_steps='',i8,'','')') n_time_steps
    write(n_out,'(''  n_tScale    ='',i3,'','')') n_tScale
    write(n_out,'(''  n_lScale    ='',i3,'','')') n_lScale
    write(n_out,'(1p,''  alpha       ='',d14.6,'','')')   alpha
    write(n_out,'(1p,''  enscale     ='',d14.6,'','')') enscale
    write(n_out,'(''  l_update_v  ='',l3,'','')') l_update_v
    write(n_out,'(''  l_update_b  ='',l3,'','')') l_update_b
    write(n_out,'(''  l_update_s  ='',l3,'','')') l_update_s
    write(n_out,'(''  l_newmap    ='',l3,'','')') l_newmap
    write(n_out,'(1p,''  alph1       ='',d14.6,'','')') alpha1
    write(n_out,'(1p,''  alph2       ='',d14.6,'','')') alpha2
    write(n_out,'(1p,''  dtstart     ='',d14.6,'','')') dtstart*tScale
    write(n_out,'(1p,''  dtMax       ='',d14.6,'','')') tScale*dtMax
    write(n_out,'(1p,''  courfac     ='',d14.6,'','')') courfac
    write(n_out,'(1p,''  alffac      ='',d14.6,'','')')  alffac
    write(n_out,'(1p,''  intfac      ='',d14.6,'','')')  intfac
    write(n_out,'(''  n_cour_step ='',i5,'','')') n_cour_step
    write(n_out,'(1p,''  difnu       ='',d14.6,'','')') difnu
    write(n_out,'(1p,''  difeta      ='',d14.6,'','')') difeta
    write(n_out,'(1p,''  difkap      ='',d14.6,'','')') difkap
    write(n_out,'(''  ldif         ='',i3,'','')') ldif
    write(n_out,'(''  ldifexp      ='',i3,'','')') ldifexp
    write(n_out,'(''  l_correct_AMe='',l3,'','')') l_correct_AMe
    write(n_out,'(''  l_correct_AMz='',l3,'','')') l_correct_AMz
    write(n_out,'(''  l_non_rot    ='',l3,'','')') l_non_rot
    write(n_out,'(''  l_runTimeLimit='',l3,'','')') l_runTimeLimit
    write(n_out,'(''  runHours     ='',i6,'','')') runTimeLimit(1)
    write(n_out,'(''  runMinutes   ='',i4,'','')') runTimeLimit(2)
    write(n_out,'(''  runSeconds   ='',i4,'','')') runTimeLimit(3)
    write(n_out,'(1p,''  tEND        ='',d14.6,'','')') tEND
    !WRITE(n_out,'(''  nThreadsRun  ='',i4,'','')') nThreadsRun
    write(n_out,*) "/"

    sq4pi=DSQRT(16.D0*DATAN(1.D0))
    write(n_out,*) "&phys_param"
    write(n_out,'(1p,''  ra          ='',d14.6,'','')') ra
    write(n_out,'(1p,''  pr          ='',d14.6,'','')') pr
    write(n_out,'(1p,''  prmag       ='',d14.6,'','')') prmag
    write(n_out,'(1p,''  ek          ='',d14.6,'','')') ek
    write(n_out,'(1p,''  epsc0       ='',d14.6,'','')') epsc/sq4pi
    write(n_out,'(1p,''  strat       ='',d14.6,'','')') strat
    write(n_out,'(1p,''  polind      ='',d14.6,'','')') polind
    write(n_out,'(1p,''  radratio    ='',d14.6,'','')') radratio
    write(n_out,'(''  g0          ='',d14.6,'','')') g0
    write(n_out,'(''  g1          ='',d14.6,'','')') g1
    write(n_out,'(''  g2          ='',d14.6,'','')') g2
    write(n_out,'(''  ktopv       ='',i3,'','')') ktopv
    write(n_out,'(''  kbotv       ='',i3,'','')') kbotv
    write(n_out,'(''  ktopb       ='',i3,'','')') ktopb
    write(n_out,'(''  kbotb       ='',i3,'','')') kbotb
    write(n_out,'("  Bottom boundary l,m,S:")')

    !--- Heat boundary condition:
    write(n_out,'(''  ktops       ='',i3,'','')') ktops
    write(n_out,'(''  kbots       ='',i3,'','')') kbots
    DO m=0,m_max,minc
        DO l=m,l_max
            IF ( bots(l,m) /= 0.D0 ) WRITE(n_out,'(1p,4x,2i4,2d14.6)') &
                 l,m,REAL(bots(l,m))/sq4pi,AIMAG(bots(l,m))/sq4pi
        END DO
    END DO
    write(n_out,'("  Top boundary l,m,S:")')
    DO m=0,m_max,minc
        DO l=m,l_max
            IF ( tops(l,m) /= 0.D0 ) WRITE(n_out,'(1p,4x,2i4,2d14.6)') &
                 l,m,REAL(tops(l,m))/sq4pi,AIMAG(tops(l,m))/sq4pi
        END DO
    END DO
    write(n_out,'(''  impS        ='',i3,'','')') impS
    rad=4.D0*DATAN(1.D0)/180.D0
    DO i=1,n_impS
       IF ( i == 1 ) THEN
          WRITE(n_out,'(A)',advance='NO') "  sCMB        ="
       ELSE
          WRITE(n_out,'(A)',advance='NO') "               "
       END IF
       WRITE(n_out,'(1p,4(D14.6,A))') &
            &peakS(i)/rad,",",&
            &thetaS(i)/rad,",",&
            &phiS(i)/rad,",",&
            &widthS(i)/rad,","
    END DO

    !----- Conductivity variation:
    write(n_out,'(''  nVarCond    ='',i3,'','')') nVarCond
    write(n_out,'(1p,''  con_DecRate    ='',d14.6,'','')') con_DecRate
    write(n_out,'(1p,''  con_RadRatio   ='',d14.6,'','')') con_RadRatio
    write(n_out,'(1p,''  con_LambdaMatch='',d14.6,'','')') con_LambdaMatch
    write(n_out,'(1p,''  con_LambdaOut  ='',d14.6,'','')') con_LambdaOut
    write(n_out,'(1p,''  con_FuncWidth  ='',d14.6,'','')') con_FuncWidth

    !----- Thermal diffusivity variation:
    write(n_out,'(1p,''  difExp    ='',d14.6,'','')') difExp
    write(n_out,'(''  nVarDiff    ='',i3,'','')') nVarDiff

    !----- Variable kinematic viscosity:
    write(n_out,'(''  nVarVisc    ='',i3,'','')') nVarVisc

    !----- Internal heating form:
    write(n_out,'(''  nVarEps     ='',i3,'','')') nVarEps

    !----- External field
    write(n_out,'(''  n_imp       ='',i3,'','')') n_imp
    write(n_out,'(1p,''  rrMP           ='',d14.6,'','')') rrMP
    write(n_out,'(1p,''  amp_imp        ='',d14.6,'','')') amp_imp
    write(n_out,'(1p,''  expo_imp       ='',d14.6,'','')') expo_imp
    write(n_out,'(1p,''  bmax_imp       ='',d14.6,'','')') bmax_imp

    write(n_out,*) "/"


    write(n_out,*) "&start_field"
    write(n_out,'(''  l_start_file='',l3,'','')') l_start_file
    length=length_to_blank(start_file)
    write(n_out,*) " start_file  = """,start_file(1:length),""","
    write(n_out,'(''  inform      ='',i3,'','')') inform
    write(n_out,'(''  l_reset_t   ='',l3,'','')') l_reset_t
    write(n_out,'(1p,''  scale_s     ='',d14.6,'','')') scale_s
    write(n_out,'(1p,''  scale_b     ='',d14.6,'','')') scale_b
    write(n_out,'(1p,''  scale_v     ='',d14.6,'','')') scale_v
    write(n_out,'(1p,''  tipdipole   ='',d14.6,'','')') tipdipole
    WRITE(n_out,'(''  init_s1     ='',i7,'','')') init_s1
    write(n_out,'(''  init_s2     ='',i3,'','')') init_s2
    write(n_out,'(''  init_v1     ='',i3,'','')') init_v1
    write(n_out,'(''  init_b1     ='',i3,'','')') init_b1
    write(n_out,'(''  imagcon     ='',i3,'','')') imagcon
    write(n_out,'(1p,''  amp_s1      ='',d14.6,'','')') amp_s1
    write(n_out,'(1p,''  amp_s2      ='',d14.6,'','')') amp_s2
    write(n_out,'(1p,''  amp_v1      ='',d14.6,'','')') amp_v1
    write(n_out,'(1p,''  amp_b1      ='',d14.6,'','')') amp_b1
    write(n_out,*) "/"

    write(n_out,*) "&output_control"
    write(n_out,'(''  n_graph_step  ='',i5,'','')') n_graph_step
    write(n_out,'(''  n_graphs      ='',i5,'','')') n_graphs
    write(n_out,'(1p,''  t_graph_start ='',d14.6,'','')') &
         t_graph_start
    write(n_out,'(1p,''  t_graph_stop  ='',d14.6,'','')') &
         t_graph_stop
    write(n_out,'(1p,''  dt_graph      ='',d14.6,'','')') &
         dt_graph
    write(n_out,'(''  n_rst_step    ='',i5,'','')') n_rst_step
    write(n_out,'(''  n_rsts        ='',i5,'','')') n_rsts
    write(n_out,'(1p,''  t_rst_start   ='',d14.6,'','')') &
         t_rst_start
    write(n_out,'(1p,''  t_rst_stop    ='',d14.6,'','')') &
         t_rst_stop
    write(n_out,'(1p,''  dt_rst        ='',d14.6,'','')') dt_rst
    write(n_out,'(''  n_stores      ='',i5,'','')') n_stores
    write(n_out,'(''  n_log_step    ='',i5,'','')') n_log_step
    write(n_out,'(''  n_logs        ='',i5,'','')') n_logs
    write(n_out,'(1p,''  t_log_start   ='',d14.6,'','')') &
         t_log_start
    write(n_out,'(1p,''  t_log_stop    ='',d14.6,'','')') &
         t_log_stop
    write(n_out,'(1p,''  dt_log        ='',d14.6,'','')') dt_log
    write(n_out,'(''  n_p_step      ='',i5,'','')') n_p_step
    write(n_out,'(''  n_ps          ='',i5,'','')') n_ps
    write(n_out,'(1p,''  t_p_start     ='',d14.6,'','')') &
         t_p_start
    write(n_out,'(1p,''  t_p_stop      ='',d14.6,'','')') &
         t_p_stop
    write(n_out,'(1p,''  dt_p          ='',d14.6,'','')') dt_p
    write(n_out,'(''  n_spec_step   ='',i5,'','')') n_spec_step
    write(n_out,'(''  n_specs       ='',i5,'','')') n_specs
    write(n_out,'(1p,''  t_spec_start  ='',d14.6,'','')') &
         t_spec_start
    write(n_out,'(1p,''  t_spec_stop   ='',d14.6,'','')') &
         t_spec_stop
    write(n_out,'(1p,''  dt_spec       ='',d14.6,'','')') dt_spec
    write(n_out,'(''  n_cmb_step    ='',i5,'','')') n_cmb_step
    write(n_out,'(''  n_cmbs        ='',i5,'','')') n_cmbs
    write(n_out,'(1p,''  t_cmb_start   ='',d14.6,'','')') t_cmb_start
    write(n_out,'(1p,''  t_cmb_stop    ='',d14.6,'','')') t_cmb_stop
    write(n_out,'(1p,''  dt_cmb        ='',d14.6,'','')') dt_cmb
    write(n_out,'(''  n_r_field_step   ='',i5,'','')') n_r_field_step
    write(n_out,'(''  n_r_fields       ='',i5,'','')') n_r_fields
    write(n_out,'(1p,''  t_r_field_start='',d14.6,'','')') t_r_field_start
    write(n_out,'(1p,''  t_r_field_stop ='',d14.6,'','')') t_r_field_stop
    write(n_out,'(1p,''  dt_r_field    ='',d14.6,'','')') dt_r_field
    write(n_out,'(''  l_movie       ='',l3,'','')') l_movie
    write(n_out,'(''  n_movie_step  ='',i5,'','')') n_movie_step
    write(n_out,'(''  n_movie_frames='',i5,'','')') n_movie_frames
    write(n_out,'(1p,''  t_movie_start ='',d14.6,'','')') &
         t_movie_start
    write(n_out,'(1p,''  t_movie_stop  ='',d14.6,'','')') &
         t_movie_stop
    write(n_out,'(1p,''  dt_movie      ='',d14.6,'','')') &
         dt_movie
    do n=1,n_movies_max
       length=len(trim(movie(n)))
       if ( length > 0 ) &
            write(n_out,'(''  movie         = '',a,'','')') &
            movie(n)(1:length)
    end do
    write(n_out,'(''  ngform        ='',i3,'','')') ngform
    write(n_out,'(''  l_average     ='',l3,'','')') l_average
    write(n_out,'(''  l_cmb_field   ='',l3,'','')') l_cmb_field
    write(n_out,'(''  l_dt_cmb_field='',l3,'','')') l_dt_cmb_field
    write(n_out,'(''  l_save_out    ='',l3,'','')') l_save_out
    write(n_out,'(''  l_true_time   ='',l3,'','')') l_true_time
    write(n_out,'(''  lVerbose      ='',l3,'','')') lVerbose
    write(n_out,'(''  l_rMagSpec    ='',l3,'','')') l_rMagSpec
    write(n_out,'(''  l_DTrMagSpec  ='',l3,'','')') l_DTrMagSpec
    write(n_out,'(''  l_max_cmb     ='',i3,'','')') l_max_cmb
    write(n_out,'(''  l_r_field     ='',l3,'','')') l_r_field
    write(n_out,'(''  l_r_fieldT    ='',l3,'','')') l_r_fieldT
    write(n_out,'(''  l_max_r       ='',i3,'','')') l_max_r
    write(n_out,'(''  n_r_step      ='',i3,'','')') n_r_step
    DO n=1,n_coeff_r_max
       write(n_out,'(''    n_coeff_r   ='',i3,'','')') n_coeff_r(n)
    END DO
    write(n_out,'(''  l_hel         ='',l3,'','')') l_hel
    write(n_out,'(''  l_AM          ='',l3,'','')') l_AM
    write(n_out,'(''  l_power       ='',l3,'','')') l_power
    write(n_out,'(''  l_viscBcCalc  ='',l3,'','')') l_viscBcCalc
    write(n_out,'(''  l_fluxProfs   ='',l3,'','')') l_fluxProfs
    write(n_out,'(''  l_drift       ='',l3,'','')') l_drift
    write(n_out,'(''  l_iner        ='',l3,'','')') l_iner
    write(n_out,'(''  l_TO          ='',l3,'','')') l_TO
    write(n_out,'(''  l_TOmovie     ='',l3,'','')') l_TOmovie
    write(n_out,'(''  l_PV          ='',l3,'','')') l_PV
    write(n_out,'(''  l_storeBpot   ='',l3,'','')') l_storeBpot
    write(n_out,'(''  l_storeVpot   ='',l3,'','')') l_storeVpot
    write(n_out,'(''  l_RMS         ='',l3,'','')') l_RMS
    write(n_out,'(''  l_RMStest     ='',l3,'','')') l_RMStest
    write(n_out,'(''  l_prms        ='',l3,'','')') l_prms
    write(n_out,'(''  l_par         ='',l3,'','')') l_par
    write(n_out,'(''  l_corrMov     ='',l3,'','')') l_corrMov
    write(n_out,'(1p,''  rCut          ='',d14.6,'','')') rCut
    write(n_out,'(1p,''  rDea          ='',d14.6,'','')') rDea
    write(n_out,*) "/"

    write(n_out,*) "&mantle"
    write(n_out,'(1p,''  conductance_ma='',d14.6,'','')') &
         conductance_ma
    write(n_out,'(1p,''  rho_ratio_ma  ='',d14.6,'','')') &
         rho_ratio_ma
    write(n_out,'(1p,''  nRotMa        ='',i4,'','')') nRotMa
    write(n_out,'(1p,''  omega_ma1     ='',d14.6,'','')') &
         omega_ma1
    write(n_out,'(1p,''  omegaOsz_ma1  ='',d14.6,'','')') &
         omegaOsz_ma1
    write(n_out,'(1p,''  tShift_ma1    ='',d14.6,'','')') &
         tShift_ma1
    write(n_out,'(1p,''  omega_ma2     ='',d14.6,'','')') &
         omega_ma2
    write(n_out,'(1p,''  omegaOsz_ma2  ='',d14.6,'','')') &
         omegaOsz_ma2
    write(n_out,'(1p,''  tShift_ma2    ='',d14.6,'','')') &
         tShift_ma2
    write(n_out,*) "/"

    write(n_out,*) "&inner_core"
    write(n_out,'(1p,''  sigma_ratio   ='',d14.6,'','')') &
         sigma_ratio
    write(n_out,'(1p,''  rho_ratio_ic  ='',d14.6,'','')') &
         rho_ratio_ic
    write(n_out,'(1p,''  nRotIc        ='',i4,'','')') nRotIc
    write(n_out,'(1p,''  omega_ic1     ='',d14.6,'','')') &
         omega_ic1
    write(n_out,'(1p,''  omegaOsz_ic1  ='',d14.6,'','')') &
         omegaOsz_ic1
    write(n_out,'(1p,''  tShift_ic1    ='',d14.6,'','')') &
         tShift_ic1
    write(n_out,'(1p,''  omega_ic2     ='',d14.6,'','')') &
         omega_ic2
    write(n_out,'(1p,''  omegaOsz_ic2  ='',d14.6,'','')') &
         omegaOsz_ic2
    write(n_out,'(1p,''  tShift_ic2    ='',d14.6,'','')') &
         tShift_ic2
    write(n_out,'(1p,''  BIC           ='',d14.6,'','')') &
         BIC
    write(n_out,*) "/"


    RETURN
  end SUBROUTINE writeNamelists

  !----------------------------------------------------------------------
  !  +-------------+----------------+------------------------------------+
  !  |                                                                   |
  !  |  Purpose of this subroutine is to set default parameters          |
  !  |  for the namelists.                                               |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+
  SUBROUTINE defaultNamelists

    !-- local:
    INTEGER n

    !-- end of declaration
    !----------------------------------------------------------------------

    !----- Namelist grid
    ! must be of form 4*integer+1
    ! Possible values for n_r_max:
    !  5,9,13,17,21,25,33,37,41,49,
    !  61,65,73,81,97,101,121,129, ...
    n_r_max       =33
    ! max degree-1 of cheb polynomia
    n_cheb_max    =31
    ! number of longitude grid points
    ! Possible valus: 
    ! 16,32,48,64,96,128,192,256,288,320,384,
    ! 400,512,576,640,768,864,1024
    n_phi_tot     =192
    ! number of grid points in inner core
    n_r_ic_max    =17
    ! number of chebs in inner core
    n_cheb_ic_max =15
    ! basic wavenumber, longitude symmetry
    minc          =1
    ! controls dealiasing in latitude and 
    ! longitude direction, no aliasing for nalias=20
    !   20 <= nalias <= 30
    nalias        =20
    !----- Namelist contrl
    mode          =0            ! self consistent dynamo !
    tag           ="default"
    n_time_steps  =100
    n_tScale      =0
    n_lScale      =0
    alpha         =0.5D0
    enscale       =1.0D0
    dtstart       =0.0D0
    dtMax         =1.0D-4
    courfac       =2.5D0
    alffac        =1.0D0
    intfac        =0.15D0
    n_cour_step   =10

    cacheblock_size_in_B=4096

    l_update_v    =.TRUE.
    l_update_b    =.TRUE.
    l_update_s    =.TRUE.
    l_correct_AMe =.FALSE.  ! Correct equatorial AM
    l_correct_AMz =.FALSE.  ! Correct axial AM
    l_non_rot     =.FALSE.  ! No Coriolis force !
    l_anel        =.FALSE.  ! Anelastic stuff !
    l_isothermal  =.FALSE.  ! Isothermal = 0 Grünesein !
    l_interior_model=.FALSE.! Interior model

    !---- Run time and number of threads:
    l_runTimeLimit=.FALSE. ! Control of absolute run time
    DO n=1,4
       runTimeLimit(n)=0
    END DO

    !nThreadsRun   =1
    tEND          =0.D0    ! numerical time where run should end

    !----- Hyperdiffusion:
    difnu         =0.D0
    difeta        =0.D0
    difkap        =0.D0
    ldif          =1
    ldifexp       =-1

    !----- Namelist phys_param:
    ra         =1.1D5
    ek         =1.D-3
    pr         =1.D0
    prmag      =5.D0
    epsc0      =0.D0
    radratio   =0.35D0
    !----- Anelatic stuff
    strat      =0.D0
    polind     =1.5D0
    r_cut_model=0.98D0 ! outer radius when using interior model
    !----- Gravity parameters: defaut value g propto r (i.e. g1=1)
    g0         =0.D0
    g1         =1.D0
    g2         =0.D0        
    !----- Boundary conditions        
    ktops      =1
    kbots      =1
    ktopv      =2
    kbotv      =2
    ktopb      =1
    kbotb      =1
    do n=1,4*n_s_bounds
       s_top(n)=0.D0
       s_bot(n)=0.D0
    end do
    impS=0
    DO n=1,n_impS_max
       peakS(n) =0.D0
       thetaS(n)=0.D0
       phiS(n)  =0.D0
       widthS(n)=0.D0
    END DO

    !----- Conductivity variation:
    nVarCond       =0
    con_DecRate    =9
    con_RadRatio   =0.75D0
    con_LambdaMatch=0.6D0
    con_LambdaOut  =0.1
    con_FuncWidth  =0.25D0

    !----- Thermal diffusivity variation:
    difExp         =-0.5D0
    nVarDiff       =0

    !----- Variable kinematic viscosity:
    nVarVisc       =0

    !----- Internal heating form:
    nVarEps        =0

    !----- Non-linear mapping parameters (Bayliss, 1992):
    l_newmap       =.FALSE.
    l_plotmap      =.FALSE.
    alph1          =2.D0
    alph2          =0.D0

    !----- External field
    n_imp          =0    ! No external field
    rrMP           =0.D0 ! r(Magnetopause)/r_cmb, used for n_imp=1
    amp_imp        =0.D0 ! amplitude of external field
    expo_imp       =0.D0 ! oscillation frequency of external field
    bmax_imp       =0.D0

    !----- Namelist start_field:
    l_start_file  =.FALSE.
    start_file    ="no_start_file"
    inform        =-1   
    runid         ="MAGIC default run"
    l_reset_t     =.false.
    scale_s       =1.D0
    scale_b       =1.D0
    scale_v       =1.D0
    tipdipole     =0.D0
    init_s1       =0
    init_s2       =0
    init_b1       =0
    init_v1       =0
    imagcon       =0
    tmagcon       =0.D0
    amp_s1        =1.D0
    amp_s2        =0.D0
    amp_v1        =0.D0
    amp_b1        =1.D0


    !----- Namelist output_control:
    l_save_out    =.FALSE.  ! Save output
    l_true_time   =.FALSE.  ! Use exact requested output times
    lVerbose      =.FALSE.  ! Tell me what you are doing
    l_average     =.FALSE.  ! Average various quantities in time

    !----- Restart files:
    n_rst_step    =0
    n_rsts        =1
    t_rst_start   =0.D0
    t_rst_stop    =0.D0
    dt_rst        =0.D0
    n_stores      =0

    !----- Log output:
    n_log_step    =50
    n_logs        =0
    t_log_start   =0.D0
    t_log_stop    =0.D0
    dt_log        =0.D0

    !----- Graphic output:
    n_graph_step  =0
    n_graphs      =1
    t_graph_start =0.D0
    t_graph_stop  =0.D0
    dt_graph      =0.D0
    ngform        =0
    l_graph_time  =.FALSE.

    !----- Spectrum files:
    n_spec_step   =0
    n_specs       =0
    t_spec_start  =0.D0
    t_spec_stop   =0.D0
    dt_spec       =0.D0

    !----- Output of poloidal magnetic field potential at CMB:
    !      also stored at times of movie frames
    l_cmb_field   =.FALSE.
    l_dt_cmb_field=.FALSE.
    l_max_cmb     =14
    n_cmb_step    =0
    n_cmbs        =0
    t_cmb_start   =0.D0
    t_cmb_stop    =0.D0
    dt_cmb        =0.D0

    !----- Output of magnetic and flow potential af five different radial levels:
    l_r_field     =.FALSE.
    l_r_fieldT    =.FALSE.
    l_max_r       =l_max
    n_r_step      =2
    DO n=1,n_coeff_r_max
       n_r_array(n)=0
    END DO
    n_r_field_step=0
    n_r_fields    =0
    t_r_field_start=0.D0
    t_r_field_stop =0.D0
    dt_r_field    =0.D0

    !----- Movie output:
    l_movie       =.FALSE.
    n_movies      =0
    n_movie_step  =0
    n_movie_frames=0
    t_movie_start =0.D0
    t_movie_stop  =0.D0
    dt_movie      =0.D0
    DO n=1,n_movies_max
       movie(n)=' '
    END DO

    !----- Output of magnetic potentials:
    l_storeBpot   =.FALSE.
    n_Bpot_step   =0
    n_Bpots       =0
    t_Bpot_start  =0.D0
    t_Bpot_stop   =0.D0
    dt_Bpot       =0.D0

    !----- Output of flow potentials:
    l_storeVpot   =.FALSE.
    n_Vpot_step   =0
    n_Vpots       =0
    t_Vpot_start  =0.D0
    t_Vpot_stop   =0.D0
    dt_Vpot       =0.D0

    !----- Output of T potential:
    l_storeTpot   =.FALSE.
    n_Tpot_step   =0
    n_Tpots       =0
    t_Tpot_start  =0.D0
    t_Tpot_stop   =0.D0
    dt_Tpot       =0.D0

    !----- Output of all potential:
    l_storePot    =.FALSE.
    n_pot_step    =0
    n_pots        =0
    t_pot_start   =0.D0
    t_pot_stop    =0.D0
    dt_pot        =0.D0

    !----- Output TOZ:
    n_TOZ_step    =0
    n_TOZs        =0
    t_TOZ_start   =0.D0
    t_TOZ_stop    =0.D0
    dt_TOZ        =0.D0

    !----- Times for different output:
    DO n=1,n_time_hits
       t_graph(n)  =-1.D0
       t_rst(n)    =-1.D0
       t_log(n)    =-1.D0
       t_p(n)      =-1.D0
       t_spec(n)   =-1.D0
       t_cmb(n)    =-1.D0
       t_r_field(n)=-1.D0
       t_movie(n)  =-1.D0
       t_Vpot(n)   =-1.D0
       t_Bpot(n)   =-1.D0
       t_Tpot(n)   =-1.D0
       t_TO(n)     =-1.D0
       t_TOZ(n)    =-1.D0
       t_TOmovie(n)=-1.D0
    END DO

    !----- Magnetic spectra for different depths
    !      at times of log output or movie frames:
    l_rMagSpec    =.FALSE.
    l_DTrMagSpec  =.FALSE.

    !----- TO output, output times same as for log outout:
    l_TO          =.FALSE. ! TO output in TOnhs.TAG, TOshs.TAG
    l_TOmovie     =.FALSE. ! TO movies 
    sDens         =1.D0    ! relative s-grid point density 
    zDens         =1.D0    ! relative z-grid point dendity 

    !----- Potential vortivity:
    l_PV          =.FALSE.

    !----- Different output, output times same as for log outout:
    l_hel         =.FALSE. ! Helicity in misc.TAG 
    l_AM          =.FALSE. ! Angular moment in AM.TAG 
    l_power       =.FALSE. ! power budget in power.TAG and dtE.TAG
    l_viscBcCalc  =.FALSE. ! dissipation layer for stress-free BCs
    l_fluxProfs   =.FALSE. ! radial profiles of flux contributions
    l_drift       =.FALSE. ! files for calculating drift rates 
    l_iner        =.FALSE. ! files for calculating inertial modes
    l_RMS         =.FALSE. ! RMS force ballance and dynamo term 
    ! ballance in dtVrms.TAG and dtBrms.TAG
    l_RMStest     =.FALSE. ! special test for RMS balance
    l_prms        =.FALSE. ! radial plot of the rms pressure (snapshot)
    l_par         =.FALSE. ! Calculate additional parameters in s_getEgeos.f
    l_corrMov     =.FALSE. ! North/south correlation movie (see s_getEgeos.f)
    rCut          =0.075D0 ! Thickness of layer to be left out at both
    ! boundaries for RMS calculation.
    ! rCut=0.075 means that 7.5% at the CMB and ICB are disregarded.
    rDea          =0.00D0  ! Controls dealiazing in  RMS calculation
    ! rDea=0.1 means that highest 10% of cheb modes are set to zero

    !----- Mantle name list:
    conductance_ma=0.D0    ! insulation mantle is default
    nRotMa        =0       ! non rotating mantle is default
    rho_ratio_ma  =1.D0    ! same density as outer core
    omega_ma1     =0.D0    ! prescribed rotation rate
    omegaOsz_ma1  =0.D0    ! oszillation frequency of mantle rotation rate
    tShift_ma1    =0.D0    ! time shift
    omega_ma2     =0.D0    ! second mantle rotation rate 
    omegaOsz_ma2  =0.D0    ! oscillation frequency of second mantle rotation
    tShift_ma2    =0.D0    ! time shift for second rotation

    !----- Inner core name list:
    sigma_ratio   =0.D0    ! no conducting inner core is default 
    nRotIc        =0       ! non rotating inner core is default
    rho_ratio_ic  =1.D0    ! same density as outer core
    omega_ic1     =0.D0    ! prescribed rotation rate, added to first one
    omegaOsz_ic1  =0.D0    ! oszillation frequency of IC rotation rate
    tShift_ic1    =0.D0    ! time shift
    omega_ic2     =0.D0    ! second prescribed rotation rate 
    omegaOsz_ic2  =0.D0    ! oszillation frequency of second IC rotation rate
    tShift_ic2    =0.D0    ! tims shift for second IC rotation
    BIC           =0.D0    ! Imposed dipole field strength at ICB


    RETURN
  END SUBROUTINE defaultNamelists

END MODULE Namelists
