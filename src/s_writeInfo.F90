!$Id$
!***********************************************************************
    SUBROUTINE writeInfo(n_out)
!***********************************************************************

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  |  Purpose of this subroutine is to write the namelist to           |
!  |  file unit n_out. This file has to be open before calling this    |
!  |  routine.                                                         |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+

    USE truncation
    USE physical_parameters
    USE num_param
    USE init_fields
    USE logic
    USE output_data
    USE const
    USE parallel_mod,only:rank

    IMPLICIT NONE

!-- input:
    INTEGER :: n_out

!-- input via common-blocks:
! include 'truncation.f'
! include 'c_phys_param.f'
! include 'c_num_param.f'
! include 'c_logic.f'
! include 'c_init_fields.f'
! include 'c_output.f'
! include 'c_const.f'

!-- end of declaration
!----------------------------------------------------------------------

    IF (rank.EQ.0) THEN
       IF ( l_save_out .AND. n_out == n_log_file ) THEN
          OPEN(n_out,FILE=log_file,STATUS='UNKNOWN',POSITION='APPEND')
       END IF

       !-- Output of mode:
       WRITE(n_out,*)
       IF ( mode == 0 ) THEN
          WRITE(n_out,*) &
               '! Self consistent dynamo integration.'
       ELSE IF ( mode == 1 ) then
          WRITE(n_out,*) &
               '! Convection integration.'
       ELSE IF ( mode == 2 ) then
          WRITE(n_out,*) &
               '! Kinematic dynamo integration.'
       ELSE IF ( mode == 3 ) then
          WRITE(n_out,*) &
               '! Magnetic decay modes.'
       ELSE IF ( mode == 4 ) then
          WRITE(n_out,*) &
               '! Magneto convection.'
       ELSE IF ( mode == 5 ) then
          WRITE(n_out,*) &
               '! Linear onset of convection.'
       ELSE IF ( mode == 6 ) then
          WRITE(n_out,*) &
               '! Self consistent dynamo integration without LF.'
       ELSE IF ( mode == 7 ) then
          WRITE(n_out,*) &
               '! Super-rotating IC, no convection, no dynamo.'
       ELSE IF ( mode == 8 ) then
          WRITE(n_out,*) &
               '! Super-rotating IC, no convection, dynamo.'
       ELSE IF ( mode == 9 ) then
          WRITE(n_out,*) &
               '! Super-rotating IC, no convection, dynamo, no LF.'
       ELSE IF ( mode == 10 ) then
          WRITE(n_out,*) &
               '! Super-rotating IC, no advection, no convection, no dynamo.'
       ELSE IF ( mode == 11 ) then
          WRITE(n_out,*) &
               '! Viscous flow, no inertia, no rotation, no dynamo.'
       END IF
    END IF
    IF (mode.GT.11) THEN
       if (rank.eq.0) WRITE(*,'(/," Mode > 11 not implemented ! ")')
       STOP
    END IF

    IF (rank.EQ.0) THEN
       !-- Output of name lists:
       WRITE(n_out, &
            &   '(/,'' ! Normalized OC moment of inertia:'',d14.6)') c_moi_oc
       WRITE(n_out, &
            &   '('' ! Normalized IC moment of inertia:'',d14.6)') c_moi_ic
       WRITE(n_out, &
            &   '('' ! Normalized MA moment of inertia:'',d14.6)') c_moi_ma
       WRITE(n_out, &
            &   '('' ! Normalized IC volume :'',d14.6)') vol_ic
       WRITE(n_out, &
            &   '('' ! Normalized OC volume :'',d14.6)') vol_oc
       WRITE(n_out, &
            &   '('' ! Normalized IC surface:'',d14.6)') surf_cmb*radratio**2
       WRITE(n_out, &
            &   '('' ! Normalized OC surface:'',d14.6)') surf_cmb
       WRITE(n_out,*)
       WRITE(n_out,*) '! Grid parameters:'
       WRITE(n_out,'(''  n_r_max      ='',i6, &
            &   '' = number of radial grid points'')') n_r_max
       WRITE(n_out,'(''  n_cheb_max   ='',i6)') n_cheb_max
       WRITE(n_out,'(''  max cheb deg.='',i6)') n_cheb_max-1
       WRITE(n_out,'(''  n_phi_max    ='',i6, &
            &   '' = no of longitude grid points'')') n_phi_max
       WRITE(n_out,'(''  n_theta_max  ='',i6, &
            &   '' = no of latitude grid points'')') n_theta_max
       WRITE(n_out,'(''  n_r_ic_max   ='',i6, &
            &   '' = number of radial grid points in IC'')') n_r_ic_max
       WRITE(n_out,'(''  n_cheb_ic_max='',i6)') &
            n_cheb_ic_max-1
       WRITE(n_out,'(''  max cheb deg ='',i6)') &
            2*(n_cheb_ic_max-1)
       WRITE(n_out,'(''  l_max        ='',i6, &
            &   '' = max degree of Plm'')') l_max
       WRITE(n_out,'(''  m_max        ='',i6, &
            &   '' = max oder of Plm'')') m_max
       WRITE(n_out,'(''  lm_max       ='',i6, &
            &   '' = no of l/m combinations'')') lm_max
       WRITE(n_out,'(''  minc         ='',i6, &
            &   '' = longitude symmetry wave no'')') minc
       WRITE(n_out,'(''  nalias       ='',i6, &
            &   '' = spher. harm. deal. factor '')') nalias

       IF ( l_save_out .AND. n_out == n_log_file ) &
            CLOSE(n_out)
    END IF

    RETURN
    end SUBROUTINE writeInfo

!----------------------------------------------------------------------
