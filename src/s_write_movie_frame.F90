!$Id$
!*************************************************************************
    SUBROUTINE write_movie_frame(n_frame,time, &
            b,db,aj,dj,b_ic,db_ic,aj_ic,dj_ic, &
                            omega_ic,omega_ma)
!*************************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!-------------------------------------------------------------------------

!  Writes different movie frames into respective output files.

!-------------------------------------------------------------------------

    USE truncation
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE blocking
    USE horizontal_data
    USE logic
    USE movie_data
    USE output_data

    IMPLICIT NONE

!-- Input of constant parameters:
! include 'truncation.f'
! include 'c_blocking.f'
! include 'c_radial.f'
! include 'c_horizontal.f'
! include 'c_phys_param.f'
! include 'c_num_param.f'
! include 'c_output.f'
! include 'c_logic.f'

!-- Input of variables:
    REAL(kind=8) :: time
    INTEGER :: n_frame
    REAL(kind=8) :: omega_ic,omega_ma

!-- Input of scalar fields:
    COMPLEX(kind=8) :: b(lm_maxMag,n_r_maxMag)
    COMPLEX(kind=8) :: db(lm_maxMag,n_r_maxMag)
    COMPLEX(kind=8) :: aj(lm_maxMag,n_r_maxMag)
    COMPLEX(kind=8) :: dj(lm_maxMag,n_r_maxMag)
    COMPLEX(kind=8) :: b_ic(lm_maxMag,n_r_ic_maxMag)
    COMPLEX(kind=8) :: db_ic(lm_maxMag,n_r_ic_maxMag)
    COMPLEX(kind=8) :: aj_ic(lm_maxMag,n_r_ic_maxMag)
    COMPLEX(kind=8) :: dj_ic(lm_maxMag,n_r_ic_maxMag)

!-- Input of frame to be written:
! include 'c_movie.f'

!-- Local:
    INTEGER :: n_fields
    INTEGER :: n_fields_ic
    INTEGER :: n_fields_oc
    INTEGER :: n_movie
    INTEGER :: n_surface
    INTEGER :: n_type
    INTEGER :: n_out
    INTEGER :: n_field,n,n_start,n_stop
    INTEGER :: n_r,n_theta,n_phi
    INTEGER :: n_r_mov_tot

    CHARACTER(len=64) :: version
            
    REAL(kind=8) :: const
    REAL(kind=8) :: r_mov_tot(n_r_max+n_r_ic_max)
    REAL(kind=4) :: dumm(n_theta_max)


!-- end of declaration
!----------------------------------------------------------------------

    t_movieS(n_frame)=time

    DO n_movie=1,n_movies

        n_type      =n_movie_type(n_movie)
        n_surface   =n_movie_surface(n_movie)
        n_fields_oc =n_movie_fields(n_movie)
        n_fields_ic =n_movie_fields_ic(n_movie)
        n_fields    =max0(n_fields_ic,n_fields_oc)
        n_out       =n_movie_file(n_movie)
        const       =movie_const(n_movie)
        IF ( n_surface == 1 ) const=const/r_cmb

    !------ Open movie file:
        IF ( l_save_out ) THEN
            OPEN(n_out,FILE=movie_file(n_movie),STATUS='UNKNOWN', &
                 FORM='UNFORMATTED',POSITION='APPEND')
        END IF

    !------ Write header if this is the first frame:
         
        IF ( n_frame == 1 ) then

        !------ Start with info about movie type:
            version='JW_Movie_Version_2'
            WRITE(n_out) version
            WRITE(n_out) FLOAT(n_type),FLOAT(n_surface), &
                             SNGL(const),FLOAT(n_fields)
            WRITE(n_out) &
                 (FLOAT(n_movie_field_type(n,n_movie)),n=1,n_fields)


        !------ Combine OC and IC radial grid points:
            n_r_mov_tot=n_r_max
            DO n_r=1,n_r_max
                r_mov_tot(n_r)=r(n_r)
            END DO
            IF ( n_r_ic_max > 0 ) THEN
                n_r_mov_tot=n_r_mov_tot+n_r_ic_max-2
                DO n_r=1,n_r_ic_max-2
                    r_mov_tot(n_r_max+n_r)=r_ic(n_r+1)
                END DO
            END IF

        !------ Now other info about grid and parameters:
            WRITE(n_out) runid          ! run identifyer (as set in namelist contrl)
            dumm( 1)=FLOAT(n_r_mov_tot)
            dumm( 2)=FLOAT(n_r_max)
            dumm( 3)=FLOAT(n_theta_max) ! no. of theta points
            dumm( 4)=FLOAT(n_phi_max)   ! no. of phi points
            dumm( 5)=FLOAT(minc)            ! imposed symmetry
            dumm( 6)=SNGL(ra)               ! control parameters
            dumm( 7)=SNGL(ek)               ! (for information only)
            dumm( 8)=SNGL(pr)               !      -"-
            dumm( 9)=SNGL(prmag)            !      -"-
            dumm(10)=SNGL(radratio)         ! ratio of inner / outer core
            dumm(11)=SNGL(tScale)           ! timescale
            WRITE(n_out) (dumm(n),n=1,11)

        !------ Write grid:
            WRITE(n_out) (SNGL(r_mov_tot(n_r)/r_cmb), &
                         n_r=1,n_r_mov_tot)
            WRITE(n_out) (SNGL(theta_ord(n_theta)), &
                         n_theta=1,n_theta_max)
            WRITE(n_out) (SNGL(phi(n_phi)), &
                         n_phi=1,n_phi_max)

        END IF  ! Write header ?
                   

    !------ Write frame number, time and IC and MA rotation rates::
        dumm(1)=FLOAT(n_frame)
        dumm(2)=SNGL(t_movieS(n_frame))
        dumm(3)=SNGL(omega_ic)
        dumm(4)=SNGL(omega_ma)
        dumm(5)=SNGL(movieDipColat)
        dumm(6)=SNGL(movieDipLon)
        dumm(7)=SNGL(movieDipStrength)
        dumm(8)=SNGL(movieDipStrengthGeo)
        WRITE(n_out) (dumm(n),n=1,8)

    !------ Write frames:
        IF ( .NOT. lStoreMov(n_movie) ) THEN
            IF ( n_type == 99 ) THEN
                WRITE(*,*) '! Use TO output for Lorentz force!'
                STOP
            ELSE
                CALL write_dtB_frame(n_movie,b,db,aj,dj, &
                                     b_ic,db_ic,aj_ic,dj_ic)
            END IF
        ELSE
            DO n_field=1,n_fields
                n_start=n_movie_field_start(n_field,n_movie)
                IF ( n_fields_oc > 0 ) THEN
                    n_stop =n_movie_field_stop(n_field,n_movie)
                ENDIF
                IF ( n_fields_ic > 0 ) THEN
                    n_stop=n_movie_field_stop(n_field+n_fields_oc,n_movie)
                ENDIF
                WRITE(n_out) (SNGL(frames(n)),n=n_start,n_stop)
            END DO
        END IF


        IF ( l_save_out ) CLOSE(n_out)


    END DO  ! Loop over movies


    RETURN
    end SUBROUTINE write_movie_frame

!--------------------------------------------------------------------------
