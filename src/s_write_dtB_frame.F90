!$Id$
!*************************************************************************
    SUBROUTINE write_dtB_frame(n_movie,b,db,aj,dj, &
                               b_ic,db_ic,aj_ic,dj_ic)
!*************************************************************************

!------------ This is release 2 level 1  --------------!
!------------ Created on 1/17/02  by JW. --------------!

!-------------------------------------------------------------------------

!  Controlls output of movie frames.

!-------------------------------------------------------------------------

    USE truncation
    USE radial_functions
    USE blocking
    USE horizontal_data
    USE dtB_mod
    USE output_data

    IMPLICIT NONE

!-- Input of constant parameters:
! include 'truncation.f'
! include 'c_output.f'
! include 'c_blocking.f'
! include 'c_radial.f'
! include 'c_horizontal.f'

!-- Input of variables:
    integer :: n_movie

!-- Input of scalar fields:
    COMPLEX(kind=8) :: b(lm_maxMag,n_r_maxMag)
    COMPLEX(kind=8) :: db(lm_maxMag,n_r_maxMag)
    COMPLEX(kind=8) :: aj(lm_maxMag,n_r_maxMag)
    COMPLEX(kind=8) :: dj(lm_maxMag,n_r_maxMag)
    COMPLEX(kind=8) :: b_ic(lm_maxMag,n_r_ic_maxMag)
    COMPLEX(kind=8) :: db_ic(lm_maxMag,n_r_ic_maxMag)
    COMPLEX(kind=8) :: aj_ic(lm_maxMag,n_r_ic_maxMag)
    COMPLEX(kind=8) :: dj_ic(lm_maxMag,n_r_ic_maxMag)

!-- Input of magnetic field generation terms:
!   Parallelization note: see s_get_dtBLMfinish.f
! include 'c_dtB.f'

!-- Output: writes movie frames

!-- Local:
    integer :: n_type
    integer :: n_out
    integer :: n_fields_oc
    integer :: n_fields_ic
    integer :: n_fields,n_field
    integer :: n_const

    integer :: n_r,n_rC,n_r_loop_max
    integer :: n_theta,n_theta_start,n_theta_block,n_theta_cal
    integer :: n_phi
    integer :: n,n_pos,n_or
    integer :: lm,l,m

    integer :: n_surface
    integer :: n_field_type
    integer :: n_field_size

    real(kind=8) :: dtB(n_theta_max)
    real(kind=8) :: dtBframe(n_r_tot*n_theta_max)
    real(kind=8) :: dtBr(nrp,nfs)
    real(kind=8) :: dtBt(nrp,nfs)
    real(kind=8) :: dtBp(nrp,nfs)
    real(kind=8) :: dtBrframe(n_r_tot*n_phi_max*n_theta_max)
    real(kind=8) :: const,rMov

    complex(kind=8) :: workA(lm_max_dtB,n_r_max_dtB)
    complex(kind=8) :: workB(lm_max_dtB,n_r_max_dtB)

    logical :: l_loop


!-- end of declaration
!----------------------------------------------------------------------


    n_type      =n_movie_type(n_movie)
    n_fields_oc =n_movie_fields(n_movie)
    n_fields_ic =n_movie_fields_ic(n_movie)
    n_fields    =max0(n_fields_ic,n_fields_oc)
    n_out       =n_movie_file(n_movie)
    n_const     =n_movie_const(n_movie)
    n_surface   =n_movie_surface(n_movie)
    const       =movie_const(n_movie)

!--- Axisymmetric dtFL or dtAB:

    if ( n_type == 31 .OR. &
         n_type == 32 .OR. n_type == 33 .OR. &
         n_type == 41 .OR. &
         n_type == 42 .OR. n_type == 43 ) then

        n_field_size=n_theta_max*(n_r_max+n_r_ic_max-2)

        do n_field=1,n_fields

            n_field_type=n_movie_field_type(n_field,n_movie)

        !--- OC contribution:
            do n_r=1,n_r_max

                if ( n_field_type == 20 ) then
                    call get_dtB(dtB,PstrLM,lm_max,n_r_max, &
                                 n_r,1,n_theta_max,.false.)
                else if ( n_field_type == 21 ) then
                    call get_dtB(dtB,PadvLM,lm_max,n_r_max, &
                                 n_r,1,n_theta_max,.false.)
                else if ( n_field_type == 22 ) then
                    call get_dtB(dtB,PdifLM,lm_max,n_r_max, &
                                 n_r,1,n_theta_max,.false.)
                else if ( n_field_type == 23 ) then
                    call get_dtB(dtB,TstrLM,lm_max,n_r_max, &
                                 n_r,1,n_theta_max,.false.)
                else if ( n_field_type == 25 ) then
                    call get_dtB(dtB,TadvLM,lm_max,n_r_max, &
                                 n_r,1,n_theta_max,.false.)
                else if ( n_field_type == 26 ) then
                    call get_dtB(dtB,TdifLM,lm_max,n_r_max, &
                                 n_r,1,n_theta_max,.false.)
                end if

            !--- Store in frame field:
                DO n_theta_cal=1,n_theta_max
                    n_theta=n_theta_cal2ord(n_theta_cal)
                    n_pos  =(n_r-1)*n_theta_max+n_theta
                    dtBframe(n_pos)=dtB(n_theta)
                END DO

            END DO

        !--- Now IC contribution:
            DO n_r=2,n_r_ic_max-1

                if ( n_field_type == 21 ) then
                    call get_dtB(dtB,PadvLMIC,lm_max,n_r_ic_max, &
                                 n_r,1,n_theta_max,.TRUE.)
                else if ( n_field_type == 22 ) then
                    call get_dtB(dtB,PdifLMIC,lm_max,n_r_ic_max, &
                                 n_r,1,n_theta_max,.TRUE.)
                else if ( n_field_type == 25 ) then
                    call get_dtB(dtB,TadvLMIC,lm_max,n_r_ic_max, &
                                 n_r,1,n_theta_max,.TRUE.)
                else if ( n_field_type == 26 ) then
                    call get_dtB(dtB,TdifLMIC,lm_max,n_r_ic_max, &
                                 n_r,1,n_theta_max,.TRUE.)
                else
                    do n_theta=1,n_theta_max
                        dtB(n_theta)=0.d0
                    end do
                end if
                                 
            !--- Store in frame field:
                DO n_theta_cal=1,n_theta_max
                    n_theta=n_theta_cal2ord(n_theta_cal)
                    n_pos  =(n_r_max+n_r-2)*n_theta_max+n_theta
                    dtBframe(n_pos)=dtB(n_theta)
                END DO

            END DO

        !--- Write frame field:
            WRITE(n_out) (REAL(dtBframe(n),4),n=1,n_field_size)

        END DO

    !--- dtBr:
    ELSE ! non-axisymmetric fields

        DO n_field=1,n_fields

            n_field_type=n_movie_field_type(n_field,n_movie)

            IF ( n_field_type == 24 ) THEN
            ! This reduces omega effect to field production of axisymm. toroidal field:
                DO n_r=1,n_r_max
                    DO lm=1,lm_max
                        m=lm2m(lm)
                        IF ( m == 0 ) THEN
                            workA(lm,n_r)=TomeLM(lm,n_r)
                        ELSE
                            workA(lm,n_r)=CMPLX(0.d0,0.d0,KIND=KIND(0d0))
                        END IF
                    END DO
                END DO
            ELSE IF ( n_field_type == 81 ) THEN
            ! This reduces poloidal field to the dipole contribution:
                DO n_r=1,n_r_max
                    DO lm=1,lm_max
                        l=lm2l(lm)
                        IF ( l == 1 ) THEN
                            workA(lm,n_r)=b(lm,n_r)
                        ELSE
                            workA(lm,n_r)=CMPLX(0.d0,0.d0,KIND=KIND(0d0))
                        END IF
                    END DO
                END DO
            END IF

        !------ Outer core contribution:

        !------ Calculate needed radial derivatives:
            if ( n_field_type == 35 ) then
                call get_drNS(PstrLM,workA,lm_max_real,1,lm_max_real, &
                                            n_r_max,n_cheb_max,workB, &
                                       i_costf_init,d_costf_init,drx)
            else if ( n_field_type == 36 ) then
                call get_drNS(PadvLM,workA,lm_max_real,1,lm_max_real, &
                                            n_r_max,n_cheb_max,workB, &
                                       i_costf_init,d_costf_init,drx)
            else if ( n_field_type == 37 ) then
                call get_drNS(PdifLM,workA,lm_max_real,1,lm_max_real, &
                                            n_r_max,n_cheb_max,workB, &
                                       i_costf_init,d_costf_init,drx)

            else if ( n_field_type == 38 ) then
                call get_drNS(TstrLM,workA,lm_max_real,1,lm_max_real, &
                                            n_r_max,n_cheb_max,workB, &
                                       i_costf_init,d_costf_init,drx)
            else if ( n_field_type == 39 ) then
                call get_drNS(TomeLM,workA,lm_max_real,1,lm_max_real, &
                                            n_r_max,n_cheb_max,workB, &
                                       i_costf_init,d_costf_init,drx)
            else if ( n_field_type == 40 ) then
                call get_drNS(TadvLM,workA,lm_max_real,1,lm_max_real, &
                                            n_r_max,n_cheb_max,workB, &
                                       i_costf_init,d_costf_init,drx)
            else if ( n_field_type == 41 ) then
                call get_drNS(TdifLM,workA,lm_max_real,1,lm_max_real, &
                                            n_r_max,n_cheb_max,workB, &
                                       i_costf_init,d_costf_init,drx)
            end if

            if ( n_surface == 0 ) then
                n_r_loop_max=n_r_max
                n_field_size=(n_r_max+n_r_ic_max-2) * &
                              n_theta_max*n_phi_max
            else if ( n_surface == 1 ) then
                n_r_loop_max=1
                n_field_size=n_theta_max*n_phi_max
            end if

            DO n_rC=1,n_r_loop_max
                n_r=n_rC
                if ( n_surface == 0 ) then
                    rMov=r(n_r)
                    n_or=(n_r-1)*n_theta_max*n_phi_max
                else if ( n_surface == 1 ) then
                    n_r=n_const
                    rMov=const
                    n_or=0
                end if

                DO n=1,nThetaBs ! theta blocks
                    n_theta_start=(n-1)*sizeThetaB+1
                !-- Br:
                    if ( n_field_type == 27 ) then
                        call get_Bpol(PstrLM(1,n_r),workA(1,n_r), &
                                             dtBr,dtBt,dtBp,rMov, &
                                n_theta_start,sizeThetaB,.false.)
                    else if ( n_field_type == 28 ) then
                        call get_Bpol(PadvLM(1,n_r),workA(1,n_r), &
                                             dtBr,dtBt,dtBp,rMov, &
                                n_theta_start,sizeThetaB,.false.)
                    else if ( n_field_type == 29 ) then
                        call get_Bpol(PdifLM(1,n_r),workA(1,n_r), &
                                             dtBr,dtBt,dtBp,rMov, &
                                n_theta_start,sizeThetaB,.false.)
                    else if ( n_field_type == 81 ) then
                    !---------- get radial field diffusion and radial dipole field:
                        call get_Bpol(PdifLM(1,n_r),workA(1,n_r), &
                                             dtBr,dtBt,dtBp,rMov, &
                                n_theta_start,sizeThetaB,.false.)
                        call get_Bpol(workA(1,n_r),workA(1,n_r), &
                                            dtBt,dtBp,dtBp,rMov, &
                               n_theta_start,sizeThetaB,.false.)
                    !-- Jr:
                    else if ( n_field_type == 30 ) then
                        call get_Bpol(aj(1,n_r),workA(1,n_r), &
                                         dtBr,dtBt,dtBp,rMov, &
                            n_theta_start,sizeThetaB,.false.)
                    else if ( n_field_type == 31 ) then
                        call get_Bpol(TstrLM(1,n_r),workA(1,n_r), &
                                             dtBr,dtBt,dtBp,rMov, &
                                n_theta_start,sizeThetaB,.false.)
                    else if ( n_field_type == 32 ) then
                        call get_Bpol(workA(1,n_r),workA(1,n_r), &
                                            dtBr,dtBt,dtBp,rMov, &
                               n_theta_start,sizeThetaB,.false.)
                    else if ( n_field_type == 33 ) then
                        call get_Bpol(TadvLM(1,n_r),workA(1,n_r), &
                                             dtBr,dtBt,dtBp,rMov, &
                                n_theta_start,sizeThetaB,.false.)
                    else if ( n_field_type == 34 ) then
                        call get_Bpol(TdifLM(1,n_r),workA(1,n_r), &
                                             dtBr,dtBt,dtBp,rMov, &
                                n_theta_start,sizeThetaB,.false.)
                    !-- Bz poloidal
                    else if ( n_field_type == 13 ) then
                        call get_Bpol(b(1,n_r),db(1,n_r), &
                                     dtBr,dtBt,dtBp,rMov, &
                        n_theta_start,sizeThetaB,.false.)
                    else if ( n_field_type == 35 ) then
                        call get_Bpol(PstrLM(1,n_r),workA(1,n_r), &
                                             dtBr,dtBt,dtBp,rMov, &
                                n_theta_start,sizeThetaB,.false.)
                    else if ( n_field_type == 36 ) then
                        call get_Bpol(PadvLM(1,n_r),workA(1,n_r), &
                                             dtBr,dtBt,dtBp,rMov, &
                                n_theta_start,sizeThetaB,.false.)
                    else if ( n_field_type == 37 ) then
                        call get_Bpol(PdifLM(1,n_r),workA(1,n_r), &
                                             dtBr,dtBt,dtBp,rMov, &
                                n_theta_start,sizeThetaB,.false.)

                    !-- Jz poloidal
                    else if ( n_field_type == 14 ) then
                        call get_Bpol(aj(1,n_r),dj(1,n_r), &
                                      dtBr,dtBt,dtBp,rMov, &
                         n_theta_start,sizeThetaB,.false.)
                    else if ( n_field_type == 38 ) then
                        call get_Bpol(TstrLM(1,n_r),workA(1,n_r), &
                                             dtBr,dtBt,dtBp,rMov, &
                                n_theta_start,sizeThetaB,.false.)
                    else if ( n_field_type == 39 ) then
                        call get_Bpol(workA(1,n_r),workA(1,n_r), &
                                            dtBr,dtBt,dtBp,rMov, &
                               n_theta_start,sizeThetaB,.false.)
                    else if ( n_field_type == 40 ) then
                        call get_Bpol(TadvLM(1,n_r),workA(1,n_r), &
                                             dtBr,dtBt,dtBp,rMov, &
                                n_theta_start,sizeThetaB,.false.)
                    else if ( n_field_type == 41 ) then
                        call get_Bpol(TdifLM(1,n_r),workA(1,n_r), &
                                             dtBr,dtBt,dtBp,rMov, &
                                n_theta_start,sizeThetaB,.false.)
                    !--- Bp toriodal:
                    else if ( n_field_type == 24 ) then
                        call get_Btor(workA(1,n_r),dtBt,dtBr,rMov, &
                                 n_theta_start,sizeThetaB,.false.)
                    else if ( n_field_type == 42 ) then
                        call get_Btor(aj(1,n_r),dtBt,dtBr,rMov, &
                              n_theta_start,sizeThetaB,.false.)
                    else if ( n_field_type == 43 ) then
                        call get_Btor(TstrLM(1,n_r),dtBt,dtBr,rMov, &
                                  n_theta_start,sizeThetaB,.false.)
                    else if ( n_field_type == 44 ) then
                        call get_Btor(workA(1,n_r),dtBt,dtBr,rMov, &
                                 n_theta_start,sizeThetaB,.false.)
                    else if ( n_field_type == 45 ) then
                        call get_Btor(TadvLM(1,n_r),dtBt,dtBr,rMov, &
                                  n_theta_start,sizeThetaB,.false.)
                    else if ( n_field_type == 46 ) then
                        call get_Btor(TdifLM(1,n_r),dtBt,dtBr,rMov, &
                                  n_theta_start,sizeThetaB,.false.)
                    else if ( n_field_type == 49 ) then
                        call get_Btor(TomeLM(1,n_r),dtBt,dtBr,rMov, &
                                  n_theta_start,sizeThetaB,.false.)
                    !--- Bt toroidal
                    else if ( n_field_type == 50 ) then
                        call get_Btor(aj(1,n_r),dtBr,dtBt,rMov, &
                              n_theta_start,sizeThetaB,.false.)
                    !--- Toroidal Potential
                    else if ( n_field_type == 51 ) then
                        call lm2pt(aj(1,n_r),dtBr,rMov, &
                         n_theta_start,.false.,.false.)
                    !--- Fieldlines for theta=const.
                    else if ( n_field_type == 52 ) then
                        call get_Btor(b(1,n_r),dtBr,dtBt,rMov, &
                             n_theta_start,sizeThetaB,.false.)
                    end if

                !--- Now store the stuff of theta block to dtBrframe:
                    DO n_theta_block=1,sizeThetaB
                        n_theta_cal=n_theta_start+n_theta_block-1
                        n_theta    =n_theta_cal2ord(n_theta_cal)
                        DO n_phi=1,n_phi_max
                            n_pos=n_or+n_phi+(n_theta-1)*n_phi_max
                            IF ( n_field_type == 13 .OR. &
                                 n_field_type == 35 .OR. &
                                 n_field_type == 36 .OR. &
                                 n_field_type == 37 .OR. &
                                 n_field_type == 14 .OR. &
                                 n_field_type == 38 .OR. &
                                 n_field_type == 39 .OR. &
                                 n_field_type == 40 .OR. &
                                 n_field_type == 41 ) then
                                dtBrframe(n_pos) = &
                                  cosTheta(n_theta_cal)*dtBr(n_phi,n_theta_block) - &
                                  sinTheta(n_theta_cal)*dtBt(n_phi,n_theta_block)
                            ELSE IF ( n_field_type == 81 ) then
                                IF ( dtBr(n_phi,n_theta_block) * &
                                     dtBt(n_phi,n_theta_block) < 0.D0 ) THEN
                                    dtBrframe(n_pos)=dtBr(n_phi,n_theta_block)
                                ELSE
                                    dtBrframe(n_pos)=0.D0
                                END IF
                            ELSE
                                dtBrframe(n_pos)=dtBr(n_phi,n_theta_block)
                            END IF ! Br (Bp) or Bz ?
                        END DO ! phi Loop
                    END DO    ! theta loop

                END DO       ! theta block loop
            END DO          ! r loop

        !------ Inner core contribution:
            l_loop=.true.
            if ( n_surface == 0 ) then
                n_r_loop_max=n_r_ic_max-1
            else if ( n_surface == 1 ) then
                n_r_loop_max=2
                if ( const >= r_icb ) l_loop= .FALSE. 
            end if
            if ( l_loop ) then

            !------ Calculate needed radial derivatives:
                if ( n_field_type == 36 ) then
                    call get_drNS_even(       PadvLMIC,workA, &
                                   lm_max_real,1,lm_max_real, &
                    n_r_ic_max,n_cheb_ic_max,dr_fac_ic,workB, &
                           i_costf1_ic_init,d_costf1_ic_init, &
                           i_costf2_ic_init,d_costf2_ic_init)
                else if ( n_field_type == 37 ) then
                    call get_drNS_even(       PdifLMIC,workA, &
                                   lm_max_real,1,lm_max_real, &
                    n_r_ic_max,n_cheb_ic_max,dr_fac_ic,workB, &
                           i_costf1_ic_init,d_costf1_ic_init, &
                           i_costf2_ic_init,d_costf2_ic_init)
                else if ( n_field_type == 40 ) then
                    call get_drNS_even(       TadvLMIC,workA, &
                                   lm_max_real,1,lm_max_real, &
                    n_r_ic_max,n_cheb_ic_max,dr_fac_ic,workB, &
                           i_costf1_ic_init,d_costf1_ic_init, &
                           i_costf2_ic_init,d_costf2_ic_init)
                else if ( n_field_type == 41 ) then
                    call get_drNS_even(       TdifLMIC,workA, &
                                   lm_max_real,1,lm_max_real, &
                    n_r_ic_max,n_cheb_ic_max,dr_fac_ic,workB, &
                           i_costf1_ic_init,d_costf1_ic_init, &
                           i_costf2_ic_init,d_costf2_ic_init)
                end if

                DO n_rC=2,n_r_loop_max
                    n_r=n_rC
                    if ( n_surface == 0 ) then
                        rMov=r_ic(n_r)
                        n_or=(n_r_max+n_r-2) * &
                              n_theta_max*n_phi_max
                    else if ( n_surface == 1 ) then
                        n_r=n_const
                        rMov=const
                        n_or=0
                    end if
                    DO n=1,nThetaBs ! theta blocks
                        n_theta_start=(n-1)*sizeThetaB+1

                        if ( n_field_type == 35 .OR. &
                             n_field_type == 43 .OR. &
                             n_field_type == 44 .OR. &
                             n_field_type == 27 .OR. &
                             n_field_type == 31 .OR. &
                             n_field_type == 32 .OR. &
                             n_field_type == 38 .OR. &
                             n_field_type == 49 .OR. &
                             n_field_type == 39 ) then

                            DO n_theta=1,sizeThetaB ! no stretching !
                                DO n_phi=1,n_phi_max
                                    dtBr(n_phi,n_theta)=0.d0
                                    dtBt(n_phi,n_theta)=0.d0
                                    dtBp(n_phi,n_theta)=0.d0
                                END DO
                            END DO
                        !--- Br:
                        else if ( n_field_type == 28 ) then
                            call get_Bpol(PadvLMIC(1,n_r),workA(1,n_r), &
                                                   dtBr,dtBt,dtBp,rMov, &
                                       n_theta_start,sizeThetaB,.true.)
                        else if ( n_field_type == 29 ) then
                            call get_Bpol(PdifLMIC(1,n_r),workA(1,n_r), &
                                                   dtBr,dtBt,dtBp,rMov, &
                                       n_theta_start,sizeThetaB,.true.)
                        !--- Jr:
                        else if ( n_field_type == 30 ) then
                            call get_Bpol(aj_ic(1,n_r),workA(1,n_r), &
                                                dtBr,dtBt,dtBp,rMov, &
                                    n_theta_start,sizeThetaB,.true.)
                        else if ( n_field_type == 33 ) then
                            call get_Bpol(TadvLMIC(1,n_r),workA(1,n_r), &
                                                   dtBr,dtBt,dtBp,rMov, &
                                       n_theta_start,sizeThetaB,.true.)
                        else if ( n_field_type == 34 ) then
                            call get_Bpol(TdifLMIC(1,n_r),workA(1,n_r), &
                                                   dtBr,dtBt,dtBp,rMov, &
                                       n_theta_start,sizeThetaB,.true.)
                        !- Bz poloidal:
                        else if ( n_field_type == 13 ) then
                            call get_Bpol(b_ic(1,n_r),db_ic(1,n_r), &
                                               dtBr,dtBt,dtBp,rMov, &
                                   n_theta_start,sizeThetaB,.true.)
                        else if ( n_field_type == 36 ) then
                            call get_Bpol(PadvLMIC(1,n_r),workA(1,n_r), &
                                                   dtBr,dtBt,dtBp,rMov, &
                                       n_theta_start,sizeThetaB,.true.)
                        else if ( n_field_type == 37 ) then
                            call get_Bpol(PdifLMIC(1,n_r),workA(1,n_r), &
                                                   dtBr,dtBt,dtBp,rMov, &
                                       n_theta_start,sizeThetaB,.true.)
                        !--- Jz poloidal:
                        else if ( n_field_type == 14 ) then
                            call get_Bpol(aj_ic(1,n_r),dj_ic(1,n_r), &
                                                dtBr,dtBt,dtBp,rMov, &
                                    n_theta_start,sizeThetaB,.true.)
                        else if ( n_field_type == 40 ) then
                            call get_Bpol(TadvLMIC(1,n_r),workA(1,n_r), &
                                                   dtBr,dtBt,dtBp,rMov, &
                                       n_theta_start,sizeThetaB,.true.)
                        else if ( n_field_type == 41 ) then
                            call get_Bpol(TdifLMIC(1,n_r),workA(1,n_r), &
                                                   dtBr,dtBt,dtBp,rMov, &
                                       n_theta_start,sizeThetaB,.true.)
                        !--- Bphi toroidal:
                        else if ( n_field_type == 42 ) then
                            call get_Btor(aj_ic(1,n_r),dtBt,dtBr,rMov, &
                                      n_theta_start,sizeThetaB,.true.)
                        else if ( n_field_type == 45 ) then
                            call get_Btor(TadvLMIC(1,n_r),dtBt,dtBr,rMov, &
                                         n_theta_start,sizeThetaB,.true.)
                        else if ( n_field_type == 46 ) then
                            call get_Btor(TdifLMIC(1,n_r),dtBt,dtBr,rMov, &
                                         n_theta_start,sizeThetaB,.true.)
                        !--- Btheta toroidal:
                        else if ( n_field_type == 50 ) then
                            call get_Btor(aj_ic(1,n_r),dtBr,dtBt,rMov, &
                                      n_theta_start,sizeThetaB,.true.)
                        !--- Toroidal Potential
                        else if ( n_field_type == 51 ) then
                            call lm2pt(aj_ic(1,n_r),dtBr,r_ic(n_r), &
                                      n_theta_start,.true.,.false.)
                        !--- Fieldlines for theta=const.
                        else if ( n_field_type == 52 ) then
                            call get_Btor(b_ic(1,n_r),dtBr,dtBt,rMov, &
                                     n_theta_start,sizeThetaB,.true.)

                        end if

                    !--- Now store the stuff:
                        DO n_theta_block=1,sizeThetaB
                            n_theta_cal=n_theta_start+n_theta_block-1
                            n_theta=    n_theta_cal2ord(n_theta_cal)
                            DO n_phi=1,n_phi_max
                                n_pos=n_or+n_phi+(n_theta-1)*n_phi_max
                                if ( n_field_type == 13 .OR. &
                                     n_field_type == 36 .OR. &
                                     n_field_type == 37 .OR. &
                                     n_field_type == 14 .OR. &
                                     n_field_type == 40 .OR. &
                                     n_field_type == 41 ) then
                                    dtBrframe(n_pos) = &
                                      cosTheta(n_theta_cal)*dtBr(n_phi,n_theta_block) - &
                                      sinTheta(n_theta_cal)*dtBt(n_phi,n_theta_block)
                                else
                                    dtBrframe(n_pos)=dtBr(n_phi,n_theta_block)
                                end if ! Br or Bz ?
                            END DO ! phi Loop
                        END DO    ! theta loop

                    END DO       ! theta block loop
                END DO          ! r loop

            END IF ! l_loop ?

        !--- Write frame field:
            WRITE(n_out) (real(dtBrframe(n),4),n=1,n_field_size)

        END DO  ! LOOP over fields in movie file

    END IF


    RETURN
    end SUBROUTINE write_dtB_frame

!--- End of subroutine  s_write_dtB_frame.f
!----------------------------------------------------------------------------
