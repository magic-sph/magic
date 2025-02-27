module init_fields
   !
   ! This module is used to construct the initial solution.
   !

   use precision_mod
   use iso_fortran_env, only: output_unit
   use parallel_mod
   use chebyshev, only: type_cheb_odd
   use finite_differences, only: type_fd
   use communications, only: lo2r_one, r2lo_one, transp_R2Phi, transp_Phi2R
   use truncation, only: n_r_max, n_r_maxMag, n_r_ic_max, m_min, l_max, &
       &                 n_phi_max, n_theta_max, n_r_tot, m_max,        &
       &                 minc, n_cheb_ic_max, lm_max, nlat_padded
   use mem_alloc, only: bytes_allocated
   use blocking, only: lo_map, st_map, llm, ulm, llmMag, ulmMag
   use horizontal_data, only: sinTheta, dLh, dTheta1S, dTheta1A, &
       &                      phi, cosTheta, hdif_B
   use logic, only: l_rot_ic, l_rot_ma, l_SRIC, l_SRMA, l_cond_ic,  &
       &            l_temperature_diff, l_chemical_conv, l_onset,   &
       &            l_anelastic_liquid, l_non_adia, l_finite_diff,  &
       &            l_non_rot, l_force_v, l_tidal
   use radial_functions, only: r_icb, r, r_cmb, r_ic, or1, jVarCon,    &
       &                       lambda, or2, dLlambda, or3, cheb_ic,    &
       &                       dcheb_ic, d2cheb_ic, cheb_norm_ic, or1, &
       &                       r_ic, orho1, chebt_ic, temp0,           &
       &                       dLtemp0, kappa, dLkappa, beta, dbeta,   &
       &                       epscProf, ddLtemp0, ddLalpha0, rgrav,   &
       &                       rho0, dLalpha0, alpha0, otemp1, ogrun,  &
       &                       rscheme_oc, or1, O_r_ic, l_R
   use radial_data, only: n_r_icb, n_r_cmb, nRstart, nRstop
   use radial_der, only: get_dr_Rloc, get_dr
   use constants, only: pi, y10_norm, c_z10_omega_ic, c_z10_omega_ma, osq4pi, &
       &                zero, one, two, three, four, third, half, sq4pi,c_moi_oc,c_moi_ic
   use useful, only: abortRun, lagrange_interp, cc2real
   use sht, only: scal_to_SH, scal_to_spat, torpol_to_spat
   use physical_parameters, only: impS, n_impS_max, n_impS, phiS, thetaS, &
       &                          peakS, widthS, radratio, imagcon, opm,  &
       &                          sigma_ratio, O_sr, kbots, ktops, opr,   &
       &                          epsc, ViscHeatFac, ThExpNb, tmelt,      &
       &                          impXi, n_impXi_max, n_impXi, phiXi,     &
       &                          thetaXi, peakXi, widthXi, osc, epscxi,  &
       &                          kbotxi, ktopxi, BuoFac, ktopp, oek,     &
       &                          epsPhase,ktopv, kbotv, LFfac,interior_model, &
       &                          ktopb, r_cut_model, tidalFac, w_orbit_th
   use algebra, only: prepare_mat, solve_mat
   use cosine_transform_odd
   use dense_matrices
   use real_matrices
   use band_matrices
   use outRot, only: get_angular_moment
   use integration, only: rInt_R
   use num_param, only: lScale, eScale
   use special, only: n_imp,amp_imp, ampForce
   use fields, only: z0v_Rloc, dz0v_Rloc, z0v_LMloc

   implicit none

   private

   !-- Initialisation of fields:
   integer, public :: init_s1,init_s2
   integer, public :: init_xi1,init_xi2
   integer, public :: init_b1,init_v1
   integer, public :: init_phi ! An integer to specify phase field initial configuration

   !----- Entropy amplitudes for initialisation:
   real(cp), public :: amp_s1,amp_s2,amp_v1,amp_b1,amp_xi1,amp_xi2
   real(cp), public :: init_b_length_min, init_b_length_max,init_b_index

   !----- Entropy at CMB and ICB (input):
   integer, public, parameter :: n_s_bounds=20
   real(cp), public :: s_bot(4*n_s_bounds)  ! input variables for tops,bots
   real(cp), public :: s_top(4*n_s_bounds)
   complex(cp), public, allocatable :: tops(:,:)
   complex(cp), public, allocatable :: bots(:,:)

   !----- Velocity at CMB (input):
   complex(cp),public, allocatable :: topv(:,:)
   !----- Velocity at ICB (input):
   complex(cp),public, allocatable :: botv(:,:)

   !----- Chemical composition
   integer, public, parameter :: n_xi_bounds=20
   real(cp), public :: xi_bot(4*n_xi_bounds)  ! input variables for topxi,botxi
   real(cp), public :: xi_top(4*n_xi_bounds)
   complex(cp), public, allocatable :: topxi(:,:)
   complex(cp), public, allocatable :: botxi(:,:)

   !---- Phase field
   real(cp), public :: phi_top ! Phase field value at the outer boundary
   real(cp), public :: phi_bot ! Phase field value at the inner boundary

   !----- Peak values for magnetic field:
   real(cp), public :: bpeakbot,bpeaktop

   !----- Initialised IC and mantle rotation rates:
   integer, public :: nRotMa,nRotIc
   real(cp), public :: omega_ma1,omegaOsz_ma1,tShift_ma1,tOmega_ma1
   real(cp), public :: omega_ma2,omegaOsz_ma2,tShift_ma2,tOmega_ma2
   real(cp), public :: omega_ic1,omegaOsz_ic1,tShift_ic1,tOmega_ic1
   real(cp), public :: omega_ic2,omegaOsz_ic2,tShift_ic2,tOmega_ic2

   !----- Rotation profile
   real(cp), public :: q_rot, norm_ome
   integer, public :: pertur_z
   integer, public :: pertur_w
   real(cp), public :: rand_num
   !----- Forcing parameters
   integer, public :: force_z_vol !Volume forcing parameter, 1= viscous compensation + drag force, 2= all (m=0) modes = 0
   real(cp), public :: tau !parameter for the drag force


   !----- About start-file:
   logical, public :: l_start_file     ! taking fields from startfile ?
   logical, public :: l_reset_t        ! reset time from startfile ?
   integer, public :: inform           ! format of start_file
   character(len=72), public :: start_file  ! name of start_file

   !-- Scales for input field:
   real(cp), public :: scale_s
   real(cp), public :: scale_xi
   real(cp), public :: scale_v
   real(cp), public :: scale_b
   real(cp), public :: tipdipole       ! adding to symetric field

   !Tidal velocity on the grid
   real(cp), public, allocatable :: vrtidal(:,:,:)
   real(cp), public, allocatable :: vttidal(:,:,:)
   real(cp), public, allocatable :: vptidal(:,:,:)
   real(cp), public, allocatable :: X(:), dX(:), d2X(:)

   public :: initialize_init_fields, initV, initS, initB, initPhi, initF, ps_cond, &
   &         pt_cond, initXi, xi_cond, finalize_init_fields, initTidal

contains

   subroutine initialize_init_fields
      !
      ! Memory allocation
      !
      tOmega_ic1=0.0_cp
      tOmega_ic2=0.0_cp

      allocate( tops(0:l_max,0:m_max), bots(0:l_max,0:m_max) )
      tops(:,:)=zero
      bots(:,:)=zero
      bots(0,0)=one
      tops(0,0)=0.0_cp
      bytes_allocated = bytes_allocated+2*(l_max+1)*(m_max+1)*SIZEOF_DEF_COMPLEX

      if ( ktopv == 3) then
         allocate( topv(0:l_max,0:m_max))
         topv(:,:)=zero
         bytes_allocated = bytes_allocated + (l_max+1)*(m_max+1)*SIZEOF_DEF_COMPLEX
      end if
      if ( kbotv == 3) then
         allocate( botv(0:l_max,0:m_max))
         botv(:,:)=zero
         bytes_allocated = bytes_allocated + (l_max+1)*(m_max+1)*SIZEOF_DEF_COMPLEX
      end if
      if ( l_chemical_conv ) then
         allocate( topxi(0:l_max,0:m_max), botxi(0:l_max,0:m_max) )
         topxi(:,:)=zero
         botxi(:,:)=zero
         botxi(0,0)=one
         topxi(0,0)=0.0_cp
         bytes_allocated = bytes_allocated+2*(l_max+1)*(m_max+1)*SIZEOF_DEF_COMPLEX
      end if

      if (l_tidal) then
         allocate(vttidal(nlat_padded,n_phi_max,nRstart:nRstop) )
         allocate(vptidal(nlat_padded,n_phi_max,nRstart:nRstop) )
         allocate(vrtidal(nlat_padded,n_phi_max,nRstart:nRstop) )
         vttidal(:,:,:)=zero
         vrtidal(:,:,:)=zero
         vptidal(:,:,:)=zero
         bytes_allocated=bytes_allocated + 3*n_phi_max*nlat_padded*(nRstop-nRstart+1)*&
              &               SIZEOF_DEF_REAL
         allocate(X(n_r_max))
         allocate(dX(n_r_max))
         allocate(d2X(n_r_max))
         X(:)=zero
         dX(:)=zero
         d2X(:)=zero
         bytes_allocated=bytes_allocated + 3*n_r_max*&
              &               SIZEOF_DEF_REAL

      end if

   end subroutine initialize_init_fields
!------------------------------------------------------------------------------
   subroutine finalize_init_fields
      !
      ! Memory deallocation
      !
      deallocate (tops, bots )
      if ( l_chemical_conv ) deallocate( topxi, botxi )
      if ( ktopv == 3 )   deallocate( topv)
      if ( kbotv == 3 )   deallocate( botv)
      if (l_tidal) deallocate (vrtidal,vttidal,vptidal,X,dX,d2X)

   end subroutine finalize_init_fields
!------------------------------------------------------------------------------
   subroutine initV(w,z,omega_ic,omega_ma)
      !
      ! Purpose of this subroutine is to initialize the velocity field
      ! So far it is only rudimentary and will be expanded later.
      ! Because s is needed for dwdt init_s has to be called before.
      !

      !-- Output variables
      complex(cp), intent(inout) :: w(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: z(llm:ulm,n_r_max)
      real(cp), intent(out) :: omega_ic,omega_ma

      !-- Local variables
      complex(cp) :: z_Rloc(lm_max,nRstart:nRstop)
      integer :: lm,l,m,l1m0
      integer :: nR,nTheta,nPhi
      real(cp) :: ra1,ra2,c_r,c_i
      real(cp) :: amp_r,rExp
      real(cp) :: rDep(n_r_max)

      !For angular momentum computation !ARS
      real(cp):: trash1(3),trash2(3), angular_moment(3)
      real(cp) :: lr_z,r_cut
      integer :: l1m1,rank_with_cmb
      integer :: sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1)
      integer :: irank,i

!      class(type_mpitransp), pointer :: r2lo_initv, lo2r_initv, r2lo_initv2, lo2r_initv2

      real(cp) :: ss,ome(nlat_padded,n_phi_max)
      complex(cp) :: omeLM(lm_max)
      complex(cp) :: dz1(llm:ulm,n_r_max)
      real(cp), allocatable :: coeffOmega(:)

! !<<<<<<< HEAD
!       allocate( type_mpiptop :: r2lo_initv )
!       allocate( type_mpiptop :: lo2r_initv )
!       allocate( type_mpiptop :: r2lo_initv2 )
!       allocate( type_mpiptop :: lo2r_initv2 )

      if (.not. l_force_v) then
         z0v_Rloc(:,:) =0.0_cp
         dz0v_Rloc(:,:) =0.0_cp
         z0v_LMloc(:,:) =0.0_cp
      end if

      !-- Initialize rotation according to
      !   given inner core and mantel rotation rate:
      if ( init_v1 == 1 .and. ( omega_ic1 /= 0.0_cp .or. omega_ma1 /= 0.0_cp ) ) then

         !-- From lo distributed to r distributed
         call lo2r_one%transp_lm2r(z, z_Rloc)

         !-- Approximating the Stewardson solution:
         do nR=nRstart,nRstop

            do nPhi=1,n_phi_max
               do nTheta=1,n_theta_max
                  ss=r(nR)*sinTheta(nTheta)
                  if ( ss <= r_icb ) then
                     ome(nTheta,nPhi)=omega_ma1+half*omega_ic1
                  else
                     ome(nTheta,nPhi)=omega_ma1
                  end if
               end do
            end do
            call scal_to_SH(ome, omeLM, l_max)

            !------- ome now in spherical harmonic space,
            !        apply operator dTheta1=1/(r sinTheta) d/ d theta sinTheta**2,
            !        additional application of r**2/(l*(l+1)) then yields
            !        the axisymmetric toriodal flow contribution:
            do lm=2,lm_max
               l   =st_map%lm2l(lm)
               m   =st_map%lm2m(lm)
               if ( l < l_max .and. l > m ) then
                  z_Rloc(lm,nR)=z_Rloc(lm,nR) + r(nR)**2/dLh(lm) * ( &
                  &          dTheta1S(lm)*omeLM(st_map%lm2lmS(lm))   &
                  &         -dTheta1A(lm)*omeLM(st_map%lm2lmA(lm)) )
               else if ( l < l_max .and. l == m ) then
                  z_Rloc(lm,nR)=z_Rloc(lm,nR) - r(nR)**2/dLh(lm) *  &
                  &          dTheta1A(lm)*omeLM(st_map%lm2lmA(lm))
               else if ( l == l_max .and. m < l ) then
                  z_Rloc(lm,nR)=z_Rloc(lm,nR) + r(nR)**2/dLh(lm) *  &
                  &          dTheta1S(lm)*omeLM(st_map%lm2lmS(lm))
               end if
            end do

         end do ! close loop over radial grid points

         !-- Transpose back to lo distributed
         call r2lo_one%transp_r2lm(z_Rloc, z)

      else if ( init_v1 == 2 ) then

         !-- From lo distributed to r distributed
         call lo2r_one%transp_lm2r(z, z_Rloc)

         !-- Approximating the Stewardson solution:
         do nR=nRstart,nRstop

            do nPhi=1,n_phi_max
               do nTheta=1,n_theta_max
                  ss=r(nR)*sinTheta(nTheta)
                  ome(nTheta,nPhi)=amp_v1/sqrt(one+ss**4)
               end do
            end do
            call scal_to_SH(ome, omeLM, l_max)

            !------------ ome now in spherical harmonic space,
            !             apply operator dTheta1=1/(r sinTheta) d/ d theta sinTheta**2,
            !             additional application of r**2/(l*(l+1)) then yields
            !             the axisymmetric toriodal flow contribution:
            do lm=2,lm_max
               l   =st_map%lm2l(lm)
               m   =st_map%lm2m(lm)
               if ( l < l_max .and. l > m ) then
                  z_Rloc(lm,nR)=z_Rloc(lm,nR) + r(nR)**2/dLh(lm) * ( &
                  &            dTheta1S(lm)*omeLM(st_map%lm2lmS(lm)) &
                  &          - dTheta1A(lm)*omeLM(st_map%lm2lmA(lm)) )
               else if ( l < l_max .and. l == m ) then
                   z_Rloc(lm,nR)=z_Rloc(lm,nR) - r(nR)**2/dLh(lm) * &
                   &           dTheta1A(lm)*omeLM(st_map%lm2lmA(lm))
               else if ( l == l_max  .and. m < l) then
                  z_Rloc(lm,nR)=z_Rloc(lm,nR) + r(nR)**2/dLh(lm) *   &
                  &            dTheta1S(lm)*omeLM(st_map%lm2lmS(lm))
               end if
            end do

         end do ! close loop over radial grid points
         !-- Transpose back to lo distributed
         call r2lo_one%transp_r2lm(z_Rloc, z)

      else if (init_v1 == 30) then
         !--- Add a rotation profile proportionnal to r^q where q is a parameter
         ! Amplitude is made to be coherent with the solid body rotation of 1/Ekman

         call lo2r_one%transp_lm2r(z, z_Rloc)

         if ( index(interior_model,'PNS_SZ_0V2S') /= 0 .OR. pertur_z == 1 ) then
            r_cut=0.25_cp
         else
            r_cut=0.6_cp
         end if
         !-- Approximating the Stewardson solution:
         !if (rank ==0) then
         do nR=nRstart,nRstop
            do nPhi=1,n_phi_max
               do nTheta=1,n_theta_max
                  !------------ start with constructing rotation rate ome:
                  if ( index(interior_model,'PNS_SZ_0V2S') /= 0 ) then
                     ss=r(nR)*sinTheta(nTheta)/(r_cut*r_cmb)
                     ome(nTheta,nPhi)=rho0(nR)*one*(one+ss**(20.0_cp*q_rot))**0.05_cp
                  else
                     ss=r(nR)*sinTheta(nTheta)/(r_cut*r_cmb)
                     ome(nTheta,nPhi)=one*(one+ss**(20.0_cp*q_rot))**0.05_cp
                  end if
               end do
            end do

            !------------ Transform to spherical hamonic space for each theta block
            call scal_to_SH(ome, omeLM, l_max)
            !------- ome now in spherical harmonic space,
            !        apply operator dTheta1=1/(r sinTheta) d/ d theta sinTheta**2,
            !        additional application of r**2/(l*(l+1)) then yields
            !        the axisymmetric toriodal flow contribution:
            do lm=2,lm_max
               l   =st_map%lm2l(lm)
               m   =st_map%lm2m(lm)
               if ( l < l_max .and. l > m ) then
                  z_Rloc(lm,nR)=z_Rloc(lm,nR) + r(nR)**2/dLh(lm) * ( &
                  &          dTheta1S(lm)*omeLM(st_map%lm2lmS(lm))   &
                  &         -dTheta1A(lm)*omeLM(st_map%lm2lmA(lm)) )
               else if ( l < l_max .and. l == m ) then
                  z_Rloc(lm,nR)=z_Rloc(lm,nR) - r(nR)**2/dLh(lm) *  &
                  &          dTheta1A(lm)*omeLM(st_map%lm2lmA(lm))
               else if ( l == l_max .and. m < l ) then
                  z_Rloc(lm,nR)=z_Rloc(lm,nR) + r(nR)**2/dLh(lm) *  &
                  &          dTheta1S(lm)*omeLM(st_map%lm2lmS(lm))
               end if
               if (nR == n_r_icb .AND. kbotv == 3) then
                  botv(l,m) = botv(l,m) + z_Rloc(lm,nR)
               end if
               if (nR == 1 .AND. ktopv == 3) then
                  topv(l,m) = topv(l,m) + z_Rloc(lm,1)
               end if
            end do
         end do ! close loop over radial grid points

         if (ktopv == 3) then !Need to broadcast topv since only the rank_with_cmb has it
#ifdef WITH_MPI
            call MPI_Bcast(topv,(l_max+1)*(m_max+1),MPI_DEF_COMPLEX,0,MPI_COMM_WORLD,ierr)
#endif
         end if
         if (kbotv == 3 ) then !Need to broadcast topv since only the rank_with_icb has it
#ifdef WITH_MPI
            call MPI_Bcast(botv,(l_max+1)*(m_max+1),MPI_DEF_COMPLEX,n_procs-1,MPI_COMM_WORLD,ierr)
#endif
         end if

         !-- Transpose back to lo distributed
         call r2lo_one%transp_r2lm(z_Rloc, z)

         l1m0 = lo_map%lm2(1,0)
         l1m1 = lo_map%lm2(1,1)
        if (l1m0 >= llm .and. l1m0 <= ulm) then
           call get_angular_moment(z(l1m0,:),z(l1m1,:),0.D0,0.D0,angular_moment,trash1,trash2)
           lr_z = angular_moment(3)
           do nR=1,n_r_max
              z(l1m0,nR)=z(l1m0,nR) - r(nR)*r(nR)*lr_z*rho0(nR)/(c_moi_oc*y10_norm) !TO REMOVE THE UNIFORM ROTATION
           end do
           norm_ome = c_moi_oc*oek/lr_z
           if (ktopv == 3) then
              write(output_unit,*) "c_z10_omega_ma",c_z10_omega_ma
              write(output_unit,*) "rscheme_oc%rnorm",rscheme_oc%rnorm
              topv(1,0)= (topv(1,0) - r(1)*r(1)*lr_z*rho0(1)/(c_moi_oc*y10_norm))*c_z10_omega_ma!*rscheme_oc%rnorm
              topv = topv*norm_ome!/rscheme_oc%rnorm
              if( index(interior_model,'PNS_SZ_0V2S') == 0) then
                 omega_ic1 = c_z10_omega_ic*real(z(l1m0,n_r_icb))*norm_ome
              end if
           end if
           if (kbotv == 3) then
              botv(1,0)= (botv(1,0) - r(n_r_icb)*r(n_r_icb)*lr_z*rho0(n_r_icb)/(c_moi_oc*y10_norm))*c_z10_omega_ic!*rscheme_oc%rnorm
              botv = botv*norm_ome!/rscheme_oc%rnorm
           end if
           write (*,*) "rho0(1)",rho0(1)
           write (*,*) "r_cut=",r_cut
           write (*,*) "ANGULAR MOMENT =", lr_z
           write (*,*) "c_moi_oc=",c_moi_oc
           write (*,*) "oek =", oek
           write (*,*) "Norm omega =",norm_ome
        end if

        !Broadcast norm_ome, omega_ic1, topv
#ifdef WITH_MPI
         call MPI_Bcast(norm_ome,1,MPI_DEF_REAL,rank_with_l1m0,MPI_COMM_WORLD,ierr)
         call MPI_Bcast(omega_ic1,1,MPI_DEF_REAL,rank_with_l1m0,MPI_COMM_WORLD,ierr)
#endif
         if (ktopv == 3 ) then !Need to broadcast topv since only the rank_with_l1m0 has it
#ifdef WITH_MPI
            call MPI_Bcast(topv,(l_max+1)*(m_max+1),MPI_DEF_COMPLEX,rank_with_l1m0,MPI_COMM_WORLD,ierr)
#endif
         end if
         if (kbotv == 3 ) then !Need to broadcast botv since only the rank_with_l1m0 has it
#ifdef WITH_MPI
            call MPI_Bcast(botv,(l_max+1)*(m_max+1),MPI_DEF_COMPLEX,rank_with_l1m0,MPI_COMM_WORLD,ierr)
#endif
         end if

         z= z*norm_ome

         if (force_z_vol > 0 .and. l_force_v) then
            z0v_LMloc(:,:) = z(:,:)
            call get_dr(z,dz1,ulm-llm+1,1,ulm-llm+1,n_r_max,rscheme_oc, &
                 &      nocopy=.true.)
            call lo2r_one%transp_lm2r(z, z0v_Rloc)
            call lo2r_one%transp_lm2r(dz1, dz0v_Rloc)
            !            call get_dr_Rloc(z0v_Rloc, dz0v_Rloc, lm_max, nRstart, nRstop, n_r_max, &
            !     &           rscheme_oc)
            z=z*0.0_cp
         end if


      else if (init_v1 == -1) then
         !--- Add a rotation profile proportionnal to r^-q where q is a parameter
         ! Amplitude is made to be coherent with the solid body rotation of 1/Ekman

         !-- From lo distributed to r distributed
         call lo2r_one%transp_lm2r(z, z_Rloc)

         if ( index(interior_model,'PNS_SZ_0V2S') /= 0 .OR. pertur_z == 1 ) then
            r_cut=0.25
         else if (index(interior_model,'HMNS_17MS_SZ') /= 0 .or. pertur_w == 1) then
            allocate(coeffOmega(15))
            coeffOmega = [ 1.92802244e+08_cp, -1.34146749e+09_cp,  3.92111322e+09_cp, -5.96788010e+09_cp, &
                 4.19854130e+09_cp,  1.00761404e+09_cp, -5.10292509e+09_cp,  5.22604331e+09_cp, &
                 -2.98578864e+09_cp,  1.05035061e+09_cp, -2.22591449e+08_cp,  2.53181706e+07_cp, &
                 -1.15262826e+06_cp,  2.03799231e+04_cp,  3.24442843e+03_cp] &
                 + [-5.80099523e-02_cp, -3.16876435e+00_cp,  1.04622507e+00_cp, -4.43496513e+00_cp, &
                 1.91004753e-01_cp, -3.42118359e+00_cp, -6.72642708e-01_cp, -1.03414631e+00_cp, &
                 -3.97218466e+00_cp, -7.47951031e-01_cp,  6.97055459e-03_cp, -9.69161466e-03_cp, &
                 -3.41950869e-03_cp, -2.82726578e-05_cp, -2.16976377e-06_cp]
            norm_ome = oek/q_rot ! to put the profile in visocus units/ in q_rot is encoded the rotation rate
            !            norm_ome =oek*1.71e-4_cp !to put the profile in viscous units
            r_cut=1.0_cp
         else
            r_cut=0.4
         end if
         !-- Approximating the Stewardson solution:
         !if (rank ==0) then
         do nR=nRstart,nRstop
            do nPhi=1,n_phi_max
               do nTheta=1,n_theta_max
                  ss=r(nR)*sinTheta(nTheta)*r_cut_model/(r_cut*r_cmb)
                  !------------ start with constructing rotation rate ome:
                  if ( index(interior_model,'PNS_SZ_0V2S') /= 0 .or.  pertur_z == 1 ) then
                     ome(nTheta,nPhi)=rho0(nR)*one/(one+ss**(20.0_cp*q_rot))**0.05_cp
                  else if (index(interior_model,'HMNS_17MS_SZ') /= 0 .or. pertur_w == 1) then
                     ome(nTheta,nPhi) = -oek !To remove the frame rotation
                     do i = 1,15
                        ome(nTheta,nPhi) = ome(nTheta,nPhi)+norm_ome*coeffOmega(i)*(ss)**(15-i)
                     end do
                     ome(nTheta,nPhi) = rho0(nR)*ome(nTheta,nPhi)
                  else
                     ome(nTheta,nPhi)=one/(one+ss**(20.0_cp*q_rot))**0.05_cp
                  end if
               end do
            end do
               !------------ Transform to spherical hamonic space for each theta block
            call scal_to_SH(ome, omeLM, l_max)

            !------- ome now in spherical harmonic space,
            !        apply operator dTheta1=1/(r sinTheta) d/ d theta sinTheta**2,
            !        additional application of r**2/(l*(l+1)) then yields
            !        the axisymmetric toriodal flow contribution:
            do lm=2,lm_max
               l   =st_map%lm2l(lm)
               m   =st_map%lm2m(lm)
               if ( l < l_max .and. l > m ) then
                  z_Rloc(lm,nR)=z_Rloc(lm,nR) + r(nR)**2/dLh(lm) * ( &
                  &          dTheta1S(lm)*omeLM(st_map%lm2lmS(lm))   &
                  &         -dTheta1A(lm)*omeLM(st_map%lm2lmA(lm)) )
               else if ( l < l_max .and. l == m ) then
                  z_Rloc(lm,nR)=z_Rloc(lm,nR) - r(nR)**2/dLh(lm) *  &
                  &          dTheta1A(lm)*omeLM(st_map%lm2lmA(lm))
               else if ( l == l_max .and. m < l ) then
                  z_Rloc(lm,nR)=z_Rloc(lm,nR) + r(nR)**2/dLh(lm) *  &
                  &          dTheta1S(lm)*omeLM(st_map%lm2lmS(lm))
               end if
               if (nR == n_r_icb .AND. kbotv == 3) then
                  botv(l,m) = botv(l,m) + z_Rloc(lm,nR)
               end if
               if (nR == 1 .AND. ktopv == 3) then
                  topv(l,m) = topv(l,m) + z_Rloc(lm,1)
               end if
            end do
         end do ! close loop over radial grid points

         if (ktopv == 3) then !Need to broadcast topv since only the rank_with_cmb has it
#ifdef WITH_MPI
            call MPI_Bcast(topv,(l_max+1)*(m_max+1),MPI_DEF_COMPLEX,0,MPI_COMM_WORLD,ierr)
#endif
         end if
         if (kbotv == 3 ) then !Need to broadcast topv since only the rank_with_icb has it
#ifdef WITH_MPI
            call MPI_Bcast(botv,(l_max+1)*(m_max+1),MPI_DEF_COMPLEX,n_procs-1,MPI_COMM_WORLD,ierr)
#endif
         end if

         !-- Transpose back to lo distributed
         call r2lo_one%transp_r2lm(z_Rloc, z)

         l1m0 = lo_map%lm2(1,0)
         l1m1 = lo_map%lm2(1,1)
        if (l1m0 >= llm .and. l1m0 <= ulm ) then
           if (index(interior_model,'HMNS_17MS_SZ') == 0 .and. pertur_w /= 1 ) then
              call get_angular_moment(z(l1m0,:),z(l1m1,:),0.D0,0.D0,angular_moment,trash1,trash2)
              lr_z = angular_moment(3)
              do nR=1,n_r_max
                 z(l1m0,nR)=z(l1m0,nR) - r(nR)*r(nR)*lr_z*rho0(nR)/(c_moi_oc*y10_norm) !TO REMOVE THE UNIFORM ROTATION
              end do
              norm_ome = c_moi_oc*oek/lr_z
              if (ktopv == 3) then
                 topv(1,0)= (topv(1,0) - r(1)*r(1)*lr_z*rho0(1)/(c_moi_oc*y10_norm))*c_z10_omega_ma!*rscheme_oc%rnorm
                 topv = topv*norm_ome!/rscheme_oc%rnorm
                 if( index(interior_model,'PNS_SZ_0V2S') == 0) then
                    omega_ic1 = c_z10_omega_ic*real(z(l1m0,n_r_icb))*norm_ome
                 end if
              end if
              if (kbotv == 3) then
                 botv(1,0)= (botv(1,0) - r(n_r_icb)*r(n_r_icb)*lr_z*rho0(n_r_icb)/(c_moi_oc*y10_norm))&
                 & *c_z10_omega_ic!*rscheme_oc%rnorm
                 botv = botv*norm_ome!/rscheme_oc%rnorm
              end if
           else
              if (ktopv ==3) then
                 topv(1,0) = topv(1,0)*c_z10_omega_ma
              end if
              if (kbotv ==2) then
                 omega_ic1 = c_z10_omega_ic*real(z(l1m0,n_r_icb))
              else if (kbotv ==3) then
                 botv(1,0) = botv(1,0)*c_z10_omega_ic
              end if
           end if

           write (output_unit,*) "rho0(1)",rho0(1),"r_cut=",r_cut,"ANGULAR MOMENT =", lr_z,&
                &"c_moi_oc=",c_moi_oc,"oek =", oek,"Norm omega =",norm_ome

        end if


        !Broadcast norm_ome, omega_ic1, topv
#ifdef WITH_MPI
         call MPI_Bcast(norm_ome,1,MPI_DEF_REAL,rank_with_l1m0,MPI_COMM_WORLD,ierr)
         call MPI_Bcast(omega_ic1,1,MPI_DEF_REAL,rank_with_l1m0,MPI_COMM_WORLD,ierr)
#endif
         if (ktopv == 3 ) then !Need to broadcast topv since only the rank_with_l1m0 has it
#ifdef WITH_MPI
            call MPI_Bcast(topv,(l_max+1)*(m_max+1),MPI_DEF_COMPLEX,rank_with_l1m0,MPI_COMM_WORLD,ierr)
#endif
         end if
         if (kbotv == 3 ) then !Need to broadcast botv since only the rank_with_l1m0 has it
#ifdef WITH_MPI
            call MPI_Bcast(botv,(l_max+1)*(m_max+1),MPI_DEF_COMPLEX,rank_with_l1m0,MPI_COMM_WORLD,ierr)
#endif
         end if

         if (index(interior_model,'HMNS_17MS_SZ') /= 0 .or. pertur_w ==1 ) then
            norm_ome =1.0_cp
            deallocate(coeffOmega)
         end if

         z= z*norm_ome

         if (force_z_vol > 0 .and. l_force_v) then
            z0v_LMloc(:,:) = z(:,:)
            call get_dr(z,dz1,ulm-llm+1,1,ulm-llm+1,n_r_max,rscheme_oc, &
                 &      nocopy=.true.)
            call lo2r_one%transp_lm2r(z, z0v_Rloc)
            call lo2r_one%transp_lm2r(dz1, dz0v_Rloc)
            !            call get_dr_Rloc(z0v_Rloc, dz0v_Rloc, lm_max, nRstart, nRstop, n_r_max, &
            !     &           rscheme_oc)
            omega_ic1=0.0_cp
            omega_ma1=0.0_cp
            omega_ic2=0.0_cp
            omega_ma2=0.0_cp
            z=z*0.0_cp
         end if

      else if ( init_v1 > 2 ) then

         !--- Add random noise toroidal field of all (l,m) modes exept (l=0,m=0):
         !    It decays likes l**(init_v1-1)
         !    Amplitude is chosen so that the (1,0) term resembles amp_v1 *
         !    the 'solid body' rotation set by inner core and mantle rotation.

         rExp=4.
         if ( omega_ic1 /= 0 ) then
            amp_r=amp_v1*omega_ic1*r_ICB**(rExp+1.)/y10_norm
         else
            amp_r=amp_v1*r_ICB**(rExp+1.)/y10_norm
         end if

         do nR=1,n_r_max
            rDep(nR)=amp_r*or1(nR)**(rExp-1.)
            !write(output_unit,"(A,I3,A,ES20.12)") "rDep(",nR,") = ",rDep(nR)
            do lm=llm,ulm
               l=lo_map%lm2l(lm)
               m=lo_map%lm2m(lm)
               if ( l /= 0 ) then
                  call random_number(ra1)
                  call random_number(ra2)
                  ra1=(-one+two*ra1)/(real(l,cp))**(init_v1-1)
                  ra2=(-one+two*ra2)/(real(l,cp))**(init_v1-1)
                  c_r=ra1*rDep(nR)
                  c_i=ra2*rDep(nR)
                  if ( m == 0 ) then  ! only axisymmetric modes
                     z(lm,nR)=z(lm,nR)+cmplx(c_r,0.0_cp,kind=cp)
                  else
                     z(lm,nR)=z(lm,nR)+cmplx(c_r,c_i,kind=cp)
                  end if
               end if
               write(output_unit,"(A,4I4,2ES20.12)") "z = ",nR,lm,l,m,z(lm,nR)
            end do
         end do

      else if ( init_v1 < -1 ) then

         !--- Add random noise poloidal field of all (l,m) modes exept (l=0,m=0):
         !    It decays likes l**(-init_v1+1)
         !    Amplitude is chosen to be comparable to amp * inner core roation speed
         !    at inner core boundary...
         if ( omega_ic1 /= 0.0_cp ) then
            amp_r=amp_v1*omega_ic1*r_icb*r_icb/(y10_norm*PI)
         else
            amp_r=amp_v1
         end if

         do nR=1,n_r_max
            rDep(nR)=-amp_r*sin( (r(nR)-r_ICB)*PI )
            do lm=llm,ulm
               l=lo_map%lm2l(lm)
               m=lo_map%lm2m(lm)
               if ( l /= 0 ) then
                  call random_number(ra1)
                  call random_number(ra2)
                  ra1=(-one+two*ra1)/(real(l,cp))**(-init_v1-1)
                  ra2=(-one+two*ra2)/(real(l,cp))**(-init_v1-1)
                  c_r=ra1*rDep(nR)
                  c_i=ra2*rDep(nR)
                  if ( m > 0 ) then  ! no axisymmetric modes
                     w(lm,nR)=w(lm,nR)+cmplx(c_r,c_i,kind=cp)
                     z(lm,nR)=z(lm,nR)+cmplx(c_r,c_i,kind=cp)
                  end if
               end if
            end do
         end do !end of the radial loop

      end if

      if (pertur_z>2) then
         !--- Add random noise toroidal field of all (l,m) modes exept (l=0,m=0):
         !    It decays likes l**(init_v1-1)
         !    Amplitude is chosen so that the (1,0) term resembles amp_v1 *
         !    the 'solid body' rotation set by inner core and mantle rotation.

         rExp=4.
         if ( omega_ic1 /= 0 ) then
            amp_r=amp_v1*omega_ic1*r_ICB**(rExp+1.)/y10_norm
         else
            amp_r=amp_v1*r_ICB**(rExp+1.)/y10_norm
         end if
         do nR=1,n_r_max
            rDep(nR)=amp_r/r(nR)**(rExp-1.)
            !write(output_unit,"(A,I3,A,ES20.12,ES20.12)") "rDep(",nR,") = ",rDep(nR)
            do lm=llm,ulm
               l=lo_map%lm2l(lm)
               m=lo_map%lm2m(lm)
               if ( l /= 0 ) then
                  call random_number(ra1)
                  call random_number(ra2)
                  ra1=(-one+two*ra1)/(real(l,cp))**(pertur_z-1)
                  ra2=(-one+two*ra2)/(real(l,cp))**(pertur_z-1)
                  c_r=ra1*rDep(nR)
                  c_i=ra2*rDep(nR)
                  if ( m == 0 ) then  ! only axisymmetric modes
                     z(lm,nR)=z(lm,nR)+cmplx(c_r,0.0_cp,kind=cp)
                  else
                     z(lm,nR)=z(lm,nR)+cmplx(c_r,c_i,kind=cp)
                  end if
               end if
               !end if
              !write(output_unit,"(A,4I4,2ES20.12)") "z = ",nR,lm,l,m,z_Rloc(lm,nR)
            end do
 !           end if
         end do !end of the radial loop
      end if

      if (pertur_z  < -2 ) then
         !--- Add random noise toroidal field of all (l,m) modes exept (l=0,m=0):
         !    It decays likes l**(init_v1-1)
         !    Amplitude is chosen so that the (1,0) term resembles amp_v1 *
         !    the 'solid body' rotation set by inner core and mantle rotation.
         rExp=4.
         write(output_unit,*) "test", rank
         if ( omega_ic1 /= 0 ) then
            amp_r=amp_v1*omega_ic1*r_ICB**(rExp+1.)/y10_norm
         else
            amp_r=amp_v1*r_ICB**(rExp+1.)/y10_norm
         end if
            !write(output_unit,"(A,I3,A,ES20.12,ES20.12)") "rDep(",nR,") = ",rDep(nR)
         do lm=llm,ulm
            l=lo_map%lm2l(lm)
            m=lo_map%lm2m(lm)
            if ( l /= 0 ) then
               if ( m == 0 ) then  ! only axisymmetric modes
                  do nR=1,n_r_max
                     call random_number(ra1)
                     call random_number(ra2)
                     ra1=(-one+two*ra1)/(real(l,cp))**(-pertur_z-1)
                     ra2=(-one+two*ra2)/(real(l,cp))**(-pertur_z-1)
                     rDep(nR)=amp_r*(r(nR)**2)*sin((r(nR)-r_icb)*PI)
                     c_r=ra1*rDep(nR)
                     c_i=ra2*rDep(nR)
                     z(lm,nR)=z(lm,nR)+cmplx(c_r,0.0_cp,kind=cp)
                     w(lm,nR)=w(lm,nR)+cmplx(c_i,0.0_cp,kind=cp)
                  end do
               end if
            end if
              !write(output_unit,"(A,4I4,2ES20.12)") "z = ",nR,lm,l,m,z_Rloc(lm,nR)
         end do
        !end of the radial loop
      end if

      if (pertur_w<-2 .or. (pertur_w ==1)) then
         !--- Add random noise poloidal field of all (l,m) modes exept (l=0,m=0):
         !    It decays likes l**(-init_v1+1)
         !    Amplitude is chosen to be comparable to amp * inner core roation speed
         !    at inner core boundary...
         !amp_r=random(rand_num*rank/n_procs+one)
         if ( omega_ic1 /= 0.0_cp .and.  r_cut_model >0.3) then
            amp_r=amp_v1*omega_ic1*r_icb*r_icb/(y10_norm*PI)
         else
            amp_r=amp_v1*0.05_cp
         end if
         if (pertur_w ==1) pertur_w =-3

         !write(output_unit,*) amp_R,"TEST"
         do nR=1,n_r_max
            rDep(nR)=-amp_r*sin( (r_cmb-r(nR))*PI)!-r_ICB)*PI ))/(1+5*sin( (r(nR)-r_ICB)*PI ))
            do lm=llm,ulm
               l=lo_map%lm2l(lm)
               m=lo_map%lm2m(lm)
               if (l/=0) then
                  call random_number(ra1)
                  call random_number(ra2)
                  ra1=(-one+two*ra1)/(real(l,cp))**(-pertur_w-1)
                  ra2=(-one+two*ra2)/(real(l,cp))**(-pertur_w-1)
                  c_r=ra1*rDep(nR)
                  c_i=ra2*rDep(nR)
                  if ( m == 0 ) then  ! no axisymmetric modes
                     w(lm,nR)=w(lm,nR)+cmplx(c_r,0.0_cp,kind=cp)
                     z(lm,nR)=z(lm,nR)+cmplx(c_r,0.0_cp,kind=cp)
                  else
                     w(lm,nR)=w(lm,nR)+cmplx(c_r,c_i,kind=cp)
                     z(lm,nR)=z(lm,nR)+cmplx(c_r,c_i,kind=cp)
                  end if
               end if
            end do
         end do !end of the radial loop
      end if


      !----- Caring for IC and mantle rotation rates if this
      !      has not been done already in read_start_file.f:
      if ( ( .not. l_start_file ) ) then

         l1m0 = lo_map%lm2(1,0)

         if ( (l1m0>=llm) .and. (l1m0<=ulm) ) then

            write(output_unit,*) '! NO STARTFILE READ, SETTING Z10!'

            if ( (l_SRIC .and. .not. (kbotv==3) ).or. l_rot_ic .and. omega_ic1 /= 0.0_cp ) then
               omega_ic=omega_ic1*cos(omegaOsz_ic1*tShift_ic1) + &
               &        omega_ic2*cos(omegaOsz_ic2*tShift_ic2)
               write(output_unit,*)
               write(output_unit,*) '! I use prescribed inner core rotation rate:'
               write(output_unit,*) '! omega_ic=',omega_ic
               z(l1m0,n_r_icb)=cmplx(omega_ic/c_z10_omega_ic,kind=cp)
            else if ( l_rot_ic .and. omega_ic1 == 0.0_cp ) then
               omega_ic=c_z10_omega_ic*real(z(l1m0,n_r_icb))
            else
               omega_ic=0.0_cp
            end if
            if ( (l_SRMA .and. .not. (ktopv==3)) .or. (l_rot_ma .and. omega_ma1 /= 0.0_cp) ) then
               omega_ma=omega_ma1*cos(omegaOsz_ma1*tShift_ma1) + &
               &        omega_ma2*cos(omegaOsz_ma2*tShift_ma2)

               write(output_unit,*)
               write(output_unit,*) '! I use prescribed mantle rotation rate:'
               write(output_unit,*) '! omega_ma=',omega_ma
               z(l1m0,n_r_cmb)=cmplx(omega_ma/c_z10_omega_ma,kind=cp)
            else if ( l_rot_ma .and. omega_ma1 == 0.0_cp ) then
               omega_ma=c_z10_omega_ma*real(z(l1m0,n_r_cmb))
            else
               omega_ma=0.0_cp
            end if
         end if

#ifdef WITH_MPI
         if ( m_min == 0 ) then
            call MPI_Bcast(omega_ic,1,MPI_DEF_REAL,rank_with_l1m0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(omega_ma,1,MPI_DEF_REAL,rank_with_l1m0,MPI_COMM_WORLD,ierr)
         end if
#endif
      else
         if ( nRotIc == 2 ) omega_ic=omega_ic1
         if ( nRotMa == 2 ) omega_ma=omega_ma1
      end if

   end subroutine initV
!--------------------------------------------------------------------
   subroutine initS(s,p)
      !
      ! Purpose of this subroutine is to initialize the entropy field
      ! according to the input control parameters.
      !
      ! +-----------------+---------------------------------------------+
      ! | Input           | value                                       |
      ! +=================+=============================================+
      ! | init_s1 < 100:  | random noise initialized                    |
      ! |                 | the noise spectrum decays as l ^ (init_s1-1)|
      ! |                 | with peak amplitude amp_s1  for l=1         |
      ! +-----------------+---------------------------------------------+
      ! | init_s1 >=100:  | a specific harmonic mode initialized        |
      ! |                 | with amplitude amp_s1.                      |
      ! |                 | init_s1 is interpreted as number llmm       |
      ! |                 | where ll: harmonic degree,                  |
      ! |                 | mm: harmonic order.                         |
      ! +-----------------+---------------------------------------------+
      ! | init_s2 >100 :  | a second harmonic mode initialized          |
      ! |                 | with amplitude amp_s2.                      |
      ! |                 | init_s2 is again interpreted as number llmm |
      ! |                 | where ll: harmonic degree,                  |
      ! |                 | mm: harmonic order.                         |
      ! +-----------------+---------------------------------------------+
      !

      !-- Output variables:
      complex(cp), intent(inout) :: s(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: p(llm:ulm,n_r_max)

      !-- Local variables:
      integer :: n_r,lm,l,m,lm00
      real(cp) :: x,c_r,c_i,s_r,s_i
      real(cp) :: ra1,ra2
      real(cp) :: s0(n_r_max),p0(n_r_max),s1(n_r_max)

      integer :: nTheta,nPhi,nS
      real(cp) :: xL,yL,zL,rH,angleL,s00,s00P
      real(cp) :: mata(n_impS_max,n_impS_max)
      real(cp) :: amp(n_impS_max)
      integer :: pivot(n_impS_max)
      real(cp) :: xS(n_impS_max),yS(n_impS_max)
      real(cp) :: zS(n_impS_max),sFac(n_impS_max)
      real(cp) :: sCMB(nlat_padded,n_phi_max)
      complex(cp) :: sLM(lm_max)
      integer :: info,i,j,filehandle
      logical :: rank_has_l0m0

      lm00=lo_map%lm2(0,0)
      rank_has_l0m0=.false.

      if ( lm00 >= llm .and. lm00 <= ulm ) then
         rank_has_l0m0=.true.
      end if

      if ( (.not. l_start_file) .and. (.not. l_non_adia) ) then

         if ( rank_has_l0m0 ) then

            open(newunit=filehandle, file='scond.dat')
            if ( l_anelastic_liquid ) then
               call pt_cond(s0,p0)
               do n_r=1,n_r_max
                  write(filehandle,'(5ES20.12)') r(n_r), osq4pi*otemp1(n_r)* &
                  &            (s0(n_r)-ViscHeatFac*ThExpNb*alpha0(n_r)*     &
                  &            temp0(n_r)*orho1(n_r)*p0(n_r)),               &
                  &            osq4pi*p0(n_r), osq4pi*s0(n_r),               &
                  &            osq4pi*alpha0(n_r)*(-rho0(n_r)*s0(n_r)+       &
                  &            ViscHeatFac*ThExpNb*(alpha0(n_r)*temp0(n_r)   &
                  &            +ogrun(n_r))*p0(n_r))
               end do
            else
               call ps_cond(s0,p0)
               do n_r=1,n_r_max
                  write(filehandle,'(5ES20.12)') r(n_r), s0(n_r)*osq4pi, &
                  &            p0(n_r)*osq4pi, osq4pi*temp0(n_r)*(       &
                  &            s0(n_r)+alpha0(n_r)*orho1(n_r)*p0(n_r)*   &
                  &            ThExpNb*ViscHeatFac), osq4pi*alpha0(n_r)* &
                  &            ThExpNb*(-rho0(n_r)*temp0(n_r)*s0(n_r)+   &
                  &            ViscHeatFac*ogrun(n_r)*p0(n_r))
               end do
            end if
            close(filehandle)
            do n_r=1,n_r_max
               s(lm00,n_r)=s0(n_r)
               p(lm00,n_r)=p0(n_r)
            end do

         end if

      end if

      !-- Radial dependence of perturbation in s1:
      do n_r=1,n_r_max
         x=two*r(n_r)-r_cmb-r_icb
         s1(n_r)=one-three*x**2+three*x**4-x**6
      end do

      if ( l_onset ) then
         !-- same amplitude on all modes
         do n_r=1,n_r_max
            do lm=llm,ulm ! all modes except l=m=0 carry the same function
               l = lo_map%lm2l(lm)
               if ( l == 0 ) cycle
               s(lm,n_r)=amp_s1*s1(n_r)
            end do
         end do

      else if ( init_s1 < 100 .and. init_s1 > 0 ) then

      !-- Random noise initialization of all (l,m) modes exept (l=0,m=0):

         do lm=llm,ulm
            l = lo_map%lm2l(lm)
            if ( l == 0 ) cycle
            m = lo_map%lm2m(lm)
            call random_number(ra1)
            call random_number(ra2)
            ra1=(-one+two*ra1)*amp_s1/(real(l,cp))**(init_s1-1)
            ra2=(-one+two*ra2)*amp_s1/(real(l,cp))**(init_s1-1)
            do n_r=1,n_r_max
               c_r=ra1*s1(n_r)
               c_i=ra2*s1(n_r)
               if ( m > 0 ) then  ! non axisymmetric modes
                  s(lm,n_r)=s(lm,n_r)+cmplx(c_r,c_i,kind=cp)
               else
                  s(lm,n_r)=s(lm,n_r)+cmplx(c_r,0.0_cp,kind=cp)
               end if
            end do
         end do

      else  if ( init_s1 >= 100 ) then

      !-- Initialize one or two modes specifically

      !----- Initialize first mode:
         l=init_s1/100
         if ( l > 99 ) l=init_s1/1000
         m=mod(init_s1,100)
         if ( l > 99 ) m=mod(init_s1,1000)
         if ( mod(m,minc) /= 0 ) then
            write(output_unit,*) '! Wave number of mode for entropy initialisation'
            write(output_unit,*) '! not compatible with phi-symmetry:',m
            call abortRun('Stop run in init')
         end if
         if ( l > l_max .or. l < m ) then
            write(output_unit,*) '! Degree of mode for entropy initialisation'
            write(output_unit,*) '! > l_max or < m !',l
            call abortRun('Stop run in init')
         end if
         lm=lo_map%lm2(l,m)
         if ( (lm>=llm) .and. (lm<=ulm) ) then
            do n_r=1,n_r_max
               c_r=s1(n_r)*amp_s1
               s(lm,n_r)=s(lm,n_r)+cmplx(c_r,0.0_cp,kind=cp)
            end do

            write(output_unit,'(/'' ! Entropy initialized at mode:'', &
            &      '' l='',i4,'' m='',i4,'' Ampl='',f8.5)') l,m,amp_s1
         end if

         !----- Initialize second mode:
         if ( init_s2 > 99 ) then
            m=mod(init_s2,100)
            if ( mod(m,minc) /= 0 ) then
               write(output_unit,*) '! Wave number of mode for entropy initialisation'
               write(output_unit,*) '! not compatible with phi-symmetry:',m
               call abortRun('Stop run in init')
            end if
            l=init_s2/100
            if ( l > l_max .or. l < m ) then
               write(output_unit,*) '! Degree of mode for entropy initialisation'
               write(output_unit,*) '! > l_max or < m !',l
               call abortRun('Stop run in init')
            end if

            lm=lo_map%lm2(l,m)
            s_r=amp_s2
            s_i=0.0_cp
            if ( (lm>=llm) .and. (lm<=ulm) ) then
               if ( amp_s2 < 0.0_cp .and. m /= 0 ) then
               !-------- Sin(phi)-mode initialized for amp_s2<0
                  s_r = 0.0_cp
                  s_i = amp_s2
               end if
               do n_r=1,n_r_max
                  c_r=s1(n_r)*s_r
                  c_i=s1(n_r)*s_i
                  s(lm,n_r)=s(lm,n_r)+cmplx(c_r,c_i,kind=cp)
               end do
               write(output_unit,'('' ! Second mode:'', &
               &       '' l='',i3,'' m='',i3,'' Ampl='',f8.5/)') l,m,amp_s2
            end if

         end if

      end if

      if ( impS == 0 ) return

      !-- Now care for the prescribed boundary condition:

      if ( minc /= 1 ) call abortRun('! impS doesnt work for minc /= 1')

      if ( abs(impS) == 1 ) then
         n_impS=2
         peakS(2)=-peakS(1)
         thetaS(2)=pi-thetaS(1)
         phiS(2)  =pi+phiS(1)
         if ( phiS(2) > 2*pi ) phiS(2)=phiS(2)-2*pi
         widthS(2)=widthS(1)
      end if

      !-- Determine the peak value vector in (xS,yS,zS) space.
      !       Then get the proportionality factors for the linear dependence
      !       of the mean (l=0,m=0) contribution on the total peak amplitude
      !       amp:
      do nS=1,n_impS

         xS(nS)=sin(thetaS(nS))*cos(phiS(nS))
         yS(nS)=sin(thetaS(nS))*sin(phiS(nS))
         zS(nS)=cos(thetaS(nS))

         do nPhi=1,n_phi_max
            do nTheta=1,n_theta_max
               xL=sinTheta(nTheta)*cos(phi(nPhi))
               yL=sinTheta(nTheta)*sin(phi(nPhi))
               zL=cosTheta(nTheta)
               rH=sqrt((xS(nS)-xL)**2 + (yS(nS)-yL)**2+(zS(nS)-zL)**2)
               !------ Opening angleL with peak value vector:
               angleL=two*abs(asin(rH/2))
               if ( angleL <= widthS(nS) ) then
                  sCMB(nTheta,nPhi)=(cos(angleL/widthS(nS)*pi)+1)/2
               else
                  sCMB(nTheta,nPhi)=0.0_cp
               end if
            end do
         end do
         call scal_to_SH(sCMB, sLM, l_max)

         !--- sFac describes the linear dependence of the (l=0,m=0) mode
         !    on the amplitude peakS, SQRT(4*pi) is a normalisation factor
         !    according to the spherical harmonic function form chosen here.
         sFac(nS)=real(sLM(st_map%lm2(0,0)))*osq4pi

      end do ! Loop over peak

      !-- Value due to prescribed (l=0,m=0) contribution
      s00P=real(tops(0,0))*osq4pi
      if ( s00P == 0.0_cp .and. impS < 0 ) then
         write(output_unit,*) '! No relative amplitudes possible!'
         write(output_unit,*) '! for impS<0 because the mean value!'
         write(output_unit,*) '! is zero! Refince s_top?'
         call abortRun('Stop run in init')
      end if
      if ( impS > 0 ) s00P=one

      !-- Determine the true amplitudes amp for the peaks by solving linear system:
      !    These amplitudes guarantee that the peak has an ampliture peakS
      !    above or below the mean (l=0,m=0)
      if ( n_impS == 1 ) then
         amp(1)=peakS(1)/(s00P*(one-sFac(1)))
      else
         do j=1,n_impS
            amp(j)=-peakS(j)/s00P
            do i=1,n_impS
               if ( i == j ) then
                  mata(i,j)=sFac(i)-1
               else
                  mata(i,j)=sFac(i)
               end if
            end do
         end do
        call prepare_mat(mata,n_impS_max,n_impS,pivot,info)
        call solve_mat(mata,n_impS_max,n_impS,pivot,amp)
      end if
      s00=0.0_cp
      do nS=1,n_impS
         s00=s00+sFac(nS)*amp(nS)
      end do

      !--- Now get the total thing so that the mean (l=0,m=0) due
      !    to the peaks is zero. The (l=0,m=0) contribution is
      !    determined (prescribed) by other means.

      do nPhi=1,n_phi_max
         do nTheta=1,n_theta_max
            xL=sinTheta(nTheta)*cos(phi(nPhi))
            yL=sinTheta(nTheta)*sin(phi(nPhi))
            zL=cosTheta(nTheta)
            sCMB(nTheta,nPhi)=-s00
            do nS=1,n_impS
               rH=sqrt((xS(nS)-xL)**2 + (yS(nS)-yL)**2+(zS(nS)-zL)**2)
               !------ Opening angle with peak value vector:
               angleL=two*abs(asin(rH/2))
               if ( angleL <= widthS(nS) )              &
               &  sCMB(nTheta,nPhi)=sCMB(nTheta,nPhi) + &
               &                    amp(nS)*(cos(angleL/widthS(nS)*pi)+1)/2
            end do
         end do
      end do

      call scal_to_SH(sCMB, sLM, l_max)

      !--- Finally store the boundary condition and care for
      !    the fact that peakS provides the relative amplitudes
      !    in comparison to the (l=0,m=0) contribution when impS<0:
      !    Note that the (l=0,m=0) has to be determined by other means
      !    for example by setting: s_top= 0 0 -1 0
      do lm=1,lm_max
         l = st_map%lm2l(l)
         m = st_map%lm2m(m)
         if ( l <= l_max .and. l > 0 ) tops(l,m)=tops(l,m)+sLM(lm)
      end do

   end subroutine initS
!---------------------------------------------------------------------------
   subroutine initXi(xi)
      !
      ! Purpose of this subroutine is to initialize the chemical composition
      ! according to the input control parameters.
      !
      ! +-----------------+---------------------------------------------+
      ! | Input           | value                                       |
      ! +=================+=============================================+
      ! | init_xi1 < 100: | random noise initialized                    |
      ! |                 | the noise spectrum decays as l^ (init_xi1-1)|
      ! |                 | with peak amplitude amp_xi1  for l=1        |
      ! +-----------------+---------------------------------------------+
      ! | init_xi1 >=100: | a specific harmonic mode initialized        |
      ! |                 | with amplitude amp_xi1.                     |
      ! |                 | init_xi1 is interpreted as number llmm      |
      ! |                 | where ll: harmonic degree,                  |
      ! |                 | mm: harmonic order.                         |
      ! +-----------------+---------------------------------------------+
      ! | init_xi2 >100 : | a second harmonic mode initialized          |
      ! |                 | with amplitude amp_xi2.                     |
      ! |                 | init_xi2 is again interpreted as number llmm|
      ! |                 | where ll: harmonic degree,                  |
      ! |                 | mm: harmonic order.                         |
      ! +-----------------+---------------------------------------------+
      !

      !-- Output variables:
      complex(cp), intent(inout) :: xi(llm:ulm,n_r_max)

      !-- Local variables:
      integer :: n_r,lm,l,m,lm00
      real(cp) :: x,c_r,c_i,xi_r,xi_i
      real(cp) :: ra1,ra2
      real(cp) :: xi0(n_r_max),xi1(n_r_max)

      integer :: nTheta,nPhi,nXi
      real(cp) :: xL,yL,zL,rH,angleL,xi00,xi00P
      real(cp) :: mata(n_impXi_max,n_impXi_max)
      real(cp) :: amp(n_impXi_max)
      integer :: pivot(n_impXi_max)
      real(cp) :: xXi(n_impXi_max),yXi(n_impXi_max)
      real(cp) :: zXi(n_impXi_max),xiFac(n_impXi_max)
      real(cp) :: xiCMB(nlat_padded,n_phi_max)
      complex(cp) :: xiLM(lm_max)
      integer :: info,i,j,fileHandle

      lm00=lo_map%lm2(0,0)

      if ( .not. l_start_file ) then

         if ( (llm <= lm00) .and. (ulm >= lm00) ) then
            call xi_cond(xi0)
            open(newunit=fileHandle, file='xicond.dat')
            do n_r=1,n_r_max
               xi(lm00,n_r)=xi0(n_r)
               write(fileHandle,*) r(n_r), xi0(n_r)*osq4pi
            end do
            close(fileHandle)
         end if

      end if

      !-- Radial dependence of perturbation in xi1:
      do n_r=1,n_r_max
         x=two*r(n_r)-r_cmb-r_icb
         xi1(n_r)=one-three*x**2+three*x**4-x**6
      end do

      if ( l_onset ) then
         !-- same amplitude on all modes
         do n_r=1,n_r_max
            do lm=llm,ulm ! all modes except l=m=0 carry the same function
               l = lo_map%lm2l(lm)
               if ( l == 0 ) cycle
               xi(lm,n_r)=amp_xi1*xi1(n_r)
            end do
         end do

      else if ( init_xi1 < 100 .and. init_xi1 > 0 ) then

      !-- Random noise initialization of all (l,m) modes exept (l=0,m=0):

         do lm=llm,ulm
            l = lo_map%lm2l(lm)
            if ( l == 0 ) cycle
            m = lo_map%lm2m(lm)
            call random_number(ra1)
            call random_number(ra2)
            ra1=(-one+two*ra1)*amp_xi1/(real(l,cp))**(init_xi1-1)
            ra2=(-one+two*ra2)*amp_xi1/(real(l,cp))**(init_xi1-1)
            do n_r=1,n_r_max
               c_r=ra1*xi1(n_r)
               c_i=ra2*xi1(n_r)
               if ( m > 0 ) then  ! non axisymmetric modes
                  xi(lm,n_r)=xi(lm,n_r)+cmplx(c_r,c_i,kind=cp)
               else
                  xi(lm,n_r)=xi(lm,n_r)+cmplx(c_r,0.0_cp,kind=cp)
               end if
            end do
         end do

      else if ( init_xi1 >= 100 ) then

      !-- Initialize one or two modes specifically

      !----- Initialize first mode:
         l=init_xi1/100
         if ( l > 99 ) l=init_xi1/1000
         m=mod(init_xi1,100)
         if ( l > 99 ) m=mod(init_xi1,1000)
         if ( mod(m,minc) /= 0 ) then
            write(output_unit,*) '! Wave number of mode for chemical composition initialisation'
            write(output_unit,*) '! not compatible with phi-symmetry:',m
            call abortRun('Stop run in init')
         end if
         if ( l > l_max .or. l < m ) then
            write(output_unit,*) '! Degree of mode for chemical composition initialisation'
            write(output_unit,*) '! > l_max or < m !',l
            call abortRun('Stop run in init')
         end if
         lm=lo_map%lm2(l,m)

         if ( (llm <= lm) .and. (ulm >= lm) ) then
            do n_r=1,n_r_max
               c_r=xi1(n_r)*amp_xi1
               xi(lm,n_r)=xi(lm,n_r)+cmplx(c_r,0.0_cp,kind=cp)
            end do

            write(output_unit,'(/'' ! Chemical composition initialized at mode:'', &
                &  '' l='',i4,'' m='',i4,'' Ampl='',f8.5)') l,m,amp_s1
         end if

      !----- Initialize second mode:
         if ( init_xi2 > 99 ) then
            m=mod(init_xi2,100)
            if ( mod(m,minc) /= 0 ) then
               write(output_unit,*) '! Wave number of mode for chemical composition initialisation'
               write(output_unit,*) '! not compatible with phi-symmetry:',m
               call abortRun('Stop run in init')
            end if
            l=init_xi2/100
            if ( l > l_max .or. l < m ) then
               write(output_unit,*) '! Degree of mode for chemical composition initialisation'
               write(output_unit,*) '! > l_max or < m !',l
               call abortRun('Stop run in init')
            end if

            lm=lo_map%lm2(l,m)
            if ( (llm <= lm) .and. (ulm >= lm) ) then
               xi_r=amp_s2
               xi_i=0.0_cp
               if ( amp_s2 < 0.0_cp .and. m /= 0 ) then
               !-------- Sin(phi)-mode initialized for amp_xi2<0
                  xi_r = 0.0_cp
                  xi_i = amp_s2
               end if
               do n_r=1,n_r_max
                  c_r=xi1(n_r)*xi_r
                  c_i=xi1(n_r)*xi_i
                  xi(lm,n_r)=xi(lm,n_r)+cmplx(c_r,c_i,kind=cp)
               end do
               write(output_unit,'('' ! Second mode:'', &
               &     '' l='',i3,'' m='',i3,'' Ampl='',f8.5/)') l,m,amp_xi2
            end if

         end if

      end if

      if ( impXi == 0 ) return

      !-- Now care for the prescribed boundary condition:
      if ( minc /= 1 ) call abortRun('! impXi doesnt work for minc /= 1')

      if ( abs(impXi) == 1 ) then
         n_impXi=2
         peakXi(2)=-peakXi(1)
         thetaXi(2)=pi-thetaXi(1)
         phiXi(2)  =pi+phiXi(1)
         if ( phiXi(2) > 2*pi ) phiXi(2)=phiXi(2)-2*pi
         widthXi(2)=widthXi(1)
      end if

      !-- Determine the peak value vector in (xXi,yXi,zXi) space.
      !     Then get the proportionality factors for the linear dependence
      !     of the mean (l=0,m=0) contribution on the total peak amplitude
      !     amp:
      do nXi=1,n_impXi

         xXi(nXi)=sin(thetaXi(nXi))*cos(phiXi(nXi))
         yXi(nXi)=sin(thetaXi(nXi))*sin(phiXi(nXi))
         zXi(nXi)=cos(thetaXi(nXi))

         do nPhi=1,n_phi_max
            do nTheta=1,n_theta_max
               xL=sinTheta(nTheta)*cos(phi(nPhi))
               yL=sinTheta(nTheta)*sin(phi(nPhi))
               zL=cosTheta(nTheta)
               rH=sqrt((xXi(nXi)-xL)**2 + (yXi(nXi)-yL)**2+(zXi(nXi)-zL)**2)
               !------ Opening angleL with peak value vector:
               angleL=two*abs(asin(rH/2))
               if ( angleL <= widthXi(nXi) ) then
                  xiCMB(nTheta,nPhi) = half*(cos(angleL/widthXi(nXi)*pi)+1)
               else
                  xiCMB(nTheta,nPhi)=0.0_cp
               end if
            end do
         end do
         call scal_to_SH(xiCMB, xiLM, l_max)

         !--- xiFac describes the linear dependence of the (l=0,m=0) mode
         !    on the amplitude peakXi, sqrt(4*pi) is a normalisation factor
         !    according to the spherical harmonic function form chosen here.
         xiFac(nXi)=real(xiLM(st_map%lm2(0,0)))*osq4pi

      end do ! Loop over peak

      !-- Value due to prescribed (l=0,m=0) contribution
      xi00P=real(topxi(0,0))*osq4pi
      if ( xi00P == 0.0_cp .and. impXi < 0 ) then
         write(output_unit,*) '! No relative amplitudes possible!'
         write(output_unit,*) '! for impXi<0 because the mean value!'
         write(output_unit,*) '! is zero! Refince xi_top?'
         call abortRun('Stop run in init')
      end if
      if ( impXi > 0 ) xi00P=one

      !-- Determine the true amplitudes amp for the peaks by solving linear system:
      !    These amplitudes guarantee that the peak as an ampliture peakXi
      !    above or below the mean (l=0,m=0)
      if ( n_impXi == 1 ) then
         amp(1)=peakXi(1)/(xi00P*(one-xiFac(1)))
      else
         do j=1,n_impXi
            amp(j)=-peakXi(j)/xi00P
            do i=1,n_impXi
               if ( i == j ) then
                  mata(i,j)=xiFac(i)-1
               else
                  mata(i,j)=xiFac(i)
               end if
            end do
         end do
        call prepare_mat(mata,n_impXi_max,n_impXi,pivot,info)
        call solve_mat(mata,n_impXi_max,n_impXi,pivot,amp)
      end if
      xi00=0.0_cp
      do nXi=1,n_impXi
         xi00=xi00+xiFac(nXi)*amp(nXi)
      end do

      !--- Now get the total thing so that the mean (l=0,m=0) due
      !    to the peaks is zero. The (l=0,m=0) contribution is
      !    determined (prescribed) by other means.
      do nPhi=1,n_phi_max
         do nTheta=1,n_theta_max
            xL=sinTheta(nTheta)*cos(phi(nPhi))
            yL=sinTheta(nTheta)*sin(phi(nPhi))
            zL=cosTheta(nTheta)
            xiCMB(nTheta,nPhi)=-xi00
            do nXi=1,n_impXi
               rH=sqrt((xXi(nXi)-xL)**2 + (yXi(nXi)-yL)**2+(zXi(nXi)-zL)**2)
               !------ Opening angle with peak value vector:
               angleL=two*abs(asin(rH/2))
               if ( angleL <= widthXi(nXi) )              &
                  xiCMB(nTheta,nPhi)=xiCMB(nTheta,nPhi) + &
                                     amp(nXi)*half*(cos(angleL/widthXi(nXi)*pi)+1)
            end do
         end do
      end do

      call scal_to_SH(xiCMB, xiLM, l_max)

      !--- Finally store the boundary condition and care for
      !    the fact that peakS provides the relative amplitudes
      !    in comparison to the (l=0,m=0) contribution when impS<0:
      !    Note that the (l=0,m=0) has to be determined by other means
      !    for example by setting: s_top= 0 0 -1 0
      do lm=1,lm_max
         l = st_map%lm2l(l)
         m = st_map%lm2m(m)
         if ( l <= l_max .and. l > 0 ) topxi(l,m)=topxi(l,m)+xiLM(lm)
      end do

   end subroutine initXi
!---------------------------------------------------------------------------
   subroutine initB(b, aj, b_ic, aj_ic)
      !
      ! Purpose of this subroutine is to initialize the magnetic field
      ! according to the control parameters imagcon and init_b1/2.
      ! In addition CMB and ICB peak values are calculated for
      ! magneto convection.
      !

      !-- Output variables:
      complex(cp), intent(inout) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(inout) :: aj(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(inout) :: b_ic(llmMag:ulmMag,n_r_ic_max)
      complex(cp), intent(inout) :: aj_ic(llmMag:ulmMag,n_r_ic_max)

      !-- Local variables:
      integer :: lm,lm0,l1,m1,l,m
      integer :: n_r, n_r_fourier
      real(cp) :: b_pol,b_tor
      complex(cp) :: aj0(n_r_max+1)
      complex(cp) :: aj0_ic(n_r_ic_max)
      real(cp) :: arg,aj_ic1,aj_ic2

      real(cp) :: b1(n_r_max)
      real(cp) :: b1_ic(n_r_ic_max)
      real(cp) :: bR,bI
      real(cp) :: aVarCon,bVarCon
      integer :: bExp

      integer :: l1m0,l2m0,l3m0,l1m1
      real(cp) :: x_smooth, l_smooth, l_corr
      real(cp) :: e_mag_p,global_e_mag_p
      real(cp) :: e_mag_p_r(n_r_max),db1(n_r_max)

      real(cp) :: alpha,beta1,c,x2,ss,r_cut
      complex(cp) :: b1_Rloc(lm_max,nRstart:nRstop)

      integer :: nTheta,n,nPhi,st_lmP

      l1m0 = lo_map%lm2(1,0)
      l2m0 = lo_map%lm2(2,0)
      l3m0 = lo_map%lm2(3,0)
      l1m1 = lo_map%lm2(1,1)

      lm0=l2m0 ! Default quadrupole field

      if ( imagcon == -1 ) then

         !----- impose l=1,m=0 poloidal field at ICB:
         lm0 = l1m0
         bpeakbot = -sqrt(third*pi)*r_icb**2*amp_b1
         bpeaktop = 0.0_cp

      else if ( imagcon == -2 ) then

         !----- impose l=1,m=0 poloidal field at CMB:
         lm0 = l1m0
         bpeakbot = 0.0_cp
         bpeaktop = -sqrt(third*pi)*r_cmb**2*amp_b1

      else if ( imagcon == 1 ) then

         !----- impose l=2,m=0 toroidal field at ICB:
         lm0 = l2m0
         bpeakbot = four*third*sqrt(pi/5.0_cp)*r_icb*amp_b1
         bpeaktop = 0.0_cp

      else if ( imagcon == 10 ) then

         !----- impose l=2,m=0 toroidal field at ICB and CMB:
         lm0 = l2m0
         bpeakbot = four*third*sqrt(pi/5.0_cp)*r_icb*amp_b1
         bpeaktop = four*third*sqrt(pi/5.0_cp)*r_cmb*amp_b1

      else if ( imagcon == 11 ) then

         !----- same as imagcon == 10 but opposite sign at CMB:
         lm0 = l2m0
         bpeakbot = four*third*sqrt(pi/5.0_cp)*r_icb*amp_b1
         bpeaktop = -four*third*sqrt(pi/5.0_cp)*r_cmb*amp_b1

      else if ( imagcon == 12 ) then

         !----- impose l=1,m=0 toroidal field at ICB and CMB:
         lm0 = l1m0
         bpeakbot = two*sqrt(third*pi)*r_icb*amp_b1
         bpeaktop = two*sqrt(third*pi)*r_cmb*amp_b1

      else if ( imagcon == 0 ) then

         lm0 = l2m0
         bpeakbot = 0.0_cp
         bpeaktop = 0.0_cp

      else if ( imagcon == -10 ) then

         !----- Test of variable conductivity case with analytical solution:
         !      Assume the magnetic diffusivity is lambda=r**5, that the aspect ratio
         !      is 0.5, and that there is no flow.
         !      The analytical stationary solution for the (l=3,m=0) toroidal field
         !      with bounday condition aj(r=r_ICB)=1, aj(r=r_CMB)=0 is then
         !      given by jVarCon(r)!
         !      A disturbed solution is used to initialize aj,
         !      the disturbance should decay with time.
         !      The solution is stored in file testVarCond.TAG at the end of the run,
         !      where the first column denotes radius, the second is aj(l=3,m=0,r) and
         !      the third is jVarCon(r). Second and third should be identical when
         !      the stationary solution has been reached.
         lm0=l3m0  ! This is l=3,m=0
         bpeakbot=one
         bpeaktop=0.0_cp
         aVarCon =-one/255.0_cp
         bVarCon =256.0_cp/255.0_cp
         if ( llm <= lm0 .and. ulm >= lm0 ) then ! select processor
            do n_r=1,n_r_max             ! Diffusive toroidal field
               jVarCon(n_r)=aVarCon*r(n_r)**2 + bVarCon/(r(n_r)**6)
               aj(lm0,n_r) =jVarCon(n_r) + 0.1_cp*sin((r(n_r)-r_ICB)*pi)
            end do
         end if

      end if

      if ( init_b1 == 1 .or. imagcon > 0 ) then
         !----- Conductive solution for toroidal field,
         !      diffusion equation solved in j_cond, amplitude defined
         !      by bpeaktop and bpeakbot respectively.
         !      bpeakbot is only used for insulating inner core !
         if ( llm <= lm0 .and. ulm >= lm0 ) then ! select processor
            call j_cond(lm0,aj0,aj0_ic)
            do n_r=1,n_r_max             ! Diffusive toroidal field
               aj(lm0,n_r)=aj0(n_r)
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  aj_ic(lm0,n_r)=aj0_ic(n_r)
               end do
            end if
         end if

      else if ( init_b1 == 2 ) then  ! l=1,m=0 analytical toroidal field
         ! with a maximum of amp_b1 at mid-radius
         ! between r_icb and r_cmb for an insulating
         ! inner core and at r_cmb/2 for a conducting
         ! inner core

         if ( llm <= l1m0 .and. ulm >= l1m0 ) then ! select processor
            b_tor=-two*amp_b1*sqrt(third*pi)  ! minus sign makes phi comp. > 0
            if ( l_cond_ic ) then
               do n_r=1,n_r_max
                  aj(l1m0,n_r)=aj(l1m0,n_r) + b_tor*r(n_r)*sin(pi*r(n_r)/r_cmb)
               end do
               do n_r=1,n_r_ic_max
                  aj_ic(l1m0,n_r)=aj_ic(l1m0,n_r) + &
                                  b_tor*r_ic(n_r)*sin(pi*r_ic(n_r)/r_cmb)
              end do
            else
               do n_r=1,n_r_max
                  aj(l1m0,n_r)=aj(l1m0,n_r) + b_tor*r(n_r)*sin(pi*(r(n_r)-r_icb))
               end do
            end if
         end if

      else if ( init_b1 == 3 ) then
         ! l=2,m=0 toroidal field and l=1,m=0 poloidal field
         ! toroidal field has again its maximum of amp_b1
         ! at mid-radius between r_icb and r_cmb for an
         ! insulating inner core and at r_cmb/2 for a
         ! conducting inner core
         ! The outer core poloidal field is defined by
         ! a homogeneous  current density, its maximum at
         ! the ICB is set to amp_b1.
         ! The inner core poloidal field is chosen accordingly.
         if ( llm <= l1m0 .and. ulm >= l1m0 ) then ! select processor
            b_tor=-four*third*amp_b1*sqrt(pi/5.0_cp)
            if ( l_cond_ic ) then
               b_pol=amp_b1*sqrt(three*pi)/(three+r_cmb)
               do n_r=1,n_r_max
                  b(l1m0,n_r)=b(l1m0,n_r) + &
                              b_pol*(r(n_r)**3 - four*third*r_cmb*r(n_r)**2)
               end do
               arg=pi*r_icb/r_cmb
               aj_ic1=(arg-two*sin(arg)*cos(arg)) / (arg+sin(arg)*cos(arg))
               aj_ic2=(one-aj_ic1)*r_icb*sin(arg)/cos(arg)
               do n_r=1,n_r_ic_max
                  b_ic(l1m0,n_r)=b_ic(l1m0,n_r)+b_pol*r_icb**2 * (    &
                  &                         half*r_ic(n_r)**2/r_icb + &
                  &                                      half*r_icb - &
                  &                                 four*third*r_cmb )
               end do
            else
               b_pol=amp_b1*sqrt(three*pi)/four
               do n_r=1,n_r_max
                  b(l1m0,n_r)=b(l1m0,n_r)+ b_pol *     (  &
                  &                          r(n_r)**3 -  &
                  &          four*third*r_cmb*r(n_r)**2 + &
                  &       third*r_icb**4*or1(n_r)    )
               end do
            end if
         end if

         if ( llm <= l2m0 .and. ulm >= l2m0 ) then ! select processor
            b_tor=-four*third*amp_b1*sqrt(pi/5.0_cp)
            if ( l_cond_ic ) then
               b_pol=amp_b1*sqrt(three*pi)/(three+r_cmb)
               do n_r=1,n_r_max
                  aj(l2m0,n_r)=aj(l2m0,n_r) + b_tor*r(n_r)*sin(pi*(r(n_r)/r_cmb))
               end do
               arg=pi*r_icb/r_cmb
               aj_ic1=(arg-two*sin(arg)*cos(arg)) / (arg+sin(arg)*cos(arg))
               aj_ic2=(one-aj_ic1)*r_icb*sin(arg)/cos(arg)
               do n_r=1,n_r_ic_max
                  aj_ic(l2m0,n_r)=aj_ic(l2m0,n_r)+b_tor*             ( &
                  &        aj_ic1*r_ic(n_r)*sin(pi*r_ic(n_r)/r_cmb) +  &
                  &                  aj_ic2*cos(pi*r_ic(n_r)/r_cmb) )
               end do
            else
               b_pol=amp_b1*sqrt(three*pi)/four
               do n_r=1,n_r_max
                  aj(l2m0,n_r)=aj(l2m0,n_r) + b_tor*r(n_r)*sin(pi*(r(n_r)-r_icb))
               end do
            end if
         end if

      else if ( init_b1 == 4 .or. imagcon == -1 ) then  ! l=1,m0 poloidal field
         ! with max field amplitude amp_b1 at r_icb
         if ( llm <= l1m0 .and. ulm >= l1m0 ) then ! select processor
            b_pol=-amp_b1*r_icb**3*sqrt(third*pi)
            do n_r=1,n_r_max
               b(l1m0,n_r)=b(l1m0,n_r)+b_pol*or1(n_r)
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  b_ic(l1m0,n_r)=b_ic(l1m0,n_r)+b_pol/r_icb* &
                                 ( -three*half + half*(r_ic(n_r)/r_icb)**2 )
               end do
            end if
         end if

      else if ( init_b1 == -3 ) then
         ! l=2,m=0 toroidal field and l=1,m=0 poloidal field
         ! toroidal field has again its maximum of amp_b1
         ! at mid-radius between r_icb and r_cmb for an
         ! insulating inner core and at r_cmb/2 for a
         ! conducting inner core
         ! The outer core poloidal field is defined by
         ! a homogeneous  current density, its maximum at
         ! the ICB is set to amp_b1.
         ! The inner core poloidal field is chosen accordingly.
         if ( llm <= l1m0 .and. ulm >= l1m0 ) then ! select processor
               b_pol=amp_b1*sqrt(three*pi)/four
               do n_r=1,n_r_max
                  b(l1m0,n_r)=b(l1m0,n_r)+ b_pol *     ( &
                                             r(n_r)**3 - &
                             four*third*r_cmb*r(n_r)**2 + &
                          third*r_icb**4/r(n_r)    )
               end do
         end if

         if ( llm <= l1m1 .and. ulm >= l1m1 ) then ! select processor
               b_pol=amp_b1*sqrt(three*pi)/four
               do n_r=1,n_r_max
                  b(l1m1,n_r)=b(l1m1,n_r)+ b_pol *     ( &
                                             r(n_r)**3 - &
                             four*third*r_cmb*r(n_r)**2 + &
                          third*r_icb**4/r(n_r)    )
               end do
        end if

      else if ( init_b1 == 5 ) then  ! l=1,m0 poloidal field
         ! constant j density, defined max field value amp_v1 at r_cmb
         if ( llm <= l1m0 .and. ulm >= l1m0 ) then ! select processor
            if ( l_cond_ic ) then
               b_pol=amp_b1*sqrt(three*pi)/r_cmb
               do n_r=1,n_r_max
                  b(l1m0,n_r)=b(l1m0,n_r)+b_pol* (   r(n_r)**3 - &
                              four*third*r_cmb * r(n_r)**2 )
               end do
               do n_r=1,n_r_ic_max
                  b_ic(l1m0,n_r)=b_ic(l1m0,n_r)+b_pol*r_icb**2 * &
                     (-5.0_cp/6.0_cp*r_icb-four*third+half*r_ic(n_r)**2/r_icb)
               end do
            else
               b_pol=amp_b1*sqrt(three*pi)/(r_cmb*(one-radratio**4))
               do n_r=1,n_r_max
                  b(l1m0,n_r)=b(l1m0,n_r)+b_pol* (   r(n_r)**3 - &
                                   four*third*r_cmb * r(n_r)**2 + &
                                   third*r_icb**4 / r(n_r)    )
               end do
            end if
         end if

      else if ( init_b1 == 6 ) then  ! l=1,m=0 poloidal field , constant in r !
         ! no potential at r_cmb but simple
         if ( llm <= l1m0 .and. ulm >= l1m0 ) then ! select processor
            b_pol=amp_b1
            do n_r=1,n_r_max
               b(l1m0,n_r)=b(l1m0,n_r)+b_pol*r(n_r)**2
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  b_ic(l1m0,n_r)=b_ic(l1m0,n_r)+b_pol*r_icb**2
               end do
            end if
         end if

      else if ( init_b1 == 7 .or. imagcon == -2 ) then  ! l=1,m0 poloidal field
         ! which is potential field at r_cmb
         if ( llm <= l1m0 .and. ulm >= l1m0 ) then ! select processor
            b_pol=amp_b1*5.0_cp*half*sqrt(third*pi)*r_icb**2
            do n_r=1,n_r_max
               b(l1m0,n_r)=b(l1m0,n_r)+b_pol*(r(n_r)/r_icb)**2 * &
                           ( one - three/5.0_cp*(r(n_r)/r_cmb)**2 )
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  b_ic(l1m0,n_r)=b_ic(l1m0,n_r)+b_pol * &
                                 ( one - three/5.0_cp*(r_ic(n_r)/r_cmb)**2 )
               end do
            end if
         end if

      else if ( init_b1 == 8 ) then  ! l=1,m0 pol. field, l=2,m=0 toroidal field
         ! which is potential field at r_cmb
         if ( llm <= l1m0 .and. ulm >= l1m0 ) then ! select processor
            b_pol=amp_b1*5.0_cp*half*sqrt(third*pi)*r_icb**2
            do n_r=1,n_r_max
               b(l1m0,n_r)=b(l1m0,n_r)+b_pol*(r(n_r)/r_cmb)**2 * &
                           ( one - three/5.0_cp*(r(n_r)/r_cmb)**2 )
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  b_ic(l1m0,n_r)=b_ic(l1m0,n_r)+b_pol * &
                                 ( one - three/5.0_cp*(r_ic(n_r)/r_cmb)**2 )
               end do
            end if
         end if

         if ( llm <= l2m0 .and. ulm >= l2m0 ) then ! select processor
            b_tor=amp_b1*three*half*sqrt(pi/5.0_cp)*r_icb**2*radratio
            do n_r=1,n_r_max
               aj(l2m0,n_r)=aj(l2m0,n_r)+b_tor*(r(n_r)/r_icb)**3 * &
                            ( one - (r(n_r)/r_cmb)**2 )
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  aj_ic(l2m0,n_r)=aj_ic(l2m0,n_r)+b_tor * &
                                  ( one - (r_ic(n_r)/r_cmb)**2 )
               end do
            end if
         end if

      else if ( init_b1 == 9 ) then  ! l=2,m0 poloidal field
         ! which is potential field at r_cmb
         if ( llm <= l2m0 .and. ulm >= l2m0 ) then ! select processor
            b_pol=amp_b1*7.0_cp/6.0_cp*sqrt(pi/5.0_cp)*r_icb**2*radratio
            do n_r=1,n_r_max
               b(l2m0,n_r)=b(l2m0,n_r)+b_pol*(r(n_r)/r_icb)**3 * &
                           ( one - 5.0_cp/7.0_cp*(r(n_r)/r_cmb)**2 )
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  b_ic(l2m0,n_r)=b_ic(l2m0,n_r)+b_pol * &
                                 ( one - 5.0_cp/7.0_cp*(r_ic(n_r)/r_cmb)**2 )
               end do
            end if
         end if

      else if ( init_b1 == 10 ) then  ! only equatorial dipole

         if ( l1m1 <= 0 ) then
            call abortRun('! Can not initialize l=1,m=1 !')
         end if

         if ( llm <= l1m1 .and. ulm >= l1m1 ) then ! select processor
            b_pol=amp_b1*5.0_cp*half*sqrt(third*pi)*r_icb**2
            do n_r=1,n_r_max
               b(l1m1,n_r)=b(l1m1,n_r)+b_pol*(r(n_r)/r_icb)**2 * &
                            ( one - three/5.0_cp*(r(n_r)/r_cmb)**2 )
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  b_ic(l1m1,n_r)=b_ic(l1m1,n_r)+b_pol * &
                                 ( one - three/5.0_cp*(r_ic(n_r)/r_cmb)**2 )
               end do
            end if
         end if

      else if ( init_b1 < 0 ) then  ! l,m mixture, random init

         bExp=abs(init_b1)
         do n_r=1,n_r_max
            b1(n_r)=(r(n_r)/r_cmb)**2 * ( one-three/5.0_cp*(r(n_r)/r_cmb)**2 )
         end do
         if ( l_cond_ic ) then
            do n_r=1,n_r_ic_max
               b1_ic(n_r)= ( one-three/5.0_cp*(r_ic(n_r)/r_cmb)**2 )
            end do
         end if

         !-- Random noise initialization of all (l,m) modes exept (l=0,m=0):
         do lm=llm,ulm
            l=lo_map%lm2l(lm)
            m=lo_map%lm2m(lm)
            if ( l > 0 ) then
               call random_number(bR)
               call random_number(bI)
               bR=(-one+two*bR)*amp_b1/(real(l,cp))**(bExp-1)
               bI=(-one+two*bI)*amp_b1/(real(l,cp))**(bExp-1)
            else
               bR=0.0_cp
               bI=0.0_cp
            end if
            if ( m == 0 ) bI=0.0_cp
            do n_r=1,n_r_max
               b(lm,n_r)=b(lm,n_r) + cmplx(bR*b1(n_r),bI*b1(n_r),kind=cp)
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  b_ic(lm,n_r)=b_ic(lm,n_r) + cmplx(bR*b1_ic(n_r), bI*b1_ic(n_r), cp)
               end do
            end if
         end do

      else if ( init_b1 == 11 ) then  ! axial and equatorial dipole

         if ( llm <= l1m0 .and. ulm >= l1m0 ) then ! select processor
            b_pol=amp_b1*5.0_cp*half*sqrt(third*pi)*r_icb**2
            do n_r=1,n_r_max
               b(l1m0,n_r)=b(l1m0,n_r)+b_pol*(r(n_r)/r_cmb)**2 * &
                           ( one - three/5.0_cp*(r(n_r)/r_cmb)**2 )
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  b_ic(l1m0,n_r)=b_ic(l1m0,n_r)+b_pol * &
                                 ( one - three/5.0_cp*(r_ic(n_r)/r_cmb)**2 )
               end do
            end if
         end if

         if ( l1m1 <= 0 ) then
            call abortRun('! Cannot initialize l=1,m=1 !')
         end if

         if ( llm <= l1m1 .and. ulm >= l1m1 ) then ! select processor
            b_pol=amp_b1*5.0_cp*half*sqrt(third*pi)*r_icb**2
            do n_r=1,n_r_max
               b(l1m1,n_r)=b(l1m1,n_r) +                   &
                           b_pol/10.0_cp*(r(n_r)/r_icb)**2 * &
                           ( one - three/5.0_cp*(r(n_r)/r_cmb)**2 )
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  b_ic(l1m1,n_r)=b_ic(l1m1,n_r)+b_pol/5.0_cp * &
                                 ( one - three/5.0_cp*(r_ic(n_r)/r_cmb)**2 )
               end do
            end if
         end if

      else if ( init_b1 == 21 ) then ! toroidal field created by inner core rotation
         ! equatorialy symmetric
         if ( llm <= l1m0 .and. ulm >= l1m0 ) then ! select processor
            do n_r=1,n_r_max
               aj0(n_r)=amp_b1*(r_icb*or1(n_r))**6
            end do
            do n_r=1,n_r_max             ! Diffusive toroidal field
               aj(l1m0,n_r)=aj(l1m0,n_r)+aj0(n_r)
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  aj_ic(l1m0,n_r)=aj_ic(l1m0,n_r)+aj0(n_r_icb)
               end do
            end if
         end if

      else if ( init_b1 == 22 ) then ! toroidal field created by inner core rotation
         ! equatorialy asymmetric
         if ( llm <= l2m0 .and. ulm >= l2m0 ) then ! select processor
            do n_r=1,n_r_max
               aj0(n_r)=amp_b1*(r_icb*or1(n_r))**6
            end do
            do n_r=1,n_r_max             ! Diffusive toroidal field
               aj(l2m0,n_r)=aj(l2m0,n_r)+aj0(n_r)
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  aj_ic(l2m0,n_r)=aj_ic(l2m0,n_r)+aj0(n_r_icb)
               end do
            end if
         end if

      else if ( init_b1 == 29 ) then ! toroidal field created by inner core rotation
         ! equatorialy asymmetric
         if ( llm <= l2m0 .and. ulm >= l2m0 ) then ! select processor
            do n_r=1,n_r_max
               aj0(n_r)=amp_b1*(r_icb*or1(n_r))
            end do
            do n_r=1,n_r_max             ! Diffusive toroidal field
               aj(l2m0,n_r)=aj(l2m0,n_r)+aj0(n_r)
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  aj_ic(l2m0,n_r)=aj_ic(l2m0,n_r)+aj0(n_r_icb)
               end do
            end if
         end if


      else if ( init_b1 == 30 ) then ! small scale random field
         ! random l,m mixture with radial structure that is a random mixture of (almost) Fourrier modes
            WRITE(OUTPUT_UNIT,*) 'enter init_b1..', rank
            !-- Random noise initialization of all (l,m) modes with wavelength between init_b_length_min and init_b_length_max
            l_smooth = 0.2_cp
            e_mag_p_r=0.0_cp
            !write(output_unit,*) 'init_length_min/max:', init_b_length_min, init_b_length_max!, '  l_smooth:',l_smooth
            if ( index(interior_model,'PNS_SZ_0V2S') /= 0 ) then
               r_cut=0.25_cp
            else if (index(interior_model,'HMNS_17MS_SZ') /= 0 ) then
               r_cut=0.2_cp
            else
               r_cut=0.40_cp
            end if
            do n_r_fourier=1,n_r_max
               if (one/n_r_fourier >= init_b_length_min .and. one/n_r_fourier <=init_b_length_max) then
                  write(output_unit,*) 'n_r_fourier:', n_r_fourier
                  do n_r=1,n_r_max
                     if (r(n_r) <= r_cut*r_cmb-half*l_smooth) then
                        b1(n_r)=0.0_cp
                        db1(n_r)=0.0_cp
                     else if (r(n_r) >=r_cut*r_cmb+half*l_smooth) then
                        b1(n_r)=cos(two*pi*n_r_fourier*(r(n_r)-r_icb))!/sqrt(rho0(n_r))!*test_hydro_test_BIS
                        db1(n_r)=-two*pi*n_r_fourier*sin(two*pi*n_r_fourier*(r(n_r)-r_icb))!/sqrt(rho0(n_r)) !- beta(n_r)*b1(n_r)/2
                     else
                        x_smooth=(r(n_r)-r_cut*r_cmb+half*l_smooth)/l_smooth
                        b1(n_r)=x_smooth*(three*x_smooth-two*x_smooth**2)*cos(two*pi*n_r_fourier*(r(n_r)-r_icb))!/sqrt(rho0(n_r))
                        db1(n_r)=(6.0_cp*x_smooth*(1.0_cp-x_smooth)/l_smooth*cos(two*pi*n_r_fourier*(r(n_r)-r_icb)) &
                        & - x_smooth*(three*x_smooth-two*x_smooth**2)*two*pi*n_r_fourier*sin(two*pi*n_r_fourier*(r(n_r)-r_icb)))
                        !/sqrt(rho0(n_r)) - beta(n_r)*b1(n_r)/2
                     end if
                     !write(output_unit,*) 'r:', r(n_r), '  b1:', b1(n_r), ' db1:', db1(n_r)
                  end do

                  do lm=llmMag,ulmMag
                     l1=lo_map%lm2l(lm)
                     m1=lo_map%lm2m(lm)
                     if (2*pi*(r_cmb+r_icb)/sqrt(l1*(l1+one)) >= init_b_length_min .and. &
                          &    2*pi*(r_cmb+r_icb)/sqrt(l1*(l1+one)) <=init_b_length_max) then
                        !write(output_unit,*) 'l:',l1,' m:',m1
                        l_corr=(l1*(l1+one)/(2*pi*(r_cmb+r_icb))**2+n_r_fourier**2)**(-0.5_cp)
                        !WRITE(OUTPUT_UNIT,*) 'l_corr:',l_corr
                        call random_number(bR)
                        call random_number(bI)
                        bR=(-one+two*bR)*l_corr**(3.0_cp+0.5_cp*init_b_index) !/D_l(st_map%lm2(l1,m1))**(bExp-1)
                        bI=(-one+two*bI)*l_corr**(3.0_cp+0.5_cp*init_b_index) !/D_l(st_map%lm2(l1,m1))**(bExp-1)
                        if ( m1 == 0 ) bI=0.0_cp
                        do n_r=1,n_r_max
                           b(lm,n_r)=b(lm,n_r) + cmplx(bR*b1(n_r),bI*b1(n_r),kind=cp)
                           e_mag_p_r(n_r)=e_mag_p_r(n_r)+dLh(st_map%lm2(l1,m1))*( dLh(st_map%lm2(l1,m1))* &
                                &     or2(n_r)*cc2real( cmplx(bR*b1(n_r),bI*b1(n_r),kind=cp),m1)  + &
                                &     cc2real(cmplx(bR*db1(n_r),bI*db1(n_r),kind=cp),m1) )
                           !WRITE(OUTPUT_UNIT,*) 'r:', r(n_r), '  b:', b(lm,n_r), ' db:', cmplx(bR*db1(n_r),bI*db1(n_r),kind=cp), ' e_mag_p_r:', e_mag_p_r(n_r)
                        end do
                     end if
                  end do
               end if
            end do

            ! renormalise the magnetic field such that the magnetic energy in the perturbed region
            ! equals that of a uniform field of amplitude amp_b1,
            ! where amp_b1 is given in units of sqrt(4*pi*rho)*d*omega
            e_mag_p=half*LFfac*eScale*rInt_R(e_mag_p_r,r,rscheme_oc) ! magnetic energy integrated over the volume
            WRITE(OUTPUT_UNIT,*) 'e_mag_p:',e_mag_p
            e_mag_p=e_mag_p/(four*third*pi*(r_cmb**3 - r_icb**3))/lScale**3 ! magnetic energy per unit volume
#ifdef WITH_MPI
            call MPI_Allreduce(e_mag_p,global_e_mag_p,1,MPI_DEF_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

            do lm=llmMag,ulmMag
               b(lm,:)=b(lm,:)*sqrt(rho0(:))*oek/sqrt(global_e_mag_p)*amp_b1
            end do
          !  b=b*oek/sqrt(global_e_mag_p)*amp_b1
            WRITE(OUTPUT_UNIT,*) 'e_mag_p/Volume:',global_e_mag_p, ' Volume:',four*third*pi*(r_cmb**3 - r_icb**3)/lScale**3
            WRITE(OUTPUT_UNIT,*) 'LFfac:',LFfac, ' oek:',oek, 'lScale:', lScale, ' eScale:', eScale

      else if ( init_b1 == 31 ) then ! small scale random field
         ! random l,m mixture with radial structure that is a random mixture of (almost) Fourrier modes respecting ktopb = 1
            WRITE(OUTPUT_UNIT,*) 'enter init_b1..', rank
            !-- Random noise initialization of all (l,m) modes with wavelength between init_b_length_min and init_b_length_max
            l_smooth = 0.2_cp
            e_mag_p_r=0.0_cp
            !write(output_unit,*) 'init_length_min/max:', init_b_length_min, init_b_length_max!, '  l_smooth:',l_smooth

            !rr=random(rand_num*rank/n_procs+one)
            do n_r_fourier=1,n_r_max
               if (one/n_r_fourier >= init_b_length_min .and. one/n_r_fourier <=init_b_length_max) then
                  write(output_unit,*) 'n_r_fourier:', n_r_fourier
                   do lm=llmMag,ulmMag
                     l1=lo_map%lm2l(lm)
                     m1=lo_map%lm2m(lm)
                     if (2*pi*(r_cmb+r_icb)/sqrt(l1*(l1+one)) >= init_b_length_min .and. &
                       &    2*pi*(r_cmb+r_icb)/sqrt(l1*(l1+one)) <=init_b_length_max) then
                        !write(output_unit,*) 'l:',l1,' m:',m1
                        l_corr=(l1*(l1+one)/(two*pi*(r_cmb+r_icb))**2+n_r_fourier**2)**(-0.5_cp)
                        !WRITE(OUTPUT_UNIT,*) 'l_corr:',l_corr
                        call random_number(bR)
                        call random_number(bI)
                        bR=(-one+two*bR)*l_corr**(3.0_cp+0.5_cp*init_b_index) !/D_l(st_map%lm2(l1,m1))**(bExp-1)
                        bI=(-one+two*bI)*l_corr**(3.0_cp+0.5_cp*init_b_index) !/D_l(st_map%lm2(l1,m1))**(bExp-1)
                        if ( m1 == 0 ) bI=0.0_cp
                        do n_r=1,n_r_max
                           if (r(n_r) <= 0.4_cp*r_cmb-half*l_smooth) then
                              b1(n_r)=0.0_cp
                              db1(n_r)=0.0_cp
                           else if (r(n_r) >=0.4_cp*r_cmb+half*l_smooth) then
                              b1(n_r)=cos(two*pi*n_r_fourier*(r(n_r)-r_icb)) - cos(two*pi*n_r_fourier*(r_cmb-r_icb))
                              db1(n_r)=-two*pi*n_r_fourier*sin(two*pi*n_r_fourier*(r(n_r)-r_icb))
                           else
                              x_smooth=(r(n_r)-0.4_cp*r_cmb+half*l_smooth)/l_smooth
                              b1(n_r)=x_smooth*(three*x_smooth-two*x_smooth**2)*(cos(two*pi*n_r_fourier*(r(n_r)-r_icb)) &
                                   & -cos(two*pi*n_r_fourier*(r_cmb-r_icb)))
                              db1(n_r)=6.0_cp*x_smooth*(1.0_cp-x_smooth)/l_smooth*(cos(two*pi*n_r_fourier*(r(n_r)-r_icb))-1)&
                                   & - x_smooth*(three*x_smooth-two*x_smooth**2)*two*pi*n_r_fourier                     &
                                   & *sin(two*pi*n_r_fourier*(r(n_r)-r_icb))
                           end if
                           !write(output_unit,*) 'r:', r(n_r), '  b1:', b1(n_r), ' db1:', db1(n_r)
                           b(lm,n_r)=b(lm,n_r) + cmplx(bR*b1(n_r),bI*b1(n_r),kind=cp)
                           e_mag_p_r(n_r)=e_mag_p_r(n_r)+dLh(st_map%lm2(l1,m1))*( dLh(st_map%lm2(l1,m1))* &
                                  &     or2(n_r)*cc2real( cmplx(bR*b1(n_r),bI*b1(n_r),kind=cp),m1)  + &
                                  &     cc2real(cmplx(bR*db1(n_r),bI*db1(n_r),kind=cp),m1) )
                           !WRITE(OUTPUT_UNIT,*) 'r:', r(n_r), '  b:', b(lm,n_r), ' db:', cmplx(bR*db1(n_r),bI*db1(n_r),kind=cp), ' e_mag_p_r:', e_mag_p_r(n_r)
                        end do
                     end if
                  end do
               end if
            end do

            ! renormalise the magnetic field such that the magnetic energy in the perturbed region
            ! equals that of a uniform field of amplitude amp_b1,
            ! where amp_b1 is given in units of sqrt(4*pi*rho)*d*omega
            e_mag_p=half*LFfac*eScale*rInt_R(e_mag_p_r,r,rscheme_oc) ! magnetic energy integrated over the volume
            WRITE(OUTPUT_UNIT,*) 'e_mag_p:',e_mag_p
            e_mag_p=e_mag_p/(four*third*pi*r_cmb**3*(one - 0.4_cp**3))/lScale**3 ! magnetic energy per unit volume
#ifdef WITH_MPI
            call MPI_Allreduce(e_mag_p,global_e_mag_p,1,MPI_DEF_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif
            b=b*oek/sqrt(global_e_mag_p)*amp_b1
            WRITE(OUTPUT_UNIT,*) 'e_mag_p/Volume:',global_e_mag_p, ' Volume:',four*third*pi*r_cmb**3*(one - 0.4_cp**3)/lScale**3
            WRITE(OUTPUT_UNIT,*) 'LFfac:',LFfac, ' oek:',oek, ' lScale:', lScale, ' eScale:', eScale

      else if ( init_b1 == 32 ) then !  Axial dipolar field with current equal to J=9*amp_b1*r0**6*r**2/(r**3+r0**3)**3
         if ( llm <= l1m0 .and. ulm >= l1m0 ) then ! select processor
            b_pol=amp_b1*sqrt(four*three*pi)*(r_cmb/2)**3
            do n_r=1,n_r_max
               b(l1m0,n_r)=b(l1m0,n_r)+b_pol*((r(n_r))**2)/((r(n_r))**3+(r_cmb/2.5)**3)
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  b_ic(l1m0,n_r)=b_ic(l1m0,n_r) + b_pol*(r_ic(n_r)**2)/(r_ic(n_r)**3+(r_cmb/2.5)**3)
               end do
            end if
         end if

       else if ( init_b1 == 33 ) then !  Axial dipolar field
         if ( llm <= l1m0 .and. ulm >= l1m0 ) then ! select processor
            do n_r=1,n_r_max
               b(l1m0,n_r)=b(l1m0,n_r)+amp_b1*sin(pi*(r(n_r)-r_icb))*r(n_r)**2/3
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  b_ic(l1m0,n_r)=b_ic(l1m0,n_r) + sin(pi*(r(n_r)-r_icb))*amp_b1*r_ic(n_r)**2/3
               end do
            end if
         end if

      else if ( init_b1 == 34 ) then ! small scale random field for anelastic models
         ! random l,m mixture with radial structure that is a random mixture of (almost) Fourrier modes
            WRITE(OUTPUT_UNIT,*) 'enter init_b1..', rank
            !-- Random noise initialization of all (l,m) modes with wavelength between init_b_length_min and init_b_length_max
            !l_smooth = 0.1_cp
            e_mag_p_r=0.0_cp
            !write(output_unit,*) 'init_length_min/max:', init_b_length_min, init_b_length_max!, '  l_smooth:',l_smooth
            !r_cut=0.25_cp
            !rr=random(rand_num*rank/n_procs+one)
            write (*,*) "Random number",rand_num*rank/n_procs+one,"rank",rank
            do n_r_fourier=1,n_r_max
               if (one/n_r_fourier >= init_b_length_min .and. one/n_r_fourier <=init_b_length_max) then
                  do n_r=1,n_r_max
                        b1(n_r)=(1-cos(two*pi*n_r_fourier*(r(n_r)-r_icb)))!/sqrt(rho0(n_r))
                        db1(n_r)=two*pi*n_r_fourier*sin(two*pi*n_r_fourier*(r(n_r)-r_icb))!&
                  end do
                  do lm=llmMag,ulmMag
                     l1=lo_map%lm2l(lm)
                     m1=lo_map%lm2m(lm)
                     if (2*pi*(r_cmb+r_icb)/sqrt(l1*(l1+one)) >= init_b_length_min .and. &
                          &    2*pi*(r_cmb+r_icb)/sqrt(l1*(l1+one)) <=init_b_length_max) then
                        !write(output_unit,*) 'l:',l1,' m:',m1
                        l_corr=(l1*(l1+one)/(2*pi*(r_cmb+r_icb))**2+n_r_fourier**2)**(-0.5_cp)
                        !WRITE(OUTPUT_UNIT,*) 'l_corr:',l_corr
                        call random_number(bR)
                        call random_number(bI)
                        bR=(-one+two*bR)*l_corr**(3.0_cp+0.5_cp*init_b_index) !/D_l(st_map%lm2(l1,m1))**(bExp-1)
                        bI=(-one+two*bI)*l_corr**(3.0_cp+0.5_cp*init_b_index) !/D_l(st_map%lm2(l1,m1))**(bExp-1)
                        if ( m1 == 0 ) bI=0.0_cp
                        do n_r=1,n_r_max
                           b(lm,n_r)=b(lm,n_r) + cmplx(bR*b1(n_r),bI*b1(n_r),kind=cp)
                           e_mag_p_r(n_r)=e_mag_p_r(n_r)+dLh(st_map%lm2(l1,m1))*( dLh(st_map%lm2(l1,m1))* &
                                &     or2(n_r)*cc2real( cmplx(bR*b1(n_r),bI*b1(n_r),kind=cp),m1)  + &
                                &     cc2real(cmplx(bR*db1(n_r),bI*db1(n_r),kind=cp),m1) )
                        !  WRITE(OUTPUT_UNIT,*) 'r:', r(n_r), '  b:', b(lm,n_r), ' db:', &
                        !        & cmplx(bR*db1(n_r),bI*db1(n_r),kind=cp), ' e_mag_p_r:', e_mag_p_r(n_r)
                        end do
                     end if
                  end do
               end if
            end do
!         end if
!            end do

            ! renormalise the magnetic field such that the magnetic energy in the perturbed region
            ! equals that of a uniform field of amplitude amp_b1,
            ! where amp_b1 is given in units of sqrt(4*pi*rho)*d*omega
            e_mag_p=half*LFfac*eScale*rInt_R(e_mag_p_r,r,rscheme_oc) ! magnetic energy integrated over the volume
            WRITE(OUTPUT_UNIT,*) 'e_mag_p:',e_mag_p
            e_mag_p=e_mag_p/(four*third*pi*(r_cmb**3 - r_icb**3))/lScale**3 ! magnetic energy per unit volume
#ifdef WITH_MPI
            call MPI_Allreduce(e_mag_p,global_e_mag_p,1,MPI_DEF_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

            do lm=llmMag,ulmMag
               b(lm,:)=b(lm,:)*sqrt(rho0(:))*oek/sqrt(global_e_mag_p)*amp_b1
            end do
          !  b=b*oek/sqrt(global_e_mag_p)*amp_b1
            WRITE(OUTPUT_UNIT,*) 'e_mag_p/Volume:',global_e_mag_p, ' Volume:',four*third*pi*(r_cmb**3 - r_icb**3)/lScale**3
            WRITE(OUTPUT_UNIT,*) 'LFfac:',LFfac, ' oek:',oek, 'lScale:', lScale, ' eScale:', eScale

       else if (init_b1 == 35) then
          if ( llm <= l1m0 .and. ulm >= l1m0 ) then ! select processor
            b_pol=amp_b1
            x2 = 0.35*r_cmb-0.3*r_cmb
            alpha = 1/(x2**3)
            beta1 = -3/(x2**2)
            c = 3/(x2)
            do n_r=1,n_r_max
               if(r(n_r) >= 0.35*r_cmb) then
                  b(l1m0,n_r)=b(l1m0,n_r)+b_pol*((r(n_r)-0.3*r_cmb))**2
               else if (r(n_r) >= 0.3*r_cmb) then
                  b(l1m0,n_r)=b(l1m0,n_r) + b_pol*((0.35*r_cmb)**2)*(alpha*((r(n_r)-0.3*r_cmb)**5)+ &
                             & beta1*((r(n_r)-0.3*r_cmb)**4)+c*((r(n_r)-0.3*r_cmb)**3))
               else
                  b(l1m0,n_r)= zero
               end if
            end do
         end if

      else if ( init_b1 == 36 ) then  ! l=1,m=0 poloidal field , constant in r !
         ! no potential at r_cmb but simple
         if ( llm <= l1m0 .and. ulm >= l1m0 ) then ! select processor
            b_pol=amp_b1
            do n_r=1,n_r_max
               b(l1m0,n_r)=b(l1m0,n_r)+b_pol*r(n_r)**2
               e_mag_p_r(n_r)=e_mag_p_r(n_r)+dLh(st_map%lm2(1,0))*( dLh(st_map%lm2(1,0))* &
                                  &     or2(n_r)*(b(l1m0,n_r)**2)+ &
                                  &     (2*r(n_r)*b_pol)**2)
            end do
            e_mag_p=half*LFfac*eScale*rInt_R(e_mag_p_r,r,rscheme_oc) ! magnetic energy integrated over the volume
            e_mag_p=e_mag_p/(four*third*pi*r_cmb**3*(one - 0.4_cp**3))/lScale**3 ! magnetic energy per unit volume
            b=b*oek*sqrt(rho0(n_r_max))/sqrt(e_mag_p)*amp_b1
            if (n_imp>=2) then
               amp_imp = oek*sqrt(rho0(n_r_max))/sqrt(e_mag_p)*amp_b1*y10_norm*b_pol*r_cmb
               WRITE(OUTPUT_UNIT,*) 'amp_imp =',amp_imp ! For restart
            end if
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  b_ic(l1m0,n_r)=b_ic(l1m0,n_r)+b_pol*oek/sqrt(e_mag_p)*amp_b1*r_icb**2
               end do
            end if
         end if
      else if ( init_b1 == 37 ) then ! small scale random toroidal field for anelastic models
         ! random l,m mixture with radial structure that is a random mixture of (almost) Fourrier modes
         WRITE(OUTPUT_UNIT,*) 'enter init_b1..', rank
         !-- Random noise initialization of all (l,m) modes with wavelength between init_b_length_min and init_b_length_max
         l_smooth = 0.1_cp
         e_mag_p_r=0.0_cp
         !write(output_unit,*) 'init_length_min/max:', init_b_length_min, init_b_length_max!, '  l_smooth:',l_smooth
         r_cut=0.25_cp
         !rr=random(rand_num*rank/n_procs+one)
         write (*,*) "Random number",rand_num*rank/n_procs+one,"rank",rank
         do n_r_fourier=1,n_r_max
            if (one/n_r_fourier >= init_b_length_min .and. one/n_r_fourier <=init_b_length_max) then
               !                  write(output_unit,*) 'n_r_fourier:', n_r_fourier
               do n_r=1,n_r_max
                     b1(n_r)=(1-cos(two*pi*n_r_fourier*(r(n_r)-r_icb)))!/sqrt(rho0(n_r))
                     db1(n_r)=two*pi*n_r_fourier*sin(two*pi*n_r_fourier*(r(n_r)-r_icb))!&
                     !& /sqrt(rho0(n_r))- beta(n_r)*b1(n_r)/2
               end do
               do lm=llmMag,ulmMag
                  l1=lo_map%lm2l(lm)
                  m1=lo_map%lm2m(lm)
                  if (2*pi*(r_cmb+r_icb)/sqrt(l1*(l1+one)) >= init_b_length_min .and. &
                       &    2*pi*(r_cmb+r_icb)/sqrt(l1*(l1+one)) <=init_b_length_max) then
                     !write(output_unit,*) 'l:',l1,' m:',m1
                     l_corr=(l1*(l1+one)/(2*pi*(r_cmb+r_icb))**2+n_r_fourier**2)**(-0.5_cp)
                     !WRITE(OUTPUT_UNIT,*) 'l_corr:',l_corr
                     call random_number(bR)
                     call random_number(bI)
                     bR=(-one+two*bR)*l_corr**(3.0_cp+0.5_cp*init_b_index) !/D_l(st_map%lm2(l1,m1))**(bExp-1)
                     bI=(-one+two*bI)*l_corr**(3.0_cp+0.5_cp*init_b_index) !/D_l(st_map%lm2(l1,m1))**(bExp-1)
                     if ( m1 == 0 ) bI=0.0_cp
                     do n_r=1,n_r_max
                        aj(lm,n_r)=aj(lm,n_r) + cmplx(bR*b1(n_r),bI*b1(n_r),kind=cp)
                        e_mag_p_r(n_r)=e_mag_p_r(n_r)+dLh(st_map%lm2(l1,m1)) * cc2real(aj(lm,n_r),m1)
                        !WRITE(OUTPUT_UNIT,*) 'r:', r(n_r), '  b:', b(lm,n_r), ' db:', cmplx(bR*db1(n_r),bI*db1(n_r),kind=cp), ' e_mag_p_r:', e_mag_p_r(n_r)
                     end do
                  end if
               end do
            end if
         end do

         ! renormalise the magnetic field such that the magnetic energy in the perturbed region
         ! equals that of a uniform field of amplitude amp_b1,
         ! where amp_b1 is given in units of sqrt(4*pi*rho)*d*omega
         e_mag_p=half*LFfac*eScale*rInt_R(e_mag_p_r,r,rscheme_oc) ! magnetic energy integrated over the volume
         WRITE(OUTPUT_UNIT,*) 'e_mag_p:',e_mag_p
         e_mag_p=e_mag_p/(four*third*pi*(r_cmb**3 - r_icb**3))/lScale**3 ! magnetic energy per unit volume
#ifdef WITH_MPI
         call MPI_Allreduce(e_mag_p,global_e_mag_p,1,MPI_DEF_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

         do lm=llmMag,ulmMag
            aj(lm,:)=aj(lm,:)*sqrt(rho0(:))*oek/sqrt(global_e_mag_p)*amp_b1
         end do
         !  b=b*oek/sqrt(global_e_mag_p)*amp_b1
         WRITE(OUTPUT_UNIT,*) 'e_mag_p/Volume:',global_e_mag_p, ' Volume:',four*third*pi*(r_cmb**3 - r_icb**3)/lScale**3
         WRITE(OUTPUT_UNIT,*) 'LFfac:',LFfac, ' oek:',oek, 'lScale:', lScale, ' eScale:', eScale

      else if ( init_b1 == 38 ) then

         if ( llm <= l2m0 .and. ulm >= l2m0 ) then ! select processor
            b_pol=-four*third*amp_b1*sqrt(pi/5.0_cp)
            if ( l_cond_ic ) then
!               b_pol=amp_b1*sqrt(three*pi)/(three+r_cmb)
               do n_r=1,n_r_max
                  aj(l2m0,n_r)=aj(l2m0,n_r) + b_pol*(r(n_r)-r_icb)**2*(r(n_r)-r_cmb)**2/r(n_r)!*sin(pi*(r(n_r)/r_cmb))
               end do
               arg=pi*r_icb/r_cmb
               aj_ic1=(arg-two*sin(arg)*cos(arg)) / (arg+sin(arg)*cos(arg))
               aj_ic2=(one-aj_ic1)*r_icb*sin(arg)/cos(arg)
               do n_r=1,n_r_ic_max
                  aj_ic(l2m0,n_r)=aj_ic(l2m0,n_r)+b_pol*             ( &
                  &        aj_ic1*r_ic(n_r)*sin(pi*r_ic(n_r)/r_cmb) +  &
                  &                  aj_ic2*cos(pi*r_ic(n_r)/r_cmb) )
               end do
            else
!               b_tor=amp_b1*sqrt(three*pi)/four
               do n_r=1,n_r_max
                  aj(l2m0,n_r)=aj(l2m0,n_r) + b_pol*(r(n_r)-r_icb)**2*(r(n_r)-r_cmb)**2/r(n_r)!*sin(pi*(r(n_r)-r_icb))
                  e_mag_p_r(n_r)=e_mag_p_r(n_r)+dLh(st_map%lm2(2,0)) * cc2real(aj(l2m0,n_r),0)
               end do
            end if
            e_mag_p=half*LFfac*eScale*rInt_R(e_mag_p_r,r,rscheme_oc)
            e_mag_p=e_mag_p/(four*third*pi*(r_cmb**3 - r_icb**3))/lScale**3 ! magnetic energy per unit volume
            if (l_non_rot) then
               aj(l2m0,:)=aj(l2m0,:)*sqrt(rho0(:))/sqrt(e_mag_p)*amp_b1
            else
               aj(l2m0,:)=aj(l2m0,:)*sqrt(rho0(:))*oek/sqrt(e_mag_p)*amp_b1
            end if
            WRITE(OUTPUT_UNIT,*) 'e_mag_p/Volume:',e_mag_p, ' Volume:',four*third*pi*(r_cmb**3 - r_icb**3)/lScale**3
            WRITE(OUTPUT_UNIT,*) 'LFfac:',LFfac, ' oek:',oek, 'lScale:', lScale, ' eScale:', eScale
         end if
      else if ( init_b1 == 39 ) then  ! l=9,m0 poloidal field
         ! which is potential field at r_cmb
         l3m0 = lo_map%lm2(8,0)
         if ( llm <= l3m0 .and. ulm >= l3m0 ) then ! select processor
            b_pol=amp_b1*13.0_cp/6.0_cp*sqrt(pi/19.0_cp)*r_icb**2*radratio !no idea why 7/6
            do n_r=1,n_r_max
               b(l3m0,n_r)=b(l3m0,n_r)+b_pol*(r(n_r)/r_icb)**3 * &
                           ( one - 11.0_cp/13.0_cp*(r(n_r)/r_cmb)**2 )
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                   b_ic(l3m0,n_r)=b_ic(l3m0,n_r)+b_pol * &
                                 ( one - 11.0_cp/13.0_cp*(r_ic(n_r)/r_cmb)**2 )
               end do
            end if
         end if
      end if

!      !-- Too lazy to calculate these:
!      lorentz_torque_ic=0.0_cp
!      lorentz_torque_ma=0.0_cp

   end subroutine initB
!-----------------------------------------------------------------------
   subroutine initPhi(s, phi)
      !
      ! This subroutine sets the initial phase field distribution. It follows
      ! a tanh function with a width of size epsPhase
      !

      !-- Input variable
      complex(cp), intent(inout) :: s(llm:ulm,n_r_max) ! Entropy/Temperature

      !-- In/out variable
      complex(cp), intent(inout) :: phi(llm:ulm,n_r_max) ! Phase field

      !-- Local variables:
      real(cp) :: temp00(n_r_max), rmelt, phi0(n_r_max), rphase, x(4), y(4)
      integer :: lm00, n_r, n_r_melt, l, lm, n_p, n_t, nPstart, nPstop
      integer :: n_r_phase, n_r_start, n_r_stop
      complex(cp), allocatable :: s_Rloc(:,:), phi_Rloc(:,:)
      real(cp), allocatable :: temp_Rloc(:,:,:), temp_Ploc(:,:,:)
      real(cp), allocatable :: phase_Rloc(:,:,:), phase_Ploc(:,:,:)
      type(load), allocatable :: phi_balance(:)

      if ( init_phi == 1 ) then
         !-- The initial phase field is set as a tanh function of width epsPhase
         !-- centered at the melting temperature
         lm00 = lo_map%lm2(0,0) ! spherically-symmetric mode
         if ( llm <= lm00 .and. ulm >= lm00 ) then
            temp00(:)=osq4pi * real(s(lm00,:))
            do n_r=2,n_r_max
               if ( temp00(n_r-1) < tmelt .and. temp00(n_r) >= tmelt ) then
                  n_r_melt=n_r
               end if
            end do
            rmelt=r(n_r_melt)
            phi0(:)=half*(one+tanh((r(:)-rmelt)/two/sqrt(two)/epsPhase))
            phi(lm00,:)=sq4pi*cmplx(phi0,0.0_cp,cp)
         end if

#ifdef WITH_MPI
         call MPI_Bcast(phi0,n_r_max,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
#endif

         !-- Make sure there is no temperature perturbation in the solid phase
         if ( init_s1 /= 0 .or. init_s2 /= 0 ) then
            do n_r=1,n_r_max
               do lm=llm,ulm
                  l = lo_map%lm2l(lm)
                  if ( l == 0 ) cycle
                  s(lm,n_r)=s(lm,n_r)*(one-phi0(n_r))
               end do
            end do
         end if

      else if ( init_phi == 2 ) then
         !-- Compute locally the temperature distribution in the physical space
         !-- and then get the corresponding phase field
         allocate(s_Rloc(lm_max,nRstart:nRstop), phi_Rloc(lm_max,nRstart:nRstop))
         allocate( phi_balance(0:n_procs-1) )
         call getBlocks(phi_balance, n_phi_max, n_procs)
         nPstart = phi_balance(rank)%nStart
         nPstop = phi_balance(rank)%nStop
         allocate( phase_Rloc(nlat_padded,n_phi_max,nRstart:nRstop) )
         allocate( temp_Rloc(nlat_padded,n_phi_max,nRstart:nRstop) )
         allocate( phase_Ploc(nlat_padded,nPstart:nPstop,n_r_max) )
         allocate( temp_Ploc(nlat_padded,nPstart:nPstop,n_r_max) )

         !-- Bring the temperature field on the R-distributed version
         call lo2r_one%transp_lm2r(s, s_Rloc)

         !-- Compute the temperature on the physical grid
         do n_r=nRstart,nRstop
            call scal_to_spat(s_Rloc(:,n_r), temp_Rloc(:,:,n_r), l_max)
         end do

         !-- Transpose the temperature on a phi-distributed grid
         call transp_R2Phi(temp_Rloc, temp_Ploc, phi_balance, nPstart, nPstop)

         !-- Determine the melt radius by 4th order Lagrangian smoothing
         do n_p=nPstart,nPstop
            do n_t=1,n_theta_max
               !-- Determine the radial level where T=Tmelt
               n_r_phase=1
               do n_r=2,n_r_max
                  if ( temp_Ploc(n_t,n_p,n_r) > tmelt .and. &
                  &    temp_Ploc(n_t,n_p,n_r-1) < tmelt ) then
                     n_r_phase=n_r
                  end if
               end do

               !-- 4th order Lagrangian interpolation of melting point
               if ( n_r_phase == n_r_cmb ) then
                  rphase=r_cmb
               else
                  if ( n_r_phase == n_r_cmb+1 ) then
                     n_r_start=n_r_phase-1
                     n_r_stop =n_r_phase+2
                  else if ( n_r_phase == n_r_icb ) then
                     n_r_start=n_r_phase-3
                     n_r_stop =n_r_phase
                  else
                     n_r_start=n_r_phase-2
                     n_r_stop =n_r_phase+1
                  end if
                  x(:)=temp_Ploc(n_t,n_p,n_r_start:n_r_stop)
                  y(:)=r(n_r_start:n_r_stop)
                  rphase=lagrange_interp(x,tmelt,y)
               end if

               phase_Ploc(n_t,n_p,:)=half*(one+tanh((r(:)-rphase)/two/epsPhase))

            end do
         end do

         !-- Transpose phase from phi-distributed to r-distributed
         call transp_Phi2R(phase_Ploc, phase_Rloc, phi_balance, nPstart, nPstop)

         !-- Now bring the phase back to spectral space
         do n_r=nRstart,nRstop
            call scal_to_SH(phase_Rloc(:,:,n_r), phi_Rloc(:,n_r), l_max)
         end do

         !-- Final transpose to get the LM-distributed array for the phase field
         call r2lo_one%transp_r2lm(phi_Rloc, phi)

         deallocate( temp_Rloc, temp_Ploc, phase_Rloc, phase_Ploc, phi_balance )
         deallocate( s_Rloc, phi_Rloc )
      end if

   end subroutine initPhi
!-----------------------------------------------------------------------
   subroutine initF(bodyForce)
      !
      ! This subroutine is used to initialize a toroidal body force
      ! of the form (-a s + b s^2) \vec{e}_\phi in lm space. This can easily
      ! be extended to other forms of body forces prescribed in physical space.
      !

      !-- In/out variables
      complex(cp), intent(inout) :: bodyForce(llm:ulm,n_r_max)

      !-- Local variables
      complex(cp) :: bf_Rloc(lm_max,nRstart:nRstop), bfLM(lm_max)
      integer :: lm,l,m,nR,nTheta,nPhi
      real(cp) :: bf_spat(nlat_padded,n_phi_max)
      real(cp) :: eta_fac, a_force, b_force

      call lo2r_one%transp_lm2r(bodyForce, bf_Rloc)

      eta_fac = 15.0_cp*pi/64.0_cp * (1.0_cp - radratio**6)/(1.0_cp-radratio**5)
      a_force = eta_fac / (r_cmb * (1.0_cp - eta_fac))
      b_force = 1.0_cp / (r_cmb**2 * ( 1.0_cp - eta_fac ))

      do nR = nRstart,nRstop
         do nPhi=1,n_phi_max
            do nTheta=1,n_theta_max
               ! This is F_{\phi}/sin(\theta) for toroidal equation
               ! where F = ( -a s + b s^2 ) \vec{e}_\phi
               ! Use F_r for poloidal equation
               bf_spat(nTheta,nPhi) = ampForce   *               &
               &                      (- a_force * r(nR)         &
               &                       + b_force * r(nR) * r(nR) &
               &                       * sinTheta(nTheta) )
            end do
         end do

         call scal_to_SH(bf_spat, bfLM, l_max)

         !------- body force is now in spherical harmonic space,
         !        For toroidal equation, get radial component of
         !        curl by applying operator
         !        dTheta1=1/(r sinTheta) d/ d theta sinTheta**2,
         !        comment out for poloidal equation
         do lm=2,lm_max
            l=st_map%lm2l(lm)
            m=st_map%lm2m(lm)
            if ( l < l_max .and. l > m ) then
               bf_Rloc(lm,nR)=dTheta1S(lm)*bfLM(st_map%lm2lmS(lm))   &
               &             -dTheta1A(lm)*bfLM(st_map%lm2lmA(lm))
            else if ( l < l_max .and. l == m ) then
               bf_Rloc(lm,nR)=dTheta1A(lm)*bfLM(st_map%lm2lmA(lm))
            else if ( l == l_max .and. m < l ) then
               bf_Rloc(lm,nR)=dTheta1S(lm)*bfLM(st_map%lm2lmS(lm))
            end if
         end do

      end do ! close loop over radial points

      call r2lo_one%transp_r2lm(bf_Rloc, bodyForce)

   end subroutine initF
!-----------------------------------------------------------------------
   subroutine initTidal(we,dwe,ddwe,wer,dwer,ddwer)
      !
      ! This subroutine is used to initialize a quasi-equilibrium tides u_e
      !
      !-- In/out variables
     complex(cp), intent(inout) :: we(llm:ulm,n_r_max)
     complex(cp), intent(inout) :: dwe(llm:ulm,n_r_max)
     complex(cp), intent(inout) :: ddwe(llm:ulm,n_r_max)
     complex(cp), intent(inout) :: wer(lm_max,nRstart:nRstop)
     complex(cp), intent(inout) :: dwer(lm_max,nRstart:nRstop)
     complex(cp), intent(inout) :: ddwer(lm_max,nRstart:nRstop)

     !-- Local variables
     complex(cp) :: work1(lm_max)
     integer :: lm,l,m,st_lmP, nR, l2m2
     !class(type_mpitransp), pointer :: r2lo_initt, lo2r_initt
     real(cp) ::  r_rcmb, alpha

      ! allocate( type_mpiptop :: r2lo_initt )
      ! allocate( type_mpiptop :: lo2r_initt )

      ! call r2lo_initt%create_comm(1)
      ! call lo2r_initt%create_comm(1)

      l2m2 = lo_map%lm2(2,2) !2

      X(:)=(r(:)**2+(2.0D0/3.0D0)*((radratio**5)*or3(:)))&
           & /(2.0D0*(1.0D0-radratio**5))
      dX(:)=2.0D0*(r(:)-((radratio**5)*or2(:)*or2(:)))&
           & /(2.0D0*(1.0D0-radratio**5))
      d2X(:)=sqrt(15.0D0/32.0D0/pi)*tidalFac*dX*or1(:)*rho0(:)
      !alpha=r(n_r_max)/r(1)
      !d2X(:)=2.0D0*(1.0D0+4.0D0*((radratio**5)*or2(:)*or3(:)))&
      !     & /(2.0D0*(1.0D0-radratio**5))
      if (l2m2 >= llm .and. l2m2 <= ulm) then
         l=2
         m=2
         do nR = 1,n_r_max
            we(l2m2,nR)= - cmplx(0.0,1.0,kind=cp)*r(nR)**2/(l*(l+1))*w_orbit_th &
                 & *tidalFac*dX(nR)*rho0(nR)/two !last divide by two is due that force is only the real part
         end do
         write(*,*) aimag(-we(l2m2,:)*l*(l+1)*or2(:))
      end if

      call get_dr(we,dwe,ulm-llm+1,1,ulm-llm+1,n_r_max,rscheme_oc, &
           &      nocopy=.true.)
      call get_dr(dwe,ddwe,ulm-llm+1,1,ulm-llm+1,n_r_max,rscheme_oc, &
           &      nocopy=.true.)

      call r2lo_one%transp_lm2r(we,wer)
      call r2lo_one%transp_lm2r(dwe,dwer)
      call r2lo_one%transp_lm2r(ddwe,ddwer)

      work1(:)=zero
      do nR = nRstart,nRstop
         call torpol_to_spat(wer(:,nR), dwer(:,nR),  work1, &
              &              vrtidal(:,:,nR), vttidal(:,:,nR),vptidal(:,:,nR), l_R(nR))
      end do

      ! call r2lo_initt%destroy_comm()
      ! call lo2r_initt%destroy_comm()

   end subroutine initTidal
!-----------------------------------------------------------------------
   subroutine j_cond(lm0, aj0, aj0_ic)
      !
      ! Purpose of this subroutine is to solve the diffusion equation
      ! for an initial toroidal magnetic field.
      !

      !-- Input variable:
      integer, intent(in) :: lm0

      !-- Output variables:
      complex(cp), intent(out) :: aj0(:)    ! aj(l=0,m=0) in the outer core
      complex(cp), intent(out) :: aj0_ic(:) ! aj(l=0,m=0) in the inner core

      !-- Local variables
      integer :: n_cheb,n_r,info,n_r_real,n_r_out, l
      real(cp) :: dL
      complex(cp) :: rhs(n_r_tot)
      complex(cp) :: work_l_ic(n_r_ic_max)
      real(cp), allocatable :: jMat(:,:)
      integer, allocatable :: jPivot(:)

      allocate( jMat(n_r_tot,n_r_tot) )
      allocate( jPivot(n_r_tot) )

      l = lo_map%lm2l(lm0)
      dL = real(l*(l+1),cp)

      n_r_real = n_r_max
      if ( l_cond_ic ) n_r_real = n_r_real+n_r_ic_max

      !----- Outer core:
      do n_r_out=1,rscheme_oc%n_max
         do n_r=2,n_r_max-1
            jMat(n_r,n_r_out)= rscheme_oc%rnorm *                   &
            &              hdif_B(l)*dL*opm*lambda(n_r)*or2(n_r) *  &
            &         (            rscheme_oc%d2rMat(n_r,n_r_out) + &
            &         dLlambda(n_r)*rscheme_oc%drMat(n_r,n_r_out) - &
            &            dL*or2(n_r)*rscheme_oc%rMat(n_r,n_r_out) )
         end do
      end do

      !----- boundary conditions:
      !----- CMB:
      do n_r_out=1,rscheme_oc%n_max     ! should be bpeaktop at CMB
         jMat(1,n_r_out)=rscheme_oc%rnorm*rscheme_oc%rMat(1,n_r_out)
      end do
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         jMat(1,n_r_out)=0.0_cp
      end do

      !----- ICB:
      if ( l_cond_ic ) then  ! matching condition at inner core:
         do n_r_out=1,rscheme_oc%n_max
            jMat(n_r_max,n_r_out)  =rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,n_r_out)
            jMat(n_r_max+1,n_r_out)=rscheme_oc%rnorm*sigma_ratio* &
            &                       rscheme_oc%drMat(n_r_max,n_r_out)
         end do
         do n_r_out=rscheme_oc%n_max+1,n_r_max
            jMat(n_r_max,n_r_out)  =0.0_cp
            jMat(n_r_max+1,n_r_out)=0.0_cp
         end do
      else
         do n_r_out=1,rscheme_oc%n_max
            jMat(n_r_max,n_r_out)=rscheme_oc%rMat(n_r_max,n_r_out)*rscheme_oc%rnorm
         end do
         do n_r_out=rscheme_oc%n_max+1,n_r_max
            jMat(n_r_max,n_r_out)=0.0_cp
         end do
      end if

      do n_r=1,n_r_max
         jMat(n_r,1)      =rscheme_oc%boundary_fac*jMat(n_r,1)
         jMat(n_r,n_r_max)=rscheme_oc%boundary_fac*jMat(n_r,n_r_max)
      end do

      !----- Inner core:
      if ( l_cond_ic ) then
         do n_cheb=1,n_r_ic_max ! counts even IC cheb modes
            do n_r=2,n_r_ic_max-1 ! counts IC radial grid point
               jMat(n_r_max+n_r,n_r_max+n_cheb) =                 &
               &        cheb_norm_ic*dL*or2(n_r_max)*opm*O_sr * ( &
               &        d2cheb_ic(n_cheb,n_r) + two*real(l+1,cp)* &
               &        O_r_ic(n_r)*dcheb_ic(n_cheb,n_r) )
            end do
            ! r=0: central point
            n_r=n_r_ic_max
            jMat(n_r_max+n_r,n_r_max+n_cheb) = cheb_norm_ic*dL*or2(n_r_max)* &
            &              opm*O_sr*(one+two*real(l+1,cp))*d2cheb_ic(n_cheb,n_r)
         end do

         !-------- boundary conditions:
         do n_cheb=1,n_cheb_ic_max
            jMat(n_r_max,n_r_max+n_cheb)=-cheb_norm_ic*cheb_ic(n_cheb,1)
            jMat(n_r_max+1,n_r_max+n_cheb)= -cheb_norm_ic * (      &
            &                                dcheb_ic(n_cheb,1)  + &
            &         real(l+1,cp)*or1(n_r_max)*cheb_ic(n_cheb,1) )
         end do
         do n_cheb=n_r_max+n_cheb_ic_max+1,n_r_tot
            jMat(n_r_max,n_cheb)  =0.0_cp
            jMat(n_r_max+1,n_cheb)=0.0_cp
         end do

         !-------- normalization for lowest Cheb mode:
         do n_r=n_r_max,n_r_tot
            jMat(n_r,n_r_max+1)=half*jMat(n_r,n_r_max+1)
            jMat(n_r,n_r_tot)  =half*jMat(n_r,n_r_tot)
         end do

         !-------- fill matrix up with zeros:
         do n_r_out=n_r_max+1,n_r_tot
            do n_r=1,n_r_max-1
               jMat(n_r,n_r_out)=0.0_cp
            end do
         end do
         do n_r_out=1,n_r_max
            do n_r=n_r_max+2,n_r_tot
               jMat(n_r,n_r_out)=0.0_cp
            end do
         end do

      end if ! conducting inner core ?

      !----- invert matrix:
      call prepare_mat(jMat,n_r_tot,n_r_real,jPivot,info)
      if ( info /= 0 ) then
         call abortRun('Singular matrix jMat in j_cond.')
      end if

      !----- zero RHS, except BC's
      do n_r=2,n_r_real-1
         rhs(n_r)=zero
      end do
      rhs(1)= bpeaktop                             ! Outer boundary
      if ( .not. l_cond_ic ) rhs(n_r_max)=bpeakbot  ! Inner boundary

      !----- solve linear system:
      call solve_mat(jMat,n_r_tot,n_r_real,jPivot,rhs)

      !----- copy result for OC:
      do n_r_out=1,rscheme_oc%n_max
         aj0(n_r_out)=rhs(n_r_out)
      end do

      do n_r_out=rscheme_oc%n_max+1,n_r_max
         aj0(n_r_out)=zero
      end do

      !----- transform to radial space:
      call rscheme_oc%costf1(aj0)

      if ( l_cond_ic ) then
         !----- copy result for IC:
         do n_cheb=1,n_cheb_ic_max
            aj0_ic(n_cheb)=rhs(n_r_max+n_cheb)
         end do
         do n_cheb=n_cheb_ic_max+1,n_r_ic_max
            aj0_ic(n_cheb)=zero
         end do

         !----- transform to radial space:
         !  Note: this is assuming that aj0_ic is an even function !
         call chebt_ic%costf1(aj0_ic,work_l_ic)
      end if

      deallocate( jMat, jPivot )

   end subroutine j_cond
!--------------------------------------------------------------------------------
   subroutine xi_cond(xi0)
      !
      ! Purpose of this subroutine is to solve the chemical composition equation
      ! for an the conductive (l=0,m=0)-mode.
      ! Output is the radial dependence of the solution in s0.
      !

      real(cp), intent(out) :: xi0(:) ! spherically-symmetric part

      !-- local variables:
      integer :: n_r_out, n_r, info, n_bands
      real(cp) :: rhs(n_r_max), dat(n_r_max,n_r_max)
      class(type_realmat), pointer :: xi0Mat

      if ( l_finite_diff ) then
         allocate( type_bandmat :: xi0Mat )
         if ( ktopxi == 1 .and. kbotxi == 1 .and. rscheme_oc%order <= 2 &
         &    .and. rscheme_oc%order_boundary <= 2 ) then
            n_bands = rscheme_oc%order+1
         else
            n_bands = max(2*rscheme_oc%order_boundary+1,rscheme_oc%order+1)
         end if
         call xi0Mat%initialize(n_bands,n_r_max,l_pivot=.true.)
      else
         allocate( type_densemat :: xi0Mat )
         call xi0Mat%initialize(n_r_max,n_r_max,l_pivot=.true.)
      end if

      !-- Set Matrix:
      do n_r_out=1,n_r_max
         do n_r=2,n_r_max-1
            dat(n_r,n_r_out)=rscheme_oc%rnorm*osc*(                       &
            &                            rscheme_oc%d2rMat(n_r,n_r_out) + &
            &  ( two*or1(n_r)+beta(n_r) )*                                &
            &                             rscheme_oc%drMat(n_r,n_r_out)  )
         end do
      end do

      !-- Set boundary conditions:
      if ( ktopxi == 1 .or. kbotxi == 2 ) then
         dat(1,:)=rscheme_oc%rMat(1,:)*rscheme_oc%rnorm
      else
         dat(1,:)=rscheme_oc%drMat(1,:)*rscheme_oc%rnorm
      end if
      if ( kbotxi == 1 ) then
         dat(n_r_max,:)=rscheme_oc%rMat(n_r_max,:)* rscheme_oc%rnorm
      else
         dat(n_r_max,:)=rscheme_oc%drMat(n_r_max,:)*rscheme_oc%rnorm
      end if

      !-- Fill with zeros:
      if ( rscheme_oc%n_max < n_r_max ) then
         do n_r_out=rscheme_oc%n_max+1,n_r_max
            dat(1,n_r_out)      =0.0_cp
            dat(n_r_max,n_r_out)=0.0_cp
         end do
      end if

      !-- Renormalize:
      do n_r=1,n_r_max
         dat(n_r,1)      =rscheme_oc%boundary_fac*dat(n_r,1)
         dat(n_r,n_r_max)=rscheme_oc%boundary_fac*dat(n_r,n_r_max)
      end do

      !-- Array copy
      call xi0Mat%set_data(dat)

      !-- Invert matrix:
      call xi0Mat%prepare(info)
      if ( info /= 0 ) call abortRun('! Singular Matrix xi0Mat in init_xi!')

      !-- Set source terms in RHS:
      do n_r=2,n_r_max-1
         rhs(n_r)=-epscxi
      end do

      !-- Set boundary values:
      if ( ktopxi == 2 .and. kbotxi == 2 ) then
         rhs(1)=0.0_cp
      else
         rhs(1)=real(topxi(0,0))
      end if
      rhs(n_r_max)=real(botxi(0,0))

      !-- Solve for xi0:
      call xi0Mat%solve(rhs)

      !-- Copy result to xi0:
      xi0(:)=rhs(:)

      !-- Set cheb-modes > rscheme_oc%n_max to zero:
      if ( rscheme_oc%n_max < n_r_max ) then
         do n_r_out=rscheme_oc%n_max+1,n_r_max
            xi0(n_r_out)=0.0_cp
         end do
      end if

      !-- Transform to radial space:
      call rscheme_oc%costf1(xi0)

      !-- Deallocate arrays
      call xi0Mat%finalize()

   end subroutine xi_cond
!--------------------------------------------------------------------------------
   subroutine pt_cond(t0,p0)
      !
      ! Purpose of this subroutine is to solve the entropy equation
      ! for an the conductive (l=0,m=0)-mode.
      ! Output is the radial dependence of the solution in t0 and p0.
      !

      real(cp), intent(out) :: t0(:) ! spherically-symmetric temperature
      real(cp), intent(out) :: p0(:) ! spherically-symmetric pressure

      !-- local variables:
      integer :: n_cheb,nCheb_p,n_r,n_r_p,info,n_cheb_in
      integer :: n_r_out, n_r_out_p
      real(cp), allocatable :: work(:), work2(:), rhs(:)
      integer, allocatable :: pt0Pivot(:)
      real(cp), allocatable :: pt0Mat_fac(:)
      real(cp), allocatable :: pt0Mat(:, :)

      allocate ( rhs(2*n_r_max), work(n_r_max), work2(n_r_max) )
      allocate ( pt0Pivot(2*n_r_max), pt0Mat_fac(2*n_r_max) )
      allocate ( pt0Mat(2*n_r_max,2*n_r_max) )

      if ( l_temperature_diff ) then

         do n_r_out=1,n_r_max
            n_r_out_p=n_r_out+n_r_max
            do n_r=1,n_r_max
               n_r_p=n_r+n_r_max

               ! Delta T = epsc
               pt0Mat(n_r,n_r_out)=rscheme_oc%rnorm*opr*kappa(n_r)* (        &
               &                          rscheme_oc%d2rMat(n_r,n_r_out) +   &
               &         ( beta(n_r)+two*or1(n_r)+dLkappa(n_r) )*            &
               &                           rscheme_oc%drMat(n_r,n_r_out) )

               pt0Mat(n_r,n_r_out_p)=0.0_cp

               ! Hydrostatic equilibrium
               pt0Mat(n_r_p,n_r_out) = -rscheme_oc%rnorm*rho0(n_r)*BuoFac*   &
               &                        rgrav(n_r)*alpha0(n_r)*              &
               &                         rscheme_oc%rMat(n_r,n_r_out)
               pt0Mat(n_r_p,n_r_out_p)= rscheme_oc%rnorm *(                  &
               &                           rscheme_oc%drMat(n_r,n_r_out)+    &
               &                      ViscHeatFac*BuoFac*(                   &
               &               ThExpNb*alpha0(n_r)*temp0(n_r)+ogrun(n_r) )*  &
               &                      alpha0(n_r)*rgrav(n_r)*                &
               &                            rscheme_oc%rMat(n_r,n_r_out) )

            end do
         end do

      else ! entropy diffusion

         do n_r_out=1,n_r_max
            n_r_out_p=n_r_out+n_r_max
            do n_r=1,n_r_max
               n_r_p=n_r+n_r_max

               ! Delta T = epsc
               pt0Mat(n_r,n_r_out)=rscheme_oc%rnorm*opr*kappa(n_r)* (           &
               &                              rscheme_oc%d2rMat(n_r,n_r_out) +  &
               & ( beta(n_r)-dLtemp0(n_r)+two*or1(n_r)+dLkappa(n_r) )*          &
               &                               rscheme_oc%drMat(n_r,n_r_out) -  &
               &         (ddLtemp0(n_r)+dLtemp0(n_r)*(dLkappa(n_r)+beta(n_r)    &
               &          +two*or1(n_r)))*     rscheme_oc%rMat(n_r,n_r_out) )
               pt0Mat(n_r,n_r_out_p)=-rscheme_oc%rnorm*opr*kappa(n_r)*          &
               &               alpha0(n_r)*temp0(n_r)*orho1(n_r)*ViscHeatFac*   &
               &          ThExpNb*(           rscheme_oc%d2rMat(n_r,n_r_out) +  &
               &          ( dLkappa(n_r)+dLtemp0(n_r)+two*or1(n_r)+             &
               &                   two*dLalpha0(n_r)-beta(n_r) ) *              &
               &                              rscheme_oc%drMat(n_r,n_r_out) +   &
               &                   ((dLalpha0(n_r)-beta(n_r))*( two*or1(n_r)+   &
               &                    dLalpha0(n_r)+dLkappa(n_r)+dLtemp0(n_r) )   &
               &                    +ddLalpha0(n_r)-dbeta(n_r) ) *              &
               &                              rscheme_oc%rMat(n_r,n_r_out))

               ! Hydrostatic equilibrium
               pt0Mat(n_r_p,n_r_out) = -rscheme_oc%rnorm*rho0(n_r)*BuoFac*    &
               &                       rgrav(n_r)*alpha0(n_r)*                &
               &                       rscheme_oc%rMat(n_r,n_r_out)
               pt0Mat(n_r_p,n_r_out_p)= rscheme_oc%rnorm *(                   &
               &                       rscheme_oc%drMat(n_r,n_r_out)+         &
               &                      ViscHeatFac*BuoFac*(                    &
               &                 ThExpNb*alpha0(n_r)*temp0(n_r)+ogrun(n_r) )* &
               &                  alpha0(n_r)*rgrav(n_r)*                     &
               &                        rscheme_oc%rMat(n_r,n_r_out) )

            end do
         end do

      end if

      !-- Set boundary conditions:
      if ( ktops == 1 .or. kbots == 2 .or. kbots == 4 ) then
         pt0Mat(1,1:n_r_max)=otemp1(1)*rscheme_oc%rnorm* &
         &                 rscheme_oc%rMat(1,1:n_r_max)
         pt0Mat(1,n_r_max+1:)=-rscheme_oc%rnorm*ViscHeatFac*ThExpNb* &
         &                   alpha0(1)*orho1(1)*rscheme_oc%rMat(1,1:n_r_max)
      else if ( ktops == 2) then ! constant entropy flux at outer boundary
         if ( rscheme_oc%version == 'cheb' ) then
            pt0Mat(1,1:n_r_max) =rscheme_oc%rnorm*otemp1(1)*(            &
            &                             rscheme_oc%drMat(1,1:n_r_max)- &
            &                dLtemp0(1)*   rscheme_oc%rMat(1,1:n_r_max) )
            pt0Mat(1,n_r_max+1:)=-rscheme_oc%rnorm*ViscHeatFac*ThExpNb*alpha0(1)* &
            &                orho1(1)*(            rscheme_oc%drMat(1,1:n_r_max)+ &
            &                 (dLalpha0(1)-beta(1))* rscheme_oc%rMat(1,1:n_r_max) )
         else
            pt0Mat(1,1:n_r_max)=-otemp1(1)*dLtemp0(1)*rscheme_oc%rMat(1,1:n_r_max)
            pt0Mat(1,1:rscheme_oc%order_boundary+1)=                      &
            &                    pt0Mat(1,1:rscheme_oc%order_boundary+1)+ &
            &                       otemp1(1)*rscheme_oc%dr_top(1,:)
            pt0Mat(1,n_r_max+1:)=-ViscHeatFac*ThExpNb*alpha0(1)*         &
            &          orho1(1)*(dLalpha0(1)-beta(1))*rscheme_oc%rMat(1,1:n_r_max)
            pt0Mat(1,n_r_max+1:n_r_max+rscheme_oc%order_boundary+1)=       &
            &    pt0Mat(1,n_r_max+1:n_r_max+rscheme_oc%order_boundary+1) - &
            &    ViscHeatFac*ThExpNb*alpha0(1)*orho1(1)*rscheme_oc%dr_top(1,:)
         end if
      else if ( ktops == 3) then ! constant temperature at outer boundary
         pt0Mat(1,1:n_r_max)  =rscheme_oc%rnorm*rscheme_oc%rMat(1,1:n_r_max)
         pt0Mat(1,n_r_max+1:) =0.0_cp
      else if ( ktops == 4) then ! constant temperature flux at outer boundary
         if ( rscheme_oc%version == 'cheb' ) then
            pt0Mat(1,1:n_r_max) =rscheme_oc%drMat(1,1:n_r_max)*rscheme_oc%rnorm
         else
            pt0Mat(1,1:rscheme_oc%order_boundary+1)=rscheme_oc%dr_top(1,:)
            pt0Mat(1,rscheme_oc%order_boundary+2:n_r_max)=0.0_cp
         end if
         pt0Mat(1,n_r_max+1:)=0.0_cp
      end if

      if ( kbots == 1 ) then        ! Constant entropy at inner boundary
         pt0Mat(n_r_max,1:n_r_max)=rscheme_oc%rMat(n_r_max,1:n_r_max)* &
         &                       rscheme_oc%rnorm*otemp1(n_r_max)
         pt0Mat(n_r_max,n_r_max+1:)=-rscheme_oc%rnorm*ViscHeatFac*ThExpNb* &
         &                        alpha0(n_r_max)*orho1(n_r_max)*          &
         &                       rscheme_oc%rMat(n_r_max,1:n_r_max)
      else if ( kbots == 2 ) then   ! Constant entropy flux at inner boundary
         if ( rscheme_oc%version == 'cheb' ) then
            pt0Mat(n_r_max,1:n_r_max) =rscheme_oc%rnorm*otemp1(n_r_max)*(   &
            &                          rscheme_oc%drMat(n_r_max,1:n_r_max)- &
            &          dLtemp0(n_r_max)*rscheme_oc%rMat(n_r_max,1:n_r_max) )
            pt0Mat(n_r_max,n_r_max+1:)=-rscheme_oc%rnorm*ViscHeatFac*ThExpNb* &
            &           alpha0(n_r_max)*orho1(n_r_max)*(                      &
            &                          rscheme_oc%drMat(n_r_max,1:n_r_max)+   &
            &          (dLalpha0(n_r_max)-beta(n_r_max))*                     &
            &                           rscheme_oc%rMat(n_r_max,1:n_r_max) )
         else
            pt0Mat(n_r_max,1:n_r_max)=-otemp1(n_r_max)*dLtemp0(n_r_max)* &
            &                         rscheme_oc%rMat(n_r_max,1:n_r_max)
            pt0Mat(n_r_max,n_r_max:n_r_max-rscheme_oc%order_boundary:-1)=      &
            &     pt0Mat(n_r_max,n_r_max:n_r_max-rscheme_oc%order_boundary:-1)+&
            &     otemp1(n_r_max)*rscheme_oc%dr_bot(1,:)
            pt0Mat(n_r_max,n_r_max+1:)=-ViscHeatFac*ThExpNb*alpha0(n_r_max)*  &
            &              orho1(n_r_max)*(dLalpha0(n_r_max)-beta(n_r_max))*  &
            &              rscheme_oc%rMat(n_r_max,1:n_r_max)
            pt0Mat(n_r_max,2*n_r_max:2*n_r_max-rscheme_oc%order_boundary:-1)=  &
            & pt0Mat(n_r_max,2*n_r_max:2*n_r_max-rscheme_oc%order_boundary:-1)-&
            &      ViscHeatFac*ThExpNb*alpha0(n_r_max)*orho1(n_r_max)*         &
            &           rscheme_oc%dr_bot(1,:)
         end if
      else if ( kbots == 3 ) then   ! Constant temperature at inner boundary
         pt0Mat(n_r_max,1:n_r_max)=rscheme_oc%rMat(n_r_max,1:n_r_max)* &
         &                         rscheme_oc%rnorm
         pt0Mat(n_r_max,n_r_max+1:)=0.0_cp
      else if ( kbots == 4 ) then   ! Constant temperature flux at inner boundary
         if ( rscheme_oc%version == 'cheb' ) then
            pt0Mat(n_r_max,1:n_r_max)=rscheme_oc%drMat(n_r_max,1:n_r_max)* &
            &                         rscheme_oc%rnorm
         else
            pt0Mat(n_r_max,1:n_r_max)=0.0_cp
            pt0Mat(n_r_max,n_r_max:n_r_max-rscheme_oc%order_boundary:-1)= &
            &        rscheme_oc%dr_bot(1,:)
         end if
         pt0Mat(n_r_max,n_r_max+1:)=0.0_cp
      end if

      if ( rscheme_oc%version == 'cheb' ) then
         pt0Mat(2*n_r_max,1:n_r_max) =0.0_cp
         pt0Mat(2*n_r_max,n_r_max+1:)=0.0_cp
      end if

      ! In case density perturbations feed back on pressure (non-Boussinesq)
      ! Impose that the integral of (rho' r^2) vanishes

      if ( ViscHeatFac*ThExpNb /= 0.0_cp .and. ktopp==1 ) then

         work(:)=ViscHeatFac*alpha0(:)*(ThExpNb*alpha0(:)*temp0(:)+&
         &       ogrun(:))*r(:)*r(:)

         work2(:)=-alpha0(:)*rho0(:)*r(:)*r(:)
         call rscheme_oc%costf1(work2)

         select type(rscheme_oc)

         type is(type_cheb_odd)
            call rscheme_oc%costf1(work)
            work(:)      =work(:)*rscheme_oc%rnorm
            work(1)      =rscheme_oc%boundary_fac*work(1)
            work(n_r_max)=rscheme_oc%boundary_fac*work(n_r_max)
            work2(:)     =work2(:)*rscheme_oc%rnorm
            work2(1)      =rscheme_oc%boundary_fac*work2(1)
            work2(n_r_max)=rscheme_oc%boundary_fac*work2(n_r_max)

            do n_cheb=1,rscheme_oc%n_max
               nCheb_p=n_cheb+n_r_max
               pt0Mat(n_r_max+1,nCheb_p)=0.0_cp
               pt0Mat(n_r_max+1,n_cheb) =0.0_cp
               do n_cheb_in=1,rscheme_oc%n_max
                  if (mod(n_cheb+n_cheb_in-2,2)==0) then
                     pt0Mat(n_r_max+1,nCheb_p)=pt0Mat(n_r_max+1,nCheb_p)+              &
                     &                       (one/(one-real(n_cheb_in-n_cheb,cp)**2)+  &
                     &                       one/(one-real(n_cheb_in+n_cheb-2,cp)**2))*&
                     &                       work(n_cheb_in)*half*rscheme_oc%rnorm
                     pt0Mat(n_r_max+1,n_cheb)=pt0Mat(n_r_max+1,n_cheb)+                &
                     &                       (one/(one-real(n_cheb_in-n_cheb,cp)**2)+  &
                     &                       one/(one-real(n_cheb_in+n_cheb-2,cp)**2))*&
                     &                       work2(n_cheb_in)*half*rscheme_oc%rnorm
                  end if
               end do
            end do

         type is(type_fd)

            !-- In the finite differences case, we restrict the integral boundary
            !-- condition to a trapezoidal rule of integration
            do n_r_out=2,rscheme_oc%n_max-1
               n_r_out_p=n_r_out+n_r_max
               pt0Mat(n_r_max+1,n_r_out)  =half*work2(n_r_out)*        &
               &                           ( r(n_r_out+1)-r(n_r_out-1) )
               pt0Mat(n_r_max+1,n_r_out_p)=half* work(n_r_out)*        &
               &                           ( r(n_r_out+1)-r(n_r_out-1) )
            end do
            pt0Mat(n_r_max+1,1)        =half*work2(1)*( r(2)-r(1) )
            pt0Mat(n_r_max+1,n_r_max+1)=half* work(1)*( r(2)-r(1) )
            pt0Mat(n_r_max+1,n_r_max)  =half*work2(n_r_max)*( r(n_r_max)-r(n_r_max-1) )
            pt0Mat(n_r_max+1,2*n_r_max)=half* work(n_r_max)*( r(n_r_max)-r(n_r_max-1) )

         end select

      else

         do n_r_out=1,rscheme_oc%n_max
            n_r_out_p=n_r_out+n_r_max
            pt0Mat(n_r_max+1,n_r_out) =0.0_cp
            pt0Mat(n_r_max+1,n_r_out_p)=rscheme_oc%rnorm*rscheme_oc%rMat(1,n_r_out)
         end do

      end if

      !-- Fill with zeros:
      if ( rscheme_oc%n_max < n_r_max ) then
         do n_r_out=rscheme_oc%n_max+1,n_r_max
            n_r_out_p=n_r_out+n_r_max
            pt0Mat(1,n_r_out)          =0.0_cp
            pt0Mat(n_r_max,n_r_out)    =0.0_cp
            pt0Mat(n_r_max+1,n_r_out)  =0.0_cp
            pt0Mat(2*n_r_max,n_r_out)  =0.0_cp
            pt0Mat(1,n_r_out_p)        =0.0_cp
            pt0Mat(n_r_max,n_r_out_p)  =0.0_cp
            pt0Mat(n_r_max+1,n_r_out_p)=0.0_cp
         end do
      end if

      !-- Renormalize:
      do n_r=1,n_r_max
         n_r_p=n_r+n_r_max
         pt0Mat(n_r,1)          =rscheme_oc%boundary_fac*pt0Mat(n_r,1)
         pt0Mat(n_r,n_r_max)    =rscheme_oc%boundary_fac*pt0Mat(n_r,n_r_max)
         pt0Mat(n_r,n_r_max+1)  =rscheme_oc%boundary_fac*pt0Mat(n_r,n_r_max+1)
         pt0Mat(n_r,2*n_r_max)  =rscheme_oc%boundary_fac*pt0Mat(n_r,2*n_r_max)
         pt0Mat(n_r_p,1)        =rscheme_oc%boundary_fac*pt0Mat(n_r_p,1)
         pt0Mat(n_r_p,n_r_max)  =rscheme_oc%boundary_fac*pt0Mat(n_r_p,n_r_max)
         pt0Mat(n_r_p,n_r_max+1)=rscheme_oc%boundary_fac*pt0Mat(n_r_p,n_r_max+1)
         pt0Mat(n_r_p,2*n_r_max)=rscheme_oc%boundary_fac*pt0Mat(n_r_p,2*n_r_max)
      end do

      ! compute the linesum of each line
      do n_r=1,2*n_r_max
         pt0Mat_fac(n_r)=one/maxval(abs(pt0Mat(n_r,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do n_r=1,2*n_r_max
         pt0Mat(n_r,:) = pt0Mat(n_r,:)*pt0Mat_fac(n_r)
      end do

      !-- Prepare matrix:
      call prepare_mat(pt0Mat,2*n_r_max,2*n_r_max,pt0Pivot,info)
      if ( info /= 0 ) call abortRun('! Singular Matrix pt0Mat in pt_cond!')

      !-- Set source terms in RHS:
      do n_r=1,n_r_max
         rhs(n_r)          =-epsc*epscProf(n_r)*orho1(n_r)
         rhs(n_r+n_r_max)  =0.0_cp
      end do

      !-- Set boundary values:
      if ( (ktops==2 .and. kbots==2) .or. (ktops==4 .and. kbots==4) ) then
         rhs(1)=0.0_cp
      else
         rhs(1)=real(tops(0,0))
      end if
      rhs(n_r_max)=real(bots(0,0))

      !-- Pressure at the top boundary
      rhs(n_r_max+1)=0.0_cp

      rhs(:) = pt0Mat_fac(:)*rhs(:)

      !-- Solve for t0 and p0
      call solve_mat(pt0Mat,2*n_r_max,2*n_r_max,pt0Pivot,rhs)

      !-- Copy result to t0 and p0:
      t0(:)=rhs(1:n_r_max)
      p0(:)=rhs(n_r_max+1:)

      !-- Set cheb-modes > rscheme_oc%n_max to zero:
      if ( rscheme_oc%n_max < n_r_max ) then
         do n_r_out=rscheme_oc%n_max+1,n_r_max
            t0(n_r_out)=0.0_cp
            p0(n_r_out)=0.0_cp
         end do
      end if

      !-- Transform to radial space:
      call rscheme_oc%costf1(t0)
      call rscheme_oc%costf1(p0)

      deallocate ( rhs, work, work2 )
      deallocate ( pt0Pivot, pt0Mat_fac )
      deallocate ( pt0Mat )

   end subroutine pt_cond
!--------------------------------------------------------------------------------
   subroutine ps_cond(s0,p0)
      !
      ! Purpose of this subroutine is to solve the entropy equation
      ! for an the conductive (l=0,m=0)-mode.
      ! Output is the radial dependence of the solution in s0 and p0.
      !

      real(cp), intent(out) :: s0(:) ! spherically-symmetric part
      real(cp), intent(out) :: p0(:) ! spherically-symmetric part

      !-- local variables:
      integer :: n_cheb,nCheb_p,n_r,n_r_p,info,n_cheb_in
      integer :: n_r_out, n_r_out_p
      real(cp), allocatable :: work(:), work2(:), rhs(:)
      integer, allocatable :: ps0Pivot(:)
      real(cp), allocatable :: ps0Mat_fac(:)
      real(cp), allocatable :: ps0Mat(:, :)

      allocate ( rhs(2*n_r_max), work(n_r_max), work2(n_r_max) )
      allocate ( ps0Pivot(2*n_r_max), ps0Mat_fac(2*n_r_max) )
      allocate ( ps0Mat(2*n_r_max,2*n_r_max) )

      if ( l_temperature_diff ) then

         do n_r_out=1,n_r_max
            n_r_out_p=n_r_out+n_r_max
            do n_r=1,n_r_max
               n_r_p=n_r+n_r_max

               ! Delta T = epsc
               ps0Mat(n_r,n_r_out)=rscheme_oc%rnorm*opr*kappa(n_r)* (          &
               &                            rscheme_oc%d2rMat(n_r,n_r_out) +   &
               &      ( beta(n_r)+two*dLtemp0(n_r)+two*or1(n_r)+dLkappa(n_r) )*&
               &                             rscheme_oc%drMat(n_r,n_r_out) +   &
               &      ( ddLtemp0(n_r)+dLtemp0(n_r)*(                           &
               &  two*or1(n_r)+dLkappa(n_r)+dLtemp0(n_r)+beta(n_r) ) ) *       &
               &                             rscheme_oc%rMat(n_r,n_r_out) )

               ps0Mat(n_r,n_r_out_p)=rscheme_oc%rnorm*opr*kappa(n_r)*          &
               &       alpha0(n_r)*orho1(n_r)*ViscHeatFac*ThExpNb*(            &
               &                             rscheme_oc%d2rMat(n_r,n_r_out) +  &
               &      ( dLkappa(n_r)+two*(dLalpha0(n_r)+dLtemp0(n_r)) -        &
               &        beta(n_r) +two*or1(n_r) ) *                            &
               &                             rscheme_oc%drMat(n_r,n_r_out) +   &
               & ( (dLkappa(n_r)+dLalpha0(n_r)+dLtemp0(n_r)+two*or1(n_r)) *    &
               &        (dLalpha0(n_r)+dLtemp0(n_r)-beta(n_r)) +               &
               &        ddLalpha0(n_r)+ddLtemp0(n_r)-dbeta(n_r) ) *            &
               &                             rscheme_oc%rMat(n_r,n_r_out) )

               ! Hydrostatic equilibrium
               ps0Mat(n_r_p,n_r_out) = -rscheme_oc%rnorm*rho0(n_r)*BuoFac* &
               &                        rgrav(n_r)*rscheme_oc%rMat(n_r,n_r_out)
               ps0Mat(n_r_p,n_r_out_p)= rscheme_oc%rnorm*(                     &
               &                                rscheme_oc%drMat(n_r,n_r_out)- &
               &                       beta(n_r)*rscheme_oc%rMat(n_r,n_r_out) )

            end do
         end do

      else ! entropy diffusion

         do n_r_out=1,n_r_max
            n_r_out_p=n_r_out+n_r_max
            do n_r=1,n_r_max
               n_r_p=n_r+n_r_max

               ! Delta T = epsc
               ps0Mat(n_r,n_r_out)=rscheme_oc%rnorm*opr*kappa(n_r)* (        &
               &                           rscheme_oc%d2rMat(n_r,n_r_out) +  &
               &      ( beta(n_r)+dLtemp0(n_r)+two*or1(n_r)+dLkappa(n_r) )*  &
               &                           rscheme_oc%drMat(n_r,n_r_out) )
               ps0Mat(n_r,n_r_out_p)=0.0_cp

               ! Hydrostatic equilibrium
               ps0Mat(n_r_p,n_r_out)=-rscheme_oc%rnorm*rho0(n_r)*BuoFac*   &
               &                     rgrav(n_r)*rscheme_oc%rMat(n_r,n_r_out)
               ps0Mat(n_r_p,n_r_out_p)= rscheme_oc%rnorm*(                 &
               &                            rscheme_oc%drMat(n_r,n_r_out)- &
               &                   beta(n_r)*rscheme_oc%rMat(n_r,n_r_out) )

            end do
         end do

      end if

      !-- Set boundary conditions:
      if ( ktops == 1 .or. kbots == 2 .or. kbots == 4 ) then
         ps0Mat(1,1:n_r_max)  =rscheme_oc%rnorm*rscheme_oc%rMat(1,:)
         ps0Mat(1,n_r_max+1:)=0.0_cp
      else if ( ktops == 2) then ! constant entropy flux at CMB
         ps0Mat(1,1:n_r_max) =rscheme_oc%drMat(1,:)*rscheme_oc%rnorm
         ps0Mat(1,n_r_max+1:)=0.0_cp
      else if ( ktops == 3) then ! constant temperature at CMB
         ps0Mat(1,1:n_r_max) =rscheme_oc%rnorm*temp0(1)*rscheme_oc%rMat(1,:)
         ps0Mat(1,n_r_max+1:)=rscheme_oc%rnorm*alpha0(1)*temp0(1)*orho1(1)* &
          &                   ViscHeatFac*ThExpNb*rscheme_oc%rMat(1,:)
      else if ( ktops == 4) then ! constant temperature flux at CMB
         ps0Mat(1,1:n_r_max) =rscheme_oc%rnorm*temp0(1)*(              &
         &                                      rscheme_oc%drMat(1,:)+ &
         &                            dLtemp0(1)*rscheme_oc%rMat(1,:) )
         ps0Mat(1,n_r_max+1:)=rscheme_oc%rnorm*orho1(1)*alpha0(1)*      &
         &                    temp0(1)*ViscHeatFac*ThExpNb*(            &
         &                       rscheme_oc%drMat(1,:)+(dLalpha0(1)+    &
         &                     dLtemp0(1)-beta(1))*rscheme_oc%rMat(1,:) )
      end if

      if ( kbots == 1 ) then        ! Constant entropy at ICB
         ps0Mat(n_r_max,1:n_r_max)=rscheme_oc%rMat(n_r_max,:)* &
         &                         rscheme_oc%rnorm
         ps0Mat(n_r_max,n_r_max+1:)=0.0_cp
      else if ( kbots == 2 ) then   ! Constant entropy flux at ICB
         ps0Mat(n_r_max,1:n_r_max)=rscheme_oc%drMat(n_r_max,:)* &
         &                         rscheme_oc%rnorm
         ps0Mat(n_r_max,n_r_max+1:)=0.0_cp
      else if ( kbots == 3 ) then   ! Constant temperature at ICB
         ps0Mat(n_r_max,1:n_r_max)=rscheme_oc%rnorm* &
         &                       rscheme_oc%rMat(n_r_max,:)*temp0(n_r_max)
         ps0Mat(n_r_max,n_r_max+1:)=rscheme_oc%rnorm*                 &
         &                         rscheme_oc%rMat(n_r_max,:)*        &
         &                           alpha0(n_r_max)*temp0(n_r_max)*  &
         &                         orho1(n_r_max)*ViscHeatFac*ThExpNb
      else if ( kbots == 4 ) then   ! Constant temperature flux at ICB
         ps0Mat(n_r_max,1:n_r_max)=rscheme_oc%rnorm*temp0(n_r_max)*(    &
         &                       rscheme_oc%drMat(n_r_max,:)+           &
         &      dLtemp0(n_r_max)*rscheme_oc%rMat(n_r_max,:) )
         ps0Mat(n_r_max,n_r_max+1:)=rscheme_oc%rnorm*orho1(n_r_max)*    &
         &                         alpha0(n_r_max)*temp0(n_r_max)*      &
         &                         ViscHeatFac*ThExpNb*(                &
         &                         rscheme_oc%drMat(n_r_max,:)+         &
         &                        (dLalpha0(n_r_max)+dLtemp0(n_r_max)-  &
         &               beta(n_r_max))*rscheme_oc%rMat(n_r_max,:) )
      end if

      if ( rscheme_oc%version == 'cheb' ) then
         ps0Mat(n_r_max+1,1:n_r_max)=0.0_cp
      end if

      ! In case density perturbations feed back on pressure (non-Boussinesq)
      ! Impose that the integral of (rho' r^2) vanishes

      if ( ViscHeatFac*ThExpNb /= 0.0_cp .and. ktopp == 1 ) then

         work(:)=ThExpNb*ViscHeatFac*ogrun(:)*alpha0(:)*r(:)*r(:)
         call rscheme_oc%costf1(work)
         work         =work*rscheme_oc%rnorm
         work(1)      =rscheme_oc%boundary_fac*work(1)
         work(n_r_max)=rscheme_oc%boundary_fac*work(n_r_max)

         work2(:)=-ThExpNb*alpha0(:)*temp0(:)*rho0(:)*r(:)*r(:)
         call rscheme_oc%costf1(work2)
         work2         =work2*rscheme_oc%rnorm
         work2(1)      =rscheme_oc%boundary_fac*work2(1)
         work2(n_r_max)=rscheme_oc%boundary_fac*work2(n_r_max)

         if ( rscheme_oc%version == 'cheb' ) then

            do n_cheb=1,rscheme_oc%n_max
               nCheb_p=n_cheb+n_r_max
               ps0Mat(n_r_max+1,nCheb_p)=0.0_cp
               ps0Mat(n_r_max+1,n_cheb)=0.0_cp
               do n_cheb_in=1,rscheme_oc%n_max
                  if (mod(n_cheb+n_cheb_in-2,2)==0) then
                  ps0Mat(n_r_max+1,nCheb_p)=ps0Mat(n_r_max+1,nCheb_p)+              &
                  &                       (one/(one-real(n_cheb_in-n_cheb,cp)**2)+  &
                  &                       one/(one-real(n_cheb_in+n_cheb-2,cp)**2))*&
                  &                       work(n_cheb_in)*half*rscheme_oc%rnorm
                  ps0Mat(n_r_max+1,n_cheb)=ps0Mat(n_r_max+1,n_cheb)+                &
                  &                       (one/(one-real(n_cheb_in-n_cheb,cp)**2)+  &
                  &                       one/(one-real(n_cheb_in+n_cheb-2,cp)**2))*&
                  &                       work2(n_cheb_in)*half*rscheme_oc%rnorm
                  end if
               end do
            end do

         else

            !-- In the finite differences case, we restrict the integral boundary
            !-- condition to a trapezoidal rule of integration
            do n_r_out=2,rscheme_oc%n_max-1
               n_r_out_p=n_r_out+n_r_max
               ps0Mat(n_r_max+1,n_r_out)  =half*work2(n_r_out)*        &
               &                           ( r(n_r_out+1)-r(n_r_out-1) )
               ps0Mat(n_r_max+1,n_r_out_p)=half*work(n_r_out)*         &
               &                           ( r(n_r_out+1)-r(n_r_out-1) )
            end do
            ps0Mat(n_r_max+1,1)        =half*work2(1)*( r(2)-r(1) )
            ps0Mat(n_r_max+1,n_r_max+1)=half* work(1)*( r(2)-r(1) )
            ps0Mat(n_r_max+1,n_r_max)  =half*work2(n_r_max)*( r(n_r_max)-r(n_r_max-1) )
            ps0Mat(n_r_max+1,2*n_r_max)=half* work(n_r_max)*( r(n_r_max)-r(n_r_max-1) )

         end if

      else

         do n_r_out=1,n_r_max
            n_r_out_p=n_r_out+n_r_max
            ps0Mat(n_r_max+1,n_r_out)  =0.0_cp
            ps0Mat(n_r_max+1,n_r_out_p)=rscheme_oc%rnorm*rscheme_oc%rMat(1,n_r_out)
         end do

      end if

      !-- Fill with zeros:
      if ( rscheme_oc%n_max < n_r_max ) then
         do n_r_out=rscheme_oc%n_max+1,n_r_max
            n_r_out_p=n_r_out+n_r_max
            ps0Mat(1,n_r_out)          =0.0_cp
            ps0Mat(n_r_max,n_r_out)    =0.0_cp
            ps0Mat(n_r_max+1,n_r_out)  =0.0_cp
            ps0Mat(2*n_r_max,n_r_out)  =0.0_cp
            ps0Mat(1,n_r_out_p)        =0.0_cp
            ps0Mat(n_r_max,n_r_out_p)  =0.0_cp
            ps0Mat(n_r_max+1,n_r_out_p)=0.0_cp
         end do
      end if

      !-- Renormalize:
      do n_r=1,n_r_max
         n_r_p=n_r+n_r_max
         ps0Mat(n_r,1)          =rscheme_oc%boundary_fac*ps0Mat(n_r,1)
         ps0Mat(n_r,n_r_max)    =rscheme_oc%boundary_fac*ps0Mat(n_r,n_r_max)
         ps0Mat(n_r,n_r_max+1)  =rscheme_oc%boundary_fac*ps0Mat(n_r,n_r_max+1)
         ps0Mat(n_r,2*n_r_max)  =rscheme_oc%boundary_fac*ps0Mat(n_r,2*n_r_max)
         ps0Mat(n_r_p,1)        =rscheme_oc%boundary_fac*ps0Mat(n_r_p,1)
         ps0Mat(n_r_p,n_r_max)  =rscheme_oc%boundary_fac*ps0Mat(n_r_p,n_r_max)
         ps0Mat(n_r_p,n_r_max+1)=rscheme_oc%boundary_fac*ps0Mat(n_r_p,n_r_max+1)
         ps0Mat(n_r_p,2*n_r_max)=rscheme_oc%boundary_fac*ps0Mat(n_r_p,2*n_r_max)
      end do

      ! compute the linesum of each line
      do n_r=1,2*n_r_max
         ps0Mat_fac(n_r)=one/maxval(abs(ps0Mat(n_r,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do n_r=1,2*n_r_max
         ps0Mat(n_r,:) = ps0Mat(n_r,:)*ps0Mat_fac(n_r)
      end do

      !-- Invert matrix:
      call prepare_mat(ps0Mat,2*n_r_max,2*n_r_max,ps0Pivot,info)
      if ( info /= 0 ) then
         call abortRun('! Singular Matrix ps0Mat in ps_cond!')
      end if

      !-- Set source terms in RHS:
      do n_r=1,n_r_max
         rhs(n_r)        =-epsc*epscProf(n_r)*orho1(n_r)
         rhs(n_r+n_r_max)=0.0_cp
      end do

      !-- Set boundary values:
      if ( (ktops==2 .and. kbots==2) .or. (ktops==4 .and. kbots==4) ) then
         rhs(1)=0.0_cp
      else
         rhs(1)=real(tops(0,0))
      end if
      rhs(n_r_max)=real(bots(0,0))

      !-- Pressure at the top boundary
      rhs(n_r_max+1)=0.0_cp

      rhs(:)=ps0Mat_fac(:)*rhs(:)

      !-- Solve for s0 and p0
      call solve_mat(ps0Mat,2*n_r_max,2*n_r_max,ps0Pivot,rhs)

      !-- Copy result to s0 and p0
      s0(:)=rhs(1:n_r_max)
      p0(:)=rhs(n_r_max+1:)

      !-- Set cheb-modes > rscheme_oc%n_max to zero:
      if ( rscheme_oc%n_max < n_r_max ) then
         do n_cheb=rscheme_oc%n_max+1,n_r_max
            s0(n_cheb)=0.0_cp
            p0(n_cheb)=0.0_cp
         end do
      end if

      !-- Transform to radial space:
      call rscheme_oc%costf1(s0)
      call rscheme_oc%costf1(p0)

      deallocate ( rhs, work, work2 )
      deallocate ( ps0Pivot, ps0Mat_fac )
      deallocate ( ps0Mat )

   end subroutine ps_cond
!--------------------------------------------------------------------------------
end module init_fields
