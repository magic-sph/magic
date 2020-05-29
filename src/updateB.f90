#include "perflib_preproc.cpp"
module updateB_mod

#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif

   use omp_lib
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, n_r_tot, n_r_ic_max, n_cheb_ic_max,    &
       &                 n_r_ic_maxMag, n_r_maxMag, n_r_totMag, n_r_cmb, &
       &                 n_r_icb, get_openmp_blocks, n_lo_loc, n_mlo_loc,&
       &                 n_mloMag_loc
   use LMmapping, only: map_mlo, map_glbl_st
   use radial_functions, only: chebt_ic,or2,r_cmb,chebt_ic_even, d2cheb_ic,    &
       &                       cheb_norm_ic,dr_fac_ic,lambda,dLlambda,o_r_ic,r,&
       &                       or1, cheb_ic, dcheb_ic, rscheme_oc, dr_top_ic
   use physical_parameters, only: n_r_LCR, opm, O_sr, kbotb, imagcon, tmagcon, &
       &                         sigma_ratio, conductance_ma, ktopb
   use init_fields, only: bpeaktop, bpeakbot
   use num_param, only: solve_counter, dct_counter
   use horizontal_data, only: hdif_B
   use logic, only: l_cond_ic, l_LCR, l_rot_ic, l_mag_nl, l_b_nl_icb, &
       &            l_b_nl_cmb, l_update_b, l_RMS, l_finite_diff,     &
       &            l_full_sphere
   use RMS, only: dtBPolLMr, dtBPol2hInt, dtBTor2hInt
   use constants, only: pi, zero, one, two, three, half
   use special
   use RMS_helpers, only: hInt2PolLM, hInt2TorLM
   use fields, only: work_LMdist
   use radial_der_even, only: get_ddr_even
   use radial_der, only: get_dr, get_ddr
   use useful, only: abortRun
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray
   use dense_matrices
   use band_matrices
   use bordered_matrices
   use real_matrices

   implicit none

   private

   !-- Local work arrays:
   complex(cp), allocatable :: workA(:,:), workB(:,:)
   complex(cp), allocatable :: work_ic_LMdist(:,:)
   real(cp), allocatable :: rhs1(:,:,:),rhs2(:,:,:)
   complex(cp), allocatable :: dtT(:), dtP(:)
   class(type_realmat), pointer :: bMat(:), jMat(:)
#ifdef WITH_PRECOND_BJ
   real(cp), allocatable :: bMat_fac(:,:)
   real(cp), allocatable :: jMat_fac(:,:)
#endif
   logical, public, allocatable :: lBmat(:)

   public :: initialize_updateB, finalize_updateB, updateB, finish_exp_mag, &
   &         get_mag_rhs_imp, get_mag_ic_rhs_imp, finish_exp_mag_ic,        &
   &         assemble_mag

contains

   subroutine initialize_updateB()
      !
      ! Purpose of this subroutine is to allocate the matrices needed
      ! to time advance the magnetic field. Depending on the
      ! radial scheme and the inner core conductivity, it can be full,
      ! bordered or band matrices.
      !

      !-- Local variables
      integer :: maxThreads, ll, n_bandsJ, n_bandsB

      if ( l_finite_diff ) then
         if ( l_cond_ic ) then
            allocate( type_bordmat :: jMat(n_lo_loc) )
            allocate( type_bordmat :: bMat(n_lo_loc) )
         else
            allocate( type_bandmat :: jMat(n_lo_loc) )
            allocate( type_bandmat :: bMat(n_lo_loc) )
         end if

         if ( kbotb == 2 .or. ktopb == 2 .or. conductance_ma /= 0 .or. &
         &    rscheme_oc%order  > 2 .or. rscheme_oc%order_boundary > 2 ) then
            !-- Perfect conductor or conducting mantle
            n_bandsJ = max(2*rscheme_oc%order_boundary+1,rscheme_oc%order+1)
         else
            n_bandsJ = rscheme_oc%order+1
         end if

         if ( conductance_ma /= 0 ) then
            !-- Second derivative on the boundary
            n_bandsB = max(2*rscheme_oc%order_boundary+3,rscheme_oc%order+1)
         else
            n_bandsB = max(2*rscheme_oc%order_boundary+1,rscheme_oc%order+1)
         end if

         if ( l_cond_ic ) then
            do ll=1,n_lo_loc
               call bMat(ll)%initialize(n_bandsB,n_r_max,.true.,n_r_ic_max)
               call jMat(ll)%initialize(n_bandsJ,n_r_max,.true.,n_r_ic_max)
            end do
         else
            do ll=1,n_lo_loc
               call bMat(ll)%initialize(n_bandsB,n_r_tot,l_pivot=.true.)
               call jMat(ll)%initialize(n_bandsJ,n_r_tot,l_pivot=.true.)
            end do
         end if
      else
         allocate( type_densemat :: jMat(n_lo_loc) )
         allocate( type_densemat :: bMat(n_lo_loc) )

         do ll=1,n_lo_loc
            call bMat(ll)%initialize(n_r_tot,n_r_tot,l_pivot=.true.)
            call jMat(ll)%initialize(n_r_tot,n_r_tot,l_pivot=.true.)
         end do
      end if

#ifdef WITH_PRECOND_BJ
      allocate(bMat_fac(n_r_tot,n_lo_loc), jMat_fac(n_r_tot,n_lo_loc))
      bytes_allocated = bytes_allocated+2*n_r_tot*n_lo_loc*SIZEOF_DEF_REAL
#endif
      allocate( lBmat(n_lo_loc) )
      bytes_allocated = bytes_allocated+n_lo_loc*SIZEOF_LOGICAL

      if ( l_RMS ) then
         allocate( workA(n_mlo_loc,n_r_max), workB(n_mlo_loc,n_r_max) )
         bytes_allocated = bytes_allocated+2*n_mlo_loc*n_r_max*SIZEOF_DEF_COMPLEX
      end if

      if ( l_cond_ic ) then
         allocate( work_ic_LMdist(n_mlo_loc,n_r_ic_max) )
         bytes_allocated = bytes_allocated+n_mlo_loc*n_r_ic_max*SIZEOF_DEF_COMPLEX
      end if

      allocate( dtT(n_mlo_loc), dtP(n_mlo_loc) )
      bytes_allocated = bytes_allocated+2*n_mlo_loc*SIZEOF_DEF_COMPLEX

#ifdef WITHOMP
      maxThreads=omp_get_max_threads()
#else
      maxThreads=1
#endif

      allocate(rhs1(n_r_tot,2*maxval(map_mlo%n_mi(:)), 1))
      allocate(rhs2(n_r_tot,2*maxval(map_mlo%n_mi(:)), 1))
      bytes_allocated=bytes_allocated+2*n_r_tot*maxThreads* &
      &               maxval(map_mlo%n_mi(:))*SIZEOF_DEF_COMPLEX

   end subroutine initialize_updateB
!-----------------------------------------------------------------------------
   subroutine finalize_updateB()
      !
      ! This subroutine deallocates the matrices involved in the time
      ! advance of the magnetic field
      !

      !-- Local variables
      integer :: ll

      do ll=1,n_lo_loc
         call jMat(ll)%finalize()
         call bMat(ll)%finalize()
      end do

      if ( l_cond_ic ) deallocate ( work_ic_LMdist )
      deallocate( lBmat )

#ifdef WITH_PRECOND_BJ
      deallocate(bMat_fac,jMat_fac)
#endif
      deallocate( dtT, dtP, rhs1, rhs2 )
      if ( l_RMS ) deallocate( workA, workB )

   end subroutine finalize_updateB
!-----------------------------------------------------------------------------
   subroutine updateB(b,db,ddb,aj,dj,ddj,dbdt,djdt,b_ic,db_ic,ddb_ic,aj_ic,  &
              &       dj_ic,ddj_ic,dbdt_ic,djdt_ic,b_nl_cmb,aj_nl_cmb,       &
              &       aj_nl_icb,time,tscheme,lRmsNext)
      !
      !
      !  Calculated update of magnetic field potential and the time
      !  stepping arrays dbdtLast, ...
      !
      !
      !  updates the magnetic field potentials b, aj and
      !  their derivatives,
      !  adds explicit part to time derivatives of b and j

      !  input:  w1 - weight for dbdt-contribution from current time step
      !               (w2=1-w1: weight for contrib. from previous step)
      !          coex - factor depending on weighting alpha of
      !                 implicit contribution
      !          mc_min,mc_max
      !               - range of mca-indices in which field is updated
      !                 (harmonic order m=(mca-1)*minc)
      !          b_nl_cmb = RHS of nonlinear BC for poloidal magnetic field
      !                      potential in the case of stress free CMB
      !          aj_nl_cmb = RHS of nonlinear BC for toroidal magnetic field
      !                      potential in the case of stress free CMB
      !          aj_nl_icb = RHS of nonlinear BC for radial derivative of
      !                      toroidal magnetic field
      !                      potential in the case of stress free ICB
      !

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme
      !@> TODO: the three following arrays could be in n_mloMag_loc space
      complex(cp),         intent(in) :: b_nl_cmb(:)  ! nonlinear BC for b at CMB
      complex(cp),         intent(in) :: aj_nl_cmb(:) ! nonlinear BC for aj at CMB
      complex(cp),         intent(in) :: aj_nl_icb(:) ! nonlinear BC for dr aj at ICB
      real(cp),            intent(in) :: time
      logical,             intent(in) :: lRmsNext

      !-- Input/output of scalar potentials and time stepping arrays:
      type(type_tarray), intent(inout) :: dbdt
      type(type_tarray), intent(inout) :: djdt
      type(type_tarray), intent(inout) :: dbdt_ic
      type(type_tarray), intent(inout) :: djdt_ic
      complex(cp),       intent(inout) :: b(n_mloMag_loc,n_r_maxMag)
      complex(cp),       intent(inout) :: aj(n_mloMag_loc,n_r_maxMag)
      complex(cp),       intent(inout) :: b_ic(n_mloMag_loc,n_r_ic_maxMag)
      complex(cp),       intent(inout) :: aj_ic(n_mloMag_loc,n_r_ic_maxMag)

      !-- Output variables:
      complex(cp), intent(out) :: db(n_mloMag_loc,n_r_maxMag)
      complex(cp), intent(out) :: ddb(n_mloMag_loc,n_r_maxMag)
      complex(cp), intent(out) :: dj(n_mloMag_loc,n_r_maxMag)
      complex(cp), intent(out) :: ddj(n_mloMag_loc,n_r_maxMag)
      complex(cp), intent(out) :: db_ic(n_mloMag_loc,n_r_ic_maxMag)
      complex(cp), intent(out) :: ddb_ic(n_mloMag_loc,n_r_ic_maxMag)
      complex(cp), intent(out) :: dj_ic(n_mloMag_loc,n_r_ic_maxMag)
      complex(cp), intent(out) :: ddj_ic(n_mloMag_loc,n_r_ic_maxMag)

      !-- Local variables:
      real(cp) :: yl0_norm,prefac!External magnetic field of general l

      integer :: l, m              ! degree and order
      integer :: lj, mi            ! l and m indices
      integer :: i                 ! global index
      integer :: n_r_out           ! No of cheb polynome (degree+1)
      integer :: nR                ! No of radial grid point

      !-- for feedback
      integer :: l1m0
      real(cp) :: ff,cimp,aimp,b10max
      real(cp), save :: direction

      if ( .not. l_update_b ) return

      l1m0=map_mlo%ml2i(0,1)

      !-- Now assemble the right hand side and store it in work_LMdist
      call tscheme%set_imex_rhs(work_LMdist, dbdt, 1, n_mloMag_loc, n_r_maxMag)
      call tscheme%set_imex_rhs(ddb, djdt, 1, n_mloMag_loc, n_r_maxMag)

      if ( l_cond_ic ) then
         call tscheme%set_imex_rhs(ddb_ic, dbdt_ic, 1, n_mloMag_loc, n_r_ic_max)
         call tscheme%set_imex_rhs(ddj_ic, djdt_ic, 1, n_mloMag_loc, n_r_ic_max)
      end if

      call solve_counter%start_count()

      if ( (n_imp == 3 .or. n_imp == 4 .or. n_imp == 7) .and. ( l_imp /= 1 ) ) then
         call abortRun('l_imp /= 1 not implemented for this imposed field setup!')
      end if

      do lj=1, n_lo_loc
         l = map_mlo%lj2l(lj)

         if ( l > 0 ) then
            if ( .not. lBmat(lj) ) then
#ifdef WITH_PRECOND_BJ
               call get_bMat(tscheme,l,hdif_B(l),bMat(lj), &
                    &        bMat_fac(:,lj),jMat(lj),jMat_fac(:,lj))
#else
               call get_bMat(tscheme,l,hdif_B(l),bMat(lj),jMat(lj))
#endif
               lBmat(lj)=.true.
            end if

            ! Loop over m corresponding to current l
            do mi=1,map_mlo%n_mi(lj)
               m = map_mlo%milj2m(mi,lj)
               i = map_mlo%milj2i(mi,lj)

               !-------- Magnetic boundary conditions, outer core:
               !         Note: the CMB condition is not correct if we assume free slip
               !         and a conducting mantle (conductance_ma>0).
               if ( l_b_nl_cmb ) then ! finitely conducting mantle
                  rhs1(1,2*mi-1,1) =  real(b_nl_cmb(map_glbl_st%lm2(l,m)))
                  rhs1(1,2*mi,1)   = aimag(b_nl_cmb(map_glbl_st%lm2(l,m)))
                  rhs2(1,2*mi-1,1) =  real(aj_nl_cmb(map_glbl_st%lm2(l,m)))
                  rhs2(1,2*mi,1)   = aimag(aj_nl_cmb(map_glbl_st%lm2(l,m)))
               else
                  rhs1(1,2*mi-1,1) = 0.0_cp
                  rhs1(1,2*mi,1)   = 0.0_cp
                  rhs2(1,2*mi-1,1) = 0.0_cp
                  rhs2(1,2*mi,1)   = 0.0_cp
               end if

               !-------- inner core
               rhs1(n_r_max,2*mi-1,1)=0.0_cp
               rhs1(n_r_max,2*mi,1)  =0.0_cp
               rhs2(n_r_max,2*mi-1,1)=0.0_cp
               rhs2(n_r_max,2*mi,1)  =0.0_cp

               if ( m == 0 ) then   ! Magnetoconvection boundary conditions
                  if ( imagcon /= 0 .and. tmagcon <= time ) then
                     if ( l_LCR ) then
                        call abortRun('LCR not compatible with imposed field!')
                     end if
                     if ( l == 2 .and. imagcon > 0 .and. imagcon  /=  12 ) then
                        rhs2(1,2*mi-1,1)      =bpeaktop
                        rhs2(1,2*mi,1)        =0.0_cp
                        rhs2(n_r_max,2*mi-1,1)=bpeakbot
                        rhs2(n_r_max,2*mi,1)  =0.0_cp
                     else if( l == 1 .and. imagcon == 12 ) then
                        rhs2(1,2*mi-1,1)      =bpeaktop
                        rhs2(1,2*mi,1)        =0.0_cp
                        rhs2(n_r_max,2*mi-1,1)=bpeakbot
                        rhs2(n_r_max,2*mi,1)  =0.0_cp
                     else if( l == 1 .and. imagcon == -1) then
                        rhs1(n_r_max,2*mi-1,1)=bpeakbot
                        rhs1(n_r_max,2*mi,1)  =0.0_cp
                     else if( l == 1 .and. imagcon == -2) then
                        rhs1(1,2*mi-1,1)      =bpeaktop
                        rhs1(1,2*mi,1)        =0.0_cp
                     else if( l == 3 .and. imagcon == -10 ) then
                        rhs2(1,2*mi-1,1)      =bpeaktop
                        rhs2(1,2*mi,1)        =0.0_cp
                        rhs2(n_r_max,2*mi-1,1)=bpeakbot
                        rhs2(n_r_max,2*mi,1)  =0.0_cp
                     end if
                  end if

                  if ( l_curr .and. (mod(l,2) /= 0) ) then ! Current carrying loop around equator of sphere, only odd harmonics

                     if ( l_LCR ) then
                        call abortRun('LCR not compatible with imposed field!')
                     end if

                     !General normalization for spherical harmonics of degree l and order 0
                     yl0_norm=half*sqrt((2*l+1)/pi)

                     !Prefactor for CMB matching condition
                     prefac = real(2*l+1,kind=cp)/real(l*(l+1),kind=cp)


                     bpeaktop=prefac*fac_loop(l)*amp_curr*r_cmb/yl0_norm

                     rhs1(1,2*mi-1,1)=bpeaktop
                     rhs1(1,2*mi,1)  =0.0_cp

                  end if

                  if ( n_imp > 1 .and. l == l_imp ) then
                      if ( l_LCR ) then
                         call abortRun('LCR not compatible with imposed field!')
                      end if
                      ! General normalization for degree l and order 0
                      yl0_norm = half*sqrt((2*l+1)/pi)
                      ! Prefactor for CMB matching condition
                      prefac = real(2*l+1,kind=cp)/real(l*(l+1),kind=cp)

                     if ( n_imp == 2 ) then
                        !  Chose external field coefficient so that amp_imp is
                        !  the amplitude of the external field:
                        bpeaktop=prefac*r_cmb/yl0_norm*amp_imp
                        rhs1(1,2*mi-1,1)=bpeaktop
                        rhs1(1,2*mi,1)  =0.0_cp
                     else if ( n_imp == 3 ) then
                        !  Chose external field coefficient so that amp_imp is
                        !  the amplitude of the external field:
                        bpeaktop=prefac*r_cmb/yl0_norm*amp_imp
                        if ( real(b(l1m0,n_r_cmb)) > 1.0e-9_cp ) &
                        &    direction=real(b(l1m0,n_r_cmb))/abs(real(b(l1m0,n_r_cmb)))
                        rhs1(1,2*mi-1,1)=bpeaktop
                        rhs1(1,2*mi,1)  =0.0_cp
                        rhs1(1,2*mi-1,1)=bpeaktop*direction
                        rhs1(1,2*mi,1)  =0.0_cp
                     else if ( n_imp == 4 ) then
                        !  I have forgotten what this was supposed to do:
                        bpeaktop=three/r_cmb*amp_imp*real(b(l1m0,n_r_cmb))**2
                        rhs1(1,2*mi-1,1)=bpeaktop/real(b(l1m0,n_r_cmb))
                        rhs1(1,2*mi,1)  =0.0_cp

                     else
                        ! Special Heyner feedback functions:
                        ! We don't provide an actual external field strength but rather its
                        ! dependence on the internal dipole via a feedback function ff:
                        !              b10_external=ff(b10_internal)
                        ! where b10_external is the external axial dipole field coeff. and
                        ! b10_internal=b(l1m0,n_r_cmb) is the internal axial dipole field coeff.
                        ! Then rhs1 is always given by:
                        !              rhs1 = (2*l+1)/r_cmb * ff
                        ! because of the special formulation of the CMB matching condition!
                        ! Note that
                        !  B_r_internal(r_cmb) = l*(l+1)/r_cmb**2 * b10_internal*y10_norm*cos(theta)
                        !  B_r_external(r)     = l*(l+1)/r_cmb**2 * b10_external*y10_norm*cos(theta)
                        ! This determines the units of ff=b10_external.
                        ! B itself is given in units sqrt(rho*omega/sigma) so that B**2 is
                        ! the Elsasser number. ff is thus given in units L**2*sqrt(rho*omega/sigma)
                        ! with L=(r_cmb-r_icb). Note that the external dipole field does not depend
                        ! on the radius.
                        ! Here amp_imp provides the relative amplitude so that the maximum of
                        ! the external field is   max(b10_external/b10_internal)=amp_imp

                        ff=0.0_cp
                        if ( n_imp == 7 ) then
                           ! Using a feedback function of the form
                           !    ff= aimp * b10_internal**expo_imp / (cimp+b10_internal**expo_imp)
                           ! Here expo_imp is an input parameter and aimp and cimp are determined
                           ! such that the maximum of b10_external/b10_internal is amp_imp
                           ! and is located at b10_internal=bmax_imp.
                           ! amp_imp and bmax_imp are two more input parameters.
                           ! Note that bmax_imp on input hat the dimensionless unit of the magnetic
                           ! field B and we convert is to dimensionless units for dipole
                           ! coefficients b10 first by multiplying with r_cmb**2
                           ! Other factors like  l*(l+1)*y10_norm*mean(cos(theta)) may also be
                           ! considered, but like r_cmb**2 they are all of order one here.
                           ! Since amp_imp is dimensionless aimp and thus ff has the dimension of
                           ! b10 as required.

                           !  I have forgotten what this was supposed to do:
                           b10max=bmax_imp*r_cmb**2
                           cimp=b10max**expo_imp / (expo_imp-1)
                           aimp=amp_imp*cimp*expo_imp / &
                           &    (cimp*(expo_imp-1))**((expo_imp-1)/expo_imp)
                           ff=aimp*real(b(l1m0,n_r_cmb))**expo_imp/ &
                           &  (cimp+real(b(l1m0,n_r_cmb))**expo_imp)
                        end if
                        rhs1(1,2*mi-1,1)=(2*l+1)/r_cmb*ff
                        rhs1(1,2*mi,1)  =0.0_cp

                     end if
                  end if ! n_imp > 1

               end if ! m == 0 special BC when imposed field

               do nR=2,n_r_max-1
                  if ( nR<=n_r_LCR ) then
                     rhs1(nR,2*mi-1,1)=0.0_cp
                     rhs1(nR,2*mi,1)  =0.0_cp
                     rhs2(nR,2*mi-1,1)=0.0_cp
                     rhs2(nR,2*mi,1)  =0.0_cp
                  else
                     rhs1(nR,2*mi-1,1)= real(work_LMdist(i,nR))
                     rhs1(nR,2*mi,1)  =aimag(work_LMdist(i,nR))
                     rhs2(nR,2*mi-1,1)= real(ddb(i,nR)) ! ddb is used as a work array here
                     rhs2(nR,2*mi,1)  =aimag(ddb(i,nR))
                  end if
               end do

               !-- Overwrite RHS when perfect conductor
               if ( ktopb == 2 ) then
                  rhs1(2,2*mi-1,1) = 0.0_cp
                  rhs1(2,2*mi,1)   = 0.0_cp
               end if
               if ( kbotb == 2 ) then
                  rhs1(n_r_max-1,2*mi-1,1) = 0.0_cp
                  rhs1(n_r_max-1,2*mi,1)   = 0.0_cp
               end if

               !-- Magnetic boundary conditions, inner core for radial derivatives
               !         of poloidal and toroidal magnetic potentials:
               if ( l_cond_ic ) then    ! inner core
                  rhs1(n_r_max+1,2*mi-1,1)=0.0_cp
                  rhs1(n_r_max+1,2*mi,1)  =0.0_cp
                  if ( l_b_nl_icb ) then
                     !@> TODO one can likely have aj_nl_icb in the right space here
                     rhs2(n_r_max+1,2*mi-1,1)=real( &
                     &                        aj_nl_icb(map_glbl_st%lm2(l,m)) )
                     rhs2(n_r_max+1,2*mi,1)  =aimag( &
                     &                        aj_nl_icb(map_glbl_st%lm2(l,m)) )
                  else
                     rhs2(n_r_max+1,2*mi-1,1)=0.0_cp
                     rhs2(n_r_max+1,2*mi,1)  =0.0_cp
                  end if


                  do nR=2,n_r_ic_max
                     rhs1(n_r_max+nR,2*mi-1,1)= real(ddb_ic(i,nR)) ! ddb_ic as work array
                     rhs1(n_r_max+nR,2*mi,1)  =aimag(ddb_ic(i,nR))
                     rhs2(n_r_max+nR,2*mi-1,1)= real(ddj_ic(i,nR))
                     rhs2(n_r_max+nR,2*mi,1)  =aimag(ddj_ic(i,nR))
                  end do
               end if

#ifdef WITH_PRECOND_BJ
               do nR=1,n_r_tot
                  rhs1(nR,2*mi-1,1)=rhs1(nR,2*mi-1,1)*bMat_fac(nR,lj)
                  rhs1(nR,2*mi,1)  =rhs1(nR,2*mi,1)*bMat_fac(nR,lj)
                  rhs2(nR,2*mi-1,1)=rhs2(nR,2*mi-1,1)*jMat_fac(nR,lj)
                  rhs2(nR,2*mi,1)  =rhs2(nR,2*mi,1)*jMat_fac(nR,lj)
               end do
#endif
            end do    ! loop over m

            call bMat(lj)%solve(rhs1(:,1:2*map_mlo%n_mi(lj),1), &
                 &                 2*map_mlo%n_mi(lj))
            call jMat(lj)%solve(rhs2(:,1:2*map_mlo%n_mi(lj),1), &
                 &                 2*map_mlo%n_mi(lj))

            if ( lRmsNext .and. tscheme%istage == 1 ) then ! Store old b,aj
               do nR=1,n_r_max
                  do mi=1,map_mlo%n_mi(lj)
                     i = map_mlo%milj2i(mi,lj)
                     workA(i,nR)= b(i,nR)
                     workB(i,nR)=aj(i,nR)
                  end do
               end do
            end if

            !----- Update magnetic field in cheb space:
            do mi=1,map_mlo%n_mi(lj)
               m = map_mlo%milj2m(mi,lj)
               i = map_mlo%milj2i(mi,lj)

               if ( m > 0 ) then ! Non-axisymmetric modes
                  do n_r_out=1,rscheme_oc%n_max  ! outer core
                     b(i,n_r_out) =cmplx(rhs1(n_r_out,2*mi-1,1), &
                     &                     rhs1(n_r_out,2*mi,1),cp)
                     aj(i,n_r_out)=cmplx(rhs2(n_r_out,2*mi-1,1), &
                     &                     rhs2(n_r_out,2*mi,1),cp)
                  end do
                  if ( l_cond_ic ) then   ! inner core
                     do n_r_out=1,n_cheb_ic_max
                        b_ic(i,n_r_out) =cmplx(rhs1(n_r_max+n_r_out,2*mi-1,  &
                        &                        1),rhs1(n_r_max+n_r_out,&
                        &                        2*mi,1),cp)
                        aj_ic(i,n_r_out)=cmplx(rhs2(n_r_max+n_r_out,2*mi-1,  &
                        &                        1),rhs2(n_r_max+n_r_out,&
                        &                        2*mi,1),cp)
                     end do
                  end if
               else ! Axisymmetric modes
                  do n_r_out=1,rscheme_oc%n_max   ! outer core
                     b(i,n_r_out) = cmplx(rhs1(n_r_out,2*mi-1,1), &
                     &                      0.0_cp,kind=cp)
                     aj(i,n_r_out)= cmplx(rhs2(n_r_out,2*mi-1,1), &
                     &                      0.0_cp,kind=cp)
                  end do
                  if ( l_cond_ic ) then    ! inner core
                     do n_r_out=1,n_cheb_ic_max
                        b_ic(i,n_r_out)= cmplx(rhs1(n_r_max+n_r_out, &
                        &                   2*mi-1,1),0.0_cp,kind=cp)
                        aj_ic(i,n_r_out)= cmplx(rhs2(n_r_max+n_r_out,2*mi-1,&
                        &                   1),0.0_cp,kind=cp)
                     end do
                  end if
               end if

            end do

         end if   ! l > 0

      end do      ! end of do loop over l'ss

      call solve_counter%stop_count(l_increment=.false.)

      !-- Set cheb modes > rscheme_oc%n_max to zero (dealiazing)
      !   for inner core modes > 2*n_cheb_ic_max = 0
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         do i=1,n_mloMag_loc
            b(i,n_r_out) =zero
            aj(i,n_r_out)=zero
         end do
      end do

      if ( l_cond_ic ) then
         do n_r_out=n_cheb_ic_max+1,n_r_ic_max
            do i=1,n_mloMag_loc
               b_ic(i,n_r_out) =zero
               aj_ic(i,n_r_out)=zero
            end do
         end do
      end if

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dbdt, 1, n_mloMag_loc, n_r_maxMag)
      call tscheme%rotate_imex(djdt, 1, n_mloMag_loc, n_r_maxMag)
      if ( l_cond_ic ) then
         call tscheme%rotate_imex(dbdt_ic, 1, n_mloMag_loc, n_r_ic_max)
         call tscheme%rotate_imex(djdt_ic, 1, n_mloMag_loc, n_r_ic_max)
      end if


      !-- Get implicit terms
      if ( tscheme%istage == tscheme%nstages ) then
         call get_mag_rhs_imp(b, db, ddb, aj, dj, ddj, dbdt, djdt, tscheme, 1, &
              &               tscheme%l_imp_calc_rhs(1), lRmsNext,             &
              &               l_in_cheb_space=.true.)

         if ( l_cond_ic ) then
            call get_mag_ic_rhs_imp(b_ic, db_ic, ddb_ic, aj_ic, dj_ic, ddj_ic,     &
                 &                  dbdt_ic, djdt_ic, 1, tscheme%l_imp_calc_rhs(1),&
                 &                  l_in_cheb_space=.true.)
         end if

      else
         call get_mag_rhs_imp(b, db, ddb, aj, dj, ddj, dbdt, djdt, tscheme, &
              &               tscheme%istage+1,                             &
              &               tscheme%l_imp_calc_rhs(tscheme%istage+1),     &
              &               lRmsNext, l_in_cheb_space=.true.)

         if ( l_cond_ic ) then
            call get_mag_ic_rhs_imp(b_ic, db_ic, ddb_ic, aj_ic, dj_ic, ddj_ic,  &
                 &                  dbdt_ic, djdt_ic, tscheme%istage+1,         &
                 &                  tscheme%l_imp_calc_rhs(tscheme%istage+1),   &
                 &                  l_in_cheb_space=.true.)
         end if

      end if

   end subroutine updateB
!-----------------------------------------------------------------------------
   subroutine finish_exp_mag_ic(b_ic, aj_ic, omega_ic, db_exp_last, dj_exp_last)

      !-- Input variables
      real(cp),    intent(in) :: omega_ic
      complex(cp), intent(in) :: aj_ic(n_mloMag_loc,n_r_ic_max)
      complex(cp), intent(in) :: b_ic(n_mloMag_loc,n_r_ic_max)

      !-- Output variables
      complex(cp), intent(inout) :: db_exp_last(n_mloMag_loc,n_r_ic_max)
      complex(cp), intent(inout) :: dj_exp_last(n_mloMag_loc,n_r_ic_max)

      !-- Local variables
      complex(cp) :: fac
      integer :: n_r, lm, l, m

      if ( omega_ic /= 0.0_cp .and. l_rot_ic .and. l_mag_nl ) then
         !$omp parallel do default(shared) private(lm,n_r,fac,l,m) collapse(2)
         do n_r=2,n_r_ic_max
            do lm=1,n_mloMag_loc
               l = map_mlo%i2l(lm)
               m = map_mlo%i2m(lm)

               fac = -omega_ic*or2(n_r_max)*cmplx(0.0_cp,real(m,cp),cp)* &
               &      real(l*(l+1),cp)
               db_exp_last(lm,n_r)=fac* b_ic(lm,n_r)
               dj_exp_last(lm,n_r)=fac*aj_ic(lm,n_r)
            end do
         end do
         !$omp end parallel do
      else
         !$omp parallel do default(shared) private(lm,n_r,fac) collapse(2)
         do n_r=2,n_r_ic_max
            do lm=1,n_mloMag_loc
               db_exp_last(lm,n_r)=zero
               dj_exp_last(lm,n_r)=zero
            end do
         end do
         !$omp end parallel do
      end if

   end subroutine finish_exp_mag_ic
!-----------------------------------------------------------------------------
   subroutine finish_exp_mag(dVxBhLM, dj_exp_last)

      !-- Input variables
      complex(cp), intent(inout) :: dVxBhLM(n_mloMag_loc,n_r_maxMag)

      !-- Output variables
      complex(cp), intent(inout) :: dj_exp_last(n_mloMag_loc,n_r_maxMag)

      !-- Local variables
      integer :: n_r, lm, start_lm, stop_lm

      !$omp parallel default(shared) private(start_lm,stop_lm)
      start_lm=1; stop_lm=n_mloMag_loc
      call get_openmp_blocks(start_lm, stop_lm)

      call get_dr( dVxBhLM, work_LMdist, n_mloMag_loc, start_lm, &
           &       stop_lm, n_r_max, rscheme_oc, nocopy=.true. )
      !$omp barrier

      !$omp do private(n_r,lm)
      do n_r=1,n_r_max
         do lm=1,n_mloMag_loc
            dj_exp_last(lm,n_r)=dj_exp_last(lm,n_r)+or2(n_r)*work_LMdist(lm,n_r)
         end do
      end do
      !$omp end do
      !$omp end parallel

   end subroutine finish_exp_mag
!-----------------------------------------------------------------------------
   subroutine get_mag_ic_rhs_imp(b_ic, db_ic, ddb_ic, aj_ic, dj_ic, ddj_ic,  &
              &                  dbdt_ic, djdt_ic, istage, l_calc_lin,       &
              &                  l_in_cheb_space)

      !-- Input variables
      integer,             intent(in) :: istage
      logical,             intent(in) :: l_calc_lin
      logical, optional,   intent(in) :: l_in_cheb_space

      !-- Output variables
      type(type_tarray), intent(inout) :: dbdt_ic
      type(type_tarray), intent(inout) :: djdt_ic
      complex(cp),       intent(inout) :: b_ic(n_mloMag_loc,n_r_ic_max)
      complex(cp),       intent(inout) :: aj_ic(n_mloMag_loc,n_r_ic_max)
      complex(cp),       intent(out) :: db_ic(n_mloMag_loc,n_r_ic_max)
      complex(cp),       intent(out) :: ddb_ic(n_mloMag_loc,n_r_ic_max)
      complex(cp),       intent(out) :: dj_ic(n_mloMag_loc,n_r_ic_max)
      complex(cp),       intent(out) :: ddj_ic(n_mloMag_loc,n_r_ic_max)

      !-- Local variables
      complex(cp) :: tmp(n_mloMag_loc,n_r_ic_max)
      real(cp) :: dL
      logical :: l_in_cheb
      integer :: l1, n_r, lm, start_lm, stop_lm

      if ( present(l_in_cheb_space) ) then
         l_in_cheb = l_in_cheb_space
      else
         l_in_cheb = .false.
      end if

      !$omp parallel default(shared)  private(start_lm, stop_lm)
      start_lm=1; stop_lm=n_mloMag_loc
      call get_openmp_blocks(start_lm,stop_lm)

      !$omp single
      call dct_counter%start_count()
      !$omp end single
      if ( l_in_cheb ) call chebt_ic%costf1( b_ic, n_mloMag_loc, &
                            &                start_lm, stop_lm, work_ic_LMdist)
      call get_ddr_even( b_ic,db_ic,ddb_ic, n_mloMag_loc, &
           &             start_lm,stop_lm, n_r_ic_max,n_cheb_ic_max, dr_fac_ic,&
           &             work_ic_LMdist, tmp, chebt_ic, chebt_ic_even )
      if ( l_in_cheb ) call chebt_ic%costf1( aj_ic, n_mloMag_loc, &
                            &                start_lm, stop_lm, work_ic_LMdist)
      call get_ddr_even( aj_ic,dj_ic,ddj_ic, n_mloMag_loc,  &
           &             start_lm, stop_lm, n_r_ic_max,n_cheb_ic_max, dr_fac_ic,  &
           &             work_ic_LMdist, tmp, chebt_ic, chebt_ic_even )
      !$omp barrier
      !$omp single
      call dct_counter%stop_count(l_increment=.false.)
      !$omp end single

      if ( istage == 1 ) then
         !$omp do private(n_r,lm,l1,dL) collapse(2)
         do n_r=1,n_r_ic_max
            do lm=1,n_mloMag_loc
               l1 = map_mlo%i2l(lm)
               dL = real(l1*(l1+1),cp)
               dbdt_ic%old(lm,n_r,istage)=dL*or2(n_r_max)* b_ic(lm,n_r)
               djdt_ic%old(lm,n_r,istage)=dL*or2(n_r_max)*aj_ic(lm,n_r)
            end do
         end do
         !$omp end do
      end if

      if ( l_calc_lin ) then
         !$omp do private(n_r,lm,l1,dL) collapse(2)
         do n_r=2,n_r_ic_max-1
            do lm=1,n_mloMag_loc
               l1 = map_mlo%i2l(lm)
               dL = real(l1*(l1+1),cp)
               dbdt_ic%impl(lm,n_r,istage)=opm*O_sr*dL*or2(n_r_max)* (  &
               &                                      ddb_ic(lm,n_r) +  &
               &         two*real(l1+1,cp)*O_r_ic(n_r)*db_ic(lm,n_r) )
               djdt_ic%impl(lm,n_r,istage)=opm*O_sr*dL*or2(n_r_max) *  (&
               &                                      ddj_ic(lm,n_r) +  &
               &         two*real(l1+1,cp)*O_r_ic(n_r)*dj_ic(lm,n_r) )
            end do
         end do
         !$omp end do
         n_r=n_r_ic_max
         !$omp do private(lm,l1,dL)
         do lm=1,n_mloMag_loc
            l1 = map_mlo%i2l(lm)
            dL = real(l1*(l1+1),cp)
            dbdt_ic%impl(lm,n_r,istage)=opm*O_sr*dL*or2(n_r_max) *  &
            &                           (one+two*real(l1+1,cp))*ddb_ic(lm,n_r)
            djdt_ic%impl(lm,n_r,istage)=opm*O_sr*dL* or2(n_r_max) *  &
            &                           (one+two*real(l1+1,cp))*ddj_ic(lm,n_r)
         end do
         !$omp end do
      end if

      !$omp end parallel

   end subroutine get_mag_ic_rhs_imp
!-----------------------------------------------------------------------------
   subroutine assemble_mag(b, db, ddb, aj, dj, ddj, b_ic, db_ic, ddb_ic, aj_ic, &
              &            dj_ic, ddj_ic, dbdt, djdt, dbdt_ic, djdt_ic,         &
              &            lRmsNext, tscheme)

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lRmsNext

      !-- Output variables
      type(type_tarray), intent(inout) :: dbdt, djdt, dbdt_ic, djdt_ic
      complex(cp),       intent(inout) :: b(n_mloMag_loc,n_r_maxMag)
      complex(cp),       intent(inout) :: aj(n_mloMag_loc,n_r_maxMag)
      complex(cp),       intent(inout) :: b_ic(n_mloMag_loc,n_r_ic_max)
      complex(cp),       intent(inout) :: aj_ic(n_mloMag_loc,n_r_ic_max)
      complex(cp),       intent(out) :: db(n_mloMag_loc,n_r_maxMag)
      complex(cp),       intent(out) :: dj(n_mloMag_loc,n_r_maxMag)
      complex(cp),       intent(out) :: ddj(n_mloMag_loc,n_r_maxMag)
      complex(cp),       intent(out) :: ddb(n_mloMag_loc,n_r_maxMag)
      complex(cp),       intent(out) :: db_ic(n_mloMag_loc,n_r_ic_max)
      complex(cp),       intent(out) :: dj_ic(n_mloMag_loc,n_r_ic_max)
      complex(cp),       intent(out) :: ddj_ic(n_mloMag_loc,n_r_ic_max)
      complex(cp),       intent(out) :: ddb_ic(n_mloMag_loc,n_r_ic_max)

      !-- Local variables
      complex(cp) :: val_bot
      real(cp) :: fac_top, fac_bot
      integer :: n_r, lm, l, m
      real(cp) :: dL

      if ( l_b_nl_cmb .or. l_b_nl_icb ) then
         call abortRun('Non linear magnetic BCs not implemented at assembly stage!')
      end if

      !-- Assemble IMEX using ddb and ddj as a work array
      call tscheme%assemble_imex(ddb, dbdt, 1, n_mloMag_loc, n_r_maxMag)
      call tscheme%assemble_imex(ddj, djdt, 1, n_mloMag_loc, n_r_maxMag)
      if ( l_cond_ic ) then
         call tscheme%assemble_imex(ddb_ic, dbdt_ic, 1, n_mloMag_loc, n_r_ic_max)
         call tscheme%assemble_imex(ddj_ic, djdt_ic, 1, n_mloMag_loc, n_r_ic_max)
      end if

      !-- Now get the toroidal potential from the assembly
      !$omp parallel default(shared)
      !$omp do private(n_r,lm,l,m,dL)
      do n_r=2,n_r_max-1
         do lm=1,n_mloMag_loc
            l = map_mlo%i2l(lm)
            m = map_mlo%i2m(lm)
            if ( l == 0 ) cycle
            dL = real(l*(l+1),cp)
            if ( m == 0 ) then
               b(lm,n_r)  = r(n_r)*r(n_r)/dL*cmplx(real(ddb(lm,n_r)),0.0_cp,cp)
               aj(lm,n_r) = r(n_r)*r(n_r)/dL*cmplx(real(ddj(lm,n_r)),0.0_cp,cp)
            else
               b(lm,n_r)  = r(n_r)*r(n_r)/dL*ddb(lm,n_r)
               aj(lm,n_r) = r(n_r)*r(n_r)/dL*ddj(lm,n_r)
            end if
         end do
      end do
      !$omp end do


      if ( l_cond_ic ) then
         !$omp do private(n_r,lm,l,m,dL)
         do n_r=2,n_r_ic_max
            do lm=1,n_mloMag_loc
               l = map_mlo%i2l(lm)
               m = map_mlo%i2m(lm)
               if ( l == 0 ) cycle
               dL = real(l*(l+1),cp)
               if ( m == 0 ) then
                  b_ic(lm,n_r)  = r(n_r_max)*r(n_r_max)/dL*cmplx(real(ddb_ic(lm,n_r)),0.0_cp,cp)
                  aj_ic(lm,n_r) = r(n_r_max)*r(n_r_max)/dL*cmplx(real(ddj_ic(lm,n_r)),0.0_cp,cp)
               else
                  b_ic(lm,n_r)  = r(n_r_max)*r(n_r_max)/dL*ddb_ic(lm,n_r)
                  aj_ic(lm,n_r) = r(n_r_max)*r(n_r_max)/dL*ddj_ic(lm,n_r)
               end if
            end do
         end do
         !$omp end do
      end if

      !-- Now handle boundary conditions !
      if ( imagcon /= 0 ) call abortRun('imagcon/=0 not implemented with assembly stage!')
      if ( conductance_ma /=0 ) call abortRun('conducting ma not implemented here!')

      !-- If conducting inner core then the solution at ICB should already be fine
      if ( l_full_sphere ) then
         if ( ktopb == 1 ) then
            !$omp do private(lm,l,fac_top)
            do lm=1,n_mloMag_loc
               l = map_mlo%i2l(lm)
               if ( l==0 ) cycle
               fac_top=real(l,cp)*or1(1)
               call rscheme_oc%robin_bc(one, fac_top, zero, 0.0_cp, one, zero, b(lm,:))
               aj(lm,1)      =zero
               aj(lm,n_r_max)=zero
            end do
            !$omp end do
         else if (ktopb == 4 ) then
            !$omp do private(lm,l)
            do lm=1,n_mloMag_loc
               l = map_mlo%i2l(lm)
               if ( l==0 ) cycle
               call rscheme_oc%robin_bc(one, 0.0_cp, zero, 0.0_cp, one, zero, b(lm,:))
               aj(lm,1)      =zero
               aj(lm,n_r_max)=zero
            end do
            !$omp end do
         else
            call abortRun('Not implemented yet!')
         end if
      else ! spherical shell
         if ( kbotb == 3 ) then ! Conducting inner core
            if ( ktopb==1 ) then ! Vacuum outside + cond. I. C.
               !$omp do private(lm,l,fac_top,fac_bot,val_bot)
               do lm=1,n_mloMag_loc
                  l = map_mlo%i2l(lm)
                  if ( l==0 ) cycle
                  fac_top=real(l,cp)*or1(1)
                  fac_bot=-dr_top_ic(1)-real(l+1,cp)*or1(n_r_max)
                  val_bot = zero
                  do n_r=2,n_r_ic_max
                     val_bot = val_bot+dr_top_ic(n_r)*b_ic(lm,n_r)
                  end do
                  call rscheme_oc%robin_bc(one, fac_top, zero, one, fac_bot, val_bot, &
                       &                   b(lm,:))
                  b_ic(lm,1)=b(lm,n_r_max)

                  val_bot = zero
                  do n_r=2,n_r_ic_max
                     val_bot = val_bot+dr_top_ic(n_r)*aj_ic(lm,n_r)
                  end do
                  call rscheme_oc%robin_bc(0.0_cp, one, zero, sigma_ratio, fac_bot, &
                       &                   val_bot, aj(lm,:))
                  aj_ic(lm,1)=aj(lm,n_r_max)
               end do
               !$omp end do
            else if ( ktopb == 4 ) then ! Pseudo-Vacuum outside + cond. I. C.
               !$omp do private(lm,l,fac_bot,val_bot)
               do lm=1,n_mloMag_loc
                  l = map_mlo%i2l(lm)
                  if ( l==0 ) cycle
                  fac_bot=-dr_top_ic(1)+real(l+1,cp)*or1(n_r_max)
                  val_bot = zero
                  do n_r=2,n_r_ic_max
                     val_bot = val_bot+dr_top_ic(n_r)*b_ic(lm,n_r)
                  end do
                  call rscheme_oc%robin_bc(one, 0.0_cp, zero, one, fac_bot, val_bot, &
                       &                   b(lm,:))
                  b_ic(lm,1)=b(lm,n_r_max)

                  val_bot = zero
                  do n_r=2,n_r_ic_max
                     val_bot = val_bot+dr_top_ic(n_r)*aj_ic(lm,n_r)
                  end do
                  call rscheme_oc%robin_bc(0.0_cp, one, zero, sigma_ratio, fac_bot, &
                       &                   val_bot, aj(lm,:))
                  aj_ic(lm,1)=aj(lm,n_r_max)
               end do
               !$omp end do
            else
               call abortRun('Not implemented yet!')
            end if
         else if ( kbotb == 4 ) then
            if ( ktopb==1 ) then
               !$omp do private(lm,l,fac_top)
               do lm=1,n_mloMag_loc
                  l = map_mlo%i2l(lm)
                  if ( l==0 ) cycle
                  fac_top=real(l,cp)*or1(1)
                  call rscheme_oc%robin_bc(one, fac_top, zero, one, 0.0_cp, zero, b(lm,:))
                  aj(lm,1)      =zero
                  aj(lm,n_r_max)=zero
               end do
               !$omp end do
            else if ( ktopb == 4 ) then
               !$omp do private(lm,l)
               do lm=1,n_mloMag_loc
                  l = map_mlo%i2l(lm)
                  if ( l==0 ) cycle
                  call rscheme_oc%robin_bc(one, 0.0_cp, zero, one, 0.0_cp, zero, b(lm,:))
                  aj(lm,1)      =zero
                  aj(lm,n_r_max)=zero
               end do
               !$omp end do
            else
               call abortRun('Not implemented yet!')
            end if
         else if ( kbotb == 1 ) then ! Insulating inner core
            if ( ktopb==1 ) then ! Vacuum on both sides
               !$omp do private(lm,l,fac_bot,fac_top)
               do lm=1,n_mloMag_loc
                  l = map_mlo%i2l(lm)
                  if ( l==0 ) cycle
                  fac_bot=-real(l+1,cp)*or1(n_r_max)
                  fac_top=real(l,cp)*or1(1)
                  call rscheme_oc%robin_bc(one, fac_top, zero, one, fac_bot, zero, b(lm,:))
                  aj(lm,1)      =zero
                  aj(lm,n_r_max)=zero
               end do
               !$omp end do
            else if ( ktopb == 4 ) then ! Pseudo-vacuum on the outer boundary
               !$omp do private(lm,l,fac_bot)
               do lm=1,n_mloMag_loc
                  l = map_mlo%i2l(lm)
                  if ( l==0 ) cycle
                  fac_bot=-real(l+1,cp)*or1(n_r_max)
                  call rscheme_oc%robin_bc(one, 0.0_cp, zero, one, fac_bot, zero, b(lm,:))
                  aj(lm,1)      =zero
                  aj(lm,n_r_max)=zero
               end do
               !$omp end do
            end if
         else
            call abortRun('Combination not implemented yet!')
         end if
      end if
      !$omp end parallel

      !-- Finally compute the required implicit stage if needed
      call get_mag_rhs_imp(b, db, ddb, aj, dj, ddj, dbdt, djdt, tscheme, 1, &
           &               tscheme%l_imp_calc_rhs(1), lRmsNext, .false.)

      if ( l_cond_ic ) then
         call get_mag_ic_rhs_imp(b_ic, db_ic, ddb_ic, aj_ic, dj_ic, ddj_ic,     &
              &                  dbdt_ic, djdt_ic, 1, tscheme%l_imp_calc_rhs(1),&
              &                  .false.)
      end if

   end subroutine assemble_mag
!-----------------------------------------------------------------------------
   subroutine get_mag_rhs_imp(b, db, ddb, aj, dj, ddj, dbdt, djdt, tscheme, &
              &               istage, l_calc_lin, lRmsNext, l_in_cheb_space)

      !-- Input variables
      integer,             intent(in) :: istage
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lRmsNext
      logical,             intent(in) :: l_calc_lin
      logical, optional,   intent(in) :: l_in_cheb_space

      !-- Output variables
      type(type_tarray), intent(inout) :: dbdt
      type(type_tarray), intent(inout) :: djdt
      complex(cp),       intent(inout) :: b(n_mloMag_loc,n_r_max)
      complex(cp),       intent(inout) :: aj(n_mloMag_loc,n_r_max)
      complex(cp),       intent(out) :: db(n_mloMag_loc,n_r_max)
      complex(cp),       intent(out) :: dj(n_mloMag_loc,n_r_max)
      complex(cp),       intent(out) :: ddj(n_mloMag_loc,n_r_max)
      complex(cp),       intent(out) :: ddb(n_mloMag_loc,n_r_max)

      !-- Local variables
      real(cp) :: dL
      logical :: l_in_cheb
      integer :: n_r_top, n_r_bot, n_r, lm, l1, start_lm, stop_lm

      if ( present(l_in_cheb_space) ) then
         l_in_cheb = l_in_cheb_space
      else
         l_in_cheb = .false.
      end if

      !$omp parallel default(shared)  private(start_lm, stop_lm)
      start_lm=1; stop_lm=n_mloMag_loc
      call get_openmp_blocks(start_lm,stop_lm)

      !$omp single
      call dct_counter%start_count()
      !$omp end single
      call get_ddr(b,db,ddb,n_mloMag_loc,start_lm,stop_lm,n_r_max,rscheme_oc,  &
           &       l_dct_in=.not. l_in_cheb)
      if ( l_in_cheb ) call rscheme_oc%costf1(b,n_mloMag_loc,start_lm,stop_lm)
      call get_ddr(aj,dj,ddj,n_mloMag_loc,start_lm,stop_lm,n_r_max,rscheme_oc, &
           &       l_dct_in=.not. l_in_cheb)
      if ( l_in_cheb ) call rscheme_oc%costf1(aj,n_mloMag_loc,start_lm,stop_lm)
      !$omp barrier
      !$omp single
      call dct_counter%stop_count(l_increment=.false.)
      !$omp end single

      if ( l_LCR ) then
         !$omp do private(n_r,lm,l1)
         do n_r=n_r_cmb,n_r_icb-1
            if ( n_r<=n_r_LCR ) then
               do lm=1,n_mloMag_loc
                  l1 = map_mlo%i2l(lm)

                  b(lm,n_r)=(r(n_r_LCR)/r(n_r))**real(l1,cp)*b(lm,n_r_LCR)
                  db(lm,n_r)=-real(l1,cp)*(r(n_r_LCR))**real(l1,cp)/  &
                  &          (r(n_r))**(real(l1,cp)+1)*b(lm,n_r_LCR)
                  ddb(lm,n_r)=real(l1,cp)*real(l1+1,cp)*(r(n_r_LCR))**real(l1,cp)/ &
                  &           (r(n_r))**(real(l1,cp)+2)*b(lm,n_r_LCR)
                  aj(lm,n_r) =zero
                  dj(lm,n_r) =zero
                  ddj(lm,n_r)=zero
               end do
            end if
         end do
         !$omp end do
      end if

      if ( istage == 1 ) then
         !$omp do private(n_r,lm,l1,dL)
         do n_r=1,n_r_max
            do lm=1,n_mloMag_loc
               l1 = map_mlo%i2l(lm)
               dL = real(l1*(l1+1),cp)
               dbdt%old(lm,n_r,istage)=dL*or2(n_r)* b(lm,n_r)
               djdt%old(lm,n_r,istage)=dL*or2(n_r)*aj(lm,n_r)
            end do
         end do
         !$omp end do
      end if

      if ( l_calc_lin .or. (tscheme%istage==tscheme%nstages .and. lRmsNext)) then
         if ( lRmsNext ) then
            n_r_top=n_r_cmb
            n_r_bot=n_r_icb
         else
            n_r_top=n_r_cmb+1
            n_r_bot=n_r_icb-1
         end if

         !$omp do private(n_r,lm,l1,dtP,dtT,dL)
         do n_r=n_r_top,n_r_bot
            do lm=1,n_mloMag_loc
               l1 = map_mlo%i2l(lm)
               dL = real(l1*(l1+1),cp)
               dbdt%impl(lm,n_r,istage)=opm*lambda(n_r)*hdif_B(l1)*     &
               &                    dL*or2(n_r)*(ddb(lm,n_r)-dL*or2(n_r)*b(lm,n_r) )
               djdt%impl(lm,n_r,istage)= opm*lambda(n_r)*hdif_B(l1)*           &
               &                    dL*or2(n_r)*( ddj(lm,n_r)+dLlambda(n_r)*   &
               &                    dj(lm,n_r)-dL*or2(n_r)*aj(lm,n_r) )
               if ( lRmsNext .and. tscheme%istage == tscheme%nstages ) then
                  dtP(lm)=dL*or2(n_r)/tscheme%dt(1) * (  b(lm,n_r)-workA(lm,n_r) )
                  dtT(lm)=dL*or2(n_r)/tscheme%dt(1) * ( aj(lm,n_r)-workB(lm,n_r) )
               end if
            end do
            if ( lRmsNext .and. tscheme%istage == tscheme%nstages ) then
               call hInt2PolLM(dtP, 1, n_mloMag_loc, n_r, dtBPolLMr(:,n_r), &
                    &          dtBPol2hInt(:,n_r))
               call hInt2TorLM(dtT, 1, n_mloMag_loc, n_r, dtBTor2hInt(:,n_r))
            end if
         end do
         !$omp end do

      end if
      !$omp end parallel

   end subroutine get_mag_rhs_imp
!-----------------------------------------------------------------------------
#ifdef WITH_PRECOND_BJ
   subroutine get_bMat(tscheme,l,hdif,bMat,bMat_fac,jMat,jMat_fac)
#else
   subroutine get_bMat(tscheme,l,hdif,bMat,jMat)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matrices
      !  bmat(i,j) and ajmat for the dynamo equations.
      !

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme        ! time step
      integer,             intent(in) :: l
      real(cp),            intent(in) :: hdif

      !-- Output variables:
      class(type_realmat), intent(inout) :: bMat
      class(type_realmat), intent(inout) :: jMat
#ifdef WITH_PRECOND_BJ
      real(cp), intent(out) :: bMat_fac(n_r_totMag),jMat_fac(n_r_totMag)
#endif

      !-- local variables:
      integer :: nR,nCheb,nR_out,nRall
      integer :: info
      real(cp) :: l_P_1
      real(cp) :: dLh
      real(cp) :: rRatio
      real(cp) :: datJmat(n_r_tot,n_r_tot)
      real(cp) :: datBmat(n_r_tot,n_r_tot)

      nRall=n_r_max
      if ( l_cond_ic ) nRall=nRall+n_r_ic_max
      dLh=real(l*(l+1),kind=cp)

      !-- matrices depend on degree l but not on order m,
      !   we thus have to construct bmat and ajmat for each l:

      l_P_1=real(l+1,kind=cp)

      do nR_out=1,n_r_max
         do nR=2,n_r_max-1
            datBmat(nR,nR_out)=                       rscheme_oc%rnorm * (  &
            &                 dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out) -      &
            &    tscheme%wimp_lin(1)*opm*lambda(nR)*hdif*dLh*or2(nR) * (    &
            &                                rscheme_oc%d2rMat(nR,nR_out) - &
            &                      dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out) ) )

            datJmat(nR,nR_out)=                       rscheme_oc%rnorm * (  &
            &                 dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out) -      &
            &   tscheme%wimp_lin(1)*opm*lambda(nR)*hdif*dLh*or2(nR) * (     &
            &                                rscheme_oc%d2rMat(nR,nR_out) + &
            &                    dLlambda(nR)*rscheme_oc%drMat(nR,nR_out) - &
            &                      dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out) ) )
         end do
      end do

      if  ( l_LCR ) then
         do nR=2,n_r_max-1
            if ( nR<=n_r_LCR ) then
               do nR_out=1,n_r_max
                  datBmat(nR,nR_out)= rscheme_oc%rnorm*(                    &
                  &                           rscheme_oc%drMat(nR,nR_out) + &
                  &    real(l,kind=cp)*or1(nR)*rscheme_oc%rMat(nR,nR_out) )

                  datJmat(nR,nR_out)= rscheme_oc%rnorm*rscheme_oc%rMat(nR,nR_out)
               end do
            end if
         end do
      end if

      !----- boundary conditions for outer core field:
      if ( ktopb == 1 .or. ktopb == 3 ) then
         !-------- at CMB (nR=1):
         !         the internal poloidal field should fit a potential
         !         field (matrix bmat) and the toroidal field has to
         !         vanish (matrix ajmat).

         datBmat(1,1:n_r_max)=    rscheme_oc%rnorm * (   &
         &                      rscheme_oc%drMat(1,:) +  &
         &     real(l,cp)*or1(1)*rscheme_oc%rMat(1,:) +  &
         &                     conductance_ma* (         &
         &                     rscheme_oc%d2rMat(1,:) -  &
         &            dLh*or2(1)*rscheme_oc%rMat(1,:) ) )

         datJmat(1,1:n_r_max)=    rscheme_oc%rnorm * (   &
         &                       rscheme_oc%rMat(1,:) +  &
         &       conductance_ma*rscheme_oc%drMat(1,:) )
      else if ( ktopb == 2 ) then
         !----- perfect conductor
         !      see Glatzmaier, JCP 55, 461-484 (1984)
         ! the (extra) condition Br=0 on Bpol is imposed just
         ! below the boundary
         datBmat(1,1:n_r_max)=rscheme_oc%rnorm*  rscheme_oc%rMat(1,:)
         datBmat(2,1:n_r_max)=rscheme_oc%rnorm*rscheme_oc%d2rMat(1,:)
         datJmat(1,1:n_r_max)=rscheme_oc%rnorm* rscheme_oc%drMat(1,:)
      else if ( ktopb == 4 ) then
         !----- pseudo vacuum condition, field has only
         !      a radial component, horizontal components
         !      vanish when aj and db are zero:
         datJmat(1,1:n_r_max)=rscheme_oc%rnorm* rscheme_oc%rMat(1,:)
         datBmat(1,1:n_r_max)=rscheme_oc%rnorm*rscheme_oc%drMat(1,:)
      end if

      !-------- at IC (nR=n_r_max):
      if ( l_full_sphere ) then
         datBmat(n_r_max,1:n_r_max)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
         datJmat(n_r_max,1:n_r_max)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
      else
         if ( kbotb == 1 ) then
            !----------- insulating IC, field has to fit a potential field:
            datJmat(n_r_max,1:n_r_max)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)

            datBmat(n_r_max,1:n_r_max)=rscheme_oc%rnorm * (      &
            &                      rscheme_oc%drMat(n_r_max,:) - &
            &    l_P_1*or1(n_r_max)*rscheme_oc%rMat(n_r_max,:) )
         else if ( kbotb == 2 ) then
            !----------- perfect conducting IC
            datBmat(n_r_max,  1:n_r_max)=rscheme_oc%rnorm*  rscheme_oc%rMat(n_r_max,:)
            datBmat(n_r_max-1,1:n_r_max)=rscheme_oc%rnorm*rscheme_oc%d2rMat(n_r_max,:)
            datJmat(n_r_max,  1:n_r_max)=rscheme_oc%rnorm* rscheme_oc%drMat(n_r_max,:)
         else if ( kbotb == 3 ) then
            !---------- finite conducting IC, four boundary conditions:
            !           continuity of b,j, (d b)/(d r) and (d j)/(d r)/sigma.
            !           note: n_r=n_r_max and n_r=n_r_max+1 stand for IC radius
            !           here we set the outer core part of the equations.
            !           the conductivity ratio sigma_ratio is used as
            !           an additional dimensionless parameter.
            datBmat(n_r_max,1:n_r_max)  =rscheme_oc%rnorm* rscheme_oc%rMat(n_r_max,:)
            datBmat(n_r_max+1,1:n_r_max)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
            datJmat(n_r_max,1:n_r_max)  =rscheme_oc%rnorm* rscheme_oc%rMat(n_r_max,:)
            datJmat(n_r_max+1,1:n_r_max)=rscheme_oc%rnorm*sigma_ratio* &
            &                      rscheme_oc%drMat(n_r_max,:)
         else if ( kbotb == 4 ) then
            !----- Pseudovacuum conduction at lower boundary:
            datJmat(n_r_max,1:n_r_max)= rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
            datBmat(n_r_max,1:n_r_max)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         end if
      end if

      !-------- Imposed fields: (overwrites above IC boundary cond.)
      if ( l == 1 .and. ( imagcon == -1 .or. imagcon == -2 ) ) then
         datBmat(n_r_max,1:n_r_max)=rscheme_oc%rnorm * rscheme_oc%rMat(n_r_max,:)
      else if ( l == 3 .and. imagcon == -10 ) then
         if ( l_LCR ) then
            call abortRun('Imposed field not compatible with weak conducting region!')
         end if
         datJmat(1,1:n_r_max)      =rscheme_oc%rnorm * rscheme_oc%rMat(1,:)
         datJmat(n_r_max,1:n_r_max)=rscheme_oc%rnorm * rscheme_oc%rMat(n_r_max,:)
      else if ( n_imp == 1 ) then
         !-- This is Uli Christensen's idea where the external field is
         !   not fixed but compensates the internal field so that the
         !   radial field component vanishes at r/r_cmb=rrMP
         if ( l_LCR ) then
            call abortRun('Imposed field not compatible with weak conducting region!')
         end if
         rRatio=rrMP**real(2*l+1,kind=cp)

         datBmat(1,1:n_r_max)=rscheme_oc%rnorm * (        &
         &                      rscheme_oc%drMat(1,:) +   &
         &     real(l,cp)*or1(1)*rscheme_oc%rMat(1,:) -   &
         &    real(2*l+1,cp)*or1(1)/(1-rRatio) +          &
         &                     conductance_ma* (          &
         &                     rscheme_oc%d2rMat(1,:) -   &
         &            dLh*or2(1)*rscheme_oc%rMat(1,:) ) )
      end if

      !----- fill up with zeros:
      do nR_out=rscheme_oc%n_max+1,n_r_max
         datBmat(1,nR_out)=0.0_cp
         datJmat(1,nR_out)=0.0_cp
         if ( ktopb == 2 ) then
            datBmat(2,nR_out)=0.0_cp
         end if

         if ( l_LCR ) then
            do nR=2,n_r_LCR
               datBmat(nR,nR_out)=0.0_cp
               datJmat(nR,nR_out)=0.0_cp
            end do
         end if

         datBmat(n_r_max,nR_out)  =0.0_cp
         datJmat(n_r_max,nR_out)  =0.0_cp
         if ( kbotb == 2 ) then
            datBmat(n_r_max-1,nR_out)=0.0_cp
         else if ( kbotb == 3 ) then
            datBmat(n_r_max+1,nR_out)=0.0_cp
            datJmat(n_r_max+1,nR_out)=0.0_cp
         end if
      end do

      !----- normalization for highest and lowest Cheb mode:
      do nR=1,n_r_max
         datBmat(nR,1)      =rscheme_oc%boundary_fac*datBmat(nR,1)
         datBmat(nR,n_r_max)=rscheme_oc%boundary_fac*datBmat(nR,n_r_max)
         datJmat(nR,1)      =rscheme_oc%boundary_fac*datJmat(nR,1)
         datJmat(nR,n_r_max)=rscheme_oc%boundary_fac*datJmat(nR,n_r_max)
      end do
      if ( kbotb == 3 ) then
         datBmat(n_r_max+1,1)      =rscheme_oc%boundary_fac*datBmat(n_r_max+1,1)
         datBmat(n_r_max+1,n_r_max)=rscheme_oc%boundary_fac* &
         &                          datBmat(n_r_max+1,n_r_max)
         datJmat(n_r_max+1,1)      =rscheme_oc%boundary_fac*datJmat(n_r_max+1,1)
         datJmat(n_r_max+1,n_r_max)=rscheme_oc%boundary_fac* &
         &                          datJmat(n_r_max+1,n_r_max)
      end if

      !----- Conducting inner core:
      if ( l_cond_ic ) then
         !----- inner core implicit time step matrices for the grid
         !      points n_r=n_r_max+1,...,n_r_max+n_r_ic
         do nCheb=1,n_r_ic_max ! counts even IC cheb modes
            do nR=2,n_r_ic_max-1 ! counts IC radial grid points
               ! n_r=1 represents ICB
               !----------- poloidal field matrix for an inner core field
               !            of the radial form: (r/r_icb)**(l+1)*cheb_ic(r)
               !            where cheb_ic are even chebs only.
               !            NOTE: no hyperdiffusion in inner core !

               datBmat(n_r_max+nR,n_r_max+nCheb) =    &
               &    cheb_norm_ic*dLh*or2(n_r_max) * ( &
               &             cheb_ic(nCheb,nR) -      &
               &     tscheme%wimp_lin(1)*opm*O_sr * ( &
               &                d2cheb_ic(nCheb,nR) + &
               &    two*l_P_1*O_r_ic(nR)*dcheb_ic(nCheb,nR) )   )

               datJmat(n_r_max+nR,n_r_max+nCheb)=datBmat(n_r_max+nR,n_r_max+nCheb)
            end do

            !----- Special treatment for r=0, asymptotic of 1/r dr
            nR=n_r_ic_max
            datBmat(n_r_max+nR,n_r_max+nCheb) =    &
            &    cheb_norm_ic*dLh*or2(n_r_max) * ( &
            &                  cheb_ic(nCheb,nR) - &
            &       tscheme%wimp_lin(1)*opm*O_sr * &
            &    (one+two*l_P_1)*d2cheb_ic(nCheb,nR) )

            datJmat(n_r_max+nR,n_r_max+nCheb)=datBmat(n_r_max+nR,n_r_max+nCheb)
         end do

         !-------- boundary condition at r_icb:
         do nCheb=1,n_cheb_ic_max
            datBmat(n_r_max,n_r_max+nCheb)=-cheb_norm_ic*cheb_ic(nCheb,1)
            datBmat(n_r_max+1,n_r_max+nCheb)=              &
            &    -cheb_norm_ic * ( dcheb_ic(nCheb,1) +     &
            &    l_P_1*or1(n_r_max)*cheb_ic(nCheb,1) )
            datJmat(n_r_max,n_r_max+nCheb)  =datBmat(n_r_max,n_r_max+nCheb)
            datJmat(n_r_max+1,n_r_max+nCheb)=datBmat(n_r_max+1,n_r_max+nCheb)
         end do ! cheb modes

         !-------- fill with zeros:
         do nCheb=n_cheb_ic_max+1,n_r_ic_max
            datBmat(n_r_max,n_r_max+nCheb)  =0.0_cp
            datBmat(n_r_max+1,n_r_max+nCheb)=0.0_cp
            datJmat(n_r_max,n_r_max+nCheb)  =0.0_cp
            datJmat(n_r_max+1,n_r_max+nCheb)=0.0_cp
         end do

         !-------- normalization for lowest Cheb mode:
         do nR=n_r_max,n_r_tot
            datBmat(nR,n_r_max+1)=half*datBmat(nR,n_r_max+1)
            datJmat(nR,n_r_max+1)=half*datJmat(nR,n_r_max+1)
            datBmat(nR,n_r_tot)  =half*datBmat(nR,n_r_tot)
            datJmat(nR,n_r_tot)  =half*datJmat(nR,n_r_tot)
         end do

         !-------- fill matrices up with zeros:
         do nCheb=n_r_max+1,n_r_tot
            do nR=1,n_r_max-1
               datBmat(nR,nCheb)=0.0_cp
               datJmat(nR,nCheb)=0.0_cp
            end do
         end do
         do nCheb=1,n_r_max
            do nR=n_r_max+2,n_r_tot
               datBmat(nR,nCheb)=0.0_cp
               datJmat(nR,nCheb)=0.0_cp
            end do
         end do

      end if  ! conducting inner core ?

#ifdef WITH_PRECOND_BJ
      ! compute the linesum of each line
      do nR=1,nRall
         bMat_fac(nR)=one/maxval(abs(datBmat(nR,1:nRall)))
         datBmat(nR,:) = datBmat(nR,:)*bMat_fac(nR)
      end do
      do nR=1,nRall
         jMat_fac(nR)=one/maxval(abs(datJmat(nR,1:nRall)))
         datJmat(nR,:) = datJmat(nR,:)*jMat_fac(nR)
      end do
#endif

#undef MATRIX_CHECK
#ifdef MATRIX_CHECK
      block

      integer :: i,j
      real(cp) :: rcond
      integer ::ipiv(n_r_tot),iwork(n_r_tot)
      real(cp) :: work(4*n_r_tot),anorm,linesum
      real(cp) :: temp_Mat(n_r_tot,n_r_tot)
      integer, save :: counter=0
      integer :: filehandle
      character(len=100) :: filename

      ! copy the bMat to a temporary variable for modification
      write(filename,"(A,I3.3,A,I3.3,A)") "bMat_",l,"_",counter,".dat"
      open(newunit=filehandle,file=trim(filename))
      counter= counter+1

      do i=1,n_r_tot
         do j=1,n_r_tot
            write(filehandle,"(2ES20.12,1X)",advance="no") datBmat(i,j)
         end do
         write(filehandle,"(A)") ""
      end do
      close(filehandle)
      temp_Mat=datBmat
      anorm = 0.0_cp
      do i=1,n_r_tot
         linesum = 0.0_cp
         do j=1,n_r_tot
            linesum = linesum + abs(temp_Mat(i,j))
         end do
         if (linesum  >  anorm) anorm=linesum
      end do
      !write(*,"(A,ES20.12)") "anorm = ",anorm
      ! LU factorization
      call dgetrf(n_r_tot,n_r_tot,temp_Mat,n_r_tot,ipiv,info)
      ! estimate the condition number
      call dgecon('I',n_r_tot,temp_Mat,n_r_tot,anorm,rcond,work,iwork,info)
      write(*,"(A,I3,A,ES11.3)") "inverse condition number of bMat for l=",l," is ",rcond

      ! The same computation for jMat.
      ! copy the jMat to a temporary variable for modification
      write(filename,"(A,I3.3,A,I3.3,A)") "jMat_",l,"_",counter,".dat"
      open(newunit=filehandle,file=trim(filename))
      counter= counter+1

      do i=1,n_r_tot
         do j=1,n_r_tot
            write(filehandle,"(2ES20.12,1X)",advance="no") datJmat(i,j)
         end do
         write(filehandle,"(A)") ""
      end do
      close(filehandle)
      temp_Mat=datJmat
      anorm = 0.0_cp
      do i=1,n_r_tot
         linesum = 0.0_cp
         do j=1,n_r_tot
            linesum = linesum + abs(temp_Mat(i,j))
         end do
         if (linesum  >  anorm) anorm=linesum
      end do
      !write(*,"(A,ES20.12)") "anorm = ",anorm
      ! LU factorization
      call dgetrf(n_r_tot,n_r_tot,temp_Mat,n_r_tot,ipiv,info)
      ! estimate the condition number
      call dgecon('I',n_r_tot,temp_Mat,n_r_tot,anorm,rcond,work,iwork,info)
      write(*,"(A,I3,A,ES11.3)") "inverse condition number of jMat for l=",l," is ",rcond

      end block
#endif

      !-- Array copy
      call bMat%set_data(datBmat)
      call jMat%set_data(datJmat)

      !----- LU decomposition:
      call bMat%prepare(info)
      if ( info /= 0 ) call abortRun('Singular matrix bMat in get_bmat')

      !----- LU decomposition:
      call jMat%prepare(info)
      if ( info /= 0 ) call abortRun('! Singular matrix jMat in get_bmat!')

   end subroutine get_bMat
!-----------------------------------------------------------------------------
end module updateB_mod
