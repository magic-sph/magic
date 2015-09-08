!$Id$
module magnetic_energy

   use parallel_mod
   use precision_mod
   use truncation, only: n_r_maxMag, n_r_ic_maxMag, n_r_max, n_r_ic_max
   use radial_data, only: n_r_cmb
   use radial_functions, only: r_icb, r_cmb, r_ic, dr_fac_ic, i_costf_init,      &
                               d_costf_init, i_costf1_ic_init, d_costf1_ic_init, &
                               sigma, orho1, r, or2, drx
   use physical_parameters, only: LFfac
   use num_param, only: eScale, tScale
   use blocking, only: st_map, lo_map, lmStartB, lmStopB
   use horizontal_data, only: dLh
   use logic, only: l_cond_ic, l_mag, l_mag_LF, l_save_out
   use movie_data, only: movieDipColat, movieDipLon, movieDipStrength, &
                         movieDipStrengthGeo
   use output_data, only: n_dipole_file, dipole_file, n_e_mag_ic_file,   &
                          e_mag_ic_file, n_e_mag_oc_file, e_mag_oc_file, &
                          tag
   use const, only: pi, zero, one, two, half, four
   use Bext, only: n_imp, rrMP
   use LMLoop_data, only: llmMag, ulmMag
   use integration, only: rInt_R,rIntIC
   use useful, only: cc2real,cc22real
 
   implicit none
 
   private
 
   real(cp), allocatable :: e_dipA(:)
   real(cp), allocatable :: e_pA(:),e_p_asA(:)
   real(cp), allocatable :: e_tA(:),e_t_asA(:)
 
   public :: initialize_magnetic_energy, get_e_mag
  
contains

   subroutine initialize_magnetic_energy

      allocate( e_dipA(n_r_max) )
      allocate( e_pA(n_r_max),e_p_asA(n_r_max) )
      allocate( e_tA(n_r_max),e_t_asA(n_r_max) )
    
   end subroutine initialize_magnetic_energy
!----------------------------------------------------------------------------
   subroutine get_e_mag(time,l_write,l_stop_time,n_e_sets,        &
        &               b,db,aj,b_ic,db_ic,aj_ic,                 &
        &               e_p,e_t,e_p_as,e_t_as,                    &
        &               e_p_ic,e_t_ic,e_p_as_ic,e_t_as_ic,        &
        &               e_p_os,e_p_as_os,e_cmb,Dip,DipCMB,        &
        &               elsAnel)
      !--------------------------------------------------------------------
      !
      !  calculates magnetic energy  = 1/2 Integral(B^2 dV)
      !  integration in theta,phi by summation over harmonic coeffs.
      !  integration in r by Chebycheff integrals
      !
      !  Output:
      !  enbp: Total poloidal        enbt: Total toroidal
      !  apome: Axisym. poloidal     atome: Axisym. toroidal
      !
      !--------------------------------------------------------------------

      !-- Input of variables:
      integer,     intent(in) :: n_e_sets
      real(cp),    intent(in) :: time
      logical,     intent(in) :: l_write
      logical,     intent(in) :: l_stop_time
      complex(cp), intent(in) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: db(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: aj(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: db_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)

      !-- Output variables:
      real(cp), intent(out) :: e_p,e_t       ! poloidal, toroidal energy
      real(cp), intent(out) :: e_p_as,e_t_as ! axisymmetric poloidal, toroidal energy
      real(cp), intent(out) :: e_p_ic,e_t_ic   
      real(cp), intent(out) :: e_p_as_ic,e_t_as_ic
      real(cp), intent(out) :: e_p_os,e_p_as_os
      real(cp), intent(out) :: e_cmb
      real(cp), intent(out) :: Dip,DipCMB
      real(cp), intent(out) :: elsAnel

      !-- local:
      integer :: nR,lm,l,m,l1m0,l1m1
      integer :: l_geo

      real(cp) :: e_p_r(n_r_max), e_p_r_global(n_r_max)
      real(cp) :: e_t_r(n_r_max), e_t_r_global(n_r_max)
      real(cp) :: els_r(n_r_max), els_r_global(n_r_max)
      real(cp) :: e_p_as_r(n_r_max), e_p_as_r_global(n_r_max)
      real(cp) :: e_t_as_r(n_r_max), e_t_as_r_global(n_r_max)
      real(cp) :: e_p_es_r(n_r_max), e_p_es_r_global(n_r_max)
      real(cp) :: e_t_es_r(n_r_max), e_t_es_r_global(n_r_max)
      real(cp) :: e_p_eas_r(n_r_max), e_p_eas_r_global(n_r_max)
      real(cp) :: e_t_eas_r(n_r_max), e_t_eas_r_global(n_r_max)
      real(cp) :: e_dipole_r(n_r_max), e_dipole_r_global(n_r_max)
      real(cp) :: e_dipole_ax_r(n_r_max), e_dipole_ax_r_global(n_r_max)

      real(cp) :: e_p_ic_r(n_r_ic_max), e_p_ic_r_global(n_r_ic_max)
      real(cp) :: e_t_ic_r(n_r_ic_max), e_t_ic_r_global(n_r_ic_max)   
      real(cp) :: e_p_as_ic_r(n_r_ic_max), e_p_as_ic_r_global(n_r_ic_max)
      real(cp) :: e_t_as_ic_r(n_r_ic_max), e_t_as_ic_r_global(n_r_ic_max)

      real(cp) :: e_geo,e_es_geo,e_as_geo,e_eas_geo
      real(cp) :: e_geo_global,e_es_geo_global,e_as_geo_global,e_eas_geo_global
      real(cp) :: e_p_ic_global, e_p_as_ic_global, e_p_os_global, e_p_as_os_global
      real(cp) :: e_p_e,e_p_as_e
      real(cp) :: e_p_e_global,e_p_as_e_global

      real(cp) :: r_ratio,fac
      real(cp) :: e_p_temp,e_t_temp
      real(cp) :: e_dipole, e_dipole_ax, e_dipole_ax_cmb
      real(cp) :: e_dipole_e, e_dipole_e_global
      real(cp) :: e_p_e_ratio
      real(cp) :: O_r_icb_E_2,rad
      real(cp) :: e_p_es,e_t_es,e_es_cmb,e_as_cmb
      real(cp) :: e_p_eas,e_t_eas,e_eas_cmb

      real(cp) :: e_dip_cmb,eTot,eDR
      real(cp) :: theta_dip,phi_dip

      complex(cp) :: r_dr_b,b10,b11

      !-- time averaging of e(r):
      character(len=80) :: filename
      real(cp) :: dt,surf
      real(cp), save :: timeLast,timeTot
      logical :: rank_has_l1m0,rank_has_l1m1
#ifdef WITH_MPI
      integer :: status(MPI_STATUS_SIZE)
#endif
      integer :: sr_tag,request1,request2

      ! some arbitrary send recv tag
      sr_tag=18654

      l_geo=11   ! max degree for geomagnetic field seen on Earth  
      ! surface

      e_p      =0.0_cp
      e_t      =0.0_cp
      e_p_as   =0.0_cp
      e_t_as   =0.0_cp
      e_p_ic   =0.0_cp
      e_t_ic   =0.0_cp
      e_p_as_ic=0.0_cp
      e_t_as_ic=0.0_cp
      e_p_os   =0.0_cp
      e_p_as_os=0.0_cp
      e_geo    =0.0_cp
      e_es_geo =0.0_cp
      e_as_geo =0.0_cp
      e_eas_geo=0.0_cp
      Dip      =0.0_cp
      DipCMB   =0.0_cp

      if ( .not.( l_mag .or. l_mag_LF ) ) return

      do nR=1,n_r_max

         e_p_r(nR)     =0.0_cp
         e_t_r(nR)     =0.0_cp
         e_p_as_r(nR)  =0.0_cp
         e_t_as_r(nR)  =0.0_cp
         e_p_es_r(nR)  =0.0_cp
         e_t_es_r(nR)  =0.0_cp
         e_p_eas_r(nR) =0.0_cp
         e_t_eas_r(nR) =0.0_cp
         e_dipole_r(nR)=0.0_cp
         e_dipole_ax_r(nR)=0.0_cp

         !do lm=2,lm_max
         do lm=max(2,llmMag),ulmMag
            l=lo_map%lm2l(lm)
            m=lo_map%lm2m(lm)

            e_p_temp= dLh(st_map%lm2(l,m))*( dLh(st_map%lm2(l,m))*    &
                 &                       or2(nR)*cc2real( b(lm,nR),m) &
                 &                             + cc2real(db(lm,nR),m) )
            e_t_temp= dLh(st_map%lm2(l,m)) * cc2real(aj(lm,nR),m)

            if ( m == 0 ) then  ! axisymmetric part 
               e_p_as_r(nR)=e_p_as_r(nR) + e_p_temp
               e_t_as_r(nR)=e_t_as_r(nR) + e_t_temp
               if ( mod(l,2) == 1 ) then
                  e_p_eas_r(nR)=e_p_eas_r(nR)+e_p_temp
               else
                  e_t_eas_r(nR)=e_t_eas_r(nR)+e_t_temp
               end if
            else
               e_p_r(nR)=e_p_r(nR) + e_p_temp
               e_t_r(nR)=e_t_r(nR) + e_t_temp
            end if
            if ( mod(l+m,2) == 1 ) then
               e_p_es_r(nR)=e_p_es_r(nR) + e_p_temp
            else
               e_t_es_r(nR)=e_t_es_r(nR) + e_t_temp
            end if
            if ( l <= l_geo .and. nR == n_r_cmb ) then     
               e_geo    =e_geo    +e_p_temp
               if ( mod(l+m,2) == 1 ) e_es_geo = e_es_geo + e_p_temp
               if ( m == 0 )          e_as_geo = e_as_geo + e_p_temp
               if ( mod(l+m,2) == 1 .and. m == 0 )                    &
                    &                 e_eas_geo=e_eas_geo+e_p_temp
            end if

            if ( l == 1 .and. m == 0 ) e_dipole_ax_r(nR)=e_p_temp
            if ( l == 1 )  e_dipole_r(nR)=e_dipole_r(nR)+e_p_temp

         end do    ! do loop over lms in block 

         e_p_r(nR)=e_p_r(nR)+e_p_as_r(nR)
         e_t_r(nR)=e_t_r(nR)+e_t_as_r(nR)

         ! In anelastic models it is also interesting to have Lambda
         els_r(nR)=(e_p_r(nR)+e_t_r(nR))*orho1(nR)*sigma(nR)

      end do    ! radial grid points 

#ifdef WITH_MPI
      ! reduce over the ranks
      call MPI_Reduce(e_p_r,    e_p_r_global,     n_r_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_t_r,    e_t_r_global,     n_r_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_p_as_r, e_p_as_r_global,  n_r_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_t_as_r, e_t_as_r_global,  n_r_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_p_es_r, e_p_es_r_global,  n_r_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_t_es_r, e_t_es_r_global,  n_r_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_p_eas_r,e_p_eas_r_global, n_r_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_t_eas_r,e_t_eas_r_global, n_r_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_dipole_ax_r,e_dipole_ax_r_global, n_r_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_dipole_r,   e_dipole_r_global,    n_r_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(els_r,   els_r_global,    n_r_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      ! reduce some scalars
      call MPI_Reduce(e_geo,     e_geo_global,    1, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_es_geo,  e_es_geo_global, 1, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_as_geo,  e_as_geo_global, 1, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_eas_geo, e_eas_geo_global,1, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
      e_p_r_global        =e_p_r
      e_t_r_global        =e_t_r
      e_p_as_r_global     =e_p_as_r
      e_t_as_r_global     =e_t_as_r
      e_p_es_r_global     =e_p_es_r
      e_t_es_r_global     =e_t_es_r
      e_p_eas_r_global    =e_p_eas_r
      e_t_eas_r_global    =e_t_eas_r
      e_dipole_ax_r_global=e_dipole_ax_r
      e_dipole_r_global   =e_dipole_r
      els_r_global        =els_r
      e_geo_global        =e_geo
      e_es_geo_global     =e_es_geo
      e_as_geo_global     =e_as_geo
      e_eas_geo_global    =e_eas_geo
#endif

      if ( rank == 0 ) then
         !-- Get Values at CMB:
         e_cmb          =e_p_r_global(n_r_cmb)+e_t_r_global(n_r_cmb)
         e_dip_cmb      =e_dipole_r_global(n_r_cmb)
         e_dipole_ax_cmb=e_dipole_ax_r_global(n_r_cmb)
         e_es_cmb       =e_p_es_r_global(n_r_cmb)
         e_as_cmb       =e_p_as_r_global(n_r_cmb)
         e_eas_cmb      =e_p_eas_r_global(n_r_cmb)

         ! NOTE: n_e_sets=0 prevents averaging
         if ( n_e_sets == 1 ) then
            timeTot=one
            do nR=1,n_r_max
               e_dipA(nR) =e_dipole_r_global(nR)
               e_pA(nR)   =e_p_r_global(nR)
               e_p_asA(nR)=e_p_r_global(nR)
               e_tA(nR)   =e_t_r_global(nR)
               e_t_asA(nR)=e_t_r_global(nR)
            end do
         else if ( n_e_sets == 2 ) then
            dt=time-timeLast
            timeTot=two*dt
            do nR=1,n_r_max
               e_dipA(nR) =dt*(e_dipA(nR) +e_dipole_r_global(nR))
               e_pA(nR)   =dt*(e_pA(nR)   +e_p_r_global(nR)     )
               e_p_asA(nR)=dt*(e_p_asA(nR)+e_p_as_r_global(nR)  )
               e_tA(nR)   =dt*(e_tA(nR)   +e_t_r_global(nR)     )
               e_t_asA(nR)=dt*(e_t_asA(nR)+e_t_as_r_global(nR)  )
            end do
         else
            dt=time-timeLast
            timeTot=timeTot+dt
            do nR=1,n_r_max
               e_dipA(nR) =e_dipA(nR) +dt*e_dipole_r_global(nR)
               e_pA(nR)   =e_pA(nR)   +dt*e_p_r_global(nR)
               e_p_asA(nR)=e_p_asA(nR)+dt*e_p_as_r_global(nR)
               e_tA(nR)   =e_tA(nR)   +dt*e_t_r_global(nR)
               e_t_asA(nR)=e_t_asA(nR)+dt*e_t_as_r_global(nR)
            end do
         end if
         if ( l_stop_time ) then
            fac=half*LFfac*eScale
            filename='eMagR.'//tag
            open(99, file=filename, status='unknown')
            do nR=1,n_r_max
               eTot=e_pA(nR)+e_tA(nR)
               if ( e_dipA(nR)  <  1.e-6_cp*eTot ) then
                  eDR=0.0_cp
               else
                  eDR=e_dipA(nR)/eTot
               end if
               surf=four*pi*r(nR)**2
               write(99,'(ES20.10,9ES15.7)') r(nR),                 &
                    &               fac*e_pA(nR)/timetot,           &
                    &               fac*e_p_asA(nR)/timetot,        &
                    &               fac*e_tA(nR)/timetot,           &
                    &               fac*e_t_asA(nR)/timetot,        &
                    &               fac*e_pA(nR)/timetot/surf,      &
                    &               fac*e_p_asA(nR)/timetot/surf,   &
                    &               fac*e_tA(nR)/timetot/surf,      &
                    &               fac*e_t_asA(nR)/timetot/surf,   &
                    &               eDR
            end do
            close(99)
         end if
         timeLast=time


         !-- Radial integrals:
         e_p        =rInt_R(e_p_r_global,n_r_max,n_r_max,drx, &
                            i_costf_init,d_costf_init)
         e_t        =rInt_R(e_t_r_global,n_r_max,n_r_max,drx, &
                            i_costf_init,d_costf_init)
         e_p_as     =rInt_R(e_p_as_r_global,n_r_max,n_r_max,drx, &
                            i_costf_init,d_costf_init)
         e_t_as     =rInt_R(e_t_as_r_global,n_r_max,n_r_max,drx, &
                            i_costf_init,d_costf_init)
         e_p_es     =rInt_R(e_p_es_r_global,n_r_max,n_r_max,drx, &
                            i_costf_init,d_costf_init)
         e_t_es     =rInt_R(e_t_es_r_global,n_r_max,n_r_max,drx, &
                            i_costf_init,d_costf_init)
         e_p_eas    =rInt_R(e_p_eas_r_global,n_r_max,n_r_max,drx, &
                            i_costf_init,d_costf_init)
         e_t_eas    =rInt_R(e_t_eas_r_global,n_r_max,n_r_max,drx, &
                            i_costf_init,d_costf_init)
         e_dipole   =rInt_R(e_dipole_r_global,n_r_max,n_r_max,drx, &
                            i_costf_init,d_costf_init)
         e_dipole_ax=rInt_R(e_dipole_ax_r_global,n_r_max,n_r_max,drx, &
                            i_costf_init,d_costf_init)
         elsAnel    =rInt_R(els_r_global,n_r_max,n_r_max,drx, &
                            i_costf_init,d_costf_init)

         fac=0.5*LFfac*eScale
         e_p        =fac*e_p
         e_t        =fac*e_t
         e_p_as     =fac*e_p_as
         e_t_as     =fac*e_t_as
         e_p_es     =fac*e_p_es
         e_t_es     =fac*e_t_es
         e_p_eas    =fac*e_p_eas
         e_t_eas    =fac*e_t_eas
         e_dipole   =fac*e_dipole
         e_dipole_ax=fac*e_dipole_ax
         elsAnel    =eScale*elsAnel

         e_cmb          =fac*e_cmb
         e_dip_cmb      =fac*e_dip_cmb
         e_dipole_ax_cmb=fac*e_dipole_ax_cmb
         e_es_cmb       =fac*e_es_cmb
         e_as_cmb       =fac*e_as_cmb
         e_eas_cmb      =fac*e_eas_cmb
         e_geo          =fac*e_geo_global
         e_es_geo       =fac*e_es_geo_global
         e_as_geo       =fac*e_as_geo_global
         e_eas_geo      =fac*e_eas_geo_global
      end if

      !-- Inner core:

      if ( l_cond_ic ) then 

         O_r_icb_E_2=one/(r(n_r_max)*r(n_r_max))

         do nR=1,n_r_ic_max

            r_ratio=r_ic(nR)/r_ic(1)

            e_p_ic_r(nR)   =0.0_cp
            e_t_ic_r(nR)   =0.0_cp
            e_p_as_ic_r(nR)=0.0_cp
            e_t_as_ic_r(nR)=0.0_cp

            !do lm=2,lm_max
            do lm=max(2,llmMag),ulmMag
               l=lo_map%lm2l(lm)
               m=lo_map%lm2m(lm)
               r_dr_b=r_ic(nR)*db_ic(lm,nR)

               e_p_temp=     dLh(st_map%lm2(l,m))*O_r_icb_E_2*r_ratio**(2*l) * (     &
                    &           real((l+1)*(2*l+1),cp)*cc2real(b_ic(lm,nR),m)     +  &
                    &           real(2*(l+1),cp)*cc22real(b_ic(lm,nR),r_dr_b,m)   +  &
                    &                                 cc2real(r_dr_b,m)            )
               e_t_temp=  dLh(st_map%lm2(l,m))*r_ratio**(2*l+2) *                  &
                    &                             cc2real(aj_ic(lm,nR),m)

               if ( m == 0 ) then  ! axisymmetric part
                  e_p_as_ic_r(nR)=e_p_as_ic_r(nR) + e_p_temp
                  e_t_as_ic_r(nR)=e_t_as_ic_r(nR) + e_t_temp
               else
                  e_p_ic_r(nR)   =e_p_ic_r(nR) + e_p_temp
                  e_t_ic_r(nR)   =e_t_ic_r(nR) + e_t_temp
               end if

            end do    ! do loop over lms in block

            e_p_ic_r(nR)=e_p_ic_r(nR)+e_p_as_ic_r(nR)
            e_t_ic_r(nR)=e_t_ic_r(nR)+e_t_as_ic_r(nR)

         end do    ! radial grid points

         ! reduce over the ranks
#ifdef WITH_MPI
         call MPI_Reduce(e_p_ic_r,    e_p_ic_r_global,    n_r_ic_max, &
              & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_Reduce(e_t_ic_r,    e_t_ic_r_global,    n_r_ic_max, &
              & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_Reduce(e_p_as_ic_r, e_p_as_ic_r_global, n_r_ic_max, &
              & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_Reduce(e_t_as_ic_r, e_t_as_ic_r_global, n_r_ic_max, &
              & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
         e_p_ic_r_global   =e_p_ic_r
         e_t_ic_r_global   =e_t_ic_r
         e_p_as_ic_r_global=e_p_as_ic_r
         e_t_as_ic_r_global=e_t_as_ic_r
#endif

         if ( rank == 0 ) then
            e_p_ic   =rIntIC(e_p_ic_r_global,n_r_ic_max,dr_fac_ic,              &
                 &                      i_costf1_ic_init,d_costf1_ic_init)
            e_t_ic   =rIntIC(e_t_ic_r_global,n_r_ic_max,dr_fac_ic,              &
                 &                      i_costf1_ic_init,d_costf1_ic_init)
            e_p_as_ic=rIntIC(e_p_as_ic_r_global,n_r_ic_max,dr_fac_ic,           &
                 &                      i_costf1_ic_init,d_costf1_ic_init)
            e_t_as_ic=rIntIC(e_t_as_ic_r_global,n_r_ic_max,dr_fac_ic,           &
                 &                      i_costf1_ic_init,d_costf1_ic_init)
            fac=half*LFfac*eScale
            e_p_ic   =fac*e_p_ic
            e_t_ic   =fac*e_t_ic
            e_p_as_ic=fac*e_p_as_ic
            e_t_as_ic=fac*e_t_as_ic
         end if

      else

         !do lm=2,lm_max
         do lm=max(2,llmMag),ulmMag
            l=lo_map%lm2l(lm)
            m=lo_map%lm2m(lm)
            fac=real(l*(l+1)*(l+1),cp)
            e_p_temp=fac*cc2real(b(lm,n_r_max),m)
            e_p_ic=e_p_ic + e_p_temp
            if ( m == 0 ) e_p_as_ic=e_p_as_ic+e_p_temp
         end do    ! do loop over lms in block

#ifdef WITH_MPI
         call MPI_Reduce(e_p_ic,    e_p_ic_global,   1, &
              & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_Reduce(e_p_as_ic, e_p_as_ic_global,1, &
              & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
         e_p_ic_global   =e_p_ic
         e_p_as_ic_global=e_p_as_ic
#endif

         if (rank == 0) then
            fac      =half*LFfac/r_icb*eScale
            e_p_ic   =fac*e_p_ic_global
            e_t_ic   =0.0_cp
            e_p_as_ic=fac*e_p_as_ic_global
            e_t_as_ic=0.0_cp
         end if

      end if  ! conducting inner core ?


      !-- Outside energy:
      nR=n_r_cmb
      e_p_os   =0.0_cp
      e_p_as_os=0.0_cp
      !do lm=2,lm_max
      do lm=max(2,llmMag),ulmMag
         l=lo_map%lm2l(lm)
         m=lo_map%lm2m(lm)
         fac=real( l*l*(l+1),cp)
         e_p_temp=fac*cc2real(b(lm,nR),m)
         e_p_os  =e_p_os + e_p_temp
         if ( m == 0 ) e_p_as_os=e_p_as_os + e_p_temp
      end do

#ifdef WITH_MPI
      call MPI_Reduce(e_p_os,    e_p_os_global,   1, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_p_as_os, e_p_as_os_global,1, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
      e_p_os_global   =e_p_os
      e_p_as_os_global=e_p_as_os
#endif

      if ( rank == 0 ) then
         fac      =half*LFfac/r_cmb*eScale
         e_p_os   =fac*e_p_os_global
         e_p_as_os=fac*e_p_as_os_global
      end if

      !-- External potential field energy in Uli case (n_imp=1)
      e_p_e     =0.0_cp
      e_p_as_e  =0.0_cp
      e_dipole_e=0.0_cp
      if ( n_imp == 1 ) then
         !do lm=2,lm_max
         do lm=max(2,llmMag),ulmMag
            l=lo_map%lm2l(lm)
            m=lo_map%lm2m(lm)
            fac=real(l*(l+1)**2*(2*l+1),cp)*one/(rrMP**(2*l+1)-one)
            e_p_temp=fac*cc2real(b(lm,nR),m)
            e_p_e   =e_p_e  + e_p_temp
            if ( m == 0 ) e_p_as_e =e_p_as_e  + e_p_temp
            if ( l == 1 ) e_dipole_e=e_dipole_e+e_p_temp
         end do

#ifdef WITH_MPI
         call MPI_Reduce(e_p_e,    e_p_e_global,   1, &
              & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_Reduce(e_p_as_e, e_p_as_e_global,1, &
              & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_Reduce(e_dipole_e, e_dipole_e_global,1, &
              & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
         e_p_e_global     =e_p_e
         e_p_as_e_global  =e_p_as_e
         e_dipole_e_global=e_dipole_e
#endif
         
         if ( rank == 0 ) then
            fac       =half*LFfac/r_cmb**2*eScale
            e_p_e     =fac*e_p_e_global 
            e_p_as_e  =fac*e_p_as_e_global
            e_dipole_e=fac*e_dipole_e_global
         end if
      end if


      !-- Output of OC and outside energies:
      if ( rank == 0 ) then
         if ( l_write ) then
            if ( l_save_out ) then
               open(n_e_mag_oc_file, file=e_mag_oc_file, status='unknown', &
                    &             position='append')
            end if
            write(n_e_mag_oc_file,'(1P,ES20.12,12ES16.8)')               &
                 &                             time*tScale,              &! 1
                 &                             e_p,e_t,                  &! 2,3
                 &                             e_p_as,e_t_as,            &! 4,5
                 &                             e_p_os,e_p_as_os,         &! 6,7
                 &                             e_p_es,e_t_es,            &! 8,9
                 &                             e_p_eas,e_t_eas,          &! 10,11
                 &                             e_p_e,e_p_as_e    ! 12,13
            if ( l_save_out ) close(n_e_mag_oc_file)
         end if

         !-- Output of IC energies:
         if ( l_write ) then
            if ( l_save_out ) then
               open(n_e_mag_ic_file, file=e_mag_ic_file, status='unknown', &
                    &             position='append')
            end if
            write(n_e_mag_ic_file,'(1P,ES20.12,4ES16.8)')             &
                 &                       time*tScale,                 &
                 &                       e_p_ic,e_t_ic,               &
                 &                       e_p_as_ic,e_t_as_ic
            if ( l_save_out ) close(n_e_mag_ic_file)
         end if
      end if

      l1m0=lo_map%lm2(1,0)
      l1m1=lo_map%lm2(1,1)
      !if ((l1m0 < lmStartB(1)).or.(l1m0 > lmStopB(1))) then
      !   if (rank == 0) PRINT*,"in get_e_mag, dipole part: l1m0 not on rank 0"
      !   stop
      !end if
      !if (  ( l1m1 > 0).and.&
      !     &( (l1m1 < lmStartB(1)).or.(l1m1 > lmStopB(1)) ) ) then
      !   if (rank == 0) PRINT*,"in get_e_mag, dipole part: l1m1 not on rank 0:"
      !   stop
      !end if
      rank_has_l1m0=.false.
      rank_has_l1m1=.false.
      !write(*,"(I5,A,2I5,A,2I5)") rank,": l1m0,l1m1 = ",l1m0,l1m1,&
      !     & ", lm block: ",lmStartB(rank+1),lmStopB(rank+1)
      if ( (l1m0 >= lmStartB(rank+1)) .and. (l1m0 <= lmStopB(rank+1)) ) then
         b10=b(l1m0,n_r_cmb)
#ifdef WITH_MPI
         if (rank /= 0) then
            call MPI_Send(b10,1,MPI_DEF_COMPLEX,0,sr_tag,MPI_COMM_WORLD,ierr)
         end if
#endif
         rank_has_l1m0=.true.
      end if
      if ( l1m1 > 0 ) then
         if ( (l1m1 >= lmStartB(rank+1)) .and. (l1m1 <= lmStopB(rank+1)) ) then
            b11=b(l1m1,n_r_cmb)
#ifdef WITH_MPI
            if (rank /= 0) then
               call MPI_Send(b11,1,MPI_DEF_COMPLEX,0,sr_tag+1,MPI_COMM_WORLD,ierr)
            end if
#endif
            rank_has_l1m1=.true.
         end if
      else
         b11=zero
         rank_has_l1m1=.true.
      end if

         
      if ( rank == 0 ) then
         !-- Calculate pole position:
         rad =180.0_cp/pi
#ifdef WITH_MPI
         if (.not.rank_has_l1m0) then
            call MPI_IRecv(b10,1,MPI_DEF_COMPLEX,MPI_ANY_SOURCE,&
                 &         sr_tag,MPI_COMM_WORLD,request1, ierr)
         end if
         if ( .not. rank_has_l1m1 ) then
            call MPI_IRecv(b11,1,MPI_DEF_COMPLEX,MPI_ANY_SOURCE,&
                 &         sr_tag+1,MPI_COMM_WORLD,request2,ierr)
         end if
         if ( .not. rank_has_l1m0 ) then
            call MPI_Wait(request1,status,ierr)
         end if
         if ( .not. rank_has_l1m1 ) then
            call MPI_Wait(request2,status,ierr)
         end if
#endif

         !print*, "------------", b10, b11
         theta_dip= rad*atan2(sqrt(two)*abs(b11),real(b10))
         if ( theta_dip < 0.0_cp ) theta_dip=180.0_cp+theta_dip
         if ( abs(b11) < 1.e-20_cp ) then
            phi_dip=0.0_cp
         else
            phi_dip=-rad*atan2(aimag(b11),real(b11))
         end if
         Dip      =e_dipole_ax/(e_p+e_t)
         DipCMB   =e_dipole_ax_cmb/e_cmb

         !-- Output of pole position:
         if ( l_write ) then
            if ( l_save_out ) then
               open(n_dipole_file, file=dipole_file, status='unknown',   &
                    &             position='append')
            end if
            if ( e_p_e == 0 ) then
               e_p_e_ratio=0.0_cp
            else
               e_p_e_ratio=e_dipole_e/e_p_e
            end if
            !write(*,"(A,4ES25.17)") "e_cmb = ",e_cmb,e_es_cmb,e_cmb-e_es_cmb,(e_cmb-e_es_cmb)/e_cmb
            ! There are still differences in field 17 of the dipole file. These
            ! differences are due to the summation for e_es_cmb and are only of the order
            ! of machine accuracy.
            write(n_dipole_file,'(1P,ES20.12,19ES16.8)')   &
                 &       time*tScale,                      &! 1
                 &       theta_dip,phi_dip,                &! 2,3
                 &       Dip,                              &! 4  
                 &       DipCMB,                           &! 5
                 &       e_dipole_ax_cmb/e_geo,            &! 6
                 &       e_dipole/(e_p+e_t),               &! 7
                 &       e_dip_cmb/e_cmb,                  &! 8
                 &       e_dip_cmb/e_geo,                  &! 9
                 &       e_dip_cmb,e_dipole_ax_cmb,        &! 10,11
                 &       e_dipole,e_dipole_ax,             &! 12,13
                 &       e_cmb,e_geo,                      &! 14,15
                 &       e_p_e_ratio,                      &! 16
                 &       (e_cmb-e_es_cmb)/e_cmb,           &! 17
                 &       (e_cmb-e_as_cmb)/e_cmb,           &! 18
                 &       (e_geo-e_es_geo)/e_geo,           &! 19
                 &       (e_geo-e_as_geo)/e_geo        ! 20
            if ( l_save_out ) close(n_dipole_file)
         end if
         ! Store values needed for movie output:
         movieDipColat      =theta_dip
         movieDipLon        =phi_dip
         movieDipStrength   =e_dip_cmb/e_cmb
         movieDipStrengthGeo=e_dip_cmb/e_geo
      end if


   end subroutine get_e_mag
!----------------------------------------------------------------------------
end module magnetic_energy
