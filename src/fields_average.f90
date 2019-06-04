module fields_average_mod
   !
   ! This module is used when one wants to store time-averaged quantities
   !

   use truncation
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use radial_data, only: n_r_cmb, n_r_icb
   use radial_functions, only: chebt_ic, chebt_ic_even, r, dr_fac_ic, &
       &                       rscheme_oc
   use blocking,only: lmStartB, lmStopB, sizeThetaB, nThetaBs, lm2, nfs
   use logic, only: l_mag, l_conv, l_save_out, l_heat, l_cond_ic, &
       &            l_chemical_conv
   use kinetic_energy, only: get_e_kin
   use magnetic_energy, only: get_e_mag
   use output_data, only: tag, n_log_file, log_file, n_graphs, l_max_cmb
   use parallel_mod, only: rank
#ifdef WITH_SHTNS
   use shtns
#else
   use horizontal_data, only: Plm, dPlm, dLh
   use fft, only: fft_thetab
   use legendre_spec_to_grid, only: legTF
#endif
   use constants, only: zero, vol_oc, vol_ic, one
   use LMLoop_data, only: llm,ulm,llmMag,ulmMag
   use communications, only: get_global_sum, gather_from_lo_to_rank0,&
       &                     gather_all_from_lo_to_rank0,gt_OC,gt_IC
   use out_coeff, only: write_Bcmb, write_Pot
   use spectra, only: spectrum, spectrum_temp
   use graphOut_mod, only: graphOut, graphOut_IC, n_graph_file
   use leg_helper_mod, only: legPrep
   use radial_der_even, only: get_drNS_even, get_ddrNS_even
   use radial_der, only: get_dr
   use fieldsLast, only: dwdtLast_LMloc, dpdtLast_LMloc, dzdtLast_lo,     &
       &                 dsdtLast_LMloc, dxidtLast_LMloc, dbdtLast_LMloc, &
       &                 djdtLast_LMloc, dbdt_icLast_LMloc,               &
       &                 djdt_icLast_LMloc
   use storeCheckPoints, only: store

   implicit none

   private

   complex(cp), allocatable :: w_ave(:,:)
   complex(cp), allocatable :: z_ave(:,:)
   complex(cp), allocatable :: s_ave(:,:)
   complex(cp), allocatable :: xi_ave(:,:)
   complex(cp), allocatable :: p_ave(:,:)
   complex(cp), allocatable :: b_ave(:,:)
   complex(cp), allocatable :: aj_ave(:,:)
   complex(cp), allocatable :: b_ic_ave(:,:)
   complex(cp), allocatable :: aj_ic_ave(:,:)
   ! on rank 0 we also allocate the following fields
   complex(cp), allocatable :: b_ave_global(:), bICB(:)
   complex(cp), allocatable :: db_ave_global(:), aj_ave_global(:)
   complex(cp), allocatable :: w_ave_global(:), dw_ave_global(:)
   complex(cp), allocatable :: z_ave_global(:), s_ave_global(:)
   complex(cp), allocatable :: p_ave_global(:), xi_ave_global(:)

   public :: initialize_fields_average_mod, fields_average, &
   &         finalize_fields_average_mod


contains

   subroutine initialize_fields_average_mod

      allocate( w_ave(llm:ulm,n_r_max) )
      allocate( z_ave(llm:ulm,n_r_max) )
      allocate( s_ave(llm:ulm,n_r_max) )
      allocate( p_ave(llm:ulm,n_r_max) )
      allocate( b_ave(llm:ulm,n_r_max) )
      allocate( aj_ave(llm:ulm,n_r_max) )
      bytes_allocated = bytes_allocated+6*(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
      allocate( b_ic_ave(llm:ulm,n_r_ic_max) )
      allocate( aj_ic_ave(llm:ulm,n_r_ic_max) )
      bytes_allocated = bytes_allocated+2*(ulm-llm+1)*n_r_ic_max*SIZEOF_DEF_COMPLEX

      if ( l_chemical_conv ) then
         allocate( xi_ave(llm:ulm,n_r_max) )
         bytes_allocated = bytes_allocated+(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
      else
         allocate( xi_ave(1,1) )
      end if

      if ( rank == 0 ) then
         allocate( bICB(1:lm_max) )
         allocate( b_ave_global(1:lm_max) )
         allocate( db_ave_global(1:lm_max) )
         allocate( aj_ave_global(1:lm_max) )
         allocate( w_ave_global(1:lm_max) )
         allocate( dw_ave_global(1:lm_max) )
         allocate( z_ave_global(1:lm_max) )
         allocate( s_ave_global(1:lm_max) )
         allocate( p_ave_global(1:lm_max) )
         bytes_allocated = bytes_allocated+9*lm_max*SIZEOF_DEF_COMPLEX
         if ( l_chemical_conv ) then
            allocate( xi_ave_global(1:lm_max) )
            bytes_allocated = bytes_allocated+lm_max*SIZEOF_DEF_COMPLEX
         end if
      else
         allocate( bICB(1) )
         allocate( b_ave_global(1) )
         allocate( db_ave_global(1) )
         allocate( aj_ave_global(1) )
         allocate( w_ave_global(1) )
         allocate( dw_ave_global(1) )
         allocate( z_ave_global(1) )
         allocate( s_ave_global(1) )
         allocate( p_ave_global(1) )
         if ( l_chemical_conv ) then
            allocate( xi_ave_global(1) )
         end if
      end if

   end subroutine initialize_fields_average_mod
!----------------------------------------------------------------------------
   subroutine finalize_fields_average_mod

      deallocate( w_ave, z_ave, s_ave, p_ave, b_ave, aj_ave, b_ic_ave )
      deallocate( aj_ic_ave, db_ave_global, aj_ave_global, w_ave_global )
      deallocate( dw_ave_global, z_ave_global, s_ave_global, p_ave_global )
      deallocate( b_ave_global, bICB )

      if ( l_chemical_conv ) deallocate( xi_ave, xi_ave_global )

   end subroutine finalize_fields_average_mod
!----------------------------------------------------------------------------
   subroutine fields_average(simtime,dt,dtNew,nAve,l_stop_time,       &
      &                      time_passed,time_norm,omega_ic,omega_ma, &
      &                      w,z,p,s,xi,b,aj,b_ic,aj_ic)
      !
      ! This subroutine averages fields b and v over time.
      !

      !-- Input of variables:
      integer,     intent(in) :: nAve         ! number for averaged time steps
      logical,     intent(in) :: l_stop_time  ! true if this is the last time step
      real(cp),    intent(in) :: time_passed  ! time passed since last log
      real(cp),    intent(in) :: time_norm    ! time passed since start of time loop
      real(cp),    intent(in) :: omega_ic,omega_ma
      real(cp),    intent(in) :: dt,dtNew
      real(cp),    intent(in) :: simtime
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(in) :: z(llm:ulm,n_r_max)
      complex(cp), intent(in) :: p(llm:ulm,n_r_max)
      complex(cp), intent(in) :: s(llm:ulm,n_r_max)
      complex(cp), intent(in) :: xi(llm:ulm,n_r_max)
      complex(cp), intent(in) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: aj(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)

      !-- Local stuff:
      ! fields for the gathering
      complex(cp) :: b_ic_ave_global(1:lm_maxMag,n_r_ic_maxMag)
      complex(cp) :: db_ic_ave_global(1:lm_maxMag,n_r_ic_maxMag)
      complex(cp) :: ddb_ic_ave_global(1:lm_maxMag,n_r_ic_maxMag)
      complex(cp) :: aj_ic_ave_global(1:lm_maxMag,n_r_ic_maxMag)
      complex(cp) :: dj_ic_ave_global(1:lm_maxMag,n_r_ic_maxMag)

      !----- Time averaged fields:
      complex(cp) :: dw_ave(llm:ulm,n_r_max)
      complex(cp) :: ds_ave(llm:ulm,n_r_max)
      complex(cp) :: dxi_ave(llm:ulm,n_r_max)
      complex(cp) :: db_ave(llm:ulm,n_r_max)
      complex(cp) :: db_ic_ave(llm:ulm,n_r_ic_max)
      complex(cp) :: ddb_ic_ave(llm:ulm,n_r_ic_max)
      complex(cp) :: dj_ic_ave(llm:ulm,n_r_ic_max)

      !----- Work array:
      complex(cp) :: workA_LMloc(llm:ulm,n_r_max)

      !----- Fields in grid space:
      real(cp) :: Br(nrp,nfs),Bt(nrp,nfs),Bp(nrp,nfs) ! B field comp.
      real(cp) :: Vr(nrp,nfs),Vt(nrp,nfs),Vp(nrp,nfs) ! B field comp.
      real(cp) :: Sr(nrp,nfs),Prer(nrp,nfs)           ! entropy,pressure
      real(cp) :: Xir(nrp,nfs)                        ! chemical composition

      !----- Help arrays for fields:
#ifndef WITH_SHTNS
      complex(cp) :: dLhb(lm_max),bhG(lm_max),bhC(lm_max)
      complex(cp) :: dLhw(lm_max),vhG(lm_max),vhC(lm_max)
      integer :: nThetaB, nThetaStart
#endif

      !----- Energies of time average field:
      real(cp) :: e_kin_p_ave,e_kin_t_ave
      real(cp) :: e_kin_p_as_ave,e_kin_t_as_ave
      real(cp) :: e_mag_p_ave,e_mag_t_ave
      real(cp) :: e_mag_p_as_ave,e_mag_t_as_ave
      real(cp) :: e_mag_p_ic_ave,e_mag_t_ic_ave
      real(cp) :: e_mag_p_as_ic_ave,e_mag_t_as_ic_ave
      real(cp) :: e_mag_os_ave,e_mag_as_os_ave
      real(cp) :: Dip,DipCMB,e_cmb,elsAnel

      integer :: lm, nR
      integer :: n_e_sets,n_spec

      character(len=72) :: graph_file
      character(len=80) :: outFile
      integer :: nOut,n_cmb_sets

      logical :: lGraphHeader

      real(cp) :: time
      real(cp) :: dt_norm

      integer :: nBpotSets,nVpotSets,nTpotSets
      integer :: lmStart,lmStop

      !-- Initialise average for first time step:

      if ( nAve == 1 ) then

         !zero=zero
         if ( n_graphs > 0 ) then
            if ( l_conv ) then
               w_ave=zero
               z_ave=zero
               p_ave=zero
            end if
            if ( l_heat ) then
               s_ave=zero
            end if
            if ( l_chemical_conv ) then
               xi_ave=zero
            end if
            if ( l_mag ) then
               b_ave=zero
               aj_ave=zero
               if ( l_cond_ic ) then
                  b_ic_ave=zero
                  aj_ic_ave=zero
               end if
            end if
         end if

      end if  ! First step

      !-- Add new time step:

      if ( l_conv ) then
         do nR=1,n_r_max
            do lm=llm,ulm
               w_ave(lm,nR)=w_ave(lm,nR) + time_passed*w(lm,nR)
               z_ave(lm,nR)=z_ave(lm,nR) + time_passed*z(lm,nR)
               p_ave(lm,nR)=p_ave(lm,nR) + time_passed*p(lm,nR)
            end do
         end do
      end if
      if ( l_heat ) then
         do nR=1,n_r_max
            do lm=llm,ulm
               s_ave(lm,nR)=s_ave(lm,nR) + time_passed*s(lm,nR)
            end do
         end do
      end if
      if ( l_chemical_conv ) then
         do nR=1,n_r_max
            do lm=llm,ulm
               xi_ave(lm,nR)=xi_ave(lm,nR) + time_passed*xi(lm,nR)
            end do
         end do
      end if
      if ( l_mag ) then
         do nR=1,n_r_max
            do lm=llm,ulm
               b_ave(lm,nR) =b_ave(lm,nR)  + time_passed*b(lm,nR)
               aj_ave(lm,nR)=aj_ave(lm,nR) + time_passed*aj(lm,nR)
            end do
         end do
         if ( l_cond_ic ) then
            do nR=1,n_r_ic_max
               do lm=llm,ulm
                  b_ic_ave(lm,nR) =b_ic_ave(lm,nR) + time_passed*b_ic(lm,nR)
                  aj_ic_ave(lm,nR)=aj_ic_ave(lm,nR)+ time_passed*aj_ic(lm,nR)
               end do
            end do
         end if
      end if

      !--- Output, intermediate output every 10th averaging to save result
      !    will be overwritten.
      if ( l_stop_time .or. mod(nAve,10) == 0 ) then

         !write(*,"(A,2ES22.15)") "w_ave = ",get_global_sum( w_ave )
         time   =-one  ! This signifies averaging in output files!
         dt_norm=one/time_norm

         if ( l_conv ) then
            do nR=1,n_r_max
               do lm=llm,ulm
                  w_ave(lm,nR)=dt_norm*w_ave(lm,nR)
                  z_ave(lm,nR)=dt_norm*z_ave(lm,nR)
                  p_ave(lm,nR)=dt_norm*p_ave(lm,nR)
               end do
            end do
         end if
         if ( l_heat ) then
            do nR=1,n_r_max
               do lm=llm,ulm
                  s_ave(lm,nR)=dt_norm*s_ave(lm,nR)
               end do
            end do
         end if
         if ( l_chemical_conv ) then
            do nR=1,n_r_max
               do lm=llm,ulm
                  xi_ave(lm,nR)=dt_norm*xi_ave(lm,nR)
               end do
            end do
         end if
         if ( l_mag ) then
            do nR=1,n_r_max
               do lm=llm,ulm
                  b_ave(lm,nR) =dt_norm*b_ave(lm,nR)
                  aj_ave(lm,nR)=dt_norm*aj_ave(lm,nR)
               end do
            end do
         end if
         if ( l_cond_ic ) then
            do nR=1,n_r_ic_max
               do lm=llm,ulm
                  b_ic_ave(lm,nR) =dt_norm*b_ic_ave(lm,nR)
                  aj_ic_ave(lm,nR)=dt_norm*aj_ic_ave(lm,nR)
               end do
            end do
         end if

         !----- Get the radial derivatives:
         lmStart=lmStartB(rank+1)
         lmStop = lmStopB(rank+1)

         call get_dr(w_ave,dw_ave,ulm-llm+1,lmStart-llm+1,lmStop-llm+1,   &
              &      n_r_max,rscheme_oc,nocopy=.true.)
         if ( l_mag ) then
            call get_dr(b_ave,db_ave,ulm-llm+1,lmStart-llm+1,lmStop-llm+1,  &
                 &      n_r_max,rscheme_oc,nocopy=.true.)
         end if
         if ( l_heat ) then
            call get_dr(s_ave,ds_ave,ulm-llm+1,lmStart-llm+1,lmStop-llm+1,  &
                 &      n_r_max,rscheme_oc,nocopy=.true.)
         end if
         if ( l_chemical_conv ) then
            call get_dr(xi_ave,dxi_ave,ulm-llm+1,lmStart-llm+1,lmStop-llm+1, &
                 &      n_r_max,rscheme_oc,nocopy=.true.)
         end if
         if ( l_cond_ic ) then
            call get_ddrNS_even(b_ic_ave,db_ic_ave,ddb_ic_ave,         &
                 &              ulm-llm+1,lmStart-llm+1,               &
                 &              lmStop-llm+1,n_r_ic_max,               &
                 &              n_cheb_ic_max,dr_fac_ic,workA_LMloc,   &
                 &              chebt_ic, chebt_ic_even)
            call get_drNS_even(aj_ic_ave,dj_ic_ave,                    &
                 &             ulm-llm+1,lmStart-llm+1,                &
                 &             lmStop-llm+1,n_r_ic_max,                &
                 &             n_cheb_ic_max,dr_fac_ic,workA_LMloc,    &
                 &             chebt_ic,chebt_ic_even)
         end if

         !----- Get averaged spectra:
         !      Note: average spectra will be in file no 0
         n_spec=0
         call spectrum(n_spec,time,.false.,nAve,l_stop_time,time_passed, &
              &        time_norm,w_ave,dw_ave,z_ave,b_ave,db_ave,        &
              &        aj_ave,b_ic_ave,db_ic_ave,aj_ic_ave)

         if ( l_heat ) then
            call spectrum_temp(n_spec,time,.false.,0,l_stop_time,     &
                 &             0.0_cp,0.0_cp,s_ave,ds_ave)
         end if
         if ( l_save_out ) then
            open(newunit=n_log_file, file=log_file, status='unknown', &
            &    position='append')
         end if

         !----- Write averaged energies into log-file at end of run:
         if ( l_stop_time ) then
            !----- Calculate energies of averaged field:
            n_e_sets=1
            call get_e_kin(time,.false.,.true.,n_e_sets, &
                 &         w_ave,dw_ave,z_ave,           &
                 &         e_kin_p_ave,e_kin_t_ave,      &
                 &         e_kin_p_as_ave,e_kin_t_as_ave)

            call get_e_mag(time,.false.,.true.,n_e_sets,                  &
                 &         b_ave,db_ave,aj_ave,                           &
                 &         b_ic_ave,db_ic_ave,aj_ic_ave,                  &
                 &         e_mag_p_ave,e_mag_t_ave,                       &
                 &         e_mag_p_as_ave,e_mag_t_as_ave,                 &
                 &         e_mag_p_ic_ave,e_mag_t_ic_ave,                 &
                 &         e_mag_p_as_ic_ave,e_mag_t_as_ic_ave,           &
                 &         e_mag_os_ave,e_mag_as_os_ave,e_cmb,Dip,DipCMB, &
                 &         elsAnel)

            if ( rank == 0 ) then
               !----- Output of energies of averaged field:
               write(n_log_file,'(/,A)')                           &
               &    ' ! ENERGIES OF TIME AVERAGED FIELD'
               write(n_log_file,                                   &
               &    '('' !  (total,poloidal,toroidal,total density)'')')
               write(n_log_file,'(1P,'' !  Kinetic energies:'',4ES16.6)') &
               &    (e_kin_p_ave+e_kin_t_ave),e_kin_p_ave,e_kin_t_ave,    &
               &    (e_kin_p_ave+e_kin_t_ave)/vol_oc
               write(n_log_file,'(1P,'' !  OC Mag  energies:'',4ES16.6)') &
               &    (e_mag_p_ave+e_mag_t_ave),e_mag_p_ave,e_mag_t_ave,    &
               &    (e_mag_p_ave+e_mag_t_ave)/vol_oc
               write(n_log_file,'(1P,'' !  IC Mag  energies:'',4ES16.6)') &
               &    (e_mag_p_ic_ave+e_mag_t_ic_ave),e_mag_p_ic_ave,       &
               &     e_mag_t_ic_ave,(e_mag_p_ic_ave+e_mag_t_ic_ave)/vol_ic
               write(n_log_file,'(1P,'' !  OS Mag  energies:'',ES16.6)')  &
               &     e_mag_os_ave
               write(n_log_file,'(/,'' !  AXISYMMETRIC PARTS:'')')
               write(n_log_file,                                          &
               &     '('' !  (total,poloidal,toroidal,total density)'')')
               write(n_log_file,'(1P,'' !  Kinetic AS energies:'',4ES16.6)') &
               &    (e_kin_p_as_ave+e_kin_t_as_ave),e_kin_p_as_ave,          &
               &     e_kin_t_as_ave,(e_kin_p_as_ave+e_kin_t_as_ave)/vol_oc
               write(n_log_file,'(1P,'' !  OC Mag  AS energies:'',4ES16.6)') &
               &    (e_mag_p_as_ave+e_mag_t_as_ave),e_mag_p_as_ave,          &
               &     e_mag_t_as_ave,(e_mag_p_as_ave+e_mag_t_as_ave)/vol_oc
               write(n_log_file,'(1P,'' !  IC Mag  AS energies:'',4ES16.6)') &
               &    (e_mag_p_as_ic_ave+e_mag_t_as_ic_ave),e_mag_p_as_ic_ave, &
               &     e_mag_t_as_ic_ave,(e_mag_p_as_ic_ave+e_mag_t_as_ic_ave) &
               &     /vol_ic
               write(n_log_file,'(1P,'' !  OC Mag  AS energies:'',ES16.6)')  &
               &     e_mag_os_ave
               write(n_log_file,'(1P,'' !  Relative ax. dip. E:'',ES16.6)')  &
               &     Dip
            end if
         end if ! End of run ?

         !----- Construct name of graphic file and open it:
         ! For the graphic file of the average fields, we gather them
         ! on rank 0 and use the old serial output routine.

         if ( rank == 0 ) then
            graph_file='G_ave.'//tag
            open(newunit=n_graph_file, file=graph_file, status='unknown', &
            &    form='unformatted')

            !----- Write header into graphic file:
            lGraphHeader=.true.
            call graphOut(time,0,Vr,Vt,Vp,Br,Bt,Bp,Sr,PreR,Xir, &
                 &        0,sizeThetaB,lGraphHeader)
         end if

         !-- This will be needed for the inner core
         if ( l_mag ) then
            call gather_from_lo_to_rank0(b_ave(llm,n_r_icb),bICB)
         end if

         !----- Outer core:
         do nR=1,n_r_max
            if ( l_mag ) then
               call gather_from_lo_to_rank0(b_ave(llm,nR),b_ave_global)
               call gather_from_lo_to_rank0(db_ave(llm,nR),db_ave_global)
               call gather_from_lo_to_rank0(aj_ave(llm,nR),aj_ave_global)
            end if
            call gather_from_lo_to_rank0(w_ave(llm,nR),w_ave_global)
            call gather_from_lo_to_rank0(dw_ave(llm,nR),dw_ave_global)
            call gather_from_lo_to_rank0(z_ave(llm,nR),z_ave_global)
            call gather_from_lo_to_rank0(p_ave(llm,nR),p_ave_global)
            if ( l_heat ) then
               call gather_from_lo_to_rank0(s_ave(llm,nR),s_ave_global)
            end if
            if ( l_chemical_conv ) then
               call gather_from_lo_to_rank0(xi_ave(llm,nR),xi_ave_global)
            end if

            if ( rank == 0 ) then
#ifdef WITH_SHTNS
               if ( l_mag ) then
                  call torpol_to_spat(b_ave_global, db_ave_global, &
                       &              aj_ave_global, Br, Bt, Bp)
               end if
               call torpol_to_spat(w_ave_global, dw_ave_global, &
                    &              z_ave_global, Vr, Vt, Vp)
               call scal_to_spat(p_ave_global, Prer)
               call scal_to_spat(s_ave_global, Sr)
               if ( l_chemical_conv ) then
                  call scal_to_spat(xi_ave_global, Xir)
               end if
               call graphOut(time, nR, Vr, Vt, Vp, Br, Bt, Bp, Sr, Prer, &
               &             Xir, 1, sizeThetaB, lGraphHeader)
#else
               if ( l_mag ) then
                  call legPrep(b_ave_global,db_ave_global,db_ave_global, &
                       &       aj_ave_global,aj_ave_global,dLh,lm_max,   &
                       &       l_max,minc,r(nR),.false.,.true.,          &
                       &       dLhb,bhG,bhC,dLhb,bhG,bhC)
               end if
               call legPrep(w_ave_global,dw_ave_global,dw_ave_global, &
                    &       z_ave_global,z_ave_global,dLh,lm_max,     &
                    &       l_max,minc,r(nR),.false.,.true.,          &
                    &       dLhw,vhG,vhC,dLhb,bhG,bhC)

               do nThetaB=1,nThetaBs
                  nThetaStart=(nThetaB-1)*sizeThetaB+1

                  !-------- Transform to grid space:
                  call legTF(dLhb,bhG,bhC,dLhw,vhG,vhC,                  &
                       &     l_max,minc,nThetaStart,sizeThetaB,          &
                       &     Plm,dPlm,.true.,.false.,                    &
                       &     Br,Bt,Bp,Br,Br,Br)
                  call legTF(dLhw,vhG,vhC,dLhw,vhG,vhC,                  &
                       &     l_max,minc,nThetaStart,sizeThetaB,          &
                       &     Plm,dPlm,.true.,.false.,                    &
                       &     Vr,Vt,Vp,Br,Br,Br)
                  call legTF(s_ave_global,vhG,vhC,dLhw,vhG,vhC,          &
                       &     l_max,minc,nThetaStart,sizeThetaB,          &
                       &     Plm,dPlm,.false.,.false.,                   &
                       &     Sr,Vt,Vp,Br,Br,Br)
                  if ( l_chemical_conv ) then
                     call legTF(xi_ave_global,vhG,vhC,dLhw,vhG,vhC,      &
                          &     l_max,minc,nThetaStart,sizeThetaB,       &
                          &     Plm,dPlm,.false.,.false.,                &
                          &     Xir,Vt,Vp,Br,Br,Br)
                     if ( .not. l_axi ) call fft_thetab(Xir,1)
                  end if
                  call legTF(p_ave_global,vhG,vhC,dLhw,vhG,vhC,          &
                       &     l_max,minc,nThetaStart,sizeThetaB,          &
                       &     Plm,dPlm,.false.,.false.,                   &
                       &     Prer,Vt,Vp,Br,Br,Br)
                  if ( .not. l_axi ) then
                     call fft_thetab(Br,1)
                     call fft_thetab(Bp,1)
                     call fft_thetab(Bt,1)
                     call fft_thetab(Vr,1)
                     call fft_thetab(Vt,1)
                     call fft_thetab(Vp,1)
                     call fft_thetab(Sr,1)
                     call fft_thetab(Prer,1)
                  end if

                  !-------- Graphic output:
                  call graphOut(time,nR,Vr,Vt,Vp,Br,Bt,Bp,Sr,Prer, &
                       &        Xir,nThetaStart,sizeThetaB,lGraphHeader)
               end do
#endif
            end if
         end do

         !----- Inner core: Transform is included in graphOut_IC!
         if ( l_mag .and. n_r_ic_max > 0 ) then
            call gather_all_from_lo_to_rank0(gt_IC,b_ic_ave,b_ic_ave_global)
            call gather_all_from_lo_to_rank0(gt_IC,db_ic_ave,db_ic_ave_global)
            call gather_all_from_lo_to_rank0(gt_IC,ddb_ic_ave,ddb_ic_ave_global)
            call gather_all_from_lo_to_rank0(gt_IC,aj_ic_ave,aj_ic_ave_global)
            call gather_all_from_lo_to_rank0(gt_IC,dj_ic_ave,dj_ic_ave_global)

            if ( rank == 0 ) then
               call graphOut_IC(b_ic_ave_global,db_ic_ave_global,   &
                    &           ddb_ic_ave_global,aj_ic_ave_global, &
                    &           dj_ic_ave_global,bICB,l_avg=.true.)
            end if
         end if

         if ( rank == 0 ) close(n_graph_file)  ! close graphic output file !

         !----- Write info about graph-file into STDOUT and log-file:
         if ( l_stop_time ) then
            if ( rank == 0 )  &
            &  write(n_log_file,'(/,'' ! WRITING AVERAGED GRAPHIC FILE !'')')
         end if

         !--- Store time averaged poloidal magnetic coeffs at cmb
         if ( l_mag) then
            outFile='B_coeff_cmb_ave.'//tag
            nOut   =93
            n_cmb_sets=-1
            !call write_Bcmb(time,b(1,n_r_cmb),lm_max,l_max,           &
            !     &           l_max_cmb,minc,lm2,n_cmb_sets,outFile,nOut)
            call write_Bcmb(time,b_ave(:,n_r_cmb),l_max_cmb,n_cmb_sets,outFile,nOut)
         end if

         !--- Store potentials of averaged field:
         !    dw_ave and db_ave used as work arrays here.
         nBpotSets=-1
         nVpotSets=-1
         nTpotSets=-1
         if ( l_mag) then
            call write_Pot(time,b_ave,aj_ave,b_ic_ave,aj_ic_ave,nBpotSets,  &
                 &        'B_lmr_ave.',omega_ma,omega_ic)
         end if
         call write_Pot(time,w_ave,z_ave,b_ic_ave,aj_ic_ave,nVpotSets,      &
              &        'V_lmr_ave.',omega_ma,omega_ic)
         if ( l_heat ) then
            call write_Pot(time,s_ave,z_ave,b_ic_ave,aj_ic_ave,nTpotSets,   &
                 &        'T_lmr_ave.',omega_ma,omega_ic)
         end if

         !--- Store checkpoint file
         call store(simtime,dt,dtNew,-1,l_stop_time,.false.,.true.,         &
              &     w_ave,z_ave,p_ave,s_ave,xi_ave,b_ave,aj_ave,b_ic_ave,   &
              &     aj_ic_ave,dwdtLast_LMloc,dzdtLast_lo,dpdtLast_LMloc,    &
              &     dsdtLast_LMloc,dxidtLast_LMloc,dbdtLast_LMloc,          &
              &     djdtLast_LMloc,dbdt_icLast_LMloc,djdt_icLast_LMloc)

         ! if ( l_chemical_conv ) then
            ! call write_Pot(time,s_ave,z_ave,b_ic_ave,aj_ic_ave,nTpotSets,   &
                 ! &        'Xi_lmr_ave.',omega_ma,omega_ic)
         ! end if

         if ( l_save_out ) close(n_log_file)

         ! now correct the stored average fields by the factor which has been
         ! applied before
         if ( l_conv ) then
            do nR=1,n_r_max
               do lm=llm,ulm
                  w_ave(lm,nR)=w_ave(lm,nR)*time_norm
                  z_ave(lm,nR)=z_ave(lm,nR)*time_norm
                  p_ave(lm,nR)=p_ave(lm,nR)*time_norm
               end do
            end do
         end if
         if ( l_heat ) then
            do nR=1,n_r_max
               do lm=llm,ulm
                  s_ave(lm,nR)=s_ave(lm,nR)*time_norm
               end do
            end do
         end if
         if ( l_chemical_conv ) then
            do nR=1,n_r_max
               do lm=llm,ulm
                  xi_ave(lm,nR)=xi_ave(lm,nR)*time_norm
               end do
            end do
         end if
         if ( l_mag ) then
            do nR=1,n_r_max
               do lm=llm,ulm
                  b_ave(lm,nR) =b_ave(lm,nR)*time_norm
                  aj_ave(lm,nR)=aj_ave(lm,nR)*time_norm
               end do
            end do
         end if
         if ( l_cond_ic ) then
            do nR=1,n_r_ic_max
               do lm=llm,ulm
                  b_ic_ave(lm,nR) =b_ic_ave(lm,nR)*time_norm
                  aj_ic_ave(lm,nR)=aj_ic_ave(lm,nR)*time_norm
               end do
            end do
         end if


      end if ! last time step ?

   end subroutine fields_average
!------------------------------------------------------------------------------
end module fields_average_mod
