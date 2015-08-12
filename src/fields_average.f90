!$Id$
module fields_average_mod

   use truncation
   use precision_mod, only: cp
   use radial_data, only: n_r_cmb
   use radial_functions, only: i_costf_init, d_costf_init, drx,    &
                               i_costf1_ic_init, d_costf1_ic_init, &
                               i_costf2_ic_init, d_costf2_ic_init, &
                               r, dr_fac_ic
   use blocking,only: lmStartB, lmStopB, sizeThetaB, nThetaBs, lm2
   use horizontal_data, only: Plm, dPlm, dLh
   use logic, only: l_mag, l_conv, l_save_out, l_heat, l_cond_ic
   use kinetic_energy, only: get_e_kin
   use magnetic_energy, only: get_e_mag
   use output_data, only: tag, graph_file, nLF, n_graph_file, ngform, &
                          log_file, n_graphs, l_max_cmb
   use parallel_mod, only: rank
#if (FFTLIB==JW)
   use fft_JW
#elif (FFTLIB==MKL)
   use fft_MKL
#endif
   use const, only: zero, vol_oc, vol_ic, one
   use LMLoop_data, only: llm,ulm,llmMag,ulmMag
   use communications, only: get_global_sum, gather_from_lo_to_rank0,&
                           & gather_all_from_lo_to_rank0,gt_OC,gt_IC
   use out_coeff, only: write_Bcmb
   use spectra, only: spectrum, spectrum_temp
   use graphOut_mod, only: graphOut, graphOut_IC
   use store_pot_mod, only: storePotW
   use leg_helper_mod, only: legPrep
   use legendre_spec_to_grid, only: legTF
   use radial_der_even, only: get_drNS_even, get_ddrNS_even
   use radial_der, only: get_drNS
 
   implicit none
 
   private
 
   complex(cp), allocatable :: w_ave(:,:)
   complex(cp), allocatable :: z_ave(:,:)
   complex(cp), allocatable :: s_ave(:,:)
   complex(cp), allocatable :: b_ave(:,:)
   complex(cp), allocatable :: aj_ave(:,:)
   complex(cp), allocatable :: b_ic_ave(:,:)
   complex(cp), allocatable :: aj_ic_ave(:,:)
 
   ! on rank 0 we also allocate the following fields
   complex(cp), allocatable :: db_ave_global(:),aj_ave_global(:)
   complex(cp), allocatable :: w_ave_global(:),dw_ave_global(:)
   complex(cp), allocatable :: z_ave_global(:), s_ave_global(:)
 
   public :: initialize_fields_average_mod, fields_average

contains

   subroutine initialize_fields_average_mod

      allocate( w_ave(llm:ulm,n_r_max) )
      allocate( z_ave(llm:ulm,n_r_max) )
      allocate( s_ave(llm:ulm,n_r_max) )
      allocate( b_ave(llm:ulm,n_r_max) )
      allocate( aj_ave(llm:ulm,n_r_max) )
      allocate( b_ic_ave(llm:ulm,n_r_ic_max) )
      allocate( aj_ic_ave(llm:ulm,n_r_ic_max) )

      if ( rank == 0 ) then
         allocate( db_ave_global(1:lm_max) )
         allocate( aj_ave_global(1:lm_max) )
         allocate( w_ave_global(1:lm_max) )
         allocate( dw_ave_global(1:lm_max) )
         allocate( z_ave_global(1:lm_max) )
         allocate( s_ave_global(1:lm_max) )
#ifdef WITH_DEBUG
      else
         allocate( db_ave_global(1) )
         allocate( aj_ave_global(1) )
         allocate( w_ave_global(1) )
         allocate( dw_ave_global(1) )
         allocate( z_ave_global(1) )
         allocate( s_ave_global(1) )
#endif
      end if

   end subroutine initialize_fields_average_mod
!----------------------------------------------------------------------------
   subroutine fields_average(nAve,l_stop_time,                        &
      &                      time_passed,time_norm,omega_ic,omega_ma, &
      &                      w,z,s,b,aj,b_ic,aj_ic)
      !  +-------------+----------------+------------------------------------+
      !  |                                                                   |
      !  |  This subroutine averages fields b and v over time.               |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+

      !-- Input of variables:
      integer,     intent(in) :: nAve         ! number for averaged time steps
      logical,     intent(in) :: l_stop_time  ! true if this is the last time step
      real(cp),    intent(in) :: time_passed  ! time passed since last log
      real(cp),    intent(in) :: time_norm    ! time passed since start of time loop
      real(cp),    intent(in) :: omega_ic,omega_ma
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(in) :: z(llm:ulm,n_r_max)
      complex(cp), intent(in) :: s(llm:ulm,n_r_max)
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
      complex(cp) :: b_ave_global(1:lm_maxMag,n_r_maxMag)

      !----- Time averaged fields:
      complex(cp) :: dw_ave(llm:ulm,n_r_max)
      complex(cp) :: ds_ave(llm:ulm,n_r_max)
      complex(cp) :: db_ave(llm:ulm,n_r_max)
      complex(cp) :: db_ic_ave(llm:ulm,n_r_ic_max)
      complex(cp) :: ddb_ic_ave(llm:ulm,n_r_ic_max)
      complex(cp) :: dj_ic_ave(llm:ulm,n_r_ic_max)

      !----- Work array:
      complex(cp) :: workA_LMloc(llm:ulm,n_r_max)

      !----- Fields in grid space:
      real(cp) :: Br(nrp,nfs),Bt(nrp,nfs),Bp(nrp,nfs) ! B field comp.
      real(cp) :: Vr(nrp,nfs),Vt(nrp,nfs),Vp(nrp,nfs) ! B field comp.
      real(cp) :: Sr(nrp,nfs)                         ! entropy

      !----- Help arrays for fields:
      complex(cp) :: dLhb(lm_max),bhG(lm_max),bhC(lm_max)
      complex(cp) :: dLhw(lm_max),vhG(lm_max),vhC(lm_max)

      !----- Energies of time average field:
      real(cp) :: ekinR(n_r_max)     ! kinetic energy w radius
      real(cp) :: e_kin_p_ave,e_kin_t_ave
      real(cp) :: e_kin_p_as_ave,e_kin_t_as_ave
      real(cp) :: e_mag_p_ave,e_mag_t_ave
      real(cp) :: e_mag_p_as_ave,e_mag_t_as_ave
      real(cp) :: e_mag_p_ic_ave,e_mag_t_ic_ave
      real(cp) :: e_mag_p_as_ic_ave,e_mag_t_as_ic_ave
      real(cp) :: e_mag_os_ave,e_mag_as_os_ave
      real(cp) :: Dip,DipCMB,e_cmb,elsAnel

      integer :: lm,nR,nThetaB,nThetaStart
      integer :: n_e_sets,n_spec

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
            end if
            if ( l_heat ) then
               s_ave=zero
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

         call get_drNS(w_ave,dw_ave,ulm-llm+1,          &
              &        lmStart-llm+1,lmStop-llm+1,      &
              &        n_r_max,n_cheb_max,workA_LMloc,  &
              &        i_costf_init,d_costf_init,drx)
         if (l_mag) then
            call get_drNS(b_ave,db_ave,ulm-llm+1,          &
                 &        lmStart-llm+1,lmStop-llm+1,      &
                 &        n_r_max,n_cheb_max,workA_LMloc,  &
                 &        i_costf_init,d_costf_init,drx)
         end if
         if ( l_heat ) then
            call get_drNS(s_ave,ds_ave,ulm-llm+1,            &
                 &        lmStart-llm+1,lmStop-llm+1,        &
                 &        n_r_max,n_cheb_max,workA_LMloc,    &
                 &        i_costf_init,d_costf_init,drx)
         end if
         if ( l_cond_ic ) then
            call get_ddrNS_even(b_ic_ave,db_ic_ave,ddb_ic_ave,         &
                 &              ulm-llm+1,lmStart-llm+1,               &
                 &              lmStop-llm+1,n_r_ic_max,               &
                 &              n_cheb_ic_max,dr_fac_ic,workA_LMloc,   &
                 &              i_costf1_ic_init,d_costf1_ic_init,     &
                 &              i_costf2_ic_init,d_costf2_ic_init)
            call get_drNS_even(aj_ic_ave,dj_ic_ave,                    &
                 &             ulm-llm+1,lmStart-llm+1,                &
                 &             lmStop-llm+1,n_r_ic_max,                &
                 &             n_cheb_ic_max,dr_fac_ic,workA_LMloc,    &
                 &             i_costf1_ic_init,d_costf1_ic_init,      &
                 &             i_costf2_ic_init,d_costf2_ic_init)
         end if

         !----- Get averaged spectra:
         !      Note: average spectra will be in file no 0
         n_spec=0
         call spectrum(time,n_spec,w_ave,dw_ave,z_ave, &
              &        b_ave,db_ave,aj_ave,            &
              &        b_ic_ave,db_ic_ave,aj_ic_ave)  

         if ( l_heat ) then
            call spectrum_temp(time,n_spec,s_ave,ds_ave)
         end if
         if ( l_save_out ) then
            open(nLF,file=log_file, status='unknown', position='append')
         end if

         !----- Write averaged energies into log-file at end of run:
         if ( l_stop_time ) then 
            !----- Calculate energies of averaged field:
            n_e_sets=1
            call get_e_kin(time,.false.,.true.,n_e_sets, &
                 &         w_ave,dw_ave,z_ave,           &
                 &         e_kin_p_ave,e_kin_t_ave,      &
                 &         e_kin_p_as_ave,e_kin_t_as_ave,&
                 &         eKinR)

            call get_e_mag(time,.false.,.true.,n_e_sets,                  &
                 &         b_ave,db_ave,aj_ave,                           &
                 &         b_ic_ave,db_ic_ave,aj_ic_ave,                  &
                 &         e_mag_p_ave,e_mag_t_ave,                       &
                 &         e_mag_p_as_ave,e_mag_t_as_ave,                 &
                 &         e_mag_p_ic_ave,e_mag_t_ic_ave,                 &
                 &         e_mag_p_as_ic_ave,e_mag_t_as_ic_ave,           &
                 &         e_mag_os_ave,e_mag_as_os_ave,e_cmb,Dip,DipCMB, &
                 &         elsAnel)

            if (rank == 0) then
               !----- Output of energies of averaged field:
               write(nLF,'(/,A)')                                        &
                    &           ' ! ENERGIES OF TIME AVERAGED FIELD'
               write(nLF,                                                &
                    &        '('' !  (total,poloidal,toroidal,total density)'')')
               write(nLF,'(1P,'' !  Kinetic energies:'',4D16.6)')                  &
                    &           (e_kin_p_ave+e_kin_t_ave),e_kin_p_ave,e_kin_t_ave, &
                    &           (e_kin_p_ave+e_kin_t_ave)/vol_oc
               write(nLF,'(1P,'' !  OC Mag  energies:'',4D16.6)')                  &
                    &           (e_mag_p_ave+e_mag_t_ave),e_mag_p_ave,e_mag_t_ave, &
                    &           (e_mag_p_ave+e_mag_t_ave)/vol_oc
               write(nLF,'(1P,'' !  IC Mag  energies:'',4D16.6)')        &
                    &           (e_mag_p_ic_ave+e_mag_t_ic_ave),         &
                    &           e_mag_p_ic_ave,e_mag_t_ic_ave,           &
                    &           (e_mag_p_ic_ave+e_mag_t_ic_ave)/vol_ic
               write(nLF,'(1P,'' !  OS Mag  energies:'',D16.6)')         &
                    &           e_mag_os_ave
               write(nLF,'(/,'' !  AXISYMMETRIC PARTS:'')')
               write(nLF,                                                &
                    &        '('' !  (total,poloidal,toroidal,total density)'')')
               write(nLF,'(1P,'' !  Kinetic AS energies:'',4D16.6)')     &
                    &           (e_kin_p_as_ave+e_kin_t_as_ave),         &
                    &           e_kin_p_as_ave,e_kin_t_as_ave,           &
                    &           (e_kin_p_as_ave+e_kin_t_as_ave)/vol_oc
               write(nLF,'(1P,'' !  OC Mag  AS energies:'',4D16.6)')     &
                    &           (e_mag_p_as_ave+e_mag_t_as_ave),         &
                    &           e_mag_p_as_ave,e_mag_t_as_ave,           &
                    &           (e_mag_p_as_ave+e_mag_t_as_ave)/vol_oc
               write(nLF,'(1P,'' !  IC Mag  AS energies:'',4D16.6)')     &
                    &           (e_mag_p_as_ic_ave+e_mag_t_as_ic_ave),   &
                    &           e_mag_p_as_ic_ave,e_mag_t_as_ic_ave,     &
                    &           (e_mag_p_as_ic_ave+e_mag_t_as_ic_ave)/vol_ic
               write(nLF,'(1P,'' !  OS Mag  AS energies:'',D16.6)')      &
                    &           e_mag_os_ave
               write(nLF,'(1P,'' !  Relative ax. dip. E:'',D16.6)')      &
                    &           Dip           
            end if
         end if ! End of run ?
            
         !----- Construct name of graphic file and open it:
         ! For the graphic file of the average fields, we gather them
         ! on rank 0 and use the old serial output routine.

         if ( rank == 0 ) then
            if ( ngform == 0 ) then
               graph_file='G_ave.'//tag
               open(n_graph_file, file=graph_file, status='unknown', form='unformatted')
            else
               graph_file='g_ave.'//tag
               open(n_graph_file, file=graph_file, status='unknown', form='formatted')
            end if

            !----- Write header into graphic file:
            lGraphHeader=.true.
            call graphOut(time,0,ngform,Vr,Vt,Vp,Br,Bt,Bp,Sr,0,sizeThetaB,lGraphHeader)
         end if

         !----- Transform and output of data:
         ! b_ave is different as it is again used later for graphOut_IC
         if ( l_mag ) then
            call gather_all_from_lo_to_rank0(gt_OC,b_ave,b_ave_global)
         end if
         !----- Outer core:
         do nR=1,n_r_max
            if ( l_mag ) then
               call gather_from_lo_to_rank0(db_ave(llm,nR),db_ave_global)
               call gather_from_lo_to_rank0(aj_ave(llm,nR),aj_ave_global)
            end if
            call gather_from_lo_to_rank0(w_ave(llm,nR),w_ave_global)
            call gather_from_lo_to_rank0(dw_ave(llm,nR),dw_ave_global)
            call gather_from_lo_to_rank0(z_ave(llm,nR),z_ave_global)
            if ( l_heat ) then
               call gather_from_lo_to_rank0(s_ave(llm,nR),s_ave_global)
            end if

            if ( rank == 0 ) then
               if ( l_mag ) then
                  call legPrep(b_ave_global(1,nR),db_ave_global,db_ave_global, &
                       &       aj_ave_global,aj_ave_global,dLh,lm_max,         &
                       &       l_max,minc,r(nR),.false.,.true.,                &
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
                  call fft_thetab(Br,1)
                  call fft_thetab(Bp,1)
                  call fft_thetab(Bt,1)
                  call fft_thetab(Vr,1)
                  call fft_thetab(Vt,1)
                  call fft_thetab(Vp,1)
                  call fft_thetab(Sr,1)

                  !-------- Graphic output:
                  call graphOut(time,nR,ngform,Vr,Vt,Vp,Br,Bt,Bp,Sr, &
                       &        nThetaStart,sizeThetaB,lGraphHeader)
               end do
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
               call graphOut_IC(ngform,b_ic_ave_global,db_ic_ave_global,&
                    &           ddb_ic_ave_global,aj_ic_ave_global,     &
                    &           dj_ic_ave_global,b_ave_global)
            end if
         end if

         if ( rank == 0 ) close(n_graph_file)  ! close graphic output file !

         !----- Write info about graph-file into STDOUT and log-file:
         if ( l_stop_time ) then
            if ( rank == 0 ) write(nLF,'(/,'' ! WRITING AVERAGED GRAPHIC FILE !'')')
         end if

         !--- Store time averaged poloidal magnetic coeffs at cmb
         if ( rank == 0 ) then
            if ( l_mag) then
               outFile='B_coeff_cmb_ave.'//tag
               nOut   =93
               n_cmb_sets=-1
               !call write_Bcmb(time,b(1,n_r_cmb),lm_max,l_max,           &
               !     &           l_max_cmb,minc,lm2,n_cmb_sets,outFile,nOut)
               call write_Bcmb(time,b_ave_global(1,n_r_cmb),1,lm_max,l_max, &
                    &          l_max_cmb,minc,lm2,n_cmb_sets,outFile,nOut)
            end if
         end if

         !--- Store potentials of averaged field:
         !    dw_ave and db_ave used as work arrays here.
         nBpotSets=-1
         nVpotSets=-1
         nTpotSets=-1
         if ( l_mag) then
            call storePotW(time,b_ave,aj_ave,b_ic_ave,aj_ic_ave,            &
                 &         workA_LMloc,dw_ave,db_ave,nBpotSets,'Bpot_ave.', &
                 &         omega_ma,omega_ic)
         end if
         call storePotW(time,w_ave,z_ave,b_ic_ave,aj_ic_ave,             &
              &         workA_LMloc,dw_ave,db_ave,nVpotSets,'Vpot_ave.', &
              &                                      omega_ma,omega_ic)
         call storePotW(time,s_ave,z_ave,b_ic_ave,aj_ic_ave,             &
              &         workA_LMloc,dw_ave,db_ave,nTpotSets,'Tpot_ave.', &
              &                                      omega_ma,omega_ic)

         if ( l_save_out ) close(nLF)

         ! now correct the stored average fields by the factor which has been
         ! applied before
         if ( l_conv ) then
            do nR=1,n_r_max
               do lm=llm,ulm
                  w_ave(lm,nR)=w_ave(lm,nR)*time_norm
                  z_ave(lm,nR)=z_ave(lm,nR)*time_norm
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
