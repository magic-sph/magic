#include "perflib_preproc.cpp"
module start_fields

#ifdef WITH_MPI
   use mpi
#endif
   use truncation
   use precision_mod
   use radial_data, only: n_r_cmb, n_r_icb
   use radial_functions, only: topcond, botcond, chebt_oc,         &
                               drx, ddrx, dr_fac_ic, chebt_ic,     &
                               chebt_ic_even, r, or1, alpha0,      &
                               dLtemp0, dLalpha0, beta, orho1,     &
                               temp0
   use physical_parameters, only: interior_model, epsS, impS, n_r_LCR,   &
                                  ktopv, kbotv, LFfac, imagcon, ThExpNb, &
                                  ViscHeatFac
   use num_param, only: dtMax, alpha
   use special, only: lGrenoble
   use blocking, only: lmStartB, lmStopB, nLMBs, lo_map
   use logic, only: l_conv, l_mag, l_cond_ic, l_heat, l_SRMA, l_SRIC,    &
                    l_mag_kin, l_mag_LF, l_rot_ic, l_z10Mat, l_LCR,      &
                    l_rot_ma, l_temperature_diff, l_single_matrix
   use init_fields, only: l_start_file, init_s1, init_b1, tops, ps_cond, &
                          initV, initS, initB, s_cond, start_file
   use fields ! The entire module is required
   use fieldsLast ! The entire module is required
   use constants, only: zero, c_lorentz_ma, c_lorentz_ic, osq4pi, &
                    one, two
   use useful, only: cc2real, logWrite
   use LMLoop_data, only: lm_per_rank, lm_on_last_rank, llm, ulm, &
                          ulmMag,llmMag
   use parallel_mod, only: rank, n_procs, nLMBs_per_rank
   use communications, only: lo2r_redist_start,lo2r_s,lo2r_z, lo2r_p,    &
                             lo2r_b, lo2r_aj, scatter_from_rank0_to_lo,  &
                             get_global_sum, lo2r_w
   use radial_der, only: get_dr, get_ddr
   use radial_der_even, only: get_ddr_even
#ifdef WITH_HDF5
   use readCheckPoints, only: readHdf5_serial,readStartFields
#else
   use readCheckPoints, only: readStartFields
#endif
    
   implicit none

   private

   public :: getStartFields

contains

   subroutine getStartFields(time,dt,dtNew,n_time_step)
      !
      !  Purpose of this subroutine is to initialize the fields and       
      !  other auxiliary parameters.                                      
      !
    
      !---- Output variables:
      real(cp), intent(out) :: time,dt,dtNew
      integer,  intent(out) :: n_time_step
    
      !-- Local variables:
      integer :: nR,l1m0,nLMB,l,m
      integer :: lm, n_r
      integer :: lmStart,lmStop
      real(cp) :: coex
      real(cp) :: d_omega_ma_dt,d_omega_ic_dt
      character(len=76) :: message
    
      real(cp) :: sEA,sES,sAA
    
      real(cp) :: s0(n_r_max),p0(n_r_max),ds0(n_r_max),dp0(n_r_max)
      real(cp) :: w1(n_r_max),w2(n_r_max)
    
      complex(cp), allocatable :: workA_LMloc(:,:),workB_LMloc(:,:)
    
      integer :: ierr
      logical :: DEBUG_OUTPUT=.false.
    
      !PERFON('getFlds')
      !print*,"Starting getStartFields"
      !write(*,"(2(A,L1))") "l_conv=",l_conv,", l_heat=",l_heat
      !---- Computations for the Nusselt number if we are anelastic
      !     Can be done before setting the fields
      if (l_heat) then
    
         if ( l_single_matrix ) then
            call ps_cond(s0,p0)
            open(unit=999, file='pscond.dat')
            do n_r=1,n_r_max
               write(999,*) r(n_r), s0(n_r)*osq4pi, p0(n_r)*osq4pi, &
                            osq4pi*temp0(n_r)*(s0(n_r)+alpha0(n_r)* &
                            orho1(n_r)*p0(n_r)*ThExpNb*ViscHeatFac)
            end do
            close(999)

            call get_dr(s0,ds0,n_r_max,n_cheb_max,w1,w2,chebt_oc,drx)

            if ( l_temperature_diff ) then
               call get_dr(p0,dp0,n_r_max,n_cheb_max,w1,w2,chebt_oc,drx)

               topcond = -osq4pi*(  temp0(1)*dLtemp0(1)*s0(1) +             &
                 &                  ViscHeatFac*ThExpNb*(                   &
                 &                  alpha0(1)*temp0(1)*orho1(1)*            &
                 &                  (dLalpha0(1)+dLtemp0(1)-beta(1))        &
                 &                  )*                  p0(1) +             &
                 &                  temp0(1)*          ds0(1) +             &
                 &                  ViscHeatFac*ThExpNb*alpha0(1)*orho1(1)* &
                 &                  temp0(1)*          dp0(1) )
               botcond = -osq4pi*( temp0(n_r_max)*dLtemp0(n_r_max)*         &
                 &                                      s0(n_r_max) +       &
                 &                  ViscHeatFac*ThExpNb*(                   &
                 &           alpha0(n_r_max)*temp0(n_r_max)*orho1(n_r_max)* &
                 & (dLtemp0(n_r_max)+dLalpha0(n_r_max)-beta(n_r_max)) )*    &
                 &                                      p0(n_r_max) +       &
                 &                  temp0(n_r_max)*    ds0(n_r_max) +       &
                 &    ViscHeatFac*ThExpNb*alpha0(n_r_max)*temp0(n_r_max)*   &
                 &                  orho1(n_r_max)*    dp0(n_r_max) )

            else ! entropy diffusion
               topcond=-osq4pi*ds0(1)
               botcond=-osq4pi*ds0(n_r_max)
            end if

         else
            if ( index(interior_model,'EARTH') /= 0 ) then
               !topcond=-one/epsS*dtemp0(1)
               !botcond=-one/epsS*dtemp0(n_r_max)
               topcond=one
               botcond=one
            else
               call s_cond(s0)
               call get_dr(s0,ds0,n_r_max,n_cheb_max,w1,w2,chebt_oc,drx)
               topcond=-osq4pi*ds0(1)
               botcond=-osq4pi*ds0(n_r_max)
            end if
         end if
      end if
    
    
      !-- Start with setting fields to zero:
      !   Touching the fields with the appropriate processor
      !   for the LM-distribute parallel region (LMLoop) makes
      !   sure that they are located close the individual
      !   processors in memory:
    
      if ( rank == 0 ) then
         if ( l_start_file ) then
            !PERFON('readFlds')
#ifdef WITH_HDF5
            if ( index(start_file,'h5_') /= 0 ) then
               call readHdf5_serial( w,dwdtLast,z,dzdtLast,p,dpdtLast,s,dsdtLast, &
                    &                b,dbdtLast,aj,djdtLast,b_ic,dbdt_icLast,     &
                    &                aj_ic,djdt_icLast,omega_ic,omega_ma,         &
                    &                lorentz_torque_icLast,lorentz_torque_maLast, &
                    &                time,dt,dtNew)
               n_time_step=0
            else
               call readStartFields( w,dwdtLast,z,dzdtLast,p,dpdtLast,s,dsdtLast, &
                    &                b,dbdtLast,aj,djdtLast,b_ic,dbdt_icLast,     &
                    &                aj_ic,djdt_icLast,omega_ic,omega_ma,         &
                    &                lorentz_torque_icLast,lorentz_torque_maLast, &
                    &                time,dt,dtNew,n_time_step)
            end if
#else
            call readStartFields( w,dwdtLast,z,dzdtLast,p,dpdtLast,s,dsdtLast, &
                 &                b,dbdtLast,aj,djdtLast,b_ic,dbdt_icLast,     &
                 &                aj_ic,djdt_icLast,omega_ic,omega_ma,         &
                 &                lorentz_torque_icLast,lorentz_torque_maLast, &
                 &                time,dt,dtNew,n_time_step)
#endif
            if ( dt > 0.0_cp ) then
               write(message,'(''! Using old time step:'',ES16.6)') dt
            else
               dt=dtMax
               write(message,'(''! Using dtMax time step:'',ES16.6)') dtMax
            end if
            !PERFOFF
         else
            ! Initialize with zero
            if ( l_conv ) then
               w       =zero
               dwdtLast=zero
               z       =zero
               dzdtLast=zero
               p       =zero
               dpdtLast=zero
            end if
            if ( l_heat ) then
               s       =zero
               dsdtLast=zero
            end if
            if ( l_mag ) then
               b       =zero
               dbdtLast=zero
               aj      =zero
               djdtLast=zero
            end if
            if ( l_cond_ic ) then
               b_ic       =zero
               dbdt_icLast=zero
               aj_ic      =zero
               djdt_icLast=zero
            end if
    
            time =0.0_cp
            dt   =dtMax
            dtNew=dtMax
            n_time_step=0
            write(message,'(''! Using dtMax time step:'',ES16.6)') dtMax
         end if
         call logWrite(message)
    
         !----- Get radial derivatives and initialize:
    
         do nLMB=1,nLMBs ! Blocking of loop over all (l,m)
            lmStart=lmStartB(nLMB)
            lmStop =lmStopB(nLMB)
    
            !----- Initialize/add magnetic field:
            if ( ( imagcon /= 0 .or. init_b1 /= 0 .or. lGrenoble ) &
                 & .and. ( l_mag .or. l_mag_LF ) ) then
               call initB(b,aj,b_ic,aj_ic, &
                    &     lorentz_torque_icLast, lorentz_torque_maLast, &
                    &     lmStart,lmStop)
            end if
    
            !----- Initialize/add velocity, set IC and ma rotation:
            if ( l_conv .or. l_mag_kin .or. l_SRIC .or. l_SRMA ) then
               call initV(w,z,omega_ic,omega_ma,lmStart,lmStop)
            end if
    
            !----- Initialize/add entropy:
            if ( ( init_s1 /= 0 .or. impS /= 0 ) .and. l_heat ) then
               call initS(s,p,lmStart,lmStop)
            end if
    
            if ( DEBUG_OUTPUT ) then
               write(*,"(A,I3,10ES22.15)") "direct after init: w,z,s,b,aj ", &
                      nLMB, sum(w), sum(z), sum(s),sum(b),sum(aj)
            end if
    
         end do ! Loop over LM blocks
    
      end if
    
      ! ========== Redistribution of the fields ============
      ! 1. Broadcast the scalars
#ifdef WITH_MPI
      call MPI_Bcast(omega_ic,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(omega_ma,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(lorentz_torque_icLast,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(lorentz_torque_maLast,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(time,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(dt,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(dtNew,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(n_time_step,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
    
      ! 2. Scatter the d?dtLast arrays, they are only used in LMLoop
      !write(*,"(4X,A)") "Start Scatter d?dtLast arrays"
      do nR=1,n_r_max
         !write(*,"(8X,A,I4)") "nR = ",nR
         call scatter_from_rank0_to_lo(dwdtLast(1,nR),dwdtLast_LMloc(llm,nR))
         call scatter_from_rank0_to_lo(dzdtLast(1,nR),dzdtLast_lo(llm,nR))
         call scatter_from_rank0_to_lo(dpdtLast(1,nR),dpdtLast_LMloc(llm,nR))
         call scatter_from_rank0_to_lo(dsdtLast(1,nR),dsdtLast_LMloc(llm,nR))
    
         if ( l_mag ) then
            call scatter_from_rank0_to_lo(dbdtLast(1,nR),dbdtLast_LMloc(llm,nR))
            call scatter_from_rank0_to_lo(djdtLast(1,nR),djdtLast_LMloc(llm,nR))
         end if
      end do
      if ( l_cond_ic ) then
         do nR=1,n_r_ic_max
            call scatter_from_rank0_to_lo(dbdt_icLast(1,nR),dbdt_icLast_LMloc(llm,nR))
            call scatter_from_rank0_to_lo(djdt_icLast(1,nR),djdt_icLast_LMloc(llm,nR))
         end do
      end if
    
      ! 3. Scatter the fields to the LMloc space
      !write(*,"(4X,A)") "Start Scatter the fields"
      if ( DEBUG_OUTPUT ) then
         if ( rank == 0 ) write(*,"(A,2ES20.12)") "init z = ",sum(z)
      end if
      do nR=1,n_r_max
         call scatter_from_rank0_to_lo(w(1,nR),w_LMloc(llm:,nR))
         call scatter_from_rank0_to_lo(z(1,nR),z_LMloc(llm:,nR))
         call scatter_from_rank0_to_lo(p(1,nR),p_LMloc(llm:,nR))
         call scatter_from_rank0_to_lo(s(1,nR),s_LMloc(llm:,nR))
         if ( l_mag ) then
            call scatter_from_rank0_to_lo(b(1,nR),b_LMloc(llmMag:,nR))
            call scatter_from_rank0_to_lo(aj(1,nR),aj_LMloc(llmMag:,nR))
         end if
         !if (DEBUG_OUTPUT) then
         !   if (rank == 0) then
         !      write(*,"(A,I4,6ES22.14)") "full arrays: ",nR, &
         !           &  SUM( s(:,nR) ),SUM( b(:,nR) ),SUM( aj(:,nR) )
         !   end if
         !   write(*,"(A,I4,6ES22.14)") "LMloc arrays: ",nR, &
         !        & SUM( s_LMloc(:,nR) ),SUM( b_LMloc(:,nR) ),SUM( aj_LMloc(:,nR) )
         !end if
    
      end do
      !write(*,"(A,2ES20.12)") "init z_LMloc = ",get_global_sum(z_LMloc)
      if ( l_cond_ic ) then
         do nR=1,n_r_ic_max
            call scatter_from_rank0_to_lo(b_ic(1,nR),b_ic_LMloc(llm,nR))
            call scatter_from_rank0_to_lo(aj_ic(1,nR),aj_ic_LMloc(llm,nR))
         end do
      end if
    
      !if (DEBUG_OUTPUT) then
      !   if (rank == 0) then
      !      write(*,"(A,4ES20.12)") "getStartFields: z,dzdtLast full = ", &
      !           & SUM( z ),SUM( dzdtLast )
      !   end if!
    
      !   write(*,"(A,4ES20.12)") "getStartFields: z,dzdtLast = ", &
      !        &    SUM( z_LMloc ),SUM( dzdtLast_lo )
      !end if
    
      allocate( workA_LMloc(llm:ulm,n_r_max) )
      allocate( workB_LMloc(llm:ulm,n_r_max) )
    
      !  print*,"Computing derivatives"
      do nLMB=1+rank*nLMBs_per_rank,MIN((rank+1)*nLMBs_per_rank,nLMBs) 
         ! Blocking of loop over all (l,m)
         lmStart=lmStartB(nLMB)
         lmStop =lmStopB(nLMB)
    
         !if (DEBUG_OUTPUT) then
         !   write(*,"(A,I3,10ES22.15)") "after init: w,z,s,b,aj ",nLMB, &
         !        & SUM(w_LMloc), SUM(z_LMloc), SUM(s_LMloc),SUM(b_LMloc),SUM(aj_LMloc)
         !end if
    
    
         if ( l_conv .or. l_mag_kin ) then
            call get_ddr( w_LMloc,dw_LMloc,ddw_LMloc,ulm-llm+1,lmStart-llm+1, &
                          lmStop-llm+1,n_r_max,n_cheb_max,workA_LMloc,        &
                          workB_LMloc,chebt_oc,drx,ddrx )
            call get_dr( z_LMloc,dz_LMloc,ulm-llm+1, lmStart-llm+1,lmStop-llm+1, &
                         n_r_max,n_cheb_max,workA_LMloc,workB_LMloc,             &
                         chebt_oc,drx )
         end if
    
         if ( l_mag .or. l_mag_kin  ) then
            call get_ddr( b_LMloc,db_LMloc,ddb_LMloc,ulmMag-llmMag+1, &
                          lmStart-llmMag+1,lmStop-llmMag+1,n_r_max,   &
                          n_cheb_max,workA_LMloc,workB_LMloc,         &
                          chebt_oc,drx,ddrx )
            call get_ddr( aj_LMloc,dj_LMloc,ddj_LMloc,ulmMag-llmMag+1, &
                          lmStart-llmMag+1,lmStop-llmMag+1,n_r_max,    &
                          n_cheb_max,workA_LMloc,workB_LMloc,          &
                          chebt_oc,drx,ddrx )
         end if
         if ( l_cond_ic ) then
            call get_ddr_even(b_ic_LMloc,db_ic_LMLoc,ddb_ic_LMloc,       &
                              ulmMag-llmMag+1,lmStart-llmMag+1,          &
                              lmStop-llmMag+1,n_r_ic_max,n_cheb_ic_max,  &
                              dr_fac_ic,workA_LMloc,workB_LMloc,         &
                              chebt_ic, chebt_ic_even)
            call get_ddr_even(aj_ic_LMloc,dj_ic_LMloc,ddj_ic_LMloc,      &
                              ulmMag-llmMag+1,lmStart-llmMag+1,          &
                              lmStop-llmMag+1,n_r_ic_max,n_cheb_ic_max,  &
                              dr_fac_ic,workA_LMloc,workB_LMloc,         &
                              chebt_ic, chebt_ic_even)
         end if
    
         if ( l_LCR ) then
            do nR=n_r_cmb,n_r_icb-1
               if ( nR<=n_r_LCR ) then
                  do lm=lmStart,lmStop
                     l=lo_map%lm2l(lm)
                     m=lo_map%lm2m(lm)
    
                     b_LMloc(lm,nR)=(r(n_r_LCR)/r(nR))**real(l,cp)* &
                                     b_LMloc(lm,n_r_LCR)
                     db_LMloc(lm,nR)=-real(l,cp)*(r(n_r_LCR))**real(l,cp)/ &
                                     (r(nR))**real(l+1,cp)*b_LMloc(lm,n_r_LCR)
                     ddb_LMloc(lm,nR)=real(l,cp)*real(l+1,cp)*    &
                                      (r(n_r_LCR))**(real(l,cp))/ &
                                      (r(nR))**real(l+2,cp)*b_LMloc(lm,n_r_LCR)
                     aj_LMloc(lm,nR)=zero
                     dj_LMloc(lm,nR)=zero
                     ddj_LMloc(lm,nR)=zero
                  end do
               end if
            end do
         end if
    
    
         if ( l_heat ) then
            !-- Get radial derivatives of entropy:
            !if (DEBUG_OUTPUT) then
             !  do nR=1,n_r_max
            !      write(*,"(A,I4)") "nR=",nR
            !      do lm=lmStart,lmStop
            !         write(*,"(4X,A,4I5,2ES22.14)") "s : ", nR,lm, &
            !              &  lo_map%lm2l(lm),lo_map%lm2m(lm),s_LMloc(lm,nR)
            !      end do
            !   end do
            !end if
            call get_dr( s_LMloc,ds_LMloc,ulm-llm+1, lmStart-llm+1,lmStop-llm+1, &
                         n_r_max,n_cheb_max,workA_LMloc,workB_LMloc,             &
                         chebt_oc,drx )
         end if
    
         if ( DEBUG_OUTPUT ) then
            !do nR=1,n_r_max
            !   write(*,"(A,I5,4ES22.14)") "Rdep: s,ds : ", nR,  &
            !        & SUM(s_LMloc(lmStart:lmStop,nR)),SUM(ds_LMloc(lmStart:lmStop,nR))
            !end do
            !do lm=lmStart,lmStop
            !   write(*,"(A,3I5,4ES22.14)") "s,ds : ", lm,lo_map%lm2l(lm), &
            !        & lo_map%lm2m(lm), SUM(s_LMloc(lm,:)),SUM(ds_LMloc(lm,:))
            !end do
            write(*,"(A,I3,10ES22.15)") "derivatives: w,z,s,b,aj ", &
                 & nLMB, SUM(dw_LMloc), SUM(dz_LMloc),              &
                 & SUM(ds_LMloc),SUM(db_LMloc),SUM(dj_LMloc)
         end if
         
      end do
    
      deallocate(workA_LMloc)
      deallocate(workB_LMloc)
      !--- Get symmetry properties of tops excluding l=m=0:
      sES=0.0_cp
      sEA=0.0_cp
      sAA=0.0_cp
      do m=0,l_max,minc
         do l=m,l_max
            if ( l > 0 ) then
               if ( MOD(l+m,2) == 0 ) then
                  sES=sES+cc2real(tops(l,m),m)
               else
                  sEA=sEA+cc2real(tops(l,m),m)
               end if
               if ( m /= 0 ) sAA=sAA+cc2real(tops(l,m),m)
            end if
         end do
      end do
      if ( sEA+sES == 0 ) then
         write(message,'(''! Only l=m=0 comp. in tops:'')')
         call logWrite(message)
      else
         sEA=sqrt(sEA/(sEA+sES))
         sAA=sqrt(sAA/(sEA+sES))
         write(message,'(''! Rel. RMS equ. asym. tops:'',ES16.6)') sEA
         call logWrite(message)
         write(message,'(''! Rel. RMS axi. asym. tops:'',ES16.6)') sAA
         call logWrite(message)
      end if
    
      !----- Get changes in mantle and ic rotation rate:
      if ( .not. l_mag_LF ) then
         lorentz_torque_icLast=0.0_cp
         lorentz_torque_maLast=0.0_cp
      end if
      if ( l_z10mat ) then
         l1m0=lo_map%lm2(1,0)
         coex=-two*(alpha-one)
         if ( ( .not. l_SRMA .and. ktopv == 2 .and. l_rot_ma ).and.&
              & (l1m0 >= llm .and.l1m0 <= ulm) ) then
            d_omega_ma_dt=LFfac*c_lorentz_ma*lorentz_torque_maLast
            d_omega_ma_dtLast=d_omega_ma_dt -           &
                 coex * ( two*or1(1)*real(z_LMloc(l1m0,1)) - &
                 real(dz_LMloc(l1m0,1)) )
         end if
         if ( ( .not. l_SRIC .and. kbotv == 2 .and. l_rot_ic ).and.&
              & (l1m0 >= llm .and. l1m0 <= ulm) ) then
            d_omega_ic_dt=LFfac*c_lorentz_ic*lorentz_torque_icLast
            d_omega_ic_dtLast= d_omega_ic_dt +                      &
                 coex * ( two*or1(n_r_max)*real(z_LMloc(l1m0,n_r_max)) - &
                 real(dz_LMloc(l1m0,n_r_max)) )
         end if
      else
         d_omega_ma_dtLast=0.0_cp
         d_omega_ic_dtLast=0.0_cp
      end if
    
    
         ! --------------- end of insertion ----------
    
      !print*,"Start redistribution in getStartfields"
      ! start the redistribution
      if (l_heat) then
         call lo2r_redist_start(lo2r_s,s_LMloc_container,s_Rloc_container)
      end if
      if (l_conv) then
         call lo2r_redist_start(lo2r_z,z_LMloc_container,z_Rloc_container)
         call lo2r_redist_start(lo2r_w,w_LMloc_container,w_Rloc_container)
         call lo2r_redist_start(lo2r_p,p_LMloc_container,p_Rloc_container)
      end if
    
      if (l_mag) then
         call lo2r_redist_start(lo2r_b,b_LMloc_container,b_Rloc_container)
         call lo2r_redist_start(lo2r_aj,aj_LMloc_container,aj_Rloc_container)
      end if
    
      !print*,"End of getStartFields"
      !PERFOFF
   end subroutine getStartFields
end module start_fields
