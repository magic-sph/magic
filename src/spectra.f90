module spectra

   use parallel_mod
   use precision_mod
   use truncation, only: n_r_max, n_r_ic_maxMag, n_r_maxMag, &
                         n_r_ic_max, l_max, minc
   use radial_data, only: n_r_cmb, n_r_icb
   use radial_functions, only: orho1, orho2, r_ic, i_costf1_ic_init, &
                               d_costf1_ic_init, r, i_costf_init,    &
                               d_costf_init, or2, r_icb, dr_fac_ic,  &
                               drx, dr_fac
   use physical_parameters, only: LFfac
   use num_param, only: eScale, tScale
   use blocking, only: lo_map, st_map
   use horizontal_data, only: dLh
   use logic, only: l_mag, l_anel, l_cond_ic, l_heat
   use output_data, only: tag, log_file, nLF, n_mag_spec_file, &
                          n_u2_spec_file, n_kin_spec_file
   use LMLoop_data,only: llm,ulm,llmMag,ulmMag
   use useful, only: cc2real, cc22real, safeOpen, safeClose
   use integration, only: rInt_R, rIntIC, rInt
   use constants, only: pi, vol_oc, half, one, four

   implicit none
  
   private
 
   real(cp), allocatable :: e_p_l_ave(:),e_p_m_ave(:)
   real(cp), allocatable :: e_p2_l_ave(:),e_p2_m_ave(:)
   real(cp), allocatable :: e_t_l_ave(:),e_t_m_ave(:)
   real(cp), allocatable :: e_t2_l_ave(:),e_t2_m_ave(:)
   real(cp), allocatable :: e_cmb_l_ave(:),e_cmb_m_ave(:)
   real(cp), allocatable :: e_cmb2_l_ave(:),e_cmb2_m_ave(:)
 
   real(cp), allocatable :: ek_p_l_ave(:),ek_p_m_ave(:)
   real(cp), allocatable :: ek_p2_l_ave(:),ek_p2_m_ave(:)
   real(cp), allocatable :: ek_t_l_ave(:),ek_t_m_ave(:)
   real(cp), allocatable :: ek_t2_l_ave(:),ek_t2_m_ave(:)

   real(cp), allocatable :: T_ave(:)
   real(cp), allocatable :: T_ICB_ave(:)
   real(cp), allocatable :: dT_ICB_ave(:)
   real(cp), allocatable :: T2_ave(:)
   real(cp), allocatable :: T_ICB2_ave(:)
   real(cp), allocatable :: dT_ICB2_ave(:)
 
   public :: initialize_spectra, spectrum, spectrum_average, &
             spectrum_temp, spectrum_temp_average

contains

   subroutine initialize_spectra

      allocate( e_p_l_ave(0:l_max),e_p_m_ave(0:l_max) )
      allocate( e_p2_l_ave(0:l_max),e_p2_m_ave(0:l_max) )
      allocate( e_t_l_ave(0:l_max),e_t_m_ave(0:l_max) )
      allocate( e_t2_l_ave(0:l_max),e_t2_m_ave(0:l_max) )
      allocate( e_cmb_l_ave(0:l_max),e_cmb_m_ave(0:l_max) )
      allocate( e_cmb2_l_ave(0:l_max),e_cmb2_m_ave(0:l_max) )

      allocate( ek_p_l_ave(0:l_max),ek_p_m_ave(0:l_max) )
      allocate( ek_p2_l_ave(0:l_max),ek_p2_m_ave(0:l_max) )
      allocate( ek_t_l_ave(0:l_max),ek_t_m_ave(0:l_max) )
      allocate( ek_t2_l_ave(0:l_max),ek_t2_m_ave(0:l_max) )

      if ( l_heat ) then
         allocate( T_ave(l_max+1) )
         allocate( T_ICB_ave(l_max+1) )
         allocate( dT_ICB_ave(l_max+1) )
         allocate( T2_ave(l_max+1) )
         allocate( T_ICB2_ave(l_max+1) )
         allocate( dT_ICB2_ave(l_max+1) )
      end if

   end subroutine initialize_spectra
!----------------------------------------------------------------------------
   subroutine spectrum_average(n_time_ave,l_stop_time,             &
       &                      time_passed,time_norm,b,aj,db,BV)

      !-- Input variables:
      integer,          intent(in) :: n_time_ave
      logical,          intent(in) :: l_stop_time
      real(cp),         intent(in) :: time_passed
      real(cp),         intent(in) :: time_norm
      complex(cp),      intent(in) :: b(llm:ulm,n_r_max)
      complex(cp),      intent(in) :: aj(llm:ulm,n_r_max)
      complex(cp),      intent(in) :: db(llm:ulm,n_r_max)
      character(len=1), intent(in) :: BV

      !-- Output variables: 
      real(cp) :: e_p_l(0:l_max),e_t_l(0:l_max)
      real(cp) :: e_cmb_l(0:l_max)
      real(cp) :: e_p_m(0:l_max),e_t_m(0:l_max)
      real(cp) :: e_cmb_m(0:l_max)

      !-- Local variables:
      character(len=85) :: outFile
      integer :: nOut
      integer :: nR,lm,l,m,ierr

      real(cp) :: e_p_temp,e_t_temp
      real(cp) :: fac
      real(cp) :: dt_norm
      real(cp) :: SDp_l,SDp_m,SDt_l,SDt_m,SDcmb_l,SDcmb_m

      real(cp) :: e_p_r_l(n_r_max,0:l_max),e_p_r_l_global(n_r_max,0:l_max)
      real(cp) :: e_t_r_l(n_r_max,0:l_max),e_t_r_l_global(n_r_max,0:l_max)
      real(cp) :: e_p_r_m(n_r_max,0:l_max),e_p_r_m_global(n_r_max,0:l_max)
      real(cp) :: e_t_r_m(n_r_max,0:l_max),e_t_r_m_global(n_r_max,0:l_max)

      if ( BV == 'V' ) then ! kinetic spectrum (correction of density)

         do nR=1,n_r_max
            do l=0,l_max
               e_p_r_l(nR,l)=0.0_cp
               e_t_r_l(nR,l)=0.0_cp
               e_p_r_m(nR,l)=0.0_cp
               e_t_r_m(nR,l)=0.0_cp
            end do
            !do lm=2,lm_max
            do lm=max(llm,2),ulm
               l =lo_map%lm2l(lm)
               m =lo_map%lm2m(lm)
               e_p_temp= orho1(nR) * dLh(st_map%lm2(l,m)) * (               &
                    &      dLh(st_map%lm2(l,m))*or2(nR)*cc2real(b(lm,nR),m) &
                    &      + cc2real(db(lm,nR),m) )
               e_t_temp=orho1(nR)*dLh(st_map%lm2(l,m))*cc2real(aj(lm,nR),m)
               e_p_r_l(nR,l)=e_p_r_l(nR,l)+e_p_temp
               e_t_r_l(nR,l)=e_t_r_l(nR,l)+e_t_temp
               e_p_r_m(nR,m)=e_p_r_m(nR,m)+e_p_temp
               e_t_r_m(nR,m)=e_t_r_m(nR,m)+e_t_temp
            end do    ! do loop over lms in block 
         end do    ! radial grid points

      else ! magnetic spectrum

         do nR=1,n_r_max
            do l=0,l_max
               e_p_r_l(nR,l)=0.0_cp
               e_t_r_l(nR,l)=0.0_cp
               e_p_r_m(nR,l)=0.0_cp
               e_t_r_m(nR,l)=0.0_cp
            end do
            do lm=max(2,llm),ulm
               l =lo_map%lm2l(lm)
               m =lo_map%lm2m(lm)
               e_p_temp=  dLh(st_map%lm2(l,m)) * (                           &
                    &       dLh(st_map%lm2(l,m))*or2(nR)*cc2real(b(lm,nR),m) &
                    &       + cc2real(db(lm,nR),m) )
               e_t_temp=dLh(st_map%lm2(l,m))*cc2real(aj(lm,nR),m)
               e_p_r_l(nR,l)=e_p_r_l(nR,l)+e_p_temp
               e_t_r_l(nR,l)=e_t_r_l(nR,l)+e_t_temp
               e_p_r_m(nR,m)=e_p_r_m(nR,m)+e_p_temp
               e_t_r_m(nR,m)=e_t_r_m(nR,m)+e_t_temp
            end do    ! do loop over lms in block 
         end do    ! radial grid points

      end if

#ifdef WITH_MPI
      call MPI_Reduce(e_p_r_l,e_p_r_l_global,n_r_max*(l_max+1),&
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_t_r_l,e_t_r_l_global,n_r_max*(l_max+1),&
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_p_r_m,e_p_r_m_global,n_r_max*(l_max+1),&
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_t_r_m,e_t_r_m_global,n_r_max*(l_max+1),&
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
      e_p_r_l_global=e_p_r_l
      e_t_r_l_global=e_t_r_l
      e_p_r_m_global=e_p_r_m
      e_t_r_m_global=e_t_r_m
#endif

      if ( rank == 0 ) then
         !-- Radial Integrals:
         fac=half*eScale
         if ( BV == 'B' ) fac=fac*LFfac
         do l=0,l_max
            e_p_l(l)  =fac*rInt_R(e_p_r_l_global(1,l),n_r_max,n_r_max,drx,   &
                 &                           i_costf_init,d_costf_init)
            e_t_l(l)  =fac*rInt_R(e_t_r_l_global(1,l),n_r_max,n_r_max,drx,   &
                 &                           i_costf_init,d_costf_init)
            e_p_m(l)  =fac*rInt_R(e_p_r_m_global(1,l),n_r_max,n_r_max,drx,   &
                 &                           i_costf_init,d_costf_init)
            e_t_m(l)  =fac*rInt_R(e_t_r_m_global(1,l),n_r_max,n_r_max,drx,   &
                 &                           i_costf_init,d_costf_init)
            if ( BV == 'B' ) then 
               e_cmb_l(l)=fac*e_p_r_l_global(1,l)
               e_cmb_m(l)=fac*e_p_r_m_global(1,l)
            end if
         end do


         !-- Averaging:
         if ( n_time_ave == 1 ) then
            do l=0,l_max
               if ( BV == 'B' ) then
                  e_p_l_ave(l)   =time_passed*e_p_l(l)
                  e_t_l_ave(l)   =time_passed*e_t_l(l)
                  e_p2_l_ave(l)  =time_passed*e_p_l(l)*e_p_l(l)
                  e_t2_l_ave(l)  =time_passed*e_t_l(l)*e_t_l(l)
                  e_p_m_ave(l)   =time_passed*e_p_m(l)
                  e_t_m_ave(l)   =time_passed*e_t_m(l)
                  e_p2_m_ave(l)  =time_passed*e_p_m(l)*e_p_m(l)
                  e_t2_m_ave(l)  =time_passed*e_t_m(l)*e_t_m(l)
                  e_cmb_l_ave(l) =time_passed*e_cmb_l(l)
                  e_cmb2_l_ave(l)=time_passed*e_cmb_l(l)*e_cmb_l(l)
                  e_cmb_m_ave(l) =time_passed*e_cmb_m(l)
                  e_cmb2_m_ave(l)=time_passed*e_cmb_m(l)*e_cmb_m(l)
               else
                  ek_p_l_ave(l)   =time_passed*e_p_l(l)
                  ek_t_l_ave(l)   =time_passed*e_t_l(l)
                  ek_p2_l_ave(l)  =time_passed*e_p_l(l)*e_p_l(l)
                  ek_t2_l_ave(l)  =time_passed*e_t_l(l)*e_t_l(l)
                  ek_p_m_ave(l)   =time_passed*e_p_m(l)
                  ek_t_m_ave(l)   =time_passed*e_t_m(l)
                  ek_p2_m_ave(l)  =time_passed*e_p_m(l)*e_p_m(l)
                  ek_t2_m_ave(l)  =time_passed*e_t_m(l)*e_t_m(l)
               end if
            end do
         else
            do l=0,l_max
               if ( BV == 'B' ) then
                  e_p_l_ave(l)   =e_p_l_ave(l)   +time_passed*e_p_l(l)
                  e_t_l_ave(l)   =e_t_l_ave(l)   +time_passed*e_t_l(l)
                  e_p2_l_ave(l)  =e_p2_l_ave(l)  +time_passed*e_p_l(l)*e_p_l(l)
                  e_t2_l_ave(l)  =e_t2_l_ave(l)  +time_passed*e_t_l(l)*e_t_l(l)
                  e_p_m_ave(l)   =e_p_m_ave(l)   +time_passed*e_p_m(l)
                  e_t_m_ave(l)   =e_t_m_ave(l)   +time_passed*e_t_m(l)
                  e_p2_m_ave(l)  =e_p2_m_ave(l)  +time_passed*e_p_m(l)*e_p_m(l)
                  e_t2_m_ave(l)  =e_t2_m_ave(l)  +time_passed*e_t_m(l)*e_t_m(l)
                  e_cmb_l_ave(l) =e_cmb_l_ave(l) +time_passed*e_cmb_l(l)
                  e_cmb2_l_ave(l)=e_cmb2_l_ave(l)+time_passed*e_cmb_l(l)*e_cmb_l(l)
                  e_cmb_m_ave(l) =e_cmb_m_ave(l) +time_passed*e_cmb_m(l)
                  e_cmb2_m_ave(l)=e_cmb2_m_ave(l)+time_passed*e_cmb_m(l)*e_cmb_m(l)
               else
                  ek_p_l_ave(l)   =ek_p_l_ave(l) +time_passed*e_p_l(l)
                  ek_t_l_ave(l)   =ek_t_l_ave(l) +time_passed*e_t_l(l)
                  ek_p2_l_ave(l)  =ek_p2_l_ave(l)+time_passed*e_p_l(l)*e_p_l(l)
                  ek_t2_l_ave(l)  =ek_t2_l_ave(l)+time_passed*e_t_l(l)*e_t_l(l)
                  ek_p_m_ave(l)   =ek_p_m_ave(l) +time_passed*e_p_m(l)
                  ek_t_m_ave(l)   =ek_t_m_ave(l) +time_passed*e_t_m(l)
                  ek_p2_m_ave(l)  =ek_p2_m_ave(l)+time_passed*e_p_m(l)*e_p_m(l)
                  ek_t2_m_ave(l)  =ek_t2_m_ave(l)+time_passed*e_t_m(l)*e_t_m(l)
               end if
            end do
         end if


         !-- Output: every 10th averaging step and at end of run
         if ( l_stop_time .or. mod(n_time_ave,10) == 0 ) then

            !------ Output:
            dt_norm=one/time_norm
            if ( BV == 'B' ) then
               outFile='mag_spec_ave.'//TAG
            else if ( BV == 'V' ) then
               outFile='kin_spec_ave.'//TAG
            else
               write(*,*) 'WRONG BV INPUT TO spectrum_average!'
               stop
            end if
            nOut   =93
            open(nOut, file=outFile, status='unknown')
            if ( BV == 'B' ) then
               do l=0,l_max
                  SDp_l = get_standard_deviation(dt_norm,e_p_l_ave(l),e_p2_l_ave(l))
                  SDp_m = get_standard_deviation(dt_norm,e_p_m_ave(l),e_p2_m_ave(l))
                  SDt_l = get_standard_deviation(dt_norm,e_t_l_ave(l),e_t2_l_ave(l))
                  SDt_m = get_standard_deviation(dt_norm,e_t_m_ave(l),e_t2_m_ave(l))
                  SDcmb_l=get_standard_deviation(dt_norm,e_cmb_l_ave(l),e_cmb2_l_ave(l))
                  SDcmb_m=get_standard_deviation(dt_norm,e_cmb_m_ave(l),e_cmb2_m_ave(l))
                  !write(*,"(A,I4,5ES22.14)") "SDcmb_m = ",l,dt_norm,SDcmb_m, &
                  !     &  dt_norm*e_cmb2_m_ave(l),    &
                  !     & (dt_norm*e_cmb_m_ave(l))**2, &
                  !     & dt_norm*e_cmb2_m_ave(l) - (dt_norm*e_cmb_m_ave(l))**2
                  write(93,'(2X,1P,I4,16ES16.8)') l,                          &
                       &  dt_norm*e_p_l_ave(l),   dt_norm*e_p_m_ave(l),       &
                       &  dt_norm*e_t_l_ave(l),   dt_norm*e_t_m_ave(l),       &
                       &  dt_norm*e_cmb_l_ave(l), dt_norm*e_cmb_m_ave(l),     &
                       &  dt_norm*e_p_l_ave(l)+SDp_l,                         &
                       &  dt_norm*e_p_l_ave(l)-SDp_l,                         &
                       &  dt_norm*e_p_m_ave(l)+SDp_m,                         &
                       &  dt_norm*e_p_m_ave(l)-SDp_m,                         &
                       &  dt_norm*e_t_l_ave(l)+SDt_l,                         &
                       &  dt_norm*e_t_l_ave(l)-SDt_l,                         &
                       &  dt_norm*e_t_m_ave(l)+SDt_m,                         &
                       &  dt_norm*e_t_m_ave(l)-SDt_m,                         &
                       &  dt_norm*e_cmb_m_ave(l)+SDcmb_m,                     &
                       &  dt_norm*e_cmb_m_ave(l)-SDcmb_m 
               end do
            else
               do l=0,l_max
                  write(93,'(2X,1P,I4,8ES16.8)') l,                           &
                       &  dt_norm*ek_p_l_ave(l), dt_norm*ek_p_m_ave(l),       &
                       &  dt_norm*ek_t_l_ave(l), dt_norm*ek_t_m_ave(l),       &
                       &  dt_norm*ek_p2_l_ave(l),dt_norm*ek_p2_m_ave(l),      &
                       &  dt_norm*ek_t2_l_ave(l),dt_norm*ek_t2_m_ave(l)
               end do
            end if
            close(nOut)

            if ( l_stop_time ) then
               call safeOpen(nLF,log_file)
               write(nLF,"(/,A,A)") ' ! TIME AVERAGED SPECTRA STORED IN FILE: ', &
                     outFile
               write(nLF,"(A,I5)")  ' !              No. of averaged spectra: ', &
                     n_time_ave
               call safeClose(nLF)
            end if

         end if
      end if

   end subroutine spectrum_average
!----------------------------------------------------------------------------
   real(cp) function get_standard_deviation(dt_norm,mean,sum_of_squares) result(stdev)

      real(cp), intent(in) :: dt_norm,mean,sum_of_squares

      real(cp) :: mean2,variance
      real(cp), parameter :: tolerance = 5.0_cp*epsilon(one)

      mean2    = (dt_norm*mean)**2
      variance = abs(dt_norm*sum_of_squares - mean2 )
      if ( mean2 /= 0.0_cp ) then
         if ( variance/mean2 < tolerance ) then
            stdev = 0.0_cp
         else
            stdev = sqrt(variance)
         end if
      else
         stdev = 0.0_cp
      end if

   end function get_standard_deviation
!----------------------------------------------------------------------------
   subroutine spectrum(time,n_spec,w,dw,z,b,db,aj,b_ic,db_ic,aj_ic)
      !
      !  calculates magnetic energy  = 1/2 Integral(B^2 dV)
      !  integration in theta,phi by summation over harmonic coeffs.
      !  integration in r by Chebycheff integrals
      !
      !  Output:
      !  enbp: Total poloidal        enbt: Total toroidal
      !  apome: Axisym. poloidal     atome: Axisym. toroidal
      !
    
      !-- Input of variables:
      integer,     intent(in) :: n_spec     ! number of spectrum/call, file
      real(cp),    intent(in) :: time
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dw(llm:ulm,n_r_max)
      complex(cp), intent(in) :: z(llm:ulm,n_r_max)
      complex(cp), intent(in) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: db(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: aj(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: db_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
    
      !-- Output:
      real(cp) :: b_rms
      real(cp) :: e_mag_p_l(l_max),e_mag_t_l(l_max)
      real(cp) :: e_kin_p_l(l_max),e_kin_t_l(l_max)
      real(cp) :: e_mag_p_ic_l(l_max),e_mag_t_ic_l(l_max)
      real(cp) :: u2_p_l(l_max),u2_t_l(l_max)
    
      real(cp) :: e_mag_p_m(l_max+1),e_mag_t_m(l_max+1)
      real(cp) :: e_kin_p_m(l_max+1),e_kin_t_m(l_max+1)
      real(cp) :: e_mag_p_ic_m(l_max+1),e_mag_t_ic_m(l_max+1)
      real(cp) :: u2_p_m(l_max+1),u2_t_m(l_max+1)
    
      real(cp) :: e_mag_cmb_l(l_max)
      real(cp) :: e_mag_cmb_m(l_max+1)
      real(cp) :: e_kin_nearSurf_l(l_max)
      real(cp) :: e_kin_nearSurf_m(l_max+1)
    
      real(cp) :: eCMB(l_max),eCMB_global(l_max)
    
      !-- local:
      character(len=14) :: string
      character(len=72) :: mag_spec_file,kin_spec_file,u2_spec_file
      integer :: n_r,lm,ml,l,mc,m,n_const
    
      real(cp) :: r_ratio,O_r_icb_E_2
      real(cp) :: e_mag_p_temp,e_mag_t_temp
      real(cp) :: e_kin_p_temp,e_kin_t_temp
      real(cp) :: u2_p_temp,u2_t_temp
      real(cp) :: O_surface
      real(cp) :: fac_mag,fac_kin
      real(cp) :: nearSurfR
    
      real(cp) :: e_mag_p_r_l(n_r_max,l_max),e_mag_p_r_l_global(n_r_max,l_max)
      real(cp) :: e_mag_t_r_l(n_r_max,l_max),e_mag_t_r_l_global(n_r_max,l_max)
      real(cp) :: e_kin_p_r_l(n_r_max,l_max),e_kin_p_r_l_global(n_r_max,l_max)
      real(cp) :: e_kin_t_r_l(n_r_max,l_max),e_kin_t_r_l_global(n_r_max,l_max)
      real(cp) :: u2_p_r_l(n_r_max,l_max),u2_p_r_l_global(n_r_max,l_max)
      real(cp) :: u2_t_r_l(n_r_max,l_max),u2_t_r_l_global(n_r_max,l_max)
      real(cp) :: e_mag_p_r_m(n_r_max,l_max+1),e_mag_p_r_m_global(n_r_max,l_max+1)
      real(cp) :: e_mag_t_r_m(n_r_max,l_max+1),e_mag_t_r_m_global(n_r_max,l_max+1)
      real(cp) :: e_kin_p_r_m(n_r_max,l_max+1),e_kin_p_r_m_global(n_r_max,l_max+1)
      real(cp) :: e_kin_t_r_m(n_r_max,l_max+1),e_kin_t_r_m_global(n_r_max,l_max+1)
      real(cp) :: u2_p_r_m(n_r_max,l_max+1),u2_p_r_m_global(n_r_max,l_max+1)
      real(cp) :: u2_t_r_m(n_r_max,l_max+1),u2_t_r_m_global(n_r_max,l_max+1)
    
      real(cp) :: e_mag_p_ic_r_l(n_r_ic_max,l_max)
      real(cp) :: e_mag_p_ic_r_l_global(n_r_ic_max,l_max)
      real(cp) :: e_mag_t_ic_r_l(n_r_ic_max,l_max)
      real(cp) :: e_mag_t_ic_r_l_global(n_r_ic_max,l_max)
      real(cp) :: e_mag_p_ic_r_m(n_r_ic_max,l_max+1)
      real(cp) :: e_mag_p_ic_r_m_global(n_r_ic_max,l_max+1)
      real(cp) :: e_mag_t_ic_r_m(n_r_ic_max,l_max+1)
      real(cp) :: e_mag_t_ic_r_m_global(n_r_ic_max,l_max+1)
    
      complex(cp) :: r_dr_b
    
    
      eCMB=0.0_cp
    
      do n_r=1,n_r_max
    
         do l=1,l_max
            if ( l_mag ) then
               e_mag_p_r_l(n_r,l)=0.0_cp
               e_mag_t_r_l(n_r,l)=0.0_cp
            end if
            if ( l_anel ) then
               u2_p_r_l(n_r,l)=0.0_cp
               u2_t_r_l(n_r,l)=0.0_cp
            end if
            e_kin_p_r_l(n_r,l)=0.0_cp
            e_kin_t_r_l(n_r,l)=0.0_cp
         end do
         do mc=1,l_max+1
            if ( l_mag ) then
               e_mag_p_r_m(n_r,mc)=0.0_cp
               e_mag_t_r_m(n_r,mc)=0.0_cp
            end if
            if ( l_anel ) then
               u2_p_r_m(n_r,mc)=0.0_cp
               u2_t_r_m(n_r,mc)=0.0_cp
            end if
            e_kin_p_r_m(n_r,mc)=0.0_cp
            e_kin_t_r_m(n_r,mc)=0.0_cp
         end do
    
         !do lm=2,lm_max
         do lm=max(llm,2),ulm
    
            l  =lo_map%lm2l(lm)
            m  =lo_map%lm2m(lm)
            mc=m+1
    
            if ( l_mag ) then
               e_mag_p_temp= dLh(st_map%lm2(l,m)) * ( &
                    &          dLh(st_map%lm2(l,m))*or2(n_r)*cc2real(b(lm,n_r),m) + &
                    &          cc2real(db(lm,n_r),m) )
               e_mag_t_temp=dLh(st_map%lm2(l,m))*cc2real(aj(lm,n_r),m)
            end if
            if ( l_anel ) then
               u2_p_temp=  orho2(n_r)*dLh(st_map%lm2(l,m)) *  ( &
                    &        dLh(st_map%lm2(l,m))*or2(n_r)*cc2real(w(lm,n_r),m) + &
                    &        cc2real(dw(lm,n_r),m) )
               u2_t_temp=orho2(n_r)*dLh(st_map%lm2(l,m))*cc2real(z(lm,n_r),m)
            end if
            e_kin_p_temp= orho1(n_r)*dLh(st_map%lm2(l,m)) *  ( &
                 &          dLh(st_map%lm2(l,m))*or2(n_r)*cc2real(w(lm,n_r),m) + &
                 &          cc2real(dw(lm,n_r),m) )
            e_kin_t_temp=orho1(n_r)*dLh(st_map%lm2(l,m))*cc2real(z(lm,n_r),m)
    
            !----- l-spectra:
            if ( l_mag ) then
               e_mag_p_r_l(n_r,l) = e_mag_p_r_l(n_r,l) + e_mag_p_temp
               e_mag_t_r_l(n_r,l) = e_mag_t_r_l(n_r,l) + e_mag_t_temp
               if ( m == 0 .and. n_r == n_r_cmb ) eCMB(l)=e_mag_p_temp
            end if
            if ( l_anel ) then
               u2_p_r_l(n_r,l) = u2_p_r_l(n_r,l) + u2_p_temp
               u2_t_r_l(n_r,l) = u2_t_r_l(n_r,l) + u2_t_temp
            end if
            e_kin_p_r_l(n_r,l) = e_kin_p_r_l(n_r,l) + e_kin_p_temp
            e_kin_t_r_l(n_r,l) = e_kin_t_r_l(n_r,l) + e_kin_t_temp
    
            !----- m-spectra:
            if ( l_mag ) then
               e_mag_p_r_m(n_r,mc) = e_mag_p_r_m(n_r,mc) + e_mag_p_temp
               e_mag_t_r_m(n_r,mc) = e_mag_t_r_m(n_r,mc) + e_mag_t_temp
            end if
            if ( l_anel ) then
               u2_p_r_m(n_r,mc) = u2_p_r_m(n_r,mc) + u2_p_temp
               u2_t_r_m(n_r,mc) = u2_t_r_m(n_r,mc) + u2_t_temp
            end if
            e_kin_p_r_m(n_r,mc)=e_kin_p_r_m(n_r,mc) + e_kin_p_temp
            e_kin_t_r_m(n_r,mc)=e_kin_t_r_m(n_r,mc) + e_kin_t_temp
    
         end do    ! do loop over lms in block
    
      end do    ! radial grid points
    
      ! ----------- We need a reduction here ----------------
      ! first the l-spectra
#ifdef WITH_MPI
      if ( l_mag ) then
         call MPI_Reduce(e_mag_p_r_l, e_mag_p_r_l_global, n_r_max*l_max,&
              &MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_Reduce(e_mag_t_r_l, e_mag_t_r_l_global, n_r_max*l_max,&
              &MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      end if
      if ( l_anel ) then
         call MPI_Reduce(u2_p_r_l, u2_p_r_l_global, n_r_max*l_max,&
              &MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_Reduce(u2_t_r_l, u2_t_r_l_global, n_r_max*l_max,&
              &MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      end if
      call MPI_Reduce(e_kin_p_r_l, e_kin_p_r_l_global, n_r_max*l_max,&
           &MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_kin_t_r_l, e_kin_t_r_l_global, n_r_max*l_max,&
           &MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    
      ! then the m-spectra
      if ( l_mag ) then
         call MPI_Reduce(e_mag_p_r_m, e_mag_p_r_m_global, n_r_max*(l_max+1),&
              &MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_Reduce(e_mag_t_r_m, e_mag_t_r_m_global, n_r_max*(l_max+1),&
              &MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_Reduce(eCMB, eCMB_global, l_max,&
              &MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      end if
      if ( l_anel ) then
         call MPI_Reduce(u2_p_r_m, u2_p_r_m_global, n_r_max*(l_max+1),&
              &MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_Reduce(u2_t_r_m, u2_t_r_m_global, n_r_max*(l_max+1),&
              &MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      end if
      call MPI_Reduce(e_kin_p_r_m, e_kin_p_r_m_global, n_r_max*(l_max+1),&
           &MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_kin_t_r_m, e_kin_t_r_m_global, n_r_max*(l_max+1),&
           &MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
      if ( l_mag ) then
         e_mag_p_r_l_global=e_mag_p_r_l
         e_mag_t_r_l_global=e_mag_t_r_l
         e_mag_p_r_m_global=e_mag_p_r_m
         e_mag_t_r_m_global=e_mag_t_r_m
         eCMB_global       =eCMB
      end if
      if ( l_anel ) then
         u2_p_r_l_global=u2_p_r_l
         u2_t_r_l_global=u2_t_r_l
         u2_p_r_m_global=u2_p_r_m
         u2_t_r_m_global=u2_t_r_m
      end if
      e_kin_p_r_l_global=e_kin_p_r_l
      e_kin_t_r_l_global=e_kin_t_r_l
      e_kin_p_r_m_global=e_kin_p_r_m
      e_kin_t_r_m_global=e_kin_t_r_m
#endif
    
      ! now switch to rank 0 for the postprocess
      
    
      ! Getting appropriate radius index for e_kin_nearSurf spectra
      nearSurfR = r_icb+0.99
      do n_r=2,n_r_max
         if ( r(n_r-1) > nearSurfR .and. r(n_r)  <= nearSurfR ) then
            if ( r(n_r-1)-nearSurfR < nearSurfR-r(n_r) ) then
               n_const=n_r-1
            else
               n_const=n_r
            end if
         end if
      end do
    
      if ( rank == 0 ) then
         !-- Save CMB energy spectra:
         O_surface=one/(four*pi*r(1)*r(1))
    
         if ( l_mag ) then
            b_rms=0.0_cp
            do l=1,l_max
               e_mag_cmb_l(l)=e_mag_p_r_l_global(1,l)
               b_rms=b_rms + e_mag_cmb_l(l)
            end do
            b_rms=sqrt(b_rms*O_surface)
            do mc=1,l_max+1
               e_mag_cmb_m(mc)=e_mag_p_r_m_global(1,mc)
            end do
         end if
    
         !-- Save nearSurf kin energy spectra:
         do l=1,l_max
            e_kin_nearSurf_l(l)=e_kin_p_r_l_global(n_const,l)
         end do
         do mc=1,l_max+1
            e_kin_nearSurf_m(mc)=e_kin_p_r_m_global(n_const,mc)
         end do
    
         !-- Radial Integrals:
         fac_mag=0.5*LFfac*eScale
         fac_kin=0.5*eScale
         do l=1,l_max
            if ( l_mag ) then
               e_mag_p_l(l)=fac_mag*rInt_R(e_mag_p_r_l_global(1,l), &
                            n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
               e_mag_t_l(l)=fac_mag*rInt_R(e_mag_t_r_l_global(1,l), &
                            n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
               e_mag_cmb_l(l)=fac_mag*e_mag_cmb_l(l)
            end if
            if ( l_anel ) then
               u2_p_l(l)  =fac_kin*rInt_R(u2_p_r_l_global(1,l),     &
                           n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
               u2_t_l(l)  =fac_kin*rInt_R(u2_t_r_l_global(1,l),     &
                           n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
            end if
            e_kin_p_l(l)  =fac_kin*rInt_R(e_kin_p_r_l_global(1,l),  &
                           n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
            e_kin_t_l(l)  =fac_kin*rInt_R(e_kin_t_r_l_global(1,l),  &
                           n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
            e_kin_nearSurf_l(l)=fac_kin*e_kin_nearSurf_l(l)
         end do
         do m=1,l_max+1 ! Note: counter m is actual order+1
            if ( l_mag )  then
               e_mag_p_m(m)=fac_mag*rInt_R(e_mag_p_r_m_global(1,m), &
                            n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
               e_mag_t_m(m)=fac_mag*rInt_R(e_mag_t_r_m_global(1,m), &
                            n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
               e_mag_cmb_m(m)=fac_mag*e_mag_cmb_m(m)
            end if
            if ( l_anel ) then
               u2_p_m(m)   =fac_kin*rInt_R(u2_p_r_m_global(1,m),    &
                            n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
               u2_t_m(m)   =fac_kin*rInt_R(u2_t_r_m_global(1,m),    &
                            n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
            end if
            e_kin_p_m(m)   =fac_kin*rInt_R(e_kin_p_r_m_global(1,m), &
                            n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
            e_kin_t_m(m)   =fac_kin*rInt_R(e_kin_t_r_m_global(1,m), &
                            n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
            e_kin_nearSurf_m(m)=fac_kin*e_kin_nearSurf_m(m)
         end do
      end if
    
      !-- inner core:
    
      if ( l_cond_ic ) then
    
         O_r_icb_E_2=one/(r_ic(1)*r_ic(1))
         do n_r=1,n_r_ic_max
            r_ratio=r_ic(n_r)/r_ic(1)
            do mc=1,l_max+1
               e_mag_p_ic_r_m(n_r,mc)=0.0_cp
               e_mag_t_ic_r_m(n_r,mc)=0.0_cp
            end do
            do l=1,l_max
               e_mag_p_ic_r_l(n_r,l)=0.0_cp
               e_mag_t_ic_r_l(n_r,l)=0.0_cp
            end do
            !do lm=2,lm_max
            do lm=max(llm,2),ulm
               l =lo_map%lm2l(lm)
               m =lo_map%lm2m(lm)
               mc=m+1
               r_dr_b=r_ic(n_r)*db_ic(lm,n_r)
    
               e_mag_p_temp=                                        &
                    dLh(st_map%lm2(l,m))*O_r_icb_E_2*r_ratio**(2*l) * ( &
                    real((2*l+1)*(l+1),cp)*cc2real(b_ic(lm,n_r),m)   + &
                    real(2*(l+1),cp)*cc22real(b_ic(lm,n_r),r_dr_b,m) + &
                    cc2real(r_dr_b,m) )
               e_mag_t_temp= dLh(st_map%lm2(l,m))*r_ratio**(2*l+2) * &
                    cc2real(aj_ic(lm,n_r),m)
    
               e_mag_p_ic_r_l(n_r,l)=e_mag_p_ic_r_l(n_r,l) + &
                    e_mag_p_temp
               e_mag_t_ic_r_l(n_r,l)=e_mag_t_ic_r_l(n_r,l) + &
                    e_mag_t_temp
               e_mag_p_ic_r_m(n_r,mc)=e_mag_p_ic_r_m(n_r,mc) + &
                    e_mag_p_temp
               e_mag_t_ic_r_m(n_r,mc)=e_mag_t_ic_r_m(n_r,mc) + &
                    e_mag_t_temp
            end do  ! loop over lm's
         end do ! loop over radial levels
    
#ifdef WITH_MPI
         call MPI_Reduce(e_mag_p_ic_r_l, e_mag_p_ic_r_l_global, n_r_ic_max*l_max,&
              &MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_Reduce(e_mag_t_ic_r_l, e_mag_t_ic_r_l_global, n_r_ic_max*l_max,&
              &MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_Reduce(e_mag_p_ic_r_m, e_mag_p_ic_r_m_global, n_r_ic_max*(l_max+1),&
              &MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_Reduce(e_mag_t_ic_r_m, e_mag_t_ic_r_m_global, n_r_ic_max*(l_max+1),&
              &MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
         e_mag_p_ic_r_l_global=e_mag_p_ic_r_l
         e_mag_t_ic_r_l_global=e_mag_t_ic_r_l
         e_mag_p_ic_r_m_global=e_mag_p_ic_r_m
         e_mag_t_ic_r_m_global=e_mag_t_ic_r_m
#endif
    
    
         if ( rank == 0 ) then
            !----- Radial Integrals:
            fac_mag=LFfac*half*eScale
            do l=1,l_max
               e_mag_p_ic_l(l)=fac_mag*rIntIC(e_mag_p_ic_r_l_global(1,l), &
                    n_r_ic_max,dr_fac_ic,i_costf1_ic_init,d_costf1_ic_init)
               e_mag_t_ic_l(l)=fac_mag*rIntIC(e_mag_t_ic_r_l_global(1,l), &
                    n_r_ic_max,dr_fac_ic,i_costf1_ic_init,d_costf1_ic_init)
            end do
            do m=1,l_max+1
               e_mag_p_ic_m(m)=fac_mag*rIntIC(e_mag_p_ic_r_m_global(1,m), &
                    n_r_ic_max,dr_fac_ic,i_costf1_ic_init,d_costf1_ic_init)
               e_mag_t_ic_m(m)=fac_mag*rIntIC(e_mag_t_ic_r_m_global(1,m), &
                    n_r_ic_max,dr_fac_ic,i_costf1_ic_init,d_costf1_ic_init)
            end do
         end if
      else
         do l=1,l_max
            e_mag_p_ic_l(l)=0.0_cp
            e_mag_t_ic_l(l)=0.0_cp
         end do
         do mc=1,l_max+1
            e_mag_p_ic_m(mc)=0.0_cp
            e_mag_t_ic_m(mc)=0.0_cp
         end do
      end if  ! conducting inner core ?
    
      if ( rank == 0 ) then
         !-- Output into files:
         if ( l_mag ) then
            write(string, *) n_spec
            mag_spec_file='mag_spec_'//trim(adjustl(string))//'.'//tag
            open(n_mag_spec_file, file=mag_spec_file, status='unknown')
            if ( n_spec == 0 ) then
               write(n_mag_spec_file,'(1x, &
                    &      ''Magnetic energy spectra of time averaged field:'')')
            else
               write(n_mag_spec_file,'(1x, &
                    &      ''Magnetic energy spectra at time:'', &
                    &      ES20.12)') time*tScale
            end if
            write(n_mag_spec_file,'(1p,i4,11ES16.8)')          &
                 0,0.0_cp,e_mag_p_m(1)   ,0.0_cp,e_mag_t_m(1), &
                 0.0_cp,e_mag_p_ic_m(1),0.0_cp,e_mag_t_ic_m(1),&
                 0.0_cp,e_mag_cmb_m(1),0.0_cp
            do ml=1,l_max
               write(n_mag_spec_file,'(1p,i4,11ES16.8)') &
                    ml,e_mag_p_l(ml),   e_mag_p_m(ml+1), &
                    e_mag_t_l(ml),   e_mag_t_m(ml+1),    &
                    e_mag_p_ic_l(ml),e_mag_p_ic_m(ml+1), &
                    e_mag_t_ic_l(ml),e_mag_t_ic_m(ml+1), &
                    e_mag_cmb_l(ml), e_mag_cmb_m(ml+1),  &
                    eCMB_global(ml)
            end do
            close(n_mag_spec_file)
    
            mag_spec_file='2D_mag_spec_'//trim(adjustl(string))//'.'//tag
            open(n_mag_spec_file, file=mag_spec_file, status='unknown', &
                 form='unformatted')
    
            write(n_mag_spec_file) time*tScale,n_r_max,l_max,minc
            write(n_mag_spec_file) r
            write(n_mag_spec_file) fac_mag*e_mag_p_r_l_global
            write(n_mag_spec_file) fac_mag*e_mag_p_r_m_global
            write(n_mag_spec_file) fac_mag*e_mag_t_r_l_global
            write(n_mag_spec_file) fac_mag*e_mag_t_r_m_global
    
            close(n_mag_spec_file)
         end if
    
         if ( l_anel ) then
            write(string, *) n_spec
            u2_spec_file='u2_spec_'//trim(adjustl(string))//'.'//tag
            open(n_u2_spec_file, file=u2_spec_file, status='unknown')
            if ( n_spec == 0 ) then
               write(n_u2_spec_file,'(1x, &
                    &     ''Velocity square spectra of time averaged field:'')')
            else
               write(n_u2_spec_file,'(1x,                &
                    &     ''Velocity square spectra at time:'', &
                    &     ES20.12)') time*tScale
            end if
            write(n_u2_spec_file,'(1p,i4,4ES16.8)') &
                 &   0,0.0_cp,u2_p_m(1),0.0_cp,u2_t_m(1)
            do ml=1,l_max
               write(n_u2_spec_file,'(1p,i4,4ES16.8)') &
                    &       ml,u2_p_l(ml),u2_p_m(ml+1),u2_t_l(ml),u2_t_m(ml+1)
            end do
            close(n_u2_spec_file)
    
            u2_spec_file='2D_u2_spec_'//trim(adjustl(string))//'.'//tag
            open(n_u2_spec_file, file=u2_spec_file, status='unknown', &
                 form='unformatted')
    
            write(n_u2_spec_file) time*tScale,n_r_max,l_max,minc
            write(n_u2_spec_file) r
            write(n_u2_spec_file) fac_kin*u2_p_r_l_global
            write(n_u2_spec_file) fac_kin*u2_p_r_m_global
            write(n_u2_spec_file) fac_kin*u2_t_r_l_global
            write(n_u2_spec_file) fac_kin*u2_t_r_m_global
    
            close(n_u2_spec_file)
    
         end if
    
         write(string, *) n_spec
         kin_spec_file='kin_spec_'//trim(adjustl(string))//'.'//tag
         open(n_kin_spec_file, file=kin_spec_file, status='unknown')
         if ( n_spec == 0 ) then
            write(n_kin_spec_file,'(1x, &
                 &      ''Kinetic energy spectra of time averaged field:'')')
         else
            write(n_kin_spec_file,'(1x,                      &
                 &      ''Kinetic energy spectra at time:'', &
                 &      ES20.12)') time*tScale
         end if
         write(n_kin_spec_file,'(1p,i4,6ES16.8)')        &
              0,0.0_cp,e_kin_p_m(1),0.0_cp,e_kin_t_m(1),    &
              0.0_cp, e_kin_nearSurf_m(1)
         do ml=1,l_max
            write(n_kin_spec_file,'(1p,i4,6ES16.8)')    &
                 ml,e_kin_p_l(ml),e_kin_p_m(ml+1),     &
                 e_kin_t_l(ml),e_kin_t_m(ml+1),        &
                 e_kin_nearSurf_l(ml), e_kin_nearSurf_m(ml+1)
         end do
         close(n_kin_spec_file)
    
         kin_spec_file='2D_kin_spec_'//trim(adjustl(string))//'.'//tag
         open(n_kin_spec_file, file=kin_spec_file, status='unknown', &
              form='unformatted')
    
         write(n_kin_spec_file) time*tScale,n_r_max,l_max,minc
         write(n_kin_spec_file) r
         write(n_kin_spec_file) fac_kin*e_kin_p_r_l_global
         write(n_kin_spec_file) fac_kin*e_kin_p_r_m_global
         write(n_kin_spec_file) fac_kin*e_kin_t_r_l_global
         write(n_kin_spec_file) fac_kin*e_kin_t_r_m_global
    
         close(n_kin_spec_file)
    
      end if
    
   end subroutine spectrum
!----------------------------------------------------------------------------
   subroutine spectrum_temp_average(n_time_ave,l_stop_time,         &
       &                       time_passed,time_norm,s,ds)

      !-- Direct input:
      integer,     intent(in) :: n_time_ave
      logical,     intent(in) :: l_stop_time
      real(cp),    intent(in) :: time_passed
      real(cp),    intent(in) :: time_norm
      complex(cp), intent(in) :: s(llm:ulm,n_r_max)
      complex(cp), intent(in) :: ds(llm:ulm,n_r_max)

      !-- Local:
      character(len=72) :: outFile
      integer :: n_r,lm,l,m,lc
      real(cp) :: T_temp
      real(cp) :: dT_temp
      real(cp) :: surf_ICB
      real(cp) :: fac,facICB
      real(cp) :: dt_norm

      real(cp) :: T_r_l(n_r_max,l_max+1),T_r_l_global(n_r_max,l_max+1)
      real(cp) :: T_l(l_max+1)
      real(cp) :: T_ICB_l(l_max+1), T_ICB_l_global(l_max+1)
      real(cp) :: dT_ICB_l(l_max+1), dT_ICB_l_global(l_max+1)

      real(cp) :: comp(l_max+1)
      !real(cp) :: y,t
      integer :: nOut,ierr

      T_ICB_l =0.0_cp
      dT_ICB_l=0.0_cp

      do n_r=1,n_r_max
         do l=1,l_max+1
            T_r_l(n_r,l)=0.0_cp
            comp(l) = 0.0_cp
         end do
         do lm=llm,ulm
            l =lo_map%lm2l(lm)
            m =lo_map%lm2m(lm)
            lc=l+1

            T_temp=sqrt(cc2real(s(lm,n_r),m))/or2(n_r)

            !local_sum = 0.0_cp
            !c = 0.0_cp          !A running compensation for lost low-order bits.
            !do i=lb,ub
            !   y = arr_local(i) - c    !So far, so good: c is zero.
            !   t = local_sum + y       !Alas, sum is big, y small, so low-order digits of y are lost.
            !   c = (t - local_sum) - y !(t - sum) recovers the high-order part of y; subtracting y recovers -(low part of y)
            !   local_sum = t           !Algebraically, c should always be zero. Beware eagerly optimising compilers!
            !   !Next time around, the lost low part will be added to y in a fresh attempt.
            !end do
#if 0
            y = T_temp - comp(lc)
            t = T_r_l(n_r,lc) + y
            comp(lc) = (t-T_r_l(n_r,lc)) - y
            T_r_l(n_r,lc) = t
#else
            T_r_l(n_r,lc) =T_r_l(n_r,lc) +  T_temp
#endif

            if ( n_r == n_r_icb ) then
               dT_temp=sqrt(cc2real(ds(lm,n_r),m))/or2(n_r)
               T_ICB_l(lc) =  T_ICB_l(lc) +  T_temp
               dT_ICB_l(lc)= dT_ICB_l(lc) + dT_temp
            end if
         end do    ! do loop over lms in block 
      end do    ! radial grid points 

      ! Reduction over all ranks
#ifdef WITH_MPI
      call MPI_Reduce(T_r_l,T_r_l_global,n_r_max*(l_max+1),&
           &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(T_ICB_l,T_ICB_l_global,l_max+1,&
           &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(dT_ICB_l,dT_ICB_l_global,l_max+1,&
           &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
      T_r_l_global   =T_r_l
      T_ICB_l_global =T_ICB_l
      dT_ICB_l_global=dT_ICB_l
#endif

      if ( rank == 0 .and. l_heat ) then
         !-- Radial Integrals:
         surf_ICB=four*pi*r_icb*r_icb
         fac      =one/vol_oc
         facICB   =one/surf_ICB
         do l=1,l_max+1
            T_l(l)=fac*rInt(T_r_l_global(1,l),n_r_max,dr_fac, &
                 &          i_costf_init,d_costf_init)
            T_ICB_l(l)=facICB*T_ICB_l_global(l)
            dT_ICB_l(l)=facICB*dT_ICB_l_global(l)
         end do

         !-- Averaging:
         if ( n_time_ave == 1 ) then
            do l=1,l_max+1
               T_ave(l)     =time_passed*T_l(l)
               T_ICB_ave(l) =time_passed*T_ICB_l(l)
               dT_ICB_ave(l)=time_passed*dT_ICB_l(l)
               T2_ave(l)    =time_passed*T_l(l)*T_l(l)
               T_ICB2_ave(l)=time_passed*T_ICB_l(l)*T_ICB_l(l)
               dT_ICB2_ave(l)=time_passed*dT_ICB_l(l)*dT_ICB_l(l)
            end do
         else
            do l=1,l_max+1
               T_ave(l)     =T_ave(l)       +time_passed*T_l(l)
               T_ICB_ave(l) =T_ICB_ave(l)   +time_passed*T_ICB_l(l)
               dT_ICB_ave(l)=dT_ICB_ave(l)  +time_passed*dT_ICB_l(l)
               T2_ave(l)    =T2_ave(l)      +time_passed*T_l(l)*T_l(l)
               T_ICB2_ave(l)=T_ICB2_ave(l)  +time_passed*T_ICB_l(l)*T_ICB_l(l)
               dT_ICB2_ave(l)=dT_ICB2_ave(l)+time_passed*dT_ICB_l(l)*dT_ICB_l(l)
            end do
         end if

         !-- Output:
         if ( l_stop_time ) then

            !------ Normalize:
            dt_norm=one/time_norm
            do l=1,l_max+1
               T2_ave(l)     =get_standard_deviation(dt_norm,T_ave(l),T2_ave(l))
               T_ICB2_ave(l) =get_standard_deviation(dt_norm,T_ICB_ave(l),T_ICB2_ave(l))
               dT_ICB2_ave(l)=get_standard_deviation(dt_norm,dT_ICB_ave(l),dT_ICB2_ave(l))
               T_ave(l)      =dt_norm*T_ave(l)
               T_ICB_ave(l)  =dt_norm*T_ICB_ave(l)
               dT_ICB_ave(l) =dt_norm*dT_ICB_ave(l)
            end do

            !------ Output:
            outFile='T_spec_ave.'//TAG
            nOut   =93
            open(nOut,file=outFile,status='unknown')
            do l=1,l_max+1
               write(93,'(2X,1P,I4,6ES16.8)') l,                   &
                    &              T_ave(l),T2_ave(l),             &
                    &              T_ICB_ave(l),T_ICB2_ave(l),     &
                    &              dT_ICB_ave(l),dT_ICB2_ave(l) 
            end do
            close(nOut)

            call safeOpen(nLF,log_file)
            write(nLF,"(/,A,A)") ' ! TIME AVERAGED T/C SPECTRA STORED IN FILE: ', outFile
            write(nLF,"(A,I5)")  ' !              No. of averaged spectra: ', n_time_ave
            call safeClose(nLF)

         end if
      end if

   end subroutine spectrum_temp_average
!----------------------------------------------------------------------------
   subroutine spectrum_temp(time,n_spec,s,ds)
      !
      !  calculates spectra of temperature and composition
      !

      !-- Input variables:
      integer,         intent(in) :: n_spec     ! number of spectrum/call, file
      real(cp),    intent(in) :: time
      complex(cp), intent(in) :: s(llm:ulm,n_r_max)
      complex(cp), intent(in) :: ds(llm:ulm,n_r_max)

      !-- Output variables
      real(cp) :: T_l(l_max+1)
      real(cp) :: T_m(l_max+1)
      real(cp) :: T_ICB_l(l_max+1),T_ICB_l_global(l_max+1)
      real(cp) :: T_ICB_m(l_max+1),T_ICB_m_global(l_max+1)
      real(cp) :: dT_ICB_l(l_max+1),dT_ICB_l_global(l_max+1)
      real(cp) :: dT_ICB_m(l_max+1),dT_ICB_m_global(l_max+1)

      !-- Local variables
      character(len=14) :: string
      character(len=72) :: spec_file
      integer :: n_r,lm,ml,l,mc,m,lc
      real(cp) :: T_temp
      real(cp) :: dT_temp
      real(cp) :: surf_ICB
      real(cp) :: fac,facICB

      real(cp) :: T_r_l(n_r_max,l_max+1),T_r_l_global(n_r_max,l_max+1)
      real(cp) :: T_r_m(n_r_max,l_max+1),T_r_m_global(n_r_max,l_max+1)


      do l=1,l_max+1
         T_l(l)=0.0_cp
         T_ICB_l(l)=0.0_cp
         dT_ICB_l(l)=0.0_cp
         T_m(l)=0.0_cp
         T_ICB_m(l)=0.0_cp
         dT_ICB_m(l)=0.0_cp
      end do

      do n_r=1,n_r_max

         do l=1,l_max+1
            T_r_l(n_r,l)=0.0_cp
            T_ICB_l(l)  =0.0_cp
            dT_ICB_l(l) =0.0_cp
            T_r_m(n_r,l)=0.0_cp
            T_ICB_m(l)  =0.0_cp
            dT_ICB_m(l) =0.0_cp
         end do

         do lm=llm,ulm
            l =lo_map%lm2l(lm)
            m =lo_map%lm2m(lm)
            lc=l+1
            mc=m+1

            T_temp=sqrt(cc2real(s(lm,n_r),m))/or2(n_r)
            dT_temp=sqrt(cc2real(ds(lm,n_r),m))/or2(n_r)
            !----- l-spectra:
            T_r_l(n_r,lc) =T_r_l(n_r,lc) +  T_temp
            !----- m-spectra:
            T_r_m(n_r,mc)=T_r_m(n_r,mc) + T_temp

            !----- ICB spectra:
            if ( n_r == n_r_icb ) then
               T_ICB_l(lc) =T_ICB_l(lc) +T_temp
               T_ICB_m(mc)=T_ICB_m(mc)+T_temp
               dT_ICB_l(lc) =dT_ICB_l(lc) +dT_temp
               dT_ICB_m(mc)=dT_ICB_m(mc)+dT_temp
            end if
         end do

      end do

      ! reduction over all ranks
#ifdef WITH_MPI
      call MPI_Reduce(T_r_l,T_r_l_global,n_r_max*(l_max+1),&
           &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(T_r_m,T_r_m_global,n_r_max*(l_max+1),&
           &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(T_ICB_l,T_ICB_l_global,l_max+1,&
           &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(T_ICB_m,T_ICB_m_global,l_max+1,&
           &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(dT_ICB_l,dT_ICB_l_global,l_max+1,&
           &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(dT_ICB_m,dT_ICB_m_global,l_max+1,&
           &          MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
      T_r_l_global   =T_r_l
      T_r_m_global   =T_r_m
      T_ICB_l_global =T_ICB_l
      T_ICB_m_global =T_ICB_m
      dT_ICB_l_global=dT_ICB_l
      dT_ICB_m_global=dT_ICB_m
#endif

      if ( rank == 0 ) then
         !-- Radial Integrals:
         surf_ICB =four*pi*r_icb*r_icb
         fac      =one/vol_oc
         facICB   =one/surf_ICB
         do l=1,l_max+1
            T_l(l)=fac*rInt(T_r_l_global(1,l),n_r_max,dr_fac, &
                            i_costf_init,d_costf_init)
            T_ICB_l(l)=facICB*T_ICB_l_global(l)
            dT_ICB_l(l)=facICB*dT_ICB_l_global(l)
         end do
         do m=1,l_max+1 ! Note: counter m is actual order+1
            T_m(m)=fac*rInt(T_r_m_global(1,m),n_r_max,dr_fac, &
                            i_costf_init,d_costf_init)
            T_ICB_m(m)=facICB*T_ICB_m_global(m)
            dT_ICB_m(m)=facICB*dT_ICB_m_global(m)
         end do

         !-- Output into files:
         write(string, *) n_spec
         spec_file='T_spec_'//trim(adjustl(string))//'.'//tag
         open(98, file=spec_file, status='unknown')
         write(98,'(1x,''TC spectra at time:'', ES20.12)') time*tScale
         do ml=1,l_max+1
            write(98,'(1P,I4,6ES12.4)')   &
                 ml-1,T_l(ml),T_m(ml),    &
                 T_ICB_l(ml),T_ICB_m(ml), &
                 dT_ICB_l(ml),dT_ICB_m(ml)
         end do
         close(98)

      end if

   end subroutine spectrum_temp
!----------------------------------------------------------------------------
end module spectra
