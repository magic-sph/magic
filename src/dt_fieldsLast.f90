module fieldsLast
   !
   ! This module contains all the work arrays of the previous time-steps
   ! needed to time advance the code.
   ! They are needed in the time-stepping scheme.
   !

   use precision_mod
   use truncation, only: n_r_max, n_r_maxMag, n_r_ic_maxMag, nRstart, nRstop,  &
       &                 nRstartMag, nRstopMag, n_lm_loc, n_lmMag_loc,         &
       &                 n_mlo_loc, n_mloMag_loc, fd_order, fd_order_bound
   use logic, only: l_chemical_conv, l_heat, l_mag, l_cond_ic, l_double_curl, &
       &            l_RMS, l_finite_diff
   use constants, only: zero
   use mem_alloc, only: bytes_allocated
   use time_array

   implicit none

   private

   type(type_tarray), public :: dsdt_dist, dwdt_dist, dpdt_dist, dzdt_dist
   type(type_tarray), public :: dxidt_dist
   type(type_tarray), public :: dbdt_dist, djdt_dist, dbdt_ic_dist, djdt_ic_dist
   type(type_tscalar), public :: domega_ma_dt, domega_ic_dt
   type(type_tscalar), public :: lorentz_torque_ic_dt, lorentz_torque_ma_dt

   complex(cp), public, allocatable, target  :: dflowdt_Rdist_container(:,:,:)
   complex(cp), public, allocatable, target  :: dsdt_Rdist_container(:,:,:)
   complex(cp), public, allocatable, target  :: dxidt_Rdist_container(:,:,:)
   complex(cp), public, allocatable, target  :: dbdt_Rdist_container(:,:,:)
   complex(cp), public, pointer :: dwdt_Rdist(:,:),dzdt_Rdist(:,:)
   complex(cp), public, pointer :: dpdt_Rdist(:,:), dsdt_Rdist(:,:), dVSrLM_Rdist(:,:)
   complex(cp), public, pointer :: dxidt_Rdist(:,:), dVXirLM_Rdist(:,:)
   complex(cp), public, pointer :: dVxVhLM_Rdist(:,:)

   complex(cp), public, pointer :: djdt_Rdist(:,:), dVxBhLM_Rdist(:,:)
   complex(cp), public, pointer :: dbdt_Rdist(:,:)

   ! The same arrays, but now the LM local part
   complex(cp), public, allocatable, target  :: dflowdt_LMdist_container(:,:,:,:)
   complex(cp), public, allocatable, target  :: dsdt_LMdist_container(:,:,:,:)
   complex(cp), public, allocatable, target  :: dxidt_LMdist_container(:,:,:,:)
   complex(cp), public, allocatable, target  :: dbdt_LMdist_container(:,:,:,:)
   complex(cp), public, pointer :: dVSrLM_LMdist(:,:,:), dVXirLM_LMdist(:,:,:)
   complex(cp), public, pointer :: dVxVhLM_LMdist(:,:,:), dVxBhLM_LMdist(:,:,:)

   complex(cp), public, allocatable :: dbdt_CMB_LMdist(:)


   public :: initialize_fieldsLast, finalize_fieldsLast

contains

   subroutine initialize_fieldsLast(nold, nexp, nimp)
      !
      ! This routine defines all the arrays needed to time advance MagIC
      !

      integer, intent(in) :: nold ! Number of storage of the old state
      integer, intent(in) :: nexp ! Number of explicit states
      integer, intent(in) :: nimp ! Number of implicit states

      integer :: n_fields

      call domega_ma_dt%initialize(nold, nexp, nimp)
      call domega_ic_dt%initialize(nold, nexp, nimp)

      call lorentz_torque_ic_dt%initialize(nold, nexp, nimp)
      call lorentz_torque_ma_dt%initialize(nold, nexp, nimp)

      if ( l_heat ) call dsdt_dist%initialize(1, n_mlo_loc, n_r_max, nold, nexp, nimp)
      if ( l_chemical_conv ) call dxidt_dist%initialize(1, n_mlo_loc, n_r_max, nold, &
                                  &                nexp, nimp)

      if ( l_mag ) then
         call dbdt_dist%initialize(1, n_mloMag_loc, n_r_maxMag, nold, nexp, nimp)
         call djdt_dist%initialize(1, n_mloMag_loc, n_r_maxMag, nold, nexp, nimp)
      end if

      if ( l_cond_ic ) then
         call dbdt_ic_dist%initialize(1, n_mloMag_loc, n_r_ic_maxMag, nold, &
              &                  nexp, nimp, l_allocate_exp=.true.)
         call djdt_ic_dist%initialize(1, n_mloMag_loc, n_r_ic_maxMag, nold, &
              &                  nexp, nimp, l_allocate_exp=.true.)
      end if

      call dwdt_dist%initialize(1, n_mlo_loc, n_r_max, nold, nexp, nimp)
      if ( (.not. l_double_curl) .or. l_RMS ) then
         call dpdt_dist%initialize(1, n_mlo_loc, n_r_max, nold, nexp, nimp)
      end if
      call dzdt_dist%initialize(1, n_mlo_loc, n_r_max, nold, nexp, nimp)


      if ( l_finite_diff .and. fd_order==2 .and. fd_order_bound==2 ) then
         n_fields=3
         if ( l_mag ) n_fields=n_fields+2
         allocate( dflowdt_Rdist_container(n_lm_loc,nRstart:nRstop,1:n_fields) )
         dwdt_Rdist(1:,nRstart:) => dflowdt_Rdist_container(1:n_lm_loc,nRstart:nRstop,1)
         dzdt_Rdist(1:,nRstart:) => dflowdt_Rdist_container(1:n_lm_loc,nRstart:nRstop,2)
         dsdt_Rdist(1:,nRstart:) => dflowdt_Rdist_container(1:n_lm_loc,nRstart:nRstop,3)
         if ( l_mag ) then
            dbdt_Rdist(1:,nRstart:) => dflowdt_Rdist_container(1:n_lm_loc,nRstart:nRstop,4)
            djdt_Rdist(1:,nRstart:) => dflowdt_Rdist_container(1:n_lm_loc,nRstart:nRstop,5)
         end if
         allocate(dpdt_Rdist(n_lm_loc,nRstart:nRstop))
         allocate(dVxVhLM_Rdist(n_lm_loc,nRstart:nRstop))
         allocate(dVSrLM_Rdist(n_lm_loc,nRstart:nRstop))
         allocate(dVxBhLM_Rdist(n_lm_loc,nRstart:nRstop))
         bytes_allocated = bytes_allocated+                                  &
         &                 6*n_lm_loc*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX+ &
         &                 3*n_lmMag_loc*(nRstopMag-nRstartMag+1)*SIZEOF_DEF_COMPLEX
      else
         if ( l_double_curl ) then
            allocate( dflowdt_Rdist_container(n_lm_loc,nRstart:nRstop,1:4) )
            dwdt_Rdist(1:,nRstart:) => dflowdt_Rdist_container(1:n_lm_loc,nRstart:nRstop,1)
            dzdt_Rdist(1:,nRstart:) => dflowdt_Rdist_container(1:n_lm_loc,nRstart:nRstop,2)
            dpdt_Rdist(1:,nRstart:) => dflowdt_Rdist_container(1:n_lm_loc,nRstart:nRstop,3)
            dVxVhLM_Rdist(1:,nRstart:) => &
            &                         dflowdt_Rdist_container(1:n_lm_loc,nRstart:nRstop,4)
            bytes_allocated = bytes_allocated+ &
            &                 4*n_lm_loc*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
         else
            allocate( dflowdt_Rdist_container(n_lm_loc,nRstart:nRstop,1:3) )
            dwdt_Rdist(1:,nRstart:) => dflowdt_Rdist_container(1:n_lm_loc,nRstart:nRstop,1)
            dzdt_Rdist(1:,nRstart:) => dflowdt_Rdist_container(1:n_lm_loc,nRstart:nRstop,2)
            dpdt_Rdist(1:,nRstart:) => dflowdt_Rdist_container(1:n_lm_loc,nRstart:nRstop,3)
            allocate( dVxVhLM_Rdist(1:1,1:1) )
            bytes_allocated = bytes_allocated+ &
            &                 3*n_lm_loc*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
         end if

         allocate( dsdt_Rdist_container(n_lm_loc,nRstart:nRstop,1:2) )
         dsdt_Rdist(1:,nRstart:)   => dsdt_Rdist_container(1:n_lm_loc,nRstart:nRstop,1)
         dVSrLM_Rdist(1:,nRstart:) => dsdt_Rdist_container(1:n_lm_loc,nRstart:nRstop,2)
         bytes_allocated = bytes_allocated+ &
         &                 2*n_lm_loc*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX

         ! the magnetic part
         allocate( dbdt_Rdist_container(n_lmMag_loc,nRstartMag:nRstopMag,1:3) )
         dbdt_Rdist(1:,nRstartMag:) => &
         &                    dbdt_Rdist_container(1:n_lmMag_loc,nRstartMag:nRstopMag,1)
         djdt_Rdist(1:,nRstartMag:) => &
         &                    dbdt_Rdist_container(1:n_lmMag_loc,nRstartMag:nRstopMag,2)
         dVxBhLM_Rdist(1:,nRstartMag:)=> &
         &                    dbdt_Rdist_container(1:n_lmMag_loc,nRstartMag:nRstopMag,3)
         bytes_allocated = bytes_allocated+ &
         &                 3*n_lmMag_loc*(nRstopMag-nRstartMag+1)*SIZEOF_DEF_COMPLEX
      end if

      if ( l_chemical_conv ) then
         allocate( dxidt_Rdist_container(n_lm_loc,nRstart:nRstop,1:2) )
         dxidt_Rdist(1:,nRstart:)   => dxidt_Rdist_container(1:n_lm_loc,nRstart:nRstop,1)
         dVXirLM_Rdist(1:,nRstart:) => dxidt_Rdist_container(1:n_lm_loc,nRstart:nRstop,2)
         bytes_allocated = bytes_allocated+ &
         &                 2*n_lm_loc*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
      else
         allocate( dxidt_Rdist_container(1,1,1:2) )
         dxidt_Rdist(1:,1:)   => dxidt_Rdist_container(1:1,1:1,1)
         dVXirLM_Rdist(1:,1:) => dxidt_Rdist_container(1:1,1:1,2)
      end if

      !-- Set the initial values to zero
      if ( l_mag ) then
         dbdt_Rdist(:,:)   =zero
         djdt_Rdist(:,:)   =zero
         dVxBhLM_Rdist(:,:)=zero
      end if
      dwdt_Rdist(:,:)=zero
      dzdt_Rdist(:,:)=zero
      dsdt_Rdist(:,:)=zero
      dpdt_Rdist(:,:)=zero
      dVSrLM_Rdist(:,:)=zero
      if ( l_double_curl ) dVxVhLM_Rdist(:,:)=zero
      if ( l_chemical_conv ) then
         dxidt_Rdist(:,:)  =zero
         dVXirLM_Rdist(:,:)=zero
      end if

      ! The same arrays, but now the LM local part
      if ( l_finite_diff .and. fd_order==2 .and. fd_order_bound==2 ) then
         n_fields=3
         if ( l_mag ) n_fields=n_fields+2
         allocate(dflowdt_LMdist_container(1:n_mlo_loc,n_r_max,1:n_fields,1:nexp))
         dwdt_dist%expl(1:,1:,1:) => dflowdt_LMdist_container(1:n_mlo_loc,1:n_r_max,1,1:nexp)
         dzdt_dist%expl(1:,1:,1:) => dflowdt_LMdist_container(1:n_mlo_loc,1:n_r_max,2,1:nexp)
         dsdt_dist%expl(1:,1:,1:) => dflowdt_LMdist_container(1:n_mlo_loc,1:n_r_max,3,1:nexp)
         bytes_allocated = bytes_allocated+3*n_mlo_loc*n_r_max*nexp*SIZEOF_DEF_COMPLEX
         if ( l_mag ) then
            dbdt_dist%expl(1:,1:,1:) => dflowdt_LMdist_container(1:n_mlo_loc,1:n_r_max,4,1:nexp)
            djdt_dist%expl(1:,1:,1:) => dflowdt_LMdist_container(1:n_mlo_loc,1:n_r_max,5,1:nexp)
            bytes_allocated = bytes_allocated+2*n_mlo_loc*n_r_max*nexp*SIZEOF_DEF_COMPLEX
         end if
         if ( (.not. l_double_curl) .or. l_RMS ) then
            allocate( dpdt_dist%expl(n_mlo_loc,n_r_max,nexp) )
            bytes_allocated = bytes_allocated+n_mlo_loc*n_r_max*nexp*SIZEOF_DEF_COMPLEX
         else
            allocate( dpdt_dist%expl(1,1,1) ) ! To avoid debug
         end if
      else
         if ( l_double_curl ) then
            allocate(dflowdt_LMdist_container(1:n_mlo_loc,n_r_max,1:4,1:nexp))
            dwdt_dist%expl(1:,1:,1:) => dflowdt_LMdist_container(1:n_mlo_loc,1:n_r_max,1,1:nexp)
            dzdt_dist%expl(1:,1:,1:) => dflowdt_LMdist_container(1:n_mlo_loc,1:n_r_max,2,1:nexp)
            dpdt_dist%expl(1:,1:,1:) => dflowdt_LMdist_container(1:n_mlo_loc,1:n_r_max,3,1:nexp)
            dVxVhLM_LMdist(1:,1:,1:) => dflowdt_LMdist_container(1:n_mlo_loc,1:n_r_max,4,1:nexp)
            bytes_allocated = bytes_allocated+4*n_mlo_loc*n_r_max*nexp*SIZEOF_DEF_COMPLEX
         else
            allocate(dflowdt_LMdist_container(1:n_mlo_loc,n_r_max,1:3,1:nexp))
            dwdt_dist%expl(1:,1:,1:) => dflowdt_LMdist_container(1:n_mlo_loc,1:n_r_max,1,1:nexp)
            dzdt_dist%expl(1:,1:,1:) => dflowdt_LMdist_container(1:n_mlo_loc,1:n_r_max,2,1:nexp)
            dpdt_dist%expl(1:,1:,1:) => dflowdt_LMdist_container(1:n_mlo_loc,1:n_r_max,3,1:nexp)
            allocate( dVxVhLM_LMdist(1:1,1:1,1:1) )
            bytes_allocated = bytes_allocated+3*n_mlo_loc*n_r_max*nexp*SIZEOF_DEF_COMPLEX
         end if

         allocate(dsdt_LMdist_container(n_mlo_loc,n_r_max,1:2,1:nexp))
         dsdt_dist%expl(1:,1:,1:) => dsdt_LMdist_container(1:n_mlo_loc,1:n_r_max,1,1:nexp)
         dVSrLM_LMdist(1:,1:,1:) => dsdt_LMdist_container(1:n_mlo_loc,1:n_r_max,2,1:nexp)
         bytes_allocated = bytes_allocated+2*n_mlo_loc*n_r_max*nexp*SIZEOF_DEF_COMPLEX

         allocate(dbdt_LMdist_container(1:n_mloMag_loc,n_r_maxMag,1:3,1:nexp))
         dbdt_dist%expl(1:,1:,1:) => dbdt_LMdist_container(1:n_mloMag_loc,1:n_r_maxMag,1,1:nexp)
         djdt_dist%expl(1:,1:,1:) => dbdt_LMdist_container(1:n_mloMag_loc,1:n_r_maxMag,2,1:nexp)
         dVxBhLM_LMdist(1:,1:,1:) => &
         &                         dbdt_LMdist_container(1:n_mloMag_loc,1:n_r_maxMag,3,1:nexp)
         bytes_allocated = bytes_allocated+ &
         &                 3*nexp*n_mloMag_loc*n_r_maxMag*SIZEOF_DEF_COMPLEX
      end if

      if ( l_chemical_conv ) then
         allocate(dxidt_LMdist_container(1:n_mlo_loc,n_r_max,1:2,1:nexp))
         dxidt_dist%expl(1:,1:,1:)   => dxidt_LMdist_container(1:n_mlo_loc,1:n_r_max,1,1:nexp)
         dVXirLM_LMdist(1:,1:,1:) => dxidt_LMdist_container(1:n_mlo_loc,1:n_r_max,2,1:nexp)
         bytes_allocated = bytes_allocated+2*n_mlo_loc*n_r_max*nexp* &
         &                 SIZEOF_DEF_COMPLEX
      else
         allocate(dxidt_LMdist_container(1,1,1:2,1))
         dxidt_dist%expl(1:,1:,1:)   => dxidt_LMdist_container(1:1,1:1,1,1:)
         dVXirLM_LMdist(1:,1:,1:) => dxidt_LMdist_container(1:1,1:1,2,1:)
      end if

      ! Only when l_dt_cmb_field is requested
      ! There might be a way to allocate only when needed
      allocate ( dbdt_CMB_LMdist(1:n_mloMag_loc) )
      bytes_allocated = bytes_allocated+n_mloMag_loc*SIZEOF_DEF_COMPLEX

   end subroutine initialize_fieldsLast
!-------------------------------------------------------------------------------
   subroutine finalize_fieldsLast
      !
      ! Memory deallocation
      !

      deallocate( dflowdt_Rdist_container, dflowdt_LMdist_container )
      if ( l_finite_diff .and. fd_order==2 .and. fd_order_bound==2 ) then
         deallocate( dpdt_Rdist, dVxVhLM_Rdist, dVxBhLM_Rdist, dVSrLM_Rdist)
      else
         deallocate( dbdt_Rdist_container, dbdt_LMdist_container )
         deallocate( dsdt_Rdist_container, dsdt_LMdist_container )
         if ( .not. l_double_curl ) deallocate( dVxVhLM_Rdist, dVxVhLM_LMdist )
      end if
      deallocate( dbdt_CMB_LMdist )
      deallocate( dxidt_Rdist_container, dxidt_LMdist_container )

      call lorentz_torque_ma_dt%finalize()
      call lorentz_torque_ic_dt%finalize()
      call domega_ma_dt%finalize()
      call domega_ic_dt%finalize()
      call dzdt_dist%finalize()
      if ( .not. l_double_curl .or. l_RMS ) call dpdt_dist%finalize()
      call dwdt_dist%finalize()
      if ( l_heat ) call dsdt_dist%finalize()
      if ( l_chemical_conv ) call dxidt_dist%finalize()
      if ( l_mag ) then
         call dbdt_dist%finalize()
         call djdt_dist%finalize()
      end if
      if ( l_cond_ic ) then
         call dbdt_ic_dist%finalize()
         call djdt_ic_dist%finalize()
      end if

   end subroutine finalize_fieldsLast
!-------------------------------------------------------------------------------
end module fieldsLast
