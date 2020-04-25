module fieldsLast
   !
   ! This module contains all the work arrays of the previous time-steps
   ! needed to time advance the code.
   ! They are needed in the time-stepping scheme.
   !

   use precision_mod
   use truncation, only: n_r_max, lm_max, n_r_maxMag, lm_maxMag, &
       &                 n_r_ic_maxMag, nRstart, nRstop,         &
       &                 nRstartMag, nRstopMag, n_lm_loc,        &
       &                 n_lmMag_loc, n_mlo_loc, n_mloMag_loc
   use blocking, only: llm, ulm, llmMag, ulmMag
   use logic, only: l_chemical_conv, l_heat, l_mag, l_cond_ic, l_double_curl, &
       &            l_RMS
   use constants, only: zero
   use mem_alloc, only: bytes_allocated
   use time_array

   implicit none

   private

   type(type_tarray), public :: dsdt, dwdt, dpdt, dzdt, dxidt
   type(type_tarray), public :: dbdt, djdt, dbdt_ic, djdt_ic
   type(type_tscalar), public :: domega_ma_dt, domega_ic_dt
   type(type_tscalar), public :: lorentz_torque_ic_dt, lorentz_torque_ma_dt

   !DIR$ ATTRIBUTES ALIGN:64 :: dwdt_Rloc,dzdt_Rloc,dpdt_Rloc,dsdt_Rloc,dVSrLM_Rloc
   complex(cp), public, allocatable, target  :: dflowdt_Rloc_container(:,:,:)
   complex(cp), public, allocatable, target  :: dsdt_Rloc_container(:,:,:)
   complex(cp), public, allocatable, target  :: dxidt_Rloc_container(:,:,:)
   complex(cp), public, allocatable, target  :: dbdt_Rloc_container(:,:,:)
   complex(cp), public, pointer :: dwdt_Rloc(:,:),dzdt_Rloc(:,:)
   complex(cp), public, pointer :: dpdt_Rloc(:,:), dsdt_Rloc(:,:), dVSrLM_Rloc(:,:)
   complex(cp), public, pointer :: dxidt_Rloc(:,:), dVXirLM_Rloc(:,:)
   complex(cp), public, pointer :: dVxVhLM_Rloc(:,:)

   !@> TODO those should supersede the previous Rloc version
   complex(cp), public, allocatable, target  :: dflowdt_Rdist_container(:,:,:)
   complex(cp), public, allocatable, target  :: dsdt_Rdist_container(:,:,:)
   complex(cp), public, allocatable, target  :: dxidt_Rdist_container(:,:,:)
   complex(cp), public, allocatable, target  :: dbdt_Rdist_container(:,:,:)
   complex(cp), public, pointer :: dwdt_Rdist(:,:),dzdt_Rdist(:,:)
   complex(cp), public, pointer :: dpdt_Rdist(:,:), dsdt_Rdist(:,:), dVSrLM_Rdist(:,:)
   complex(cp), public, pointer :: dxidt_Rdist(:,:), dVXirLM_Rdist(:,:)
   complex(cp), public, pointer :: dVxVhLM_Rdist(:,:)

   !DIR$ ATTRIBUTES ALIGN:64 :: djdt_Rloc,dbdt_Rloc,dVxBhLM_Rloc
   complex(cp), public, pointer :: djdt_Rloc(:,:), dVxBhLM_Rloc(:,:)
   complex(cp), public, pointer :: dbdt_Rloc(:,:)

   !@> TODO those should supersede the previous Rloc version
   complex(cp), public, pointer :: djdt_Rdist(:,:), dVxBhLM_Rdist(:,:)
   complex(cp), public, pointer :: dbdt_Rdist(:,:)

   ! The same arrays, but now the LM local part
   complex(cp), public, allocatable, target  :: dflowdt_LMloc_container(:,:,:,:)
   complex(cp), public, allocatable, target  :: dsdt_LMloc_container(:,:,:,:)
   complex(cp), public, allocatable, target  :: dxidt_LMloc_container(:,:,:,:)
   complex(cp), public, allocatable, target  :: dbdt_LMloc_container(:,:,:,:)
   complex(cp), public, pointer :: dVSrLM_LMloc(:,:,:), dVXirLM_LMloc(:,:,:)
   complex(cp), public, pointer :: dVxVhLM_LMloc(:,:,:), dVxBhLM_LMloc(:,:,:)

   complex(cp), public, allocatable :: dbdt_CMB_LMloc(:)

   ! The same arrays, but now the LM local part
   complex(cp), public, allocatable, target  :: dflowdt_LMdist_container(:,:,:,:)
   complex(cp), public, allocatable, target  :: dsdt_LMdist_container(:,:,:,:)
   complex(cp), public, allocatable, target  :: dxidt_LMdist_container(:,:,:,:)
   complex(cp), public, allocatable, target  :: dbdt_LMdist_container(:,:,:,:)
   complex(cp), public, pointer :: dVSrLM_LMdist(:,:,:), dVXirLM_LMdist(:,:,:)
   complex(cp), public, pointer :: dVxVhLM_LMdist(:,:,:), dVxBhLM_LMdist(:,:,:)

   complex(cp), public, allocatable :: dbdt_CMB_LMdist(:)



   public :: initialize_fieldsLast, finalize_fieldsLast, gather_dt_fields

contains

   subroutine initialize_fieldsLast(nold, nexp, nimp)
      !
      ! This routine defines all the arrays needed to time advance MagIC
      !

      integer, intent(in) :: nold ! Number of storage of the old state
      integer, intent(in) :: nexp ! Number of explicit states
      integer, intent(in) :: nimp ! Number of implicit states

      call domega_ma_dt%initialize(nold, nexp, nimp)
      call domega_ic_dt%initialize(nold, nexp, nimp)

      call lorentz_torque_ic_dt%initialize(nold, nexp, nimp)
      call lorentz_torque_ma_dt%initialize(nold, nexp, nimp)

      if ( l_heat ) call dsdt%initialize(llm, ulm, n_r_max, nold, nexp, nimp)
      if ( l_chemical_conv ) call dxidt%initialize(llm, ulm, n_r_max, nold, &
                                  &                nexp, nimp)

      if ( l_mag ) then
         call dbdt%initialize(llmMag, ulmMag, n_r_maxMag, nold, nexp, nimp)
         call djdt%initialize(llmMag, ulmMag, n_r_maxMag, nold, nexp, nimp)
      end if

      if ( l_cond_ic ) then
         call dbdt_ic%initialize(llmMag, ulmMag, n_r_ic_maxMag, nold, &
              &                  nexp, nimp, l_allocate_exp=.true.)
         call djdt_ic%initialize(llmMag, ulmMag, n_r_ic_maxMag, nold, &
              &                  nexp, nimp, l_allocate_exp=.true.)
      end if

      call dwdt%initialize(llm, ulm, n_r_max, nold, nexp, nimp)
      if ( (.not. l_double_curl) .or. l_RMS ) then
         call dpdt%initialize(llm, ulm, n_r_max, nold, nexp, nimp)
      end if
      call dzdt%initialize(llm, ulm, n_r_max, nold, nexp, nimp)

      if ( l_double_curl ) then
         allocate( dflowdt_Rloc_container(lm_max,nRstart:nRstop,1:4) )
         dwdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,1)
         dzdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,2)
         dpdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,3)
         dVxVhLM_Rloc(1:,nRstart:) => &
         &                         dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,4)
         bytes_allocated = bytes_allocated+ &
         &                 4*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX

         allocate( dflowdt_Rdist_container(n_lm_loc,nRstart:nRstop,1:4) )
         dwdt_Rdist(1:,nRstart:) => dflowdt_Rdist_container(1:n_lm_loc,nRstart:nRstop,1)
         dzdt_Rdist(1:,nRstart:) => dflowdt_Rdist_container(1:n_lm_loc,nRstart:nRstop,2)
         dpdt_Rdist(1:,nRstart:) => dflowdt_Rdist_container(1:n_lm_loc,nRstart:nRstop,3)
         dVxVhLM_Rdist(1:,nRstart:) => &
         &                         dflowdt_Rdist_container(1:n_lm_loc,nRstart:nRstop,4)
         bytes_allocated = bytes_allocated+ &
         &                 4*n_lm_loc*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
      else
         allocate( dflowdt_Rloc_container(lm_max,nRstart:nRstop,1:3) )
         dwdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,1)
         dzdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,2)
         dpdt_Rloc(1:,nRstart:) => dflowdt_Rloc_container(1:lm_max,nRstart:nRstop,3)
         allocate( dVxVhLM_Rloc(1:1,1:1) )
         bytes_allocated = bytes_allocated+ &
         &                 3*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX

         allocate( dflowdt_Rdist_container(n_lm_loc,nRstart:nRstop,1:3) )
         dwdt_Rdist(1:,nRstart:) => dflowdt_Rdist_container(1:n_lm_loc,nRstart:nRstop,1)
         dzdt_Rdist(1:,nRstart:) => dflowdt_Rdist_container(1:n_lm_loc,nRstart:nRstop,2)
         dpdt_Rdist(1:,nRstart:) => dflowdt_Rdist_container(1:n_lm_loc,nRstart:nRstop,3)
         allocate( dVxVhLM_Rdist(1:1,1:1) )
         bytes_allocated = bytes_allocated+ &
         &                 3*n_lm_loc*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
      end if

      allocate( dsdt_Rloc_container(lm_max,nRstart:nRstop,1:2) )
      dsdt_Rloc(1:,nRstart:)   => dsdt_Rloc_container(1:lm_max,nRstart:nRstop,1)
      dVSrLM_Rloc(1:,nRstart:) => dsdt_Rloc_container(1:lm_max,nRstart:nRstop,2)
      bytes_allocated = bytes_allocated+ &
      &                 2*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX

      allocate( dsdt_Rdist_container(n_lm_loc,nRstart:nRstop,1:2) )
      dsdt_Rdist(1:,nRstart:)   => dsdt_Rdist_container(1:n_lm_loc,nRstart:nRstop,1)
      dVSrLM_Rdist(1:,nRstart:) => dsdt_Rdist_container(1:n_lm_loc,nRstart:nRstop,2)
      bytes_allocated = bytes_allocated+ &
      &                 2*n_lm_loc*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX

      if ( l_chemical_conv ) then
         allocate( dxidt_Rloc_container(lm_max,nRstart:nRstop,1:2) )
         dxidt_Rloc(1:,nRstart:)   => dxidt_Rloc_container(1:lm_max,nRstart:nRstop,1)
         dVXirLM_Rloc(1:,nRstart:) => dxidt_Rloc_container(1:lm_max,nRstart:nRstop,2)
         bytes_allocated = bytes_allocated+ &
         &                 2*lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX

         allocate( dxidt_Rdist_container(n_lm_loc,nRstart:nRstop,1:2) )
         dxidt_Rdist(1:,nRstart:)   => dxidt_Rdist_container(1:n_lm_loc,nRstart:nRstop,1)
         dVXirLM_Rdist(1:,nRstart:) => dxidt_Rdist_container(1:n_lm_loc,nRstart:nRstop,2)
         bytes_allocated = bytes_allocated+ &
         &                 2*n_lm_loc*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
      else
         allocate( dxidt_Rloc_container(1,1,1:2) )
         dxidt_Rloc(1:,1:)   => dxidt_Rloc_container(1:1,1:1,1)
         dVXirLM_Rloc(1:,1:) => dxidt_Rloc_container(1:1,1:1,2)

         allocate( dxidt_Rdist_container(1,1,1:2) )
         dxidt_Rdist(1:,1:)   => dxidt_Rdist_container(1:1,1:1,1)
         dVXirLM_Rdist(1:,1:) => dxidt_Rdist_container(1:1,1:1,2)
      end if

      ! the magnetic part
      allocate( dbdt_Rloc_container(lm_maxMag,nRstartMag:nRstopMag,1:3) )
      dbdt_Rloc(1:,nRstartMag:) => &
      &                    dbdt_Rloc_container(1:lm_maxMag,nRstartMag:nRstopMag,1)
      djdt_Rloc(1:,nRstartMag:) => &
      &                    dbdt_Rloc_container(1:lm_maxMag,nRstartMag:nRstopMag,2)
      dVxBhLM_Rloc(1:,nRstartMag:)=> &
      &                    dbdt_Rloc_container(1:lm_maxMag,nRstartMag:nRstopMag,3)
      bytes_allocated = bytes_allocated+ &
      &                 3*lm_maxMag*(nRstopMag-nRstartMag+1)*SIZEOF_DEF_COMPLEX

      allocate( dbdt_Rdist_container(n_lmMag_loc,nRstartMag:nRstopMag,1:3) )
      dbdt_Rdist(1:,nRstartMag:) => &
      &                    dbdt_Rdist_container(1:n_lmMag_loc,nRstartMag:nRstopMag,1)
      djdt_Rdist(1:,nRstartMag:) => &
      &                    dbdt_Rdist_container(1:n_lmMag_loc,nRstartMag:nRstopMag,2)
      dVxBhLM_Rdist(1:,nRstartMag:)=> &
      &                    dbdt_Rdist_container(1:n_lmMag_loc,nRstartMag:nRstopMag,3)
      bytes_allocated = bytes_allocated+ &
      &                 3*n_lmMag_loc*(nRstopMag-nRstartMag+1)*SIZEOF_DEF_COMPLEX

      !-- Set the initial values to zero
      if ( l_mag ) then
         dbdt_Rloc(:,:)   =zero
         djdt_Rloc(:,:)   =zero
         dVxBhLM_Rloc(:,:)=zero
      end if
      dwdt_Rloc(:,:)=zero
      dzdt_Rloc(:,:)=zero
      dsdt_Rloc(:,:)=zero
      dpdt_Rloc(:,:)=zero
      dVSrLM_Rloc(:,:)=zero
      if ( l_double_curl ) dVxVhLM_Rloc(:,:)=zero
      if ( l_chemical_conv ) then
         dxidt_Rloc(:,:)  =zero
         dVXirLM_Rloc(:,:)=zero
      end if

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
      if ( l_double_curl ) then
         allocate(dflowdt_LMloc_container(llm:ulm,n_r_max,1:4,1:nexp))
         dwdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,1,1:nexp)
         dzdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,2,1:nexp)
         dpdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,3,1:nexp)
         dVxVhLM_LMloc(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,4,1:nexp)
         bytes_allocated = bytes_allocated+4*(ulm-llm+1)*n_r_max*nexp* &
         &                 SIZEOF_DEF_COMPLEX
      else
         allocate(dflowdt_LMloc_container(llm:ulm,n_r_max,1:3,1:nexp))
         dwdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,1,1:nexp)
         dzdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,2,1:nexp)
         dpdt%expl(llm:,1:,1:) => dflowdt_LMloc_container(llm:ulm,1:n_r_max,3,1:nexp)
         allocate( dVxVhLM_LMloc(1:1,1:1,1:1) )
         bytes_allocated = bytes_allocated+3*(ulm-llm+1)*n_r_max*nexp* &
         &                 SIZEOF_DEF_COMPLEX
      end if

      allocate(dsdt_LMloc_container(llm:ulm,n_r_max,1:2,1:nexp))
      dsdt%expl(llm:,1:,1:) => dsdt_LMloc_container(llm:ulm,1:n_r_max,1,1:nexp)
      dVSrLM_LMloc(llm:,1:,1:) => dsdt_LMloc_container(llm:ulm,1:n_r_max,2,1:nexp)
      bytes_allocated = bytes_allocated+2*(ulm-llm+1)*n_r_max*nexp* &
      &                 SIZEOF_DEF_COMPLEX

      if ( l_chemical_conv ) then
         allocate(dxidt_LMloc_container(llm:ulm,n_r_max,1:2,1:nexp))
         dxidt%expl(llm:,1:,1:)   => dxidt_LMloc_container(llm:ulm,1:n_r_max,1,1:nexp)
         dVXirLM_LMloc(llm:,1:,1:) => dxidt_LMloc_container(llm:ulm,1:n_r_max,2,1:nexp)
         bytes_allocated = bytes_allocated+2*(ulm-llm+1)*n_r_max*nexp* &
         &                 SIZEOF_DEF_COMPLEX
      else
         allocate(dxidt_LMloc_container(1,1,1:2,1))
         dxidt%expl(1:,1:,1:)   => dxidt_LMloc_container(1:1,1:1,1,1:)
         dVXirLM_LMloc(1:,1:,1:) => dxidt_LMloc_container(1:1,1:1,2,1:)
      end if

      allocate(dbdt_LMloc_container(llmMag:ulmMag,n_r_maxMag,1:3,1:nexp))
      dbdt%expl(llmMag:,1:,1:) => dbdt_LMloc_container(llmMag:ulmMag,1:n_r_maxMag,1,1:nexp)
      djdt%expl(llmMag:,1:,1:) => dbdt_LMloc_container(llmMag:ulmMag,1:n_r_maxMag,2,1:nexp)
      dVxBhLM_LMloc(llmMag:,1:,1:) => &
      &                         dbdt_LMloc_container(llmMag:ulmMag,1:n_r_maxMag,3,1:nexp)
      bytes_allocated = bytes_allocated+ &
      &                 3*nexp*(ulmMag-llmMag+1)*n_r_maxMag*SIZEOF_DEF_COMPLEX

      !-- Set the initial values to zero
      if ( l_mag ) dVxBhLM_LMloc(:,:,:)=zero
      dVSrLM_LMloc(:,:,:)=zero
      if ( l_double_curl ) dVxVhLM_LMloc(:,:,:)=zero
      if ( l_chemical_conv ) dVXirLM_LMloc(:,:,:)=zero

      ! Only when l_dt_cmb_field is requested
      ! There might be a way to allocate only when needed
      allocate ( dbdt_CMB_LMloc(llmMag:ulmMag) )
      bytes_allocated = bytes_allocated+(ulmMag-llmMag+1)*SIZEOF_DEF_COMPLEX

      ! The same arrays, but now the LM local part
      if ( l_double_curl ) then
         allocate(dflowdt_LMdist_container(1:n_mlo_loc,n_r_max,1:4,1:nexp))
         dwdt%expl_dist(1:,1:,1:) => dflowdt_LMdist_container(1:n_mlo_loc,1:n_r_max,1,1:nexp)
         dzdt%expl_dist(1:,1:,1:) => dflowdt_LMdist_container(1:n_mlo_loc,1:n_r_max,2,1:nexp)
         dpdt%expl_dist(1:,1:,1:) => dflowdt_LMdist_container(1:n_mlo_loc,1:n_r_max,3,1:nexp)
         dVxVhLM_LMdist(1:,1:,1:) => dflowdt_LMdist_container(1:n_mlo_loc,1:n_r_max,4,1:nexp)
         bytes_allocated = bytes_allocated+4*n_mlo_loc*n_r_max*nexp*SIZEOF_DEF_COMPLEX
      else
         allocate(dflowdt_LMdist_container(1:n_mlo_loc,n_r_max,1:3,1:nexp))
         dwdt%expl_dist(1:,1:,1:) => dflowdt_LMdist_container(1:n_mlo_loc,1:n_r_max,1,1:nexp)
         dzdt%expl_dist(1:,1:,1:) => dflowdt_LMdist_container(1:n_mlo_loc,1:n_r_max,2,1:nexp)
         dpdt%expl_dist(1:,1:,1:) => dflowdt_LMdist_container(1:n_mlo_loc,1:n_r_max,3,1:nexp)
         allocate( dVxVhLM_LMdist(1:1,1:1,1:1) )
         bytes_allocated = bytes_allocated+3*n_mlo_loc*n_r_max*nexp*SIZEOF_DEF_COMPLEX
      end if

      allocate(dsdt_LMdist_container(n_mlo_loc,n_r_max,1:2,1:nexp))
      dsdt%expl_dist(1:,1:,1:) => dsdt_LMdist_container(1:n_mlo_loc,1:n_r_max,1,1:nexp)
      dVSrLM_LMdist(1:,1:,1:) => dsdt_LMdist_container(1:n_mlo_loc,1:n_r_max,2,1:nexp)
      bytes_allocated = bytes_allocated+2*n_mlo_loc*n_r_max*nexp*SIZEOF_DEF_COMPLEX

      if ( l_chemical_conv ) then
         allocate(dxidt_LMdist_container(1:n_mlo_loc,n_r_max,1:2,1:nexp))
         dxidt%expl_dist(1:,1:,1:)   => dxidt_LMdist_container(1:n_mlo_loc,1:n_r_max,1,1:nexp)
         dVXirLM_LMdist(1:,1:,1:) => dxidt_LMdist_container(1:n_mlo_loc,1:n_r_max,2,1:nexp)
         bytes_allocated = bytes_allocated+2*n_mlo_loc*n_r_max*nexp* &
         &                 SIZEOF_DEF_COMPLEX
      else
         allocate(dxidt_LMdist_container(1,1,1:2,1))
         dxidt%expl_dist(1:,1:,1:)   => dxidt_LMdist_container(1:1,1:1,1,1:)
         dVXirLM_LMdist(1:,1:,1:) => dxidt_LMdist_container(1:1,1:1,2,1:)
      end if

      allocate(dbdt_LMdist_container(1:n_mloMag_loc,n_r_maxMag,1:3,1:nexp))
      dbdt%expl_dist(1:,1:,1:) => dbdt_LMdist_container(1:n_mloMag_loc,1:n_r_maxMag,1,1:nexp)
      djdt%expl_dist(1:,1:,1:) => dbdt_LMdist_container(1:n_mloMag_loc,1:n_r_maxMag,2,1:nexp)
      dVxBhLM_LMdist(1:,1:,1:) => &
      &                         dbdt_LMdist_container(1:n_mloMag_loc,1:n_r_maxMag,3,1:nexp)
      bytes_allocated = bytes_allocated+ &
      &                 3*nexp*n_mloMag_loc*n_r_maxMag*SIZEOF_DEF_COMPLEX

      !-- Set the initial values to zero
      if ( l_mag ) dVxBhLM_LMdist(:,:,:)=zero
      dVSrLM_LMdist(:,:,:)=zero
      if ( l_double_curl ) dVxVhLM_LMdist(:,:,:)=zero
      if ( l_chemical_conv ) dVXirLM_LMdist(:,:,:)=zero

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

      deallocate( dflowdt_Rloc_container, dsdt_Rloc_container )
      deallocate( dflowdt_Rdist_container, dsdt_Rdist_container )
      deallocate( dbdt_Rloc_container, dflowdt_LMloc_container )
      deallocate( dsdt_LMloc_container, dbdt_LMloc_container )
      deallocate( dbdt_CMB_LMloc )
      deallocate( dxidt_Rloc_container, dxidt_LMloc_container )
      deallocate( dbdt_Rdist_container, dxidt_Rdist_container )

      deallocate( dflowdt_LMdist_container, dbdt_CMB_LMdist )
      deallocate( dsdt_LMdist_container, dbdt_LMdist_container )
      deallocate( dxidt_LMdist_container )

      if ( .not. l_double_curl ) deallocate( dVxVhLM_Rloc, dVxVhLM_LMloc )
      if ( .not. l_double_curl ) deallocate( dVxVhLM_Rdist, dVxVhLM_LMdist )

      call lorentz_torque_ma_dt%finalize()
      call lorentz_torque_ic_dt%finalize()
      call domega_ma_dt%finalize()
      call domega_ic_dt%finalize()
      call dzdt%finalize()
      if ( .not. l_double_curl .or. l_RMS ) call dpdt%finalize()
      call dwdt%finalize()
      if ( l_heat ) call dsdt%finalize()
      if ( l_chemical_conv ) call dxidt%finalize()
      if ( l_mag ) then
         call dbdt%finalize()
         call djdt%finalize()
      end if
      if ( l_cond_ic ) then
         call dbdt_ic%finalize()
         call djdt_ic%finalize()
      end if

   end subroutine finalize_fieldsLast
!-------------------------------------------------------------------------------
   !@> TODO: delete this subroutine once merging is done
   subroutine gather_dt_fields
      use communications, only: gather_Flm

      integer :: nR

      do nR=nRstart,nRstop
         call gather_Flm(dwdt_Rdist(:,nR), dwdt_Rloc(:,nR))
         call gather_Flm(dzdt_Rdist(:,nR), dzdt_Rloc(:,nR))
         call gather_Flm(dpdt_Rdist(:,nR), dpdt_Rloc(:,nR))
         if ( l_double_curl ) call gather_Flm(dVxVhLM_Rdist(:,nR), dVxVhLM_Rloc(:,nR))
         call gather_Flm(dsdt_Rdist(:,nR), dsdt_Rloc(:,nR))
         call gather_Flm(dVSrLM_Rdist(:,nR), dVSrLM_Rloc(:,nR))
         if ( l_chemical_conv ) then
            call gather_Flm(dxidt_Rdist(:,nR), dxidt_Rloc(:,nR))
            call gather_Flm(dVXirLM_Rdist(:,nR), dVXirLM_Rloc(:,nR))
         end if
         if ( l_mag ) then
            call gather_Flm(dbdt_Rdist(:,nR), dbdt_Rloc(:,nR))
            call gather_Flm(djdt_Rdist(:,nR), djdt_Rloc(:,nR))
            call gather_Flm(dVxBhLM_Rdist(:,nR), dVxBhLM_Rloc(:,nR))
         end if
      end do

   end subroutine gather_dt_fields
!-------------------------------------------------------------------------------
end module fieldsLast
