module fields
   !
   !  This module contains the potential fields and their radial
   !  derivatives
   !
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_mloMag_loc, n_mlo_loc, nRstart, nRstop, n_r_max, &
       &                 n_lm_loc, n_lmMag_loc, n_r_ic_max, lm_maxMag,      &
       &                 n_r_ic_maxMag, n_r_loc, n_r_maxMag, nRstartMag,    &
       &                 nRstopMag, fd_order, fd_order_bound
   use logic, only: l_chemical_conv, l_finite_diff, l_mag
   use parallel_mod, only: l_master_rank

   implicit none

   private

   !-- Velocity potentials:
   complex(cp), public, allocatable, target :: flow_LMdist_container(:,:,:)
   complex(cp), public, allocatable, target :: flow_Rdist_container(:,:,:)
   complex(cp), public, allocatable, target :: press_LMdist_container(:,:,:)
   complex(cp), public, allocatable, target :: press_Rdist_container(:,:,:)
   complex(cp), public, pointer :: w_LMdist(:,:),dw_LMdist(:,:),ddw_LMdist(:,:)
   complex(cp), public, pointer :: w_Rdist(:,:), dw_Rdist(:,:), ddw_Rdist(:,:)

   complex(cp), public, pointer :: z_LMdist(:,:),dz_LMdist(:,:)
   complex(cp), public, pointer :: z_Rdist(:,:), dz_Rdist(:,:)

   !-- Entropy:
   complex(cp), public, allocatable, target :: s_LMdist_container(:,:,:)
   complex(cp), public, allocatable, target :: s_Rdist_container(:,:,:)
   complex(cp), public, pointer :: s_LMdist(:,:), ds_LMdist(:,:)
   complex(cp), public, pointer :: s_Rdist(:,:), ds_Rdist(:,:)

   !-- Chemical composition:
   complex(cp), public, allocatable, target :: xi_LMdist_container(:,:,:)
   complex(cp), public, allocatable, target :: xi_Rdist_container(:,:,:)
   complex(cp), public, pointer :: xi_LMdist(:,:), dxi_LMdist(:,:)
   complex(cp), public, pointer :: xi_Rdist(:,:), dxi_Rdist(:,:)

   !-- Pressure:
   complex(cp), public, pointer :: p_LMdist(:,:), dp_LMdist(:,:)
   complex(cp), public, pointer :: p_Rdist(:,:), dp_Rdist(:,:)

   !-- Magnetic field potentials:
   complex(cp), public, allocatable :: bICB(:)
   complex(cp), public, allocatable, target :: field_LMdist_container(:,:,:)
   complex(cp), public, allocatable, target :: field_Rdist_container(:,:,:)
   complex(cp), public, pointer :: b_LMdist(:,:), db_LMdist(:,:), ddb_LMdist(:,:)
   complex(cp), public, pointer :: b_Rdist(:,:), db_Rdist(:,:), ddb_Rdist(:,:)
   complex(cp), public, pointer :: aj_LMdist(:,:), dj_LMdist(:,:), ddj_LMdist(:,:)
   complex(cp), public, pointer :: aj_Rdist(:,:), dj_Rdist(:,:)

   !-- Magnetic field potentials in inner core:
   complex(cp), public, allocatable :: b_ic(:,:)
   complex(cp), public, allocatable :: db_ic(:,:)
   complex(cp), public, allocatable :: ddb_ic(:,:)
   complex(cp), public, allocatable :: aj_ic(:,:)
   complex(cp), public, allocatable :: dj_ic(:,:)
   complex(cp), public, allocatable :: ddj_ic(:,:)
   complex(cp), public, allocatable :: b_ic_LMdist(:,:)
   complex(cp), public, allocatable :: db_ic_LMdist(:,:)
   complex(cp), public, allocatable :: ddb_ic_LMdist(:,:)
   complex(cp), public, allocatable :: aj_ic_LMdist(:,:)
   complex(cp), public, allocatable :: dj_ic_LMdist(:,:)
   complex(cp), public, allocatable :: ddj_ic_LMdist(:,:)

   complex(cp), public, allocatable :: work_LMdist(:,:) ! Needed in update routines
   
   !-- Rotation rates:
   real(cp), public :: omega_ic,omega_ma

   public :: initialize_fields, finalize_fields

contains

   subroutine initialize_fields

      integer :: n_fields

      !-- Velocity potentials:
      if ( l_master_rank ) then
         allocate( bICB(lm_maxMag) )
         bytes_allocated = bytes_allocated + lm_maxMag*SIZEOF_DEF_COMPLEX
         allocate( b_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( db_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( ddb_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( aj_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( dj_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( ddj_ic(lm_maxMag,n_r_ic_maxMag) )
         bytes_allocated = bytes_allocated + &
         &                 6*lm_maxMag*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX

      else
         allocate( bICB(1) )
         bytes_allocated = bytes_allocated + SIZEOF_DEF_COMPLEX
         allocate( b_ic(1,n_r_ic_maxMag) )
         allocate( db_ic(1,n_r_ic_maxMag) )
         allocate( ddb_ic(1,n_r_ic_maxMag) )
         allocate( aj_ic(1,n_r_ic_maxMag) )
         allocate( dj_ic(1,n_r_ic_maxMag) )
         allocate( ddj_ic(1,n_r_ic_maxMag) )
         bytes_allocated = bytes_allocated + 6*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX
      end if

      if ( l_finite_diff .and. fd_order==2 .and. fd_order_bound==2 ) then
         n_fields = 3
         if ( l_mag ) n_fields = n_fields+2
         allocate( flow_LMdist_container(n_mlo_loc,n_r_max,1:n_fields) )
         w_LMdist(1:,1:)   => flow_LMdist_container(1:n_mlo_loc,1:n_r_max,1)
         z_LMdist(1:,1:)   => flow_LMdist_container(1:n_mlo_loc,1:n_r_max,2)
         s_LMdist(1:,1:)   => flow_LMdist_container(1:n_mlo_loc,1:n_r_max,3)
         if ( l_mag ) then
            b_LMdist(1:,1:)   => flow_LMdist_container(1:n_mlo_loc,1:n_r_max,4)
            aj_LMdist(1:,1:)   => flow_LMdist_container(1:n_mlo_loc,1:n_r_max,5)
         end if
         allocate( dw_LMdist(n_mlo_loc,n_r_max), ddw_LMdist(n_mlo_loc,n_r_max) )
         allocate( dz_LMdist(n_mlo_loc,n_r_max), ds_LMdist(n_mlo_loc,n_r_max) )
         allocate( db_LMdist(n_mloMag_loc,n_r_maxMag) )
         allocate( ddb_LMdist(n_mloMag_loc,n_r_maxMag) )
         allocate( dj_LMdist(n_mloMag_loc,n_r_maxMag) )
         allocate( ddj_LMdist(n_mloMag_loc,n_r_maxMag) )

         allocate( flow_Rdist_container(n_lm_loc,nRstart:nRstop,1:n_fields) )
         w_Rdist(1:,nRstart:) => flow_Rdist_container(1:n_lm_loc,nRstart:nRstop,1)
         z_Rdist(1:,nRstart:) => flow_Rdist_container(1:n_lm_loc,nRstart:nRstop,2)
         s_Rdist(1:,nRstart:) => flow_Rdist_container(1:n_lm_loc,nRstart:nRstop,3)
         if ( l_mag ) then
            b_Rdist(1:,nRstart:) => flow_Rdist_container(1:n_lm_loc,nRstart:nRstop,4)
            aj_Rdist(1:,nRstart:)=> flow_Rdist_container(1:n_lm_loc,nRstart:nRstop,5)
         end if
         allocate(dw_Rdist(n_lm_loc,nRstart:nRstop), ddw_Rdist(n_lm_loc,nRstart:nRstop))
         allocate(dz_Rdist(n_lm_loc,nRstart:nRstop), ds_Rdist(n_lm_loc,nRstart:nRstop))
         allocate(db_Rdist(n_lmMag_loc,nRstartMag:nRstopMag))
         allocate(ddb_Rdist(n_lmMag_loc,nRstartMag:nRstopMag))
         allocate(dj_Rdist(n_lmMag_loc,nRstartMag:nRstopMag))
      else
         allocate( flow_LMdist_container(n_mlo_loc,n_r_max,1:5) )
         w_LMdist(1:,1:)   => flow_LMdist_container(1:n_mlo_loc,1:n_r_max,1)
         dw_LMdist(1:,1:)  => flow_LMdist_container(1:n_mlo_loc,1:n_r_max,2)
         ddw_LMdist(1:,1:) => flow_LMdist_container(1:n_mlo_loc,1:n_r_max,3)
         z_LMdist(1:,1:)   => flow_LMdist_container(1:n_mlo_loc,1:n_r_max,4)
         dz_LMdist(1:,1:)  => flow_LMdist_container(1:n_mlo_loc,1:n_r_max,5)

         allocate( flow_Rdist_container(n_lm_loc,nRstart:nRstop,1:5) )
         w_Rdist(1:,nRstart:)   => flow_Rdist_container(1:n_lm_loc,nRstart:nRstop,1)
         dw_Rdist(1:,nRstart:)  => flow_Rdist_container(1:n_lm_loc,nRstart:nRstop,2)
         ddw_Rdist(1:,nRstart:) => flow_Rdist_container(1:n_lm_loc,nRstart:nRstop,3)
         z_Rdist(1:,nRstart:)   => flow_Rdist_container(1:n_lm_loc,nRstart:nRstop,4)
         dz_Rdist(1:,nRstart:)  => flow_Rdist_container(1:n_lm_loc,nRstart:nRstop,5)

         !-- Entropy:
         allocate( s_LMdist_container(n_mlo_loc,n_r_max,1:2) )
         s_LMdist(1:,1:)  => s_LMdist_container(1:n_mlo_loc,1:n_r_max,1)
         ds_LMdist(1:,1:) => s_LMdist_container(1:n_mlo_loc,1:n_r_max,2)
         allocate( s_Rdist_container(n_lm_loc,nRstart:nRstop,1:2) )
         s_Rdist(1:,nRstart:)  => s_Rdist_container(1:n_lm_loc,nRstart:nRstop,1)
         ds_Rdist(1:,nRstart:) => s_Rdist_container(1:n_lm_loc,nRstart:nRstop,2)

         !-- Magnetic field potentials:
         allocate( field_LMdist_container(n_mloMag_loc,n_r_maxMag,1:6) )
         b_LMdist(1:,1:)   => field_LMdist_container(1:n_mloMag_loc,1:n_r_maxMag,1)
         db_LMdist(1:,1:)  => field_LMdist_container(1:n_mloMag_loc,1:n_r_maxMag,2)
         ddb_LMdist(1:,1:) => field_LMdist_container(1:n_mloMag_loc,1:n_r_maxMag,3)
         aj_LMdist(1:,1:)  => field_LMdist_container(1:n_mloMag_loc,1:n_r_maxMag,4)
         dj_LMdist(1:,1:)  => field_LMdist_container(1:n_mloMag_loc,1:n_r_maxMag,5)
         ddj_LMdist(1:,1:) => field_LMdist_container(1:n_mloMag_loc,1:n_r_maxMag,6)

         allocate( field_Rdist_container(n_lmMag_loc,nRstart:nRstop,1:5) )
         b_Rdist(1:,nRstart:)   => field_Rdist_container(1:n_lmMag_loc,nRstart:nRstop,1)
         db_Rdist(1:,nRstart:)  => field_Rdist_container(1:n_lmMag_loc,nRstart:nRstop,2)
         ddb_Rdist(1:,nRstart:) => field_Rdist_container(1:n_lmMag_loc,nRstart:nRstop,3)
         aj_Rdist(1:,nRstart:)  => field_Rdist_container(1:n_lmMag_loc,nRstart:nRstop,4)
         dj_Rdist(1:,nRstart:)  => field_Rdist_container(1:n_lmMag_loc,nRstart:nRstop,5)
      end if

      allocate( press_LMdist_container(n_mlo_loc,n_r_max,1:2) )
      p_LMdist(1:,1:)   => press_LMdist_container(1:n_mlo_loc,1:n_r_max,1)
      dp_LMdist(1:,1:)  => press_LMdist_container(1:n_mlo_loc,1:n_r_max,2)

      allocate( press_Rdist_container(n_lm_loc,nRstart:nRstop,1:2) )
      p_Rdist(1:,nRstart:)   => press_Rdist_container(1:n_lm_loc,nRstart:nRstop,1)
      dp_Rdist(1:,nRstart:)  => press_Rdist_container(1:n_lm_loc,nRstart:nRstop,2)


      bytes_allocated = bytes_allocated + 9*n_mlo_loc*n_r_max*SIZEOF_DEF_COMPLEX
      bytes_allocated = bytes_allocated + 6*n_mloMag_loc*n_r_maxMag*SIZEOF_DEF_COMPLEX
      bytes_allocated = bytes_allocated + 9*n_lm_loc*n_r_loc*SIZEOF_DEF_COMPLEX
      bytes_allocated = bytes_allocated + 5*n_lmMag_loc*n_r_loc*SIZEOF_DEF_COMPLEX

      !-- Chemical composition:
      if ( l_chemical_conv ) then
         allocate( xi_LMdist_container(1:n_mlo_loc,n_r_max,1:2) )
         xi_LMdist(1:,1:)  => xi_LMdist_container(1:n_mlo_loc,1:n_r_max,1)
         dxi_LMdist(1:,1:) => xi_LMdist_container(1:n_mlo_loc,1:n_r_max,2)
         allocate( xi_Rdist_container(n_lm_loc,nRstart:nRstop,1:2) )
         xi_Rdist(1:,nRstart:)  => xi_Rdist_container(1:n_lm_loc,nRstart:nRstop,1)
         dxi_Rdist(1:,nRstart:) => xi_Rdist_container(1:n_lm_loc,nRstart:nRstop,2)
         bytes_allocated = bytes_allocated + &
         &                 2*(n_lm_loc)*n_r_max*SIZEOF_DEF_COMPLEX
         bytes_allocated = bytes_allocated + &
         &                 2*n_lm_loc*n_r_loc*SIZEOF_DEF_COMPLEX
      else
         allocate( xi_LMdist_container(1,1,2) ) ! For debugging
         xi_LMdist(1:,1:)  => xi_LMdist_container(1:1,1:1,1)
         dxi_LMdist(1:,1:) => xi_LMdist_container(1:1,1:1,2)
         allocate( xi_Rdist_container(1,1,2) )
         xi_Rdist(1:,1:)   => xi_Rdist_container(1:1,1:1,1)
         dxi_Rdist(1:,1:)  => xi_Rdist_container(1:1,1:1,2)
      end if

      !-- Magnetic field potentials in inner core:
      !   NOTE: n_r-dimension may be smaller once CHEBFT is adopted
      !         for even chebs
      allocate( b_ic_LMdist(1:n_mloMag_loc,n_r_ic_maxMag) )
      allocate( db_ic_LMdist(1:n_mloMag_loc,n_r_ic_maxMag) )
      allocate( ddb_ic_LMdist(1:n_mloMag_loc,n_r_ic_maxMag) )
      allocate( aj_ic_LMdist(1:n_mloMag_loc,n_r_ic_maxMag) )
      allocate( dj_ic_LMdist(1:n_mloMag_loc,n_r_ic_maxMag) )
      allocate( ddj_ic_LMdist(1:n_mloMag_loc,n_r_ic_maxMag) )

      bytes_allocated = bytes_allocated + &
      &                 6*n_mloMag_loc*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX

      allocate( work_LMdist(1:n_mlo_loc,1:n_r_max) )
      bytes_allocated = bytes_allocated + n_mlo_loc*n_r_max*SIZEOF_DEF_COMPLEX

   end subroutine initialize_fields
!----------------------------------------------------------------------------
   subroutine finalize_fields

      deallocate( bICB, b_ic, db_ic, ddb_ic, aj_ic, dj_ic, ddj_ic )
      deallocate( press_LMdist_container, press_Rdist_container )
      deallocate( flow_Rdist_container, flow_LMdist_container )
      if ( l_finite_diff .and. fd_order==2 .and. fd_order_bound==2 ) then
         deallocate( dw_LMdist, ddw_LMdist, dz_LMdist, ds_LMdist)
         deallocate( db_LMdist, ddb_LMdist, dj_LMdist, ddj_LMdist)
         deallocate( dw_Rdist, ddw_Rdist, dz_Rdist, ds_Rdist)
         deallocate( db_Rdist, ddb_Rdist, dj_Rdist)
      else
         deallocate( s_LMdist_container, s_Rdist_container )
         deallocate( field_LMdist_container, field_Rdist_container )
      end if
      deallocate( b_ic_LMdist, db_ic_LMdist, ddb_ic_LMdist, aj_ic_LMdist )
      deallocate( dj_ic_LMdist, ddj_ic_LMdist )
      deallocate( xi_LMdist_container, xi_Rdist_container )
      deallocate( work_LMdist )

   end subroutine finalize_fields
!----------------------------------------------------------------------------
end module fields
