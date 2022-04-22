module grid_blocking

   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: nlat_padded, n_phi_max, lm_max, lmP_max
   use radial_data, only: nRstart, nRstop

   implicit none
   
   private

   logical :: l_batch

   integer, public :: n_phys_space   ! Number of indices in physical space
   integer, public :: n_spec_space   ! Number of indices in spectral space
   integer, public :: n_spec_space_lmP   ! Number of indices in spectral space

   integer, public, allocatable :: spat2lat(:), spat2lon(:)
   integer, public, allocatable :: spat2rad(:), spec2rad(:)
   integer, public, allocatable :: radlatlon2spat(:,:,:)
   integer, public, allocatable :: spec2lm(:)

   public :: initialize_grid_blocking, finalize_grid_blocking

contains

   subroutine initialize_grid_blocking(l_batched)

      !-- Input variables
      logical, intent(in) :: l_batched

      !-- Local variables
      integer :: nelem, n_phi, n_theta, n_r, lm

      l_batch = l_batched


      if ( .not. l_batched ) then
         n_phys_space = nlat_padded * n_phi_max
         n_spec_space = lm_max !* n_r_max
         n_spec_space_lmP = lmP_max
      else
         n_phys_space = nlat_padded * n_phi_max * (nRstop-nRstart+1)
         n_spec_space = lm_max * (nRstop-nRstart+1)
         n_spec_space_lmP = lmP_max * (nRstop-nRstart+1)
      end if

      allocate( spat2lat(n_phys_space), spat2lon(n_phys_space) )
      allocate( spat2rad(n_phys_space) )
      bytes_allocated = bytes_allocated+(3*n_phys_space)*SIZEOF_INTEGER
      allocate( radlatlon2spat(nlat_padded,n_phi_max,nRstart:nRstop) )
      bytes_allocated = bytes_allocated+nlat_padded*n_phi_max*(nRstop-nRstart+1)* &
      &                 SIZEOF_INTEGER
      allocate( spec2lm(n_spec_space), spec2rad(n_spec_space) )
      bytes_allocated = bytes_allocated+2*n_spec_space*SIZEOF_INTEGER

      if ( .not. l_batched ) then
         do n_r=nRstart,nRstop
            nelem = 1
            do n_phi=1,n_phi_max
               do n_theta=1,nlat_padded
                  spat2lat(nelem)=n_theta
                  spat2lon(nelem)=n_phi
                  !spat2rad(nelem,n_r)=n_r
                  radlatlon2spat(n_theta,n_phi,n_r)=nelem
                  nelem = nelem+1
               end do
            end do
         end do

         do n_r=nRstart,nRstop
            nelem=1
            do lm=1,lm_max
               spec2lm(nelem)=lm
               !spec2rad(nelem,n_r)=n_r
               nelem=nelem+1
            end do
         end do

      else ! Batched transforms with theta layout

         nelem = 1
         do n_phi=1,n_phi_max
            do n_r=nRstart,nRstop
               do n_theta=1,nlat_padded
                  spat2lat(nelem)=n_theta
                  spat2lon(nelem)=n_phi
                  spat2rad(nelem)=n_r
                  radlatlon2spat(n_theta,n_phi,n_r)=nelem
                  nelem = nelem+1
               end do
            end do
         end do

         nelem=1
         do n_r=nRstart,nRstop
            do lm=1,lm_max
               spec2lm(nelem) =lm
               spec2rad(nelem)=n_r
               nelem=nelem+1
            end do
         end do
      end if

   end subroutine initialize_grid_blocking
!---------------------------------------------------------------------------------
   subroutine finalize_grid_blocking

      deallocate( spat2lat, spat2lon, spat2rad, radlatlon2spat )
      deallocate( spec2lm, spec2rad )

   end subroutine finalize_grid_blocking
!---------------------------------------------------------------------------------
end module grid_blocking
