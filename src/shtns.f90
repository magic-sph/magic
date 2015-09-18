module shtns

#ifdef WITH_SHTNS

   use truncation, only: ncp, nrp, m_max, l_max, n_theta_max, n_phi_max, minc, lm_max
   use blocking, only: nfs
   use horizontal_data, only: dLh, gauss, theta_ord
   use radial_functions, only: r
   use parallel_mod
   include 'shtns.f'

contains

   subroutine init_shtns()
      use truncation, only: lm_max
      implicit none

      integer it, ip
      integer nlm
      integer norm

      if (rank .eq. 0) then
          call shtns_verbose(1)
      end if
      call shtns_use_threads(0)

      norm = SHT_ORTHONORMAL + SHT_NO_CS_PHASE

      call shtns_set_size(l_max, m_max/minc, minc, norm)
      call shtns_precompute(SHT_GAUSS, SHT_PHI_CONTIGUOUS, 1.e-12, &
          n_theta_max, n_phi_max)

      call shtns_calc_nlm(nlm, l_max, m_max/minc, minc)
      if (lm_max .ne. nlm) then
          print*, "error: nlm /= lm_max", nlm, lm_max
          call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      end if

   end subroutine

   ! transform a spherical harmonic field into grid space
   subroutine scal_to_spat(Slm, fieldc)
      implicit none
      complex(kind=8), intent(in) :: Slm(:)
      real(kind=8), intent(out) :: fieldc(:)

      call shtns_SH_to_spat(Slm, fieldc)

   end subroutine

   ! transform a scalar spherical harmonic field into it's gradient
   ! on the grid
   subroutine scal_to_grad_spat(Slm, gradtc, gradpc)
      implicit none
      complex(kind=8), intent(in) :: Slm(:)
      real(kind=8), intent(out) :: gradtc(:), gradpc(:)

      call shtns_sph_to_spat(Slm, gradtc, gradpc)

   end subroutine

   subroutine pol_to_grad_spat(Slm, gradtc, gradpc)
      use horizontal_data

      implicit none
      complex(kind=8), intent(in) :: Slm(:)
      real(kind=8), intent(out) :: gradtc(:), gradpc(:)

      ! local
      complex(kind=8), allocatable :: Qlm(:)
      integer :: ip, it, lm

      allocate(Qlm(lm_max))

      do lm = 1, lm_max
          Qlm(lm) = dLh(lm) * Slm(lm)
      end do

      call shtns_sph_to_spat(Qlm, gradtc, gradpc)
      
      deallocate(Qlm)

   end subroutine

   subroutine torpol_to_spat(Wlm, dWlm, Zlm, nR, vrc, vtc, vpc)
      use truncation, only: lm_max
      use horizontal_data

      implicit none
      complex(kind=8), intent(in) :: Wlm(:), dWlm(:), Zlm(:)
      real(kind=8), intent(out) :: vrc(:), vtc(:), vpc(:)
      integer, intent(in) :: nR

      ! local
      complex(kind=8), allocatable :: Qlm(:)
      integer :: ip, it, lm

      allocate(Qlm(lm_max))

      do lm = 1, lm_max
          Qlm(lm) = dLh(lm) * Wlm(lm)
      end do

      call shtns_qst_to_spat(Qlm, dWlm, Zlm, vrc, vtc, vpc)

      deallocate(Qlm)
   end subroutine

   subroutine pol_to_curlr_spat(Qlm, cvrc)
      use truncation, only: lm_max

      implicit none
      complex(kind=8), intent(in) :: Qlm(:)
      real(kind=8), intent(out) :: cvrc(:)

      ! local
      complex(kind=8), allocatable :: dQlm(:)
      integer :: lm

      allocate(dQlm(lm_max))

      do lm = 1, lm_max
          dQlm(lm) = dLh(lm) * Qlm(lm)
      end do

      call shtns_SH_to_spat(dQlm, cvrc)

      deallocate(dQlm)
      end subroutine

   subroutine torpol_to_curl_spat(Blm, dBlm, ddBlm, Jlm, dJlm, nR, &
                                 cvrc, cvtc, cvpc)
      use truncation, only: lm_max
      use horizontal_data, only: sinTheta

      implicit none
      complex(kind=8), intent(in) :: Blm(:), dBlm(:), ddBlm(:)
      complex(kind=8), intent(in) :: Jlm(:), dJlm(:)
      integer, intent(in) :: nR
      real(kind=8), intent(out) :: cvrc(:), cvtc(:), cvpc(:)

      ! local
      complex(kind=8), allocatable :: Qlm(:), Tlm(:)
      integer :: it, ip, lm

      allocate(Qlm(lm_max), Tlm(lm_max))

      do lm = 1, lm_max
          Qlm(lm) = dLh(lm) * Jlm(lm)
          Tlm(lm) = 1/r(nR)**2 * dLh(lm) * Blm(lm) - ddBlm(lm)
      end do

      call shtns_qst_to_spat(Qlm, dJlm, Tlm, cvrc, cvtc, cvpc)

      deallocate(Qlm, Tlm)
   end subroutine

#endif

end module shtns

