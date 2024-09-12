module out_movie_IC
   !
   ! This module handles the computation of the inner core contribution to
   ! the movie files. Relevant outputs are stored in frames(*) from the
   ! movie.f90 module.
   !

   use precision_mod
   use parallel_mod
   use communications, only: gather_from_lo_to_rank0
   use truncation, only: minc, lm_maxMag, n_r_ic_maxMag, n_theta_max,  &
       &                 nlat_padded, n_phi_max, lm_max, n_r_ic_max, l_max
   use radial_functions, only: r_ic, r_ICB, O_r_ic2, O_r_ic
   use physical_parameters, only: LFfac
   use horizontal_data, only: n_theta_cal2ord, O_sin_theta
   use logic, only: l_cond_ic
   use movie_data, only: frames, n_movie_field_stop, n_movie_field_start, &
       &                 n_movie_const, n_movie_fields_ic,                &
       &                 n_movie_surface, n_movies, n_movie_field_type,   &
       &                 n_movie_fields
   use out_movie, only: get_fl
   use constants, only: zero, one
   use blocking, only: llmMag, ulmMag
   use sht, only: torpol_to_spat_IC, torpol_to_curl_spat_IC

   implicit none

   private

   public :: store_movie_frame_IC

contains

   subroutine store_movie_frame_IC(bICB,b_ic,db_ic,ddb_ic,aj_ic,dj_ic)
      !
      ! Controls storage of IC magnetic field in movie frame.
      !

      !-- Input of scalar fields:
      complex(cp), intent(in) :: bICB(lm_maxMag)
      complex(cp), intent(in) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: db_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: ddb_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: dj_ic(llmMag:ulmMag,n_r_ic_maxMag)

      !-- Local variables:
      integer :: n_movie, n_field, n_surface, n_const
      integer :: nR, nTheta
      integer :: n_field_type, n_store_last
      integer :: n_field_size, n_fields, n_fields_oc, n_fields_ic

      do n_movie=1,n_movies

         n_fields_oc=n_movie_fields(n_movie)
         n_fields_ic=n_movie_fields_ic(n_movie)
         n_fields   =n_fields_oc+n_fields_ic
         n_surface  =n_movie_surface(n_movie)
         n_const    =n_movie_const(n_movie)

         select case(n_surface)

            case(0) ! 3d
               do n_field=n_fields_oc+1,n_fields
                  n_field_type= n_movie_field_type(n_field,n_movie)
                  n_store_last= n_movie_field_start(n_field,n_movie)-1
                  if ( n_store_last >= 0 ) then
                     call store_fields_3d(bICB, b_ic, db_ic, ddb_ic, aj_ic, dj_ic, &
                          &               n_store_last, n_field_type)
                  end if
               end do  ! Loop over fields

            case(1)  ! r-slice
               if ( n_fields_ic > 0 ) then
                  nR=abs(n_const)
                  do n_field=n_fields_oc+1,n_fields
                     n_field_type=n_movie_field_type(n_field,n_movie)
                     n_store_last=n_movie_field_start(n_field,n_movie)-1
                     if ( n_store_last >= 0 ) then
                        call store_fields_r(bICB, b_ic, db_ic, aj_ic, &
                             &              n_store_last, n_field_type, nR)
                     end if
                  end do  ! Loop over fields
               end if ! if there is anything in the inner core

            case(2) ! Theta-slice or equat
               nTheta =n_const

               do n_field=n_fields_oc+1,n_fields
                  n_field_type=n_movie_field_type(n_field,n_movie)
                  n_store_last=n_movie_field_start(n_field,n_movie)-1
                  if ( n_store_last >= 0 ) then
                     call store_fields_t(bICB, b_ic, db_ic, ddb_ic, aj_ic, dj_ic, &
                          &              n_store_last, n_field_type, nTheta)
                  end if
               end do  ! Loop over fields

            case(3) ! Phi-slice or Phi-avg
               do n_field=n_fields_oc+1,n_fields
                  n_field_type=n_movie_field_type(n_field,n_movie)
                  n_store_last=n_movie_field_start(n_field,n_movie)-1
                  n_field_size=( n_movie_field_stop(n_field,n_movie)-n_store_last )/2
                  if ( n_store_last >= 0 ) then
                     call store_fields_p(bICB, b_ic, db_ic, ddb_ic, aj_ic, dj_ic, &
                          &              n_store_last, n_field_type, n_const,     &
                          &              n_field_size)
                  end if
               end do    ! Loop over fields

         end select  ! Which surface ?

      end do     ! Loop over movies

   end subroutine store_movie_frame_IC
!----------------------------------------------------------------------------
   subroutine store_fields_r(bICB, b_ic, db_ic, aj_ic, n_store_last,  &
              &              n_field_type, nR)
      !
      ! This subroutine stores the frames corresponding to radial cuts movie files
      ! for inner core field.
      !

      !-- Input variables:
      complex(cp), intent(in) :: bICB(lm_maxMag)
      complex(cp), intent(in) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: db_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      integer,     intent(in) :: n_store_last
      integer,     intent(in) :: n_field_type
      integer,     intent(in) :: nR

      !-- Local variables
      integer :: nTheta, nThetaR, nPhi, n_o, n_o_r
      real(cp) :: help
      real(cp) :: BrB(nlat_padded,n_phi_max), BtB(nlat_padded,n_phi_max)
      real(cp) :: BpB(nlat_padded,n_phi_max)
      complex(cp) :: b_ic_r(lm_maxMag), db_ic_r(lm_maxMag), aj_ic_r(lm_maxMag)

      if ( l_cond_ic ) then
         call gather_from_lo_to_rank0(b_ic(:,nR), b_ic_r)
         call gather_from_lo_to_rank0(db_ic(:,nR), db_ic_r)
         call gather_from_lo_to_rank0(aj_ic(:,nR), aj_ic_r)
      else
         db_ic_r(:)=zero
         aj_ic_r(:)=zero
      end if

      if ( rank == 0 ) then
         if ( l_cond_ic ) then
            call torpol_to_spat_IC(r_ic(nR),r_ICB,b_ic_r,db_ic_r,aj_ic_r, &
                 &                 BrB,BtB,BpB)
         else
            call torpol_to_spat_IC(r_ic(nR),r_ICB,bICB,db_ic_r,aj_ic_r, &
                 &                 BrB,BtB,BpB)
         end if

         n_o_r=n_store_last
         if ( n_field_type == 1 ) then
            do nPhi=1,n_phi_max
               do nThetaR=1,n_theta_max
                  nTheta=n_theta_cal2ord(nThetaR)
                  n_o=n_o_r+(nTheta-1)*n_phi_max
                  frames(nPhi+n_o)=BrB(nThetaR,nPhi)*O_r_ic2(nR)
               end do
            end do
         else if ( n_field_type == 2 ) then
            do nPhi=1,n_phi_max
               do nThetaR=1,n_theta_max
                  nTheta=n_theta_cal2ord(nThetaR)
                  n_o=n_o_r+(nTheta-1)*n_phi_max
                  help=O_r_ic(nR)*O_sin_theta(nThetaR)
                  frames(nPhi+n_o)=help*BtB(nThetaR,nPhi)
               end do
            end do
         else if ( n_field_type == 3 ) then
            do nPhi=1,n_phi_max
               do nThetaR=1,n_theta_max
                  nTheta=n_theta_cal2ord(nThetaR)
                  n_o=n_o_r+(nTheta-1)*n_phi_max
                  help=O_r_ic(nR)*O_sin_theta(nThetaR)
                  frames(nPhi+n_o)=help*BpB(nThetaR,nPhi)
               end do
            end do
         end if

      end if ! rank==0 ?

   end subroutine store_fields_r
!----------------------------------------------------------------------------
   subroutine store_fields_p(bICB, b_ic, db_ic, ddb_ic, aj_ic, dj_ic, &
              &              n_store_last, n_field_type, n_phi_const, &
              &              n_field_size)
      !
      ! This subroutine stores the frames corresponding to phi cuts movie files
      ! for inner core field.
      !

      !-- Input variables:
      complex(cp), intent(in) :: bICB(lm_maxMag)
      complex(cp), intent(in) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: db_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: ddb_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: dj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      integer,     intent(in) :: n_store_last
      integer,     intent(in) :: n_field_type
      integer,     intent(in) :: n_phi_const
      integer,     intent(in) :: n_field_size

      !-- Local variables
      integer :: nR, nTheta, nThetaR, nPhi, n_o, n_o_r, nPhi0, nPhi180
      real(cp) :: help, fl(n_theta_max), phi_norm
      real(cp) :: BrB(nlat_padded,n_phi_max), BtB(nlat_padded,n_phi_max)
      real(cp) :: BpB(nlat_padded,n_phi_max), cBrB(nlat_padded,n_phi_max)
      real(cp) :: cBtB(nlat_padded,n_phi_max), cBpB(nlat_padded,n_phi_max)
      complex(cp) :: b_ic_r(lm_maxMag), db_ic_r(lm_maxMag), ddb_ic_r(lm_maxMag)
      complex(cp) :: aj_ic_r(lm_maxMag), dj_ic_r(lm_maxMag)

      phi_norm=one/n_phi_max ! 1 /n_phi_max

      do nR=1,n_r_ic_max

         if ( l_cond_ic ) then
            call gather_from_lo_to_rank0(b_ic(:,nR), b_ic_r)
            call gather_from_lo_to_rank0(db_ic(:,nR), db_ic_r)
            call gather_from_lo_to_rank0(aj_ic(:,nR), aj_ic_r)

            if ( n_field_type == 54 ) then
               call gather_from_lo_to_rank0(dj_ic(:,nR), dj_ic_r)
               call gather_from_lo_to_rank0(ddb_ic(:,nR), ddb_ic_r)
            end if
         else
            db_ic_r(:)=zero
            aj_ic_r(:)=zero
         end if

         if ( rank == 0 ) then
            if ( l_cond_ic ) then
               call torpol_to_spat_IC(r_ic(nR),r_ICB,b_ic_r,db_ic_r,aj_ic_r, &
                    &                 BrB,BtB,BpB)
               if ( n_field_type == 54 ) then
                  call torpol_to_curl_spat_IC(r_ic(nR),r_ICB,db_ic_r,ddb_ic_r, &
                       &                      aj_ic_r,dj_ic_r,cBrB,cBtB,cBpB)
               end if
            else
               call torpol_to_spat_IC(r_ic(nR),r_ICB,bICB,db_ic_r,aj_ic_r, &
                    &                 BrB,BtB,BpB)
               if ( n_field_type == 54 ) then
                  cBrB(:,:)=0.0_cp
                  cBtB(:,:)=0.0_cp
                  cBpB(:,:)=0.0_cp
               end if
            end if

            !------ Get phi no. for left and righty halfspheres:
            nPhi0=n_phi_const
            if ( mod(minc,2) == 1 ) then
               nPhi180=n_phi_max/2+nPhi0
            else
               nPhi180=nPhi0
            end if

            n_o_r=n_store_last+(nR-1)*n_theta_max

            if ( n_field_type == 1 ) then
               do nThetaR=1,n_theta_max
                  nTheta=n_theta_cal2ord(nThetaR)
                  n_o=n_o_r+nTheta
                  frames(n_o)=BrB(nThetaR,nPhi0)*O_r_ic2(nR)
                  frames(n_o+n_field_size)=O_r_ic2(nR)*BrB(nThetaR,nPhi180)
               end do
            else if ( n_field_type == 2 ) then
               do nThetaR=1,n_theta_max
                  nTheta=n_theta_cal2ord(nThetaR)
                  n_o=n_o_r+nTheta
                  frames(n_o)=BtB(nThetaR,nPhi0) *                &
                  &           O_r_ic(nR)*O_sin_theta(nThetaR)
                  frames(n_o+n_field_size)=BtB(nThetaR,nPhi180) * &
                  &             O_r_ic(nR)*O_sin_theta(nThetaR)
              end do
            else if ( n_field_type == 3 ) then
                do nThetaR=1,n_theta_max
                   nTheta=n_theta_cal2ord(nThetaR)
                   n_o=n_o_r+nTheta
                   frames(n_o)=BpB(nThetaR,nPhi0) *                 &
                   &           O_r_ic(nR)*O_sin_theta(nThetaR)
                   frames(n_o+n_field_size)= BpB(nThetaR,nPhi180) * &
                   &                 O_r_ic(nR)*O_sin_theta(nThetaR)
               end do
            else if ( n_field_type == 8 ) then
               call get_fl(fl,nR,.true.)
               do nThetaR=1,n_theta_max
                  nTheta=n_theta_cal2ord(nThetaR)
                  n_o=n_o_r+nTheta
                  !write(*,"(A,I5,A)") "store_movie_IC: frames(",n_o,")"
                  frames(n_o)=fl(nThetaR)
               end do
            else if ( n_field_type == 9 ) then
               do nThetaR=1,n_theta_max
                  nTheta=n_theta_cal2ord(nThetaR)
                  n_o=n_o_r+nTheta
                  help=0.0_cp
                  do nPhi=1,n_phi_max
                     help=help+BpB(nThetaR,nPhi)
                  end do
                  frames(n_o)=phi_norm*help*O_r_ic(nR)* &
                  &           O_sin_theta(nThetaR)
               end do
            else if ( n_field_type == 54 ) then
               help=LFfac*O_r_ic(nR)*O_r_ic2(nR)
               do nThetaR=1,n_theta_max
                  nTheta=n_theta_cal2ord(nThetaR)
                  n_o=n_o_r+nTheta
                  frames(n_o)=        help*O_sin_theta(nThetaR) * &
                  &    ( cBrB(nThetaR,nPhi0)*BtB(nThetaR,nPhi0) - &
                  &      cBtB(nThetaR,nPhi0)*BrB(nThetaR,nPhi0) )
                  frames(n_o+n_field_size)=help*O_sin_theta(nThetaR) * &
                  &    ( cBrB(nThetaR,nPhi180)*BtB(nThetaR,nPhi180) -  &
                  &      cBtB(nThetaR,nPhi180)*BrB(nThetaR,nPhi180) )
               end do
            end if

         end if ! rank == 0?

      end do     ! Loop over radius


   end subroutine store_fields_p
!----------------------------------------------------------------------------
   subroutine store_fields_t(bICB, b_ic, db_ic, ddb_ic, aj_ic, dj_ic, &
              &              n_store_last, n_field_type, nTheta)
      !
      ! This subroutine stores the frames corresponding to theta cuts movie files
      ! for inner core field.
      !

      !-- Input variables:
      complex(cp), intent(in) :: bICB(lm_maxMag)
      complex(cp), intent(in) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: db_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: ddb_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: dj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      integer,     intent(in) :: n_store_last
      integer,     intent(in) :: n_field_type
      integer,     intent(in) :: nTheta

      !-- Local variables
      integer :: nR, nPhi, n_o
      real(cp) :: help
      real(cp) :: BrB(nlat_padded,n_phi_max), BtB(nlat_padded,n_phi_max)
      real(cp) :: BpB(nlat_padded,n_phi_max), cBrB(nlat_padded,n_phi_max)
      real(cp) :: cBtB(nlat_padded,n_phi_max), cBpB(nlat_padded,n_phi_max)
      complex(cp) :: b_ic_r(lm_maxMag), db_ic_r(lm_maxMag), ddb_ic_r(lm_maxMag)
      complex(cp) :: aj_ic_r(lm_maxMag), dj_ic_r(lm_maxMag)

      do nR=1,n_r_ic_max

         if ( l_cond_ic ) then
            call gather_from_lo_to_rank0(b_ic(:,nR), b_ic_r)
            call gather_from_lo_to_rank0(db_ic(:,nR), db_ic_r)
            call gather_from_lo_to_rank0(aj_ic(:,nR), aj_ic_r)

            if ( n_field_type == 14 ) then
               call gather_from_lo_to_rank0(dj_ic(:,nR), dj_ic_r)
               call gather_from_lo_to_rank0(ddb_ic(:,nR), ddb_ic_r)
            end if
         else
            db_ic_r(:)=zero
            aj_ic_r(:)=zero
         end if

         if ( rank == 0 ) then
            if ( l_cond_ic ) then
               call torpol_to_spat_IC(r_ic(nR),r_ICB,b_ic_r,db_ic_r,aj_ic_r, &
                    &                 BrB,BtB,BpB)
               if ( n_field_type == 14 ) then
                  call torpol_to_curl_spat_IC(r_ic(nR),r_ICB,db_ic_r,ddb_ic_r, &
                       &                      aj_ic_r,dj_ic_r,cBrB,cBtB,cBpB)
               end if
            else
               call torpol_to_spat_IC(r_ic(nR),r_ICB,bICB,db_ic_r,aj_ic_r, &
                    &                 BrB,BtB,BpB)
               if ( n_field_type == 14 ) then
                  cBrB(:,:)=0.0_cp
                  cBtB(:,:)=0.0_cp
                  cBpB(:,:)=0.0_cp
               end if
            end if

            n_o=n_store_last+(nR-1)*n_phi_max
            if ( n_field_type == 1 ) then !-- Br
               do nPhi=1,n_phi_max
                  frames(nPhi+n_o)=BrB(nTheta,nPhi)*O_r_ic2(nR)
               end do
            else if ( n_field_type == 2 ) then !-- Btheta
               help=O_r_ic(nR)*O_sin_theta(nTheta)
               do nPhi=1,n_phi_max
                  frames(nPhi+n_o)=help*BtB(nTheta,nPhi)
               end do
            else if ( n_field_type == 3 ) then !-- Bphi
               help=O_r_ic(nR)*O_sin_theta(nTheta)
               do nPhi=1,n_phi_max
                  frames(nPhi+n_o)=help*BpB(nTheta,nPhi)
               end do
            else if ( n_field_type == 13 ) then
               help=-O_r_ic(nR)*O_sin_theta(nTheta)
               do nPhi=1,n_phi_max
                  frames(nPhi+n_o)=help*BtB(nTheta,nPhi)
               end do
            else if ( n_field_type == 14 ) then !-- jtheta
               help=-O_r_ic(nR)*O_sin_theta(nTheta)
               do nPhi=1,n_phi_max
                  frames(nPhi+n_o)=help*cBtB(nTheta,nPhi)
               end do
            end if

         end if ! rank == 0?

      end do     ! Loop over radius

   end subroutine store_fields_t
!----------------------------------------------------------------------------
   subroutine store_fields_3d(bICB, b_ic, db_ic, ddb_ic, aj_ic, dj_ic, &
              &               n_store_last, n_field_type)
      !
      ! This subroutine stores the frames corresponding to 3-D movie files
      ! for inner core field.
      !

      !-- Input variables:
      complex(cp), intent(in) :: bICB(lm_maxMag)
      complex(cp), intent(in) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: db_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: ddb_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: dj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      integer,     intent(in) :: n_store_last
      integer,     intent(in) :: n_field_type

      !-- Local variables
      integer :: nR, nPhi, nThetaR, nTheta, n_o, n_o_r
      real(cp) :: help
      real(cp) :: BrB(nlat_padded,n_phi_max), BtB(nlat_padded,n_phi_max)
      real(cp) :: BpB(nlat_padded,n_phi_max), cBrB(nlat_padded,n_phi_max)
      real(cp) :: cBtB(nlat_padded,n_phi_max), cBpB(nlat_padded,n_phi_max)
      complex(cp) :: b_ic_r(lm_maxMag), db_ic_r(lm_maxMag), ddb_ic_r(lm_maxMag)
      complex(cp) :: aj_ic_r(lm_maxMag), dj_ic_r(lm_maxMag)

      do nR=1,n_r_ic_max

         if ( l_cond_ic ) then
            call gather_from_lo_to_rank0(b_ic(:,nR), b_ic_r)
            call gather_from_lo_to_rank0(db_ic(:,nR), db_ic_r)
            call gather_from_lo_to_rank0(aj_ic(:,nR), aj_ic_r)

            if ( n_field_type == 54 ) then
               call gather_from_lo_to_rank0(dj_ic(:,nR), dj_ic_r)
               call gather_from_lo_to_rank0(ddb_ic(:,nR), ddb_ic_r)
            end if
         else
            db_ic_r(:)=zero
            aj_ic_r(:)=zero
         end if

         if ( rank == 0 ) then
            if ( l_cond_ic ) then
               call torpol_to_spat_IC(r_ic(nR),r_ICB,b_ic_r,db_ic_r,aj_ic_r, &
                    &                 BrB,BtB,BpB)
               if ( n_field_type == 54 ) then
                  call torpol_to_curl_spat_IC(r_ic(nR),r_ICB,db_ic_r,ddb_ic_r, &
                       &                      aj_ic_r,dj_ic_r,cBrB,cBtB,cBpB)
               end if
            else
               call torpol_to_spat_IC(r_ic(nR),r_ICB,bICB,db_ic_r,aj_ic_r, &
                    &                 BrB,BtB,BpB)
               if ( n_field_type == 54 ) then
                  cBrB(:,:)=0.0_cp
                  cBtB(:,:)=0.0_cp
                  cBpB(:,:)=0.0_cp
               end if
            end if

            !------ Calculate magnetic field on grid points:
            n_o_r=n_store_last + (nR-1)*n_theta_max*n_phi_max

            if ( n_field_type == 1 ) then
               do nPhi=1,n_phi_max
                  do nThetaR=1,n_theta_max
                     nTheta=n_theta_cal2ord(nThetaR)
                     n_o=n_o_r+(nTheta-1)*n_phi_max
                     frames(nPhi+n_o)=BrB(nThetaR,nPhi)*O_r_ic2(nR)
                  end do
               end do
            else if ( n_field_type == 2 ) then
               do nPhi=1,n_phi_max
                  do nThetaR=1,n_theta_max
                     nTheta=n_theta_cal2ord(nThetaR)
                     n_o=n_o_r+(nTheta-1)*n_phi_max
                     help=O_r_ic(nR)*O_sin_theta(nThetaR)
                     frames(nPhi+n_o)=help*BtB(nThetaR,nPhi)
                  end do
               end do
            else if ( n_field_type == 3 ) then
               do nPhi=1,n_phi_max
                  do nThetaR=1,n_theta_max
                     nTheta=n_theta_cal2ord(nThetaR)
                     n_o=n_o_r+(nTheta-1)*n_phi_max
                     help=O_r_ic(nR)*O_sin_theta(nThetaR)
                     frames(nPhi+n_o)=help*BpB(nThetaR,nPhi)
                  end do
               end do
            else if ( n_field_type == 54 ) then
               help=LFfac*O_r_ic(nR)*O_r_ic2(nR)
               do nPhi=1,n_phi_max
                  do nThetaR=1,n_theta_max
                     nTheta=n_theta_cal2ord(nThetaR)
                     n_o=n_o_r+(nTheta-1)*n_phi_max
                     frames(nPhi+n_o)= help*O_sin_theta(nThetaR) * &
                     &    ( cBrB(nThetaR,nPhi)*BtB(nThetaR,nPhi) - &
                     &      cBtB(nThetaR,nPhi)*BrB(nThetaR,nPhi) )
                  end do
               end do
            end if

         end if ! rank == 0?

      end do        ! Loop over radius

   end subroutine store_fields_3d
!----------------------------------------------------------------------------
end module out_movie_IC
