module out_movie_IC

   use precision_mod
   use truncation, only: minc, lm_maxMag, n_r_maxMag, n_r_ic_maxMag, &
       &                 n_phi_max, lm_max, n_r_ic_max, l_max,       &
       &                 n_theta_max, l_axi, n_r_icb, nlat_padded
   use radial_functions, only: r_ic, r_ICB, O_r_ic2, O_r_ic
   use physical_parameters, only: LFfac
   use horizontal_data, only: n_theta_cal2ord, O_sin_theta
   use logic, only: l_cond_ic
   use movie_data, only: frames, n_movie_field_stop, n_movie_field_start, &
       &                 n_movie_type, n_movie_const, n_movie_fields_ic,  &
       &                 n_movie_surface, n_movies, n_movie_field_type,   &
       &                 n_movie_fields
   use out_movie, only: get_fl
   use constants, only: one
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
      complex(cp), intent(in) :: b_ic(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: db_ic(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: ddb_ic(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: aj_ic(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: dj_ic(lm_maxMag,n_r_ic_maxMag)

      !-- Local variables:
      integer :: n_movie        ! No. of movie
      integer :: n_field        ! No. of field
      integer :: n_type         ! Movie type
      integer :: n_surface      ! Surface (1=r,2=theta,3=phi)
      integer :: n_const        ! Gives surface
      integer :: nR
      integer :: n_field_type   ! Numbers field types
      integer :: n_store_last   ! Position i in frame(i) were field starts
      integer :: nTheta,nThetaR
      integer :: nPhi,nPhi0,nPhi180
      integer :: n_field_size,n_fields,n_fields_oc,n_fields_ic
      integer :: n_o,n_o_r

      real(cp) :: BrB(nlat_padded,n_phi_max), BtB(nlat_padded,n_phi_max)
      real(cp) :: BpB(nlat_padded,n_phi_max), cBrB(nlat_padded,n_phi_max)
      real(cp) :: cBtB(nlat_padded,n_phi_max), cBpB(nlat_padded,n_phi_max)
      real(cp) :: fl(n_theta_max),help

      real(cp) ::  phi_norm

      phi_norm=one/n_phi_max ! 2 pi /n_phi_max

      do n_movie=1,n_movies

         n_fields_oc=n_movie_fields(n_movie)
         n_fields_ic=n_movie_fields_ic(n_movie)
         n_fields   =n_fields_oc+n_fields_ic
         n_surface  =n_movie_surface(n_movie)
         n_const    =n_movie_const(n_movie)
         n_type     =n_movie_type(n_movie)

         if ( n_surface == 0 ) then

            do nR=1,n_r_ic_max

               if ( l_cond_ic ) then
                  call torpol_to_spat_IC(r_ic(nR),r_ICB,b_ic(:,nR),  &
                       &                 db_ic(:,nR),aj_ic(:,nR),BrB,BtB,BpB)
                  call torpol_to_curl_spat_IC(r_ic(nR),r_ICB,db_ic(:,nR),  &
                       &                      ddb_ic(:,nR),aj_ic(:,nR),    &
                       &                      dj_ic(:,nR),cBrB,cBtB,cBpB)
               else
                  call torpol_to_spat_IC(r_ic(nR),r_ICB,bICB(:),db_ic(:,1),&
                       &                 aj_ic(:,1), BrB, BtB, BpB)
                  call torpol_to_curl_spat_IC(r_ic(nR),r_ICB,db_ic(:,1),   &
                       &                      ddb_ic(:,1),aj_ic(:,1),      &
                       &                      dj_ic(:,1),cBrB,cBtB,cBpB)
               end if

               !------ Calculate magnetic field on grid points:
               do n_field=n_fields_oc+1,n_fields
                  n_field_type= n_movie_field_type(n_field,n_movie)
                  n_store_last= n_movie_field_start(n_field,n_movie)-1

                  if ( n_store_last >= 0 ) then
                     n_o_r=n_store_last + (nR-2)*n_theta_max*n_phi_max

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

                  end if

               end do  ! Loop over fields

            end do        ! Loop over radius

         else if ( n_surface == 1 ) then  ! r-slice

            if ( n_fields_ic > 0 ) then
               nR=abs(n_const)

               if ( l_cond_ic ) then
                  call torpol_to_spat_IC(r_ic(nR),r_ICB,b_ic(:, nR),db_ic(:, nR),&
                       &                 aj_ic(:, nR),BrB,BtB,BpB)
               else
                  call torpol_to_spat_IC(r_ic(nR),r_ICB,bICB(:),db_ic(:,1),  &
                       &                 aj_ic(:,1),BrB,BtB,BpB)
               end if

               !------ Calculate magnetic field on grid points:
               do n_field=n_fields_oc+1,n_fields
                  n_field_type=n_movie_field_type(n_field,n_movie)
                  n_store_last=n_movie_field_start(n_field,n_movie)-1

                  if ( n_store_last >= 0 ) then
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
                  end if

               end do  ! Loop over fields

            end if ! if there anything in the inner core

         else if ( n_surface == 2 ) then ! Theta-slice or equat

            do nR=1,n_r_ic_max

               if ( l_cond_ic ) then
                  call torpol_to_spat_IC(r_ic(nR),r_ICB,b_ic(:,nR),db_ic(:,nR),&
                       &                 aj_ic(:,nR),BrB,BtB,BpB)
                  call torpol_to_curl_spat_IC(r_ic(nR),r_ICB,db_ic(:,nR),  &
                       &                      ddb_ic(:,nR),aj_ic(:,nR),    &
                       &                      dj_ic(:,nR),cBrB,cBtB,cBpB)
               else
                  call torpol_to_spat_IC(r_ic(nR),r_ICB,bICB(:),db_ic(:,1),&
                       &                 aj_ic(:,1), BrB, BtB, BpB)
                  call torpol_to_curl_spat_IC(r_ic(nR),r_ICB,db_ic(:,1),   &
                       &                      ddb_ic(:,1),aj_ic(:,1),      &
                       &                      dj_ic(:,1),cBrB,cBtB,cBpB)
               end if

               nTheta =n_const
               nThetaR=n_const

               do n_field=n_fields_oc+1,n_fields

                  n_field_type=n_movie_field_type(n_field,n_movie)
                  n_store_last=n_movie_field_start(n_field,n_movie)-1

                  if ( n_store_last >= 0 ) then
                     n_o=n_store_last+(nR-2)*n_phi_max
                     if ( n_field_type == 1 ) then !-- Br
                        do nPhi=1,n_phi_max
                           frames(nPhi+n_o)=BrB(nThetaR,nPhi)*O_r_ic2(nR)
                        end do
                     else if ( n_field_type == 2 ) then !-- Btheta
                        help=O_r_ic(nR)*O_sin_theta(nTheta)
                        do nPhi=1,n_phi_max
                           frames(nPhi+n_o)=help*BtB(nThetaR,nPhi)
                        end do
                     else if ( n_field_type == 3 ) then !-- Bphi
                        help=O_r_ic(nR)*O_sin_theta(nTheta)
                        do nPhi=1,n_phi_max
                           frames(nPhi+n_o)=help*BpB(nThetaR,nPhi)
                        end do
                     else if ( n_field_type == 13 ) then
                        help=-O_r_ic(nR)*O_sin_theta(nTheta)
                        do nPhi=1,n_phi_max
                           frames(nPhi+n_o)=help*BtB(nThetaR,nPhi)
                        end do
                     else if ( n_field_type == 14 ) then !-- jtheta
                        help=-O_r_ic(nR)*O_sin_theta(nTheta)
                        do nPhi=1,n_phi_max
                           frames(nPhi+n_o)=help*cBtB(nThetaR,nPhi)
                        end do
                     end if

                  end if

               end do  ! Loop over fields

            end do     ! Loop over radius


         else if ( abs(n_surface) == 3 ) then ! Phi-slice or Phi-avg

            do nR=1,n_r_ic_max

               if ( l_cond_ic ) then
                  call torpol_to_spat_IC(r_ic(nR),r_ICB,b_ic(:,nR),  &
                       &                 db_ic(:,nR),aj_ic(:,nR),BrB,BtB,BpB)
                  call torpol_to_curl_spat_IC(r_ic(nR),r_ICB,db_ic(:,nR),  &
                       &                      ddb_ic(:,nR),aj_ic(:,nR),    &
                       &                      dj_ic(:,nR),cBrB,cBtB,cBpB)
               else
                  call torpol_to_spat_IC(r_ic(nR),r_ICB,bICB(:),db_ic(:,1),&
                       &                 aj_ic(:,1), BrB, BtB, BpB)
                  call torpol_to_curl_spat_IC(r_ic(nR),r_ICB,db_ic(:,1),   &
                       &                      ddb_ic(:,1),aj_ic(:,1),      &
                       &                      dj_ic(:,1),cBrB,cBtB,cBpB)
               end if

               !------ Get phi no. for left and righty halfspheres:
               nPhi0=n_const
               if ( mod(minc,2) == 1 ) then
                  nPhi180=n_phi_max/2+nPhi0
               else
                  nPhi180=nPhi0
               end if

               !------ Calculate magnetic field on grid points:
               if ( n_type == 30 ) then
                  !------ get_fl returns field for field line plot:
                  call get_fl(fl,nR,.true.)
               end if

               do n_field=n_fields_oc+1,n_fields

                  n_field_type=n_movie_field_type(n_field,n_movie)
                  n_store_last=n_movie_field_start(n_field,n_movie)-1
                  n_field_size=( n_movie_field_stop(n_field,n_movie)-n_store_last )/2

                  if ( n_store_last >= 0 ) then
                     n_o_r=n_store_last+(nR-2)*n_theta_max

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
                           frames(n_o)=            help*O_sin_theta(nThetaR) * &
                           &    ( cBrB(nThetaR,nPhi180)*BtB(nThetaR,nPhi180) - &
                           &      cBtB(nThetaR,nPhi180)*BrB(nThetaR,nPhi180) )
                        end do
                     else
                        do nThetaR=1,n_theta_max
                           nTheta=n_theta_cal2ord(nThetaR)
                           n_o=n_o_r+nTheta
                           frames(n_o)=0.0_cp
                           frames(n_o+n_field_size)=0.0_cp
                        end do
                     end if
                  end if

               end do    ! Loop over fields

            end do          ! Loop over r

         end if  ! Which surface ?

      end do     ! Loop over movies

   end subroutine store_movie_frame_IC
!----------------------------------------------------------------------------
end module out_movie_IC
