#define XSH_COURANT 0
module courant_mod
   !
   ! This module handles the computation of Courant factors on grid space.
   ! It then checks whether the timestep size needs to be modified.
   !

   use parallel_mod
   use precision_mod
   use truncation, only: n_phi_max, n_theta_max, nlat_padded
   use radial_data, only: nRstart, nRstop
   use radial_functions, only: orho1, orho2, or4, or2
   use physical_parameters, only: LFfac, opm
   use num_param, only: delxr2, delxh2
   use horizontal_data, only: osn2
   use logic, only: l_mag, l_mag_LF, l_mag_kin, l_cour_alf_damp
   use useful, only: logWrite
   use constants, only: half, one, two

   implicit none

   private

   integer :: file_handle

   public :: courant, dt_courant, initialize_courant, finalize_courant, &
   &         courant_batch

contains

   subroutine initialize_courant(time, dt, tag)
      !
      ! This subroutine opens the timestep.TAG file which stores the time step
      ! changes of MagIC.
      !

      !-- Input variables
      real(cp),         intent(in) :: time ! time
      real(cp),         intent(in) :: dt   ! time step
      character(len=*), intent(in) :: tag ! trailing of the fime

      if ( rank == 0 ) then
         open(newunit=file_handle, file='timestep.'//tag, status='new')
         write(file_handle, '(1p, es20.12, es16.8)')  time, dt
      end if

   end subroutine initialize_courant
!------------------------------------------------------------------------------
   subroutine finalize_courant()

      if ( rank == 0 ) close(file_handle)

   end subroutine finalize_courant
!------------------------------------------------------------------------------
   subroutine courant(n_r,dtrkc,dthkc,vr,vt,vp,br,bt,bp,courfac,alffac)
      !
      !  courant condition check: calculates Courant
      !  advection lengths in radial direction dtrkc
      !  and in horizontal direction dthkc
      !  on the local radial level n_r
      !
      !  for the effective velocity, the abs. sum of fluid
      !  velocity and Alfven velocity is taken
      !
      !  instead of the full Alfven velocity
      !  a modified Alfven velocity is employed that takes
      !  viscous and Joule damping into account. Different
      !  Courant factors are used for the fluid velocity and
      !  the such modified Alfven velocity
      !
      !

      !-- Input variable:
      integer,  intent(in) :: n_r       ! radial level
      real(cp), intent(in) :: vr(:,:)   ! radial velocity
      real(cp), intent(in) :: vt(:,:)   ! longitudinal velocity
      real(cp), intent(in) :: vp(:,:)   ! azimuthal velocity
      real(cp), intent(in) :: br(:,:)   ! radial magnetic field
      real(cp), intent(in) :: bt(:,:)   ! longitudinal magnetic field
      real(cp), intent(in) :: bp(:,:)   ! azimuthal magnetic field
      real(cp), intent(in) :: courfac
      real(cp), intent(in) :: alffac

      !-- Output:
      real(cp), intent(inout) :: dtrkc    ! Courant step (based on radial advection)
                                          ! for the range of points covered
      real(cp), intent(inout) :: dthkc    ! Courant step based on horizontal advection

      !-- Local  variables:
      integer :: n_theta, n_theta_nhs, n_phi
      real(cp) :: valri2,valhi2,valh2,valh2m,vr2max,vh2max
      real(cp) :: valr,valr2,vflr2,vflh2,cf2,af2

#if (XSH_COURANT==1)
      real(cp) :: dx2_ua, dr2, dh2, dh2n, vflr2max, vflh2max, valr2max, valh2max
      real(cp) :: dtrkc_new,dthkc_new

      vflr2max=0.0_cp
      valr2max=0.0_cp
      vflh2max=0.0_cp
      valh2max=0.0_cp
      cf2=courfac*courfac

      if ( l_mag .and. l_mag_LF .and. .not. l_mag_kin ) then

         af2=alffac*alffac

         !$omp parallel do default(shared) &
         !$omp private(n_theta,n_theta_nhs,n_phi) &
         !$omp private(vflr2,valr,valr2,vflh2,valh2,valh2m) &
         !$omp reduction(max:vflr2max,valr2max,vflh2max,valh2max)
         do n_phi=1,n_phi_max
            do n_theta=1,n_theta_max
               n_theta_nhs=(n_theta+1)/2 ! northern hemisphere=odd n_theta

               vflr2=orho2(n_r)*vr(n_theta,n_phi)*vr(n_theta,n_phi)
               valr =br(n_theta,n_phi)*br(n_theta,n_phi) * &
               &     LFfac*orho1(n_r)
               valr2=valr*valr/(valr+valri2)
               vflr2max=max(vflr2max,or4(n_r)*cf2*vflr2)
               valr2max=max(valr2max,or4(n_r)*af2*valr2)


               vflh2= ( vt(n_theta,n_phi)*vt(n_theta,n_phi) +  &
               &        vp(n_theta,n_phi)*vp(n_theta,n_phi) )* &
               &        osn2(n_theta_nhs)*orho2(n_r)
               valh2= ( bt(n_theta,n_phi)*bt(n_theta,n_phi) +  &
               &        bp(n_theta,n_phi)*bp(n_theta,n_phi) )* &
               &        LFfac*osn2(n_theta_nhs)*orho1(n_r)
               valh2m=valh2*valh2/(valh2+valhi2)
               vflh2max=max(vflh2max,or2(n_r)*cf2*vflh2)
               valh2max=max(valh2max,or2(n_r)*af2*valh2)
            end do
         end do
         !$omp end parallel do

         !-- We must resolve the shortest period of Alfven waves
         if ( l_cour_alf_damp ) then
            dx2_ua=(half*(one+opm))**2/(valr2max+valh2max)
            if ( dx2_ua > delxr2(n_r) )  then
               dr2 = dx2_ua
            else
               dr2 = delxr2(n_r)
            end if

            if ( dx2_ua > delxh2(n_r) )  then
               dh2 = dx2_ua
            else
               dh2 = delxh2(n_r)
            end if
         else
            dr2 = delxr2(n_r)
            dh2 = delxh2(n_r)
         end if

         if ( vflr2max /= 0.0_cp .and. valr2max /= 0.0_cp ) then
            dtrkc_new = min(sqrt(delxr2(n_r)/vflr2max),sqrt(dr2/valr2max))
         else
            dtrkc_new = dtrkc
         end if

         if ( vflh2max /= 0.0_cp .and. valh2max /= 0.0_cp ) then
            dthkc_new = min(sqrt(delxh2(n_r)/vflh2max),sqrt(dh2/valh2max))
         else
            dthkc_new = dthkc
         end if

      else   ! Magnetic field ?

         !$omp parallel do default(shared) &
         !$omp private(n_theta,n_theta_nhs,n_phi,vflr2,vflh2) &
         !$omp reduction(max:vflr2max,vflh2max)
         do n_phi=1,n_phi_max
            do n_theta=1,n_theta_max
               n_theta_nhs=(n_theta+1)/2 ! northern hemisphere=odd n_theta

               vflr2=orho2(n_r)*vr(n_theta,n_phi)*vr(n_theta,n_phi)
               vflr2max=max(vflr2max,cf2*or4(n_r)*vflr2)

               vflh2= ( vt(n_theta,n_phi)*vt(n_theta,n_phi) +  &
               &        vp(n_theta,n_phi)*vp(n_theta,n_phi) )* &
               &        osn2(n_theta_nhs)*orho2(n_r)
               vflh2max=max(vflh2max,cf2*or2(n_r)*vflh2)
            end do
         end do
         !$omp end parallel do

         if ( vflr2max /= 0.0_cp ) then
            dtrkc_new = delxr2(n_r)/vflr2max
         else
            dtrkc_new = dtrkc
         end if

         if ( vflh2max /= 0.0_cp ) then
            dthkc_new = delxh2(n_r)/vflh2max
         else
            dthkc_new = dthkc
         end if

      end if   ! Magnetic field ?

      !$omp critical
      dtrkc=min(dtrkc,dtrkc_new)
      dthkc=min(dthkc,dthkc_new)
      !$omp end critical

#elif ( XSH_COURANT==0)
      if ( l_cour_alf_damp ) then
         valri2=(half*(one+opm))**2/delxr2(n_r)
         valhi2=(half*(one+opm))**2/delxh2(n_r)
      else
         valri2=0.0_cp
         valhi2=0.0_cp
      end if

      vr2max=0.0_cp
      vh2max=0.0_cp
      cf2=courfac*courfac

      if ( l_mag .and. l_mag_LF .and. .not. l_mag_kin ) then

         af2=alffac*alffac

         !$omp parallel do default(shared) &
         !$omp private(n_theta,n_theta_nhs,n_phi) &
         !$omp private(vflr2,valr,valr2,vflh2,valh2,valh2m) &
         !$omp reduction(max:vr2max,vh2max)
         do n_phi=1,n_phi_max
            do n_theta=1,n_theta_max
               n_theta_nhs=(n_theta+1)/2 ! northern hemisphere=odd n_theta

               vflr2=orho2(n_r)*vr(n_theta,n_phi)*vr(n_theta,n_phi)
               valr =br(n_theta,n_phi)*br(n_theta,n_phi)*LFfac*orho1(n_r)
               valr2=valr*valr/(valr+valri2)
               vr2max=max(vr2max,or4(n_r)*(cf2*vflr2+af2*valr2))

               vflh2= ( vt(n_theta,n_phi)*vt(n_theta,n_phi) +  &
               &        vp(n_theta,n_phi)*vp(n_theta,n_phi) )* &
               &        osn2(n_theta_nhs)*orho2(n_r)
               valh2= ( bt(n_theta,n_phi)*bt(n_theta,n_phi) +  &
               &        bp(n_theta,n_phi)*bp(n_theta,n_phi) )* &
               &        LFfac*osn2(n_theta_nhs)*orho1(n_r)
               valh2m=valh2*valh2/(valh2+valhi2)
               vh2max=max(vh2max,or2(n_r)*(cf2*vflh2+af2*valh2m))

            end do
         end do
         !$omp end parallel do

      else   ! Magnetic field ?

         !$omp parallel do default(shared) &
         !$omp private(n_theta,n_theta_nhs,n_phi,vflr2,vflh2) &
         !$omp reduction(max:vr2max,vh2max)
         do n_phi=1,n_phi_max
            do n_theta=1,n_theta_max
               n_theta_nhs=(n_theta+1)/2 ! northern hemisphere=odd n_theta

               vflr2=orho2(n_r)*vr(n_theta,n_phi)*vr(n_theta,n_phi)
               vr2max=max(vr2max,cf2*or4(n_r)*vflr2)

               vflh2= ( vt(n_theta,n_phi)*vt(n_theta,n_phi) +  &
               &        vp(n_theta,n_phi)*vp(n_theta,n_phi) )* &
               &        osn2(n_theta_nhs)*orho2(n_r)
               vh2max=max(vh2max,cf2*or2(n_r)*vflh2)

            end do

         end do
         !$omp end parallel do
      end if   ! Magnetic field ?

      !$omp critical
      if ( vr2max /= 0.0_cp ) dtrkc=min(dtrkc,sqrt(delxr2(n_r)/vr2max))
      if ( vh2max /= 0.0_cp ) dthkc=min(dthkc,sqrt(delxh2(n_r)/vh2max))
      !$omp end critical
#endif

   end subroutine courant
!------------------------------------------------------------------------------
   subroutine courant_batch(n_r_l,n_r_u,dtrkc,dthkc,vr,vt,vp,br,bt,bp,courfac,alffac)
      !
      !  courant condition check: calculates Courant
      !  advection lengths in radial direction dtrkc
      !  and in horizontal direction dthkc
      !  on the local radial level n_r
      !
      !  for the effective velocity, the abs. sum of fluid
      !  velocity and Alfven velocity is taken
      !
      !  instead of the full Alfven velocity
      !  a modified Alfven velocity is employed that takes
      !  viscous and Joule damping into account. Different
      !  Courant factors are used for the fluid velocity and
      !  the such modified Alfven velocity
      !

      !-- Input variable:
      integer,  intent(in) :: n_r_l       ! Starting radial level
      integer,  intent(in) :: n_r_u       ! Stopping radial level
      real(cp), intent(in) :: vr(nlat_padded,n_r_l:n_r_u,*)   ! radial velocity
      real(cp), intent(in) :: vt(nlat_padded,n_r_l:n_r_u,*)   ! longitudinal velocity
      real(cp), intent(in) :: vp(nlat_padded,n_r_l:n_r_u,*)   ! azimuthal velocity
      real(cp), intent(in) :: br(nlat_padded,n_r_l:n_r_u,*)   ! radial magnetic field
      real(cp), intent(in) :: bt(nlat_padded,n_r_l:n_r_u,*)   ! longitudinal magnetic field
      real(cp), intent(in) :: bp(nlat_padded,n_r_l:n_r_u,*)   ! azimuthal magnetic field
      real(cp), intent(in) :: courfac
      real(cp), intent(in) :: alffac

      !-- Output:
      real(cp), intent(inout) :: dtrkc(n_r_l:n_r_u)  ! Courant step (radial)
      real(cp), intent(inout) :: dthkc(n_r_l:n_r_u)  ! Courant step (horizontal)

      !-- Local  variables:
      integer :: n_theta, n_phi, n_theta_nhs, n_r
      real(cp) :: valri2,valhi2,valh2,valh2m
      real(cp) :: vr2max(n_r_l:n_r_u),vh2max(n_r_l:n_r_u)
      real(cp) :: valr,valr2,vflr2,vflh2,cf2,af2

      vr2max=0.0_cp
      vh2max=0.0_cp
      cf2=courfac*courfac

      if ( l_mag .and. l_mag_LF .and. .not. l_mag_kin ) then

         af2=alffac*alffac

         !$omp parallel do default(shared) &
         !$omp private(n_theta,n_r,n_theta_nhs,n_phi) &
         !$omp private(vflr2,valr,valr2,vflh2,valh2,valh2m) &
         !$omp reduction(max:vr2max,vh2max)
         do n_phi=1,n_phi_max
            do n_r=n_r_l,n_r_u
               if ( l_cour_alf_damp ) then
                  valri2=(half*(one+opm))**2/delxr2(n_r)
                  valhi2=(half*(one+opm))**2/delxh2(n_r)
               else
                  valri2=0.0_cp
                  valhi2=0.0_cp
               end if
               do n_theta=1,n_theta_max
                  n_theta_nhs=(n_theta+1)/2 ! northern hemisphere=odd n_theta

                  vflr2=orho2(n_r)*vr(n_theta,n_r,n_phi)*vr(n_theta,n_r,n_phi)
                  valr =br(n_theta,n_r,n_phi)*br(n_theta,n_r,n_phi)*LFfac*orho1(n_r)
                  valr2=valr*valr/(valr+valri2)
                  vr2max(n_r)=max(vr2max(n_r),or4(n_r)*(cf2*vflr2+af2*valr2))

                  vflh2= ( vt(n_theta,n_r,n_phi)*vt(n_theta,n_r,n_phi) +  &
                  &        vp(n_theta,n_r,n_phi)*vp(n_theta,n_r,n_phi) )* &
                  &        osn2(n_theta_nhs)*orho2(n_r)
                  valh2= ( bt(n_theta,n_r,n_phi)*bt(n_theta,n_r,n_phi) +  &
                  &        bp(n_theta,n_r,n_phi)*bp(n_theta,n_r,n_phi) )* &
                  &        LFfac*osn2(n_theta_nhs)*orho1(n_r)
                  valh2m=valh2*valh2/(valh2+valhi2)
                  vh2max(n_r)=max(vh2max(n_r),or2(n_r)*(cf2*vflh2+af2*valh2m))

               end do
            end do
         end do
         !$omp end parallel do

      else   ! Magnetic field ?

         !$omp parallel do default(shared) &
         !$omp private(n_theta,n_r,n_theta_nhs,n_phi,vflr2,vflh2) &
         !$omp reduction(max:vr2max,vh2max)
         do n_phi=1,n_phi_max
            do n_r=n_r_l,n_r_u
               do n_theta=1,n_theta_max
                  n_theta_nhs=(n_theta+1)/2 ! northern hemisphere=odd n_theta

                  vflr2=orho2(n_r)*vr(n_theta,n_r,n_phi)*vr(n_theta,n_r,n_phi)
                  vr2max(n_r)=max(vr2max(n_r),cf2*or4(n_r)*vflr2)

                  vflh2= ( vt(n_theta,n_r,n_phi)*vt(n_theta,n_r,n_phi) +  &
                  &        vp(n_theta,n_r,n_phi)*vp(n_theta,n_r,n_phi) )* &
                  &        osn2(n_theta_nhs)*orho2(n_r)
                  vh2max(n_r)=max(vh2max(n_r),cf2*or2(n_r)*vflh2)

               end do
            end do
         end do
         !$omp end parallel do
      end if   ! Magnetic field ?

      !$omp critical
      do n_r=n_r_l,n_r_u
         if ( vr2max(n_r) /= 0.0_cp ) then
            dtrkc(n_r)=min(dtrkc(n_r),sqrt(delxr2(n_r)/vr2max(n_r)))
         end if
         if ( vh2max(n_r) /= 0.0_cp ) then
            dthkc(n_r)=min(dthkc(n_r),sqrt(delxh2(n_r)/vh2max(n_r)))
         end if
      end do
      !$omp end critical

   end subroutine courant_batch
!------------------------------------------------------------------------------
   subroutine dt_courant(l_new_dt,dt,dt_new,dtMax,dtrkc,dthkc,time)
      !
      !  This subroutine checks if the Courant criterion based on combined
      !  fluid and Alfven velocity is satisfied. It returns a new value for
      !  the time step size.
      !

      !-- Input variables:
      real(cp), intent(in) :: dt ! old time step size
      real(cp), intent(in) :: dtMax ! Maximum time step size
      real(cp), intent(in) :: dtrkc(nRstart:nRstop) ! radial Courant time step as function of radial level
      real(cp), intent(in) :: dthkc(nRstart:nRstop) ! horizontal Courant time step as function of radial level
      real(cp), intent(in) :: time ! Current time

      !-- Output variables:
      logical,  intent(out) :: l_new_dt ! flag indicating that time step is changed (=1) or not (=0)
      real(cp), intent(out) :: dt_new ! new time step size

      !-- Local variables:
      integer :: n_r
      real(cp) :: dtMin, dt_2, dt_fac, dt_rh(2)
      character(len=200) :: message

      dt_fac=two
      dt_rh(:)  =1000.0_cp*dtMax
      do n_r=nRstart,nRstop
         dt_rh(1)=min(dtrkc(n_r),dt_rh(1))
         dt_rh(2)=min(dthkc(n_r),dt_rh(2))
      end do
#ifdef WITH_MPI
      call MPI_Allreduce(MPI_IN_PLACE,dt_rh,2,MPI_DEF_REAL,MPI_MIN, &
           &             MPI_COMM_WORLD,ierr)
#endif

      dtMin=min(dt_rh(1),dt_rh(2))
      dt_2 =min(half*(one/dt_fac+one)*dtMin,dtMax)

      if ( dt > dtMax ) then
         l_new_dt=.true.
         dt_new=dtMax
         write(message,'(1P," ! COURANT: dt=dtMax =",ES12.4,A)') dtMax,&
         &     " ! Think about changing dtMax !"
         call logWrite(message)
      else if ( dt > dtMin ) then
         l_new_dt=.true.
         dt_new  =dt_2
         write(message,'(1P," ! COURANT: dt=",ES11.4," > dt_r=",ES12.4, &
         &            " and dt_h=",ES12.4)') dt,dt_rh(1),dt_rh(2)
         call logWrite(message)
         if ( rank == 0 ) then
            write(file_handle, '(1p, es20.12, es16.8)')  time, dt_new
         end if
      else if ( dt_fac*dt < dtMin .and. dt < dtMax ) then
         l_new_dt=.true.
         dt_new=dt_2
         write(message,'(" ! COURANT: ",F4.1,1P,"*dt=",ES11.4, &
         &          " < dt_r=",ES12.4," and dt_h=",ES12.4)')   &
         &          dt_fac,dt_fac*dt,dt_rh(1),dt_rh(2)
         call logWrite(message)
         if ( rank == 0 ) then
            write(file_handle, '(1p, es20.12, es16.8)')  time, dt_new
         end if
      else
         l_new_dt = .false.
         dt_new = dt
      end if

   end subroutine dt_courant
!------------------------------------------------------------------------------
end module courant_mod
