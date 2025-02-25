module radial_spectra

   use precision_mod
   use parallel_mod
   use communications, only: reduce_radial
   use truncation, only: n_r_max, n_r_ic_max
   use radial_data, only: n_r_icb
   use radial_functions, only: or2, r_icb, r_ic
   use num_param, only: eScale
   use blocking, only: llm, ulm
   use logic, only: l_cond_ic
   use output_data, only: tag
   use useful, only: cc2real
   use LMmapping, only: mappings
   use constants, only: pi, one, four, half

   implicit none

   private

   integer :: fileHandle

   public :: rBrSpec, rBpSpec

contains

   subroutine rBrSpec(time,Pol,PolIC,fileRoot,lIC,map)

      !-- Input variables
      real(cp),         intent(in) :: time
      complex(cp),      intent(in) :: Pol(llm:ulm,n_r_max)
      complex(cp),      intent(in) :: PolIC(llm:ulm,n_r_ic_max)
      character(len=*), intent(in) :: fileRoot
      logical,          intent(in) :: lIC
      type(mappings),   intent(in) :: map

      !-- Output to file:
      real(cp) :: e_p_AS(6,n_r_max+n_r_ic_max), e_p_AS_global(6,n_r_max+n_r_ic_max)
      real(cp) :: e_p(6,n_r_max+n_r_ic_max), e_p_global(6,n_r_max+n_r_ic_max)

      !-- Local:
      character(len=72) :: specFile
      integer :: n_r,lm,l,m,n_r_tot
      real(cp) :: fac,O_r_icb_E_2,rRatio,amp,e_p_temp
      logical :: lAS


      fac=half*eScale/(four*pi)
      n_r_tot=n_r_ic_max+n_r_max
      e_p(:,:)   =0.0_cp
      e_p_AS(:,:)=0.0_cp

      do n_r=1,n_r_max
         do lm=llm,ulm
            l=map%lm2l(lm)
            if ( l > 0 .and. l <= 6 ) then
               m=map%lm2m(lm)
               amp=real(Pol(lm,n_r))
               e_p_temp=real(l*(l+1),cp)**2 *or2(n_r)*cc2real(Pol(lm,n_r),m)
               if ( m == 0 ) then
                  if ( abs(amp)/=0.0_cp ) then
                     e_p_AS(l,n_r)=fac*amp/abs(amp)*e_p_temp
                  end if
               end if
               e_p(l,n_r)=e_p(l,n_r)+fac*e_p_temp
            end if
         end do    ! do loop over lms in block
      end do    ! radial grid points

      !-- Inner core:
      if ( lIC ) then

         lAS=.true.
         if ( trim(adjustl(fileRoot)) == 'rBrAdvSpec' ) lAS= .false.

         O_r_icb_E_2=one/r_icb**2

         do n_r=2,n_r_ic_max
            rRatio=r_ic(n_r)/r_ic(1)
            do l=1,6
               e_p(l,n_r_max-1+n_r)=0.0_cp
               e_p_AS(l,n_r_max-1+n_r)=0.0_cp
            end do
            do lm=llm,ulm
               l=map%lm2l(lm)
               if ( l > 0 .and. l <= 6 ) then
                  m=map%lm2m(lm)
                  if ( m /= 0 .or. lAS ) then
                     if ( l_cond_ic ) then
                        e_p_temp=real(l*(l+1),cp)*rRatio**(2*l) * &
                        &        real(l*(l+1),cp)*O_r_icb_E_2*    &
                        &        cc2real(PolIC(lm,n_r),m)
                        amp=real(PolIC(lm,n_r))
                     else
                        e_p_temp=real(l*(l+1),cp)*O_r_icb_E_2*rRatio**(2*l) * &
                        &        real(l*(l+1),cp)*cc2real(Pol(lm,n_r_icb),m)
                        amp=real(Pol(lm,n_r_icb))
                     end if
                     if ( m == 0 ) then
                        if ( abs(amp) /= 0.0_cp) then
                           e_p_AS(l,n_r_max-1+n_r)= fac*amp/abs(amp)*e_p_temp
                        end if
                     end if
                     e_p(l,n_r_max-1+n_r)=e_p(l,n_r_max-1+n_r) + fac*e_p_temp
                  end if
               end if
            end do
         end do
      else
         do n_r=2,n_r_ic_max
            do l=1,6
               e_p_AS(l,n_r_max-1+n_r)=0.0_cp
               e_p(l,n_r_max-1+n_r)   =0.0_cp
            end do
         end do
      end if

      call reduce_radial(e_p, e_p_global, 0)
      call reduce_radial(e_p_AS, e_p_AS_global, 0)

      if ( rank == 0 ) then

         !-- Output into file:
         !     writing l=0/1/2 magnetic energy
         specFile=trim(adjustl(fileRoot))//'.'//tag
         open(newunit=fileHandle, file=specFile, form='unformatted', &
         &    status='unknown', position='append')

         write(fileHandle) real(time,kind=outp),                                &
         &                (real(e_p_global(1,n_r),kind=outp),n_r=1,n_r_tot-1),  &
         &                (real(e_p_global(2,n_r),kind=outp),n_r=1,n_r_tot-1),  &
         &                (real(e_p_global(3,n_r),kind=outp),n_r=1,n_r_tot-1),  &
         &                (real(e_p_global(4,n_r),kind=outp),n_r=1,n_r_tot-1),  &
         &                (real(e_p_global(5,n_r),kind=outp),n_r=1,n_r_tot-1),  &
         &                (real(e_p_global(6,n_r),kind=outp),n_r=1,n_r_tot-1)
         write(fileHandle) real(time,kind=outp),                                  &
         &                (real(e_p_AS_global(1,n_r),kind=outp),n_r=1,n_r_tot-1), &
         &                (real(e_p_AS_global(2,n_r),kind=outp),n_r=1,n_r_tot-1), &
         &                (real(e_p_AS_global(3,n_r),kind=outp),n_r=1,n_r_tot-1), &
         &                (real(e_p_AS_global(4,n_r),kind=outp),n_r=1,n_r_tot-1), &
         &                (real(e_p_AS_global(5,n_r),kind=outp),n_r=1,n_r_tot-1), &
         &                (real(e_p_AS_global(6,n_r),kind=outp),n_r=1,n_r_tot-1)

         close(fileHandle)

      end if

   end subroutine rBrSpec
!----------------------------------------------------------------------------
   subroutine rBpSpec(time,Tor,TorIC,fileRoot,lIC,map)
      !
      !  Called from rank0, map gives the lm order of Tor and TorIC
      !

      !-- Input variables:
      real(cp),         intent(in) :: time
      complex(cp),      intent(in) :: Tor(llm:ulm,n_r_max)
      complex(cp),      intent(in) :: TorIC(llm:ulm,n_r_ic_max)
      character(len=*), intent(in) :: fileRoot
      logical,          intent(in) :: lIC
      type(mappings),   intent(in) :: map

      !-- Output:
      real(cp) :: e_t_AS(6,n_r_max+n_r_ic_max), e_t_AS_global(6,n_r_max+n_r_ic_max)
      real(cp) :: e_t(6,n_r_max+n_r_ic_max), e_t_global(6,n_r_max+n_r_ic_max)

      !-- Local:
      character(len=72) :: specFile
      integer :: n_r,lm,l,m,n_r_tot
      real(cp) :: fac,rRatio,amp,e_t_temp
      logical :: lAS

      fac=half*eScale/(four*pi)
      n_r_tot=n_r_max+n_r_ic_max
      e_t(:,:)   =0.0_cp
      e_t_AS(:,:)=0.0_cp

      do n_r=1,n_r_max
         do lm=llm,ulm
            l=map%lm2l(lm)
            if ( l > 0 .and. l <= 6 ) then
               m=map%lm2m(lm)
               amp=real(Tor(lm,n_r))
               e_t_temp=real(l*(l+1),cp)*cc2real(Tor(lm,n_r),m)
               if ( abs(amp)/=0.0_cp ) then
                  if ( m == 0 ) e_t_AS(l,n_r)=fac*amp/abs(amp)*e_t_temp
               end if
               e_t(l,n_r)=e_t(l,n_r)+fac*e_t_temp
            end if
         end do    ! do loop over lms in block
      end do    ! radial grid points

      !-- Inner core:
      do n_r=2,n_r_ic_max
         do l=1,6
            e_t_AS(l,n_r_max-1+n_r)=0.0_cp
            e_t(l,n_r_max-1+n_r)   =0.0_cp
         end do
      end do
      if ( lIC .and. l_cond_ic ) then

         lAS=.true.
         if ( trim(adjustl(fileRoot)) == 'rBrAdvSpec' ) lAS= .false.

         do n_r=2,n_r_ic_max
            rRatio=r_ic(n_r)/r_ic(1)
            do lm=llm,ulm
               l=map%lm2l(lm)
               if ( l > 0 .and. l <= 6 ) then
                  m=map%lm2m(lm)
                  if ( m /= 0 .or. lAS ) then
                     e_t_temp=real(l*(l+1),cp)*rRatio**(2*l+2) &
                     &        * cc2real(TorIC(lm,n_r),m)
                     amp=real(TorIC(lm,n_r))
                     if ( abs(amp)/=0.0_cp ) then
                        if ( m == 0 ) e_t_AS(l,n_r_max-1+n_r)=fac*amp/abs(amp)*e_t_temp
                     end if
                     e_t(l,n_r_max-1+n_r)=e_t(l,n_r_max-1+n_r)+fac*e_t_temp
                  end if
               end if
            end do
         end do

      end if

      call reduce_radial(e_t, e_t_global, 0)
      call reduce_radial(e_t_AS, e_t_AS_global, 0)

      if ( rank == 0 ) then

         !-- Output into file:
         !     writing l=0/1/2 magnetic energy
         specFile=trim(adjustl(fileRoot))//'.'//tag
         open(newunit=fileHandle, file=specFile, form='unformatted', &
         &    status='unknown', position='append')

         write(fileHandle) real(time,kind=outp),                                  &
         &                (real(e_t_global(1,n_r),kind=outp),n_r=1,n_r_tot-1),    &
         &                (real(e_t_global(2,n_r),kind=outp),n_r=1,n_r_tot-1),    &
         &                (real(e_t_global(3,n_r),kind=outp),n_r=1,n_r_tot-1),    &
         &                (real(e_t_global(4,n_r),kind=outp),n_r=1,n_r_tot-1),    &
         &                (real(e_t_global(5,n_r),kind=outp),n_r=1,n_r_tot-1),    &
         &                (real(e_t_global(6,n_r),kind=outp),n_r=1,n_r_tot-1)
         write(fileHandle) real(time,kind=outp),                                  &
         &                (real(e_t_AS_global(1,n_r),kind=outp),n_r=1,n_r_tot-1), &
         &                (real(e_t_AS_global(2,n_r),kind=outp),n_r=1,n_r_tot-1), &
         &                (real(e_t_AS_global(3,n_r),kind=outp),n_r=1,n_r_tot-1), &
         &                (real(e_t_AS_global(4,n_r),kind=outp),n_r=1,n_r_tot-1), &
         &                (real(e_t_AS_global(5,n_r),kind=outp),n_r=1,n_r_tot-1), &
         &                (real(e_t_AS_global(6,n_r),kind=outp),n_r=1,n_r_tot-1)

         close(fileHandle)

      end if

   end subroutine rBpSpec
!----------------------------------------------------------------------------
end module radial_spectra
