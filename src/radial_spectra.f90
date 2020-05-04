module radial_spectra

   use precision_mod
   use parallel_mod
   use communications, only: reduce_radial
   use LMmapping, only: map_mlo
   use truncation, only: n_r_max, n_r_ic_max, n_mlo_loc, n_r_tot, n_r_icb
   use radial_functions, only: or2, r_icb, r_ic
   use num_param, only: eScale
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

   subroutine rBrSpec(time,Pol,PolIC,fileRoot,lIC)

      !-- Input variables
      real(cp),         intent(in) :: time
      complex(cp),      intent(in) :: Pol(n_mlo_loc,n_r_max)
      complex(cp),      intent(in) :: PolIC(n_mlo_loc,n_r_ic_max)
      character(len=*), intent(in) :: fileRoot
      logical,          intent(in) :: lIC
    
      !-- Output to file:
      real(cp) :: e_p_AS(6,n_r_tot), e_p_AS_global(6,n_r_tot)
      real(cp) :: e_p(6,n_r_tot), e_p_global(6,n_r_tot)
    
      !-- Local:
      character(len=72) :: specFile
      integer :: n_r,lm,l,m
      real(cp) :: fac,O_r_icb_E_2,rRatio,amp
      real(cp) :: e_p_temp,dLh
      logical :: lAS
    

      fac=half*eScale/(four*pi)
    
      do n_r=1,n_r_max
         ! setting zero
         e_p(1:6,n_r)   =0.0_cp
         e_p_AS(1:6,n_r)=0.0_cp
    
         do lm=1,n_mlo_loc
            l=map_mlo%i2l(lm)
            if ( l == 0 ) cycle
            if ( l <= 6 ) then
               m=map_mlo%i2m(lm)
               dLh = real(l*(l+1),cp)
               amp=real(Pol(lm,n_r))
               e_p_temp=dLh**2 *or2(n_r)*cc2real(Pol(lm,n_r),m)
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
            do lm=1,n_mlo_loc
               l=map_mlo%i2l(lm)
               if ( l == 0 ) cycle
               if ( l <= 6 ) then
                  dLh = real(l*(l+1),cp)
                  m=map_mlo%i2m(lm)
                  if ( m /= 0 .or. lAS ) then
                     if ( l_cond_ic ) then
                        e_p_temp=dLh*rRatio**(2*l)*dLh*O_r_icb_E_2*    &
                        &        cc2real(PolIC(lm,n_r),m)
                        amp=real(PolIC(lm,n_r))
                     else
                        e_p_temp=dLh*O_r_icb_E_2*rRatio**(2*l) * &
                        &        dLh*cc2real(PolIC(lm,n_r_ICB),m)
                        amp=real(Pol(lm,n_r_ICB))
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
      
      if ( l_master_rank ) then

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
   subroutine rBpSpec(time,Tor,TorIC,fileRoot,lIC)
      !
      !  Called from rank0, map gives the lm order of Tor and TorIC
      !

      !-- Input variables:
      real(cp),         intent(in) :: time
      complex(cp),      intent(in) :: Tor(n_mlo_loc,n_r_max)
      complex(cp),      intent(in) :: TorIC(n_mlo_loc,n_r_ic_max)
      character(len=*), intent(in) :: fileRoot
      logical,          intent(in) :: lIC
    
      !-- Output:
      real(cp) :: e_t_AS(6,n_r_tot), e_t_AS_global(6,n_r_tot)
      real(cp) :: e_t(6,n_r_tot), e_t_global(6,n_r_tot)
    
      !-- Local:
      character(len=72) :: specFile
      integer :: n_r,lm,l,m
      real(cp) :: fac,rRatio,amp,e_t_temp,dLh
      LOGICAl :: lAS
    
      fac=half*eScale/(four*pi)
    
      do n_r=1,n_r_max
         do l=1,6
            e_t(l,n_r)   =0.0_cp
            e_t_AS(l,n_r)=0.0_cp
         end do
         do lm=1,n_mlo_loc
            l=map_mlo%i2l(lm)
            if ( l == 0 ) cycle
            if ( l <= 6 ) then
               m=map_mlo%i2m(lm)
               dLh = real(l*(l+1),cp)
               amp=real(Tor(lm,n_r))
               e_t_temp=dLh*cc2real(Tor(lm,n_r),m)
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
            do lm=1,n_mlo_loc
               l=map_mlo%i2l(lm)
               if ( l == 0 ) cycle
               if ( l <= 6 ) then
                  m=map_mlo%i2m(lm)
                  dLh = real(l*(l+1),cp)
                  if ( m /= 0 .or. lAS ) then
                     e_t_temp= dLh*rRatio**(2*l+2)*cc2real(TorIC(lm,n_r),m)
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
      
      if ( l_master_rank ) then
    
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
