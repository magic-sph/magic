!$Id$
!********************************************************************
subroutine rBpSpec(time,Tor,TorIC,fileRoot,lIC,map)
   !--------------------------------------------------------------------
   !  Called from rank0, map gives the lm order of Tor and TorIC
   !--------------------------------------------------------------------

   use truncation
   use radial_functions
   use physical_parameters
   use num_param
   use blocking
   use horizontal_data
   use logic
   use output_data
   use useful, only: cc2real
   use LMmapping,only:mappings
 
   implicit none
 
   !-- Input variables:
   real(kind=8),     intent(in) :: time
   complex(kind=8),  intent(in) :: Tor(lm_max,n_r_max)
   complex(kind=8),  intent(in) :: TorIC(lm_max,n_r_ic_max)
   character(len=*), intent(in) :: fileRoot
   logical,          intent(in) :: lIC
   type(mappings),   intent(IN) :: map
 
   !-- Output:
   real(kind=8) :: e_t_AS(l_max,n_r_tot)
   real(kind=8) :: e_t(l_max,n_r_tot)
 
   !-- Local:
   character(len=72) :: specFile
   integer :: n_r,lm,l,m
   real(kind=8) :: fac,rRatio,amp
   real(kind=8) :: e_t_temp
   LOGICAl :: lAS
 
   fac=0.5D0*eScale/(16.D0*datan(1.D0))
 
   do n_r=1,n_r_max
      do l=1,6
         e_t(l,n_r)=0.D0
         e_t_AS(l,n_r) = 0.0D0
      end do
      do lm=2,lm_max
         l=map%lm2l(lm)
         if ( l <= 6 ) then
            m=map%lm2m(lm)
            amp=real(Tor(lm,n_r))
            e_t_temp=dLh(st_map%lm2(l,m))*cc2real(Tor(lm,n_r),m)
            if ( abs(amp)/=0d0 ) then
               if ( m == 0 ) e_t_AS(l,n_r)=fac*amp/abs(amp)*e_t_temp
            end if
            e_t(l,n_r)=e_t(l,n_r)+fac*e_t_temp
         end if
      end do    ! do loop over lms in block
   end do    ! radial grid points
 
   !-- Inner core:
   do n_r=2,n_r_ic_max
      do l=1,6
         e_t_AS(l,n_r_max-1+n_r)=0.D0
         e_t(l,n_r_max-1+n_r)   =0.D0
      end do
   end do
   if ( lIC .and. l_cond_ic ) then
 
      lAS=.true.
      if ( trim(adjustl(fileRoot)) == 'rBrAdvSpec' ) lAS= .false. 
 
      do n_r=2,n_r_ic_max
         rRatio=r_ic(n_r)/r_ic(1)
         do lm=2,lm_max
            l=map%lm2l(lm)
            if ( l <= 6 ) then
               m=map%lm2m(lm)
               if ( m /= 0 .or. lAS ) then
                  e_t_temp= dLh(st_map%lm2(l,m))*rRatio**(2*l+2) &
                       &    * cc2real(TorIC(lm,n_r),m)
                  amp=real(TorIC(lm,n_r))
                  if ( ABS(amp)/=0d0 ) then
                     if ( m == 0 ) e_t_AS(l,n_r_max-1+n_r)= &
                          fac*amp/ABS(amp)*e_t_temp
                  end if
                  e_t(l,n_r_max-1+n_r)=e_t(l,n_r_max-1+n_r)+fac*e_t_temp
               end if
            end if
         end do
      end do
 
   end if
 
   !-- Output into file:
   !     writing l=0/1/2 magnetic energy
   specFile=trim(adjustl(fileRoot))//'.'//tag
   open(91, file=specFile, form='unformatted', status='unknown', &
        position='append')
 
   write(91) sngl(time), (sngl(e_t(1,n_r)),n_r=1,n_r_tot-1),    &
                         (sngl(e_t(2,n_r)),n_r=1,n_r_tot-1),    &
                         (sngl(e_t(3,n_r)),n_r=1,n_r_tot-1),    &
                         (sngl(e_t(4,n_r)),n_r=1,n_r_tot-1),    &
                         (sngl(e_t(5,n_r)),n_r=1,n_r_tot-1),    &
                         (sngl(e_t(6,n_r)),n_r=1,n_r_tot-1)
   write(91) sngl(time), (sngl(e_t_AS(1,n_r)),n_r=1,n_r_tot-1), &
                         (sngl(e_t_AS(2,n_r)),n_r=1,n_r_tot-1), &
                         (sngl(e_t_AS(3,n_r)),n_r=1,n_r_tot-1), &
                         (sngl(e_t_AS(4,n_r)),n_r=1,n_r_tot-1), &
                         (sngl(e_t_AS(5,n_r)),n_r=1,n_r_tot-1), &
                         (sngl(e_t_AS(6,n_r)),n_r=1,n_r_tot-1)
 
   close(91)
 
end subroutine rBpSpec
!----------------------------------------------------------------------------
