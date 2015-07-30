!$Id$
!********************************************************************
subroutine rBrSpec(time,Pol,PolIC,fileRoot,lIC,map)

   use truncation
   use radial_functions
   use physical_parameters
   use num_param,only: eScale
   use blocking
   use horizontal_data
   use logic
   use output_data
   use useful, only: cc2real
   use LMmapping,only:mappings

   implicit none
 
   real(kind=8),     intent(in) :: time
   complex(kind=8),  intent(in) :: Pol(lm_max,n_r_max)
   complex(kind=8),  intent(in) :: PolIC(lm_max,n_r_ic_max)
   character(len=*), intent(in) :: fileRoot
   logical,          intent(in) :: lIC
   type(mappings),   intent(in) :: map
 
   !-- Output to file:
   real(kind=8) :: e_p_AS(l_max,n_r_tot)
   real(kind=8) :: e_p(l_max,n_r_tot)
 
   !-- Local:
   character(len=72) :: specFile
   integer :: n_r,lm,l,m
   real(kind=8) :: fac,O_r_icb_E_2,rRatio,amp
   real(kind=8) :: e_p_temp
   logical :: lAS
 

   fac=0.5D0*eScale/(16.D0*datan(1.D0))
 
   do n_r=1,n_r_max
      ! setting zero
      e_p(1:6,n_r)=0.0D0
      e_p_AS(1:6,n_r)=0.0D0
 
      do lm=2,lm_max
         l=map%lm2l(lm)
         if ( l <= 6 ) then
            m=map%lm2m(lm)
            amp=real(Pol(lm,n_r))
            e_p_temp=dLh(st_map%lm2(l,m))**2 *or2(n_r)*cc2real(Pol(lm,n_r),m)
            if ( m == 0 ) then
               if ( ABS(amp)/=0.d0 ) then
                  e_p_AS(l,n_r)=fac*amp/ABS(amp)*e_p_temp
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
 
      O_r_icb_E_2=1.d0/r_icb**2
 
      do n_r=2,n_r_ic_max
         rRatio=r_ic(n_r)/r_ic(1)
         do l=1,6
            e_p(l,n_r_max-1+n_r)=0.D0
            e_p_AS(l,n_r_max-1+n_r)=0.D0
         end do
         do lm=2,lm_max
            l=map%lm2l(lm)
            if ( l <= 6 ) then
               m=map%lm2m(lm)
               if ( m /= 0 .or. lAS ) then
                  IF( l_cond_ic ) then
                     e_p_temp=dLh(st_map%lm2(l,m))*rRatio**(2*l) * &
                              dLh(st_map%lm2(l,m))*O_r_icb_E_2*    &
                              cc2real(PolIC(lm,n_r),m)
                     amp=real(PolIC(lm,n_r))
                  else
                     e_p_temp=dLh(st_map%lm2(l,m))*O_r_icb_E_2*rRatio**(2*l) * &
                              dLh(st_map%lm2(l,m))*cc2real(PolIC(lm,n_r_ICB),m)
                     amp=real(Pol(lm,n_r_ICB))
                  end if
                  if ( m == 0 ) then
                     if ( abs(amp) /= 0d0) then
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
            e_p_AS(l,n_r_max-1+n_r)=0.d0
            e_p(l,n_r_max-1+n_r)   =0.d0
         end do
      end do
   end if
   
   !-- Output into file:
   !     writing l=0/1/2 magnetic energy
   specFile=trim(adjustl(fileRoot))//'.'//tag
   open(91, file=specFile, form='unformatted', status='unknown', &
        position='append')
 
   write(91) real(time,kind=4),                       &
        (real(e_p(1,n_r),kind=4),n_r=1,n_r_tot-1),    &
        (real(e_p(2,n_r),kind=4),n_r=1,n_r_tot-1),    &
        (real(e_p(3,n_r),kind=4),n_r=1,n_r_tot-1),    &
        (real(e_p(4,n_r),kind=4),n_r=1,n_r_tot-1),    &
        (real(e_p(5,n_r),kind=4),n_r=1,n_r_tot-1),    &
        (real(e_p(6,n_r),kind=4),n_r=1,n_r_tot-1)
   write(91) real(time,kind=4),                       &
        (real(e_p_AS(1,n_r),kind=4),n_r=1,n_r_tot-1), &
        (real(e_p_AS(2,n_r),kind=4),n_r=1,n_r_tot-1), &
        (real(e_p_AS(3,n_r),kind=4),n_r=1,n_r_tot-1), &
        (real(e_p_AS(4,n_r),kind=4),n_r=1,n_r_tot-1), &
        (real(e_p_AS(5,n_r),kind=4),n_r=1,n_r_tot-1), &
        (real(e_p_AS(6,n_r),kind=4),n_r=1,n_r_tot-1)
 
   close(91)
 
end subroutine rBrSpec
!----------------------------------------------------------------------------
