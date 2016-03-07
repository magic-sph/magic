module readCheckPoints_helper
   use precision_mod
   use truncation, only: n_r_max,lm_max,lm_maxMag, &
                         l_max,m_max,     &
                         minc
   use cosine_transform_odd
   use blocking, only: lmStartB,lmStopB,nLMBs,lm2l,lm2m
   use logic, only: l_heat
   use radial_functions, only: chebt_oc, cheb_norm, chebt_ic, cheb_norm_ic
   use init_fields, only: scale_b, scale_v, scale_s
   implicit none

   integer(lip) :: bytes_allocated=0

   contains

!------------------------------------------------------------------------------
   subroutine getLm2lmO(n_r_max,n_r_max_old,l_max,l_max_old, &
                        m_max,minc,minc_old,inform,lm_max,   &
                        lm_max_old,n_data_oldP,lm2lmo)

      !--- Input variables
      integer, intent(in) :: n_r_max,l_max,m_max,minc
      integer, intent(in) :: n_r_max_old,l_max_old,minc_old
      integer, intent(in) :: inform,lm_max

      !--- Output variables
      integer,intent(out) :: lm2lmo(lm_max)
      integer,intent(out) :: n_data_oldP
      integer,intent(out) :: lm_max_old

      !--- Local variables
      integer :: n_data,n_data_old
      integer :: m_max_old
      integer :: l,m,lm,lmo,lo,mo

      !-- Outer core fields:
      n_data  = lm_max*n_r_max
      !-- This allows to increase the number of grid points by 10!

      if ( l_max==l_max_old .and. minc==minc_old .and. n_r_max==n_r_max_old ) then

         !----- Direct reading of fields, grid not changed:
         write(*,'(/,'' ! Reading fields directly.'')')

         n_data_old=n_data
         if ( inform>2 ) then
            n_data_oldP=n_data
         else
            !----- In the past an 'extra' radial grid point has been
            !      stored which was not really necessary
            n_data_oldP=lm_max*(n_r_max+1)
         end if

         lm_max_old=lm_max
      else

         !----- Mapping onto new grid !
         write(*,'(/,'' ! Mapping onto new grid.'')')

         if ( MOD(minc_old,minc) /= 0 )                                &
              &     write(6,'('' ! Warning: Incompatible old/new minc= '',2i3)')

         m_max_old =(l_max_old/minc_old)*minc_old
         lm_max_old=m_max_old*(l_max_old+1)/minc_old -                &
              &     m_max_old*(m_max_old-minc_old)/(2*minc_old) +     &
              &     l_max_old-m_max_old+1

         n_data_old=lm_max_old*n_r_max_old
         if ( inform>2 ) then
            n_data_oldP=n_data_old
         else
            n_data_oldP=lm_max_old*(n_r_max_old+1)
         end if

         !-- Write info to STdoUT:
         write(*,'('' ! Old/New  l_max= '',2I4,''  m_max= '',2I4,     &
              &       ''  minc= '',2I3,''  lm_max= '',2I5/)')         &
              &           l_max_old,l_max,m_max_old,m_max,            &
              &           minc_old,minc,lm_max_old,lm_max
         if ( n_r_max_old /= n_r_max )                                &
              &        write(*,'('' ! Old/New n_r_max='',2i4)')       &
             &              n_r_max_old,n_r_max

      end if

      do lm=1,lm_max
         l=lm2l(lm)
         m=lm2m(lm)
         lm2lmo(lm)=-1 ! -1 means that there is no data in the startfile
         lmo=0
         do mo=0,l_max_old,minc_old
            do lo=mo,l_max_old
               lmo=lmo+1
               if ( lo==l .and. mo==m ) then
                  lm2lmo(lm)=lmo ! data found in startfile
                  cycle
               end if
            end do
         end do
      end do

   end subroutine getLm2lmO
!------------------------------------------------------------------------------
   subroutine mapDataHydro(wo,zo,po,so,n_data_oldP,lm2lmo, &
                          n_r_max_old,lm_max_old,n_r_maxL, &
                          lbc1,lbc2,lbc3,lbc4,w,z,p,s )

      !--- Input variables
      integer,         intent(in) :: n_r_max_old,lm_max_old
      integer,         intent(in) :: n_r_maxL,n_data_oldP
      logical,         intent(in) :: lbc1,lbc2,lbc3,lbc4
      integer,         intent(in) :: lm2lmo(lm_max)
      complex(cp), intent(in) :: wo(n_data_oldP),zo(n_data_oldP)
      complex(cp), intent(in) :: po(n_data_oldP),so(n_data_oldP)

      !--- Output variables
      complex(cp),intent(out) :: w(lm_max,n_r_max),z(lm_max,n_r_max)
      complex(cp),intent(out) :: p(lm_max,n_r_max),s(lm_max,n_r_max)

      !--- Local variables
      integer :: lm,lmo,n,nR,lmStart,lmStop,nLMB
      complex(cp),allocatable :: woR(:),zoR(:)
      complex(cp),allocatable :: poR(:),soR(:)

      !PRINT*,omp_get_thread_num(),": Before nLMB loop, nLMBs=",nLMBs
      allocate( woR(n_r_maxL),zoR(n_r_maxL) )
      allocate( poR(n_r_maxL),soR(n_r_maxL) )
      bytes_allocated = bytes_allocated + 4*n_r_maxL*SIZEOF_DEF_COMPLEX
      write(*,"(A,I12)") "maximal allocated bytes in mapData are ",bytes_allocated

      !PERFON('mD_map')
      do nLMB=1,nLMBs ! Blocking of loop over all (l,m)
         lmStart=lmStartB(nLMB)
         lmStop =lmStopB(nLMB)

         !PRINT*,nLMB,lmStart,lmStop
         do lm=lmStart,lmStop
            lmo=lm2lmo(lm)
            if ( lmo > 0 ) then
               if ( n_r_max /= n_r_max_old ) then
                  do nR=1,n_r_max_old  ! copy on help arrays
                     n=lmo+(nR-1)*lm_max_old
                     woR(nR)=wo(n)
                     zoR(nR)=zo(n)
                     poR(nR)=po(n)
                     if ( l_heat ) soR(nR)=so(n)
                  end do
                  call mapDataR(woR,n_r_max,n_r_max_old,n_r_maxL,lBc1,.FALSE.)
                  call mapDataR(zoR,n_r_max,n_r_max_old,n_r_maxL,lBc2,.FALSE.)
                  call mapDataR(poR,n_r_max,n_r_max_old,n_r_maxL,lBc3,.FALSE.)
                  if ( l_heat ) call mapDataR(soR,n_r_max,n_r_max_old, & 
                                              n_r_maxL,lBc4,.FALSE.)
                  do nR=1,n_r_max
                     if ( lm > 1 ) then
                        w(lm,nR)=scale_v*woR(nR)
                        z(lm,nR)=scale_v*zoR(nR)
                        p(lm,nR)=scale_v*poR(nR)
                     end if
                     if ( l_heat ) s(lm,nR)=scale_s*soR(nR)
                  end do
               else
                  do nR=1,n_r_max
                     n=lmo+(nR-1)*lm_max_old
                     if ( lm > 1 ) then
                        w(lm,nR)=scale_v*wo(n)
                        z(lm,nR)=scale_v*zo(n)
                        p(lm,nR)=scale_v*po(n)
                     end if
                     if ( l_heat ) s(lm,nR)=scale_s*so(n)
                  end do
               end if
            end if
         end do
      end do
      !PERFOFF
      !PRINT*,omp_get_thread_num(),": After nLMB loop"
      deallocate(woR,zoR,poR,soR)
      bytes_allocated = bytes_allocated - 4*n_r_maxL*SIZEOF_DEF_COMPLEX

   end subroutine mapDataHydro
!------------------------------------------------------------------------------
   subroutine mapDataMag( wo,zo,po,so,n_data_oldP,n_rad_tot,n_r_max_old, &
                          lm_max_old,n_r_maxL,lm2lmo,dim1,l_IC,          & 
                          w,z,p,s )

      !--- Input variables
      integer,     intent(in) :: n_rad_tot,n_r_max_old,lm_max_old
      integer,     intent(in) :: n_r_maxL,n_data_oldP,dim1
      integer,     intent(in) :: lm2lmo(lm_max)
      logical,     intent(in) :: l_IC
      complex(cp), intent(in) :: wo(n_data_oldP),zo(n_data_oldP)
      complex(cp), intent(in) :: po(n_data_oldP),so(n_data_oldP)

      !--- Output variables
      complex(cp), intent(out) :: w(lm_maxMag,dim1),z(lm_maxMag,dim1)
      complex(cp), intent(out) :: p(lm_maxMag,dim1),s(lm_maxMag,dim1)

      !--- Local variables
      integer :: lm,lmo,n,nR,lmStart,lmStop,nLMB
      complex(cp), allocatable :: woR(:),zoR(:),poR(:),soR(:)

      !PRINT*,omp_get_thread_num(),": Before nLMB loop, nLMBs=",nLMBs
      allocate( woR(n_r_maxL),zoR(n_r_maxL) )
      allocate( poR(n_r_maxL),soR(n_r_maxL) )
      bytes_allocated = bytes_allocated + 4*n_r_maxL*SIZEOF_DEF_COMPLEX
      write(*,"(A,I12)") "maximal allocated bytes in mapData are ",bytes_allocated

      !PERFON('mD_map')
      do nLMB=1,nLMBs ! Blocking of loop over all (l,m)
         lmStart=lmStartB(nLMB)
         lmStop =lmStopB(nLMB)
         lmStart=max(2,lmStart)
         do lm=lmStart,lmStop
            lmo=lm2lmo(lm)
            if ( lmo > 0 ) then
               if ( n_rad_tot /= n_r_max_old ) then
                  do nR=1,n_r_max_old  ! copy on help arrays
                     n=lmo+(nR-1)*lm_max_old
                     woR(nR)=wo(n)
                     zoR(nR)=zo(n)
                     poR(nR)=po(n)
                     soR(nR)=so(n)
                  end do
                  call mapDataR(woR,dim1,n_r_max_old,n_r_maxL,.FALSE.,l_IC)
                  call mapDataR(zoR,dim1,n_r_max_old,n_r_maxL,.TRUE.,l_IC)
                  call mapDataR(poR,dim1,n_r_max_old,n_r_maxL,.TRUE.,l_IC)
                  call mapDataR(soR,dim1,n_r_max_old,n_r_maxL,.FALSE.,l_IC)
                  do nR=1,n_rad_tot
                     w(lm,nR)=scale_b*woR(nR)
                     z(lm,nR)=scale_b*zoR(nR)
                     p(lm,nR)=scale_b*poR(nR)
                     s(lm,nR)=scale_b*soR(nR)
                  end do
               else
                  do nR=1,n_rad_tot
                     n=lmo+(nR-1)*lm_max_old
                     w(lm,nR)=scale_b*wo(n)
                     z(lm,nR)=scale_b*zo(n)
                     p(lm,nR)=scale_b*po(n)
                     s(lm,nR)=scale_b*so(n)
                  end do
               end if
            end if
         end do
      end do
      !PERFOFF
      !PRINT*,omp_get_thread_num(),": After nLMB loop"
      deallocate(woR,zoR,poR,soR)
      bytes_allocated = bytes_allocated - 4*n_r_maxL*SIZEOF_DEF_COMPLEX

   end subroutine mapDataMag
!------------------------------------------------------------------------------
   subroutine mapDataR(dataR,n_rad_tot,n_r_max_old,n_r_maxL,lBc,l_IC)
      !
      !
      !  Copy (interpolate) data (read from disc file) from old grid structure
      !  to new grid. Linear interploation is used in r if the radial grid
      !  structure differs
      !
      !  called in mapdata
      !
      !

      !--- Input variables
      integer, intent(in) :: n_r_max_old
      integer, intent(in) :: n_r_maxL,n_rad_tot
      logical, intent(in) :: lBc,l_IC

      !--- Output variables
      complex(cp), intent(inout) :: dataR(:)  ! old data

      !-- Local variables
      integer :: nR, n_r_index_start
      type(costf_odd_t) :: chebt_oc_old
      complex(cp) :: work(n_r_maxL)
      real(cp) :: cheb_norm_old,scale

      !----- Initialize transform to cheb space:
      call chebt_oc_old%initialize(n_r_max_old, 2*n_r_maxL+2,2*n_r_maxL+5)

      !-- Guess the boundary values, since they have not been stored:
      if ( .not. l_IC .and. lBc ) then
         dataR(1)=2.0_cp*dataR(2)-dataR(3)
         dataR(n_r_max_old)=2.0_cp*dataR(n_r_max_old-1)-dataR(n_r_max_old-2)
      end if

      !----- Transform old data to cheb space:
      call chebt_oc_old%costf1(dataR,work)

      !----- Fill up cheb polynomial with zeros:
      if ( n_rad_tot>n_r_max_old ) then
         if ( l_IC) then
            n_r_index_start=n_r_max_old
         else 
            n_r_index_start=n_r_max_old+1
         end if
         do nR=n_r_index_start,n_rad_tot
            dataR(nR)=0.0_cp
         end do
      end if
    
      !----- Now transform to new radial grid points:
      if ( l_IC ) then
         call chebt_ic%costf1(dataR,work)
         !----- Rescale :
         cheb_norm_old=sqrt(2.0_cp/real(n_r_max_old-1,kind=cp))
         scale=cheb_norm_old/cheb_norm_ic
      else
         call chebt_oc%costf1(dataR,work)
         !----- Rescale :
         cheb_norm_old=sqrt(2.0_cp/real(n_r_max_old-1,kind=cp))
         scale=cheb_norm_old/cheb_norm
      end if
      do nR=1,n_rad_tot
         dataR(nR)=scale*dataR(nR)
      end do

      call chebt_oc_old%finalize()

   end subroutine mapDataR
!---------------------------------------------------------------------
   function mapDataR_copy(data_old,n_rad_new,n_rad_old,lBc,l_IC) result(data_new)
      ! A version of mapDataR() that does not write the result in-place of the
      ! input array

      !--- Input variables
      integer, intent(in) :: n_rad_old
      integer, intent(in) :: n_rad_new
      logical, intent(in) :: lBc,l_IC
      complex(cp), intent(in) :: data_old(:)

      !--- Result
      complex(cp) :: data_new(1:n_rad_new)

      ! intermediate
      complex(cp) :: dataR(1:max(n_rad_old, n_rad_new))

      if (.not. ubound(data_old, dim=1) == n_rad_old) then
        stop "Wrong bounds for data_old in function mapDataR"
      endif

      dataR(:) = 0.0_cp
      dataR(1:n_rad_old) = data_old(1:n_rad_old)
      call mapDataR(dataR,n_rad_new,n_rad_old,max(n_rad_new, n_rad_old),lBc,l_IC)
      data_new(1:n_rad_new) = dataR(1:n_rad_new)

    end function mapDataR_copy

end module
