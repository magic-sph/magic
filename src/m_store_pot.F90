!$Id$
module store_pot_mod

   use truncation
   use radial_functions
   use physical_parameters
   use horizontal_data
   use logic
   use output_data
   use LMLoop_data, only: llm, ulm, llm_real, ulm_real
   use parallel_mod, only: rank
   use communications, only: gather_from_lo_to_rank0
    
   implicit none

   private

   public :: storePot, storePotW

contains

   subroutine storePot(time,b,aj,b_ic,aj_ic,nPotSets,root,omega_ma,omega_ic)

      !-- Input of variables:
      real(kind=8),     intent(in) :: time
      complex(kind=8),  intent(in) :: b(lm_max,n_r_max)
      complex(kind=8),  intent(in) :: aj(lm_max,n_r_max)
      complex(kind=8),  intent(in) :: b_ic(lm_max,n_r_ic_max)
      complex(kind=8),  intent(in) :: aj_ic(lm_max,n_r_ic_max)
      character(len=*), intent(in) :: root
      real(kind=8),     intent(in) :: omega_ma,omega_ic

      integer,         intent(inout) :: nPotSets
    
      !-- Work arrays:
      complex(kind=8) :: workA(lm_max,n_r_max)
      complex(kind=8) :: workB(lm_max,n_r_max)
      complex(kind=8) :: workC(lm_max,n_r_max)
    
      character(80) :: string
      integer :: n_r,lm,n_cheb
      character(80) :: fileName
      logical :: lVB
      real(kind=8) :: chebNorm
    
      nPotSets=nPotSets+1
      lVB=.false.
      if ( root(1:1) /= 'T' ) lVB= .true. 
    
      !--- Copy:
      do n_r=1,n_r_max
         do lm =1,lm_max
            workA(lm,n_r)= b(lm,n_r)
            if ( lVB ) workB(lm,n_r)=aj(lm,n_r)
         end do
      end do
    
      !--- Transform to Cheb-space:
      call costf1(workA,lm_max_real,1,lm_max_real, &
                  workC,i_costf_init,d_costf_init)
      if ( lVB ) &
           call costf1(workB,lm_max_real,1,lm_max_real, &
                       workC,i_costf_init,d_costf_init)
    
      !--- Correct amplitude:
      chebNorm=dsqrt(2.D0/(n_r_max-1))
      do n_cheb=1,n_cheb_max
         do lm=1,lm_max
            if ( n_cheb == 1 .or. n_cheb == n_r_max ) then
               workA(lm,n_cheb)=chebNorm*0.5D0*workA(lm,n_cheb)
               if ( lVB ) workB(lm,n_cheb)=chebNorm*0.5D0*workB(lm,n_cheb)
            else
               workA(lm,n_cheb)=chebNorm*workA(lm,n_cheb)
               if ( lVB ) workB(lm,n_cheb)=chebNorm*workB(lm,n_cheb)
            end if
         end do
      end do
    
      !--- Write:
      if ( nPotSets == 0 ) then ! nPotSets=-1 on call
         fileName=trim(adjustl(root))//'.'//tag
      else
         !------ Names including the time:
         !           if ( l_graph_time ) then
         !              call dble2string(time,'_',6,string,length)
         !              fileName=root(1:lengthR)//'_t='//string(1:length)//'.'//tag
         !           else
         !------ Numbered names:
         write(string, *) nPotSets
         fileName=trim(adjustl(root))//'_'//trim(adjustl(string))//'.'//tag
         !           end if
      end if
    
      open(99, file=fileName, form='unformatted', status='unknown')
    
      write(99) l_max,n_cheb_max,n_cheb_ic_max,minc,lm_max
      write(99) sngl(ra),sngl(ek),sngl(pr),sngl(prmag), &
                      sngl(radratio),sngl(sigma_ratio), &
                         sngl(omega_ma),sngl(omega_ic)
      write(99) sngl(time), ( (cmplx(   real(workA(lm,n_cheb)), &
                                       aimag(workA(lm,n_cheb)), &
                                       kind=kind(0.D0) ),       &
                              lm=1,lm_max),n_cheb=1,n_cheb_max )
      if ( lVB ) then
         write(99) sngl(time), ( (cmplx(   real(workB(lm,n_cheb)), &
                                          aimag(workB(lm,n_cheb)), &
                                          kind=kind(0.D0) ),       &
                                 lm=1,lm_max),n_cheb=1,n_cheb_max )
      end if
    
      !-- Now inner core field
      if ( root(1:1) == 'B' .and. l_cond_ic ) then
    
         write(*,*) 'WRITING IC DATA INTO FILE:',fileName
    
         do n_r=1,n_r_ic_max
            do lm =1,lm_max
               workA(lm,n_r)= b_ic(lm,n_r)
               workB(lm,n_r)=aj_ic(lm,n_r)
            end do
         end do
    
         call costf1(workA,lm_max_real,1,lm_max_real, &
                     workC,i_costf1_ic_init,d_costf1_ic_init)
         call costf1(workB,lm_max_real,1,lm_max_real, &
                     workC,i_costf1_ic_init,d_costf1_ic_init)
    
         chebNorm=dsqrt(2.D0/(n_r_ic_max-1))
         do n_cheb=1,n_cheb_ic_max
            do lm=1,lm_max
               if ( n_cheb == 1 .or. n_cheb == n_r_ic_max ) then
                  workA(lm,n_cheb)=chebNorm*0.5D0*workA(lm,n_cheb)
                  workB(lm,n_cheb)=chebNorm*0.5D0*workB(lm,n_cheb)
               else
                  workA(lm,n_cheb)=chebNorm*workA(lm,n_cheb)
                  workB(lm,n_cheb)=chebNorm*workB(lm,n_cheb)
               end if
            end do
         end do
         write(99) sngl(time), ( (cmplx(   real(workA(lm,n_cheb)), &
                                          aimag(workA(lm,n_cheb)), &
                                          kind=kind(0d0)),         &
                                lm=1,lm_max),n_cheb=1,n_cheb_ic_max )
         write(99) sngl(time), ( (cmplx(   real(workB(lm,n_cheb)), &
                                          aimag(workB(lm,n_cheb)), &
                                          kind=kind(0d0)),         &
                                lm=1,lm_max),n_cheb=1,n_cheb_ic_max )
    
      end if
    
      close(99)

   end subroutine storePot
!----------------------------------------------------------------------
   subroutine storePotW(time,b,aj,b_ic,aj_ic,workA,workB,workC, &
                        nPotSets,root,omega_ma,omega_ic)

      !-- Input of variables:
      real(kind=8),     intent(in) :: time
      complex(kind=8),  intent(in) :: b(llm:ulm,n_r_max)
      complex(kind=8),  intent(in) :: aj(llm:ulm,n_r_max)
      complex(kind=8),  intent(in) :: b_ic(llm:ulm,n_r_ic_max)
      complex(kind=8),  intent(in) :: aj_ic(llm:ulm,n_r_ic_max)
      character(len=9), intent(in) :: root
      real(kind=8),     intent(in) :: omega_ma,omega_ic

      integer,          intent(inout) :: nPotSets

      !-- Work arrays:
      complex(kind=8) :: workA(llm:ulm,n_r_max)
      complex(kind=8) :: workB(llm:ulm,n_r_max)
      complex(kind=8) :: workC(llm:ulm,n_r_max)

      !-- Local variables
      complex(kind=8), allocatable :: workA_global(:,:),workB_global(:,:)
      character(len=80) :: string
      integer :: n_r,lm,n_cheb
      character(len=80) :: fileName
      logical :: lVB
      real(kind=8) :: chebNorm
       
      nPotSets=nPotSets+1
      lVB=.false.
      if ( root(1:1) /= 'T' ) lVB= .true. 

      !--- Copy:
      do n_r=1,n_r_max
         do lm =llm,ulm
            workA(lm,n_r)= b(lm,n_r)
            if ( lVB ) workB(lm,n_r)=aj(lm,n_r)
         end do
      end do

      !--- Transform to Cheb-space:
      call costf1(workA,ulm_real-llm_real+1,1,ulm_real-llm_real+1, &
                  workC,i_costf_init,d_costf_init)
      if ( lVB ) call costf1(workB,ulm_real-llm_real+1,1,ulm_real-llm_real+1, &
           &                 workC,i_costf_init,d_costf_init)

      !--- Correct amplitude:
      chebNorm=dsqrt(2.D0/(n_r_max-1))
      do n_cheb=1,n_cheb_max
         do lm=llm,ulm
            if ( n_cheb == 1 .or. n_cheb == n_r_max ) then
               workA(lm,n_cheb)=chebNorm*0.5D0*workA(lm,n_cheb)
               if ( lVB) workB(lm,n_cheb)=chebNorm*0.5D0*workB(lm,n_cheb)
            else
               workA(lm,n_cheb)=chebNorm*workA(lm,n_cheb)
               if ( lVB) workB(lm,n_cheb)=chebNorm*workB(lm,n_cheb)
            end if
         end do
      end do

      ! now gather the fields on rank 0 and write them to file
      ! it would be nicer to write the fields with MPI IO in parallel
      ! but then presumably the file format will change
      if ( rank == 0 ) then
         allocate(workA_global(1:lm_max,1:n_cheb_max))
         allocate(workB_global(1:lm_max,1:n_cheb_max))
#ifdef WITH_DEBUG
      else
         allocate(workA_global(1,1:n_cheb_max))
         allocate(workB_global(1,1:n_cheb_max))
#endif
      end if

      do n_cheb=1,n_cheb_max
         call gather_from_lo_to_rank0(workA(llm,n_cheb),workA_global(1,n_cheb))
         call gather_from_lo_to_rank0(workB(llm,n_cheb),workB_global(1,n_cheb))
      end do


      if ( rank == 0 ) then
         !--- Write:
         if ( nPotSets == 0 ) then ! nPotSets=-1 on call
            fileName=trim(adjustl(root))//'.'//tag
         else
            !------ Names including the time:
            !           if ( l_graph_time ) then
            !              call dble2string(time,'_',6,string,length)
            !              fileName=root(1:lengthR)//'_t='//string(1:length)//'.'//tag
            !           else
            !------ Numbered names:
            write(string, *) nPotSets
            fileName=trim(adjustl(root))//'_'//trim(adjustl(string))//'.'//tag
            !         end if
         end if

         open(99, file=fileName, form='unformatted', status='unknown')

         write(99) l_max,n_cheb_max,n_cheb_ic_max,minc,lm_max
         write(99) sngl(ra),sngl(ek),sngl(pr),sngl(prmag), &
                   sngl(radratio),sngl(sigma_ratio),       &
                   sngl(omega_ma),sngl(omega_ic)
         write(99) sngl(time), ( (cmplx(  real(workA_global(lm,n_cheb)), &
                                         aimag(workA_global(lm,n_cheb)), &
                                         kind=kind(0d0) ),               &
                                 lm=1,lm_max),n_cheb=1,n_cheb_max )
         if ( lVB ) then
            write(99) sngl(time), ( (cmplx(  real(workB_global(lm,n_cheb)), &
                                            aimag(workB_global(lm,n_cheb)), &
                                            kind=kind(0d0) ),               &
                                   lm=1,lm_max),n_cheb=1,n_cheb_max)
         end if
      end if
      
      !-- Now inner core field
      if ( root(1:1) == 'B' .and. l_cond_ic ) then

         do n_r=1,n_r_ic_max
            do lm =llm,ulm
               workA(lm,n_r)= b_ic(lm,n_r)
               workB(lm,n_r)=aj_ic(lm,n_r)
            end do
         end do

         call costf1(workA,ulm_real-llm_real+1,1,ulm_real-llm_real+1, &
                     workC,i_costf1_ic_init,d_costf1_ic_init)
         call costf1(workB,ulm_real-llm_real+1,1,ulm_real-llm_real+1, &
                     workC,i_costf1_ic_init,d_costf1_ic_init)

         chebNorm=dsqrt(2.D0/(n_r_ic_max-1))
         do n_cheb=1,n_cheb_ic_max
            do lm=llm,ulm
               if ( n_cheb == 1 .or. n_cheb == n_r_ic_max ) then
                  workA(lm,n_cheb)=chebNorm*0.5D0*workA(lm,n_cheb)
                  workB(lm,n_cheb)=chebNorm*0.5D0*workB(lm,n_cheb)
               else
                  workA(lm,n_cheb)=chebNorm*workA(lm,n_cheb)
                  workB(lm,n_cheb)=chebNorm*workB(lm,n_cheb)
               end if
            end do
         end do

         do n_cheb=1,n_cheb_ic_max
            call gather_from_lo_to_rank0(workA(llm,n_cheb),workA_global(1,n_cheb))
            call gather_from_lo_to_rank0(workB(llm,n_cheb),workB_global(1,n_cheb))
         end do

         if ( rank == 0 ) then
            write(*,*) 'WRITING IC DATA INTO FILE:',fileName

            write(99) sngl(time), ( (cmplx(   real(workA_global(lm,n_cheb)), &
                                             aimag(workA_global(lm,n_cheb)), &
                                             kind=kind(0d0) ),               &
                                    lm=1,lm_max),n_cheb=1,n_cheb_ic_max )
            write(99) sngl(time), ( (cmplx(   real(workB_global(lm,n_cheb)), &
                                             aimag(workB_global(lm,n_cheb)), &
                                             kind=kind(0d0)),                &
                                    lm=1,lm_max),n_cheb=1,n_cheb_ic_max )
         end if

      end if

      close(99)

      if ( rank == 0 ) deallocate( workA_global,workB_global ) 

   end subroutine storePotW
!----------------------------------------------------------------------
end module store_pot_mod
