module store_pot_mod
   !
   ! This module contains the subroutines that can be used to write unformatted
   ! fortran files that contain the flow/magnetic field potentials (in both
   ! Chebyshev and spectral space)
   !

   use precision_mod
   use truncation, only: n_r_max, n_r_ic_max, lm_max, &
       &                 n_cheb_ic_max, minc, l_max
   use radial_functions, only: chebt_ic, rscheme_oc, cheb_norm_ic
   use physical_parameters, only: ra, ek, pr, prmag, radratio, sigma_ratio
   use logic, only: l_cond_ic
   use output_data, only: tag
   use LMLoop_data, only: llm, ulm
   use parallel_mod, only: rank
   use communications, only: gather_from_lo_to_rank0
   use constants, only: two, half
    
   implicit none

   private

   integer :: fileHandle

   public :: storePot, storePotW

contains

   subroutine storePot(time,b,aj,b_ic,aj_ic,nPotSets,root,omega_ma,omega_ic)

      !-- Input of variables:
      real(cp),         intent(in) :: time
      complex(cp),      intent(in) :: b(llm:ulm,n_r_max)
      complex(cp),      intent(in) :: aj(llm:ulm,n_r_max)
      complex(cp),      intent(in) :: b_ic(llm:ulm,n_r_ic_max)
      complex(cp),      intent(in) :: aj_ic(llm:ulm,n_r_ic_max)
      character(len=*), intent(in) :: root
      real(cp),         intent(in) :: omega_ma,omega_ic

      integer,          intent(inout) :: nPotSets
    
      !-- Work arrays:
      complex(cp) :: workA(llm:ulm,n_r_max)
      complex(cp) :: workB(llm:ulm,n_r_max)
      complex(cp) :: workC(llm:ulm,n_r_max)
      complex(cp), allocatable :: workA_global(:,:)
      complex(cp), allocatable :: workB_global(:,:)
    
      character(80) :: string
      character(:), allocatable :: head
      integer :: n_r,lm,n_cheb, n_r_out
      character(80) :: fileName
      logical :: lVB
    
      head = trim(adjustl(root))
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
      call rscheme_oc%costf1(workA,ulm-llm+1,1,ulm-llm+1)
      if ( lVB ) call rscheme_oc%costf1(workB,ulm-llm+1,1,ulm-llm+1)
    
      !--- Correct amplitude:
      do n_r_out=1,rscheme_oc%n_max
         do lm=llm,ulm
            if ( n_r_out == 1 .or. n_r_out == n_r_max ) then
               workA(lm,n_r_out)=rscheme_oc%rnorm*rscheme_oc%boundary_fac*&
               &                 workA(lm,n_r_out)
               if ( lVB ) workB(lm,n_r_out)=rscheme_oc%rnorm* &
               &                rscheme_oc%boundary_fac*workB(lm,n_r_out)
            else
               workA(lm,n_r_out)=rscheme_oc%rnorm*workA(lm,n_r_out)
               if ( lVB ) workB(lm,n_r_out)=rscheme_oc%rnorm*workB(lm,n_r_out)
            end if
         end do
      end do

      ! now gather the fields on rank 0 and write them to file
      ! it would be nicer to write the fields with MPI IO in parallel
      ! but then presumably the file format will change
      if ( rank == 0 ) then
         allocate(workA_global(1:lm_max,1:rscheme_oc%n_max))
         allocate(workB_global(1:lm_max,1:rscheme_oc%n_max))
      else
         allocate(workA_global(1,1:rscheme_oc%n_max))
         allocate(workB_global(1,1:rscheme_oc%n_max))
      end if

      do n_r_out=1,rscheme_oc%n_max
         call gather_from_lo_to_rank0(workA(llm,n_r_out),workA_global(1:,n_r_out))
         call gather_from_lo_to_rank0(workB(llm,n_r_out),workB_global(1:,n_r_out))
      end do

      if ( rank == 0 ) then
         !--- Write:
         if ( nPotSets == 0 ) then ! nPotSets=-1 on call
            fileName=head//tag
         else
            !------ Names including the time:
            !           if ( l_graph_time ) then
            !              call dble2string(time,'_',6,string,length)
            !              fileName=root(1:lengthR)//'_t='//string(1:length)//'.'//tag
            !           else
            !------ Numbered names:
            write(string, *) nPotSets
            fileName=head(1:len(head)-1)//'_'//trim(adjustl(string))//'.'//tag
            !           end if
         end if
       
         open(newunit=fileHandle, file=fileName, form='unformatted', &
         &    status='unknown')
       
         write(fileHandle) l_max,rscheme_oc%n_max,n_cheb_ic_max,minc,lm_max
         write(fileHandle) real(ra,kind=outp), real(ek,kind=outp),               &
         &                 real(pr,kind=outp), real(prmag,kind=outp),            &
         &                 real(radratio,kind=outp), real(sigma_ratio,kind=outp),&
         &                 real(omega_ma,kind=outp), real(omega_ic,kind=outp)
         write(fileHandle) real(time,kind=outp), ( (cmplx(real(workA(lm,n_cheb)),&
         &                 aimag(workA(lm,n_cheb)),kind=outp ),lm=1,lm_max),     &
         &                 n_cheb=1,rscheme_oc%n_max )
         if ( lVB ) then
            write(fileHandle) real(time,kind=outp),((cmplx(real(workB(lm,n_cheb)),&
            &                 aimag(workB(lm,n_cheb)),kind=outp),lm=1,lm_max),    &
            &                 n_cheb=1,rscheme_oc%n_max )
         end if
      end if
    
      !-- Now inner core field
      if ( root(1:1) == 'B' .and. l_cond_ic ) then
    
         write(*,*) 'WRITING IC DATA INTO FILE:',fileName
    
         do n_r=1,n_r_ic_max
            do lm =llm,ulm
               workA(lm,n_r)= b_ic(lm,n_r)
               workB(lm,n_r)=aj_ic(lm,n_r)
            end do
         end do
    
         call chebt_ic%costf1(workA,ulm-llm+1,1,ulm-llm+1,workC)
         call chebt_ic%costf1(workB,ulm-llm+1,1,ulm-llm+1,workC)

         do n_cheb=1,n_cheb_ic_max
            do lm=llm,ulm
               if ( n_cheb == 1 .or. n_cheb == n_r_ic_max ) then
                  workA(lm,n_cheb)=cheb_norm_ic*half*workA(lm,n_cheb)
                  workB(lm,n_cheb)=cheb_norm_ic*half*workB(lm,n_cheb)
               else
                  workA(lm,n_cheb)=cheb_norm_ic*workA(lm,n_cheb)
                  workB(lm,n_cheb)=cheb_norm_ic*workB(lm,n_cheb)
               end if
            end do
         end do

         do n_cheb=1,n_cheb_ic_max
            call gather_from_lo_to_rank0(workA(llm,n_cheb),workA_global(1,n_cheb))
            call gather_from_lo_to_rank0(workB(llm,n_cheb),workB_global(1,n_cheb))
         end do

         if ( rank == 0 ) then

            write(fileHandle) real(time,kind=outp),((cmplx(real(workA(lm,n_cheb)),&
            &                 aimag(workA(lm,n_cheb)),kind=outp),lm=1,lm_max),    &
            &                 n_cheb=1,n_cheb_ic_max )
            write(fileHandle) real(time,kind=outp),((cmplx(real(workB(lm,n_cheb)),&
            &                 aimag(workB(lm,n_cheb)),kind=outp),lm=1,lm_max),    &
            &                 n_cheb=1,n_cheb_ic_max )

         end if
    
      end if

      deallocate( workA_global, workB_global )
    
      if ( rank == 0 ) then
         close(fileHandle)
      end if

   end subroutine storePot
!----------------------------------------------------------------------
   subroutine storePotW(time,b,aj,b_ic,aj_ic,workA,workB,workC, &
              &         nPotSets,root,omega_ma,omega_ic)
      !
      ! This routine stores the fields in spectral and Chebyshev space
      !

      !-- Input of variables:
      real(cp),         intent(in) :: time
      complex(cp),      intent(in) :: b(llm:ulm,n_r_max)
      complex(cp),      intent(in) :: aj(llm:ulm,n_r_max)
      complex(cp),      intent(in) :: b_ic(llm:ulm,n_r_ic_max)
      complex(cp),      intent(in) :: aj_ic(llm:ulm,n_r_ic_max)
      character(len=9), intent(in) :: root
      real(cp),         intent(in) :: omega_ma,omega_ic

      integer,          intent(inout) :: nPotSets

      !-- Work arrays:
      complex(cp) :: workA(llm:ulm,n_r_max)
      complex(cp) :: workB(llm:ulm,n_r_max)
      complex(cp) :: workC(llm:ulm,n_r_max)

      !-- Local variables
      complex(cp), allocatable :: workA_global(:,:),workB_global(:,:)
      character(len=80) :: string
      character(:), allocatable :: head
      integer :: n_r,lm,n_r_out,n_cheb
      character(len=80) :: fileName
      logical :: lVB
       
      head = trim(adjustl(root))
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
      call rscheme_oc%costf1(workA,ulm-llm+1,1,ulm-llm+1)
      if ( lVB ) call rscheme_oc%costf1(workB,ulm-llm+1,1,ulm-llm+1)

      !--- Correct amplitude:
      do n_r_out=1,rscheme_oc%n_max
         do lm=llm,ulm
            if ( n_r_out == 1 .or. n_r_out == n_r_max ) then
               workA(lm,n_r_out)=rscheme_oc%rnorm*rscheme_oc%boundary_fac* &
               &                 workA(lm,n_r_out)
               if ( lVB) workB(lm,n_r_out)=rscheme_oc%rnorm* &
               &                 rscheme_oc%boundary_fac*workB(lm,n_r_out)
            else
               workA(lm,n_r_out)=rscheme_oc%rnorm*workA(lm,n_r_out)
               if ( lVB) workB(lm,n_r_out)=rscheme_oc%rnorm*workB(lm,n_r_out)
            end if
         end do
      end do

      ! now gather the fields on rank 0 and write them to file
      ! it would be nicer to write the fields with MPI IO in parallel
      ! but then presumably the file format will change
      if ( rank == 0 ) then
         allocate(workA_global(1:lm_max,1:rscheme_oc%n_max))
         allocate(workB_global(1:lm_max,1:rscheme_oc%n_max))
      else
         allocate(workA_global(1,1:rscheme_oc%n_max))
         allocate(workB_global(1,1:rscheme_oc%n_max))
      end if

      do n_r_out=1,rscheme_oc%n_max
         call gather_from_lo_to_rank0(workA(llm,n_r_out),workA_global(1:,n_r_out))
         call gather_from_lo_to_rank0(workB(llm,n_r_out),workB_global(1:,n_r_out))
      end do

      if ( rank == 0 ) then
         !--- Write:
         if ( nPotSets == 0 ) then ! nPotSets=-1 on call
            fileName=head//tag
         else
            !------ Names including the time:
            !           if ( l_graph_time ) then
            !              call dble2string(time,'_',6,string,length)
            !              fileName=root(1:lengthR)//'_t='//string(1:length)//'.'//tag
            !           else
            !------ Numbered names:
            write(string, *) nPotSets
            fileName=head(1:len(head)-1)//'_'//trim(adjustl(string))//'.'//tag
            !         end if
         end if

         open(newunit=fileHandle, file=fileName, form='unformatted', &
         &    status='unknown')

         write(fileHandle) l_max,rscheme_oc%n_max,n_cheb_ic_max,minc,lm_max
         write(fileHandle) real(ra,kind=outp), real(ek,kind=outp),       &
         &                 real(pr,kind=outp), real(prmag,kind=outp),    &
         &                 real(radratio,kind=outp),                     &
         &                 real(sigma_ratio,kind=outp),                  &
         &                 real(omega_ma,kind=outp),                     &
         &    real(omega_ic,kind=outp)
         write(fileHandle) real(time,kind=outp),                         &
         &                 ((cmplx(real(workA_global(lm,n_r_out)),       &
         &                 aimag(workA_global(lm,n_r_out)),kind=outp ),  &
         &                 lm=1,lm_max),n_r_out=1,rscheme_oc%n_max )
         if ( lVB ) then
            write(fileHandle) real(time,kind=outp),                      &
            &                 ((cmplx(real(workB_global(lm,n_r_out)),     &
            &                 aimag(workB_global(lm,n_r_out)),kind=outp ),&
            &                 lm=1,lm_max),n_r_out=1,rscheme_oc%n_max)
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

         call chebt_ic%costf1(workA,ulm-llm+1,1,ulm-llm+1,workC)
         call chebt_ic%costf1(workB,ulm-llm+1,1,ulm-llm+1,workC)

         do n_cheb=1,n_cheb_ic_max
            do lm=llm,ulm
               if ( n_cheb == 1 .or. n_cheb == n_r_ic_max ) then
                  workA(lm,n_cheb)=cheb_norm_ic*half*workA(lm,n_cheb)
                  workB(lm,n_cheb)=cheb_norm_ic*half*workB(lm,n_cheb)
               else
                  workA(lm,n_cheb)=cheb_norm_ic*workA(lm,n_cheb)
                  workB(lm,n_cheb)=cheb_norm_ic*workB(lm,n_cheb)
               end if
            end do
         end do

         do n_cheb=1,n_cheb_ic_max
            call gather_from_lo_to_rank0(workA(llm,n_cheb),workA_global(1,n_cheb))
            call gather_from_lo_to_rank0(workB(llm,n_cheb),workB_global(1,n_cheb))
         end do

         if ( rank == 0 ) then
            write(*,*) 'WRITING IC DATA INTO FILE:',fileName

            write(fileHandle) real(time,kind=outp),                       &
                 &   ( (cmplx( real(workA_global(lm,n_cheb)),             &
                 &            aimag(workA_global(lm,n_cheb)), kind=outp ),&
                 &     lm=1,lm_max),n_cheb=1,n_cheb_ic_max )
            write(fileHandle) real(time,kind=outp),                       &
                 &   ( (cmplx( real(workB_global(lm,n_cheb)),             &
                 &            aimag(workB_global(lm,n_cheb)), kind=outp), &
                 &     lm=1,lm_max),n_cheb=1,n_cheb_ic_max )
         end if

      end if

      deallocate( workA_global,workB_global ) 

      if ( rank == 0 ) then
         close(fileHandle)
      end if

   end subroutine storePotW
!----------------------------------------------------------------------
end module store_pot_mod
