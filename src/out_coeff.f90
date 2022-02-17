module out_coeff
   !
   ! This module contains the subroutines that calculate the Bcmb files,
   ! the [B|V|T]_coeff_r files and the [B|V|T]_lmr files
   !

   use precision_mod
   use parallel_mod
   use mem_alloc, only: bytes_allocated
   use logic, only: l_r_field, l_cmb_field, l_save_out, l_average, &
       &            l_cond_ic, l_r_fieldT, l_r_fieldXi, l_mag
   use radial_functions, only: r, rho0
   use radial_data, only: nRstart, nRstop
   use physical_parameters, only: ra, ek, pr, prmag, radratio, sigma_ratio, &
       &                          raxi, sc
   use num_param, only: tScale
   use blocking, only: lm2, llm, ulm
   use truncation, only: lm_max, l_max, minc, n_r_max, n_r_ic_max, minc
   use communications, only: gather_from_lo_to_rank0, gather_all_from_lo_to_rank0,&
       &                     gt_IC, gt_OC
   use output_data, only: tag, n_coeff_r, n_r_array, n_r_step, l_max_r, n_coeff_r_max
   use constants, only: two, half

   implicit none

   private

   integer :: fileHandle
   integer, allocatable :: n_v_r_sets(:), n_b_r_sets(:)
   integer, allocatable :: n_T_r_sets(:), n_Xi_r_sets(:)
   integer, allocatable :: n_v_r_file(:)
   integer, allocatable :: n_t_r_file(:)
   integer, allocatable :: n_xi_r_file(:)
   integer, allocatable:: n_b_r_file(:)
   character(len=72), allocatable :: v_r_file(:)
   character(len=72), allocatable :: t_r_file(:)
   character(len=72), allocatable :: xi_r_file(:)
   character(len=72), allocatable :: b_r_file(:)


   public :: write_Bcmb, write_coeffs, write_Pot, initialize_coeff, finalize_coeff
#ifdef WITH_MPI
   public :: write_Pot_mpi
#endif

contains

   subroutine initialize_coeff

      integer :: n
      character(len=72) :: string

      !-- Coeffs at radial levels:
      allocate ( n_coeff_r(n_coeff_r_max))
      allocate ( n_v_r_file(n_coeff_r_max), v_r_file(n_coeff_r_max) )
      allocate ( n_v_r_sets(n_coeff_r_max) )
      n_v_r_sets=0

      if ( l_mag ) then
         allocate ( n_b_r_file(n_coeff_r_max), b_r_file(n_coeff_r_max) )
         allocate ( n_b_r_sets(n_coeff_r_max) )
         n_b_r_sets=0
      end if

      do n=1,n_coeff_r_max
         write(string,*) n
         v_r_file(n)='V_coeff_r'//trim(adjustl(string))//'.'//tag
         if ( l_mag ) then
            B_r_file(n)='B_coeff_r'//trim(adjustl(string))//'.'//tag
         end if
      end do

      if ( l_r_fieldT ) then
         allocate ( n_t_r_file(n_coeff_r_max), t_r_file(n_coeff_r_max) )
         allocate ( n_t_r_sets(n_coeff_r_max) )
         n_T_r_sets=0

         do n=1,n_coeff_r_max
            write(string,*) n
            t_r_file(n)='T_coeff_r'//trim(adjustl(string))//'.'//tag
         end do
      end if

      if ( l_r_fieldXi ) then
         allocate ( n_xi_r_file(n_coeff_r_max), xi_r_file(n_coeff_r_max) )
         allocate ( n_xi_r_sets(n_coeff_r_max) )
         n_Xi_r_sets=0

         do n=1,n_coeff_r_max
            write(string,*) n
            xi_r_file(n)='Xi_coeff_r'//trim(adjustl(string))//'.'//tag
         end do
      end if

      if ( count(n_r_array>0)> 0 ) then
         n_coeff_r=n_r_array(1:n_coeff_r_max)
      else
         n_r_step=max(n_r_step,1)
         do n=1,n_coeff_r_max
            n_coeff_r(n)=n*n_r_step  ! used every n_r_step point !
         end do
      end if

      if ( rank == 0 .and. ( .not. l_save_out ) ) then

         if ( l_r_field ) then
            do n=1,n_coeff_r_max
               open(newunit=n_v_r_file(n), file=v_r_file(n), &
               &    status='new', form='unformatted')
               if ( l_mag ) then
                  open(newunit=n_b_r_file(n), file=b_r_file(n), &
                  &    status='new', form='unformatted')
               end if
            end do
         end if
         if ( l_r_fieldT ) then
            do n=1,n_coeff_r_max
               open(newunit=n_t_r_file(n), file=t_r_file(n), &
               &    status='new', form='unformatted')
            end do
         end if
         if ( l_r_fieldXi ) then
            do n=1,n_coeff_r_max
               open(newunit=n_xi_r_file(n), file=xi_r_file(n), &
               &    status='new', form='unformatted')
            end do
         end if

      end if

   end subroutine initialize_coeff
!----------------------------------------------------------------------
   subroutine finalize_coeff

      integer :: n

      if ( rank == 0 .and. ( .not. l_save_out ) ) then
         if ( l_r_field ) then
            do n=1,n_coeff_r_max
               close(n_v_r_file(n))
               if ( l_mag ) then
                  close(n_b_r_file(n))
               end if
            end do
         end if

         if ( l_r_fieldT ) then
            do n=1,n_coeff_r_max
               close(n_t_r_file(n))
            end do
         end if

         if ( l_r_fieldXi ) then
            do n=1,n_coeff_r_max
               close(n_xi_r_file(n))
            end do
         end if
      end if

      deallocate ( n_coeff_r, n_v_r_file, v_r_file, n_v_r_sets )
      if ( l_mag ) deallocate ( n_b_r_file, b_r_file, n_b_r_sets )
      if ( l_r_fieldT ) deallocate ( n_t_r_file, t_r_file, n_t_r_sets )
      if ( l_r_fieldXi ) deallocate ( n_xi_r_file, xi_r_file, n_xi_r_sets )

   end subroutine finalize_coeff
!----------------------------------------------------------------------
   subroutine write_Bcmb(time,b_LMloc,l_max_cmb,n_cmb_sets,cmb_file,n_cmb_file)
      !
      ! Each call of this subroutine writes time and the poloidal magnetic
      ! potential coefficients b at the CMB up to degree and order l_max_cmb
      ! at the end of output file cmb_file. The parameters l_max_cmb, minc and
      ! the number of stored coeffs are written into the header of cmb_file.
      ! Each further set contains:
      !
      ! Real and imaginary part of b(*) for all orders m<=l are written
      ! for a specific degree l, then for the degrees l+1, l+2, l_max_cmb.
      !

      !-- Input variables:
      real(cp),         intent(in) :: time             ! Time
      complex(cp),      intent(in) :: b_LMloc(llm:ulm) ! Poloidal field potential
      character(len=*), intent(in) :: cmb_file         ! Name of output file

      !-- Output variables:
      integer, intent(inout) :: n_cmb_file    ! Output unit for $cmb_file
      integer, intent(inout) :: l_max_cmb     ! Max degree of output
      integer, intent(inout) :: n_cmb_sets    ! Total no. of cmb sets,
      ! should be set to zero before first call

      !-- Local variables:
      integer :: n,n_out        ! counters
      integer :: l,m            ! degree and order
      integer :: lm             ! position of (l,m) in b(*,n_r_cmb)
      integer :: m_max_cmb      ! Max order of output
      integer :: lm_max_cmb     ! Max no of combinations l,m for output
      integer :: n_data         ! No of output data

      complex(cp) :: work(lm_max)
      real(cp), allocatable ::  output(:) ! Output array

      call gather_from_lo_to_rank0(b_LMloc, work)

      if ( rank == 0 ) then

         !--- Definition of max degree for output
         if ( l_max < l_max_cmb ) l_max_cmb=l_max

         !--- Calculate no of data for l_max_cmb:
         m_max_cmb=(l_max_cmb/minc)*minc
         lm_max_cmb= m_max_cmb*(l_max_cmb+1)/minc        &
         &          -m_max_cmb*(m_max_cmb-minc)/(2*minc) &
         &          +l_max_cmb-m_max_cmb+1
         n_data=2*lm_max_cmb-l_max_cmb-2

         allocate(output(n_data))

         !--- Increase no. of cmb_sets:
         n_cmb_sets=n_cmb_sets+1

         !--- Open output file name:
         if ( l_save_out .or. n_cmb_sets == 0 ) then
            open(unit=n_cmb_file, file=cmb_file, position='append', &
            &    form='unformatted')
         end if

         !--- If this is the first set write l_max_cmb and minc into file:
         if ( n_cmb_sets <= 1 ) write(n_cmb_file) l_max_cmb,minc,n_data

         !--- Write b(*) into output array output(*):
         n_out=0

         !--- Axisymmetric part: (m=0) only real part stored
         do l=1,l_max_cmb
            lm=lm2(l,0)
            n_out=n_out+1
            output(n_out)=real(work(lm))
         end do

         !--- Non-axisymmetric part: store real and imag part
         do m=minc,l_max_cmb,minc
            do l=m,l_max_cmb
               lm=lm2(l,m)
               n_out=n_out+1
               output(n_out)=real(work(lm))
               n_out=n_out+1
               output(n_out)=aimag(work(lm))
            end do
         end do

         !--- Finally write output array output(*) into cmb_file:
         write(n_cmb_file) time,(output(n),n=1,n_out)

         !--- Close cmb_file
         if ( l_save_out .or. n_cmb_sets == 0 ) close(n_cmb_file)

         deallocate(output)

      end if

   end subroutine write_Bcmb
!----------------------------------------------------------------------
   subroutine write_coeffs(w_LMloc, dw_LMloc, ddw_LMloc, z_LMLoc, b_LMLoc,  &
              &            db_LMloc, ddb_LMloc, aj_LMloc, dj_LMloc, s_LMloc,&
              &            xi_LMloc, timeScaled)
      !
      ! This routine handles the writing of coefficients at a given depth
      !

      real(cp),    intent(in) :: timeScaled
      complex(cp), intent(in) :: w_LMloc(:,:)
      complex(cp), intent(in) :: dw_LMloc(:,:)
      complex(cp), intent(in) :: ddw_LMloc(:,:)
      complex(cp), intent(in) :: z_LMloc(:,:)
      complex(cp), intent(in) :: b_LMloc(:,:)
      complex(cp), intent(in) :: db_LMloc(:,:)
      complex(cp), intent(in) :: ddb_LMloc(:,:)
      complex(cp), intent(in) :: aj_LMloc(:,:)
      complex(cp), intent(in) :: dj_LMloc(:,:)
      complex(cp), intent(in) :: s_LMloc(:,:)
      complex(cp), intent(in) :: xi_LMloc(:,:)

      integer :: n, nR

      !--- Store potential coeffs for velocity fields and magnetic fields
      do n=1,n_coeff_r_max
         nR=n_coeff_r(n)
         call write_coeff_r(timeScaled,w_LMloc(:,nR),dw_LMloc(:,nR),  &
              &             ddw_LMloc(:,nR),z_LMloc(:,nR),r(nR),      &
              &             l_max_r,n_v_r_sets(n),v_r_file(n),        &
              &             n_v_r_file(n),1)
         if ( l_mag )                                                   &
            call write_coeff_r(timeScaled,b_LMloc(:,nR),db_LMloc(:,nR), &
                 &             ddb_LMloc(:,nR),aj_LMloc(:,nR),r(nR),    &
                 &             l_max_r,n_b_r_sets(n),b_r_file(n),       &
                 &             n_b_r_file(n),2)
         if ( l_r_fieldT )                                              &
            call write_coeff_r(timeScaled,s_LMloc(:,nR),db_LMloc(:,nR), &
                 &             ddb_LMloc(:,nR),aj_LMloc(:,nR),r(nR),    &
                 &             l_max_r,n_T_r_sets(n),T_r_file(n),       &
                 &             n_t_r_file(n),3)
         if ( l_r_fieldXi )                                              &
            call write_coeff_r(timeScaled,xi_LMloc(:,nR),db_LMloc(:,nR), &
                 &             ddb_LMloc(:,nR),aj_LMloc(:,nR),r(nR),     &
                 &             l_max_r,n_Xi_r_sets(n),Xi_r_file(n),      &
                 &             n_xi_r_file(n),3)
      end do

   end subroutine write_coeffs
!----------------------------------------------------------------------
   subroutine write_coeff_r(time,w_LMloc,dw_LMloc,ddw_LMloc,z_LMloc,r,  &
              &             l_max_r,n_sets,file,n_file,nVBS)
      !
      ! Each call of this subroutine writes time and the poloidal and
      ! toroidal coeffitients w,dw,z at a specific radius up to degree
      ! and order l_max_r at the end of output file $file.
      ! If the input is magnetic field (nVBS=2) ddw is stored as well.
      ! If the input is entropy field (nVBS=3) only ddw=s is stored.
      ! The parameters l_max_r, minc, the number of stored coeffs and
      ! radius in the outer core are written into the first line of $file.
      ! Each further set contains:
      !
      ! Real and imaginary part of w(*) for all orders m<=l are written
      ! for a specific degree l, then for the degrees l+1, l+2, l_max_r.
      !

      !-- Input variables:
      real(cp),         intent(in) :: r                 ! radius of coeffs
      real(cp),         intent(in) :: time              ! Time
      complex(cp),      intent(in) :: w_LMloc(llm:ulm)  ! Poloidal field potential
      complex(cp),      intent(in) :: dw_LMloc(llm:ulm) ! dr of Poloidal field potential
      complex(cp),      intent(in) :: ddw_LMloc(llm:ulm)! dr^2 of Poloidal field potential
      complex(cp),      intent(in) :: z_LMloc(llm:ulm)   ! Toroidal field potential
      character(len=*), intent(in) :: file         ! Name of output file
      integer,          intent(in) :: nVBS   ! 1 if flow, 2 if magnetic field and 3 if temperature or composition

      !-- Output:
      integer, intent(inout) :: n_file      ! Output unit for $file
      integer, intent(inout) :: l_max_r     ! Max degree of output
      integer, intent(inout) :: n_sets      ! Total no. of cmb sets,
      ! should be set to zero before first call

      !-- Local variables:
      integer :: n,n_out        ! counter
      integer :: l,m            ! degree and order
      integer :: lm             ! position of (l,m) in v(*),..
      integer :: m_max_r        ! Max order of output
      integer :: lm_max_r       ! Max no of combinations l,m for output
      integer :: n_data         ! No of output data

      complex(cp) :: work(lm_max)
      real(cp), allocatable ::  output(:)! Output array

      !--- Definition of max degree for output
      if ( l_max < l_max_r ) l_max_r=l_max

      !--- Calculate no of data for l_max_r:
      m_max_r=(l_max_r/minc)*minc
      lm_max_r=m_max_r*(l_max_r+1)/minc- m_max_r*(m_max_r-minc)/(2*minc) + &
      &        l_max_r-m_max_r+1
      n_data=2*lm_max_r-l_max_r-2

      if ( rank == 0 ) then
         if ( nVBS == 1 ) then
            allocate(output(3*n_data))
         else if ( nVBS == 2 ) then
            allocate(output(4*n_data))
         else if ( nVBS == 3 ) then
            allocate(output(n_data+1))
         end if
      end if

      !--- Increase no. of sets:
      n_sets=n_sets+1

      !--- Write b(*) into output array output(*):
      n_out=0

      call gather_from_lo_to_rank0(w_LMloc, work)

      if ( rank == 0 ) then
         if ( nVBS == 3 ) then
            !--- Axisymmetric part of s: (m=0) only real part stored
            do l=0,l_max ! start with l=0
               lm=lm2(l,0)
               if ( l <= l_max_r ) then
                  n_out=n_out+1
                  output(n_out)=real(work(lm))
               end if
            end do
         else
            !--- Axisymmetric part of w: (m=0) only real part stored
            do l=1,l_max ! start with l=1
               lm=lm2(l,0)
               if ( l <= l_max_r ) then
                  n_out=n_out+1
                  output(n_out)=real(work(lm))
               end if
            end do
         end if

         !--- Non-axisymmetric part of w: store real and imag part
         do m=minc,l_max_r,minc
            do l=m,l_max
               lm=lm2(l,m)
               if ( l <= l_max_r ) then
                  n_out=n_out+1
                  output(n_out)=real(work(lm))
                  n_out=n_out+1
                  output(n_out)=aimag(work(lm))
               end if
            end do
         end do
      end if

      if ( nVBS /= 3 ) then

         call gather_from_lo_to_rank0(dw_LMloc, work)

         if ( rank == 0 ) then
            !-- Now output for flow or magnetic field only:
            !--- Axisymmetric part of dw: (m=0) only real part stored
            do l=1,l_max
               lm=lm2(l,0)
               if ( l <= l_max_r ) then
                  n_out=n_out+1
                  output(n_out)=real(work(lm))
               end if
            end do
            !--- Non-axisymmetric part of dv: store real and imag part
            do m=minc,l_max_r,minc
               do l=m,l_max
                  lm=lm2(l,m)
                  if ( l <= l_max_r ) then
                     n_out=n_out+1
                     output(n_out)=real(work(lm))
                     n_out=n_out+1
                     output(n_out)=aimag(work(lm))
                  end if
               end do
            end do
         end if

         call gather_from_lo_to_rank0(z_LMloc, work)

         if ( rank == 0 ) then
            !--- Axisymmetric part of z: (m=0) only real part stored
            do l=1,l_max
               lm=lm2(l,0)
               if ( l <= l_max_r ) then
                  n_out=n_out+1
                  output(n_out)=real(work(lm))
               end if
            end do
            !--- Non-axisymmetric part of z: store real and imag part
            do m=minc,l_max_r,minc
               do l=m,l_max
                  lm=lm2(l,m)
                  if ( l <= l_max_r ) then
                     n_out=n_out+1
                     output(n_out)=real(work(lm))
                     n_out=n_out+1
                     output(n_out)=aimag(work(lm))
                  end if
               end do
            end do
         end if
      end if

      !--- If this is a magnetic field I also store the second radial derivative
      !    of the poloidal potential to caluclate diffusion:
      if ( nVBS == 2 ) then

         call gather_from_lo_to_rank0(ddw_LMloc, work)

         if ( rank == 0 ) then
            !--- Axisymmetric part of ddw: (m=0) only real part stored
            do l=1,l_max
               lm=lm2(l,0)
               if ( l <= l_max_r ) then
                  n_out=n_out+1
                  output(n_out)=real(work(lm))
               end if
            end do
            !--- Non-axisymmetric part of ddw: store real and imag part
            do m=minc,l_max_r,minc
               do l=m,l_max
                  lm=lm2(l,m)
                  if ( l <= l_max_r ) then
                     n_out=n_out+1
                     output(n_out)=real(work(lm))
                     n_out=n_out+1
                     output(n_out)=aimag(work(lm))
                  end if
               end do
            end do
         end if
      end if

      if ( rank == 0 ) then

         !--- Open output file with name $file:
         if ( l_save_out ) then
            open(unit=n_file, file=file, form='unformatted', status='unknown', &
            &    position='append')
         end if

         !--- If this is the first set write, l_max_r and minc into first line:
         if ( n_sets == 1 ) then
            write(n_file) l_max_r,minc,n_data,r
         end if

         !--- Finally write output array output(*) into file:
         write(n_file) time,(output(n),n=1,n_out)

         !--- Close file
         if ( l_save_out ) close(n_file)

         deallocate(output)

      end if

   end subroutine write_coeff_r
!-------------------------------------------------------------------------------
#ifdef WITH_MPI
   subroutine write_Pot_mpi(time,b,aj,b_ic,aj_ic,nPotSets,root,omega_ma,omega_ic)
      !
      ! This routine stores the fields in (lm,r) space using MPI-IO
      !

      !-- Input of variables:
      real(cp),         intent(in) :: time ! Time
      complex(cp),      intent(in) :: b(lm_max,nRstart:nRstop) ! Poloidal potential
      complex(cp),      intent(in) :: aj(lm_max,nRstart:nRstop)! Toroidal potential
      complex(cp),      intent(in) :: b_ic(llm:ulm,n_r_ic_max)
      complex(cp),      intent(in) :: aj_ic(llm:ulm,n_r_ic_max)
      character(len=*), intent(in) :: root
      real(cp),         intent(in) :: omega_ma,omega_ic
      integer,          intent(in) :: nPotSets

      !-- Local variables
      complex(outp), allocatable :: tmp(:,:)
      complex(cp), allocatable :: work(:,:)
      integer :: info, fh, version, istat(MPI_STATUS_SIZE), datatype
      integer :: arr_size(2), arr_loc_size(2), arr_start(2)
      integer(lip) :: disp, offset, size_tmp
      character(80) :: string
      character(:), allocatable :: head
      character(80) :: fileName
      logical :: lVB

      version = 1 ! file version


      allocate( tmp(lm_max,nRstart:nRstop) )

      head = trim(adjustl(root))
      lVB=.false.
      if ( root(1:1) /= 'T' .and. root(1:2) /= 'Xi' ) lVB= .true.

      if ( nPotSets == 0 ) then
         fileName=head//tag
      else
         write(string, *) nPotSets
         fileName=head(1:len(head)-1)//'_'//trim(adjustl(string))//'.'//tag
         !         end if
      end if

      !--  MPI-IO setup
      call mpiio_setup(info)

      !-- Open file
      call MPI_File_Open(MPI_COMM_WORLD, fileName, ior(MPI_MODE_WRONLY, &
           &             MPI_MODE_CREATE), info, fh, ierr)

      !-- Set the first view
      disp = 0
      call MPI_File_Set_View(fh, disp, MPI_BYTE, MPI_BYTE, "native", &
           &                 info, ierr)

      if ( rank == 0 ) then
         !-- Write the header of the file
         call MPI_File_Write(fh, version, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, real(time*tScale,outp), 1, MPI_OUT_REAL, &
              &              istat, ierr)
         call MPI_File_Write(fh, real(ra,outp), 1, MPI_OUT_REAL, istat, ierr)
         call MPI_File_Write(fh, real(pr,outp), 1, MPI_OUT_REAL, istat, ierr)
         call MPI_File_Write(fh, real(raxi,outp), 1, MPI_OUT_REAL, istat, ierr)
         call MPI_File_Write(fh, real(sc,outp), 1, MPI_OUT_REAL, istat, ierr)
         call MPI_File_Write(fh, real(prmag,outp), 1, MPI_OUT_REAL, istat, ierr)
         call MPI_File_Write(fh, real(ek,outp), 1, MPI_OUT_REAL, istat, ierr)
         call MPI_File_Write(fh, real(radratio,outp), 1, MPI_OUT_REAL, istat, &
              &              ierr)
         call MPI_File_Write(fh, real(sigma_ratio, outp), 1, MPI_OUT_REAL, &
              &              istat, ierr)
         call MPI_File_Write(fh, n_r_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, n_r_ic_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, l_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, minc, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, lm_max, 1, MPI_INTEGER, istat, ierr)

         call MPI_File_Write(fh, real(omega_ic,outp), 1, MPI_OUT_REAL, &
              &              istat, ierr)
         call MPI_File_Write(fh, real(omega_ma,outp), 1, MPI_OUT_REAL, &
              &              istat, ierr)

         call MPI_File_Write(fh, real(r,outp), n_r_max, MPI_OUT_REAL, &
              &              istat, ierr)
         call MPI_File_Write(fh, real(rho0,outp), n_r_max, MPI_OUT_REAL, &
              &              istat, ierr)

         !-- Rank 0 gets the displacement
         call MPI_File_get_position(fh, offset, ierr)
         call MPI_File_get_byte_offset(fh, offset, disp, ierr)
      end if

      !-- Broadcast the displacement
      call MPI_Bcast(disp, 1, MPI_OFFSET, 0, MPI_COMM_WORLD, ierr)

      arr_size(1) = lm_max
      arr_size(2) = n_r_max
      arr_loc_size(1) = lm_max
      arr_loc_size(2) = nR_per_rank
      arr_start(1) = 0
      arr_start(2) = nRstart-1
      call MPI_Type_Create_Subarray(2, arr_size, arr_loc_size, arr_start, &
           &                        MPI_ORDER_FORTRAN, MPI_COMPLEX8,      &
           &                        datatype, ierr)
      call MPI_Type_Commit(datatype, ierr)

      !-- Copy into a single precision array
      tmp(:,:) = cmplx(b(:,:), kind=outp)

      !-- Set the view after the header
      call MPI_File_Set_View(fh, disp, MPI_COMPLEX8, datatype, "native", &
           &                 info, ierr)

      size_tmp=int(lm_max,kind=lip)*int(n_r_max,kind=lip)* &
      &        int(2*SIZEOF_OUT_REAL,kind=lip)

      !-- Poloidal potential
      call MPI_File_Write_all(fh, tmp, lm_max*nR_per_rank, MPI_COMPLEX8, &
           &                  istat, ierr)
      disp = disp+size_tmp
      call MPI_File_Set_View(fh, disp, MPI_COMPLEX8, datatype, "native", &
           &                 info, ierr)

      !-- Toroidal potential
      if ( lVB ) then
         tmp(:,:) = cmplx(aj(:,:), kind=outp)
         call MPI_File_Write_all(fh, tmp, lm_max*nR_per_rank, MPI_COMPLEX8, &
              &                  istat, ierr)
         disp = disp+size_tmp
         call MPI_File_Set_View(fh, disp, MPI_COMPLEX8, datatype, "native", &
              &                 info, ierr)
      end if


      !-- Displacement at the end of the file
      offset = 0
      call MPI_File_Seek(fh, offset, MPI_SEEK_END, ierr)
      call MPI_File_get_byte_offset(fh, offset, disp, ierr)
      call MPI_File_Set_View(fh, disp, MPI_BYTE, MPI_BYTE, "native", &
           &                 info, ierr)

      call MPI_Type_Free(datatype, ierr)
      deallocate( tmp )

      !-- Now inner core field
      if ( root(1:1) == 'B' .and. l_cond_ic ) then

         if ( rank == 0 ) then
            allocate ( work(lm_max, n_r_ic_max), tmp(lm_max, n_r_ic_max) )
         else
            allocate ( work(1,1), tmp(1,1) )
         end if

         call gather_all_from_lo_to_rank0(gt_IC, b_ic, work)
         if ( rank == 0 ) then
            tmp(:,:)=cmplx(work(:,:), kind=outp)
            call MPI_File_Write(fh, tmp, lm_max*n_r_ic_max, MPI_COMPLEX8, &
                 &              istat, ierr)
         end if

         call gather_all_from_lo_to_rank0(gt_IC, aj_ic, work)
         if ( rank == 0 ) then
            tmp(:,:)=cmplx(work(:,:), kind=outp)
            call MPI_File_Write(fh, tmp, lm_max*n_r_ic_max, MPI_COMPLEX8, &
                 &              istat, ierr)
         end if

         deallocate( tmp, work )

      end if

      call MPI_Info_Free(info, ierr)
      call MPI_File_Close(fh, ierr)

   end subroutine write_Pot_mpi
!-------------------------------------------------------------------------------
#endif
   subroutine write_Pot(time,b,aj,b_ic,aj_ic,nPotSets,root,omega_ma,omega_ic)
      !
      ! This routine stores the fields in spectral and radial space
      !

      !-- Input of variables:
      real(cp),         intent(in) :: time ! Output time
      complex(cp),      intent(in) :: b(llm:ulm,n_r_max)
      complex(cp),      intent(in) :: aj(llm:ulm,n_r_max)
      complex(cp),      intent(in) :: b_ic(llm:ulm,n_r_ic_max)
      complex(cp),      intent(in) :: aj_ic(llm:ulm,n_r_ic_max)
      character(len=*), intent(in) :: root ! File prefix
      real(cp),         intent(in) :: omega_ma,omega_ic
      integer,          intent(in) :: nPotSets

      !-- Work arrays:
      complex(cp), allocatable :: workA_global(:,:)
      complex(cp), allocatable :: workB_global(:,:)

      !-- File outputs:
      character(80) :: string
      character(:), allocatable :: head
      integer :: n_r, lm, version
      character(80) :: fileName
      logical :: lVB

      version = 1

      head = trim(adjustl(root))
      lVB=.false.
      if ( root(1:1) /= 'T' .and. root(1:2) /= 'Xi' ) lVB= .true.


      ! now gather the fields on rank 0 and write them to file
      ! it would be nicer to write the fields with MPI IO in parallel
      ! but then presumably the file format will change
      if ( rank == 0 ) then
         allocate(workA_global(lm_max,n_r_max))
         allocate(workB_global(lm_max,n_r_max))
      else
         allocate(workA_global(1,n_r_max))
         allocate(workB_global(1,n_r_max))
      end if

      call gather_all_from_lo_to_rank0(gt_OC, b, workA_global)
      call gather_all_from_lo_to_rank0(gt_OC, aj, workB_global)

      if ( rank == 0 ) then
         !--- Write:
         if ( nPotSets == 0 ) then ! nPotSets=-1 on call
            fileName=head//tag
         else
            write(string, *) nPotSets
            fileName=head(1:len(head)-1)//'_'//trim(adjustl(string))//'.'//tag
            !         end if
         end if

         open(newunit=fileHandle, file=fileName, form='unformatted', &
         &    status='unknown', access='stream')

         write(fileHandle) version, real(time*tScale,kind=outp)
         write(fileHandle) real(ra,kind=outp), real(pr,kind=outp),     &
         &                 real(raxi,kind=outp), real(sc,kind=outp),   &
         &                 real(prmag,kind=outp), real(ek,kind=outp),  &
         &                 real(radratio,kind=outp),                   &
         &                 real(sigma_ratio,kind=outp)

         write(fileHandle) n_r_max,n_r_ic_max,l_max,minc,lm_max

         write(fileHandle) real(omega_ic,kind=outp), real(omega_ma,kind=outp)

         write(fileHandle) real(r,kind=outp), real(rho0, kind=outp)

         write(fileHandle) ((cmplx(real(workA_global(lm,n_r)),         &
         &                 aimag(workA_global(lm,n_r)),kind=outp ),    &
         &                 lm=1,lm_max),n_r=1,n_r_max)
         if ( lVB ) then
            write(fileHandle) ((cmplx(real(workB_global(lm,n_r)),      &
            &                 aimag(workB_global(lm,n_r)),kind=outp ), &
            &                 lm=1,lm_max),n_r=1,n_r_max)
         end if
      end if


      !-- Now inner core field
      if ( root(1:1) == 'B' .and. l_cond_ic ) then

         call gather_all_from_lo_to_rank0(gt_IC, b_ic, workA_global)
         call gather_all_from_lo_to_rank0(gt_IC, aj_ic, workB_global)

         if ( rank == 0 ) then
            write(fileHandle) ( (cmplx( real(workA_global(lm,n_r)),    &
            &                 aimag(workA_global(lm,n_r)), kind=outp ),&
            &          lm=1,lm_max),n_r=1,n_r_ic_max )
            write(fileHandle) ( (cmplx( real(workB_global(lm,n_r)),    &
            &                 aimag(workB_global(lm,n_r)), kind=outp), &
            &          lm=1,lm_max),n_r=1,n_r_ic_max )
         end if

      end if

      if ( rank == 0 ) then
         close(fileHandle)
      end if

      deallocate( workA_global, workB_global )

   end subroutine write_Pot
!------------------------------------------------------------------------------
end module out_coeff
