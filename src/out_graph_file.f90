module graphOut_mod
   !
   ! This module contains the subroutines that store the 3-D graphic files.
   !

   use parallel_mod
   use precision_mod
   use constants, only: one
   use truncation, only: lm_maxMag, n_r_maxMag, n_r_ic_maxMag, lm_max, &
       &                 n_theta_max, n_phi_tot, n_r_max, l_max, minc, &
       &                 n_phi_max, n_r_ic_max, l_axi, nlat_padded,    &
       &                 nRstart, n_r_cmb, nThetaStart, nThetaStop
   use radial_functions, only: r_cmb, orho1, or1, or2, r, r_icb, r_ic, &
       &                       O_r_ic, O_r_ic2
   use physical_parameters, only: ra, ek, pr, prmag, radratio, sigma_ratio, raxi, &
       &                          sc
   use num_param, only: vScale
   use horizontal_data, only: theta_ord, O_sin_theta, n_theta_cal2ord
   use logic, only: l_mag, l_cond_ic, l_PressGraph, l_chemical_conv, l_heat, &
       &            l_save_out
   use output_data, only: runid, n_log_file, log_file, tag
   use sht, only: torpol_to_spat_IC
   use communications, only: gather_f

   implicit none

   private

   integer :: n_graph = 0
   integer :: info
   integer, public :: n_graph_file
#ifdef WITH_MPI
   integer :: graph_mpi_fh
   integer(kind=MPI_OFFSET_KIND) :: size_of_header, n_fields
   public :: graphOut, graphOut_mpi, graphOut_IC, graphOut_mpi_header, &
   &         open_graph_file, close_graph_file, graphOut_header
#else
   public :: graphOut, graphOut_IC, graphOut_header, open_graph_file, &
   &         close_graph_file
#endif

contains

   subroutine open_graph_file(n_time_step, timeScaled)

      !-- Input variables
      integer,  intent(in) :: n_time_step
      real(cp), intent(in) :: timeScaled

      !-- Local variables
      character(len=72) :: graph_file
      character(len=20) :: string

      n_graph = n_graph+1
      write(string, *) n_graph
      graph_file='G_'//trim(adjustl(string))//'.'//tag

      if ( l_master_rank ) then
         write(*,'(1p,/,A,/,A,ES20.10,/,A,i15,/,A,A)')&
         &    " ! Storing graphic file:",             &
         &    "             at time=",timeScaled,     &
         &    "            step no.=",n_time_step,    &
         &    "           into file=",graph_file
         if ( l_save_out ) then
            open(newunit=n_log_file, file=log_file, status='unknown', &
            &    position='append')
         end if
         write(n_log_file,'(1p,/,A,/,A,ES20.10,/,A,i15,/,A,A)') &
         &    " ! Storing graphic file:",                       &
         &    "             at time=",timeScaled,               &
         &    "            step no.=",n_time_step,              &
         &    "           into file=",graph_file
         if ( l_save_out ) close(n_log_file)
      end if

      !-- Setup MPI/IO
      call mpiio_setup(info)

#ifdef WITH_MPI
      call MPI_File_open(comm_r,graph_file,ior(MPI_MODE_WRONLY,MPI_MODE_CREATE),  &
           &             info,graph_mpi_fh,ierr)
#else
      open(newunit=n_graph_file, file=graph_file, status='new',  &
      &    form='unformatted', access='stream')
#endif

   end subroutine open_graph_file
!--------------------------------------------------------------------------------
   subroutine close_graph_file

#ifdef WITH_MPI
         call MPI_File_close(graph_mpi_fh,ierr)
#else
         close(n_graph_file)
#endif

   end subroutine close_graph_file
!--------------------------------------------------------------------------------
   subroutine graphOut(n_r,vr,vt,vp,br,bt,bp,sr,prer,xir)
      !
      !  Output of components of velocity, magnetic field vector, entropy
      !  and composition for graphic outputs.
      !

      !-- Input variables
      integer,  intent(in) :: n_r                    ! radial grod point no.
      real(cp), intent(in) :: vr(:,:),vt(:,:),vp(:,:)
      real(cp), intent(in) :: br(:,:),bt(:,:),bp(:,:)
      real(cp), intent(in) :: sr(:,:),prer(:,:),xir(:,:)

      !-- Local variables:
      integer :: n_phi, n_theta, n_theta_cal, version
      real(cp) :: fac, fac_r
      real(outp) :: dummy(n_theta_max,n_phi_max)

      !-- Write header & colatitudes for n_r=0:
      version = 13

      !-- Calculate and write radial velocity:
      fac=or2(n_r)*vScale*orho1(n_r)
      do n_phi=1,n_phi_max ! do loop over phis
         do n_theta_cal=1,n_theta_max
            n_theta =n_theta_cal2ord(n_theta_cal)
            dummy(n_theta,n_phi)=real(fac*vr(n_theta_cal,n_phi),kind=outp)
         end do
      end do
      write(n_graph_file) dummy(:,:)

      !-- Calculate and write latitudinal velocity:
      fac_r=or1(n_r)*vScale*orho1(n_r)
      do n_phi=1,n_phi_max
         do n_theta_cal=1,n_theta_max
            fac=fac_r*O_sin_theta(n_theta_cal)
            n_theta =n_theta_cal2ord(n_theta_cal)
            dummy(n_theta,n_phi)=real(fac*vt(n_theta_cal,n_phi),kind=outp)
         end do
      end do
      write(n_graph_file) dummy(:,:)

      !-- Calculate and write longitudinal velocity:
      fac_r=or1(n_r)*vScale*orho1(n_r)
      do n_phi=1,n_phi_max
         do n_theta_cal=1,n_theta_max
            fac=fac_r*O_sin_theta(n_theta_cal)
            n_theta =n_theta_cal2ord(n_theta_cal)
            dummy(n_theta,n_phi)=real(fac*vp(n_theta_cal,n_phi),kind=outp)
         end do
      end do
      write(n_graph_file) dummy(:,:)

      if ( l_heat ) then
         !-- Write entropy:
         do n_phi=1,n_phi_max ! do loop over phis
            do n_theta_cal=1,n_theta_max
               n_theta =n_theta_cal2ord(n_theta_cal)
               dummy(n_theta,n_phi)=real(sr(n_theta_cal,n_phi),kind=outp)
            end do
         end do
         write(n_graph_file) dummy(:,:)
      end if

      if ( l_chemical_conv ) then
         !-- Write entropy:
         do n_phi=1,n_phi_max ! do loop over phis
            do n_theta_cal=1,n_theta_max
               n_theta =n_theta_cal2ord(n_theta_cal)
               dummy(n_theta,n_phi)=real(xir(n_theta_cal,n_phi),kind=outp)
            end do
         end do
         write(n_graph_file) dummy(:,:)
      end if

      if ( l_PressGraph ) then
         !-- Write entropy:
         do n_phi=1,n_phi_max ! do loop over phis
            do n_theta_cal=1,n_theta_max
               n_theta =n_theta_cal2ord(n_theta_cal)
               dummy(n_theta,n_phi)=real(prer(n_theta_cal,n_phi),kind=outp)
            end do
         end do
         write(n_graph_file) dummy(:,:)
      end if

      if ( l_mag ) then
         !-- Calculate and write radial magnetic field:
         fac=or2(n_r)*vScale*orho1(n_r)
         do n_phi=1,n_phi_max ! do loop over phis
            do n_theta_cal=1,n_theta_max
               n_theta =n_theta_cal2ord(n_theta_cal)
               dummy(n_theta,n_phi)=real(fac*br(n_theta_cal,n_phi),kind=outp)
            end do
         end do
         write(n_graph_file) dummy(:,:)

         !-- Calculate and write latitudinal magnetic field:
         fac_r=or1(n_r)*vScale*orho1(n_r)
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               fac=fac_r*O_sin_theta(n_theta_cal)
               n_theta =n_theta_cal2ord(n_theta_cal)
               dummy(n_theta,n_phi)=real(fac*bt(n_theta_cal,n_phi),kind=outp)
            end do
         end do
         write(n_graph_file) dummy(:,:)

         !-- Calculate and write longitudinal magnetic field:
         fac_r=or1(n_r)*vScale*orho1(n_r)
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               fac=fac_r*O_sin_theta(n_theta_cal)
               n_theta =n_theta_cal2ord(n_theta_cal)
               dummy(n_theta,n_phi)=real(fac*bp(n_theta_cal,n_phi),kind=outp)
            end do
         end do
         write(n_graph_file) dummy(:,:)
      end if ! l_mag ?

   end subroutine graphOut
!-----------------------------------------------------------------------
   subroutine graphOut_header(time)

      !-- Input variables
      real(cp), intent(in) :: time

      !-- Local variables:
      integer :: n_theta, version

      version = 13

      !-- Write parameters:
      write(n_graph_file) version
      write(n_graph_file) runid
      write(n_graph_file) real(time,outp)
      write(n_graph_file) real(ra,outp), real(pr,outp), real(raxi,outp),    &
      &                   real(sc,outp), real(ek,outp), real(prmag,outp),   &
      &                   real(radratio,outp), real(sigma_ratio,outp)
      write(n_graph_file) n_r_max, n_theta_max, n_phi_tot, minc, n_r_ic_max

      write(n_graph_file) l_heat,l_chemical_conv, l_mag, l_PressGraph, &
      &                   l_cond_ic

      !-- Write colatitudes:
      write(n_graph_file) (real(theta_ord(n_theta),outp), n_theta=1,n_theta_max)

      !-- Write radius:
      write(n_graph_file) real(r,outp)
      if ( l_mag .and. n_r_ic_max > 1 ) then
         write(n_graph_file) real(r_ic,outp)
      end if

   end subroutine graphOut_header
!-------------------------------------------------------------------------------
#ifdef WITH_MPI
   subroutine graphOut_mpi(n_r,vr,vt,vp,br,bt,bp,sr,prer,xir)
      !
      ! MPI version of the graphOut subroutine (use of MPI_IO)
      !

      !-- Input variables:
      integer,  intent(in) :: n_r                      ! radial grod point no.
      real(cp), intent(in) :: vr(nThetaStart:nThetaStop,n_phi_max)
      real(cp), intent(in) :: vt(nThetaStart:nThetaStop,n_phi_max)
      real(cp), intent(in) :: vp(nThetaStart:nThetaStop,n_phi_max)
      real(cp), intent(in) :: br(nThetaStart:nThetaStop,n_phi_max)
      real(cp), intent(in) :: bt(nThetaStart:nThetaStop,n_phi_max)
      real(cp), intent(in) :: bp(nThetaStart:nThetaStop,n_phi_max)
      real(cp), intent(in) :: sr(nThetaStart:nThetaStop,n_phi_max)
      real(cp), intent(in) :: prer(nThetaStart:nThetaStop,n_phi_max)
      real(cp), intent(in) :: xir(nThetaStart:nThetaStop,n_phi_max)

      !-- Local variables:
      integer :: n_phi, n_theta, n_theta_cal
      real(cp) :: fac, fac_r
      real(cp) :: dat(n_theta_max,n_phi_max)
      real(outp) :: dummy(n_theta_max,n_phi_max)

      !$omp critical
      !-- Calculate and write radial velocity:
      if ( n_ranks_theta > 1 ) call gather_f(vr, dat)
      fac=or2(n_r)*vScale*orho1(n_r)
      do n_phi=1,n_phi_max
         do n_theta_cal=1,n_theta_max
            n_theta =n_theta_cal2ord(n_theta_cal)
            dummy(n_theta,n_phi)=real(fac*dat(n_theta_cal,n_phi),kind=outp)
         end do
      end do
      call write_one_field(dummy, graph_mpi_fh, n_phi_max, n_theta_max)

      !-- Calculate and write latitudinal velocity:
      if ( n_ranks_theta > 1 ) call gather_f(vt, dat)
      fac_r=or1(n_r)*vScale*orho1(n_r)
      do n_phi=1,n_phi_max
         do n_theta_cal=1,n_theta_max
            n_theta =n_theta_cal2ord(n_theta_cal)
            fac=fac_r*O_sin_theta(n_theta_cal)
            dummy(n_theta,n_phi)=real(fac*dat(n_theta_cal,n_phi),kind=outp)
         end do
      end do
      call write_one_field(dummy, graph_mpi_fh, n_phi_max, n_theta_max)

      !-- Calculate and write longitudinal velocity:
      if ( n_ranks_theta > 1 ) call gather_f(vp, dat)
      fac_r=or1(n_r)*vScale*orho1(n_r)
      do n_phi=1,n_phi_max
         do n_theta_cal=1,n_theta_max
            n_theta =n_theta_cal2ord(n_theta_cal)
            fac=fac_r*O_sin_theta(n_theta_cal)
            dummy(n_theta,n_phi)=real(fac*dat(n_theta_cal,n_phi),kind=outp)
         end do
      end do
      call write_one_field(dummy, graph_mpi_fh, n_phi_max, n_theta_max)

      !-- Write entropy:
      if ( l_heat ) then
         if ( n_ranks_theta > 1 ) call gather_f(sr, dat)
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta =n_theta_cal2ord(n_theta_cal)
               dummy(n_theta,n_phi)=real(dat(n_theta_cal,n_phi),kind=outp)
            end do
         end do
         call write_one_field(dummy, graph_mpi_fh, n_phi_max, n_theta_max)
      end if

      !-- Write composition:
      if ( l_chemical_conv ) then
         if ( n_ranks_theta > 1 ) call gather_f(xir, dat)
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta =n_theta_cal2ord(n_theta_cal)
               dummy(n_theta,n_phi)=real(dat(n_theta_cal,n_phi),kind=outp)
            end do
         end do
         call write_one_field(dummy, graph_mpi_fh, n_phi_max, n_theta_max)
      end if

      !-- Write pressure:
      if ( l_PressGraph ) then
         if ( n_ranks_theta > 1 ) call gather_f(prer, dat)
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta =n_theta_cal2ord(n_theta_cal)
               dummy(n_theta,n_phi)=real(dat(n_theta_cal,n_phi),kind=outp)
            end do
         end do
         call write_one_field(dummy, graph_mpi_fh, n_phi_max, n_theta_max)
      end if

      if ( l_mag ) then

         !-- Calculate and write radial magnetic field:
         if ( n_ranks_theta > 1 ) call gather_f(br, dat)
         fac=or2(n_r)
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta =n_theta_cal2ord(n_theta_cal)
               dummy(n_theta,n_phi)=real(fac*dat(n_theta_cal,n_phi),kind=outp)
            end do
         end do
         call write_one_field(dummy, graph_mpi_fh, n_phi_max, n_theta_max)

         !-- Calculate and write latitudinal magnetic field:
         if ( n_ranks_theta > 1 ) call gather_f(bt, dat)
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta =n_theta_cal2ord(n_theta_cal)
               fac=or1(n_r)*O_sin_theta(n_theta_cal)
               dummy(n_theta,n_phi)=real(fac*dat(n_theta_cal,n_phi),kind=outp)
            end do
         end do
         call write_one_field(dummy, graph_mpi_fh, n_phi_max, n_theta_max)

         !-- Calculate and write longitudinal magnetic field:
         if ( n_ranks_theta > 1 ) call gather_f(bp, dat)
         do n_phi=1,n_phi_max
            do n_theta_cal=1,n_theta_max
               n_theta =n_theta_cal2ord(n_theta_cal)
               fac=or1(n_r)*O_sin_theta(n_theta_cal)
               dummy(n_theta,n_phi)=real(fac*dat(n_theta_cal,n_phi),kind=outp)
            end do
         end do
         call write_one_field(dummy, graph_mpi_fh, n_phi_max, n_theta_max)

      end if ! l_mag ?
      !$omp end critical

   end subroutine graphOut_mpi
!----------------------------------------------------------------------------
   subroutine graphOut_mpi_header(time)
      !
      ! Writes the header of the G file (MPI version)
      !

      !-- Input variables:
      real(cp), intent(in) :: time

      !-- Local variables:
      integer :: n_theta, version
      integer :: st(MPI_STATUS_SIZE)
      integer(kind=MPI_OFFSET_kind) :: disp, offset

      version = 13
      n_fields = 3
      if ( l_mag ) n_fields = n_fields+3
      if ( l_heat ) n_fields = n_fields+1
      if ( l_PressGraph ) n_fields = n_fields+1
      if ( l_chemical_conv ) n_fields = n_fields+1

      if ( l_master_rank ) then
         !-------- Write parameters:
         call MPI_File_Write(graph_mpi_fh,version,1,MPI_INTEGER,st,ierr)
         call MPI_File_Write(graph_mpi_fh,runid,len(runid),MPI_CHARACTER,st,ierr)
         call MPI_File_Write(graph_mpi_fh,real(time,outp),1,MPI_OUT_REAL,st,ierr)
         call MPI_File_Write(graph_mpi_fh,real(ra,outp),1,MPI_OUT_REAL,st,ierr)
         call MPI_File_Write(graph_mpi_fh,real(pr,outp),1,MPI_OUT_REAL,st,ierr)
         call MPI_File_Write(graph_mpi_fh,real(raxi,outp),1,MPI_OUT_REAL,st,ierr)
         call MPI_File_Write(graph_mpi_fh,real(sc,outp),1,MPI_OUT_REAL,st,ierr)
         call MPI_File_Write(graph_mpi_fh,real(ek,outp),1,MPI_OUT_REAL,st,ierr)
         call MPI_File_Write(graph_mpi_fh,real(prmag,outp),1,MPI_OUT_REAL,st,ierr)
         call MPI_File_Write(graph_mpi_fh,real(radratio,outp),1,MPI_OUT_REAL,st,ierr)
         call MPI_File_Write(graph_mpi_fh,real(sigma_ratio,outp),1,MPI_OUT_REAL,st,ierr)

         call MPI_File_Write(graph_mpi_fh,n_r_max,1,MPI_INTEGER,st,ierr)
         call MPI_File_Write(graph_mpi_fh,n_theta_max,1,MPI_INTEGER,st,ierr)
         call MPI_File_Write(graph_mpi_fh,n_phi_tot,1,MPI_INTEGER,st,ierr)
         call MPI_File_Write(graph_mpi_fh,minc,1,MPI_INTEGER,st,ierr)
         call MPI_File_Write(graph_mpi_fh,n_r_ic_max,1,MPI_INTEGER,st,ierr)

         call MPI_File_Write(graph_mpi_fh,l_heat,1,MPI_LOGICAL,st,ierr)
         call MPI_File_Write(graph_mpi_fh,l_chemical_conv,1,MPI_LOGICAL,st,ierr)
         call MPI_File_Write(graph_mpi_fh,l_mag,1,MPI_LOGICAL,st,ierr)
         call MPI_File_Write(graph_mpi_fh,l_PressGraph,1,MPI_LOGICAL,st,ierr)
         call MPI_File_Write(graph_mpi_fh,l_cond_ic,1,MPI_LOGICAL,st,ierr)

         !-------- Write colatitudes:
         do n_theta=1,n_theta_max
            call MPI_File_Write(graph_mpi_fh,real(theta_ord(n_theta),outp),1, &
                 &              MPI_OUT_REAL,st,ierr)
         end do

         !-------- Write radius:
         call MPI_File_Write(graph_mpi_fh,real(r,outp),n_r_max,MPI_OUT_REAL,st,ierr)
         if ( l_mag .and. n_r_ic_max > 1 ) then
            call MPI_File_Write(graph_mpi_fh,real(r_ic,outp),n_r_ic_max,MPI_OUT_REAL, &
                 &              st,ierr)
         end if

         !-- master gets the displacement
         call MPI_File_get_position(graph_mpi_fh, offset, ierr)
         call MPI_File_get_byte_offset(graph_mpi_fh, offset, size_of_header, ierr)
      end if ! rank 0

      !-- Broadcast the displacement
      call MPI_Bcast(size_of_header, 1, MPI_OFFSET, 0, MPI_COMM_WORLD, ierr)

      !-- Add the size of the arrays
      disp = size_of_header+(nRStart-1)*n_phi_max*n_theta_max*n_fields*SIZEOF_OUT_REAL

      call MPI_File_Set_View(graph_mpi_fh, disp, MPI_OUT_REAL, MPI_OUT_REAL, &
           &                 "native", info, ierr)

   end subroutine graphOut_mpi_header
!----------------------------------------------------------------------------
   subroutine write_one_field(dummy, graph_mpi_fh, n_phis, n_thetas)

      !-- Input variables:
      integer,    intent(in) :: n_thetas       ! number of first colatitude value
      integer,    intent(in) :: n_phis         ! number of logitudes to be printed
      real(outp), intent(in) :: dummy(:,:)     ! data
      integer,    intent(in) :: graph_mpi_fh   ! mpi handle of the mpi file

      !-- MPI related variables
      integer :: st(MPI_STATUS_SIZE), n_counts
      integer(kind=MPI_OFFSET_KIND) :: offset

      n_counts = 0
      do while (n_phis*n_thetas /= n_counts)
          offset = -n_counts*SIZEOF_OUT_REAL
          if (n_counts /= 0 ) call MPI_File_Seek(graph_mpi_fh, offset, MPI_SEEK_CUR, ierr)
          call MPI_File_Write(graph_mpi_fh, dummy, n_phis*n_thetas, &
               &              MPI_OUT_REAL, st, ierr)
          call MPI_Get_Count(st, MPI_OUT_REAL, n_counts, ierr)
      enddo

   end subroutine write_one_field
#endif
!----------------------------------------------------------------------------
   subroutine graphOut_IC(b_ic,db_ic,aj_ic,bICB,l_avg)
      !
      !  Purpose of this subroutine is to write inner core magnetic
      !  field onto graphic output file. If the inner core is
      !  insulating (l_cond_ic=false) the potential field is calculated
      !  from the outer core field at r=r_cmb.
      !  This version assumes that the fields are fully local on the rank
      !  which is calling this routine (usually rank 0).
      !

      !-- Input variables:
      complex(cp), intent(in) :: b_ic(:,:), db_ic(:,:), aj_ic(:,:)
      complex(cp), intent(in) :: bICB(:)
      logical, optional, intent(in) :: l_avg

      !-- Local variables:
      logical :: l_avg_loc
      integer :: nR, nPhi, nTheta, nTheta_cal

      real(cp) :: BrB(nlat_padded,n_phi_max), BtB(nlat_padded,n_phi_max)
      real(cp) :: BpB(nlat_padded,n_phi_max)
      real(outp) :: Br(n_theta_max,n_phi_max),Bt(n_theta_max,n_phi_max)
      real(outp) :: Bp(n_theta_max,n_phi_max)

#ifdef WITH_MPI
      !-- MPI specific variables
      integer(kind=MPI_OFFSET_KIND) :: disp
#endif

      if ( present(l_avg) ) then
         l_avg_loc = l_avg
      else
         l_avg_loc = .false.
      end if

#ifdef WITH_MPI
      !-- One has to bring rank=0 to the end of the file
      if ( .not. l_avg_loc ) then
         disp = size_of_header+n_fields*n_r_max*n_phi_max*n_theta_max* &
         &      SIZEOF_OUT_REAL
         call MPI_File_Set_View(graph_mpi_fh, disp, MPI_OUT_REAL, MPI_OUT_REAL, &
              &                 "native", info, ierr)
      end if
#endif

      if ( l_master_rank ) then
         !-- Loop over all radial levels:
         do nR=1,n_r_ic_max  ! nR=1 is ICB

            if ( l_cond_ic ) then
               call torpol_to_spat_IC(r_ic(nR), r_ICB, b_ic(:, nR), db_ic(:, nR), &
                    &                 aj_ic(:, nR), BrB, BtB, BpB)
            else
               call torpol_to_spat_IC(r_ic(nR), r_ICB, bICB(:), db_ic(:,1), &
                    &                 aj_ic(:,1), BrB, BtB, BpB)
            end if

            do nPhi=1,n_phi_max
               do nTheta_cal=1,n_theta_max
                  nTheta=n_theta_cal2ord(nTheta_cal)
                  Br(nTheta,nPhi)=real(BrB(nTheta_cal,nPhi)*O_r_ic2(nR),kind=outp)
                  Bt(nTheta,nPhi)=real(BtB(nTheta_cal,nPhi)*O_r_ic(nR) * &
                  &                    O_sin_theta(nTheta_cal),kind=outp)
                  Bp(nTheta,nPhi)=real(BpB(nTheta_cal,nPhi)*O_r_ic(nR) * &
                  &                    O_sin_theta(nTheta_cal),kind=outp)
               end do
            end do

            !-- Write radial magnetic field:
#ifdef WITH_MPI
            if ( .not. l_avg_loc ) then
               call write_one_field(Br, graph_mpi_fh, n_phi_max, n_theta_max)
            else
               write(n_graph_file) Br(:,:)
            end if
#else
            write(n_graph_file) Br(:,:)
#endif

            !-- Write latitudinal magnetic field:
#ifdef WITH_MPI
            if ( .not. l_avg_loc ) then
               call write_one_field(Bt, graph_mpi_fh, n_phi_max, n_theta_max)
            else
               write(n_graph_file) Bt(:,:)
            end if
#else
            write(n_graph_file) Bt(:,:)
#endif

            !-- Write longitudinal magnetic field:
#ifdef WITH_MPI
            if ( .not. l_avg_loc ) then
               call write_one_field(Bp, graph_mpi_fh, n_phi_max, n_theta_max)
            else
               write(n_graph_file) Bp(:,:)
            end if
#else
            write(n_graph_file) Bp(:,:)
#endif

         end do  ! Do loop over radial levels nR
      end if ! Only rank==0 writes

   end subroutine graphOut_IC
!----------------------------------------------------------------------------
end module graphOut_mod
