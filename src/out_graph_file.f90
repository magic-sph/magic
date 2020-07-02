#include "perflib_preproc.cpp"
#define ONE_LARGE_BLOCK

module graphOut_mod
   !
   ! This module contains the subroutines that store the 3-D graphic files.
   !

   use parallel_mod
   use precision_mod
   use constants, only: one
   use truncation, only: lm_maxMag, n_r_maxMag, n_r_ic_maxMag, lm_max, &
       &                 n_theta_max, n_phi_tot, n_r_max, l_max, minc, &
       &                 n_phi_max, n_r_ic_max, l_axi, nlat_padded
   use radial_functions, only: r_cmb, orho1, or1, or2, r, r_icb, r_ic, &
       &                       O_r_ic, O_r_ic2
   use radial_data, only: nRstart
   use physical_parameters, only: ra, ek, pr, prmag, radratio, sigma_ratio
   use num_param, only: vScale
   use horizontal_data, only: theta_ord, O_sin_theta
   use logic, only: l_mag, l_cond_ic, l_PressGraph, l_chemical_conv,  &
       &            l_save_out
   use output_data, only: runid, n_log_file, log_file, tag
   use sht, only: torpol_to_spat_IC

   implicit none

   private

   integer :: n_graph = 0
   integer, public :: n_graph_file
#ifdef WITH_MPI
   integer :: graph_mpi_fh
#endif

#ifdef WITH_MPI
   public :: graphOut, graphOut_mpi, graphOut_IC, graphOut_mpi_header, &
   &         open_graph_file, close_graph_file
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
      integer :: info

      n_graph = n_graph+1
      write(string, *) n_graph
      graph_file='G_'//trim(adjustl(string))//'.'//tag

      if ( rank == 0 ) then
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
      call MPI_File_open(MPI_COMM_WORLD,graph_file,             &
           &             IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),  &
           &             MPI_INFO_NULL,graph_mpi_fh,ierr)
#else
      open(newunit=n_graph_file,file=graph_file,status='new',  &
      &    form='unformatted')
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
   subroutine graphOut(time,n_r,vr,vt,vp,br,bt,bp,sr,prer,xir,lGraphHeader)
      !
      !
      !  Output of components of velocity, magnetic field vector and
      !  entropy for graphics.
      !
      !  * n_r: (input) for n_r = 0 a header is written.
      !    for n_r > 0 values at radial level n_r are written
      !
      !  * vr...sr: (input) arrays with grid-point values
      !
      !  * n_theta_start : (input) values are written for theta-points :
      !    ``n_theta_start <= n_theta <= n_theta_start-1+n_theta_block``
      !

      !-- Input variables
      real(cp), intent(in) :: time
      integer,  intent(in) :: n_r                    ! radial grod point no.
      real(cp), intent(in) :: vr(:,:),vt(:,:),vp(:,:)
      real(cp), intent(in) :: br(:,:),bt(:,:),bp(:,:)
      real(cp), intent(in) :: sr(:,:),prer(:,:),xir(:,:)

      logical,  intent(inout) :: lGraphHeader

      !-- Local variables:
      integer :: n_phi         ! counter for longitude
      integer :: n_theta       ! counter for colatitude

      real(cp) :: fac,fac_r
      real(outp) :: dummy(n_phi_max,n_theta_max)

      character(len=20) :: version


      !-- Write header & colatitudes for n_r=0:
      if ( l_chemical_conv ) then
         if ( l_PressGraph ) then
            version='Graphout_Version_6'
         else
            version='Graphout_Version_5'
         end if
      else
         if ( l_PressGraph ) then
            version='Graphout_Version_8'
         else
            version='Graphout_Version_7'
         end if
      end if

      if ( lGraphHeader ) then

         !-------- Write parameters:
         write(n_graph_file) version
         write(n_graph_file) runid
         write(n_graph_file) real(time,outp), real(n_r_max,outp),          &
         &                   real(n_theta_max,outp), real(n_phi_tot,outp), &
         &                   real(n_r_ic_max-1,outp), real(minc,outp),     &
         &                   real(one,outp), real(ra,outp),                &
         &                   real(ek,outp), real(pr,outp),                 &
         &                   real(prmag,outp), real(radratio,outp),        &
         &                   real(sigma_ratio,outp)

         !-------- Write colatitudes:
         write(n_graph_file) (real(theta_ord(n_theta),outp), n_theta=1,n_theta_max)

         lGraphHeader=.false.

      else  ! Call not for writing header

         !*******************************************************************
         !  define CRITICAL section, so that the write statements do not
         !  get mixed up
         !*******************************************************************

         write(n_graph_file) real(n_r-1,outp),real(r(n_r)/r(1),outp), &
         &                   real(1,outp),real(n_theta_max,outp)

         !-- Write entropy:
         do n_phi=1,n_phi_max ! do loop over phis
            do n_theta=1,n_theta_max,2
               dummy(n_phi,n_theta)  =real(sr(n_theta,n_phi),kind=outp)   ! NHS
               dummy(n_phi,n_theta+1)=real(sr(n_theta+1,n_phi),kind=outp) ! SHS
            end do
         end do
         call graph_write(n_phi_max,n_theta_max,dummy,n_graph_file)

         !-- Calculate and write radial velocity:
         fac=or2(n_r)*vScale*orho1(n_r)
         do n_phi=1,n_phi_max ! do loop over phis
            do n_theta=1,n_theta_max,2
               dummy(n_phi,n_theta)  =real(fac*vr(n_theta,n_phi),kind=outp)
               dummy(n_phi,n_theta+1)=real(fac*vr(n_theta+1,n_phi),kind=outp)
            end do
         end do
         call graph_write(n_phi_max,n_theta_max,dummy,n_graph_file)


         !-- Calculate and write latitudinal velocity:
         fac_r=or1(n_r)*vScale*orho1(n_r)
         do n_phi=1,n_phi_max
            do n_theta=1,n_theta_max,2
               fac=fac_r*O_sin_theta(n_theta)
               dummy(n_phi,n_theta)  =real(fac*vt(n_theta,n_phi),kind=outp)
               dummy(n_phi,n_theta+1)=real(fac*vt(n_theta+1,n_phi),kind=outp)
            end do
         end do
         call graph_write(n_phi_max,n_theta_max,dummy,n_graph_file)

         !-- Calculate and write longitudinal velocity:
         fac_r=or1(n_r)*vScale*orho1(n_r)
         do n_phi=1,n_phi_max
            do n_theta=1,n_theta_max,2
               fac=fac_r*O_sin_theta(n_theta)
               dummy(n_phi,n_theta)  =real(fac*vp(n_theta,n_phi),kind=outp)
               dummy(n_phi,n_theta+1)=real(fac*vp(n_theta+1,n_phi),kind=outp)
            end do
         end do
         call graph_write(n_phi_max,n_theta_max,dummy,n_graph_file)

         if ( version == 'Graphout_Version_5' .or. version == 'Graphout_Version_6') then
            !-- Write composition:
            do n_phi=1,n_phi_max ! do loop over phis
               do n_theta=1,n_theta_max,2
                  dummy(n_phi,n_theta)  =real(xir(n_theta,n_phi),kind=outp)   ! NHS
                  dummy(n_phi,n_theta+1)=real(xir(n_theta+1,n_phi),kind=outp) ! SHS
               end do
            end do
            call graph_write(n_phi_max,n_theta_max,dummy,n_graph_file)
         end if

         if ( version == 'Graphout_Version_6' .or. version == 'Graphout_Version_8') then
            !-- Write pressure:
            do n_phi=1,n_phi_max ! do loop over phis
               do n_theta=1,n_theta_max,2
                  dummy(n_phi,n_theta)  =real(prer(n_theta,n_phi),kind=outp)   ! NHS
                  dummy(n_phi,n_theta+1)=real(prer(n_theta+1,n_phi),kind=outp) ! SHS
               end do
            end do
            call graph_write(n_phi_max,n_theta_max,dummy,n_graph_file)
         end if

         if ( l_mag ) then

            !-- Calculate and write radial magnetic field:
            fac=or2(n_r)
            do n_phi=1,n_phi_max
               do n_theta=1,n_theta_max,2
                  dummy(n_phi,n_theta)  =real(fac*br(n_theta,n_phi),kind=outp)
                  dummy(n_phi,n_theta+1)=real(fac*br(n_theta+1,n_phi),kind=outp)
               end do
            end do
            call graph_write(n_phi_max,n_theta_max,dummy,n_graph_file)

            !-- Calculate and write latitudinal magnetic field:
            do n_phi=1,n_phi_max
               do n_theta=1,n_theta_max,2
               fac=or1(n_r)*O_sin_theta(n_theta)
                  dummy(n_phi,n_theta)  =real(fac*bt(n_theta,n_phi),kind=outp)
                  dummy(n_phi,n_theta+1)=real(fac*bt(n_theta+1,n_phi),kind=outp)
               end do
            end do
            call graph_write(n_phi_max,n_theta_max,dummy,n_graph_file)

            !-- Calculate and write longitudinal magnetic field:
            do n_phi=1,n_phi_max
               do n_theta=1,n_theta_max,2
                  fac=or1(n_r)*O_sin_theta(n_theta)
                  dummy(n_phi,n_theta)  =real(fac*bp(n_theta,n_phi),kind=outp)
                  dummy(n_phi,n_theta+1)=real(fac*bp(n_theta+1,n_phi),kind=outp)
               end do
            end do
            call graph_write(n_phi_max,n_theta_max,dummy,n_graph_file)

         end if ! l_mag ?

         !-- End of CRITICAL section ********************************************

      end if

   end subroutine graphOut
!-----------------------------------------------------------------------
   subroutine graphOut_header(time)

      !-- Input variables
      real(cp), intent(in) :: time

      !-- Local variables:
      character(len=20) :: version
      integer :: n_theta


      !-- Write header & colatitudes for n_r=0:
      if ( l_chemical_conv ) then
         if ( l_PressGraph ) then
            version='Graphout_Version_6'
         else
            version='Graphout_Version_5'
         end if
      else
         if ( l_PressGraph ) then
            version='Graphout_Version_8'
         else
            version='Graphout_Version_7'
         end if
      end if

      !-------- Write parameters:
      write(n_graph_file) version
      write(n_graph_file) runid
      write(n_graph_file) real(time,outp), real(n_r_max,outp),          &
      &                   real(n_theta_max,outp), real(n_phi_tot,outp), &
      &                   real(n_r_ic_max-1,outp), real(minc,outp),     &
      &                   real(one,outp), real(ra,outp),                &
      &                   real(ek,outp), real(pr,outp),                 &
      &                   real(prmag,outp), real(radratio,outp),        &
      &                   real(sigma_ratio,outp)

      !-------- Write colatitudes:
      write(n_graph_file) (real(theta_ord(n_theta),outp), n_theta=1,n_theta_max)

   end subroutine graphOut_header
!-------------------------------------------------------------------------------
#ifdef WITH_MPI
   subroutine graphOut_mpi(time,n_r,vr,vt,vp,br,bt,bp,sr,prer,xir,lGraphHeader)
      !
      ! MPI version of the graphOut subroutine (use of MPI_IO)
      !

      !-- Input variables:
      real(cp), intent(in) :: time
      integer,  intent(in) :: n_r                      ! radial grod point no.
      real(cp), intent(in) :: vr(:,:),vt(:,:),vp(:,:)
      real(cp), intent(in) :: br(:,:),bt(:,:),bp(:,:)
      real(cp), intent(in) :: sr(:,:),prer(:,:),xir(:,:)

      logical, intent(inout) :: lGraphHeader

      !-- Local variables:
      integer :: n_phi         ! counter for longitude
      integer :: n_theta       ! counter for colatitude

      real(cp) :: fac,fac_r
      real(outp) :: dummy(n_phi_max,n_theta_max)

      character(len=20) :: version

      ! MPI related variables
      !integer :: info
      integer :: status(MPI_STATUS_SIZE)
      !character(len=MPI_MAX_ERROR_STRING) :: error_string
      !integer :: count
      integer :: bytes_written!,length_of_error
      integer :: size_of_header, size_of_data_per_r
      integer :: size_of_data_per_thetaB
      integer(kind=MPI_OFFSET_kind) :: disp
      integer :: etype,filetype
      character(len=MPI_MAX_DATAREP_STRING) :: datarep
      ! end of MPI related variables

      if ( l_chemical_conv ) then
         if ( l_PressGraph ) then
            version='Graphout_Version_12'
         else
            version='Graphout_Version_11'
         end if
      else
         if ( l_PressGraph ) then
            version='Graphout_Version_10'
         else
            version='Graphout_Version_9'
         end if
      end if

      !$OMP CRITICAL
      if ( lGraphHeader ) then
         size_of_header = 8+len(version)+8+len(runid)+8+ &
         &                13*SIZEOF_OUT_REAL+8+n_theta_max*SIZEOF_OUT_REAL

#ifdef ONE_LARGE_BLOCK
         size_of_data_per_thetaB = 8+4*SIZEOF_OUT_REAL+3* &
         &                         (8+n_phi_max*SIZEOF_OUT_REAL*n_theta_max)
         if ( version=='Graphout_Version_10' .or. version=='Graphout_Version_11') then
            size_of_data_per_thetaB = size_of_data_per_thetaB + &
            &                      (8+n_phi_max*SIZEOF_OUT_REAL*n_theta_max)
         else if ( version=='Graphout_Version_12') then
            size_of_data_per_thetaB = size_of_data_per_thetaB + &
                                   2*(8+n_phi_max*SIZEOF_OUT_REAL*n_theta_max)
         end if

         if ( l_mag ) size_of_data_per_thetaB = size_of_data_per_thetaB + &
                      &     3*(8+n_phi_max*SIZEOF_OUT_REAL*n_theta_max)
#else
         size_of_data_per_thetaB = 8+4*SIZEOF_OUT_REAL+ &
         &                         3*(8+n_phi_max*SIZEOF_OUT_REAL)*n_theta_max
         if ( version=='Graphout_Version_10' .or. version=='Graphout_Version_11') then
            size_of_data_per_thetaB = size_of_data_per_thetaB + &
            &                      (8+n_phi_max*SIZEOF_OUT_REAL)*n_theta_max
         else if ( version=='Graphout_Version_12') then
            size_of_data_per_thetaB = size_of_data_per_thetaB + &
            &                      2*(8+n_phi_max*SIZEOF_OUT_REAL)*n_theta_max
         end if

         if ( l_mag ) size_of_data_per_thetaB = size_of_data_per_thetaB + &
                            3*(8+n_phi_max*SIZEOF_OUT_REAL)*n_theta_max
#endif
         size_of_data_per_r = size_of_data_per_thetaB

         if ( rank == 0 ) then
            ! rank zero writes the Header
            disp = 0
            call MPI_FILE_SET_VIEW(graph_mpi_fh,disp,MPI_CHARACTER, &
                 &                 MPI_CHARACTER,"native",MPI_INFO_NULL,ierr)
         else
            disp = size_of_header+(nRstart-1)*size_of_data_per_r
            call MPI_FILE_SET_VIEW(graph_mpi_fh,disp,MPI_CHARACTER, &
                 &                 MPI_CHARACTER,"native",MPI_INFO_NULL,ierr)
         end if

         call MPI_FILE_GET_VIEW(graph_mpi_fh,disp,etype,filetype,datarep,ierr)

         bytes_written = 0
         !-- Write header & colatitudes for n_r=0:
         if ( rank == 0 ) then
            !-------- Write parameters:
            call MPI_FILE_WRITE(graph_mpi_fh,len(version),1,MPI_INTEGER, &
                 &              status,ierr)
            call MPI_FILE_WRITE(graph_mpi_fh,version,len(version), &
                 &              MPI_CHARACTER,status,ierr)
            call MPI_FILE_WRITE(graph_mpi_fh,len(version),1,MPI_INTEGER, &
                 &              status,ierr)

            call MPI_FILE_WRITE(graph_mpi_fh,len(runid),1,MPI_INTEGER,  &
                 &              status,ierr)
            call MPI_FILE_WRITE(graph_mpi_fh,runid,len(runid), &
                 &              MPI_CHARACTER,status,ierr)
            call MPI_FILE_WRITE(graph_mpi_fh,len(runid),1,MPI_INTEGER, &
                 &              status,ierr)

            call MPI_FILE_WRITE(graph_mpi_fh,13*SIZEOF_OUT_REAL,1, &
                 &              MPI_INTEGER,status,ierr)
            call MPI_FILE_WRITE(graph_mpi_fh,real(time,outp),1,MPI_OUT_REAL, &
                 &              status,ierr)
            call MPI_FILE_WRITE(graph_mpi_fh,real(n_r_max,outp),1, &
                 &              MPI_OUT_REAL,status,ierr)
            call MPI_FILE_WRITE(graph_mpi_fh,real(n_theta_max,outp),1, &
                 &              MPI_OUT_REAL,status,ierr)
            call MPI_FILE_WRITE(graph_mpi_fh,real(n_phi_tot,outp),1, &
                 &              MPI_OUT_REAL,status,ierr)
            call MPI_FILE_WRITE(graph_mpi_fh,real(n_r_ic_max-1,outp),1, &
                 &              MPI_OUT_REAL,status,ierr)
            call MPI_FILE_WRITE(graph_mpi_fh,real(minc,outp),1,MPI_OUT_REAL, &
                 &              status,ierr)
            call MPI_FILE_WRITE(graph_mpi_fh,real(one,outp),1, &
                 &              MPI_OUT_REAL,status,ierr)
            call MPI_FILE_WRITE(graph_mpi_fh,real(ra,outp),1,MPI_OUT_REAL, &
                 &              status,ierr)
            call MPI_FILE_WRITE(graph_mpi_fh,real(ek,outp),1,MPI_OUT_REAL, &
                 &              status,ierr)
            call MPI_FILE_WRITE(graph_mpi_fh,real(pr,outp),1,MPI_OUT_REAL, &
                 &              status,ierr)
            call MPI_FILE_WRITE(graph_mpi_fh,real(prmag,outp),1,MPI_OUT_REAL, &
                 &              status,ierr)
            call MPI_FILE_WRITE(graph_mpi_fh,real(radratio,outp),1, &
                 &              MPI_OUT_REAL,status,ierr)
            call MPI_FILE_WRITE(graph_mpi_fh,real(sigma_ratio,outp),1, &
                 &              MPI_OUT_REAL,status,ierr)
            call MPI_FILE_WRITE(graph_mpi_fh,13*SIZEOF_OUT_REAL,1, &
                 &              MPI_INTEGER,status,ierr)

            !-------- Write colatitudes:
            call MPI_FILE_WRITE(graph_mpi_fh,n_theta_max*SIZEOF_OUT_REAL,1, &
                 &              MPI_INTEGER,status,ierr)
            do n_theta=1,n_theta_max
               call MPI_FILE_WRITE(graph_mpi_fh,real(theta_ord(n_theta),outp),1,&
                    &              MPI_OUT_REAL,status,ierr)
            end do
            call MPI_FILE_WRITE(graph_mpi_fh,n_theta_max*SIZEOF_OUT_REAL,1, &
                 &              MPI_INTEGER,status,ierr)

         end if
         lGraphHeader=.false.
      else  ! Call not for writing header

         !PERFON('mw_data')
         bytes_written=0

         !-- Determine radius and thetas in this block:
         call MPI_FILE_WRITE(graph_mpi_fh,4*SIZEOF_OUT_REAL,1,MPI_INTEGER, &
              &              status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,real(n_r-1,outp),1,MPI_OUT_REAL, &
              &              status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,real(r(n_r)/r(1),outp),1, &
              &              MPI_OUT_REAL,status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,real(1,outp),1, &
              &              MPI_OUT_REAL,status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,real(n_theta_max,outp),1, &
              &              MPI_OUT_REAL,status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,4*SIZEOF_OUT_REAL,1,MPI_INTEGER,&
              &              status,ierr)

         !-- Write entropy:
         do n_phi=1,n_phi_max
            do n_theta=1,n_theta_max,2
               dummy(n_phi,n_theta)  =real(sr(n_theta,n_phi),kind=outp)   ! NHS
               dummy(n_phi,n_theta+1)=real(sr(n_theta+1,n_phi),kind=outp) ! SHS
            end do
         end do
         call graph_write_mpi(n_phi_max,n_theta_max,dummy,graph_mpi_fh)

         !-- Calculate and write radial velocity:
         fac=or2(n_r)*vScale*orho1(n_r)
         do n_phi=1,n_phi_max
            do n_theta=1,n_theta_max,2
               dummy(n_phi,n_theta)  =real(fac*vr(n_theta,n_phi),kind=outp)
               dummy(n_phi,n_theta+1)=real(fac*vr(n_theta+1,n_phi),kind=outp)
            end do
         end do
         call graph_write_mpi(n_phi_max,n_theta_max,dummy,graph_mpi_fh)

         !-- Calculate and write latitudinal velocity:
         fac_r=or1(n_r)*vScale*orho1(n_r)
         do n_phi=1,n_phi_max
            do n_theta=1,n_theta_max,2
               fac=fac_r*O_sin_theta(n_theta)
               dummy(n_phi,n_theta)  =real(fac*vt(n_theta,n_phi),kind=outp)
               dummy(n_phi,n_theta+1)=real(fac*vt(n_theta+1,n_phi),kind=outp)
            end do
         end do
         call graph_write_mpi(n_phi_max,n_theta_max,dummy,graph_mpi_fh)

         !-- Calculate and write longitudinal velocity:
         fac_r=or1(n_r)*vScale*orho1(n_r)
         do n_phi=1,n_phi_max
            do n_theta=1,n_theta_max,2
               fac=fac_r*O_sin_theta(n_theta)
               dummy(n_phi,n_theta)  =real(fac*vp(n_theta,n_phi),kind=outp)
               dummy(n_phi,n_theta+1)=real(fac*vp(n_theta+1,n_phi),kind=outp)
            end do
         end do
         call graph_write_mpi(n_phi_max,n_theta_max,dummy,graph_mpi_fh)

         !-- Write composition:
         if ( version == 'Graphout_Version_11' .or. version == 'Graphout_Version_12' ) then
            do n_phi=1,n_phi_max
               do n_theta=1,n_theta_max,2
                  dummy(n_phi,n_theta)  =real(xir(n_theta,n_phi),kind=outp)   ! NHS
                  dummy(n_phi,n_theta+1)=real(xir(n_theta+1,n_phi),kind=outp) ! SHS
               end do
            end do
            call graph_write_mpi(n_phi_max,n_theta_max,dummy,graph_mpi_fh)
         end if

         !-- Write pressure:
         if ( version == 'Graphout_Version_10' .or. version == 'Graphout_Version_12' ) then
            do n_phi=1,n_phi_max
               do n_theta=1,n_theta_max,2
                  dummy(n_phi,n_theta)  =real(prer(n_theta,n_phi),kind=outp)   ! NHS
                  dummy(n_phi,n_theta+1)=real(prer(n_theta+1,n_phi),kind=outp) ! SHS
               end do
            end do
            call graph_write_mpi(n_phi_max,n_theta_max,dummy,graph_mpi_fh)
         end if

         if ( l_mag ) then

            !-- Calculate and write radial magnetic field:
            fac=or2(n_r)
            do n_phi=1,n_phi_max
               do n_theta=1,n_theta_max,2
                  dummy(n_phi,n_theta)  =real(fac*br(n_theta,n_phi),kind=outp)
                  dummy(n_phi,n_theta+1)=real(fac*br(n_theta+1,n_phi),kind=outp)
               end do
            end do
            call graph_write_mpi(n_phi_max,n_theta_max,dummy,graph_mpi_fh)

            !-- Calculate and write latitudinal magnetic field:
            do n_phi=1,n_phi_max
               do n_theta=1,n_theta_max,2
                  fac=or1(n_r)*O_sin_theta(n_theta)
                  dummy(n_phi,n_theta)  =real(fac*bt(n_theta,n_phi),kind=outp)
                  dummy(n_phi,n_theta+1)=real(fac*bt(n_theta+1,n_phi),kind=outp)
               end do
            end do
            call graph_write_mpi(n_phi_max,n_theta_max,dummy,graph_mpi_fh)

            !-- Calculate and write longitudinal magnetic field:
            do n_phi=1,n_phi_max
               do n_theta=1,n_theta_max,2
                  fac=or1(n_r)*O_sin_theta(n_theta)
                  dummy(n_phi,n_theta)  =real(fac*bp(n_theta,n_phi),kind=outp)
                  dummy(n_phi,n_theta+1)=real(fac*bp(n_theta+1,n_phi),kind=outp)
               end do
            end do
            call graph_write_mpi(n_phi_max,n_theta_max,dummy,graph_mpi_fh)

         end if ! l_mag ?

         !PERFOFF
      end if
      !$OMP END CRITICAL

   end subroutine graphOut_mpi
!----------------------------------------------------------------------------
   subroutine graphOut_mpi_header(time)
      !
      ! Writes the header (MPI version)
      !

      !-- Input variables:
      real(cp), intent(in) :: time

      !-- Local variables:
      integer :: n_theta       ! counter for colatitude
      character(len=20) :: version

      !-- MPI related variables
      integer :: status(MPI_STATUS_SIZE)
      integer :: bytes_written
      integer :: size_of_header, size_of_data_per_r
      integer :: size_of_data_per_thetaB
      integer(kind=MPI_OFFSET_kind) :: disp
      integer :: etype,filetype
      character(len=MPI_MAX_DATAREP_STRING) :: datarep
      ! end of MPI related variables

      !----- Unformatted output:
      if ( l_chemical_conv ) then
         if ( l_PressGraph ) then
            version='Graphout_Version_12'
         else
            version='Graphout_Version_11'
         end if
      else
         if ( l_PressGraph ) then
            version='Graphout_Version_10'
         else
            version='Graphout_Version_9'
         end if
      end if

      size_of_header = 8+len(version)+8+len(runid)+8+13*SIZEOF_OUT_REAL+8+ &
      &                n_theta_max*SIZEOF_OUT_REAL

#ifdef ONE_LARGE_BLOCK
      size_of_data_per_thetaB = 8+4*SIZEOF_OUT_REAL+4* &
      &                         (8+n_phi_max*SIZEOF_OUT_REAL*n_theta_max)
      if ( version=='Graphout_Version_10' .or. version=='Graphout_Version_11') then
         size_of_data_per_thetaB = size_of_data_per_thetaB + &
         &                         (8+n_phi_max*SIZEOF_OUT_REAL*n_theta_max)
      else if ( version=='Graphout_Version_12') then
         size_of_data_per_thetaB = size_of_data_per_thetaB + &
         &                         2*(8+n_phi_max*SIZEOF_OUT_REAL*n_theta_max)
      end if


      if ( l_mag ) size_of_data_per_thetaB = size_of_data_per_thetaB + &
                   &               3*(8+n_phi_max*SIZEOF_OUT_REAL*n_theta_max)
#else
      size_of_data_per_thetaB = 8+4*SIZEOF_OUT_REAL+4* &
      &                        (8+n_phi_max*SIZEOF_OUT_REAL)*n_theta_max
      if ( version=='Graphout_Version_10' .or. version=='Graphout_Version_11') then
         size_of_data_per_thetaB = size_of_data_per_thetaB + &
         &                         (8+n_phi_max*SIZEOF_OUT_REAL)*n_theta_max
      else if ( version=='Graphout_Version_12') then
         size_of_data_per_thetaB = size_of_data_per_thetaB + &
         &                         2*(8+n_phi_max*SIZEOF_OUT_REAL)*n_theta_max
      end if
      if ( l_mag ) size_of_data_per_thetaB = size_of_data_per_thetaB + &
                   & 3*(8+n_phi_max*SIZEOF_OUT_REAL)*n_theta_max
#endif
      size_of_data_per_r = size_of_data_per_thetaB

      if ( rank == 0 ) then
         ! rank zero writes the Header
         disp = 0
         call MPI_FILE_SET_VIEW(graph_mpi_fh,disp,MPI_CHARACTER, &
                                MPI_CHARACTER,"native",MPI_INFO_NULL,ierr)
      else
         disp = size_of_header+(nRstart-1)*size_of_data_per_r
         call MPI_FILE_SET_VIEW(graph_mpi_fh,disp,&
              & MPI_CHARACTER,MPI_CHARACTER,"native",MPI_INFO_NULL,ierr)
      end if

      call mpi_file_get_view(graph_mpi_fh,disp,etype,filetype,datarep,ierr)

      bytes_written = 0
      !-- Write header & colatitudes for n_r=0:
      if ( rank == 0 ) then

         !-------- Write parameters:
         call MPI_FILE_WRITE(graph_mpi_fh,len(version),1,MPI_INTEGER,status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,version,len(version), &
              &              MPI_CHARACTER,status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,len(version),1,MPI_INTEGER,status,ierr)

         call MPI_FILE_WRITE(graph_mpi_fh,len(runid),1,MPI_INTEGER,status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,runid,len(runid), &
              &              MPI_CHARACTER,status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,len(runid),1,MPI_INTEGER,status,ierr)

         call MPI_FILE_WRITE(graph_mpi_fh,13*SIZEOF_OUT_REAL,1, &
              &              MPI_INTEGER,status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,real(time,outp),1,MPI_OUT_REAL, &
              &              status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,real(n_r_max,outp),1,MPI_OUT_REAL,&
              &              status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,real(n_theta_max,outp),1, &
              &              MPI_OUT_REAL,status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,real(n_phi_tot,outp),1,MPI_OUT_REAL, &
              &              status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,real(n_r_ic_max-1,outp),1, &
              &              MPI_OUT_REAL,status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,real(minc,outp),1,MPI_OUT_REAL, &
              &              status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,real(one,outp),1,MPI_OUT_REAL,&
              &              status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,real(ra,outp),1,MPI_OUT_REAL, &
              &              status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,real(ek,outp),1,MPI_OUT_REAL, &
              &              status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,real(pr,outp),1,MPI_OUT_REAL, &
              &              status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,real(prmag,outp),1,MPI_OUT_REAL, &
              &              status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,real(radratio,outp),1,MPI_OUT_REAL,&
              &              status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,real(sigma_ratio,outp),1, &
              &              MPI_OUT_REAL,status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,13*SIZEOF_OUT_REAL,1, &
              &              MPI_INTEGER,status,ierr)

         !-------- Write colatitudes:
         call MPI_FILE_WRITE(graph_mpi_fh,n_theta_max*SIZEOF_OUT_REAL,1, &
              &              MPI_INTEGER,status,ierr)
         do n_theta=1,n_theta_max
            call MPI_FILE_WRITE(graph_mpi_fh,real(theta_ord(n_theta),outp),1, &
                 &              MPI_OUT_REAL,status,ierr)
         end do
         call MPI_FILE_WRITE(graph_mpi_fh,n_theta_max*SIZEOF_OUT_REAL,1, &
              &              MPI_INTEGER,status,ierr)

      end if

   end subroutine graphOut_mpi_header
#endif
!----------------------------------------------------------------------------
   subroutine graphOut_IC(b_ic,db_ic,ddb_ic,aj_ic,dj_ic,bICB,l_avg)
      !
      !  Purpose of this subroutine is to write inner core magnetic
      !  field onto graphic output file. If the inner core is
      !  insulating (l_cond_ic=false) the potential field is calculated
      !  from the outer core field at r=r_cmb.
      !  This version assumes that the fields are fully local on the rank
      !  which is calling this routine (usually rank 0).
      !

      !-- Input variables:
      complex(cp), intent(in) :: b_ic(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: db_ic(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: ddb_ic(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: aj_ic(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: dj_ic(lm_maxMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: bICB(lm_maxMag)
      logical, optional, intent(in) :: l_avg

      !-- Local variables:
      logical :: l_avg_loc
      integer :: nR, nPhi, nTheta

      real(cp) :: BrB(nlat_padded,n_phi_max), BtB(nlat_padded,n_phi_max)
      real(cp) :: BpB(nlat_padded,n_phi_max)
      real(outp) :: Br(n_phi_max,n_theta_max)
      real(outp) :: Bt(n_phi_max,n_theta_max)
      real(outp) :: Bp(n_phi_max,n_theta_max)

#ifdef WITH_MPI
      ! MPI specific variables
      integer :: status(MPI_STATUS_SIZE)
      integer(kind=MPI_OFFSET_KIND) :: offset
      ! end MPI variables
#endif

      if ( present(l_avg) ) then
         l_avg_loc = l_avg
      else
         l_avg_loc = .false.
      end if

#ifdef WITH_MPI
      !-- One has to bring rank=0 to the end of the file
      if ( .not. l_avg_loc ) then
         offset = 0
         call MPI_File_Seek(graph_mpi_fh, offset, MPI_SEEK_END, ierr)
      end if
#endif

      !-- Loop over all radial levels:

      do nR=2,n_r_ic_max  ! nR=1 is ICB

         if ( l_cond_ic ) then
            call torpol_to_spat_IC(r_ic(nR), r_ICB, b_ic(:, nR), db_ic(:, nR), &
                 &                 aj_ic(:, nR), BrB, BtB, BpB)
         else
            call torpol_to_spat_IC(r_ic(nR), r_ICB, bICB(:),db_ic(:,1), &
                 &                 aj_ic(:,1), BrB, BtB, BpB)
         end if

         do nPhi=1,n_phi_max
            do nTheta=1,n_theta_max
               Br(nPhi,nTheta)=real(BrB(nTheta,nPhi)*O_r_ic2(nR),kind=outp)
               Bt(nPhi,nTheta)=real(BtB(nTheta,nPhi)*O_r_ic(nR) * &
               &                    O_sin_theta(nTheta),kind=outp)
               Bp(nPhi,nTheta)=real(BpB(nTheta,nPhi)*O_r_ic(nR) * &
               &                    O_sin_theta(nTheta),kind=outp)
            end do
         end do

#ifdef WITH_MPI
         ! in process n_procs-1 the last oc fields have been written,
         ! Now just append on this process.
         if ( .not. l_avg_loc ) then
            call MPI_FILE_WRITE(graph_mpi_fh,4*SIZEOF_OUT_REAL,1,  &
                 &              MPI_INTEGER,status,ierr)
            call MPI_FILE_WRITE(graph_mpi_fh,real(n_r_max+nR-2,outp),1, &
                 &              MPI_OUT_REAL,status,ierr)
            call MPI_FILE_WRITE(graph_mpi_fh,real(r_ic(nR)/r_cmb,outp),1, &
                 &              MPI_OUT_REAL,status,ierr)
            call MPI_FILE_WRITE(graph_mpi_fh,1.e0_outp,1,MPI_OUT_REAL, &
                 &              status,ierr)
            call MPI_FILE_WRITE(graph_mpi_fh,real(n_theta_max,outp),1, &
                 &              MPI_OUT_REAL,status,ierr)
            call MPI_FILE_WRITE(graph_mpi_fh,4*SIZEOF_OUT_REAL,1, &
                 &              MPI_INTEGER,status,ierr)
         else
            write(n_graph_file) real(n_r_max+nR-2,outp),real(r_ic(nR)/r_cmb,outp),&
                 &              1.e0_outp,real(n_theta_max,outp)
         end if
#else
         write(n_graph_file) real(n_r_max+nR-2,outp),real(r_ic(nR)/r_cmb,outp),&
              &              1.e0_outp,real(n_theta_max,outp)
#endif


         !-- Write radial magnetic field:
#ifdef WITH_MPI
         if ( .not. l_avg_loc ) then
            call graph_write_mpi(n_phi_max,n_theta_max,Br,graph_mpi_fh)
         else
            call graph_write(n_phi_max,n_theta_max,Br,n_graph_file)
         end if
#else
         call graph_write(n_phi_max,n_theta_max,Br,n_graph_file)
#endif

         !-- Write latitudinal magnetic field:
#ifdef WITH_MPI
         if ( .not. l_avg_loc ) then
            call graph_write_mpi(n_phi_max,n_theta_max,Bt,graph_mpi_fh)
         else
            call graph_write(n_phi_max,n_theta_max,Bt,n_graph_file)
         end if
#else
         call graph_write(n_phi_max,n_theta_max,Bt,n_graph_file)
#endif

         !-- Write longitudinal magnetic field:
#ifdef WITH_MPI
         if ( .not. l_avg_loc ) then
            call graph_write_mpi(n_phi_max,n_theta_max,Bp,graph_mpi_fh)
         else
            call graph_write(n_phi_max,n_theta_max,Bp,n_graph_file)
         end if
#else
         call graph_write(n_phi_max,n_theta_max,Bp,n_graph_file)
#endif

      end do  ! Do loop over radial levels nR

   end subroutine graphOut_IC
!-----------------------------------------------------------------------
   subroutine graph_write(n_phis,n_thetas,dummy,n_graph_file)
      !
      !  This subroutine writes the data for one theta-band
      !  (stored in 'dummy'). Version May, 5, 2000.
      !

      !-- Input variables:
      integer,    intent(in) :: n_thetas            ! number of first colatitude value
      integer,    intent(in) :: n_phis              ! number of logitudes to be printed
      real(outp), intent(in) :: dummy(n_phi_max,*)  ! data
      integer,    intent(in) :: n_graph_file        ! output unit

      !-- Local variables:
      integer :: n_phi,n_theta


      !PERFON('gwrite')
      do n_theta=1,n_thetas
         write(n_graph_file) (dummy(n_phi,n_theta),n_phi=1,n_phis)
      end do
      !PERFOFF

   end subroutine graph_write
!------------------------------------------------------------------------------
#ifdef WITH_MPI
   subroutine graph_write_mpi(n_phis,n_thetas,dummy,graph_mpi_fh)

      !-- Input variables
      integer,    intent(in) :: n_thetas          ! number of first colatitude value
      integer,    intent(in) :: n_phis            ! number of logitudes to be printed
      real(outp), intent(in) :: dummy(n_phi_max,*)! data
      integer,    intent(in) :: graph_mpi_fh      ! mpi handle of the mpi file

      !-- Local variables:
#ifndef ONE_LARGE_BLOCK
      integer :: n_theta
#endif

      !-- MPI related variables
      integer :: status(MPI_STATUS_SIZE), count
      integer(kind=MPI_OFFSET_KIND) :: offset

#ifdef ONE_LARGE_BLOCK
      call MPI_FILE_WRITE(graph_mpi_fh,n_phis*n_thetas*SIZEOF_OUT_REAL,1, &
           &              MPI_INTEGER,status,ierr)
      ! call MPI_FILE_WRITE(graph_mpi_fh,dummy(:,1:n_thetas),n_phis*n_thetas, &
      !                     MPI_OUT_REAL,status,ierr)
      count = 0
      do while (n_phis*n_thetas /= count)
          offset = -count*SIZEOF_OUT_REAL
          if (count /= 0 ) call MPI_File_seek(graph_mpi_fh, offset, MPI_SEEK_CUR, ierr)
          call MPI_File_write(graph_mpi_fh,dummy(:,1:n_thetas),n_phis*n_thetas, &
               &              MPI_OUT_REAL,status,ierr)
          call MPI_Get_count(status, MPI_OUT_REAL, count, ierr)
      enddo
      call MPI_FILE_WRITE(graph_mpi_fh,n_phis*n_thetas*SIZEOF_OUT_REAL,1, &
           &              MPI_INTEGER,status,ierr)
#else
      !PERFON('gwrite_M')
      do n_theta=1,n_thetas

         call MPI_FILE_WRITE(graph_mpi_fh,n_phis*SIZEOF_OUT_REAL,1, &
              &              MPI_INTEGER,status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,dummy(:,n_theta),n_phis, &
              &              MPI_OUT_REAL,status,ierr)
         call MPI_FILE_WRITE(graph_mpi_fh,n_phis*SIZEOF_OUT_REAL,1, &
              &              MPI_INTEGER,status,ierr)

      end do
      !PERFOFF

#endif
   end subroutine graph_write_mpi
#endif
!----------------------------------------------------------------------------
end module graphOut_mod
