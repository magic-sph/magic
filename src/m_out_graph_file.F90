!$Id$
#include "perflib_preproc.cpp"
#include "intrinsic_sizes.h"
#define ONE_LARGE_BLOCK

module graphOut_mod

   use parallel_mod
   use truncation, only: lm_maxMag, n_r_maxMag, n_r_ic_maxMag, lm_max, &
                         n_theta_max, n_phi_tot, n_r_max, l_max, minc, &
                         n_phi_max, nrp, n_r_ic_max
   use radial_data, only: n_r_icb
   use radial_functions, only: r_cmb, orho1, or1, or2, r, r_icb, r_ic, &
                               O_r_ic, O_r_ic2
   use physical_parameters, only: ra, ek, pr, prmag, radratio, sigma_ratio
   use num_param, only: vScale
   use blocking, only: nThetaBs, sizeThetaB
   use horizontal_data, only: theta_ord, dLh, Plm, dPlm, O_sin_theta
   use logic, only: l_mag, l_cond_ic
   use output_data, only: n_graph_file, runid, graph_mpi_fh
#if (FFTLIB==JW)
   use fft_JW
#elif (FFTLIB==MKL)
   use fft_MKL
#endif

   implicit none

   private

#ifdef WITH_MPI
   public :: graphOut, graphOut_mpi, graphOut_IC, graphOut_mpi_header
#else
   public :: graphOut, graphOut_IC
#endif

contains

   subroutine graphOut(time,n_r,format,vr,vt,vp,br,bt,bp,sr, &
     &              n_theta_start,n_theta_block_size,lGraphHeader)
      !-------------------------------------------------------------------------
    
      !  called in amhd
    
      !  output of components of velocity, magnetic field vector and
      !  entropy for graphics. Version May, 23, 2000.
    
      !  n_r: (input) for n_r = 0 a header is written
      !               for n_r > 0 values at radial level n_r are written
    
      !  format : (input) = 0: unformatted
      !                   = 1: formatted writing
      !                   =-1: comment lines are included into file for
      !                        easier reading (cannot be used for graphics
      !                        processing in this form)
    
      !  vr...sr: (input) arrays with grid-point values
    
      !  n_theta_start : (input) values are written for theta-points :
      !                  n_theta_start <= n_theta <= n_theta_start-1+n_theta_block
    
      !-------------------------------------------------------------------------
    
      !-- Input variables
      real(kind=8), intent(in) :: time
      integer,      intent(in) :: n_r                    ! radial grod point no.
      integer,      intent(in) :: format                 ! determines format
      integer,      intent(in) :: n_theta_start          ! start theta no.
      integer,      intent(in) :: n_theta_block_size     ! size of theta block
      real(kind=8), intent(in) :: vr(nrp,*),vt(nrp,*),vp(nrp,*)
      real(kind=8), intent(in) :: br(nrp,*),bt(nrp,*),bp(nrp,*)
      real(kind=8), intent(in) :: sr(nrp,*)

      logical,      intent(inout) :: lGraphHeader
    
      !-- Local variables:
      integer :: n_theta_stop           ! end theta no.
      integer :: n_phi   ! counter for longitude
      integer :: n_theta ! counter for colatitude
      integer :: n_theta_cal                  ! position of block colat in all colats
    
      integer :: prec
    
      real(kind=8) :: fac,fac_r
      real(kind=4) :: dummy(n_phi_max,nfs)
    
      character(len=20) :: version
    
    
      !-- Write header & colatitudes for n_r=0:
    
      if ( lGraphHeader ) then
    
         if ( format /= 0 ) then
    
            !----- Formatted output:
    
            version='Graphout_Version_6'
    
            !-------- Write parameters:
            write(n_graph_file,'(a)')  version ! identifying header
            write(n_graph_file,'(A64)') runid
            if ( format < 0 ) write(n_graph_file,'(/ &
                 &  "Time,truncation,outputgrid:")')
            write(n_graph_file,'(e12.5,1x,10i6)') &
                 time,n_r_max,n_theta_max,n_phi_tot, &
                 n_r_ic_max-1,minc,nThetaBs
            if ( format < 0 ) write(n_graph_file, &
                 & '(/"Ra, Ek, Pr, Pm, radratio, sigma_ratio:")')
            write(n_graph_file,'(6e14.6)') &
                 ra,ek,pr,prmag,radratio,sigma_ratio
    
            !-------- Write colatitudes:
            if ( format < 0 ) then
               write(n_graph_file, &
                    &   '(/"Colatitudes:", /"point#, theta")')
               do n_theta=1,n_theta_max/2   ! NHS
                  write(n_graph_file,'(i4, f9.4)') &
                       n_theta,theta_ord(n_theta)
               end do
               write(n_graph_file,'("--- equator ---")')
               do n_theta=n_theta_max/2+1,n_theta_max   ! SHS
                  write(n_graph_file,'(i4, f9.4)') &
                       n_theta,theta_ord(n_theta)
               end do
            else if ( format >= 0 ) then
               write(n_graph_file,'(129(1x,f8.5))') &
                    (theta_ord(n_theta),n_theta=1,n_theta_max)
            end if
    
         else
    
            !----- Unformatted output:
    
            version='Graphout_Version_7'
    
            !-------- Write parameters:
            write(n_graph_file) version
            write(n_graph_file) runid
            write(n_graph_file)              sngl(time), &
                      float(n_r_max),float(n_theta_max), &
                   float(n_phi_tot),float(n_r_ic_max-1), &
                            float(minc),float(nThetaBs), &
                 sngl(ra),sngl(ek),sngl(pr),sngl(prmag), &
                 sngl(radratio),sngl(sigma_ratio)
    
            !-------- Write colatitudes:
            write(n_graph_file) (sngl(theta_ord(n_theta)), &
                                 n_theta=1,n_theta_max)
    
         end if
    
         lGraphHeader=.false.
    
      else  ! Call not for writing header
    
         !*******************************************************************
         !  define CRITICAL section, so that the write statements do not
         !  get mixed up
         !*******************************************************************
    
         !-- Determine radius and thetas in this block:
         n_theta_stop=n_theta_start+n_theta_block_size-1
         if ( format < 0 ) write(n_graph_file, &
              &   '(/"--- Radial level, radius, i1, i2 ---")')
    
         if ( format /= 0 ) then
            write(n_graph_file,'(I4, F9.5, 2I4)') &
                 n_r-1,r(n_r)/r(1),n_theta_start,n_theta_stop
         else
            write(n_graph_file) float(n_r-1),sngl(r(n_r)/r(1)), &
                                float(n_theta_start),float(n_theta_stop)
         end if
    
    
         !-- Write entropy:
         prec=1
         if ( format < 0 ) write(n_graph_file,'(/" S: ")')
         do n_theta=1,n_theta_block_size,2
            do n_phi=1,n_phi_max ! do loop over phis
               dummy(n_phi,n_theta)  =sngl(sr(n_phi,n_theta))   ! NHS
               dummy(n_phi,n_theta+1)=sngl(sr(n_phi,n_theta+1)) ! SHS
            end do
         end do
         call graph_write(n_phi_max,n_theta_block_size,dummy, &
              prec,format,n_graph_file)
    
         prec=0
    
         !-- Calculate and write radial velocity:
         if ( format < 0 ) write(n_graph_file,'(/" Vr: ")')
         fac=or2(n_r)*vScale*orho1(n_r)
         do n_theta=1,n_theta_block_size,2
            do n_phi=1,n_phi_max
               dummy(n_phi,n_theta)  =sngl(fac*vr(n_phi,n_theta))
               dummy(n_phi,n_theta+1)=sngl(fac*vr(n_phi,n_theta+1))
            end do
         end do
         call graph_write(n_phi_max,n_theta_block_size,dummy, &
              prec,format,n_graph_file)
    
    
         !-- Calculate and write latitudinal velocity:
         if ( format < 0 ) write(n_graph_file,'(/" Vt: ")')
         fac_r=or1(n_r)*vScale*orho1(n_r)
         do n_theta=1,n_theta_block_size,2
            n_theta_cal=n_theta_start+n_theta-1
            fac=fac_r*O_sin_theta(n_theta_cal)
            do n_phi=1,n_phi_max
               dummy(n_phi,n_theta)  =sngl(fac*vt(n_phi,n_theta))
               dummy(n_phi,n_theta+1)=sngl(fac*vt(n_phi,n_theta+1))
            end do
         end do
         call graph_write(n_phi_max,n_theta_block_size,dummy, &
              prec,format,n_graph_file)
    
         !-- Calculate and write longitudinal velocity:
         if ( format < 0 ) write(n_graph_file,'(/" Vp: ")')
         fac_r=or1(n_r)*vScale*orho1(n_r)
         do n_theta=1,n_theta_block_size,2
            n_theta_cal=n_theta_start+n_theta-1
            fac=fac_r*O_sin_theta(n_theta_cal)
            do n_phi=1,n_phi_max
               dummy(n_phi,n_theta)  =sngl(fac*vp(n_phi,n_theta))
               dummy(n_phi,n_theta+1)=sngl(fac*vp(n_phi,n_theta+1))
            end do
         end do
         call graph_write(n_phi_max,n_theta_block_size,dummy, &
              prec,format,n_graph_file)
    
         if ( l_mag ) then
    
            !-- Calculate and write radial magnetic field:
            if ( format < 0 ) write(n_graph_file,'(/'' Br: '')')
            fac=or2(n_r)
            do n_theta=1,n_theta_block_size,2
               do n_phi=1,n_phi_max
                  dummy(n_phi,n_theta)  =sngl(fac*br(n_phi,n_theta))
                  dummy(n_phi,n_theta+1)=sngl(fac*br(n_phi,n_theta+1))
               end do
            end do
            call graph_write(n_phi_max,n_theta_block_size,dummy, &
                 prec,format,n_graph_file)
    
            !-- Calculate and write latitudinal magnetic field:
            if ( format < 0 ) write(n_graph_file,'(/'' Bt: '')')
            do n_theta=1,n_theta_block_size,2
               n_theta_cal=n_theta_start+n_theta-1
               fac=or1(n_r)*O_sin_theta(n_theta_cal)
               do n_phi=1,n_phi_max
                  dummy(n_phi,n_theta)  =sngl(fac*bt(n_phi,n_theta))
                  dummy(n_phi,n_theta+1)=sngl(fac*bt(n_phi,n_theta+1))
               end do
            end do
            call graph_write(n_phi_max,n_theta_block_size,dummy, &
                 prec,format,n_graph_file)
    
            !-- Calculate and write longitudinal magnetic field:
            if ( format < 0 ) write(n_graph_file,'(/'' Bp: '')')
            do n_theta=1,n_theta_block_size,2
               n_theta_cal=n_theta_start+n_theta-1
               fac=or1(n_r)*O_sin_theta(n_theta_cal)
               do n_phi=1,n_phi_max
                  dummy(n_phi,n_theta)  =sngl(fac*bp(n_phi,n_theta))
                  dummy(n_phi,n_theta+1)=sngl(fac*bp(n_phi,n_theta+1))
               end do
            end do
            call graph_write(n_phi_max,n_theta_block_size,dummy, &
                 prec,format,n_graph_file)
    
         end if ! l_mag ?
    
         !-- End of CRITICAL section ********************************************
    
      end if

   end subroutine graphOut
!-----------------------------------------------------------------------
#ifdef WITH_MPI
   subroutine graphOut_mpi(time,n_r,which_form,vr,vt,vp,br,bt,bp,sr, &
            &              n_theta_start,n_theta_block_size,lGraphHeader)

      !-- Input variables:
      real(kind=8), intent(in) :: time
      integer,      intent(in) :: n_r                      ! radial grod point no.
      integer,      intent(in) :: which_form               ! determins format
      integer,      intent(in) :: n_theta_start            ! start theta no.
      integer,      intent(in) :: n_theta_block_size       ! size of theta block
      real(kind=8), intent(in) :: vr(nrp,*),vt(nrp,*),vp(nrp,*)
      real(kind=8), intent(in) :: br(nrp,*),bt(nrp,*),bp(nrp,*)
      real(kind=8), intent(in) :: sr(nrp,*)

      logical, intent(inout) :: lGraphHeader

      !-- Local variables:
      integer :: n_phi         ! counter for longitude
      integer :: n_theta       ! counter for colatitude
      integer :: n_theta_cal   ! position of block colat in all colats
      integer :: n_theta_stop  ! end theta no.

      integer :: prec

      real(kind=8) :: fac,fac_r
      real(kind=4) :: dummy(n_phi_max,nfs)

      character(len=20) :: version

      ! MPI related variables
      !integer :: info
      integer :: status(MPI_STATUS_SIZE)
      !character(len=MPI_MAX_ERROR_STRING) :: error_string
      !integer :: count
      integer :: bytes_written!,length_of_error
      integer :: size_of_header, size_of_data_per_rank, size_of_data_per_r
      integer :: size_of_data_per_thetaB
      integer(kind=MPI_OFFSET_kind) :: disp
      integer :: etype,filetype
      character(len=MPI_MAX_DATAREP_STRING) :: datarep
      ! end of MPI related variables

      !$OMP CRITICAL
      if ( lGraphHeader ) then
         size_of_header = 8+LEN(version)+8+LEN(runid)+8+ &
                          13*SIZEOF_INTEGER+8+n_theta_max*SIZEOF_REAL

#ifdef ONE_LARGE_BLOCK
         size_of_data_per_thetaB = 8+4*SIZEOF_REAL+4* &
                                   (8+n_phi_max*SIZEOF_REAL*n_theta_block_size)
         if ( l_mag ) size_of_data_per_thetaB = size_of_data_per_thetaB + &
                            3*(8+n_phi_max*SIZEOF_REAL*n_theta_block_size)
#else
         size_of_data_per_thetaB = 8+4*SIZEOF_REAL+ &
                                   4*(8+n_phi_max*SIZEOF_REAL)*n_theta_block_size
         if ( l_mag ) size_of_data_per_thetaB = size_of_data_per_thetaB + &
                            3*(8+n_phi_max*SIZEOF_REAL)*n_theta_block_size
#endif
         size_of_data_per_r = size_of_data_per_thetaB * nThetaBs
         size_of_data_per_rank = size_of_data_per_r * nr_per_rank

         if ( rank == 0 ) then
            ! rank zero writes the Header
            disp = 0
            call MPI_FILE_SET_VIEW(graph_mpi_fh,disp,MPI_CHARACTER, &
                                   MPI_CHARACTER,"external32",MPI_INFO_NULL,ierr)
         else
            disp = size_of_header+rank*size_of_data_per_rank
            call MPI_FILE_SET_VIEW(graph_mpi_fh,disp,&
                 & MPI_CHARACTER,MPI_CHARACTER,"external32",MPI_INFO_NULL,ierr)
         end if

         call mpi_file_get_view(graph_mpi_fh,disp,etype,filetype,datarep,ierr)

         bytes_written = 0
         !-- Write header & colatitudes for n_r=0:
         if ( rank == 0 ) then
            if ( which_form /= 0 ) then
#if 0
               !----- Formatted output:
               version='Graphout_Version_6'

               !-------- Write parameters:
               write(n_graph_file,'(a)')  version ! identifying header
               ! for IDL magsym5 routine
               write(n_graph_file,'(A64)') runid
               if ( which_form < 0 ) write(n_graph_file,'(/ &
                    &  "Time,truncation,outputgrid:")')
               write(n_graph_file,'(e12.5,1x,10i6)') &
                    time,n_r_max,n_theta_max,n_phi_tot, &
                    n_r_ic_max-1,minc,nThetaBs
               if ( which_form < 0 ) write(n_graph_file, &
                    & '(/"Ra, Ek, Pr, Pm, radratio, sigma_ratio:")')
               write(n_graph_file,'(6e14.6)') &
                    ra,ek,pr,prmag,radratio,sigma_ratio

               !-------- Write colatitudes:
               if ( which_form < 0 ) then
                  write(n_graph_file, &
                       &   '(/"Colatitudes:", /"point#, theta")')
                  do n_theta=1,n_theta_max/2   ! NHS
                     write(n_graph_file,'(i4, f9.4)') &
                          n_theta,theta_ord(n_theta)
                  end do
                  write(n_graph_file,'("--- equator ---")')
                  do n_theta=n_theta_max/2+1,n_theta_max   ! SHS
                     write(n_graph_file,'(i4, f9.4)') &
                          n_theta,theta_ord(n_theta)
                  end do
               else if ( which_form >= 0 ) then
                  write(n_graph_file,'(129(1x,f8.5))') &
                       (theta_ord(n_theta),n_theta=1,n_theta_max)
               end if
#endif
            else

               !----- Unformatted output:
               version='Graphout_Version_9'

               !-------- Write parameters:
               call MPI_FILE_write(graph_mpi_fh,len(version),1,MPI_INTEGER,status,ierr)
               !call mpi_get_count(status,MPI_INTEGER,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_INTEGER
               call MPI_FILE_write(graph_mpi_fh,version,LEN(version), &
                                   MPI_CHARACTER,status,ierr)
               !call mpi_get_count(status,MPI_CHARACTER,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_CHARACTER
               call MPI_FILE_write(graph_mpi_fh,len(version),1,MPI_INTEGER,status,ierr)
               !call mpi_get_count(status,MPI_INTEGER,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_INTEGER

               call MPI_FILE_write(graph_mpi_fh,LEN(runid),1,MPI_INTEGER,status,ierr)
               !call mpi_get_count(status,MPI_INTEGER,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_INTEGER
               call MPI_FILE_write(graph_mpi_fh,runid,len(runid), &
                                   MPI_CHARACTER,status,ierr)
               !call mpi_get_count(status,MPI_CHARACTER,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_CHARACTER
               call MPI_FILE_write(graph_mpi_fh,len(runid),1,MPI_INTEGER,status,ierr)
               !call mpi_get_count(status,MPI_INTEGER,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_INTEGER

               call MPI_FILE_write(graph_mpi_fh,13*4,1,MPI_INTEGER,status,ierr)
               !call mpi_get_count(status,MPI_INTEGER,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_INTEGER
               call MPI_FILE_write(graph_mpi_fh,sngl(time),1,MPI_REAL,status,ierr)
               !call mpi_get_count(status,MPI_REAL,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_REAL
               call MPI_FILE_write(graph_mpi_fh,float(n_r_max),1,MPI_REAL,status,ierr)
               !call mpi_get_count(status,MPI_REAL,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_REAL
               call MPI_FILE_write(graph_mpi_fh,float(n_theta_max),1, &
                                   MPI_REAL,status,ierr)
               !call mpi_get_count(status,MPI_REAL,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_REAL
               call MPI_FILE_write(graph_mpi_fh,float(n_phi_tot),1, &
                                   MPI_REAL,status,ierr)
               !call mpi_get_count(status,MPI_REAL,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_REAL
               call MPI_FILE_write(graph_mpi_fh,float(n_r_ic_max-1),1, &
                                   MPI_REAL,status,ierr)
               !call mpi_get_count(status,MPI_REAL,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_REAL
               call MPI_FILE_write(graph_mpi_fh,float(minc),1,MPI_REAL,status,ierr)
               !call mpi_get_count(status,MPI_REAL,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_REAL
               call MPI_FILE_write(graph_mpi_fh,float(nThetaBs),1, &
                                   MPI_REAL,status,ierr)
               !call mpi_get_count(status,MPI_REAL,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_REAL
               call MPI_FILE_write(graph_mpi_fh,sngl(ra),1,MPI_REAL,status,ierr)
               !call mpi_get_count(status,MPI_REAL,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_REAL
               call MPI_FILE_write(graph_mpi_fh,sngl(ek),1,MPI_REAL,status,ierr)
               !call mpi_get_count(status,MPI_REAL,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_REAL
               call MPI_FILE_write(graph_mpi_fh,sngl(pr),1,MPI_REAL,status,ierr)
               !call mpi_get_count(status,MPI_REAL,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_REAL
               call MPI_FILE_write(graph_mpi_fh,sngl(prmag),1,MPI_REAL,status,ierr)
               !call mpi_get_count(status,MPI_REAL,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_REAL
               call MPI_FILE_write(graph_mpi_fh,sngl(radratio),1, &
                                   MPI_REAL,status,ierr)
               !call mpi_get_count(status,MPI_REAL,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_REAL
               call MPI_FILE_write(graph_mpi_fh,sngl(sigma_ratio),1, &
                                   MPI_REAL,status,ierr)
               !call mpi_get_count(status,MPI_REAL,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_REAL
               call MPI_FILE_write(graph_mpi_fh,13*4,1,MPI_INTEGER,status,ierr)
               !call mpi_get_count(status,MPI_INTEGER,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_INTEGER

               !-------- Write colatitudes:
               call MPI_FILE_write(graph_mpi_fh,n_theta_max*SIZEOF_REAL,1, &
                                   MPI_INTEGER,status,ierr)
               !call mpi_get_count(status,MPI_INTEGER,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_INTEGER
               do n_theta=1,n_theta_max
                  call MPI_FILE_write(graph_mpi_fh,sngl(theta_ord(n_theta)),1, &
                                      MPI_REAL,status,ierr)
                  !call mpi_get_count(status,MPI_REAL,count,ierr)
                  !bytes_written = bytes_written + count*SIZEOF_REAL
               end do
               call MPI_FILE_write(graph_mpi_fh,n_theta_max*SIZEOF_REAL,1, &
                                   MPI_INTEGER,status,ierr)
               !call mpi_get_count(status,MPI_INTEGER,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_INTEGER

            end if
         end if
         lGraphHeader=.false.
         !PRINT*,"For the header, we wrote ",bytes_written," bytes."
      else  ! Call not for writing header

         !PERFON('mw_data')
         bytes_written=0

         !-- Determine radius and thetas in this block:
         n_theta_stop=n_theta_start+n_theta_block_size-1
#if 0
         if ( which_form < 0 ) write(n_graph_file, &
              &   '(/"--- Radial level, radius, i1, i2 ---")')

         if ( which_form /= 0 ) then
            write(n_graph_file,'(I4, F9.5, 2I4)') &
                 n_r-1, r(n_r)/r(1), n_theta_start, n_theta_stop
         else
#endif
            call MPI_FILE_write(graph_mpi_fh,4*SIZEOF_REAL,1,MPI_INTEGER,status,ierr)
            !call mpi_get_count(status,MPI_INTEGER,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_INTEGER
            call MPI_FILE_write(graph_mpi_fh,float(n_r-1),1,MPI_REAL,status,ierr)
            !call mpi_get_count(status,MPI_REAL,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_REAL
            call MPI_FILE_write(graph_mpi_fh,sngl(r(n_r)/r(1)),1,MPI_REAL,status,ierr)
            !call mpi_get_count(status,MPI_REAL,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_REAL
            call MPI_FILE_write(graph_mpi_fh,float(n_theta_start),1,MPI_REAL,status,ierr)
            !call mpi_get_count(status,MPI_REAL,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_REAL
            call MPI_FILE_write(graph_mpi_fh,float(n_theta_stop),1,MPI_REAL,status,ierr)
            !call mpi_get_count(status,MPI_REAL,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_REAL
            call MPI_FILE_write(graph_mpi_fh,4*SIZEOF_REAL,1,MPI_INTEGER,status,ierr)
            !call mpi_get_count(status,MPI_INTEGER,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_INTEGER
#if 0
         end if
#endif

         !-- Write entropy:
         prec=1
         !if ( which_form < 0 ) write(n_graph_file,'(/" S: ")')
         do n_theta=1,n_theta_block_size,2
            do n_phi=1,n_phi_max ! do loop over phis
               dummy(n_phi,n_theta)  =sngl(sr(n_phi,n_theta))   ! NHS
               dummy(n_phi,n_theta+1)=sngl(sr(n_phi,n_theta+1)) ! SHS
            end do
         end do
         call graph_write_mpi(n_phi_max,n_theta_block_size,dummy, &
              prec,which_form,graph_mpi_fh)

         prec=0

         !-- Calculate and write radial velocity:
         !if ( which_form < 0 ) write(n_graph_file,'(/" Vr: ")')
         fac=or2(n_r)*vScale*orho1(n_r)
         do n_theta=1,n_theta_block_size,2
            do n_phi=1,n_phi_max
               dummy(n_phi,n_theta)  =sngl(fac*vr(n_phi,n_theta))
               dummy(n_phi,n_theta+1)=sngl(fac*vr(n_phi,n_theta+1))
            end do
         end do
         call graph_write_mpi(n_phi_max,n_theta_block_size,dummy, &
              prec,which_form,graph_mpi_fh)


         !-- Calculate and write latitudinal velocity:
         !if ( which_form < 0 ) write(n_graph_file,'(/" Vt: ")')
         fac_r=or1(n_r)*vScale*orho1(n_r)
         do n_theta=1,n_theta_block_size,2
            n_theta_cal=n_theta_start+n_theta-1
            fac=fac_r*O_sin_theta(n_theta_cal)
            do n_phi=1,n_phi_max
               dummy(n_phi,n_theta)  =sngl(fac*vt(n_phi,n_theta))
               dummy(n_phi,n_theta+1)=sngl(fac*vt(n_phi,n_theta+1))
            end do
         end do
         call graph_write_mpi(n_phi_max,n_theta_block_size,dummy, &
              prec,which_form,graph_mpi_fh)

         !-- Calculate and write longitudinal velocity:
         !if ( which_form < 0 ) write(n_graph_file,'(/" Vp: ")')
         fac_r=or1(n_r)*vScale*orho1(n_r)
         do n_theta=1,n_theta_block_size,2
            n_theta_cal=n_theta_start+n_theta-1
            fac=fac_r*O_sin_theta(n_theta_cal)
            do n_phi=1,n_phi_max
               dummy(n_phi,n_theta)  =sngl(fac*vp(n_phi,n_theta))
               dummy(n_phi,n_theta+1)=sngl(fac*vp(n_phi,n_theta+1))
            end do
         end do
         call graph_write_mpi(n_phi_max,n_theta_block_size,dummy, &
              prec,which_form,graph_mpi_fh)

         if ( l_mag ) then

            !-- Calculate and write radial magnetic field:
            !if ( which_form < 0 ) write(n_graph_file,'(/'' Br: '')')
            fac=or2(n_r)
            do n_theta=1,n_theta_block_size,2
               do n_phi=1,n_phi_max
                  dummy(n_phi,n_theta)  =sngl(fac*br(n_phi,n_theta))
                  dummy(n_phi,n_theta+1)=sngl(fac*br(n_phi,n_theta+1))
               end do
            end do
            call graph_write_mpi(n_phi_max,n_theta_block_size,dummy, &
                 prec,which_form,graph_mpi_fh)

            !-- Calculate and write latitudinal magnetic field:
            !if ( which_form < 0 ) write(n_graph_file,'(/'' Bt: '')')
            do n_theta=1,n_theta_block_size,2
               n_theta_cal=n_theta_start+n_theta-1
               fac=or1(n_r)*O_sin_theta(n_theta_cal)
               do n_phi=1,n_phi_max
                  dummy(n_phi,n_theta)  =sngl(fac*bt(n_phi,n_theta))
                  dummy(n_phi,n_theta+1)=sngl(fac*bt(n_phi,n_theta+1))
               end do
            end do
            call graph_write_mpi(n_phi_max,n_theta_block_size,dummy, &
                 prec,which_form,graph_mpi_fh)

            !-- Calculate and write longitudinal magnetic field:
            !if ( which_form < 0 ) write(n_graph_file,'(/'' Bp: '')')
            do n_theta=1,n_theta_block_size,2
               n_theta_cal=n_theta_start+n_theta-1
               fac=or1(n_r)*O_sin_theta(n_theta_cal)
               do n_phi=1,n_phi_max
                  dummy(n_phi,n_theta)  =sngl(fac*bp(n_phi,n_theta))
                  dummy(n_phi,n_theta+1)=sngl(fac*bp(n_phi,n_theta+1))
               end do
            end do
            call graph_write_mpi(n_phi_max,n_theta_block_size,dummy, &
                 prec,which_form,graph_mpi_fh)

         end if ! l_mag ?

         !write(*,"(A,I8)") "bytes_written = ",bytes_written

         !PERFOFF
      end if
      !$OMP END CRITICAL
   end subroutine graphOut_mpi
!----------------------------------------------------------------------------
   subroutine graphOut_mpi_header(time,n_r,which_form, &
        &              n_theta_start,n_theta_block_size)

      !-- Input variables:
      real(kind=8), intent(in) :: time

      integer,      intent(in) :: n_r                    ! radial grod point no.
      integer,      intent(in) :: which_form             ! determines format
      integer,      intent(in) :: n_theta_start          ! start theta no.
      integer,      intent(in) :: n_theta_block_size     ! size of theta block

      !-- Local variables:
      integer :: n_theta       ! counter for colatitude
      character(len=20) :: version

      !-- MPI related variables
      integer :: status(MPI_STATUS_SIZE)
      integer :: bytes_written
      integer :: size_of_header, size_of_data_per_rank, size_of_data_per_r
      integer :: size_of_data_per_thetaB
      integer(kind=MPI_OFFSET_kind) :: disp
      integer :: etype,filetype
      character(len=MPI_MAX_DATAREP_STRING) :: datarep
      ! end of MPI related variables

      size_of_header = 8+LEN(version)+8+LEN(runid)+8+13*SIZEOF_INTEGER+8+ &
                       n_theta_max*SIZEOF_REAL

#ifdef ONE_LARGE_BLOCK
      size_of_data_per_thetaB = 8+4*SIZEOF_REAL+4* &
                                (8+n_phi_max*SIZEOF_REAL*n_theta_block_size)
      if ( l_mag ) size_of_data_per_thetaB = size_of_data_per_thetaB + &
                                           3*(8+n_phi_max*SIZEOF_REAL*n_theta_block_size)
#else
      size_of_data_per_thetaB = 8+4*SIZEOF_REAL+4* &
                               (8+n_phi_max*SIZEOF_REAL)*n_theta_block_size
      if ( l_mag ) size_of_data_per_thetaB = size_of_data_per_thetaB + &
                                           3*(8+n_phi_max*SIZEOF_REAL)*n_theta_block_size
#endif
      size_of_data_per_r = size_of_data_per_thetaB * nThetaBs
      size_of_data_per_rank = size_of_data_per_r * nr_per_rank

      if ( rank == 0 ) then
         ! rank zero writes the Header
         disp = 0
         call MPI_FILE_SET_VIEW(graph_mpi_fh,disp,MPI_CHARACTER, &
                                MPI_CHARACTER,"external32",MPI_INFO_NULL,ierr)
      else
         disp = size_of_header+rank*size_of_data_per_rank
         call MPI_FILE_SET_VIEW(graph_mpi_fh,disp,&
              & MPI_CHARACTER,MPI_CHARACTER,"external32",MPI_INFO_NULL,ierr)
      end if

      call mpi_file_get_view(graph_mpi_fh,disp,etype,filetype,datarep,ierr)

      bytes_written = 0
      !-- Write header & colatitudes for n_r=0:
      if (rank == 0) then
         if ( which_form /= 0 ) then
         else

            !----- Unformatted output:
            version='Graphout_Version_9'

            !-------- Write parameters:
            call MPI_FILE_write(graph_mpi_fh,len(version),1,MPI_INTEGER,status,ierr)
            !call mpi_get_count(status,MPI_INTEGER,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_INTEGER
            call MPI_FILE_write(graph_mpi_fh,version,LEN(version), &
                                MPI_CHARACTER,status,ierr)
            !call mpi_get_count(status,MPI_CHARACTER,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_CHARACTER
            call MPI_FILE_write(graph_mpi_fh,len(version),1,MPI_INTEGER,status,ierr)
            !call mpi_get_count(status,MPI_INTEGER,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_INTEGER


            call MPI_FILE_write(graph_mpi_fh,LEN(runid),1,MPI_INTEGER,status,ierr)
            !call mpi_get_count(status,MPI_INTEGER,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_INTEGER
            call MPI_FILE_write(graph_mpi_fh,runid,len(runid), &
                                MPI_CHARACTER,status,ierr)
            !call mpi_get_count(status,MPI_CHARACTER,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_CHARACTER
            call MPI_FILE_write(graph_mpi_fh,len(runid),1,MPI_INTEGER,status,ierr)
            !call mpi_get_count(status,MPI_INTEGER,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_INTEGER

            call MPI_FILE_write(graph_mpi_fh,13*4,1,MPI_INTEGER,status,ierr)
            !call mpi_get_count(status,MPI_INTEGER,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_INTEGER
            call MPI_FILE_write(graph_mpi_fh,sngl(time),1,MPI_REAL,status,ierr)
            !call mpi_get_count(status,MPI_REAL,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_REAL
            call MPI_FILE_write(graph_mpi_fh,float(n_r_max),1,MPI_REAL,status,ierr)
            !call mpi_get_count(status,MPI_REAL,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_REAL
            call MPI_FILE_write(graph_mpi_fh,float(n_theta_max),1, &
                                MPI_REAL,status,ierr)
            !call mpi_get_count(status,MPI_REAL,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_REAL
            call MPI_FILE_write(graph_mpi_fh,float(n_phi_tot),1,MPI_REAL,status,ierr)
            !call mpi_get_count(status,MPI_REAL,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_REAL
            call MPI_FILE_write(graph_mpi_fh,float(n_r_ic_max-1),1, &
                                MPI_REAL,status,ierr)
            !call mpi_get_count(status,MPI_REAL,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_REAL
            call MPI_FILE_write(graph_mpi_fh,float(minc),1,MPI_REAL,status,ierr)
            !call mpi_get_count(status,MPI_REAL,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_REAL
            call MPI_FILE_write(graph_mpi_fh,float(nThetaBs),1,MPI_REAL,status,ierr)
            !call mpi_get_count(status,MPI_REAL,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_REAL
            call MPI_FILE_write(graph_mpi_fh,sngl(ra),1,MPI_REAL,status,ierr)
            !call mpi_get_count(status,MPI_REAL,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_REAL
            call MPI_FILE_write(graph_mpi_fh,sngl(ek),1,MPI_REAL,status,ierr)
            !call mpi_get_count(status,MPI_REAL,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_REAL
            call MPI_FILE_write(graph_mpi_fh,sngl(pr),1,MPI_REAL,status,ierr)
            !call mpi_get_count(status,MPI_REAL,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_REAL
            call MPI_FILE_write(graph_mpi_fh,sngl(prmag),1,MPI_REAL,status,ierr)
            !call mpi_get_count(status,MPI_REAL,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_REAL
            call MPI_FILE_write(graph_mpi_fh,sngl(radratio),1,MPI_REAL,status,ierr)
            !call mpi_get_count(status,MPI_REAL,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_REAL
            call MPI_FILE_write(graph_mpi_fh,sngl(sigma_ratio),1,MPI_REAL,status,ierr)
            !call mpi_get_count(status,MPI_REAL,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_REAL
            call MPI_FILE_write(graph_mpi_fh,13*4,1,MPI_INTEGER,status,ierr)
            !call mpi_get_count(status,MPI_INTEGER,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_INTEGER

            !-------- Write colatitudes:
            call MPI_FILE_write(graph_mpi_fh,n_theta_max*SIZEOF_REAL,1, &
                                MPI_INTEGER,status,ierr)
            !call mpi_get_count(status,MPI_INTEGER,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_INTEGER
            do n_theta=1,n_theta_max
               call MPI_FILE_write(graph_mpi_fh,sngl(theta_ord(n_theta)),1, &
                                   MPI_REAL,status,ierr)
               !call mpi_get_count(status,MPI_REAL,count,ierr)
               !bytes_written = bytes_written + count*SIZEOF_REAL
            end do
            call MPI_FILE_write(graph_mpi_fh,n_theta_max*SIZEOF_REAL,1, &
                                MPI_INTEGER,status,ierr)
            !call mpi_get_count(status,MPI_INTEGER,count,ierr)
            !bytes_written = bytes_written + count*SIZEOF_INTEGER

         end if
      end if

   end subroutine graphOut_mpi_header
#endif
!----------------------------------------------------------------------------
   subroutine graphOut_IC(format,b_ic,db_ic,ddb_ic,aj_ic,dj_ic,b)
      !  +-------------------------------------------------------------------+
      !  |                                                                   |
      !  |  Purpose of this subroutine is to write inner core magnetic       |
      !  |  field onto graphic output file. If the inner core is             |
      !  |  insulating (l_cond_ic=false) the potential field is calculated   |
      !  |  from the outer core field at r=r_cmb.                            |
      !  |  This version assumes that the fields are fully local on the rank |
      !  |  which is calling this routine (usually rank 0).                  |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+

      !-- Input variables:
      integer,         intent(in) :: format    ! controls output format
      complex(kind=8), intent(in) :: b_ic(lm_maxMag,n_r_ic_maxMag)
      complex(kind=8), intent(in) :: db_ic(lm_maxMag,n_r_ic_maxMag)
      complex(kind=8), intent(in) :: ddb_ic(lm_maxMag,n_r_ic_maxMag)
      complex(kind=8), intent(in) :: aj_ic(lm_maxMag,n_r_ic_maxMag)
      complex(kind=8), intent(in) :: dj_ic(lm_maxMag,n_r_ic_maxMag)
      complex(kind=8), intent(in) :: b(lm_maxMag,n_r_maxMag)
    
      !-- Local variables:
      integer :: nR
      integer :: nThetaB,nTheta,nThetaStart,nThetaC
      integer :: nPhi
    
      complex(kind=8) :: dLhb(lm_max)
      complex(kind=8) :: bhG(lm_max)
      complex(kind=8) :: bhC(lm_max)
      complex(kind=8) :: dLhj(lm_max)
      complex(kind=8) :: cbhG(lm_max)
      complex(kind=8) :: cbhC(lm_max)
      real(kind=8) :: BrB(nrp,nfs)
      real(kind=8) :: BtB(nrp,nfs)
      real(kind=8) :: BpB(nrp,nfs)
      real(kind=4) :: Br(n_phi_max,n_theta_max)
      real(kind=4) :: Bt(n_phi_max,n_theta_max)
      real(kind=4) :: Bp(n_phi_max,n_theta_max)
    
    
      integer :: prec  ! controls output precision in graph_write
    
#ifdef WITH_MPI
      ! MPI specific variables
      integer :: status(MPI_STATUS_SIZE)
      ! end MPI variables
#endif

      prec=1
    
      !-- Loop over all radial levels:
    
      do nR=2,n_r_ic_max  ! nR=1 is ICB
    
         if ( l_cond_ic ) then
            call legPrep_IC(b_ic(1,nR),db_ic(1,nR),ddb_ic(1,nR), &
                 aj_ic(1,nR),dj_ic(1,nR), &
                 dLh,lm_max,l_max,minc, &
                 r_ic(nR),r_ICB,.false.,.true.,l_cond_ic, &
                 dLhb,bhG,bhC,dLhj,cbhG,cbhC)
         else
            call legPrep_IC(b(1,n_r_icb),db_ic(1,1), &
                 ddb_ic(1,1),aj_ic(1,1), &
                 dj_ic(1,1),dLh,lm_max,l_max,minc, &
                 r_ic(nR),r_ICB,.false.,.true.,l_cond_ic, &
                 dLhb,bhG,bhC,dLhj,cbhG,cbhC)
         end if
    
         do nThetaB=1,nThetaBs
            nThetaStart=(nThetaB-1)*sizeThetaB+1
    
            !------ Preform Legendre transform:
            call legTF(dLhb,bhG,bhC,dLhj,cbhG,cbhC, &
                 l_max,minc,nThetaStart,sizeThetaB, &
                 Plm,dPlm,.true.,.false., &
                 BrB,BtB,BpB,BrB,BrB,BrB)
    
            call fft_thetab(BrB,1)
            call fft_thetab(BtB,1)
            call fft_thetab(BpB,1)
    
            !------ Copy theta block and calculate real components:
            do nTheta=1,sizeThetaB
               nThetaC=nThetaStart-1+nTheta
               do nPhi=1,n_phi_max
                  Br(nPhi,nThetaC)=sngl(BrB(nPhi,nTheta)*O_r_ic2(nR))
                  Bt(nPhi,nThetaC)=sngl(BtB(nPhi,nTheta)*O_r_ic(nR) * &
                       O_sin_theta(nThetaC))
                  Bp(nPhi,nThetaC)=sngl(BpB(nPhi,nTheta)*O_r_ic(nR) * &
                       O_sin_theta(nThetaC))
               end do
            end do
    
         end do
    
    
         !-- Write radius and theta information (all thetas at one go)
         if ( format == -1 ) write(n_graph_file, &
              &   '(/"--- Radial level IC, radius, i1, i2 ---")')
    
         if ( format /= 0 ) then
            write(n_graph_file,'(I4, F9.5, 2I4)') &
                 n_r_max+nR-2,r_ic(nR)/r_cmb, &
                 1,n_theta_max
         else if ( format == 0 ) then
#ifdef WITH_MPI
            ! in process n_procs-1 the last oc fields have been written,
            ! Now just append on this process.
            if (rank == n_procs-1) then
               call MPI_FILE_write(graph_mpi_fh,4*4,1,MPI_INTEGER,status,ierr)
               call MPI_FILE_write(graph_mpi_fh,float(n_r_max+nR-2),1, &
                                   MPI_REAL,status,ierr)
               call MPI_FILE_write(graph_mpi_fh,sngl(r_ic(nR)/r_cmb),1, &
                                   MPI_REAL,status,ierr)
               call MPI_FILE_write(graph_mpi_fh,1.E0,1,MPI_REAL,status,ierr)
               call MPI_FILE_write(graph_mpi_fh,float(n_theta_max),1, &
                                   MPI_REAL,status,ierr)
               call MPI_FILE_write(graph_mpi_fh,4*4,1,MPI_INTEGER,status,ierr)
            end if
#else
            write(n_graph_file) float(n_r_max+nR-2),sngl(r_ic(nR)/r_cmb), &
                 &              1.E0,float(n_theta_max)
#endif
         end if


         !-- Write radial magnetic field:
         if ( format == -1 ) write(n_graph_file,'(/'' Br IC: '')')
#ifdef WITH_MPI
         if (rank == n_procs-1) then
            call graph_write_mpi(n_phi_max,n_theta_max,Br, &
                 prec,FORMAT,graph_mpi_fh)
         end if
#else
         call graph_write(n_phi_max,n_theta_max,Br, &
              prec,format,n_graph_file)
#endif

         !-- Write latitudinal magnetic field:
         if ( format == -1 ) write(n_graph_file,'(/'' Bt IC: '')')
#ifdef WITH_MPI
         if (rank == n_procs-1) then
            call graph_write_mpi(n_phi_max,n_theta_max,Bt, &
                 prec,format,graph_mpi_fh)
         end if
#else
         call graph_write(n_phi_max,n_theta_max,Bt, &
              prec,format,n_graph_file)
#endif
  
         !-- Write longitudinal magnetic field:
         if ( format == -1 ) write(n_graph_file,'(/'' Bp IC: '')')
#ifdef WITH_MPI
         if (rank == n_procs-1) then
            call graph_write_mpi(n_phi_max,n_theta_max,Bp, &
                 prec,format,graph_mpi_fh)
         end if
#else
         call graph_write(n_phi_max,n_theta_max,Bp, &
                          prec,format,n_graph_file)
#endif

      end do  ! Do loop over radial levels nR


   end subroutine graphOut_IC
!-----------------------------------------------------------------------
   subroutine graph_write(n_phis,n_thetas,dummy,prec,FORMAT,n_graph_file)
      !------------------------------------------------------------------------------
      !  This subroutine writes the data for one theta-band
      !  (stored in 'dummy'). Version May, 5, 2000.
      !------------------------------------------------------------------------------

      !-- Input variables:
      integer,      intent(in) :: n_thetas            ! number of first colatitude value
      integer,      intent(in) :: n_phis              ! number of logitudes to be printed
      real(kind=4), intent(in) :: dummy(n_phi_max,*)  ! data
      integer,      intent(in) :: prec                ! determines precision if output
      integer,      intent(in) :: format              ! formatted/unformatted output
      integer,      intent(in) :: n_graph_file        ! output unit

      !-- Local variables:
      integer :: n_phi,n_theta


      !PERFON('gwrite')
      do n_theta=1,n_thetas
           
         if ( format == 0 ) then ! unformatted output
            write(n_graph_file) (dummy(n_phi,n_theta),n_phi=1,n_phis)
         else                ! formatted output
            if ( prec == 0 ) then
               write(n_graph_file,900) &
                    (dummy(n_phi,n_theta),n_phi=1,n_phis)
            else if ( prec == 1 ) then
               write(n_graph_file,901) &
                    (dummy(n_phi,n_theta),n_phi=1,n_phis)
            end if
         end if

      end do

      900 FORMAT(512(1X,F7.2))
      901 FORMAT(512(1X,F7.3))

      !PERFOFF

   end subroutine graph_write
!------------------------------------------------------------------------------
#ifdef WITH_MPI
   subroutine graph_write_mpi(n_phis,n_thetas,dummy,prec,which_form,graph_mpi_fh)

      !-- Input variables
      integer,      intent(in) :: n_thetas          ! number of first colatitude value
      integer,      intent(in) :: n_phis            ! number of logitudes to be printed
      real(kind=4), intent(in) :: dummy(n_phi_max,*)! data
      integer,      intent(in) :: prec              ! determines precision if output
      integer,      intent(in) :: which_form        ! formatted/unformatted output
      integer,      intent(in) :: graph_mpi_fh      ! mpi handle of the mpi file

      !-- Local variables:
      integer :: n_phi,n_theta

      !-- MPI related variables
      integer :: status(MPI_STATUS_SIZE)!,count

#ifdef ONE_LARGE_BLOCK
      call MPI_FILE_write(graph_mpi_fh,n_phis*n_thetas*SIZEOF_REAL,1, &
                          MPI_INTEGER,status,ierr)
      call MPI_FILE_write(graph_mpi_fh,dummy(:,1:n_thetas),n_phis*n_thetas, &
                          MPI_REAL,status,ierr)
      call MPI_FILE_write(graph_mpi_fh,n_phis*n_thetas*SIZEOF_REAL,1, &
                          MPI_INTEGER,status,ierr)
#else
      !PERFON('gwrite_M')
      do n_theta=1,n_thetas

         if ( which_form == 0 ) then ! unformatted output
            call MPI_FILE_write(graph_mpi_fh,n_phis*SIZEOF_REAL,1, &
                                MPI_INTEGER,status,ierr)
            call MPI_FILE_write(graph_mpi_fh,dummy(1,n_theta),n_phis, &
                                MPI_REAL,status,ierr)
            call MPI_FILE_write(graph_mpi_fh,n_phis*SIZEOF_REAL,1, &
                                MPI_INTEGER,status,ierr)

            !call MPI_FILE_write(n_graph_file) &
            !      (dummy(n_phi,n_theta),n_phi=1,n_phis)
         else                ! formatted output
#if 0  
            if ( prec == 0 ) then
               write(n_graph_file,900) &
                    (dummy(n_phi,n_theta),n_phi=1,n_phis)
            else if ( prec == 1 ) then
               write(n_graph_file,901) &
                    (dummy(n_phi,n_theta),n_phi=1,n_phis)
            end if
#endif
         end if

      end do

      900 FORMAT(512(1X,F7.2))
      901 FORMAT(512(1X,F7.3))

    !PERFOFF
#endif
   end subroutine graph_write_mpi
#endif
!----------------------------------------------------------------------------
end module graphOut_mod
