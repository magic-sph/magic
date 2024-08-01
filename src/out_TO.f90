module outTO_mod
   !
   ! This module handles the writing of TO-related outputs: zonal force balance and
   ! z-integrated terms. This is a re-implementation of the spectral method used up
   ! to MagIC 5.10, which formerly relies on calculation of Plm on the cylindrical grid.
   ! This was quite costly and not portable on large truncations. As such, a local
   ! 4th order method is preferred here.
   !
   use precision_mod
   use parallel_mod
   use constants, only: one, two, pi, vol_oc
   use truncation, only: n_r_max, n_theta_max, n_phi_max, minc
   use mem_alloc, only: bytes_allocated
   use num_param, only: tScale
   use output_data, only: sDens, zDens, tag, runid, log_file, n_log_file
   use radial_data, only: radial_balance, nRstart, nRstop
   use radial_functions, only: r_ICB, r_CMB, r
   use horizontal_data, only: theta_ord, phi
   use logic, only: l_save_out, l_TOmovie, l_full_sphere, l_phase_field, l_mag
   use physical_parameters, only: ra, ek, pr,prmag, radratio, LFfac
   use torsional_oscillations, only: dzCorAS_Rloc, dzdVpAS_Rloc, dzddVpAS_Rloc,  &
       &                             dzRstrAS_Rloc, dzAstrAS_Rloc, dzStrAS_Rloc, &
       &                             dzLFAS_Rloc, V2AS_Rloc, Bs2AS_Rloc,         &
       &                             BspdAS_Rloc, BpsdAS_Rloc, BzpdAS_Rloc,      &
       &                             BpzdAS_Rloc, BpzAS_Rloc, BspAS_Rloc,        &
       &                             BszAS_Rloc, VAS_Rloc, dzPenAS_Rloc
   use useful, only: logWrite
   use integration, only: cylmean_otc, cylmean_itc, simps

   implicit none

   private

   real(cp), allocatable :: cyl(:) ! Cylindrical grid
   real(cp), allocatable :: h(:)   ! height
   real(cp), allocatable :: Oh(:)  ! 1/h
   real(cp), allocatable :: dzCorAS(:,:), dzdVpAS(:,:), dzddVpAS(:,:), BszAS(:,:)
   real(cp), allocatable :: dzRstrAS(:,:), dzAstrAS(:,:), dzStrAS(:,:), BpzAS(:,:)
   real(cp), allocatable :: dzLFAS(:,:), V2AS(:,:), Bs2AS(:,:), BspAS(:,:), VAS(:,:)
   real(cp), allocatable :: BspdAS(:,:), BpsdAS(:,:), BzpdAS(:,:), BpzdAS(:,:)
   real(cp), allocatable :: dzPenAS(:,:)
   integer :: n_s_otc, n_s_max, n_NHS_file, n_SHS_file, n_TOmov_file, n_TO_file
   real(cp) :: volcyl_oc
   character(len=64) :: movFile, TOFile
   integer :: nTOmovSets ! Number of TO_mov frames

   public :: initialize_outTO_mod, finalize_outTO_mod, outTO


contains

   subroutine initialize_outTO_mod()
      !
      ! Memory allocation of arrays needed for TO outputs
      !

      !-- Local variables
      integer :: n_s, nFields, n
      real(cp) :: smin, smax, ds, dumm(12)
      character(len=64) :: TOfileNhs,TOfileShs
      character(len=64) :: version

      if ( rank == 0 ) then
         allocate( dzCorAS(n_theta_max,n_r_max), dzdVpAS(n_theta_max,n_r_max) )
         allocate( dzddVpAS(n_theta_max,n_r_max), dzRstrAS(n_theta_max,n_r_max) )
         allocate( dzAstrAS(n_theta_max,n_r_max), dzStrAS(n_theta_max,n_r_max) )
         allocate( V2AS(n_theta_max,n_r_max), VAS(n_theta_max,n_r_max) )
         bytes_allocated=bytes_allocated+8*n_theta_max*n_r_max*SIZEOF_DEF_REAL
         if ( l_mag ) then
            allocate( dzLFAS(n_theta_max,n_r_max) )
            allocate( Bs2AS(n_theta_max,n_r_max), BspAS(n_theta_max,n_r_max) )
            allocate( BszAS(n_theta_max,n_r_max), BpzAS(n_theta_max,n_r_max) )
            allocate( BspdAS(n_theta_max,n_r_max), BpsdAS(n_theta_max,n_r_max) )
            allocate( BzpdAS(n_theta_max,n_r_max), BpzdAS(n_theta_max,n_r_max) )
            bytes_allocated=bytes_allocated+9*n_theta_max*n_r_max*SIZEOF_DEF_REAL
         end if
         if ( l_phase_field ) then
            allocate(dzPenAS(n_theta_max,n_r_max))
            bytes_allocated=bytes_allocated+n_theta_max*n_r_max*SIZEOF_DEF_REAL
         end if
      else
         allocate( dzCorAS(1,1), dzdVpAS(1,1), dzddVpAS(1,1), dzRstrAS(1,1) )
         allocate( dzAstrAS(1,1), dzStrAS(1,1), VAS(1,1), V2AS(1,1) )
         if ( l_mag ) then
            allocate( BspAS(1,1), BpzAS(1,1), BszAS(1,1), Bs2AS(1,1), dzLFAS(1,1) )
            allocate( BspdAS(1,1), BpsdAS(1,1), BzpdAS(1,1), BpzdAS(1,1) )
         end if
         if ( l_phase_field ) allocate(dzPenAS(1,1))
      end if

      !-- Cylindrical radius
      n_s_max = n_r_max+int(r_ICB*n_r_max)
      n_s_max = int(sDens*n_s_max)
      allocate( cyl(n_s_max) )
      bytes_allocated=bytes_allocated+n_s_max*SIZEOF_DEF_REAL

      ! Minimum and maximum cylindrical radii
      smin = r_CMB*sin(theta_ord(1))
      !smin = 0.0_cp
      smax = r_CMB

      !-- Grid spacing
      ds = (smax-smin)/(n_s_max-1)

      !-- Cylindrical grid
      do n_s=1,n_s_max
         cyl(n_s)=r_cmb-(n_s-1)*ds
      end do

      !-- Height
      allocate( h(n_s_max), Oh(n_s_max) )
      bytes_allocated=bytes_allocated+2*n_s_max*SIZEOF_DEF_REAL
      do n_s=1,n_s_max
         if ( cyl(n_s) >= r_ICB ) then
            h(n_s)=two*sqrt(r_CMB**2-cyl(n_s)**2)
            n_s_otc=n_s ! Last index outside TC
         else
            h(n_s)=sqrt(r_CMB**2-cyl(n_s)**2)-sqrt(r_ICB**2-cyl(n_s)**2)
         end if
      end do

      Oh(1)=0.0_cp
      Oh(2:n_s_max)=one/h(2:n_s_max)

      if ( rank == 0 ) then
         !-- Handle file openings and headers
         TOfileNhs='TOnhs.'//tag
         TOfileShs='TOshs.'//tag
         TOFile   ='Tay.'//tag

         open(newunit=n_NHS_file, file=TOfileNhs, status='unknown', form='unformatted', &
         &    position='append')
         write(n_NHS_file) real(n_s_max,kind=outp)                    ! 1
         write(n_NHS_file) (real(cyl(n_s),kind=outp),n_s=1,n_s_max)   ! 2

         open(newunit=n_SHS_file, file=TOfileShs, status='unknown', form='unformatted', &
         &    position='append')
         write(n_SHS_file) real(n_s_max,kind=outp)                    ! 1
         write(n_SHS_file) (real(cyl(n_s),kind=outp),n_s=1,n_s_max)   ! 2

         if (.not. l_save_out ) open(newunit=n_TO_file, file=TOFile, status='new')
      end if

      !-- Determine the volume of the spherical shell interpolated on the cylindrical grid
      if ( .not. l_full_sphere ) then
         volcyl_oc = two * pi * (simps(h*cyl, cyl)+simps(h(n_s_otc+1:)*cyl(n_s_otc+1:), &
         &           cyl(n_s_otc+1:)))
      else
         volcyl_oc = two * pi * simps(h*cyl, cyl)
      end if

      !-- TO movie related outputs
      if ( l_TOmovie ) then
         nTOmovSets=0
         movFile  ='TO_mov.'//tag

         !-- Write the header of the TO_mov.TAG file
         if ( rank == 0 ) then
            open(newunit=n_TOmov_file, file=movFile, status='unknown',  &
            &    form='unformatted', position='append')

            nFields=7
            if ( l_phase_field ) nFields=nFields+1
            version='JW_Movie_Version_2'
            write(n_TOmov_file) version
            dumm(1)=102           ! type of input
            dumm(2)=3             ! marker for constant phi plane
            dumm(3)=0.0_cp        ! surface constant
            dumm(4)=nFields       ! no of fields
            write(n_TOmov_file) (real(dumm(n),kind=outp),n=1,4)

            dumm(1)=11.0          ! Field marker for AS vPhi
            dumm(2)=61.0          ! Field marker for Reynolds Force
            dumm(3)=62.0          ! Field marker for Advective Force
            dumm(4)=63.0          ! Field marker for Viscous Force
            dumm(5)=64.0          ! Field marker for Lorentz Force
            dumm(6)=65.0          ! Field marker for Coriolis force
            dumm(7)=66.0          ! Field marker for dtVp
            write(n_TOmov_file) (real(dumm(n),kind=outp),n=1,nFields)

            !------ Now other info about grid and parameters:
            write(n_TOmov_file) runid ! run identifier
            dumm( 1)=n_r_max          ! total number of radial points
            dumm( 2)=n_r_max          ! no of radial point in outer core
            dumm( 3)=n_theta_max      ! no. of theta points
            dumm( 4)=n_phi_max        ! no. of phi points
            dumm( 5)=minc             ! imposed symmetry
            dumm( 6)=ra               ! control parameters
            dumm( 7)=ek               ! (for information only)
            dumm( 8)=pr               !      -"-
            dumm( 9)=prmag            !      -"-
            dumm(10)=radratio         ! ratio of inner / outer core
            dumm(11)=tScale           ! timescale
            write(n_TOmov_file) (real(dumm(n),kind=outp),     n=1,11)
            write(n_TOmov_file) (real(r(n)/r_CMB,kind=outp),  n=1,n_r_max)
            write(n_TOmov_file) (real(theta_ord(n),kind=outp),n=1,n_theta_max)
            write(n_TOmov_file) (real(phi(n),kind=outp),      n=1,n_phi_max)
         end if
      end if

   end subroutine initialize_outTO_mod
!------------------------------------------------------------------------------
   subroutine finalize_outTO_mod()
      !
      ! Memory de-allocation of arrays needed for TO outputs
      !

      !-- Close files
      if ( rank == 0 ) then
         close(n_NHS_file)
         close(n_SHS_file)
         if ( l_TOmovie ) close(n_TOmov_file)
         if ( .not. l_save_out ) close(n_TO_file)
      end if

      !-- Deallocate arrays
      deallocate(cyl, h, Oh, VAS)
      deallocate(dzCorAS, dzdVpAS, dzddVpAS, dzRstrAS, dzAstrAS, dzStrAS, V2AS)
      if ( l_mag ) then
         deallocate(dzLFAS, BspAS, BpzAS, BszAS, Bs2AS, BspdAS)
         deallocate(BpsdAS, BzpdAS, BpzdAS)
      end if
      if ( l_phase_field ) deallocate(dzPenAS)

   end subroutine finalize_outTO_mod
!------------------------------------------------------------------------------
   subroutine outTO(time,n_time_step,eKin,eKinTAS,lTOmov)
      !
      !   Output of axisymmetric zonal flow, its relative strength,
      !   its time variation, and all forces acting on it.
      !

      !-- Input variables
      real(cp), intent(in) :: time ! Time
      real(cp), intent(in) :: eKin ! Kinetic energy
      real(cp), intent(in) :: eKinTAS ! Toroidal axisymmetric energy
      integer,  intent(in) :: n_time_step ! Iteration number
      logical,  intent(in) :: lTOmov ! Do we need to store the movie files as well

      !-- Local variables
      character(len=255) :: message
      logical :: lTC
      integer :: n_s, n
      real(cp) :: V2IntS(n_s_max),V2IntN(n_s_max),Bs2IntS(n_s_max),Bs2IntN(n_s_max)
      real(cp) :: BspIntS(n_s_max),BspIntN(n_s_max),BspdIntS(n_s_max),BspdIntN(n_s_max)
      real(cp) :: BpsdIntS(n_s_max),BpsdIntN(n_s_max),dVpIntS(n_s_max),dVpIntN(n_s_max)
      real(cp) :: ddVpIntS(n_s_max),ddVpIntN(n_s_max),LFIntS(n_s_max),LFIntN(n_s_max)
      real(cp) :: TayIntS(n_s_max),TayIntN(n_s_max),TayVIntS(n_s_max),TayVIntN(n_s_max)
      real(cp) :: TayRIntS(n_s_max),TayRIntN(n_s_max),RstrIntS(n_s_max),RstrIntN(n_s_max)
      real(cp) :: AstrIntS(n_s_max),AstrIntN(n_s_max),StrIntS(n_s_max),StrIntN(n_s_max)
      real(cp) :: VpIntN(n_s_max),VpIntS(n_s_max),VpRIntN(n_s_max),VpRIntS(n_s_max)
      real(cp) :: TauBS(n_s_max),dTauBS(n_s_max),dTTauBS(n_s_max),SBspIntN(n_s_max)
      real(cp) :: TauBN(n_s_max),dTauBN(n_s_max),dTTauBN(n_s_max),SBspIntS(n_s_max)
      real(cp) :: SVpIntN(n_s_max),SVpIntS(n_s_max),SBs2IntN(n_s_max),SBs2IntS(n_s_max)
      real(cp) :: dSVpIntN(n_s_max),dSVpIntS(n_s_max),d2SVpIntN(n_s_max)
      real(cp) :: d2SVpIntS(n_s_max),dSBspIntN(n_s_max),dSBspIntS(n_s_max)
      real(cp) :: dSBs2IntN(n_s_max),dSBs2IntS(n_s_max),TauN(n_s_max),TauS(n_s_max)
      real(cp) :: dTTauS(n_s_max),dTTauN(n_s_max)
      real(cp) :: VgRMS, VpRMS, VRMS, TayRMS, TayRRMS, TayVRMS

      !-- For boundaries:
      real(cp) :: BspB(2),BspdB(2),BpsdB(2),zMin,zMax,dumm(8)
      real(cp) :: Bs2B(2),BszB(2),BpzB(2),BzpdB(2),BpzdB(2)

      !-- Gather R-distributed arrays on rank == 0
      call gather_from_Rloc_to_rank0(dzCorAS_Rloc, dzCorAS)
      if ( l_phase_field ) call gather_from_Rloc_to_rank0(dzPenAS_Rloc, dzPenAS)
      call gather_from_Rloc_to_rank0(dzdVpAS_Rloc, dzdVpAS)
      call gather_from_Rloc_to_rank0(dzddVpAS_Rloc, dzddVpAS)
      call gather_from_Rloc_to_rank0(dzRstrAS_Rloc, dzRstrAS)
      call gather_from_Rloc_to_rank0(dzAstrAS_Rloc, dzAstrAS)
      call gather_from_Rloc_to_rank0(dzStrAS_Rloc, dzStrAS)

      call gather_from_Rloc_to_rank0(VAS_Rloc, VAS)
      call gather_from_Rloc_to_rank0(V2AS_Rloc, V2AS)

      if ( l_mag ) then
         call gather_from_Rloc_to_rank0(dzLFAS_Rloc, dzLFAS)
         call gather_from_Rloc_to_rank0(Bs2AS_Rloc, Bs2AS)
         call gather_from_Rloc_to_rank0(BszAS_Rloc, BszAS)
         call gather_from_Rloc_to_rank0(BspAS_Rloc, BspAS)
         call gather_from_Rloc_to_rank0(BpzAS_Rloc, BpzAS)
         call gather_from_Rloc_to_rank0(BspdAS_Rloc, BspdAS)
         call gather_from_Rloc_to_rank0(BpsdAS_Rloc, BpsdAS)
         call gather_from_Rloc_to_rank0(BzpdAS_Rloc, BzpdAS)
         call gather_from_Rloc_to_rank0(BpzdAS_Rloc, BpzdAS)
      end if

      !-- Starting from here only rank=0 will do the cylindrical integration
      if ( rank == 0 ) then
         !-- Z-integrals
         call cylmean(VAS, VpIntN, VpIntS)
         call cylmean(dzdVpAS, dVpIntN, dVpIntS)
         call cylmean(dzddVpAS, ddVpIntN, ddVpIntS)
         call cylmean(abs(dzRstrAS), TayRIntN, TayRIntS)
         call cylmean(abs(dzStrAS), TayVIntN, TayVIntS)
         call cylmean(dzRstrAS, RstrIntN, RstrIntS)
         call cylmean(dzAstrAS, AstrIntN, AstrIntS)
         call cylmean(dzStrAS, StrIntN, StrIntS)
         call cylmean(V2AS, V2IntN, V2IntS)
         if ( l_mag ) then
            call cylmean(dzLFAS, LFIntN, LFIntS)
            call cylmean(abs(dzLFAS), TayIntN, TayIntS)
            call cylmean(Bs2AS, Bs2IntN, Bs2IntS)
            call cylmean(BspAS, BspIntN, BspIntS)
            call cylmean(BspdAS, BspdIntN, BspdIntS)
            call cylmean(BpsdAS, BpsdIntN, BpsdIntS)
         else
            LFIntN(:)  =0.0_cp
            LFIntS(:)  =0.0_cp
            TayIntN(:) =0.0_cp
            TayIntS(:) =0.0_cp
            Bs2IntN(:) =0.0_cp
            Bs2IntS(:) =0.0_cp
            BspIntN(:) =0.0_cp
            BspIntS(:) =0.0_cp
            BspdIntN(:)=0.0_cp
            BspdIntS(:)=0.0_cp
            BpsdIntN(:)=0.0_cp
            BpsdIntS(:)=0.0_cp
            TauBN(:)   =0.0_cp
            TauBS(:)   =0.0_cp
            dTauBN(:)  =0.0_cp
            dTauBS(:)  =0.0_cp
            dTTauBN(:) =0.0_cp
            dTTauBS(:) =0.0_cp
         end if

         do n_s=1,n_s_max
            if ( cyl(n_s) < r_ICB ) then
               lTC = .true.
            else
               lTC = .false.
            end if

            if ( V2IntN(n_s) < 0.0_cp ) then
               VpRIntN(n_s)=one
            else if ( abs(V2IntN(n_s)) <= 10.0_cp*epsilon(one) ) then
               VpRIntN(n_s)=0.0_cp
            else
               VpRIntN(n_s)=abs(VpIntN(n_s))/sqrt(V2IntN(n_s))
            end if
            if ( V2IntS(n_s) < 0.0_cp ) then
               VpRIntS(n_s)=one
            else if ( abs(V2IntS(n_s)) <= 10.0_cp*epsilon(one) ) then
               VpRIntS(n_s)=0.0_cp
            else
               VpRIntS(n_s)=abs(VpIntS(n_s))/sqrt(V2IntS(n_s))
            end if
            VpRIntN(n_s)=min(one,VpRIntN(n_s))
            VpRIntS(n_s)=min(one,VpRIntS(n_s))
            if ( abs(TayIntN(n_s)) > 0.0_cp ) then
               TayIntN(n_s) =LFIntN(n_s)/TayIntN(n_s)
               TayIntS(n_s) =LFIntS(n_s)/TayIntS(n_s)
            end if
            if ( abs(TayRIntN(n_s)) > 0.0_cp ) then
               TayRIntN(n_s)=RstrIntN(n_s)/TayRIntN(n_s)
               TayRIntS(n_s)=RstrIntS(n_s)/TayRIntS(n_s)
            end if
            if ( abs(TayVIntN(n_s)) > 0.0_cp ) then
               TayVIntN(n_s)=StrIntN(n_s)/TayVIntN(n_s)
               TayVIntS(n_s)=StrIntS(n_s)/TayVIntS(n_s)
            end if

            !-- Boundary Values:
            if ( l_mag ) then
               if ( lTC ) then ! Inside TC
                  zMax = sqrt(r_CMB*r_CMB-cyl(n_s)*cyl(n_s))
                  zMin = -zMax
                  call interp_theta(BspAS(:,1),BspB,r_CMB,cyl(n_s),theta_ord)
                  call interp_theta(BspdAS(:,1),BspdB,r_CMB,cyl(n_s),theta_ord)
                  call interp_theta(BpsdAS(:,1),BpsdB,r_CMB,cyl(n_s),theta_ord)
                  call interp_theta(Bs2AS(:,1),Bs2B,r_CMB,cyl(n_s),theta_ord)
                  call interp_theta(BszAS(:,1),BszB,r_CMB,cyl(n_s),theta_ord)
                  call interp_theta(BpzAS(:,1),BpzB,r_CMB,cyl(n_s),theta_ord)
                  call interp_theta(BzpdAS(:,1),BzpdB,r_CMB,cyl(n_s),theta_ord)
                  call interp_theta(BpzdAS(:,1),BpzdB,r_CMB,cyl(n_s),theta_ord)

                  TauBS(n_s)  =-(BpzB(2)+cyl(n_s)/zMin*BspB(2))
                  dTauBS(n_s) =-(BpzdB(2)+BzpdB(2) + cyl(n_s)/zMin*(BspdB(2)+BpsdB(2)))
                  dTTauBS(n_s)=-(BszB(2)+cyl(n_s)/zMin*Bs2B(2))
                  TauBN(n_s)  =  BpzB(1)+cyl(n_s)/zMax*BspB(1)
                  dTauBN(n_s) =  BpzdB(1)+BzpdB(1) + cyl(n_s)/zMax*(BspdB(1)+BpsdB(1))
                  dTTauBN(n_s)=  BszB(1)+cyl(n_s)/zMax*Bs2B(1)

                  zMax = sqrt(r_ICB*r_ICB-cyl(n_s)*cyl(n_s))
                  zMin = -zMax
                  call interp_theta(BspAS(:,n_r_max),BspB,r_ICB,cyl(n_s),theta_ord)
                  call interp_theta(BspdAS(:,n_r_max),BspdB,r_ICB,cyl(n_s),theta_ord)
                  call interp_theta(BpsdAS(:,n_r_max),BpsdB,r_ICB,cyl(n_s),theta_ord)
                  call interp_theta(Bs2AS(:,n_r_max),Bs2B,r_ICB,cyl(n_s),theta_ord)
                  call interp_theta(BszAS(:,n_r_max),BszB,r_ICB,cyl(n_s),theta_ord)
                  call interp_theta(BpzAS(:,n_r_max),BpzB,r_ICB,cyl(n_s),theta_ord)
                  call interp_theta(BzpdAS(:,n_r_max),BzpdB,r_ICB,cyl(n_s),theta_ord)
                  call interp_theta(BpzdAS(:,n_r_max),BpzdB,r_ICB,cyl(n_s),theta_ord)

                  if ( n_s > 1 ) then
                     TauBS(n_s)  =TauBS(n_s)  +(BpzB(2)+cyl(n_s)/zMin*BspB(2))
                     dTauBS(n_s) =dTauBS(n_s) +(BpzdB(2)+BzpdB(2) +              &
                     &                         cyl(n_s)/zMin*(BspdB(2)+BpsdB(2)))
                     dTTauBS(n_s)=dTTauBS(n_s)+(BszB(2)+cyl(n_s)/zMin*Bs2B(2))
                     TauBN(n_s)  =TauBN(n_s)  -(BpzB(1) +cyl(n_s)/zMax*BspB(1))
                     dTauBN(n_s) =dTauBN(n_s) -(BpzdB(1)+BzpdB(1) +              &
                     &                         cyl(n_s)/zMax*(BspdB(1)+BpsdB(1)))
                     dTTauBN(n_s)=dTTauBN(n_s)-(BszB(1)+cyl(n_s)/zMax*Bs2B(1))
                  else ! n_s=1
                     TauBS(n_s)  =0.0_cp
                     dTauBS(n_s) =0.0_cp
                     dTTauBS(n_s)=0.0_cp
                     TauBN(n_s)  =0.0_cp
                     dTauBN(n_s) =0.0_cp
                     dTTauBN(n_s)=0.0_cp
                  end if
               else
                  zMax = sqrt(r_CMB*r_CMB-cyl(n_s)*cyl(n_s))
                  zMin = -zMax

                  call interp_theta(BspAS(:,1),BspB,r_CMB,cyl(n_s),theta_ord)
                  call interp_theta(BspdAS(:,1),BspdB,r_CMB,cyl(n_s),theta_ord)
                  call interp_theta(BpsdAS(:,1),BpsdB,r_CMB,cyl(n_s),theta_ord)
                  call interp_theta(Bs2AS(:,1),Bs2B,r_CMB,cyl(n_s),theta_ord)
                  call interp_theta(BszAS(:,1),BszB,r_CMB,cyl(n_s),theta_ord)
                  call interp_theta(BpzAS(:,1),BpzB,r_CMB,cyl(n_s),theta_ord)
                  call interp_theta(BzpdAS(:,1),BzpdB,r_CMB,cyl(n_s),theta_ord)
                  call interp_theta(BpzdAS(:,1),BpzdB,r_CMB,cyl(n_s),theta_ord)

                  if ( n_s > 1 ) then
                     TauBS(n_s)  =BpzB(1)+cyl(n_s)/zMax*BspB(1) - &
                     &            BpzB(2)-cyl(n_s)/zMin*BspB(2)
                     dTauBS(n_s) =BpzdB(1)+BzpdB(1) + cyl(n_s)/zMax*(BspdB(1)+BpsdB(1)) - &
                     &            BpzdB(2)-BzpdB(2) - cyl(n_s)/zMin*(BspdB(2)+BpsdB(2))
                     dTTauBS(n_s)=BszB(1)+cyl(n_s)/zMax*Bs2B(1) - &
                     &            BszB(2)-cyl(n_s)/zMin*Bs2B(2)
                     TauBN(n_s)  =TauBS(n_s)
                     dTauBN(n_s) =dTauBS(n_s)
                     dTTauBN(n_s)=dTTauBS(n_s)
                  else ! n_s=1
                     TauBS(n_s)  =0.0_cp
                     dTauBS(n_s) =0.0_cp
                     dTTauBS(n_s)=0.0_cp
                     TauBN(n_s)  =0.0_cp
                     dTauBN(n_s) =0.0_cp
                     dTTauBN(n_s)=0.0_cp
                  end if
               end if
            end if

            !-- s-derivatives:
            !-- Create arrays to be differentiated:
            SVpIntN(n_s) =VpIntN(n_s)/cyl(n_s)
            SVpIntS(n_s) =VpIntS(n_s)/cyl(n_s)
            if ( l_mag ) then
               SBspIntN(n_s)=h(n_s)*cyl(n_s)*cyl(n_s)*BspIntN(n_s)
               SBs2IntN(n_s)=h(n_s)*cyl(n_s)**3*Bs2IntN(n_s)
               SBspIntS(n_s)=h(n_s)*cyl(n_s)*cyl(n_s)*BspIntS(n_s)
               SBs2IntS(n_s)=h(n_s)*cyl(n_s)**3*Bs2IntS(n_s)
            end if
         end do

         !-- Calculate s-derivatives
         call get_ds(SVpIntN, dSVpIntN, cyl)
         call get_dds(SVpIntN, d2SVpIntN, cyl)
         call get_ds(SVpIntS, dSVpIntS, cyl)
         call get_dds(SVpIntS, d2SVpIntS, cyl)
         if ( l_mag ) then
            call get_ds(SBspIntN, dSBspIntN, cyl)
            call get_ds(SBspIntS, dSBspIntS, cyl)
            call get_ds(SBs2IntN, dSBs2IntN, cyl)
            call get_ds(SBs2IntS, dSBs2IntS, cyl)

            TauN(:)=Oh(:)*(dSBspIntN(:)/cyl(:)**2+TauBN(:))
            dTTauN(:)=cyl(:)*d2SVpIntN(:)*Bs2IntN(:) +                       &
            &                    Oh(:)*(dSVpIntN(:)*dSBs2IntN(:)/cyl(:)**2 + &
            &                             cyl(:)*dSVpIntN(:)*dTTauBN(:) )
            TauS(:)=Oh(:)*(dSBspIntS(:)/cyl(:)**2+TauBS(:))
            dTTauS(:)=cyl(:)*d2SVpIntS(:)*Bs2IntS(:) +                       &
            &                    Oh(:)*(dSVpIntS(:)*dSBs2IntS(:)/cyl(:)**2 + &
            &                             cyl(:)*dSVpIntS(:)*dTTauBS(:) )
         else
            dSBspIntN(:)=0.0_cp
            dSBspIntS(:)=0.0_cp
            dSBs2IntN(:)=0.0_cp
            dSBs2IntS(:)=0.0_cp
            TauN(:)     =0.0_cp
            TauS(:)     =0.0_cp
            dTTauN(:)   =0.0_cp
            dTTauS(:)   =0.0_cp
         end if

         !--- Output of z-integral:
         write(n_NHS_file)  real(time,kind=outp),                          &! 3
         &      (real(VpIntN(n_s),kind=outp)             ,n_s=1,n_s_max),  &! 4
         &      (real(dVpIntN(n_s),kind=outp)            ,n_s=1,n_s_max),  &! 5
         &      (real(ddVpIntN(n_s),kind=outp)           ,n_s=1,n_s_max),  &! 6
         &      (real(VpRIntN(n_s),kind=outp)            ,n_s=1,n_s_max),  &! 7
         &      (real(RstrIntN(n_s),kind=outp)           ,n_s=1,n_s_max),  &! 8
         &      (real(AstrIntN(n_s),kind=outp)           ,n_s=1,n_s_max),  &! 9
         &      (real(LFfac*LFIntN(n_s),kind=outp)       ,n_s=1,n_s_max),  &! 10
         &      (real(StrIntN(n_s),kind=outp)            ,n_s=1,n_s_max),  &! 11
         &      (real(TayIntN(n_s),kind=outp)            ,n_s=1,n_s_max),  &! 12
         &      (real(LFfac*TauN(n_s),kind=outp)         ,n_s=1,n_s_max),  &! 13
         &      (real(LFfac*TauBN(n_s)*Oh(n_s),kind=outp) ,n_s=1,n_s_max), &! 14
         &      (real(LFfac*BspdIntN(n_s),kind=outp)     ,n_s=1,n_s_max),  &! 15 For first part of dTau
         &      (real(LFfac*BpsdIntN(n_s),kind=outp)     ,n_s=1,n_s_max),  &! 16 For second part of dTau
         &      (real(LFfac*dTauBN(n_s)*Oh(n_s),kind=outp) ,n_s=1,n_s_max),&! 17 Boundary contribution !
         &      (real(LFfac*dTTauN(n_s),kind=outp)       ,n_s=1,n_s_max),  &! 18
         &      (real(LFfac*dTTauBN(n_s)*Oh(n_s),kind=outp),n_s=1,n_s_max),&! 19
         &      (real(LFfac*Bs2IntN(n_s),kind=outp)      ,n_s=1,n_s_max)    ! 20

         write(n_SHS_file)  real(time,kind=outp),                          &! 3
         &      (real(VpIntS(n_s),kind=outp)             ,n_s=1,n_s_max),  &! 4
         &      (real(dVpIntS(n_s),kind=outp)            ,n_s=1,n_s_max),  &! 5
         &      (real(ddVpIntS(n_s),kind=outp)           ,n_s=1,n_s_max),  &! 6
         &      (real(VpRIntS(n_s),kind=outp)            ,n_s=1,n_s_max),  &! 7
         &      (real(RstrIntS(n_s),kind=outp)           ,n_s=1,n_s_max),  &! 8
         &      (real(AstrIntS(n_s),kind=outp)           ,n_s=1,n_s_max),  &! 9
         &      (real(LFfac*LFIntS(n_s),kind=outp)       ,n_s=1,n_s_max),  &! 10
         &      (real(StrIntS(n_s),kind=outp)            ,n_s=1,n_s_max),  &! 11
         &      (real(TayIntS(n_s),kind=outp)            ,n_s=1,n_s_max),  &! 12
         &      (real(LFfac*TauS(n_s),kind=outp)         ,n_s=1,n_s_max),  &! 13
         &      (real(LFfac*TauBS(n_s)*Oh(n_s),kind=outp),n_s=1,n_s_max),  &! 14
         &      (real(LFfac*BspdIntS(n_s),kind=outp)     ,n_s=1,n_s_max),  &! 15 For first part of dTau
         &      (real(LFfac*BpsdIntS(n_s),kind=outp)     ,n_s=1,n_s_max),  &! 16 For second part of dTau
         &      (real(LFfac*dTauBS(n_s)*Oh(n_s),kind=outp),n_s=1,n_s_max), &! 17 Boundary contribution !
         &      (real(LFfac*dTTauS(n_s),kind=outp)       ,n_s=1,n_s_max),  &! 18
         &      (real(LFfac*dTTauBS(n_s)*Oh(n_s),kind=outp),n_s=1,n_s_max),&! 19
         &      (real(LFfac*Bs2IntS(n_s),kind=outp)      ,n_s=1,n_s_max)    ! 20

         !-- Integrate Geostrophic azimuthal flow energy and Taylor measure:
         VgRMS = simps(VpIntN*VpIntN*cyl*h, cyl)
         TayRMS = simps(abs(TayIntN)*cyl*h, cyl)
         TayRRMS = simps(abs(TayRIntN)*cyl*h, cyl)
         TayVRMS = simps(abs(TayVIntN)*cyl*h, cyl)
         if ( .not. l_full_sphere ) then
            VgRMS = VgRMS+simps(VpIntS(n_s_otc+1:n_s_max)*VpIntS(n_s_otc+1:n_s_max)* &
            &                   cyl(n_s_otc+1:n_s_max)*h(n_s_otc+1:n_s_max),         &
            &                   cyl(n_s_otc+1:n_s_max))
            TayRMS = TayRMS+simps(abs(TayIntS(n_s_otc+1:n_s_max))*cyl(n_s_otc+1:n_s_max)* &
            &                     h(n_s_otc+1:n_s_max), cyl(n_s_otc+1:n_s_max))
            TayRRMS = TayRRMS+simps(abs(TayRIntS(n_s_otc+1:n_s_max))*cyl(n_s_otc+1:n_s_max)*&
            &                       h(n_s_otc+1:n_s_max), cyl(n_s_otc+1:n_s_max))
            TayVRMS = TayVRMS+simps(abs(TayVIntS(n_s_otc+1:n_s_max))*cyl(n_s_otc+1:n_s_max)*&
            &                       h(n_s_otc+1:n_s_max), cyl(n_s_otc+1:n_s_max))
         end if

         !-- Normalisation factors
         VgRMS  =sqrt(two*pi*VgRMS/volcyl_oc)
         TayRMS =two*pi*TayRMS/volcyl_oc
         TayRRMS= two*pi*TayRRMS/volcyl_oc
         TayVRMS= two*pi*TayVRMS/volcyl_oc

         !--- Relative importance of azimuthal and geostrophic flow
         VRMS =sqrt(two*eKin/vol_oc)
         VpRMS=sqrt(two*eKinTAS/vol_oc)

         if ( VRMS /= 0.0_cp) then
            VpRMS=VpRMS/VRMS
            VgRMS=VgRMS/VRMS
         end if

         !-- Time series Tay.TAG
         if ( l_save_out ) then
            open(newunit=n_TO_file, file=TOFile, status='unknown', &
            &    position='append')
         end if
         write(n_TO_file, '(1P,ES20.12,6ES16.8)') time, VpRMS**2, VgRMS**2, TayRMS, &
         &                                        TayRRMS, TayVRMS, eKin
         if ( l_save_out ) close(n_TO_file)

         !-- To movie outputs
         if ( lTOmov ) then
            nTOmovSets=nTOmovSets+1

            !-- Write time and frame number
            dumm(1)=nTOmovSets        ! time frame number for movie
            dumm(2)=time              ! time
            dumm(3)=0.0_cp
            dumm(4)=0.0_cp
            dumm(5)=0.0_cp
            dumm(6)=0.0_cp
            dumm(7)=0.0_cp
            dumm(8)=0.0_cp
            write(n_TOmov_file) (real(dumm(n),kind=outp),n=1,8)

            !-- Write fields
            write(n_TOmov_file) real(VAS,kind=outp)
            write(n_TOmov_file) real(dzRstrAS,kind=outp)
            write(n_TOmov_file) real(dzAstrAS,kind=outp)
            write(n_TOmov_file) real(dzStrAS,kind=outp)
            if ( l_mag ) then
               write(n_TOmov_file) real(LFfac*dzLFAS,kind=outp)
            else
               write(n_TOmov_file) real(0.0_cp*VAS,kind=outp)
            end if
            write(n_TOmov_file) real(dzCorAS,kind=outp)
            write(n_TOmov_file) real(dzdVpAS,kind=outp)
            if ( l_phase_field ) write(n_TOmov_file) real(dzPenAS,kind=outp)
         end if

         if ( l_save_out ) then
            open(newunit=n_log_file, file=log_file, status='unknown', &
            &    position='append')
         end if
         if ( lTOmov ) then
            call logWrite(' ')
            write(message,'(1p,A,I8,A,ES16.6,I8)')              &
            &      "! WRITING TO MOVIE FRAME NO ",nTOmovSets,   &
            &      "       at time/step",time*tScale,n_time_step
            call logWrite(message)
         end if
         if ( l_save_out ) close(n_log_file)

      end if

   end subroutine outTO
!------------------------------------------------------------------------------
   subroutine cylmean(dat, datN, datS)
      !
      ! This routine computes the z-average inside and outside TC
      !

      !-- Input variables
      real(cp), intent(in) :: dat(:,:) ! input data

      !-- Output variables
      real(cp), intent(out) :: datN(n_s_max)   ! z-average outside T.C. + N.H.
      real(cp), intent(out) :: datS(n_s_max)   ! z-average outside T.C. + S.H.

      !-- Local variables
      real(cp) :: dat_OTC(n_s_max), dat_ITC_N(n_s_max), dat_ITC_S(n_s_max)

      !-- z-averaging
      call cylmean_otc(dat,dat_OTC,n_s_max,n_s_otc,r,cyl,theta_ord,zDens)
      call cylmean_itc(dat,dat_ITC_N,dat_ITC_S,n_s_max,n_s_otc,r,cyl,theta_ord,zDens)

      !-- rearange
      datN(1:n_s_otc)=dat_OTC(1:n_s_otc)
      datS(1:n_s_otc)=dat_OTC(1:n_s_otc)
      datN(n_s_otc+1:n_s_max)=dat_ITC_N(n_s_otc+1:n_s_max)
      datS(n_s_otc+1:n_s_max)=dat_ITC_S(n_s_otc+1:n_s_max)

   end subroutine cylmean
!------------------------------------------------------------------------------
   subroutine interp_theta(a,ac,rr,cyl,theta)
      !
      ! This routine computes the interpolation of a value at the surface of a
      ! spherical shell onto the cylindrical grid. This is only a theta interpolation
      ! using the cylindrical theta's.
      !

      !-- Input variables
      real(cp), intent(in) :: rr       ! Radius at which we compute the extrapolation (can be ri or ro)
      real(cp), intent(in) :: cyl      ! Cylindrical radius
      real(cp), intent(in) :: theta(:) ! Colatitude
      real(cp), intent(in) :: a(:)     ! Field at the outer radius

      !-- Output variable
      real(cp), intent(out) :: ac(2)   ! Surface values for NH and SH as a function of s

      !-- Local variables
      integer :: n_hs, n_th
      integer :: n_th0, n_th1, n_th2, n_th3
      real(cp) :: eps, thet, tt0, tt1, tt2, tt3, t10, t20, t30, t21, t31, t32
      real(cp) :: a01, a12, a23, a012, a123

      eps=10.0_cp*epsilon(1.0_cp)
      ac(:)=0.0_cp

      !-- Loop over hemispheres
      do n_hs=1,2

         !-- Get the polar angle in cylindrical coordinates
         if (n_hs == 1) then ! Northern Hemisphere
            thet=asin(cyl/rr) ! theta_c(NH) = arcsin(s/rr)
         else ! Southern Hemisphere
            thet=pi-asin(cyl/rr) ! theta_c(SH)=pi-arcsin(s/rr)
         end if

         !-- Find indices of angular grid levels that bracket thet
         n_th1=n_theta_max
         tbracket: do n_th=n_theta_max,1,-1
            if ( theta(n_th) <= thet) then
               n_th1=n_th
               exit tbracket
            end if
         end do tbracket
         if ( n_th1 == n_theta_max ) n_th1=n_theta_max-2
         if ( n_th1 == n_theta_max-1 ) n_th1=n_theta_max-2
         if ( n_th1 == 1 ) n_th1=2
         n_th2=n_th1+1
         n_th3=n_th1+2
         n_th0=n_th1-1

         !--  Calculate differences in theta for 4th-order interpolation
         tt0=thet-theta(n_th0)
         tt1=thet-theta(n_th1)
         tt2=thet-theta(n_th2)
         tt3=thet-theta(n_th3)
         t10=1.0_cp/(theta(n_th1)-theta(n_th0))
         t20=1.0_cp/(theta(n_th2)-theta(n_th0))
         t30=1.0_cp/(theta(n_th3)-theta(n_th0))
         t21=1.0_cp/(theta(n_th2)-theta(n_th1))
         t31=1.0_cp/(theta(n_th3)-theta(n_th1))
         t32=1.0_cp/(theta(n_th3)-theta(n_th2))

         !-- Interpolation in theta-direction
         a01=(tt0*a(n_th1)-tt1*a(n_th0))*t10
         a12=(tt1*a(n_th2)-tt2*a(n_th1))*t21
         a23=(tt2*a(n_th3)-tt3*a(n_th2))*t32

         a012=(tt0*a12-tt2*a01)*t20
         a123=(tt1*a23-tt3*a12)*t31
         ac(n_hs)=ac(n_hs)+(tt0*a123-tt3*a012)*t30
      end do

      !-- Set the boundary points
      if ( abs(rr*sin(theta(1))-cyl) <= eps ) then
         ac(1)=a(n_theta_max)
         ac(2)=a(1)
      end if

   end subroutine interp_theta
!------------------------------------------------------------------------------------
   subroutine get_ds(arr, darr, cyl)
      !
      ! This subroutine is used to compute the 4th order accurate first s-derivative
      ! on the regularly spaced grid
      !

      !-- Input variables
      real(cp), intent(in) :: arr(:) ! Array to be differentiated
      real(cp), intent(in) :: cyl(:) ! Cylindrical grid

      !-- Output variables
      real(cp), intent(out) :: darr(:) ! s-derivative of the input array

      !-- Local variables
      integer :: nsMax,n_s
      real(cp) :: ds, f1

      nsMax = size(cyl)
      ds = (cyl(nsMax)-cyl(1))/real(nsMax-1,cp)

      !-- Bulk points
      f1 = one/(12.0_cp*ds)
      do n_s=3,nsMax-2
         darr(n_s)=f1*(arr(n_s-2)-8.0_cp*arr(n_s-1)+8.0_cp*arr(n_s+1)-arr(n_s+2))
      end do

      !-- Boundary points
      darr(2)=f1*(-3.0_cp*arr(1)-10.0_cp*arr(2)+18.0_cp*arr(3)-6.0_cp*arr(4)+arr(5))
      darr(1)=f1*(-25.0_cp*arr(1)+48.0_cp*arr(2)-36.0_cp*arr(3)+16.0_cp*arr(4)- &
      &           3.0_cp*arr(5))
      darr(nSMax-1)=f1*(3.0_cp*arr(nSMax)+10.0_cp*arr(nSMax-1)-18.0_cp*arr(nSMax-2) &
      &                 +6.0_cp*arr(nSMax-3)-arr(nSMax-4))
      darr(nSmax)=f1*(25.0_cp*arr(nSmax)-48.0_cp*arr(nSmax-1)+36.0_cp*arr(nSmax-2)-&
      &               16.0_cp*arr(nSmax-3)+3.0_cp*arr(nSmax-4))

   end subroutine get_ds
!------------------------------------------------------------------------------------
   subroutine get_dds(arr, ddarr, cyl)
      !
      ! This subroutine is used to compute the 4th order accurate 2nd s-derivative 
      ! on the regularly spaced grid
      !
      ! https://bellaard.com/tools/Finite%20difference%20coefficient%20calculator/
      !

      !-- Input variables
      real(cp), intent(in) :: arr(:) ! Array to be differentiated
      real(cp), intent(in) :: cyl(:) ! Cylindrical grid

      !-- Output variables
      real(cp), intent(out) :: ddarr(:) ! s-derivative of the input array

      !-- Local variables
      integer :: nsMax,n_s
      real(cp) :: ds, f1

      nsMax = size(cyl)
      ds = (cyl(nsMax)-cyl(1))/real(nsMax-1,cp)

      !-- Bulk points
      f1 = one/(12.0_cp*ds*ds)
      do n_s=3,nsMax-2
         ddarr(n_s)=f1*(-arr(n_s-2)+16.0_cp*arr(n_s-1)-30.0_cp*arr(n_s)+ &
         &              16.0_cp*arr(n_s+1)-arr(n_s+2))
      end do

      !-- Boundary points
      ddarr(1)=f1*(45.0_cp*arr(1)-154.0_cp*arr(2)+214.0_cp*arr(3)-156.0_cp*arr(4)+&
      &            61.0_cp*arr(5)-10.0_cp*arr(6))
      ddarr(2)=f1*(10.0_cp*arr(1)-15.0_cp*arr(2)-4.0_cp*arr(3)+14.0_cp*arr(4)-&
      &            6.0_cp*arr(5)+arr(6))
      ddarr(nSmax)=f1*(45.0_cp*arr(nSmax)-154.0_cp*arr(nSmax-1)+214.0_cp*arr(nSmax-2)-&
      &                156.0_cp*arr(nSmax-3)+61.0_cp*arr(nSmax-4)-10.0_cp*arr(nSmax-5))
      ddarr(nSmax-1)=f1*(10.0_cp*arr(nSmax)-15.0_cp*arr(nSmax-1)-4.0_cp*arr(nSmax-2)+ &
      &                  14.0_cp*arr(nSmax-3)-6.0_cp*arr(nSmax-4)+arr(nSmax-5))

   end subroutine get_dds
!------------------------------------------------------------------------------------
   subroutine gather_from_Rloc_to_rank0(arr_Rloc, arr)
      !
      ! This subroutine gathers the r-distributed array
      !
      
      !-- Input variable
      real(cp), intent(in) :: arr_Rloc(n_theta_max,nRstart:nRstop)

      !-- Output variable
      real(cp), intent(out) :: arr(n_theta_max,n_r_max)

#ifdef WITH_MPI
      !-- Local variables:
      integer :: sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1)
      integer :: p

      sendcount  = nR_per_rank*n_theta_max
      do p=0,n_procs-1
         recvcounts(p)=radial_balance(p)%n_per_rank*n_theta_max
      end do
      displs(0)=0
      do p=1,n_procs-1
         displs(p)=displs(p-1)+recvcounts(p-1)
      end do

      call MPI_GatherV(arr_Rloc, sendcount, MPI_DEF_REAL, arr, recvcounts, &
           &           displs, MPI_DEF_REAL, 0, MPI_COMM_WORLD,ierr)

#else
      arr(:,:)=arr_Rloc(:,:)
#endif

   end subroutine gather_from_Rloc_to_rank0
!------------------------------------------------------------------------------
end module outTO_mod
