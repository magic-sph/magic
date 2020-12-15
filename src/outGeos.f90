module geos
   !
   ! This module is used to compute z-integrated diagnostics such as the degree
   ! of geostrophy or the separation of energies between inside and outside the
   ! tangent cylinder. This makes use of a local Simpson's method.
   !

   use precision_mod
   use parallel_mod
   use constants, only: half, two, pi, one, four, third
   use mem_alloc, only: bytes_allocated
   use radial_data, only: radial_balance, nRstart, nRstop
   use radial_functions, only: or1, or2, r_ICB, r_CMB, r, orho1, orho2, beta
   use output_data, only: sDens, tag
   use horizontal_data, only: n_theta_cal2ord, O_sin_theta_E2, theta_ord, &
       &                      O_sin_theta, cosTheta, sinTheta
   use truncation, only: n_phi_max, n_theta_max, n_r_max
   use integration, only: simps
   use logic, only: l_save_out

   implicit none

   private

   real(cp), allocatable :: cyl(:) ! Cylindrical grid
   real(cp), allocatable :: h(:)   ! h(s)
   real(cp), allocatable :: us_Rloc(:,:,:), up_Rloc(:,:,:), uz_Rloc(:,:,:)
   real(cp), allocatable :: us_Ploc(:,:,:), up_Ploc(:,:,:), uz_Ploc(:,:,:)
   real(cp), allocatable :: wz_Rloc(:,:,:), wz_Ploc(:,:,:)
   real(cp) :: vol_otc ! volume outside tangent cylinder

   integer :: nPstart ! Starting nPhi index when MPI distributed
   integer :: nPstop  ! Stoping nPhi index when MPI distributed
   type(load), allocatable :: phi_balance(:) ! phi-distributed balance
   integer :: n_geos_file ! file unit for geos.TAG
   integer :: n_s_max ! Number of cylindrical points
   integer :: n_s_otc ! Index for last point outside TC
   character(len=72) :: geos_file ! file name

   public :: initialize_geos, finalize_geos, calcGeos, outGeos

contains

   subroutine initialize_geos()
      !
      ! Memory allocation and definition of the cylindrical grid
      !

      integer :: n_s
      real(cp) :: smin, smax, ds

      !-- Volume outside tangent cylinder
      vol_otc = four*third*pi*(r_CMB**2-r_ICB**2)**(1.5_cp)

      !-- R-distributed arrays
      allocate( us_Rloc(n_theta_max,n_phi_max,nRstart:nRstop) )
      allocate( up_Rloc(n_theta_max,n_phi_max,nRstart:nRstop) )
      allocate( uz_Rloc(n_theta_max,n_phi_max,nRstart:nRstop) )
      allocate( wz_Rloc(n_theta_max,n_phi_max,nRstart:nRstop) )
      bytes_allocated=bytes_allocated+4*n_phi_max*n_theta_max*(nRstop-nRstart+1)*&
      &               SIZEOF_DEF_REAL

      !-- Distribute over the ranks
      allocate(phi_balance(0:n_procs-1))
      call getBlocks(phi_balance, n_phi_max, n_procs)
      nPstart = phi_balance(rank)%nStart
      nPstop = phi_balance(rank)%nStop

      !-- Phi-distributed arrays
      allocate( us_Ploc(n_theta_max,n_r_max,nPstart:nPstop) )
      allocate( up_Ploc(n_theta_max,n_r_max,nPstart:nPstop) )
      allocate( uz_Ploc(n_theta_max,n_r_max,nPstart:nPstop) )
      allocate( wz_Ploc(n_theta_max,n_r_max,nPstart:nPstop) )
      bytes_allocated=bytes_allocated+4*n_r_max*n_theta_max*(nPstop-nPstart+1)*&
      &               SIZEOF_DEF_REAL

      !-- Cylindrical radius
      n_s_max = n_r_max+int(r_ICB*n_r_max)
      n_s_max = int(sDens*n_s_max)
      allocate( cyl(n_s_max) )
      bytes_allocated=bytes_allocated+n_s_max*SIZEOF_DEF_REAL

      ! Minimum and maximum cylindrical radii
      smin = r_CMB*sin(theta_ord(1))
      smax = r_CMB

      !-- Grid spacing
      ds = (smax-smin)/(n_s_max-1)

      !-- Cylindrical grid
      do n_s=1,n_s_max
         cyl(n_s)=r_cmb-(n_s-1)*ds
      end do

      !-- Height
      allocate( h(n_s_max) )
      do n_s=1,n_s_max
         if ( cyl(n_s) >= r_ICB ) then
            h(n_s)=two*sqrt(r_CMB**2-cyl(n_s)**2)
            n_s_otc=n_s ! Last index outside TC
         else
            h(n_s)=sqrt(r_CMB**2-cyl(n_s)**2)-sqrt(r_ICB**2-cyl(n_s)**2)
         end if
      end do

      !-- Open geos.TAG file
      if ( rank == 0 ) then
         geos_file='geos.'//tag
         if ( (.not. l_save_out) ) then
            open(newunit=n_geos_file, file=geos_file, status='new')
         end if
      end if

   end subroutine initialize_geos
!------------------------------------------------------------------------------------
   subroutine finalize_geos()
      !
      ! Memory deallocation
      !

      if ( rank == 0 .and. (.not. l_save_out) ) close(n_geos_file)
      deallocate( us_Ploc, up_Ploc, uz_Ploc )
      deallocate( us_Rloc, up_Rloc, uz_Rloc )
      deallocate( cyl, h )

   end subroutine finalize_geos
!------------------------------------------------------------------------------------
   subroutine calcGeos(vr,vt,vp,cvr,dvrdp,dvpdr,nR)
      !
      ! This routine computes the term needed for geos.TAG outputs in physical
      ! space.
      !

      !-- Input variables
      integer,  intent(in) :: nR ! Radial grid point
      real(cp), intent(in) :: vr(:,:), vt(:,:), vp(:,:)
      real(cp), intent(in) :: cvr(:,:), dvrdp(:,:), dvpdr(:,:)

      !-- Local variables
      integer :: nPhi, nTheta, nTheta1

      !$omp parallel do default(shared) private(nTheta,nPhi,nTheta1)
      do nPhi=1,n_phi_max
         do nTheta=1,n_theta_max
            nTheta1=n_theta_cal2ord(nTheta)
            !- us=ur*sin(theta)+ut*cos(theta)
            us_Rloc(nTheta1,nPhi,nR)=or2(nR)*orho1(nR)*sinTheta(nTheta)*vr(nTheta,nPhi) &
            &                        +or1(nR)*orho1(nR)*cosTheta(nTheta)*               &
            &                        O_sin_theta(nTheta)*vt(nTheta,nPhi)

            !- uz=ur*cos(theta)-ut*sin(theta)
            uz_Rloc(nTheta1,nPhi,nR)=or2(nR)*orho1(nR)*cosTheta(nTheta)*vr(nTheta,nPhi) &
            &                        -or1(nR)*orho1(nR)*vt(nTheta,nPhi)

            !- uphi
            up_Rloc(nTheta1,nPhi,nR)=or1(nR)*O_sin_theta(nTheta)*orho1(nR)* &
            &                        vp(nTheta,nPhi)

            !-- z-vorticity
            wz_Rloc(nTheta1,nPhi,nR)=or1(nR)*orho1(nR)*(                 &
            &                cosTheta(nTheta)*or1(nR)*cvr(nTheta,nPhi) - &
            &            or2(nR)*dvrdp(nTheta,nPhi)+dvpdr(nTheta,nPhi) - &
            &                       beta(nR)*vp(nTheta,nPhi)   )
         end do
      end do
      !$omp end parallel do

   end subroutine calcGeos
!------------------------------------------------------------------------------------
   subroutine outGeos(time,Geos,GeosA,GeosZ,GeosM,GeosNAP,Ekin)
      !
      ! This routine handles the output of geos.TAG
      !

      !-- Input variable
      real(cp), intent(in) :: time

      !-- Output variables
      real(cp), intent(out) :: Geos, GeosA, GeosZ, GeosM, GeosNAP, Ekin

      !-- Local variables
      real(cp) :: phiNorm
      real(cp) :: tmp(n_theta_max,n_r_max)
      real(cp) :: us_axi(n_theta_max,n_r_max),us_axi_dist(n_theta_max,n_r_max)
      real(cp) :: up_axi(n_theta_max,n_r_max),up_axi_dist(n_theta_max,n_r_max)
      real(cp) :: uz_axi(n_theta_max,n_r_max),uz_axi_dist(n_theta_max,n_r_max)
      real(cp) :: CVZ_dist(n_s_max), CVZ_full(n_s_max), CVZ
      real(cp) :: CVor_dist(n_s_max), CVor_full(n_s_max), CVor
      real(cp) :: CHel_dist(n_s_max), CHel_full(n_s_max), CHel
      real(cp) :: Ekin_dist(n_s_max), Ekin_full(n_s_max)
      real(cp) :: Ek_nAP_dist(n_s_max), Ek_nAP_full(n_s_max), Ek_nAP
      real(cp) :: Ek_NTC_dist(n_s_max), Ek_NTC_full(n_s_max), EkNTC
      real(cp) :: Ek_STC_dist(n_s_max), Ek_STC_full(n_s_max), EkSTC
      real(cp) :: Egeos_dist(n_s_max), Egeos_full(n_s_max), Egeos
      real(cp) :: Eg_nAP_dist(n_s_max), Eg_nAP_full(n_s_max), Eg_nAP
      real(cp) :: uzSN_OTC(n_s_max), uzNN_OTC(n_s_max), wzNN_OTC(n_s_max)
      real(cp) :: Ek_ITC_N(n_s_max), Ek_ITC_S(n_s_max), Ek_OTC(n_s_max)
      real(cp) :: us_ITC_N(n_s_max), us_ITC_S(n_s_max), us_OTC(n_s_max)
      real(cp) :: up_ITC_N(n_s_max), up_ITC_S(n_s_max), up_OTC(n_s_max)
      real(cp) :: uz_ITC_N(n_s_max), uz_ITC_S(n_s_max), uz_OTC(n_s_max)
      real(cp) :: EkS, EkN, EA, EM, EZ, EgeosA, EgeosM, EgeosZ
      integer :: n_p, n_s, n_t

      phiNorm = one/n_phi_max

      !-- MPI transpositions for  us, uz and up
      call transp_R2Phi(us_Rloc, us_Ploc)
      call transp_R2Phi(uz_Rloc, uz_Ploc)
      call transp_R2Phi(up_Rloc, up_Ploc)
      call transp_R2Phi(wz_Rloc, wz_Ploc)

      !-- Get axisymmetric flows for partial geos business
      us_axi_dist(:,:)=0.0_cp
      up_axi_dist(:,:)=0.0_cp
      uz_axi_dist(:,:)=0.0_cp
      do n_p=nPstart,nPstop
         us_axi_dist(:,:)=us_axi_dist(:,:)+phiNorm*us_Ploc(:,:,n_p)
         up_axi_dist(:,:)=up_axi_dist(:,:)+phiNorm*up_Ploc(:,:,n_p)
         uz_axi_dist(:,:)=uz_axi_dist(:,:)+phiNorm*uz_Ploc(:,:,n_p)
      end do

#ifdef WITH_MPI
      call MPI_AllReduce(us_axi_dist,us_axi,n_r_max*n_theta_max,MPI_DEF_REAL,MPI_SUM, &
           &             MPI_COMM_WORLD,ierr)
      call MPI_AllReduce(up_axi_dist,up_axi,n_r_max*n_theta_max,MPI_DEF_REAL,MPI_SUM, &
           &             MPI_COMM_WORLD,ierr)
      call MPI_AllReduce(uz_axi_dist,uz_axi,n_r_max*n_theta_max,MPI_DEF_REAL,MPI_SUM, &
           &             MPI_COMM_WORLD,ierr)
#else
      us_axi(:,:)=us_axi_dist(:,:)
      up_axi(:,:)=up_axi_dist(:,:)
      uz_axi(:,:)=uz_axi_dist(:,:)
#endif

      if ( rank == 0 ) then
         !-- Axisymmetric kinetic energy
         tmp(:,:) = half*(us_axi(:,:)**2+up_axi(:,:)**2+uz_axi(:,:)**2)
         call cylmean(tmp,Ek_OTC,Ek_ITC_N,Ek_ITC_S)
         Ekin_full(:)=Ek_ITC_N(:)+Ek_ITC_S(:)+Ek_OTC(:)
         EA=two*pi*simps(Ekin_full*cyl*h,cyl)

         !-- Meridional kinetic energy
         tmp(:,:) = half*(us_axi(:,:)**2+uz_axi(:,:)**2)
         call cylmean(tmp,Ek_OTC,Ek_ITC_N,Ek_ITC_S)
         Ekin_full(:)=Ek_ITC_N(:)+Ek_ITC_S(:)+Ek_OTC(:)
         EM=two*pi*simps(Ekin_full*cyl*h,cyl)

         !-- Zonal kinetic energy
         tmp(:,:) = half*up_axi(:,:)**2
         call cylmean(tmp,Ek_OTC,Ek_ITC_N,Ek_ITC_S)
         Ekin_full(:)=Ek_ITC_N(:)+Ek_ITC_S(:)+Ek_OTC(:)
         EZ=two*pi*simps(Ekin_full*cyl*h,cyl)

         !-- Axisymmetric geostrophic energy
         call cylmean(us_axi,us_OTC,us_ITC_N,us_ITC_S)
         call cylmean(up_axi,up_OTC,up_ITC_N,up_ITC_S)
         call cylmean(uz_axi,uz_OTC,uz_ITC_N,uz_ITC_S)
         Egeos_full(:)= half*(us_ITC_N(:)**2+us_ITC_S(:)**2+us_OTC(:)**2+ &
         &                    up_ITC_N(:)**2+up_ITC_S(:)**2+up_OTC(:)**2+ &
         &                    uz_ITC_N(:)**2+uz_ITC_S(:)**2+uz_OTC(:)**2)
         EgeosA=two*pi*simps(Egeos_full*cyl*h,cyl)

         !-- Meridional geostrophic energy
         Egeos_full(:)=half*(us_ITC_N(:)**2+us_ITC_S(:)**2+us_OTC(:)**2+ &
         &                   uz_ITC_N(:)**2+uz_ITC_S(:)**2+uz_OTC(:)**2)
         EgeosM=two*pi*simps(Egeos_full*cyl*h,cyl)

         !-- Zonal geostrophic energy
         Egeos_full(:)=half*(up_ITC_N(:)**2+up_ITC_S(:)**2+up_OTC(:)**2)
         EgeosZ=two*pi*simps(Egeos_full*cyl*h,cyl)
      end if

      !-- Loop over phi points
      Egeos_dist(:) =0.0_cp
      Eg_nAP_dist(:)=0.0_cp
      Ekin_dist(:)  =0.0_cp
      Ek_nAP_dist(:)=0.0_cp
      Ek_NTC_dist(:)=0.0_cp
      Ek_STC_dist(:)=0.0_cp
      CVZ_dist(:)   =0.0_cp
      CVor_dist(:)  =0.0_cp
      CHel_dist(:)  =0.0_cp
      do n_p=nPstart,nPstop

         !-- Total kinetic energy
         tmp(:,:) = half*(us_Ploc(:,:,n_p)**2+up_Ploc(:,:,n_p)**2+uz_Ploc(:,:,n_p)**2)
         call cylmean(tmp,Ek_OTC,Ek_ITC_N,Ek_ITC_S)
         Ekin_dist(:)=Ekin_dist(:)+Ek_ITC_N(:)+Ek_ITC_S(:)+Ek_OTC(:)
         Ek_NTC_dist(:)=Ek_NTC_dist(:)+Ek_ITC_N(:)
         Ek_STC_dist(:)=Ek_STC_dist(:)+Ek_ITC_S(:)

         !-- Non-axisymmetric kinetic energy and perpendicular to rotation axis
         tmp(:,:) = half*((us_Ploc(:,:,n_p)-us_axi(:,:))**2+ &
         &                (up_Ploc(:,:,n_p)-up_axi(:,:))**2)
         call cylmean(tmp,Ek_OTC,Ek_ITC_N,Ek_ITC_S)
         Ek_nAP_dist(:)=Ek_nAP_dist(:)+Ek_ITC_N(:)+Ek_ITC_S(:)+Ek_OTC(:)

         !-- Geostrophic energy
         call cylmean(us_Ploc(:,:,n_p),us_OTC,us_ITC_N,us_ITC_S)
         call cylmean(up_Ploc(:,:,n_p),up_OTC,up_ITC_N,up_ITC_S)
         call cylmean(uz_Ploc(:,:,n_p),uz_OTC,uz_ITC_N,uz_ITC_S)

         Egeos_dist(:)=Egeos_dist(:)+                                        &
         &                 half*(us_ITC_N(:)**2+us_ITC_S(:)**2+us_OTC(:)**2+ &
         &                       up_ITC_N(:)**2+up_ITC_S(:)**2+up_OTC(:)**2+ &
         &                       uz_ITC_N(:)**2+uz_ITC_S(:)**2+uz_OTC(:)**2)

         !-- Geostrophic non-axisymmetric energy
         call cylmean(us_Ploc(:,:,n_p)-us_axi,us_OTC,us_ITC_N,us_ITC_S)
         call cylmean(up_Ploc(:,:,n_p)-up_axi,up_OTC,up_ITC_N,up_ITC_S)
         Eg_nAP_dist(:)=Eg_nAP_dist(:)+                                      &
         &                 half*(us_ITC_N(:)**2+us_ITC_S(:)**2+us_OTC(:)**2+ &
         &                       up_ITC_N(:)**2+up_ITC_S(:)**2+up_OTC(:)**2)

         !-- Correlations
         ! Calculate the Pearson correlation coeff between
         ! z-integrated stuff in northern and southern HS.
         ! Considered are Vz,z-vorticity Vor, and axial helicity Vz*Vor.
         ! The correlation is averaged over the shell.
         !-- Reuse tmp as a work array

         !-- Uz North-South correlation
         do n_t=1,n_theta_max
            tmp(n_t,:)=uz_Ploc(n_t,:,n_p)*uz_Ploc(n_theta_max+1-n_t,:,n_p)
         end do
         call cylmean_otc(tmp,uzSN_OTC,n_s_max,n_s_otc,r,cyl,theta_ord)
         tmp(:,:)=uz_Ploc(:,:,n_p)*uz_Ploc(:,:,n_p)
         call cylmean_otc(tmp,uzNN_OTC,n_s_max,n_s_otc,r,cyl,theta_ord)

         do n_s=1,n_s_otc
            if ( uzNN_OTC(n_s) > 0.0_cp ) then
               CVZ_dist(n_s)=CVZ_dist(n_s)+uzSN_OTC(n_s)/uzNN_OTC(n_s)
            else
               CVZ_dist(n_s)=CVZ_dist(n_s)+one
            end if
         end do

         !-- z-vorticity
         do n_t=1,n_theta_max
            tmp(n_t,:)=wz_Ploc(n_t,:,n_p)*wz_Ploc(n_theta_max+1-n_t,:,n_p)
         end do
         call cylmean_otc(tmp,uzSN_OTC,n_s_max,n_s_otc,r,cyl,theta_ord)
         tmp(:,:)=wz_Ploc(:,:,n_p)*wz_Ploc(:,:,n_p)
         call cylmean_otc(tmp,wzNN_OTC,n_s_max,n_s_otc,r,cyl,theta_ord)

         do n_s=1,n_s_otc
            if ( wzNN_OTC(n_s) > 0.0_cp ) then
               CVor_dist(n_s)=CVor_dist(n_s)+uzSN_OTC(n_s)/wzNN_OTC(n_s)
            else
               CVor_dist(n_s)=CVor_dist(n_s)+one
            end if
         end do

         !-- Helicity z: u_z*w_z
         do n_t=1,n_theta_max/2
            tmp(n_t,:)=uz_Ploc(n_t,:,n_p)*wz_Ploc(n_t,:,n_p)- &
            &          uz_Ploc(n_theta_max+1-n_t,:,n_p)*wz_Ploc(n_theta_max+1-n_t,:,n_p)
         end do
         do n_t=n_theta_max/2+1,n_theta_max
            tmp(n_t,:)=-uz_Ploc(n_t,:,n_p)*wz_Ploc(n_t,:,n_p)+ &
            &          uz_Ploc(n_theta_max+1-n_t,:,n_p)*wz_Ploc(n_theta_max+1-n_t,:,n_p)
         end do
         call cylmean_otc(tmp,uzSN_OTC,n_s_max,n_s_otc,r,cyl,theta_ord)

         do n_s=1,n_s_max
            if ( uzNN_OTC(n_s) > 0.0_cp .and. wzNN_OTC(n_s) > 0.0_cp ) then
               CHel_dist(n_s)=CHel_dist(n_s)+uzSN_OTC(n_s)/sqrt(uzNN_OTC(n_s))/ &
               &              sqrt(wzNN_OTC(n_s))
            else
               CHel_dist(n_s)=CHel_dist(n_s)+one
            end if
         end do

      end do

#ifdef WITH_MPI
      call MPI_Reduce(Egeos_dist,Egeos_full,n_s_max,MPI_DEF_REAL,MPI_SUM,0, &
           &          MPI_COMM_WORLD,ierr)
      call MPI_Reduce(Ekin_dist,Ekin_full,n_s_max,MPI_DEF_REAL,MPI_SUM,0, &
           &          MPI_COMM_WORLD,ierr)
      call MPI_Reduce(Ek_nAP_dist,Ek_nAP_full,n_s_max,MPI_DEF_REAL,MPI_SUM,0, &
           &          MPI_COMM_WORLD,ierr)
      call MPI_Reduce(Eg_nAP_dist,Eg_nAP_full,n_s_max,MPI_DEF_REAL,MPI_SUM,0, &
           &          MPI_COMM_WORLD,ierr)
      call MPI_Reduce(Ek_NTC_dist,Ek_NTC_full,n_s_max,MPI_DEF_REAL,MPI_SUM,0, &
           &          MPI_COMM_WORLD,ierr)
      call MPI_Reduce(Ek_STC_dist,Ek_STC_full,n_s_max,MPI_DEF_REAL,MPI_SUM,0, &
           &          MPI_COMM_WORLD,ierr)
      call MPI_Reduce(CVZ_dist,CVZ_full,n_s_max,MPI_DEF_REAL,MPI_SUM,0, &
           &          MPI_COMM_WORLD,ierr)
      call MPI_Reduce(CVor_dist,CVor_full,n_s_max,MPI_DEF_REAL,MPI_SUM,0, &
           &          MPI_COMM_WORLD,ierr)
      call MPI_Reduce(CHel_dist,CHel_full,n_s_max,MPI_DEF_REAL,MPI_SUM,0, &
           &          MPI_COMM_WORLD,ierr)
#else
      Egeos_full(:) =Egeos_dist(:)
      Eg_nAP_full(:)=Eg_nAP_dist(:)
      Ekin_full(:)  =Ekin_dist(:)
      Ek_nAP_full(:)=Ek_nAP_dist(:)
      Ek_NTC_full(:)=Ek_NTC_dist(:)
      Ek_STC_full(:)=Ek_STC_dist(:)
      CVZ_full(:)   =CVZ_dist(:)
      CVor_full(:)  =CVor_dist(:)
      CHel_full(:)  =CHel_dist(:)
#endif

      if ( rank == 0 ) then

         phiNorm=two*pi/n_phi_max

         !-- Cylindrical integration using Simpson's rule
         Egeos =phiNorm*simps(Egeos_full*cyl*h,cyl)
         Eg_nAP=phiNorm*simps(Eg_nAP_full*cyl*h,cyl)
         Ekin  =phiNorm*simps(Ekin_full*cyl*h,cyl)
         Ek_nAP=phiNorm*simps(Ek_nAP_full*cyl*h,cyl)

         EkNTC =phiNorm*simps(Ek_NTC_full*cyl*h,cyl)
         EkSTC =phiNorm*simps(Ek_STC_full*cyl*h,cyl)

         !-- Correlations outside tangent cylinder
         CVZ=phiNorm*simps(CVZ_full(1:n_s_otc)*cyl(1:n_s_otc)*h(1:n_s_otc),cyl(1:n_s_otc))
         CVZ=abs(CVZ)/vol_otc
         CVor=phiNorm*simps(CVor_full(1:n_s_otc)*cyl(1:n_s_otc)*h(1:n_s_otc), &
         &                  cyl(1:n_s_otc))
         CVor=abs(CVor)/vol_otc
         CHel=phiNorm*simps(CHel_full(1:n_s_otc)*cyl(1:n_s_otc)*h(1:n_s_otc), &
         &                  cyl(1:n_s_otc))
         CHel=abs(CHel)/vol_otc

         !-- Relative geostrophic energy
         if ( Ekin>0.0_cp ) then
            Geos = Egeos/Ekin
            EkN  = EkNTC/Ekin
            EkS  = EkSTC/Ekin
         else
            Geos = 0.0_cp
            EkN  = 0.0_cp
            EkS  = 0.0_cp
         end if

         !-- Relative geostrophic energy for non-axisymmetric flows perpendicular to
         !-- rotation axis
         if ( Ek_nAP>0.0_cp ) then
            GeosNAP = Eg_nAP/Ek_nAP
         else
            GeosNAP = 0.0_cp
         end if

         !-- Relative axisymmetric geostrophic energy
         if ( EA > 0.0_cp ) then
            GeosA = EgeosA/EA
         else
            GeosA = 0.0_cp
         end if

         !-- Relative meridional geostrophic energy
         if ( EM > 0.0_cp ) then
            GeosM = EgeosM/EM
         else
            GeosM = 0.0_cp
         end if

         !-- Relative zonal geostrophic energy
         if ( EZ > 0.0_cp ) then
            GeosZ = EgeosZ/EZ
         else
            GeosZ = 0.0_cp
         end if

         if ( l_save_out ) then
            open(newunit=n_geos_file, file=geos_file, status='unknown', &
            &    position='append')
         end if

         write(n_geos_file,'(1P,ES20.12,11ES16.8)') time, Geos, EkN, EkS, Ekin, CVZ, &
         &                                          CVor, CHel, GeosA, GeosZ, GeosM, &
         &                                          GeosNAP
         if ( l_save_out ) close(n_geos_file)
      end if

   end subroutine outGeos
!------------------------------------------------------------------------------------
   subroutine transp_R2Phi(arr_Rloc, arr_Ploc)
      !
      ! This subroutine is used to compute a MPI transpose between a R-distributed
      ! array and a Phi-distributed array
      !

      !-- Input array
      real(cp), intent(in) :: arr_Rloc(n_theta_max,n_phi_max,nRstart:nRstop)

      !-- Output array
      real(cp), intent(out) :: arr_Ploc(n_theta_max,n_r_max,nPstart:nPstop)

      !-- Local variables
      integer :: n_r, n_t, n_p
#ifdef WITH_MPI
      integer, allocatable :: rcounts(:), scounts(:), rdisp(:), sdisp(:)
      real(cp), allocatable :: sbuff(:), rbuff(:)
      integer :: p, ii, my_phi_counts

      !-- Set displacements vectors and buffer sizes
      allocate( rcounts(0:n_procs-1), scounts(0:n_procs-1) )
      allocate( rdisp(0:n_procs-1), sdisp(0:n_procs-1) )
      do p=0,n_procs-1
         my_phi_counts=phi_balance(p)%n_per_rank
         scounts(p)=nR_per_rank*my_phi_counts*n_theta_max
         rcounts(p)=radial_balance(p)%n_per_rank*(nPStop-nPStart+1)*n_theta_max
      end do

      rdisp(0)=0
      sdisp(0)=0
      do p=1,n_procs-1
         sdisp(p)=sdisp(p-1)+scounts(p-1)
         rdisp(p)=rdisp(p-1)+rcounts(p-1)
      end do
      allocate( sbuff(sum(scounts)), rbuff(sum(rcounts)) )

      !-- Prepare buffer
      do p=0,n_procs-1
         ii=sdisp(p)+1
         do n_r=nRstart,nRstop
            do n_p=phi_balance(p)%nStart,phi_balance(p)%nStop
               do n_t=1,n_theta_max
                  sbuff(ii)=arr_Rloc(n_t,n_p,n_r)
                  ii=ii+1
               end do
            end do
         end do
      end do

      !-- All to all
      call MPI_Alltoallv(sbuff, scounts, sdisp, MPI_DEF_REAL, &
           &             rbuff, rcounts, rdisp, MPI_DEF_REAL, &
           &             MPI_COMM_WORLD, ierr)

      !-- Reassemble array
      do p=0,n_procs-1
         ii=rdisp(p)+1
         do n_r=radial_balance(p)%nStart,radial_balance(p)%nStop
            do n_p=nPstart,nPstop
               do n_t=1,n_theta_max
                  arr_Ploc(n_t,n_r,n_p)=rbuff(ii)
                  ii=ii+1
               end do
            end do
         end do
      end do

      !-- Clear memory from temporary arrays
      deallocate( rcounts, scounts, rdisp, sdisp, rbuff, sbuff )
#else
      do n_r=1,n_r_max
         do n_p=1,n_phi_max
            do n_t=1,n_theta_max
               arr_Ploc(n_t,n_r,n_p)=arr_Rloc(n_t,n_p,n_r)
            end do
         end do
      end do
#endif

   end subroutine transp_R2Phi
!------------------------------------------------------------------------------------
   subroutine cylmean(dat, dat_OTC, dat_ITC_N, dat_ITC_S)
      !
      ! This routine computes the z-average inside and outside TC
      !

      !-- Input variables
      real(cp), intent(in) :: dat(n_theta_max,n_r_max) ! input data

      !-- Output variables
      real(cp), intent(out) :: dat_OTC(n_s_max)   ! z-average outside T.C.
      real(cp), intent(out) :: dat_ITC_N(n_s_max) ! z-average inside T.C. N.H.
      real(cp), intent(out) :: dat_ITC_S(n_s_max) ! z-average inside T.C. S.H.

      call cylmean_otc(dat,dat_OTC,n_s_max,n_s_otc,r,cyl,theta_ord)
      call cylmean_itc(dat,dat_ITC_N,dat_ITC_S,n_s_max,n_s_otc,r,cyl,theta_ord)

   end subroutine cylmean
!------------------------------------------------------------------------------------
   subroutine cylmean_otc(a,v,n_s_max,n_s_otc,r,s,theta)
      !
      ! This routine computes a z-averaging using Simpsons's rule outside T.C.
      !

      !-- Input variables
      integer,  intent(in) :: n_s_max
      integer,  intent(in) :: n_s_otc
      real(cp), intent(in) :: r(:)         ! Spherical radius
      real(cp), intent(in) :: s(n_s_max)   ! Cylindrical radius
      real(cp), intent(in) :: theta(:)     ! Colatitude
      real(cp), intent(in) :: a(:,:)

      !-- Output variable
      real(cp), intent(out) :: v(n_s_max)

      !-- Local variables
      integer :: n_z, nz, n_s, n_th, n_r, itr
      integer :: n_r0, n_r1, n_r2, n_r3, n_th0, n_th1, n_th2, n_th3
      real(cp) :: zmin, zmax, dz, z, eps, r_cmb, rc, thet
      real(cp) :: rr0, rr1, rr2, rr3, r10, r20, r30, r21, r31, r32
      real(cp) :: tt0, tt1, tt2, tt3, t10, t20, t30, t21, t31, t32
      real(cp) :: a01, a12, a23, a012, a123, tot
      real(cp) :: ac(-n_s_max:n_s_max,n_s_max), ait(0:3)

      eps=10.0_cp*epsilon(1.0_cp)
      r_cmb=r(1)
      v(:) =0.0_cp

      !-- Loop over axial cylinders starts here
      !$omp parallel do default(shared) &
      !$omp private(n_s,zmax,zmin,nz,dz,n_z,z,rc,thet,n_r2,n_r,n_r3,n_r1,n_r0) &
      !$omp private(n_th,n_th1,n_th2,n_th3,n_th0,rr0,rr1,rr2,rr3,r10,r20,r30)  &
      !$omp private(r21,r31,r32,tt0,tt1,tt2,tt3,t10,t20,t30,t21,t31,t32,itr)   &
      !$omp private(a01,a12,a23,a012,a123,tot,ait)
      sLoop: do n_s=1,n_s_otc

         zmax = sqrt(r_cmb*r_cmb-s(n_s)*s(n_s)) ! zmax
         zmin = 0.0_cp
         nz = 2*int(n_s_max*(zmax-zmin)/(2.0_cp*r_cmb)) ! Number of z points (one HS)
         nz = max(nz, 4)  ! Minimum to 4 for Simpson integration
         dz = (zmax-zmin)/real(nz,cp)

         !-- Loop over z starts
         do n_z=-nz,nz
            z=zmin+dz*n_z
            rc=sqrt(s(n_s)*s(n_s)+z*z)   ! radius from center
            if (rc >= r_cmb) rc=r_cmb-eps
            thet=0.5_cp*pi-atan(z/s(n_s))  ! polar angle of point (rax,z)
            ac(n_z,n_s)=0.0_cp
            !
            !  **** Interpolate values from (theta,r)-grid onto equidistant
            !  **** (z,rax)-grid using a fourth-order Lagrangian scheme
            !
            !--  Find indices of radial grid levels that bracket rc
            n_r2=n_r_max-1
            rbracket: do n_r=n_r_max-1,1,-1
               if ( r(n_r) >= rc ) then
                  n_r2 = n_r
                  exit rbracket
               end if
            end do rbracket
            if(n_r2 == n_r_max-1) n_r2=n_r_max-2
            if(n_r2 == 1 ) n_r2=2
            n_r3=n_r2-1
            n_r1=n_r2+1
            n_r0=n_r2+2

            !-- Find indices of angular grid levels that bracket thet
            n_th1=n_theta_max
            tbracket: do n_th=n_theta_max,1,-1
               if( theta(n_th) <= thet) then
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

            !--  Calculate differences in r for 4th-order interpolation
            rr0=rc-r(n_r0)
            rr1=rc-r(n_r1)
            rr2=rc-r(n_r2)
            rr3=rc-r(n_r3)
            r10= 1.0_cp/(r(n_r1)-r(n_r0))
            r20= 1.0_cp/(r(n_r2)-r(n_r0))
            r30= 1.0_cp/(r(n_r3)-r(n_r0))
            r21= 1.0_cp/(r(n_r2)-r(n_r1))
            r31= 1.0_cp/(r(n_r3)-r(n_r1))
            r32= 1.0_cp/(r(n_r3)-r(n_r2))

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

            !-- Loop over 4 neighboring grid angles
            do itr=0,3
               n_th=n_th0+itr

               !-- Interpolation in r-direction
               a01=(rr0*a(n_th,n_r1)-rr1*a(n_th,n_r0))*r10
               a12=(rr1*a(n_th,n_r2)-rr2*a(n_th,n_r1))*r21
               a23=(rr2*a(n_th,n_r3)-rr3*a(n_th,n_r2))*r32

               a012=(rr0*a12-rr2*a01)*r20
               a123=(rr1*a23-rr3*a12)*r31

               ait(itr)=(rr0*a123-rr3*a012)*r30
            end do

            !-- Interpolation in theta-direction
            a01=(tt0*ait(1)-tt1*ait(0))*t10
            a12=(tt1*ait(2)-tt2*ait(1))*t21
            a23=(tt2*ait(3)-tt3*ait(2))*t32

            a012=(tt0*a12-tt2*a01)*t20
            a123=(tt1*a23-tt3*a12)*t31
            ac(n_z,n_s)=ac(n_z,n_s)+(tt0*a123-tt3*a012)*t30
         end do

         !-- Simpson integration
         tot=ac(-nz,n_s)+ac(nz,n_s)
         do n_z=-nz+1,nz-1,2
            tot=tot+4.0_cp*ac(n_z,n_s)
         enddo
         do n_z=-nz+2,nz-2,2
            tot=tot+2.0_cp*ac(n_z,n_s)
         enddo
         v(n_s)=tot/(6.0_cp*nz)
      end do sLoop
      !$omp end parallel do

      !-- Equatorial point
      v(1)=0.5_cp*(a(n_theta_max/2,1)+a(n_theta_max/2+1,1))

   end subroutine cylmean_otc
!------------------------------------------------------------------------------------
   subroutine cylmean_itc(a,vn,vs,n_s_max,n_s_otc,r,s,theta)
      !
      ! This routine computes a z-averaging using Simpsons's rule inside T.C.
      !

      !-- Input variables
      integer,  intent(in) :: n_s_max
      integer,  intent(in) :: n_s_otc
      real(cp), intent(in) :: r(:)         ! Spherical radius
      real(cp), intent(in) :: s(n_s_max)         ! Cylindrical radius
      real(cp), intent(in) :: theta(:) ! Colatitude
      real(cp), intent(in) :: a(n_theta_max, n_r_max)

      !-- Output variable
      real(cp), intent(out) :: vn(n_s_max)
      real(cp), intent(out) :: vs(n_s_max)

      !-- Local variables
      integer :: n_z, nz, n_s, n_th, n_r, itr, n_hs
      integer :: n_r0, n_r1, n_r2, n_r3, n_th0, n_th1, n_th2, n_th3
      real(cp) :: zmin, zmax, dz, z, eps, r_cmb, r_icb, rc, thet
      real(cp) :: rr0, rr1, rr2, rr3, r10, r20, r30, r21, r31, r32
      real(cp) :: tt0, tt1, tt2, tt3, t10, t20, t30, t21, t31, t32
      real(cp) :: a01, a12, a23, a012, a123, tot1, tot2
      real(cp) :: ac(0:n_s_max,n_s_max,2), ait(0:3)

      eps=10.0_cp*epsilon(1.0_cp)
      r_cmb=r(1)
      r_icb=r(size(r))

      vn(:)=0.0_cp
      vs(:)=0.0_cp

      !-- Loop over axial cylinders starts here
      !$omp parallel do default(shared) &
      !$omp private(n_s,zmax,zmin,nz,dz,n_z,z,rc,thet,n_r2,n_r,n_r3,n_r1,n_r0) &
      !$omp private(n_th,n_th1,n_th2,n_th3,n_th0,rr0,rr1,rr2,rr3,r10,r20,r30)  &
      !$omp private(r21,r31,r32,tt0,tt1,tt2,tt3,t10,t20,t30,t21,t31,t32,itr)   &
      !$omp private(a01,a12,a23,a012,a123,tot1,tot2,n_hs,ait)
      sLoop: do n_s=n_s_otc+1,n_s_max

         zmax = sqrt(r_cmb*r_cmb-s(n_s)*s(n_s)) ! zmax
         zmin = sqrt(r_icb*r_icb-s(n_s)*s(n_s))
         nz = 2*int(n_s_max*(zmax-zmin)/(2.0_cp*r_cmb)) ! Number of z points (one HS)
         nz = max(nz, 4)  ! Minimum to 4 for Simpson integration
         dz = (zmax-zmin)/real(nz,cp)

         !-- Loop over z starts
         do n_z=0,nz
            z=zmin+dz*n_z
            rc=sqrt(s(n_s)*s(n_s)+z*z)   ! radius from center
            if (rc >= r_cmb) rc=r_cmb-eps
            if (rc <= r_icb) rc=r_icb+eps
            thet=0.5_cp*pi-atan(z/s(n_s))  ! polar angle of point (rax,z)
            ac(n_z,n_s,1)=0.0_cp
            ac(n_z,n_s,2)=0.0_cp
            !
            !  **** Interpolate values from (theta,r)-grid onto equidistant
            !  **** (z,rax)-grid using a fourth-order Lagrangian scheme
            !
            !--  Find indices of radial grid levels that bracket rc
            n_r2=n_r_max-1
            rbracket: do n_r=n_r_max-1,1,-1
               if ( r(n_r) >= rc ) then
                  n_r2 = n_r
                  exit rbracket
               end if
            end do rbracket
            if(n_r2 == n_r_max-1) n_r2=n_r_max-2
            if(n_r2 == 1) n_r2=2
            n_r3=n_r2-1
            n_r1=n_r2+1
            n_r0=n_r2+2

            !-- Find indices of angular grid levels that bracket thet
            n_th1=n_theta_max
            tbracket: do n_th=n_theta_max,1,-1
               if( theta(n_th) <= thet) then
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

            !--  Calculate differences in r for 4th-order interpolation
            rr0=rc-r(n_r0)
            rr1=rc-r(n_r1)
            rr2=rc-r(n_r2)
            rr3=rc-r(n_r3)
            r10= 1.0_cp/(r(n_r1)-r(n_r0))
            r20= 1.0_cp/(r(n_r2)-r(n_r0))
            r30= 1.0_cp/(r(n_r3)-r(n_r0))
            r21= 1.0_cp/(r(n_r2)-r(n_r1))
            r31= 1.0_cp/(r(n_r3)-r(n_r1))
            r32= 1.0_cp/(r(n_r3)-r(n_r2))

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

            !-- Loop over North/Sooth
            do n_hs=1,2

               !-- Loop over 4 neighboring grid angles
               do itr=0,3
                  n_th=n_th0+itr

                  if ( n_hs == 2 ) n_th = n_theta_max+1-n_th ! Southern hemisphere

                  !-- Interpolation in r-direction
                  a01=(rr0*a(n_th,n_r1)-rr1*a(n_th,n_r0))*r10
                  a12=(rr1*a(n_th,n_r2)-rr2*a(n_th,n_r1))*r21
                  a23=(rr2*a(n_th,n_r3)-rr3*a(n_th,n_r2))*r32

                  a012=(rr0*a12-rr2*a01)*r20
                  a123=(rr1*a23-rr3*a12)*r31

                  ait(itr)=(rr0*a123-rr3*a012)*r30
               end do

               !-- Interpolation in theta-direction
               a01=(tt0*ait(1)-tt1*ait(0))*t10
               a12=(tt1*ait(2)-tt2*ait(1))*t21
               a23=(tt2*ait(3)-tt3*ait(2))*t32

               a012=(tt0*a12-tt2*a01)*t20
               a123=(tt1*a23-tt3*a12)*t31
               ac(n_z,n_s,n_hs)=ac(n_z,n_s,n_hs)+(tt0*a123-tt3*a012)*t30
            end do

         end do

         !-- Simpson integration
         tot1=ac(0,n_s,1)+ac(nz,n_s,1)
         tot2=ac(0,n_s,2)+ac(nz,n_s,2)
         do n_z=1,nz-1,2
            tot1=tot1+4.0_cp*ac(n_z,n_s,1)
            tot2=tot2+4.0_cp*ac(n_z,n_s,2)
         enddo
         do n_z=2,nz-2,2
            tot1=tot1+2.0_cp*ac(n_z,n_s,1)
            tot2=tot2+2.0_cp*ac(n_z,n_s,2)
         enddo
         vn(n_s)=tot1/(3.0_cp*nz)
         vs(n_s)=tot2/(3.0_cp*nz)
      end do sLoop
      !$omp end parallel do

   end subroutine cylmean_itc
!------------------------------------------------------------------------------------
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
            if( theta(n_th) <= thet) then
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

   end subroutine interp_theta
!------------------------------------------------------------------------------------
end module geos
