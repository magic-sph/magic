module outMRI_mod

   use parallel_mod
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, n_r_maxStr, n_theta_maxStr, l_max, &
       &                 n_theta_max, n_phi_max, minc, lStressMem,   &
       &                 lm_max
   use radial_functions, only: r_ICB, rscheme_oc, r, r_CMB, orho1, rscheme_oc
   use radial_data, only: nRstart, nRstop, n_r_cmb, radial_balance
   use physical_parameters, only: ra, ek, pr, prmag, radratio, LFfac
   use mri
   use num_param, only: tScale
   use blocking, only: lo_map, llm, ulm
   use horizontal_data, only: phi, sinTheta, theta_ord, gauss
   use logic, only: lVerbose, l_save_out, l_MRICalc, l_mag
   use output_data, only: sDens, zDens, tag, log_file, runid, n_log_file
   use constants, only: pi, vol_oc, one, two, half, four
   use integration, only: rInt_R
   use useful, only: logWrite, abortRun
   use cosine_transform_odd
   use init_fields, only: q_rot, norm_ome
   use sht, only: spat_to_SH_axi
   use integration, only: simps,cylmean_otc,cylmean_itc
   
   implicit none

   private

!   integer :: nSmax, nZmaxA,ierrtest
!   integer :: nSstart, nSstop, nS_per_rank, nS_on_last_rank
!   real(cp) :: timeLast!,tNorm
!   real(outp) :: timeAve
   !-- s-distributed arrays

   !-- global arrays for outputs
!   integer, allocatable :: nZmaxS(:)
!   real(cp), allocatable :: zZ(:,:)

   !-- (l,r) Representation of the different contributions
   real(cp), allocatable :: cyl(:) ! Cylindrical grid                                              
   real(cp), allocatable :: h(:)   ! height                                     
   real(cp), allocatable :: Oh(:)  ! 1/h    
   real(cp), allocatable :: BpsAS(:,:)
   real(cp), allocatable :: VpsAS(:,:)
   real(cp), allocatable :: BpAS(:,:)
   real(cp), allocatable :: VpAS(:,:)
   real(cp), allocatable :: VsAS(:,:)
   real(cp), allocatable :: VzAS(:,:)
   real(cp), allocatable :: BsAS(:,:)
   real(cp), allocatable :: BzAS(:,:)
   real(cp), allocatable :: Bs2AS(:,:)
   real(cp), allocatable :: Bz2AS(:,:)
   real(cp), allocatable :: Bp2AS(:,:)
   real(cp), allocatable :: Vs2AS(:,:)
   real(cp), allocatable :: Vz2AS(:,:)
   real(cp), allocatable :: Vp2AS(:,:)
   integer :: n_s_otc, n_s_max
   real(cp) :: volcyl_oc
   
!   real(cp), allocatable :: EkinLMr_Rloc(:,:)
!   real(cp), allocatable :: EmagLMr_Rloc(:,:)
!   real(cp), allocatable :: EkinLMr(:,:)
!   real(cp), allocatable :: EmagLMr(:,:)
  

   !-- Output files
   character(len=64) :: MRIfilemean,MRI2DFile
   character(len=64) :: SMRIfileSprofile,NMRIfileSprofile
   public :: initialize_outMRI_mod, finalize_outMRI_mod, outMRI

 contains

   subroutine initialize_outMRI_mod

      !-- Local variables                                                         
      integer :: n_s, nFields, n
      real(cp) :: smin, smax, ds
      if ( lVerbose ) write(*,*) '! Before Initialisation outMRI!'

      if ( rank == 0 ) then 
         allocate( VpAS(n_theta_max,n_r_max), VsAS(n_theta_max,n_r_max) )
         allocate( VzAS(n_theta_max,n_r_max), BpAS(n_theta_max,n_r_max) )
         allocate( BsAS(n_theta_max,n_r_max), BzAS(n_theta_max,n_r_max) )
         allocate( Vp2AS(n_theta_max,n_r_max), Vz2AS(n_theta_max,n_r_max) )
         allocate( Bs2AS(n_theta_max,n_r_max), Vs2AS(n_theta_max,n_r_max) )
         allocate( Bz2AS(n_theta_max,n_r_max), Bp2AS(n_theta_max,n_r_max) )
         allocate( BpsAS(n_theta_max,n_r_max), VpsAS(n_theta_max,n_r_max) )
         bytes_allocated=bytes_allocated+14*n_theta_max*n_r_max*SIZEOF_DEF_REAL
      else
         allocate( VpAS(1,1), VzAS(1,1), VsAS(1,1), BsAS(1,1) )
         allocate( BzAS(1,1), BpAS(1,1), Vp2AS(1,1), Vz2AS(1,1) )
         allocate( Vs2AS(1,1), Bp2AS(1,1), Bs2AS(1,1), Bz2AS(1,1) )
         allocate( BpsAS(1,1), VpsAS(1,1))
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
      
      !-- R-distributed arrays
!      nSmax=n_r_max+int(r_ICB*real(n_r_max,cp))
!      nSmax=int(sDens*nSmax)
!      nZmaxA=2*nSmax+1
      !-- Determine the volume of the spherical shell interpolated on the cylindrical grid                          
      volcyl_oc = two * pi * (simps(h*cyl, cyl)+simps(h(n_s_otc+1:)*cyl(n_s_otc+1:), &
      &           cyl(n_s_otc+1:)))

      MRIfilemean='MaxReyavgMRI.'//tag
      NMRIfileSprofile='NMaxReySMRI.'//tag
      SMRIfileSprofile='SMaxReySMRI.'//tag
!      MRI2Dfile = 'BpsVps2DMRI.'//tag
      !movFile  ='MRImov.'//tag
    end subroutine initialize_outMRI_mod
!----------------------------------------------------------------------------
   subroutine finalize_outMRI_mod
      !
      ! Memory deallocation
     !
      deallocate(cyl, h, Oh)
      deallocate(VpsAS)!,EkinAS_Rloc,EmagAS_Rloc)
      deallocate(VpAS,Vp2AS,Vz2AS,Vs2AS,VzAS,VsAS)

      if (l_mag) then
         deallocate(BpsAS)!,EkinAS,EmagAS)
         deallocate(BpAS,Bp2AS,BzAS,BsAS,Bz2AS,Bs2AS)
     end if

!      deallocate( nZmaxS), zZ)
   end subroutine finalize_outMRI_mod
!----------------------------------------------------------------------------
   subroutine outMRI(time,n_time_step,nMRIsets)
      !
      !   Output of transport coefficients, Maxwell stress and Reynolds stress for the MRI.
      !   The slowest part in the MRI process is the repetitious calculation
      !   of plms by subroutine plm_theta. They are needed in getAStr and
      !   getPAStr when I transform on the cylindrical grid.
      !   The necessary plms could simply be calculated one and then
      !   be stored for later use!
      !
      !-- Input of variables:
      integer,          intent(in) :: n_time_step
      real(cp), intent(in)         :: time
      integer, intent(inout) ::nMRIsets

      !---- Local variables
      logical :: lStopRun,lTC
      !integer :: lm,l,m ! counter for degree and order
      integer :: nOutFile, nOutFile2, nOutFile3
      integer :: n_s,nSI, nS
      integer :: n      ! counter for theta blocks

      !-- Global arrays
      real(cp) :: BpsintN(n_s_max), VpsintN(n_s_max)
      real(cp) :: BsintN(n_s_max), VsintN(n_s_max)
      real(cp) :: BpintN(n_s_max), VpintN(n_s_max)
      real(cp) :: BzintN(n_s_max), VzintN(n_s_max)
      real(cp) :: Bp2intN(n_s_max), Vp2intN(n_s_max)
      real(cp) :: Bs2intN(n_s_max), Vs2intN(n_s_max)
      real(cp) :: Bz2intN(n_s_max), Vz2intN(n_s_max)
      real(cp) :: BpsintS(n_s_max), VpsintS(n_s_max)
      real(cp) :: BsintS(n_s_max), VsintS(n_s_max)
      real(cp) :: BpintS(n_s_max), VpintS(n_s_max)
      real(cp) :: BzintS(n_s_max), VzintS(n_s_max)
      real(cp) :: Bp2intS(n_s_max), Vp2intS(n_s_max)
      real(cp) :: Bs2intS(n_s_max), Vs2intS(n_s_max)
      real(cp) :: Bz2intS(n_s_max), Vz2intS(n_s_max)
      
      real(cp) :: zMax_global(n_s_max)
      real(cp):: PowerVol_fR(nRstart:nRstop)
      real(cp):: PowerVol_global(n_r_max)
      real(cp):: V0_th,V_temp
      real(cp) :: fvisc_nr_cmb
      real(cp) :: fvisc_nr_1,zMax,zMin
      
      !-- mean values
      real(cp) :: Maxwellstress, Reynoldsstress,Bs2,Bz2,Bp2,Vs2,Vz2,Vp2,Ekin,Emag,PowerVol
      real(cp) :: Bs,Bz,Bp,Vs,Vz,Vp

      !-- For boundaries:                  
!      real(cp) :: BpsB(2),BpsdB(2),BpsdB(2),zMin,zMax,dumm(8)
!      real(cp) :: Bs2B(2),BszB(2),BpzB(2),BzpdB(2),BpzdB(2)
      
      character(len=255) :: message
      character(len=64) :: version,fileName

#ifdef WITH_MPI
      real(cp) :: global_sum
      integer :: i,sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1)
#endif

      if ( lVerbose ) write(*,*) '! Starting outMRI!'

      l_MRICalc=.FALSE.

      nMRIsets=nMRIsets+1

      !if ( lVerbose ) write(*,*) '! Before LGTransform outMRI!'

      !-- Gather R-distributed arrays on rank == 0                                                  
      call gather_from_Rloc_to_rank0(BspAS_Rloc, BpsAS)
      call gather_from_Rloc_to_rank0(VspAS_Rloc, VpsAS)
      call gather_from_Rloc_to_rank0(BpAS_Rloc, BpAS)
      call gather_from_Rloc_to_rank0(VpAS_Rloc, VpAS)
      call gather_from_Rloc_to_rank0(BsAS_Rloc, BsAS)
      call gather_from_Rloc_to_rank0(VsAS_Rloc, VsAS)
      call gather_from_Rloc_to_rank0(BzAS_Rloc, BzAS)
      call gather_from_Rloc_to_rank0(VzAS_Rloc, VzAS)
      call gather_from_Rloc_to_rank0(Bp2AS_Rloc, Bp2AS)
      call gather_from_Rloc_to_rank0(Vp2AS_Rloc, Vp2AS)
      call gather_from_Rloc_to_rank0(Bs2AS_Rloc, Bs2AS)
      call gather_from_Rloc_to_rank0(Vs2AS_Rloc, Vs2AS)
      call gather_from_Rloc_to_rank0(Bz2AS_Rloc, Bz2AS)
      call gather_from_Rloc_to_rank0(Vz2AS_Rloc, Vz2AS)

      if (rank == 0) then
         call cylmean(VpAS, VpIntN, VpIntS)
         call cylmean(VsAS, VsIntN, VsIntS)
         call cylmean(VzAS, VzIntN, VzIntS)
         call cylmean(BpAS, BpIntN, BpIntS)
         call cylmean(BsAS, BsIntN, BsIntS)
         call cylmean(BzAS, BzIntN, BzIntS)
         call cylmean(Vp2AS, Vp2IntN, Vp2IntS)
         call cylmean(Vs2AS, Vs2IntN, Vs2IntS)
         call cylmean(Vz2AS, Vz2IntN, Vz2IntS)
         call cylmean(Bp2AS, Bp2IntN, Bp2IntS)
         call cylmean(Bs2AS, Bs2IntN, Bs2IntS)
         call cylmean(Bz2AS, Bz2IntN, Bz2IntS)
         call cylmean(VpsAS, VpsIntN, VpsIntS)
         call cylmean(BpsAS, BpsIntN, BpsIntS)
         
         if  ( nMRIsets == 1 ) then
            do n_s=1,n_s_max
               zMax = sqrt(r_CMB*r_CMB-cyl(n_s)*cyl(n_s))
               zMin = -zMax
               zMax_global(n_s)= zMax
            end do
         end if
 
         ! Integration finished
         
         if (lVerbose) write (*,*) '!Before Volume Force Work '
         
!   if (nRstop > n_r_cmb +2 .AND. nRstart < n_r_cmb+1) then
!      fvisc_nr_cmb=0.0_cp
!      fvisc_nr_1 = 0.0_cp
!      do nTheta=1,n_theta_maxStr ! Loop over theta blocks
!         nThetaNHS=(nTheta+1)/2
!         fvisc_nr_cmb=fvisc_nr_cmb+gauss(nThetaNHS)*fvisc_nr_cmbAS(nThetaNHS)
!         fvisc_nr_1=fvisc_nr_1+gauss(nThetaNHS)*fvisc_nr_1AS(nThetaNHS)
!      end do
!   end if

!   do nR=nRstart,nRstop
!      PowerVol_fR(nR)=0.0_cp
!      do n=1,n_theta_maxStr ! Loop over theta blocks
!         nThetaNHS=(nTheta+1)/2
!         ss=r(nR)*sinTheta(nTheta)/(0.4*r_cmb)
!         V0_th = r(nR)*norm_ome/(one+ss**(20.0_cp*q_rot))**0.05_cp - r(nR)/ek
!         V_temp= VzAS_Rloc(nThetaNHS,nR) + VpAS_Rloc(nThetaNHS,nR) + VsAS_Rloc(nThetaNHS,nR)
!         PowerVol_fR(nR)=PowerVol_fR(nR)-gauss(nThetaNHS)*tau*(VpAS_Rloc(nThetaNHS,nR)-V0_th)*V_temp
!      end do
!   end do


         if ( lVerbose ) write(*,*) '! Before Final output outMRI!'


         !--- Output of z-integral:
         open(newunit=nOutFile, file=NMRIfileSprofile, status='unknown',    &
              &       form='unformatted', position='append')
         if ( nMRIsets == 1 ) then
            
            write(message,'(" ! MRI: No. of s-values:",i4)') int(n_s_max/r_cmb)
            call logWrite(message)
            
            write(nOutFile) real(n_s_max,kind=outp)                        ! 1
            write(nOutFile) (real(cyl(nS),kind=outp),nS=1,n_s_max)          ! 2
            write(nOutFile) (real(zMax_global(nS),kind=outp),ns=1,n_s_max) ! 3 Zmax
         end if
         write(nOutFile)  real(time,kind=outp),                          &! 4
              &  (real(BpsIntN(nS)*LFfac,kind=outp)     ,nS=1,n_s_max),       &! 5  Bps
              &  (real(VpsIntN(nS),kind=outp)     ,nS=1,n_s_max),             &! 6  Vps
              &  (real(BsIntN(nS)*LFfac,kind=outp)     ,nS=1,n_s_max),        &! 7  Bs
              &  (real(VsIntN(nS),kind=outp)     ,nS=1,n_s_max),              &! 8  Vs
              &  (real(BpIntN(nS)*LFfac,kind=outp)     ,nS=1,n_s_max),        &! 9  Bp
              &  (real(VpIntN(nS),kind=outp)     ,nS=1,n_s_max),              &! 10  Vp
              &  (real(BzIntN(nS)*LFfac,kind=outp)     ,nS=1,n_s_max),        &! 11 Bz
              &  (real(VzIntN(nS),kind=outp)     ,nS=1,n_s_max),              &! 12 Vz
              &  (real(Bs2IntN(nS)*LFfac,kind=outp)     ,nS=1,n_s_max),       &! 13  Bs2
              &  (real(Vs2IntN(nS),kind=outp)     ,nS=1,n_s_max),             &! 14  Vs2
              &  (real(Bp2IntN(nS)*LFfac,kind=outp)     ,nS=1,n_s_max),       &! 15  Bp2
              &  (real(Vp2IntN(nS),kind=outp)     ,nS=1,n_s_max),             &! 16  Vp2
              &  (real(Bz2IntN(nS)*LFfac,kind=outp)     ,nS=1,n_s_max),       &! 17 Bz2
              &  (real(Vz2IntN(nS),kind=outp)     ,nS=1,n_s_max)               ! 18 Vz2

         
         !           &  (real(EkinS(nS),kind=outp)     ,nS=1,n_s_max),            &! 12 Ekin
         !           &  (real(EmagS(nS),kind=outp)     ,nS=1,n_s_max)              ! 13 Emag
         close(nOutFile)
         
         !--- Output of z-integral southern hemisphere: 
         open(newunit=nOutFile3, file=SMRIfileSprofile, status='unknown',    &
              &       form='unformatted', position='append')
         if ( nMRIsets == 1 ) then
            
            write(message,'(" ! MRI: No. of s-values:",i4)') int(n_s_max/r_cmb)
            call logWrite(message)

            write(nOutFile3) real(n_s_max,kind=outp)                        ! 1                                    
            write(nOutFile3) (real(cyl(nS),kind=outp),nS=1,n_s_max)          ! 2                                    
            write(nOutFile3) (real(zMax_global(nS),kind=outp),ns=1,n_s_max) ! 3 Zmax                               
         end if
         write(nOutFile3)  real(time,kind=outp),                          &! 4                                     
              &  (real(BpsIntS(nS)*LFfac,kind=outp)     ,nS=1,n_s_max),       &! 5  Bps                              
              &  (real(VpsIntS(nS),kind=outp)     ,nS=1,n_s_max),             &! 6  Vps                              
              &  (real(BsIntS(nS)*LFfac,kind=outp)     ,nS=1,n_s_max),        &! 7  Bs                               
              &  (real(VsIntS(nS),kind=outp)     ,nS=1,n_s_max),              &! 8  Vs                               
              &  (real(BpIntS(nS)*LFfac,kind=outp)     ,nS=1,n_s_max),        &! 9  Bp                               
              &  (real(VpIntS(nS),kind=outp)     ,nS=1,n_s_max),              &! 10  Vp                              
              &  (real(BzIntS(nS)*LFfac,kind=outp)     ,nS=1,n_s_max),        &! 11 Bz                               
              &  (real(VzIntS(nS),kind=outp)     ,nS=1,n_s_max),              &! 12 Vz                               
              &  (real(Bs2IntS(nS)*LFfac,kind=outp)     ,nS=1,n_s_max),       &! 13  Bs2                             
              &  (real(Vs2IntS(nS),kind=outp)     ,nS=1,n_s_max),             &! 14  Vs2                             
              &  (real(Bp2IntS(nS)*LFfac,kind=outp)     ,nS=1,n_s_max),       &! 15  Bp2                             
              &  (real(Vp2IntS(nS),kind=outp)     ,nS=1,n_s_max),             &! 16  Vp2                             
              &  (real(Bz2IntS(nS)*LFfac,kind=outp)     ,nS=1,n_s_max),       &! 17 Bz2                              
              &  (real(Vz2IntS(nS),kind=outp)     ,nS=1,n_s_max)               ! 18 Vz2                              
         !           &  (real(EkinS(nS),kind=outp)     ,nS=1,n_s_max),            &! 12 Ekin                           
         !           &  (real(EmagS(nS),kind=outp)     ,nS=1,n_s_max)              ! 13 Emag                          
         close(nOutFile3)
      
         !      PowerVol=rInt_R(PowerVol_global,r,rscheme_oc)

         ! S-Integration
         Maxwellstress = simps(BpsIntN*cyl*h, cyl)
         Maxwellstress = Maxwellstress + simps(BpsIntS(n_s_otc+1:n_s_max)*cyl(n_s_otc+1:n_s_max)* &
              &                        h(n_s_otc+1:n_s_max),cyl(n_s_otc+1:n_s_max))
         Reynoldsstress = simps(VpsIntN*cyl*h, cyl)
         Reynoldsstress = Reynoldsstress + simps(VpsIntS(n_s_otc+1:n_s_max)*cyl(n_s_otc+1:n_s_max)* &
              &                        h(n_s_otc+1:n_s_max),cyl(n_s_otc+1:n_s_max))
         Vs = simps(VsIntN*cyl*h, cyl)
         Vs = Vs+simps(VsIntS(n_s_otc+1:n_s_max)*cyl(n_s_otc+1:n_s_max)* &
              &                        h(n_s_otc+1:n_s_max),cyl(n_s_otc+1:n_s_max))
         Vp = simps(VpIntN*cyl*h, cyl)
         Vp = Vp+simps(VpIntS(n_s_otc+1:n_s_max)*cyl(n_s_otc+1:n_s_max)* &
              &                        h(n_s_otc+1:n_s_max),cyl(n_s_otc+1:n_s_max))
         Vz = simps(VzIntN*cyl*h, cyl)
         Vz = Vz+simps(VzIntS(n_s_otc+1:n_s_max)*cyl(n_s_otc+1:n_s_max)* &
              &                        h(n_s_otc+1:n_s_max),cyl(n_s_otc+1:n_s_max))
         Vs2 = simps(Vs2IntN*cyl*h, cyl)
         Vs2 = Vs2+simps(Vs2IntS(n_s_otc+1:n_s_max)*cyl(n_s_otc+1:n_s_max)* &
              &                        h(n_s_otc+1:n_s_max),cyl(n_s_otc+1:n_s_max))
         Vp2 = simps(Vp2IntN*cyl*h, cyl)
         Vp2 = Vp2+simps(Vp2IntS(n_s_otc+1:n_s_max)*cyl(n_s_otc+1:n_s_max)* &
              &                        h(n_s_otc+1:n_s_max),cyl(n_s_otc+1:n_s_max))
         Vz2 = simps(Vz2IntN*cyl*h, cyl)
         Vz2 = Vz2+simps(Vz2IntS(n_s_otc+1:n_s_max)*cyl(n_s_otc+1:n_s_max)* &
              &                        h(n_s_otc+1:n_s_max),cyl(n_s_otc+1:n_s_max))
         Bs = simps(BsIntN*cyl*h, cyl)
         Bs = Bs+simps(BsIntS(n_s_otc+1:n_s_max)*cyl(n_s_otc+1:n_s_max)* &
              &                        h(n_s_otc+1:n_s_max),cyl(n_s_otc+1:n_s_max))
         Bp = simps(BpIntN*cyl*h, cyl)
         Bp = Bp+simps(BpIntS(n_s_otc+1:n_s_max)*cyl(n_s_otc+1:n_s_max)* &
              &                        h(n_s_otc+1:n_s_max),cyl(n_s_otc+1:n_s_max))
         Bz = simps(BzIntN*cyl*h, cyl)
         Bz = Bz+simps(BzIntS(n_s_otc+1:n_s_max)*cyl(n_s_otc+1:n_s_max)* &
              &                        h(n_s_otc+1:n_s_max),cyl(n_s_otc+1:n_s_max))
         Bs2 = simps(Bs2IntN*cyl*h, cyl)
         Bs2 = Bs2+simps(Bs2IntS(n_s_otc+1:n_s_max)*cyl(n_s_otc+1:n_s_max)* &
              &                        h(n_s_otc+1:n_s_max),cyl(n_s_otc+1:n_s_max))
         Bp2 = simps(Bp2IntN*cyl*h, cyl)
         Bp2 = Bp2+simps(Bp2IntS(n_s_otc+1:n_s_max)*cyl(n_s_otc+1:n_s_max)* &
              &                        h(n_s_otc+1:n_s_max),cyl(n_s_otc+1:n_s_max))
         Bz2 = simps(Bz2IntN*cyl*h, cyl)
         Bz2 = Bz2+simps(Bz2IntS(n_s_otc+1:n_s_max)*cyl(n_s_otc+1:n_s_max)* &
              &                        h(n_s_otc+1:n_s_max),cyl(n_s_otc+1:n_s_max))
         
         Maxwellstress = two*pi*Maxwellstress/volcyl_oc
         Reynoldsstress = two*pi*Reynoldsstress/volcyl_oc
         Vs = two*pi*Vs/volcyl_oc
         Vp = two*pi*Vp/volcyl_oc
         Vz = two*pi*Vz/volcyl_oc
         Vs2 = two*pi*Vs2/volcyl_oc
         Vp2 = two*pi*Vp2/volcyl_oc
         Vz2 = two*pi*Vz2/volcyl_oc
         Bs = two*pi*Bs/volcyl_oc
         Bp = two*pi*Bp/volcyl_oc
         Bz = two*pi*Bz/volcyl_oc
         Bs2 = two*pi*Bs2/volcyl_oc
         Bp2 = two*pi*Bp2/volcyl_oc
         Bz2 = two*pi*Bz2/volcyl_oc
         Ekin =half*( Vs2 +Vp2 +Vz2)
         Emag = half*(Bs2 +Bp2 +Bz2)
         open(newunit=nOutFile2, file=MRIfilemean, status='unknown',      &
              &        position='append')
         write(nOutFile2,'(ES20.10,13ES15.7)')  real(time,kind=outp),     &! 1 time
              & real(Maxwellstress,kind=outp),                           &! 2 Maxwell stress
              & real(Reynoldsstress,kind=outp),                           &! 3 Reynolds stress
              & real(Vs2,kind=outp),                                     &! 4 Vs2
              & real(Vp2,kind=outp),                                      &! 5 Vp2
              & real(Vz2,kind=outp),                                     &! 6 Vz2
              & real(Bs2,kind=outp),                                      &! 7 Bs2
              & real(Bp2,kind=outp),                                     &! 8 Bp2
              & real(Bz2,kind=outp),                                      &! 9 Bz2
              & real(Ekin,kind=outp),                                    &! 10 Ekin
              & real(Emag,kind=outp)                                      ! 11 Emag
         !           & real(PowerVol,kind=outp),                                  &!Power of the volume force
         !           & real(fvisc_nr_cmb,kind = outp),                                        &!Visc Flux at the boundary
         !           & real(fvisc_nr_1,kind = outp)
         close(nOutFile2)
         
!      open(newunit=nOutFile3,file=MRI2DFile,status='unknown',    &
!           &       form='unformatted', position='append')
!      if ( nMRIsets == 1 ) then
!         write(nOutFile3) real(nZmaxA,kind=outp),real(n_s_max,kind=outp)   !  Numbers of elements
!         write(nOutFile3) (real(sZ(nS),kind=outp),nS=1,n_s_max)            !  Cylindrical radius
!         !      do nS=1,n_s_max
!         !        write(nOutFile3) (real(zZ(nZ,nS),kind=outp) ,nZ=1,nZmaxA)
!         !      end do
!      end if
!      write(nOutFile3) real(time,kind=outp)                              ! 1 time
!      do nS=1,n_s_max
!         if ( sZ(nS) < r_ICB ) then
!            lTC=.true.
!         else
!            lTC=.false.
!         end if
!
!         nZmax=nZmaxS(nS)
!         if ( lTC ) then
!            nZmaxNS=2*nZmax
!            do nZ=1,nZmax
!               zALL(nZ)=zZ(nZ,nS)
!               zALL(nZmaxNS-nZ+1)=-zZ(nZ,nS)
!            end do
!         else
!            nZmaxNS=nZmax
!            do nZ=1,nZmax
!               zALL(nZ)=zZ(nZ,nS)
!            end do
!         end if
!         write(nOutFile3) real(nZmaxNS,kind=outp)                          ! 2 Number of elements
!         write(nOutFile3) (real(zALL(nZ),kind=outp),nZ=1,nZmaxNS),           & ! 3 z grid
!              &           (real(VpsZS(nZ,nS),kind=outp) ,nZ=1,nZmaxNS),      & ! 4 Vps
!              &           (real(BpsZS(nZ,nS),kind=outp) ,nZ=1,nZmaxNS),      & ! 5 Bps
!              &           (real(VsZS(nZ,nS),kind=outp) ,nZ=1,nZmaxNS),       & ! 6 Vs
!              &           (real(BsZS(nZ,nS),kind=outp) ,nZ=1,nZmaxNS),       & ! 7 Bs
!              &           (real(VpZS(nZ,nS),kind=outp) ,nZ=1,nZmaxNS),       & ! 8 Vp
!              &           (real(BpZS(nZ,nS),kind=outp) ,nZ=1,nZmaxNS),       & ! 9 Bp
!              &           (real(VzZS(nZ,nS),kind=outp) ,nZ=1,nZmaxNS),       & ! 10 Vz
!              &           (real(BzZS(nZ,nS),kind=outp) ,nZ=1,nZmaxNS),       & ! 11 Bz
!              &           (real(EkinZS(nZ,nS),kind=outp) ,nZ=1,nZmaxNS),     & ! 12 kinetic energy
!              &           (real(EmagZS(nZ,nS),kind=outp) ,nZ=1,nZmaxNS)        ! 13 magnetic energy

!      end do
         close(nOutFile3)
      end if
  
    end subroutine outMRI

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

 !----------------------------------------------------------------------------
end module outMRI_mod
