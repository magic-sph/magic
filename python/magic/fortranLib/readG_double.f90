module greader_double

   implicit none

   real(kind=8) :: ra,ek,pr,prmag,radratio,sigma
   real(kind=8) :: time
   integer :: nr,nt,np,minc,nric,nThetasBs
   real(kind=8), allocatable :: radius(:),colat(:),radius_ic(:)
   real(kind=8), allocatable :: entropy(:,:,:),vr(:,:,:),vt(:,:,:),vp(:,:,:)
   real(kind=8), allocatable :: Br(:,:,:),Bt(:,:,:),Bp(:,:,:),pre(:,:,:)
   real(kind=8), allocatable :: Br_ic(:,:,:),Bt_ic(:,:,:),Bp_ic(:,:,:)
   real(kind=8), allocatable :: xi(:,:,:)

contains

   subroutine readG(filename,endian)

      !-- Input variables
      character(len=*), intent(in) :: filename
      character(len=1), intent(in) :: endian
   
      !-- Local variables
      integer :: i,j,nth_loc,nn,n_th,read_ok
      character(len=20) :: version
      character(len=64) :: runid
      real(kind=8) :: ir,rad,ilat1,ilat2,dumm
      real(kind=8) :: nrF,ntF,npF,mincF,nricF,nThetasBsF
      real(kind=8), allocatable :: dummy(:,:)
   
      if ( endian == 'B' ) then
         open(unit=10, file=filename, form='unformatted', convert='big_endian')
      else
         open(unit=10, file=filename, form='unformatted', convert='little_endian')
      end if

      read(10) version
      read(10) runid
      read(10) time,nrF,ntF,npF,nricF,mincF,nThetasBsF,ra,ek,pr,prmag, &
      &        radratio,sigma
   
      nr=int(nrF)
      nt=int(ntF)
      np=int(npF)
      minc=int(mincF)
      np=np/minc
      nThetasBs=int(nThetasBsF)
      nric=int(nricF)
   
      if ( allocated(colat) ) then
         deallocate( colat )
         deallocate( radius )
         deallocate( entropy )
         deallocate( vr )
         deallocate( vt )
         deallocate( vp )
         if ( version == 'Graphout_Version_8' .or. version == 'Graphout_Version_10') then
            deallocate( pre )
         end if
         if ( prmag /= 0. ) then
            deallocate( Br )
            deallocate( Bt )
            deallocate( Bp )
         end if
         if ( (prmag /= 0.) .and. (nric > 1) ) then
            deallocate( radius_ic )
            deallocate( Br_ic, Bt_ic, Bp_ic )
         end if
      end if
   
      allocate( colat(1:nt) )
      allocate( radius(1:nr) )
   
      allocate( dummy(1:np,1:nt) )
   
      allocate( entropy(1:np,1:nt,1:nr) )
      if ( version=='Graphout_Version_6' .or. version=='Graphout_Version_8'  &
          & .or. version=='Graphout_Version_10' .or. version=='Graphout_Version_12') then
         allocate( pre(1:np,1:nt,1:nr) )
      end if
      if ( version=='Graphout_Version_5' .or. version=='Graphout_Version_6' &
          & .or. version=='Graphout_Version_11'.or. version=='Graphout_Version_12') then
         allocate( xi(1:np,1:nt,1:nr) )
      end if
      allocate( vr(1:np,1:nt,1:nr) )
      allocate( vt(1:np,1:nt,1:nr) )
      allocate( vp(1:np,1:nt,1:nr) )
      if ( prmag /= 0. ) then
         allocate( Br(1:np,1:nt,1:nr) )
         allocate( Bt(1:np,1:nt,1:nr) )
         allocate( Bp(1:np,1:nt,1:nr) )
      end if

      if ( (prmag /= 0.) .and. (nric > 1) ) then
         allocate( radius_ic(1:nric) )
         allocate( Br_ic(1:np,1:nt,1:nric) )
         allocate( Bt_ic(1:np,1:nt,1:nric) )
         allocate( Bp_ic(1:np,1:nt,1:nric) )
      end if

      read(10) colat
   
      !reading
      if ( version=='Graphout_Version_9' .or. version=='Graphout_Version_10' &
           .or. version=='Graphout_Version_11' .or. version=='Graphout_Version_12') then
         do i=1,nr*nThetasBs
            read(10) ir, rad, ilat1, ilat2
            radius(int(ir)+1) = rad
            nth_loc=int(ilat2)-int(ilat1)+1
            read(10) entropy(:,int(ilat1):int(ilat2),int(ir+1))
            read(10) vr(:,int(ilat1):int(ilat2),int(ir+1))
            read(10) vt(:,int(ilat1):int(ilat2),int(ir+1))
            read(10) vp(:,int(ilat1):int(ilat2),int(ir+1))
            if ( version=='Graphout_Version_11' .or. version=='Graphout_Version_12') then
               read(10) xi(:,int(ilat1):int(ilat2),int(ir+1))
            end if
            if ( version=='Graphout_Version_10' .or. version=='Graphout_Version_12') then
               read(10) pre(:,int(ilat1):int(ilat2),int(ir+1))
            end if
            if ( prmag /= 0 ) then
               read(10) Br(:,int(ilat1):int(ilat2),int(ir+1))
               read(10) Bt(:,int(ilat1):int(ilat2),int(ir+1))
               read(10) Bp(:,int(ilat1):int(ilat2),int(ir+1))
            end if
         end do

         if ( (prmag /= 0.) .and. (nric > 1) ) then
            ic_loop1: do i=1,nric
               read(10, iostat=read_ok) ir, rad, ilat1, ilat2
               if ( read_ok /= 0 ) then
                  exit ic_loop1
               else
                  radius_ic(int(ir)+1-nr) = rad
                  read(10, iostat=read_ok) Br_ic(:,int(ilat1):int(ilat2),int(ir+1-nr))
                  if ( read_ok /= 0 ) exit ic_loop1
                  read(10, iostat=read_ok) Bt_ic(:,int(ilat1):int(ilat2),int(ir+1-nr))
                  if ( read_ok /= 0 ) exit ic_loop1
                  read(10, iostat=read_ok) Bp_ic(:,int(ilat1):int(ilat2),int(ir+1-nr))
                  if ( read_ok /= 0 ) exit ic_loop1
               end if
            end do ic_loop1
         end if
      else
         do i=1,nr*nThetasBs
            read(10) ir, rad, ilat1, ilat2
            radius(int(ir)+1) = rad
            nth_loc=int(ilat2)-int(ilat1)+1
            do j=int(ilat1),int(ilat2)
              read(10) entropy(:,j,int(ir+1))
            end do
            do j=int(ilat1),int(ilat2)
              read(10) vr(:,j,int(ir+1))
            end do
            do j=int(ilat1),int(ilat2)
              read(10) vt(:,j,int(ir+1))
            end do
            do j=int(ilat1),int(ilat2)
              read(10) vp(:,j,int(ir+1))
            end do
            if ( version=='Graphout_Version_5' .or. version=='Graphout_Version_6' ) then
               do j=int(ilat1),int(ilat2)
                 read(10) xi(:,j,int(ir+1))
               end do
            end if
            if ( version=='Graphout_Version_6' .or. version=='Graphout_Version_8' ) then
               do j=int(ilat1),int(ilat2)
                 read(10) pre(:,j,int(ir+1))
               end do
            end if
            if ( prmag /= 0 ) then
               do j=int(ilat1),int(ilat2)
                  read(10) Br(:,j,int(ir+1))
               end do
               do j=int(ilat1),int(ilat2)
                  read(10) Bt(:,j,int(ir+1))
               end do
               do j=int(ilat1),int(ilat2)
                  read(10) Bp(:,j,int(ir+1))
               end do
            end if
         end do

         if ( (prmag /= 0.) .and. (nric > 1) ) then
            ic_loop: do i=1,nric
               read(10, iostat=read_ok) ir, rad, ilat1, ilat2
               if ( read_ok /= 0 ) then
                  exit ic_loop
               else
                  radius_ic(int(ir)+1-nr) = rad
                  do j=int(ilat1),int(ilat2)
                     read(10) Br_ic(:,j,int(ir+1-nr))
                  end do
                  do j=int(ilat1),int(ilat2)
                     read(10) Bt_ic(:,j,int(ir+1-nr))
                  end do
                  do j=int(ilat1),int(ilat2)
                     read(10) Bp_ic(:,j,int(ir+1-nr))
                  end do
               end if
            end do ic_loop
         end if
   
      end if
   
      close(10)
   
      !rearanging hemispherical data
      do i=1,nr
         dummy(:,:) = entropy(:,:,i)
         do j=1,nt/2
            entropy(:,j,i)=dummy(:,2*(j-1)+1)
            entropy(:,j+nt/2,i)=dummy(:,nt-1-2*(j-1)+1)
         end do
         dummy(:,:) = vr(:,:,i)
         do j=1,nt/2
            vr(:,j,i)=dummy(:,2*(j-1)+1)
            vr(:,j+nt/2,i)=dummy(:,nt-1-2*(j-1)+1)
         end do
         dummy(:,:) = vt(:,:,i)
         do j=1,nt/2
            vt(:,j,i)=dummy(:,2*(j-1)+1)
            vt(:,j+nt/2,i)=dummy(:,nt-1-2*(j-1)+1)
         end do
         dummy(:,:) = vp(:,:,i)
         do j=1,nt/2
            vp(:,j,i)=dummy(:,2*(j-1)+1)
            vp(:,j+nt/2,i)=dummy(:,nt-1-2*(j-1)+1)
         end do
         if ( version=='Graphout_Version_11' .or. version=='Graphout_Version_12'&
              .or. version=='Graphout_Version_5' .or. version=='Graphout_Version_6' ) then
            dummy(:,:) = xi(:,:,i)
            do j=1,nt/2
               xi(:,j,i)=dummy(:,2*(j-1)+1)
               xi(:,j+nt/2,i)=dummy(:,nt-1-2*(j-1)+1)
            end do
         end if
         if ( version=='Graphout_Version_10' .or. version=='Graphout_Version_8'&
            & .or. version=='Graphout_Version_12' .or. version=='Graphout_Version_6' ) then
            dummy(:,:) = pre(:,:,i)
            do j=1,nt/2
               pre(:,j,i)=dummy(:,2*(j-1)+1)
               pre(:,j+nt/2,i)=dummy(:,nt-1-2*(j-1)+1)
            end do
         end if
         if ( prmag /= 0 ) then
            dummy(:,:) = Br(:,:,i)
            do j=1,nt/2
               Br(:,j,i)=dummy(:,2*(j-1)+1)
               Br(:,j+nt/2,i)=dummy(:,nt-1-2*(j-1)+1)
            end do
            dummy(:,:) = Bt(:,:,i)
            do j=1,nt/2
               Bt(:,j,i)=dummy(:,2*(j-1)+1)
               Bt(:,j+nt/2,i)=dummy(:,nt-1-2*(j-1)+1)
            end do
            dummy(:,:) = Bp(:,:,i)
            do j=1,nt/2
               Bp(:,j,i)=dummy(:,2*(j-1)+1)
               Bp(:,j+nt/2,i)=dummy(:,nt-1-2*(j-1)+1)
            end do
         end if
      end do

      !rearanging hemispherical data
      if ( (prmag /= 0) .and. (nric > 1)  .and. (read_ok == 0) ) then
         do i=1,nric
            dummy(:,:)= Br_ic(:,:,i)
            do j=1,nt/2
               Br_ic(:,j,i)     =dummy(:,2*(j-1)+1)
               Br_ic(:,j+nt/2,i)=dummy(:,nt-1-2*(j-1)+1)
            end do
            dummy(:,:) = Bt_ic(:,:,i)
            do j=1,nt/2
               Bt_ic(:,j,i)     =dummy(:,2*(j-1)+1)
               Bt_ic(:,j+nt/2,i)=dummy(:,nt-1-2*(j-1)+1)
            end do
            dummy(:,:) = Bp_ic(:,:,i)
            do j=1,nt/2
               Bp_ic(:,j,i)     =dummy(:,2*(j-1)+1)
               Bp_ic(:,j+nt/2,i)=dummy(:,nt-1-2*(j-1)+1)
            end do
         end do
      end if

      deallocate(dummy)
   
      radius(:) = radius(:)/(1.-radratio)
      if ( (prmag /= 0.) .and. (nric > 1) .and. (read_ok==0) ) then
         radius_ic(:) = radius_ic(:)/(1.-radratio)
      end if
   
   end subroutine readG
!----------------------------------------------------------------------------
end module greader_double
