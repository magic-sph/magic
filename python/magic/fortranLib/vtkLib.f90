subroutine pvts_scal(filename,radius,scalars,scalcodes,nfiles,minc,np,nt,nr,nscals)

   implicit none

   ! input variables
   integer :: nr,nt,np,nscals
   integer,           intent(in) :: minc
   character(len=64), intent(in) :: filename
   real(kind=4),      intent(in) :: radius(nr)
   real(kind=4),      intent(in) :: scalars(nscals,np,nt,nr)
   integer,           intent(in) :: scalcodes(nscals)
   integer,           intent(in) :: nfiles

   ! local variables
   character(len=20) :: scalnames(nscals)
   character(len=20) :: vecnames(0)
   integer :: veccodes(0)
   real(kind=4) :: xyz(3,np*(nt/nfiles+1)*nr)
   real(kind=4) :: scals(nscals,np*(nt/nfiles+1)*nr)
   real(kind=4) :: vecs(0,0,0)
   real(kind=4) :: theta(nt)
   real(kind=4) :: phi(np)
   real(kind=4) :: pi
   integer(kind=4) :: nnos
   character(len=200) :: buffer
   character(len=64) :: fName
   character(len=36) :: str1 
   character(len=1) :: lf 
   integer(kind=4) :: ivtk=10
   integer :: h,i,j,k,l,m,jinit,jend,iscal
   integer :: nvecs=0

   lf = char(10) ! line feed character

   pi = acos(-1.)
   nnos = nr*(nt/nfiles+1)*np
   nvecs = 0

   phi=(/(2.*pi/minc/(np-1)*(i-1),i=1,np)/)
   theta=(/(pi/(nt-1)*(i-1),i=1,nt)/)

   call getCode(scalcodes,scalnames,nscals,veccodes,vecnames,nvecs)

   open(unit=ivtk,file=trim(filename)//'.pvts',form='unformatted', access='stream', &
   &    convert='big_endian')

   buffer = '<?xml version="1.0"?>'//lf 
   write(ivtk) trim(buffer)
   buffer = '<VTKFile type="PStructuredGrid" version="0.1" &
    &   byte_order="BigEndian">'//lf
   write(ivtk) trim(buffer)
   write(str1(1:36), '(6i6)') 0, nr-1, 0, nt-1, 0, np-1
   buffer = '  <PStructuredGrid WholeExtent="'//str1//'">'//lf
   write(ivtk) trim(buffer)
   buffer = '      <PPointData>'
   write(ivtk) trim(buffer)
   do h=1,nscals
      buffer = '         <DataArray type="Float32" Name="'//trim(scalnames(h))//'"/>'//lf
      write(ivtk) trim(buffer)
   end do
   buffer = '      </PPointData>'//lf
   write(ivtk) trim(buffer)
   buffer = '      <PCellData>  </PCellData>'//lf
   write(ivtk) trim(buffer)
   buffer = '      <PPoints>'//lf   
   write(ivtk) trim(buffer)
   buffer = '        <DataArray type="Float32" Name="coordinates" &
   & NumberOfComponents="3"/>'//lf
   write(ivtk) trim(buffer)
   buffer = '      </PPoints>'//lf
   write(ivtk) trim(buffer)

   do m=1,nfiles
      write(str1(1:3), '(i0.3)') m-1
      !fName = '2ScalPart_'//trim(str1(1:3))//'.vts'
      fName = 'part_'//trim(str1(1:3))//'.vts'
      if (m==1) then
         jinit=1+(m-1)*nt/nfiles
         jend=m*nt/nfiles+1
      else
         jinit=(m-1)*nt/nfiles
         jend=m*nt/nfiles
      end if
      l=1
      do k=1,np
         do j=jinit,jend
            do i=1,nr
               xyz(1,l)=radius(i)*cos(phi(k))*sin(theta(j))
               xyz(2,l)=radius(i)*sin(phi(k))*sin(theta(j))
               xyz(3,l)=radius(i)*cos(theta(j))
  
               l=l+1
            enddo
         enddo
      enddo

      do iscal=1,nscals
         l=1
         do k=1,np
            do j=jinit,jend
               do i=1,nr
                  scals(iscal,l)=scalars(iscal,k,j,i)
                  l=l+1
               enddo
            enddo
         enddo
      enddo

      call partXMLVts(fName,xyz,scals,scalnames,nscals,vecs,vecnames,nvecs, &
             &        nnos,nr,jinit,jend,np)
      write(str1(1:36), '(6i6)') 0, nr-1, jinit-1, jend-1, 0, np-1
      buffer = '    <Piece Extent="'//str1//'" Source="' &
               //trim(fName)//'">'//lf
      write(ivtk) trim(buffer)
      buffer = '    </Piece>'//lf
      write(ivtk) trim(buffer)
   enddo
   buffer = '  </PStructuredGrid>'//lf  
   write(ivtk) trim(buffer)

   buffer = '</VTKFile>'//lf
   write(ivtk) trim(buffer)
   close(ivtk)

end subroutine pvts_scal
!-------------------------------------------------------------------------
subroutine pvts(filename,radius,vecr,vect,vecp,scalars,scalcodes,veccodes,&
                nfiles,minc,nvecs,nscals,np,nt,nr)

   implicit none

   ! input variables
   integer :: nr,nt,np,nscals,nvecs
   character(len=64), intent(in) :: filename
   integer,           intent(in) :: minc
   real(kind=4),      intent(in) :: radius(nr)
   integer,           intent(in) :: scalcodes(nscals)
   integer,           intent(in) :: veccodes(nvecs)
   real(kind=4),      intent(in) :: scalars(nscals,np,nt,nr)
   real(kind=4),      intent(in) :: vecr(nvecs,np,nt,nr)
   real(kind=4),      intent(in) :: vect(nvecs,np,nt,nr)
   real(kind=4),      intent(in) :: vecp(nvecs,np,nt,nr)
   integer,           intent(in) :: nfiles

   ! local variables
   character(len=20) :: scalnames(nscals), vecnames(nvecs)
   real(kind=4) :: xyz(3,np*(nt/nfiles+1)*nr)
   real(kind=4) :: scals(nscals,np*(nt/nfiles+1)*nr)
   real(kind=4) :: vecs(nvecs,3,np*(nt/nfiles+1)*nr)
   real(kind=4) :: theta(nt)
   real(kind=4) :: phi(np)
   real(kind=4) :: pi,vs
   integer(kind=4) :: nnos
   character(len=200) :: buffer
   character(len=64) :: fName
   character(len=36) :: str1 
   character(len=1) :: lf 
   integer(kind=4) :: ivtk=10
   integer :: h,i,j,k,l,m,jinit,jend,ivec,iscal

   lf = char(10) ! line feed character

   pi = acos(-1.)
   nnos = nr*(nt/nfiles+1)*np

   phi=(/(2.*pi/minc/(np-1)*(i-1),i=1,np)/)
   theta=(/(pi/(nt-1)*(i-1),i=1,nt)/)

   call getCode(scalcodes,scalnames,nscals,veccodes,vecnames,nvecs)

   open(unit=ivtk,file=trim(filename)//'.pvts',form='unformatted', access='stream', &
   &    convert='big_endian')

   buffer = '<?xml version="1.0"?>'//lf 
   write(ivtk) trim(buffer)
   buffer = '<VTKFile type="PStructuredGrid" version="0.1" &
    &   byte_order="BigEndian">'//lf
   write(ivtk) trim(buffer)
   write(str1(1:36), '(6i6)') 0, nr-1, 0, nt-1, 0, np-1
   buffer = '  <PStructuredGrid WholeExtent="'//str1//'">'//lf
   write(ivtk) trim(buffer)
   buffer = '      <PPointData>'
   write(ivtk) trim(buffer)
   do h=1,nscals
      buffer = '         <DataArray type="Float32" Name="'//trim(scalnames(h))//'"/>'//lf
      write(ivtk) trim(buffer)
   end do
   do h=1,nvecs
      buffer = '         <DataArray type="Float32" NumberOfComponents="3" Name="'//trim(vecnames(h))//'"/>'//lf
      write(ivtk) trim(buffer)
   end do
   buffer = '      </PPointData>'//lf
   write(ivtk) trim(buffer)
   buffer = '      <PCellData>  </PCellData>'//lf
   write(ivtk) trim(buffer)
   buffer = '      <PPoints>'//lf   
   write(ivtk) trim(buffer)
   buffer = '        <DataArray type="Float32" Name="coordinates" &
   & NumberOfComponents="3"/>'//lf
   write(ivtk) trim(buffer)
   buffer = '      </PPoints>'//lf
   write(ivtk) trim(buffer)

   do m=1,nfiles
      write(str1(1:3), '(i0.3)') m-1
      fName = 'vecpart_'//trim(str1(1:3))//'.vts'
      if (m==1) then
         jinit=1+(m-1)*nt/nfiles
         jend=m*nt/nfiles+1
      else
         jinit=(m-1)*nt/nfiles
         jend=m*nt/nfiles
      end if
      do ivec=1,nvecs
         l=1
         do k=1,np
            do j=jinit,jend
               do i=1,nr
                  if (ivec == 1) then
                     xyz(1,l)=radius(i)*cos(phi(k))*sin(theta(j))
                     xyz(2,l)=radius(i)*sin(phi(k))*sin(theta(j))
                     xyz(3,l)=radius(i)*cos(theta(j))
                  end if
      
                  vs = vecr(ivec,k,j,i)*sin(theta(j))+vect(ivec,k,j,i)*cos(theta(j))
                  vecs(ivec,1,l)=vs*cos(phi(k))-vecp(ivec,k,j,i)*sin(phi(k))
                  vecs(ivec,2,l)=vs*sin(phi(k))+vecp(ivec,k,j,i)*cos(phi(k))
                  vecs(ivec,3,l)=vecr(ivec,k,j,i)*cos(theta(j))-vect(ivec,k,j,i)*sin(theta(j))
      
                  ! enddo
                  l=l+1
               end do
            end do
         end do
      end do

      do iscal=1,nscals
         l=1
         do k=1,np
            do j=jinit,jend
               do i=1,nr
                  scals(iscal,l)=scalars(iscal,k,j,i)
                  l=l+1
               enddo
            enddo
         enddo
      enddo

      call partXMLVts(fName,xyz,scals,scalnames,nscals,vecs,vecnames,nvecs, &
             &        nnos,nr,jinit,jend,np)
      write(str1(1:36), '(6i6)') 0, nr-1, jinit-1, jend-1, 0, np-1
      buffer = '    <Piece Extent="'//str1//'" Source="' &
               //trim(fName)//'">'//lf
      write(ivtk) trim(buffer)
      buffer = '    </Piece>'//lf
      write(ivtk) trim(buffer)
   enddo
   buffer = '  </PStructuredGrid>'//lf  
   write(ivtk) trim(buffer)

   buffer = '</VTKFile>'//lf
   write(ivtk) trim(buffer)
   close(ivtk)

end subroutine pvts
!-------------------------------------------------------------------------
subroutine partXMLVts(filename,xyz,scals,scalnames,nscals,vecs,vecnames,nvecs,&
                 &    nnos,nr,thinit,thend,nphi)

   implicit none

   integer(kind=4) :: ioff
   integer(kind=4) :: nbytes_scal,nbytes_vec
   integer(kind=4) :: nnos,i,j,k,intSize,nscals,nvecs
   integer(kind=4) :: nr,thinit,thend,nphi

   character(len=*), intent(in) :: filename
   real(kind=4),     intent(in) :: xyz(3,nnos)
   real(kind=4),     intent(in) :: scals(nscals,nnos)
   real(kind=4),     intent(in) :: vecs(nvecs,3,nnos)
   character(len=*), intent(in) :: scalnames(nscals)
   character(len=*), intent(in) :: vecnames(nvecs)

   character(len=200) :: buffer
   character(len=1) :: lf 
   character(len=12) :: offset
   character(len=36) :: str1 
   integer(kind=4) :: ivtk=9
   real(kind=4) :: floatSize

   lf = char(10) ! line feed character

   nbytes_scal   = nnos*int(sizeof(floatSize),kind=4)
   nbytes_vec    = 3*nnos*int(sizeof(floatSize),kind=4)

   open(unit=ivtk,file=filename,form='unformatted', access='stream', &
   &    convert='big_endian')

   buffer = '<?xml version="1.0"?>'//lf 
   write(ivtk) trim(buffer)
   buffer = '<VTKFile type="StructuredGrid" version="0.1" &
    &   byte_order="BigEndian">'//lf
   write(ivtk) trim(buffer)
   write(str1(1:36), '(6i6)') 0, nr-1, thinit-1, thend-1, 0, nphi-1
   buffer = '  <StructuredGrid WholeExtent="'//str1//'">'//lf
   write(ivtk) trim(buffer)
   buffer = '    <Piece Extent="'//str1//'">'//lf
   write(ivtk) trim(buffer)
   buffer = '      <PointData> '//lf
   write(ivtk) trim(buffer)
   ioff = 0
   do i=1,nscals
      write(offset(1:12),'(i12)') ioff
      buffer = '         <DataArray type="Float32" Name="'&
      & //trim(scalnames(i))//&
      & '" format="appended" offset="'//offset(1:12)//'"       />'//lf
      write(ivtk) trim(buffer)
      ioff=ioff+int(sizeof(intSize),kind=4)+nbytes_scal
   end do
   do i=1,nvecs
      write(offset(1:12),'(i12)') ioff
      buffer = '         <DataArray type="Float32" Name="'&
      & //trim(vecnames(i))//&
      & '" NumberOfComponents="3" format="appended" offset="'//offset(1:12)//'"       />'//lf
      write(ivtk) trim(buffer)
      ioff=ioff+int(sizeof(intSize),kind=4)+nbytes_vec
   end do
   buffer = '      </PointData>'//lf
   write(ivtk) trim(buffer)
   buffer = '      <CellData>  </CellData>'//lf
   write(ivtk) trim(buffer)
   buffer = '      <Points>'//lf   
   write(ivtk) trim(buffer)
   write(offset(1:12),'(i12)') ioff
   buffer = '        <DataArray type="Float32" Name="coordinates" &
   & NumberOfComponents="3" format="appended" offset="            &
   & '//trim(offset(1:12))//'" />'//lf
   write(ivtk) trim(buffer)
   buffer = '      </Points>'//lf
   write(ivtk) trim(buffer)
   buffer = '    </Piece>'//lf
   write(ivtk) trim(buffer)
   buffer = '  </StructuredGrid>'//lf  
   write(ivtk) trim(buffer)
   buffer = '  <AppendedData encoding="raw">'//lf
   write(ivtk) trim(buffer)
   buffer = '_'
   write(ivtk) trim(buffer)
   do j=1,nscals
      write(ivtk) nbytes_scal, (scals(j,i),i=1,nnos)
   end do
   do k=1,nvecs
      write(ivtk) nbytes_vec, ((vecs(k,i,j),i=1,3),j=1,nnos)
   end do
   write(ivtk) nbytes_vec, ((xyz(i,j),i=1,3),j=1,nnos)
   buffer = lf//'  </AppendedData>'//lf
   write(ivtk) trim(buffer)
   buffer = '</VTKFile>'//lf
   write(ivtk) trim(buffer)

   close(ivtk)

end subroutine partXMLVts
!-------------------------------------------------------------------------
subroutine getCode(scalcodes,scalnames,nscals,veccodes,vecnames,nvecs)

   implicit none

   !-- Input variables
   integer,           intent(in) :: nscals
   integer,           intent(in) :: nvecs
   integer,           intent(in) :: scalcodes(nscals)
   integer,           intent(in) :: veccodes(nvecs)

   !-- Output variables
   character(len=20), intent(out) :: scalnames(nscals)
   character(len=20), intent(out) :: vecnames(nvecs)

   !-- Local variables
   integer :: i

   do i=1,nscals
      if ( scalcodes(i) == -1) scalnames(i) = 'Radius'
      if ( scalcodes(i) == 1) scalnames(i) = 'Entropy'
      if ( scalcodes(i) == 2) scalnames(i) = 'Magnetic energy'
      if ( scalcodes(i) == 3) scalnames(i) = 'z vorticity'
      if ( scalcodes(i) == 4) scalnames(i) = 'Radial velocity'
      if ( scalcodes(i) == 5) scalnames(i) = 'Zonal velocity'
      if ( scalcodes(i) == 6) scalnames(i) = 'Fluct. entropy'
      if ( scalcodes(i) == 7) scalnames(i) = 'Fluct. z vorticity'
      if ( scalcodes(i) == 8) scalnames(i) = 'Radial vorticity'
      if ( scalcodes(i) == 9) scalnames(i) = 'Fluct. temp.'
      if ( scalcodes(i) == 10) scalnames(i) = 'Kinetic energy'
      if ( scalcodes(i) == 11) scalnames(i) = 'Radial mag. field'
      if ( scalcodes(i) == 12) scalnames(i) = 'Cyl Radial Velocity'
      if ( scalcodes(i) == 13) scalnames(i) = 'Colatitude'
      if ( scalcodes(i) == 14) scalnames(i) = 'Composition'
      if ( scalcodes(i) == 15) scalnames(i) = 'Fluct. composition'
   end do

   do i=1,nvecs
      if ( veccodes(i) == 1) vecnames(i) = 'Velocity'
      if ( veccodes(i) == 2) vecnames(i) = 'Fluct. velocity'
      if ( veccodes(i) == 3) vecnames(i) = 'Magnetic field'
      if ( veccodes(i) == 4) vecnames(i) = 'Fluct. mag. field'
   end do

end subroutine getCode
!-------------------------------------------------------------------------
subroutine vts(filename,radius,vecr,vect,vecp,scalars,scalcodes,veccodes,&
                minc,nvecs,nscals,np,nt,nr)

   implicit none

   !-- Input variables
   integer :: nr,nt,np,nvecs,nscals
   character(len=64), intent(in) :: filename
   integer,           intent(in) :: minc
   real(kind=4),      intent(in) :: radius(nr)
   real(kind=4),      intent(in) :: vecr(nvecs,np,nt,nr)
   real(kind=4),      intent(in) :: vect(nvecs,np,nt,nr)
   real(kind=4),      intent(in) :: vecp(nvecs,np,nt,nr)
   real(kind=4),      intent(in) :: scalars(nscals,np,nt,nr)
   integer,           intent(in) :: scalcodes(nscals)
   integer,           intent(in) :: veccodes(nvecs)

   !-- Local variables
   character(len=20) :: scalnames(nscals),vecnames(nvecs)
   real(kind=4) :: xyz(3,np*nt*nr)
   real(kind=4) :: scals(nscals,np*nt*nr)
   real(kind=4) :: vecs(nvecs,3,np*nt*nr)
   real(kind=4) :: theta(nt)
   real(kind=4) :: phi(np)
   real(kind=4) :: pi,Bs
   integer(kind=8) :: nnos
   integer :: i,j,k,l,iscal,ivec

   pi = acos(-1.)
   nnos = nr*nt*np

   phi=(/(2.*pi/minc/(np-1)*(i-1),i=1,np)/)
   theta=(/(pi/(nt-1)*(i-1),i=1,nt)/)

   call getCode(scalcodes,scalnames,nscals,veccodes,vecnames,nvecs)

   do ivec=1,nvecs
      l=1
      do k=1,np
         do j=1,nt
            do i=1,nr
               if (ivec == 1) then ! Calculate the coordinates on the first time
                  xyz(1,l)=radius(i)*cos(phi(k))*sin(theta(j))
                  xyz(2,l)=radius(i)*sin(phi(k))*sin(theta(j))
                  xyz(3,l)=radius(i)*cos(theta(j))
               end if
     
               Bs = vecr(ivec,k,j,i)*sin(theta(j))+vect(ivec,k,j,i)*cos(theta(j))
               vecs(ivec,1,l)=Bs*cos(phi(k))-vecp(ivec,k,j,i)*sin(phi(k))
               vecs(ivec,2,l)=Bs*sin(phi(k))+vecp(ivec,k,j,i)*cos(phi(k))
               vecs(ivec,3,l)=vecr(ivec,k,j,i)*cos(theta(j))-vect(ivec,k,j,i)*sin(theta(j))
     
               l=l+1
            enddo
         enddo
      enddo
   enddo

   do iscal=1,nscals
      l=1
      do k=1,np
         do j=1,nt
            do i=1,nr
               scals(iscal,l)=scalars(iscal,k,j,i)
               l=l+1
            enddo
         enddo
      enddo
   enddo

   call WriteXMLFormat(filename,xyz,scals,scalnames,vecs, &
                       vecnames,nnos,nscals,nvecs,nr,nt,np)

end subroutine vts
!-------------------------------------------------------------------------
subroutine vts_scal(filename,radius,scalars,scalcodes,minc,nscals,np,nt,nr)

   implicit none

   !-- Input variables
   integer :: nr,nt,np,nscals
   character(len=64), intent(in) :: filename
   integer,           intent(in) :: minc
   real(kind=4),      intent(in) :: radius(nr)
   real(kind=4),      intent(in) :: scalars(nscals,np,nt,nr)
   integer,           intent(in) :: scalcodes(nscals)

   !-- Local variables
   integer :: nvecs=0
   integer :: veccodes(0)
   character(len=20) :: scalnames(nscals)
   character(len=20) :: vecnames(0)
   real(kind=4) :: xyz(3,np*nt*nr)
   real(kind=4) :: scals(nscals,np*nt*nr)
   real(kind=4) :: vecs(0,0,0)
   real(kind=4) :: theta(nt)
   real(kind=4) :: phi(np)
   real(kind=4) :: pi
   integer(kind=8) :: nnos
   integer :: i,j,k,l,iscal

   pi = acos(-1.)
   nnos = nr*nt*np

   phi=(/(2.*pi/minc/(np-1)*(i-1),i=1,np)/)
   theta=(/(pi/(nt-1)*(i-1),i=1,nt)/)

   call getCode(scalcodes,scalnames,nscals,veccodes,vecnames,nvecs)

   l=1
   do k=1,np
      do j=1,nt
         do i=1,nr
            xyz(1,l)=radius(i)*cos(phi(k))*sin(theta(j))
            xyz(2,l)=radius(i)*sin(phi(k))*sin(theta(j))
            xyz(3,l)=radius(i)*cos(theta(j))
  
            l=l+1
         enddo
      enddo
   enddo

   do iscal=1,nscals
      l=1
      do k=1,np
         do j=1,nt
            do i=1,nr
               scals(iscal,l)=scalars(iscal,k,j,i)
               l=l+1
            enddo
         enddo
      enddo
   enddo

   call WriteXMLFormat(filename,xyz,scals,scalnames,vecs, &
                       vecnames,nnos,nscals,nvecs,nr,nt,np)

end subroutine vts_scal
!-------------------------------------------------------------------------
subroutine WriteXMLFormat(filename,xyz,scals,scalnames,vecs, &
                      vecnames,nnos,nscals,nvecs,nr,ntheta,nphi)

   implicit none

   integer(kind=8) :: ioff
   integer(kind=8) :: nbytes_scal,nbytes_vec
   integer(kind=8) :: nnos,i,j,k,intSize
   integer(kind=4) :: nr,ntheta,nphi
   integer(kind=4) :: nscals,nvecs

   character(len=*), intent(in) :: filename
   real(kind=4),     intent(in) :: xyz(3,nnos)
   real(kind=4),     intent(in) :: scals(nscals,nnos)
   real(kind=4),     intent(in) :: vecs(nvecs,3,nnos)
   character(len=*), intent(in) :: scalnames(nscals)
   character(len=*), intent(in) :: vecnames(nvecs)

   character(len=200) :: buffer
   character(len=1) :: lf 
   character(len=12) :: offset
   character(len=36) :: str1 
   integer(kind=4) :: ivtk=9
   real(kind=4) :: floatSize

   lf = char(10) ! line feed character

   nbytes_scal   = nnos * sizeof(floatSize)
   nbytes_vec    = 3  * nnos * sizeof(floatSize)

   open(unit=ivtk,file=trim(filename)//'.vts',form='unformatted', access='stream', &
   &    convert='big_endian')

   buffer = '<?xml version="1.0"?>'//lf 
   write(ivtk) trim(buffer)
   buffer = '<VTKFile type="StructuredGrid" version="0.1" &
    &   byte_order="BigEndian">'//lf
   write(ivtk) trim(buffer)
   write(str1(1:36), '(6i6)') 0, nr-1, 0, ntheta-1, 0, nphi-1
   buffer = '  <StructuredGrid WholeExtent="'//str1//'">'//lf
   write(ivtk) trim(buffer)
   buffer = '    <Piece Extent="'//str1//'">'//lf
   write(ivtk) trim(buffer)
   buffer = '      <PointData> '//lf
   write(ivtk) trim(buffer)
   ioff = 4
   do i=1,nscals
      write(offset(1:12),'(i12)') ioff
      buffer = '         <DataArray type="Float32" Name="'&
       & //trim(scalnames(i))//&
       & '" format="appended" offset="'//offset(1:12)//'"       />'//lf
      write(ivtk) trim(buffer)
      ioff=ioff+sizeof(intSize)+nbytes_scal
   enddo
   do i=1,nvecs
      write(offset(1:12),'(i12)') ioff
      buffer = '         <DataArray type="Float32" Name="'&
      & //trim(vecnames(i))//&
      &  '" NumberOfComponents="3" format="appended" offset="&
      & '//offset(1:12)//'" />'//lf
      write(ivtk) trim(buffer)
      ioff=ioff+sizeof(intSize)+nbytes_vec
   enddo
   buffer = '      </PointData>'//lf
   write(ivtk) trim(buffer)
   buffer = '      <CellData>  </CellData>'//lf
   write(ivtk) trim(buffer)
   buffer = '      <Points>'//lf   
   write(ivtk) trim(buffer)
   write(offset(1:12),'(i12)') ioff
   buffer = '        <DataArray type="Float32" Name="coordinates" &
   & NumberOfComponents="3" format="appended" offset="            &
   & '//trim(offset(1:12))//'" />'//lf
   write(ivtk) trim(buffer)
   buffer = '      </Points>'//lf
   write(ivtk) trim(buffer)
   buffer = '    </Piece>'//lf
   write(ivtk) trim(buffer)
   buffer = '  </StructuredGrid>'//lf  
   write(ivtk) trim(buffer)
   buffer = '  <AppendedData encoding="raw">'//lf
   write(ivtk) trim(buffer)
   buffer = '_'
   write(ivtk) trim(buffer)
   do j=1,nscals
      write(ivtk) nbytes_scal  , (scals(j,i),i=1,nnos)
   end do
   do k=1,nvecs
      write(ivtk) nbytes_vec   , ((vecs(k,i,j),i=1,3),j=1,nnos)
   end do
   write(ivtk) nbytes_vec   , ((xyz(i,j),i=1,3),j=1,nnos)
   buffer = lf//'  </AppendedData>'//lf
   write(ivtk) trim(buffer)
   buffer = '</VTKFile>'//lf
   write(ivtk) trim(buffer)

   close(ivtk)

end subroutine WriteXMLFormat
!-------------------------------------------------------------------------
subroutine vti(filename,vecx,vecy,vecz,scalars,scalcodes,veccodes,&
               minc,gridMax,spacng,nvecs,nscals,nz,ny,nx)

   implicit none

   !-- Input variables
   integer :: nx,ny,nz,nvecs,nscals
   character(len=64), intent(in) :: filename
   real(kind=4),      intent(in) :: gridMax,spacng
   integer,           intent(in) :: minc
   real(kind=4),      intent(in) :: vecx(nvecs,nz,ny,nx)
   real(kind=4),      intent(in) :: vecy(nvecs,nz,ny,nx)
   real(kind=4),      intent(in) :: vecz(nvecs,nz,ny,nx)
   real(kind=4),      intent(in) :: scalars(nscals,nz,ny,nx)
   integer,           intent(in) :: scalcodes(nscals)
   integer,           intent(in) :: veccodes(nvecs)

   !-- Local variables
   character(len=20) :: scalnames(nscals),vecnames(nvecs)
   real(kind=4) :: scals(nscals,nz*ny*nx)
   real(kind=4) :: vecs(nvecs,3,nz*ny*nx)
   integer(kind=8) :: nnos
   integer :: i,j,k,l,iscal,ivec

   nnos = int(nx*ny*nz,kind=8)

   call getCode(scalcodes,scalnames,nscals,veccodes,vecnames,nvecs)

   do ivec=1,nvecs
      l=1
      do k=1,nz
         do j=1,ny
            do i=1,nx
               vecs(ivec,1,l)=vecx(ivec,k,j,i)
               vecs(ivec,2,l)=vecy(ivec,k,j,i)
               vecs(ivec,3,l)=vecz(ivec,k,j,i)
     
               l=l+1
            enddo
         enddo
      enddo
   enddo

   do iscal=1,nscals
      l=1
      do k=1,nz
         do j=1,ny
            do i=1,nx
               scals(iscal,l)=scalars(iscal,k,j,i)
               l=l+1
            enddo
         enddo
      enddo
   enddo

   call WriteXmlVTI(filename,scals,scalnames,vecs, &
                    vecnames,nnos,nscals,nvecs,gridMax,spacng,nx,ny,nz)

end subroutine vti
!-------------------------------------------------------------------------
subroutine vti_scal(filename,scalars,scalcodes,minc,gridMax,spacng,nscals,nz,ny,nx)

   implicit none

   !-- Input variables
   integer :: nx,ny,nz,nscals
   character(len=64), intent(in) :: filename
   real(kind=4),      intent(in) :: gridMax,spacng
   integer,           intent(in) :: minc
   real(kind=4),      intent(in) :: scalars(nscals,nz,ny,nx)
   integer,           intent(in) :: scalcodes(nscals)

   !-- Local variables
   integer :: nvecs=0
   character(len=20) :: scalnames(nscals)
   real(kind=4) :: scals(nscals,nz*ny*nx)
   character(len=20) :: vecnames(0)
   integer :: veccodes(0)
   real(kind=4) :: vecs(0,0,0)
   integer(kind=8) :: nnos
   integer :: i,j,k,l,iscal

   nnos = nx*ny*nz

   call getCode(scalcodes,scalnames,nscals,veccodes,vecnames,nvecs)

   do iscal=1,nscals
      l=1
      do k=1,nz
         do j=1,ny
            do i=1,nx
               scals(iscal,l)=scalars(iscal,k,j,i)
               l=l+1
            enddo
         enddo
      enddo
   enddo

   call WriteXmlVTI(filename,scals,scalnames,vecs, &
                    vecnames,nnos,nscals,nvecs,gridMax,spacng,nx,ny,nz)

end subroutine vti_scal
!-------------------------------------------------------------------------
subroutine WriteXmlVTI(filename,scals,scalnames,vecs,vecnames, &
                       nnos,nscals,nvecs,gridMax,spacng,nx,ny,nz)

   implicit none

   !-- Input variables
   integer,          intent(in) :: nx,ny,nz
   real(kind=4),     intent(in) :: scals(nscals,nnos)
   real(kind=4),     intent(in) :: vecs(nvecs,3,nnos)
   real(kind=4),     intent(in) :: gridMax,spacng
   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: scalnames(nscals)
   character(len=*), intent(in) :: vecnames(nvecs)

   !-- Local variables
   integer(kind=8) :: nbytes_scal, nbytes_vec, nnos, ioff
   integer(kind=8) :: i, j, k
   integer :: nscals,nvecs
   character(len=200) :: buffer
   character(len=1) :: lf
   character(len=12) :: offset
   character(len=36) :: str1
   character(len=30) :: str2,str3
   integer   :: ivtk=9
   real(kind=4) :: floatSize

   lf = char(10) ! line feed character

   nbytes_scal   = nnos * sizeof(floatSize) 
   nbytes_vec    = 3  * nnos * sizeof(floatSize)

   open(unit=ivtk,file=trim(filename)//'.vti',form='unformatted',access='stream', &
   &    convert='big_endian')

   buffer = '<?xml version="1.0"?>'//lf 
   write(ivtk) trim(buffer)
   buffer = '<VTKFile type="ImageData" version="0.1" &
            & byte_order="BigEndian">'//lf
   write(ivtk) trim(buffer)
   write(str1(1:36), '(6i6)') 0, nx-1, 0, ny-1, 0, nz-1
   write(str2(1:30), '(3f9.5)') -gridMax,-gridMax,-gridMax
   write(str3(1:30), '(3f9.5)') spacng,spacng,spacng
   buffer = '  <ImageData WholeExtent="'//str1//'" Origin="'&
           &  //str2//'" Spacing="'//str3//'">'//lf
   write(ivtk) trim(buffer)
   buffer = '    <Piece Extent="'//str1//'">'//lf
   write(ivtk) trim(buffer)
   buffer = '      <PointData> '//lf
   write(ivtk) trim(buffer)
   ioff = 0
   do i=1,nscals
      write(offset(1:12),'(i12)') ioff
      buffer = '         <DataArray type="Float32" Name="'&
               & //trim(scalnames(i))//&
               & '" format="appended" offset="'//offset(1:12)//'"       />'//lf
      write(ivtk) trim(buffer)
      ioff=ioff+sizeof(ivtk)+nbytes_scal
   enddo
   do i=1,nvecs
      write(offset(1:12),'(i12)') ioff
      buffer = '         <DataArray type="Float32" Name="'&
               & //trim(vecnames(i))//&
               & '" NumberOfComponents="3" format="appended" offset="&
               & '//offset(1:12)//'" />'//lf
      write(ivtk) trim(buffer)
      ioff=ioff+sizeof(ivtk)+nbytes_vec
   enddo
   buffer = '      </PointData>'//lf
   write(ivtk) trim(buffer)
   buffer = '      <CellData>  </CellData>'//lf
   write(ivtk) trim(buffer)
   buffer = '    </Piece>'//lf
   write(ivtk) trim(buffer)
   buffer = '  </ImageData>'//lf  
   write(ivtk) trim(buffer)
   buffer = '  <AppendedData encoding="raw">'//lf
   write(ivtk) trim(buffer)
   buffer = '_'
   write(ivtk) trim(buffer)
   ! write scalars
   do j=1,nscals
      write(ivtk) nbytes_scal  , (scals(j,i),i=1,nnos)
   end do
   ! write vectors
   do k=1,nvecs
      write(ivtk) nbytes_vec   , ((vecs(k,i,j),i=1,3),j=1,nnos)
   end do
   buffer = lf//'  </AppendedData>'//lf
   write(ivtk) trim(buffer)
   buffer = '</VTKFile>'//lf
   write(ivtk) trim(buffer)

   close(ivtk)

end subroutine WriteXmlVTI
