subroutine pvts(filename,radius,scalin,np,nt,nr,nfiles,minc)

   implicit none

   ! input variables
   integer :: nr,nt,np,nscals,nvecs,h,i,j,k,l,m,jinit,jend

   integer,           intent(in) :: minc
   character(len=64), intent(in) :: filename
   real(kind=4),      intent(in) :: radius(nr)
   real(kind=4),      intent(in) :: scalin(np,nt,nr)
   integer,           intent(in) :: nfiles

   ! local variables
   character(len=20) :: scalnames(2)
   character(len=20) :: vecnames(0)
   real(kind=4) :: xyz(3,np*(nt/nfiles+1)*nr)
   real(kind=4) :: scals(2,np*(nt/nfiles+1)*nr)
   real(kind=4) :: vecs(0,0,0)
   !real(kind=4) :: scals(1,np*(nt/nfiles+1)*nr)
   real(kind=4) :: theta(nt)
   real(kind=4) :: phi(np)
   real(kind=4) :: pi
   integer(kind=4) :: nnos
   character(len=200) :: buffer
   character(len=64) :: fName
   character(len=36) :: str1 
   character(len=1) :: lf 
   integer(kind=4) :: ivtk=10

   lf = char(10) ! line feed character

   pi = acos(-1.)
   nnos = nr*(nt/nfiles+1)*np
   nscals = 2
   nvecs = 0

   phi=(/(2.*pi/minc/(np-1)*(i-1),i=1,np)/)
   theta=(/(pi/(nt-1)*(i-1),i=1,nt)/)

   scalnames(1)='Temperature'
   scalnames(2)='Radius'

   open(unit=ivtk,file=trim(filename)//'.pvts',form='unformatted', access='stream')
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
  
               ! do h=1,nscals
               scals(1,l)=scalin(k,j,i)
               scals(2,l)=radius(i)
               ! enddo
               l=l+1
            enddo
         enddo
      enddo
      !call partXMLVts(fName,xyz,scals,scalnames,nscals,nnos,nr,jinit,jend,np)
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
subroutine pvtsVec(filename,radius,vr,vt,vp,np,nt,nr,nfiles,minc)

   implicit none

   ! input variables
   integer :: nr,nt,np,nscals,nvecs,h,i,j,k,l,m,jinit,jend
   character(len=64), intent(in) :: filename
   integer,           intent(in) :: minc
   real(kind=4),      intent(in) :: radius(nr)
   real(kind=4),      intent(in) :: vr(np,nt,nr),vt(np,nt,nr),vp(np,nt,nr)
   integer,           intent(in) :: nfiles

   ! local variables
   character(len=20) :: scalnames(1), vecnames(1)
   real(kind=4) :: xyz(3,np*(nt/nfiles+1)*nr)
   real(kind=4) :: scals(1,np*(nt/nfiles+1)*nr)
   real(kind=4) :: vecs(1,3,np*(nt/nfiles+1)*nr)
   real(kind=4) :: theta(nt)
   real(kind=4) :: phi(np)
   real(kind=4) :: pi,vs
   integer(kind=4) :: nnos
   character(len=200) :: buffer
   character(len=64) :: fName
   character(len=36) :: str1 
   character(len=1) :: lf 
   integer(kind=4) :: ivtk=10

   lf = char(10) ! line feed character

   pi = acos(-1.)
   nnos = nr*(nt/nfiles+1)*np
   nscals = 1
   nvecs = 1

   phi=(/(2.*pi/minc/(np-1)*(i-1),i=1,np)/)
   theta=(/(pi/(nt-1)*(i-1),i=1,nt)/)

   scalnames(1)='Radius'
   vecnames(1)='Velocity'

   open(unit=ivtk,file=trim(filename)//'.pvts',form='unformatted', access='stream')
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
      l=1
      do k=1,np
         do j=jinit,jend
            do i=1,nr
               xyz(1,l)=radius(i)*cos(phi(k))*sin(theta(j))
               xyz(2,l)=radius(i)*sin(phi(k))*sin(theta(j))
               xyz(3,l)=radius(i)*cos(theta(j))
   
               ! do h=1,nscals
               scals(1,l)=radius(i)
   
               vs = vr(k,j,i)*sin(theta(j))+vt(k,j,i)*cos(theta(j))
               vecs(1,1,l)=vs*cos(phi(k))-vp(k,j,i)*sin(phi(k))
               vecs(1,2,l)=vs*sin(phi(k))+vp(k,j,i)*cos(phi(k))
               vecs(1,3,l)=vr(k,j,i)*cos(theta(j))-vt(k,j,i)*sin(theta(j))
   
               ! enddo
               l=l+1
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

end subroutine pvtsVec
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

   nbytes_scal   = nnos*sizeof(floatSize)
   nbytes_vec    = 3*nnos*sizeof(floatSize)

   open(unit=ivtk,file=filename,form='unformatted', access='stream')

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
      ioff=ioff+sizeof(intSize)+nbytes_scal
   end do
   do i=1,nvecs
      write(offset(1:12),'(i12)') ioff
      buffer = '         <DataArray type="Float32" Name="'&
      & //trim(vecnames(i))//&
      & '" NumberOfComponents="3" format="appended" offset="'//offset(1:12)//'"       />'//lf
      write(ivtk) trim(buffer)
      ioff=ioff+sizeof(intSize)+nbytes_vec
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
subroutine vts(filename,radius,Br,Bt,Bp,scal1,scal2,minc,np,nt,nr)

   implicit none

   ! input variables
   integer :: nr,nt,np,i,j,k,l
   character(len=64), intent(in) :: filename
   integer,           intent(in) :: minc
   real(kind=4),      intent(in) :: radius(nr)
   real(kind=4),      intent(in) :: Br(np,nt,nr),Bt(np,nt,nr),Bp(np,nt,nr)
   real(kind=4),      intent(in) :: scal1(np,nt,nr),scal2(np,nt,nr)

   ! local variables
   character(len=20) :: scalnames(4),vecnames(1)
   real(kind=4) :: xyz(3,np*nt*nr)
   real(kind=4) :: scals(4,np*nt*nr)
   real(kind=4) :: vecs(1,3,np*nt*nr)
   real(kind=4) :: theta(nt)
   real(kind=4) :: phi(np)
   real(kind=4) :: pi,Bs
   integer(kind=8) :: nnos

   pi = acos(-1.)
   nnos = nr*nt*np

   phi=(/(2.*pi/minc/(np-1)*(i-1),i=1,np)/)
   theta=(/(pi/(nt-1)*(i-1),i=1,nt)/)

   scalnames(1)='Radius'
   scalnames(2)='vr'
   scalnames(3)='vortr'
   scalnames(4)='entropy'
   vecnames(1)='B'
   l=1
   do k=1,np
      do j=1,nt
         do i=1,nr
            xyz(1,l)=radius(i)*cos(phi(k))*sin(theta(j))
            xyz(2,l)=radius(i)*sin(phi(k))*sin(theta(j))
            xyz(3,l)=radius(i)*cos(theta(j))
  
            Bs = Br(k,j,i)*sin(theta(j))+Bt(k,j,i)*cos(theta(j))
            vecs(1,1,l)=Bs*cos(phi(k))-Bp(k,j,i)*sin(phi(k))
            vecs(1,2,l)=Bs*sin(phi(k))+Bp(k,j,i)*cos(phi(k))
            vecs(1,3,l)= Br(k,j,i)*cos(theta(j))-Bt(k,j,i)*sin(theta(j))
  
            scals(1,l)=radius(i)
            scals(2,l)=Br(k,j,i)
            scals(3,l)=scal1(k,j,i)
            scals(4,l)=scal2(k,j,i)
  
            l=l+1
         enddo
      enddo
   enddo

   call WriteXMLFormat(filename,xyz,scals,scalnames,vecs, &
                       vecnames,nnos,4,1,nr,nt,np)

end subroutine vts
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

   open(unit=ivtk,file=trim(filename)//'.vts',form='unformatted', access='stream')

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
subroutine WriteVTR(filename,scals,vecs, &
                    nnos,nscals,nvecs,gridMax,spacng,nx,ny,nz)

   implicit none

   integer :: ioff
   integer :: nbytes_scal, nbytes_vec
   integer :: nnos
   integer :: nscals,nvecs,i,j,k

   integer,          intent(in) :: nx,ny,nz
   real(kind=4),     intent(in) :: scals(nscals,nnos)
   real(kind=4),     intent(in) :: vecs(nvecs,3,nnos)
   real(kind=4),     intent(in) :: gridMax,spacng
   character(len=*), intent(in) :: filename

   character(len=20) :: scalnames(nscals)
   character(len=20) :: vecnames(nvecs)
   character(len=200) :: buffer
   character(len=1) :: lf
   character(len=12) :: offset
   character(len=36) :: str1
   character(len=30) :: str2,str3
   integer   :: ivtk=9
   real(kind=4) :: floatSize

   lf = char(10) ! line feed character

   scalnames(1)='Radius'
   scalnames(2)='Br'
   scalnames(3)='Energy'
   vecnames(1) ='B'

   nbytes_scal   = int( nnos * sizeof(floatSize)     )
   nbytes_vec    = int( 3  * nnos * sizeof(floatSize))

   open(unit=ivtk,file=trim(filename)//'.vti',form='unformatted',access='stream')

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
      ioff=int(ioff+sizeof(ivtk)+nbytes_scal)
   enddo
   do i=1,nvecs
      write(offset(1:12),'(i12)') ioff
      buffer = '         <DataArray type="Float32" Name="'&
               & //trim(vecnames(i))//&
               & '" NumberOfComponents="3" format="appended" offset="&
               & '//offset(1:12)//'" />'//lf
      write(ivtk) trim(buffer)
      ioff=int(ioff+sizeof(ivtk)+nbytes_vec)
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

end subroutine WriteVTR
