    subroutine vts(filename, radius,Br,Bt,Bp,np,nt,nr)

    implicit none

    ! input variables
    character(len=64),intent(in) :: filename
    integer :: nr,nt,np
    real(kind=4),dimension(nr),intent(in) :: radius
    real(kind=4),dimension(np,nt,nr),intent(in) :: Br,Bt,Bp

    ! local variables
    character(len=20) :: scalnames(3),vecnames(1)
    real(kind=4),dimension(3,np*nt*nr) :: xyz,scals
    real(kind=4),dimension(1,3,np*nt*nr) :: vecs
    real(kind=4),dimension(nt) :: theta
    real(kind=4),dimension(np) :: phi
    real(kind=4) :: pi,Bs
    integer :: nnos,i,j,k,l

    pi = acos(-1.)
    nnos = nr*nt*np

    phi=(/(2.*pi/(np-1)*(i-1),i=1,np)/)
    theta=(/(pi/(nt-1)*(i-1),i=1,nt)/)

    scalnames(1)='Radius'
    scalnames(2)='Br'
    scalnames(3)='Energy'
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
        scals(3,l)=vecs(1,1,l)**2+vecs(1,2,l)**2+vecs(1,3,l)**2

        l=l+1
        enddo
      enddo
    enddo

    call WriteXMLFormat(filename,xyz,scals,scalnames,vecs, &
                        vecnames,nnos,3,1,nr,nt,np)

    end subroutine
!-------------------------------------------------------------------------
    subroutine WriteXMLFormat(filename,xyz,scals,scalnames,vecs, &
                      vecnames,nnos,nscals,nvecs,nr,ntheta,nphi)

    implicit none

    integer :: ioff
    integer :: nbytes_scal, nbytes_vec
    integer :: nr, ntheta, nphi
    integer :: nnos
    integer :: nscals,nvecs,i,j,k
    character(len=*),intent(in) :: filename
    real(kind=4),intent(in),dimension(3,nnos) :: xyz
    real(kind=4),intent(in),dimension(nscals,nnos) :: scals
    real(kind=4),intent(in),dimension(nvecs,3,nnos) :: vecs

    character(len=*),intent(in),dimension(nscals) :: scalnames
    character(len=*),intent(in),dimension(nvecs) :: vecnames
    character(len=200) :: buffer
    character(len=1) :: lf 
    character(len=12) :: offset
    character(len=24) :: str1 
    integer :: ivtk=9,intSize
    real(kind=4) :: floatSize

    lf = char(10) ! line feed character

    nbytes_scal   = int(nnos * sizeof(floatSize))
    nbytes_vec    = int(3  * nnos * sizeof(floatSize))

    open(unit=ivtk,file=trim(filename)//'.vts',form='unformatted', access='stream')

    buffer = '<?xml version="1.0"?>'//lf 
    write(ivtk) trim(buffer)
    buffer = '<VTKFile type="StructuredGrid" version="0.1" &
     &   byte_order="BigEndian">'//lf
    write(ivtk) trim(buffer)
    write(str1(1:24), '(6i4)') 0, nr-1, 0, ntheta-1, 0, nphi-1
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
      ioff=int(ioff+sizeof(intSize)+nbytes_scal)
    enddo
    do i=1,nvecs
      write(offset(1:12),'(i12)') ioff
      buffer = '         <DataArray type="Float32" Name="'&
      & //trim(vecnames(i))//&
      &  '" NumberOfComponents="3" format="appended" offset="&
      & '//offset(1:12)//'" />'//lf
      write(ivtk) trim(buffer)
      ioff=int(ioff+sizeof(intSize)+nbytes_vec)
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

    end subroutine
!-------------------------------------------------------------------------
    subroutine WriteVTR(filename,scals,vecs, &
                        nnos,nscals,nvecs,gridMax,spacng,nx,ny,nz)

    implicit none

    integer :: ioff
    integer :: nbytes_scal, nbytes_vec
    integer :: nnos
    integer :: nscals,nvecs,i,j,k
    integer,intent(in) :: nx,ny,nz
    real(kind=4),intent(in),dimension(nscals,nnos) :: scals
    real(kind=4),intent(in),dimension(nvecs,3,nnos) :: vecs
    real(kind=4),intent(in) :: gridMax,spacng
    character(len=*),intent(in) :: filename

    character(len=20),dimension(nscals) :: scalnames
    character(len=20),dimension(nvecs) :: vecnames
    character(len=200) :: buffer
    character(len=1) :: lf
    character(len=12) :: offset
    character(len=24) :: str1
    character(len=30) :: str2,str3
    integer   :: ivtk=9, intSize
    real(kind=4) :: floatSize

    lf = char(10) ! line feed character

    scalnames(1)='Radius'
    scalnames(2)='Br'
    scalnames(3)='Energy'
    vecnames(1)='B'

    nbytes_scal   = int( nnos * sizeof(floatSize)     )
    nbytes_vec    = int( 3  * nnos * sizeof(floatSize))

    open(unit=ivtk,file=trim(filename)//'.vti',form='unformatted',access='stream')

    buffer = '<?xml version="1.0"?>'//lf 
    write(ivtk) trim(buffer)
    buffer = '<VTKFile type="ImageData" version="0.1" &
             & byte_order="BigEndian">'//lf
    write(ivtk) trim(buffer)
    write(str1(1:24), '(6i4)') 0, nx-1, 0, ny-1, 0, nz-1
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
      ioff=int(ioff+sizeof(intSize)+nbytes_scal)
    enddo
    do i=1,nvecs
      write(offset(1:12),'(i12)') ioff
      buffer = '         <DataArray type="Float32" Name="'&
               & //trim(vecnames(i))//&
               & '" NumberOfComponents="3" format="appended" offset="&
               & '//offset(1:12)//'" />'//lf
      write(ivtk) trim(buffer)
      ioff=int(ioff+sizeof(intSize)+nbytes_vec)
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

    end subroutine
