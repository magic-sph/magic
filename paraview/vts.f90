!$Id$
!----------------------------------------------------------------------------
        program main

        implicit none

        integer   :: nscals,nvecs,mode
        integer   :: nr,nt,np,minc,npF,nnos
        integer   :: i,j,k,l=1
        real      :: time
        real      :: pi=acos(-1.d0), vs, Bs
        real*4, dimension(:), allocatable :: radius, colat
        real*4, dimension(:), allocatable :: theta, phi
        real*4, dimension(:,:,:), allocatable :: entropy,vr,vt,vp
        real*4, dimension(:,:,:), allocatable :: Br,Bt,Bp
        real*4, dimension(:,:), allocatable :: xyz
        real*4, dimension(:,:,:), allocatable :: vecs
        real*4, dimension(:,:), allocatable   :: scals
        character*20 :: scalnames(10),vecnames(10)
        character*20 ::filename
        character*2  ::str
        logical      :: vrSelected=.FALSE.,vtSelected=.FALSE.,&
                        vpSelected=.FALSE.,enSelected=.FALSE.,&
                        brSelected=.FALSE.,btSelected=.FALSE.,&
                        bpSelected=.FALSE.
        logical      :: velSelected=.FALSE.,magnSelected=.FALSE.
        integer      :: vrind=1,vtind=1,vpind=1,entind=1
        integer      :: brind=1,btind=1,bpind=1
        integer      :: velind=1,magnind=1

        namelist /input/ mode,filename,nscals,scalnames,nvecs,vecnames

        open(1,file='converter.nml',status='unknown')
        read(1,input)
        close(1)

        do i=1,nscals
           str = scalnames(i)
           if (str.eq.'vr') then
             vrSelected=.TRUE.
             vrind=i
           end if
           if (str.eq.'vt') then
             vtSelected=.TRUE.
             vtind=i
           end if
           if (str.eq.'vp') then
             vpSelected=.TRUE.
             vpind=i
           end if
           if (str.eq.'en') then
             enSelected=.TRUE.
             entind=i
           end if
           if (str.eq.'br') then
             brSelected=.TRUE.
             brind=i
           end if
           if (str.eq.'bt') then
             btSelected=.TRUE.
             btind=i
           end if
           if (str.eq.'bp') then
             bpSelected=.TRUE.
             bpind=i
           end if
        enddo

        do i=1,nvecs
           str = vecnames(i)
           if (str.eq.'ve') then
             velSelected=.TRUE.
             velind=i
           end if
           if (str.eq.'ma') then
             magnSelected=.TRUE.
             magnind=i
           end if
        enddo

        call header(filename,nr,nt,np,minc,time)

        if (int(np).eq.int(minc)*int(nt)) then
           np = np/minc
        end if

        npF = minc * np + 1 ! symmetry
        nnos = nr*npF*nt
        allocate(entropy(1:npF,1:nt,1:nr))
        allocate(vr(1:npF,1:nt,1:nr))
        allocate(vt(1:npF,1:nt,1:nr))
        allocate(vp(1:npF,1:nt,1:nr))
        allocate(Br(1:npF,1:nt,1:nr))
        allocate(Bt(1:npF,1:nt,1:nr))
        allocate(Bp(1:npF,1:nt,1:nr))
        allocate(radius(1:nr))
        allocate(colat(1:nt))
        if (mode.eq.0) then ! hydro case
          call readG(filename,entropy,vr,vt,vp,radius,colat,nr,nt,npF)
        else ! dynamo case
          call readGmagn(filename,entropy,vr,vt,vp,Br,Bt,Bp, &
                       radius,colat,nr,nt,npF)
        endif

        allocate(phi(1:npF))
        allocate(theta(1:nt))
        allocate(scals(1:nscals,1:nnos))
        allocate(vecs(1:nvecs,1:3,1:nnos))
        allocate(xyz(1:3,1:nnos))
        phi=(/(2.*pi/(npF-1)*(i-1),i=1,npF)/)
        theta=(/(pi/(nt-1)*(i-1),i=1,nt)/)

        radius = radius/maxval(radius)
        DO  k=1,npF
          DO j=1,nt
            DO i=1,nr
              xyz(1, l) = radius(i)*cos(phi(k))*sin(theta(j))
              xyz(2, l) = radius(i)*sin(phi(k))*sin(theta(j))
              xyz(3, l) = radius(i)*cos(theta(j))
              if (vrSelected) then
                 scals(vrind,l) = vr(k,j,i)
              end if
              if (vtSelected) then
                 scals(vtind,l) = vt(k,j,i)
              end if
              if (vpSelected) then
                 scals(vpind,l) = vp(k,j,i)
              end if
              if (enSelected) then
                 scals(entind,l) = entropy(k,j,i)
              end if
              if (brSelected) then
                 scals(brind,l) = br(k,j,i)
              end if
              if (btSelected) then
                 scals(btind,l) = bt(k,j,i)
              end if
              if (bpSelected) then
                 scals(bpind,l) = bp(k,j,i)
              end if
              if (velSelected) then
                 vs = vr(k,j,i)*sin(theta(j))+vt(k,j,i)*cos(theta(j))
                 vecs(velind,1,l)= vs*cos(phi(k))-vp(k,j,i)*sin(phi(k))
                 vecs(velind,2,l)= vs*sin(phi(k))+vp(k,j,i)*cos(phi(k))
                 vecs(velind,3,l)= vr(k,j,i)*cos(theta(j))- &
                                   vt(k,j,i)*sin(theta(j))
              end if
              if (magnSelected) then
                 Bs = Br(k,j,i)*sin(theta(j))+Bt(k,j,i)*cos(theta(j))
                 vecs(magnind,1,l)= Bs*cos(phi(k))-Bp(k,j,i)*sin(phi(k))
                 vecs(magnind,2,l)= Bs*sin(phi(k))+Bp(k,j,i)*cos(phi(k))
                 vecs(magnind,3,l)= Br(k,j,i)*cos(theta(j))- &
                                    Bt(k,j,i)*sin(theta(j))
              end if
              l = l+1
            END DO
          END DO
        END DO

        call WriteXMLFormat(filename,xyz,scals,scalnames,vecs,vecnames,&
                            nnos,nscals,nvecs,nr,nt,npF)

        end program
!----------------------------------------------------------------------------
        subroutine header(filename,nr,nt,np,minc,time)

        implicit none

        character*20 :: str1
        character*64 :: runid
        character(len=*), intent(in) :: filename
        integer, intent(out):: nr,nt,np,minc
        real, intent(out):: time
        real         :: nrf,ntf,npf,nric,mincf,nThetasBsf
        real         :: ra,ek,pr,prmag,radratio,sigma

        open(unit=11, file=filename,form='unformatted')
        read(11) str1
        read(11) runid
        read(11) time,nrf,ntf,npf,nric,mincf,nThetasBsf,ra,ek,pr,prmag, &
                 radratio,sigma
        nr = int(nrf)
        nt = int(ntf)
        np = int(npf)
        minc = int(mincf)

        close(11)
        end subroutine
!----------------------------------------------------------------------------
        subroutine readG(filename,entropy,vr,vt,vp,radius,colat,nri,nti,npF)

        implicit none

        character(len=*), intent(in) :: filename
        integer, intent(in)  :: nri,nti,npF
        integer      :: i,j, nth_loc
        character*20 :: str1
        character*64 :: runid
        real         :: time, nr,nt,np,nric,minc,nThetasBs
        real         :: ra,ek,pr,prmag,radratio,sigma
        real         :: ir, rad, ilat1, ilat2
        real, dimension(:,:), allocatable :: dummy
        real, dimension(:,:,:), allocatable :: entropyI,vrI,vtI,vpI
        real, dimension(1:nri), intent(out) :: radius
        real, dimension(1:nti), intent(out) :: colat
        real, dimension(1:npF,1:nti,1:nri), intent(out) :: entropy,vr,vt,vp

        open(unit=10, file=filename,form='unformatted')
        read(10) str1
        read(10) runid
        read(10) time,nr,nt,np,nric,minc,nThetasBs,ra,ek,pr,prmag, &
                 radratio,sigma
        read(10) colat

        if (int(np).eq.int(minc)*int(nt)) then
           np = np/minc
        end if

        allocate(dummy(1:int(np),1:int(nt)))
        allocate(entropyI(1:int(np),1:int(nt),1:int(nr)))
        allocate(vrI(1:int(np),1:int(nt),1:int(nr)))
        allocate(vtI(1:int(np),1:int(nt),1:int(nr)))
        allocate(vpI(1:int(np),1:int(nt),1:int(nr)))

        !reading
        do i=1,int(nr*nThetasBs)
           read(10) ir, rad, ilat1, ilat2
           radius(int(ir+1)) = rad
           nth_loc=int(ilat2)-int(ilat1)+1
           do j=int(ilat1),int(ilat2)
             read(10) entropyI(:,j,int(ir+1))
           end do
           do j=int(ilat1),int(ilat2)
             read(10) vrI(:,j,int(ir+1))
           end do
           do j=int(ilat1),int(ilat2)
             read(10) vtI(:,j,int(ir+1))
           end do
           do j=int(ilat1),int(ilat2)
             read(10) vpI(:,j,int(ir+1))
           end do
        end do

        close(10)

        !rearanging hemispherical data
        do i=1,int(nr)
           dummy = entropyI(:,:,i)
           do j=1,int(nt)/2
              entropyI(:,j,i)=dummy(:,2*(j-1)+1)
              entropyI(:,j+int(nt)/2,i)=dummy(:,int(nt)-1-2*(j-1)+1)
           end do
           dummy = vrI(:,:,i)
           do j=1,int(nt)/2
              vrI(:,j,i)=dummy(:,2*(j-1)+1)
              vrI(:,j+int(nt)/2,i)=dummy(:,int(nt)-1-2*(j-1)+1)
           end do
           dummy = vtI(:,:,i)
           do j=1,int(nt)/2
              vtI(:,j,i)=dummy(:,2*(j-1)+1)
              vtI(:,j+int(nt)/2,i)=dummy(:,int(nt)-1-2*(j-1)+1)
           end do
           dummy = vpI(:,:,i)
           do j=1,int(nt)/2
              vpI(:,j,i)=dummy(:,2*(j-1)+1)
              vpI(:,j+int(nt)/2,i)=dummy(:,int(nt)-1-2*(j-1)+1)
           end do
        end do

        do i=1,int(minc)
          entropy((i-1)*int(np)+1:i*int(np),:,:) = entropyI
          vr((i-1)*int(np)+1:i*int(np),:,:) = vrI
          vt((i-1)*int(np)+1:i*int(np),:,:) = vtI
          vp((i-1)*int(np)+1:i*int(np),:,:) = vpI
        end do

        deallocate(entropyI)
        deallocate(vrI)
        deallocate(vtI)
        deallocate(vpI)

        entropy(npF,:,:)=entropy(1,:,:)
        vr(npF,:,:)=vr(1,:,:)
        vt(npF,:,:)=vt(1,:,:)
        vp(npF,:,:)=vp(1,:,:)

        radius = radius/(1.-radratio)

        end subroutine
!----------------------------------------------------------------------------
        subroutine readGmagn(filename,entropy,vr,vt,vp,Br,Bt,Bp, &
                             radius,colat,nri,nti,npF)

        implicit none

        character(len=*), intent(in) :: filename
        integer, intent(in)  :: nri,nti,npF
        integer      :: i,j, nth_loc
        character*20 :: str1
        character*64 :: runid
        real         :: time, nr,nt,np,nric,minc,nThetasBs
        real         :: ra,ek,pr,prmag,radratio,sigma
        real         :: ir, rad, ilat1, ilat2
        real, dimension(:,:), allocatable :: dummy
        real, dimension(:,:,:), allocatable :: entropyI,vrI,vtI,vpI
        real, dimension(:,:,:), allocatable :: BrI,BtI,BpI
        real, dimension(1:nri), intent(out) :: radius
        real, dimension(1:nti), intent(out) :: colat
        real, dimension(1:npF,1:nti,1:nri), intent(out) :: entropy,vr,vt,vp
        real, dimension(1:npF,1:nti,1:nri), intent(out) :: Br,Bt,Bp

        open(unit=10, file=filename,form='unformatted')
        read(10) str1
        read(10) runid
        read(10) time,nr,nt,np,nric,minc,nThetasBs,ra,ek,pr,prmag, &
                 radratio,sigma
        read(10) colat

        if (int(np).eq.int(minc)*int(nt)) then
           np = np/minc
        end if

        allocate(dummy(1:int(np),1:int(nt)))
        allocate(entropyI(1:int(np),1:int(nt),1:int(nr)))
        allocate(vrI(1:int(np),1:int(nt),1:int(nr)))
        allocate(vtI(1:int(np),1:int(nt),1:int(nr)))
        allocate(vpI(1:int(np),1:int(nt),1:int(nr)))
        allocate(BrI(1:int(np),1:int(nt),1:int(nr)))
        allocate(BtI(1:int(np),1:int(nt),1:int(nr)))
        allocate(BpI(1:int(np),1:int(nt),1:int(nr)))

        !reading
        do i=1,int(nr*nThetasBs)
           read(10) ir, rad, ilat1, ilat2
           radius(int(ir+1)) = rad
           nth_loc=int(ilat2)-int(ilat1)+1
           do j=int(ilat1),int(ilat2)
             read(10) entropyI(:,j,int(ir+1))
           end do
           do j=int(ilat1),int(ilat2)
             read(10) vrI(:,j,int(ir+1))
           end do
           do j=int(ilat1),int(ilat2)
             read(10) vtI(:,j,int(ir+1))
           end do
           do j=int(ilat1),int(ilat2)
             read(10) vpI(:,j,int(ir+1))
           end do
           do j=int(ilat1),int(ilat2)
             read(10) BrI(:,j,int(ir+1))
           end do
           do j=int(ilat1),int(ilat2)
             read(10) BtI(:,j,int(ir+1))
           end do
           do j=int(ilat1),int(ilat2)
             read(10) BpI(:,j,int(ir+1))
           end do
        end do

        close(10)

        !rearanging hemispherical data
        do i=1,int(nr)
           dummy = entropyI(:,:,i)
           do j=1,int(nt)/2
              entropyI(:,j,i)=dummy(:,2*(j-1)+1)
              entropyI(:,j+int(nt)/2,i)=dummy(:,int(nt)-1-2*(j-1)+1)
           end do
           dummy = vrI(:,:,i)
           do j=1,int(nt)/2
              vrI(:,j,i)=dummy(:,2*(j-1)+1)
              vrI(:,j+int(nt)/2,i)=dummy(:,int(nt)-1-2*(j-1)+1)
           end do
           dummy = vtI(:,:,i)
           do j=1,int(nt)/2
              vtI(:,j,i)=dummy(:,2*(j-1)+1)
              vtI(:,j+int(nt)/2,i)=dummy(:,int(nt)-1-2*(j-1)+1)
           end do
           dummy = vpI(:,:,i)
           do j=1,int(nt)/2
              vpI(:,j,i)=dummy(:,2*(j-1)+1)
              vpI(:,j+int(nt)/2,i)=dummy(:,int(nt)-1-2*(j-1)+1)
           end do
           dummy = BrI(:,:,i)
           do j=1,int(nt)/2
              BrI(:,j,i)=dummy(:,2*(j-1)+1)
              BrI(:,j+int(nt)/2,i)=dummy(:,int(nt)-1-2*(j-1)+1)
           end do
           dummy = BtI(:,:,i)
           do j=1,int(nt)/2
              BtI(:,j,i)=dummy(:,2*(j-1)+1)
              BtI(:,j+int(nt)/2,i)=dummy(:,int(nt)-1-2*(j-1)+1)
           end do
           dummy = BpI(:,:,i)
           do j=1,int(nt)/2
              BpI(:,j,i)=dummy(:,2*(j-1)+1)
              BpI(:,j+int(nt)/2,i)=dummy(:,int(nt)-1-2*(j-1)+1)
           end do
        end do

        do i=1,int(minc)
          entropy((i-1)*int(np)+1:i*int(np),:,:) = entropyI
          vr((i-1)*int(np)+1:i*int(np),:,:) = vrI
          vt((i-1)*int(np)+1:i*int(np),:,:) = vtI
          vp((i-1)*int(np)+1:i*int(np),:,:) = vpI
          Br((i-1)*int(np)+1:i*int(np),:,:) = BrI
          Bt((i-1)*int(np)+1:i*int(np),:,:) = BtI
          Bp((i-1)*int(np)+1:i*int(np),:,:) = BpI
        end do

        deallocate(entropyI)
        deallocate(vrI)
        deallocate(vtI)
        deallocate(vpI)
        deallocate(BrI)
        deallocate(BtI)
        deallocate(BpI)

        entropy(npF,:,:)=entropy(1,:,:)
        vr(npF,:,:)=vr(1,:,:)
        vt(npF,:,:)=vt(1,:,:)
        vp(npF,:,:)=vp(1,:,:)
        Br(npF,:,:)=Br(1,:,:)
        Bt(npF,:,:)=Bt(1,:,:)
        Bp(npF,:,:)=Bp(1,:,:)

        radius = radius/(1.-radratio)

        end subroutine
!----------------------------------------------------------------------------
        subroutine WriteXMLFormat(filename,xyz,scals,scalnames,vecs, &
                          vecnames,nnos,nscals,nvecs,nr,ntheta,nphi)

        implicit none

        integer   :: ioff
        integer   :: nbytes_scal, nbytes_vec
        integer   :: nr, ntheta, nphi
        integer   :: nnos
        integer   :: nscals,nvecs,i,j,k
        character(len=*), intent(in) :: filename
        real*4, intent(in)    :: xyz (3,nnos)
        real*4, intent(in)    :: scals(nscals,nnos)
        real*4, intent(in)    :: vecs(nvecs,3,nnos)

        character(len=*), intent(in) :: scalnames(nscals)
        character(len=*), intent(in) :: vecnames(nvecs)
        character :: buffer*200, lf*1, offset*12, str1*24, str2*30
        integer   :: ivtk = 9, int
        real*4    :: float

        lf = char(10) ! line feed character

        nbytes_scal   = nnos * sizeof(float)
        nbytes_vec    = 3  * nnos * sizeof(float)

        open(unit=ivtk,file=trim(filename)//'.vts',form='binary')

        buffer = '<?xml version="1.0"?>'//lf 
        write(ivtk) trim(buffer)
        buffer = '<VTKFile type="StructuredGrid" version="0.1" &
        byte_order="BigEndian">'//lf
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
          //trim(scalnames(i))//&
          '" format="appended" offset="'//offset(1:12)//'"       />'//lf
          write(ivtk) trim(buffer)
          ioff=ioff+sizeof(int)+nbytes_scal
        enddo
        do i=1,nvecs
          write(offset(1:12),'(i12)') ioff
          buffer = '         <DataArray type="Float32" Name="'&
          //trim(vecnames(i))//&
          '" NumberOfComponents="3" format="appended" offset="&
          '//offset(1:12)//'" />'//lf
          write(ivtk) trim(buffer)
          ioff=ioff+sizeof(int)+nbytes_vec
        enddo
        buffer = '      </PointData>'//lf
        write(ivtk) trim(buffer)
        buffer = '      <CellData>  </CellData>'//lf
        write(ivtk) trim(buffer)
        buffer = '      <Points>'//lf   
        write(ivtk) trim(buffer)
        write(offset(1:12),'(i12)') ioff
        buffer = '        <DataArray type="Float32" Name="coordinates" &
        NumberOfComponents="3" format="appended" offset="&
        '//trim(offset(1:12))//'" />'//lf
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
!----------------------------------------------------------------------------
