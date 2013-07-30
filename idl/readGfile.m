%    Routine to load magic data into matlab
%    Works for binary data !
%    ok so let's scan the data first
%    Julien Aubert 11/02, jaubert@gwdg.de


% open the file
     fileName=uigetfile('G_*');
     fid  =fopen(fileName);

% get the various outputs

     version=char(fread(fid,32,'char')')
     tit=char(fread(fid,64,'char')')
     gridpar=fread(fid,9,'float');
     nrOc  =gridpar(4);
     nt    =gridpar(5);
     npI   =gridpar(6);
     nrIc  =gridpar(7);
     azsym =gridpar(8);
     nblock=gridpar(9);
     
     nr=nrOc+nrIc;

% Zero out input fields, define dimentions

     VrI=zeros(npI,nt,nrOc);
     VpI=zeros(npI,nt,nrOc);
     VtI=zeros(npI,nt,nrOc);
     BrI=zeros(npI,nt,nr);
     BpI=zeros(npI,nt,nr);
     BtI=zeros(npI,nt,nr);
     TI =zeros(npI,nt,nrOc);

% Get the physical parameters and the grid

     phypar=fread(fid,6,'float');
     dummy =fread(fid,2,'float');
     theta =fread(fid,nt,'float');
     cost  =cos(theta);
     sint  =sin(theta);

     r=zeros(nr,1);

% read the data for outer core

     for ir=1:nrOc*nblock

% Lines are separated from each other with two meaningless floats. I read
% them into the dummy variable

     dummy=fread(fid,2,'float');

     raddata=fread(fid,4,'float');
     rlevel=raddata(1);
     r(rlevel+1)=raddata(2);

     for i=raddata(3):raddata(4)
     dummy=fread(fid,2,'float');
     dummy=fread(fid,npI,'float');
     TI(:,i,rlevel+1)=dummy;
     end

     for i=raddata(3):raddata(4)
     dummy=fread(fid,2,'float');
     dummy=fread(fid,npI,'float');
     VrI(:,i,rlevel+1)=dummy;
     end

     for i=raddata(3):raddata(4)
     dummy=fread(fid,2,'float');
     dummy=fread(fid,npI,'float');
     VtI(:,i,rlevel+1)=dummy;
     end

     for i=raddata(3):raddata(4)
     dummy=fread(fid,2,'float');
     dummy=fread(fid,npI,'float');
     VpI(:,i,rlevel+1)=dummy;
     end

     for i=raddata(3):raddata(4)
     dummy=fread(fid,2,'float');
     dummy=fread(fid,npI,'float');
     BrI(:,i,rlevel+1)=dummy;
     end

     for i=raddata(3):raddata(4)
     dummy=fread(fid,2,'float');
     dummy=fread(fid,npI,'float');
     BtI(:,i,rlevel+1)=dummy;
     end

     for i=raddata(3):raddata(4)
     dummy=fread(fid,2,'float');
     dummy=fread(fid,npI,'float');
     BpI(:,i,rlevel+1)=dummy;
     end

     end

% Now inner core magnetic field

     for ir=1:nrIc

     dummy=fread(fid,2,'float');
     raddata=fread(fid,4,'float');
     rlevel=raddata(1);
     r(rlevel+1)=raddata(2);

     for i=raddata(3):raddata(4)
     dummy=fread(fid,2,'float');
     dummy=fread(fid,npI,'float');
     BrI(:,i,rlevel+1)=dummy;
     end

     for i=raddata(3):raddata(4)
     dummy=fread(fid,2,'float');
     dummy=fread(fid,npI,'float');
     BtI(:,i,rlevel+1)=dummy;
     end

     for i=raddata(3):raddata(4)
     dummy=fread(fid,2,'float');
     dummy=fread(fid,npI,'float');
     BpI(:,i,rlevel+1)=dummy;
     end

     end

     fclose(fid);
     
% Input done !
     disp('Input finished')
     disp('Number of radial,theta,phi point')
     disp([nr,npI,nt])

%    now for sorting the data out of this strange hemispherical
%    way that the magic code uses

     for ir=1:nrOc
         
     dummy2=TI(:,:,ir);
     for i=1:nt/2
     TI(:,i,ir)=dummy2(:,2*(i-1)+1);
     TI(:,i+nt/2,ir)=dummy2(:,nt-1-2*(i-1)+1);
     end

     dummy2=VrI(:,:,ir);
     for i=1:nt/2
     VrI(:,i,ir)=dummy2(:,2*(i-1)+1);
     VrI(:,i+nt/2,ir)=dummy2(:,nt-1-2*(i-1)+1);
     end

     dummy2=VpI(:,:,ir);
     for i=1:nt/2
     VpI(:,i,ir)=dummy2(:,2*(i-1)+1);
     VpI(:,i+nt/2,ir)=dummy2(:,nt-1-2*(i-1)+1);
     end
     
     dummy2=VtI(:,:,ir);
     for i=1:nt/2
     VtI(:,i,ir)=dummy2(:,2*(i-1)+1);
     VtI(:,i+nt/2,ir)=dummy2(:,nt-1-2*(i-1)+1);
     end
     
     end

     for ir=1:nr

     dummy2=BrI(:,:,ir);
     for i=1:nt/2
     BrI(:,i,ir)=dummy2(:,2*(i-1)+1);
     BrI(:,i+nt/2,ir)=dummy2(:,nt-1-2*(i-1)+1);
     end

     dummy2=BtI(:,:,ir);
     for i=1:nt/2
     BtI(:,i,ir)=dummy2(:,2*(i-1)+1);
     BtI(:,i+nt/2,ir)=dummy2(:,nt-1-2*(i-1)+1);
     end

     dummy2=BpI(:,:,ir);
     for i=1:nt/2
     BpI(:,i,ir)=dummy2(:,2*(i-1)+1);
     BpI(:,i+nt/2,ir)=dummy2(:,nt-1-2*(i-1)+1);
     end

     end

% RENORMALIZE r!

     r=r/(1-radRatio);
     r_cmb=1.0/(1-radRatio);
     r_icb=r_cmb-1;


% Now get the full solution by repeating ms times the structure (in case
% ms is not one)
     np=azsym*npI+1;
     
     Vr=zeros(np,nt,nrOc);
     Vp=zeros(np,nt,nrOc);
     Vt=zeros(np,nt,nrOc);
     Br=zeros(np,nt,nr);
     Bp=zeros(np,nt,nr);
     Bt=zeros(np,nt,nr);
     T =zeros(np,nt,nrOc);

     for i=1:azsym
     Vr((i-1)*npI+1:i*npI,:,:)=VrI;
     Vt((i-1)*npI+1:i*npI,:,:)=VtI;
     Vp((i-1)*npI+1:i*npI,:,:)=VpI;
     Br((i-1)*npI+1:i*npI,:,:)=BrI;
     Bt((i-1)*npI+1:i*npI,:,:)=BtI;
     Bp((i-1)*npI+1:i*npI,:,:)=BpI;
     T((i-1)*npI+1:i*npI,:,:) =TI;
     end
     Vr(np,:,:)=Vr(1,:,:);
     Vt(np,:,:)=Vt(1,:,:);
     Vp(np,:,:)=Vp(1,:,:);
     Br(np,:,:)=Br(1,:,:);
     Bt(np,:,:)=Bt(1,:,:);
     Bp(np,:,:)=Bp(1,:,:);
     T(np,:,:) =T(1,:,:);
    
     phi=(0:np)*2*pi/(np-1);
     rS=zeros(azsym*np+1,nt,nr);
     pS=rS;
     tS=rS;
     for i=1:nr
         rS(:,:,i)=r(i);
     end
     for i=1:np
         pS(i,:,:)=phi(i);
     end
     for i=1:nt
         tS(:,i,:)=theta(i);
     end

     disp('Fields ready for use')
     disp('Field Vr,size,min,max')
     disp([size(Vr),min(min(min(Vr))),max(max(max(Vr)))])
     disp('Field Vt,size,min,max')
     disp([size(Vt),min(min(min(Vt))),max(max(max(Vt)))])
     disp('Field Vp,size,min,max')
     disp([size(Vp),min(min(min(Vp))),max(max(max(Vp)))])
     disp('Field Br,size,min,max')
     disp([size(Br),min(min(min(Br))),max(max(max(Br)))])
     disp('Field Bt,size,min,max')
     disp([size(Bt),min(min(min(Bt))),max(max(max(Bt)))])
     disp('Field Bp,size,min,max')
     disp([size(Bp),min(min(min(Bp))),max(max(max(Bp)))])
     disp('Field T,size,min,max')
     disp([size(T),min(min(min(T))),max(max(max(T)))])


%---------------------------------------------------------------------------

% For temperature isosurfaces

% first define cartesian locations of the points in the spherical mesh
% These are the x,y,z arrays

     phi=[0:azsym*npI]'/(azsym*npI)*2*pi;
     x1=cos(phi)*sint';
     y1=sin(phi)*sint';
     z1=ones(np,1)*cost';

     xS=zeros(np,nt,nr);
     yS=zeros(np,nt,nr);
     zS=zeros(np,nt,nr);
     for i=1:nr
         xS(:,:,i)=r(i)*x1(:,:);
         yS(:,:,i)=r(i)*y1(:,:);
         zS(:,:,i)=r(i)*z1(:,:);
     end

     rMin=r(1);
     rMax=r(nRM);
     pMin=pS(1,1,1);
     pMax=pS(nPM,1,1);
     tMin=theta(1);
     tMax=theta(nTM);
     rAMin=r(1);
     rAMax=r(nRM);
     pAMin=pMin;
     pAMax=pMax;
     tAMin=theta(1);
     tAMax=theta(nTM);

