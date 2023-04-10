clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e+22; % cm^-3 58.0e18, 6.12e+18

minch=-4;
maxch=0;

% mycmap=colormap(cmap('Jet',1024));

qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;

% len_fr=length(frm);
% scrsz = get(0,'ScreenSize');
% hfig=figure('Position',[50 10 1200 650]);
% set(gcf,'Renderer','OpenGL');

% mycmap=colormap(cmap('blue',1024));
mycmap=colormap(jet(2048));


fl_nm1='b3-000150.h5'; %charge density
[xg,yg,zg,dsetb3,x1lt,x2lt,x3lt,time]=AshReadHDF5DenDat3D(fl_nm1);


%{
xlen=100; ylen=100; zlen=100;
tlen=xlen*ylen*zlen;
b3vec=zeros(tlen,1);
xvec=zeros(tlen,1);
yvec=zeros(tlen,1);
zvec=zeros(tlen,1);


cnt=1;
for ii=500:(xlen+500)
    for jj=500:(ylen+500)
        for kk=500:(zlen+500)
            xvec(cnt)=xg(ii);
            yvec(cnt)=yg(jj);
            zvec(cnt)=zg(kk);
            b3vec(cnt)=dsetb3(ii,jj,kk);
            cnt=cnt+1;
        end
    end
end
fscatter3(xvec,yvec,zvec,b3vec);
%}

seg = dsetb3;
segRange = [min(seg(:)) max(seg(:))];
sizeSeg = size(seg); %will never need to change.
sn = [300 200 0];
visSeg = seg;
visSeg = permute(visSeg, [3 1 2]);

xyview = seg;
xzview = seg;
yzview = seg;

xyview = flipdim(xyview, 1);
xzview = flipdim(xzview, 1);
xzview = permute(xzview, [3 2 1]);
%xzview = flipdim(xzview, 2);
yzview = permute(yzview, [1 3 2]);
yzview = flipdim(yzview, 1);

imHandles1 = imagesc(squeeze(xyview(:,:,200)), segRange); 
colormap(mycmap); 
axis image; 
axis xy;

figure;
imHandles2 = imagesc(squeeze(yzview(sn(1),:,:)), segRange); 
colormap(mycmap); 
axis image; 
axis xy;

figure;
imHandles3 = imagesc(squeeze(xzview(:,sn(2),:)), segRange); 
colormap(mycmap); 
axis image; 
axis xy;

figure;
hslc=slice(visSeg,sn(1),sn(2),sn(3)); shading('interp');
colormap(mycmap);
colorbar;
axis equal; 
axis vis3d; 
% set(hslc(1:3),'LineStyle','none');