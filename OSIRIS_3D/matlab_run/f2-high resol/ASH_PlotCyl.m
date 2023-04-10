
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

fl_nm1='/Users/repulsion/Desktop/OSIRIS_3D/Plot_cycl/FLD/b1/b1-000250.h5'; %charge density
fl_nm2='/Users/repulsion/Desktop/OSIRIS_3D/Plot_cycl/FLD/b2/b2-000250.h5'; %charge density
fl_nm3='/Users/repulsion/Desktop/OSIRIS_3D/Plot_cycl/FLD/b3/b3-000250.h5'; %charge density


[x1cart,x2cart,x3cart,b1cart,x1lt,x2lt,x3lt,time]=...
    AshReadHDF5DenDat3D(fl_nm1);
[~,~,~,b2cart,~,~,~,~]=AshReadHDF5DenDat3D(fl_nm2);
[~,~,~,b3cart,~,~,~,~]=AshReadHDF5DenDat3D(fl_nm3);

Ly=x2cart(end);
Lz=x3cart(end);
x2cart=x2cart-Ly/2;
x3cart=x3cart-Lz/2;
[b1cyl,b2cyl,b3cyl]=Cart2Cyl(b3cart,b2cart,b1cart,x3cart,x2cart,x1cart);

mycmap=colormap(jet(2048));

seg = b2cyl;
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



%{
% xele=length(xg); xx=(1:xele)';
% yele=length(yg); yy=(1:yele)';
% zele=length(zg); zz=(1:zele)';
xele=50; xx=(1:xele)';
yele=50; yy=(1:yele)';
zele=50; zz=(1:zele)';
yeze=yele*zele;
tele=xele*yele*zele;

zvec=mod((1:tele)',zele);
zvec(zvec==0)=zele;

yvec=zeros(tele,1);
cnt1=1; cnt2=zele;
for jj=1:yele
    yvec(cnt1:cnt2)=yy(jj);
    cnt1=cnt2+1; cnt2=jj*zele;
end

xvec=zeros(tele,1);
cnt1=1; cnt2=yeze;
for ii=1:xele
    yvec(cnt1:cnt2)=xx(ii);
    cnt1=cnt2+1; cnt2=ii*yeze;
end

dataB3=zeros(tele,1);
szB3=size(dset_b3);
for dd=1:tele
    [ii,jj,kk]=ind2sub(szB3,dd);
    dataB3(dd)=dset_b3(ii,jj,kk);
end


scatter3(xvec,yvec,zvec,[],dataB3);


% pcolor3(dset_b3(1:20,1:20,1:20));

%}