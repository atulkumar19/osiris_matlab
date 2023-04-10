
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

fl_nm1='b3-000250.h5'; %charge density
[xg,yg,zg,dset_b3,x1lt,x2lt,x3lt,time]=AshReadHDF5DenDat3D(fl_nm1);

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

