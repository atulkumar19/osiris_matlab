
clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e+22; % cm^-3 58.0e18, 6.12e+18
frm=250;

minch=-4;
maxch=0;

pth0='/Users/repulsion/Desktop/OSIRIS_3D/atul_osiris/f2-finite-beam';
pth1='/MS/FLD/b1/';
pth2='/MS/FLD/b2/';
pth3='/MS/FLD/b3/';

str_h1='b1-';
str_h2='b2-';
str_h3='b3-';
str_ext='.h5';
str_num='000000';
len_str_num=length(str_num);

% mycmap=colormap(cmap('Jet',1024));

qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;

len_fr=length(frm);
% scrsz = get(0,'ScreenSize');
% hfig=figure('Position',[50 10 1200 650]);
% set(gcf,'Renderer','OpenGL');
for ii=1:length(frm)    
    
    frno=frm(ii);
    str_frno=num2str(frno);
    len_frno=length(str_frno);
    str_num((len_str_num-len_frno+1):end)=str_frno;

    fl_nm1=strcat(pth0,pth1,str_h1,str_num,str_ext); 
    fl_nm2=strcat(pth0,pth2,str_h2,str_num,str_ext); 
    fl_nm3=strcat(pth0,pth3,str_h3,str_num,str_ext); 
    
    
%     [xg,yg,zg,dset_b1,x1lt,x2lt,x3lt,time]=AshReadHDF5DenDat3D(fl_nm1);
%     [xg,yg,zg,dset_b2,x1lt,x2lt,x3lt,time]=AshReadHDF5DenDat3D(fl_nm2);
%     [xg,yg,zg,dset_b3,x1lt,x2lt,x3lt,time]=AshReadHDF5DenDat3D(fl_nm3);

% fl_nm1='/Users/repulsion/Desktop/OSIRIS_3D/Plot_cycl/FLD/b1/b1-000250.h5'; %charge density
% fl_nm2='/Users/repulsion/Desktop/OSIRIS_3D/Plot_cycl/FLD/b2/b2-000250.h5'; %charge density
% fl_nm3='/Users/repulsion/Desktop/OSIRIS_3D/Plot_cycl/FLD/b3/b3-000250.h5'; %charge density


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


subplot(2,2,1)
imHandles1 = imagesc(squeeze(xyview(:,:,200)),...
    'XData',[x1cart x1cart],'YData',[x2cart+10 x2cart+10],...
    segRange); 
set(gca,'YDir','normal');
set(gca,'FontSize',24,'FontWeight','Bold');
shading('interp');
colormap(mycmap); 
ylabel('x1(c/\omega_{pe})'); xlabel('x2(c/\omega_{pe})');
colorbar; 
title(strcat('B-field, XY-View, t=',num2str(time)));
axis equal;
axis image; 
axis xy;


subplot(2,2,2);
imHandles2 = imagesc(squeeze(yzview(sn(1),:,:)),...
    'XData',[x2cart+10 x2cart+10],'YData',[x3cart+10 x3cart+10],...
    segRange);  
set(gca,'YDir','normal');
set(gca,'FontSize',24,'FontWeight','Bold');
shading('interp');
colormap(mycmap); 
ylabel('x3(c/\omega_{pe})'); xlabel('x2(c/\omega_{pe})');
colorbar; 
title(strcat('B-field, YZ-View, t=',num2str(time)));
axis equal;
axis image; 
axis xy;

subplot(2,2,3);
imHandles3 = imagesc(squeeze(xzview(:,sn(2),:)),...
    'XData',[x1cart x1cart],'YData',[x3cart+10 x3cart+10],...
    segRange);  
set(gca,'YDir','normal');
set(gca,'FontSize',24,'FontWeight','Bold');
shading('interp');
colormap(mycmap); 
ylabel('x3(c/\omega_{pe})'); xlabel('x1(c/\omega_{pe})');
colorbar; 
title(strcat('B-field, XZ-View, t=',num2str(time)));
axis equal;
axis image; 
axis xy;

subplot(2,2,4);
hslc=slice(visSeg,sn(1),sn(2),sn(3),'cubic'); 
set(gca,'FontSize',24,'FontWeight','Bold');
shading('interp');
colormap(mycmap); 
ylabel('x2(c/20\omega_{pe})'); xlabel('x1(c/20\omega_{pe})');zlabel('x3(c/20\omega_{pe})');
colorbar; 
title(strcat('B-field, t=',num2str(time)));
axis equal; 
axis vis3d; 

end

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