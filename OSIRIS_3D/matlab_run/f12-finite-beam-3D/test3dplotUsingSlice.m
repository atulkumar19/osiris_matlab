clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e+22; % cm^-3 58.0e18, 6.12e+18
frm=500;

minch=-4;
maxch=0;

pth0='/home/bhavesh/Desktop/PIC_KH';
pth1='/MS/FLD/b3/';

str_h1='b3-';

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

    fl_nm1=strcat(pth0,pth1,str_h1,str_num,str_ext); %charge density
    
    
    [xg,yg,zg,dset_b3,x1lt,x2lt,x3lt,time]=AshReadHDF5DenDat3D(fl_nm1);
        
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

 
mycmap=colormap(jet(4098));
subplot(2,2,1)
imHandles1 = imagesc(squeeze(xyview(:,:,200)),...
    'XData',[yg yg],'YData',[xg xg],...
    segRange); 
set(gca,'YDir','normal');
set(gca,'FontSize',18,'FontWeight','Bold');
shading('interp');
colormap(mycmap); 
ylabel('x1(c/\omega_{pe})'); xlabel('x2(c/\omega_{pe})');
colorbar; 
title(strcat('e-Density, XY-View, t=',num2str(time)));
axis equal;
axis image; 
axis xy;


subplot(2,2,2);
imHandles2 = imagesc(squeeze(yzview(sn(1),:,:)),...
    'XData',[yg yg],'YData',[zg zg],...
    segRange);  
set(gca,'YDir','normal');
set(gca,'FontSize',18,'FontWeight','Bold');
shading('interp');
colormap(mycmap); 
ylabel('x3(c/\omega_{pe})'); xlabel('x2(c/\omega_{pe})');
colorbar; 
title(strcat('e-Density, YZ-View, t=',num2str(time)));
axis equal;
axis image; 
axis xy;

subplot(2,2,3);
imHandles3 = imagesc(squeeze(xzview(:,sn(2),:)),...
    'XData',[xg xg],'YData',[zg zg],...
    segRange);  
set(gca,'YDir','normal');
set(gca,'FontSize',18,'FontWeight','Bold');
shading('interp');
colormap(mycmap); 
ylabel('x3(c/\omega_{pe})'); xlabel('x1(c/\omega_{pe})');
colorbar; 
title(strcat('e-Density, XZ-View, t=',num2str(time)));
axis equal;
axis image; 
axis xy;

subplot(2,2,4);
hslc=slice(visSeg,sn(1),sn(2),sn(3),'cubic'); 
set(gca,'FontSize',18,'FontWeight','Bold');
shading('interp');
colormap(mycmap); 
ylabel('x2(c/\omega_{pe})'); xlabel('x1(c/\omega_{pe})');zlabel('x3(c/\omega_{pe})');
colorbar; 
title(strcat('e-Density, t=',num2str(time)));
axis equal; 
axis vis3d; 

% set(hslc(1:3),'LineStyle','none');