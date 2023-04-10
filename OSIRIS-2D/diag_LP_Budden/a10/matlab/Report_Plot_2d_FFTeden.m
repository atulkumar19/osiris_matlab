clc;
clear all;
close all;

addpath('UtilFun');

ne=3.1e20; % cm^-3
frm=1270:1270;

% pth0='/home/testuser/OSIRIS/Budden_plane_test';
% pth0='/home/testuser/OSIRIS/Budden_m1600';
% pth0='/home/testuser/OSIRIS/LaserPlasma';
% pth0='/home/testuser/OSIRIS/Budden_m25_3pulse';
% pth0='/home/testuser/OSIRIS/laser_test_tilt';
pth0='/home/testuser/OSIRIS/oblique_laser_LH';
% pth0='/home/testuser/OSIRIS/MagSonic_Soliton';
% pth0='/home/testuser/OSIRIS/dipole';
pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/ion/charge/';
% pth3='/MS/DENSITY/Elec3/charge/';
% pth4='/MS/DENSITY/Elec4/charge/';
% pth5='/MS/DENSITY/Elec5/charge/';
pth4='/MS/FLD/b3/';

str_h1='charge-Elec1-';
str_h2='charge-ion-';
% str_h3='charge-Elec3-';
% str_h4='charge-Elec4-';
% str_h5='charge-Elec5-';
str_h4='b3-';
str_ext='.h5';
str_num='000000';
len_str_num=length(str_num);

qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;

minch=-5;
maxch=0;
scrsz=get(0,'ScreenSize');
figure('Position',[50 10 1200 650]);
% Get the width and height of the figure
%============AVI-FILE-NAME==========================
writerObj = VideoWriter ('eden-Bfield.avi');
writerObj.FrameRate=5;
open(writerObj);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for ii=1:length(frm)   
    
    frno=frm(ii);
    
    str_frno=num2str(frno);
    len_frno=length(str_frno);
    str_num((len_str_num-len_frno+1):end)=str_frno;
    fl_nm1=strcat(pth0,pth1,str_h1,str_num,str_ext);
    fl_nm2=strcat(pth0,pth2,str_h2,str_num,str_ext);
%     fl_nm3=strcat(pth0,pth3,str_h3,str_num,str_ext);
    fl_nm4=strcat(pth0,pth4,str_h4,str_num,str_ext);
   
    
    [xg_qe,yg_qe,dset_qe1,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm1);
    [xg_qe,yg_qe,dset_qe2,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm2);
%     [xg_qe,yg_qe,dset_qe3,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm3);
   
    [xg_b3,yg_b3,dset_b3,x1lt_b3,x2lt_b3,time_b3]=...
        ReadB313Aug2014(fl_nm4);
    
    
    dset_qe=-dset_qe1+dset_qe2;
    
    x1minb3=x1lt_b3(1); x1maxb3=x1lt_b3(2);
    x2minb3=x2lt_b3(1); x2maxb3=x2lt_b3(2);
    
    dx_qe=abs(xg_qe(2)-xg_qe(1));
    dy_qe=abs(yg_qe(2)-yg_qe(1));
    x1minqe=x1lt_qe(1); x1maxqe=x1lt_qe(2);
    x2minqe=x2lt_qe(1); x2maxqe=x2lt_qe(2);
    [Ngx_LOqe,Ngy_LOqe]=size(dset_qe);
    Ngyby2qe=round(Ngy_LOqe/2);
    LOqe_dat=dset_qe(:,Ngyby2qe);
    %-------fft----------------%
%     nfft=1024;
%    
    subplot(121)
%{    E=zeros(Ngx_LOb3,1);
%    E(:)=squeeze(dset_b3(45,:));
%   F3=fftshift(fft( E ,nfft));
% %    kx = (-nfft/2:nfft/2-1)/(dx_b3*nfft);%wave Vector
%   Pows=F3.*conj(F3)/(nfft*nfft); %computing power with proper scaling
%   %gama1=smoothn(abs(F3)/(Ngx_LOb3*dx_b3),'0.1');
% 
% %      E(:,1)=mean(gama,2);%mean along y axix
% %      plot((kx),abs(F3)/(Ngx_LOb3*dx_b3),'k','LineWidth',2);hold on;
%    %----1sided------%
%     kx=(0:nfft/2-1)/(dx_b3*nfft);	 	 
%  loglog((kx),abs( F3(nfft/2+1:nfft))/(Ngx_LOb3*dx_b3),'k','LineWidth',2);
%  print('autoExample', '-depsc2', '-r300');
%   set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
%} 

  
  
  
  nfft=Ngx_LOqe;nfft1=Ngy_LOqe;E=zeros(Ngx_LOqe,Ngy_LOqe);
  F3=fftshift(fft2( dset_qe ,nfft,nfft1));
   Pows=F3.*conj(F3)/(nfft*nfft1);
   
 kx=zeros(Ngx_LOqe/2,1);
 ky=zeros(Ngy_LOqe/2,1);
 kx(:,1) =(0:nfft/2-1)*1/(dx_qe*nfft);
 ky(:,1) = (0:nfft1/2-1)*1/(dy_qe*nfft1);
 [kx,ky]=meshgrid(kx,ky);
 E=abs(Pows(nfft/2+1:nfft,nfft1/2+1:nfft1))/(Ngx_LOqe*dx_qe*Ngy_LOqe*dy_qe);
 Pxx=abs(F3).^2;
 Pxx=Pxx(nfft/2+1:nfft,nfft1/2+1:nfft1);
 pcolor(kx,ky,(E'));colorbar;shading('interp'); 
 set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
%  xlim([0 0.04]);ylim([0 0.5]);
  
  
  
  
  
  
    %-----------end fft------------------%
subplot(122)
    
    imagesc(dset_qe','XData',[x1minb3 x1maxb3],...
        'YData',[x2minb3 x2maxb3]);
    set(gca,'YDir','normal');
    colormap('jet'); colorbar('south')
    shading('interp');
    set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
    print('autoExample', '-depsc2', '-r300');
    xlabel('x1(c/\omega_{pe})');
    
    ylabel('x2(c/\omega_{pe})');
   
    xlim([0 25]); ylim([0 25]);
    
%     title(num2str(ii))
     
    drawnow;
    
%     saveas(gcf,'CurrentFig.png');
%     img4=imread('CurrentFig.png');
%     writeVideo(writerObj,img4);
% gif_add_frame(gcf,'finite_beam_spectra.gif',2);

end

% close(writerObj);
