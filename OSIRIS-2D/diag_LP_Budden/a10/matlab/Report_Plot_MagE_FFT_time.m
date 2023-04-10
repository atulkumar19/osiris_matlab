clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=250:250;

% pth0='/home/testuser/OSIRIS/Budden_plane_test';
% pth0='/home/testuser/OSIRIS/Budden_m1600';
% pth0='/home/testuser/OSIRIS/LaserPlasma';
% pth0='/home/testuser/OSIRIS/Budden_m25_3pulse';
% pth0='/home/testuser/OSIRIS/laser_test_tilt';
% pth0='/home/testuser/OSIRIS/oblique_laser';
pth0='/home/testuser/OSIRIS/oblique_laser_gaussian';
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
    
   
    dset_qe=dset_qe1+dset_qe2;
    
    
    dx_b3=abs(xg_b3(2)-xg_b3(1));
    x1minb3=x1lt_b3(1); x1maxb3=x1lt_b3(2);
    x2minb3=x2lt_b3(1); x2maxb3=x2lt_b3(2);
    [Ngx_LOb3,Ngy_LOb3]=size(dset_b3);
    Ngyby2b3=round(Ngy_LOb3/2);
    LOb3_dat=dset_b3(:,Ngyby2b3);
    
   B3(ii)=sum(squeeze(dset_b3(:,2500)));
end

  nfft=10240;
  dt=0.0707*30;
  Fs=1/dt;
  t=(1/Fs)*(1:nfft);
%   B3=detrend(B3);
B3=B3-mean(B3);
  f_fft=abs(fftshift(fft(B3(0:50:length(B3)),nfft)));
  y_fft2=f_fft(nfft/2+1:nfft);
  f2=Fs*(0:nfft/2-1)/nfft;
  plot(2*pi*f2,y_fft2,'-','linewidth',2.0,'color',[0.7 0 0]);
 xlabel('Frequency ');
ylabel('Amplitude');
  

    

%      set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
%     title('Growth of Magnetic Enery');
%       xlabel('time(\omega_{pe})');
%     
%     ylabel('log|B_z|^2');
% %  aviobj = close ( aviobj );
% % close(writerObj);