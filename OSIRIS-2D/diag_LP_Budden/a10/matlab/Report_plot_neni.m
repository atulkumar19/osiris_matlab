clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=0:1100;

pth0='/home/testuser/OSIRIS/Budden_plane_test';
pth1='/MS_m25/DENSITY/Elec1/charge/';
pth2='/MS_m25/DENSITY/ion/charge/';
% pth3='/MS/DENSITY/Elec3/charge/';
pth4='/MS_m25/DENSITY/Elec1/j1/';
pth5='/MS_m25/DENSITY/ion/j1/';
% pth6='/MS/DENSITY/Elec3/j1/';
pth7='/MS_m25/DENSITY/Elec1/j2/';
pth8='/MS_m25/DENSITY/ion/j2/';
% pth9='/MS/DENSITY/Elec3/j2/';
pth10='/MS_m25/FLD/e1/';
pth11='/MS_m25/FLD/e2/';
pth12='/MS_m25/FLD/e3/';

pth13='/MS_m25/FLD/b3/';

str_h1='charge-Elec1-';
str_h2='charge-ion-';
% str_h3='charge-Elec3-';
str_h4='j1-Elec1-';
str_h5='j1-ion-';
% str_h6='j1-Elec3-';
str_h7='j2-Elec1-';
str_h8='j2-ion-';
% str_h9='j2-Elec3-';
str_h10='e1-';
str_h11='e2-';
str_h12='e3-';

str_h13='b3-';

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
writerObj = VideoWriter ('neni_m9.avi');
writerObj.FrameRate=5;
open(writerObj);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for ii=1:5:length(frm)   
    
    frno=frm(ii);
    
    str_frno=num2str(frno);
    len_frno=length(str_frno);
    str_num((len_str_num-len_frno+1):end)=str_frno;
    fl_nm1=strcat(pth0,pth1,str_h1,str_num,str_ext);
    fl_nm2=strcat(pth0,pth2,str_h2,str_num,str_ext);
%     fl_nm3=strcat(pth0,pth3,str_h3,str_num,str_ext);
    fl_nm4=strcat(pth0,pth4,str_h4,str_num,str_ext);
    fl_nm5=strcat(pth0,pth5,str_h5,str_num,str_ext);
%     fl_nm6=strcat(pth0,pth6,str_h6,str_num,str_ext);
    fl_nm7=strcat(pth0,pth7,str_h7,str_num,str_ext);
    fl_nm8=strcat(pth0,pth8,str_h8,str_num,str_ext);
%     fl_nm9=strcat(pth0,pth9,str_h9,str_num,str_ext);
    fl_nm10=strcat(pth0,pth10,str_h10,str_num,str_ext);
    fl_nm11=strcat(pth0,pth11,str_h11,str_num,str_ext);
    fl_nm12=strcat(pth0,pth12,str_h12,str_num,str_ext);
    fl_nm13=strcat(pth0,pth13,str_h13,str_num,str_ext);
    
    
    
    
    [xg_qe,yg_qe,dset_qe1,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm1);
    [xg_qe,yg_qe,dset_qe2,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm2);
%     [xg_qe,yg_qe,dset_qe3,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm3);
    [xg_qe,yg_qe,dset_J1,x1lt_J1,x2lt_J1,time_J1]=...
        ReadECharge13Aug2014(fl_nm4);
    [xg_qe,yg_qe,dset_J2,x1lt_J2,x2lt_J2,time_J2]=...
        ReadECharge13Aug2014(fl_nm5);
%     [xg_qe,yg_qe,dset_J3,x1lt_J3,x2lt_J3,time_J3]=...
%         ReadECharge13Aug2014(fl_nm6);
    [xg_qe,yg_qe,dset_J4,x1lt_J4,x2lt_J4,time_J4]=...
        ReadECharge13Aug2014(fl_nm7);
    [xg_qe,yg_qe,dset_J5,x1lt_J5,x2lt_J5,time_J5]=...
        ReadECharge13Aug2014(fl_nm8);
%     [xg_qe,yg_qe,dset_J6,x1lt_J6,x2lt_J6,time_J6]=...
%         ReadECharge13Aug2014(fl_nm9);
%     
    [xg_e1,yg_e1,dset_e1,x1lt_e1,x2lt_e1,time_e1]=...
        ReadE113Aug2014(fl_nm10);
    [xg_e2,yg_e2,dset_e2,x1lt_e2,x2lt_e2,time_e2]=...
        ReadE213Aug2014(fl_nm11);
    [xg_e3,yg_e3,dset_e3,x1lt_e3,x2lt_e3,time_e3]=...
        ReadE213Aug2014(fl_nm12);
    
    [xg_b3,yg_b3,dset_b3,x1lt_b3,x2lt_b3,time_b3]=...
        ReadB313Aug2014(fl_nm13); 
    
    dset_qe=dset_qe1+dset_qe2;
    dset_J1=dset_J1+dset_J2;
    dset_J2=dset_J4+dset_J5;
     mean1_qe1=mean(dset_qe1,2);
      mean1_qe2=mean(dset_qe2,2);
    mean1_qe=mean(dset_qe,2);
    mean1_J1=mean(dset_J1,2);
    mean1_J2=mean(dset_J2,2);
    mean1_E1=mean(dset_e1,2);
    mean1_E2=mean(dset_e2,2);
    mean1_b3=mean(dset_b3,2);
    
    
%     subplot(2,3,1)
%     plot(xg_qe,mean1_J1,'b','LineWidth',2);
%     set(gca,'linewidth',2.0,'fontsize',14,'fontweight','b');
%     title(strcat(' t= ',num2str(time_qe)));
%     xlim([0 500]);
%     xlabel('x1(c/\omega_{pe})');
%     ylabel('J_x(n_0ec)');
%     subplot(2,3,2)
%     plot(xg_qe,mean1_J2,'b','LineWidth',2);
%      set(gca,'linewidth',2.0,'fontsize',14,'fontweight','b');
%     title(strcat(' t= ',num2str(time_qe)));
%     xlim([0 500]);
%       xlabel('x1(c/\omega_{pe})');
%     
%    ylabel('J_y(n_0ec)');
%     subplot(2,3,3)
%     plot(xg_b3,mean1_b3,'b','LineWidth',2);
%      set(gca,'linewidth',2.0,'fontsize',14,'fontweight','b');
%      title(strcat(' t= ',num2str(time_qe)));
%      xlim([0 500]);
%       xlabel('x1(c/\omega_{pe})');
%     
%    ylabel('B_z(mc\omega_{pe}/e)');
%    
%    subplot(2,3,4)
%     plot(xg_qe,mean1_E1,'b','LineWidth',2);
%      set(gca,'linewidth',2.0,'fontsize',14,'fontweight','b');
%       title(strcat(' t= ',num2str(time_qe)));
%       xlim([0 500]);
%       xlabel('x1(c/\omega_{pe})');
%     
%     ylabel('E_x(m\omega_{pe}/e)');
%     subplot(2,3,5)
%     plot(xg_e2,mean1_E2,'b','LineWidth',2);
%      set(gca,'linewidth',2.0,'fontsize',14,'fontweight','b');
%     title(strcat(' t= ',num2str(time_qe)));
%     xlim([0 500]);
%       xlabel('x1(c/\omega_{pe})');
%     
%     ylabel('E_y(m\omega_{pe}/e)');
    subplot(1,2,1)
    plot(xg_e3,mean1_qe1,'b','LineWidth',2);
     set(gca,'linewidth',2.0,'fontsize',14,'fontweight','b');
      title(strcat(' t= ',num2str(time_qe)));
      xlim([0 500]);
      xlabel('x1(c/\omega_{pe})')
       ylabel('e-den(n_0)');
      subplot(1,2,2)
    plot(xg_e3,mean1_qe2,'b','LineWidth',2);
     set(gca,'linewidth',2.0,'fontsize',14,'fontweight','b');
      title(strcat(' t= ',num2str(time_qe)));
      xlim([0 500]);
      
      xlabel('x1(c/\omega_{pe})')
       ylabel('i-den(n_0)');
     drawnow
    saveas(gcf,'CurrentFig.png');
    img4=imread('CurrentFig.png');
    writeVideo(writerObj,img4);
     L = graphicsversion('handlegraphics'); 
     currFrame = getframe;


 end
close(writerObj)