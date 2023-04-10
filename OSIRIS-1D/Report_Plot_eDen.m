clc;
clear al
close all;

addpath('UtilFun');

ne=1.0e22; % cm^-3
frm=3:3;

pth0='/Users/repulsion/Desktop/OSIRIS/OSIRIS4.0/OSIRIS_RUN/OSIRIS-1D/MS1D-7.0';
pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/ion/charge/';
%  pth3='/MS/DENSITY/Elec3/charge/';
%   

pth4='/MS/FLD/b3/';
pth5='/MS/FLD/e1/';
pth6='/MS/FLD/e2/';


pth7='/MS/DENSITY/Elec1/j1/';
% pth8='/MS/DENSITY/Elec2/j2/';
%  pth9='/MS/DENSITY/Elec3/j2/';
  
str_h1='charge-Elec1-';
str_h2='charge-ion-';
%  str_h3='charge-Elec3-';


str_h4='b3-';
str_h5='e1-';
str_h6='e2-';

str_h7='j1-Elec1-';
% % str_h8='j2-Elec2-';
% % str_h9='j2-Elec3-';
% 

str_ext='.h5';
str_num='000000';
len_str_num=length(str_num);

qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;

% minch=-5;
% maxch=0;
% scrsz=get(0,'ScreenSize');
% figure('Position',[50 10 1200 650]);
% Get the width and height of the figure
%============AVI-FILE-NAME==========================
writerObj = VideoWriter ('f2-1D-vth-0.0004.avi');
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
%    fl_nm3=strcat(pth0,pth3,str_h3,str_num,str_ext);
   
    fl_nm4=strcat(pth0,pth4,str_h4,str_num,str_ext);
    
    fl_nm5=strcat(pth0,pth5,str_h5,str_num,str_ext);
    fl_nm6=strcat(pth0,pth6,str_h6,str_num,str_ext);
    
    fl_nm7=strcat(pth0,pth7,str_h7,str_num,str_ext);
%     
%     fl_nm8=strcat(pth0,pth8,str_h8,str_num,str_ext);
%   fl_nm9=strcat(pth0,pth9,str_h9,str_num,str_ext);
%     
    
    [xg_qe,dset_qe1,x1lt_qe,time_qe1,]=...
        ReadCharge10June2016(fl_nm1);
    [xg_qe,dset_qe2,x1lt_qe,time_qe2,]=...
        ReadCharge10June2016(fl_nm2);
%     [xg_qe,dset_qe3,x1lt_qe,time_qe3,]=...
%         ReadCharge10June2016(fl_nm3);
%   
   
    [xg_b3,dset_b3,x1lt_b3,time_b3]=...
       ReadB310June2016(fl_nm4);
    [xg_b3,dset_e1,x1lt_b3,time_b3]=...
       ReadB310June2016(fl_nm5);
    
[xg_b3,dset_e2,x1lt_b3,time_b3]=...
        ReadB310June2016(fl_nm6);
    [xg_qe,dset_j1,x1lt_qe,time_qe1,]=...
        ReadCharge10June2016(fl_nm7);
%     [xg_qe,dset_j2,x1lt_qe,time_qe2,]=...
%         ReadCharge10June2016(fl_nm8);
%     [xg_qe,dset_j3,x1lt_qe,time_qe3,]=...
%         ReadCharge10June2016(fl_nm9);
  
%     
%     
%     
    
    dset_qe=(dset_qe1+dset_qe2);
     dset_jy=(dset_j1);
    
    x1minb3=x1lt_b3(1); x1maxb3=x1lt_b3(2);

    dx_qe=abs(xg_qe(2)-xg_qe(1));
    x1minqe=x1lt_qe(1); x1maxqe=x1lt_qe(2);

    [Ngx_LOqe]=length(dset_qe);
    Ngxby2qe=round(Ngx_LOqe/2);
    LOqe_dat=dset_qe(Ngxby2qe);
     
   
    
  
    ax1=subplot(221);
    plot(xg_qe,(dset_qe2),'g');
    set(gca,'FontSize',16,'FontWeight','Bold');
    xlabel('x2(c/\omega_{pe})'); ylabel('J_y');
    title(strcat('e-Density, time:',num2str(time_qe1)));
%     ylim([0.9989 1.001]);
%     xlim([0 25]);
    
    
    
 ax2=subplot(222);
    
   plot(xg_qe,(dset_b3),'k');
    set(gca,'FontSize',16,'FontWeight','Bold');
    xlabel('x2(c/\omega_{pe})'); ylabel('B_z');
    title(strcat('B-Field, time:',num2str(time_b3)));
    
%    xlim([0 25]);
    

     ax3=subplot(223);
    
   plot(xg_qe,(dset_e1),'r');
    set(gca,'FontSize',16,'FontWeight','Bold');
    xlabel('x2(c/\omega_{pe})'); ylabel('E_{1x}');
    title(strcat('E_{x}-Field, time:',num2str(time_b3)));
%     xlim([0 25]);
      ax4=subplot(224);
    
   plot(xg_qe,(dset_e2),'g');
    set(gca,'FontSize',16,'FontWeight','Bold');
    xlabel('x2(c/\omega_{pe})'); ylabel('E_{1y}');
    title(strcat('E_y-Field, time:',num2str(time_b3)));
% xlim([0 25]);
    drawnow
    saveas(gcf,'CurrentFig.png');
    img4=imread('CurrentFig.png');
    writeVideo(writerObj,img4);
     L = graphicsversion('handlegraphics'); 
     currFrame = getframe;

 end
  close(writerObj);