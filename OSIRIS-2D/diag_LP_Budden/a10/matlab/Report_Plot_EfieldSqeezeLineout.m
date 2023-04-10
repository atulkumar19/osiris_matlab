clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=0:500;

pth0='/Users/repulsion/Desktop/OSIRIS/OSIRIS_2D/atul/counter_finite_0.9c';
pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/Elec2/charge/';
pth3='/MS/DENSITY/Elec3/charge/';
pth4='/MS/FLD/e1/';
pth5='/MS/FLD/e2/';
pth6='/MS/FLD/e3/';


str_h1='charge-Elec1-';
str_h2='charge-Elec2-';
str_h3='charge-Elec3-';
str_h4='e1-';
str_h5='e2-';
str_h6='e3-';


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

ExUL=zeros(length(frm),1);
ExLL=zeros(length(frm),1);

EyUL=zeros(length(frm),1);
EyLL=zeros(length(frm),1);
timE=zeros(length(frm),1);
cnt=1;
%============AVI-FILE-NAME==========================
writerObj = VideoWriter ('eden-b-field.avi');
writerObj.FrameRate=5;
open(writerObj);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for ii=1:10:length(frm)   
    
    frno=frm(ii);
    
    str_frno=num2str(frno);
    len_frno=length(str_frno);
    str_num((len_str_num-len_frno+1):end)=str_frno;
    fl_nm1=strcat(pth0,pth1,str_h1,str_num,str_ext);
    fl_nm2=strcat(pth0,pth2,str_h2,str_num,str_ext);
    fl_nm3=strcat(pth0,pth3,str_h3,str_num,str_ext);
    fl_nm4=strcat(pth0,pth4,str_h4,str_num,str_ext);
    fl_nm5=strcat(pth0,pth5,str_h5,str_num,str_ext);
    fl_nm6=strcat(pth0,pth6,str_h6,str_num,str_ext);
    
    
    
    [xg_qe,yg_qe,dset_qe1,x1lt_qe,x2lt_qe,time_qe]=...
        ReadECharge13Aug2014(fl_nm1);
    [xg_qe,yg_qe,dset_qe2,x1lt_qe,x2lt_qe,time_qe]=...
        ReadECharge13Aug2014(fl_nm2);
    [xg_qe,yg_qe,dset_qe3,x1lt_qe,x2lt_qe,time_qe]=...
        ReadECharge13Aug2014(fl_nm3);
    [xg_e1,yg_e1,dset_e1,x1lt_e1,x2lt_e1,time_e1]=...
        ReadECharge13Aug2014(fl_nm4);
    [xg_e2,yg_e2,dset_e2,x1lt_e2,x2lt_e2,time_e2]=...
        ReadECharge13Aug2014(fl_nm5);
    [xg_e3,yg_e3,dset_e3,x1lt_e3,x2lt_e3,time_e3]=...
        ReadECharge13Aug2014(fl_nm6);
    
    dset_qe=dset_qe1+dset_qe2+dset_qe3;
%     dset_J=dset_J1+dset_J2+dset_J3;

 dx_e3=abs(xg_e3(2)-xg_e3(1));
    x1mine3=x1lt_e3(1); x1maxe3=x1lt_e3(2);
    x2mine3=x2lt_e3(1); x2maxe3=x2lt_e3(2);
    [Ngx_LOe3,Ngy_LOe3]=size(dset_e3);
    Ngyby2e3=round(Ngy_LOe3/2);
    LOe3_dat=dset_e3(:,Ngyby2e3);
   E1=zeros(Ngx_LOe3,1);
   E2=zeros(Ngx_LOe3,1);
   E3=zeros(Ngx_LOe3,1);
   E4=zeros(Ngx_LOe3,1);
   E1(:)=(squeeze(dset_e1(:,750)));
   E2(:)=(squeeze(dset_e1(:,500)));
   E3(:)=(squeeze(dset_e2(:,750)));
   E4(:)=(squeeze(dset_e2(:,500)));
   
   ExUL(cnt)=max(E1(:));
   ExLL(cnt)=max(E2(:));
   
   EyUL(cnt)=max(E3(:));
   EyLL(cnt)=max(E4(:));
    

    timE(cnt)=time_e3;
    cnt=cnt+1;
    
%     subplot(2,2,1)
%     plot(yg_e1,E1(:),'b','LineWidth',2);
%     set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
%      title(strcat(' t= ',num2str(time_qe)));
%       xlabel('x2(c/\omega_{pe})');
%     
%     ylabel('E_x(mc\omega_{pe}/e)');
%     subplot(2,2,2)
%     plot(yg_e2,E2(:),'b','LineWidth',2);
%      set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
%     title(strcat('t= ',num2str(time_qe)));
%       xlabel('x2(c/\omega_{pe})');
%     
%     ylabel('E_y(mc\omega_{pe}/e)');
%     subplot(2,2,[3 4])
%     plot(yg_e3,E3(:),'b','LineWidth',2);
%      set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
%       title(strcat('t= ',num2str(time_qe)));
%       xlabel('x2(c/\omega_{pe})');
%     
%     ylabel('E_z(mc\omega_{pe}/e)');
%     drawnow
%     saveas(gcf,'CurrentFig.png');
%     img4=imread('CurrentFig.png');
%     writeVideo(writerObj,img4);
%      L = graphicsversion('handlegraphics'); 
%      currFrame = getframe;


 end
% close(writerObj);
figure(1)
plot(timE,ExUL,'b','LineWidth',2);hold on;
plot(timE,ExLL,'r','LineWidth',2);

     set(gca,'linewidth',2.0,'fontsize',14,'fontweight','b');
    title('Value of Efield at the boundary');
      xlabel('time(\omega_{pe})');
    
    ylabel('E_{1x}');
figure(2)
plot(timE,EyUL,'b','LineWidth',2);hold on;
plot(timE,EyLL,'r','LineWidth',2);

     set(gca,'linewidth',2.0,'fontsize',14,'fontweight','b');
    title('Value of Efield at the boundary');
      xlabel('time(\omega_{pe})');
    
    ylabel('E_{1y}');