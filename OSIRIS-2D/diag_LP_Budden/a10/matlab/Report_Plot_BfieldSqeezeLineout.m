clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=1:100;

pth0='/Users/repulsion/Desktop/OSIRIS/OSIRIS_2D/atul/counter_finite_0.9c';
pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/Elec2/charge/';
pth3='/MS/DENSITY/Elec3/charge/';
pth4='/MS/FLD/b1/';
pth5='/MS/FLD/b2/';
pth6='/MS/FLD/b3/';


str_h1='charge-Elec1-';
str_h2='charge-Elec2-';
str_h3='charge-Elec3-';
str_h4='b1-';
str_h5='b2-';
str_h6='b3-';


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

BzUL=zeros(length(frm),1);
BzLL=zeros(length(frm),1);


timE=zeros(length(frm),1);
cnt=1;
% Get the width and height of the figure
%============AVI-FILE-NAME==========================
writerObj = VideoWriter ('eden-b-field.avi');
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
    [xg_b3,yg_b1,dset_b1,x1lt_b1,x2lt_b1,time_b1]=...
        ReadECharge13Aug2014(fl_nm4);
    [xg_b2,yg_b2,dset_b2,x1lt_b2,x2lt_b2,time_b2]=...
        ReadECharge13Aug2014(fl_nm5);
    [xg_b3,yg_b3,dset_b3,x1lt_b3,x2lt_b3,time_b3]=...
        ReadECharge13Aug2014(fl_nm6);
    
    dset_qe=dset_qe1+dset_qe2+dset_qe3;
%     dset_J=dset_J1+dset_J2+dset_J3;

    dx_b3=abs(xg_b3(2)-xg_b3(1));
    x1minb3=x1lt_b3(1); x1maxb3=x1lt_b3(2);
    x2minb3=x2lt_b3(1); x2maxb3=x2lt_b3(2);
    [Ngx_LOb3,Ngy_LOb3]=size(dset_b3);
    Ngyby2b3=round(Ngy_LOb3/2);
    LOb3_dat=dset_b3(:,Ngyby2b3);
   B1=zeros(Ngx_LOb3,1);
   B2=zeros(Ngx_LOb3,1);
   B3=zeros(Ngx_LOb3,1);
   B1(:)=(squeeze(dset_b3(:,750)));
   B2(:)=(squeeze(dset_b3(:,500)));
   BzUL(cnt)=max(B1(:));
   BzLL(cnt)=max(B2(:));
   

    timE(cnt)=time_b3;
    cnt=cnt+1;
    


 end
% close(writerObj);
figure(1)
plot(timE,(BzUL),'b','LineWidth',2);hold on;
plot(timE,BzLL,'r','LineWidth',2);

     set(gca,'linewidth',2.0,'fontsize',14,'fontweight','b');
    title('Value of Efield at the boundary');
      xlabel('time(\omega_{pe})');
    
    ylabel('B_{1z}');

    
    
   