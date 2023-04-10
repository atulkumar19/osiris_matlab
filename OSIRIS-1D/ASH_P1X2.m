clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=00:500;

minch=-4;
maxch=0;

pth0='/Users/repulsion/Desktop/OSIRIS/OSIRIS4.0/OSIRIS_RUN/OSIRIS-1D/f2-1d';
pth1='/MS/DENSITY/Elec1/charge/';
% pth2='/MS/DENSITY/ion/charge/';
% pth3='/MS/DENSITY/Elec3/charge/';
pth4='/MS/RAW/Elec1/';

str_h1='charge-Elec1-';
% % str_h2='charge-ion-';
% str_h3='charge-Elec3-';
str_h4='RAW-Elec1-';

str_ext='.h5';
str_num='000000';
len_str_num=length(str_num);


qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;

len_fr=length(frm);
scrsz = get(0,'ScreenSize');
hfig=figure('Position',[50 10 1200 650]);
% Get the width and height of the figure
TotE=zeros(length(frm),1);

cnt=1;
% %============AVI-FILE-NAME==========================
% writerObj = VideoWriter ('phase_space.avi');
% writerObj.FrameRate=5;
% open(writerObj);
% % aviobj = avifile ( 'eden-b-field', 'fps',5); 
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for ii=1:length(frm)    
    
    frno=frm(ii);
    str_frno=num2str(frno);
    len_frno=length(str_frno);
    str_num((len_str_num-len_frno+1):end)=str_frno;
    fl_nm1=strcat(pth0,pth1,str_h1,str_num,str_ext);
%     fl_nm2=strcat(pth0,pth2,str_h2,str_num,str_ext);
%     fl_nm3=strcat(pth0,pth3,str_h3,str_num,str_ext);

    fl_nm4=strcat(pth0,pth4,str_h4,str_num,str_ext); %raw data
    
   [xg_qe,dset_qe1,x1lt_qe,time_qe1]=...
        ReadCharge10June2016(fl_nm1);
%     [xg_qe,yg_qe,dset_qe2,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm2);
%     [xg_qe,yg_qe,dset_qe3,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm3);

    [ene_raw,x1_raw,p1_raw,tag_raw,xlt_in,xlt_fi]=...
       ReadRaw4May2017_1D(fl_nm4);
   
    dx_qe=abs(xg_qe(2)-xg_qe(1));
    
    
    TotE(cnt)=sum(ene_raw)*dx_qe;
    timE(cnt)=time_qe1;
    cnt=cnt+1;
    
end
    semilogy(timE,TotE,'m','LineWidth',2);
set(gca,'fontsize',20,'fontweight','b');
xlabel('time(\omega_{pe})');
  
 
% aviobj = close ( aviobj );
%   close(writerObj);