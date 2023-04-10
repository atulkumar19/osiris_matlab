% clc;
% clear all;
% close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=595:595;

% pth0='/home/atul/no-RR-0.9999999-ion';
pth0='/home/atul/OSIRIS/Sym_Ion_Recon';
pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/Elec2/charge/';
%  pth3='/MS/DENSITY/Elec3/charge/';


pth5='/MS/FLD/b3/';

str_h1='charge-Elec1-';
str_h2='charge-Elec2-';
%  str_h3='charge-Elec3-';





str_h5='b3-';

str_ext='.h5';
str_num='000000';
len_str_num=length(str_num);

qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;

minch=-5;
maxch=0;
% scrsz=get(0,'ScreenSize');
% figure('Position',[50 10 1200 650]);
% Get the width and height of the figure
MagE=zeros(length(frm),1);
timE=zeros(length(frm),1);
cnt=1;
% %============AVI-FILE-NAME==========================
%  aviobj = avifile ( 'squeeze-eden-b-cur', 'fps',5); 
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for ii=1:length(frm)   
    
    frno=frm(ii);
    
    str_frno=num2str(frno);
    len_frno=length(str_frno);
    str_num((len_str_num-len_frno+1):end)=str_frno;
    fl_nm1=strcat(pth0,pth1,str_h1,str_num,str_ext);
    fl_nm2=strcat(pth0,pth2,str_h2,str_num,str_ext);
%     fl_nm3=strcat(pth0,pth3,str_h3,str_num,str_ext);
    
    fl_nm5=strcat(pth0,pth5,str_h5,str_num,str_ext);
    
    
    [xg_qe,yg_qe,dset_qe1,x1lt_qe,x2lt_qe,time_qe]=...
        ReadECharge13Aug2014(fl_nm1);
    [xg_qe,yg_qe,dset_qe2,x1lt_qe,x2lt_qe,time_qe]=...
        ReadECharge13Aug2014(fl_nm2);
   
%     [xg_qe,yg_qe,dset_qe3,x1lt_qe,x2lt_qe,time_qe]=...
%         ReadECharge13Aug2014(fl_nm3);
% %    
%     
    [xg_b3,yg_b3,dset_b3,x1lt_b3,x2lt_b3,time_b3]=...
        ReadB313Aug2014(fl_nm5);    
%     dset_qe=dset_qe1+dset_qe2+dset_qe3;
dset_qe=dset_qe1+dset_qe2;
    
    
    dx_b3=abs(xg_b3(2)-xg_b3(1));
    x1minb3=x1lt_b3(1); x1maxb3=x1lt_b3(2);
    x2minb3=x2lt_b3(1); x2maxb3=x2lt_b3(2);
    [Ngx_LOb3,Ngy_LOb3]=size(dset_b3);
    Ngyby2b3=round(Ngy_LOb3/2);
    LOb3_dat=dset_b3(:,Ngyby2b3);
    
    MagE(cnt)=sum(sum(dset_b3))*dx_b3.^2;
    timE(cnt)=time_b3;
    cnt=cnt+1;
    
    
   
end

semilogy(timE,MagE,'b','LineWidth',2);hold on;
     set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
    title('Growth of Magnetic Enery');
      xlabel('time(\omega_{pe})');
      
    
    ylabel('log|B_z|^2');
%     xlim([0 140]);
%  aviobj = close ( aviobj );
% close(writerObj);