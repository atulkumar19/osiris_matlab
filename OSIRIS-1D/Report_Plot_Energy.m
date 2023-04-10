% clc;
% clear all;
% close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=00:500;

pth0='/Users/repulsion/Desktop/OSIRIS/OSIRIS4.0/OSIRIS_RUN/OSIRIS-1D/23.7.18/f2-sin-1-low';
pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/Elec2/charge/';
% pth3='/MS/DENSITY/Elec3/charge/';
pth4='/MS/DENSITY/Elec1/ene/';
pth5='/MS/DENSITY/Elec2/ene/';
% pth6='/MS/DENSITY/Elec3/ene/';
pth7='/MS/FLD/ene_emf/';
pth8='/MS/FLD/ene_e/';
pth9='/MS/FLD/ene_b/';

pth10='/MS/FLD/b3/';
pth11='/MS/FLD/e1/';
pth12='/MS/FLD/e2/';





str_h1='charge-Elec1-';
str_h2='charge-Elec2-';
% str_h3='charge-Elec3-';

str_h4='ene-Elec1-';
str_h5='ene-Elec2-';
% str_h6='ene-Elec3-';

str_h7='ene_emf-';
str_h8='ene_e-';
str_h9='ene_b-';

str_h10='b3-';
str_h11='e1-';
str_h12='e2-';



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
% % Get the width and height of the figure
MagE=zeros(length(frm),1);
ElecE=zeros(length(frm),1);
EME=zeros(length(frm),1);
KinE=zeros(length(frm),1);
timE=zeros(length(frm),1);
TotE=zeros(length(frm),1);
KinE0=zeros(length(frm),1);
KinE_per=zeros(length(frm),1);
E1E=zeros(length(frm),1);
E2E=zeros(length(frm),1);
B3E=zeros(length(frm),1);
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
    fl_nm3=strcat(pth0,pth3,str_h3,str_num,str_ext);
    fl_nm4=strcat(pth0,pth4,str_h4,str_num,str_ext);
    fl_nm5=strcat(pth0,pth5,str_h5,str_num,str_ext);
    fl_nm6=strcat(pth0,pth6,str_h6,str_num,str_ext);
    fl_nm7=strcat(pth0,pth7,str_h7,str_num,str_ext);
    fl_nm8=strcat(pth0,pth8,str_h8,str_num,str_ext);
    fl_nm9=strcat(pth0,pth9,str_h9,str_num,str_ext);
      fl_nm10=strcat(pth0,pth10,str_h10,str_num,str_ext);
    
    fl_nm11=strcat(pth0,pth11,str_h11,str_num,str_ext);
     fl_nm12=strcat(pth0,pth12,str_h12,str_num,str_ext);
    
    
    
   [xg_qe,dset_qe1,x1lt_qe,time_qe1,]=...
        ReadCharge10June2016(fl_nm1);
    [xg_qe,dset_qe2,x1lt_qe,time_qe2,]=...
        ReadCharge10June2016(fl_nm2);
%     [xg_qe,dset_qe3,x1lt_qe,time_qe3,]=...
%         ReadCharge10June2016(fl_nm3);
  
   
    [xg_b3,dset_b3,x1lt_b3,time_b3]=...
          ReadCharge10June2016(fl_nm4);
%     
   
    [xg_ene1,dset_ene1,x1lt_ene1,time_ene1]=...
        ReadCharge10June2016(fl_nm4);
    [xg_ene1,dset_ene2,x1lt_ene1,time_ene1]=...
        ReadCharge10June2016(fl_nm5);
%     [xg_ene1,dset_ene3,x1lt_ene1,time_ene1]=...
%         ReadCharge10June2016(fl_nm6);
%     
    [xg_ene_EMF,dset_ene_EMF,x1lt_ene_EMF,time_ene_EMF]=...
        ReadCharge10June2016(fl_nm7);
    
    [xg_ene_e,dset_ene_e,x1lt_ene_e,time_ene_e]=...
        ReadCharge10June2016(fl_nm8);
    
    [xg_ene_b,dset_ene_b,x1lt_ene_b,time_ene_b]=...
       ReadCharge10June2016(fl_nm9);
   
   
     [xg_b3,dset_b3,x1lt_b3,time_b3]=...
       ReadB310June2016(fl_nm10);
     [xg_b3,dset_e1,x1lt_b3,time_b3]=...
       ReadB310June2016(fl_nm11);
    
     [xg_b3,dset_e2,x1lt_b3,time_b3]=...
        ReadB310June2016(fl_nm12);
    
    
    dset_qe=dset_qe1+dset_qe2+dset_qe3;
    dset_ene=dset_ene1+dset_ene2+dset_ene3;
    dset_tot=dset_ene+dset_ene_EMF;
   
    dx_qe=abs(xg_qe(2)-xg_qe(1));
    
   
    E1E(cnt)=sum((dset_e1).*conj(dset_e1))*dx_qe;
    E2E(cnt)=sum((dset_e2).*conj(dset_e2))*dx_qe;
    B3E(cnt)=sum((dset_b3).*conj(dset_b3))*dx_qe;
    KinE(cnt)=sum((dset_ene))*dx_qe;
    MagE(cnt)=sum((dset_ene_b))*dx_qe;
    ElecE(cnt)=sum((dset_ene_e))*dx_qe;
    EMFE(cnt)=sum((dset_ene_EMF))*dx_qe;
    TotE(cnt)=sum(((1/2)*dset_ene_EMF))*dx_qe+sum((dset_ene))*dx_qe-0.1953;
    timE(cnt)=time_ene1;
   
    cnt=cnt+1;
  

end
% figure(1)
% plot(timE,KinE,'b','LineWidth',2);
% set(gca,'fontsize',20,'fontweight','b');
% xlabel('time(\omega_{pe})');
% ylabel('Energy');
% 
% hold on;
% plot(timE,MagE,'r','LineWidth',2);
% set(gca,'fontsize',20,'fontweight','b');
% xlabel('time(\omega_{pe})');
% ylabel('Energy');
% hold on;
% 
% plot(timE,ElecE,'g','LineWidth',2);
% set(gca,'fontsize',20,'fontweight','b');
% xlabel('time(\omega_{pe})');
% 
% hold on;
figure(1);
semilogy(timE,MagE,'m','LineWidth',2);
set(gca,'fontsize',20,'fontweight','b');
xlabel('time(\omega_{pe})');
figure(2);
semilogy(timE,E2E,'m','LineWidth',2);
set(gca,'fontsize',20,'fontweight','b');
xlabel('time(\omega_{pe})');

% hold off;
% 
% figure(2)
% 
% 
% plot(timE,MagE,'r','LineWidth',2);
% set(gca,'fontsize',20,'fontweight','b');
% xlabel('time(\omega_{pe})');
% ylabel('Energy');
% hold on;
% 
% plot(timE,ElecE,'g','LineWidth',2);
% set(gca,'fontsize',20,'fontweight','b');
% xlabel('time(\omega_{pe})');
% plot(timE,log(MagE),'g','LineWidth',2);hold on;
%      set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
%     title('Growth of Magnetic Enery');
%       xlabel('time(\omega_{pe})');
%     
%     ylabel('log|B_z|^2');




       


       
       
