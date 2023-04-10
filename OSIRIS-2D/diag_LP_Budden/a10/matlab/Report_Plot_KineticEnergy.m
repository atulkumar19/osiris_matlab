clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=1:500;

pth0='/Users/repulsion/Desktop/OSIRIS/OSIRIS_2D/atul/counter_finite_0.9c';
pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/Elec2/charge/';
pth3='/MS/DENSITY/Elec3/charge/';
pth4='/MS/DENSITY/Elec1/ene/';
pth5='/MS/DENSITY/Elec2/ene/';
pth6='/MS/DENSITY/Elec3/ene/';
pth7='/MS/FLD/ene_emf/';
pth8='/MS/FLD/ene_e/';
pth9='/MS/FLD/ene_b/';






str_h1='charge-Elec1-';
str_h2='charge-Elec2-';
str_h3='charge-Elec3-';

str_h4='ene-Elec1-';
str_h5='ene-Elec2-';
str_h6='ene-Elec3-';

str_h7='ene_emf-';
str_h8='ene_e-';
str_h9='ene_b-';

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
MagE=zeros(length(frm),1);
ElecE=zeros(length(frm),1);
EME=zeros(length(frm),1);
KinE=zeros(length(frm),1);
timE=zeros(length(frm),1);
TotE=zeros(length(frm),1);
KinE0=zeros(length(frm),1);
KinE_per=zeros(length(frm),1);
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
    
    
    
    
    [xg_qe,yg_qe,dset_qe1,x1lt_qe,x2lt_qe,time_qe]=...
        ReadECharge13Aug2014(fl_nm1);
    [xg_qe,yg_qe,dset_qe2,x1lt_qe,x2lt_qe,time_qe]=...
        ReadECharge13Aug2014(fl_nm2);
    [xg_qe,yg_qe,dset_qe3,x1lt_qe,x2lt_qe,time_qe]=...
        ReadECharge13Aug2014(fl_nm3);
   
    [xg_ene1,yg_ene1,dset_ene1,x1lt_ene1,x2lt_ene1,time_ene1]=...
        ReadECharge13Aug2014(fl_nm4);
    [xg_ene2,yg_ene2,dset_ene2,x1lt_ene2,x2lt_ene2,time_ene2]=...
        ReadECharge13Aug2014(fl_nm5);
    [xg_ene3,yg_ene3,dset_ene3,x1lt_ene3,x2lt_ene3,time_ene3]=...
        ReadECharge13Aug2014(fl_nm6);
    
    [xg_ene_EMF,yg_ene_EMF,dset_ene_EMF,x1lt_ene_EMF,x2lt_ene_EMF,time_ene_EMF]=...
        ReadB313Aug2014(fl_nm7);
    
    [xg_ene_e,yg_ene_e,dset_ene_e,x1lt_ene_e,x2lt_ene_e,time_ene_e]=...
        ReadE113Aug2014(fl_nm8);
    
    [xg_ene_b,yg_ene_b,dset_ene_b,x1lt_ene_b,x2lt_ene_b,time_ene_b]=...
        ReadB313Aug2014(fl_nm9);
    
    
    dset_qe=dset_qe1+dset_qe2+dset_qe3;
    dset_ene=dset_ene1+dset_ene2+dset_ene3;
    dset_tot=dset_ene+dset_ene_EMF;
   
   
    
    
    KinE(cnt)=sum(sum(dset_ene));
    MagE(cnt)=sum(sum(dset_ene_b));
    ElecE(cnt)=sum(sum(dset_ene_e));
    
    TotE(cnt)=sum(sum((1/2)*dset_ene_EMF+dset_ene));
    timE(cnt)=time_ene1;
   
    cnt=cnt+1;
   

end
plot(timE,KinE,'b','LineWidth',2);
set(gca,'fontsize',20,'fontweight','b');
xlabel('time(\omega_{pe})');
ylabel('Energy');
ylim([0 20000]);
% xlim([0 55]);
hold on;
plot(timE,MagE,'r','LineWidth',2);
set(gca,'fontsize',20,'fontweight','b');
xlabel('time(\omega_{pe})');
ylabel('Energy');
hold on;

plot(timE,ElecE,'g','LineWidth',2);
set(gca,'fontsize',20,'fontweight','b');
xlabel('time(\omega_{pe})');

hold on;
plot(timE,TotE,'m','LineWidth',2);
set(gca,'fontsize',20,'fontweight','b');
xlabel('time(\omega_{pe})');

hold off;

figure(2)

semilogy(timE,MagE,'r','LineWidth',2);
set(gca,'fontsize',20,'fontweight','b');
xlabel('time(\omega_{pe})');
ylabel('log|B|^2');



       


       
       

   
%  aviobj = close ( aviobj );
% close(writerObj);