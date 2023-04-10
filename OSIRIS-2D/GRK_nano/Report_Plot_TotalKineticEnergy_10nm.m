% clc;
% clear all;
% close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=0:220;



pth0='/home/testuser/GRK_plane';
pth1='/MS/DENSITY/Elec1/ene/';



pth37='/MS/FLD/ene_emf/';
pth38='/MS/FLD/ene_e/';
pth39='/MS/FLD/ene_b/';


str_h1='ene-Elec1-';



str_h37='ene_emf-';
str_h38='ene_e-';
str_h39='ene_b-';

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
ElecE=zeros(length(frm),1);
EME=zeros(length(frm),1);
KinE=zeros(length(frm),1);
timE=zeros(length(frm),1);
TotE=zeros(length(frm),1);
KinE0=zeros(length(frm),1);
KinE_per=zeros(length(frm),1);
Absorption=zeros(length(frm),1);
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
    
    
    fl_nm37=strcat(pth0,pth37,str_h37,str_num,str_ext);
    fl_nm38=strcat(pth0,pth38,str_h38,str_num,str_ext);
    fl_nm39=strcat(pth0,pth39,str_h39,str_num,str_ext);
    
    
    
%     [xg_qe,yg_qe,dset_qe1,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm1);
%     [xg_qe,yg_qe,dset_qe2,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm2);
%     [xg_qe,yg_qe,dset_qe3,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm3);
%     [xg_qe,yg_qe,dset_qe4,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm4);
%     [xg_qe,yg_qe,dset_qe5,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm5);
%      [xg_qe,yg_qe,dset_qe6,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm6);
%     [xg_qe,yg_qe,dset_qe7,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm7);
%     [xg_qe,yg_qe,dset_qe8,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm8);
%     [xg_qe,yg_qe,dset_qe9,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm9);
%     [xg_qe,yg_qe,dset_qe10,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm10);
%      [xg_qe,yg_qe,dset_qe11,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm11);
%     [xg_qe,yg_qe,dset_qe12,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm12);
%     [xg_qe,yg_qe,dset_qe13,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm13);
%     [xg_qe,yg_qe,dset_qe14,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm14);
%     [xg_qe,yg_qe,dset_qe15,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm15);
%     [xg_qe,yg_qe,dset_qe16,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm16);
%     [xg_qe,yg_qe,dset_qe17,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm17);
%     [xg_qe,yg_qe,dset_qe18,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm18);
    
   
    [xg_ene1,yg_ene1,dset_ene1,x1lt_ene1,x2lt_ene1,time_ene1]=...
        ReadECharge13Aug2014(fl_nm1);
    
    
    
    [xg_ene_EMF,yg_ene_EMF,dset_ene_EMF,x1lt_ene_EMF,x2lt_ene_EMF,time_ene_EMF]=...
        ReadB313Aug2014(fl_nm37);
    
    [xg_ene_e,yg_ene_e,dset_ene_e,x1lt_ene_e,x2lt_ene_e,time_ene_e]=...
        ReadE113Aug2014(fl_nm38);
    
    [xg_ene_b,yg_ene_b,dset_ene_b,x1lt_ene_b,x2lt_ene_b,time_ene_b]=...
        ReadB313Aug2014(fl_nm39);
    
    
    
    dset_ene=dset_ene1;
    dset_tot=dset_ene+dset_ene_EMF;
   
    dx_ene=abs(xg_ene1(2)-xg_ene1(1));
    
   
    
    KinE(cnt)=sum(sum(dset_ene))*dx_ene.^2;
%     MagE(cnt)=sum(sum(dset_ene_b))*dx_qe.^2;
%     ElecE(cnt)=sum(sum(dset_ene_e))*dx_qe.^2;
%     
%     TotE(cnt)=sum(sum((1/2)*dset_ene_EMF+dset_ene))*dx_qe.^2;
    timE(cnt)=time_ene1;
   
    cnt=cnt+1;


end
% 
% for i=1:length(frm)-1
% 
% Absorption(i)=KinE(i+1)-KinE(i);
% 
% end

plot(timE, KinE,'*-g','LineWidth',2);
hold on;
set(gca,'fontsize',20,'fontweight','b');
xlabel('time(\omega_{pe})');
ylabel('KinEnergy');